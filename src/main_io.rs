mod kmer;

use crayfish::collective;
use crayfish::finish;
use crayfish::place;
use crayfish::place::world_size;
use crayfish::logging::*;
use crayfish::place::Place;
use crayfish::shared::PlaceLocal;
use crayfish::shared::PlaceLocalWeak;
use std::fs::File;
use std::sync::Mutex;

use kmer::AbstractKMer;
use kmer::KMeru64;
use kmer::DNA;

const KMER_LEN: usize = 31;

type CountBin = Vec<u64>;
type KMer = KMeru64<DNA, KMER_LEN>;

#[crayfish::activity]
async fn update_kmer(kmers: Vec<u64>, final_ptr: PlaceLocalWeak<Mutex<CountBin>>) {
    let ptr = final_ptr.upgrade().unwrap();
    let mut h = ptr.lock().unwrap();
    h.extend_from_slice(&kmers[..]);
}

fn get_partition(kmer: &KMer) -> usize {
    let mut key = kmer.data;
    key = !key + (key << 21);
	key = key ^ key >> 24;
	key = (key + (key << 3)) + (key << 8);
	key = key ^ key >> 14;
	key = (key + (key << 2)) + (key << 4);
	key = key ^ key >> 28;
	key = key + (key << 31);
    (key % place::world_size() as u64) as usize
    // let key = kmer.data;
    // let mut hs = rustc_hash::FxHasher::default();
    // key.hash(&mut hs);
    // (hs.finish() % place::world_size() as u64) as usize
}

// TODO: stupid fasta/fastq reader
struct SeqReader<I>
{
    lines: I,
}

impl<T> SeqReader<T>
{
    pub fn new(lines: T) -> Self {
        SeqReader { lines }
    }
}

impl<'a, T> Iterator for SeqReader<T>
where
    T: Iterator<Item = String> + 'a,
{
    type Item = Vec<u8>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut skip_q = false;
        loop {
            let line = self.lines.next()?;
            if !line.is_empty() {
                let line = line.into_bytes();
                match line[0] {
                    b'@' => {
                        skip_q = false;
                    }
                    b'>' => (),
                    b'+' => {
                        skip_q = true;
                    }
                    _ => {
                        if !skip_q {
                            return Some(line);
                        }
                    }
                }
            }
        }
    }
}

struct Lines<'a>{
    data: &'a[u8]
}
impl<'a> Lines<'a>{
    fn new(data:&'a [u8]) -> Self{
        Lines{
            data
        }
    }
}

impl<'a> Iterator for Lines<'a>{
    type Item = &'a[u8];
    fn next(&mut self) -> Option<Self::Item> {
        let pos = memchr::memchr(b'\n', self.data)?;
        let ret = &self.data[..pos];
        if pos == self.data.len() - 1{ // is last
            self.data = &self.data[pos..pos]; // empty
        }else{
            self.data = &self.data[pos + 1..];
        }
        Some(ret)
    }
}

// desugered finish
#[crayfish::main]
async fn inner_main() {
    let count_bin = PlaceLocal::new(Mutex::new(CountBin::default()));
    collective::barrier().await;
    // ctx contains a new finish id now
        let mut kmers = vec![vec![]; place::world_size()];
        let args = std::env::args().collect::<Vec<_>>();
        let filename = &args[1];
        let file = File::open(filename).unwrap();
        use std::io::BufReader;
        use std::io::BufRead;
        let bf = BufReader::new(file);
        let lines = SeqReader::new(bf.lines().map(|lr|lr.unwrap()));
        // let lines = Lines

        let world_size = world_size();
        let here = place::here();

        let chunk_size = 40960usize;

        finish! {
        for (l_num, read) in lines.enumerate() {
            if here as usize != l_num % world_size as usize{
                continue
            }
            if read.len() < KMer::kmer_len() {
                continue;
            }

            let mut next_pos = 0;
            let mut start = true;
            let end = read.len();
            let mut current_kmer = KMer::new(0); // fake start, won't be extended
            while next_pos < end {
                if start {
                    match KMer::from_bytes(&read[next_pos..]) {
                        Some(k) => {
                            let k = k.get_canonical();
                            // TODO should depends on trait. struct field k.data used here
                            kmers[get_partition(&k)].push(k.data);
                            current_kmer = k;
                            next_pos += KMer::kmer_len();
                            start = false;
                        }
                        None => {
                            next_pos += 1;
                        }
                    }
                } else {
                    match current_kmer.extend(read[next_pos]) {
                        Some(k) => {
                            let k = k.get_canonical();
                            kmers[get_partition(&k)].push(k.data);
                            current_kmer = k;
                        }
                        None => {
                            start = true;
                        }
                    }
                    next_pos += 1;
                }
            }

            // interleave communication and computing
            if l_num / world_size  % chunk_size == 0 {
                let mut new_kmers = vec![vec![]; place::world_size()];
                std::mem::swap(&mut new_kmers, &mut kmers);
                for (dst, kmer_list) in new_kmers.into_iter().enumerate() {
                    crayfish::ff!(dst as Place, update_kmer(kmer_list, count_bin.downgrade()));
                }
            }
        }

        for (dst, kmer_list) in kmers.into_iter().enumerate() {
            crayfish::ff!(dst as Place, update_kmer(kmer_list, count_bin.downgrade()));
        }
        info!("k-mer gen done");
            
        }
    collective::barrier().await;
    info!("start counting");

    let mut sorted_bin = count_bin.lock().unwrap();
    use crate::kmer::radix::RadixSort;
    sorted_bin.voracious_sort();

    let mut hist = vec![0usize; 1024];
    let hist_len = hist.len();

    if sorted_bin.is_empty() {
        return
    }
    let mut current = sorted_bin[0];
    let mut count = 0;

    for kmer in sorted_bin.iter() {
        if current != *kmer {
            if count <= hist_len {
                hist[count - 1] += 1;
            }
            current = *kmer;
            count = 1;
        }else{
            count += 1;
        }
    }
    info!("{:?}", hist);
}
