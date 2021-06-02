mod kmer;

use crayfish::collective;
use crayfish::finish;
use crayfish::global_id;
use crayfish::global_id::world_size;
use crayfish::logging::*;
use crayfish::place::Place;
use crayfish::shared::PlaceLocal;
use crayfish::shared::PlaceLocalWeak;
use rustc_hash::FxHashMap;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::Hash;
use std::hash::Hasher;
use std::io;
use std::io::BufRead;
use std::io::BufReader;
use std::sync::Mutex;

use kmer::AbstractKMer;
use kmer::KMeru64;
use kmer::DNA;

const KMER_LEN: usize = 31;

type Reads = Vec<Vec<u8>>;
type CountTable = FxHashMap<KMer, usize>;
type KMer = KMeru64<DNA, KMER_LEN>;

#[crayfish::activity]
async fn update_kmer(kmers: Vec<u64>, final_ptr: PlaceLocalWeak<Mutex<CountTable>>) {
    // info!("Got {} kmers. Counting", kmers.len());
    let ptr = final_ptr.upgrade().unwrap();
    let mut h = ptr.lock().unwrap();
    for kmer in kmers {
        *h.entry(KMer::new(kmer)).or_insert(0) += 1;
    }
}

fn get_partition(kmer: &KMer) -> usize {
    let mut hasher = DefaultHasher::new();
    kmer.hash(&mut hasher);
    (hasher.finish() % global_id::world_size() as u64) as usize
}

#[crayfish::activity]
async fn kmer_counting(reads: Reads, final_ptr: PlaceLocalWeak<Mutex<CountTable>>) {
    info!("Got {} reads. Spliting into Kmers", reads.len());

    let mut kmers = vec![vec![]; global_id::world_size()];
    for read in reads {
        // drop too short read
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
    }

    info!("Sending kmers to destination");
    for (dst, kmer_list) in kmers.into_iter().enumerate() {
        crayfish::ff!(dst as Place, update_kmer(kmer_list, final_ptr.clone()));
    }
}

// TODO: stupid fasta reader
struct SeqReader<I>
where
    I: Iterator<Item = io::Result<String>>,
{
    lines: I,
}

impl<T> SeqReader<T>
where
    T: Iterator<Item = io::Result<String>>,
{
    pub fn new(lines: T) -> Self {
        SeqReader { lines }
    }
}

impl<T> Iterator for SeqReader<T>
where
    T: Iterator<Item = io::Result<String>>,
{
    type Item = Vec<u8>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut skip_q = false;
        loop {
            let line = self.lines.next()?.unwrap().into_bytes();
            if !line.is_empty() {
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

// desugered finish
#[crayfish::main]
async fn inner_main() {
    let count_table_ptr = PlaceLocal::new(Mutex::new(CountTable::default()));
    collective::barrier().await;
    if global_id::here() == 0 {
        // ctx contains a new finish id now
        let chunk_size = 81920;
        let args = std::env::args().collect::<Vec<_>>();
        let filename = &args[1];
        let file = File::open(filename).unwrap();
        let lines = BufReader::new(file).lines();
        let lines = SeqReader::new(lines.into_iter());

        let world_size = world_size();
        let mut next_place: Place = 0;
        let mut buffer: Reads = vec![];

        finish! {
        for (l_num, line) in lines.enumerate() {
                if buffer.len() == chunk_size {
                    info!(
                        "Sending {}~{} reads to {}",
                        l_num,
                        l_num + chunk_size - 1,
                        next_place
                    );
                    let mut new_read = vec![];
                    std::mem::swap(&mut new_read, &mut buffer);
                    crayfish::ff!(next_place, kmer_counting(new_read, count_table_ptr.downgrade()));
                    next_place = (next_place + 1) % (world_size as Place);
                }
                buffer.push(line);
        }
        }
    }
    collective::barrier().await;

    let global_table = count_table_ptr.lock().unwrap();

    let mut hist = vec![0usize; 1024];
    for (_, count) in global_table.iter() {
        if *count <= hist.len() {
            hist[*count - 1] += 1;
        }
    }
    println!("{:?}", hist);
}
