use std::cmp::Ordering;
use std::fmt;
use std::fmt::Debug;
use std::fmt::Formatter;
use std::hash::Hash;
use std::hash::Hasher;
use std::marker::PhantomData;
use std::mem::size_of;

pub extern crate voracious_radix_sort as radix;

pub trait Alphabet {
    // size_of_usize % unit_len must be 0
    const UNIT_LEN: usize;
    fn to_unit(base: u8) -> Option<u8>;
    fn to_u8(unit: u8) -> Option<u8>;
}

pub struct DNA {}

impl Alphabet for DNA {
    const UNIT_LEN: usize = 2;
    fn to_unit(base: u8) -> Option<u8> {
        Some(match base {
            b'A' => 0b00,
            b'T' => 0b11,
            b'G' => 0b10,
            b'C' => 0b01,
            _ => return None,
        })
    }
    fn to_u8(unit: u8) -> Option<u8> {
        Some(match unit {
            0b00 => b'A',
            0b11 => b'T',
            0b10 => b'G',
            0b01 => b'C',
            _ => return None,
        })
    }
}

pub trait AbstractKMer
where
    Self: Sized + Ord + Copy + Hash,
{
    fn kmer_len() -> usize;
    fn extend(&self, base: u8) -> Option<Self>;
    fn from_bytes(bytes: &[u8]) -> Option<Self>;
    fn reverse(&self) -> Self;
    fn complement(&self) -> Self;
    fn get_canonical(&self) -> Self {
        let rc = self.reverse().complement();
        if self > &rc {
            rc
        } else {
            *self
        }
    }
}

pub struct KMeru64<A, const KMERLEN: usize>
where
    A: Alphabet,
{
    pub data: u64,
    _mark: PhantomData<A>,
}

impl<A, const N:usize> radix::Radixable<u64> for KMeru64<A, N> where A: Alphabet{
    type Key = u64;
    fn key(&self) -> Self::Key{
        self.data
    }
}
impl<A, const N: usize> Copy for KMeru64<A, N> where A: Alphabet {}
impl<A, const N: usize> Clone for KMeru64<A, N>
where
    A: Alphabet,
{
    fn clone(&self) -> Self {
        Self::new(self.data)
    }
}

impl<A, const N: usize> Hash for KMeru64<A, N>
where
    A: Alphabet,
{
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.data.hash(state);
    }
}

impl<A, const N: usize> Eq for KMeru64<A, N> where A: Alphabet {}

impl<A, const N: usize> PartialEq for KMeru64<A, N>
where
    A: Alphabet,
{
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl<A, const N: usize> Ord for KMeru64<A, N>
where
    A: Alphabet,
{
    fn cmp(&self, other: &Self) -> Ordering {
        self.data.cmp(&other.data)
    }
}

impl<A, const N: usize> PartialOrd for KMeru64<A, N>
where
    A: Alphabet,
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.data.partial_cmp(&other.data)
    }
}

impl<A, const N: usize> Debug for KMeru64<A, N>
where
    A: Alphabet,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        f.write_str(&self.to_string())
    }
}

impl<A, const N: usize> Default for KMeru64<A, N>
where
    A: Alphabet,
{
    fn default() -> Self {
        Self::new(0)
    }
}

impl<A, const N: usize> KMeru64<A, N>
where
    A: Alphabet,
{
    pub fn new(data: u64) -> Self {
        KMeru64 {
            data,
            _mark: PhantomData,
        }
    }
    fn set_unit(&mut self, at: usize, unit: u8) {
        // this won't rewrite if unit is set
        let unit = unit as u64;
        debug_assert!(at < Self::unit_num());
        debug_assert!(unit < 1 << A::UNIT_LEN);
        self.data |= unit << (Self::data_len() - A::UNIT_LEN * (at + 1))
    }

    fn get_unit(&self, at: usize) -> u8 {
        debug_assert!(at < size_of::<u64>() * 8 / A::UNIT_LEN);
        let mut ret = self.data >> (Self::data_len() - A::UNIT_LEN * (at + 1));
        ret &= (1 << A::UNIT_LEN) - 1;
        ret as u8
    }

    fn unused_bits() -> usize {
        size_of::<u64>() * 8 - Self::used_bits()
    }

    fn used_bits() -> usize {
        A::UNIT_LEN * Self::kmer_len()
    }

    fn data_len() -> usize {
        size_of::<u64>() * 8
    }

    fn unit_num() -> usize {
        Self::data_len() / A::UNIT_LEN
    }
}

impl<A, const N: usize> AbstractKMer for KMeru64<A, N>
where
    A: Alphabet,
{
    fn extend(&self, base: u8) -> Option<Self> {
        let mut next = *self;
        let unit = A::to_unit(base)? as u64;
        next.data = self.data << A::UNIT_LEN | unit;
        next.data &= (u64::MAX) >> Self::unused_bits();
        Some(next)
    }
    fn kmer_len() -> usize {
        N
    }
    fn complement(&self) -> Self {
        let data = !self.data & u64::MAX >> Self::unused_bits();
        Self::new(data)
    }
    fn reverse(&self) -> Self {
        #[cfg(not(target_feature = "sse"))]
        compile_error!("should support ssse3");
        #[cfg(target_feature = "sse")]
        let data = ssse3::reverse_u64_pack_2(self.data) >> Self::unused_bits();
        Self::new(data)
    }

    fn from_bytes(s: &[u8]) -> Option<Self> {
        let mut kmer = Self::default();

        if s.len() < Self::kmer_len() {
            return None;
        }
        let mut iter = s.iter().take(Self::kmer_len()).enumerate();
        while let Some((i, c)) = iter.next() {
            kmer.set_unit(i, A::to_unit(*c)?)
        }
        kmer.data >>= Self::unused_bits();
        Some(kmer)
    }
}

impl<A, const N: usize> std::str::FromStr for KMeru64<A, N>
where
    A: Alphabet,
{
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_bytes(s.as_bytes()).ok_or(())
    }
}

impl<A, const N: usize> std::string::ToString for KMeru64<A, N>
where
    A: Alphabet,
{
    fn to_string(&self) -> String {
        let mut s = vec![];
        for i in Self::unused_bits() / A::UNIT_LEN..size_of::<u64>() * 8 / A::UNIT_LEN {
            s.push(A::to_u8(self.get_unit(i)).unwrap())
        }
        String::from_utf8(s).unwrap()
    }
}

#[cfg(target_feature = "sse")]
mod ssse3 {
    #[cfg(target_arch = "x86")]
    use std::arch::x86::*;
    #[cfg(target_arch = "x86_64")]
    use std::arch::x86_64::*;

    unsafe fn reverse_m128i_pack_2(mut v: __m128i) -> __m128i {
        // from https://github.com/ParBLiSS/kmerind/blob/0062fe91fdeef66fce4d1e897c15318241130277/src/common/test/kmer_reverse_helper.hpp#L269
        // reverse byte

        // println!("hapi {:X?}", v);
        v = _mm_shuffle_epi8(
            v,
            _mm_set_epi8(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15),
        );
        // println!("hapi {:X?}", v);
        let lo_mask = _mm_set_epi8(
            0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F, 0x0F,
            0x0F, 0x0F,
        );
        // lower part of each byte
        let lo = _mm_and_si128(lo_mask, v);
        let hi = _mm_srli_epi16(_mm_andnot_si128(lo_mask, v), 4);

        // println!("lo {:X?}", lo);
        #[allow(overflowing_literals)]
        let slo = _mm_shuffle_epi8(
            _mm_setr_epi8(
                0x00, 0x40, 0x80, 0xc0, 0x10, 0x50, 0x90, 0xd0, 0x20, 0x60, 0xa0, 0xe0, 0x30, 0x70,
                0xb0, 0xf0,
            ),
            lo,
        );
        // println!("slo {:X?}", slo);
        // println!("ho {:X?}", hi);
        let shi = _mm_shuffle_epi8(
            _mm_setr_epi8(
                0x00, 0x04, 0x08, 0x0c, 0x01, 0x05, 0x09, 0x0d, 0x02, 0x06, 0x0a, 0x0e, 0x03, 0x07,
                0x0b, 0x0f,
            ),
            hi,
        );
        // println!("shi {:X?}", shi);
        _mm_or_si128(slo, shi)
    }

    pub fn reverse_u64_pack_2(data: u64) -> u64 {
        unsafe {
            let v = _mm_set_epi64x(data as i64, data as i64);
            _mm_extract_epi64(reverse_m128i_pack_2(v), 1) as u64
        }
    }

    #[cfg(test)]
    mod test {
        use super::*;
        #[test]
        pub fn test_reverse_u64_pack_2() {
            let data: u64 = 0x0123456789ABCDEF;
            assert_eq!(reverse_u64_pack_2(reverse_u64_pack_2(data)), data)
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    type KMer31 = KMeru64<DNA, 31>;
    #[test]
    pub fn test_parse() {
        let read = "TCGCGTAGCTAGCATATATTCGCGGCTAGTAC";
        let kmer = read.parse::<KMer31>().unwrap();
        assert_eq!(&kmer.to_string(), &read[..KMer31::kmer_len()]);
    }

    #[test]
    pub fn test_get_unit() {
        let read = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mut kmer = read.parse::<KMer31>().unwrap();
        kmer.set_unit(1, DNA::to_unit(b'C').unwrap());
        assert_eq!(DNA::to_u8(kmer.get_unit(1)).unwrap(), b'C');
    }

    #[test]
    pub fn test_complement() {
        let kmer = "AAAAAAAAAAAAAAAAAAAAAAAAAAACCCC".parse::<KMer31>().unwrap();
        let comp = "TTTTTTTTTTTTTTTTTTTTTTTTTTTGGGG".parse::<KMer31>().unwrap();
        assert_eq!(kmer.complement(), comp);
    }

    #[test]
    pub fn test_reverse() {
        let kmer = "AAAAAAAAAAAAAAAAAAAAAAAAAAACCCC".parse::<KMer31>().unwrap();
        let reverse = "CCCCAAAAAAAAAAAAAAAAAAAAAAAAAAA".parse::<KMer31>().unwrap();
        assert_eq!(kmer.reverse(), reverse);
    }

    #[test]
    pub fn test_canonical() {
        let kmer = "CCCCAAAAAAAAAAAAAAAAAAAAAAAAAAA".parse::<KMer31>().unwrap();
        assert_eq!(kmer, kmer.get_canonical());

        let kmer = "GGGGGGGGGGGGGGGGGGGGGGGGGGGCCCC".parse::<KMer31>().unwrap();
        let c = "GGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
            .parse::<KMer31>()
            .unwrap();
        assert_eq!(c, kmer.get_canonical());
    }

    #[test]
    pub fn test_extend() {
        let kmer = "CCCCAAAAAAAAAAAAAAAAAAAAAAAAAAA".parse::<KMer31>().unwrap();
        let kmer_e = "CCCAAAAAAAAAAAAAAAAAAAAAAAAAAAT".parse::<KMer31>().unwrap();
        assert_eq!(kmer.extend(b'T').unwrap(), kmer_e);
    }
}
