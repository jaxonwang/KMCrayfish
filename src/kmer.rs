use std::fmt;
use std::fmt::Debug;
use std::fmt::Formatter;
use std::marker::PhantomData;
use std::mem::size_of;

trait Alphabet: Debug + Copy + Ord {
    // size_of_usize % unit_len must be 0
    const UNIT_LEN: usize;
    fn to_unit(base: char) -> Option<u8>;
    fn to_char(unit: u8) -> Option<char>;
}

#[derive(Debug, Copy, Clone, Ord, PartialOrd, PartialEq, Eq)]
struct DNA {}

impl Alphabet for DNA {
    const UNIT_LEN: usize = 2;
    fn to_unit(base: char) -> Option<u8> {
        Some(match base {
            'A' => 0b00,
            'T' => 0b11,
            'G' => 0b10,
            'C' => 0b01,
            _ => return None,
        })
    }
    fn to_char(unit: u8) -> Option<char> {
        Some(match unit {
            0b00 => 'A',
            0b11 => 'T',
            0b10 => 'G',
            0b01 => 'C',
            _ => return None,
        })
    }
}

const KMER_LEN: usize = 31;
//
trait KMer
where
    Self: Sized + Ord + Copy,
{
    fn kmer_len() -> usize;
    fn extend(&self, base: char) -> Option<Self>;
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

#[derive(Copy, Clone, Ord, PartialOrd)]
struct KMer64<A>
where
    A: Alphabet,
{
    data: u64,
    _mark: PhantomData<A>,
}

impl<A> Eq for KMer64<A> where A: Alphabet {}

impl<A> PartialEq for KMer64<A>
where
    A: Alphabet,
{
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl<A> Debug for KMer64<A>
where
    A: Alphabet,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> Result<(), fmt::Error> {
        f.write_str(&self.to_string())
    }
}

impl<A> Default for KMer64<A>
where
    A: Alphabet,
{
    fn default() -> Self {
        Self::new(0)
    }
}

impl<A> KMer64<A>
where
    A: Alphabet,
{
    fn new(data: u64) -> Self {
        KMer64 {
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
        A::UNIT_LEN * KMER_LEN
    }

    fn data_len() -> usize {
        size_of::<u64>() * 8
    }

    fn unit_num() -> usize {
        Self::data_len() / A::UNIT_LEN
    }
}

impl<A> KMer for KMer64<A>
where
    A: Alphabet,
{
    fn extend(&self, base: char) -> Option<Self> {
        let mut next = *self;
        let unit = A::to_unit(base)? as u64;
        next.data = self.data << A::UNIT_LEN | unit;
        next.data &= (u64::MAX) >> Self::unused_bits();
        Some(next)
    }
    fn kmer_len() -> usize {
        KMER_LEN
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
}

impl<A> std::str::FromStr for KMer64<A>
where
    A: Alphabet,
{
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut kmer = Self::default();
        debug_assert!(s.len() >= KMER_LEN);
        let mut iter = s.chars().take(KMER_LEN).enumerate();
        while let Some((i, c)) = iter.next() {
            kmer.set_unit(i, A::to_unit(c).ok_or(())?)
        }
        kmer.data >>= Self::unused_bits();
        Ok(kmer)
    }
}

impl<A> std::string::ToString for KMer64<A>
where
    A: Alphabet,
{
    fn to_string(&self) -> String {
        let mut s = String::default();
        for i in Self::unused_bits() / A::UNIT_LEN..size_of::<u64>() * 8 / A::UNIT_LEN {
            s.push(A::to_char(self.get_unit(i)).unwrap())
        }
        s
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
                // 0x00, 0x80, 0x40, 0xc0, 0x20, 0xa0, 0x60, 0xe0, 0x10, 0x90, 0x50, 0xd0, 0x30, 0xb0,
                // 0x70, 0xf0,
                0x00, 0x40, 0x80, 0xc0, 0x10, 0x50, 0x90, 0xd0, 0x20, 0x60, 0xa0, 0xe0, 0x30, 0x70,
                0xb0, 0xf0,
            ),
            lo,
        );
        // println!("slo {:X?}", slo);
        // println!("ho {:X?}", hi);
        let shi = _mm_shuffle_epi8(
            _mm_setr_epi8(
                // 0x00, 0x08, 0x04, 0x0c, 0x02, 0x0a, 0x06, 0x0e, 0x01, 0x09, 0x05, 0x0d, 0x03, 0x0b,
                // 0x07, 0x0f,
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
    #[test]
    pub fn test_parse() {
        let read = "TCGCGTAGCTAGCATATATTCGCGGCTAGTAC";
        let kmer = read.parse::<KMer64<DNA>>().unwrap();
        assert_eq!(&kmer.to_string(), &read[..KMER_LEN]);
    }

    #[test]
    pub fn test_get_unit() {
        let read = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mut kmer = read.parse::<KMer64<DNA>>().unwrap();
        kmer.set_unit(1, DNA::to_unit('C').unwrap());
        assert_eq!(DNA::to_char(kmer.get_unit(1)).unwrap(), 'C');
    }

    #[test]
    pub fn test_complement() {
        let kmer = "AAAAAAAAAAAAAAAAAAAAAAAAAAACCCC"
            .parse::<KMer64<DNA>>()
            .unwrap();
        let comp = "TTTTTTTTTTTTTTTTTTTTTTTTTTTGGGG"
            .parse::<KMer64<DNA>>()
            .unwrap();
        assert_eq!(kmer.complement(), comp);
    }

    #[test]
    pub fn test_reverse() {
        let kmer = "AAAAAAAAAAAAAAAAAAAAAAAAAAACCCC"
            .parse::<KMer64<DNA>>()
            .unwrap();
        let reverse = "CCCCAAAAAAAAAAAAAAAAAAAAAAAAAAA"
            .parse::<KMer64<DNA>>()
            .unwrap();
        assert_eq!(kmer.reverse(), reverse);
    }

    #[test]
    pub fn test_canonical() {
        let kmer = "CCCCAAAAAAAAAAAAAAAAAAAAAAAAAAA"
            .parse::<KMer64<DNA>>()
            .unwrap();
        assert_eq!(kmer, kmer.get_canonical());

        let kmer = "GGGGGGGGGGGGGGGGGGGGGGGGGGGCCCC"
            .parse::<KMer64<DNA>>()
            .unwrap();
        let c = "GGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCC"
            .parse::<KMer64<DNA>>()
            .unwrap();
        assert_eq!(c, kmer.get_canonical());
    }
}
