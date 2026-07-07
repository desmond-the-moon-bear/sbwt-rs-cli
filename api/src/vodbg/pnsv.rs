// Code by Martin Kostadinov.

pub mod balanced_parenthesis;
pub mod matrix;
pub mod ranges;
pub mod scan;
pub mod wavelet;

use balanced_parenthesis as bp;

use crate::ContractLeft;
use crate::ExtendRight;
use crate::LcsArray;

pub use balanced_parenthesis::LcsPnsvBp;
pub use balanced_parenthesis::PnsvBp;
pub use matrix::Matrix as PnsvMatrix;
pub use ranges::Ranges;
pub use scan::AugmentedBoundedScan as ABS;
pub use scan::LcsSimd;
pub use wavelet::WindowedWaveletTree as WWT;

use simple_sds_sbwt::serialize::Serialize;

/// Previous/Next Smaller value.
pub trait Pnsv {
    fn previous(&self, index: usize, target_length: usize) -> usize;
    fn next(&self, index: usize, target_length: usize) -> usize;
    fn override_contract_left(&self) -> bool { false }
    #[allow(unused_variables)]
    fn overriden_contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> { 0..0 }
    fn max_target(&self) -> usize { 0 }
}

impl<T: ?Sized + Pnsv> ContractLeft for T {
    fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {
        if self.override_contract_left() {
            // note(mk): Godbolt says that if the method is "overriden" by changing the default
            // implementation of Pnsv::override_contract_left the compiler will optimise this
            // condition away and just include the code of the overriden method.
            return self.overriden_contract_left(I, target_len);
        }
        let new_start = self.previous(I.start, target_len);
        let new_end = self.next(I.end, target_len);
        new_start..new_end
    }
}

// note(mk): Probably better to implement this with enum_dispatch.
pub struct PnsvDynOwned {
    pub structures: Vec<Box<dyn Pnsv>>,
}

impl Pnsv for PnsvDynOwned {
    fn previous(&self, index: usize, target_length: usize) -> usize {
        for i in 0..self.structures.len() - 1 {
            if target_length <= self.structures[i].max_target() {
                return self.structures[i].previous(index, target_length);
            }
        }
        self.structures[self.structures.len() - 1].previous(index, target_length)
    }

    fn next(&self, index: usize, target_length: usize) -> usize {
        for i in 0..self.structures.len() - 1 {
            if target_length <= self.structures[i].max_target() {
                return self.structures[i].next(index, target_length);
            }
        }
        self.structures[self.structures.len() - 1].next(index, target_length)
    }

    fn override_contract_left(&self) -> bool {
        true
    }

    #[allow(non_snake_case)]
    fn overriden_contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {
        for i in 0..self.structures.len() - 1 {
            if target_len <= self.structures[i].max_target() {
                return self.structures[i].contract_left(I, target_len);
            }
        }
        self.structures[self.structures.len() - 1].contract_left(I, target_len)
    }
}

// Experimentally the scan is fastest if the average length of the ranges it searches is below
// around 200 i.e. the target length. This value is equal to approx log_4(200).
const TARGET_LENGTH_LOG_4_FLOOR: usize = 3;

pub fn pnsv_abs_simd(extend: &impl ExtendRight, lcs: &LcsArray) -> PnsvDynOwned {
    let count = lcs.len();

    log::info!("[pnsv_abs_simd] creating ranges...");
    let mut ranges_upper_bound = 0;
    let mut bits_in_current_level_of_ranges = usize::BITS as usize * 4;
    while bits_in_current_level_of_ranges < count {
        ranges_upper_bound += 1;
        bits_in_current_level_of_ranges *= 4;
    }
    ranges_upper_bound = ranges_upper_bound.min(Ranges::MAX_K);
    let ranges = Ranges::new(extend, count, ranges_upper_bound);
    let ranges_box = Box::new(ranges);

    let iterator = (0..count).map(|index| lcs.access(index) as u8);

    let log_4 = (usize::BITS - count.leading_zeros()).div_ceil(2) as usize;

    // log_4(count / 200) == log_4(count) - log_4(200)
    let matrix_upper_bound = log_4 - TARGET_LENGTH_LOG_4_FLOOR; 
    
    log::info!("[pnsv_abs_simd] creating lcs simd...");
    let lcs_simd = LcsSimd::from_iterator(iterator.clone(), count);

    log::info!("[pnsv_abs_simd] creating augmented bounded scan...");
    let abs = ABS::from_iterator(lcs_simd, iterator, 8, ranges_upper_bound + 1, matrix_upper_bound);
    let abs_box = Box::new(abs);

    log::info!("[pnsv_abs_simd] target length ranges: 1:{}:{}:..", ranges_upper_bound, matrix_upper_bound);

    PnsvDynOwned {
        structures: vec![ranges_box, abs_box],
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PnsvTuned {
    pub ranges: Ranges,
    pub matrix: PnsvMatrix,
    pub lcs_simd: LcsSimd,
    pub scan_bound: usize,
    pub fallback_scan_overlap: usize,
}

impl PnsvTuned {
    pub const DEFAULT_SCAN_BOUND: usize = 16;
    pub const DEFAULT_FALLBACK_OVERLAP: usize = 2;

    pub fn new_with_default_values(extend: &impl ExtendRight, lcs: &LcsArray, max_k: usize) -> Self {
        Self::new(extend, lcs, max_k, Self::DEFAULT_SCAN_BOUND, Self::DEFAULT_FALLBACK_OVERLAP)
    }

    pub fn new(extend: &impl ExtendRight, lcs: &LcsArray, max_k: usize, scan_bound: usize, mut fallback_scan_overlap: usize) -> Self {
        let count = lcs.len();

        log::info!("[PnsvTuned::new] creating ranges...");
        let ranges = make_ranges(extend, count, max_k);
        let ranges_upper_bound = ranges.max_target();

        let log_4 = (usize::BITS - count.leading_zeros()).div_ceil(2) as usize;
        let mut matrix_upper_bound = log_4.saturating_sub(TARGET_LENGTH_LOG_4_FLOOR);
        matrix_upper_bound = matrix_upper_bound.min(ranges_upper_bound + 1 + matrix::MAX_ROWS);
        matrix_upper_bound = matrix_upper_bound.min(max_k);
        
        let mut matrix_option = None;
        let mut lcs_simd_option = None;

        let matrix_iterator = (0..count).map(|index| lcs.access(index) as u8);
        let lcs_simd_iterator = (0..count).map(|index| lcs.access(index) as u8);

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(2).build().unwrap();
        thread_pool.scope(|s| {
            if matrix_upper_bound > ranges_upper_bound {
                s.spawn(|_| {
                    log::info!("[PnsvTuned::new] creating matrix...");
                    let matrix = PnsvMatrix::from_iterator(matrix_iterator, count, ranges_upper_bound + 1, matrix_upper_bound);
                    fallback_scan_overlap = fallback_scan_overlap.min(matrix.max_target());
                    matrix_option = Some(matrix);
                });
            }
            s.spawn(|_| {
                log::info!("[PnsvTuned::new] creating lcs simd...");
                let lcs_simd = LcsSimd::from_iterator(lcs_simd_iterator, count);
                lcs_simd_option = Some(lcs_simd);
            });
        });

        let matrix = match matrix_option {
            Some(value) => value,
            None => PnsvMatrix::empty(),
        };
        let lcs_simd = lcs_simd_option.expect("Creating lcs simd should not fail.");

        log::info!("[PnsvTuned::new] target length ranges: 1:{}:{}:..", ranges_upper_bound, matrix_upper_bound);

        Self {
            ranges,
            matrix,
            lcs_simd,
            scan_bound,
            fallback_scan_overlap,
        }
    }

    pub fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut written = 0;
        out.write_all(&(self.scan_bound as u64).to_le_bytes())?;
        out.write_all(&(self.fallback_scan_overlap as u64).to_le_bytes())?;
        written += 2 * size_of::<u64>();
        written += self.ranges.serialize(out)?;
        written += self.matrix.serialize(out)?;
        written += self.lcs_simd.serialize(out)?;
        Ok(written)
    }

    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
        let scan_bound = u64::from_le(u64::load(input)?) as usize;
        let fallback_scan_overlap = u64::from_le(u64::load(input)?) as usize;
        let ranges = Ranges::load(input)?;
        let matrix = PnsvMatrix::load(input)?;
        let lcs_simd = LcsSimd::load(input)?;
        let result = Self {
            ranges,
            matrix,
            lcs_simd,
            scan_bound,
            fallback_scan_overlap,
        };
        Ok(result)
    }

}

impl Pnsv for PnsvTuned {
    fn previous(&self, index: usize, target_length: usize) -> usize {
        if target_length <= self.ranges.max_target() {
            return self.ranges.previous(index, target_length);
        }

        if !self.matrix.is_empty() {
            if target_length <= self.matrix.max_target() - self.fallback_scan_overlap {
                return self.matrix.previous(index, target_length);
            }

            if target_length <= self.matrix.max_target() {
                let result = self.lcs_simd.scan_left_bounded(index, target_length as u8, self.scan_bound);
                return match result {
                    Ok(index) => index,
                    Err(continue_search_index) => {
                        self.matrix.previous(continue_search_index, target_length)
                    }
                };
            }
        }

        self.lcs_simd.scan_left(index, target_length as u8)
    }

    fn next(&self, index: usize, target_length: usize) -> usize {
        if target_length <= self.ranges.max_target() {
            return self.ranges.next(index, target_length);
        }

        if !self.matrix.is_empty() {
            if target_length <= self.matrix.max_target() - self.fallback_scan_overlap {
                return self.matrix.next(index, target_length);
            }

            if target_length <= self.matrix.max_target() {
                let result = self.lcs_simd.scan_right_bounded(index, target_length as u8, self.scan_bound);
                return match result {
                    Ok(index) => index,
                    Err(continue_search_index) => {
                        self.matrix.next(continue_search_index, target_length)
                    }
                };
            }
        }

        self.lcs_simd.scan_right(index, target_length as u8)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PnsvSafe {
    pub ranges: Ranges,
    pub wwt: WWT,
    pub lcs_simd: LcsSimd,
    pub scan_bound: usize,
}

impl PnsvSafe {
    pub const DEFAULT_SCAN_BOUND: usize = 16;

    pub fn new_with_default_values(extend: &impl ExtendRight, lcs: &LcsArray, max_k: usize) -> Self {
        Self::new(extend, lcs, max_k, Self::DEFAULT_SCAN_BOUND)
    }

    pub fn new(extend: &impl ExtendRight, lcs: &LcsArray, max_k: usize, scan_bound: usize) -> Self {
        let count = lcs.len();

        log::info!("[PnsvSafe::new] creating ranges...");
        let ranges = make_ranges(extend, count, max_k);
        let ranges_upper_bound = ranges.max_target();

        let mut wwt_option = None;
        let mut lcs_simd_option = None;

        let wwt_iterator = (0..count).map(|index| lcs.access(index));
        let lcs_simd_iterator = (0..count).map(|index| lcs.access(index) as u8);

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(2).build().unwrap();
        thread_pool.scope(|s| {
            s.spawn(|_| {
                log::info!("[PnsvSafe::new] creating windowed wavelet tree...");
                // The last few target lengths will have a small maximum range length, therefore,
                // it is better to not include them in the range of the windowed wavelet tree to
                // save on memory.
                let guaranteed_scan_hits_last_lengths = (usize::BITS - (scan_bound * LcsSimd::LANES - 1).leading_zeros()) as usize / 2;
                let wwt_upper_bound = max_k.saturating_sub(guaranteed_scan_hits_last_lengths);
                let wwt = if wwt_upper_bound > ranges_upper_bound {
                    let window_size = wwt_upper_bound - ranges_upper_bound + 1;
                    WWT::from_iterator(wwt_iterator, count, ranges_upper_bound, window_size)
                } else {
                    WWT::empty()
                };
                wwt_option = Some(wwt);
            });
            s.spawn(|_| {
                log::info!("[PnsvSafe::new] creating lcs simd...");
                let lcs_simd = LcsSimd::from_iterator(lcs_simd_iterator, count);
                lcs_simd_option = Some(lcs_simd);
            });
        });

        let wwt = wwt_option.expect("Creating windowed wavelet tree should not fail.");
        let lcs_simd = lcs_simd_option.expect("Creating lcs simd should not fail.");
        
        log::info!("[PnsvSafe::new] target length ranges: 1:{}:{}:..", ranges_upper_bound, wwt.max_target());

        Self {
            ranges,
            wwt,
            lcs_simd,
            scan_bound,
        }
    }

    pub fn serialize<W: std::io::Write>(&self, out: &mut W) -> std::io::Result<usize> {
        let mut written = 0;
        out.write_all(&(self.scan_bound as u64).to_le_bytes())?;
        written += size_of::<u64>();
        written += self.ranges.serialize(out)?;
        written += self.wwt.serialize(out)?;
        written += self.lcs_simd.serialize(out)?;
        Ok(written)
    }

    pub fn load<R: std::io::Read>(input: &mut R) -> std::io::Result<Self> {
        let scan_bound = u64::from_le(u64::load(input)?) as usize;
        let ranges = Ranges::load(input)?;
        let wwt = WWT::load(input)?;
        let lcs_simd = LcsSimd::load(input)?;
        let result = Self {
            ranges,
            wwt,
            lcs_simd,
            scan_bound,
        };
        Ok(result)
    }
}

impl Pnsv for PnsvSafe {
    fn previous(&self, index: usize, target_length: usize) -> usize {
        if target_length <= self.ranges.max_target() {
            return self.ranges.previous(index, target_length);
        }
        if !self.wwt.is_empty() {
            let result = self.lcs_simd.scan_left_bounded(index, target_length as u8, self.scan_bound);
            match result {
                Ok(index) => index,
                Err(continue_search_index) => {
                    self.wwt.previous(continue_search_index, target_length)
                }
            }
        } else {
            self.lcs_simd.previous(index, target_length)
        }
    }

    fn next(&self, index: usize, target_length: usize) -> usize {
        if target_length <= self.ranges.max_target() {
            return self.ranges.next(index, target_length);
        }
        if !self.wwt.is_empty() {
            let result = self.lcs_simd.scan_right_bounded(index, target_length as u8, self.scan_bound);
            match result {
                Ok(index) => index,
                Err(continue_search_index) => {
                    self.wwt.next(continue_search_index, target_length)
                }
            }
        } else {
            self.lcs_simd.next(index, target_length)
        }
    }
}

pub fn make_ranges(extend: &impl ExtendRight, count: usize, max_k: usize) -> Ranges {
    let mut ranges_upper_bound = 0;
    let mut bits_in_current_level_of_ranges = usize::BITS as usize * 4;
    let mut total_bits = bits_in_current_level_of_ranges;
    while total_bits < count {
        ranges_upper_bound += 1;
        bits_in_current_level_of_ranges *= 4;
        total_bits += bits_in_current_level_of_ranges;
    }
    ranges_upper_bound = ranges_upper_bound.min(Ranges::MAX_K);
    ranges_upper_bound = ranges_upper_bound.min(max_k);
    Ranges::new(extend, count, ranges_upper_bound)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BitPackedKmerSortingMem, SbwtIndexBuilder};
    use crate::{SbwtIndex, SubsetMatrix, LcsArray};

    fn setup(kmer_count: usize, max_k: usize) -> (SbwtIndex<SubsetMatrix>, LcsArray) {
        use rand_chacha::ChaCha20Rng;
        use rand_chacha::rand_core::SeedableRng;
        use rand_chacha::rand_core::RngCore;
        let mut rng = ChaCha20Rng::from_seed([53; 32]);

        let mut seqs = Vec::<Vec<u8>>::new();
        for _ in 0..kmer_count {
            let kmer: Vec<u8> = (0..max_k).map(|_| match rng.next_u32() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }).collect();
            seqs.push(kmer);
        }

        seqs.sort();
        seqs.dedup();

        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(max_k).build_lcs(true)
            .build_select_support(true)
            .run_from_vecs(seqs.as_slice());

        (sbwt, lcs.unwrap())
    }

    #[test]
    fn serialize_and_load_pnsv_safe() {
        let kmer_count: usize = u16::MAX as usize;
        let max_k: usize = 20;
        let (sbwt, lcs) = setup(kmer_count, max_k);
        let pnsv_safe = PnsvSafe::new(&sbwt, &lcs, max_k, 1);
        let mut buffer = Vec::<u8>::new();
        let written = pnsv_safe.serialize(&mut buffer).unwrap();
        assert_eq!(buffer.len(), written);
        let pnsv_safe_loaded = PnsvSafe::load(&mut buffer.as_slice()).unwrap();
        assert_eq!(pnsv_safe, pnsv_safe_loaded);
    }

    #[test]
    fn serialize_and_load_pnsv_tuned() {
        let kmer_count: usize = u16::MAX as usize;
        let max_k: usize = 20;
        let (sbwt, lcs) = setup(kmer_count, max_k);
        let pnsv_tuned = PnsvTuned::new_with_default_values(&sbwt, &lcs, max_k);
        let mut buffer = Vec::<u8>::new();
        let written = pnsv_tuned.serialize(&mut buffer).unwrap();
        assert_eq!(buffer.len(), written);
        let pnsv_tuned_loaded = PnsvTuned::load(&mut buffer.as_slice()).unwrap();
        assert_eq!(pnsv_tuned, pnsv_tuned_loaded);
    }

    #[test]
    fn randomised_kmers() {
        let kmer_count: usize = u16::MAX as usize;
        let min_k: usize = 3;
        let max_k: usize = 8;
        let (sbwt, lcs) = setup(kmer_count, max_k);
        let pnsv_safe = PnsvSafe::new(&sbwt, &lcs, max_k, 1);
        println!("pnsv_safe.ranges.max_target(): {}", pnsv_safe.ranges.max_target());
        println!("pnsv_safe.wwt.max_target():    {}", pnsv_safe.wwt.max_target());
        let pnsv_tuned = PnsvTuned::new_with_default_values(&sbwt, &lcs, max_k);
        println!("pnsv_tuned.ranges.max_target(): {}", pnsv_tuned.ranges.max_target());
        println!("pnsv_tuned.matrix.max_target(): {}", pnsv_tuned.matrix.max_target());
        for target_length in min_k..=max_k {
            for i in 0..lcs.len() {
                assert_eq!(pnsv_safe.lcs_simd.previous(i, target_length), pnsv_safe.previous(i, target_length));
                assert_eq!(pnsv_safe.lcs_simd.previous(i, target_length), pnsv_tuned.previous(i, target_length));
            }
        }
    }

}
