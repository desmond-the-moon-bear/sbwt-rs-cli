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
pub use matrix::MatrixSux as PnsvMatrixSux;
pub use ranges::Ranges;
pub use scan::AugmentedBoundedScan as ABS;
pub use scan::LcsSimd;
pub use scan::ScanWithFallback;
pub use wavelet::WindowedWaveletTree as WWT;

/// Previous/Next Smaller value.
pub trait Pnsv {
    fn previous(&self, index: usize, target_length: usize) -> usize;
    fn next(&self, index: usize, target_length: usize) -> usize;
    fn override_contract_left(&self) -> bool { false }
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

pub struct PnsvDyn<'a, const LEVELS: usize> {
    pub structures: [&'a dyn Pnsv; LEVELS],
}

impl<'a, const LEVELS: usize> ContractLeft for PnsvDyn<'a, LEVELS> {
    #[allow(non_snake_case)]
    fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {
        for i in 0..self.structures.len() - 1 {
            if target_len <= self.structures[i].max_target() {
                return self.structures[i].contract_left(I, target_len);
            }
        }
        self.structures[self.structures.len() - 1].contract_left(I, target_len)
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

pub fn pnsv_simd_fallback_matrix(extend: &impl ExtendRight, lcs: &LcsArray, scan_bound: usize) -> PnsvDynOwned {
    let count = lcs.len();

    let mut structures: Vec<Box<dyn Pnsv>> = vec![];

    log::info!("[pnsv_simd_fallback_matrix] creating ranges...");

    let mut ranges_upper_bound = 0;
    let mut bits_in_current_level_of_ranges = usize::BITS as usize * 4;
    while bits_in_current_level_of_ranges < count {
        ranges_upper_bound += 1;
        bits_in_current_level_of_ranges *= 4;
    }
    ranges_upper_bound = ranges_upper_bound.min(Ranges::MAX_K);
    let ranges = Ranges::new(extend, count, ranges_upper_bound);
    let ranges_box = Box::new(ranges);
    structures.push(ranges_box);

    let iterator = (0..count).map(|index| lcs.access(index) as u8);

    // Experimentally determined, the ranges shrink 4-fold initially and afterwards the ratio
    // between two consecutive average region lengths (e.g. region lengths for k=7 and k=8) becomes
    // less than 4. Around that time as well the average range length becomes around 200 i.e. the
    // target length. Therefore a simple logarithm should find the target length upper bound for
    // the matrix solution.
    let log_4 = (usize::BITS - count.leading_zeros()).div_ceil(2) as usize;

    // log_4(count / 200) == log_4(count) - log_4(200)
    let mut matrix_upper_bound = log_4 - TARGET_LENGTH_LOG_4_FLOOR; 
    matrix_upper_bound = matrix_upper_bound.min(ranges_upper_bound + 1 + matrix::MAX_ROWS);

    if matrix_upper_bound > ranges_upper_bound {
        log::info!("[pnsv_simd_fallback_matrix] creating matrix...");
        // let matrix = PnsvMatrix::from_iterator(iterator.clone(), count, ranges_upper_bound + 1, matrix_upper_bound);
        let matrix = PnsvMatrixSux::from_iterator(iterator.clone(), count, ranges_upper_bound + 1, matrix_upper_bound);

        log::info!("[pnsv_simd_fallback_matrix] creating lcs simd...");
        let lcs_simd = LcsSimd::from_iterator(iterator, count);

        let swf = ScanWithFallback::new(lcs_simd, scan_bound, matrix);
        let swf_box = Box::new(swf);
        structures.push(swf_box);
    } else {
        log::info!("[pnsv_simd_fallback_matrix] creating lcs simd...");
        let lcs_simd = LcsSimd::from_iterator(iterator, count);
        let lcs_simd_box = Box::new(lcs_simd);
        structures.push(lcs_simd_box);
    }

    if matrix_upper_bound > ranges_upper_bound {
        log::info!("[pnsv_simd_fallback_matrix] target length ranges: 1:{}:{}:..", ranges_upper_bound, matrix_upper_bound);
    } else {
        log::info!("[pnsv_simd_fallback_matrix] target length ranges: 1:{}:..", ranges_upper_bound);
    }

    PnsvDynOwned {
        structures,
    }
}

pub fn pnsv_matrix_simd(extend: &impl ExtendRight, lcs: &LcsArray) -> PnsvDynOwned {
    let count = lcs.len();

    let mut structures: Vec<Box<dyn Pnsv>> = vec![];

    log::info!("[pnsv_matrix_simd] creating ranges...");
    let mut ranges_upper_bound = 0;
    let mut bits_in_current_level_of_ranges = usize::BITS as usize * 4;
    while bits_in_current_level_of_ranges < count {
        ranges_upper_bound += 1;
        bits_in_current_level_of_ranges *= 4;
    }
    ranges_upper_bound = ranges_upper_bound.min(Ranges::MAX_K);
    let ranges = Ranges::new(extend, count, ranges_upper_bound);
    let ranges_box = Box::new(ranges);
    structures.push(ranges_box);

    let iterator = (0..count).map(|index| lcs.access(index) as u8);

    // Experimentally determined, the ranges shrink 4-fold initially and afterwards the ratio
    // between two consecutive average region lengths (e.g. region lengths for k=7 and k=8) becomes
    // less than 4. Around that time as well the average range length becomes around 200 i.e. the
    // target length. Therefore a simple logarithm should find the target length upper bound for
    // the matrix solution.
    let log_4 = (usize::BITS - count.leading_zeros()).div_ceil(2) as usize;

    // log_4(count / 200) == log_4(count) - log_4(200)
    let mut matrix_upper_bound = log_4 - TARGET_LENGTH_LOG_4_FLOOR; 
    matrix_upper_bound = matrix_upper_bound.min(ranges_upper_bound + 1 + matrix::MAX_ROWS);

    if matrix_upper_bound > ranges_upper_bound {
        log::info!("[pnsv_matrix_simd] creating matrix...");
        // let matrix = PnsvMatrix::from_iterator(iterator.clone(), count, ranges_upper_bound + 1, matrix_upper_bound);
        let matrix = PnsvMatrixSux::from_iterator(iterator.clone(), count, ranges_upper_bound + 1, matrix_upper_bound);
        let matrix_box = Box::new(matrix);
        structures.push(matrix_box);
    }

    log::info!("[pnsv_matrix_simd] creating lcs simd...");
    let lcs_simd = LcsSimd::from_iterator(iterator, count);
    let lcs_simd_box = Box::new(lcs_simd);
    structures.push(lcs_simd_box);

    log::info!("[pnsv_matrix_simd] target length ranges: 1:{}:{}:..", ranges_upper_bound, matrix_upper_bound);

    PnsvDynOwned {
        structures,
    }
}

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

pub struct PnsvTuned {
    pub ranges: Ranges,
    pub matrix: PnsvMatrixSux,
    pub lcs_simd: LcsSimd,
    pub scan_bound: usize,
    pub fallback_scan_overlap: usize,
}

impl PnsvTuned {
    pub const DEFAULT_SCAN_BOUND: usize = 4;
    pub const DEFAULT_FALLBACK_OVERLAP: usize = 2;

    pub fn new_with_default_values(extend: &impl ExtendRight, lcs: &LcsArray) -> Self {
        Self::new(extend, lcs, Self::DEFAULT_SCAN_BOUND, Self::DEFAULT_FALLBACK_OVERLAP)
    }

    pub fn new(extend: &impl ExtendRight, lcs: &LcsArray, scan_bound: usize, fallback_scan_overlap: usize) -> Self {
        let count = lcs.len();

        log::info!("[PnsvTuned::new] creating ranges...");
        let mut ranges_upper_bound = 0;
        let mut bits_in_current_level_of_ranges = usize::BITS as usize * 4;
        while bits_in_current_level_of_ranges < count {
            ranges_upper_bound += 1;
            bits_in_current_level_of_ranges *= 4;
        }
        ranges_upper_bound = ranges_upper_bound.min(Ranges::MAX_K);
        let ranges = Ranges::new(extend, count, ranges_upper_bound);

        let iterator = (0..count).map(|index| lcs.access(index) as u8);

        let log_4 = (usize::BITS - count.leading_zeros()).div_ceil(2) as usize;
        let mut matrix_upper_bound = log_4 - TARGET_LENGTH_LOG_4_FLOOR; 
        matrix_upper_bound = matrix_upper_bound.min(ranges_upper_bound + 1 + matrix::MAX_ROWS);

        log::info!("[PnsvTuned::new] creating matrix...");
        let matrix = PnsvMatrixSux::from_iterator(iterator.clone(), count, ranges_upper_bound + 1, matrix_upper_bound);
        assert!(matrix.max_target() >= fallback_scan_overlap, "Make sure the fallback scan overlap has a reasonably small value (i.e. most 5).");

        log::info!("[PnsvTuned::new] creating lcs simd...");
        let lcs_simd = LcsSimd::from_iterator(iterator, count);

        log::info!("[PnsvTuned::new] target length ranges: 1:{}:{}:..", ranges_upper_bound, matrix_upper_bound);

        Self {
            ranges,
            matrix,
            lcs_simd,
            scan_bound,
            fallback_scan_overlap,
        }
    }
}

impl Pnsv for PnsvTuned {
    fn previous(&self, index: usize, target_length: usize) -> usize {
        if target_length <= self.ranges.max_target() {
            return self.ranges.previous(index, target_length);
        }

        if target_length <= self.matrix.max_target() - self.fallback_scan_overlap {
            return self.ranges.previous(index, target_length);
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

        self.lcs_simd.scan_left(index, target_length as u8)
    }

    fn next(&self, index: usize, target_length: usize) -> usize {
        if target_length <= self.ranges.max_target() {
            return self.ranges.next(index, target_length);
        }

        if target_length <= self.matrix.max_target() - self.fallback_scan_overlap {
            return self.ranges.next(index, target_length);
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

        self.lcs_simd.scan_right(index, target_length as u8)
    }
}

#[cfg(test)]
mod tests {
    // todo(mk): test ScanWithFallback...
}
