// Code by Martin Kostadinov.

use crate::ExtendRight;
use crate::util::DNA_ALPHABET;

/// Precompute the ranges in the SBWT index for each suffix of length k up to a given number. Use
/// those range borders to perform ContractLeft. Each level has at most O(4^k) border values.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ranges {
    // note(mk): Maybe SIMD can be used here as well.
    pub levels: Vec<Vec<usize>>,
    /// The number of k-mers in the SBWT (including the dummy ones).
    pub n: usize,
}

impl Ranges {
    pub const MAX_K: usize = 8;

    pub fn new<E: ExtendRight>(index: &E, n: usize, max_k: usize) -> Self {
        let max_k = if max_k > Self::MAX_K {
            log::warn!("This structure should be used only for small values of k. Provided was: {}, max is: {}.", max_k, Self::MAX_K);
            Self::MAX_K
        } else {
            max_k
        };

        let max_k = max_k.max(1);
        let mut levels: Vec<Vec<usize>> = vec![];

        let entire_range = 0..n;
        let mut current_level: Vec<std::ops::Range<usize>> = vec![entire_range];
        let mut next_level: Vec<std::ops::Range<usize>> = vec![];

        let mut buffer = vec![];
        for i in 0..=max_k {
            for range in &current_level {
                if range.start != 0 {
                    // The other procedures will assume there is an imaginary 0 at the beginning of
                    // the level vectors.
                    buffer.push(range.start);
                }
                for c in DNA_ALPHABET {
                    let new_range = index.extend_right(range.clone(), c);
                    if !new_range.is_empty() {
                        next_level.push(new_range);
                    }
                }
            }

            // The other procedures will also assume that there is an imaginary 'n' element at the
            // end of the level vectors.

            // This sort should not take that much time given max_k is not very big.
            buffer.sort();
            if i == max_k {
                break;
            }
            levels.push(buffer.clone());

            std::mem::swap(&mut current_level, &mut next_level);
            next_level.clear();
        }

        levels.push(buffer);

        Self {
            levels,
            n,
        }
    }

    const SCAN_UPPER_BOUND: usize = 64;

    /// Given an index and a target value finds the index of the previous element in the LCS array
    /// which has value smaller than the target one. This is equivalent to finding a the maximum
    /// border index which is smaller than target_len in the corresponding level.
    pub fn previous(&self, index: usize, target_length: usize) -> usize {
        assert!(target_length < self.levels.len());
        if self.levels[target_length].len() >= Self::SCAN_UPPER_BOUND {
            let result = self.levels[target_length].binary_search(&index);
            // If the value is found, then the value under the index is smaller than the target
            // length. Otherwise we want the previous element.
            match result {
                Ok(_) => { return index; }
                Err(insertion_index) => {
                    if insertion_index == 0 {
                        return 0;
                    }
                    return self.levels[target_length][insertion_index - 1];
                }
            }
        }
        let mut answer = 0;
        for &item in &self.levels[target_length] {
            if item > index { break; }
            answer = item;
        }
        answer
    }

    /// Similar to the find_previous_range method. Expands the upper bound of the range for the
    /// ContractLeft operation.
    pub fn next(&self, index: usize, target_length: usize) -> usize {
        assert!(target_length < self.levels.len());
        if self.levels[target_length].len() >= Self::SCAN_UPPER_BOUND {
            let result = self.levels[target_length].binary_search(&index);
            // If the value is found, then the value under the index is smaller than the target
            // length. Otherwise, the slot where the index can be inserted contains the next
            // smaller value if the found index is in the bounds of the array.
            match result {
                Ok(_) => { return index; }
                Err(insertion_index) => {
                    if insertion_index >= self.levels[target_length].len() {
                        return self.n;
                    }
                    return self.levels[target_length][insertion_index];
                },
            };
        }
        for &item in &self.levels[target_length] {
            if item >= index { 
                return item;
            }
        }
        self.n
    }
}

impl super::Pnsv for Ranges {
    #[inline]
    fn previous(&self, index: usize, target_length: usize) -> usize {
        self.previous(index, target_length)
    }

    #[inline]
    fn next(&self, index: usize, target_length: usize) -> usize {
        self.next(index, target_length)
    }

    #[inline]
    fn max_target(&self) -> usize { self.levels.len() - 1 }
}

