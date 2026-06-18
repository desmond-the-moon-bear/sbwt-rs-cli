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
    pub const MAX_K: usize = 9;

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
        for _ in 0..max_k {
            buffer = vec![];

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
            levels.push(buffer);

            std::mem::swap(&mut current_level, &mut next_level);
            next_level.clear();
        }

        Self {
            levels,
            n,
        }
    }

    const SCAN_UPPER_BOUND: usize = 64;

    /// Given an index and a target value finds the index of the previous element in the LCS array
    /// which has value smaller than the target one. This is equivalent to finding a the maximum
    /// border index which is smaller than target_len in the corresponding level.
    pub fn previous(&self, index: usize, target_len: usize) -> usize {
        let mut answer = 0;
        // To be consistent with the ContractLeft semantics this method accepts a target value for
        // which to expand the interval. However, we want to expand the range to a level (i.e.
        // value of k) which is smaller than it.
        let level_to_search = target_len - 1;

        if self.levels[level_to_search].len() >= Self::SCAN_UPPER_BOUND {
            let result = self.levels[level_to_search].binary_search(&index);
            // If the value is found we want the previous element (if it exists, otherwise 0). If
            // the value is not found, then the value at the found index is greater than the
            // parameter index and the previous value is smaller. We want the previous value either
            // way.
            let found_index = match result {
                Ok(value) => value,
                Err(value) => value,
            };
            if found_index == 0 {
                return 0;
            }
            return self.levels[level_to_search][found_index - 1];
        }

        for &item in &self.levels[level_to_search] {
            if item >= index { break; }
            answer = item;
        }
        answer
    }

    /// Similar to the find_previous_range method. Expands the upper bound of the range for the
    /// ContractLeft operation.
    pub fn next(&self, index: usize, target_len: usize) -> usize {
        let level_to_search = target_len - 1;

        if self.levels[level_to_search].len() >= Self::SCAN_UPPER_BOUND {
            let result = self.levels[level_to_search].binary_search(&index);
            // If the value is found we want the next element (if it exists, otherwise n). If
            // the value is not found, then the value at the found index is greater than the
            // parameter index so we return that value.
            match result {
                Ok(matched_index) => {
                    if matched_index + 1 >= self.levels[level_to_search].len() {
                        return self.n;
                    }
                    return self.levels[level_to_search][matched_index + 1];
                }
                Err(insertion_index) => {
                    if insertion_index >= self.levels[level_to_search].len() {
                        return self.n;
                    }
                    return self.levels[level_to_search][insertion_index];
                },
            };
        }

        for &item in &self.levels[level_to_search] {
            if item > index { 
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
    fn max_target(&self) -> usize { self.levels.len() }
}

