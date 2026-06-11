// Code by Martin Kostadinov.

use crate::ContractLeft;

pub type Word = wide::u8x32;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LcsSimd {
    pub words: Vec<Word>,
    pub n: usize,
}

impl LcsSimd {
    const ZERO: [u8; Word::LANES as usize] = [0; Word::LANES as usize];
    const LANES: usize = Word::LANES as usize;

    pub fn from_iterator<T, I>(input: I, n: usize) -> Self
    where 
        T: Into<u8>,
        I: Iterator<Item = T>
    {
        #[allow(clippy::manual_div_ceil)]
        let mut words = Vec::with_capacity((n + Self::LANES - 1) / Self::LANES);
        let mut array = Self::ZERO;
        for (i, item) in input.enumerate() {
            if i != 0 && i % Self::LANES == 0 {
                words.push(Word::new(array));
                array = Self::ZERO;
            }
            array[i % Self::LANES] = item.into();
        }
        words.push(Word::new(array));

        Self { words, n }
    }

    pub fn scan_left(&self, index: usize, target_length: u8) -> usize {
        if index >= self.n {
            return self.n;
        }

        let word_index = index / Self::LANES;
        let index_in_word = index % Self::LANES;

        // Scan the values in the SIMD word the index is located in individually.
        let near_word = self.words[word_index].as_array();
        for i in (0..=index_in_word).rev() {
            if near_word[i] < target_length {
                return word_index * Self::LANES + i; 
            }
        }

        // Scan the rest of the words using SIMD operations.
        for w in (0..word_index).rev() {
            let comparison_result = self.words[w].simd_lt(target_length);
            if !comparison_result.any() {
                continue;
            }
            let bitmask = comparison_result.to_bitmask();
            let rightmost_smaller_element = Self::LANES - 1 - bitmask.leading_zeros() as usize;
            return w * Self::LANES + rightmost_smaller_element;
        }

        0
    }

    pub fn scan_right(&self, index: usize, target_length: u8) -> usize {
        if index >= self.n {
            return self.n;
        }

        let word_index = index / Self::LANES;
        let index_in_word = index % Self::LANES;

        // Similarly to the scan_left procedure, first scan values in the word the index is located
        // in individually.
        let near_word = self.words[word_index].as_array();
        for i in index_in_word..Self::LANES {
            if near_word[i] < target_length {
                return word_index * Self::LANES + i; 
            }
        }

        // Then scan the rest of the words using SIMD.
        for w in (word_index + 1)..self.words.len() {
            let comparison_result = self.words[w].simd_lt(target_length);
            if !comparison_result.any() {
                continue;
            }
            let bitmask = comparison_result.to_bitmask();
            let leftmost_smaller_element = bitmask.trailing_zeros() as usize;
            let result = w * Self::LANES + leftmost_smaller_element;
            return result.min(self.n);
        }

        self.n
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.n
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
}

impl ContractLeft for LcsSimd {
    fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {
        let new_start = self.scan_left(I.start, target_len as u8);
        let new_end = self.scan_right(I.end, target_len as u8);
        new_start..new_end
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lcs_simd_words_are_correct() {
        let items: &[u8] = &[
             0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
            10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
            20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
            20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
        ];

        let lcs_simd = LcsSimd::from_iterator(items.iter().cloned(), items.len());
        for i in 0..items.len() {
            let item_in_word = lcs_simd.words[i / LcsSimd::LANES].as_array()[i % LcsSimd::LANES];
            assert_eq!(items[i], item_in_word);
        }
    }

    #[test]
    fn lcs_simd_scan_left() {
        let items: &[u8] = &[
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 5, 5, 5, 5, 5, // 25: 1

            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, // 57: 2

            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, // 89: 3

            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, // 121: 4
        ];

        let lcs_simd = LcsSimd::from_iterator(items.iter().cloned(), items.len());
        assert_eq!(lcs_simd.scan_left(125, 5), 121, "Result in the same SIMD word failed.");
        assert_eq!(lcs_simd.scan_left(100, 5), 89, "Result in the previous SIMD word failed.");
        assert_eq!(lcs_simd.scan_left(100, 4), 89, "Smaller than 4 failed.");
        assert_eq!(lcs_simd.scan_left(100, 3), 57, "Smaller than 3 failed.");
        assert_eq!(lcs_simd.scan_left(100, 2), 25, "Smaller than 2 failed.");
        assert_eq!(lcs_simd.scan_left(100, 1), 0, "Scan to beginning failed.");
    }

    #[test]
    fn lcs_simd_scan_right() {
        let items: &[u8] = &[
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, // 25: 4

            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 5, 5, 5, 5, 5, 5, // 57: 3

            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 5, 5, 5, 5, 5, 5, // 89: 2

            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 5, 5, 5, 5, 5, 5, // 121: 1
        ];

        let lcs_simd = LcsSimd::from_iterator(items.iter().cloned(), items.len());
        assert_eq!(lcs_simd.scan_right(20, 5), 25, "Result in the same SIMD word failed.");
        assert_eq!(lcs_simd.scan_right(30, 5), 57, "Result in the next SIMD word failed.");
        assert_eq!(lcs_simd.scan_right(30, 4), 57, "Smaller than 4 failed.");
        assert_eq!(lcs_simd.scan_right(30, 3), 89, "Smaller than 3 failed.");
        assert_eq!(lcs_simd.scan_right(30, 2), 121, "Smaller than 2 failed.");
        assert_eq!(lcs_simd.scan_right(30, 1), items.len(), "Scan to the end failed.");
    }
}

