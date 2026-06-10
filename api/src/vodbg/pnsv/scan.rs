// Code by Martin Kostadinov.

use crate::ContractLeft;

pub type Word = wide::u8x32;

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
        let mut words = vec![];
        let mut array = Self::ZERO;
        for (i, item) in input.enumerate() {
            if (i + 1) % Self::LANES == 0 {
                words.push(Word::new(array));
                array = Self::ZERO;
            }
            array[i % Self::LANES] = item.into();
        }
        words.push(Word::new(array));

        Self { words, n }
    }

    pub fn scan_left(&self, index: usize, target_length: u8) -> usize {
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
        let word_index = index / Self::LANES;
        let index_in_word = (index % Self::LANES) + 1;

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

