// Code by Martin Kostadinov.

use super::bp::nearest_neighbor_dictionary::NearestNeighbourDictionary as NND;
use super::Pnsv;

pub type Word = wide::u8x32;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LcsSimd {
    pub words: Vec<Word>,
    pub n: usize,
}

impl LcsSimd {
    const ZERO: [u8; Word::LANES as usize] = [0; Word::LANES as usize];
    pub const LANES: usize = Word::LANES as usize;

    pub fn from_iterator<T, I>(input: I, n: usize) -> Self
    where
        T: Into<u8>,
        I: Iterator<Item = T>,
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
            return 0;
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

    /// Bound gives the maximum number of words to be scanned in addition to the one the index is located.
    pub fn scan_left_bounded(&self, index: usize, target_length: u8, bound: usize) -> Result<usize, usize> {
        if index >= self.n {
            return Err(self.n);
        }

        let word_index = index / Self::LANES;
        let index_in_word = index % Self::LANES;

        // Scan the values in the SIMD word the index is located in individually.
        let near_word = self.words[word_index].as_array();
        for i in (0..=index_in_word).rev() {
            if near_word[i] < target_length {
                return Ok(word_index * Self::LANES + i);
            }
        }

        let lower_bound_word_index = word_index.saturating_sub(bound);

        // Scan the rest of the words using SIMD operations.
        for w in (lower_bound_word_index..word_index).rev() {
            let comparison_result = self.words[w].simd_lt(target_length);
            if !comparison_result.any() {
                continue;
            }
            let bitmask = comparison_result.to_bitmask();
            let rightmost_smaller_element = Self::LANES - 1 - bitmask.leading_zeros() as usize;
            return Ok(w * Self::LANES + rightmost_smaller_element);
        }

        Err(lower_bound_word_index * Self::LANES)
    }

    /// Bound gives the maximum number of words to be scanned in addition to the one the index is located.
    pub fn scan_right_bounded(&self, index: usize, target_length: u8, bound: usize) -> Result<usize, usize> {
        if index >= self.n {
            return Err(self.n);
        }

        let word_index = index / Self::LANES;
        let index_in_word = index % Self::LANES;

        // Similarly to the scan_left procedure, first scan values in the word the index is located
        // in individually.
        let near_word = self.words[word_index].as_array();
        for i in index_in_word..Self::LANES {
            if near_word[i] < target_length {
                return Ok(word_index * Self::LANES + i);
            }
        }

        let mut upper_bound_word_index = (word_index + bound).min(self.words.len());
        if upper_bound_word_index < self.words.len() {
            upper_bound_word_index += 1;
        }

        // Then scan the rest of the words using SIMD.
        for w in (word_index + 1)..upper_bound_word_index {
            let comparison_result = self.words[w].simd_lt(target_length);
            if !comparison_result.any() {
                continue;
            }
            let bitmask = comparison_result.to_bitmask();
            let leftmost_smaller_element = bitmask.trailing_zeros() as usize;
            let result = w * Self::LANES + leftmost_smaller_element;
            return Ok(result.min(self.n));
        }

        Err(upper_bound_word_index * Self::LANES - 1)
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

impl Pnsv for LcsSimd {
    #[inline]
    fn previous(&self, index: usize, target_length: usize) -> usize {
        self.scan_left(index, target_length as u8)
    }

    #[inline]
    fn next(&self, index: usize, target_length: usize) -> usize {
        self.scan_right(index, target_length as u8)
    }
}

pub struct AugmentedBoundedScan<'a> {
    pub lcs_simd: &'a LcsSimd,
    pub target_length_lower: usize,
    pub target_length_upper: usize,
    pub scan_word_bound: usize,
    pub levels: Vec<NND<16>>,
}

impl<'a> AugmentedBoundedScan<'a> {
    pub fn from_iterator<T, I>(lcs_simd: &'a LcsSimd, input: I, n: usize, scan_word_bound: usize, target_length_lower: usize, target_length_upper: usize) -> Self
    where
        T: Into<usize>,
        I: Iterator<Item = T> + Clone,
    {
        let mut levels = Vec::with_capacity(target_length_upper - target_length_lower + 1);

        for target_length in target_length_lower..=target_length_upper {
            let select = SelectFromSimdWord::new(input.clone(), n, target_length);
            let nnd: NND<16> = NND::new(select, n);
            levels.push(nnd);
        }

        Self {
            lcs_simd,
            target_length_lower,
            target_length_upper,
            scan_word_bound,
            levels
        }
    }
}

impl<'a> Pnsv for AugmentedBoundedScan<'a> {
    fn previous(&self, index: usize, target_length: usize) -> usize {
        let result = self.lcs_simd.scan_left_bounded(index, target_length as u8, self.scan_word_bound);
        match result {
            Ok(index) => index,
            Err(continue_search_index) => {
                let nnd_index = target_length - self.target_length_lower;
                if nnd_index < self.levels.len() {
                    self.levels[nnd_index].previous(continue_search_index)
                } else {
                    index
                }
            }
        }
    }

    fn next(&self, index: usize, target_length: usize) -> usize {
        let result = self.lcs_simd.scan_right_bounded(index, target_length as u8, self.scan_word_bound);
        match result {
            Ok(index) => index,
            Err(continue_search_index) => {
                let nnd_index = target_length - self.target_length_lower;
                if nnd_index < self.levels.len() {
                    self.levels[nnd_index].next(continue_search_index)
                } else {
                    index
                }
            }
        }
    }
}

// An iterator over the indices of the first and last value in each simd word smaller than a given
// target length.
struct SelectFromSimdWord<T, I>
where 
    T: Into<usize>,
    I: Iterator<Item = T> + Clone
{
    target_length: usize,
    index: usize,
    previous_index_of_smaller: usize,
    last_yielded: usize,
    count: usize,
    pending: Option<usize>,
    underlying_values: I,
}

impl<T, I> SelectFromSimdWord<T, I> 
where 
    T: Into<usize>,
    I: Iterator<Item = T> + Clone
{
    fn new(underlying_values: I, count: usize, target_length: usize) -> Self {
        Self {
            target_length,
            index: 0,
            previous_index_of_smaller: count,
            last_yielded: count + 1,
            count,
            pending: None,
            underlying_values,
        }
    }
}

impl<T, I> Iterator for SelectFromSimdWord<T, I>
where 
    T: Into<usize>,
    I: Iterator<Item = T> + Clone
{
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        if self.pending.is_some() {
            self.last_yielded = self.pending.unwrap();
            return self.pending.take();
        }

        if self.index == self.count {
            return None;
        }

        #[allow(clippy::while_let_on_iterator)]
        while let Some(value) = self.underlying_values.next() {
            let value = value.into();

            let index = self.index;
            self.index += 1;

            if value >= self.target_length {
                continue;
            }

            let previous_index_of_smaller = self.previous_index_of_smaller;
            self.previous_index_of_smaller = index;

            if previous_index_of_smaller < self.count {
                let previous_word = previous_index_of_smaller / LcsSimd::LANES;
                let current_word = index / LcsSimd::LANES;

                if previous_word != current_word {
                    if self.last_yielded == previous_index_of_smaller {
                        self.last_yielded = index;
                        return Some(index);
                    } else {
                        self.pending = Some(index);
                        // note(mk): This next line is probably not necessary.
                        self.last_yielded = previous_index_of_smaller;
                        return Some(previous_index_of_smaller);
                    }
                }
            } else {
                return Some(index);
            }
        }

        if self.previous_index_of_smaller < self.count && self.previous_index_of_smaller != self.last_yielded {
            return Some(self.previous_index_of_smaller);
        }

        None
    }
}

// note(mk): There were some errors with the automatically derived trait. 
impl<T, I> Clone for SelectFromSimdWord<T, I>
where 
    T: Into<usize>,
    I: Iterator<Item = T> + Clone
{
    fn clone(&self) -> Self {
        Self {
            target_length: self.target_length,
            index: self.index,
            previous_index_of_smaller: self.previous_index_of_smaller,
            last_yielded: self.last_yielded,
            count: self.count,
            pending: self.pending,
            underlying_values: self.underlying_values.clone()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lcs_simd_words_are_correct() {
        let items: &[u8] = &[
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
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
        assert_eq!(
            lcs_simd.previous(125, 5),
            121,
            "Result in the same SIMD word failed."
        );
        assert_eq!(
            lcs_simd.previous(100, 5),
            89,
            "Result in the previous SIMD word failed."
        );
        assert_eq!(lcs_simd.previous(100, 4), 89, "Smaller than 4 failed.");
        assert_eq!(lcs_simd.previous(100, 3), 57, "Smaller than 3 failed.");
        assert_eq!(lcs_simd.previous(100, 2), 25, "Smaller than 2 failed.");
        assert_eq!(lcs_simd.previous(100, 1), 0, "Scan to beginning failed.");
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
        assert_eq!(
            lcs_simd.next(20, 5),
            25,
            "Result in the same SIMD word failed."
        );
        assert_eq!(
            lcs_simd.next(30, 5),
            57,
            "Result in the next SIMD word failed."
        );
        assert_eq!(lcs_simd.next(30, 4), 57, "Smaller than 4 failed.");
        assert_eq!(lcs_simd.next(30, 3), 89, "Smaller than 3 failed.");
        assert_eq!(lcs_simd.next(30, 2), 121, "Smaller than 2 failed.");
        assert_eq!(lcs_simd.next(30, 1), items.len(), "Scan to the end failed.");
    }

    #[test]
    fn lcs_simd_scan_left_bounded() {
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
        let max_bound = lcs_simd.words.len();
        assert_eq!(
            lcs_simd.scan_left_bounded(125, 5, max_bound),
            Ok(121),
            "Result in the same SIMD word failed."
        );
        assert_eq!(
            lcs_simd.scan_left_bounded(100, 5, max_bound),
            Ok(89),
            "Result in the previous SIMD word failed."
        );
        assert_eq!(
            lcs_simd.scan_left_bounded(100, 4, max_bound),
            Ok(89),
            "Smaller than 4 failed."
        );
        assert_eq!(
            lcs_simd.scan_left_bounded(100, 3, max_bound),
            Ok(57),
            "Smaller than 3 failed."
        );
        assert_eq!(
            lcs_simd.scan_left_bounded(100, 2, max_bound),
            Ok(25),
            "Smaller than 2 failed."
        );
        assert!(
            lcs_simd.scan_left_bounded(100, 1, max_bound).is_err(),
            "Scan to beginning failed."
        );

        assert!(lcs_simd.scan_left_bounded(100, 5, 0).is_err());
        assert!(lcs_simd.scan_left_bounded(100, 4, 0).is_err());
        assert!(lcs_simd.scan_left_bounded(100, 3, 1).is_err());
        assert!(lcs_simd.scan_left_bounded(100, 2, 2).is_err());

        assert_eq!(lcs_simd.scan_left_bounded(125, 5, 0), Ok(121));
        assert_eq!(lcs_simd.scan_left_bounded(100, 5, 1), Ok(89));
        assert_eq!(lcs_simd.scan_left_bounded(100, 4, 1), Ok(89));
        assert_eq!(lcs_simd.scan_left_bounded(100, 3, 2), Ok(57));
        assert_eq!(lcs_simd.scan_left_bounded(100, 2, 3), Ok(25));
    }

    #[test]
    fn lcs_simd_scan_right_bounded() {
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
        let max_bound = lcs_simd.words.len();
        assert_eq!(
            lcs_simd.scan_right_bounded(20, 5, max_bound),
            Ok(25),
            "Result in the same SIMD word failed."
        );
        assert_eq!(
            lcs_simd.scan_right_bounded(30, 5, max_bound),
            Ok(57),
            "Result in the next SIMD word failed."
        );
        assert_eq!(
            lcs_simd.scan_right_bounded(30, 4, max_bound),
            Ok(57),
            "Smaller than 4 failed."
        );
        assert_eq!(
            lcs_simd.scan_right_bounded(30, 3, max_bound),
            Ok(89),
            "Smaller than 3 failed."
        );
        assert_eq!(
            lcs_simd.scan_right_bounded(30, 2, max_bound),
            Ok(121),
            "Smaller than 2 failed."
        );
        assert!(
            lcs_simd.scan_right_bounded(30, 1, max_bound).is_err(),
            "Scan to the end failed."
        );

        assert!(lcs_simd.scan_right_bounded(30, 5, 0).is_err());
        assert!(lcs_simd.scan_right_bounded(30, 4, 0).is_err());
        assert!(lcs_simd.scan_right_bounded(30, 3, 1).is_err());
        assert!(lcs_simd.scan_right_bounded(30, 2, 2).is_err());

        assert_eq!(lcs_simd.scan_right_bounded(20, 5, 0), Ok(25));
        assert_eq!(lcs_simd.scan_right_bounded(30, 5, 1), Ok(57));
        assert_eq!(lcs_simd.scan_right_bounded(30, 4, 1), Ok(57));
        assert_eq!(lcs_simd.scan_right_bounded(30, 3, 2), Ok(89));
        assert_eq!(lcs_simd.scan_right_bounded(30, 2, 3), Ok(121));
    }

    #[test]
    fn select_from_simd_word() {
        let items: &[u8] = &[
            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,

            1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,

            1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        ];

        let less_than_5: Vec<usize> = SelectFromSimdWord::new(items.iter().cloned(), items.len(), 5).collect();
        let less_than_4: Vec<usize> = SelectFromSimdWord::new(items.iter().cloned(), items.len(), 4).collect();
        let less_than_3: Vec<usize> = SelectFromSimdWord::new(items.iter().cloned(), items.len(), 3).collect();
        let less_than_2: Vec<usize> = SelectFromSimdWord::new(items.iter().cloned(), items.len(), 2).collect();

        assert_eq!(less_than_5, &[0, 31, 32, 63, 64, 95, 96, 127, 160, 192]);
        assert_eq!(less_than_4, &[1, 30, 33, 62, 65, 94, 97, 126, 160, 192]);
        assert_eq!(less_than_3, &[2, 29, 34, 61, 66, 93, 98, 125, 160, 192]);
        assert_eq!(less_than_2, &[3, 28, 35, 60, 67, 92, 99, 124, 160, 192]);
    }

    #[test]
    fn bounded_scan_with_fallback() {
        let items: &[u8] = &[
            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            4, 3, 2, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 3, 4,

            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,

            1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,

            1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
        ];

        let target_length_lower = 2;
        let target_length_upper = 5;

        let lcs_simd = LcsSimd::from_iterator(items.iter().cloned(), items.len());
        let abs = AugmentedBoundedScan::from_iterator(&lcs_simd, items.iter().cloned(), lcs_simd.len(), 0, target_length_lower, target_length_upper);

        for target_length in target_length_lower..=target_length_upper {
            for i in 0..items.len() {
                assert_eq!(
                    abs.previous(i, target_length),
                    lcs_simd.previous(i, target_length),
                    "previous; i: {}; target_length: {}",
                    i, target_length
                );
                assert_eq!(
                    abs.next(i, target_length),
                    lcs_simd.next(i, target_length),
                    "next; i: {}; target_length: {}",
                    i, target_length
                );
            }
        }

    }
}
