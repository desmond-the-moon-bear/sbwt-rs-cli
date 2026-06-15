// Code by Martin Kostadinov.

use super::Pnsv;

use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::ops::{BitVec, Rank, Select};
use simple_sds_sbwt::raw_vector::{AccessRaw, RawVector};

/// Supports previous smaller value for a range of target lengths.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Matrix {
    pub lower_bound: usize,
    pub upper_bound: usize,
    pub width: usize,
    pub rows: Vec<BitVector>,
}

impl Matrix {
    pub fn from_iterator<T, I>(input: I, width: usize, lower_bound: usize, upper_bound: usize) -> Self
    where 
        T: Into<usize>,
        I: Iterator<Item = T>,
    {
        let mut rows_raw: Vec<RawVector> = Vec::with_capacity(upper_bound - lower_bound + 1);

        for _ in lower_bound..=upper_bound {
            rows_raw.push(RawVector::with_len(width, false));
        }

        for (index, item) in input.enumerate() {
            let value = item.into();
            for target_length in lower_bound..=upper_bound {
                let row_index = target_length - lower_bound;
                if value < target_length {
                    rows_raw[row_index].set_bit(index, true);
                }
            }
        }
        
        let mut rows: Vec<BitVector> = Vec::with_capacity(rows_raw.len());
        for row in rows_raw {
            let mut bit_vector: BitVector = row.into();
            bit_vector.enable_rank();
            bit_vector.enable_select();
            rows.push(bit_vector);
        }

        Self {
            lower_bound,
            upper_bound,
            width,
            rows,
        }
    }

    pub fn previous(&self, index: usize, target_length: usize) -> usize {
        if index >= self.width {
            return self.width;
        }
        if target_length < self.lower_bound || target_length > self.upper_bound {
            return 0;
        }
        let row_index = target_length - self.lower_bound;
        let is_one = self.rows[row_index].get(index);
        if is_one {
            return index;
        }
        let one_rank = self.rows[row_index].rank(index);
        if one_rank == 0 {
            return 0;
        }
        // There is at least one smaller value before.
        self.rows[row_index].select(one_rank - 1).unwrap()
    }

    pub fn next(&self, index: usize, target_length: usize) -> usize {
        if index >= self.width {
            return self.width;
        }
        if target_length < self.lower_bound || target_length > self.upper_bound {
            return self.width;
        }
        let row_index = target_length - self.lower_bound;
        let is_one = self.rows[row_index].get(index);
        if is_one {
            return index;
        }
        let one_rank = self.rows[row_index].rank(index);
        if one_rank == self.rows[row_index].count_ones() {
            return self.width;
        }
        // There is at least one smaller value after.
        self.rows[row_index].select(one_rank).unwrap()
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.width
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.width == 0
    }
}

impl Pnsv for Matrix {
    #[inline]
    fn previous(&self, index: usize, target_length: usize) -> usize {
        self.previous(index, target_length)
    }

    #[inline]
    fn next(&self, index: usize, target_length: usize) -> usize {
        self.next(index, target_length)
    }

    fn max_target(&self) -> usize {
        self.upper_bound
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn previous(items: &[usize], index: usize, target_length: usize, lower_bound: usize, upper_bound: usize) -> usize {
        for i in (0..=index).rev() {
            let item = items[i].clamp(lower_bound - 1, upper_bound);
            if item < target_length {
                return i;
            }
        }
        0
    }

    fn next(items: &[usize], index: usize, target_length: usize, lower_bound: usize, upper_bound: usize) -> usize {
        for i in index..items.len() {
            let item = items[i].clamp(lower_bound - 1, upper_bound);
            if item < target_length {
                return i;
            }
        }
        items.len()
    }

    fn test_with_parameters(items: &[usize], lower_bound: usize, upper_bound: usize) {
        let matrix = Matrix::from_iterator(items.iter().cloned(), items.len(), lower_bound, upper_bound);
        for i in 0..items.len() {
            for target_length in lower_bound..=upper_bound {
                assert_eq!(
                    previous(items, i, target_length, lower_bound, upper_bound),
                    matrix.previous(i, target_length),
                    "previous; i: {}, target_length: {}",
                    i,
                    target_length
                );
                assert_eq!(
                    next(items, i, target_length, lower_bound, upper_bound),
                    matrix.next(i, target_length),
                    "next; i: {}, target_length: {}",
                    i,
                    target_length
                );
            }
        }
    }

    #[test]
    fn pnsv_matrix_all_01() {
        let items: &[usize] = &[
            2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9,
            9, 9, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2,
        ];
        test_with_parameters(items, 4, 6);
        test_with_parameters(items, 5, 7);
        test_with_parameters(items, 2, 7);
    }

    #[test]
    fn pnsv_matrix_all_02() {
        let items: &[usize] = &[
            2, 2, 2, 3, 3, 4, 3, 3, 1, 4, 6, 4, 5, 8, 8, 5, 6, 6, 7, 6, 7, 3, 2, 7, 8, 6, 6, 6, 9, 9, 9, 9,
            9, 2, 9, 9, 3, 8, 8, 8, 4, 4, 7, 7, 6, 5, 6, 6, 8, 8, 2, 5, 4, 4, 8, 8, 9, 3, 3, 3, 3, 2, 2, 2,
        ];
        test_with_parameters(items, 4, 6);
        test_with_parameters(items, 5, 7);
        test_with_parameters(items, 2, 7);
    }
}

