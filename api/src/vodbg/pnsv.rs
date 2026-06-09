use balanced_parenthesis as bp;

use crate::{ContractLeft, LcsArray};
use simple_sds_sbwt::{
    raw_vector::{
        AccessRaw,
        RawVector
    }
};

pub mod balanced_parenthesis;

/// Previous and Next Smaller Value using Balanced Parenthesis.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PnsvBp {
    pub previous: bp::Bp,
    pub next: bp::Bp,
    pub count: usize,
}

impl PnsvBp {
    /// Constructs a data structure which supports Previous and Next Smaller Value queries on the
    /// values in the original array.
    pub fn from_iterator<T, I>(values: I, count: usize, block_size: usize) -> Self
    where
        T: Ord + Copy,
        I: Iterator<Item = T> + DoubleEndedIterator + Clone,
    {
        let mut stack = vec![];
        let previous_smaller_bp_vector = Self::make_psv_bp_vector(values.clone(), count, &mut stack);
        let next_smaller_bp_vector = Self::make_psv_bp_vector(values.rev(), count, &mut stack);

        let previous = bp::Bp::new(previous_smaller_bp_vector, block_size);
        let next = bp::Bp::new(next_smaller_bp_vector, block_size);

        Self { previous, next, count }
    }

    /// Given an index in the original array finds the index of the previous smaller value.
    #[inline]
    pub fn previous(&self, index: usize) -> usize {
        let parenthesis_index = self.previous.select(index);
        let parent_parenthesis_index = self.previous.enclose(parenthesis_index);
        self.previous.rank(parent_parenthesis_index)
    }

    /// Given an index in the original array finds the index of the next smaller value.
    #[inline]
    pub fn next(&self, index: usize) -> usize {
        let reverse_index = self.count - 1 - index;
        let parenthesis_index = self.next.select(reverse_index);
        let parent_parenthesis_index = self.next.enclose(parenthesis_index);
        let reverse_result_index = self.next.rank(parent_parenthesis_index);
        self.count - 1 - reverse_result_index
    }

    /// Given an iterator of values constructs the balanced parenthesis representation of the
    /// Previous Smaller Value tree. A PSV tree is such that each vertex's parent is the previous
    /// smaller value in the array.
    pub fn make_psv_bp_vector<T, I>(values: I, count: usize, stack: &mut Vec<T>) -> RawVector
    where
        T: Ord + Copy,
        I: Iterator<Item = T>,
    {
        // The maximum height of the stack is the number of unique values that the array contains.
        // In the case of the Longest Common Suffix array that would be k.
        stack.clear();

        let mut result = RawVector::with_len(count * 2, false);
        let mut bit_iterator = 0;

        for value in values {
            loop {
                if stack.is_empty() {
                    break;
                }

                let previous_value = *stack.last().unwrap();
                if previous_value < value {
                    // The previous value is smaller, therefore it is this value's parent.
                    break;
                }

                // The previous value is either greater than this value in which case it is a sibling
                // or a child of a sibling of this value, or it is equal to this value in which case it
                // is a sibling of this value in the PSV tree.
                stack.pop();

                // By moving the iterator forward we leave a 0 in the resulting bitvector which is
                // equivalent to closing the parenthesis.
                bit_iterator += 1;
            }

            result.set_bit(bit_iterator, true);
            bit_iterator += 1;
            stack.push(value);
        }

        result
    }
}

pub struct LcsPnsvBp {
    pub lcs: LcsArray,
    pub pnsv: PnsvBp,
}

impl LcsPnsvBp {
    pub fn new(lcs: LcsArray, block_size: usize) -> Self {
        let iterator = (0..lcs.len()).map(|index| lcs.access(index));
        let pnsv = PnsvBp::from_iterator(iterator, lcs.len(), block_size);
        Self {
            lcs,
            pnsv,
        }
    }
}

impl ContractLeft for LcsPnsvBp {
    #[allow(non_snake_case)]
    fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {
        let mut new_start = I.start;
        while self.lcs.access(new_start) >= target_len {
            new_start = self.pnsv.previous(new_start);
        }
        let mut new_end = I.end;
        while self.lcs.access(new_end) >= target_len {
            new_end = self.pnsv.next(new_start);
        }
        new_start..new_end
    }
}

#[cfg(test)]
mod tests {
    use super::balanced_parenthesis as bp;
    use super::*;

    fn make_psv_tree<T, I>(values: I, count: usize) -> Vec<usize>
    where
        T: Ord + Copy,
        I: Iterator<Item = T>,
    {
        let mut stack: Vec<(usize, T)> = vec![];

        let mut result = vec![count; count];

        for (index, value) in values.enumerate() {
            loop {
                if stack.is_empty() {
                    break;
                }
                let (previous_index, previous_value) = *stack.last().unwrap();
                if previous_value < value {
                    result[index] = previous_index;
                    break;
                }
                stack.pop();
            }

            stack.push((index, value));
        }

        result
    }

    #[test]
    fn bp_method_random_sequence() {
        use rand::distributions::uniform::Uniform;
        use rand::distributions::Distribution;
        use rand::rngs::StdRng;
        use rand::SeedableRng;

        let n = 50000;
        let k = 31u8;
        let mut rng = StdRng::from_seed([42; 32]);
        let uniform = Uniform::from(0..k);
        let numbers = uniform.sample_iter(&mut rng).take(n).collect::<Vec<_>>();

        let answers = make_psv_tree(numbers.iter(), numbers.len());

        let mut bp_value_stack = vec![];
        let bp_vector = PnsvBp::make_psv_bp_vector(numbers.iter(), numbers.len(), &mut bp_value_stack);

        for block_size in (10..=100).step_by(10) {
            let psv_tree = bp::Bp::new(bp_vector.clone(), block_size);

            for i in 0..n {
                let parenthesis_index = psv_tree.select(i);
                let parenthesis_index_of_parent = psv_tree.enclose(parenthesis_index);
                if answers[i] == numbers.len() {
                    assert!(
                        parenthesis_index_of_parent == parenthesis_index,
                        "Failed at i: {}, parenthesis_index: {}",
                        i,
                        parenthesis_index
                    );
                    println!("{i} -> null");
                    continue;
                }
                let parent_index = psv_tree.rank(parenthesis_index_of_parent);
                assert_eq!(
                    answers[i], parent_index,
                    "Failed at i: {}, parenthesis_index: {}",
                    i, parenthesis_index
                );
                println!("{i} -> {}", parent_index);
            }
        }
    }

    // TODO(mk): write tests whether LcsPnsvBp works correctly.
}
