use balanced_parenthesis as bp;

use crate::ContractLeft;
use simple_sds_sbwt::raw_vector::{AccessRaw, RawVector};

pub mod balanced_parenthesis;

/// Previous and Next Smaller Value using Balanced Parenthesis.
pub struct PnsvBp {
    pub previous: bp::Bp,
    pub next: bp::Bp,
}

impl PnsvBp {
    pub fn new<T, I>(values: I, count: usize, block_size: usize) -> Self
    where
        T: Ord + Copy,
        I: Iterator<Item = T> + DoubleEndedIterator + Clone,
    {
        let mut stack = vec![];
        let previous_smaller_bp_vector = Self::make_psv_bp_vector(values.clone(), count, &mut stack);
        let next_smaller_bp_vector = Self::make_psv_bp_vector(values.rev(), count, &mut stack);

        let previous = bp::Bp::new(previous_smaller_bp_vector, block_size);
        let next = bp::Bp::new(next_smaller_bp_vector, block_size);

        Self { previous, next }
    }

    pub fn make_psv_bp_vector<T, I>(values: I, count: usize, stack: &mut Vec<T>) -> RawVector
    where
        T: Ord + Copy,
        I: Iterator<Item = T>,
    {
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
                bit_iterator += 1;
            }

            result.set_bit(bit_iterator, true);
            bit_iterator += 1;
            stack.push(value);
        }

        result
    }
}

impl ContractLeft for PnsvBp {
    #[allow(non_snake_case)]
    #[allow(unused)]
    fn contract_left(&self, I: std::ops::Range<usize>, target_len: usize) -> std::ops::Range<usize> {
        todo!()
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

        let psv_tree = bp::Bp::new(bp_vector, 20);

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
