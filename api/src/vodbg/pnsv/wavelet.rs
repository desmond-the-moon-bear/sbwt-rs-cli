use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::ops::{BitVec, Rank, Select, SelectZero};
use simple_sds_sbwt::raw_vector::{PushRaw, RawVector};

use std::collections::vec_deque::VecDeque;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct WindowedWaveletTree {
    pub lower_bound: usize,
    pub window_size: usize,
    pub tree: Vec<Node>,
    pub data: Vec<BitVector>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Node {
    pub parent: usize,
    pub lower_bound: usize,
    pub window_size: usize,
    pub left_child: usize,
    pub right_child: usize,
    pub data: usize,
}

impl Node {
    pub const NULL: usize = 0;
}

impl WindowedWaveletTree {
    pub fn from_iterator<T, I>(input: I, lower_bound: usize, window_size: usize) -> Self
    where
        T: Into<usize>,
        I: Iterator<Item = T> + Clone,
    {
        let mut queue: VecDeque<(usize, usize, usize)> = VecDeque::new();
        let mut tree: Vec<Node> = vec![];
        let mut data_index: usize = 0;

        queue.push_back((0, lower_bound, window_size));

        // Create the nodes based on the lower bound and the window size.
        while !queue.is_empty() {
            let current_node_index = tree.len();
            let (parent, lower_bound, window_size) = queue.pop_front().unwrap();

            let right_child_window_size = window_size >> 1;
            let left_child_window_size = window_size - right_child_window_size;

            let mut left_child = 0;
            let mut right_child = 0;

            // Add 1 for the node that will be pushed now.
            let mut last_node_index = tree.len() + queue.len() + 1;
            if left_child_window_size > 1 {
                left_child = last_node_index;
                last_node_index += 1;
                queue.push_back((current_node_index, lower_bound, left_child_window_size));
            }

            if right_child_window_size > 1 {
                right_child = last_node_index;
                queue.push_back((
                    current_node_index,
                    lower_bound + left_child_window_size,
                    right_child_window_size,
                ));
            }

            let data = data_index;
            data_index += 1;

            tree.push(Node {
                parent,
                lower_bound,
                window_size,
                left_child,
                right_child,
                data,
            });
        }

        // println!("{:#?}", tree);

        // First pass over the items in order to get the number of bits per bitvector. When
        // allocating the bitvectors this count will be used to allocate enough space so that extra
        // allocations from the expansion of the bitvector when the bits are pushed one by one are
        // not needed.
        let upper_bound = lower_bound + window_size - 1;
        let mut bits_per_node = vec![0_usize; tree.len()];
        let mut current_node;
        for item in input.clone() {
            let clamped_item = item.into().clamp(lower_bound, upper_bound);

            current_node = 0;

            loop {
                bits_per_node[current_node] += 1;
                let left_child_window_size =
                    (tree[current_node].window_size & 1) + (tree[current_node].window_size >> 1);

                if clamped_item < tree[current_node].lower_bound + left_child_window_size {
                    // Go to the left child of the current node.
                    current_node = tree[current_node].left_child;
                    if current_node == Node::NULL {
                        break;
                    }
                    continue;
                }

                // Go to the right child of the current node.
                current_node = tree[current_node].right_child;
                if current_node == Node::NULL {
                    break;
                }
            }
        }

        // Allocate the raw vectors with the needed capacity.
        let mut data_raw: Vec<RawVector> = vec![];
        for bit_count in bits_per_node {
            data_raw.push(RawVector::with_capacity(bit_count));
        }

        // Second pass to populate the bits in the bitvectors.
        for item in input.clone() {
            let clamped_item = item.into().clamp(lower_bound, upper_bound);

            current_node = 0;

            loop {
                let left_child_window_size =
                    (tree[current_node].window_size & 1) + (tree[current_node].window_size >> 1);
                let data_index = tree[current_node].data;
                if clamped_item < tree[current_node].lower_bound + left_child_window_size {
                    data_raw[data_index].push_bit(false);
                    current_node = tree[current_node].left_child;
                    if current_node == Node::NULL {
                        break;
                    }
                    continue;
                }

                data_raw[data_index].push_bit(true);
                current_node = tree[current_node].right_child;
                if current_node == Node::NULL {
                    break;
                }
            }
        }

        // Transform the raw vectors into bitvectors and enable rank and select support.
        let mut data: Vec<BitVector> = Vec::with_capacity(data_raw.len());
        for raw_vector in data_raw {
            let mut bit_vector: BitVector = raw_vector.into();
            bit_vector.enable_rank();
            bit_vector.enable_select();
            bit_vector.enable_select_zero();
            data.push(bit_vector);
        }

        Self {
            lower_bound,
            window_size,
            tree,
            data,
        }
    }

    pub fn absolute_index(&self, mut node: usize, mut index_in_node: usize) -> usize {
        let mut parent;
        while node != 0 {
            parent = self.tree[node].parent;
            let parent_data_index = self.tree[parent].data;
            index_in_node = if self.tree[parent].left_child == node {
                // If the node is the left chlid of its parent, then the bit representing the
                // current index in the parent has a value of 0.
                match self.data[parent_data_index].select_zero(index_in_node) {
                    Some(value) => value,
                    None => {
                        return self.data[0].len();
                    }
                }
            } else {
                // Otherwise it has a value of 1.
                match self.data[parent_data_index].select(index_in_node) {
                    Some(value) => value,
                    None => {
                        return self.data[0].len();
                    }
                }
            };
            node = parent;
        }
        index_in_node
    }

    // Find the previous smaller value
    pub fn previous(&self, mut index: usize, target_length: usize) -> usize {
        if target_length <= self.lower_bound {
            return 0;
        }

        let mut result = 0;
        let mut current_node = 0;

        loop {
            if self.tree[current_node].lower_bound >= target_length {
                break;
            }

            let max_value = self.tree[current_node].lower_bound + self.tree[current_node].window_size - 1;
            if max_value < target_length {
                result = result.max(self.absolute_index(current_node, index));
                break;
            }

            let left_child_window_size = (self.tree[current_node].window_size & 1)
                + (self.tree[current_node].window_size >> 1);
            let min_value_in_right_child =
                self.tree[current_node].lower_bound + left_child_window_size;
            let max_value_in_left_child = min_value_in_right_child - 1;

            let data_index = self.tree[current_node].data;
            let is_one = self.data[data_index].get(index);
            if !is_one && max_value_in_left_child < target_length {
                result = result.max(self.absolute_index(current_node, index));
            }

            if min_value_in_right_child >= target_length {
                // Move to the left child as there are no possible smaller values in the right
                // child.
                current_node = self.tree[current_node].left_child;
                if current_node == Node::NULL {
                    break;
                }
                // Update the index so that it is valid in the left child.
                index = self.data[data_index].rank_zero(index);
                if is_one {
                    if index == 0 {
                        // There are no more zeroes to the left. Therefore, there are no more smaller
                        // values that can be found.
                        break;
                    }
                    index -= 1;
                }
                continue;
            }

            // Try to move to the right child if there are any elements which are a part of it.
            // current_node = self.tree[current_node].right_child;
            // if current_node == Node::NULL {
            //     break;
            // }
            // // Update the index so that it is valid in the right child.
            // index = self.data[data_index].rank(index); - if !is_one { 1 } else { 0 };
            todo!()
        }

        result
    }

    pub fn next(&self, mut index: usize, target_value: usize) -> usize {
        if target_value <= self.lower_bound {
            return self.len();
        }

        let mut result = self.len();
        let mut current_node = 0;

        loop {

            todo!();

            // let data_index = self.tree[current_node].data;
            // let left_child_window_size = (self.tree[current_node].window_size & 1)
            //     + (self.tree[current_node].window_size >> 1);
            // let max_value_in_left_child =
            //     self.tree[current_node].lower_bound + left_child_window_size - 1;
            //
            // if !self.data[data_index].get(index) { 
            //     if max_value_in_left_child < target_value {
            //         result = result.min(self.absolute_index(current_node, index));
            //     }
            //
            //     // max_value_in_left_child >= target_value
            //     // There are no elements smaller than the target_value in the right child. Proceed
            //     // to the left child.
            //     current_node = self.tree[current_node].left_child;
            //     if current_node == Node::NULL {
            //         break;
            //     }
            //     index = self.data[data_index].rank_zero(index);
            //     continue;
            // }
            //
            // // At the current index the bit is set to 1.
            // let number_of_smaller_elements_before = self.data[data_index].rank_zero(index);
            // if max_value_in_left_child < target_value {
            //     if number_of_smaller_elements_before < self.data[data_index].count_zeros() {
            //         // There is at least one smaller element afterwards.
            //         let position_of_leftmost_smaller_element = self.data[data_index]
            //             .select_zero(number_of_smaller_elements_before)
            //             .unwrap();
            //         result = result.min(self.absolute_index(current_node, position_of_leftmost_smaller_element));
            //     }
            //
            //     // Move to the right child to search for further candidates.
            //     if max_value_in_left_child + 1 < target_value {
            //         // If there could exists a value smaller than the target value in the range of
            //         // the right child, we have to check.
            //         current_node = self.tree[current_node].right_child;
            //         if current_node == Node::NULL {
            //             break;
            //         }
            //         index = self.data[data_index].rank(index);
            //     } else {
            //         // max_value_in_left_child + 1 == min_value_in_right_child >= target_value
            //         //
            //         // Otherwise, all values in the right child are greater than or equal to the
            //         // target value so there is no need to explore the tree further.
            //         break;
            //     }
            // } else {
            //     // max_value_in_left_child == min_value_in_right_child - 1 >= target_value
            //     // There are no possible values that are smaller than the target value in the right
            //     // child, so the only possibility is the left child.
            //
            //     if number_of_smaller_elements_before == self.data[data_index].count_zeros() {
            //         // There aren't any possible smaller elements from the left child to the right
            //         // of this index, therefore, we can stop the search.
            //         break;
            //     }
            //
            //     // Move to the left child.
            //     current_node = self.tree[current_node].right_child;
            //     if current_node == 0 {
            //         break;
            //     }
            //     index = number_of_smaller_elements_before - 1;
            // }
        }

        result
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.data[0].len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.data[0].is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn previous(items: &[usize], index: usize, target_length: usize, lower_bound: usize) -> usize {
        for i in (0..=index).rev() {
            if target_length <= lower_bound && items[i] < lower_bound {
                continue;
            }
            if items[i] < target_length {
                return i;
            }
        }
        0
    }

    fn next(items: &[usize], index: usize, target_length: usize, lower_bound: usize) -> usize {
        for i in index..items.len() {
            if target_length <= lower_bound && items[i] < lower_bound {
                continue;
            }
            if items[i] < target_length {
                return i;
            }
        }
        items.len()
    }

    #[test]
    fn windowed_wavelet_tree_all() {
        let items: &[usize] = &[
            4, 4, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9,
            9, 9, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4,
        ];
        let lower_bound = 5;
        let wavelet = WindowedWaveletTree::from_iterator(items.iter().cloned(), lower_bound, 3);
        for i in 0..items.len() {
            for target_length in 4..=9 {
                assert_eq!(
                    previous(items, i, target_length, lower_bound),
                    wavelet.previous(i, target_length),
                    "ws: 3, previous; i: {}, target_length: {}",
                    i,
                    target_length
                );
                assert_eq!(
                    next(items, i, target_length, lower_bound),
                    wavelet.next(i, target_length),
                    "ws: 3, next; i: {}, target_length: {}",
                    i,
                    target_length
                );
            }
        }
        let wavelet = WindowedWaveletTree::from_iterator(items.iter().cloned(), 5, 4);
        for i in 0..items.len() {
            for target_length in 4..=9 {
                assert_eq!(
                    previous(items, i, target_length, lower_bound),
                    wavelet.previous(i, target_length),
                    "ws: 4, previous; i: {}, target_length: {}",
                    i,
                    target_length
                );
                assert_eq!(
                    next(items, i, target_length, lower_bound),
                    wavelet.next(i, target_length),
                    "ws: 4, next; i: {}, target_length: {}",
                    i,
                    target_length
                );
            }
        }
    }
}
