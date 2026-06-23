//! An important observation for some operations on the VoDbg is that any correct range of k-mers
//! in the SBWT will contain at most one dummy k-mer which is comprised entirely of $ except for
//! the shared suffix for the given range. Another important note is that this k-mer node will
//! always be located at the beginning of the range.

// Module and submodule contributions by Martin Kostadinov.

#![allow(clippy::ptr_arg)]

use crate::{ContractLeft, ExtendRight, SbwtIndex};
use crate::subsetseq::SubsetSeq;
use pnsv::Pnsv;

use sux::bits::{BitVec, BitFieldVec};
use sux::traits::{BitVecOpsMut, Rank};
use value_traits::slices::{SliceByValue, SliceByValueMut};
use sux::rank_sel::Rank9;

/// Module for Previous and Next Smaller Value support.
pub mod pnsv;
pub mod benchmark;

#[derive(Clone, Debug)]
pub struct VoDbg<'a, SS: SubsetSeq + Send + Sync, P: Pnsv + Send + Sync> {
    sbwt: &'a SbwtIndex<SS>,
    pnsv: &'a P,
    dummy_marks: Rank9,
    dummy_lengths: BitFieldVec,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, Hash)]
pub struct Node {
    pub start: usize,
    pub end: usize,
    pub k: usize,
}

impl<'a, SS: SubsetSeq + Send + Sync, P: Pnsv + Send + Sync> VoDbg<'a, SS, P> {
    /// Marks the dummy k-mers and records the number of $ each has.
    pub fn compute_auxiliary_data_about_dummies(sbwt: &SbwtIndex<SS>) -> (Rank9, BitFieldVec) {
        // Node, depth.
        let mut dfs_stack = Vec::<(usize, usize)>::new(); 
        let mut outlabels = Vec::<u8>::new();

        let mut dummy_count = 0;
        let mut dummy_marks = BitVec::new(sbwt.n_sets());

        // First pass to calculate the dummies.
        dfs_stack.push((0, 0)); // Colex rank of $, depth of $
        while let Some((node, depth)) = dfs_stack.pop() { 
            if !dummy_marks[node] {
                dummy_count += 1;
            }

            dummy_marks.set(node, true);

            if depth + 1 < sbwt.k() {
                outlabels.clear();
                sbwt.sbwt.append_set_to_buf(node, &mut outlabels);
                for &c_idx in outlabels.iter() {
                    let u = sbwt.lf_step(node, c_idx as usize);
                    dfs_stack.push((u, depth + 1));
                }
            }
        }
        
        let dummy_marks = Rank9::new(dummy_marks);

        let bit_width = (64 - sbwt.k().leading_zeros()) as usize;
        let mut dummy_lengths = BitFieldVec::new(bit_width, dummy_count);

        // Second pass to calculate their depths now that we have their positions.
        dfs_stack.push((0, 0)); // Colex rank of $, depth of $
        while let Some((node, depth)) = dfs_stack.pop() { 
            let dummy_index = dummy_marks.rank(node);
            
            // note(mk): I had to introduce another dependency just to set a value at a given index
            // in this packed array... This is extremely disappointing and perhaps I should just
            // use a Vec<usize>.
            dummy_lengths.set_value(dummy_index, depth);

            if depth + 1 < sbwt.k() {
                outlabels.clear();
                sbwt.sbwt.append_set_to_buf(node, &mut outlabels);
                for &c_idx in outlabels.iter() {
                    let u = sbwt.lf_step(node, c_idx as usize);
                    dfs_stack.push((u, depth + 1));
                }
            }
        }

        (dummy_marks, dummy_lengths)
    }

    /// Initializes supports for de Bruijn graph operation based on the given [SbwtIndex].
    /// If the Lcs array of the SBWT is available, it can be given to significantly speed up construction.
    /// IMPORTANT: [select support][SbwtIndex::build_select()] must be built before calling this function. 
    pub fn new(sbwt: &'a SbwtIndex<SS>, pnsv: &'a P) -> Self
    {
        assert!(sbwt.sbwt.has_select_support());
        log::info!("[VoDbg::new] marking dummy nodes...");
        let (dummy_marks, dummy_lengths) = Self::compute_auxiliary_data_about_dummies(sbwt);
        Self {
            sbwt,
            pnsv,
            dummy_marks,
            dummy_lengths,
        }
    }

    /// Push the k-mer string of the node to the given buffer.
    pub fn push_node_kmer(&self, node: Node, buf: &mut Vec<u8>) {
        assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);
        let mut colex_rank = node.start;
        let buf_start = buf.len();
        for _ in 0..node.k {
            match self.sbwt.inlabel(colex_rank) {
                Some(c) => {
                    buf.push(c);
                    // Can unwrap because c != None.
                    colex_rank = self.sbwt.inverse_lf_step(colex_rank).unwrap(); 
                },
                None => {
                    buf.push(b'$');
                },
            }
        }
        buf[buf_start..buf_start+node.k].reverse();
    }

    /// Get the k-mer string label of a node. To avoid memory allocation, check
    /// [VoDbg::push_node_kmer].
    pub fn get_kmer(&self, node: Node) -> Vec<u8> {
        assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);
        let mut buf = Vec::<u8>::with_capacity(node.k);
        self.push_node_kmer(node, &mut buf);
        buf
    }

    /// Get a handle to the node corresponding to the given k-mer, if exists in the graph.
    pub fn get_node(&self, kmer: &[u8]) -> Option<Node> {
        assert!(kmer.len() <= self.sbwt.k());
        self.sbwt.search(kmer).map(|range| Node {
            start: range.start,
            end: range.end,
            k: kmer.len(),
        })
    }

    /// Climb to a "lower" level of the graph where the strings in the nodes are shorter by
    /// removing characters from the left.
    pub fn contract_left(&self, node: Node, target_length: usize) -> Node {
        // note(mk): Think about asserting here...
        // assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);
        assert!(node.k > target_length);
        let range = self.pnsv.contract_left(node.start..node.end, target_length);
        Node {
            start: range.start,
            end: range.end,
            k: target_length,
        }
    }

    /// Climb to a "lower" level of the graph where the strings in the nodes are shorter by
    /// removing characters from the right.
    pub fn contract_right(&self, node: Node, target_length: usize) -> Option<Node> {
        // note(mk): Think about asserting here...
        // assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);
        assert!(node.k > target_length);
        let representative = self.get_representative(node);
        let mut current_length = node.k;
        let mut start = representative;
        while current_length > target_length {
            start = self.sbwt.inverse_lf_step(start)?;
            current_length -= 1;
        }
        let end = self.pnsv.next(start, target_length);
        let node = Node {
            start,
            end,
            k: target_length,
        };
        Some(node)
    }

    // note(mk): idea - use Pnsv trait to find the subdivisions of this node's range...
    // pub fn extend_left(&self, node: Node, character: u8) -> Option<Node> {
    //     assert!(node.k < self.sbwt.k());
    //     todo!()
    // }

    /// Climb to an "upper" level of the graph where the strings in the nodes are one longer.
    pub fn extend_right(&self, node: Node, character: u8) -> Option<Node> {
        let result = self.sbwt.extend_right(node.start..node.end, character);
        let length_increase = if node.k < self.sbwt.k() { 1 } else { 0 };
        if result.is_empty() {
            return None;
        }
        let node = Node {
            start: result.start,
            end: result.end,
            k: node.k + length_increase,
        };
        Some(node)
    }

    /// Returns the number of outgoing edges from the given node.
    pub fn outdegree(&self, node: Node) -> usize {
        assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);
        todo!("Fix this for node.k == sbwt.k");
        // let mut outdegree = 0;
        // for character_index in 0..self.sbwt.alphabet().len() {
        //     let after  = self.sbwt.sbwt.rank(character_index as u8, node.end);
        //     let before = self.sbwt.sbwt.rank(character_index as u8, node.start);
        //     let sets_containing_character_in_range = after - before;
        //     if sets_containing_character_in_range != 0 {
        //         outdegree += 1;
        //     }
        // }
        // outdegree
    }

    /// Returns the number of incoming edges to the given node.
    pub fn indegree(&self, node: Node) -> usize {
        assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);

        //
        // The overall idea is to count the number of ranges for suffixes of length node.k in the
        // range of the suffix which is equal to the (node.k-1) prefix of the node's k-mer.
        //

        let in_neighbours_whole_range = self.contract_right(node, node.k - 1);
        let (mut in_neighbour_start, end) = if let Some(node) = in_neighbours_whole_range {
            (node.start, node.end)
        } else {
            return 0;
        };

        //
        // Let's assume the given node's k-mer is ---ACG i.e. the max k is 6, but the node
        // is part of a lower order de Bruijn graph. After the contract right operation the maximal
        // k-mers in the new range will be of the form ----AC. We want to find the number of slices
        // of the range which have the same (node.k)-suffix. For example, if the range contains:
        //
        // ---AAC
        // ---AAC
        // ---CAC
        // ---GAC
        //
        // then the indegree of the node should be 3. As per the observation about the ranges,
        // there might be a single k-mer at the beginning of the range whose (node.k)-suffix
        // contains at most 1 '$' symbol. We should skip it.
        //

        if self.dummy_marks[in_neighbour_start] {
            let length = self.get_dummy_length(in_neighbour_start);
            if length < node.k {
                in_neighbour_start += 1;
            }
        }

        let mut count = 0;
        let target_length = node.k;
        while in_neighbour_start < end {
            count += 1;
            in_neighbour_start = self.pnsv.next(in_neighbour_start, target_length);
        }

        count
    }
    
    /// The k-mers in the SBWT which are colexicographically the smallest with a given (k-1) suffix
    /// will be referred to as "representatives" i.e. those k-mers in the SBWT which (should) have
    /// a non-empty set of outgoing edges (where the value k in (k-1) is equal to the k of the
    /// SBWT). This method returns the reprsentative for the k-mer in the SBWT at the start of the
    /// range of the given node.
    pub fn get_representative(&self, node: Node) -> usize {
        // For reference I will refer to nodes and k-mers whose length is equal to the k-mers
        // in the SBWT as "maximal" and shorter k-mers and their corresponding nodes as "non-maximal".
        //
        // If the node is non-maximal, then the start of its range is guaranteed to be a
        // representative k-mer. This is because if the start of the range is not the first element
        // with a given (k-1) suffix, then there exists a k-mer immediately before it (in the
        // colexicographic order) with the same (k-1) suffix. This means that this previous k-mer
        // shares any suffix with length less than or equal to (k-1) with the k-mer at the start of
        // the range which means that the range is not correct. Assuming that every node has a
        // correct range, this is a contradiction.
        //
        // If, however, the node is maximal, then its range contains only 1 k-mer which is not
        // guaranteed to be a representative and thus we must find that representative.
        if node.k == self.sbwt.k() {
            // Returns the colex rank of the smallest k-mer (possibly dummy) that has the same suffix of
            // length (k-1) as the given colex position (possibly dummy).
            self.pnsv.previous(node.start, self.sbwt.k() - 1)
        } else {
            node.start
        }
    }

    // For each outgoing edge from the given node, pushes to the output vector a pair
    // (v, c), where v is the target node and c is the edge label.
    pub fn push_out_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>) {
        assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);
        todo!("Fix this for node.k == sbwt.k");
        // for character_index in 0..self.sbwt.alphabet().len() {
        //     let after  = self.sbwt.sbwt.rank(character_index as u8, node.end);
        //     let before = self.sbwt.sbwt.rank(character_index as u8, node.start);
        //     let sets_containing_character_in_range = after - before;
        //     if sets_containing_character_in_range != 0 {
        //         outdegree += 1;
        //     }
        // }
        // let representative = self.get_representative(node);
        // for (i, &c) in self.sbwt.alphabet().iter().enumerate() {
        //     if self.has_outlabel(node, edge_label)
        //     if self.sbwt.sbwt.set_contains(representative, i as u8) {
        //         let outnode = self.follow_outedge(node, c).expect("The given outnode should exist.");
        //         output.push((outnode, c));
        //     }
        // }
    }

    // For each incoming edge to the given node, pushes to the output vector a pair
    // (v, c), where v is the source node and c is the edge label. The edge label
    // will be the same for all in-neighbors because it has to be equal to the
    // last character of the destination k-mer.
    pub fn push_in_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>) {
        assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);
        let inlabel = self.get_last_character(node);
        let in_neighbours_whole_range = self.contract_right(node, node.k - 1);
        let (mut in_neighbour_start, end) = if let Some(node) = in_neighbours_whole_range {
            (node.start, node.end)
        } else {
            return;
        };
        let target_length = node.k;
        loop {
            let in_neighbour_end = self.pnsv.next(in_neighbour_start + 1, target_length);
            let innode = Node {
                start: in_neighbour_start,
                end: in_neighbour_end,
                k: target_length,
            };
            output.push((innode, inlabel));
            if in_neighbour_end >= end {
                break;
            }
            in_neighbour_start = in_neighbour_end;
        }
    }

    // Gets the last character of the k-mer string of the given node.
    pub fn get_last_character(&self, node: Node) -> u8 {
        assert!(0 < node.k && (node.k < self.sbwt.k() || !self.dummy_marks[node.start]));
        // Can unwrap because this is not a dummy node.
        // note(mk): Think about whether the previous statement is true...
        self.sbwt.inlabel(node.start).unwrap() 
    }

    // Returns whether the given node has an outgoing edge labeled with `edge_label`.
    pub fn has_outlabel(&self, node: Node, edge_label: u8) -> bool {
        assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);
        let character_index = self.sbwt.char_idx(edge_label) as u8;
        let after  = self.sbwt.sbwt.rank(character_index, node.end);
        let before = self.sbwt.sbwt.rank(character_index, node.start);
        let sets_containing_character_in_range = after - before;
        sets_containing_character_in_range != 0
    }

    // Pushes the labels of all outgoing edges from the given node to the output vector.
    pub fn push_outlabels(&self, node: Node, output: &mut Vec<u8>) {
        assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);
        for character_index in 0..self.sbwt.alphabet().len() {
            let after  = self.sbwt.sbwt.rank(character_index as u8, node.end);
            let before = self.sbwt.sbwt.rank(character_index as u8, node.start);
            let sets_containing_character_in_range = after - before;
            if sets_containing_character_in_range != 0 {
                let character = self.sbwt.alphabet()[character_index];
                output.push(character);
            }
        }
    }

    // Follows the outgoing edge labeled with edge_label from the given node.
    // Returns None if the edge does not exist.
    pub fn follow_outedge(&self, node: Node, edge_label: u8) -> Option<Node>{
        assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);
        // The SBWT is based on a node-centric de Bruijn graph and thus it is not guaranteed that
        // for each two nodes there is a corresponding (k+1)-mer.
        let contracted = self.contract_left(node, node.k);
        self.extend_right(contracted, edge_label)
    }

    // Follows backward the incoming edge that comes from the i-th smallest k-mer
    // (i ∈ [0, indegree(node)) in colexicographic order that has an outgoing edge to `node`. 
    // Returns None if i ≥ indegree(node). 
    pub fn follow_inedge(&self, node: Node, i: usize) -> Option<Node>{
        assert!(node.k < self.sbwt.k() || !self.dummy_marks[node.start]);
        let in_neighbours_whole_range = self.contract_right(node, node.k - 1);
        let (mut in_neighbour_start, end) = if let Some(node) = in_neighbours_whole_range {
            (node.start, node.end)
        } else {
            return None;
        };
        let target_length = node.k;
        let mut current_index = 0;
        loop {
            let in_neighbour_end = self.pnsv.next(in_neighbour_start + 1, target_length);
            if current_index == i {
                let innode = Node {
                    start: in_neighbour_start,
                    end: in_neighbour_end,
                    k: target_length,
                };
                return Some(innode);
            }
            if in_neighbour_end >= end {
                break;
            }
            current_index += 1;
            in_neighbour_start = in_neighbour_end;
        }
        None
    }

    pub fn is_dummy_colex_position(&self, position: usize) -> bool {
        self.dummy_marks[position]
    }

    pub fn get_dummy_length(&self, position: usize) -> usize {
        let dummy_index = self.dummy_marks.rank(position);
        self.dummy_lengths.get_value(dummy_index)
            .expect("For each dummy there must be a corresponding length of its suffix.")
    }
}
