// Module and submodule contributions by Martin Kostadinov.

#![allow(unused)]
#![allow(clippy::ptr_arg)]

use crate::{ContractLeft, ExtendRight, SbwtIndex};
use crate::subsetseq::SubsetSeq;
use pnsv::Pnsv;

/// Module for Previous and Next Smaller Value support.
pub mod pnsv;
pub mod benchmark;

#[derive(Clone, Debug)]
pub struct VoDbg<'a, SS: SubsetSeq + Send + Sync, P: Pnsv + Send + Sync> {
    sbwt: &'a SbwtIndex<SS>,
    pnsv: &'a P,
    dummy_marks: bitvec::vec::BitVec,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, Hash)]
pub struct Node {
    pub start: usize,
    pub end: usize,
    pub k: usize,
}

impl<'a, SS: SubsetSeq + Send + Sync, P: Pnsv + Send + Sync> VoDbg<'a, SS, P> {
    // Some code adapted from Dbg.
    
    // An internal function for marking the dummy nodes in the SBWT.
    fn mark_dummies(sbwt: &SbwtIndex<SS>) -> bitvec::vec::BitVec {
        sbwt.compute_dummy_node_marks()
    } 

    // An internal function marking for each (k-1)-mer the smallest k-mer that has that (k-1)-mer as a suffix.
    // fn mark_k_minus_1_mers(lcs: &crate::streaming_index::LcsArray, k: usize) -> bitvec::vec::BitVec {
    //     let mut k_minus_1_marks = bitvec![0; lcs.len()];
    //     for i in 0..lcs.len(){
    //         let len = lcs.access(i);
    //         if len < k - 1 {
    //             k_minus_1_marks.set(i,true);
    //         }
    //     }
    //     k_minus_1_marks
    // }

    // Returns next 1-bit to the right of i, or bv.len() if does not exist
    // fn next_1_bit(bv: &bitvec::vec::BitVec, mut i: usize) -> usize {
    //     while i < bv.len() && !bv[i] {
    //         i += 1;
    //     }
    //     i
    // }

    // Initializes supports for de Bruijn graph operation based on the given [SbwtIndex].
    // If the Lcs array of the SBWT is available, it can be given to significantly speed up construction.
    // IMPORTANT: [select support][SbwtIndex::build_select()] must be built before calling this function. 
    pub fn new(sbwt: &'a SbwtIndex<SS>, pnsv: &'a P) -> Self
    {
        assert!(sbwt.sbwt.has_select_support());
        // let k_minus_1_marks = match lcs {
        //     Some(lcs) => {
        //         log::info!("Building (k-1)-mer marks from LCS array");
        //         Self::mark_k_minus_1_mers(lcs, sbwt.k())
        //     }
        //     None => {
        //         log::info!("No LCS-array given. Building (k-1)-mer marks with column inversion.");
        //         sbwt.mark_k_minus_1_mers(n_threads)
        //     }
        // };

        log::info!("[VoDbg::new] marking dummy nodes...");
        let dummy_marks = Self::mark_dummies(sbwt);
        Self {
            sbwt,
            pnsv,
            dummy_marks
        }
    }

    // Push the k-mer string of the node to the given buffer.
    pub fn push_node_kmer(&self, node: Node, buf: &mut Vec<u8>) {
        // assert!(!self.dummy_marks[node.start]);
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

    // Get the k-mer string label of a node. To avoid memory allocation, check
    // [VoDbg::push_node_kmer].
    pub fn get_kmer(&self, node: Node) -> Vec<u8> {
        // assert!(!self.dummy_marks[node.id]);
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
        // assert!(node.k < self.sbwt.k());
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

    // Returns the number of outgoing edges from the given node.
    pub fn outdegree(&self, node: Node) -> usize {
        // assert!(!self.dummy_marks[node.id]);
        // self.sbwt.sbwt.subset_size(self.get_suffix_group_start(node.id))
        let representative = self.get_representative(node);
        self.sbwt.sbwt.subset_size(representative)
    }

    // Returns the number of incoming edges to the given node.
    pub fn indegree(&self, node: Node) -> usize {
        // assert!(!self.dummy_marks[node.id]);
        // match self.follow_inedge(node, 0) {
        //     Some(v) => {
        //         let s = v.id; // This is the smallest non-dummy in-neighbor
        //         let e = Self::next_1_bit(&self.k_minus_1_marks, s + 1);
        //         e - s
        //     },
        //     None => 0,
        // }

        // The overall idea is to count the number of ranges for suffixes of length node.k in the
        // range of the suffix which is equal to the (node.k-1) prefix of the node's k-mer.
        let in_neighbours_whole_range = self.contract_right(node, node.k - 1);
        let (mut in_neighbour_start, end) = if let Some(node) = in_neighbours_whole_range {
            (node.start + 1, node.end)
        } else {
            return 0;
        };
        let mut count = 1;
        let target_length = node.k;
        loop {
            in_neighbour_start = self.pnsv.next(in_neighbour_start + 1, target_length);
            if in_neighbour_start >= end {
                break;
            }
            count += 1;
        }
        count
    }
    
    // Returns the colex rank of the smallest k-mer (possibly dummy) that has the same suffix of
    // length (k-1) as the given colex position (possibly dummy).
    fn get_suffix_group_start(&self, colex: usize) -> usize {
        // while !self.k_minus_1_marks[colex] { // index 0 is always marked so we're good
        //     colex -= 1;
        // }
        // colex
        self.pnsv.previous(colex, self.sbwt.k() - 1)
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
            self.get_suffix_group_start(node.start)
        } else {
            node.start
        }
    }

    // For each outgoing edge from the given node, pushes to the output vector a pair
    // (v, c), where v is the target node and c is the edge label.
    pub fn push_out_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>) {
        // assert!(!self.dummy_marks[node.id]);
        // let rep = self.get_suffix_group_start(node.id);
        // for (i, &c) in self.sbwt.alphabet().iter().enumerate() {
        //     if self.sbwt.sbwt.set_contains(rep, i as u8) {
        //         let outnode = Node{id: self.sbwt.lf_step(rep, i)};
        //         output.push((outnode, c));
        //     }
        // }
        let representative = self.get_representative(node);
        for (i, &c) in self.sbwt.alphabet().iter().enumerate() {
            if self.sbwt.sbwt.set_contains(representative, i as u8) {
                let outnode = self.follow_outedge(node, c).expect("The given outnode should exist.");
                output.push((outnode, c));
            }
        }
    }

    // For each incoming edge to the given node, pushes to the output vector a pair
    // (v, c), where v is the source node and c is the edge label. The edge label
    // will be the same for all in-neighbors because it has to be equal to the
    // last character of the destination k-mer.
    pub fn push_in_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>) {
        // assert!(!self.dummy_marks[node.id]);
        // let inlabel = self.get_last_character(node);
        // if let Some(v) = self.sbwt.inverse_lf_step(node.id) { // Predecessor
        //     let vrep = self.get_suffix_group_start(v);
        //     let end = Self::next_1_bit(&self.k_minus_1_marks, vrep+1);
        //     (vrep..end).filter(|&i| !self.dummy_marks[i]).for_each(|i|{
        //         output.push((Node{id: i}, inlabel));
        //     });
        // }

        // The overall idea is similar to the indegree method i.e. to subdivide the range of the
        // suffix which is equal to the (node.k-1) prefix of the node's k-mer into ranges for nodes
        // of length (node.k).
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
        assert!(!self.dummy_marks[node.start]);
        self.sbwt.inlabel(node.start).unwrap() // Can unwrap because this is not a dummy node
    }

    // Returns whether the given node has an outgoing edge labeled with `edge_label`.
    pub fn has_outlabel(&self, node: Node, edge_label: u8) -> bool {
        // assert!(!self.dummy_marks[node.id]);
        let representative = self.get_representative(node);
        let c_idx = self.sbwt.char_idx(edge_label) as u8;
        self.sbwt.sbwt.set_contains(representative, c_idx)
    }

    // Pushes the labels of all outgoing edges from the given node to the output vector.
    pub fn push_outlabels(&self, node: Node, output: &mut Vec<u8>) {
        // assert!(!self.dummy_marks[node.id]);
        let representative = self.get_representative(node);
        self.sbwt.sbwt.append_set_to_buf(representative, output);
        for c in output.iter_mut() { // Map from 0123 to ACGT
            *c = self.sbwt.alphabet()[*c as usize];
        }
    }

    // Follows the outgoing edge labeled with edge_label from the given node.
    // Returns None if the edge does not exist.
    pub fn follow_outedge(&self, node: Node, edge_label: u8) -> Option<Node>{
        // assert!(!self.dummy_marks[node.id]);
        // if !self.has_outlabel(node, edge_label) {
        //     return None;
        // }
        // let rep = self.get_suffix_group_start(node.id);
        // Some(Node{id: self.sbwt.lf_step(rep, self.sbwt.char_idx(edge_label))})
        let extended = self.extend_right(node, edge_label)?;
        let node = self.contract_left(extended, node.k);
        Some(node)
    }

    // Follows backward the incoming edge that comes from the i-th smallest k-mer
    // (i ∈ [0, indegree(node)) in colexicographic order that has an outgoing edge to `node`. 
    // Returns None if i ≥ indegree(node). 
    pub fn follow_inedge(&self, node: Node, i: usize) -> Option<Node>{
        // assert!(!self.dummy_marks[node.id]);
        // let v = self.sbwt.inverse_lf_step(node.id)?;
        //
        // // Inverse lf step always takes us to a suffix group start
        // let vrep = v;
        //
        // let end = Self::next_1_bit(&self.k_minus_1_marks, vrep+1);
        // let mut non_dummies = 0_usize;
        //
        // // Return the position of the 0-bit with rank in_edge_number in the range
        // for j in vrep..end {
        //     if !self.dummy_marks[j] {
        //         if non_dummies == i {
        //             return Some(Node{id: j});
        //         }
        //         non_dummies += 1; 
        //     }
        // }
        //
        // None
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

    pub fn is_dummy_colex_position(&self, pos: usize) -> bool {
        self.dummy_marks[pos]
    }
}
