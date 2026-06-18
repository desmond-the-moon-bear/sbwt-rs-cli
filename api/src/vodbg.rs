// Module and submodule contributions by Martin Kostadinov.

#![allow(unused)]
#![allow(clippy::ptr_arg)]

use crate::{ContractLeft, SbwtIndex};
use crate::subsetseq::SubsetSeq;
use pnsv::Pnsv;

/// Module for Previous and Next Smaller Value support.
pub mod pnsv;
pub mod benchmark;

#[derive(Clone, Debug)]
pub struct VoDbg<'a, SS: SubsetSeq + Send + Sync, P: Pnsv> {
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

impl<'a, SS: SubsetSeq + Send + Sync, P: Pnsv> VoDbg<'a, SS, P> {
    // Code adapted from Dbg.
    
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

    /// Climb to a "lower" level of the graph where the strings in the nodes are one shorter.
    pub fn contract_left(&self, node: Node) -> Node {
        assert!(node.k > 0);
        let range = self.pnsv.contract_left(node.start..node.end, node.k - 1);
        Node {
            start: range.start,
            end: range.end,
            k: node.k - 1,
        }
    }

    // pub fn extend_right(&self, node: Node, character: u8) -> Option<Node> {
    //     assert!(node.k < self.sbwt.k());
    //     todo!()
    // }

    // note(mk): idea - use Pnsv trait to find the next and previous value...
    // pub fn extend_left(&self, node: Node, character: u8) -> Option<Node> {
    //     assert!(node.k < self.sbwt.k());
    //     todo!()
    // }

    // Returns the number of outgoing edges from the given node.
    pub fn outdegree(&self, node: Node) -> usize {
        // assert!(!self.dummy_marks[node.id]);
        // self.sbwt.sbwt.subset_size(self.get_suffix_group_start(node.id))
        todo!()
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
        todo!()
    }
    
    // Returns the colex rank of the smallest k-mer (possibly dummy)
    // that has the same suffix of length (k-1) as the given colex position (possibly dummy)
    fn get_suffix_group_start(&self, mut colex: usize) -> usize {
        // while !self.k_minus_1_marks[colex] { // index 0 is always marked so we're good
        //     colex -= 1;
        // }
        // colex
        todo!()
    }

    // For each outgoing edge from the given node, pushes to the output vector a pair
    // (v, c), where v is the target node and c is the edge label.
    pub fn push_out_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>) {
        // assert!(!self.dummy_marks[node.id]);
        //
        // let rep = self.get_suffix_group_start(node.id);
        //
        // for (i, &c) in self.sbwt.alphabet().iter().enumerate() {
        //     if self.sbwt.sbwt.set_contains(rep, i as u8) {
        //         let outnode = Node{id: self.sbwt.lf_step(rep, i)};
        //         output.push((outnode, c));
        //     }
        // }
        todo!()
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
        todo!()
    }

    // Gets the last character of the k-mer string of the given node.
    pub fn get_last_character(&self, node: Node) -> u8 {
        // assert!(!self.dummy_marks[node.id]);
        // self.sbwt.inlabel(node.id).unwrap() // Can unwrap because this is not a dummy node
        todo!()
    }

    // Returns whether the given node has an outgoing edge labeled with `edge_label`.
    pub fn has_outlabel(&self, node: Node, edge_label: u8) -> bool {
        // assert!(!self.dummy_marks[node.id]);
        // let rep = self.get_suffix_group_start(node.id);
        // let c_idx = self.sbwt.char_idx(edge_label) as u8;
        // self.sbwt.sbwt.set_contains(rep, c_idx)
        todo!()
    }

    // Pushes the labels of all outgoing edges from the given node to the output vector.
    pub fn push_outlabels(&self, node: Node, output: &mut Vec<u8>) {
        // assert!(!self.dummy_marks[node.id]);
        // let rep = self.get_suffix_group_start(node.id);
        // self.sbwt.sbwt.append_set_to_buf(rep, output);
        // for c in output.iter_mut() { // Map from 0123 to ACGT
        //     *c = self.sbwt.alphabet()[*c as usize];
        // }
        todo!()
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
        todo!()
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
        todo!()
    }

    pub fn is_dummy_colex_position(&self, pos: usize) -> bool {
        self.dummy_marks[pos]
    }
}
