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
    /// A bitvector which marks the dummy nodes in the SBWT.
    dummy_marks: Rank9,
    /// A packed integer array with the maximum length of a suffix of a dummy node which does not
    /// contain $ characters. These lengths appear in the same order as the marks in the
    /// [VoDbg::dummy_marks] bitvector and together with this packed integer array supports random
    /// access to these values given the position of a dummy node.
    dummy_lengths: BitFieldVec,
}

#[derive(Copy, Clone, Eq, PartialEq, Debug, Hash)]
pub struct Node {
    pub start: usize,
    pub end: usize,
    pub k: usize,
    // Limit the construction of this structure to this module.
    // note(mk): Think about whether this really is a good idea...
    _phantom: std::marker::PhantomData<()>,
}

#[inline]
fn new_node(start: usize, end: usize, k: usize) -> Node {
    Node {
        start,
        end,
        k,
        _phantom: Default::default(),
    }
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
        log::info!("[VoDbg::new] computing auxiliary data about dummy k-mers...");
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
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
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
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        let mut buf = Vec::<u8>::with_capacity(node.k);
        self.push_node_kmer(node, &mut buf);
        buf
    }

    /// Get a handle to the node corresponding to the given k-mer, if exists in the graph.
    pub fn get_node(&self, kmer: &[u8]) -> Option<Node> {
        assert!(kmer.len() <= self.sbwt.k());
        self.sbwt.search(kmer).map(|range| new_node(
            range.start,
            range.end,
            kmer.len(),
        ))
    }

    /// Climb to a "lower" level of the graph where the strings in the nodes are shorter by
    /// removing characters from the left.
    pub fn contract_left(&self, node: Node, target_length: usize) -> Node {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        assert!(node.k > target_length);
        let range = self.pnsv.contract_left(node.start..node.end, target_length);
        new_node(range.start, range.end, target_length)
    }

    /// Climb to a "lower" level of the graph where the strings in the nodes are shorter by
    /// removing characters from the right.
    pub fn contract_right(&self, node: Node, target_length: usize) -> Option<Node> {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        assert!(node.k > target_length);
        let representative = self.get_representative(node);
        let mut current_length = node.k;
        let mut start = representative;
        while current_length > target_length {
            start = self.sbwt.inverse_lf_step(start)?;
            current_length -= 1;
        }
        let end = self.pnsv.next(start + 1, target_length);
        let start = self.pnsv.previous(start, target_length);
        let node = new_node(start, end, target_length);
        Some(node)
    }

    /// Climb to an "upper" level of the graph where the strings in the nodes are one longer. Given
    /// an index returns the corresponding node in colexicographic order if it exists.
    pub fn extend_left_with_index(&self, node: Node, index: usize) -> Option<Node> {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));

        let mut start = node.start;
        // Skip the dummy node which has only $ except for the suffix of the range.
        if self.is_dummy(start) && self.get_dummy_length(start) == node.k {
            start += 1;
        }

        let target_length = node.k + 1;
        let mut end;
        let mut current_index = 0;
        while start < node.end {
            end = self.pnsv.next(start + 1, target_length);
            if current_index == index {
                let extended_node = new_node(start, end, target_length);
                return Some(extended_node);
            }
            start = end;
            current_index += 1;
        }

        None
    }

    /// Climb to an "upper" level of the graph where the strings in the nodes are one longer.
    /// Returns a node only if it exists.
    pub fn extend_left_with_character(&self, node: Node, character: u8) -> Option<Node> {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        let mut kmer_buffer = Vec::<u8>::with_capacity(node.k + 1);

        let mut start = node.start;
        // Skip the dummy node which has only $ except for the suffix of the range.
        if self.is_dummy(start) && self.get_dummy_length(start) == node.k {
            start += 1;
        }

        let target_length = node.k + 1;
        let mut end; 
        while start < node.end {
            end = self.pnsv.next(start + 1, target_length);
            let extended_node = new_node(start, end, target_length);
            kmer_buffer.clear();
            self.push_node_kmer(extended_node, &mut kmer_buffer);
            
            // It is guaranteed that there is at least one character in the suffix.
            if kmer_buffer[0] == character {
                return Some(extended_node);
            }
            start = end;
        }

        None
    }

    /// Climb to an "upper" level of the graph where the strings in the nodes are one longer.
    pub fn extend_right(&self, node: Node, character: u8) -> Option<Node> {
        let result = self.sbwt.extend_right(node.start..node.end, character);
        let length_increase = if node.k < self.sbwt.k() { 1 } else { 0 };
        if result.is_empty() {
            return None;
        }
        let node = new_node(result.start, result.end, node.k + length_increase);
        Some(node)
    }

    /// Returns the number of outgoing edges from the given node.
    pub fn outdegree(&self, node: Node) -> usize {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        // note(mk): We have to do a contract left and an extend right to find whether the given
        // neighbour exists (node centric de Bruijn graph). This is a costly operation. This note
        // is also valid for every other location where information about the out neighbours is
        // needed.
        let mut outdegree = 0;
        for &c in self.sbwt.alphabet() {
            if self.follow_outedge(node, c).is_some() {
                outdegree += 1;
            }
        }
        outdegree
    }

    /// Returns the number of incoming edges to the given node.
    pub fn indegree(&self, node: Node) -> usize {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));

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
        // contains at most 1 $ symbol. We should skip it.
        //

        if self.is_dummy(in_neighbour_start) {
            let length = self.get_dummy_length(in_neighbour_start);
            if length < node.k {
                in_neighbour_start += 1;
            }
        }

        let mut count = 0;
        let target_length = node.k;
        while in_neighbour_start < end {
            count += 1;
            in_neighbour_start = self.pnsv.next(in_neighbour_start + 1, target_length);
        }

        count
    }
 
    /// The k-mers in the SBWT which are colexicographically the smallest with a given (k-1) suffix
    /// will be referred to as "representatives" i.e. those k-mers in the SBWT which (should) have
    /// a non-empty set of outgoing edges (where the value k in (k-1) is equal to the k of the
    /// SBWT). This method returns the reprsentative for the k-mer in the SBWT at the start of the
    /// range of the given node.
    fn get_representative(&self, node: Node) -> usize {
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
            self.get_representative_of_maximal_node(node)
        } else {
            node.start
        }
    }

    /// Returns the colex rank of the smallest k-mer (possibly dummy) that has the same suffix of
    /// length (k-1) as the given colex position (possibly dummy). Serves a similar purpose as
    /// [super::dbg::Dbg::get_suffix_group_start()] from the Dbg structure.
    #[inline]
    fn get_representative_of_maximal_node(&self, node: Node) -> usize {
        self.pnsv.previous(node.start, self.sbwt.k() - 1)
    }

    /// For each outgoing edge from the given node to nodes in the same order de Bruijn graph,
    /// pushes to the output vector a pair (v, c), where v is the target node and c is the edge
    /// label.
    pub fn push_out_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>) {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        for &c in self.sbwt.alphabet() {
            if let Some(neighbor) = self.follow_outedge(node, c) {
                output.push((neighbor, c));
            }
        }
    }

    /// For each incoming edge to the given node in the same order de Bruijn graph, pushes to the
    /// output vector a pair (v, c), where v is the source node and c is the edge label. The edge
    /// label will be the same for all in-neighbors because it has to be equal to the last character
    /// of the destination k-mer.
    pub fn push_in_neighbors(&self, node: Node, output: &mut Vec<(Node, u8)>) {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        let inlabel = self.get_last_character(node);
        let in_neighbours_whole_range = self.contract_right(node, node.k - 1);
        let (mut in_neighbour_start, end) = if let Some(node) = in_neighbours_whole_range {
            (node.start, node.end)
        } else {
            return;
        };

        if self.is_dummy(in_neighbour_start) {
            let length = self.get_dummy_length(in_neighbour_start);
            if length < node.k {
                in_neighbour_start += 1;
            }
        }

        let target_length = node.k;
        while in_neighbour_start < end {
            let in_neighbour_end = self.pnsv.next(in_neighbour_start + 1, target_length);
            let innode = new_node(in_neighbour_start, in_neighbour_end, target_length);
            output.push((innode, inlabel));
            in_neighbour_start = in_neighbour_end;
        }
    }

    /// Gets the last character of the k-mer string of the given node. Panics if the k-mer the node
    /// represents is empty.
    pub fn get_last_character(&self, node: Node) -> u8 {
        assert!(0 < node.k);
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        // If the length of the suffix of the node is greater than 0, then the last character of
        // the k-mer this node represents is not a $. That is, node.start == 0 if and only if
        // node.k == 0.
        self.sbwt.inlabel(node.start).unwrap() 
    }

    /// Returns whether the given node has an outgoing edge labeled with `edge_label` in the same
    /// order de Bruijn graph.
    pub fn has_outlabel(&self, node: Node, edge_label: u8) -> bool {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        self.follow_outedge(node, edge_label).is_some()
    }

    /// Pushes the labels of all outgoing edges in the same order de Bruijn graph from the given
    /// node to the output vector.
    pub fn push_outlabels(&self, node: Node, output: &mut Vec<u8>) {
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        for &c in self.sbwt.alphabet() {
            if self.follow_outedge(node, c).is_some() {
                output.push(c);
            }
        }
    }

    /// Follows the outgoing edge labeled with `edge_label` from the given node in the same order de
    /// Bruijn graph. Returns None if the edge does not exist.
    pub fn follow_outedge(&self, node: Node, edge_label: u8) -> Option<Node>{
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        // The SBWT is based on a node-centric de Bruijn graph and thus it is not guaranteed that
        // for each two nodes there is a corresponding (k+1)-mer.
        let contracted = self.contract_left(node, node.k - 1);
        self.extend_right(contracted, edge_label)
    }

    /// Follows backward the incoming edge that comes from the i-th smallest k-mer
    /// (i ∈ [0, indegree(node)) in colexicographic order that has an outgoing edge to `node` in the
    /// same order de Bruijn graph. Returns None if i ≥ indegree(node).
    pub fn follow_inedge(&self, node: Node, i: usize) -> Option<Node>{
        assert!(node.k < self.sbwt.k() || !self.is_dummy(node.start));
        let in_neighbours_whole_range = self.contract_right(node, node.k - 1);
        let (mut in_neighbour_start, end) = if let Some(node) = in_neighbours_whole_range {
            (node.start, node.end)
        } else {
            return None;
        };

        if self.is_dummy(in_neighbour_start) {
            let length = self.get_dummy_length(in_neighbour_start);
            if length < node.k {
                in_neighbour_start += 1;
            }
        }

        let target_length = node.k;
        let mut current_index = 0;
        loop {
            let in_neighbour_end = self.pnsv.next(in_neighbour_start + 1, target_length);
            if current_index == i {
                let innode = new_node(in_neighbour_start, in_neighbour_end, target_length);
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

    #[inline]
    pub fn is_dummy(&self, position: usize) -> bool {
        self.dummy_marks[position]
    }

    pub fn get_dummy_length(&self, position: usize) -> usize {
        assert!(self.is_dummy(position));
        let dummy_index = self.dummy_marks.rank(position);
        self.dummy_lengths.get_value(dummy_index)
            .expect("For each dummy there must be a corresponding length of its suffix.")
    }
}

#[cfg(test)]
mod tests {
    #![allow(unused)]
    use super::*;
    use crate::{BitPackedKmerSortingMem, LcsArray, SbwtIndexBuilder, SubsetMatrix};
    use crate::dbg::{Dbg, Node as DbgNode};
    use pnsv::PnsvTuned;

    #[test]
    fn all() {
        use rand_chacha::ChaCha20Rng;
        use rand_chacha::rand_core::SeedableRng;
        use rand_chacha::rand_core::RngCore;

        const MIN_K: usize = 3;
        let k: usize = 20;
        let kmer_count = 256;
        let mut rng = ChaCha20Rng::from_seed([42; 32]);

        let mut seqs = Vec::<Vec<u8>>::new();
        for _ in 0..kmer_count {
            let kmer: Vec<u8> = (0..k).map(|_| match rng.next_u32() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }).collect();
            seqs.push(kmer);
            // seqs.push(kmer.clone());
            // seqs_hashset.insert(kmer);
        }

        seqs.sort();
        seqs.dedup();

        let mut sbwt_indices: Vec<(SbwtIndex<SubsetMatrix>, Option<LcsArray>)> = Vec::with_capacity(k);
        let mut graphs = Vec::with_capacity(k);

        for i in MIN_K..=k {
            let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
                .k(i).build_lcs(true)
                .build_select_support(true)
                .run_from_vecs(seqs.as_slice());
            let last = sbwt_indices.len();
            sbwt_indices.push((sbwt, lcs));
        }

        for i in 0..sbwt_indices.len() {
            let dbg = Dbg::new(&sbwt_indices[i].0, sbwt_indices[i].1.as_ref(), 3);
            graphs.push(dbg);
        }

        let vodbg_sbwt = &sbwt_indices[k - MIN_K].0;
        let vodbg_lcs = sbwt_indices[k - MIN_K].1.as_ref().unwrap();
        let pnsv_tuned = PnsvTuned::new_with_default_values(vodbg_sbwt, vodbg_lcs);
        let vodbg = VoDbg::new(vodbg_sbwt, &pnsv_tuned);

        let alphabet = sbwt_indices[k - MIN_K].0.alphabet();

        let mut dbg_buffer = vec![];
        let mut vodbg_buffer = vec![];

        let mut dbg_outlabels = vec![];
        let mut vodbg_outlabels = vec![];

        for current_k in MIN_K..=k {
            let dbg_index = current_k - MIN_K;
            let dbg = &graphs[dbg_index];
            for sequence_index in 0..seqs.len() {
                for sequence_start in 0..k-current_k {
                    let sequence = &seqs[sequence_index][sequence_start..sequence_start + current_k];

                    // get_node
                    let dbg_node = dbg.get_node(sequence).expect("Should exist.");
                    let vodbg_node = vodbg.get_node(sequence).expect("Should exist.");
 
                    // get_kmer
                    assert_eq!(dbg.get_kmer(dbg_node), sequence);
                    assert_eq!(vodbg.get_kmer(vodbg_node), sequence);

                    // indegree, outdegree
                    assert_eq!(dbg.indegree(dbg_node), vodbg.indegree(vodbg_node));
                    assert_eq!(dbg.outdegree(dbg_node), vodbg.outdegree(vodbg_node));

                    // push_in_neighbors
                    dbg_buffer.clear();
                    dbg.push_in_neighbors(dbg_node, &mut dbg_buffer);
                    vodbg_buffer.clear();
                    vodbg.push_in_neighbors(vodbg_node, &mut vodbg_buffer);

                    assert_eq!(dbg_buffer.len(), vodbg_buffer.len());
                    for i in 0..dbg_buffer.len() {
                        assert_eq!(dbg_buffer[i].1, vodbg_buffer[i].1);
                        let dbg_in_neighbor = dbg_buffer[i].0;
                        let vodbg_in_neighbor = vodbg_buffer[i].0;

                        let dbg_in_neighbour_kmer = dbg.get_kmer(dbg_in_neighbor);
                        let vodbg_in_neighbour_kmer = vodbg.get_kmer(vodbg_in_neighbor);
                        assert_eq!(dbg_in_neighbour_kmer, vodbg_in_neighbour_kmer);

                        // follow_inedge
                        let dbg_from_inedge = dbg.follow_inedge(dbg_node, i).expect("Should exist.");
                        let vodbg_from_inedge = vodbg.follow_inedge(vodbg_node, i).expect("Should exist.");
                        assert_eq!(dbg_in_neighbour_kmer, dbg.get_kmer(dbg_from_inedge));
                        assert_eq!(vodbg_in_neighbour_kmer, vodbg.get_kmer(vodbg_from_inedge));
                    }

                    // extend_left_with_character
                    let vodbg_contracted = vodbg.contract_left(vodbg_node, vodbg_node.k - 1);
                    let vodbg_extended = vodbg.extend_left_with_character(vodbg_contracted, sequence[0])
                        .expect("Should exist.");
                    assert_eq!(vodbg_node, vodbg_extended);

                    // extend_left_with_index
                    {
                        let mut index = 0;
                        loop {
                            let extended_op = vodbg.extend_left_with_index(vodbg_node, index);
                            match extended_op {
                                Some(extended) => {
                                    let mut contracted = vodbg.contract_left(extended, vodbg_node.k);
                                    assert_eq!(contracted, vodbg_node);
                                },
                                None => {
                                    break;
                                }
                            }
                            index += 1;
                        }
                    }

                    // push_out_neighbors
                    dbg_buffer.clear();
                    dbg.push_out_neighbors(dbg_node, &mut dbg_buffer);
                    vodbg_buffer.clear();
                    vodbg.push_out_neighbors(vodbg_node, &mut vodbg_buffer);

                    assert_eq!(dbg_buffer.len(), vodbg_buffer.len());
                    for i in 0..dbg_buffer.len() {
                        assert_eq!(dbg_buffer[i].1, vodbg_buffer[i].1);
                        let dbg_out_neighbor = dbg_buffer[i].0;
                        let vodbg_out_neighbor = vodbg_buffer[i].0;
                        assert_eq!(dbg.get_kmer(dbg_out_neighbor), vodbg.get_kmer(vodbg_out_neighbor));
                    }

                    // get_last_character
                    assert_eq!(dbg.get_last_character(dbg_node), vodbg.get_last_character(vodbg_node));

                    // has_outlabel
                    for &c in alphabet {
                        assert_eq!(dbg.has_outlabel(dbg_node, c), vodbg.has_outlabel(vodbg_node, c));
                    }

                    // push_outlabels
                    dbg_outlabels.clear();
                    dbg.push_outlabels(dbg_node, &mut dbg_outlabels);
                    vodbg_outlabels.clear();
                    vodbg.push_outlabels(vodbg_node, &mut vodbg_outlabels);
                    assert_eq!(dbg_outlabels, vodbg_outlabels);
                }
            }
        }
    }
}
