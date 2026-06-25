
use crate::SubsetSeq;

use super::VoDbg;
use super::pnsv::Pnsv;

use sux::bits::BitFieldVec;
use sux::traits::Rank;
use sux::rank_sel::Rank9;
use value_traits::slices::SliceByValue;

pub fn node_iterator_with_k<'a, SS, P>(vodbg: &'a VoDbg<'a, SS, P>, k: usize) -> NodeIterator<'a, P>
where 
    SS: SubsetSeq + Send + Sync,
    P: Pnsv + Send + Sync,
{
    NodeIterator {
        start: 1,
        end: vodbg.sbwt.n_sets(),
        k,
        pnsv: vodbg.pnsv,
        dummy_marks: &vodbg.dummy_marks,
        dummy_lengths: &vodbg.dummy_lengths,
    }
}

#[derive(Clone, Debug)]
pub struct NodeIterator<'a, P: Pnsv + Send + Sync> {
    start: usize,
    end: usize,
    k: usize,
    pnsv: &'a P,
    dummy_marks: &'a Rank9,
    dummy_lengths: &'a BitFieldVec,
}

impl<P> NodeIterator<'_, P>
where P: Pnsv + Send + Sync
{
    pub fn reset(&mut self, k: usize) {
        self.start = 1;
        self.k = k;
    }
}

impl<P> Iterator for NodeIterator<'_, P>
where P: Pnsv + Send + Sync
{
    type Item = super::Node;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start >= self.end {
            return None;
        }

        let mut dummy_index;
        if self.dummy_marks[self.start] {
            dummy_index = self.dummy_marks.rank(self.start);
            loop {
                let invalid_start =
                    self.start < self.end
                    && self.dummy_marks[self.start]
                    && self.dummy_lengths.get_value(dummy_index).expect("Dummy index should be valid.") < self.k;
                if !invalid_start {
                    break;
                }
                self.start += 1;
                dummy_index += 1;
            }
        }

        if self.start >= self.end {
            return None;
        }
        
        let region_end = self.pnsv.next(self.start + 1, self.k);
        let node = super::new_node(self.start, region_end, self.k);
        self.start = region_end;
        Some(node)
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::vodbg::pnsv::PnsvTuned;
    use crate::{BitPackedKmerSortingMem, LcsArray, SbwtIndex, SbwtIndexBuilder, SubsetMatrix};
    use crate::dbg::Dbg;

    #[test]
    pub fn randomised_kmers() {
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

        let mut vodbg_node_iterator = node_iterator_with_k(&vodbg, 0);

        for current_k in MIN_K..=k {
            let dbg_index = current_k - MIN_K;
            let dbg = &graphs[dbg_index];

            vodbg_node_iterator.reset(current_k);
            let mut dbg_node_iterator = dbg.node_iterator();
            
            loop {
                let dbg_node_op = dbg_node_iterator.next();
                let vodbg_node_op = vodbg_node_iterator.next();
                if dbg_node_op.is_none() {
                    assert!(vodbg_node_op.is_none());
                    break;
                }
                let dbg_node = dbg_node_op.unwrap();
                let vodbg_node = vodbg_node_op.expect("Should exist.");

                let dbg_kmer = dbg.get_kmer(dbg_node);
                let vodbg_kmer = vodbg.get_kmer(vodbg_node);
                assert_eq!(dbg_kmer, vodbg_kmer);
            }
        }
    }
}

