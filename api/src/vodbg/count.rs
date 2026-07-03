use crate::{ContractLeft, ExtendRight, SeqStream, StreamingIndex};

use super::DummyInfo;

pub struct Counts {
    pub individual_counts: Vec<u8>,
    pub sample_distance: usize,
    pub sampled_counts: Vec<u64>,
    pub large_counts_up_to_sample: Vec<usize>,
    pub large_counts: Vec<u64>,
}

impl Counts {
    pub const DEFAULT_SAMPLE_DISTANCE: usize = 16;

    pub fn try_new_with_default_values<SS, E, C>(
        sequence_stream: SS,
        streaming_index: StreamingIndex<'_, E, C>,
        dummy_info: &impl DummyInfo
    ) -> Option<Self>
    where 
        SS: SeqStream + Send,
        E: ExtendRight,
        C: ContractLeft,
    {
        Self::try_new(sequence_stream, streaming_index, dummy_info, Self::DEFAULT_SAMPLE_DISTANCE)
    }

    pub fn try_new<SS, E, C>(
        mut sequence_stream: SS,
        streaming_index: StreamingIndex<'_, E, C>,
        dummy_info: &impl DummyInfo,
        sample_distance: usize
    ) -> Option<Self>
    where 
        SS: SeqStream + Send,
        E: ExtendRight,
        C: ContractLeft,
    {
        let mut individual_counts: Vec<u8> = vec![0; streaming_index.n];
        let number_of_samples = streaming_index.n / sample_distance + 1;
        let mut large_counts_up_to_sample: Vec<usize> = vec![0; number_of_samples];

        // note(mk): Think about whether this is efficient enough...
        let mut large_counts = std::collections::BTreeMap::<usize, u64>::new();

        while let Some(sequence) = sequence_stream.stream_next() {
            for (length, range) in streaming_index.matching_statistics_iter(sequence) {
                if length == 0 {
                    continue;
                }

                let representative = range.start;

                let sbwt_was_not_built_with_all_dummies = length < streaming_index.k
                    && (!dummy_info.is_dummy(representative) || dummy_info.get_dummy_length(representative) != length);
                if sbwt_was_not_built_with_all_dummies {
                    return None;
                }

                let sample = representative / sample_distance + 1;
                if individual_counts[representative] == u8::MAX - 1 {
                    large_counts_up_to_sample[sample] += 1;
                    large_counts.insert(representative, 0);
                }

                if individual_counts[representative] == u8::MAX {
                    let extra_count = large_counts.entry(representative).or_default();
                    *extra_count += 1;
                } else {
                    individual_counts[representative] += 1;
                }

            }
        }

        // The values are in the BTreeMap are sorted by key i.e. by index of the sum in the array.
        // Using the array large_counts_before_sample and while scanning we can find the position
        // of the extra sum for the corresponding item in the array without the need of the key.
        let large_counts: Vec<u64> = large_counts.into_values().collect();

        let mut individual_index = 0;
        let mut sampled_counts: Vec<u64> = vec![0; number_of_samples];
        let mut large_count_index = 0;
        // The first sample is "before" the beginning of the array and will have a value of 0.
        for i in 1..number_of_samples {
            sampled_counts[i] = sampled_counts[i - 1];
            large_counts_up_to_sample[i] += large_counts_up_to_sample[i - 1];
            for _ in 0..sample_distance {
                sampled_counts[i] += individual_counts[individual_index] as u64;
                if individual_counts[individual_index] == u8::MAX {
                    sampled_counts[i] += large_counts[large_count_index];
                    large_count_index += 1;
                }
                individual_index += 1;
            }
        }

        let result = Self {
            individual_counts,
            sample_distance,
            sampled_counts,
            large_counts_up_to_sample,
            large_counts,
        };

        Some(result)
    }

    pub fn range_sum(&self, start: usize, end: usize) -> u64 {
        self.prefix_sum(end) - self.prefix_sum(start)
    }

    /// Returns the prefix sum given an end index of an individual count which is not included in
    /// that sum.
    pub fn prefix_sum(&self, end: usize) -> u64 {
        let previous_sample = end / self.sample_distance;
        let mut scan_index = previous_sample * self.sample_distance;
        let mut sum = self.sampled_counts[previous_sample];
        let mut large_count_index = self.large_counts_up_to_sample[previous_sample];
        while scan_index < end {
            sum += self.individual_counts[scan_index] as u64;
            if self.individual_counts[scan_index] == u8::MAX {
                sum += self.large_counts[large_count_index];
                large_count_index += 1;
            }
            scan_index += 1;
        }
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BitPackedKmerSortingMem, SbwtIndexBuilder};
    use crate::vodbg;
    use crate::vodbg::pnsv::PnsvTuned;

    #[test]
    fn sbwt_was_not_built_with_all_dummies() {
        let max_k: usize = 4;
        let seqs: Vec<Vec<u8>> = vec![
            b"AAAAAAAA".to_vec(),
            b"ACACACAC".to_vec(),
            b"ACGTACGT".to_vec()
        ];

        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(max_k).build_lcs(true)
            .build_select_support(true)
            .run_from_vecs(&seqs);
        let lcs = lcs.unwrap();

        let pnsv_tuned = PnsvTuned::new_with_default_values(&sbwt, &lcs, max_k);

        let sequence_stream = crate::util::VecSeqStream::new(&seqs);
        let streaming_index = StreamingIndex {
            extend_right: &sbwt,
            contract_left: &pnsv_tuned,
            n: sbwt.n_sets(),
            k: max_k,
        };

        let vodbg = vodbg::VoDbg::new(&sbwt, &pnsv_tuned);
        let counts = Counts::try_new_with_default_values(sequence_stream, streaming_index, &vodbg);
        assert!(counts.is_none())
    }

    #[test]
    fn randomised_kmers() {
        use rand_chacha::ChaCha20Rng;
        use rand_chacha::rand_core::SeedableRng;
        use rand_chacha::rand_core::RngCore;

        let max_k: usize = 16;
        let kmer_count = 256;
        let mut rng = ChaCha20Rng::from_seed([35; 32]);

        let mut seqs = Vec::<Vec<u8>>::new();
        for _ in 0..kmer_count {
            let kmer: Vec<u8> = (0..max_k).map(|_| match rng.next_u32() % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            }).collect();
            seqs.push(kmer);
        }

        seqs.push(b"AAAAAAAAAAAAAAAAAA!AAAAAAAAAAA".to_vec());
        seqs.push(b"AAAAAAAA!CCCCCCCCCCCCC!AAAAAAA".to_vec());
        seqs.push(vec![b'A'; 1024]);

        seqs.sort();
        seqs.dedup();

        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(max_k).build_lcs(true)
            .add_all_dummy_paths(true)
            .build_select_support(true)
            .run_from_vecs(&seqs);
        let lcs = lcs.unwrap();
        
        let pnsv_tuned = PnsvTuned::new_with_default_values(&sbwt, &lcs, max_k);

        let sequence_stream = crate::util::VecSeqStream::new(&seqs);
        let streaming_index = StreamingIndex {
            extend_right: &sbwt,
            contract_left: &pnsv_tuned,
            n: sbwt.n_sets(),
            k: max_k,
        };

        let vodbg = vodbg::VoDbg::new(&sbwt, &pnsv_tuned);
        let counts = Counts::try_new_with_default_values(sequence_stream, streaming_index, &vodbg).unwrap();

        for current_k in 1..=max_k {
            for node in vodbg::iter::node_iterator_with_k(&vodbg, current_k) {
                let kmer = vodbg.get_kmer(node);
                let true_count = count(&seqs, &kmer);
                let range_count = counts.range_sum(node.start, node.end);
                assert_eq!(true_count, range_count);
            }
        }
    }

    fn count(input: &[Vec<u8>], sequence: &[u8]) -> u64 {
        let mut count = 0;
        for input_sequence in input {
            for window in input_sequence.windows(sequence.len()) {
                if window == sequence {
                    count += 1;
                }
            }
        }
        count
    }
}

