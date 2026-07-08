use crate::{ContractLeft, ExtendRight, SeqStream, StreamingIndex};

use super::DummyInfo;

use std::sync::atomic::{AtomicU8, AtomicU64, AtomicUsize};
use crossbeam::channel::{Sender, Receiver};
use dashmap::DashMap;

#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct Sample {
    pub count: u64,
    pub large_counts_up_to_sample: usize,
}

struct AtomicSample {
    count: AtomicU64,
    large_counts_up_to_sample: AtomicUsize,
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct CountsSingleThreaded {
    pub individual_counts: Vec<u8>,
    pub sample_distance: usize,
    // Store the sample count and the number of large counts up to the sample interleaved to
    // reduce the number of cache misses.
    pub sample_information: Vec<Sample>,
    pub large_counts: Vec<u64>,
}

impl CountsSingleThreaded {
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
        #[allow(unused)] dummy_info: &impl DummyInfo,
        sample_distance: usize
    ) -> Option<Self>
    where 
        SS: SeqStream + Send,
        E: ExtendRight,
        C: ContractLeft,
    {
        let mut individual_counts: Vec<u8> = vec![0; streaming_index.n];
        let sample_count = streaming_index.n / sample_distance + 1;
        let mut sample_information: Vec<Sample> = vec![Sample::default(); sample_count];

        // note(mk): Check other hash maps and/or think for other solutions...
        // let mut large_counts = std::collections::HashMap::<usize, u64>::new();
        let mut large_counts = ahash::AHashMap::<usize, u64>::new();

        let sequence_step = 5000;
        let mut sequence_index: usize = 0;
        let mut progress = 0;

        while let Some(sequence) = sequence_stream.stream_next() {
            for (length, range) in streaming_index.matching_statistics_iter(sequence) {
                if length == 0 {
                    continue;
                }

                let representative = range.start;

                // let sbwt_was_not_built_with_all_dummies = length < streaming_index.k
                //     && (!dummy_info.is_dummy(representative) || dummy_info.get_dummy_length(representative) != length);
                // if sbwt_was_not_built_with_all_dummies {
                //     return None;
                // }

                let sample = representative / sample_distance + 1;
                if individual_counts[representative] == u8::MAX - 1 {
                    sample_information[sample].large_counts_up_to_sample += 1;
                    large_counts.insert(representative, 0);
                }

                if individual_counts[representative] == u8::MAX {
                    let extra_count = large_counts.entry(representative).or_default();
                    *extra_count += 1;
                } else {
                    individual_counts[representative] += 1;
                }
            }
            sequence_index += 1;
            if sequence_index > sequence_step {
                sequence_index = 0;
                progress += sequence_step;
                log::info!("[Counts::new] progress ({}/?)...", progress);
            }
        }

        let mut pairs: Vec<_> = large_counts.into_iter().collect();
        pairs.sort();
        let large_counts: Vec<u64> = pairs.into_iter().map(|(_, count)| count).collect();

        let mut individual_index = 0;
        let mut large_count_index = 0;
        // The first sample is "before" the beginning of the array and will have a value of 0.
        for i in 1..sample_count {
            sample_information[i].count = sample_information[i - 1].count;
            sample_information[i].large_counts_up_to_sample += sample_information[i - 1].large_counts_up_to_sample;
            for _ in 0..sample_distance {
                sample_information[i].count += individual_counts[individual_index] as u64;
                if individual_counts[individual_index] == u8::MAX {
                    sample_information[i].count += large_counts[large_count_index];
                    large_count_index += 1;
                }
                individual_index += 1;
            }
        }

        let result = Self {
            individual_counts,
            sample_distance,
            sample_information,
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
        let mut sum = self.sample_information[previous_sample].count;
        let mut large_count_index = self.sample_information[previous_sample].large_counts_up_to_sample;
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

    // pub fn iter<'a>(&'a self) -> Iter<'a> {
    //     Iter { index: 0, large_count_index: 0, counts: self }
    // }
}

// #[derive(Debug, Clone)]
// pub struct Iter<'a> {
//     index: usize,
//     large_count_index: usize,
//     counts: &'a CountsSingleThreaded,
// }
//
// impl Iterator for Iter<'_> {
//     type Item = u64;
//
//     fn next(&mut self) -> Option<Self::Item> {
//         if self.index >= self.counts.individual_counts.len() {
//             return None;
//         }
//
//         let mut sum = self.counts.individual_counts[self.index] as u64;
//         self.index += 1;
//         if sum == u8::MAX as u64 {
//             sum += self.counts.large_counts[self.large_count_index];
//             self.large_count_index += 1;
//         }
//
//         Some(sum)
//     }
// }

type LargeCountMap = DashMap<usize, u64>;

struct Batch {
    buffer: Vec<u8>,
    bounds: Vec<usize>,
}

impl Batch {
    fn iter<'a>(&'a self) -> BatchIterator<'a> {
        BatchIterator { batch: self, index: 0 }
    }
}

struct BatchIterator<'a> {
    batch: &'a Batch,
    index: usize,
}

impl<'a> Iterator for BatchIterator<'a> {
    type Item = &'a [u8];
    fn next(&mut self) -> Option<Self::Item> {
        if self.index + 1 < self.batch.bounds.len() {
            let start = self.batch.bounds[self.index];
            let end = self.batch.bounds[self.index + 1];
            let seq = &self.batch.buffer[start..end];
            self.index += 1;
            Some(seq)
        } else {
            None
        }
    }
}

fn sequence_producer_thread<SS>(mut sequence_stream: SS, buffer_capacity: usize, output: Sender<Batch>)
where SS: SeqStream + Send,
{
    let mut buffer = Vec::<u8>::with_capacity(buffer_capacity);
    let mut bounds = Vec::<usize>::new();

    while let Some(sequence) = sequence_stream.stream_next() {
        bounds.push(buffer.len());
        buffer.extend(sequence);

        if buffer.len() >= buffer_capacity {
            // End sentinel.
            bounds.push(buffer.len());
            let batch = Batch { buffer, bounds };
            output.send(batch).unwrap();

            buffer = Vec::<u8>::with_capacity(buffer_capacity);
            bounds = Vec::<usize>::new();
        }
    }

    if !buffer.is_empty() {
        bounds.push(buffer.len());
        let batch = Batch { buffer, bounds };
        output.send(batch).unwrap();
    }
}

fn sequence_consumer_thread_with_hash_map<E, C, D>(
    input: Receiver<Batch>,
    streaming_index: StreamingIndex<'_, E, C>,
    dummy_info: &D,
    sample_distance: usize,
    individual_counts: &[AtomicU8],
    sample_info: &[AtomicSample],
    large_counts: &LargeCountMap,
)
where 
    E: ExtendRight + Send + Sync,
    C: ContractLeft + Send + Sync,
    D: DummyInfo + Send + Sync
{
    use std::sync::atomic::Ordering::*;
    let _ = dummy_info;
    while let Ok(batch) = input.recv() {
        for sequence in batch.iter() {
            for (length, range) in streaming_index.matching_statistics_iter(sequence) {
                if length == 0 {
                    continue;
                }

                let representative = range.start;

                // let sbwt_was_not_built_with_all_dummies = length < streaming_index.k
                //     && (!dummy_info.is_dummy(representative) || dummy_info.get_dummy_length(representative) != length);
                // if sbwt_was_not_built_with_all_dummies {
                //     return None;
                // }

                let sample = representative / sample_distance + 1;
                let count = individual_counts[representative].load(Acquire);
                // individual_counts[representative].compare_exchange(current, new, success, failure)
                if count == u8::MAX - 1 {
                    sample_info[sample].large_counts_up_to_sample.fetch_add(1, AcqRel);
                }

                if count == u8::MAX {
                    let mut entry = large_counts.entry(representative).or_default();
                    *entry += 1;
                } else {
                    todo!();
                }
            }
        }
    }
}

pub struct CountsConcurrentHashMap {
    pub individual_counts: Vec<u8>,
    pub sample_distance: usize,
    // Store the sample count and the number of large counts up to the sample interleaved to
    // reduce the number of cache misses.
    pub sample_information: Vec<Sample>,
    pub large_counts: Vec<u64>,
}

impl CountsConcurrentHashMap {
    #[allow(unused)]
    pub fn try_new<SS, E, C, D>(
        mut sequence_stream: SS,
        streaming_index: StreamingIndex<'_, E, C>,
        dummy_info: &D,
        sample_distance: usize,
        mut thread_count: usize,
    ) -> Option<Self>
    where 
        SS: SeqStream + Send,
        E: ExtendRight + Send + Sync,
        C: ContractLeft + Send + Sync,
        D: DummyInfo + Send + Sync
    {
        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(thread_count).build().unwrap();

        let hash_map: DashMap<usize, u64> = DashMap::new();

        let individual_counts: Vec<AtomicU8> = vec![0u8; streaming_index.n]
            .into_iter()
            .map(AtomicU8::new)
            .collect();

        let sample_count = streaming_index.n / sample_distance + 1;
        let sample_information: Vec<AtomicSample> = vec![Sample::default(); sample_count]
            .into_iter()
            .map(|sample| AtomicSample {
                count: AtomicU64::new(sample.count),
                large_counts_up_to_sample: AtomicUsize::new(sample.large_counts_up_to_sample),
            })
            .collect();

        use crossbeam::channel::bounded;

        thread_pool.scope(|s| {
            let hash_map = &hash_map;

            let consumer_thread_count = (thread_count - 1).max(1);
            let (batches_in, batches_out) = bounded(consumer_thread_count);
            let buffer_capacity = 1_usize << 20;

            s.spawn(move |_| {
                sequence_producer_thread(sequence_stream, buffer_capacity, batches_in);
            });

            todo!();
        });

        // let individual_counts: Vec<u8> = individual_counts.into_iter().map(|i| i.load(Ordering::Relaxed)).collect();

        todo!()
    }
}

#[cfg(test)]
mod tests {
    // use super::*;
    // use crate::{BitPackedKmerSortingMem, SbwtIndexBuilder};
    // use crate::vodbg;
    // use crate::vodbg::pnsv::PnsvTuned;

    // #[test]
    // fn sbwt_was_not_built_with_all_dummies() {
    //     let max_k: usize = 4;
    //     let seqs: Vec<Vec<u8>> = vec![
    //         b"AAAAAAAA".to_vec(),
    //         b"ACACACAC".to_vec(),
    //         b"ACGTACGT".to_vec()
    //     ];
    //
    //     let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
    //         .k(max_k).build_lcs(true)
    //         .build_select_support(true)
    //         .run_from_vecs(&seqs);
    //     let lcs = lcs.unwrap();
    //
    //     let pnsv_tuned = PnsvTuned::new_with_default_values(&sbwt, &lcs, max_k);
    //
    //     let sequence_stream = crate::util::VecSeqStream::new(&seqs);
    //     let streaming_index = StreamingIndex {
    //         extend_right: &sbwt,
    //         contract_left: &pnsv_tuned,
    //         n: sbwt.n_sets(),
    //         k: max_k,
    //     };
    //
    //     let vodbg = vodbg::VoDbg::new(&sbwt, &pnsv_tuned);
    //     let counts = CountsSingleThreaded::try_new_with_default_values(sequence_stream, streaming_index, &vodbg);
    //     assert!(counts.is_none())
    // }
    //
    // #[test]
    // fn randomised_kmers() {
    //     use rand_chacha::ChaCha20Rng;
    //     use rand_chacha::rand_core::SeedableRng;
    //     use rand_chacha::rand_core::RngCore;
    //
    //     let max_k: usize = 16;
    //     let kmer_count = 256;
    //     let mut rng = ChaCha20Rng::from_seed([35; 32]);
    //
    //     let mut seqs = Vec::<Vec<u8>>::new();
    //     for _ in 0..kmer_count {
    //         let kmer: Vec<u8> = (0..max_k).map(|_| match rng.next_u32() % 4 {
    //             0 => b'A',
    //             1 => b'C',
    //             2 => b'G',
    //             _ => b'T',
    //         }).collect();
    //         seqs.push(kmer);
    //     }
    //
    //     seqs.push(b"AAAAAAAAAAAAAAAAAA!AAAAAAAAAAA".to_vec());
    //     seqs.push(b"AAAAAAAA!CCCCCCCCCCCCC!AAAAAAA".to_vec());
    //     seqs.push(vec![b'A'; 1024]);
    //     seqs.push(vec![b'C'; 1024]);
    //
    //     seqs.sort();
    //     seqs.dedup();
    //
    //     let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
    //         .k(max_k).build_lcs(true)
    //         .add_all_dummy_paths(true)
    //         .build_select_support(true)
    //         .run_from_vecs(&seqs);
    //     let lcs = lcs.unwrap();
    //
    //     let pnsv_tuned = PnsvTuned::new_with_default_values(&sbwt, &lcs, max_k);
    //
    //     let sequence_stream = crate::util::VecSeqStream::new(&seqs);
    //     let streaming_index = StreamingIndex {
    //         extend_right: &sbwt,
    //         contract_left: &pnsv_tuned,
    //         n: sbwt.n_sets(),
    //         k: max_k,
    //     };
    //
    //     let vodbg = vodbg::VoDbg::new(&sbwt, &pnsv_tuned);
    //     let counts = CountsSingleThreaded::try_new_with_default_values(sequence_stream, streaming_index, &vodbg).unwrap();
    //
    //     for current_k in 1..max_k {
    //         for node in vodbg::iter::node_iterator_with_k(&vodbg, current_k) {
    //             let kmer = vodbg.get_kmer(node);
    //             let true_count = count(&seqs, &kmer);
    //             let range_count = counts.range_sum(node.start, node.end);
    //             assert_eq!(true_count, range_count);
    //         }
    //     }
    //
    //     let mut counts_iter = counts.iter();
    //     let mut count_from_iterator;
    //     for node in vodbg::iter::node_iterator_with_k(&vodbg, max_k) {
    //         let kmer = vodbg.get_kmer(node);
    //         let true_count = count(&seqs, &kmer);
    //         let range_count = counts.range_sum(node.start, node.end);
    //         count_from_iterator = counts_iter.next().unwrap();
    //         while counts_iter.index <= node.start {
    //             count_from_iterator = counts_iter.next().unwrap();
    //         }
    //         assert_eq!(true_count, range_count);
    //         assert_eq!(true_count, count_from_iterator);
    //     }
    // }

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

