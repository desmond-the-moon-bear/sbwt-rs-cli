use crate::{ContractLeft, ExtendRight, SeqStream, StreamingIndex};

use super::DummyInfo;

use std::sync::atomic::{AtomicBool, AtomicU8, AtomicUsize};
use crossbeam::channel::{Sender, Receiver};
use dashmap::DashMap;

type LargeCountMap = DashMap<usize, u64>;

#[derive(Debug, Default, Copy, Clone, PartialEq, Eq)]
pub struct Sample {
    pub count: u64,
    pub large_counts_up_to_sample: usize,
}

struct PartialAtomicSample {
    count: u64,
    large_counts_up_to_sample: AtomicUsize,
}

#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct Counts {
    pub individual_counts: Vec<u8>,
    pub sample_distance: usize,
    // Store the sample count and the number of large counts up to the sample interleaved to
    // reduce the number of cache misses.
    pub sample_information: Vec<Sample>,
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

        Self::process_sample_information(
            sample_distance,
            sample_count, 
            &individual_counts,
            &large_counts,
            &mut sample_information
        );

        let result = Self {
            individual_counts,
            sample_distance,
            sample_information,
            large_counts,
        };

        Some(result)
    }

    /// Counts the number of occurrences of each k-mer in the input sequence concurrently. There
    /// must always be one sequence producer thread and at least one sequence consumer thread i.e.
    /// the procedure will panic if thread_count is < 1.
    pub fn try_new_concurrent<SS, E, C, D>(
        sequence_stream: SS,
        streaming_index: &StreamingIndex<'_, E, C>,
        dummy_info: &D,
        sample_distance: usize,
        additional_memory_bound_gb: usize,
        thread_count: usize,
    ) -> Option<Self>
    where 
        SS: SeqStream + Send,
        E: ExtendRight + Sync,
        C: ContractLeft + Sync,
        D: DummyInfo + Sync
    {
        assert!(thread_count > 1);

        let individual_counts: Vec<AtomicU8> = vec![0u8; streaming_index.n]
            .into_iter()
            .map(AtomicU8::new)
            .collect();

        let sample_count = streaming_index.n / sample_distance + 1;
        let sample_information: Vec<PartialAtomicSample> = vec![Sample::default(); sample_count]
            .into_iter()
            .map(|sample| PartialAtomicSample {
                count: sample.count,
                large_counts_up_to_sample: AtomicUsize::new(sample.large_counts_up_to_sample),
            })
            .collect();

        let large_counts = LargeCountMap::new();

        let success = AtomicBool::new(true);
        let consumer_thread_count = thread_count.saturating_sub(1).max(1);

        use crossbeam::channel::bounded;
        // To ensure that the additional memory bound is respected, bound the channel size.
        let (batches_in, batches_out) = bounded(consumer_thread_count);

        // Divide the memory evenly between the producer thread and the consumer threads.
        let buffer_capacity = additional_memory_bound_gb * ((1_usize << 30) / thread_count);
        std::thread::scope(|s| {
            s.spawn(move || {
                sequence_producer_thread(sequence_stream, buffer_capacity, batches_in);
            });
        });

        let thread_pool = rayon::ThreadPoolBuilder::new().num_threads(consumer_thread_count).build().unwrap();
        thread_pool.scope(|s| {

            let individual_counts = &individual_counts;
            let sample_information = &sample_information;
            let large_counts = &large_counts;
            let success = &success;
            for _ in 0..consumer_thread_count {
                let batches_out_cloned = batches_out.clone();
                s.spawn(move |_| {
                    let result = sequence_consumer_thread_with_hash_map(
                        batches_out_cloned,
                        streaming_index,
                        dummy_info,
                        sample_distance,
                        individual_counts,
                        sample_information,
                        large_counts
                    );
                    if result.is_ok() {
                        success.store(false, Ordering::SeqCst);
                    }
                });
            }
        });

        if !success.load(Ordering::SeqCst) {
            // return None;
        }

        let mut pairs: Vec<_> = large_counts.into_iter().collect();
        pairs.sort();
        let large_counts: Vec<u64> = pairs.into_iter().map(|(_, count)| count).collect();

        // The following transformations from the atomic should ultimately become NOPs.
        use std::sync::atomic::Ordering;
        let individual_counts: Vec<u8> = individual_counts
            .into_iter()
            .map(|i| i.load(Ordering::Relaxed))
            .collect();
        let mut sample_information: Vec<Sample> = sample_information
            .into_iter()
            .map(|sample| Sample {
                count: sample.count,
                large_counts_up_to_sample: sample.large_counts_up_to_sample.load(Ordering::Relaxed)
            }).collect();

        Self::process_sample_information(
            sample_distance,
            sample_count, 
            &individual_counts,
            &large_counts,
            &mut sample_information
        );

        let result = Self {
            individual_counts,
            sample_distance,
            sample_information,
            large_counts,
        };
        Some(result)
    }

    fn process_sample_information(
        sample_distance: usize,
        sample_count: usize,
        individual_counts: &[u8],
        large_counts: &[u64],
        sample_information: &mut [Sample]
    ) {
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
    }

    /// Returns the sum of the counts in a given range.
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

    pub fn iter<'a>(&'a self) -> Iter<'a> {
        Iter { index: 0, large_count_index: 0, counts: self }
    }
}

#[derive(Debug, Clone)]
pub struct Iter<'a> {
    index: usize,
    large_count_index: usize,
    counts: &'a Counts,
}

impl Iterator for Iter<'_> {
    type Item = u64;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.counts.individual_counts.len() {
            return None;
        }

        let mut sum = self.counts.individual_counts[self.index] as u64;
        self.index += 1;
        if sum == u8::MAX as u64 {
            sum += self.counts.large_counts[self.large_count_index];
            self.large_count_index += 1;
        }

        Some(sum)
    }
}

/// A batch of input sequences. Idea borrowed from
/// [`crate::bitpacked_kmer_sorting::kmer_splitter`].
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

    let mut progress: f64 = 0.0;

    while let Some(sequence) = sequence_stream.stream_next() {
        bounds.push(buffer.len());
        buffer.extend(sequence);

        if buffer.len() >= buffer_capacity {
            progress += buffer.len() as f64;
            log::info!("[reading sequencess] {} done.", human_bytes::human_bytes(progress));

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

    drop(output);
}

fn sequence_consumer_thread_with_hash_map<E, C, D>(
    input: Receiver<Batch>,
    streaming_index: &StreamingIndex<'_, E, C>,
    dummy_info: &D,
    sample_distance: usize,
    individual_counts: &[AtomicU8],
    sample_information: &[PartialAtomicSample],
    large_counts: &LargeCountMap,
) -> Result<(), ()>
where 
    E: ExtendRight,
    C: ContractLeft,
    D: DummyInfo
{
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
                //     return Err(());
                // }

                // note(mk):
                //  * fetch_update will be renamed to try_update in a future version of Rust.
                //  * the ordering could be more relaxed probably?
                use std::sync::atomic::Ordering;
                let result = individual_counts[representative].fetch_update(Ordering::SeqCst, Ordering::SeqCst, |value| {
                    if value < u8::MAX {
                        Some(value + 1)
                    } else {
                        None
                    }
                });
                let previous = match result {
                    Ok(value) => value,
                    Err(value) => value,
                };

                let sample = representative / sample_distance + 1;
                if previous == u8::MAX - 1 {
                    sample_information[sample].large_counts_up_to_sample.fetch_add(1, Ordering::AcqRel);
                    // It is possible that another thread gets to execute right here and to update
                    // the count of the same representative. Using insert(representative, 0) would
                    // overwrite that count.
                    large_counts.entry(representative).or_default();
                }

                if previous == u8::MAX {
                    let mut entry = large_counts.entry(representative).or_default();
                    *entry += 1;
                }
            }
        }
    }
    Ok(())
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
        seqs.push(vec![b'C'; 1024]);

        seqs.sort();
        seqs.dedup();

        let (sbwt, lcs) = SbwtIndexBuilder::<BitPackedKmerSortingMem>::new()
            .k(max_k).build_lcs(true)
            .add_all_dummy_paths(true)
            .build_select_support(true)
            .run_from_vecs(&seqs);
        let lcs = lcs.unwrap();

        let pnsv_tuned = PnsvTuned::new_with_default_values(&sbwt, &lcs, max_k);

        let vodbg = vodbg::VoDbg::new(&sbwt, &pnsv_tuned);

        let streaming_index = StreamingIndex {
            extend_right: &sbwt,
            contract_left: &pnsv_tuned,
            n: sbwt.n_sets(),
            k: max_k,
        };

        // let counts = Counts::try_new_with_default_values(sequence_stream, streaming_index, &vodbg).unwrap();
        let sequence_stream = crate::util::VecSeqStream::new(&seqs);
        let counts = Counts::try_new_concurrent(
            sequence_stream,
            &streaming_index,
            &vodbg,
            Counts::DEFAULT_SAMPLE_DISTANCE, 1, 4
        ).unwrap();

        for current_k in 1..max_k {
            for node in vodbg::iter::node_iterator_with_k(&vodbg, current_k) {
                let kmer = vodbg.get_kmer(node);
                let true_count = count(&seqs, &kmer);
                let range_count_many_threads  = counts.range_sum(node.start, node.end);
                assert_eq!(true_count, range_count_many_threads);
            }
        }

        let mut counts_iter = counts.iter();
        let mut count_from_iterator;
        for node in vodbg::iter::node_iterator_with_k(&vodbg, max_k) {
            let kmer = vodbg.get_kmer(node);
            let true_count = count(&seqs, &kmer);
            let range_count  = counts.range_sum(node.start, node.end);

            count_from_iterator = counts_iter.next().unwrap();
            while counts_iter.index <= node.start {
                count_from_iterator = counts_iter.next().unwrap();
            }

            assert_eq!(true_count, range_count);
            assert_eq!(true_count, count_from_iterator);
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

