
pub mod preprocessing;
pub mod input_structures;

use crate::vodbg::count::Counts;
use crate::{LcsArray, SbwtIndex, SubsetSeq};
use std::io::Read;

use input_structures::{Bwt, Lcp};
use bitvec::vec::BitVec;
use simple_sds_sbwt::bit_vector::BitVector;
use simple_sds_sbwt::int_vector::IntVector;
use simple_sds_sbwt::ops::{BitVec as BitVecTrait, Push, Rank, Select};
use simple_sds_sbwt::raw_vector::{AccessRaw, RawVector};

const FULL_SET: u8 = 0b00001111;

pub struct Output<SS: SubsetSeq + Send> {
    pub sbwt: SbwtIndex<SS>,
    pub lcs: Option<LcsArray>,
    pub counts: Option<Counts>,
}

type Row = BitVec<u64>;

pub fn from_input_build_without_redundant_dummies<RLcp, RBwt, SS>(
    bwt_input: &mut RBwt,
    lcp_input: &mut RLcp,
    length: usize,
    k: usize,
    build_lcs: bool,
) -> std::io::Result<Output<SS>>
where
    RLcp: Read,
    RBwt: Read,
    SS: SubsetSeq + Send
{
    let bwt = preprocessing::ascii_to_bwt(bwt_input, length)?;
    let lcp = preprocessing::truncate_lcp::<_, false>(lcp_input, length, k)?;
    let result = if build_lcs {
        build_without_redundant_dummies::<_, true>(&bwt, &lcp, k)
    } else {
        build_without_redundant_dummies::<_, false>(&bwt, &lcp, k)
    };
    Ok(result)
}

pub fn build_without_redundant_dummies<SS: SubsetSeq + Send, const BUILD_LCS: bool>(
    bwt: &Bwt,
    lcp: &Lcp,
    k: usize
) -> Output<SS> {
    let aux = build_full_auxiliary_bitvectors(bwt, lcp, k);
    let dummy_marks = build_dummy_marks(bwt, k, &aux);

    #[inline]
    fn push_set(rows: &mut [Row], set: u8) {
        for i in 0..4 {
            rows[i].push(set & (1 << i) != 0);
        }
    }

    #[inline]
    fn include_letter(bwt: &Bwt, index: usize, current_set: u8) -> u8 {
        (1 << (bwt.get_char_index(index) - 1)) as u8 | current_set
    }

    let mut rows = Vec::<BitVec<u64>>::new();
    for _ in 0..4 {
        rows.push(BitVec::with_capacity(bwt.len()));
    }

    let bit_width = (usize::BITS - k.leading_zeros()) as usize;
    let mut lcs: Option<IntVector> = if BUILD_LCS {
        let mut value = IntVector::with_capacity(bwt.len(), bit_width).unwrap();
        value.push(0); // '$...$' dummy k-mer
        Some(value)
    } else {
        None
    };

    let mut current_set: u8 = 0;
    let separator_count = bwt.counts[1];

    for index in 1..separator_count {
        if dummy_marks.keep_dummy.bit(index) {
            current_set = include_letter(bwt, index, current_set);
            if current_set == FULL_SET {
                break;
            }
        }
    }
    push_set(&mut rows, current_set);
    log::info!("[build_sets_and_lcs] done with $ range");

    current_set = 0;
    let mut current_lcs_value  = k - 1;
    let mut include_dummy_kmer = false;
    let mut has_dummy_kmer     = false;
    let mut k_range_count = 0;
    for index in separator_count..bwt.len() {
        if aux.k_minus_one_ranges.get(index) {
            if has_dummy_kmer && !include_dummy_kmer {
                k_range_count -= 1;
            }
            while k_range_count > 0 {
                push_set(&mut rows, current_set);
                current_set = 0;
                k_range_count -= 1;
            }

            current_set = 0;
            has_dummy_kmer = false;
            include_dummy_kmer = false;
            k_range_count = 0;
        }

        current_lcs_value = current_lcs_value.min(lcp.get(index));

        let is_start_of_k_range = aux.k_ranges.bit(index);
        if is_start_of_k_range {
            k_range_count += 1;
        }

        if aux.shorter_than_k.get(index) {
            has_dummy_kmer = true;
            if dummy_marks.keep_dummy.bit(index) {
                if !include_dummy_kmer {
                    if BUILD_LCS {
                        lcs.as_mut().unwrap().push(current_lcs_value as u64);
                    }
                    current_lcs_value = k - 1;
                }
                include_dummy_kmer = true;
                current_set = include_letter(bwt, index, current_set);
            }
            if dummy_marks.keep_outedge.bit(index) {
                current_set = include_letter(bwt, index, current_set);
            }
        } else {
            if is_start_of_k_range {
                if BUILD_LCS {
                    lcs.as_mut().unwrap().push(current_lcs_value as u64);
                }
                current_lcs_value = k - 1;
            }
            current_set = include_letter(bwt, index, current_set);
        }
    }

    if has_dummy_kmer && !include_dummy_kmer {
        k_range_count -= 1;
    }
    while k_range_count > 0 {
        push_set(&mut rows, current_set);
        current_set = 0;
        k_range_count -= 1;
    }
    log::info!("[build_sets_and_lcs] done with other ranges");

    let C: Vec<usize> = crate::util::get_C_array(&rows);
    let mut subset_rank = SS::new_from_bit_vectors(rows);
    subset_rank.build_rank();
    let n_sets  = subset_rank.len();
    let n_kmers = aux.kmer_count;
    let mut index = SbwtIndex::<SS>::from_components(
        subset_rank,
        n_kmers,
        k,
        C,
        crate::PrefixLookupTable::new_empty(n_sets)
    );
    let prefix_lookup_table = crate::PrefixLookupTable::new(&index, 8);
    index.set_lookup_table(prefix_lookup_table);
    let lcs = lcs.map(LcsArray::new);

    Output {
        sbwt: index,
        lcs,
        counts: None,
    }
}

struct FullAuxiliaryBitVectors {
    kmer_count: usize,
    shorter_than_k: BitVector,
    equal_to_k: RawVector,
    k_minus_one_ranges: BitVector,
    k_ranges: RawVector,
}

fn build_full_auxiliary_bitvectors(bwt: &Bwt, lcp: &Lcp, k: usize) -> FullAuxiliaryBitVectors {
    // The prefix of a suffix up to the first '$' will be referred to as the true prefix and
    // its length as the true length.
    //
    // An N-range is a contiguous range of suffixes which have the same true prefix after being
    // truncated from the right to a length of N (or if the true length of the suffix is less than
    // N, they are padded with imaginary '$').
    //
    // A (k-1)-range which contains suffixes with true lengths less than k-1 will be referred
    // to as a small range. A (k-1)-range which contains suffixes with true length equal to k
    // will be referred to as a big range.
    //
    // A big range can be further divided into k-ranges.

    log::info!("[build_full_auxiliary_bitvectors] begin");
    let len = bwt.len();
    let mut shorter_than_k     = RawVector::with_len(len, false);
    let mut equal_to_k         = RawVector::with_len(len, false);
    let mut k_minus_one_ranges = RawVector::with_len(len, false);
    let mut k_ranges           = RawVector::with_len(len, false);
    let mut kmer_count = 0;
    let mut order = 1;
    let mut current_length = 0;
    for _ in 1..len {
        let (next_order, character) = bwt.lf_step(order);
        order = next_order;
        if character == b'$' {
            current_length = 0;
            shorter_than_k.set_bit(order, true);
        } else {
            current_length += 1;
            if current_length < k {
                shorter_than_k.set_bit(order, true);
            } else if current_length == k {
                equal_to_k.set_bit(order, true);
            } else {
                current_length = k;
            }
        }
        let lcp_value = lcp.get(order);
        if lcp_value < current_length {
            k_ranges.set_bit(order, true);
            if current_length >= k {
                kmer_count += 1;
            }
            if current_length < k || lcp_value < k - 1 {
                k_minus_one_ranges.set_bit(order, true);
            }
        }
    }

    // Skip the '#' region at the beginning.
    k_minus_one_ranges.set_bit(1, true);
    k_ranges.set_bit(1, true);

    log::info!("[build_full_auxiliary_bitvectors] rank for shorter than k k-mers bitvector");
    let mut shorter_than_k = BitVector::from(shorter_than_k);
    shorter_than_k.enable_rank();
    log::info!("[build_full_auxiliary_bitvectors] rank and select for (k-1)-ranges bitvector");
    let mut k_minus_one_ranges = BitVector::from(k_minus_one_ranges);
    k_minus_one_ranges.enable_rank();
    k_minus_one_ranges.enable_select();
    FullAuxiliaryBitVectors {
        kmer_count,
        shorter_than_k,
        equal_to_k,
        k_minus_one_ranges,
        k_ranges
    }
}

struct DummyMarks {
    keep_dummy: RawVector,
    keep_outedge: RawVector,
}

fn build_dummy_marks(bwt: &Bwt, k: usize, aux: &FullAuxiliaryBitVectors) -> DummyMarks {
    log::info!("[build_dummy_marks] begin");
    let len = bwt.len();
    let mut keep_dummy   = RawVector::with_len(len, false);
    let mut keep_outedge = RawVector::with_len(len, false);

    let start = bwt.counts[1];
    let mut predecessor_confirmed = false;
    for index in start..bwt.len() {
        if aux.k_ranges.bit(index) {
            predecessor_confirmed = false;
        }

        if aux.equal_to_k.bit(index) {
            let predecessor = bwt.inverse_lf_step(index);
            // If we haven't found a full k-mer as a predecessor for this k-range, search for it.
            if !predecessor_confirmed {
                predecessor_confirmed |= has_full_kmer_predecessor(
                    predecessor, bwt, &aux.k_minus_one_ranges, &aux.shorter_than_k
                );
            }

            if predecessor_confirmed {
                keep_outedge.set_bit(predecessor, true);
            } else {
                keep_predecessors(predecessor, bwt, k, &mut keep_dummy);
                predecessor_confirmed = true;
            }
        } else if !aux.shorter_than_k.get(index) {
            // If the true length of the prefix of this suffix is not equal to k and it is not
            // shorter than k, this means that it is longer than k. If this is the case, this means
            // that this k-range has a predecessor.
            predecessor_confirmed = true;
        }
    }

    DummyMarks {
        keep_dummy,
        keep_outedge,
    }
}

fn has_full_kmer_predecessor(
    predecessor: usize,
    bwt: &Bwt,
    k_minus_one_ranges: &BitVector, 
    shorter_than_k: &BitVector
) -> bool {
    let range_start = predecessor;
    let one_index = k_minus_one_ranges.rank(range_start + 1);
    let range_end = if one_index == k_minus_one_ranges.count_ones() {
        bwt.len()
    } else {
        // There is at least one 1 after the current position.
        k_minus_one_ranges.select(one_index).unwrap()
    };
    let range_length = range_end - range_start;
    let number_of_prefixes_with_true_length_smaller_than_k =
        shorter_than_k.rank(range_end) - shorter_than_k.rank(range_start);
    number_of_prefixes_with_true_length_smaller_than_k < range_length
}

fn keep_predecessors(mut predecessor: usize, bwt: &Bwt, mut k: usize, keep_suffix: &mut RawVector) {
    while k > 0 {
        keep_suffix.set_bit(predecessor, true);
        predecessor = bwt.inverse_lf_step(predecessor);
        k -= 1;
    }
}

struct PartialAuxiliaryBitVectors {
    kmer_count: usize,
    k_ranges: RawVector,
}

fn build_parital_auxiliary_bitvectors(bwt: &Bwt, lcp: &Lcp, k: usize) -> PartialAuxiliaryBitVectors {
    log::info!("[build_partial_auxiliary_bitvectors] begin");
    let len = bwt.len();
    let mut k_ranges = RawVector::with_len(len, false);
    let mut kmer_count = 0;
    let mut order = 1;
    let mut current_length = 0;
    for _ in 1..len {
        let (next_order, character) = bwt.lf_step(order);
        order = next_order;
        if character == b'$' {
            current_length = 0;
        } else {
            current_length += 1;
            if current_length > k {
                current_length = k;
            }
        }
        let lcp_value = lcp.get(order);
        if lcp_value < current_length {
            k_ranges.set_bit(order, true);
            if current_length >= k {
                kmer_count += 1;
            }
        }
    }

    // Skip the '#' region at the beginning.
    k_ranges.set_bit(1, true);

    PartialAuxiliaryBitVectors {
        kmer_count,
        k_ranges
    }
}

