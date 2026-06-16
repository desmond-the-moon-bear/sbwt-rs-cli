// Code by Martin Kostadinov.

#![allow(unused)]

use crate::streaming_index::{ContractLeft, ExtendRight, StreamingIndex};
use crate::SbwtIndexVariant;
use crate::streaming_index::LcsArray;

use super::pnsv::Pnsv;

// note(mk): Temporary solution for reading input data.
pub fn read_index_and_lcs(arguments_start: usize) -> (SbwtIndexVariant, LcsArray) {
    let mut args = std::env::args().skip(arguments_start);
    let sbwt_path = args.next().expect("expected sbwt index path");
    let lcs_path = args.next().expect("expected lcs index path");

    let mut index_reader = std::io::BufReader::new(std::fs::File::open(sbwt_path).unwrap());
    let index = crate::load_sbwt_index_variant(&mut index_reader).unwrap();

    let mut lcs_reader = std::io::BufReader::new(std::fs::File::open(lcs_path).unwrap());
    let lcs = LcsArray::load(&mut lcs_reader).unwrap();

    (index, lcs)
}

pub fn read_query(arguments_start: usize) -> Vec<Vec<u8>> {
    use crate::SeqStream;

    let mut args = std::env::args().skip(arguments_start);
    let query_path = args.next().expect("expected query path");
    let query_file = std::fs::File::open(query_path).unwrap();
    let buf_reader = std::io::BufReader::new(query_file);
    let mut stream = crate::JSeqIOSeqStreamWrapper {
        inner: jseqio::reader::DynamicFastXReader::new(buf_reader).unwrap(),
    };
    let mut result: Vec<Vec<u8>> = vec![];
    while let Some(sequence) = stream.stream_next() {
        result.push(sequence.into());
    }
    result
}

pub fn benchmark_bms_separate_queries<E, C>(
    index: &StreamingIndex<'_, E, C>,
    queries: &[Vec<u8>],
    bound: usize,
) where
    E: ExtendRight,
    C: ContractLeft,
{
    let start_time = std::time::Instant::now();
    let mut checksum = 0_usize;
    let mut n_kmers_queried = 0_usize;
    for query in queries {
        for x in index.bounded_matching_statistics(query, bound).iter() {
            checksum += x.0;
        }
        n_kmers_queried += std::cmp::max(query.len() as isize - index.k as isize + 1, 0) as usize;
    }
    let end_time = std::time::Instant::now();
    print!("{:.2},", (end_time - start_time).as_secs_f64());
    let nanos_per_kmer = (end_time - start_time).as_nanos() as f64 / n_kmers_queried as f64;
    print!("{:.2},{}", nanos_per_kmer, checksum);
    // nanos_per_kmer.round() as usize
}

// fn benchmark_bms_joined_queries_lcs(bound: usize) {}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vodbg::pnsv::{
        ABS,
        LcsPnsvBp,
        LcsSimd,
        PnsvDyn,
        PnsvHybrid,
        PnsvMatrix,
        Ranges,
        WWT,
    };

    // note(mk): Temporary solution for running the benchmark. Command to run these benchmarks:
    // cargo t (--release) vodbg::test::comparison -- --ignored --nocapture sbwt_path lcs_path query_path
    #[ignore]
    #[test]
    fn comparison() {
        println!("loading data...");
        let (index, lcs) = read_index_and_lcs(4);
        let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
        println!("lcs.len: {}", lcs.len());
        println!("index.n_sets: {}", sbwt.n_sets());
        let queries = read_query(6);

        // println!("creating standard bp structure...");
        // let lcs_pnsv = LcsPnsvBp::new(&lcs, 2048);

        println!("creating hybrid data structure...");
        println!("creating ranges...");
        let ranges = Ranges::new(&sbwt, sbwt.n_sets(), 7);
        let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);
        println!("creating lcs_simd...");
        let lcs_simd = LcsSimd::from_iterator(iterator.clone(), lcs.len());
        // println!("creating wavelet...");
        // let wavelet = WWT::from_iterator(iterator, 7, 2);
        // let pnsv_hybrid = PnsvHybrid {
        //     ranges,
        //     wavelet,
        //     lcs_simd,
        // };

        println!("creating matrix...");
        let matrix = PnsvMatrix::from_iterator(iterator, lcs.len(), 8, 10);
        
        let pnsv_dyn = PnsvDyn {
            structures: [&ranges, &matrix, &lcs_simd]
        };

        // let lcs_index = StreamingIndex {
        //     extend_right: &sbwt,
        //     contract_left: &lcs,
        //     n: sbwt.n_sets(),
        //     k: sbwt.k(),
        // };
        //
        // let lcs_pnsv_bp_index = StreamingIndex {
        //     extend_right: &sbwt,
        //     contract_left: &lcs_pnsv,
        //     n: sbwt.n_sets(),
        //     k: sbwt.k(),
        // };
        // let pnsv_hybrid_index = StreamingIndex {
        //     extend_right: &sbwt,
        //     contract_left: &pnsv_hybrid,
        //     n: sbwt.n_sets(),
        //     k: sbwt.k(),
        // };
        let pnsv_dyn_index = StreamingIndex {
            extend_right: &sbwt,
            contract_left: &pnsv_dyn,
            n: sbwt.n_sets(),
            k: sbwt.k(),
        };

        println!("running benchmarks...");

        let lower = 8;
        let upper = 31;

        for bound in 1..upper {
            print!("dyn,{},", bound);
            benchmark_bms_separate_queries(&pnsv_dyn_index, &queries, bound);
            println!();
        }
        // for bound in 1..upper {
        //     print!("hyb,{},", bound);
        //     benchmark_bms_separate_queries(&pnsv_hybrid_index, &queries, bound);
        //     println!();
        // }
        // for bound in 1..upper {
        //     print!("bp,{},", bound);
        //     benchmark_bms_separate_queries(&lcs_pnsv_bp_index, &queries, bound);
        //     println!();
        // }
        //
        // for bound in lower..upper {
        //     print!("scan,{},", bound);
        //     benchmark_bms_separate_queries(&lcs_index, &queries, bound);
        //     println!();
        // }
    }

    #[ignore]
    #[test]
    fn correctness() {
        println!("loading data...");
        let (index, lcs) = read_index_and_lcs(4);
        let SbwtIndexVariant::SubsetMatrix(sbwt) = index;

        let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);
        let lcs_simd = LcsSimd::from_iterator(iterator.clone(), lcs.len());
        let wavelet = WWT::from_iterator(iterator, 7, 2);

        let ten_percent = lcs.len() / 10;

        for i in 0..lcs.len() {
            for target_length in wavelet.lower_bound + 1..wavelet.lower_bound + wavelet.window_size {
                let lcs_answer = lcs_simd.previous(i, target_length);
                let wavelet_answer = wavelet.previous(i, target_length);
                assert_eq!(lcs_answer, wavelet_answer, "p; i: {}, target_length: {}", i, target_length);
                let lcs_answer = lcs_simd.next(i, target_length);
                let wavelet_answer = wavelet.next(i, target_length);
                assert_eq!(lcs_answer, wavelet_answer, "n; i: {}, target_length: {}", i, target_length);
            }
            if i % ten_percent == ten_percent - 1 {
                println!("{}0%", 1 + i / ten_percent);
            }
        }
    }

    #[ignore]
    #[test]
    fn simd_scan_compare() {
        let (index, lcs) = read_index_and_lcs(4);
        let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
        let queries = read_query(6);

        let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);
        let lcs_simd = LcsSimd::from_iterator(iterator, lcs.len());

        let lcs_index = StreamingIndex {
            extend_right: &sbwt,
            contract_left: &lcs,
            n: sbwt.n_sets(),
            k: sbwt.k(),
        };

        let lcs_simd_index = StreamingIndex {
            extend_right: &sbwt,
            contract_left: &lcs_simd,
            n: sbwt.n_sets(),
            k: sbwt.k(),
        };

        let bound = 10;

        print!("scan,{},", bound);
        benchmark_bms_separate_queries(&lcs_index, &queries, bound);
        println!();

        print!("simd,{},", bound);
        benchmark_bms_separate_queries(&lcs_simd_index, &queries, bound);
        println!();
    }

    fn statistics_lcs_simd(lcs_simd: &LcsSimd, target_length: usize, bound: usize) -> (f64, f64, f64, f64) {
        let n = lcs_simd.n;
        let target_length = target_length as u8;

        let mut successful_previous = 0;
        let start_time = std::time::Instant::now();
        for i in 0..n {
            successful_previous += if lcs_simd.scan_left_bounded(i, target_length, bound).is_ok() {
                1
            } else {
                0
            };
        }
        let end_time = std::time::Instant::now();
        let nanos_per_previous = (end_time - start_time).as_nanos() as f64 / n as f64;
        let percentage_previous = successful_previous as f64 / n as f64;

        let mut successful_next = 0;
        let start_time = std::time::Instant::now();
        for i in 0..n {
            successful_next += if lcs_simd.scan_left_bounded(i, target_length, bound).is_ok() {
                1
            } else {
                0
            };
        }
        let end_time = std::time::Instant::now();
        let nanos_per_next = (end_time - start_time).as_nanos() as f64 / n as f64;
        let percentage_next = successful_next as f64 / n as f64;

        (percentage_previous, percentage_next, nanos_per_next, nanos_per_previous)
    }

    fn statistics_pnsv_matrix(matrix: &PnsvMatrix, target_length: usize) -> (f64, f64) {
        let n = matrix.width;

        let start_time = std::time::Instant::now();
        for i in 0..n {
            let _ = std::hint::black_box(matrix.previous(i, target_length));
        }
        let end_time = std::time::Instant::now();
        let nanos_per_previous = (end_time - start_time).as_nanos() as f64 / n as f64;

        let start_time = std::time::Instant::now();
        for i in 0..n {
            let _ = std::hint::black_box(matrix.next(i, target_length));
        }
        let end_time = std::time::Instant::now();
        let nanos_per_next = (end_time - start_time).as_nanos() as f64 / n as f64;

        (nanos_per_previous, nanos_per_next)
    }

    #[ignore]
    #[test]
    fn simd_bounded_scan_time() {
        let (index, lcs) = read_index_and_lcs(4);
        let SbwtIndexVariant::SubsetMatrix(sbwt) = index;

        let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);

        let lower_bound = 8;
        let upper_bound = 10;

        println!("creating lcs_simd...");
        let lcs_simd = LcsSimd::from_iterator(iterator.clone(), lcs.len());

        println!("creating matrix...");
        let matrix = PnsvMatrix::from_iterator(iterator, lcs.len(), lower_bound, upper_bound);

        println!("timing...");
        let item_bound: usize = 1000;
        let word_bound = item_bound.div_ceil(LcsSimd::LANES);

        for target_length in lower_bound..=upper_bound {
            let (
                percentage_previous,
                percentage_next,
                nanos_per_next_scan,
                nanos_per_previous_scan,
            ) = statistics_lcs_simd(&lcs_simd, target_length, word_bound);

            let (
                nanos_per_previous_matrix,
                nanos_per_next_matrix
            ) = statistics_pnsv_matrix(&matrix, target_length);

            println!("target_length: {}", target_length);
            println!(
                "%previous: {:.3} <> t_scan/t_bitvector: {:.3} ({:.3}/{:.3})",
                percentage_previous,
                nanos_per_previous_scan / nanos_per_previous_matrix,
                nanos_per_previous_scan,
                nanos_per_previous_matrix
            );

            println!(
                "%next: {:.3} <> t_scan/t_bitvector: {:.3} ({:.3}/{:.3})",
                percentage_next,
                nanos_per_next_scan / nanos_per_next_matrix,
                nanos_per_next_scan,
                nanos_per_next_matrix
            );
            println!();
        }
    }

    fn statistics_lcs_simd_with_matrix_fallback(lcs_simd: &LcsSimd, matrix: &PnsvMatrix, target_length: usize, bound: usize) -> (f64, f64) {
        let n = lcs_simd.n;
        let target_length_u8 = target_length as u8;

        let start_time = std::time::Instant::now();
        for i in 0..n {
            if lcs_simd.scan_left_bounded(i, target_length_u8, bound).is_err() {
                let i = i.saturating_sub(LcsSimd::LANES * bound);
                matrix.previous(i, target_length);
            }
        }
        let end_time = std::time::Instant::now();
        let nanos_per_previous = (end_time - start_time).as_nanos() as f64 / n as f64;

        let start_time = std::time::Instant::now();
        for i in 0..n {
            if lcs_simd.scan_right_bounded(i, target_length_u8, bound).is_err() {
                let i = i + LcsSimd::LANES * bound;
                matrix.next(i, target_length);
            }
        }
        let end_time = std::time::Instant::now();
        let nanos_per_next = (end_time - start_time).as_nanos() as f64 / n as f64;

        (nanos_per_next, nanos_per_previous)
    }

    fn statistics_augmented_bounded_scan(abs: &ABS, target_length: usize) -> (f64, f64) {
        let n = abs.lcs_simd.len();
        let target_length_u8 = target_length as u8;

        let start_time = std::time::Instant::now();
        for i in 0..n {
            let _ = abs.previous(i, target_length);
        }
        let end_time = std::time::Instant::now();
        let nanos_per_previous = (end_time - start_time).as_nanos() as f64 / n as f64;

        let start_time = std::time::Instant::now();
        for i in 0..n {
            let _ = abs.next(i, target_length);
        }
        let end_time = std::time::Instant::now();
        let nanos_per_next = (end_time - start_time).as_nanos() as f64 / n as f64;

        (nanos_per_next, nanos_per_previous)
    }

    #[ignore]
    #[test]
    fn simd_bounded_scan_with_fallback_time() {
        let (index, lcs) = read_index_and_lcs(4);
        let SbwtIndexVariant::SubsetMatrix(sbwt) = index;

        let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);

        let target_length_lower = 8;
        let target_length_upper = 10;

        println!("creating lcs_simd...");
        let lcs_simd = LcsSimd::from_iterator(iterator.clone(), lcs.len());

        println!("creating matrix...");
        let matrix = PnsvMatrix::from_iterator(iterator.clone(), lcs.len(), target_length_lower, target_length_upper);

        let item_bound: usize = 256;
        let word_bound = item_bound.div_ceil(LcsSimd::LANES);

        println!("creating augmented bounded scan...");
        let abs = ABS::from_iterator(&lcs_simd, iterator, lcs.len(), word_bound, target_length_lower, target_length_upper);

        println!("timing...");
        for target_length in target_length_lower..=target_length_upper {
            // let (
            //     nanos_per_previous_abs,
            //     nanos_per_next_abs
            // ) = statistics_augmented_bounded_scan(&abs, target_length);

            let (
                nanos_per_previous_matrix_fallback,
                nanos_per_next_matrix_fallback,
            ) = statistics_lcs_simd_with_matrix_fallback(&lcs_simd, &matrix, target_length, word_bound);

            let (
                nanos_per_previous_matrix,
                nanos_per_next_matrix
            ) = statistics_pnsv_matrix(&matrix, target_length);

            println!("target_length: {}", target_length);
            // println!(
            //     "augmented bounded scan    | previous: {:.3} | next: {:.3}",
            //     nanos_per_previous_abs,
            //     nanos_per_next_abs
            // );
            println!(
                "scan with matrix fallback | previous: {:.3} | next: {:.3}",
                nanos_per_previous_matrix_fallback,
                nanos_per_next_matrix_fallback
            );
            println!(
                "matrix only               | previous: {:.3} | next: {:.3}",
                nanos_per_previous_matrix,
                nanos_per_next_matrix
            );
            println!();
        }
    }
}
