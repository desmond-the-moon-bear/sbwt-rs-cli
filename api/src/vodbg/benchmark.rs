#![allow(unused)]

use crate::streaming_index::{ContractLeft, ExtendRight, StreamingIndex};
use crate::SbwtIndexVariant;
use crate::streaming_index::LcsArray;

// note(mk): Temporary solution for reading input data.
fn read_index_and_lcs() -> (SbwtIndexVariant, LcsArray) {
    let mut args = std::env::args().skip(4);
    let sbwt_path = args.next().expect("expected sbwt index path");
    let lcs_path = args.next().expect("expected lcs index path");

    let mut index_reader = std::io::BufReader::new(std::fs::File::open(sbwt_path).unwrap());
    let index = crate::load_sbwt_index_variant(&mut index_reader).unwrap();

    let mut lcs_reader = std::io::BufReader::new(std::fs::File::open(lcs_path).unwrap());
    let lcs = LcsArray::load(&mut lcs_reader).unwrap();

    (index, lcs)
}

fn read_query() -> Vec<Vec<u8>> {
    use crate::SeqStream;

    let mut args = std::env::args().skip(6);
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

fn benchmark_bms_separate_queries<E, C>(
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

fn benchmark_bms_joined_queries_lcs(bound: usize) {}

#[cfg(test)]
mod test {
    #![allow(unused)]

    use super::*;
    use crate::vodbg::pnsv::{
        LcsPnsvBp,
        LcsSimd,
        PnsvHybrid,
        Ranges,
    };

    // note(mk): Temporary solution for running the benchmark. Command to run these benchmarks:
    // cargo t (--release) vodbg::benchmark -- --ignored --nocapture sbwt_path lcs_path query_path
    #[ignore]
    #[test]
    fn comparison() {
        let (index, lcs) = read_index_and_lcs();
        let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
        let queries = read_query();

        // let lcs_pnsv = LcsPnsvBp::new(&lcs, 2048);

        let ranges = Ranges::new(&sbwt, sbwt.n_sets(), 7);
        let iterator = (0..lcs.len()).map(|index| lcs.access(index) as u8);
        let lcs_simd = LcsSimd::from_iterator(iterator, lcs.len());
        let pnsv_hybrid = PnsvHybrid {
            ranges,
            lcs_simd,
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
        
        let pnsv_hybrid_index = StreamingIndex {
            extend_right: &sbwt,
            contract_left: &pnsv_hybrid,
            n: sbwt.n_sets(),
            k: sbwt.k(),
        };

        let lower = 8;
        let upper = 31;

        for bound in 1..upper {
            print!("hyb,{},", bound);
            benchmark_bms_separate_queries(&pnsv_hybrid_index, &queries, bound);
            println!();
        }

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
    fn simd_scan_compare() {
        let (index, lcs) = read_index_and_lcs();
        let SbwtIndexVariant::SubsetMatrix(sbwt) = index;
        let queries = read_query();

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
}
