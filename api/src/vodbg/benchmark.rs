#![allow(unused)]

use crate::SbwtIndexVariant;
use crate::streaming_index::StreamingIndex;
use crate::{sbwt::SbwtIndex, streaming_index::LcsArray, subsetseq::SubsetSeq};

// note(mk): Temporary solution for reading input data.
pub fn read_index_and_lcs() -> (SbwtIndexVariant, LcsArray){
    let mut args = std::env::args().skip(4);
    let sbwt_path = args.next().expect("expected sbwt index path");
    let lcs_path = args.next().expect("expected lcs index path");

    let mut index_reader = std::io::BufReader::new(std::fs::File::open(sbwt_path).unwrap());
    let index = crate::load_sbwt_index_variant(&mut index_reader).unwrap();

    let mut lcs_reader = std::io::BufReader::new(std::fs::File::open(lcs_path).unwrap());
    let lcs = LcsArray::load(&mut lcs_reader).unwrap();

    (index, lcs)
}

pub fn read_query() {
    let mut args = std::env::args().skip(6);
    let query_path = args.next().expect("query path");
    let query_file = std::fs::File::open(query_path).unwrap();
    let input = std::io::BufReader::new(query_file);
    let input = crate::JSeqIOSeqStreamWrapper{ inner: jseqio::reader::DynamicFastXReader::new(input).unwrap() };
    // input.stream_next()
    todo!()
}

fn benchmark_bound_matching_statistics_lcs_only() {
//     log::info!("Querying {} random sequences of length {}", n_queries, query_len);
//     let start_time = std::time::Instant::now();
//     let mut checksum = 0_usize;
//     let mut n_kmers_queried = 0_usize;
//     for query in queries.iter(){
//         for x in index.matching_statistics(query).iter(){
//             checksum += x.0;
//         }
//         n_kmers_queried += std::cmp::max(query.len() as isize - sbwt.k() as isize + 1, 0) as usize;
//     }
//     let end_time = std::time::Instant::now();
//     log::info!("Sum of answers: {}", checksum);
//     println!("Elapsed time: {:.2} seconds", (end_time - start_time).as_secs_f64());
//     let nanos_per_kmer = (end_time - start_time).as_nanos() as f64 / n_kmers_queried as f64;
//     println!("{:.2} nanoseconds / k-mer", nanos_per_kmer);
//
//     nanos_per_kmer.round() as usize
    todo!()
}


#[cfg(test)]
mod test {
    use super::*;

    /// note(mk): Temporary solution for running the benchmark. Command to run these benchmarks:
    /// cargo t vodbg::benchmark -- --ignored --nocapture sbwt_path lcs_path query_path
    ///
    #[ignore]
    #[test]
    fn run_benchmark() {
        // read_index_and_lcs();
        todo!()
    }
}

