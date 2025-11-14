use biometal::alignment::{smith_waterman_naive, smith_waterman, ScoringMatrix};

#[cfg(feature = "gpu")]
use biometal::alignment::smith_waterman_gpu;

fn main() {
    // Test the exact failing batch from property test
    let queries = vec![
        b"TGTCTTTAGGGTCAAC".as_slice(),
        b"T".as_slice(),
    ];
    let references = vec![
        b"AGAGGAACCCTCCCATTTCTACAA".as_slice(),
        b"ACTACACAGGATGACAGGGAGCTATACGCTGACGTGT".as_slice(),
    ];
    let scoring = ScoringMatrix::default();

    println!("Batch alignment test (failing property test case):\n");

    // Test each individually with CPU
    for (i, (query, reference)) in queries.iter().zip(references.iter()).enumerate() {
        let result_cpu = smith_waterman_naive(query, reference, &scoring);
        println!("Alignment {}:", i);
        println!("  Query: {}", String::from_utf8_lossy(query));
        println!("  Reference: {}", String::from_utf8_lossy(reference));
        println!("  CPU Score: {}", result_cpu.score);
    }

    #[cfg(feature = "gpu")]
    {
        use biometal::alignment::gpu::smith_waterman_batch_gpu;

        println!("\nBatch GPU processing:");
        match smith_waterman_batch_gpu(&queries, &references, &scoring) {
            Ok(results) => {
                for (i, result) in results.iter().enumerate() {
                    println!("  Alignment {}: GPU Score: {}", i, result.score);

                    // Compare with CPU
                    let cpu_result = smith_waterman_naive(queries[i], references[i], &scoring);
                    if cpu_result.score != result.score {
                        println!("    ❌ MISMATCH: CPU={}, GPU={}", cpu_result.score, result.score);
                    }
                }
            }
            Err(e) => {
                println!("  Error: {}", e);
            }
        }
    }
}
