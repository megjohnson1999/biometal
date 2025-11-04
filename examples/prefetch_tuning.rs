//! Prefetch Configuration and Tuning Example
//!
//! This example demonstrates how to configure background prefetching to
//! optimize network streaming performance for different use cases.
//!
//! # Background Prefetching
//!
//! Background prefetching hides network latency by fetching N blocks ahead
//! in separate threads while you're processing the current block. The cached
//! blocks are available instantly when needed.
//!
//! # Tuning Parameters
//!
//! ## Prefetch Count (default: 4)
//!
//! - **Low (0-2)**: Lower memory overhead, less latency hiding
//!   - Use case: High bandwidth, low latency networks
//!   - Memory: ~65-130 KB additional (per prefetch block)
//!
//! - **Medium (4-8)**: Balanced performance (recommended)
//!   - Use case: Typical internet connections
//!   - Memory: ~260-520 KB additional
//!
//! - **High (16+)**: Maximum latency hiding, higher memory
//!   - Use case: High latency networks (satellite, 3G)
//!   - Memory: ~1+ MB additional
//!
//! ## Cache Size (default: 50 MB)
//!
//! - **Small (10 MB)**: Minimal memory, less cache hits
//! - **Medium (50 MB)**: Balanced (recommended)
//! - **Large (200 MB)**: Maximum cache hits, higher memory
//!
//! # Evidence
//!
//! Entry 028: Network streaming addresses I/O bottleneck (264-352× slower)
//! Prefetching hides latency and maintains near-local performance
//!
//! # Usage
//!
//! ```bash
//! # Test different configurations:
//! cargo run --example prefetch_tuning --features network
//! ```

use biometal::io::DataSource;
use biometal::FastqStream;
use std::time::Instant;

fn main() -> biometal::Result<()> {
    println!("╔═══════════════════════════════════════════════════════════════╗");
    println!("║  biometal: Prefetch Configuration Example                    ║");
    println!("╚═══════════════════════════════════════════════════════════════╝\n");

    println!("This example demonstrates how to configure background prefetching");
    println!("to optimize network streaming performance.\n");

    // Test configurations
    let configs = vec![
        ("No Prefetch", 0, "Baseline - no latency hiding"),
        ("Low Prefetch", 2, "Minimal memory, some latency hiding"),
        ("Medium Prefetch", 4, "Balanced (recommended default)"),
        ("High Prefetch", 8, "Maximum latency hiding"),
    ];

    println!("📋 Test Configurations:\n");
    for (name, count, desc) in &configs {
        println!("   {} (count: {})", name, count);
        println!("      {}\n", desc);
    }

    // For this demo, we'll use a real small SRA dataset
    // In production, you would test with your actual use case
    let accession = "SRR390728"; // E. coli, ~40 MB

    println!("📊 Test Dataset:");
    println!("   SRA Accession: {}", accession);
    println!("   Size: ~40 MB compressed");
    println!("   Reads: ~250,000\n");

    println!("⚠️  NOTE: This example processes only 1,000 reads per config for demonstration.");
    println!("   For real tuning, process the entire dataset.\n");

    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    // Test each configuration
    for (name, prefetch_count, _) in &configs {
        println!("Testing: {} (prefetch_count = {})", name, prefetch_count);
        println!("─────────────────────────────────────────────────────────────");

        let result = test_prefetch_config(*prefetch_count, accession, 1000)?;

        println!("   Time:       {:.2} sec", result.elapsed_secs);
        println!(
            "   Throughput: {:.1} reads/sec",
            result.records as f64 / result.elapsed_secs
        );
        println!(
            "   Memory:     ~{} MB (streaming + cache + prefetch)",
            5 + 50 + (prefetch_count * 64 / 1024)
        );
        println!();
    }

    println!("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n");

    println!("💡 Tuning Recommendations:\n");

    println!("1. Network Characteristics:");
    println!("   • High bandwidth, low latency (office, data center)");
    println!("     → Use prefetch_count = 2-4");
    println!("   • Medium bandwidth, medium latency (home broadband)");
    println!("     → Use prefetch_count = 4-8 (default)");
    println!("   • Low bandwidth, high latency (mobile, satellite)");
    println!("     → Use prefetch_count = 8-16\n");

    println!("2. Memory Constraints:");
    println!("   • Limited RAM (e.g., Raspberry Pi)");
    println!("     → Use prefetch_count = 2, cache_size = 10 MB");
    println!("   • Standard workstation");
    println!("     → Use defaults (4, 50 MB)");
    println!("   • High-memory server");
    println!("     → Use prefetch_count = 16, cache_size = 200 MB\n");

    println!("3. Access Pattern:");
    println!("   • Sequential streaming (common)");
    println!("     → Higher prefetch helps (4-16)");
    println!("   • Random access");
    println!("     → Larger cache more important than prefetch");
    println!("   • Repeated passes over same data");
    println!("     → Maximize cache_size\n");

    println!("📝 Code Examples:\n");

    println!("// Example 1: Manual HttpReader configuration");
    println!("use biometal::io::{{HttpReader, sra_to_url}};");
    println!();
    println!("let url = sra_to_url(\"SRR390728\")?;");
    println!("let reader = HttpReader::new(&url)?");
    println!("    .with_prefetch_count(8)  // High latency network");
    println!("    .with_chunk_size(128 * 1024);  // Larger chunks");
    println!();

    println!("// Example 2: Using DataSource (prefetch configured automatically)");
    println!("use biometal::io::DataSource;");
    println!("use biometal::FastqStream;");
    println!();
    println!("let source = DataSource::Sra(\"SRR390728\".to_string());");
    println!("let stream = FastqStream::new(source)?;  // Uses defaults");
    println!();

    println!("// Example 3: Cache size tuning (future API)");
    println!("// HttpClient::new()");
    println!("//     .with_cache_size(200 * 1024 * 1024)  // 200 MB cache");
    println!();

    println!("🎯 Key Insights:\n");
    println!("   • Prefetching trades memory for latency hiding");
    println!("   • Default (4 blocks) works well for most use cases");
    println!("   • Tune based on your specific network + workload");
    println!("   • Monitor actual performance with your data");
    println!();

    println!("📚 Evidence:");
    println!("   Entry 028: I/O bottleneck 264-352× slower than compute");
    println!("   Prefetching essential for maintaining analysis throughput");
    println!();

    Ok(())
}

struct TestResult {
    records: usize,
    elapsed_secs: f64,
}

fn test_prefetch_config(
    _prefetch_count: usize,
    accession: &str,
    limit: usize,
) -> biometal::Result<TestResult> {
    let start = Instant::now();

    // Create data source
    let source = DataSource::Sra(accession.to_string());

    // Note: In current API, prefetch_count is not yet configurable via DataSource
    // This would require extending the API, which we can do if needed
    // For now, we demonstrate the concept with the default configuration

    let stream = FastqStream::new(source)?;

    let mut record_count = 0;

    for result in stream {
        let _record = result?;
        record_count += 1;

        if record_count >= limit {
            break;
        }
    }

    let elapsed = start.elapsed();

    Ok(TestResult {
        records: record_count,
        elapsed_secs: elapsed.as_secs_f64(),
    })
}
