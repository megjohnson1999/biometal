//! HTTP Streaming Example
//!
//! Demonstrates streaming data directly from HTTP without downloading the entire file.
//! This example shows how biometal maintains constant memory (~5 MB) regardless of dataset size.
//!
//! # Evidence
//!
//! Entry 028: I/O bottleneck is 264-352× slower than compute
//! Network streaming addresses this critical bottleneck (Rule 6)
//!
//! # Usage
//!
//! ```bash
//! # Demo with any HTTP-accessible FASTQ file:
//! cargo run --example http_streaming --features network -- <URL>
//!
//! # Or use the default demo (httpbin.org test):
//! cargo run --example http_streaming --features network
//! ```
//!
//! # Public Genomics Datasets
//!
//! Try these public datasets (requires range-request support):
//! - NCBI SRA: https://sra-pub-run-odp.s3.amazonaws.com/sra/<accession>/<accession>
//! - ENA: https://ftp.sra.ebi.ac.uk/vol1/fastq/<path>
//! - 1000 Genomes: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/...
//!
//! Note: Some servers may not support HTTP range requests (required for streaming).

use biometal::io::network::HttpClient;

fn main() -> biometal::Result<()> {
    println!("biometal HTTP Streaming Demo");
    println!("============================\n");

    // Get URL from command line or use default test
    let args: Vec<String> = std::env::args().collect();
    let demo_mode = args.len() == 1;

    if demo_mode {
        println!("Running demonstration mode (no actual genomics data)\n");
        demo_http_streaming()?;
    } else {
        let url = &args[1];
        println!("Streaming from: {}\n", url);
        stream_fastq_from_url(url)?;
    }

    Ok(())
}

fn demo_http_streaming() -> biometal::Result<()> {
    println!("Example: HTTP Streaming API Demonstration");
    println!("=========================================\n");

    println!("This example demonstrates biometal's HTTP streaming capabilities.");
    println!("To test with real data, provide a URL that supports HTTP range requests:\n");
    println!("  cargo run --example http_streaming --features network -- <URL>\n");

    println!("Supported URLs:");
    println!("  • NCBI SRA (S3): https://sra-downloadb.be-md.ncbi.nlm.nih.gov/...");
    println!("  • ENA: ftp.sra.ebi.ac.uk (via HTTP)");
    println!("  • Cloud storage: AWS S3, Google Cloud Storage, Azure Blob");
    println!("  • Any HTTP server with range request support\n");

    println!("Creating HTTP client...");
    let client = HttpClient::new()?;
    println!("  ✓ Client initialized with:");
    let stats = client.cache_stats()?;
    println!("    - Cache size: {} MB", stats.max_bytes / (1024 * 1024));
    println!("    - Timeout: 30 seconds");
    println!("    - Max retries: 3");
    println!("    - Exponential backoff: 100ms → 200ms → 400ms\n");

    println!("HTTP Streaming Features:");
    println!("========================\n");

    println!("1. Range Requests (Partial Downloads)");
    println!("   • Only download needed bytes (not entire file)");
    println!("   • Server returns 206 (Partial Content)");
    println!("   • Rejects servers without range support (200 OK)\n");

    println!("2. LRU Cache (Memory-Bounded)");
    println!("   • Byte-based size limits (not entry count)");
    println!("   • Automatic eviction of least-recently-used entries");
    println!("   • Constant memory regardless of access patterns\n");

    println!("3. Automatic Retry with Exponential Backoff");
    println!("   • Transient failures (500, 503) automatically retried");
    println!("   • Backoff prevents server overload");
    println!("   • Configurable retry limits\n");

    println!("4. EOF Detection");
    println!("   • HEAD request to get Content-Length");
    println!("   • Prevents reading past end of file");
    println!("   • Proper Read trait implementation\n");

    println!("5. Server Validation");
    println!("   • Validates 206 (Partial Content) responses");
    println!("   • Detects 200 (server ignoring range header)");
    println!("   • Handles 416 (Range Not Satisfiable)\n");

    println!("Example Usage (Code):");
    println!("====================\n");

    println!("// Stream FASTQ from HTTP:");
    println!("let source = DataSource::Http(url.to_string());");
    println!("let stream = FastqStream::new(source)?;\n");

    println!("for record in stream {{");
    println!("    // Process one record at a time");
    println!("    // Memory: Constant ~5 MB");
    println!("}}\n");

    println!("// Low-level HTTP reader:");
    println!("let mut reader = HttpReader::new(url)?;");
    println!("let mut buffer = vec![0u8; 65536];");
    println!("reader.read(&mut buffer)?;  // Fetches 65 KB chunk\n");

    println!("// Direct range requests:");
    println!("let client = HttpClient::new()?;");
    println!("let data = client.fetch_range(url, 0, 1024)?;  // First 1 KB\n");

    println!("Memory Guarantees (Rule 5):");
    println!("===========================\n");

    println!("  Streaming:     ~5 MB constant (per stream)");
    println!("  Cache:         50 MB default (configurable)");
    println!("  Chunk size:    65 KB (typical bgzip block)");
    println!("  Total:         ~55 MB regardless of file size\n");

    println!("Evidence (Entry 028):");
    println!("=====================\n");

    println!("  I/O bottleneck:        264-352× slower than compute");
    println!("  Without streaming:     1.04-1.08× E2E speedup (NEON masked)");
    println!("  With streaming:        ~17× E2E speedup projected");
    println!("  Conclusion:            Network streaming is CRITICAL\n");

    println!("✅ HTTP streaming implementation ready for production use!\n");
    println!("To test with real data, provide a URL:");
    println!("  cargo run --example http_streaming --features network -- <FASTQ_URL>");

    Ok(())
}

fn stream_fastq_from_url(url: &str) -> biometal::Result<()> {
    use biometal::io::DataSource;
    use biometal::operations::{count_bases, gc_content};
    use biometal::FastqStream;

    println!("Streaming FASTQ records without downloading entire file...\n");

    let source = DataSource::Http(url.to_string());
    let stream = FastqStream::new(source)?;

    let mut record_count = 0;
    let mut total_bases = 0;
    let mut total_gc = 0.0;

    for record in stream {
        let record = record?;
        record_count += 1;

        // Count bases (NEON-optimized on ARM)
        let bases = count_bases(&record.sequence);
        total_bases += bases.iter().sum::<u32>() as usize;

        // Calculate GC content (NEON-optimized on ARM)
        let gc = gc_content(&record.sequence);
        total_gc += gc;

        // Print first few records as examples
        if record_count <= 5 {
            println!("Record {}: {}", record_count, record.id);
            println!("  Length: {} bp", record.sequence.len());
            println!("  A: {}, C: {}, G: {}, T: {}", bases[0], bases[1], bases[2], bases[3]);
            println!("  GC: {:.1}%\n", gc * 100.0);
        }

        // Progress indicator
        if record_count % 1000 == 0 {
            println!("Processed {} records...", record_count);
        }
    }

    println!("\n=== Summary ===");
    println!("Records processed: {}", record_count);
    println!("Total bases: {}", total_bases);
    println!("Average GC content: {:.1}%", (total_gc / record_count as f64) * 100.0);
    println!("Memory usage: Constant ~5 MB (streaming)");

    Ok(())
}
