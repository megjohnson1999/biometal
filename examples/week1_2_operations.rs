//! Week 1-2 Operations Demo
//!
//! Demonstrates all core features implemented in Week 1-2:
//! - FASTQ streaming with quality filtering (Rule 1: 25.1× speedup)
//! - FASTA streaming with GC content (Rule 1: 20.3× speedup)
//! - Paired-end read streaming
//! - NEON-optimized operations

use biometal::operations::{count_bases, gc_content, mean_quality, passes_quality_filter};
use biometal::{FastaStream, FastqStream, PairedFastqStream};
use std::io::Cursor;

fn main() -> biometal::Result<()> {
    println!("biometal Week 1-2 Operations Demo");
    println!("===================================\n");

    demo_fastq_quality_filter()?;
    println!();
    demo_fasta_gc_content()?;
    println!();
    demo_paired_end()?;

    Ok(())
}

/// Demo: FASTQ streaming with quality filtering
fn demo_fastq_quality_filter() -> biometal::Result<()> {
    println!("1. FASTQ Quality Filtering (Rule 1: NEON 25.1× speedup)");
    println!("--------------------------------------------------------");

    let fastq_data = b"@read1 high quality
GATTACA
+
IIIIIII
@read2 low quality
ACGTACGT
+
!!!!!!!!
@read3 medium quality
TGCATGCA
+
55555555
";

    let stream = FastqStream::from_reader(Cursor::new(fastq_data));
    let min_quality = 20.0;

    println!("Filtering reads with mean quality >= {}\n", min_quality);

    let mut passed = 0;
    let mut filtered = 0;

    for record in stream {
        let record = record?;
        let mean_q = mean_quality(&record.quality);
        let passes = passes_quality_filter(&record.quality, min_quality);

        println!("  {}: mean_q={:.1} -> {}",
                 record.id,
                 mean_q,
                 if passes { "PASS ✓" } else { "FILTERED ✗" });

        if passes {
            passed += 1;
        } else {
            filtered += 1;
        }
    }

    println!("\nSummary: {} passed, {} filtered", passed, filtered);
    Ok(())
}

/// Demo: FASTA streaming with GC content calculation
fn demo_fasta_gc_content() -> biometal::Result<()> {
    println!("2. FASTA GC Content (Rule 1: NEON 20.3× speedup)");
    println!("-------------------------------------------------");

    let fasta_data = b">seq1 AT-rich region
ATATATATATATAT
>seq2 GC-rich region
GCGCGCGCGCGCGC
>seq3 balanced
ACGTACGTACGTACGT
>seq4 with N bases
ACGTNNNACGT
";

    let stream = FastaStream::from_reader(Cursor::new(fasta_data));

    println!("Calculating GC content for each sequence:\n");

    for record in stream {
        let record = record?;
        let gc = gc_content(&record.sequence);
        let bases = count_bases(&record.sequence);

        println!("  {}: GC={:.1}% (A={}, C={}, G={}, T={})",
                 record.id,
                 gc * 100.0,
                 bases[0], bases[1], bases[2], bases[3]);
    }

    Ok(())
}

/// Demo: Paired-end read streaming
fn demo_paired_end() -> biometal::Result<()> {
    println!("3. Paired-End Streaming (Rule 5: Constant Memory)");
    println!("--------------------------------------------------");

    let r1_data = b"@read1/1
GATTACA
+
IIIIIII
@read2/1
ACGTACGT
+
IIIIIIII
";

    let r2_data = b"@read1/2
TGTAATC
+
IIIIIII
@read2/2
ACGTACGT
+
IIIIIIII
";

    let stream1 = FastqStream::from_reader(Cursor::new(r1_data));
    let stream2 = FastqStream::from_reader(Cursor::new(r2_data));
    let paired = PairedFastqStream::from_streams(stream1, stream2);

    println!("Processing paired-end reads:\n");

    for pair in paired {
        let (r1, r2) = pair?;

        // Calculate quality and GC for both reads
        let r1_quality = mean_quality(&r1.quality);
        let r2_quality = mean_quality(&r2.quality);
        let r1_gc = gc_content(&r1.sequence);
        let r2_gc = gc_content(&r2.sequence);

        println!("  Pair: {} / {}",
                 r1.id.split_whitespace().next().unwrap_or(&r1.id),
                 r2.id.split_whitespace().next().unwrap_or(&r2.id));
        println!("    R1: Q={:.1}, GC={:.1}%",
                 r1_quality,
                 r1_gc * 100.0);
        println!("    R2: Q={:.1}, GC={:.1}%",
                 r2_quality,
                 r2_gc * 100.0);
    }

    Ok(())
}
