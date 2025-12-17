//! Biometal CLI - ARM-native bioinformatics toolkit
//!
//! A command-line interface for biometal's SIMD-accelerated bioinformatics primitives.
//! Provides Unix-style tools for sequence processing, quality analysis, and format conversion.
//!
//! # Performance
//!
//! - ARM NEON SIMD: 16.7-25.1× speedup on ARM64 platforms
//! - Streaming architecture: ~5MB constant memory usage
//! - Evidence-based optimization: 40,710+ measurements across 1,357 experiments
//!
//! # Usage
//!
//! ```bash
//! # Core statistics (NEON-optimized)
//! biometal count-bases input.fq
//! biometal gc-content input.fq
//! biometal mean-quality input.fq
//!
//! # Sequence operations
//! biometal reverse-complement input.fa
//! biometal trim-quality --end input.fq
//!
//! # Pattern matching
//! biometal find-pattern --pattern "AGATCGGAAGAGC" input.fq
//! ```

use std::env;
use std::process;

mod cli;

fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        print_usage();
        process::exit(1);
    }

    let command = &args[1];
    let command_args = &args[2..];

    match command.as_str() {
        // Tier 1 - Core Statistics (NEON-optimized 16.7-25.1× speedups)
        "count-bases" => cli::statistics::count_bases(command_args),
        "gc-content" => cli::statistics::gc_content(command_args),
        "mean-quality" => cli::statistics::mean_quality(command_args),
        "complexity-score" => cli::statistics::complexity_score(command_args),

        // Tier 1 - Essential Sequence Operations
        "reverse-complement" => cli::sequence::reverse_complement(command_args),
        "complement" => cli::sequence::complement(command_args),
        "reverse" => cli::sequence::reverse(command_args),
        "validate-dna" => cli::sequence::validate_dna(command_args),

        // Tier 2 - Quality Management
        "trim-quality" => cli::quality::trim_quality(command_args),
        "mask-low-quality" => cli::quality::mask_low_quality(command_args),
        "extract-region" => cli::quality::extract_region(command_args),

        // Tier 2 - Format Conversion & I/O
        "fastq-to-fasta" => cli::conversion::fastq_to_fasta(command_args),
        "count-reads" => cli::conversion::count_reads(command_args),

        // Tier 2 - Pattern Matching (NEW - extracted from SIMD work)
        "find-pattern" => cli::pattern::find_pattern(command_args),
        "count-pattern" => cli::pattern::count_pattern(command_args),
        "find-adapters" => cli::pattern::find_adapters(command_args),


        // Utilities
        "help" | "--help" | "-h" => print_help(),
        "version" | "--version" | "-V" => print_version(),

        _ => {
            eprintln!("Error: Unknown command '{}'", command);

            // Suggest similar commands
            let suggestions = suggest_similar_command(command);
            if !suggestions.is_empty() {
                eprintln!();
                eprintln!("Did you mean:");
                for suggestion in suggestions {
                    eprintln!("    biometal {}", suggestion);
                }
            }

            eprintln!();
            eprintln!("Run 'biometal --help' to see all available commands.");
            eprintln!("Run 'biometal <COMMAND> --help' for command-specific help.");
            process::exit(1);
        }
    }
}

fn print_usage() {
    println!("biometal {}", env!("CARGO_PKG_VERSION"));
    println!("{}", env!("CARGO_PKG_DESCRIPTION"));
    println!();
    println!("USAGE:");
    println!("    biometal <COMMAND> [OPTIONS] [INPUT]");
    println!();
    println!("COMMANDS:");
    println!("    Core Statistics (NEON-optimized):");
    println!("        count-bases        Count A/T/G/C frequencies");
    println!("        gc-content         Calculate GC percentage");
    println!("        mean-quality       Calculate mean quality scores");
    println!("        complexity-score   Calculate sequence complexity (Shannon entropy)");
    println!();
    println!("    Sequence Operations:");
    println!("        reverse-complement Reverse complement DNA sequences");
    println!("        complement         DNA complement only");
    println!("        reverse            Reverse sequences only");
    println!("        validate-dna       Validate DNA sequence characters");
    println!();
    println!("    Quality Management:");
    println!("        trim-quality       Trim low-quality bases");
    println!("        mask-low-quality   Replace low-quality bases with 'N'");
    println!("        extract-region     Extract subsequence with coordinates");
    println!();
    println!("    Format Conversion:");
    println!("        fastq-to-fasta     Convert FASTQ to FASTA format");
    println!("        count-reads        Count reads in file");
    println!();
    println!("    Pattern Matching:");
    println!("        find-pattern       Find pattern occurrences");
    println!("        count-pattern      Count pattern occurrences");
    println!("        find-adapters      Multi-pattern adapter detection");
    println!();
    println!("    Help:");
    println!("        help               Show this help message");
    println!("        version            Show version information");
    println!();
    println!("For command-specific help: biometal <COMMAND> --help");
    println!();
    println!("I/O Support: files, stdin/stdout, HTTP/HTTPS URLs, SRA accessions");
}

fn print_help() {
    print_usage();
}

fn print_version() {
    println!("biometal {}", env!("CARGO_PKG_VERSION"));
    println!("ARM-native bioinformatics library with streaming architecture");
    println!("Performance: 16.7-25.1× NEON speedup on ARM64 platforms");
    println!("Repository: {}", env!("CARGO_PKG_REPOSITORY"));
    println!("License: {}", env!("CARGO_PKG_LICENSE"));
}

/// Suggest similar commands based on common typos and alternative names
fn suggest_similar_command(input: &str) -> Vec<&'static str> {
    let all_commands = vec![
        "count-bases", "gc-content", "mean-quality", "complexity-score",
        "reverse-complement", "complement", "reverse", "validate-dna",
        "trim-quality", "mask-low-quality", "extract-region",
        "fastq-to-fasta", "count-reads",
        "find-pattern", "count-pattern", "find-adapters"
    ];

    let mut suggestions = Vec::new();

    // Direct matches for common alternative names
    match input.to_lowercase().as_str() {
        "count" | "wc" => suggestions.push("count-reads"),
        "bases" | "composition" => suggestions.push("count-bases"),
        "gc" | "gc-ratio" | "gc-percent" => suggestions.push("gc-content"),
        "quality" | "qual" => suggestions.push("mean-quality"),
        "entropy" | "complexity" => suggestions.push("complexity-score"),
        "revcomp" | "rc" => suggestions.push("reverse-complement"),
        "rev" => suggestions.push("reverse"),
        "comp" => suggestions.push("complement"),
        "validate" | "check" => suggestions.push("validate-dna"),
        "trim" => suggestions.push("trim-quality"),
        "mask" => suggestions.push("mask-low-quality"),
        "extract" | "substr" | "substring" => suggestions.push("extract-region"),
        "convert" | "fq2fa" | "fastq2fasta" => suggestions.push("fastq-to-fasta"),
        "search" | "grep" | "find" => suggestions.push("find-pattern"),
        "adapters" | "contamination" => suggestions.push("find-adapters"),
        _ => {}
    }

    // Fuzzy matching for typos (simple edit distance)
    if suggestions.is_empty() {
        for &command in &all_commands {
            if edit_distance(input, command) <= 2 && input.len() > 3 {
                suggestions.push(command);
            }
        }
    }

    // Limit to 3 suggestions to avoid overwhelming users
    suggestions.truncate(3);
    suggestions
}

/// Simple edit distance calculation for fuzzy matching
fn edit_distance(s1: &str, s2: &str) -> usize {
    let len1 = s1.len();
    let len2 = s2.len();
    let mut matrix = vec![vec![0; len2 + 1]; len1 + 1];

    for i in 0..=len1 {
        matrix[i][0] = i;
    }
    for j in 0..=len2 {
        matrix[0][j] = j;
    }

    for (i, c1) in s1.chars().enumerate() {
        for (j, c2) in s2.chars().enumerate() {
            let cost = if c1 == c2 { 0 } else { 1 };
            matrix[i + 1][j + 1] = std::cmp::min(
                std::cmp::min(
                    matrix[i][j + 1] + 1,      // deletion
                    matrix[i + 1][j] + 1,      // insertion
                ),
                matrix[i][j] + cost,           // substitution
            );
        }
    }

    matrix[len1][len2]
}