//! K-mer operations for sequence analysis
//!
//! K-mer operations included pragmatically for contamination detection and alignment preprocessing.
//! These commands provide scalar-only implementations as k-mer operations are data-structure-bound.

use std::env;
use std::process;

/// Extract minimizers for long-read alignment preprocessing
///
/// Usage: biometal extract-minimizers [OPTIONS] [INPUT]
///
/// OPTIONS:
///     --kmer K           K-mer size (default: 21)
///     --window W         Window size for minimizer extraction (default: 10)
///     --output FILE      Output file (default: stdout)
///     --format FORMAT    Output format: positions (default), sequences, hash
///     --help             Show help message
///
/// INPUT:
///     FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified
pub fn extract_minimizers(args: &[String]) {
    let mut input_file = None;
    let mut output_file = None;
    let mut kmer_size = 21;
    let mut window_size = 10;
    let mut output_format = "positions";
    let mut i = 0;

    while i < args.len() {
        match args[i].as_str() {
            "--kmer" | "-k" => {
                if i + 1 < args.len() {
                    match args[i + 1].parse::<usize>() {
                        Ok(k) => {
                            if k == 0 || k > 32 {
                                eprintln!("Error: K-mer size must be between 1 and 32");
                                process::exit(1);
                            }
                            kmer_size = k;
                        }
                        Err(_) => {
                            eprintln!("Error: Invalid k-mer size");
                            process::exit(1);
                        }
                    }
                    i += 2;
                } else {
                    eprintln!("Error: --kmer requires a value");
                    process::exit(1);
                }
            }
            "--window" | "-w" => {
                if i + 1 < args.len() {
                    match args[i + 1].parse::<usize>() {
                        Ok(w) => {
                            if w == 0 {
                                eprintln!("Error: Window size must be > 0");
                                process::exit(1);
                            }
                            window_size = w;
                        }
                        Err(_) => {
                            eprintln!("Error: Invalid window size");
                            process::exit(1);
                        }
                    }
                    i += 2;
                } else {
                    eprintln!("Error: --window requires a value");
                    process::exit(1);
                }
            }
            "--output" | "-o" => {
                if i + 1 < args.len() {
                    output_file = Some(&args[i + 1]);
                    i += 2;
                } else {
                    eprintln!("Error: --output requires a value");
                    process::exit(1);
                }
            }
            "--format" => {
                if i + 1 < args.len() {
                    output_format = &args[i + 1];
                    i += 2;
                } else {
                    eprintln!("Error: --format requires a value");
                    process::exit(1);
                }
            }
            "--help" | "-h" => {
                print_extract_minimizers_help();
                return;
            }
            arg if !arg.starts_with('-') => {
                if input_file.is_none() {
                    input_file = Some(arg);
                    i += 1;
                } else {
                    eprintln!("Error: Multiple input files specified");
                    process::exit(1);
                }
            }
            _ => {
                eprintln!("Error: Unknown option '{}'", args[i]);
                process::exit(1);
            }
        }
    }

    match output_format {
        "positions" | "sequences" | "hash" => {}
        _ => {
            eprintln!("Error: Invalid format '{}'. Supported: positions, sequences, hash", output_format);
            process::exit(1);
        }
    }

    if window_size > kmer_size {
        eprintln!("Warning: Window size ({}) larger than k-mer size ({})", window_size, kmer_size);
    }

    println!("TODO: Implement extract_minimizers with biometal::operations::extract_minimizers_fast");
    println!("Input: {:?}", input_file.map_or("<stdin>", |f| f));
    println!("Output: {:?}", output_file.map_or("<stdout>", |f| f));
    println!("K-mer size: {}", kmer_size);
    println!("Window size: {}", window_size);
    println!("Format: {}", output_format);
}

/// Generate k-mer frequency spectrum for contamination detection
///
/// Usage: biometal kmer-spectrum [OPTIONS] [INPUT]
///
/// OPTIONS:
///     --kmer K           K-mer size (default: 21)
///     --output FILE      Output file (default: stdout)
///     --format FORMAT    Output format: histogram (default), counts, json
///     --min-count N      Minimum count to include (default: 1)
///     --help             Show help message
///
/// INPUT:
///     FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified
pub fn kmer_spectrum(args: &[String]) {
    let mut input_file = None;
    let mut output_file = None;
    let mut kmer_size = 21;
    let mut output_format = "histogram";
    let mut min_count = 1;
    let mut i = 0;

    while i < args.len() {
        match args[i].as_str() {
            "--kmer" | "-k" => {
                if i + 1 < args.len() {
                    match args[i + 1].parse::<usize>() {
                        Ok(k) => {
                            if k == 0 || k > 32 {
                                eprintln!("Error: K-mer size must be between 1 and 32");
                                process::exit(1);
                            }
                            kmer_size = k;
                        }
                        Err(_) => {
                            eprintln!("Error: Invalid k-mer size");
                            process::exit(1);
                        }
                    }
                    i += 2;
                } else {
                    eprintln!("Error: --kmer requires a value");
                    process::exit(1);
                }
            }
            "--output" | "-o" => {
                if i + 1 < args.len() {
                    output_file = Some(&args[i + 1]);
                    i += 2;
                } else {
                    eprintln!("Error: --output requires a value");
                    process::exit(1);
                }
            }
            "--format" => {
                if i + 1 < args.len() {
                    output_format = &args[i + 1];
                    i += 2;
                } else {
                    eprintln!("Error: --format requires a value");
                    process::exit(1);
                }
            }
            "--min-count" => {
                if i + 1 < args.len() {
                    match args[i + 1].parse::<u64>() {
                        Ok(c) => min_count = c,
                        Err(_) => {
                            eprintln!("Error: Invalid min-count value");
                            process::exit(1);
                        }
                    }
                    i += 2;
                } else {
                    eprintln!("Error: --min-count requires a value");
                    process::exit(1);
                }
            }
            "--help" | "-h" => {
                print_kmer_spectrum_help();
                return;
            }
            arg if !arg.starts_with('-') => {
                if input_file.is_none() {
                    input_file = Some(arg);
                    i += 1;
                } else {
                    eprintln!("Error: Multiple input files specified");
                    process::exit(1);
                }
            }
            _ => {
                eprintln!("Error: Unknown option '{}'", args[i]);
                process::exit(1);
            }
        }
    }

    match output_format {
        "histogram" | "counts" | "json" => {}
        _ => {
            eprintln!("Error: Invalid format '{}'. Supported: histogram, counts, json", output_format);
            process::exit(1);
        }
    }

    println!("TODO: Implement kmer_spectrum with biometal::operations::kmer_spectrum");
    println!("Input: {:?}", input_file.map_or("<stdin>", |f| f));
    println!("Output: {:?}", output_file.map_or("<stdout>", |f| f));
    println!("K-mer size: {}", kmer_size);
    println!("Format: {}", output_format);
    println!("Min count: {}", min_count);
}

fn print_extract_minimizers_help() {
    println!("biometal extract-minimizers - Extract minimizers for long-read alignment");
    println!();
    println!("USAGE:");
    println!("    biometal extract-minimizers [OPTIONS] [INPUT]");
    println!();
    println!("OPTIONS:");
    println!("    --kmer K, -k       K-mer size (default: 21, range: 1-32)");
    println!("    --window W, -w     Window size for minimizer extraction (default: 10)");
    println!("    --output FILE, -o  Output file (default: stdout)");
    println!("    --format FORMAT    Output format: positions (default), sequences, hash");
    println!("    --help, -h         Show this help message");
    println!();
    println!("INPUT:");
    println!("    FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified");
    println!();
    println!("PERFORMANCE:");
    println!("    Scalar implementation (k-mer operations are data-structure-bound)");
    println!("    Optimized sliding window with monotonic deque");
    println!();
    println!("EXAMPLES:");
    println!("    biometal extract-minimizers --kmer 21 --window 10 reads.fq");
    println!("    biometal extract-minimizers -k 15 -w 5 --format hash sequences.fa");
    println!("    biometal extract-minimizers --format sequences long_reads.fq");
}

fn print_kmer_spectrum_help() {
    println!("biometal kmer-spectrum - Generate k-mer frequency spectrum");
    println!();
    println!("USAGE:");
    println!("    biometal kmer-spectrum [OPTIONS] [INPUT]");
    println!();
    println!("OPTIONS:");
    println!("    --kmer K, -k       K-mer size (default: 21, range: 1-32)");
    println!("    --output FILE, -o  Output file (default: stdout)");
    println!("    --format FORMAT    Output format: histogram (default), counts, json");
    println!("    --min-count N      Minimum count to include (default: 1)");
    println!("    --help, -h         Show this help message");
    println!();
    println!("INPUT:");
    println!("    FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified");
    println!();
    println!("USE CASES:");
    println!("    - Contamination detection (unusual k-mer frequency patterns)");
    println!("    - Assembly quality assessment");
    println!("    - Genome size estimation");
    println!();
    println!("EXAMPLES:");
    println!("    biometal kmer-spectrum --kmer 21 sample.fq");
    println!("    biometal kmer-spectrum -k 17 --min-count 2 --format json reads.fq");
    println!("    biometal kmer-spectrum --format counts assembly.fa");
}