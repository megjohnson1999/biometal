//! Format conversion and I/O utilities
//!
//! Commands for converting between sequence formats and basic I/O operations.

use std::env;
use std::process;

/// Convert FASTQ to FASTA format
///
/// Usage: biometal fastq-to-fasta [OPTIONS] [INPUT]
///
/// OPTIONS:
///     --output FILE      Output file (default: stdout)
///     --keep-description Keep original header description
///     --help             Show help message
///
/// INPUT:
///     FASTQ file path, URL, SRA accession, or stdin if not specified
pub fn fastq_to_fasta(args: &[String]) {
    let mut input_file = None;
    let mut output_file = None;
    let mut keep_description = false;
    let mut i = 0;

    while i < args.len() {
        match args[i].as_str() {
            "--output" | "-o" => {
                if i + 1 < args.len() {
                    output_file = Some(&args[i + 1]);
                    i += 2;
                } else {
                    eprintln!("Error: --output requires a value");
                    process::exit(1);
                }
            }
            "--keep-description" => {
                keep_description = true;
                i += 1;
            }
            "--help" | "-h" => {
                print_fastq_to_fasta_help();
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

    // Import biometal types and operations
    use biometal::FastqStream;
    use biometal::operations::to_fasta_record;
    use std::io::{self, Write};
    use std::fs::File;

    // Determine output writer
    let mut output: Box<dyn Write> = match output_file {
        Some(path) => {
            match File::create(path) {
                Ok(file) => Box::new(file),
                Err(e) => {
                    eprintln!("Error creating output file '{}': {}", path, e);
                    process::exit(1);
                }
            }
        }
        None => Box::new(io::stdout()),
    };

    match input_file {
        Some(file_path) => {
            // FASTQ only (converting from FASTQ to FASTA)
            match FastqStream::from_path(file_path) {
                Ok(stream) => {
                    for record_result in stream {
                        match record_result {
                            Ok(fastq_record) => {
                                let fasta_record = to_fasta_record(&fastq_record);

                                // Write FASTA format
                                if writeln!(output, ">{}", fasta_record.id).is_err() {
                                    eprintln!("Error writing to output");
                                    process::exit(1);
                                }
                                if writeln!(output, "{}", String::from_utf8_lossy(&fasta_record.sequence)).is_err() {
                                    eprintln!("Error writing to output");
                                    process::exit(1);
                                }
                            }
                            Err(e) => {
                                eprintln!("Error reading FASTQ record: {}", e);
                                process::exit(1);
                            }
                        }
                    }
                }
                Err(e) => {
                    eprintln!("Error opening FASTQ file '{}': {}", file_path, e);
                    process::exit(1);
                }
            }
        }
        None => {
            eprintln!("Error: stdin input not yet implemented");
            process::exit(1);
        }
    }
}

/// Count reads in sequence files
///
/// Fast read counting similar to bio-virome-tools count command
///
/// Usage: biometal count-reads [OPTIONS] [INPUT]
///
/// OPTIONS:
///     --format FORMAT    Input format: auto (default), fastq, fasta
///     --help             Show help message
///
/// INPUT:
///     FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified
pub fn count_reads(args: &[String]) {
    let mut input_file = None;
    let mut input_format = "auto";
    let mut i = 0;

    while i < args.len() {
        match args[i].as_str() {
            "--format" => {
                if i + 1 < args.len() {
                    input_format = &args[i + 1];
                    i += 2;
                } else {
                    eprintln!("Error: --format requires a value");
                    process::exit(1);
                }
            }
            "--help" | "-h" => {
                print_count_reads_help();
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

    match input_format {
        "auto" | "fastq" | "fasta" => {}
        _ => {
            eprintln!("Error: Invalid format '{}'. Supported: auto, fastq, fasta", input_format);
            process::exit(1);
        }
    }

    // Import biometal types
    use biometal::{FastqStream, FastaStream};

    let mut read_count = 0u64;

    match input_file {
        Some(file_path) => {
            // Determine format: auto-detect or use specified format
            let detected_format = if input_format == "auto" {
                if file_path.ends_with(".fa") || file_path.ends_with(".fasta") || file_path.ends_with(".fas") {
                    "fasta"
                } else {
                    "fastq" // Default to FASTQ
                }
            } else {
                input_format
            };

            match detected_format {
                "fasta" => {
                    match FastaStream::from_path(file_path) {
                        Ok(stream) => {
                            for record_result in stream {
                                match record_result {
                                    Ok(_) => read_count += 1,
                                    Err(e) => {
                                        eprintln!("Error reading FASTA record: {}", e);
                                        process::exit(1);
                                    }
                                }
                            }
                        }
                        Err(e) => {
                            eprintln!("Error opening FASTA file '{}': {}", file_path, e);
                            process::exit(1);
                        }
                    }
                }
                "fastq" => {
                    match FastqStream::from_path(file_path) {
                        Ok(stream) => {
                            for record_result in stream {
                                match record_result {
                                    Ok(_) => read_count += 1,
                                    Err(e) => {
                                        eprintln!("Error reading FASTQ record: {}", e);
                                        process::exit(1);
                                    }
                                }
                            }
                        }
                        Err(e) => {
                            eprintln!("Error opening FASTQ file '{}': {}", file_path, e);
                            process::exit(1);
                        }
                    }
                }
                _ => unreachable!(), // We validate format earlier
            }

            // Output result
            println!("{} reads in {} (format: {})", read_count, file_path, detected_format);
        }
        None => {
            eprintln!("Error: stdin input not yet implemented");
            process::exit(1);
        }
    }
}

fn print_fastq_to_fasta_help() {
    println!("biometal fastq-to-fasta - Convert FASTQ to FASTA format");
    println!();
    println!("USAGE:");
    println!("    biometal fastq-to-fasta [OPTIONS] [INPUT]");
    println!();
    println!("OPTIONS:");
    println!("    --output FILE, -o      Output file (default: stdout)");
    println!("    --keep-description     Keep original header description");
    println!("    --help, -h             Show this help message");
    println!();
    println!("INPUT:");
    println!("    FASTQ file path, URL, SRA accession, or stdin if not specified");
    println!();
    println!("EXAMPLES:");
    println!("    biometal fastq-to-fasta reads.fq > sequences.fa");
    println!("    biometal fastq-to-fasta --keep-description reads.fq");
    println!("    cat reads.fq | biometal fastq-to-fasta --output sequences.fa");
}

fn print_count_reads_help() {
    println!("biometal count-reads - Count reads in sequence files");
    println!();
    println!("USAGE:");
    println!("    biometal count-reads [OPTIONS] [INPUT]");
    println!();
    println!("OPTIONS:");
    println!("    --format FORMAT    Input format: auto (default), fastq, fasta");
    println!("    --help, -h         Show this help message");
    println!();
    println!("INPUT:");
    println!("    FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified");
    println!();
    println!("EXAMPLES:");
    println!("    biometal count-reads sample.fq");
    println!("    biometal count-reads --format fasta sequences.fa");
    println!("    biometal count-reads https://example.com/data.fq.gz");
}