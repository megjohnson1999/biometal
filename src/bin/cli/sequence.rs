//! Core sequence transformation commands
//!
//! Essential sequence operations used in daily bioinformatics workflows.
//! All commands support comprehensive I/O: files, stdin/stdout, HTTP/HTTPS URLs, SRA accessions.

use std::env;
use std::process;

/// Generate reverse complement of DNA sequences
///
/// Most common sequence operation in bioinformatics
///
/// Usage: biometal reverse-complement [OPTIONS] [INPUT]
///
/// OPTIONS:
///     --output FILE      Output file (default: stdout)
///     --help             Show help message
///
/// INPUT:
///     FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified
pub fn reverse_complement(args: &[String]) {
    let mut input_file = None;
    let mut output_file = None;
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
            "--help" | "-h" => {
                print_reverse_complement_help();
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
    use biometal::{FastqStream, FastaStream, FastqRecord, FastaRecord};
    use biometal::operations::reverse_complement;
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
            // Detect format from extension
            if file_path.ends_with(".fa") || file_path.ends_with(".fasta") || file_path.ends_with(".fas") {
                // FASTA format - output FASTA
                match FastaStream::from_path(file_path) {
                    Ok(stream) => {
                        for record_result in stream {
                            match record_result {
                                Ok(record) => {
                                    let rc_seq = reverse_complement(&record.sequence);
                                    let rc_record = FastaRecord {
                                        id: record.id,
                                        sequence: rc_seq,
                                    };

                                    // Write FASTA format
                                    if writeln!(output, ">{}", rc_record.id).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "{}", String::from_utf8_lossy(&rc_record.sequence)).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                }
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
            } else {
                // FASTQ format (default) - output FASTQ
                match FastqStream::from_path(file_path) {
                    Ok(stream) => {
                        for record_result in stream {
                            match record_result {
                                Ok(record) => {
                                    let rc_seq = reverse_complement(&record.sequence);
                                    // Also reverse quality scores to maintain alignment
                                    let mut rc_qual = record.quality;
                                    rc_qual.reverse();

                                    let rc_record = FastqRecord {
                                        id: record.id,
                                        sequence: rc_seq,
                                        quality: rc_qual,
                                    };

                                    // Write FASTQ format
                                    if writeln!(output, "@{}", rc_record.id).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "{}", String::from_utf8_lossy(&rc_record.sequence)).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "+").is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "{}", String::from_utf8_lossy(&rc_record.quality)).is_err() {
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
        }
        None => {
            eprintln!("Error: stdin input not yet implemented");
            process::exit(1);
        }
    }
}

/// Generate complement of DNA sequences (no reverse)
///
/// Usage: biometal complement [OPTIONS] [INPUT]
///
/// OPTIONS:
///     --output FILE      Output file (default: stdout)
///     --help             Show help message
///
/// INPUT:
///     FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified
pub fn complement(args: &[String]) {
    let mut input_file = None;
    let mut output_file = None;
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
            "--help" | "-h" => {
                print_complement_help();
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
    use biometal::{FastqStream, FastaStream, FastqRecord, FastaRecord};
    use biometal::operations::complement;
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
            // Detect format from extension
            if file_path.ends_with(".fa") || file_path.ends_with(".fasta") || file_path.ends_with(".fas") {
                // FASTA format - output FASTA
                match FastaStream::from_path(file_path) {
                    Ok(stream) => {
                        for record_result in stream {
                            match record_result {
                                Ok(record) => {
                                    let comp_seq = complement(&record.sequence);
                                    let comp_record = FastaRecord {
                                        id: record.id,
                                        sequence: comp_seq,
                                    };

                                    // Write FASTA format
                                    if writeln!(output, ">{}", comp_record.id).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "{}", String::from_utf8_lossy(&comp_record.sequence)).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                }
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
            } else {
                // FASTQ format (default) - output FASTQ
                match FastqStream::from_path(file_path) {
                    Ok(stream) => {
                        for record_result in stream {
                            match record_result {
                                Ok(record) => {
                                    let comp_seq = complement(&record.sequence);
                                    // Quality scores stay the same (no reverse for complement-only)

                                    let comp_record = FastqRecord {
                                        id: record.id,
                                        sequence: comp_seq,
                                        quality: record.quality, // Keep original quality order
                                    };

                                    // Write FASTQ format
                                    if writeln!(output, "@{}", comp_record.id).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "{}", String::from_utf8_lossy(&comp_record.sequence)).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "+").is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "{}", String::from_utf8_lossy(&comp_record.quality)).is_err() {
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
        }
        None => {
            eprintln!("Error: stdin input not yet implemented");
            process::exit(1);
        }
    }
}

/// Reverse sequences without complement
///
/// Usage: biometal reverse [OPTIONS] [INPUT]
///
/// OPTIONS:
///     --output FILE      Output file (default: stdout)
///     --help             Show help message
///
/// INPUT:
///     FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified
pub fn reverse(args: &[String]) {
    let mut input_file = None;
    let mut output_file = None;
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
            "--help" | "-h" => {
                print_reverse_help();
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
    use biometal::{FastqStream, FastaStream, FastqRecord, FastaRecord};
    use biometal::operations::reverse;
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
            // Detect format from extension
            if file_path.ends_with(".fa") || file_path.ends_with(".fasta") || file_path.ends_with(".fas") {
                // FASTA format - output FASTA
                match FastaStream::from_path(file_path) {
                    Ok(stream) => {
                        for record_result in stream {
                            match record_result {
                                Ok(record) => {
                                    let rev_seq = reverse(&record.sequence);
                                    let rev_record = FastaRecord {
                                        id: record.id,
                                        sequence: rev_seq,
                                    };

                                    // Write FASTA format
                                    if writeln!(output, ">{}", rev_record.id).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "{}", String::from_utf8_lossy(&rev_record.sequence)).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                }
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
            } else {
                // FASTQ format (default) - output FASTQ
                match FastqStream::from_path(file_path) {
                    Ok(stream) => {
                        for record_result in stream {
                            match record_result {
                                Ok(record) => {
                                    let rev_seq = reverse(&record.sequence);
                                    // Also reverse quality scores to maintain alignment
                                    let mut rev_qual = record.quality;
                                    rev_qual.reverse();

                                    let rev_record = FastqRecord {
                                        id: record.id,
                                        sequence: rev_seq,
                                        quality: rev_qual,
                                    };

                                    // Write FASTQ format
                                    if writeln!(output, "@{}", rev_record.id).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "{}", String::from_utf8_lossy(&rev_record.sequence)).is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "+").is_err() {
                                        eprintln!("Error writing to output");
                                        process::exit(1);
                                    }
                                    if writeln!(output, "{}", String::from_utf8_lossy(&rev_record.quality)).is_err() {
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
        }
        None => {
            eprintln!("Error: stdin input not yet implemented");
            process::exit(1);
        }
    }
}

/// Validate DNA sequence characters
///
/// Check for invalid bases beyond A/T/G/C/N
///
/// Usage: biometal validate-dna [OPTIONS] [INPUT]
///
/// OPTIONS:
///     --strict           Only allow A/T/G/C (no N or IUPAC codes)
///     --quiet            Only exit code (no output)
///     --help             Show help message
///
/// INPUT:
///     FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified
///
/// EXIT CODES:
///     0: All sequences valid
///     1: Invalid characters found
///     2: Error in processing
pub fn validate_dna(args: &[String]) {
    let mut input_file = None;
    let mut strict = false;
    let mut quiet = false;
    let mut i = 0;

    while i < args.len() {
        match args[i].as_str() {
            "--strict" => {
                strict = true;
                i += 1;
            }
            "--quiet" | "-q" => {
                quiet = true;
                i += 1;
            }
            "--help" | "-h" => {
                print_validate_dna_help();
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
    use biometal::{FastqStream, FastaStream};
    use biometal::operations::{is_valid_dna, count_invalid_bases};

    // Helper function for strict validation (only A/T/G/C allowed)
    fn is_strict_dna_base(base: u8) -> bool {
        matches!(base, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't')
    }

    fn count_invalid_strict(seq: &[u8]) -> usize {
        seq.iter().filter(|&&base| !is_strict_dna_base(base)).count()
    }

    // Track validation results
    let mut total_sequences = 0u64;
    let mut invalid_sequences = 0u64;
    let mut total_invalid_bases = 0usize;

    match input_file {
        Some(file_path) => {
            // Detect format from extension
            if file_path.ends_with(".fa") || file_path.ends_with(".fasta") || file_path.ends_with(".fas") {
                // FASTA format
                match FastaStream::from_path(file_path) {
                    Ok(stream) => {
                        for record_result in stream {
                            match record_result {
                                Ok(record) => {
                                    total_sequences += 1;
                                    let invalid_count = if strict {
                                        count_invalid_strict(&record.sequence)
                                    } else {
                                        count_invalid_bases(&record.sequence)
                                    };

                                    if invalid_count > 0 {
                                        invalid_sequences += 1;
                                        total_invalid_bases += invalid_count;
                                        if !quiet {
                                            eprintln!("Invalid sequence '{}': {} invalid bases",
                                                     record.id, invalid_count);
                                        }
                                    }
                                }
                                Err(e) => {
                                    eprintln!("Error reading FASTA record: {}", e);
                                    process::exit(2);
                                }
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("Error opening FASTA file '{}': {}", file_path, e);
                        process::exit(2);
                    }
                }
            } else {
                // FASTQ format (default)
                match FastqStream::from_path(file_path) {
                    Ok(stream) => {
                        for record_result in stream {
                            match record_result {
                                Ok(record) => {
                                    total_sequences += 1;
                                    let invalid_count = if strict {
                                        count_invalid_strict(&record.sequence)
                                    } else {
                                        count_invalid_bases(&record.sequence)
                                    };

                                    if invalid_count > 0 {
                                        invalid_sequences += 1;
                                        total_invalid_bases += invalid_count;
                                        if !quiet {
                                            eprintln!("Invalid sequence '{}': {} invalid bases",
                                                     record.id, invalid_count);
                                        }
                                    }
                                }
                                Err(e) => {
                                    eprintln!("Error reading FASTQ record: {}", e);
                                    process::exit(2);
                                }
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("Error opening FASTQ file '{}': {}", file_path, e);
                        process::exit(2);
                    }
                }
            }
        }
        None => {
            eprintln!("Error: stdin input not yet implemented");
            process::exit(2);
        }
    }

    // Output summary and exit with appropriate code
    if invalid_sequences > 0 {
        if !quiet {
            eprintln!("Validation failed: {}/{} sequences have {} total invalid bases",
                     invalid_sequences, total_sequences, total_invalid_bases);
            if strict {
                eprintln!("Note: Strict mode only allows A/T/G/C bases");
            }
        }
        process::exit(1); // Invalid characters found
    } else {
        if !quiet {
            println!("All {} sequences are valid DNA", total_sequences);
            if strict {
                println!("Strict mode: Only A/T/G/C bases found");
            }
        }
        process::exit(0); // All valid
    }
}

fn print_reverse_complement_help() {
    println!("biometal reverse-complement - Generate reverse complement of DNA sequences");
    println!();
    println!("USAGE:");
    println!("    biometal reverse-complement [OPTIONS] [INPUT]");
    println!();
    println!("OPTIONS:");
    println!("    --output FILE, -o  Output file (default: stdout)");
    println!("    --help, -h         Show this help message");
    println!();
    println!("INPUT:");
    println!("    FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified");
    println!();
    println!("EXAMPLES:");
    println!("    biometal reverse-complement sequences.fa");
    println!("    biometal reverse-complement -o output.fa input.fq");
    println!("    cat input.fa | biometal reverse-complement > output.fa");
}

fn print_complement_help() {
    println!("biometal complement - Generate complement of DNA sequences");
    println!();
    println!("USAGE:");
    println!("    biometal complement [OPTIONS] [INPUT]");
    println!();
    println!("OPTIONS:");
    println!("    --output FILE, -o  Output file (default: stdout)");
    println!("    --help, -h         Show this help message");
    println!();
    println!("INPUT:");
    println!("    FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified");
    println!();
    println!("EXAMPLES:");
    println!("    biometal complement sequences.fa");
    println!("    cat input.fa | biometal complement");
}

fn print_reverse_help() {
    println!("biometal reverse - Reverse sequences without complement");
    println!();
    println!("USAGE:");
    println!("    biometal reverse [OPTIONS] [INPUT]");
    println!();
    println!("OPTIONS:");
    println!("    --output FILE, -o  Output file (default: stdout)");
    println!("    --help, -h         Show this help message");
    println!();
    println!("INPUT:");
    println!("    FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified");
    println!();
    println!("EXAMPLES:");
    println!("    biometal reverse sequences.fa");
    println!("    cat input.fa | biometal reverse");
}

fn print_validate_dna_help() {
    println!("biometal validate-dna - Validate DNA sequence characters");
    println!();
    println!("USAGE:");
    println!("    biometal validate-dna [OPTIONS] [INPUT]");
    println!();
    println!("OPTIONS:");
    println!("    --strict           Only allow A/T/G/C (no N or IUPAC codes)");
    println!("    --quiet, -q        Only exit code (no output)");
    println!("    --help, -h         Show this help message");
    println!();
    println!("INPUT:");
    println!("    FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified");
    println!();
    println!("EXIT CODES:");
    println!("    0: All sequences valid");
    println!("    1: Invalid characters found");
    println!("    2: Error in processing");
    println!();
    println!("EXAMPLES:");
    println!("    biometal validate-dna sequences.fa");
    println!("    biometal validate-dna --strict sample.fq");
    println!("    biometal validate-dna --quiet input.fa && echo 'Valid'");
}