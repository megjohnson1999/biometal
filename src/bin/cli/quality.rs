//! Quality-based sequence operations
//!
//! Commands for quality-based trimming, masking, and region extraction.
//! Essential for daily preprocessing workflows in genomics.

use std::process;

/// Trim low-quality bases from sequences
///
/// Usage: biometal trim-quality [OPTIONS] [INPUT]
///
/// OPTIONS:
///     --start            Trim from start only
///     --end              Trim from end only (default)
///     --both             Trim from both ends
///     --window SIZE      Sliding window size (default: 4)
///     --threshold QUAL   Quality threshold (default: 20)
///     --output FILE      Output file (default: stdout)
///     --help             Show help message
///
/// INPUT:
///     FASTQ file path, URL, SRA accession, or stdin if not specified
pub fn trim_quality(args: &[String]) {
    let mut input_file = None;
    let mut output_file = None;
    let mut trim_mode = "end";
    let mut window_size = 4;
    let mut threshold = 20;
    let mut i = 0;

    while i < args.len() {
        match args[i].as_str() {
            "--start" => {
                trim_mode = "start";
                i += 1;
            }
            "--end" => {
                trim_mode = "end";
                i += 1;
            }
            "--both" => {
                trim_mode = "both";
                i += 1;
            }
            "--window" => {
                if i + 1 < args.len() {
                    match args[i + 1].parse::<usize>() {
                        Ok(w) => window_size = w,
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
            "--threshold" => {
                if i + 1 < args.len() {
                    match args[i + 1].parse::<u8>() {
                        Ok(t) => threshold = t,
                        Err(_) => {
                            eprintln!("Error: Invalid threshold value");
                            process::exit(1);
                        }
                    }
                    i += 2;
                } else {
                    eprintln!("Error: --threshold requires a value");
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
            "--help" | "-h" => {
                print_trim_quality_help();
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
    use biometal::operations::{trim_quality_start, trim_quality_end, trim_quality_both, trim_quality_window};
    use biometal::io::FastqWriter;
    
    use std::path::PathBuf;

    match input_file {
        Some(file_path) => {
            // FASTQ only (quality trimming requires quality scores)
            let mut writer = match output_file {
                Some(path) => {
                    match FastqWriter::create(&PathBuf::from(path)) {
                        Ok(w) => w,
                        Err(e) => {
                            eprintln!("Error creating output file '{}': {}", path, e);
                            process::exit(1);
                        }
                    }
                }
                None => {
                    match FastqWriter::stdout() {
                        Ok(w) => w,
                        Err(e) => {
                            eprintln!("Error creating stdout writer: {}", e);
                            process::exit(1);
                        }
                    }
                }
            };

            match FastqStream::from_path_streaming(file_path) {
                Ok(stream) => {
                    for record_result in stream {
                        match record_result {
                            Ok(record) => {
                                let trimmed_result = match trim_mode {
                                    "start" => trim_quality_start(&record, threshold),
                                    "end" => trim_quality_end(&record, threshold),
                                    "both" => trim_quality_both(&record, threshold),
                                    "window" => trim_quality_window(&record, threshold, window_size),
                                    _ => unreachable!(),
                                };

                                match trimmed_result {
                                    Ok(trimmed_record) => {
                                        if let Err(e) = writer.write_record(&trimmed_record) {
                                            eprintln!("Error writing FASTQ record: {}", e);
                                            process::exit(1);
                                        }
                                    }
                                    Err(e) => {
                                        eprintln!("Error trimming record '{}': {}", record.id, e);
                                        process::exit(1);
                                    }
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

/// Mask low-quality bases with 'N'
///
/// Usage: biometal mask-low-quality [OPTIONS] [INPUT]
///
/// OPTIONS:
///     --threshold QUAL   Quality threshold (default: 15)
///     --output FILE      Output file (default: stdout)
///     --help             Show help message
///
/// INPUT:
///     FASTQ file path, URL, SRA accession, or stdin if not specified
pub fn mask_low_quality(args: &[String]) {
    let mut input_file = None;
    let mut output_file = None;
    let mut threshold = 15;
    let mut i = 0;

    while i < args.len() {
        match args[i].as_str() {
            "--threshold" => {
                if i + 1 < args.len() {
                    match args[i + 1].parse::<u8>() {
                        Ok(t) => threshold = t,
                        Err(_) => {
                            eprintln!("Error: Invalid threshold value");
                            process::exit(1);
                        }
                    }
                    i += 2;
                } else {
                    eprintln!("Error: --threshold requires a value");
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
            "--help" | "-h" => {
                print_mask_low_quality_help();
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
    use biometal::operations::mask_low_quality_copy;
    use biometal::io::FastqWriter;
    
    use std::path::PathBuf;

    match input_file {
        Some(file_path) => {
            // FASTQ only (masking requires quality scores)
            let mut writer = match output_file {
                Some(path) => {
                    match FastqWriter::create(&PathBuf::from(path)) {
                        Ok(w) => w,
                        Err(e) => {
                            eprintln!("Error creating output file '{}': {}", path, e);
                            process::exit(1);
                        }
                    }
                }
                None => {
                    match FastqWriter::stdout() {
                        Ok(w) => w,
                        Err(e) => {
                            eprintln!("Error creating stdout writer: {}", e);
                            process::exit(1);
                        }
                    }
                }
            };

            match FastqStream::from_path_streaming(file_path) {
                Ok(stream) => {
                    for record_result in stream {
                        match record_result {
                            Ok(record) => {
                                match mask_low_quality_copy(&record, threshold) {
                                    Ok(masked_record) => {
                                        if let Err(e) = writer.write_record(&masked_record) {
                                            eprintln!("Error writing FASTQ record: {}", e);
                                            process::exit(1);
                                        }
                                    }
                                    Err(e) => {
                                        eprintln!("Error masking record '{}': {}", record.id, e);
                                        process::exit(1);
                                    }
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

/// Extract subsequence with coordinates
///
/// Usage: biometal extract-region [OPTIONS] [INPUT]
///
/// OPTIONS:
///     --start N          Start position (1-based, required)
///     --end N            End position (1-based, required)
///     --output FILE      Output file (default: stdout)
///     --help             Show help message
///
/// INPUT:
///     FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified
pub fn extract_region(args: &[String]) {
    let mut input_file = None;
    let mut output_file = None;
    let mut start_pos: Option<usize> = None;
    let mut end_pos: Option<usize> = None;
    let mut i = 0;

    while i < args.len() {
        match args[i].as_str() {
            "--start" => {
                if i + 1 < args.len() {
                    match args[i + 1].parse::<usize>() {
                        Ok(s) => start_pos = Some(s),
                        Err(_) => {
                            eprintln!("Error: Invalid start position");
                            process::exit(1);
                        }
                    }
                    i += 2;
                } else {
                    eprintln!("Error: --start requires a value");
                    process::exit(1);
                }
            }
            "--end" => {
                if i + 1 < args.len() {
                    match args[i + 1].parse::<usize>() {
                        Ok(e) => end_pos = Some(e),
                        Err(_) => {
                            eprintln!("Error: Invalid end position");
                            process::exit(1);
                        }
                    }
                    i += 2;
                } else {
                    eprintln!("Error: --end requires a value");
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
            "--help" | "-h" => {
                print_extract_region_help();
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

    if start_pos.is_none() {
        eprintln!("Error: --start is required");
        process::exit(1);
    }
    if end_pos.is_none() {
        eprintln!("Error: --end is required");
        process::exit(1);
    }

    if let (Some(start), Some(end)) = (start_pos, end_pos) {
        if start >= end {
            eprintln!("Error: Start position must be less than end position");
            process::exit(1);
        }
        if start == 0 {
            eprintln!("Error: Positions are 1-based (start must be ≥ 1)");
            process::exit(1);
        }
    }

    // Import biometal types and operations
    use biometal::{FastqStream, FastaStream, FastaRecord};
    use biometal::operations::extract_region;
    use biometal::io::{FastqWriter, fasta::FastaWriter};
    
    use std::path::PathBuf;

    let start = start_pos.unwrap();
    let end = end_pos.unwrap();

    // Convert to 0-based indexing for internal operations
    let start_idx = start - 1;
    let end_idx = end - 1;

    match input_file {
        Some(file_path) => {
            // Detect format from extension
            if file_path.ends_with(".fa") || file_path.ends_with(".fasta") || file_path.ends_with(".fas") {
                // FASTA format - manual extraction
                let mut writer = match output_file {
                    Some(path) => {
                        match FastaWriter::create(&PathBuf::from(path)) {
                            Ok(w) => w,
                            Err(e) => {
                                eprintln!("Error creating output file '{}': {}", path, e);
                                process::exit(1);
                            }
                        }
                    }
                    None => {
                        match FastaWriter::stdout() {
                            Ok(w) => w,
                            Err(e) => {
                                eprintln!("Error creating stdout writer: {}", e);
                                process::exit(1);
                            }
                        }
                    }
                };

                match FastaStream::from_path_streaming(file_path) {
                    Ok(stream) => {
                        for record_result in stream {
                            match record_result {
                                Ok(record) => {
                                    if end_idx >= record.sequence.len() {
                                        eprintln!("Error: End position {} exceeds sequence length {} for '{}'",
                                                 end, record.sequence.len(), record.id);
                                        process::exit(1);
                                    }

                                    let extracted_seq = &record.sequence[start_idx..=end_idx];
                                    let extracted_record = FastaRecord {
                                        id: record.id,
                                        sequence: extracted_seq.to_vec(),
                                    };

                                    if let Err(e) = writer.write_record(&extracted_record) {
                                        eprintln!("Error writing FASTA record: {}", e);
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
                // FASTQ format (default) - use biometal operation
                let mut writer = match output_file {
                    Some(path) => {
                        match FastqWriter::create(&PathBuf::from(path)) {
                            Ok(w) => w,
                            Err(e) => {
                                eprintln!("Error creating output file '{}': {}", path, e);
                                process::exit(1);
                            }
                        }
                    }
                    None => {
                        match FastqWriter::stdout() {
                            Ok(w) => w,
                            Err(e) => {
                                eprintln!("Error creating stdout writer: {}", e);
                                process::exit(1);
                            }
                        }
                    }
                };

                match FastqStream::from_path_streaming(file_path) {
                    Ok(stream) => {
                        for record_result in stream {
                            match record_result {
                                Ok(record) => {
                                    match extract_region(&record, start_idx, end_idx + 1) {
                                        Ok(extracted_record) => {
                                            if let Err(e) = writer.write_record(&extracted_record) {
                                                eprintln!("Error writing FASTQ record: {}", e);
                                                process::exit(1);
                                            }
                                        }
                                        Err(e) => {
                                            eprintln!("Error extracting region from '{}': {}", record.id, e);
                                            process::exit(1);
                                        }
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

fn print_trim_quality_help() {
    println!("biometal trim-quality - Trim low-quality bases from sequences");
    println!();
    println!("USAGE:");
    println!("    biometal trim-quality [OPTIONS] [INPUT]");
    println!();
    println!("OPTIONS:");
    println!("    --start            Trim from start only");
    println!("    --end              Trim from end only (default)");
    println!("    --both             Trim from both ends");
    println!("    --window SIZE      Sliding window size (default: 4)");
    println!("    --threshold QUAL   Quality threshold (default: 20)");
    println!("    --output FILE, -o  Output file (default: stdout)");
    println!("    --help, -h         Show this help message");
    println!();
    println!("INPUT:");
    println!("    FASTQ file path, URL, SRA accession, or stdin if not specified");
    println!();
    println!("EXAMPLES:");
    println!("    biometal trim-quality --end reads.fq");
    println!("    biometal trim-quality --both --threshold 25 reads.fq");
    println!("    biometal trim-quality --window 6 --threshold 15 reads.fq");
}

fn print_mask_low_quality_help() {
    println!("biometal mask-low-quality - Replace low-quality bases with 'N'");
    println!();
    println!("USAGE:");
    println!("    biometal mask-low-quality [OPTIONS] [INPUT]");
    println!();
    println!("OPTIONS:");
    println!("    --threshold QUAL   Quality threshold (default: 15)");
    println!("    --output FILE, -o  Output file (default: stdout)");
    println!("    --help, -h         Show this help message");
    println!();
    println!("INPUT:");
    println!("    FASTQ file path, URL, SRA accession, or stdin if not specified");
    println!();
    println!("EXAMPLES:");
    println!("    biometal mask-low-quality reads.fq");
    println!("    biometal mask-low-quality --threshold 20 reads.fq");
}

fn print_extract_region_help() {
    println!("biometal extract-region - Extract subsequence with coordinates");
    println!();
    println!("USAGE:");
    println!("    biometal extract-region [OPTIONS] [INPUT]");
    println!();
    println!("OPTIONS:");
    println!("    --start N          Start position (1-based, required)");
    println!("    --end N            End position (1-based, required)");
    println!("    --output FILE, -o  Output file (default: stdout)");
    println!("    --help, -h         Show this help message");
    println!();
    println!("INPUT:");
    println!("    FASTA/FASTQ file path, URL, SRA accession, or stdin if not specified");
    println!();
    println!("EXAMPLES:");
    println!("    biometal extract-region --start 10 --end 100 sequences.fa");
    println!("    biometal extract-region --start 1 --end 50 reads.fq");
}