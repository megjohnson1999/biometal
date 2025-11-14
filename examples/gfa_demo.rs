//! GFA format parsing examples.
//!
//! This example demonstrates:
//! - Parsing GFA assembly graphs
//! - Path extraction and manipulation
//! - Graph statistics and analysis

use biometal::formats::gfa::{GfaParser, GfaRecord};
use std::collections::HashMap;
use std::io::Cursor;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== GFA Format Parsing Demo ===\n");

    // Sample assembly graph
    let gfa_data = "\
H\tVN:Z:1.0
S\tcontig1\tACGTACGTACGT\tLN:i:12\tRC:i:100
S\tcontig2\tTGCAGCTGCAGC\tLN:i:12\tRC:i:150
S\tcontig3\tGGATCCGATCCA\tLN:i:12\tRC:i:80
S\tcontig4\tCCATGGCATGCA\tLN:i:12\tRC:i:120
L\tcontig1\t+\tcontig2\t+\t4M
L\tcontig2\t+\tcontig3\t+\t3M
L\tcontig1\t+\tcontig4\t-\t2M
P\tpath1\tcontig1+,contig2+,contig3+\t4M,3M
P\tpath2\tcontig1+,contig4-\t2M
";

    // Parse all records
    println!("--- Parsing GFA Records ---");
    let parser = GfaParser::new(Cursor::new(gfa_data.as_bytes()));

    let mut segments = Vec::new();
    let mut links = Vec::new();
    let mut paths = Vec::new();
    let mut header = None;

    for result in parser {
        let record = result?;
        match record {
            GfaRecord::Header(h) => {
                println!("Header: Version {}", h.tags.get("VN").unwrap_or(&"unknown".to_string()));
                header = Some(h);
            }
            GfaRecord::Segment(s) => {
                println!(
                    "  Segment: {} (len={}, RC={})",
                    s.name,
                    s.length(),
                    s.tags.get("RC").unwrap_or(&"?".to_string())
                );
                segments.push(s);
            }
            GfaRecord::Link(l) => {
                println!(
                    "  Link: {}{} -> {}{}  (overlap: {})",
                    l.from_segment, l.from_orient, l.to_segment, l.to_orient, l.overlap
                );
                links.push(l);
            }
            GfaRecord::Path(p) => {
                println!(
                    "  Path: {} ({} segments: {})",
                    p.name,
                    p.segments.len(),
                    p.segments.join(" -> ")
                );
                paths.push(p);
            }
            _ => {}
        }
    }

    // Graph statistics
    println!("\n--- Graph Statistics ---");
    println!("  Total segments: {}", segments.len());
    println!("  Total links: {}", links.len());
    println!("  Total paths: {}", paths.len());

    // Calculate total graph length
    let total_length: usize = segments.iter().map(|s| s.length()).sum();
    println!("  Total sequence length: {} bp", total_length);

    // Segment degree analysis
    println!("\n--- Segment Connectivity ---");
    let mut in_degree: HashMap<String, usize> = HashMap::new();
    let mut out_degree: HashMap<String, usize> = HashMap::new();

    for link in &links {
        *out_degree.entry(link.from_segment.clone()).or_insert(0) += 1;
        *in_degree.entry(link.to_segment.clone()).or_insert(0) += 1;
    }

    for seg in &segments {
        let in_deg = in_degree.get(&seg.name).unwrap_or(&0);
        let out_deg = out_degree.get(&seg.name).unwrap_or(&0);
        println!(
            "  {}: in={}, out={} (total={})",
            seg.name,
            in_deg,
            out_deg,
            in_deg + out_deg
        );
    }

    // Path extraction and analysis
    println!("\n--- Path Analysis ---");
    for path in &paths {
        println!("  Path: {}", path.name);
        println!("    Segments: {}", path.segments.join(" -> "));

        // Calculate path length
        let mut path_length = 0;
        for seg_name in &path.segments {
            // Remove orientation suffix (+/-)
            let seg_id = seg_name.trim_end_matches(|c| c == '+' || c == '-');
            if let Some(seg) = segments.iter().find(|s| s.name == seg_id) {
                path_length += seg.length();
            }
        }
        println!("    Total length: {} bp", path_length);
        println!("    Overlaps: {}", path.overlaps.join(", "));
    }

    // Identify terminal segments (potential start/end points)
    println!("\n--- Terminal Segments ---");
    for seg in &segments {
        let in_deg = in_degree.get(&seg.name).unwrap_or(&0);
        let out_deg = out_degree.get(&seg.name).unwrap_or(&0);

        if *in_deg == 0 {
            println!("  {} is a source (no incoming links)", seg.name);
        }
        if *out_deg == 0 {
            println!("  {} is a sink (no outgoing links)", seg.name);
        }
    }

    // Identify branch points
    println!("\n--- Branch Points ---");
    for seg in &segments {
        let in_deg = *in_degree.get(&seg.name).unwrap_or(&0);
        let out_deg = *out_degree.get(&seg.name).unwrap_or(&0);

        if in_deg > 1 {
            println!("  {} has {} incoming branches (merge point)", seg.name, in_deg);
        }
        if out_deg > 1 {
            println!(
                "  {} has {} outgoing branches (divergence point)",
                seg.name, out_deg
            );
        }
    }

    // Coverage statistics
    println!("\n--- Coverage Statistics ---");
    let read_counts: Vec<usize> = segments
        .iter()
        .filter_map(|s| {
            s.tags
                .get("RC")
                .and_then(|v| v.split(':').last())
                .and_then(|v| v.parse::<usize>().ok())
        })
        .collect();

    if !read_counts.is_empty() {
        let min_cov = read_counts.iter().min().unwrap();
        let max_cov = read_counts.iter().max().unwrap();
        let avg_cov: usize = read_counts.iter().sum::<usize>() / read_counts.len();
        println!("  Read count range: {} - {} reads", min_cov, max_cov);
        println!("  Average coverage: {} reads", avg_cov);
    }

    println!("\n✓ GFA analysis complete!");

    Ok(())
}
