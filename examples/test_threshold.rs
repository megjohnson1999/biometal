//! Test the corrected threshold logic with real files

use std::env;

fn main() -> biometal::Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: cargo run --example test_threshold <file.fastq.gz>");
        std::process::exit(1);
    }

    let file_path = &args[1];
    
    println!("Testing threshold logic with: {}", file_path);
    
    // Check file size
    let metadata = std::fs::metadata(file_path)?;
    let file_size_mb = metadata.len() / (1024 * 1024);
    println!("File size: {} MB", file_size_mb);
    
    // With new logic:
    // - Files ≤200MB should use memory mapping  
    // - Files >200MB should use streaming
    if file_size_mb > 200 {
        println!("Expected: Streaming I/O (file > 200MB threshold)");
    } else {
        println!("Expected: Memory mapping (file ≤ 200MB threshold)");
    }
    
    // Test with regular from_path (should use threshold logic)
    let stream = biometal::FastqStream::from_path(file_path)?;
    
    // Count a few records to test functionality  
    let mut count = 0;
    for (i, record_result) in stream.enumerate() {
        let _record = record_result?;
        count += 1;
        if i >= 1000 { // Just test first 1000 records
            break;
        }
    }
    
    println!("✅ Successfully processed {} records using threshold logic", count);
    
    Ok(())
}
