use biometal::FastqStream;

fn main() -> biometal::Result<()> {
    let path = "/Users/scotthandley/Code/apple-silicon-bio-bench/datasets/huge_10m_150bp.fq.gz";

    let mut total_records = 0u64;
    let mut total_bases = 0u64;

    for record in FastqStream::from_path(path)? {
        let record = record?;
        total_bases += record.sequence.len() as u64;
        total_records += 1;
    }

    println!("Processed {} records ({} bases)", total_records, total_bases);
    Ok(())
}
