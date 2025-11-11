use caf::types::CafHeader;

fn main() {
    // Create header like writer does
    let header1 = CafHeader::new(10000, vec![]);
    let bytes1 = bincode::serialize(&header1).unwrap();
    println!("Empty header size: {} bytes", bytes1.len());

    // Add references like we do in set_references
    let mut header2 = CafHeader::new(10000, vec![]);
    header2.num_refs = 3;
    header2.ref_names = vec!["chr1".to_string(), "chr2".to_string(), "chr22".to_string()];
    header2.ref_lengths = vec![248956422, 242193529, 50818468];
    let sam_header = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n@SQ\tSN:chr2\tLN:242193529\n@SQ\tSN:chr22\tLN:50818468\n@PG\tID:generate_test_bam\tPN:generate_test_bam\tVN:1.0\tCL:generate_test_bam.sh\n@PG\tID:samtools\tPN:samtools\tPP:generate_test_bam\tVN:1.22.1\tCL:samtools view -b -o ../test-data/synthetic_100000.bam ../test-data/synthetic_100000.sam\n@PG\tID:samtools.1\tPN:samtools\tPP:samtools\tVN:1.22.1\tCL:samtools sort -o ../test-data/synthetic_100000.bam.sorted ../test-data/synthetic_100000.bam\n".to_vec();
    header2.sam_header = sam_header;
    let bytes2 = bincode::serialize(&header2).unwrap();
    println!("Header with refs size: {} bytes", bytes2.len());

    println!("\nMagic (4) + Header = {} bytes total", 4 + bytes2.len());
    println!("First block should start at offset: {}", 4 + bytes2.len());
}
