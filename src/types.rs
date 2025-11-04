//! Common types used throughout biometal

/// A FASTQ record
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    /// Sequence identifier (without '@' prefix)
    pub id: String,
    /// DNA/RNA sequence
    pub sequence: Vec<u8>,
    /// Quality scores (Phred+33)
    pub quality: Vec<u8>,
}

impl FastqRecord {
    /// Create a new FASTQ record
    pub fn new(id: String, sequence: Vec<u8>, quality: Vec<u8>) -> Self {
        Self { id, sequence, quality }
    }
}

/// A FASTA record
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    /// Sequence identifier (without '>' prefix)
    pub id: String,
    /// DNA/RNA/protein sequence
    pub sequence: Vec<u8>,
}

impl FastaRecord {
    /// Create a new FASTA record
    pub fn new(id: String, sequence: Vec<u8>) -> Self {
        Self { id, sequence }
    }
}
