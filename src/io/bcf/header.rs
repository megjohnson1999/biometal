//! BCF header parsing and field dictionaries.

use crate::error::{BiometalError, Result};
use std::collections::HashMap;
use std::io::Read;

/// BCF header with field dictionaries.
///
/// Contains the VCF text header plus dictionaries mapping field names
/// to integer indices for efficient binary encoding.
#[derive(Debug, Clone)]
pub struct BcfHeader {
    /// Raw VCF header text
    text: String,
    /// Sample names
    samples: Vec<String>,
    /// FILTER name to index mapping (0 = PASS)
    filter_dict: HashMap<String, usize>,
    /// INFO field name to index mapping
    info_dict: HashMap<String, usize>,
    /// FORMAT field name to index mapping
    format_dict: HashMap<String, usize>,
    /// Contig names
    contigs: Vec<String>,
}

impl BcfHeader {
    /// Parse BCF header from reader.
    ///
    /// Reads the length-prefixed text header and builds field dictionaries.
    pub fn parse<R: Read>(reader: &mut R) -> Result<Self> {
        // Read header text length (u32 little-endian)
        let mut len_buf = [0u8; 4];
        reader.read_exact(&mut len_buf)?;
        let text_len = u32::from_le_bytes(len_buf) as usize;

        // Read header text
        let mut text_buf = vec![0u8; text_len];
        reader.read_exact(&mut text_buf)?;
        let text = String::from_utf8(text_buf).map_err(|e| BiometalError::InvalidInput {
            msg: format!("Invalid UTF-8 in BCF header: {}", e),
        })?;

        // Parse header to extract field definitions
        let mut samples = Vec::new();
        let mut filter_dict = HashMap::new();
        let mut info_dict = HashMap::new();
        let mut format_dict = HashMap::new();
        let mut contigs = Vec::new();

        // PASS is always index 0
        filter_dict.insert("PASS".to_string(), 0);
        let mut filter_idx = 1;
        let mut info_idx = 0;
        let mut format_idx = 0;

        for line in text.lines() {
            if line.starts_with("##FILTER=") {
                if let Some(id) = extract_id(line) {
                    let idx = extract_idx(line).unwrap_or_else(|| {
                        if id == "PASS" {
                            0
                        } else {
                            let next_idx = filter_idx;
                            filter_idx += 1;
                            next_idx
                        }
                    });
                    filter_dict.insert(id, idx);
                }
            } else if line.starts_with("##INFO=") {
                if let Some(id) = extract_id(line) {
                    let idx = extract_idx(line).unwrap_or_else(|| {
                        let next_idx = info_idx;
                        info_idx += 1;
                        next_idx
                    });
                    info_dict.insert(id, idx);
                }
            } else if line.starts_with("##FORMAT=") {
                if let Some(id) = extract_id(line) {
                    let idx = extract_idx(line).unwrap_or_else(|| {
                        let next_idx = format_idx;
                        format_idx += 1;
                        next_idx
                    });
                    format_dict.insert(id, idx);
                }
            } else if line.starts_with("##contig=") {
                if let Some(id) = extract_id(line) {
                    contigs.push(id);
                }
            } else if line.starts_with("#CHROM") {
                // Parse sample names from column header
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() > 9 {
                    samples = fields[9..].iter().map(|s| s.to_string()).collect();
                }
            }
        }

        Ok(BcfHeader {
            text,
            samples,
            filter_dict,
            info_dict,
            format_dict,
            contigs,
        })
    }

    /// Get header text.
    pub fn text(&self) -> &str {
        &self.text
    }

    /// Get sample names.
    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    /// Get contig names.
    pub fn contigs(&self) -> &[String] {
        &self.contigs
    }

    /// Get FILTER index by name.
    pub fn filter_index(&self, name: &str) -> Option<usize> {
        self.filter_dict.get(name).copied()
    }

    /// Get FILTER name by index.
    pub fn filter_name(&self, index: usize) -> Option<&str> {
        self.filter_dict
            .iter()
            .find(|(_, &idx)| idx == index)
            .map(|(name, _)| name.as_str())
    }

    /// Get INFO field index by name.
    pub fn info_index(&self, name: &str) -> Option<usize> {
        self.info_dict.get(name).copied()
    }

    /// Get INFO field name by index.
    pub fn info_name(&self, index: usize) -> Option<&str> {
        self.info_dict
            .iter()
            .find(|(_, &idx)| idx == index)
            .map(|(name, _)| name.as_str())
    }

    /// Get FORMAT field index by name.
    pub fn format_index(&self, name: &str) -> Option<usize> {
        self.format_dict.get(name).copied()
    }

    /// Get FORMAT field name by index.
    pub fn format_name(&self, index: usize) -> Option<&str> {
        self.format_dict
            .iter()
            .find(|(_, &idx)| idx == index)
            .map(|(name, _)| name.as_str())
    }

    /// Get number of samples.
    pub fn n_samples(&self) -> usize {
        self.samples.len()
    }
}

/// Extract ID from VCF header line.
///
/// Parses lines like: `##INFO=<ID=DP,Number=1,Type=Integer,...>`
fn extract_id(line: &str) -> Option<String> {
    line.split("ID=")
        .nth(1)?
        .split(',')
        .next()
        .map(|s| s.to_string())
}

/// Extract IDX from BCF header line.
///
/// Parses lines like: `##INFO=<ID=DP,...,IDX=2>`
fn extract_idx(line: &str) -> Option<usize> {
    line.split("IDX=")
        .nth(1)?
        .split('>')
        .next()?
        .parse()
        .ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn create_test_header() -> Vec<u8> {
        let header_text = "\
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description=\"All filters passed\">
##FILTER=<ID=LowQual,Description=\"Low quality\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2
";

        let mut data = Vec::new();
        // Write length (u32 little-endian)
        data.extend_from_slice(&(header_text.len() as u32).to_le_bytes());
        // Write text
        data.extend_from_slice(header_text.as_bytes());
        data
    }

    #[test]
    fn test_parse_header() {
        let data = create_test_header();
        let mut cursor = Cursor::new(data);

        let header = BcfHeader::parse(&mut cursor).unwrap();

        // Check samples
        assert_eq!(header.samples(), &["Sample1", "Sample2"]);
        assert_eq!(header.n_samples(), 2);

        // Check contigs
        assert_eq!(header.contigs(), &["chr1", "chr2"]);

        // Check FILTER dictionary
        assert_eq!(header.filter_index("PASS"), Some(0));
        assert_eq!(header.filter_index("LowQual"), Some(1));
        assert_eq!(header.filter_name(0), Some("PASS"));
        assert_eq!(header.filter_name(1), Some("LowQual"));

        // Check INFO dictionary
        assert_eq!(header.info_index("DP"), Some(0));
        assert_eq!(header.info_index("AF"), Some(1));
        assert_eq!(header.info_name(0), Some("DP"));
        assert_eq!(header.info_name(1), Some("AF"));

        // Check FORMAT dictionary
        assert_eq!(header.format_index("GT"), Some(0));
        assert_eq!(header.format_index("GQ"), Some(1));
        assert_eq!(header.format_name(0), Some("GT"));
        assert_eq!(header.format_name(1), Some("GQ"));
    }

    #[test]
    fn test_extract_id() {
        assert_eq!(
            extract_id("##INFO=<ID=DP,Number=1,Type=Integer>"),
            Some("DP".to_string())
        );
        assert_eq!(
            extract_id("##FILTER=<ID=PASS,Description=\"Pass\">"),
            Some("PASS".to_string())
        );
        assert_eq!(extract_id("##fileformat=VCFv4.2"), None);
    }
}
