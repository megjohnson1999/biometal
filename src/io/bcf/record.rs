//! BCF record parsing.

use super::header::BcfHeader;
use super::typed_value::TypedValue;
use crate::error::{BiometalError, Result};
use std::collections::HashMap;
use std::io::Read;

/// A BCF record representing a variant.
///
/// Contains both shared data (applies to all samples) and
/// individual genotype data (sample-specific).
#[derive(Debug, Clone)]
pub struct BcfRecord {
    /// Shared data
    shared: SharedData,
    /// Individual/genotype data
    individual: Vec<HashMap<String, TypedValue>>,
}

/// Shared variant data (applies to all samples).
#[derive(Debug, Clone)]
struct SharedData {
    /// Chromosome/contig index
    chrom_id: i32,
    /// Position (0-based in BCF, convert to 1-based for VCF)
    pos: i32,
    /// Reference sequence end position (0-based exclusive)
    rlen: i32,
    /// Quality score
    qual: f32,
    /// Number of alleles (including reference)
    n_allele: u16,
    /// Number of INFO fields
    n_info: u16,
    /// Number of genotype/FORMAT fields
    n_fmt: u8,
    /// Number of samples (for validation)
    n_sample: u32,
    /// Variant ID
    id: Option<String>,
    /// Reference allele
    reference: String,
    /// Alternate alleles
    alternate: Vec<String>,
    /// FILTER indices
    filter: Vec<usize>,
    /// INFO fields (name -> value)
    info: HashMap<String, TypedValue>,
}

impl BcfRecord {
    /// Parse a BCF record from reader.
    ///
    /// Reads the complete BCF record including shared and individual data.
    pub fn parse<R: Read>(reader: &mut R, header: &BcfHeader) -> Result<Self> {
        // Read record size (u32 little-endian)
        let mut size_buf = [0u8; 4];
        reader.read_exact(&mut size_buf)?;
        let l_shared = u32::from_le_bytes(size_buf);

        let mut size_buf = [0u8; 4];
        reader.read_exact(&mut size_buf)?;
        let l_indiv = u32::from_le_bytes(size_buf);

        // Read shared data
        let shared = Self::parse_shared(reader, header, l_shared as usize)?;

        // Read individual data
        let individual = Self::parse_individual(reader, header, &shared, l_indiv as usize)?;

        Ok(BcfRecord { shared, individual })
    }

    /// Parse shared data section.
    fn parse_shared<R: Read>(
        reader: &mut R,
        header: &BcfHeader,
        size: usize,
    ) -> Result<SharedData> {
        // Read into buffer for easier processing
        let mut buf = vec![0u8; size];
        reader.read_exact(&mut buf)?;
        let mut cursor = std::io::Cursor::new(buf);

        // Read fixed fields (24 bytes)
        let mut fixed = [0u8; 24];
        cursor.read_exact(&mut fixed)?;

        let chrom_id = i32::from_le_bytes([fixed[0], fixed[1], fixed[2], fixed[3]]);
        let pos = i32::from_le_bytes([fixed[4], fixed[5], fixed[6], fixed[7]]);
        let rlen = i32::from_le_bytes([fixed[8], fixed[9], fixed[10], fixed[11]]);
        let qual = f32::from_le_bytes([fixed[12], fixed[13], fixed[14], fixed[15]]);
        // n_allele (16 bits) and n_info (16 bits) packed into 32 bits
        let n_allele_info = u32::from_le_bytes([fixed[16], fixed[17], fixed[18], fixed[19]]);
        let n_info = (n_allele_info & 0xFFFF) as u16;  // Lower 16 bits
        let n_allele = ((n_allele_info >> 16) & 0xFFFF) as u16;  // Upper 16 bits
        let n_fmt_sample = [fixed[20], fixed[21], fixed[22], fixed[23]];
        let n_fmt = n_fmt_sample[0];
        let n_sample = u32::from_le_bytes([n_fmt_sample[1], n_fmt_sample[2], n_fmt_sample[3], 0]);

        // Read ID
        let id_value = TypedValue::read(&mut cursor)?;
        let id = match id_value {
            TypedValue::String(s) if s == "." => None,
            TypedValue::String(s) => Some(s),
            TypedValue::Missing => None,
            _ => None,
        };

        // Read alleles (REF + ALT)
        // BCF2 stores each allele as a separate typed value
        let mut alleles = Vec::new();
        for _i in 0..n_allele {
            let allele_value = TypedValue::read(&mut cursor)?;
            match allele_value {
                TypedValue::String(s) => alleles.push(s),
                _ => {
                    return Err(BiometalError::InvalidInput {
                        msg: format!("Expected string for allele, got {:?}", allele_value),
                    })
                }
            }
        }

        if alleles.is_empty() {
            return Err(BiometalError::InvalidInput {
                msg: "No alleles found".to_string(),
            });
        }

        let reference = alleles.get(0).cloned().ok_or_else(|| BiometalError::InvalidInput {
            msg: "Missing reference allele".to_string(),
        })?;
        let alternate = alleles[1..].to_vec();

        // Read FILTER
        let filter_value = TypedValue::read(&mut cursor)?;
        let filter = match filter_value {
            TypedValue::Int8(v) => vec![v as usize],
            TypedValue::Int8Array(v) => v.iter().map(|&x| x as usize).collect(),
            TypedValue::Int16(v) => vec![v as usize],
            TypedValue::Int16Array(v) => v.iter().map(|&x| x as usize).collect(),
            TypedValue::Missing => vec![],
            _ => vec![],
        };

        // Read INFO fields
        let mut info = HashMap::new();
        for _ in 0..n_info {
            // Read INFO key (index into header dictionary)
            let key_value = TypedValue::read(&mut cursor)?;
            let key_idx = key_value.as_int().ok_or_else(|| BiometalError::InvalidInput {
                msg: "Expected integer for INFO key".to_string(),
            })? as usize;

            let key_name = header.info_name(key_idx).ok_or_else(|| BiometalError::InvalidInput {
                msg: format!("Unknown INFO index: {}", key_idx),
            })?;

            // Read INFO value
            let value = TypedValue::read(&mut cursor)?;
            info.insert(key_name.to_string(), value);
        }

        Ok(SharedData {
            chrom_id,
            pos,
            rlen,
            qual,
            n_allele,
            n_info,
            n_fmt,
            n_sample,
            id,
            reference,
            alternate,
            filter,
            info,
        })
    }

    /// Parse individual/genotype data section.
    fn parse_individual<R: Read>(
        reader: &mut R,
        header: &BcfHeader,
        shared: &SharedData,
        size: usize,
    ) -> Result<Vec<HashMap<String, TypedValue>>> {
        if size == 0 || shared.n_fmt == 0 {
            return Ok(vec![HashMap::new(); header.n_samples()]);
        }

        // Read into buffer
        let mut buf = vec![0u8; size];
        reader.read_exact(&mut buf)?;
        let mut cursor = std::io::Cursor::new(buf);

        let mut samples = vec![HashMap::new(); header.n_samples()];
        // Read FORMAT fields
        // BCF encodes FORMAT as: [key1, all_samples_data1, key2, all_samples_data2, ...]
        for _i in 0..shared.n_fmt {
            // Read FORMAT key (index)
            let key_value = TypedValue::read(&mut cursor)?;
            let key_idx = key_value.as_int().ok_or_else(|| BiometalError::InvalidInput {
                msg: format!("Expected integer for FORMAT key, got {:?}", key_value),
            })? as usize;

            let key_name =
                header
                    .format_name(key_idx)
                    .ok_or_else(|| BiometalError::InvalidInput {
                        msg: format!("Unknown FORMAT index: {}", key_idx),
                    })?;

            // Read type byte for this FORMAT field
            // For FORMAT fields, the length encodes values PER SAMPLE, not total values
            let mut type_byte_buf = [0u8; 1];
            cursor.read_exact(&mut type_byte_buf)?;
            let type_byte = type_byte_buf[0];
            let value_type = crate::io::bcf::ValueType::from_u8(type_byte)?;
            let length_per_sample = ((type_byte >> 4) & 0x0F) as usize;

            // Total values = length_per_sample * n_samples
            let total_values = length_per_sample * header.n_samples();

            // Read all values as a single array
            let mut all_values_i8 = Vec::new();
            let mut all_values_i16 = Vec::new();
            let mut all_values_i32 = Vec::new();
            let mut all_values_f32 = Vec::new();

            match value_type {
                crate::io::bcf::ValueType::Int8 => {
                    for _ in 0..total_values {
                        let mut buf = [0u8; 1];
                        cursor.read_exact(&mut buf)?;
                        all_values_i8.push(i8::from_le_bytes(buf));
                    }
                }
                crate::io::bcf::ValueType::Int16 => {
                    for _ in 0..total_values {
                        let mut buf = [0u8; 2];
                        cursor.read_exact(&mut buf)?;
                        all_values_i16.push(i16::from_le_bytes(buf));
                    }
                }
                crate::io::bcf::ValueType::Int32 => {
                    for _ in 0..total_values {
                        let mut buf = [0u8; 4];
                        cursor.read_exact(&mut buf)?;
                        all_values_i32.push(i32::from_le_bytes(buf));
                    }
                }
                crate::io::bcf::ValueType::Float => {
                    for _ in 0..total_values {
                        let mut buf = [0u8; 4];
                        cursor.read_exact(&mut buf)?;
                        all_values_f32.push(f32::from_le_bytes(buf));
                    }
                }
                crate::io::bcf::ValueType::Char => {
                    // String handling would go here
                    return Err(BiometalError::InvalidInput {
                        msg: "Char type not yet supported for FORMAT fields".to_string(),
                    });
                }
            }

            // Distribute the data to individual samples
            if !all_values_i8.is_empty() {
                // Distribute Int8 values
                if length_per_sample == 1 {
                    // Simple case: one value per sample
                    for (sample_idx, &val) in all_values_i8.iter().enumerate() {
                        samples[sample_idx].insert(key_name.to_string(), TypedValue::Int8(val));
                    }
                } else {
                    // Multi-value case (e.g., GT with 2 alleles): store as array for each sample
                    for sample_idx in 0..header.n_samples() {
                        let start = sample_idx * length_per_sample;
                        let end = start + length_per_sample;
                        let sample_values: Vec<i8> = all_values_i8[start..end].to_vec();
                        samples[sample_idx].insert(key_name.to_string(), TypedValue::Int8Array(sample_values));
                    }
                }
            } else if !all_values_i16.is_empty() {
                if length_per_sample == 1 {
                    for (sample_idx, &val) in all_values_i16.iter().enumerate() {
                        samples[sample_idx].insert(key_name.to_string(), TypedValue::Int16(val));
                    }
                } else {
                    for sample_idx in 0..header.n_samples() {
                        let start = sample_idx * length_per_sample;
                        let end = start + length_per_sample;
                        let sample_values: Vec<i16> = all_values_i16[start..end].to_vec();
                        samples[sample_idx].insert(key_name.to_string(), TypedValue::Int16Array(sample_values));
                    }
                }
            } else if !all_values_i32.is_empty() {
                if length_per_sample == 1 {
                    for (sample_idx, &val) in all_values_i32.iter().enumerate() {
                        samples[sample_idx].insert(key_name.to_string(), TypedValue::Int32(val));
                    }
                } else {
                    for sample_idx in 0..header.n_samples() {
                        let start = sample_idx * length_per_sample;
                        let end = start + length_per_sample;
                        let sample_values: Vec<i32> = all_values_i32[start..end].to_vec();
                        samples[sample_idx].insert(key_name.to_string(), TypedValue::Int32Array(sample_values));
                    }
                }
            } else if !all_values_f32.is_empty() {
                if length_per_sample == 1 {
                    for (sample_idx, &val) in all_values_f32.iter().enumerate() {
                        samples[sample_idx].insert(key_name.to_string(), TypedValue::Float(val));
                    }
                } else {
                    for sample_idx in 0..header.n_samples() {
                        let start = sample_idx * length_per_sample;
                        let end = start + length_per_sample;
                        let sample_values: Vec<f32> = all_values_f32[start..end].to_vec();
                        samples[sample_idx].insert(key_name.to_string(), TypedValue::FloatArray(sample_values));
                    }
                }
            }
        }

        Ok(samples)
    }

    /// Get chromosome/contig name from header.
    pub fn chrom<'a>(&self, header: &'a BcfHeader) -> &'a str {
        header
            .contigs()
            .get(self.shared.chrom_id as usize)
            .map(|s| s.as_str())
            .unwrap_or(".")
    }

    /// Get chromosome ID (index).
    pub fn chrom_id(&self) -> i32 {
        self.shared.chrom_id
    }

    /// Get position (1-based, VCF format).
    ///
    /// BCF stores 0-based positions, this converts to 1-based.
    pub fn pos(&self) -> i32 {
        self.shared.pos + 1
    }

    /// Get position (0-based, BCF format).
    pub fn pos_0based(&self) -> i32 {
        self.shared.pos
    }

    /// Get reference sequence end position (0-based exclusive).
    pub fn end_pos(&self) -> i32 {
        self.shared.pos + self.shared.rlen
    }

    /// Get reference allele.
    pub fn reference(&self) -> &str {
        &self.shared.reference
    }

    /// Get alternate alleles.
    pub fn alternate(&self) -> &[String] {
        &self.shared.alternate
    }

    /// Get quality score.
    pub fn qual(&self) -> Option<f32> {
        if self.shared.qual.is_nan() || self.shared.qual.to_bits() == 0x7F800001 {
            None
        } else {
            Some(self.shared.qual)
        }
    }

    /// Get variant ID.
    pub fn id(&self) -> Option<&str> {
        self.shared.id.as_deref()
    }

    /// Get FILTER status.
    ///
    /// Returns filter names from the header.
    pub fn filter<'a>(&self, header: &'a BcfHeader) -> Vec<&'a str> {
        self.shared
            .filter
            .iter()
            .filter_map(|&idx| header.filter_name(idx))
            .collect()
    }

    /// Check if variant passed all filters.
    pub fn is_pass(&self) -> bool {
        self.shared.filter.len() == 1 && self.shared.filter[0] == 0
    }

    /// Get INFO field value by name.
    pub fn info(&self, name: &str) -> Option<&TypedValue> {
        self.shared.info.get(name)
    }

    /// Get all INFO fields.
    pub fn info_fields(&self) -> &HashMap<String, TypedValue> {
        &self.shared.info
    }

    /// Get genotype data for a sample by index.
    pub fn genotype(&self, sample_idx: usize) -> Option<&HashMap<String, TypedValue>> {
        self.individual.get(sample_idx)
    }

    /// Get FORMAT field value for a specific sample.
    pub fn format(&self, sample_idx: usize, field: &str) -> Option<&TypedValue> {
        self.individual.get(sample_idx)?.get(field)
    }

    /// Get number of samples in this record.
    pub fn n_samples(&self) -> usize {
        self.individual.len()
    }

    /// Get number of alternate alleles.
    pub fn n_alleles(&self) -> usize {
        self.shared.alternate.len() + 1 // +1 for reference
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shared_data_creation() {
        let shared = SharedData {
            chrom_id: 0,
            pos: 1000,
            rlen: 1,
            qual: 30.0,
            n_allele: 2,
            n_info: 1,
            n_fmt: 1,
            n_sample: 2,
            id: Some("rs123".to_string()),
            reference: "A".to_string(),
            alternate: vec!["T".to_string()],
            filter: vec![0],
            info: HashMap::new(),
        };

        assert_eq!(shared.chrom_id, 0);
        assert_eq!(shared.pos, 1000);
        assert_eq!(shared.reference, "A");
        assert_eq!(shared.alternate, vec!["T"]);
    }
}
