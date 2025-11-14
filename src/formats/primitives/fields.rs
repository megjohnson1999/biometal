//! Field parsing utilities for tab-delimited formats.
//!
//! This module provides helper functions for parsing individual fields
//! from tab-delimited records. These utilities handle common patterns:
//! - Required vs optional fields
//! - Type conversion with error handling
//! - Missing value markers (`.`, `*`, empty string)
//! - Key=value attribute parsing (GFF, VCF INFO fields)
//!
//! # Examples
//!
//! ```
//! use biometal::formats::primitives::fields::{parse_required, parse_optional, split_fields};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Parse required field
//! let start: u64 = parse_required("12345", "start", 1)?;
//! assert_eq!(start, 12345);
//!
//! // Parse optional field (. means missing)
//! let name: Option<String> = parse_optional(".", "name", 1)?;
//! assert_eq!(name, None);
//!
//! let name: Option<String> = parse_optional("gene1", "name", 1)?;
//! assert_eq!(name, Some("gene1".to_string()));
//!
//! // Split tab-delimited line
//! let fields = split_fields("chr1\t100\t200", Some(3), 1)?;
//! assert_eq!(fields.len(), 3);
//! # Ok(())
//! # }
//! ```

use crate::formats::primitives::{FormatError, Result};
use std::collections::HashMap;
use std::str::FromStr;

/// Parses a required field with type conversion.
///
/// # Arguments
///
/// * `field` - The field string to parse
/// * `field_name` - Name of the field (for error messages)
/// * `line` - Line number (for error messages)
///
/// # Errors
///
/// Returns [`FormatError::InvalidField`] if parsing fails.
///
/// # Examples
///
/// ```
/// use biometal::formats::primitives::fields::parse_required;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let start: u64 = parse_required("12345", "start", 1)?;
/// assert_eq!(start, 12345);
///
/// let score: f64 = parse_required("3.14", "score", 1)?;
/// assert_eq!(score, 3.14);
/// # Ok(())
/// # }
/// ```
pub fn parse_required<T: FromStr>(field: &str, field_name: &str, line: usize) -> Result<T>
where
    T::Err: std::error::Error + Send + Sync + 'static,
{
    field.parse().map_err(|e: T::Err| FormatError::InvalidField {
        field: field_name.to_string(),
        line,
        reason: e.to_string(),
    })
}

/// Parses an optional field with type conversion.
///
/// Returns `None` if the field is a missing value marker:
/// - `.` (standard missing value in BED, GFF)
/// - `*` (missing value in some formats)
/// - Empty string
///
/// # Arguments
///
/// * `field` - The field string to parse
/// * `field_name` - Name of the field (for error messages)
/// * `line` - Line number (for error messages)
///
/// # Errors
///
/// Returns [`FormatError::InvalidField`] if the field is not a missing marker
/// and parsing fails.
///
/// # Examples
///
/// ```
/// use biometal::formats::primitives::fields::parse_optional;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// // Missing value markers
/// let name: Option<String> = parse_optional(".", "name", 1)?;
/// assert_eq!(name, None);
///
/// let name: Option<String> = parse_optional("*", "name", 1)?;
/// assert_eq!(name, None);
///
/// let name: Option<String> = parse_optional("", "name", 1)?;
/// assert_eq!(name, None);
///
/// // Present value
/// let score: Option<f64> = parse_optional("3.14", "score", 1)?;
/// assert_eq!(score, Some(3.14));
/// # Ok(())
/// # }
/// ```
pub fn parse_optional<T: FromStr>(field: &str, field_name: &str, line: usize) -> Result<Option<T>>
where
    T::Err: std::error::Error + Send + Sync + 'static,
{
    // Check for missing value markers
    if field == "." || field == "*" || field.is_empty() {
        return Ok(None);
    }

    field
        .parse()
        .map(Some)
        .map_err(|e: T::Err| FormatError::InvalidField {
            field: field_name.to_string(),
            line,
            reason: e.to_string(),
        })
}

/// Splits a tab-delimited line into fields with optional validation.
///
/// # Arguments
///
/// * `line` - The line to split
/// * `expected` - Expected number of fields (None = any number allowed)
/// * `line_number` - Line number (for error messages)
///
/// # Errors
///
/// Returns [`FormatError::FieldCount`] if the number of fields doesn't match
/// the expected count.
///
/// # Examples
///
/// ```
/// use biometal::formats::primitives::fields::split_fields;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// // With validation
/// let fields = split_fields("chr1\t100\t200", Some(3), 1)?;
/// assert_eq!(fields.len(), 3);
/// assert_eq!(fields[0], "chr1");
///
/// // Without validation
/// let fields = split_fields("chr1\t100\t200\tgene", None, 1)?;
/// assert_eq!(fields.len(), 4);
/// # Ok(())
/// # }
/// ```
pub fn split_fields(line: &str, expected: Option<usize>, line_number: usize) -> Result<Vec<&str>> {
    let fields: Vec<&str> = line.split('\t').collect();

    if let Some(expected_count) = expected {
        if fields.len() < expected_count {
            return Err(FormatError::FieldCount {
                expected: expected_count,
                actual: fields.len(),
                line: line_number,
            });
        }
    }

    Ok(fields)
}

/// Parses key=value attribute pairs.
///
/// Used for:
/// - GFF/GTF attributes (9th column): `ID=gene1;Name=ABC1;biotype=protein_coding`
/// - VCF INFO field: `DP=100;AF=0.5;DB`
///
/// # Format
///
/// Attributes are semicolon-separated `key=value` pairs:
/// - `key=value;key2=value2`
/// - Keys without values are stored with empty string value
/// - Whitespace is trimmed from keys and values
///
/// # Examples
///
/// ```
/// use biometal::formats::primitives::fields::parse_attributes;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let attrs = parse_attributes("ID=gene1;Name=ABC1;biotype=protein_coding", 1)?;
/// assert_eq!(attrs.get("ID"), Some(&"gene1".to_string()));
/// assert_eq!(attrs.get("Name"), Some(&"ABC1".to_string()));
/// assert_eq!(attrs.get("biotype"), Some(&"protein_coding".to_string()));
///
/// // Flag attributes (no value)
/// let attrs = parse_attributes("DP=100;DB", 1)?;
/// assert_eq!(attrs.get("DP"), Some(&"100".to_string()));
/// assert_eq!(attrs.get("DB"), Some(&"".to_string()));
/// # Ok(())
/// # }
/// ```
pub fn parse_attributes(attr_str: &str, line: usize) -> Result<HashMap<String, String>> {
    let mut attributes = HashMap::new();

    // Handle empty or missing attributes
    if attr_str.is_empty() || attr_str == "." {
        return Ok(attributes);
    }

    for pair in attr_str.split(';') {
        let pair = pair.trim();
        if pair.is_empty() {
            continue;
        }

        // Split on first '=' only
        if let Some(pos) = pair.find('=') {
            let key = pair[..pos].trim().to_string();
            let value = pair[pos + 1..].trim().to_string();
            attributes.insert(key, value);
        } else {
            // Flag attribute (no value, e.g., "DB" in VCF)
            attributes.insert(pair.to_string(), String::new());
        }
    }

    Ok(attributes)
}

/// Parses a comma-separated list of values.
///
/// Used for VCF ALT alleles, GFF parent lists, etc.
///
/// # Examples
///
/// ```
/// use biometal::formats::primitives::fields::parse_comma_list;
///
/// let values = parse_comma_list("A,T,G");
/// assert_eq!(values, vec!["A", "T", "G"]);
///
/// let values = parse_comma_list("gene1");
/// assert_eq!(values, vec!["gene1"]);
///
/// let values = parse_comma_list("");
/// assert_eq!(values.len(), 0);
/// ```
pub fn parse_comma_list(list_str: &str) -> Vec<&str> {
    if list_str.is_empty() || list_str == "." {
        return Vec::new();
    }

    list_str.split(',').map(|s| s.trim()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_required_valid() {
        let result: u64 = parse_required("12345", "test", 1).unwrap();
        assert_eq!(result, 12345);

        let result: f64 = parse_required("3.14", "test", 1).unwrap();
        assert_eq!(result, 3.14);

        let result: String = parse_required("hello", "test", 1).unwrap();
        assert_eq!(result, "hello");
    }

    #[test]
    fn test_parse_required_invalid() {
        let result: Result<u64> = parse_required("not_a_number", "test", 1);
        assert!(result.is_err());
        match result {
            Err(FormatError::InvalidField { field, line, .. }) => {
                assert_eq!(field, "test");
                assert_eq!(line, 1);
            }
            _ => panic!("Expected InvalidField error"),
        }
    }

    #[test]
    fn test_parse_optional_missing() {
        let result: Option<String> = parse_optional(".", "test", 1).unwrap();
        assert_eq!(result, None);

        let result: Option<String> = parse_optional("*", "test", 1).unwrap();
        assert_eq!(result, None);

        let result: Option<String> = parse_optional("", "test", 1).unwrap();
        assert_eq!(result, None);
    }

    #[test]
    fn test_parse_optional_present() {
        let result: Option<u64> = parse_optional("12345", "test", 1).unwrap();
        assert_eq!(result, Some(12345));

        let result: Option<String> = parse_optional("hello", "test", 1).unwrap();
        assert_eq!(result, Some("hello".to_string()));
    }

    #[test]
    fn test_split_fields_exact() {
        let fields = split_fields("chr1\t100\t200", Some(3), 1).unwrap();
        assert_eq!(fields.len(), 3);
        assert_eq!(fields[0], "chr1");
        assert_eq!(fields[1], "100");
        assert_eq!(fields[2], "200");
    }

    #[test]
    fn test_split_fields_more_than_expected() {
        // More fields than expected is OK (only checks minimum)
        let fields = split_fields("chr1\t100\t200\textra", Some(3), 1).unwrap();
        assert_eq!(fields.len(), 4);
    }

    #[test]
    fn test_split_fields_fewer_than_expected() {
        let result = split_fields("chr1\t100", Some(3), 1);
        assert!(result.is_err());
        match result {
            Err(FormatError::FieldCount {
                expected,
                actual,
                line,
            }) => {
                assert_eq!(expected, 3);
                assert_eq!(actual, 2);
                assert_eq!(line, 1);
            }
            _ => panic!("Expected FieldCount error"),
        }
    }

    #[test]
    fn test_split_fields_no_validation() {
        let fields = split_fields("chr1\t100", None, 1).unwrap();
        assert_eq!(fields.len(), 2);

        let fields = split_fields("chr1\t100\t200\t300", None, 1).unwrap();
        assert_eq!(fields.len(), 4);
    }

    #[test]
    fn test_parse_attributes_basic() {
        let attrs = parse_attributes("ID=gene1;Name=ABC1", 1).unwrap();
        assert_eq!(attrs.get("ID"), Some(&"gene1".to_string()));
        assert_eq!(attrs.get("Name"), Some(&"ABC1".to_string()));
    }

    #[test]
    fn test_parse_attributes_flag() {
        let attrs = parse_attributes("DP=100;DB", 1).unwrap();
        assert_eq!(attrs.get("DP"), Some(&"100".to_string()));
        assert_eq!(attrs.get("DB"), Some(&"".to_string()));
    }

    #[test]
    fn test_parse_attributes_whitespace() {
        let attrs = parse_attributes("ID = gene1 ; Name = ABC1", 1).unwrap();
        assert_eq!(attrs.get("ID"), Some(&"gene1".to_string()));
        assert_eq!(attrs.get("Name"), Some(&"ABC1".to_string()));
    }

    #[test]
    fn test_parse_attributes_empty() {
        let attrs = parse_attributes("", 1).unwrap();
        assert_eq!(attrs.len(), 0);

        let attrs = parse_attributes(".", 1).unwrap();
        assert_eq!(attrs.len(), 0);
    }

    #[test]
    fn test_parse_comma_list() {
        let values = parse_comma_list("A,T,G");
        assert_eq!(values, vec!["A", "T", "G"]);

        let values = parse_comma_list("single");
        assert_eq!(values, vec!["single"]);

        let values = parse_comma_list("");
        assert_eq!(values.len(), 0);

        let values = parse_comma_list(".");
        assert_eq!(values.len(), 0);
    }
}
