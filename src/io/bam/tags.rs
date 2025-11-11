//! BAM optional tags (auxiliary data).
//!
//! Optional tags store additional information about alignments such as:
//! - Edit distance (NM:i)
//! - Alignment score (AS:i)
//! - Secondary alignment status (SA:Z)
//! - Read group (RG:Z)
//! - MD string (MD:Z)
//! - Many others
//!
//! # Format
//!
//! Each tag is encoded as:
//! - 2 bytes: tag name (e.g., "NM")
//! - 1 byte: value type (A, i, f, Z, H, B)
//! - N bytes: value (format depends on type)
//!
//! # Tag Types
//!
//! - **A**: Printable character (1 byte)
//! - **i**: Signed integer (c, C, s, S, i, I - variable width)
//! - **f**: Float (4 bytes, IEEE 754)
//! - **Z**: Null-terminated string
//! - **H**: Hex string (stored as string)
//! - **B**: Array of numeric values
//!
//! # Phase 4 Implementation
//!
//! Phase 4 completes tag parsing with type-safe accessors, enabling:
//! - Individual tag extraction by name
//! - Type-safe value access
//! - Tag iteration
//! - Common tag convenience methods

use std::io;
use std::fmt;
use super::error::BamDecodeError;

/// Tag value types in BAM format.
///
/// Each tag has a type code and corresponding value encoding.
#[derive(Debug, Clone, PartialEq)]
pub enum TagValue {
    /// Character (A): Single printable character
    Char(u8),
    /// Integer (i): Signed integer (variable width)
    Int(i64),
    /// Float (f): IEEE 754 single-precision float
    Float(f32),
    /// String (Z): Null-terminated string
    String(String),
    /// Hex string (H): Hex-encoded data
    Hex(String),
    /// Array (B): Typed array of numbers
    Array(ArrayValue),
}

/// Array value types for tag arrays (B type).
#[derive(Debug, Clone, PartialEq)]
pub enum ArrayValue {
    /// Array of signed 8-bit integers
    Int8(Vec<i8>),
    /// Array of unsigned 8-bit integers
    UInt8(Vec<u8>),
    /// Array of signed 16-bit integers
    Int16(Vec<i16>),
    /// Array of unsigned 16-bit integers
    UInt16(Vec<u16>),
    /// Array of signed 32-bit integers
    Int32(Vec<i32>),
    /// Array of unsigned 32-bit integers
    UInt32(Vec<u32>),
    /// Array of 32-bit floats
    Float(Vec<f32>),
}

/// A single BAM tag with name and value.
#[derive(Debug, Clone, PartialEq)]
pub struct Tag {
    /// Two-character tag name (e.g., "NM", "AS", "RG")
    pub name: [u8; 2],
    /// Tag value
    pub value: TagValue,
}

impl Tag {
    /// Get tag name as a string slice.
    ///
    /// # Example
    ///
    /// ```
    /// # use biometal::io::bam::Tag;
    /// # use biometal::io::bam::TagValue;
    /// let tag = Tag {
    ///     name: *b"NM",
    ///     value: TagValue::Int(5),
    /// };
    /// assert_eq!(tag.name_str(), "NM");
    /// ```
    pub fn name_str(&self) -> &str {
        std::str::from_utf8(&self.name).unwrap_or("??")
    }
}

impl fmt::Display for Tag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:", self.name_str())?;
        match &self.value {
            TagValue::Char(c) => write!(f, "A:{}", *c as char),
            TagValue::Int(i) => write!(f, "i:{}", i),
            TagValue::Float(fl) => write!(f, "f:{}", fl),
            TagValue::String(s) => write!(f, "Z:{}", s),
            TagValue::Hex(h) => write!(f, "H:{}", h),
            TagValue::Array(arr) => {
                write!(f, "B:")?;
                match arr {
                    ArrayValue::Int8(v) => write!(f, "c,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::UInt8(v) => write!(f, "C,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::Int16(v) => write!(f, "s,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::UInt16(v) => write!(f, "S,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::Int32(v) => write!(f, "i,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::UInt32(v) => write!(f, "I,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                    ArrayValue::Float(v) => write!(f, "f,{}", v.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(",")),
                }
            }
        }
    }
}

/// Container for BAM optional tags.
///
/// # Phase 4 Implementation
///
/// Now provides full tag parsing with type-safe accessors:
/// - Individual tag extraction by name
/// - Tag iteration
/// - Type-safe value access
/// - Common tag convenience methods
///
/// # Example
///
/// ```
/// # use biometal::io::bam::Tags;
/// # use biometal::io::bam::TagValue;
/// // Create tags from raw BAM data
/// let tags = Tags::from_raw(vec![
///     // NM:i:5 (edit distance)
///     b'N', b'M', b'i', 5, 0, 0, 0,
/// ]);
///
/// // Access tag by name
/// if let Some(tag) = tags.get(b"NM").unwrap() {
///     if let TagValue::Int(nm) = tag.value {
///         assert_eq!(nm, 5);
///     }
/// }
/// ```
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Tags {
    /// Raw tag data from BAM record.
    ///
    /// Preserves all tag information for round-trip validation.
    /// Tags are parsed on-demand when accessed.
    data: Vec<u8>,
}

impl Tags {
    /// Create empty tags.
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }

    /// Create tags from raw BAM data.
    ///
    /// # Arguments
    ///
    /// * `data` - Raw tag bytes from BAM record
    ///
    /// # Phase 1 Note
    ///
    /// Currently just stores the raw bytes. Full parsing will be
    /// added in Phase 6 when we need tag-level access.
    pub fn from_raw(data: Vec<u8>) -> Self {
        Self { data }
    }

    /// Get the raw tag data.
    ///
    /// Useful for:
    /// - Differential testing (compare raw bytes)
    /// - Round-trip testing (preserve exact encoding)
    /// - Custom tag parsers
    pub fn as_raw(&self) -> &[u8] {
        &self.data
    }

    /// Check if tags are empty.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get the total size of tag data in bytes.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Get a specific tag by name.
    ///
    /// Returns `Ok(Some(tag))` if found, `Ok(None)` if not found, or `Err` if parsing fails.
    ///
    /// # Example
    ///
    /// ```
    /// # use biometal::io::bam::{Tags, TagValue};
    /// let tags = Tags::from_raw(vec![b'N', b'M', b'i', 5, 0, 0, 0]);
    ///
    /// let nm_tag = tags.get(b"NM").unwrap().unwrap();
    /// assert_eq!(nm_tag.name_str(), "NM");
    /// if let TagValue::Int(val) = nm_tag.value {
    ///     assert_eq!(val, 5);
    /// }
    /// ```
    pub fn get(&self, name: &[u8; 2]) -> io::Result<Option<Tag>> {
        let mut cursor = 0;

        while cursor < self.data.len() {
            // Read tag name (2 bytes)
            if cursor + 2 > self.data.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Incomplete tag name at offset {}", cursor),
                ));
            }
            let tag_name = [self.data[cursor], self.data[cursor + 1]];
            cursor += 2;

            // Read tag type (1 byte)
            if cursor >= self.data.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Missing tag type at offset {}", cursor),
                ));
            }
            let tag_type = self.data[cursor];
            cursor += 1;

            // Parse value based on type
            let (value, value_size) = parse_tag_value(&self.data[cursor..], tag_type)?;

            // Check if this is the tag we're looking for
            if &tag_name == name {
                return Ok(Some(Tag {
                    name: tag_name,
                    value,
                }));
            }

            cursor += value_size;
        }

        Ok(None)
    }

    /// Iterate over all tags.
    ///
    /// # Example
    ///
    /// ```
    /// # use biometal::io::bam::Tags;
    /// let tags = Tags::from_raw(vec![
    ///     b'N', b'M', b'i', 5, 0, 0, 0,
    ///     b'A', b'S', b'i', 100, 0, 0, 0,
    /// ]);
    ///
    /// for tag in tags.iter().unwrap() {
    ///     println!("{}: {:?}", tag.name_str(), tag.value);
    /// }
    /// ```
    pub fn iter(&self) -> io::Result<Vec<Tag>> {
        let mut tags = Vec::new();
        let mut seen_tags = std::collections::HashSet::new();
        let mut cursor = 0;

        while cursor < self.data.len() {
            // Read tag name (2 bytes)
            if cursor + 2 > self.data.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Incomplete tag name at offset {}", cursor),
                ));
            }
            let tag_name = [self.data[cursor], self.data[cursor + 1]];
            cursor += 2;

            // Check for duplicate tags (BAM spec violation)
            if !seen_tags.insert(tag_name) {
                return Err(BamDecodeError::DuplicateTag { tag: tag_name }.into());
            }

            // Read tag type (1 byte)
            if cursor >= self.data.len() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("Missing tag type at offset {}", cursor),
                ));
            }
            let tag_type = self.data[cursor];
            cursor += 1;

            // Parse value
            let (value, value_size) = parse_tag_value(&self.data[cursor..], tag_type)?;
            cursor += value_size;

            tags.push(Tag {
                name: tag_name,
                value,
            });
        }

        Ok(tags)
    }

    /// Get an integer tag value.
    ///
    /// Convenience method that returns the value if the tag exists and is an integer.
    pub fn get_i32(&self, name: &[u8; 2]) -> io::Result<Option<i32>> {
        match self.get(name)? {
            Some(tag) => match tag.value {
                TagValue::Int(i) => Ok(Some(i as i32)),
                _ => Ok(None),
            },
            None => Ok(None),
        }
    }

    /// Get a string tag value.
    ///
    /// Convenience method that returns the value if the tag exists and is a string.
    pub fn get_string(&self, name: &[u8; 2]) -> io::Result<Option<String>> {
        match self.get(name)? {
            Some(tag) => match tag.value {
                TagValue::String(s) => Ok(Some(s)),
                _ => Ok(None),
            },
            None => Ok(None),
        }
    }

    /// Get a character tag value.
    ///
    /// Convenience method that returns the value if the tag exists and is a character.
    ///
    /// # Example
    ///
    /// ```
    /// # use biometal::io::bam::Tags;
    /// let data = vec![b'C', b'C', b'A', b'X']; // CC:A:X
    /// let tags = Tags::from_raw(data);
    /// assert_eq!(tags.get_char(b"CC").unwrap(), Some(b'X'));
    /// ```
    pub fn get_char(&self, name: &[u8; 2]) -> io::Result<Option<u8>> {
        match self.get(name)? {
            Some(tag) => match tag.value {
                TagValue::Char(c) => Ok(Some(c)),
                _ => Ok(None),
            },
            None => Ok(None),
        }
    }

    /// Get a float tag value.
    ///
    /// Convenience method that returns the value if the tag exists and is a float.
    ///
    /// # Example
    ///
    /// ```
    /// # use biometal::io::bam::Tags;
    /// let data = vec![b'F', b'V', b'f', 0x00, 0x00, 0x80, 0x3F]; // FV:f:1.0
    /// let tags = Tags::from_raw(data);
    /// assert_eq!(tags.get_float(b"FV").unwrap(), Some(1.0));
    /// ```
    pub fn get_float(&self, name: &[u8; 2]) -> io::Result<Option<f32>> {
        match self.get(name)? {
            Some(tag) => match tag.value {
                TagValue::Float(f) => Ok(Some(f)),
                _ => Ok(None),
            },
            None => Ok(None),
        }
    }

    /// Get a hex string tag value.
    ///
    /// Convenience method that returns the value if the tag exists and is a hex string.
    pub fn get_hex(&self, name: &[u8; 2]) -> io::Result<Option<String>> {
        match self.get(name)? {
            Some(tag) => match tag.value {
                TagValue::Hex(h) => Ok(Some(h)),
                _ => Ok(None),
            },
            None => Ok(None),
        }
    }

    /// Get an array tag value.
    ///
    /// Convenience method that returns the value if the tag exists and is an array.
    ///
    /// # Example
    ///
    /// ```
    /// # use biometal::io::bam::{Tags, ArrayValue};
    /// // Create array tag: TS:B:I:3,1,2,3 (UInt32 array)
    /// let data = vec![
    ///     b'T', b'S', b'B', b'I',      // Tag name, type, array type
    ///     0x03, 0x00, 0x00, 0x00,      // count = 3
    ///     0x01, 0x00, 0x00, 0x00,      // [0] = 1
    ///     0x02, 0x00, 0x00, 0x00,      // [1] = 2
    ///     0x03, 0x00, 0x00, 0x00,      // [2] = 3
    /// ];
    /// let tags = Tags::from_raw(data);
    /// if let Some(ArrayValue::UInt32(values)) = tags.get_array(b"TS").unwrap() {
    ///     assert_eq!(values, vec![1, 2, 3]);
    /// }
    /// ```
    pub fn get_array(&self, name: &[u8; 2]) -> io::Result<Option<ArrayValue>> {
        match self.get(name)? {
            Some(tag) => match tag.value {
                TagValue::Array(arr) => Ok(Some(arr)),
                _ => Ok(None),
            },
            None => Ok(None),
        }
    }

    // ==================================================================
    // Common BAM Tag Convenience Methods
    // ==================================================================

    /// Get edit distance (NM tag).
    ///
    /// Returns the number of differences between the reference and query sequence.
    /// This is a standard SAM optional tag.
    ///
    /// # Example
    ///
    /// ```
    /// # use biometal::io::bam::Tags;
    /// let data = vec![b'N', b'M', b'i', 0x05, 0x00, 0x00, 0x00]; // NM:i:5
    /// let tags = Tags::from_raw(data);
    /// assert_eq!(tags.edit_distance().unwrap(), Some(5));
    /// ```
    pub fn edit_distance(&self) -> io::Result<Option<i32>> {
        self.get_i32(b"NM")
    }

    /// Get alignment score (AS tag).
    ///
    /// Returns the aligner-specific alignment score.
    pub fn alignment_score(&self) -> io::Result<Option<i32>> {
        self.get_i32(b"AS")
    }

    /// Get mapping quality (MQ tag).
    ///
    /// Returns the mapping quality of the mate/next segment.
    pub fn mate_mapping_quality(&self) -> io::Result<Option<i32>> {
        self.get_i32(b"MQ")
    }

    /// Get read group (RG tag).
    ///
    /// Returns the read group identifier.
    ///
    /// # Example
    ///
    /// ```
    /// # use biometal::io::bam::Tags;
    /// let data = vec![b'R', b'G', b'Z', b'r', b'g', b'1', 0x00]; // RG:Z:rg1
    /// let tags = Tags::from_raw(data);
    /// assert_eq!(tags.read_group().unwrap(), Some("rg1".to_string()));
    /// ```
    pub fn read_group(&self) -> io::Result<Option<String>> {
        self.get_string(b"RG")
    }

    /// Get MD string (MD tag).
    ///
    /// Returns the string encoding mismatched and deleted reference bases.
    pub fn md_string(&self) -> io::Result<Option<String>> {
        self.get_string(b"MD")
    }

    /// Get barcode sequence (BC tag).
    ///
    /// Returns the barcode sequence identifying the sample.
    pub fn barcode(&self) -> io::Result<Option<String>> {
        self.get_string(b"BC")
    }

    /// Get UMI (unique molecular identifier) sequence (RX tag).
    ///
    /// Returns the sequence bases of the UMI.
    pub fn umi(&self) -> io::Result<Option<String>> {
        self.get_string(b"RX")
    }

    /// Get comment (CO tag).
    ///
    /// Returns any free-text comment.
    pub fn comment(&self) -> io::Result<Option<String>> {
        self.get_string(b"CO")
    }
}

/// Parse tags from BAM record data.
///
/// # Phase 1 Implementation
///
/// Currently just copies the raw bytes. Full tag parsing will be
/// added in Phase 6.
///
/// # Arguments
///
/// * `data` - Raw tag data from end of BAM record
///
/// # Returns
///
/// Tags structure containing the raw data.
///
/// # Future Work
///
/// - Validate tag format (2-char name, type code, value)
/// - Parse individual tags into structured types
/// - Provide type-safe accessors
///
/// # Example
///
/// ```
/// use biometal::io::bam::Tags;
///
/// // Empty tags
/// let data = vec![];
/// let tags = Tags::from_raw(data);
/// assert!(tags.is_empty());
/// ```
pub fn parse_tags(data: &[u8]) -> io::Result<Tags> {
    // Phase 4: Store raw bytes (parsing happens on-demand via get/iter)
    Ok(Tags::from_raw(data.to_vec()))
}

/// Parse a single tag value from raw bytes.
///
/// Returns `(TagValue, bytes_consumed)`.
///
/// # Tag Type Codes
///
/// - `A` (65): Character (1 byte)
/// - `c` (99): Int8 (1 byte signed)
/// - `C` (67): UInt8 (1 byte unsigned)
/// - `s` (115): Int16 (2 bytes signed)
/// - `S` (83): UInt16 (2 bytes unsigned)
/// - `i` (105): Int32 (4 bytes signed)
/// - `I` (73): UInt32 (4 bytes unsigned)
/// - `f` (102): Float (4 bytes)
/// - `Z` (90): Null-terminated string
/// - `H` (72): Hex string (null-terminated)
/// - `B` (66): Typed array
fn parse_tag_value(data: &[u8], type_code: u8) -> io::Result<(TagValue, usize)> {
    match type_code {
        b'A' => {
            // Character (1 byte)
            if data.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for character tag",
                ));
            }
            Ok((TagValue::Char(data[0]), 1))
        }
        b'c' => {
            // Int8 (1 byte signed)
            if data.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for int8 tag",
                ));
            }
            Ok((TagValue::Int(data[0] as i8 as i64), 1))
        }
        b'C' => {
            // UInt8 (1 byte unsigned)
            if data.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for uint8 tag",
                ));
            }
            Ok((TagValue::Int(data[0] as i64), 1))
        }
        b's' => {
            // Int16 (2 bytes signed, little-endian)
            if data.len() < 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for int16 tag",
                ));
            }
            let value = i16::from_le_bytes([data[0], data[1]]);
            Ok((TagValue::Int(value as i64), 2))
        }
        b'S' => {
            // UInt16 (2 bytes unsigned, little-endian)
            if data.len() < 2 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for uint16 tag",
                ));
            }
            let value = u16::from_le_bytes([data[0], data[1]]);
            Ok((TagValue::Int(value as i64), 2))
        }
        b'i' => {
            // Int32 (4 bytes signed, little-endian)
            if data.len() < 4 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for int32 tag",
                ));
            }
            let value = i32::from_le_bytes([data[0], data[1], data[2], data[3]]);
            Ok((TagValue::Int(value as i64), 4))
        }
        b'I' => {
            // UInt32 (4 bytes unsigned, little-endian)
            if data.len() < 4 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for uint32 tag",
                ));
            }
            let value = u32::from_le_bytes([data[0], data[1], data[2], data[3]]);
            Ok((TagValue::Int(value as i64), 4))
        }
        b'f' => {
            // Float (4 bytes, little-endian)
            if data.len() < 4 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for float tag",
                ));
            }
            let value = f32::from_le_bytes([data[0], data[1], data[2], data[3]]);
            Ok((TagValue::Float(value), 4))
        }
        b'Z' | b'H' => {
            // Null-terminated string or hex string
            let null_pos = data.iter().position(|&b| b == 0).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Missing null terminator in string tag",
                )
            })?;

            let string = String::from_utf8(data[..null_pos].to_vec()).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "Invalid UTF-8 in string tag")
            })?;

            let value = if type_code == b'Z' {
                TagValue::String(string)
            } else {
                TagValue::Hex(string)
            };

            Ok((value, null_pos + 1))
        }
        b'B' => {
            // Typed array
            if data.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for array tag type",
                ));
            }

            let array_type = data[0];
            if data.len() < 5 {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "Insufficient data for array tag count",
                ));
            }

            let count_u32 = u32::from_le_bytes([data[1], data[2], data[3], data[4]]);
            let count = usize::try_from(count_u32)
                .map_err(|_| BamDecodeError::ArrayCountOverflow { count: count_u32 })?;
            let mut offset = 5;

            let array_value = match array_type {
                b'c' => {
                    // Int8 array
                    if data.len() < offset + count {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for int8 array",
                        ));
                    }
                    let values: Vec<i8> = data[offset..offset + count]
                        .iter()
                        .map(|&b| b as i8)
                        .collect();
                    offset += count;
                    ArrayValue::Int8(values)
                }
                b'C' => {
                    // UInt8 array
                    if data.len() < offset + count {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for uint8 array",
                        ));
                    }
                    let values = data[offset..offset + count].to_vec();
                    offset += count;
                    ArrayValue::UInt8(values)
                }
                b's' => {
                    // Int16 array
                    let bytes_needed = count.checked_mul(2).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Array size overflow: count={}", count),
                        )
                    })?;
                    if data.len() < offset + bytes_needed {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for int16 array",
                        ));
                    }
                    let values: Vec<i16> = (0..count)
                        .map(|i| {
                            let idx = offset + i * 2;
                            i16::from_le_bytes([data[idx], data[idx + 1]])
                        })
                        .collect();
                    offset += bytes_needed;
                    ArrayValue::Int16(values)
                }
                b'S' => {
                    // UInt16 array
                    let bytes_needed = count.checked_mul(2).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Array size overflow: count={}", count),
                        )
                    })?;
                    if data.len() < offset + bytes_needed {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for uint16 array",
                        ));
                    }
                    let values: Vec<u16> = (0..count)
                        .map(|i| {
                            let idx = offset + i * 2;
                            u16::from_le_bytes([data[idx], data[idx + 1]])
                        })
                        .collect();
                    offset += bytes_needed;
                    ArrayValue::UInt16(values)
                }
                b'i' => {
                    // Int32 array
                    let bytes_needed = count.checked_mul(4).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Array size overflow: count={}", count),
                        )
                    })?;
                    if data.len() < offset + bytes_needed {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for int32 array",
                        ));
                    }
                    let values: Vec<i32> = (0..count)
                        .map(|i| {
                            let idx = offset + i * 4;
                            i32::from_le_bytes([data[idx], data[idx + 1], data[idx + 2], data[idx + 3]])
                        })
                        .collect();
                    offset += bytes_needed;
                    ArrayValue::Int32(values)
                }
                b'I' => {
                    // UInt32 array
                    let bytes_needed = count.checked_mul(4).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Array size overflow: count={}", count),
                        )
                    })?;
                    if data.len() < offset + bytes_needed {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for uint32 array",
                        ));
                    }
                    let values: Vec<u32> = (0..count)
                        .map(|i| {
                            let idx = offset + i * 4;
                            u32::from_le_bytes([data[idx], data[idx + 1], data[idx + 2], data[idx + 3]])
                        })
                        .collect();
                    offset += bytes_needed;
                    ArrayValue::UInt32(values)
                }
                b'f' => {
                    // Float array
                    let bytes_needed = count.checked_mul(4).ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("Array size overflow: count={}", count),
                        )
                    })?;
                    if data.len() < offset + bytes_needed {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            "Insufficient data for float array",
                        ));
                    }
                    let values: Vec<f32> = (0..count)
                        .map(|i| {
                            let idx = offset + i * 4;
                            f32::from_le_bytes([data[idx], data[idx + 1], data[idx + 2], data[idx + 3]])
                        })
                        .collect();
                    offset += bytes_needed;
                    ArrayValue::Float(values)
                }
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("Invalid array type code: {}", array_type),
                    ));
                }
            };

            Ok((TagValue::Array(array_value), offset))
        }
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("Invalid tag type code: {}", type_code),
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_tags() {
        let tags = Tags::new();
        assert!(tags.is_empty());
        assert_eq!(tags.len(), 0);
        assert_eq!(tags.as_raw(), &[] as &[u8]);
    }

    #[test]
    fn test_from_raw() {
        let data = vec![0x4E, 0x4D, 0x69, 0x01, 0x00, 0x00, 0x00]; // NM:i:1
        let tags = Tags::from_raw(data.clone());
        assert!(!tags.is_empty());
        assert_eq!(tags.len(), 7);
        assert_eq!(tags.as_raw(), &data);
    }

    #[test]
    fn test_parse_tags_empty() {
        let data = vec![];
        let tags = parse_tags(&data).unwrap();
        assert!(tags.is_empty());
    }

    #[test]
    fn test_parse_tags_preserves_data() {
        let data = vec![0x4E, 0x4D, 0x69, 0x01, 0x00, 0x00, 0x00]; // NM:i:1
        let tags = parse_tags(&data).unwrap();
        assert_eq!(tags.as_raw(), &data);
    }

    #[test]
    fn test_tags_clone() {
        let tags1 = Tags::from_raw(vec![1, 2, 3]);
        let tags2 = tags1.clone();
        assert_eq!(tags1, tags2);
    }

    #[test]
    fn test_tags_equality() {
        let tags1 = Tags::from_raw(vec![1, 2, 3]);
        let tags2 = Tags::from_raw(vec![1, 2, 3]);
        let tags3 = Tags::from_raw(vec![1, 2, 4]);

        assert_eq!(tags1, tags2);
        assert_ne!(tags1, tags3);
    }

    // Noodles-inspired robustness tests

    #[test]
    fn test_duplicate_tags_rejected() {
        let data = vec![
            b'N', b'M', b'i', 5, 0, 0, 0,  // NM:i:5
            b'N', b'M', b'i', 3, 0, 0, 0,  // NM:i:3 (duplicate!)
        ];
        let tags = Tags::from_raw(data);
        let result = tags.iter();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Duplicate tag"));
    }

    #[test]
    fn test_array_count_overflow() {
        // This test is primarily for 32-bit platforms, but good for completeness
        #[cfg(target_pointer_width = "32")]
        {
            let data = vec![
                b'T', b'S', b'B', b'I',  // TS:B:I (uint32 array)
                0xFF, 0xFF, 0xFF, 0xFF,  // count = u32::MAX (overflows usize on 32-bit)
            ];
            let tags = Tags::from_raw(data);
            let result = tags.iter();
            assert!(result.is_err());
            assert!(result.unwrap_err().to_string().contains("too large"));
        }
    }

    #[test]
    fn test_string_tag_missing_nul() {
        let data = vec![
            b'R', b'G', b'Z',  // RG:Z:...
            b'r', b'g', b'0',  // "rg0" (no NUL terminator!)
        ];
        let tags = Tags::from_raw(data);
        let result = tags.iter();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("null terminator"));
    }

    #[test]
    fn test_truncated_array() {
        let data = vec![
            b'T', b'S', b'B', b'I',  // TS:B:I (int32 array)
            0x03, 0x00, 0x00, 0x00,  // count = 3
            0x01, 0x00, 0x00, 0x00,  // [0] = 1
            0x02, 0x00, // [1] = incomplete! (only 2 bytes of 4)
        ];
        let tags = Tags::from_raw(data);
        let result = tags.iter();
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Insufficient"));
    }

    #[test]
    fn test_array_size_overflow() {
        // Test that count * element_size overflow is detected
        // Use a count that when multiplied by 4 (for i32) would overflow usize
        // On 64-bit: usize::MAX/4 + 1 would overflow when multiplied by 4
        // On 32-bit: use smaller value that still overflows
        let huge_count = if cfg!(target_pointer_width = "64") {
            ((usize::MAX / 4) as u64 + 1).min(u32::MAX as u64) as u32
        } else {
            u32::MAX / 2
        };

        let mut data = vec![
            b'A', b'R', b'B', b'I',  // AR:B:I (int32 array)
        ];
        // Huge count that would overflow when multiplied by 4
        data.extend_from_slice(&huge_count.to_le_bytes());

        let tags = Tags::from_raw(data);
        let result = tags.iter();
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("overflow") || err_msg.contains("Insufficient"),
                "Expected overflow or insufficient data error, got: {}", err_msg);
    }

    #[test]
    fn test_array_size_overflow_int16() {
        // Test overflow for Int16 arrays (2 bytes per element)
        // Use a count that when multiplied by 2 would overflow usize
        let huge_count = if cfg!(target_pointer_width = "64") {
            ((usize::MAX / 2) as u64 + 1).min(u32::MAX as u64) as u32
        } else {
            u32::MAX / 2 + 1
        };

        let mut data = vec![
            b'A', b'R', b'B', b's',  // AR:B:s (int16 array)
        ];
        data.extend_from_slice(&huge_count.to_le_bytes());

        let tags = Tags::from_raw(data);
        let result = tags.iter();
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("overflow") || err_msg.contains("Insufficient"),
                "Expected overflow or insufficient data error, got: {}", err_msg);
    }

    #[test]
    fn test_array_size_overflow_float() {
        // Test overflow for Float arrays (4 bytes per element)
        // Use a count that when multiplied by 4 would overflow usize
        let huge_count = if cfg!(target_pointer_width = "64") {
            ((usize::MAX / 4) as u64 + 1).min(u32::MAX as u64) as u32
        } else {
            u32::MAX / 2
        };

        let mut data = vec![
            b'F', b'L', b'B', b'f',  // FL:B:f (float array)
        ];
        data.extend_from_slice(&huge_count.to_le_bytes());

        let tags = Tags::from_raw(data);
        let result = tags.iter();
        assert!(result.is_err());
        let err_msg = result.unwrap_err().to_string();
        assert!(err_msg.contains("overflow") || err_msg.contains("Insufficient"),
                "Expected overflow or insufficient data error, got: {}", err_msg);
    }

    // Tests for new accessor methods

    #[test]
    fn test_get_char() {
        let data = vec![b'C', b'C', b'A', b'X']; // CC:A:X
        let tags = Tags::from_raw(data);
        assert_eq!(tags.get_char(b"CC").unwrap(), Some(b'X'));
        assert_eq!(tags.get_char(b"XX").unwrap(), None); // Non-existent tag
    }

    #[test]
    fn test_get_float() {
        let data = vec![b'F', b'V', b'f', 0x00, 0x00, 0x80, 0x3F]; // FV:f:1.0
        let tags = Tags::from_raw(data);
        assert_eq!(tags.get_float(b"FV").unwrap(), Some(1.0));
        assert_eq!(tags.get_float(b"XX").unwrap(), None);
    }

    #[test]
    fn test_get_hex() {
        let data = vec![b'H', b'X', b'H', b'A', b'B', b'C', b'D', 0x00]; // HX:H:ABCD
        let tags = Tags::from_raw(data);
        assert_eq!(tags.get_hex(b"HX").unwrap(), Some("ABCD".to_string()));
    }

    #[test]
    fn test_get_array_uint32() {
        // TS:B:I:3,1,2,3 (UInt32 array)
        let data = vec![
            b'T', b'S', b'B', b'I',      // Tag name, type, array type
            0x03, 0x00, 0x00, 0x00,      // count = 3
            0x01, 0x00, 0x00, 0x00,      // [0] = 1
            0x02, 0x00, 0x00, 0x00,      // [1] = 2
            0x03, 0x00, 0x00, 0x00,      // [2] = 3
        ];
        let tags = Tags::from_raw(data);
        if let Some(ArrayValue::UInt32(values)) = tags.get_array(b"TS").unwrap() {
            assert_eq!(values, vec![1, 2, 3]);
        } else {
            panic!("Expected UInt32 array");
        }
    }

    #[test]
    fn test_get_array_float() {
        // FL:B:f:2,1.5,2.5 (Float array)
        let data = vec![
            b'F', b'L', b'B', b'f',      // Tag name, type, array type
            0x02, 0x00, 0x00, 0x00,      // count = 2
            0x00, 0x00, 0xC0, 0x3F,      // [0] = 1.5
            0x00, 0x00, 0x20, 0x40,      // [1] = 2.5
        ];
        let tags = Tags::from_raw(data);
        if let Some(ArrayValue::Float(values)) = tags.get_array(b"FL").unwrap() {
            assert_eq!(values.len(), 2);
            assert!((values[0] - 1.5).abs() < 0.001);
            assert!((values[1] - 2.5).abs() < 0.001);
        } else {
            panic!("Expected Float array");
        }
    }

    #[test]
    fn test_edit_distance() {
        let data = vec![b'N', b'M', b'i', 0x05, 0x00, 0x00, 0x00]; // NM:i:5
        let tags = Tags::from_raw(data);
        assert_eq!(tags.edit_distance().unwrap(), Some(5));
    }

    #[test]
    fn test_alignment_score() {
        let data = vec![b'A', b'S', b'i', 0x64, 0x00, 0x00, 0x00]; // AS:i:100
        let tags = Tags::from_raw(data);
        assert_eq!(tags.alignment_score().unwrap(), Some(100));
    }

    #[test]
    fn test_read_group() {
        let data = vec![b'R', b'G', b'Z', b'r', b'g', b'1', 0x00]; // RG:Z:rg1
        let tags = Tags::from_raw(data);
        assert_eq!(tags.read_group().unwrap(), Some("rg1".to_string()));
    }

    #[test]
    fn test_md_string() {
        let data = vec![b'M', b'D', b'Z', b'1', b'0', b'0', 0x00]; // MD:Z:100
        let tags = Tags::from_raw(data);
        assert_eq!(tags.md_string().unwrap(), Some("100".to_string()));
    }

    #[test]
    fn test_barcode() {
        let data = vec![b'B', b'C', b'Z', b'A', b'C', b'G', b'T', 0x00]; // BC:Z:ACGT
        let tags = Tags::from_raw(data);
        assert_eq!(tags.barcode().unwrap(), Some("ACGT".to_string()));
    }

    #[test]
    fn test_umi() {
        let data = vec![b'R', b'X', b'Z', b'G', b'A', b'T', b'C', 0x00]; // RX:Z:GATC
        let tags = Tags::from_raw(data);
        assert_eq!(tags.umi().unwrap(), Some("GATC".to_string()));
    }

    #[test]
    fn test_multiple_tags_with_accessors() {
        // Multiple tags: NM:i:5, AS:i:100, RG:Z:rg1
        let data = vec![
            b'N', b'M', b'i', 0x05, 0x00, 0x00, 0x00,  // NM:i:5
            b'A', b'S', b'i', 0x64, 0x00, 0x00, 0x00,  // AS:i:100
            b'R', b'G', b'Z', b'r', b'g', b'1', 0x00,  // RG:Z:rg1
        ];
        let tags = Tags::from_raw(data);

        assert_eq!(tags.edit_distance().unwrap(), Some(5));
        assert_eq!(tags.alignment_score().unwrap(), Some(100));
        assert_eq!(tags.read_group().unwrap(), Some("rg1".to_string()));
    }

    #[test]
    fn test_accessor_returns_none_for_wrong_type() {
        // Tag exists but is wrong type
        let data = vec![b'N', b'M', b'Z', b'x', b'y', b'z', 0x00]; // NM:Z:xyz (should be int)
        let tags = Tags::from_raw(data);

        // get_i32 should return None for string tag
        assert_eq!(tags.get_i32(b"NM").unwrap(), None);

        // get_string should work
        assert_eq!(tags.get_string(b"NM").unwrap(), Some("xyz".to_string()));
    }

    #[test]
    fn test_all_array_types() {
        // Test Int8 array
        {
            let data = vec![
                b'A', b'1', b'B', b'c',  // A1:B:c (Int8 array)
                0x02, 0x00, 0x00, 0x00,  // count = 2
                0xFF, 0x01,               // [-1, 1]
            ];
            let tags = Tags::from_raw(data);
            assert!(matches!(tags.get_array(b"A1").unwrap(), Some(ArrayValue::Int8(_))));
        }

        // Test UInt8 array
        {
            let data = vec![
                b'A', b'2', b'B', b'C',  // A2:B:C (UInt8 array)
                0x02, 0x00, 0x00, 0x00,  // count = 2
                0x01, 0x02,               // [1, 2]
            ];
            let tags = Tags::from_raw(data);
            assert!(matches!(tags.get_array(b"A2").unwrap(), Some(ArrayValue::UInt8(_))));
        }

        // Test Int16 array
        {
            let data = vec![
                b'A', b'3', b'B', b's',  // A3:B:s (Int16 array)
                0x02, 0x00, 0x00, 0x00,  // count = 2
                0xFF, 0xFF, 0x01, 0x00,   // [-1, 1]
            ];
            let tags = Tags::from_raw(data);
            assert!(matches!(tags.get_array(b"A3").unwrap(), Some(ArrayValue::Int16(_))));
        }

        // Test UInt16 array
        {
            let data = vec![
                b'A', b'4', b'B', b'S',  // A4:B:S (UInt16 array)
                0x02, 0x00, 0x00, 0x00,  // count = 2
                0x01, 0x00, 0x02, 0x00,   // [1, 2]
            ];
            let tags = Tags::from_raw(data);
            assert!(matches!(tags.get_array(b"A4").unwrap(), Some(ArrayValue::UInt16(_))));
        }

        // Test Int32 array
        {
            let data = vec![
                b'A', b'5', b'B', b'i',  // A5:B:i (Int32 array)
                0x02, 0x00, 0x00, 0x00,  // count = 2
                0xFF, 0xFF, 0xFF, 0xFF,   // -1
                0x01, 0x00, 0x00, 0x00,   // 1
            ];
            let tags = Tags::from_raw(data);
            assert!(matches!(tags.get_array(b"A5").unwrap(), Some(ArrayValue::Int32(_))));
        }

        // Test UInt32 array
        {
            let data = vec![
                b'A', b'6', b'B', b'I',  // A6:B:I (UInt32 array)
                0x02, 0x00, 0x00, 0x00,  // count = 2
                0x01, 0x00, 0x00, 0x00,   // 1
                0x02, 0x00, 0x00, 0x00,   // 2
            ];
            let tags = Tags::from_raw(data);
            assert!(matches!(tags.get_array(b"A6").unwrap(), Some(ArrayValue::UInt32(_))));
        }

        // Test Float array
        {
            let data = vec![
                b'A', b'7', b'B', b'f',  // A7:B:f (Float array)
                0x02, 0x00, 0x00, 0x00,  // count = 2
                0x00, 0x00, 0x80, 0x3F,   // 1.0
                0x00, 0x00, 0x00, 0x40,   // 2.0
            ];
            let tags = Tags::from_raw(data);
            assert!(matches!(tags.get_array(b"A7").unwrap(), Some(ArrayValue::Float(_))));
        }
    }
}
