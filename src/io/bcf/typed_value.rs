//! BCF typed value encoding and decoding.
//!
//! BCF uses a sophisticated type system where each value has a type byte
//! followed by the actual data. This enables compact binary encoding while
//! maintaining type safety.

use crate::error::{BiometalError, Result};
use std::io::Read;

/// BCF value types.
///
/// These correspond to the lower 4 bits of the type byte.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum ValueType {
    /// 8-bit integer
    Int8 = 1,
    /// 16-bit integer
    Int16 = 2,
    /// 32-bit integer
    Int32 = 3,
    /// 32-bit float
    Float = 5,
    /// Character/string
    Char = 7,
}

impl ValueType {
    /// Parse value type from type byte lower 4 bits.
    pub fn from_u8(value: u8) -> Result<Self> {
        match value & 0x0F {
            1 => Ok(ValueType::Int8),
            2 => Ok(ValueType::Int16),
            3 => Ok(ValueType::Int32),
            5 => Ok(ValueType::Float),
            7 => Ok(ValueType::Char),
            other => Err(BiometalError::InvalidInput {
                msg: format!("Invalid BCF value type: {}", other),
            }),
        }
    }

    /// Get size in bytes for this type.
    pub fn size(&self) -> usize {
        match self {
            ValueType::Int8 => 1,
            ValueType::Int16 => 2,
            ValueType::Int32 => 4,
            ValueType::Float => 4,
            ValueType::Char => 1,
        }
    }
}

/// A typed value in BCF format.
///
/// Values can be scalars or arrays of primitive types.
#[derive(Debug, Clone, PartialEq)]
pub enum TypedValue {
    /// Missing value
    Missing,
    /// Single 8-bit integer
    Int8(i8),
    /// Array of 8-bit integers
    Int8Array(Vec<i8>),
    /// Single 16-bit integer
    Int16(i16),
    /// Array of 16-bit integers
    Int16Array(Vec<i16>),
    /// Single 32-bit integer
    Int32(i32),
    /// Array of 32-bit integers
    Int32Array(Vec<i32>),
    /// Single float
    Float(f32),
    /// Array of floats
    FloatArray(Vec<f32>),
    /// String
    String(String),
}

impl TypedValue {
    /// Read a typed value from a reader.
    ///
    /// The type byte encodes both the value type and array length:
    /// - Lower 4 bits: type
    /// - Upper 4 bits: length (15 = length follows as typed int)
    pub fn read<R: Read>(reader: &mut R) -> Result<Self> {
        // Read type byte
        let mut type_byte = [0u8; 1];
        reader.read_exact(&mut type_byte)?;
        let type_byte = type_byte[0];

        // Extract type and length
        let value_type = ValueType::from_u8(type_byte)?;
        let length_bits = (type_byte >> 4) & 0x0F;

        // Decode length
        let length = if length_bits == 15 {
            // Length follows as a typed integer
            let len_value = Self::read(reader)?;
            match len_value {
                TypedValue::Int8(n) => n as usize,
                TypedValue::Int16(n) => n as usize,
                TypedValue::Int32(n) => n as usize,
                _ => {
                    return Err(BiometalError::InvalidInput {
                        msg: "Invalid length encoding in BCF typed value".to_string(),
                    })
                }
            }
        } else {
            length_bits as usize
        };

        // Read value(s)
        Self::read_values(reader, value_type, length)
    }

    /// Read values of a specific type and length.
    fn read_values<R: Read>(reader: &mut R, value_type: ValueType, length: usize) -> Result<Self> {
        if length == 0 {
            return Ok(TypedValue::Missing);
        }

        match value_type {
            ValueType::Int8 => {
                if length == 1 {
                    let mut buf = [0u8; 1];
                    reader.read_exact(&mut buf)?;
                    let value = i8::from_le_bytes(buf);
                    if value == i8::MIN {
                        Ok(TypedValue::Missing)
                    } else {
                        Ok(TypedValue::Int8(value))
                    }
                } else {
                    let mut values = Vec::with_capacity(length);
                    for _ in 0..length {
                        let mut buf = [0u8; 1];
                        reader.read_exact(&mut buf)?;
                        values.push(i8::from_le_bytes(buf));
                    }
                    Ok(TypedValue::Int8Array(values))
                }
            }
            ValueType::Int16 => {
                if length == 1 {
                    let mut buf = [0u8; 2];
                    reader.read_exact(&mut buf)?;
                    let value = i16::from_le_bytes(buf);
                    if value == i16::MIN {
                        Ok(TypedValue::Missing)
                    } else {
                        Ok(TypedValue::Int16(value))
                    }
                } else {
                    let mut values = Vec::with_capacity(length);
                    for _ in 0..length {
                        let mut buf = [0u8; 2];
                        reader.read_exact(&mut buf)?;
                        values.push(i16::from_le_bytes(buf));
                    }
                    Ok(TypedValue::Int16Array(values))
                }
            }
            ValueType::Int32 => {
                if length == 1 {
                    let mut buf = [0u8; 4];
                    reader.read_exact(&mut buf)?;
                    let value = i32::from_le_bytes(buf);
                    if value == i32::MIN {
                        Ok(TypedValue::Missing)
                    } else {
                        Ok(TypedValue::Int32(value))
                    }
                } else {
                    let mut values = Vec::with_capacity(length);
                    for _ in 0..length {
                        let mut buf = [0u8; 4];
                        reader.read_exact(&mut buf)?;
                        values.push(i32::from_le_bytes(buf));
                    }
                    Ok(TypedValue::Int32Array(values))
                }
            }
            ValueType::Float => {
                if length == 1 {
                    let mut buf = [0u8; 4];
                    reader.read_exact(&mut buf)?;
                    let value = f32::from_le_bytes(buf);
                    // Check for missing value sentinel (0x7F800001)
                    if value.to_bits() == 0x7F800001 {
                        Ok(TypedValue::Missing)
                    } else {
                        Ok(TypedValue::Float(value))
                    }
                } else {
                    let mut values = Vec::with_capacity(length);
                    for _ in 0..length {
                        let mut buf = [0u8; 4];
                        reader.read_exact(&mut buf)?;
                        values.push(f32::from_le_bytes(buf));
                    }
                    Ok(TypedValue::FloatArray(values))
                }
            }
            ValueType::Char => {
                let mut buf = vec![0u8; length];
                reader.read_exact(&mut buf)?;

                // Remove null terminator if present
                if buf.last() == Some(&0) {
                    buf.pop();
                }

                String::from_utf8(buf)
                    .map(TypedValue::String)
                    .map_err(|e| BiometalError::InvalidInput {
                        msg: format!("Invalid UTF-8 in BCF string: {}", e),
                    })
            }
        }
    }

    /// Convert to integer if possible.
    pub fn as_int(&self) -> Option<i32> {
        match self {
            TypedValue::Int8(v) => Some(*v as i32),
            TypedValue::Int16(v) => Some(*v as i32),
            TypedValue::Int32(v) => Some(*v),
            _ => None,
        }
    }

    /// Convert to float if possible.
    pub fn as_float(&self) -> Option<f32> {
        match self {
            TypedValue::Float(v) => Some(*v),
            TypedValue::Int8(v) => Some(*v as f32),
            TypedValue::Int16(v) => Some(*v as f32),
            TypedValue::Int32(v) => Some(*v as f32),
            _ => None,
        }
    }

    /// Convert to string if possible.
    pub fn as_string(&self) -> Option<&str> {
        match self {
            TypedValue::String(s) => Some(s),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_value_type_from_u8() {
        assert_eq!(ValueType::from_u8(0x01).unwrap(), ValueType::Int8);
        assert_eq!(ValueType::from_u8(0x12).unwrap(), ValueType::Int16); // Upper bits ignored
        assert_eq!(ValueType::from_u8(0x03).unwrap(), ValueType::Int32);
        assert_eq!(ValueType::from_u8(0x05).unwrap(), ValueType::Float);
        assert_eq!(ValueType::from_u8(0x07).unwrap(), ValueType::Char);
        assert!(ValueType::from_u8(0x04).is_err());
    }

    #[test]
    fn test_value_type_size() {
        assert_eq!(ValueType::Int8.size(), 1);
        assert_eq!(ValueType::Int16.size(), 2);
        assert_eq!(ValueType::Int32.size(), 4);
        assert_eq!(ValueType::Float.size(), 4);
        assert_eq!(ValueType::Char.size(), 1);
    }

    #[test]
    fn test_read_int8_scalar() {
        // Type byte: 0x11 = Int8 (0x01) + length 1 (0x10)
        let data = vec![0x11, 42];
        let mut cursor = Cursor::new(data);
        let value = TypedValue::read(&mut cursor).unwrap();
        assert_eq!(value, TypedValue::Int8(42));
    }

    #[test]
    fn test_read_int16_scalar() {
        // Type byte: 0x12 = Int16 (0x02) + length 1 (0x10)
        let data = vec![0x12, 0x00, 0x01]; // 256 in little-endian
        let mut cursor = Cursor::new(data);
        let value = TypedValue::read(&mut cursor).unwrap();
        assert_eq!(value, TypedValue::Int16(256));
    }

    #[test]
    fn test_read_string() {
        // Type byte: 0x57 = Char (0x07) + length 5 (0x50)
        let data = vec![0x57, b'h', b'e', b'l', b'l', b'o'];
        let mut cursor = Cursor::new(data);
        let value = TypedValue::read(&mut cursor).unwrap();
        assert_eq!(value, TypedValue::String("hello".to_string()));
    }

    #[test]
    fn test_read_int_array() {
        // Type byte: 0x31 = Int8 (0x01) + length 3 (0x30)
        let data = vec![0x31, 1, 2, 3];
        let mut cursor = Cursor::new(data);
        let value = TypedValue::read(&mut cursor).unwrap();
        assert_eq!(value, TypedValue::Int8Array(vec![1, 2, 3]));
    }

    #[test]
    fn test_missing_value_int8() {
        // Type byte: 0x11 = Int8 + length 1
        // Value: 0x80 (i8::MIN) = missing
        let data = vec![0x11, 0x80];
        let mut cursor = Cursor::new(data);
        let value = TypedValue::read(&mut cursor).unwrap();
        assert_eq!(value, TypedValue::Missing);
    }
}
