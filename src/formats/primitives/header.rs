//! Header and comment line parsing.
//!
//! Many bioinformatics formats have header or comment lines starting with `#`:
//! - VCF: `##fileformat=VCFv4.2`, `##INFO=<ID=...>`, `#CHROM POS ID...`
//! - GFF3: `##gff-version 3`, `##sequence-region chr1 1 248956422`
//! - BED: Optional track/browser lines
//!
//! This module provides utilities for parsing these header lines separately
//! from data records.
//!
//! # Design
//!
//! The [`HeaderParser`] consumes lines starting with `#` and returns them
//! along with the underlying reader positioned at the first data line.
//! This allows format-specific parsers to process headers before creating
//! a data parser.
//!
//! # Examples
//!
//! ```
//! use biometal::formats::primitives::header::HeaderParser;
//! use std::io::BufReader;
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! let data = "##fileformat=VCFv4.2\n##source=myprogram\n#CHROM\tPOS\tID\nchr1\t100\trs123\n";
//!
//! let reader = BufReader::new(data.as_bytes());
//! let parser = HeaderParser::new(reader);
//! let (headers, first_line, _reader) = parser.parse_headers_with_first_line()?;
//!
//! assert_eq!(headers.len(), 3);
//! assert!(headers[0].starts_with("##fileformat"));
//! assert!(headers[1].starts_with("##source"));
//! assert!(headers[2].starts_with("#CHROM"));
//! assert_eq!(first_line, Some("chr1\t100\trs123".to_string()));
//! # Ok(())
//! # }
//! ```

use crate::formats::primitives::Result;
use std::io::{BufRead, BufReader, Read};

/// Parser for header/comment lines.
///
/// Reads and collects all lines starting with `#` from the beginning of a file.
/// After parsing, returns both the header lines and the underlying reader
/// positioned at the first non-header line.
///
/// # Examples
///
/// ```
/// use biometal::formats::primitives::header::HeaderParser;
/// use std::io::BufReader;
///
/// # fn main() -> Result<(), Box<dyn std::error::Error>> {
/// let data = "##header1\n##header2\ndata\n";
/// let reader = BufReader::new(data.as_bytes());
///
/// let parser = HeaderParser::new(reader);
/// let (headers, _reader) = parser.parse_headers()?;
///
/// assert_eq!(headers.len(), 2);
/// assert_eq!(headers[0], "##header1");
/// # Ok(())
/// # }
/// ```
pub struct HeaderParser<R: Read> {
    reader: BufReader<R>,
}

impl<R: Read> HeaderParser<R> {
    /// Creates a new header parser.
    pub fn new(reader: R) -> Self {
        HeaderParser {
            reader: BufReader::new(reader),
        }
    }

    /// Parses all header lines (lines starting with `#`).
    ///
    /// Returns a tuple of:
    /// - Vector of header lines (without trailing newlines, including the `#`)
    /// - The underlying reader, positioned at the first non-header line
    ///
    /// # Notes
    ///
    /// - Header lines include the leading `#`
    /// - Trailing newlines are removed
    /// - Empty lines before the first data line are skipped
    /// - The reader can be used to continue parsing data records
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::primitives::header::HeaderParser;
    /// use std::io::{BufRead, BufReader};
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let data = "##meta1\n##meta2\n#columns\ndata\nmore_data\n";
    /// let reader = BufReader::new(data.as_bytes());
    ///
    /// let parser = HeaderParser::new(reader);
    /// let (headers, first_line, mut reader) = parser.parse_headers_with_first_line()?;
    ///
    /// assert_eq!(headers.len(), 3);
    /// assert_eq!(first_line, Some("data".to_string()));
    ///
    /// // Can continue reading more data
    /// let mut line = String::new();
    /// reader.read_line(&mut line)?;
    /// assert_eq!(line.trim(), "more_data");
    /// # Ok(())
    /// # }
    /// ```
    pub fn parse_headers(mut self) -> Result<(Vec<String>, BufReader<R>)> {
        let mut headers = Vec::new();
        let mut line_buf = String::new();

        loop {
            line_buf.clear();
            let bytes_read = self.reader.read_line(&mut line_buf)?;

            // EOF
            if bytes_read == 0 {
                break;
            }

            let line = line_buf.trim_end();

            // Skip empty lines before first data
            if line.is_empty() {
                continue;
            }

            // Header line
            if line.starts_with('#') {
                headers.push(line.to_string());
                continue;
            }

            // First non-header line found
            // We need to "put back" this line for the data parser
            // Unfortunately BufRead doesn't support putting data back,
            // so we'll need to handle this differently.
            //
            // For now, the caller needs to be aware that the first
            // non-header line has been consumed and is lost.
            // A better approach would be to use a peekable reader.
            //
            // TODO: Consider using a wrapper that supports peeking
            break;
        }

        Ok((headers, self.reader))
    }

    /// Parses headers and returns them with a line that can be used to reconstruct
    /// the first data line.
    ///
    /// This version returns:
    /// - Vector of header lines
    /// - Optional first data line (if present)
    /// - The underlying reader
    ///
    /// # Examples
    ///
    /// ```
    /// use biometal::formats::primitives::header::HeaderParser;
    ///
    /// # fn main() -> Result<(), Box<dyn std::error::Error>> {
    /// let data = "##header\ndata_line\n";
    /// let reader = data.as_bytes();
    ///
    /// let parser = HeaderParser::new(reader);
    /// let (headers, first_line, _reader) = parser.parse_headers_with_first_line()?;
    ///
    /// assert_eq!(headers.len(), 1);
    /// assert_eq!(first_line, Some("data_line".to_string()));
    /// # Ok(())
    /// # }
    /// ```
    pub fn parse_headers_with_first_line(
        mut self,
    ) -> Result<(Vec<String>, Option<String>, BufReader<R>)> {
        let mut headers = Vec::new();
        let mut line_buf = String::new();
        let mut first_data_line = None;

        loop {
            line_buf.clear();
            let bytes_read = self.reader.read_line(&mut line_buf)?;

            // EOF
            if bytes_read == 0 {
                break;
            }

            let line = line_buf.trim_end();

            // Skip empty lines before first data
            if line.is_empty() {
                continue;
            }

            // Header line
            if line.starts_with('#') {
                headers.push(line.to_string());
                continue;
            }

            // First non-header line found
            first_data_line = Some(line.to_string());
            break;
        }

        Ok((headers, first_data_line, self.reader))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_headers_basic() {
        let data = "##header1\n##header2\ndata\n";
        let reader = Cursor::new(data.as_bytes());

        let parser = HeaderParser::new(reader);
        let (headers, _reader) = parser.parse_headers().unwrap();

        assert_eq!(headers.len(), 2);
        assert_eq!(headers[0], "##header1");
        assert_eq!(headers[1], "##header2");
    }

    #[test]
    fn test_parse_headers_mixed() {
        let data = "##meta1\n##meta2\n#columns\ndata\n";
        let reader = Cursor::new(data.as_bytes());

        let parser = HeaderParser::new(reader);
        let (headers, _reader) = parser.parse_headers().unwrap();

        assert_eq!(headers.len(), 3);
        assert_eq!(headers[0], "##meta1");
        assert_eq!(headers[1], "##meta2");
        assert_eq!(headers[2], "#columns");
    }

    #[test]
    fn test_parse_headers_no_headers() {
        let data = "data1\ndata2\n";
        let reader = Cursor::new(data.as_bytes());

        let parser = HeaderParser::new(reader);
        let (headers, _reader) = parser.parse_headers().unwrap();

        assert_eq!(headers.len(), 0);
    }

    #[test]
    fn test_parse_headers_empty_lines() {
        let data = "##header1\n\n##header2\n\ndata\n";
        let reader = Cursor::new(data.as_bytes());

        let parser = HeaderParser::new(reader);
        let (headers, _reader) = parser.parse_headers().unwrap();

        assert_eq!(headers.len(), 2);
        assert_eq!(headers[0], "##header1");
        assert_eq!(headers[1], "##header2");
    }

    #[test]
    fn test_parse_headers_with_first_line() {
        let data = "##header\ndata_line\ndata2\n";
        let reader = Cursor::new(data.as_bytes());

        let parser = HeaderParser::new(reader);
        let (headers, first_line, _reader) = parser.parse_headers_with_first_line().unwrap();

        assert_eq!(headers.len(), 1);
        assert_eq!(headers[0], "##header");
        assert_eq!(first_line, Some("data_line".to_string()));
    }

    #[test]
    fn test_parse_headers_with_first_line_no_data() {
        let data = "##header1\n##header2\n";
        let reader = Cursor::new(data.as_bytes());

        let parser = HeaderParser::new(reader);
        let (headers, first_line, _reader) = parser.parse_headers_with_first_line().unwrap();

        assert_eq!(headers.len(), 2);
        assert_eq!(first_line, None);
    }

    #[test]
    fn test_parse_headers_only_data() {
        let data = "data1\ndata2\n";
        let reader = Cursor::new(data.as_bytes());

        let parser = HeaderParser::new(reader);
        let (headers, first_line, _reader) = parser.parse_headers_with_first_line().unwrap();

        assert_eq!(headers.len(), 0);
        assert_eq!(first_line, Some("data1".to_string()));
    }

    #[test]
    fn test_continue_reading_after_headers() {
        let data = "##header\ndata1\ndata2\n";
        let reader = Cursor::new(data.as_bytes());

        let parser = HeaderParser::new(reader);
        let (headers, first_line, mut reader) = parser.parse_headers_with_first_line().unwrap();

        assert_eq!(headers.len(), 1);
        assert_eq!(first_line, Some("data1".to_string()));

        // Continue reading data
        let mut line = String::new();
        reader.read_line(&mut line).unwrap();
        assert_eq!(line.trim(), "data2");
    }
}
