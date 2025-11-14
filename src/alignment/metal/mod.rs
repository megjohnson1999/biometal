//! Metal GPU acceleration module
//!
//! This module provides Metal GPU compute shader support for Smith-Waterman
//! alignment on Apple Silicon.

pub mod gpu;

pub use gpu::{is_metal_available, MetalContext};
