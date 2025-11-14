//! Machine Learning operations using Apple Neural Engine
//!
//! This module provides ML-powered bioinformatics operations that leverage
//! Apple's Neural Engine via ONNX Runtime with CoreML backend.
//!
//! # Architecture
//!
//! ```text
//! PyTorch/TensorFlow → ONNX Export → ort crate → CoreML → Neural Engine
//! ```
//!
//! # Platform Support
//!
//! - **macOS ARM64 (Apple Silicon)**: Neural Engine acceleration
//! - **macOS x86_64**: CPU fallback via ONNX Runtime
//! - **Other platforms**: Feature disabled (compile-time)
//!
//! # Feature Flag
//!
//! This module requires the `neural-engine` feature flag:
//! ```toml
//! [dependencies]
//! biometal = { version = "1.8", features = ["neural-engine"] }
//! ```

#[cfg(feature = "neural-engine")]
pub mod neural_engine;

#[cfg(feature = "neural-engine")]
pub mod quality;

#[cfg(feature = "neural-engine")]
pub use neural_engine::{is_neural_engine_available, NeuralEngineContext};
