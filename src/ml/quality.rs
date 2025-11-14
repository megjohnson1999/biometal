//! ML-powered read quality prediction
//!
//! Uses Neural Engine to predict if a FASTQ read passes quality thresholds.
//!
//! # Model Architecture
//!
//! - Input: Concatenated sequence + quality scores (300 features for 150bp read)
//! - Hidden: 2-3 dense layers (64-128 neurons)
//! - Output: Sigmoid activation → probability (PASS/FAIL)
//!
//! # Example
//!
//! ```no_run
//! use biometal::ml::quality::QualityPredictor;
//!
//! let predictor = QualityPredictor::new("quality_model.onnx").unwrap();
//! let passes = predictor.predict_quality(b"ACGT...", &[30, 35, 40, ...]).unwrap();
//! ```

use super::neural_engine::NeuralEngineContext;
use std::path::Path;

/// ML-powered quality predictor
///
/// Uses Neural Engine to predict read quality in <1ms.
pub struct QualityPredictor {
    context: NeuralEngineContext,
    sequence_length: usize,
}

impl QualityPredictor {
    /// Create a new quality predictor from ONNX model
    ///
    /// # Arguments
    ///
    /// * `model_path` - Path to trained ONNX model
    ///
    /// # Returns
    ///
    /// Quality predictor ready for inference
    pub fn new<P: AsRef<Path>>(model_path: P) -> Result<Self, String> {
        let context = NeuralEngineContext::new(model_path)?;

        // Infer sequence length from model metadata
        let metadata = context.metadata();
        let sequence_length = if let Some(input) = metadata.inputs.first() {
            // Assuming input shape is [batch, features]
            // features = sequence_length * 2 (seq + quality)
            if input.shape.len() >= 2 {
                (input.shape[1] / 2) as usize
            } else {
                150  // Default for 150bp reads
            }
        } else {
            150
        };

        Ok(Self {
            context,
            sequence_length,
        })
    }

    /// Predict if a read passes quality thresholds
    ///
    /// # Arguments
    ///
    /// * `sequence` - DNA sequence (ACGT)
    /// * `quality` - Phred quality scores (0-93)
    ///
    /// # Returns
    ///
    /// - `true` if read passes quality
    /// - `false` if read fails quality
    ///
    /// # Example
    ///
    /// ```no_run
    /// use biometal::ml::quality::QualityPredictor;
    ///
    /// let predictor = QualityPredictor::new("model.onnx").unwrap();
    /// let sequence = b"ACGTACGTACGT...";
    /// let quality = vec![30, 35, 40, 38, ...];
    ///
    /// if predictor.predict_quality(sequence, &quality).unwrap() {
    ///     println!("Read passes quality");
    /// }
    /// ```
    pub fn predict_quality(&mut self, sequence: &[u8], quality: &[u8]) -> Result<bool, String> {
        // Encode sequence and quality as features
        let features = self.encode_features(sequence, quality)?;

        // Run inference
        let output = self.context.infer(&features, &[1, features.len() as i64])?;

        // Extract probability (sigmoid output)
        let probability = output.get(0).ok_or("No output returned")?;

        // Threshold at 0.5
        Ok(*probability > 0.5)
    }

    /// Predict quality with probability score
    ///
    /// Returns (passes, probability) tuple
    pub fn predict_quality_with_score(
        &mut self,
        sequence: &[u8],
        quality: &[u8],
    ) -> Result<(bool, f32), String> {
        let features = self.encode_features(sequence, quality)?;
        let output = self.context.infer(&features, &[1, features.len() as i64])?;
        let probability = *output.get(0).ok_or("No output returned")?;
        Ok((probability > 0.5, probability))
    }

    /// Encode sequence and quality as feature vector
    fn encode_features(&self, sequence: &[u8], quality: &[u8]) -> Result<Vec<f32>, String> {
        if sequence.len() != quality.len() {
            return Err(format!(
                "Sequence and quality length mismatch: {} vs {}",
                sequence.len(),
                quality.len()
            ));
        }

        if sequence.len() > self.sequence_length {
            return Err(format!(
                "Sequence too long: {} > {}",
                sequence.len(),
                self.sequence_length
            ));
        }

        // Encode sequence: A=0, C=1, G=2, T=3, N=4
        let mut features = Vec::with_capacity(self.sequence_length * 2);

        for &base in sequence {
            let encoded = match base {
                b'A' | b'a' => 0.0,
                b'C' | b'c' => 1.0,
                b'G' | b'g' => 2.0,
                b'T' | b't' => 3.0,
                _ => 4.0,  // N or other
            };
            features.push(encoded / 4.0);  // Normalize to [0, 1]
        }

        // Pad sequence if needed
        while features.len() < self.sequence_length {
            features.push(0.0);
        }

        // Encode quality scores (normalize Phred scores)
        for &qual in quality {
            features.push((qual as f32) / 93.0);  // Normalize to [0, 1]
        }

        // Pad quality if needed
        while features.len() < self.sequence_length * 2 {
            features.push(0.0);
        }

        Ok(features)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_features() {
        // This will be tested once we have a model
        // For now, just verify the module compiles
    }
}
