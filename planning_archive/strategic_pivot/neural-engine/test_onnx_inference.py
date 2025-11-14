#!/usr/bin/env python3
"""
Test ONNX Model Inference

Verifies that the exported ONNX model produces correct predictions
and is ready for Neural Engine deployment.
"""

import argparse
from pathlib import Path

import numpy as np
import onnxruntime as ort


def encode_features(sequence: str, quality: str, max_length: int = 150) -> np.ndarray:
    """
    Encode DNA sequence and quality scores as feature vector
    (Same encoding as training script)

    Args:
        sequence: DNA sequence
        quality: Phred quality scores
        max_length: Maximum sequence length

    Returns:
        Feature vector of shape (1, max_length * 2)
    """
    seq = sequence[:max_length].upper()
    qual = quality[:max_length]

    # Encode sequence
    base_map = {'A': 0.0, 'C': 1.0, 'G': 2.0, 'T': 3.0}
    seq_features = []
    for base in seq:
        encoded = base_map.get(base, 4.0)
        seq_features.append(encoded / 4.0)

    while len(seq_features) < max_length:
        seq_features.append(0.0)

    # Encode quality scores
    qual_features = []
    for q in qual:
        phred = ord(q) - 33
        qual_features.append(phred / 93.0)

    while len(qual_features) < max_length:
        qual_features.append(0.0)

    # Concatenate and reshape to (1, features)
    features = np.array([seq_features + qual_features], dtype=np.float32)

    return features


def test_onnx_model(model_path: Path):
    """
    Test ONNX model with sample inputs

    Args:
        model_path: Path to ONNX model
    """
    print("=" * 60)
    print("ONNX Model Inference Test")
    print("=" * 60)
    print(f"Model: {model_path}")

    # Create ONNX Runtime session
    print("\nCreating ONNX Runtime session...")

    # Try CoreML execution provider first (macOS only)
    try:
        session = ort.InferenceSession(
            str(model_path),
            providers=['CoreMLExecutionProvider', 'CPUExecutionProvider']
        )
        provider = session.get_providers()[0]
        print(f"Execution provider: {provider}")

        if provider == 'CoreMLExecutionProvider':
            print("✅ Neural Engine available!")
        else:
            print("⚠️  Using CPU fallback (CoreML not available)")
    except Exception as e:
        print(f"⚠️  CoreML provider failed, using CPU: {e}")
        session = ort.InferenceSession(str(model_path), providers=['CPUExecutionProvider'])

    # Print model info
    print(f"\nModel inputs:")
    for input in session.get_inputs():
        print(f"  - {input.name}: {input.shape} ({input.type})")

    print(f"\nModel outputs:")
    for output in session.get_outputs():
        print(f"  - {output.name}: {output.shape} ({output.type})")

    # Test cases
    test_cases = [
        {
            'name': 'High quality read (should PASS)',
            'sequence': 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',
            'quality': 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII',  # Q40
            'expected': True
        },
        {
            'name': 'Low quality read (should FAIL)',
            'sequence': 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',
            'quality': '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!',  # Q0
            'expected': False
        },
        {
            'name': 'Medium quality read (borderline)',
            'sequence': 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',
            'quality': '5555555555555555555555555555555555555555555555555555555555555555',  # Q20
            'expected': True
        },
    ]

    print("\n" + "=" * 60)
    print("Running Test Cases")
    print("=" * 60)

    for i, test in enumerate(test_cases, 1):
        print(f"\nTest {i}: {test['name']}")
        print(f"Sequence: {test['sequence'][:40]}...")
        print(f"Quality:  {test['quality'][:40]}...")

        # Encode features
        features = encode_features(test['sequence'], test['quality'])

        # Run inference
        input_name = session.get_inputs()[0].name
        output_name = session.get_outputs()[0].name

        outputs = session.run([output_name], {input_name: features})
        probability = outputs[0][0][0]

        # Prediction (threshold = 0.5)
        prediction = probability > 0.5
        status = "✅ PASS" if prediction else "❌ FAIL"

        print(f"Probability: {probability:.4f}")
        print(f"Prediction:  {status}")
        print(f"Expected:    {'✅ PASS' if test['expected'] else '❌ FAIL'}")

        if prediction == test['expected']:
            print("Result: ✅ CORRECT")
        else:
            print("Result: ⚠️  MISMATCH (may need more training)")

    print("\n" + "=" * 60)
    print("Test Complete")
    print("=" * 60)
    print("\nNext steps:")
    print("1. If tests pass, the model is ready for Rust integration")
    print("2. Test with Rust Neural Engine code:")
    print("   cargo test --features neural-engine test_quality_predictor")
    print("3. Benchmark Neural Engine performance:")
    print("   cargo bench --features neural-engine neural_engine")


def main():
    parser = argparse.ArgumentParser(description='Test ONNX model inference')
    parser.add_argument('--model', type=Path, required=True, help='Path to ONNX model')
    args = parser.parse_args()

    if not args.model.exists():
        print(f"Error: Model file not found: {args.model}")
        return

    test_onnx_model(args.model)


if __name__ == '__main__':
    main()
