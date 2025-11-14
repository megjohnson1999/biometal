#!/usr/bin/env python3
"""
Create Minimal ONNX Model for Testing

This creates a simple dummy model for testing Neural Engine integration
while the full PyTorch training completes.

The model is a simple linear classifier that approximates quality filtering:
- Average the quality features (last 150 inputs)
- Apply threshold (if avg > 0.215, predict PASS ~= Q20)
"""

import numpy as np
from pathlib import Path


def create_minimal_onnx_model(output_path: Path):
    """
    Create minimal ONNX model without PyTorch dependency

    This is a simple linear model that averages quality scores
    and applies a threshold - approximating quality_filter logic.
    """
    try:
        import onnx
        from onnx import helper, TensorProto
    except ImportError:
        print("Error: onnx package required")
        print("Install with: pip install onnx")
        return False

    # Input: [batch_size, 300] features
    # - First 150: sequence features (ignored by this simple model)
    # - Last 150: quality features (averaged)
    # Output: [batch_size, 1] probability

    # Create a simple model: average quality features and threshold
    # quality features are at indices 150-299
    # avg_quality = sum(features[150:300]) / 150
    # if avg_quality > 0.215 (≈ Q20/93), predict PASS

    input_size = 300

    # Define input
    input_tensor = helper.make_tensor_value_info(
        'input', TensorProto.FLOAT, [1, input_size]
    )

    # Define output
    output_tensor = helper.make_tensor_value_info(
        'output', TensorProto.FLOAT, [1, 1]
    )

    # Create weights that extract quality features (last 150)
    # Shape: [input_size, 1]
    weights = np.zeros((input_size, 1), dtype=np.float32)
    # Only quality features contribute (indices 150-299)
    weights[150:300, 0] = 1.0 / 150.0  # Average quality scores

    # Bias: threshold at -0.215 (so prediction > 0.5 when avg_quality > 0.215)
    bias = np.array([-0.215], dtype=np.float32)

    # Create weight initializers
    weight_init = helper.make_tensor(
        'weights',
        TensorProto.FLOAT,
        [input_size, 1],
        weights.flatten().tolist()
    )

    bias_init = helper.make_tensor(
        'bias',
        TensorProto.FLOAT,
        [1],
        bias.tolist()
    )

    # Create nodes
    # 1. Linear transformation: input @ weights + bias
    matmul_node = helper.make_node(
        'MatMul',
        inputs=['input', 'weights'],
        outputs=['matmul_out']
    )

    add_node = helper.make_node(
        'Add',
        inputs=['matmul_out', 'bias'],
        outputs=['linear_out']
    )

    # 2. Sigmoid activation
    sigmoid_node = helper.make_node(
        'Sigmoid',
        inputs=['linear_out'],
        outputs=['output']
    )

    # Create graph
    graph = helper.make_graph(
        nodes=[matmul_node, add_node, sigmoid_node],
        name='QualityPredictor',
        inputs=[input_tensor],
        outputs=[output_tensor],
        initializer=[weight_init, bias_init]
    )

    # Create model with IR version 9 (compatible with ort 2.0)
    model = helper.make_model(graph, producer_name='biometal', ir_version=9)
    model.opset_import[0].version = 13  # CoreML compatible

    # Save model
    onnx.save(model, str(output_path))

    print(f"✅ Minimal ONNX model created: {output_path}")
    print(f"Model size: {output_path.stat().st_size / 1024:.2f} KB")
    print(f"\nThis is a simple linear model for testing.")
    print(f"It averages quality scores and thresholds at Q20.")
    print(f"\nReplace with trained model for better performance.")

    return True


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Create minimal ONNX model for testing')
    parser.add_argument('--output', type=Path, default='minimal_model.onnx',
                        help='Output ONNX model path')
    args = parser.parse_args()

    success = create_minimal_onnx_model(args.output)

    if success:
        print("\nTest with:")
        print(f"  python test_onnx_inference.py --model {args.output}")
        print(f"  cargo run --example neural_quality --features neural-engine")


if __name__ == '__main__':
    main()
