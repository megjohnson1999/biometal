#!/usr/bin/env python3
"""
Train Read Quality Prediction Model for Neural Engine

This script trains a simple MLP to predict if a FASTQ read passes quality thresholds.
The trained model is exported to ONNX format for Apple Neural Engine inference.

Model Architecture:
- Input: Concatenated sequence + quality scores (300 features for 150bp read)
- Hidden: 2-3 dense layers (64-128 neurons)
- Output: Sigmoid activation → probability (PASS/FAIL)

Training Data:
- 10,000-100,000 FASTQ reads from existing test data
- Labels from biometal's quality_filter() as ground truth
- Features: Encoded sequences (A=0, C=1, G=2, T=3, N=4) + normalized Phred scores

Export:
- ONNX format (opset 13) for maximum CoreML compatibility
- Output: quality_model.onnx (ready for Neural Engine)
"""

import argparse
import gzip
import random
from pathlib import Path
from typing import List, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score


class FastqRecord:
    """Simple FASTQ record representation"""
    def __init__(self, name: str, sequence: str, quality: str):
        self.name = name
        self.sequence = sequence
        self.quality = quality

    def __repr__(self):
        return f"FastqRecord(name={self.name}, len={len(self.sequence)})"


def parse_fastq(fastq_path: Path, max_records: int = None) -> List[FastqRecord]:
    """
    Parse FASTQ file (gzipped or plain text)

    Args:
        fastq_path: Path to FASTQ file (.fq, .fastq, .fq.gz, .fastq.gz)
        max_records: Maximum number of records to read (None = all)

    Returns:
        List of FastqRecord objects
    """
    records = []

    # Determine if file is gzipped
    open_fn = gzip.open if str(fastq_path).endswith('.gz') else open

    with open_fn(fastq_path, 'rt') as f:
        while True:
            # Read 4 lines per record
            name_line = f.readline()
            if not name_line:
                break

            seq_line = f.readline()
            plus_line = f.readline()
            qual_line = f.readline()

            if not qual_line:
                break

            # Parse record
            name = name_line.strip()[1:]  # Remove '@'
            sequence = seq_line.strip()
            quality = qual_line.strip()

            records.append(FastqRecord(name, sequence, quality))

            if max_records and len(records) >= max_records:
                break

    return records


def quality_filter_python(sequence: str, quality: str, min_quality: int = 20, min_length: int = 50) -> bool:
    """
    Python implementation of biometal's quality_filter logic

    This is the ground truth labeling function.

    Args:
        sequence: DNA sequence
        quality: Phred quality scores (ASCII-33 encoded)
        min_quality: Minimum average quality score
        min_length: Minimum sequence length

    Returns:
        True if read passes quality, False otherwise
    """
    # Check length
    if len(sequence) < min_length:
        return False

    # Calculate average quality
    if not quality:
        return False

    # Convert ASCII to Phred scores (ASCII - 33)
    phred_scores = [ord(q) - 33 for q in quality]
    avg_quality = sum(phred_scores) / len(phred_scores)

    return avg_quality >= min_quality


def encode_features(sequence: str, quality: str, max_length: int = 150) -> np.ndarray:
    """
    Encode DNA sequence and quality scores as feature vector

    Encoding:
    - Sequence: A=0, C=1, G=2, T=3, N=4 (normalized to [0, 1])
    - Quality: Phred scores (ASCII-33) normalized to [0, 1]

    Args:
        sequence: DNA sequence
        quality: Phred quality scores
        max_length: Maximum sequence length (padding/truncation)

    Returns:
        Feature vector of shape (max_length * 2,)
    """
    # Truncate or pad sequence
    seq = sequence[:max_length].upper()
    qual = quality[:max_length]

    # Encode sequence
    base_map = {'A': 0.0, 'C': 1.0, 'G': 2.0, 'T': 3.0}
    seq_features = []
    for base in seq:
        encoded = base_map.get(base, 4.0)  # N or unknown → 4.0
        seq_features.append(encoded / 4.0)  # Normalize to [0, 1]

    # Pad sequence features if needed
    while len(seq_features) < max_length:
        seq_features.append(0.0)

    # Encode quality scores
    qual_features = []
    for q in qual:
        phred = ord(q) - 33  # ASCII to Phred
        qual_features.append(phred / 93.0)  # Normalize to [0, 1] (max Phred = 93)

    # Pad quality features if needed
    while len(qual_features) < max_length:
        qual_features.append(0.0)

    # Concatenate sequence + quality
    features = np.array(seq_features + qual_features, dtype=np.float32)

    return features


class ReadQualityDataset(Dataset):
    """PyTorch dataset for read quality prediction"""

    def __init__(self, records: List[FastqRecord], max_length: int = 150):
        """
        Args:
            records: List of FastqRecord objects
            max_length: Maximum sequence length for encoding
        """
        self.records = records
        self.max_length = max_length

        # Pre-compute features and labels
        self.features = []
        self.labels = []

        for record in records:
            # Encode features
            features = encode_features(record.sequence, record.quality, max_length)
            self.features.append(features)

            # Generate label using quality_filter
            label = 1.0 if quality_filter_python(record.sequence, record.quality) else 0.0
            self.labels.append(label)

        self.features = np.array(self.features)
        self.labels = np.array(self.labels)

    def __len__(self):
        return len(self.records)

    def __getitem__(self, idx):
        return torch.tensor(self.features[idx]), torch.tensor(self.labels[idx])


class QualityPredictorMLP(nn.Module):
    """
    Simple MLP for read quality prediction

    Architecture:
    - Input: 300 features (150bp * 2 for seq + quality)
    - Hidden1: 128 neurons + ReLU + Dropout(0.2)
    - Hidden2: 64 neurons + ReLU + Dropout(0.2)
    - Output: 1 neuron + Sigmoid → probability
    """

    def __init__(self, input_size: int = 300, hidden1: int = 128, hidden2: int = 64):
        super(QualityPredictorMLP, self).__init__()

        self.network = nn.Sequential(
            nn.Linear(input_size, hidden1),
            nn.ReLU(),
            nn.Dropout(0.2),

            nn.Linear(hidden1, hidden2),
            nn.ReLU(),
            nn.Dropout(0.2),

            nn.Linear(hidden2, 1),
            nn.Sigmoid()
        )

    def forward(self, x):
        return self.network(x)


def train_model(
    train_loader: DataLoader,
    val_loader: DataLoader,
    model: nn.Module,
    epochs: int = 20,
    lr: float = 0.001,
    device: str = 'cpu'
) -> Tuple[nn.Module, List[float], List[float]]:
    """
    Train the quality prediction model

    Args:
        train_loader: Training data loader
        val_loader: Validation data loader
        model: PyTorch model
        epochs: Number of training epochs
        lr: Learning rate
        device: Device to train on ('cpu' or 'cuda')

    Returns:
        Trained model, training losses, validation losses
    """
    model = model.to(device)
    criterion = nn.BCELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr)

    train_losses = []
    val_losses = []

    print(f"\nTraining on {device}...")
    print(f"Epochs: {epochs}, Learning rate: {lr}")
    print("-" * 60)

    for epoch in range(epochs):
        # Training phase
        model.train()
        train_loss = 0.0

        for features, labels in train_loader:
            features = features.to(device)
            labels = labels.to(device).unsqueeze(1)

            # Forward pass
            optimizer.zero_grad()
            outputs = model(features)
            loss = criterion(outputs, labels)

            # Backward pass
            loss.backward()
            optimizer.step()

            train_loss += loss.item()

        train_loss /= len(train_loader)
        train_losses.append(train_loss)

        # Validation phase
        model.eval()
        val_loss = 0.0

        with torch.no_grad():
            for features, labels in val_loader:
                features = features.to(device)
                labels = labels.to(device).unsqueeze(1)

                outputs = model(features)
                loss = criterion(outputs, labels)

                val_loss += loss.item()

        val_loss /= len(val_loader)
        val_losses.append(val_loss)

        # Print progress
        print(f"Epoch {epoch+1:2d}/{epochs} | Train Loss: {train_loss:.4f} | Val Loss: {val_loss:.4f}")

    return model, train_losses, val_losses


def evaluate_model(model: nn.Module, test_loader: DataLoader, device: str = 'cpu') -> dict:
    """
    Evaluate model performance on test set

    Args:
        model: Trained PyTorch model
        test_loader: Test data loader
        device: Device to evaluate on

    Returns:
        Dictionary with accuracy, precision, recall, F1 score
    """
    model.eval()
    model = model.to(device)

    all_predictions = []
    all_labels = []

    with torch.no_grad():
        for features, labels in test_loader:
            features = features.to(device)
            outputs = model(features)

            # Convert probabilities to binary predictions (threshold = 0.5)
            predictions = (outputs.cpu().numpy() > 0.5).astype(int).flatten()

            all_predictions.extend(predictions)
            all_labels.extend(labels.numpy().astype(int))

    # Calculate metrics
    accuracy = accuracy_score(all_labels, all_predictions)
    precision = precision_score(all_labels, all_predictions, zero_division=0)
    recall = recall_score(all_labels, all_predictions, zero_division=0)
    f1 = f1_score(all_labels, all_predictions, zero_division=0)

    return {
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1': f1
    }


def export_to_onnx(model: nn.Module, output_path: Path, input_size: int = 300):
    """
    Export PyTorch model to ONNX format

    Args:
        model: Trained PyTorch model
        output_path: Path to save ONNX model
        input_size: Input feature size
    """
    model.eval()

    # Create dummy input (batch_size=1, features=input_size)
    dummy_input = torch.randn(1, input_size)

    # Export to ONNX
    torch.onnx.export(
        model,
        dummy_input,
        str(output_path),
        export_params=True,
        opset_version=13,  # CoreML supports opset 13
        do_constant_folding=True,
        input_names=['input'],
        output_names=['output'],
        dynamic_axes={
            'input': {0: 'batch_size'},
            'output': {0: 'batch_size'}
        }
    )

    print(f"\nModel exported to ONNX: {output_path}")
    print(f"Input shape: [batch_size, {input_size}]")
    print(f"Output shape: [batch_size, 1]")


def main():
    parser = argparse.ArgumentParser(description='Train read quality prediction model for Neural Engine')
    parser.add_argument('--fastq', type=Path, required=True, help='Path to FASTQ file (.fq.gz or .fq)')
    parser.add_argument('--max-records', type=int, default=50000, help='Maximum records to use (default: 50000)')
    parser.add_argument('--max-length', type=int, default=150, help='Maximum sequence length (default: 150)')
    parser.add_argument('--epochs', type=int, default=20, help='Training epochs (default: 20)')
    parser.add_argument('--batch-size', type=int, default=64, help='Batch size (default: 64)')
    parser.add_argument('--lr', type=float, default=0.001, help='Learning rate (default: 0.001)')
    parser.add_argument('--output', type=Path, default='quality_model.onnx', help='Output ONNX model path')
    parser.add_argument('--seed', type=int, default=42, help='Random seed (default: 42)')

    args = parser.parse_args()

    # Set random seeds for reproducibility
    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    print("=" * 60)
    print("Read Quality Prediction Model Training")
    print("=" * 60)
    print(f"FASTQ file: {args.fastq}")
    print(f"Max records: {args.max_records}")
    print(f"Max length: {args.max_length}")
    print(f"Epochs: {args.epochs}")
    print(f"Batch size: {args.batch_size}")
    print(f"Learning rate: {args.lr}")
    print(f"Output: {args.output}")
    print("=" * 60)

    # Parse FASTQ file
    print(f"\nParsing FASTQ file...")
    records = parse_fastq(args.fastq, max_records=args.max_records)
    print(f"Loaded {len(records)} records")

    # Create dataset
    print(f"\nEncoding features...")
    dataset = ReadQualityDataset(records, max_length=args.max_length)

    # Print class distribution
    num_pass = int(dataset.labels.sum())
    num_fail = len(dataset.labels) - num_pass
    print(f"Class distribution:")
    print(f"  PASS: {num_pass} ({num_pass/len(dataset.labels)*100:.1f}%)")
    print(f"  FAIL: {num_fail} ({num_fail/len(dataset.labels)*100:.1f}%)")

    # Split into train/val/test (70/15/15)
    train_idx, temp_idx = train_test_split(
        range(len(dataset)),
        test_size=0.3,
        random_state=args.seed,
        stratify=dataset.labels
    )
    val_idx, test_idx = train_test_split(
        temp_idx,
        test_size=0.5,
        random_state=args.seed,
        stratify=dataset.labels[temp_idx]
    )

    train_subset = torch.utils.data.Subset(dataset, train_idx)
    val_subset = torch.utils.data.Subset(dataset, val_idx)
    test_subset = torch.utils.data.Subset(dataset, test_idx)

    print(f"\nDataset split:")
    print(f"  Train: {len(train_subset)} records")
    print(f"  Val:   {len(val_subset)} records")
    print(f"  Test:  {len(test_subset)} records")

    # Create data loaders
    train_loader = DataLoader(train_subset, batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(val_subset, batch_size=args.batch_size, shuffle=False)
    test_loader = DataLoader(test_subset, batch_size=args.batch_size, shuffle=False)

    # Create model
    input_size = args.max_length * 2  # sequence + quality
    model = QualityPredictorMLP(input_size=input_size)

    print(f"\nModel architecture:")
    print(model)
    print(f"\nTotal parameters: {sum(p.numel() for p in model.parameters()):,}")

    # Train model
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model, train_losses, val_losses = train_model(
        train_loader, val_loader, model,
        epochs=args.epochs, lr=args.lr, device=device
    )

    # Evaluate on test set
    print("\n" + "=" * 60)
    print("Test Set Evaluation")
    print("=" * 60)
    metrics = evaluate_model(model, test_loader, device=device)
    print(f"Accuracy:  {metrics['accuracy']:.4f}")
    print(f"Precision: {metrics['precision']:.4f}")
    print(f"Recall:    {metrics['recall']:.4f}")
    print(f"F1 Score:  {metrics['f1']:.4f}")

    # Export to ONNX
    print("\n" + "=" * 60)
    print("Exporting to ONNX")
    print("=" * 60)
    export_to_onnx(model, args.output, input_size=input_size)

    print("\n" + "=" * 60)
    print("Training Complete!")
    print("=" * 60)
    print(f"\nNext steps:")
    print(f"1. Test ONNX model with ort in Python:")
    print(f"   python test_onnx_inference.py --model {args.output}")
    print(f"2. Test with Rust Neural Engine integration:")
    print(f"   cargo test --features neural-engine")
    print(f"3. Benchmark Neural Engine performance")


if __name__ == '__main__':
    main()
