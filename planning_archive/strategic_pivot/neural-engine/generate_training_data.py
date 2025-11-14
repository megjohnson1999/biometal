#!/usr/bin/env python3
"""
Generate Synthetic FASTQ Data for Neural Engine Training

Creates synthetic FASTQ reads with varying quality scores for training
the read quality prediction model.
"""

import gzip
import random
from pathlib import Path


def generate_sequence(length: int = 150) -> str:
    """Generate random DNA sequence"""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))


def generate_quality(length: int = 150, avg_quality: int = 30, variation: int = 10) -> str:
    """
    Generate quality string with specified average

    Args:
        length: Length of quality string
        avg_quality: Target average Phred score
        variation: Standard deviation of quality scores

    Returns:
        Quality string (ASCII-33 encoded)
    """
    qualities = []
    for _ in range(length):
        # Generate Phred score with normal distribution
        phred = int(random.gauss(avg_quality, variation))
        # Clamp to valid range [0, 93]
        phred = max(0, min(93, phred))
        # Convert to ASCII (Phred + 33)
        qualities.append(chr(phred + 33))

    return ''.join(qualities)


def generate_fastq_record(record_id: int, avg_quality: int, length: int = 150) -> str:
    """
    Generate a single FASTQ record

    Format:
    @read_id
    SEQUENCE
    +
    QUALITY

    Args:
        record_id: Unique record identifier
        avg_quality: Target average quality score
        length: Sequence length

    Returns:
        FASTQ record as string
    """
    name = f"@read_{record_id:08d} avg_q={avg_quality}"
    sequence = generate_sequence(length)
    plus = "+"
    quality = generate_quality(length, avg_quality, variation=10)

    return f"{name}\n{sequence}\n{plus}\n{quality}\n"


def generate_training_dataset(
    output_path: Path,
    num_records: int = 50000,
    length: int = 150,
    seed: int = 42
):
    """
    Generate balanced training dataset with varying quality

    Creates reads with quality scores distributed across range:
    - Low quality (Q5-Q15): 25% of reads → FAIL
    - Medium quality (Q15-Q25): 25% of reads → Mixed
    - High quality (Q25-Q40): 50% of reads → PASS

    Args:
        output_path: Path to output FASTQ file (.fq.gz)
        num_records: Number of records to generate
        length: Sequence length
        seed: Random seed for reproducibility
    """
    random.seed(seed)

    print(f"Generating {num_records} synthetic FASTQ records...")
    print(f"Output: {output_path}")
    print(f"Sequence length: {length}")
    print(f"Random seed: {seed}")

    # Create quality distribution
    quality_ranges = [
        (5, 15, 0.25),    # Low quality: 25% (FAIL)
        (15, 25, 0.25),   # Medium quality: 25% (Mixed)
        (25, 40, 0.50),   # High quality: 50% (PASS)
    ]

    # Calculate records per range
    records_per_range = []
    for min_q, max_q, proportion in quality_ranges:
        count = int(num_records * proportion)
        records_per_range.append((min_q, max_q, count))
        print(f"  Q{min_q}-Q{max_q}: {count} records ({proportion*100:.0f}%)")

    # Generate records
    record_id = 0
    records_generated = 0

    with gzip.open(output_path, 'wt') as f:
        for min_q, max_q, count in records_per_range:
            for _ in range(count):
                # Random quality in range
                avg_quality = random.randint(min_q, max_q)

                # Generate record
                record = generate_fastq_record(record_id, avg_quality, length)
                f.write(record)

                record_id += 1
                records_generated += 1

                # Progress indicator
                if records_generated % 10000 == 0:
                    print(f"  Generated {records_generated}/{num_records} records...")

    print(f"\n✅ Dataset generation complete!")
    print(f"Total records: {records_generated}")
    print(f"File size: {output_path.stat().st_size / 1024 / 1024:.2f} MB")


def main():
    import argparse

    parser = argparse.ArgumentParser(description='Generate synthetic FASTQ training data')
    parser.add_argument('--output', type=Path, default='training_data.fq.gz',
                        help='Output FASTQ file (.fq.gz)')
    parser.add_argument('--num-records', type=int, default=50000,
                        help='Number of records (default: 50000)')
    parser.add_argument('--length', type=int, default=150,
                        help='Sequence length (default: 150)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed (default: 42)')

    args = parser.parse_args()

    generate_training_dataset(
        output_path=args.output,
        num_records=args.num_records,
        length=args.length,
        seed=args.seed
    )

    print("\nNext step:")
    print(f"  python train_quality_model.py --fastq {args.output} --output quality_model.onnx")


if __name__ == '__main__':
    main()
