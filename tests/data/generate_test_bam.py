#!/usr/bin/env python3
"""
Generate small test BAM files for integration testing.

Creates BAM files with various edge cases:
- Different CIGAR operations
- Various tags (NM, AS, RG, MD)
- Paired and unpaired reads
- Mapped and unmapped reads
- Multiple references
- Edge cases (long CIGARs, many tags, etc.)
"""

import struct
import gzip
import sys
from typing import List, Tuple

def write_bgzf_block(data: bytes) -> bytes:
    """Write a BGZF block (gzip with extra fields)"""
    # BGZF uses gzip with specific extra fields
    # For simplicity, we'll use standard gzip which is compatible
    compressed = gzip.compress(data, compresslevel=6)
    return compressed

def write_bam_header(references: List[Tuple[str, int]]) -> bytes:
    """Write BAM header"""
    # Magic number
    header = b'BAM\x01'

    # SAM header text
    sam_text = "@HD\tVN:1.6\tSO:coordinate\n"
    for ref_name, ref_len in references:
        sam_text += f"@SQ\tSN:{ref_name}\tLN:{ref_len}\n"
    sam_text += "@RG\tID:RG001\tSM:SAMPLE1\n"
    sam_text_bytes = sam_text.encode('ascii')

    # Header length
    header += struct.pack('<I', len(sam_text_bytes))
    header += sam_text_bytes

    # Number of references
    header += struct.pack('<I', len(references))

    # Reference sequences
    for ref_name, ref_len in references:
        name_bytes = ref_name.encode('ascii') + b'\x00'
        header += struct.pack('<I', len(name_bytes))
        header += name_bytes
        header += struct.pack('<I', ref_len)

    return header

def encode_sequence(seq: str) -> bytes:
    """Encode sequence in 4-bit format"""
    seq_encoded = bytearray()
    base_map = {'=': 0, 'A': 1, 'C': 2, 'M': 3, 'G': 4, 'R': 5, 'S': 6, 'V': 7,
                'T': 8, 'W': 9, 'Y': 10, 'H': 11, 'K': 12, 'D': 13, 'B': 14, 'N': 15}

    for i in range(0, len(seq), 2):
        b1 = base_map.get(seq[i], 15)
        b2 = base_map.get(seq[i+1], 15) if i+1 < len(seq) else 0
        seq_encoded.append((b1 << 4) | b2)

    return bytes(seq_encoded)

def encode_cigar(cigar_ops: List[Tuple[str, int]]) -> bytes:
    """Encode CIGAR operations"""
    op_map = {'M': 0, 'I': 1, 'D': 2, 'N': 3, 'S': 4, 'H': 5, 'P': 6, '=': 7, 'X': 8}
    cigar_bytes = bytearray()

    for op, length in cigar_ops:
        value = (length << 4) | op_map[op]
        cigar_bytes.extend(struct.pack('<I', value))

    return bytes(cigar_bytes)

def encode_tags(tags: List[Tuple[str, str, any]]) -> bytes:
    """Encode optional tags"""
    tag_bytes = bytearray()

    for tag_name, tag_type, tag_value in tags:
        tag_bytes.extend(tag_name.encode('ascii'))
        tag_bytes.extend(tag_type.encode('ascii'))

        if tag_type == 'i':
            tag_bytes.extend(struct.pack('<i', tag_value))
        elif tag_type == 'Z':
            tag_bytes.extend(tag_value.encode('ascii') + b'\x00')
        elif tag_type == 'A':
            tag_bytes.append(ord(tag_value))

    return bytes(tag_bytes)

def write_bam_record(
    name: str,
    flag: int,
    ref_id: int,
    pos: int,
    mapq: int,
    cigar_ops: List[Tuple[str, int]],
    mate_ref_id: int,
    mate_pos: int,
    tlen: int,
    seq: str,
    qual: bytes,
    tags: List[Tuple[str, str, any]]
) -> bytes:
    """Write a single BAM record"""

    # Encode components
    name_bytes = name.encode('ascii') + b'\x00'
    cigar_bytes = encode_cigar(cigar_ops)
    seq_bytes = encode_sequence(seq)
    tag_bytes = encode_tags(tags)

    # Calculate block size (excluding the size field itself)
    block_size = (
        4 +  # ref_id
        4 +  # pos
        1 +  # l_read_name
        1 +  # mapq
        2 +  # bin
        2 +  # n_cigar_op
        2 +  # flag
        4 +  # l_seq
        4 +  # next_ref_id
        4 +  # next_pos
        4 +  # tlen
        len(name_bytes) +
        len(cigar_bytes) +
        len(seq_bytes) +
        len(qual) +
        len(tag_bytes)
    )

    record = bytearray()

    # Block size
    record.extend(struct.pack('<I', block_size))

    # Core fields
    record.extend(struct.pack('<i', ref_id))
    record.extend(struct.pack('<i', pos))
    record.append(len(name_bytes))
    record.append(mapq)
    record.extend(struct.pack('<H', 0))  # bin (can be 0)
    record.extend(struct.pack('<H', len(cigar_ops)))
    record.extend(struct.pack('<H', flag))
    record.extend(struct.pack('<I', len(seq)))
    record.extend(struct.pack('<i', mate_ref_id))
    record.extend(struct.pack('<i', mate_pos))
    record.extend(struct.pack('<i', tlen))

    # Variable-length fields
    record.extend(name_bytes)
    record.extend(cigar_bytes)
    record.extend(seq_bytes)
    record.extend(qual)
    record.extend(tag_bytes)

    return bytes(record)

def generate_test_bam():
    """Generate comprehensive test BAM file"""

    # Define references
    references = [
        ("chr1", 10000),
        ("chr2", 5000),
        ("chrM", 1000),
    ]

    # Start with header
    bam_data = write_bam_header(references)

    # Generate various test records
    records = []

    # 1. Simple mapped read with basic CIGAR
    records.append(write_bam_record(
        name="READ001",
        flag=99,  # paired, mapped, mate mapped, first in pair
        ref_id=0,
        pos=1000,
        mapq=60,
        cigar_ops=[('M', 100)],
        mate_ref_id=0,
        mate_pos=1200,
        tlen=300,
        seq="A" * 100,
        qual=b'~' * 100,  # Q~63
        tags=[('NM', 'i', 0), ('AS', 'i', 100), ('RG', 'Z', 'RG001')]
    ))

    # 2. Read with insertion
    records.append(write_bam_record(
        name="READ002",
        flag=147,  # paired, mapped, mate mapped, reverse, second in pair
        ref_id=0,
        pos=1200,
        mapq=60,
        cigar_ops=[('M', 50), ('I', 5), ('M', 45)],
        mate_ref_id=0,
        mate_pos=1000,
        tlen=-300,
        seq="C" * 100,
        qual=b'~' * 100,
        tags=[('NM', 'i', 5), ('AS', 'i', 90), ('RG', 'Z', 'RG001')]
    ))

    # 3. Read with deletion
    records.append(write_bam_record(
        name="READ003",
        flag=0,  # unmapped, single-end
        ref_id=0,
        pos=2000,
        mapq=40,
        cigar_ops=[('M', 40), ('D', 10), ('M', 50)],
        mate_ref_id=-1,
        mate_pos=-1,
        tlen=0,
        seq="G" * 90,
        qual=b'I' * 90,  # Q40
        tags=[('NM', 'i', 10), ('AS', 'i', 80), ('MD', 'Z', '40^AAAAAAAAAA50')]
    ))

    # 4. Read with soft clipping
    records.append(write_bam_record(
        name="READ004",
        flag=16,  # reverse strand
        ref_id=0,
        pos=3000,
        mapq=30,
        cigar_ops=[('S', 10), ('M', 80), ('S', 10)],
        mate_ref_id=-1,
        mate_pos=-1,
        tlen=0,
        seq="T" * 100,
        qual=b'?' * 100,  # Q30
        tags=[('NM', 'i', 2), ('AS', 'i', 78)]
    ))

    # 5. Read with N operation (RNA-seq)
    records.append(write_bam_record(
        name="READ005",
        flag=0,
        ref_id=0,
        pos=4000,
        mapq=50,
        cigar_ops=[('M', 45), ('N', 1000), ('M', 55)],
        mate_ref_id=-1,
        mate_pos=-1,
        tlen=0,
        seq="A" * 100,
        qual=b'~' * 100,
        tags=[('NM', 'i', 1), ('AS', 'i', 95), ('NH', 'i', 1)]
    ))

    # 6. Unmapped read
    records.append(write_bam_record(
        name="READ006",
        flag=4,  # unmapped
        ref_id=-1,
        pos=-1,
        mapq=0,
        cigar_ops=[],
        mate_ref_id=-1,
        mate_pos=-1,
        tlen=0,
        seq="N" * 100,
        qual=b'!' * 100,  # Q0
        tags=[]
    ))

    # 7. Duplicate read (for QC)
    records.append(write_bam_record(
        name="READ007",
        flag=1024,  # duplicate
        ref_id=0,
        pos=5000,
        mapq=60,
        cigar_ops=[('M', 100)],
        mate_ref_id=-1,
        mate_pos=-1,
        tlen=0,
        seq="C" * 100,
        qual=b'~' * 100,
        tags=[('NM', 'i', 0), ('AS', 'i', 100)]
    ))

    # 8. Secondary alignment
    records.append(write_bam_record(
        name="READ008",
        flag=256,  # secondary alignment
        ref_id=1,
        pos=1000,
        mapq=20,
        cigar_ops=[('M', 100)],
        mate_ref_id=-1,
        mate_pos=-1,
        tlen=0,
        seq="G" * 100,
        qual=b'5' * 100,  # Q20
        tags=[('NM', 'i', 3), ('AS', 'i', 85)]
    ))

    # 9. Read with hard clipping
    records.append(write_bam_record(
        name="READ009",
        flag=0,
        ref_id=0,
        pos=6000,
        mapq=40,
        cigar_ops=[('H', 20), ('M', 80), ('H', 20)],
        mate_ref_id=-1,
        mate_pos=-1,
        tlen=0,
        seq="T" * 80,
        qual=b'I' * 80,
        tags=[('NM', 'i', 1), ('AS', 'i', 79)]
    ))

    # 10. Complex CIGAR (edge case)
    records.append(write_bam_record(
        name="READ010",
        flag=0,
        ref_id=0,
        pos=7000,
        mapq=35,
        cigar_ops=[('S', 5), ('M', 30), ('I', 2), ('M', 20), ('D', 3), ('M', 38), ('S', 5)],
        mate_ref_id=-1,
        mate_pos=-1,
        tlen=0,
        seq="A" * 100,
        qual=b'F' * 100,  # Q35
        tags=[('NM', 'i', 5), ('AS', 'i', 85), ('MD', 'Z', '50^AAA38')]
    ))

    # Add all records to BAM data
    for record in records:
        bam_data += record

    # Add EOF marker (empty BGZF block)
    eof_marker = bytes.fromhex('1f8b08040000000000ff06004243020000')
    bam_data += eof_marker

    return bam_data

def main():
    """Generate test BAM file"""
    print("Generating test BAM file...")

    bam_data = generate_test_bam()

    # Write uncompressed (our BAM reader handles compression)
    output_path = "tests/data/test.bam"
    with open(output_path, 'wb') as f:
        f.write(bam_data)

    print(f"✅ Generated {output_path} ({len(bam_data)} bytes)")
    print(f"   Contains 10 test records with various edge cases")
    print(f"   - Different CIGAR operations (M, I, D, S, H, N)")
    print(f"   - Various tags (NM, AS, RG, MD)")
    print(f"   - Mapped and unmapped reads")
    print(f"   - Paired and unpaired reads")
    print(f"   - Duplicate and secondary alignments")

if __name__ == '__main__':
    main()
