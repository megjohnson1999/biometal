//! Property-based tests for format parsers.
//!
//! Tests invariants and round-trip parsing for BED, GFA, VCF, and GFF3 formats.
//! Uses proptest for randomized testing with thousands of generated inputs.

use biometal::formats::bed::{Bed3Record, Bed6Record};
use biometal::formats::gfa::{GfaSegment, GfaLink, Orientation};
use biometal::formats::gff::Gff3Record;
use biometal::formats::vcf::VcfRecord;
use biometal::formats::primitives::{GenomicInterval, Strand};
use biometal::formats::TabDelimitedRecord;
use proptest::prelude::*;
use std::collections::HashMap;

// ============================================================================
// BED Format Property Tests
// ============================================================================

mod bed_properties {
    use super::*;

    /// Generate arbitrary chromosome names
    fn arb_chrom() -> impl Strategy<Value = String> {
        prop_oneof![
            Just("chr1".to_string()),
            Just("chr2".to_string()),
            Just("chrX".to_string()),
            Just("chrM".to_string()),
        ]
    }

    /// Generate arbitrary strand (BED always has a strand field)
    fn arb_strand() -> impl Strategy<Value = Option<Strand>> {
        prop_oneof![
            Just(Some(Strand::Forward)),
            Just(Some(Strand::Reverse)),
            Just(Some(Strand::Unknown)),
        ]
    }

    /// Strategy for generating valid BED3 records
    fn arb_bed3_record() -> impl Strategy<Value = Bed3Record> {
        (arb_chrom(), 0u64..1_000_000u64, 1u64..1_000_000u64)
            .prop_filter("end must be > start", |(_, start, end)| end > start)
            .prop_map(|(chrom, start, end)| Bed3Record {
                interval: GenomicInterval { chrom, start, end },
            })
    }

    /// Strategy for generating valid BED6 records
    fn arb_bed6_record() -> impl Strategy<Value = Bed6Record> {
        (
            arb_bed3_record(),
            prop::option::of("[a-zA-Z0-9_]+"),
            prop::option::of(0u32..1000u32),
            arb_strand(),
        )
            .prop_map(|(bed3, name, score, strand)| Bed6Record {
                bed3,
                name,
                score,
                strand,
            })
    }

    proptest! {
        #[test]
        fn bed3_round_trip(record in arb_bed3_record()) {
            let line = record.to_line();
            let parsed = Bed3Record::from_line(&line).unwrap();
            prop_assert_eq!(record, parsed);
        }

        #[test]
        fn bed3_valid_coordinates(record in arb_bed3_record()) {
            prop_assert!(record.interval.end > record.interval.start);
            prop_assert!(record.interval.length() > 0);
        }

        #[test]
        fn bed3_length_invariant(record in arb_bed3_record()) {
            let length = record.interval.length();
            prop_assert_eq!(length, record.interval.end - record.interval.start);
        }

        #[test]
        fn bed6_round_trip(record in arb_bed6_record()) {
            let line = record.to_line();
            let parsed = Bed6Record::from_line(&line).unwrap();
            prop_assert_eq!(record, parsed);
        }

        #[test]
        fn bed6_inherits_bed3_invariants(record in arb_bed6_record()) {
            prop_assert!(record.bed3.interval.end > record.bed3.interval.start);
            prop_assert!(record.bed3.interval.length() > 0);
        }

        #[test]
        fn bed6_score_range(record in arb_bed6_record()) {
            if let Some(score) = record.score {
                prop_assert!(score <= 1000, "BED score should be 0-1000");
            }
        }
    }
}

// ============================================================================
// GFA Format Property Tests
// ============================================================================

mod gfa_properties {
    use super::*;

    /// Generate arbitrary segment names
    fn arb_seg_name() -> impl Strategy<Value = String> {
        "[a-zA-Z][a-zA-Z0-9_]{0,10}"
    }

    /// Generate arbitrary DNA sequence
    fn arb_dna_sequence() -> impl Strategy<Value = String> {
        prop::collection::vec(prop::sample::select(vec!['A', 'C', 'G', 'T']), 10..100)
            .prop_map(|chars| chars.into_iter().collect())
    }

    /// Generate arbitrary orientation
    fn arb_orientation() -> impl Strategy<Value = Orientation> {
        prop_oneof![
            Just(Orientation::Forward),
            Just(Orientation::Reverse),
        ]
    }

    /// Strategy for generating valid GFA segments
    fn arb_gfa_segment() -> impl Strategy<Value = GfaSegment> {
        (arb_seg_name(), arb_dna_sequence()).prop_map(|(name, sequence)| {
            let mut tags = HashMap::new();
            tags.insert("LN".to_string(), format!("i:{}", sequence.len()));
            GfaSegment {
                name,
                sequence,
                tags,
            }
        })
    }

    /// Strategy for generating valid GFA links
    fn arb_gfa_link() -> impl Strategy<Value = GfaLink> {
        (
            arb_seg_name(),
            arb_orientation(),
            arb_seg_name(),
            arb_orientation(),
        )
            .prop_map(|(from_seg, from_orient, to_seg, to_orient)| GfaLink {
                from_segment: from_seg,
                from_orient,
                to_segment: to_seg,
                to_orient,
                overlap: "8M".to_string(),
                tags: HashMap::new(),
            })
    }

    proptest! {
        #[test]
        fn gfa_segment_sequence_valid(seg in arb_gfa_segment()) {
            // All characters should be valid DNA bases
            prop_assert!(seg.sequence.chars().all(|c| "ACGT".contains(c)));
        }

        #[test]
        fn gfa_segment_length_tag_consistent(seg in arb_gfa_segment()) {
            let length = seg.length();
            prop_assert_eq!(length, seg.sequence.len());
        }

        #[test]
        fn gfa_link_segments_nonempty(link in arb_gfa_link()) {
            prop_assert!(!link.from_segment.is_empty());
            prop_assert!(!link.to_segment.is_empty());
        }

        #[test]
        fn gfa_link_overlap_nonempty(link in arb_gfa_link()) {
            prop_assert!(!link.overlap.is_empty());
        }
    }
}

// ============================================================================
// VCF Format Property Tests
// ============================================================================

mod vcf_properties {
    use super::*;

    /// Generate arbitrary chromosome
    fn arb_chrom() -> impl Strategy<Value = String> {
        prop_oneof![
            Just("chr1".to_string()),
            Just("chr21".to_string()),
            Just("chrX".to_string()),
        ]
    }

    /// Generate arbitrary nucleotide
    fn arb_base() -> impl Strategy<Value = String> {
        prop::sample::select(vec!["A", "C", "G", "T"]).prop_map(|s| s.to_string())
    }

    /// Generate arbitrary alleles (ref + alts)
    fn arb_alleles() -> impl Strategy<Value = (String, Vec<String>)> {
        (arb_base(), prop::collection::vec(arb_base(), 1..3))
    }

    /// Strategy for generating valid VCF records
    fn arb_vcf_record() -> impl Strategy<Value = VcfRecord> {
        (
            arb_chrom(),
            1u64..100_000_000u64,
            prop::option::of("[a-z0-9]+"),
            arb_alleles(),
            prop::option::of(0.0f64..100.0f64),
            prop::option::of(prop_oneof![Just("PASS"), Just("LowQual")]),
        )
            .prop_map(|(chrom, pos, id, (reference, alternate), quality, filter)| {
                let mut info = HashMap::new();
                if let Some(q) = quality {
                    info.insert("DP".to_string(), format!("{}", (q * 10.0) as u32));
                }

                VcfRecord {
                    chrom,
                    pos,
                    id: id.map(|s| s.to_string()),
                    reference,
                    alternate,
                    quality,
                    filter: filter.map(|s| s.to_string()),
                    info,
                    format: None,
                    samples: vec![],
                }
            })
    }

    proptest! {
        #[test]
        fn vcf_round_trip(record in arb_vcf_record()) {
            let line = record.to_line();
            let parsed = VcfRecord::from_line(&line).unwrap();
            prop_assert_eq!(record, parsed);
        }

        #[test]
        fn vcf_position_positive(record in arb_vcf_record()) {
            prop_assert!(record.pos > 0, "VCF position must be >= 1");
        }

        #[test]
        fn vcf_alleles_nonempty(record in arb_vcf_record()) {
            prop_assert!(!record.reference.is_empty());
            prop_assert!(!record.alternate.is_empty());
        }

        #[test]
        fn vcf_quality_nonnegative(record in arb_vcf_record()) {
            if let Some(qual) = record.quality {
                prop_assert!(qual >= 0.0);
            }
        }

        #[test]
        fn vcf_snp_classification(record in arb_vcf_record()) {
            // SNP = single base ref and single base alt
            let is_single_base = record.reference.len() == 1
                && record.alternate.iter().all(|a| a.len() == 1);

            if is_single_base {
                // This should be classified as SNP (when methods exist)
                prop_assert_eq!(record.reference.len(), 1);
            }
        }
    }
}

// ============================================================================
// GFF3 Format Property Tests
// ============================================================================

mod gff3_properties {
    use super::*;

    /// Generate arbitrary sequence ID
    fn arb_seqid() -> impl Strategy<Value = String> {
        prop_oneof![
            Just("chr1".to_string()),
            Just("chr2".to_string()),
            Just("chrX".to_string()),
        ]
    }

    /// Generate arbitrary feature type
    fn arb_feature_type() -> impl Strategy<Value = String> {
        prop_oneof![
            Just("gene".to_string()),
            Just("mRNA".to_string()),
            Just("exon".to_string()),
            Just("CDS".to_string()),
        ]
    }

    /// Generate arbitrary strand
    fn arb_strand() -> impl Strategy<Value = Strand> {
        prop_oneof![
            Just(Strand::Forward),
            Just(Strand::Reverse),
            Just(Strand::Unknown),
        ]
    }

    /// Strategy for generating valid GFF3 records
    fn arb_gff3_record() -> impl Strategy<Value = Gff3Record> {
        (
            arb_seqid(),
            arb_feature_type(),
            1u64..100_000u64,
            1u64..100_000u64,
            arb_strand(),
            prop::option::of(0u8..3u8),
        )
            .prop_filter("end must be >= start", |(_, _, start, end, _, _)| end >= start)
            .prop_map(|(seqid, feature_type, start, end, strand, phase)| {
                let mut attributes = HashMap::new();
                attributes.insert("ID".to_string(), format!("{}1", feature_type));

                if feature_type != "gene" {
                    attributes.insert("Parent".to_string(), "gene1".to_string());
                }

                Gff3Record {
                    seqid,
                    source: "test".to_string(),
                    feature_type,
                    start,
                    end,
                    score: None,
                    strand,
                    phase,
                    attributes,
                }
            })
    }

    proptest! {
        #[test]
        fn gff3_round_trip(record in arb_gff3_record()) {
            let line = record.to_line();
            let parsed = Gff3Record::from_line(&line).unwrap();
            prop_assert_eq!(record, parsed);
        }

        #[test]
        fn gff3_coordinates_ordered(record in arb_gff3_record()) {
            // GFF3 uses 1-based inclusive coordinates
            prop_assert!(record.end >= record.start);
        }

        #[test]
        fn gff3_length_calculation(record in arb_gff3_record()) {
            let length = record.length();
            // GFF3 length = end - start + 1 (inclusive)
            prop_assert_eq!(length, record.end - record.start + 1);
        }

        #[test]
        fn gff3_phase_range(record in arb_gff3_record()) {
            if let Some(phase) = record.phase {
                prop_assert!(phase < 3, "Phase must be 0, 1, or 2");
            }
        }

        #[test]
        fn gff3_has_id_attribute(record in arb_gff3_record()) {
            prop_assert!(record.get_id().is_some(), "All features should have ID");
        }

        #[test]
        fn gff3_to_0based_conversion(record in arb_gff3_record()) {
            let interval = record.interval().unwrap();

            // GFF3 1-based [start, end] -> BED 0-based [start, end)
            prop_assert_eq!(interval.start, record.start - 1);
            prop_assert_eq!(interval.end, record.end);
        }
    }
}

// ============================================================================
// Cross-Format Coordinate Tests
// ============================================================================

mod coordinate_properties {
    use super::*;

    proptest! {
        #[test]
        fn bed_to_gff_coordinate_conversion(
            chrom in "[a-z]+",
            start in 0u64..1_000_000u64,
            end in 1u64..1_000_000u64
        ) {
            prop_assume!(end > start);

            // BED: 0-based half-open [start, end)
            let bed_interval = GenomicInterval {
                chrom: chrom.clone(),
                start,
                end,
            };

            // Convert to GFF3: 1-based inclusive [start, end]
            let gff_start = bed_interval.start + 1;
            let gff_end = bed_interval.end;

            // Length should be preserved
            let bed_length = bed_interval.length();
            let gff_length = gff_end - gff_start + 1;
            prop_assert_eq!(bed_length, gff_length);
        }

        #[test]
        fn interval_length_always_positive(
            start in 0u64..1_000_000u64,
            end in 1u64..1_000_000u64
        ) {
            prop_assume!(end > start);

            let interval = GenomicInterval {
                chrom: "chr1".to_string(),
                start,
                end,
            };

            prop_assert!(interval.length() > 0);
            prop_assert_eq!(interval.length(), end - start);
        }
    }
}
