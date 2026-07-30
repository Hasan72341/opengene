from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "packages" / "alignment"))
sys.path.insert(0, str(ROOT / "packages" / "sequencing"))

from opengene_sequencing import (  # noqa: E402
    Ab1ParseError,
    Ab1UnavailableError,
    SequencingError,
    SequencingRead,
    ab1_reader_status,
    align_read_to_reference,
    build_consensus,
    decode_phred,
    parse_fastq,
    read_ab1,
    write_fastq,
)


class SequencingReadAndConsensusTests(unittest.TestCase):
    def test_fastq_round_trip_preserves_real_quality_values(self) -> None:
        parsed = parse_fastq("@read one\nACGT\n+\n!I5?\n")

        self.assertEqual(parsed[0].name, "read one")
        self.assertEqual(parsed[0].qualities, (0, 40, 20, 30))
        self.assertEqual(write_fastq(parsed), "@read one\nACGT\n+\n!I5?\n")
        self.assertEqual(parsed[0].reverse_complement().sequence, "ACGT")
        self.assertEqual(parsed[0].reverse_complement().qualities, (30, 20, 40, 0))
        self.assertEqual(decode_phred("I"), (40,))

    def test_fastq_validation_and_missing_quality_export_are_clear(self) -> None:
        with self.assertRaises(SequencingError):
            parse_fastq("@r\nAC\n+\n!\n")
        with self.assertRaises(SequencingError):
            SequencingRead("no_quality", "AC").fastq_quality()
        with self.assertRaises(SequencingError):
            SequencingRead("bad", "AC", qualities=(40,))

    def test_reverse_complement_mapping_uses_oriented_quality_coordinates(self) -> None:
        reference = "TTTAACCGGGTAAA"
        read = SequencingRead("reverse", "CCCGGTT", qualities=(10, 11, 12, 13, 14, 15, 16))

        mapped = align_read_to_reference(reference, read)

        self.assertEqual(mapped.strand, -1)
        self.assertEqual(mapped.oriented_sequence, "AACCGGG")
        self.assertEqual(mapped.oriented_qualities, (16, 15, 14, 13, 12, 11, 10))
        self.assertEqual((mapped.alignment.reference_start, mapped.alignment.reference_end), (3, 10))

    def test_consensus_uses_iupac_for_conflicting_high_quality_observations(self) -> None:
        reads = (
            SequencingRead("c", "ACGT", qualities=(40, 40, 40, 40)),
            SequencingRead("t", "ATGT", qualities=(40, 40, 40, 40)),
        )

        consensus = build_consensus(
            "ACGT",
            reads,
            minimum_base_fraction=0.75,
            ambiguity_fraction=0.2,
        )

        self.assertEqual(consensus.sequence, "AYGT")
        self.assertEqual(consensus.coverage, (2, 2, 2, 2))
        self.assertEqual(consensus.calls[1].base, "Y")
        self.assertIn(1, consensus.ambiguous_positions)

    def test_quality_weighting_prefers_high_confidence_base_and_marks_uncovered_bases(self) -> None:
        reads = (
            SequencingRead("high", "ACGT", qualities=(40, 40, 40, 40)),
            SequencingRead("low", "ATGT", qualities=(40, 0, 40, 40)),
        )

        weighted = build_consensus("ACGTA", reads, minimum_base_fraction=0.75)

        self.assertEqual(weighted.sequence, "ACGTN")
        self.assertEqual(weighted.coverage, (2, 2, 2, 2, 0))
        self.assertGreater(weighted.calls[1].weight_for("C"), weighted.calls[1].weight_for("T"))
        self.assertTrue(any("insufficient" in warning for warning in weighted.warnings))

    def test_ambiguous_read_base_distributes_evidence_without_claiming_a_single_base(self) -> None:
        consensus = build_consensus(
            "A",
            (SequencingRead("ambiguous", "R", qualities=(40,)),),
            minimum_base_fraction=0.75,
            ambiguity_fraction=0.2,
        )

        self.assertEqual(consensus.sequence, "R")
        self.assertAlmostEqual(consensus.calls[0].weight_for("A"), consensus.calls[0].weight_for("G"))

    def test_ab1_boundary_never_fakes_a_parser(self) -> None:
        status = ab1_reader_status()
        with tempfile.TemporaryDirectory() as directory:
            invalid_path = Path(directory) / "not-an-ab1.ab1"
            invalid_path.write_bytes(b"not a binary AB1 chromatogram")
            if status.available:
                with self.assertRaises(Ab1ParseError):
                    read_ab1(invalid_path)
            else:
                with self.assertRaises(Ab1UnavailableError):
                    read_ab1(invalid_path)


if __name__ == "__main__":
    unittest.main()
