from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "packages" / "alignment"))

from opengene_alignment import (  # noqa: E402
    AlignmentError,
    AlignmentMode,
    AlignmentResult,
    AlignmentVariant,
    VariantKind,
    align_global,
    align_local,
    needleman_wunsch,
    smith_waterman,
    variants_from_alignment,
)


class PairwiseAlignmentTests(unittest.TestCase):
    def test_global_alignment_preserves_full_spans_and_deletion_coordinates(self) -> None:
        result = needleman_wunsch("ACTG", "ACG")

        self.assertEqual(result.mode, AlignmentMode.GLOBAL)
        self.assertEqual((result.reference_start, result.reference_end), (0, 4))
        self.assertEqual((result.query_start, result.query_end), (0, 3))
        self.assertEqual(result.aligned_reference, "ACTG")
        self.assertEqual(result.aligned_query, "AC-G")
        self.assertEqual(result.cigar, "2=1D1=")
        self.assertEqual(
            result.variants,
            (
                AlignmentVariant(
                    VariantKind.DELETION,
                    reference_start=2,
                    reference_end=3,
                    reference_bases="T",
                    alternate_bases="",
                    query_start=2,
                    query_end=2,
                ),
            ),
        )

    def test_local_alignment_reports_only_best_source_intervals(self) -> None:
        result = smith_waterman("TTGATTACAGG", "CCGATTACATT")

        self.assertEqual(result.mode, AlignmentMode.LOCAL)
        self.assertEqual(result.aligned_reference, "GATTACA")
        self.assertEqual(result.aligned_query, "GATTACA")
        self.assertEqual((result.reference_start, result.reference_end), (2, 9))
        self.assertEqual((result.query_start, result.query_end), (2, 9))
        self.assertEqual(result.score, 14)
        self.assertEqual(result.variants, ())

    def test_variant_calling_groups_adjacent_events_and_keeps_indels_separate(self) -> None:
        result = AlignmentResult(
            mode="global",
            score=0,
            aligned_reference="AC--GTTAC",
            aligned_query="ATTTG--AT",
            reference_start=10,
            reference_end=17,
            query_start=20,
            query_end=27,
            matches=3,
            mismatches=2,
            gap_columns=4,
        )

        self.assertEqual(
            variants_from_alignment(result),
            (
                AlignmentVariant(
                    "substitution",
                    reference_start=11,
                    reference_end=12,
                    reference_bases="C",
                    alternate_bases="T",
                    query_start=21,
                    query_end=22,
                ),
                AlignmentVariant(
                    "insertion",
                    reference_start=12,
                    reference_end=12,
                    reference_bases="",
                    alternate_bases="TT",
                    query_start=22,
                    query_end=24,
                ),
                AlignmentVariant(
                    "deletion",
                    reference_start=13,
                    reference_end=15,
                    reference_bases="TT",
                    alternate_bases="",
                    query_start=25,
                    query_end=25,
                ),
                AlignmentVariant(
                    "substitution",
                    reference_start=16,
                    reference_end=17,
                    reference_bases="C",
                    alternate_bases="T",
                    query_start=26,
                    query_end=27,
                ),
            ),
        )

    def test_ambiguity_matching_is_explicit_and_affects_statistics(self) -> None:
        exact = needleman_wunsch("AN", "AG")
        compatible = needleman_wunsch("AN", "AG", ambiguity_matches=True)

        self.assertEqual(exact.score, 1)
        self.assertEqual((exact.matches, exact.mismatches), (1, 1))
        self.assertEqual(compatible.score, 4)
        self.assertEqual((compatible.matches, compatible.mismatches), (2, 0))
        self.assertEqual(compatible.variants, ())

    def test_empty_alignment_and_function_aliases_are_valid(self) -> None:
        global_result = align_global("", "AC")
        local_result = align_local("AAAA", "TTTT")

        self.assertEqual(global_result.aligned_reference, "--")
        self.assertEqual(global_result.aligned_query, "AC")
        self.assertEqual(global_result.cigar, "2I")
        self.assertEqual(local_result.alignment_length, 0)
        self.assertEqual(local_result.score, 0)
        self.assertEqual(local_result.variants, ())

    def test_result_validation_rejects_inconsistent_coordinates_and_statistics(self) -> None:
        with self.assertRaises(AlignmentError):
            AlignmentResult(
                mode="global",
                score=0,
                aligned_reference="AC",
                aligned_query="AT",
                reference_start=0,
                reference_end=2,
                query_start=0,
                query_end=2,
                matches=2,
                mismatches=0,
                gap_columns=0,
            )
        with self.assertRaises(AlignmentError):
            needleman_wunsch("A-C", "AC")
        with self.assertRaises(AlignmentError):
            smith_waterman("A", "T", gap_score=1)


if __name__ == "__main__":
    unittest.main()
