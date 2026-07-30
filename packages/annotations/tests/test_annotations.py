from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "packages" / "annotations"))

from opengene_annotations import (  # noqa: E402
    AnnotationError,
    AnnotationLocation,
    GENERIC_CODON_USAGE,
    detect_orf_features,
    find_orfs,
    optimize_codon_usage,
    reverse_translate,
    translate,
    translate_dna,
)


class TranslationTests(unittest.TestCase):
    def test_translation_supports_frames_reverse_strand_and_ambiguity(self) -> None:
        result = translate_dna("AATGAAATAGC", frame=1, stop_at_stop=True)
        self.assertEqual(result.protein_sequence, "MK*")
        self.assertTrue(result.terminated)
        self.assertEqual(result.trailing_bases, "C")
        self.assertEqual(translate("TTATTTCAT", strand=-1), "MK*")
        self.assertEqual(translate("ATH"), "I")
        self.assertEqual(translate("ATN"), "X")

    def test_translation_can_hide_and_stop_at_termination(self) -> None:
        result = translate_dna("ATGTAATGG", include_stop=False, stop_at_stop=True)
        self.assertEqual(result.protein_sequence, "M")
        self.assertEqual(result.codons, ("ATG", "TAA"))
        self.assertEqual(result.trailing_bases, "TGG")


class ORFTests(unittest.TestCase):
    def test_linear_orfs_are_found_on_both_strands(self) -> None:
        forward = find_orfs("ATGAAATAA")
        self.assertEqual(len(forward), 1)
        self.assertEqual(forward[0].location, AnnotationLocation(0, 9, strand=1))
        self.assertEqual(forward[0].protein_sequence, "MK")
        self.assertEqual(forward[0].stop_codon, "TAA")
        self.assertEqual(forward[0].reading_frame, 1)

        reverse = find_orfs("TTATTTCAT")
        self.assertEqual(len(reverse), 1)
        self.assertEqual(reverse[0].location, AnnotationLocation(0, 9, strand=-1))
        self.assertEqual(reverse[0].protein_sequence, "MK")
        self.assertEqual(reverse[0].reading_frame, -1)

    def test_circular_orf_can_cross_origin_without_linear_false_positive(self) -> None:
        sequence = "TAACCCATGAAA"
        self.assertEqual(find_orfs(sequence), ())
        circular = find_orfs(sequence, topology="circular")
        self.assertEqual(len(circular), 1)
        orf = circular[0]
        self.assertEqual(orf.location, AnnotationLocation(6, 3, strand=1, wraps_origin=True))
        self.assertEqual(orf.nucleotide_sequence, "ATGAAATAA")
        self.assertEqual(orf.protein_sequence, "MK")

    def test_partial_orfs_are_opt_in_and_feature_results_keep_evidence(self) -> None:
        self.assertEqual(find_orfs("ATGCCC"), ())
        partial = find_orfs("ATGCCC", include_partial=True)
        self.assertEqual(len(partial), 1)
        self.assertFalse(partial[0].has_terminal_stop)
        detected = detect_orf_features("ATGCCC", include_partial=True)
        self.assertEqual(detected.features[0].feature_type, "CDS")
        self.assertEqual(detected.features[0].qualifiers["partial"], True)
        self.assertLess(detected.features[0].confidence, 1.0)

    def test_circular_whole_molecule_orfs_are_conservatively_excluded(self) -> None:
        self.assertEqual(find_orfs("ATGAAATAA", topology="circular"), ())


class CodonOptimisationTests(unittest.TestCase):
    def test_generic_reverse_translation_is_explicit_and_deterministic(self) -> None:
        self.assertEqual(reverse_translate("MAK*"), "ATGGCTAAATAA")
        result = optimize_codon_usage("MAK*")
        self.assertEqual(result.codons, ("ATG", "GCT", "AAA", "TAA"))
        self.assertEqual(result.usage_table_name, "generic-balanced-v1")
        self.assertEqual(GENERIC_CODON_USAGE["A"][0].codon, "GCT")

    def test_forbidden_motif_can_select_a_lower_weight_synonym_or_report_warning(self) -> None:
        result = optimize_codon_usage("AA", forbidden_motifs=("GCTGCT",))
        self.assertEqual(result.dna_sequence, "GCTGCC")
        impossible = optimize_codon_usage("M", forbidden_motifs=("ATG",))
        self.assertEqual(impossible.dna_sequence, "ATG")
        self.assertEqual(len(impossible.warnings), 1)

    def test_invalid_input_is_rejected_rather_than_silently_guessed(self) -> None:
        with self.assertRaises(AnnotationError):
            reverse_translate("MX")
        with self.assertRaises(AnnotationError):
            find_orfs("ATGAAA", min_protein_length=-1)


if __name__ == "__main__":
    unittest.main()
