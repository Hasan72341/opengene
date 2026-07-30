from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "packages" / "io"))
sys.path.insert(0, str(ROOT / "packages" / "core"))

from opengene_core import Feature, Location, SequenceDocument, Topology  # noqa: E402
from opengene_io import parse_fasta, parse_genbank, write_fasta, write_genbank  # noqa: E402


class InterchangeTests(unittest.TestCase):
    def test_fasta_round_trip(self) -> None:
        records = parse_fasta(
            ">alpha useful construct\nAACCGG\n>beta\nTTAA\n",
            topology=Topology.LINEAR,
        )
        self.assertEqual([record.name for record in records], ["alpha", "beta"])
        self.assertEqual(records[0].metadata["description"], "useful construct")
        encoded = write_fasta(records, line_width=4)
        self.assertIn(">alpha useful construct\nAACC\nGG", encoded)
        self.assertEqual([record.sequence for record in parse_fasta(encoded)], ["AACCGG", "TTAA"])

    def test_genbank_circular_feature_round_trip(self) -> None:
        document = SequenceDocument.create(
            "plasmid",
            "ATGCGTACGT",
            topology=Topology.CIRCULAR,
            features=[
                Feature.create(
                    "CDS",
                    "origin-spanning gene",
                    [Location(8, 2, strand=-1, wraps_origin=True)],
                    {"label": "origin-spanning gene", "note": "round trip"},
                )
            ],
            metadata={"description": "example circular record", "accession": "OGS0001"},
        )
        encoded = write_genbank(document)
        self.assertIn("complement(join(9..10,1..2))", encoded)
        recovered = parse_genbank(encoded)[0]
        self.assertEqual(recovered.sequence, document.sequence)
        self.assertEqual(recovered.topology, Topology.CIRCULAR)
        feature = recovered.features[0]
        self.assertEqual(feature.locations, (Location(8, 2, strand=-1, wraps_origin=True),))
        self.assertEqual(feature.qualifiers["label"], "origin-spanning gene")

    def test_parse_common_genbank_locations_and_multiline_qualifier(self) -> None:
        source = """LOCUS       example                 12 bp    DNA    linear    01-JAN-1980
DEFINITION  location parser test.
FEATURES             Location/Qualifiers
     gene            complement(2..6)
                     /gene=\"reverse_gene\"
     misc_feature    join(8..10,11..12)
                     /note=\"joined annotation that is intentionally long and
                     continues here\"
ORIGIN
        1 aaccggttaacc
//
"""
        document = parse_genbank(source)[0]
        self.assertEqual(document.features[0].locations, (Location(1, 6, strand=-1),))
        self.assertEqual(document.features[1].locations, (Location(7, 10), Location(10, 12)))
        self.assertEqual(
            document.features[1].qualifiers["note"],
            "joined annotation that is intentionally long andcontinues here",
        )


if __name__ == "__main__":
    unittest.main()
