from __future__ import annotations

import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "packages" / "core"))

from opengene_core import (  # noqa: E402
    AddFeature,
    CommandHistory,
    DeleteBases,
    Feature,
    InsertBases,
    Location,
    Project,
    ProvenanceNode,
    ReplaceRegion,
    ReverseComplement,
    RotateOrigin,
    SequenceDocument,
    Topology,
)


class ModelAndCommandTests(unittest.TestCase):
    def test_circular_location_has_explicit_segments(self) -> None:
        location = Location(8, 3, strand=-1, wraps_origin=True)
        self.assertEqual(location.segments(10), ((8, 10), (0, 3)))
        self.assertEqual(location.length(10), 5)
        self.assertTrue(location.contains(9, 10))
        self.assertTrue(location.contains(1, 10))
        self.assertFalse(location.contains(5, 10))

    def test_insert_delete_and_undo_remap_annotations(self) -> None:
        feature = Feature.create("CDS", "reporter", [Location(2, 6)])
        document = SequenceDocument.create(
            "vector", "AACCGGTT", topology=Topology.LINEAR, features=[feature]
        )
        history = CommandHistory(document)
        inserted = history.execute(InsertBases(4, "TA"))
        self.assertEqual(inserted.sequence, "AACCTAGGTT")
        self.assertEqual(inserted.features[0].locations, (Location(2, 8),))
        deleted = history.execute(DeleteBases(3, 7))
        self.assertEqual(deleted.sequence, "AACGTT")
        self.assertEqual(deleted.features[0].locations, (Location(2, 4),))
        self.assertEqual(history.undo().sequence, inserted.sequence)
        self.assertEqual(history.undo(), document)
        self.assertEqual(history.redo().sequence, inserted.sequence)

    def test_replace_is_one_revision_and_removes_fully_deleted_feature(self) -> None:
        feature = Feature.create("misc_feature", "temporary", [Location(2, 4)])
        document = SequenceDocument.create("linear", "AACCGG", topology=Topology.LINEAR, features=[feature])
        result = ReplaceRegion(2, 4, "T").apply(document)
        self.assertEqual(result.sequence, "AATGG")
        self.assertEqual(result.revision, document.revision + 1)
        self.assertEqual(result.features, ())

    def test_reverse_complement_and_rotate_preserve_circular_arc(self) -> None:
        feature = Feature.create("CDS", "wrap", [Location(8, 2, strand=1, wraps_origin=True)])
        document = SequenceDocument.create(
            "circle", "ATGCGTACGT", topology=Topology.CIRCULAR, features=[feature]
        )
        reversed_document = ReverseComplement().apply(document)
        self.assertEqual(reversed_document.sequence, "ACGTACGCAT")
        self.assertEqual(
            reversed_document.features[0].locations,
            (Location(8, 2, strand=-1, wraps_origin=True),),
        )
        rotated = RotateOrigin(3).apply(document)
        self.assertEqual(rotated.sequence, "CGTACGTATG")
        self.assertEqual(
            rotated.features[0].locations,
            (Location(5, 9, strand=1),),
        )

    def test_deleting_one_side_of_a_wrapped_feature_makes_it_linear(self) -> None:
        feature = Feature.create("gene", "wrap", [Location(8, 2, wraps_origin=True)])
        document = SequenceDocument.create(
            "circle", "AAAAAAAAAA", topology=Topology.CIRCULAR, features=[feature]
        )
        result = DeleteBases(0, 2).apply(document)
        self.assertEqual(result.sequence, "AAAAAAAA")
        self.assertEqual(result.features[0].locations, (Location(6, 8),))

    def test_json_project_round_trip_preserves_provenance(self) -> None:
        source = SequenceDocument.create("source", "AACCGG", topology=Topology.LINEAR)
        product = SequenceDocument.create("product", "AACCTT", topology=Topology.LINEAR)
        operation = ProvenanceNode.create(
            "replace_region",
            [source.id],
            product.id,
            parameters={"start": 4, "end": 6, "bases": "TT"},
            software_version="test",
        )
        project = Project.create("demo", [source, product], metadata={"owner": "lab"}).with_provenance(operation)
        recovered = Project.from_json(project.to_json())
        self.assertEqual(recovered.to_dict(), project.to_dict())
        self.assertEqual(recovered.document("product").sequence, "AACCTT")

    def test_add_feature_is_undoable(self) -> None:
        document = SequenceDocument.create("doc", "AACCGG", topology=Topology.LINEAR)
        feature = Feature.create("gene", "x", [Location(1, 4)])
        history = CommandHistory(document)
        history.execute(AddFeature(feature))
        self.assertEqual(history.document.features, (feature,))
        self.assertEqual(history.undo().features, ())


if __name__ == "__main__":
    unittest.main()
