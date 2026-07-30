from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from opengene import Project
from opengene.cloning import gibson_assembly
from opengene_core import SequenceDocument, Topology


class AutomationApiTests(unittest.TestCase):
    def test_gibson_product_becomes_a_project_provenance_node(self) -> None:
        left = SequenceDocument.create("left", "AAAAACCCCCGGGGGTTTTT", topology=Topology.LINEAR)
        right = SequenceDocument.create("right", "GGGGGTTTTTCCCCCAAAAA", topology=Topology.LINEAR)
        project = Project.create("automation", [left, right])

        product = gibson_assembly(
            [project.sequence("left"), project.sequence("right")],
            minimum_overlap=10,
            circular_product=False,
            product_name="assembled",
        )
        project.add(product)

        self.assertEqual(project.sequence("assembled").sequence, "AAAAACCCCCGGGGGTTTTTCCCCCAAAAA")
        self.assertEqual(len(project.provenance), 1)
        self.assertEqual(project.provenance[0].operation_type, "gibson_assembly")

    def test_ogs_save_and_open_are_local_and_validated(self) -> None:
        project = Project.create("empty")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "automation.ogs"
            self.assertEqual(project.save(path), path)
            reopened = Project.open(path)
            self.assertEqual(reopened.core.name, "empty")


if __name__ == "__main__":
    unittest.main()
