from __future__ import annotations

from hashlib import sha256
import json
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path
import zipfile


ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "packages" / "io"))
sys.path.insert(0, str(ROOT / "packages" / "core"))

from opengene_core import Feature, Location, Project, ProvenanceNode, SequenceDocument, Topology  # noqa: E402
from opengene_io import OGS_FORMAT, OgsError, read_ogs, write_ogs  # noqa: E402


def _canonical_json(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def _rewrite_archive(
    source: Path,
    destination: Path,
    replacements: dict[str, bytes],
    *,
    refresh_checksums: bool,
) -> None:
    with zipfile.ZipFile(source) as archive:
        contents = {info.filename: archive.read(info.filename) for info in archive.infolist() if not info.is_dir()}
    contents.update(replacements)
    if refresh_checksums:
        checksums = json.loads(contents["checksums.json"])
        checksums["members"] = {
            name: sha256(data).hexdigest()
            for name, data in sorted(contents.items())
            if name != "checksums.json"
        }
        contents["checksums.json"] = _canonical_json(checksums)
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for name, data in contents.items():
            archive.writestr(name, data)


class OgsContainerTests(unittest.TestCase):
    def _project(self) -> Project:
        source = SequenceDocument.create(
            "source_vector",
            "AACCGGTT",
            topology=Topology.LINEAR,
            features=[Feature.create("CDS", "reporter", [Location(1, 7)], {"label": "reporter"})],
            metadata={"owner": "lab", "tags": ["test", "source"]},
        )
        duplicate = SequenceDocument.create("same_sequence", "AACCGGTT", topology=Topology.LINEAR)
        product = SequenceDocument.create("product", "AACCTT", topology=Topology.LINEAR)
        provenance = ProvenanceNode.create(
            "replace_region",
            [source.id],
            product.id,
            parameters={"start": 4, "end": 6, "bases": "TT"},
            software_version="test-suite",
            engine_version="core-1",
        )
        return Project.create(
            "integrity fixture",
            [source, duplicate, product],
            metadata={"purpose": "OGS round trip"},
        ).with_provenance(provenance)

    def test_round_trip_writes_real_sqlite_and_content_addressed_blobs(self) -> None:
        project = self._project()
        with tempfile.TemporaryDirectory() as directory:
            archive_path = Path(directory) / "fixture.ogs"
            write_ogs(archive_path, project, application_version="test-1", created_at="2026-07-23T00:00:00Z")
            recovered = read_ogs(archive_path)
            self.assertEqual(recovered.to_dict(), project.to_dict())

            with zipfile.ZipFile(archive_path) as archive:
                names = archive.namelist()
                self.assertEqual(names[0:2], ["manifest.json", "project.sqlite"])
                self.assertEqual(names[-1], "checksums.json")
                manifest = json.loads(archive.read("manifest.json"))
                self.assertEqual(manifest["format"], OGS_FORMAT)
                self.assertEqual(manifest["documentIds"], [str(document.id) for document in project.documents])
                sequence_members = [name for name in names if name.startswith("sequences/")]
                self.assertEqual(len(sequence_members), 2)
                checksums = json.loads(archive.read("checksums.json"))
                self.assertEqual(set(checksums["members"]), set(names) - {"checksums.json"})

                connection = sqlite3.connect(":memory:")
                try:
                    connection.deserialize(archive.read("project.sqlite"))
                    self.assertEqual(connection.execute("SELECT COUNT(*) FROM documents").fetchone()[0], 3)
                    self.assertEqual(connection.execute("SELECT COUNT(*) FROM provenance").fetchone()[0], 1)
                finally:
                    connection.close()

    def test_same_project_has_byte_stable_default_archive_output(self) -> None:
        project = self._project()
        with tempfile.TemporaryDirectory() as directory:
            first = Path(directory) / "one.ogs"
            second = Path(directory) / "two.ogs"
            write_ogs(first, project)
            write_ogs(second, project)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_checksum_tampering_is_rejected_before_project_reconstruction(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            original = Path(directory) / "original.ogs"
            tampered = Path(directory) / "tampered.ogs"
            write_ogs(original, self._project())
            with zipfile.ZipFile(original) as archive:
                sequence_member = next(name for name in archive.namelist() if name.startswith("sequences/"))
                changed = archive.read(sequence_member) + b"A"
            _rewrite_archive(original, tampered, {sequence_member: changed}, refresh_checksums=False)
            with self.assertRaisesRegex(OgsError, "checksum mismatch"):
                read_ogs(tampered)

    def test_unknown_manifest_version_is_rejected_even_when_checksums_are_valid(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            original = Path(directory) / "original.ogs"
            incompatible = Path(directory) / "incompatible.ogs"
            write_ogs(original, self._project())
            with zipfile.ZipFile(original) as archive:
                manifest = json.loads(archive.read("manifest.json"))
            manifest["formatVersion"] = 999
            _rewrite_archive(
                original,
                incompatible,
                {"manifest.json": _canonical_json(manifest)},
                refresh_checksums=True,
            )
            with self.assertRaisesRegex(OgsError, "unsupported formatVersion"):
                read_ogs(incompatible)


if __name__ == "__main__":
    unittest.main()
