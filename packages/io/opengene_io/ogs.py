from __future__ import annotations

from collections.abc import Mapping
from datetime import datetime
from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath
import re
import sqlite3
import tempfile
from typing import Any
from uuid import UUID
import zipfile

from opengene_core import Project


OGS_FORMAT = "OpenGene Studio project"
OGS_FORMAT_VERSION = 1
OGS_HASH_ALGORITHM = "sha256"
OGS_MANIFEST_MEMBER = "manifest.json"
OGS_DATABASE_MEMBER = "project.sqlite"
OGS_CHECKSUM_MEMBER = "checksums.json"
_CANONICAL_CREATED_AT = "1980-01-01T00:00:00Z"
_ZIP_EPOCH = (1980, 1, 1, 0, 0, 0)
_SHA256_HEX = re.compile(r"^[0-9a-f]{64}$")


class OgsError(ValueError):
    pass


def _canonical_json(value: Any, *, label: str) -> bytes:

    try:
        return json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise OgsError(f"{label} is not strict JSON serializable: {error}") from error


def _json_object(data: bytes, *, label: str) -> dict[str, Any]:
    try:
        decoded = data.decode("utf-8")
        value = json.loads(decoded)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise OgsError(f"{label} is not valid UTF-8 JSON: {error}") from error
    if not isinstance(value, dict):
        raise OgsError(f"{label} must contain a JSON object")
    return value


def _sha256(data: bytes) -> str:
    return sha256(data).hexdigest()


def _sequence_member(sequence_hash: str) -> str:
    return f"sequences/sha256-{sequence_hash}.seq"


def _require_string(value: Any, *, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise OgsError(f"{label} must be a non-empty string")
    return value


def _require_exact_int(value: Any, *, label: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise OgsError(f"{label} must be an integer")
    return value


def _validate_member_name(name: str) -> None:

    if not isinstance(name, str) or not name or "\\" in name:
        raise OgsError("archive contains an invalid member name")
    path = PurePosixPath(name)
    if path.is_absolute() or ".." in path.parts or name.startswith("/"):
        raise OgsError(f"archive contains unsafe member name {name!r}")


def _archive_members(archive: zipfile.ZipFile) -> dict[str, bytes]:
    infos = archive.infolist()
    names: list[str] = []
    for info in infos:
        _validate_member_name(info.filename)
        if info.flag_bits & 0x1:
            raise OgsError(f"encrypted archive member {info.filename!r} is not supported")
        if not info.is_dir():
            names.append(info.filename)
    if len(names) != len(set(names)):
        raise OgsError("archive contains duplicate member names")
    required = {OGS_MANIFEST_MEMBER, OGS_DATABASE_MEMBER, OGS_CHECKSUM_MEMBER}
    missing = required - set(names)
    if missing:
        raise OgsError(f"archive is missing required member(s): {', '.join(sorted(missing))}")
    try:
        return {name: archive.read(name) for name in names}
    except (KeyError, RuntimeError, zipfile.BadZipFile, OSError) as error:
        raise OgsError(f"could not read .ogs archive members: {error}") from error


def _validate_checksums(members: Mapping[str, bytes]) -> None:
    checksums = _json_object(members[OGS_CHECKSUM_MEMBER], label=OGS_CHECKSUM_MEMBER)
    if checksums.get("format") != "opengene.checksums":
        raise OgsError("checksums.json has an unsupported format")
    if _require_exact_int(checksums.get("formatVersion"), label="checksums formatVersion") != OGS_FORMAT_VERSION:
        raise OgsError("checksums.json has an unsupported formatVersion")
    if checksums.get("hashAlgorithm") != OGS_HASH_ALGORITHM:
        raise OgsError("checksums.json must use sha256")
    declared = checksums.get("members")
    if not isinstance(declared, Mapping):
        raise OgsError("checksums.json members must be an object")
    actual_names = set(members) - {OGS_CHECKSUM_MEMBER}
    declared_names = set(declared)
    if declared_names != actual_names:
        missing = actual_names - declared_names
        extra = declared_names - actual_names
        detail: list[str] = []
        if missing:
            detail.append(f"missing: {', '.join(sorted(missing))}")
        if extra:
            detail.append(f"unknown: {', '.join(sorted(extra))}")
        raise OgsError("checksums.json member list does not match archive (" + "; ".join(detail) + ")")
    for name in sorted(actual_names):
        expected = declared[name]
        if not isinstance(expected, str) or not _SHA256_HEX.fullmatch(expected):
            raise OgsError(f"checksums.json has an invalid sha256 for {name!r}")
        actual = _sha256(members[name])
        if actual != expected:
            raise OgsError(f"checksum mismatch for archive member {name!r}")


def _validate_manifest(manifest: Mapping[str, Any]) -> None:
    if manifest.get("format") != OGS_FORMAT:
        raise OgsError("manifest.json is not an OpenGene Studio project")
    if _require_exact_int(manifest.get("formatVersion"), label="manifest formatVersion") != OGS_FORMAT_VERSION:
        raise OgsError("manifest.json has an unsupported formatVersion")
    if manifest.get("hashAlgorithm") != OGS_HASH_ALGORITHM:
        raise OgsError("manifest.json must use sha256")
    _require_string(manifest.get("applicationVersion"), label="manifest applicationVersion")
    created_at = _require_string(manifest.get("createdAt"), label="manifest createdAt")
    try:
        datetime.fromisoformat(created_at.replace("Z", "+00:00"))
    except ValueError as error:
        raise OgsError("manifest createdAt must be an ISO-8601 timestamp") from error
    try:
        UUID(_require_string(manifest.get("projectId"), label="manifest projectId"))
    except ValueError as error:
        raise OgsError("manifest projectId must be a UUID") from error
    document_ids = manifest.get("documentIds")
    if not isinstance(document_ids, list) or not all(isinstance(item, str) for item in document_ids):
        raise OgsError("manifest documentIds must be a unique list")
    if len(document_ids) != len(set(document_ids)):
        raise OgsError("manifest documentIds must be a unique list")
    try:
        for document_id in document_ids:
            UUID(str(document_id))
    except ValueError as error:
        raise OgsError("manifest documentIds must contain UUID values") from error
    if not isinstance(manifest.get("pluginRequirements"), list):
        raise OgsError("manifest pluginRequirements must be a list")
    if not isinstance(manifest.get("databaseVersions"), Mapping):
        raise OgsError("manifest databaseVersions must be an object")
    _require_exact_int(manifest.get("projectSchemaVersion"), label="manifest projectSchemaVersion")
    sequence_members = manifest.get("sequenceMembers")
    if not isinstance(sequence_members, Mapping):
        raise OgsError("manifest sequenceMembers must be an object")
    if set(sequence_members) != set(document_ids):
        raise OgsError("manifest sequenceMembers must map every document id")
    for document_id, member in sequence_members.items():
        _validate_member_name(str(member))
        if not str(member).startswith("sequences/sha256-"):
            raise OgsError(f"manifest sequence member for {document_id!r} is not content addressed")


def _archive_payload(project: Project) -> tuple[dict[str, Any], dict[str, bytes]]:

    if not isinstance(project, Project):
        raise TypeError("project must be an opengene_core.Project")
    payload = _json_object(_canonical_json(project.to_dict(), label="project"), label="project")
    documents = payload.get("documents")
    if not isinstance(documents, list):
        raise OgsError("project serialization has no documents list")
    blobs: dict[str, bytes] = {}
    for document in documents:
        if not isinstance(document, dict):
            raise OgsError("project serialization contains an invalid document")
        try:
            sequence = document.pop("sequence")
            sequence_hash = document["sequence_sha256"]
            document_id = document["id"]
        except KeyError as error:
            raise OgsError(f"project document is missing {error.args[0]!r}") from error
        if not isinstance(sequence, str) or not isinstance(sequence_hash, str):
            raise OgsError("project document sequence data is invalid")
        try:
            sequence_bytes = sequence.encode("ascii")
        except UnicodeEncodeError as error:
            raise OgsError(f"sequence for document {document_id!r} is not ASCII") from error
        actual_hash = _sha256(sequence_bytes)
        if actual_hash != sequence_hash:
            raise OgsError(f"project sequence hash mismatch for document {document_id!r}")
        member = _sequence_member(sequence_hash)
        document["sequence_member"] = member
        blobs.setdefault(member, sequence_bytes)
    return payload, blobs


def _build_sqlite(payload: Mapping[str, Any]) -> bytes:

    documents = payload.get("documents")
    provenance = payload.get("provenance")
    if not isinstance(documents, list) or not isinstance(provenance, list):
        raise OgsError("project payload is missing documents or provenance")
    project_id = _require_string(payload.get("id"), label="project id")
    project_name = _require_string(payload.get("name"), label="project name")
    schema_version = _require_exact_int(payload.get("schema_version"), label="project schema_version")
    payload_json = _canonical_json(payload, label="project payload").decode("utf-8")
    connection = sqlite3.connect(":memory:")
    try:
        connection.execute("PRAGMA page_size = 4096")
        connection.execute("PRAGMA auto_vacuum = NONE")
        connection.execute("PRAGMA application_id = 1330667347")
        connection.execute(f"PRAGMA user_version = {OGS_FORMAT_VERSION}")
        connection.executescript(
            """
            CREATE TABLE project (
                singleton INTEGER PRIMARY KEY CHECK (singleton = 1),
                project_id TEXT NOT NULL,
                name TEXT NOT NULL,
                schema_version INTEGER NOT NULL,
                project_json TEXT NOT NULL
            );
            CREATE TABLE documents (
                ordinal INTEGER PRIMARY KEY,
                document_id TEXT NOT NULL UNIQUE,
                name TEXT NOT NULL,
                molecule_type TEXT NOT NULL,
                topology TEXT NOT NULL,
                sequence_sha256 TEXT NOT NULL,
                sequence_member TEXT NOT NULL,
                revision INTEGER NOT NULL,
                document_json TEXT NOT NULL
            );
            CREATE TABLE provenance (
                ordinal INTEGER PRIMARY KEY,
                operation_id TEXT NOT NULL UNIQUE,
                output_document_id TEXT NOT NULL,
                operation_type TEXT NOT NULL,
                node_json TEXT NOT NULL
            );
            """
        )
        connection.execute(
            "INSERT INTO project VALUES (1, ?, ?, ?, ?)",
            (project_id, project_name, schema_version, payload_json),
        )
        for ordinal, document in enumerate(documents):
            if not isinstance(document, Mapping):
                raise OgsError("project payload contains an invalid document")
            try:
                values = (
                    ordinal,
                    _require_string(document["id"], label="document id"),
                    _require_string(document["name"], label="document name"),
                    _require_string(document["molecule_type"], label="document molecule_type"),
                    _require_string(document["topology"], label="document topology"),
                    _require_string(document["sequence_sha256"], label="document sequence_sha256"),
                    _require_string(document["sequence_member"], label="document sequence_member"),
                    _require_exact_int(document["revision"], label="document revision"),
                    _canonical_json(document, label="document payload").decode("utf-8"),
                )
            except KeyError as error:
                raise OgsError(f"document payload is missing {error.args[0]!r}") from error
            connection.execute(
                "INSERT INTO documents VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                values,
            )
        for ordinal, node in enumerate(provenance):
            if not isinstance(node, Mapping):
                raise OgsError("project payload contains an invalid provenance node")
            try:
                values = (
                    ordinal,
                    _require_string(node["operation_id"], label="operation id"),
                    _require_string(node["output_document_id"], label="output document id"),
                    _require_string(node["operation_type"], label="operation type"),
                    _canonical_json(node, label="provenance payload").decode("utf-8"),
                )
            except KeyError as error:
                raise OgsError(f"provenance payload is missing {error.args[0]!r}") from error
            connection.execute(
                "INSERT INTO provenance VALUES (?, ?, ?, ?, ?)",
                values,
            )
        connection.commit()
        return bytes(connection.serialize())
    except (sqlite3.DatabaseError, sqlite3.ProgrammingError) as error:
        raise OgsError(f"could not build project.sqlite: {error}") from error
    finally:
        connection.close()


def _read_sqlite(data: bytes) -> dict[str, Any]:

    connection = sqlite3.connect(":memory:")
    try:
        connection.deserialize(data)
        connection.execute("PRAGMA query_only = ON")
        required_tables = {"project", "documents", "provenance"}
        tables = {
            row[0]
            for row in connection.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
        }
        missing = required_tables - tables
        if missing:
            raise OgsError(f"project.sqlite is missing table(s): {', '.join(sorted(missing))}")
        version = connection.execute("PRAGMA user_version").fetchone()[0]
        if version != OGS_FORMAT_VERSION:
            raise OgsError("project.sqlite has an unsupported user_version")
        project_rows = connection.execute(
            "SELECT singleton, project_id, name, schema_version, project_json FROM project ORDER BY singleton"
        ).fetchall()
        if len(project_rows) != 1:
            raise OgsError("project.sqlite must contain exactly one project row")
        singleton, project_id, project_name, schema_version, project_json = project_rows[0]
        if singleton != 1:
            raise OgsError("project.sqlite project row has an invalid singleton key")
        payload = _json_object(str(project_json).encode("utf-8"), label="project.sqlite project_json")
        if (
            payload.get("id") != project_id
            or payload.get("name") != project_name
            or payload.get("schema_version") != schema_version
        ):
            raise OgsError("project.sqlite project row does not match project_json")
        documents = payload.get("documents")
        provenance = payload.get("provenance")
        if not isinstance(documents, list) or not isinstance(provenance, list):
            raise OgsError("project.sqlite project_json is missing documents or provenance")

        document_rows = connection.execute(
            "SELECT ordinal, document_id, name, molecule_type, topology, sequence_sha256, "
            "sequence_member, revision, document_json FROM documents ORDER BY ordinal"
        ).fetchall()
        if len(document_rows) != len(documents):
            raise OgsError("project.sqlite document row count does not match project_json")
        for ordinal, (row, document) in enumerate(zip(document_rows, documents)):
            if not isinstance(document, Mapping):
                raise OgsError("project.sqlite project_json contains an invalid document")
            row_ordinal, document_id, name, molecule_type, topology, sequence_hash, sequence_member, revision, raw_json = row
            if row_ordinal != ordinal:
                raise OgsError("project.sqlite document ordinals are not contiguous")
            expected = (
                document.get("id"),
                document.get("name"),
                document.get("molecule_type"),
                document.get("topology"),
                document.get("sequence_sha256"),
                document.get("sequence_member"),
                document.get("revision"),
            )
            if (document_id, name, molecule_type, topology, sequence_hash, sequence_member, revision) != expected:
                raise OgsError("project.sqlite document row does not match project_json")
            if _json_object(str(raw_json).encode("utf-8"), label="project.sqlite document_json") != document:
                raise OgsError("project.sqlite document_json does not match project_json")

        provenance_rows = connection.execute(
            "SELECT ordinal, operation_id, output_document_id, operation_type, node_json "
            "FROM provenance ORDER BY ordinal"
        ).fetchall()
        if len(provenance_rows) != len(provenance):
            raise OgsError("project.sqlite provenance row count does not match project_json")
        for ordinal, (row, node) in enumerate(zip(provenance_rows, provenance)):
            if not isinstance(node, Mapping):
                raise OgsError("project.sqlite project_json contains an invalid provenance node")
            row_ordinal, operation_id, output_document_id, operation_type, raw_json = row
            if row_ordinal != ordinal:
                raise OgsError("project.sqlite provenance ordinals are not contiguous")
            expected = (node.get("operation_id"), node.get("output_document_id"), node.get("operation_type"))
            if (operation_id, output_document_id, operation_type) != expected:
                raise OgsError("project.sqlite provenance row does not match project_json")
            if _json_object(str(raw_json).encode("utf-8"), label="project.sqlite node_json") != node:
                raise OgsError("project.sqlite node_json does not match project_json")
        return payload
    except OgsError:
        raise
    except (sqlite3.DatabaseError, sqlite3.ProgrammingError, UnicodeEncodeError) as error:
        raise OgsError(f"could not read project.sqlite: {error}") from error
    finally:
        connection.close()


def _hydrate_project(
    payload: Mapping[str, Any], members: Mapping[str, bytes], manifest: Mapping[str, Any]
) -> Project:

    documents = payload.get("documents")
    if not isinstance(documents, list):
        raise OgsError("project payload has no documents list")
    manifest_document_ids = manifest["documentIds"]
    manifest_members = manifest["sequenceMembers"]
    document_ids = [document.get("id") if isinstance(document, Mapping) else None for document in documents]
    if document_ids != manifest_document_ids:
        raise OgsError("manifest documentIds do not match project.sqlite")
    restored_documents: list[dict[str, Any]] = []
    expected_blob_members: set[str] = set()
    for document in documents:
        if not isinstance(document, Mapping):
            raise OgsError("project payload contains an invalid document")
        restored = dict(document)
        document_id = _require_string(restored.get("id"), label="document id")
        sequence_hash = _require_string(restored.get("sequence_sha256"), label="document sequence_sha256")
        if not _SHA256_HEX.fullmatch(sequence_hash):
            raise OgsError(f"document {document_id!r} has an invalid sequence sha256")
        member = _require_string(restored.pop("sequence_member", None), label="document sequence_member")
        expected_member = _sequence_member(sequence_hash)
        if member != expected_member or manifest_members.get(document_id) != member:
            raise OgsError(f"document {document_id!r} has a non-canonical sequence member")
        if member not in members:
            raise OgsError(f"document {document_id!r} is missing its sequence blob")
        sequence_bytes = members[member]
        if _sha256(sequence_bytes) != sequence_hash:
            raise OgsError(f"sequence blob checksum mismatch for document {document_id!r}")
        try:
            restored["sequence"] = sequence_bytes.decode("ascii")
        except UnicodeDecodeError as error:
            raise OgsError(f"sequence blob for document {document_id!r} is not ASCII") from error
        expected_blob_members.add(member)
        restored_documents.append(restored)
    archive_blob_members = {name for name in members if name.startswith("sequences/")}
    if archive_blob_members != expected_blob_members:
        raise OgsError("archive sequence blobs do not exactly match project documents")
    restored_payload = dict(payload)
    restored_payload["documents"] = restored_documents
    try:
        project = Project.from_dict(restored_payload)
    except ValueError as error:
        raise OgsError(f"project data is invalid: {error}") from error
    return project


def _zip_info(name: str) -> zipfile.ZipInfo:
    info = zipfile.ZipInfo(name, date_time=_ZIP_EPOCH)
    info.compress_type = zipfile.ZIP_DEFLATED
    info.create_system = 3
    info.external_attr = 0o600 << 16
    return info


def _write_member(archive: zipfile.ZipFile, name: str, data: bytes) -> None:
    archive.writestr(_zip_info(name), data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)


def _normalise_created_at(created_at: str | None) -> str:
    if created_at is None:
        return _CANONICAL_CREATED_AT
    if not isinstance(created_at, str) or not created_at.strip():
        raise OgsError("created_at must be a non-empty ISO-8601 string")
    try:
        datetime.fromisoformat(created_at.replace("Z", "+00:00"))
    except ValueError as error:
        raise OgsError("created_at must be an ISO-8601 timestamp") from error
    return created_at


def write_ogs(
    path: str | Path,
    project: Project,
    *,
    application_version: str = "0.1.0",
    created_at: str | None = None,
) -> None:

    application_version = _require_string(application_version, label="application_version")
    created_at = _normalise_created_at(created_at)
    payload, blobs = _archive_payload(project)
    document_ids = [document["id"] for document in payload["documents"]]
    sequence_members = {document["id"]: document["sequence_member"] for document in payload["documents"]}
    manifest = {
        "format": OGS_FORMAT,
        "formatVersion": OGS_FORMAT_VERSION,
        "applicationVersion": application_version,
        "createdAt": created_at,
        "hashAlgorithm": OGS_HASH_ALGORITHM,
        "projectId": str(project.id),
        "projectSchemaVersion": project.schema_version,
        "documentIds": document_ids,
        "sequenceMembers": sequence_members,
        "pluginRequirements": [],
        "databaseVersions": {},
    }
    manifest_bytes = _canonical_json(manifest, label="manifest")
    database_bytes = _build_sqlite(payload)
    members: dict[str, bytes] = {
        OGS_MANIFEST_MEMBER: manifest_bytes,
        OGS_DATABASE_MEMBER: database_bytes,
    }
    members.update(blobs)
    checksums = {
        "format": "opengene.checksums",
        "formatVersion": OGS_FORMAT_VERSION,
        "hashAlgorithm": OGS_HASH_ALGORITHM,
        "members": {name: _sha256(data) for name, data in sorted(members.items())},
    }
    checksum_bytes = _canonical_json(checksums, label="checksums")
    destination = Path(path)
    file_descriptor: int | None = None
    temporary_name: str | None = None
    try:
        file_descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
        )
        os.close(file_descriptor)
        file_descriptor = None
        with zipfile.ZipFile(
            temporary_name,
            mode="w",
            compression=zipfile.ZIP_DEFLATED,
            compresslevel=9,
            strict_timestamps=True,
        ) as archive:
            _write_member(archive, OGS_MANIFEST_MEMBER, members[OGS_MANIFEST_MEMBER])
            _write_member(archive, OGS_DATABASE_MEMBER, members[OGS_DATABASE_MEMBER])
            for name in sorted(blobs):
                _write_member(archive, name, blobs[name])
            _write_member(archive, OGS_CHECKSUM_MEMBER, checksum_bytes)
        os.replace(temporary_name, destination)
        temporary_name = None
    except (OSError, zipfile.BadZipFile) as error:
        raise OgsError(f"could not write .ogs archive {destination}: {error}") from error
    finally:
        if file_descriptor is not None:
            os.close(file_descriptor)
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass


def read_ogs(path: str | Path) -> Project:

    source = Path(path)
    try:
        with zipfile.ZipFile(source, mode="r") as archive:
            members = _archive_members(archive)
    except (OSError, zipfile.BadZipFile) as error:
        raise OgsError(f"could not open .ogs archive {source}: {error}") from error
    _validate_checksums(members)
    manifest = _json_object(members[OGS_MANIFEST_MEMBER], label=OGS_MANIFEST_MEMBER)
    _validate_manifest(manifest)
    payload = _read_sqlite(members[OGS_DATABASE_MEMBER])
    if payload.get("id") != manifest["projectId"]:
        raise OgsError("manifest projectId does not match project.sqlite")
    if payload.get("schema_version") != manifest.get("projectSchemaVersion"):
        raise OgsError("manifest projectSchemaVersion does not match project.sqlite")
    return _hydrate_project(payload, members, manifest)


__all__ = [
    "OGS_FORMAT",
    "OGS_FORMAT_VERSION",
    "OgsError",
    "read_ogs",
    "write_ogs",
]
