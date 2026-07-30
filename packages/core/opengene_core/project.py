from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from enum import Enum
from typing import Any
from uuid import UUID, uuid4

from .models import (
    Feature,
    Location,
    ModelValidationError,
    MoleculeType,
    SequenceDocument,
    Topology,
    _freeze_value,
    thaw_value,
)


PROJECT_SCHEMA_VERSION = 1


class ProjectValidationError(ValueError):
    pass


def _validate_provenance_dag(nodes: tuple["ProvenanceNode", ...]) -> None:

    edges: dict[UUID, set[UUID]] = {}
    for node in nodes:
        for input_id in node.input_document_ids:
            if input_id == node.output_document_id:
                raise ProjectValidationError(
                    f"provenance {node.operation_id} cannot consume and produce the same document"
                )
            edges.setdefault(input_id, set()).add(node.output_document_id)

    visiting: set[UUID] = set()
    visited: set[UUID] = set()

    def visit(document_id: UUID) -> None:
        if document_id in visiting:
            raise ProjectValidationError("provenance graph contains a cycle")
        if document_id in visited:
            return
        visiting.add(document_id)
        for output_id in edges.get(document_id, ()):
            visit(output_id)
        visiting.remove(document_id)
        visited.add(document_id)

    for document_id in tuple(edges):
        visit(document_id)


def _json_safe(value: Any, *, path: str = "value") -> Any:

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, UUID):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item, path=f"{path}.{key}") for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_safe(item, path=f"{path}[{index}]") for index, item in enumerate(value)]
    raise ProjectValidationError(
        f"{path} contains {type(value).__name__}, which is not JSON serializable"
    )


def location_to_dict(location: Location) -> dict[str, Any]:
    return {
        "start": location.start,
        "end": location.end,
        "strand": location.strand,
        "wraps_origin": location.wraps_origin,
    }


def location_from_dict(payload: Mapping[str, Any]) -> Location:
    try:
        return Location(
            start=payload["start"],
            end=payload["end"],
            strand=payload.get("strand", 1),
            wraps_origin=payload.get("wraps_origin", False),
        )
    except (KeyError, TypeError, ModelValidationError) as error:
        raise ProjectValidationError(f"invalid location: {error}") from error


def feature_to_dict(feature: Feature) -> dict[str, Any]:
    return {
        "id": str(feature.id),
        "feature_type": feature.feature_type,
        "name": feature.name,
        "locations": [location_to_dict(location) for location in feature.locations],
        "qualifiers": _json_safe(thaw_value(feature.qualifiers), path="feature.qualifiers"),
    }


def feature_from_dict(payload: Mapping[str, Any]) -> Feature:
    try:
        return Feature(
            id=UUID(str(payload["id"])),
            feature_type=str(payload["feature_type"]),
            name=str(payload["name"]),
            locations=tuple(location_from_dict(item) for item in payload["locations"]),
            qualifiers=payload.get("qualifiers", {}),
        )
    except (KeyError, TypeError, ValueError, ModelValidationError) as error:
        raise ProjectValidationError(f"invalid feature: {error}") from error


def document_to_dict(document: SequenceDocument) -> dict[str, Any]:
    return {
        "id": str(document.id),
        "name": document.name,
        "molecule_type": document.molecule_type.value,
        "topology": document.topology.value,
        "sequence": document.sequence,
        "features": [feature_to_dict(feature) for feature in document.features],
        "metadata": _json_safe(thaw_value(document.metadata), path="document.metadata"),
        "revision": document.revision,
        "sequence_sha256": document.sequence_hash,
    }


def document_from_dict(payload: Mapping[str, Any]) -> SequenceDocument:
    try:
        document = SequenceDocument(
            id=UUID(str(payload["id"])),
            name=str(payload["name"]),
            molecule_type=MoleculeType(payload["molecule_type"]),
            topology=Topology(payload["topology"]),
            sequence=str(payload["sequence"]),
            features=tuple(feature_from_dict(item) for item in payload.get("features", [])),
            metadata=payload.get("metadata", {}),
            revision=payload.get("revision", 1),
        )
    except (KeyError, TypeError, ValueError, ModelValidationError) as error:
        raise ProjectValidationError(f"invalid sequence document: {error}") from error
    supplied_hash = payload.get("sequence_sha256")
    if supplied_hash is not None and supplied_hash != document.sequence_hash:
        raise ProjectValidationError(
            f"sequence_sha256 does not match sequence for document {document.id}"
        )
    return document


@dataclass(frozen=True, slots=True)
class ProvenanceNode:

    operation_id: UUID
    operation_type: str
    input_document_ids: tuple[UUID, ...]
    output_document_id: UUID
    parameters: Mapping[str, Any] = field(default_factory=dict)
    software_version: str = "unknown"
    warnings: tuple[str, ...] = ()
    engine_version: str | None = None
    reproducibility_status: str = "reproducible"

    def __post_init__(self) -> None:
        if not isinstance(self.operation_id, UUID) or not isinstance(self.output_document_id, UUID):
            raise ProjectValidationError("provenance operation and output ids must be UUIDs")
        inputs = tuple(self.input_document_ids)
        if not all(isinstance(identifier, UUID) for identifier in inputs):
            raise ProjectValidationError("provenance input ids must be UUIDs")
        if not self.operation_type.strip():
            raise ProjectValidationError("provenance operation_type cannot be empty")
        if not self.software_version.strip():
            raise ProjectValidationError("provenance software_version cannot be empty")
        object.__setattr__(self, "input_document_ids", inputs)
        object.__setattr__(self, "parameters", _freeze_value(self.parameters))
        object.__setattr__(self, "warnings", tuple(str(warning) for warning in self.warnings))

    @classmethod
    def create(
        cls,
        operation_type: str,
        input_document_ids: tuple[UUID, ...] | list[UUID],
        output_document_id: UUID,
        *,
        parameters: Mapping[str, Any] | None = None,
        software_version: str = "unknown",
        warnings: tuple[str, ...] | list[str] = (),
        engine_version: str | None = None,
        reproducibility_status: str = "reproducible",
    ) -> "ProvenanceNode":
        return cls(
            operation_id=uuid4(),
            operation_type=operation_type,
            input_document_ids=tuple(input_document_ids),
            output_document_id=output_document_id,
            parameters=parameters or {},
            software_version=software_version,
            warnings=tuple(warnings),
            engine_version=engine_version,
            reproducibility_status=reproducibility_status,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "operation_id": str(self.operation_id),
            "operation_type": self.operation_type,
            "input_document_ids": [str(identifier) for identifier in self.input_document_ids],
            "output_document_id": str(self.output_document_id),
            "parameters": _json_safe(thaw_value(self.parameters), path="provenance.parameters"),
            "software_version": self.software_version,
            "engine_version": self.engine_version,
            "warnings": list(self.warnings),
            "reproducibility_status": self.reproducibility_status,
        }

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "ProvenanceNode":
        try:
            return cls(
                operation_id=UUID(str(payload["operation_id"])),
                operation_type=str(payload["operation_type"]),
                input_document_ids=tuple(
                    UUID(str(identifier)) for identifier in payload.get("input_document_ids", [])
                ),
                output_document_id=UUID(str(payload["output_document_id"])),
                parameters=payload.get("parameters", {}),
                software_version=str(payload.get("software_version", "unknown")),
                engine_version=payload.get("engine_version"),
                warnings=tuple(payload.get("warnings", [])),
                reproducibility_status=str(
                    payload.get("reproducibility_status", "reproducible")
                ),
            )
        except (KeyError, TypeError, ValueError, ProjectValidationError) as error:
            raise ProjectValidationError(f"invalid provenance node: {error}") from error


@dataclass(frozen=True, slots=True)
class Project:

    id: UUID
    name: str
    documents: tuple[SequenceDocument, ...] = ()
    provenance: tuple[ProvenanceNode, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)
    schema_version: int = PROJECT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if not isinstance(self.id, UUID):
            raise ProjectValidationError("project id must be a UUID")
        if not self.name.strip():
            raise ProjectValidationError("project name cannot be empty")
        documents = tuple(self.documents)
        if not all(isinstance(document, SequenceDocument) for document in documents):
            raise ProjectValidationError("project documents must be SequenceDocument values")
        identifiers = [document.id for document in documents]
        if len(set(identifiers)) != len(identifiers):
            raise ProjectValidationError("document ids must be unique within a project")
        provenance = tuple(self.provenance)
        if not all(isinstance(node, ProvenanceNode) for node in provenance):
            raise ProjectValidationError("project provenance must contain ProvenanceNode values")
        operation_ids = [node.operation_id for node in provenance]
        if len(set(operation_ids)) != len(operation_ids):
            raise ProjectValidationError("provenance operation ids must be unique")
        document_ids = set(identifiers)
        for node in provenance:
            referenced = set(node.input_document_ids) | {node.output_document_id}
            missing = referenced - document_ids
            if missing:
                raise ProjectValidationError(
                    f"provenance {node.operation_id} references missing documents: {sorted(map(str, missing))}"
                )
        _validate_provenance_dag(provenance)
        if self.schema_version != PROJECT_SCHEMA_VERSION:
            raise ProjectValidationError(
                f"unsupported project schema {self.schema_version}; expected {PROJECT_SCHEMA_VERSION}"
            )
        object.__setattr__(self, "documents", documents)
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "metadata", _freeze_value(self.metadata))

    @classmethod
    def create(
        cls,
        name: str,
        documents: tuple[SequenceDocument, ...] | list[SequenceDocument] = (),
        *,
        metadata: Mapping[str, Any] | None = None,
    ) -> "Project":
        return cls(id=uuid4(), name=name, documents=tuple(documents), metadata=metadata or {})

    def document(self, identifier: UUID | str) -> SequenceDocument:
        if isinstance(identifier, UUID):
            for document in self.documents:
                if document.id == identifier:
                    return document
            raise KeyError(f"document {identifier} does not exist")
        matches = [document for document in self.documents if document.name == identifier]
        if not matches:
            raise KeyError(f"document named {identifier!r} does not exist")
        if len(matches) > 1:
            raise KeyError(f"document name {identifier!r} is ambiguous")
        return matches[0]

    def with_document(self, document: SequenceDocument) -> "Project":

        found = False
        documents: list[SequenceDocument] = []
        for existing in self.documents:
            if existing.id == document.id:
                documents.append(document)
                found = True
            else:
                documents.append(existing)
        if not found:
            documents.append(document)
        return replace(self, documents=tuple(documents))

    def without_document(self, document_id: UUID) -> "Project":
        if any(
            document_id in node.input_document_ids or document_id == node.output_document_id
            for node in self.provenance
        ):
            raise ProjectValidationError(
                "cannot remove a document referenced by provenance; remove/rewrite history first"
            )
        documents = tuple(document for document in self.documents if document.id != document_id)
        if len(documents) == len(self.documents):
            raise KeyError(f"document {document_id} does not exist")
        return replace(self, documents=documents)

    def with_provenance(self, node: ProvenanceNode) -> "Project":
        if any(existing.operation_id == node.operation_id for existing in self.provenance):
            raise ProjectValidationError(f"operation {node.operation_id} already exists")
        return replace(self, provenance=(*self.provenance, node))

    def to_dict(self) -> dict[str, Any]:
        return {
            "format": "opengene.project",
            "schema_version": self.schema_version,
            "id": str(self.id),
            "name": self.name,
            "metadata": _json_safe(thaw_value(self.metadata), path="project.metadata"),
            "documents": [document_to_dict(document) for document in self.documents],
            "provenance": [node.to_dict() for node in self.provenance],
        }

    def to_json(self, *, indent: int | None = 2) -> str:
        return json.dumps(
            self.to_dict(),
            ensure_ascii=False,
            indent=indent,
            sort_keys=True,
            separators=None if indent is not None else (",", ":"),
        ) + "\n"

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> "Project":
        if payload.get("format") != "opengene.project":
            raise ProjectValidationError("not an opengene.project document")
        try:
            return cls(
                id=UUID(str(payload["id"])),
                name=str(payload["name"]),
                documents=tuple(document_from_dict(item) for item in payload.get("documents", [])),
                provenance=tuple(
                    ProvenanceNode.from_dict(item) for item in payload.get("provenance", [])
                ),
                metadata=payload.get("metadata", {}),
                schema_version=payload.get("schema_version", PROJECT_SCHEMA_VERSION),
            )
        except (KeyError, TypeError, ValueError, ProjectValidationError) as error:
            if isinstance(error, ProjectValidationError):
                raise
            raise ProjectValidationError(f"invalid project: {error}") from error

    @classmethod
    def from_json(cls, source: str | bytes | bytearray) -> "Project":
        try:
            payload = json.loads(source)
        except (TypeError, json.JSONDecodeError) as error:
            raise ProjectValidationError(f"invalid project JSON: {error}") from error
        if not isinstance(payload, Mapping):
            raise ProjectValidationError("project JSON root must be an object")
        return cls.from_dict(payload)


def project_to_json(project: Project, *, indent: int | None = 2) -> str:
    return project.to_json(indent=indent)


def project_from_json(source: str | bytes | bytearray) -> Project:
    return Project.from_json(source)
