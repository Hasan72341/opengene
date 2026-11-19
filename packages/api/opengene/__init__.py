from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import Iterable
from uuid import UUID

from opengene_core import Project as CoreProject
from opengene_core import ProvenanceNode, SequenceDocument
from opengene_io import read_ogs, read_project_json, write_ogs, write_project_json


class Project:

    def __init__(self, project: CoreProject, *, path: str | Path | None = None) -> None:
        if not isinstance(project, CoreProject):
            raise TypeError("project must be an opengene_core.Project")
        self._project = project
        self._path = Path(path) if path is not None else None

    @classmethod
    def create(
        cls,
        name: str,
        documents: Iterable[SequenceDocument] = (),
        *,
        metadata: dict[str, object] | None = None,
    ) -> "Project":
        return cls(CoreProject.create(name, tuple(documents), metadata=metadata))

    @classmethod
    def open(cls, path: str | Path) -> "Project":

        source = Path(path)
        if source.suffix.lower() == ".ogs":
            project = read_ogs(source)
        else:
            project = read_project_json(source)
        return cls(project, path=source)

    @property
    def core(self) -> CoreProject:

        return self._project

    @property
    def documents(self) -> tuple[SequenceDocument, ...]:
        return self._project.documents

    @property
    def provenance(self) -> tuple[ProvenanceNode, ...]:
        return self._project.provenance

    def sequence(self, identifier: UUID | str) -> SequenceDocument:

        return self._project.document(identifier)

    def add(self, document: SequenceDocument, *, name: str | None = None) -> SequenceDocument:

        if not isinstance(document, SequenceDocument):
            raise TypeError("document must be an opengene_core.SequenceDocument")
        if name is not None:
            document = document.with_updates(name=name)
        self._project = self._project.with_document(document)
        metadata = document.metadata
        recipe = metadata.get("opengene.provenance") if hasattr(metadata, "get") else None
        if isinstance(recipe, Mapping):
            raw_inputs = recipe.get("input_document_ids", ())
            try:
                input_ids = tuple(UUID(str(identifier)) for identifier in raw_inputs)
            except (TypeError, ValueError) as error:
                raise ValueError("automation provenance has invalid input document IDs") from error
            known_ids = {item.id for item in self._project.documents}
            if set(input_ids).issubset(known_ids):
                node = ProvenanceNode.create(
                    str(recipe.get("operation_type", "automation")),
                    input_ids,
                    document.id,
                    parameters=dict(recipe.get("parameters", {})),
                    software_version=str(recipe.get("software_version", "opengene automation API")),
                    engine_version=str(recipe.get("engine_version", "unknown")),
                    warnings=tuple(str(warning) for warning in recipe.get("warnings", ())),
                    reproducibility_status=str(recipe.get("reproducibility_status", "reproducible")),
                )
                self._project = self._project.with_provenance(node)
        return document

    def save(self, path: str | Path | None = None) -> Path:

        target = Path(path) if path is not None else self._path
        if target is None:
            raise ValueError("provide a path the first time this project is saved")
        if target.suffix.lower() == ".ogs":
            write_ogs(target, self._project)
        else:
            write_project_json(target, self._project)
        self._path = target
        return target


# Miscellaneous and automation submodules
from . import misc
from .misc import assgn_two, target_gene, AssignmentResult

__all__ = ["Project", "SequenceDocument", "misc", "assgn_two", "target_gene", "AssignmentResult"]
