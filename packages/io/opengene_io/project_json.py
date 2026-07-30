from __future__ import annotations

from pathlib import Path

from opengene_core import Project


def read_project_json(path: str | Path) -> Project:
    return Project.from_json(Path(path).read_text(encoding="utf-8"))


def write_project_json(path: str | Path, project: Project, *, indent: int | None = 2) -> None:
    Path(path).write_text(project.to_json(indent=indent), encoding="utf-8")
