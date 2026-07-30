from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from opengene_core import MoleculeType, SequenceDocument, Topology


class FastaError(ValueError):
    pass


def parse_fasta(
    source: str | bytes,
    *,
    molecule_type: MoleculeType = MoleculeType.DNA,
    topology: Topology = Topology.LINEAR,
) -> list[SequenceDocument]:

    if isinstance(source, bytes):
        try:
            source = source.decode("utf-8")
        except UnicodeDecodeError as error:
            raise FastaError("FASTA must be UTF-8 text") from error
    if not isinstance(source, str):
        raise FastaError("FASTA source must be text or bytes")

    records: list[SequenceDocument] = []
    header: str | None = None
    sequence_lines: list[str] = []

    def finish_record() -> None:
        nonlocal header, sequence_lines
        if header is None:
            return
        parts = header.split(maxsplit=1)
        name = parts[0] if parts else "untitled"
        metadata: dict[str, str] = {"source_format": "fasta"}
        if len(parts) == 2:
            metadata["description"] = parts[1]
        try:
            records.append(
                SequenceDocument.create(
                    name,
                    "".join(sequence_lines),
                    molecule_type=molecule_type,
                    topology=topology,
                    metadata=metadata,
                )
            )
        except ValueError as error:
            raise FastaError(f"invalid FASTA record {name!r}: {error}") from error
        header = None
        sequence_lines = []

    for line_number, raw_line in enumerate(source.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            finish_record()
            header = line[1:].strip()
            if not header:
                raise FastaError(f"FASTA header on line {line_number} is empty")
            continue
        if header is None:
            raise FastaError(f"sequence appears before a FASTA header on line {line_number}")
        sequence_lines.append(line)
    finish_record()
    return records


def read_fasta(
    path: str | Path,
    *,
    molecule_type: MoleculeType = MoleculeType.DNA,
    topology: Topology = Topology.LINEAR,
) -> list[SequenceDocument]:
    return parse_fasta(
        Path(path).read_text(encoding="utf-8"),
        molecule_type=molecule_type,
        topology=topology,
    )


def write_fasta(
    documents: SequenceDocument | Iterable[SequenceDocument], *, line_width: int = 80
) -> str:

    if line_width < 1:
        raise FastaError("line_width must be positive")
    if isinstance(documents, SequenceDocument):
        documents = (documents,)
    lines: list[str] = []
    for document in documents:
        if not isinstance(document, SequenceDocument):
            raise FastaError("FASTA export requires SequenceDocument values")
        name = document.name.replace("\n", " ").replace("\r", " ").strip()
        description = document.metadata.get("description")
        header = name
        if isinstance(description, str) and description.strip():
            header = f"{header} {description.replace(chr(10), ' ').replace(chr(13), ' ').strip()}"
        lines.append(f">{header}")
        if document.sequence:
            lines.extend(
                document.sequence[offset : offset + line_width]
                for offset in range(0, document.length, line_width)
            )
        else:
            lines.append("")
    return "\n".join(lines) + ("\n" if lines else "")


def write_fasta_file(
    path: str | Path,
    documents: SequenceDocument | Iterable[SequenceDocument],
    *,
    line_width: int = 80,
) -> None:
    Path(path).write_text(write_fasta(documents, line_width=line_width), encoding="utf-8")
