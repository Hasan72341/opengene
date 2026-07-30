from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from opengene_core import Feature, Location, MoleculeType, SequenceDocument, Topology
from opengene_core.models import thaw_value


class GenBankError(ValueError):
    pass


_LOCATION_RANGE = re.compile(r"^[<>]?(\d+)[<>]?\.\.[<>]?(\d+)[<>]?$")
_LOCATION_POINT = re.compile(r"^[<>]?(\d+)[<>]?$")
_LOCATION_BETWEEN = re.compile(r"^(\d+)\^(\d+)$")


def _split_top_level(value: str) -> list[str]:
    parts: list[str] = []
    depth = 0
    begin = 0
    for index, char in enumerate(value):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth < 0:
                raise GenBankError(f"unbalanced feature location {value!r}")
        elif char == "," and depth == 0:
            parts.append(value[begin:index])
            begin = index + 1
    if depth != 0:
        raise GenBankError(f"unbalanced feature location {value!r}")
    parts.append(value[begin:])
    return [part for part in parts if part]


def _outer_function(value: str) -> tuple[str, str] | None:
    match = re.match(r"^([A-Za-z_]+)\(", value)
    if not match or not value.endswith(")"):
        return None
    depth = 0
    for index, char in enumerate(value[match.end() - 1 :], start=match.end() - 1):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth == 0:
                if index != len(value) - 1:
                    return None
                return match.group(1).lower(), value[match.end() : -1]
    return None


def _parse_location(expression: str, strand: int = 1) -> list[Location]:
    value = re.sub(r"\s+", "", expression)
    if not value:
        raise GenBankError("feature location is empty")
    outer = _outer_function(value)
    if outer is not None:
        function, inner = outer
        if function == "complement":
            return _parse_location(inner, -strand)
        if function in {"join", "order"}:
            locations: list[Location] = []
            for part in _split_top_level(inner):
                locations.extend(_parse_location(part, strand))
            return locations
        raise GenBankError(f"unsupported GenBank location operator {function!r}")
    if ":" in value:
        raise GenBankError("remote GenBank feature locations are not supported")
    between = _LOCATION_BETWEEN.match(value)
    if between:
        left, right = (int(number) for number in between.groups())
        if right != left + 1:
            raise GenBankError(f"invalid between-base location {expression!r}")
        return [Location(left, left, strand)]
    range_match = _LOCATION_RANGE.match(value)
    if range_match:
        start, end = (int(number) for number in range_match.groups())
        if start < 1 or end < start:
            raise GenBankError(f"invalid GenBank interval {expression!r}")
        return [Location(start - 1, end, strand)]
    point_match = _LOCATION_POINT.match(value)
    if point_match:
        point = int(point_match.group(1))
        if point < 1:
            raise GenBankError(f"invalid GenBank point {expression!r}")
        return [Location(point - 1, point, strand)]
    raise GenBankError(f"unsupported GenBank location {expression!r}")


def _merge_origin_wrap(
    locations: list[Location], topology: Topology, sequence_length: int
) -> tuple[Location, ...]:

    if (
        topology is Topology.CIRCULAR
        and len(locations) == 2
        and not locations[0].wraps_origin
        and not locations[1].wraps_origin
        and locations[0].end == sequence_length
        and locations[1].start == 0
        and locations[0].start > locations[1].end
        and locations[0].strand == locations[1].strand
    ):
        return (
            Location(
                locations[0].start,
                locations[1].end,
                locations[0].strand,
                wraps_origin=True,
            ),
        )
    return tuple(locations)


def _append_qualifier(target: dict[str, Any], key: str, value: Any) -> None:
    if key not in target:
        target[key] = value
    elif isinstance(target[key], list):
        target[key].append(value)
    else:
        target[key] = [target[key], value]


def _parse_feature_lines(lines: list[str]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    pending_key: str | None = None
    pending_value: str = ""
    pending_quoted = False

    def finish_pending() -> None:
        nonlocal pending_key, pending_value, pending_quoted
        if current is not None and pending_key is not None:
            value = pending_value
            if pending_quoted and value.endswith('"'):
                value = value[:-1]
            _append_qualifier(current["qualifiers"], pending_key, value.replace('\\"', '"'))
        pending_key = None
        pending_value = ""
        pending_quoted = False

    for raw_line in lines:
        line = raw_line.rstrip("\n\r")
        key_column = line[5:21].strip() if len(line) >= 6 else ""
        if line.startswith("     ") and key_column:
            finish_pending()
            current = {"feature_type": key_column, "location": line[21:].strip(), "qualifiers": {}}
            records.append(current)
            continue
        if current is None:
            continue
        content = line[21:].strip() if len(line) > 21 else line.strip()
        if content.startswith("/"):
            finish_pending()
            qualifier = content[1:]
            if "=" not in qualifier:
                _append_qualifier(current["qualifiers"], qualifier, True)
                continue
            qualifier_key, qualifier_value = qualifier.split("=", 1)
            pending_key = qualifier_key
            pending_quoted = qualifier_value.startswith('"')
            pending_value = qualifier_value[1:] if pending_quoted else qualifier_value
            if not pending_quoted or pending_value.endswith('"'):
                finish_pending()
        elif pending_key is not None:
            pending_value += content
            if pending_quoted and pending_value.endswith('"'):
                finish_pending()
        else:
            current["location"] += content
    finish_pending()
    return records


def _header_values(lines: list[str]) -> dict[str, str]:
    values: dict[str, list[str]] = {}
    current_key: str | None = None
    for line in lines:
        key = line[:12].strip()
        value = line[12:].strip() if len(line) > 12 else ""
        if key:
            current_key = key
            values.setdefault(key, []).append(value)
        elif current_key is not None and value:
            values[current_key].append(value)
    return {key: " ".join(value).strip() for key, value in values.items()}


def _locus_metadata(locus: str) -> tuple[str, MoleculeType, Topology, int | None]:
    tokens = locus.split()
    if not tokens:
        raise GenBankError("LOCUS line has no name")
    name = tokens[0]
    declared_length: int | None = None
    if len(tokens) > 1 and tokens[1].isdigit():
        declared_length = int(tokens[1])
    lowered = {token.lower() for token in tokens}
    topology = Topology.CIRCULAR if "circular" in lowered else Topology.LINEAR
    if "rna" in lowered or "mrna" in lowered or "trna" in lowered or "rrna" in lowered:
        molecule_type = MoleculeType.RNA
    elif "aa" in lowered or "protein" in lowered:
        molecule_type = MoleculeType.PROTEIN
    else:
        molecule_type = MoleculeType.DNA
    return name, molecule_type, topology, declared_length


def _parse_record(record: str) -> SequenceDocument:
    lines = [line.rstrip("\r") for line in record.splitlines() if line.strip() != "//"]
    locus_index = next((index for index, line in enumerate(lines) if line.startswith("LOCUS")), None)
    origin_index = next((index for index, line in enumerate(lines) if line.startswith("ORIGIN")), None)
    if locus_index is None:
        raise GenBankError("record has no LOCUS line")
    if origin_index is None:
        raise GenBankError("record has no ORIGIN section")
    locus = lines[locus_index][12:].strip() if len(lines[locus_index]) > 12 else ""
    name, molecule_type, topology, declared_length = _locus_metadata(locus)
    features_index = next(
        (index for index, line in enumerate(lines) if line.startswith("FEATURES")), origin_index
    )
    header = _header_values(lines[locus_index + 1 : features_index])
    sequence = "".join(
        re.sub(r"[^A-Za-z*\-]", "", line)
        for line in lines[origin_index + 1 :]
    ).upper()
    if declared_length is not None and declared_length != len(sequence):
        header["declared_length"] = str(declared_length)

    features: list[Feature] = []
    feature_lines = lines[features_index + 1 : origin_index] if features_index < origin_index else []
    for raw_feature in _parse_feature_lines(feature_lines):
        locations = _merge_origin_wrap(
            _parse_location(raw_feature["location"]), topology, len(sequence)
        )
        qualifiers = raw_feature["qualifiers"]
        qualifiers["_genbank_location"] = raw_feature["location"]
        name_value: str | None = None
        for field in ("label", "gene", "standard_name", "product", "note"):
            candidate = qualifiers.get(field)
            if isinstance(candidate, list):
                candidate = candidate[0] if candidate else None
            if isinstance(candidate, str) and candidate.strip():
                name_value = candidate.strip()
                break
        features.append(
            Feature.create(
                raw_feature["feature_type"],
                name_value or raw_feature["feature_type"],
                locations,
                qualifiers,
            )
        )
    metadata: dict[str, Any] = {"source_format": "genbank"}
    metadata.update(
        {
            key.lower(): value
            for key, value in header.items()
            if key in {"DEFINITION", "ACCESSION", "VERSION", "KEYWORDS", "SOURCE", "COMMENT"}
        }
    )
    return SequenceDocument.create(
        name,
        sequence,
        molecule_type=molecule_type,
        topology=topology,
        features=features,
        metadata=metadata,
    )


def parse_genbank(source: str | bytes) -> list[SequenceDocument]:

    if isinstance(source, bytes):
        try:
            source = source.decode("utf-8")
        except UnicodeDecodeError as error:
            raise GenBankError("GenBank input must be UTF-8 text") from error
    if not isinstance(source, str):
        raise GenBankError("GenBank source must be text or bytes")
    records = [record for record in re.split(r"(?m)^//\s*$", source) if record.strip()]
    return [_parse_record(record) for record in records]


def read_genbank(path: str | Path) -> list[SequenceDocument]:
    return parse_genbank(Path(path).read_text(encoding="utf-8"))


def _format_interval(location: Location, sequence_length: int) -> str:
    if location.wraps_origin:
        return f"join({location.start + 1}..{sequence_length},1..{location.end})"
    if location.start == location.end:
        if location.start < sequence_length:
            return f"{location.start}^{location.start + 1}"
        return str(sequence_length)
    if location.end == location.start + 1:
        return str(location.start + 1)
    return f"{location.start + 1}..{location.end}"


def _format_feature_location(feature: Feature, sequence_length: int) -> str:
    tokens = [_format_interval(location, sequence_length) for location in feature.locations]
    strands = {location.strand for location in feature.locations}
    if len(tokens) == 1:
        result = tokens[0]
    else:
        result = f"join({','.join(tokens)})"
    if strands == {-1}:
        return f"complement({result})"
    if strands == {1}:
        return result
    mixed_tokens = [
        f"complement({token})" if location.strand == -1 else token
        for location, token in zip(feature.locations, tokens, strict=True)
    ]
    return f"join({','.join(mixed_tokens)})"


def _wrap_feature_value(prefix: str, value: str, *, width: int = 79) -> list[str]:
    room = max(1, width - len(prefix))
    chunks = [value[offset : offset + room] for offset in range(0, len(value), room)] or [""]
    return [prefix + chunks[0], *(" " * 21 + chunk for chunk in chunks[1:])]


def _format_qualifier(key: str, value: Any) -> str:
    if value is True:
        return f"/{key}"
    text = str(value).replace('"', '\\"')
    return f'/{key}="{text}"'


def _feature_output_lines(feature: Feature, sequence_length: int) -> list[str]:
    feature_type = feature.feature_type.replace("\n", " ")[:15]
    location = _format_feature_location(feature, sequence_length)
    lines = _wrap_feature_value(f"     {feature_type:<15} ", location)
    qualifiers = thaw_value(feature.qualifiers)
    for key, value in qualifiers.items():
        if str(key).startswith("_"):
            continue
        values = value if isinstance(value, list) else [value]
        for item in values:
            lines.extend(_wrap_feature_value(" " * 21, _format_qualifier(str(key), item)))
    return lines


def _header_line(key: str, value: str) -> list[str]:
    return _wrap_feature_value(f"{key:<12}", value)


def write_genbank(documents: SequenceDocument | Iterable[SequenceDocument]) -> str:

    if isinstance(documents, SequenceDocument):
        documents = (documents,)
    records: list[str] = []
    for document in documents:
        if not isinstance(document, SequenceDocument):
            raise GenBankError("GenBank export requires SequenceDocument values")
        molecule = "aa" if document.molecule_type is MoleculeType.PROTEIN else document.molecule_type.value
        topology = document.topology.value
        safe_name = re.sub(r"\s+", "_", document.name).strip("_") or "untitled"
        lines = [
            f"LOCUS       {safe_name:<16} {document.length:>7} bp    {molecule:<6} {topology:<8}  01-JAN-1980"
        ]
        metadata = thaw_value(document.metadata)
        definition = metadata.get("definition") or metadata.get("description")
        if definition:
            lines.extend(_header_line("DEFINITION", str(definition)))
        if metadata.get("accession"):
            lines.extend(_header_line("ACCESSION", str(metadata["accession"])))
        if metadata.get("version"):
            lines.extend(_header_line("VERSION", str(metadata["version"])))
        lines.append("FEATURES             Location/Qualifiers")
        for feature in document.features:
            lines.extend(_feature_output_lines(feature, document.length))
        lines.append("ORIGIN")
        for offset in range(0, document.length, 60):
            segment = document.sequence[offset : offset + 60].lower()
            groups = " ".join(segment[index : index + 10] for index in range(0, len(segment), 10))
            lines.append(f"{offset + 1:>9} {groups}")
        lines.append("//")
        records.append("\n".join(lines))
    return "\n".join(records) + ("\n" if records else "")


def write_genbank_file(
    path: str | Path, documents: SequenceDocument | Iterable[SequenceDocument]
) -> None:
    Path(path).write_text(write_genbank(documents), encoding="utf-8")
