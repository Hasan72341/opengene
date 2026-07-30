from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from threading import RLock
from typing import Any
from uuid import UUID

from .models import (
    Feature,
    Location,
    ModelValidationError,
    MoleculeType,
    SequenceDocument,
    Topology,
    thaw_value,
)


class CommandError(ValueError):
    pass


def _check_boundary(position: int, length: int, *, label: str = "position") -> None:
    if isinstance(position, bool) or not isinstance(position, int):
        raise CommandError(f"{label} must be an integer")
    if position < 0 or position > length:
        raise CommandError(f"{label} {position} is outside 0..{length}")


def _check_region(start: int, end: int, length: int) -> None:
    _check_boundary(start, length, label="start")
    _check_boundary(end, length, label="end")
    if start > end:
        raise CommandError("editing a region across the circular origin is ambiguous; rotate first")


def _clean_edit_bases(bases: str) -> str:
    if not isinstance(bases, str):
        raise CommandError("bases must be a string")
    return "".join(bases.split()).upper()


def _location_from_arc(
    start: int,
    end: int,
    arc_length: int,
    sequence_length: int,
    strand: int,
) -> Location | None:

    if arc_length <= 0 or sequence_length <= 0:
        return None
    if arc_length >= sequence_length:
        return Location(0, sequence_length, strand)
    if start < end:
        return Location(start, end, strand)
    if start > end:
        return Location(start, end, strand, wraps_origin=True)
    return None


def _insert_location(location: Location, position: int, old_length: int, amount: int) -> Location | None:

    new_length = old_length + amount
    if location.length(old_length) == 0:
        boundary = location.start + amount if position <= location.start else location.start
        return Location(boundary, boundary, location.strand)
    new_start = location.start + amount if position <= location.start else location.start
    new_end = location.end + amount if position < location.end else location.end
    if location.wraps_origin:
        inserted_inside = position < location.end or position > location.start
    else:
        inserted_inside = location.start < position < location.end
    arc_length = location.length(old_length) + (amount if inserted_inside else 0)
    return _location_from_arc(
        new_start, new_end, arc_length, new_length, location.strand
    )


def _delete_location(location: Location, start: int, end: int, old_length: int) -> Location | None:

    amount = end - start
    new_length = old_length - amount
    if new_length == 0:
        return None

    def map_retained_boundary(value: int) -> int:
        return value if value <= start else value - amount

    if location.length(old_length) == 0:
        boundary = location.start
        if start < boundary < end:
            boundary = start
        elif boundary >= end:
            boundary -= amount
        return Location(boundary, boundary, location.strand)

    retained: list[tuple[int, int]] = []
    for segment_start, segment_end in location.segments(old_length):
        left_end = min(segment_end, start)
        if segment_start < left_end:
            retained.append(
                (
                    map_retained_boundary(segment_start),
                    map_retained_boundary(left_end),
                )
            )
        right_start = max(segment_start, end)
        if right_start < segment_end:
            retained.append(
                (
                    map_retained_boundary(right_start),
                    map_retained_boundary(segment_end),
                )
            )
    if not retained:
        return None

    retained.sort()
    merged: list[list[int]] = []
    for segment_start, segment_end in retained:
        if merged and segment_start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], segment_end)
        else:
            merged.append([segment_start, segment_end])
    if len(merged) == 1:
        return Location(merged[0][0], merged[0][1], location.strand)
    if merged[0][0] == 0 and merged[-1][1] == new_length:
        return Location(
            merged[-1][0],
            merged[0][1],
            location.strand,
            wraps_origin=True,
        )
    return None


def _map_features_insert(
    document: SequenceDocument, position: int, amount: int
) -> tuple[Feature, ...]:
    mapped: list[Feature] = []
    for feature in document.features:
        locations = tuple(
            location
            for source in feature.locations
            if (location := _insert_location(source, position, document.length, amount))
            is not None
        )
        if locations:
            mapped.append(feature.with_updates(locations=locations))
    return tuple(mapped)


def _map_features_delete(
    document: SequenceDocument, start: int, end: int
) -> tuple[Feature, ...]:
    mapped: list[Feature] = []
    for feature in document.features:
        locations = tuple(
            location
            for source in feature.locations
            if (location := _delete_location(source, start, end, document.length))
            is not None
        )
        if locations:
            mapped.append(feature.with_updates(locations=locations))
    return tuple(mapped)


def _insert_document(document: SequenceDocument, position: int, bases: str) -> SequenceDocument:
    if not bases:
        return document
    sequence = document.sequence[:position] + bases + document.sequence[position:]
    return document.revised(
        sequence=sequence,
        features=_map_features_insert(document, position, len(bases)),
    )


def _delete_document(document: SequenceDocument, start: int, end: int) -> SequenceDocument:
    if start == end:
        return document
    sequence = document.sequence[:start] + document.sequence[end:]
    return document.revised(
        sequence=sequence,
        features=_map_features_delete(document, start, end),
    )


class EditCommand(ABC):

    @abstractmethod
    def apply(self, document: SequenceDocument) -> SequenceDocument:
        pass

    def inverse(self, before: SequenceDocument) -> "EditCommand":

        return RestoreDocument(before)


@dataclass(frozen=True, slots=True)
class RestoreDocument(EditCommand):

    document: SequenceDocument

    def apply(self, document: SequenceDocument) -> SequenceDocument:  # noqa: ARG002
        return self.document


@dataclass(frozen=True, slots=True)
class InsertBases(EditCommand):
    position: int
    bases: str

    def apply(self, document: SequenceDocument) -> SequenceDocument:
        _check_boundary(self.position, document.length)
        bases = _clean_edit_bases(self.bases)
        return _insert_document(document, self.position, bases)


@dataclass(frozen=True, slots=True)
class DeleteBases(EditCommand):
    start: int
    end: int

    def apply(self, document: SequenceDocument) -> SequenceDocument:
        _check_region(self.start, self.end, document.length)
        return _delete_document(document, self.start, self.end)


@dataclass(frozen=True, slots=True)
class ReplaceRegion(EditCommand):
    start: int
    end: int
    bases: str

    def apply(self, document: SequenceDocument) -> SequenceDocument:
        _check_region(self.start, self.end, document.length)
        bases = _clean_edit_bases(self.bases)
        deleted = _delete_document(document, self.start, self.end)
        if deleted is not document:
            deleted = deleted.with_updates(revision=document.revision)
        result = _insert_document(deleted, self.start, bases)
        if result is document:
            return document
        return result.with_updates(revision=document.revision + 1)


_DNA_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN-", "TGCAYRSWMKVHDBN-")
_RNA_COMPLEMENT = str.maketrans("ACGURYSWKMBDHVN-", "UGCAYRSWMKVHDBN-")


def _reverse_complement(sequence: str, molecule_type: MoleculeType) -> str:
    if molecule_type is MoleculeType.DNA:
        return sequence.translate(_DNA_COMPLEMENT)[::-1]
    if molecule_type is MoleculeType.RNA:
        return sequence.translate(_RNA_COMPLEMENT)[::-1]
    raise CommandError("reverse complement is only defined for DNA or RNA")


@dataclass(frozen=True, slots=True)
class ReverseComplement(EditCommand):

    def apply(self, document: SequenceDocument) -> SequenceDocument:
        sequence = _reverse_complement(document.sequence, document.molecule_type)
        length = document.length
        features: list[Feature] = []
        for feature in document.features:
            locations: list[Location] = []
            for location in reversed(feature.locations):
                if location.wraps_origin:
                    locations.append(
                        Location(
                            length - location.end,
                            length - location.start,
                            -location.strand,
                            wraps_origin=True,
                        )
                    )
                else:
                    locations.append(
                        Location(
                            length - location.end,
                            length - location.start,
                            -location.strand,
                        )
                    )
            features.append(feature.with_updates(locations=tuple(locations)))
        return document.revised(sequence=sequence, features=tuple(features))


@dataclass(frozen=True, slots=True)
class RotateOrigin(EditCommand):

    origin: int

    def apply(self, document: SequenceDocument) -> SequenceDocument:
        if document.topology is not Topology.CIRCULAR:
            raise CommandError("only circular documents can be rotated")
        if document.length == 0:
            raise CommandError("an empty document has no origin to rotate")
        _check_boundary(self.origin, document.length)
        if self.origin == document.length:
            raise CommandError("origin must identify a base, not the terminal boundary")
        if self.origin == 0:
            return document
        length = document.length
        features: list[Feature] = []
        for feature in document.features:
            locations: list[Location] = []
            for location in feature.locations:
                arc_length = location.length(length)
                start = (location.start - self.origin) % length
                end = (location.end - self.origin) % length
                mapped = _location_from_arc(start, end, arc_length, length, location.strand)
                if mapped is not None:
                    locations.append(mapped)
            if locations:
                features.append(feature.with_updates(locations=tuple(locations)))
        return document.revised(
            sequence=document.sequence[self.origin :] + document.sequence[: self.origin],
            features=tuple(features),
        )


@dataclass(frozen=True, slots=True)
class ChangeTopology(EditCommand):
    topology: Topology

    def apply(self, document: SequenceDocument) -> SequenceDocument:
        topology = Topology(self.topology)
        if topology is document.topology:
            return document
        if topology is Topology.LINEAR and any(
            location.wraps_origin
            for feature in document.features
            for location in feature.locations
        ):
            raise CommandError(
                "rotate or split origin-spanning features before converting to linear"
            )
        return document.revised(topology=topology)


@dataclass(frozen=True, slots=True)
class AddFeature(EditCommand):
    feature: Feature

    def apply(self, document: SequenceDocument) -> SequenceDocument:
        if any(existing.id == self.feature.id for existing in document.features):
            raise CommandError(f"feature id {self.feature.id} already exists")
        try:
            self.feature.validate(document.length, document.topology)
        except ModelValidationError as error:
            raise CommandError(str(error)) from error
        return document.revised(features=(*document.features, self.feature))


@dataclass(frozen=True, slots=True)
class RemoveFeature(EditCommand):
    feature_id: UUID

    def apply(self, document: SequenceDocument) -> SequenceDocument:
        remaining = tuple(
            feature for feature in document.features if feature.id != self.feature_id
        )
        if len(remaining) == len(document.features):
            raise CommandError(f"feature {self.feature_id} does not exist")
        return document.revised(features=remaining)


@dataclass(frozen=True, slots=True)
class UpdateFeature(EditCommand):

    feature_id: UUID
    replacement: Feature

    def apply(self, document: SequenceDocument) -> SequenceDocument:
        if self.replacement.id != self.feature_id:
            raise CommandError("feature replacement must preserve the feature id")
        try:
            self.replacement.validate(document.length, document.topology)
        except ModelValidationError as error:
            raise CommandError(str(error)) from error
        found = False
        features: list[Feature] = []
        for feature in document.features:
            if feature.id == self.feature_id:
                features.append(self.replacement)
                found = True
            else:
                features.append(feature)
        if not found:
            raise CommandError(f"feature {self.feature_id} does not exist")
        return document.revised(features=tuple(features))


_STANDARD_GENETIC_CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass(frozen=True, slots=True)
class TranslateSequence(EditCommand):

    frame: int = 0
    stop_at_stop: bool = False

    def apply(self, document: SequenceDocument) -> SequenceDocument:
        if document.molecule_type not in (MoleculeType.DNA, MoleculeType.RNA):
            raise CommandError("only DNA or RNA can be translated")
        if self.frame not in (0, 1, 2):
            raise CommandError("translation frame must be 0, 1, or 2")
        template = document.sequence.replace("U", "T")
        amino_acids: list[str] = []
        for offset in range(self.frame, len(template) - 2, 3):
            amino_acid = _STANDARD_GENETIC_CODE.get(template[offset : offset + 3], "X")
            if amino_acid == "*" and self.stop_at_stop:
                break
            amino_acids.append(amino_acid)
        metadata = thaw_value(document.metadata)
        metadata["translated_from_hash"] = document.sequence_hash
        metadata["translation_frame"] = self.frame
        return document.revised(
            name=f"{document.name} translation",
            molecule_type=MoleculeType.PROTEIN,
            topology=Topology.LINEAR,
            sequence="".join(amino_acids),
            features=(),
            metadata=metadata,
        )


class CommandHistory:

    def __init__(self, document: SequenceDocument) -> None:
        self._document = document
        self._undo: list[tuple[EditCommand, SequenceDocument]] = []
        self._redo: list[tuple[EditCommand, SequenceDocument]] = []
        self._lock = RLock()

    @property
    def document(self) -> SequenceDocument:
        return self._document

    @property
    def can_undo(self) -> bool:
        return bool(self._undo)

    @property
    def can_redo(self) -> bool:
        return bool(self._redo)

    @property
    def undo_depth(self) -> int:
        return len(self._undo)

    @property
    def redo_depth(self) -> int:
        return len(self._redo)

    def execute(self, command: EditCommand) -> SequenceDocument:
        if not isinstance(command, EditCommand):
            raise TypeError("command must implement EditCommand")
        with self._lock:
            before = self._document
            after = command.apply(before)
            if after == before:
                return after
            self._undo.append((command, before))
            self._redo.clear()
            self._document = after
            return after

    def undo(self) -> SequenceDocument:
        with self._lock:
            if not self._undo:
                raise CommandError("there is no command to undo")
            command, before = self._undo.pop()
            self._redo.append((command, self._document))
            self._document = before
            return before

    def redo(self) -> SequenceDocument:
        with self._lock:
            if not self._redo:
                raise CommandError("there is no command to redo")
            command, after = self._redo.pop()
            self._undo.append((command, self._document))
            self._document = after
            return after

    def clear(self) -> None:
        with self._lock:
            self._undo.clear()
            self._redo.clear()
