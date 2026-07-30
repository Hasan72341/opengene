from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from enum import Enum
from hashlib import sha256
from types import MappingProxyType
from typing import Any, TypeAlias
from uuid import UUID, uuid4


class ModelValidationError(ValueError):
    pass


class MoleculeType(str, Enum):
    DNA = "DNA"
    RNA = "RNA"
    PROTEIN = "PROTEIN"


class Topology(str, Enum):
    LINEAR = "linear"
    CIRCULAR = "circular"


JsonScalar: TypeAlias = str | int | float | bool | None
JsonValue: TypeAlias = JsonScalar | tuple["JsonValue", ...] | Mapping[str, "JsonValue"]

_DNA_ALPHABET = frozenset("ACGTRYSWKMBDHVN-")
_RNA_ALPHABET = frozenset("ACGURYSWKMBDHVN-")
_PROTEIN_ALPHABET = frozenset("ABCDEFGHIJKLMNOPQRSTUVWXYZ*-")


def _freeze_value(value: Any) -> Any:

    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze_value(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_value(item) for item in value)
    if isinstance(value, set):
        return tuple(sorted((_freeze_value(item) for item in value), key=repr))
    return value


def thaw_value(value: Any) -> Any:

    if isinstance(value, Mapping):
        return {str(key): thaw_value(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [thaw_value(item) for item in value]
    return value


def _normalise_sequence(sequence: str, molecule_type: MoleculeType) -> str:
    if not isinstance(sequence, str):
        raise ModelValidationError("sequence must be a string")
    normalised = "".join(sequence.split()).upper()
    alphabet = {
        MoleculeType.DNA: _DNA_ALPHABET,
        MoleculeType.RNA: _RNA_ALPHABET,
        MoleculeType.PROTEIN: _PROTEIN_ALPHABET,
    }[molecule_type]
    invalid = sorted(set(normalised) - alphabet)
    if invalid:
        raise ModelValidationError(
            f"invalid {molecule_type.value} character(s): {''.join(invalid)}"
        )
    return normalised


@dataclass(frozen=True, slots=True)
class Location:

    start: int
    end: int
    strand: int = 1
    wraps_origin: bool = False

    def __post_init__(self) -> None:
        if isinstance(self.start, bool) or isinstance(self.end, bool):
            raise ModelValidationError("location coordinates must be integers")
        if not isinstance(self.start, int) or not isinstance(self.end, int):
            raise ModelValidationError("location coordinates must be integers")
        if self.start < 0 or self.end < 0:
            raise ModelValidationError("location coordinates cannot be negative")
        if self.strand not in (-1, 1):
            raise ModelValidationError("strand must be 1 or -1")
        if self.wraps_origin and self.start <= self.end:
            raise ModelValidationError(
                "an origin-spanning location requires start > end"
            )
        if not self.wraps_origin and self.start > self.end:
            raise ModelValidationError(
                "non-wrapping locations require start <= end"
            )

    def validate(self, sequence_length: int, topology: Topology) -> None:
        if sequence_length < 0:
            raise ModelValidationError("sequence length cannot be negative")
        if self.start > sequence_length or self.end > sequence_length:
            raise ModelValidationError(
                f"location [{self.start}, {self.end}) exceeds sequence length {sequence_length}"
            )
        if self.wraps_origin:
            if topology is not Topology.CIRCULAR:
                raise ModelValidationError("origin-spanning locations require circular topology")
            if sequence_length == 0 or self.start >= sequence_length or self.end == 0:
                raise ModelValidationError("invalid origin-spanning location for sequence length")

    def segments(self, sequence_length: int) -> tuple[tuple[int, int], ...]:

        self.validate(sequence_length, Topology.CIRCULAR if self.wraps_origin else Topology.LINEAR)
        if not self.wraps_origin:
            return ((self.start, self.end),)
        return ((self.start, sequence_length), (0, self.end))

    def length(self, sequence_length: int) -> int:
        return sum(end - start for start, end in self.segments(sequence_length))

    def contains(self, offset: int, sequence_length: int) -> bool:
        if offset < 0 or offset >= sequence_length:
            return False
        return any(start <= offset < end for start, end in self.segments(sequence_length))


@dataclass(frozen=True, slots=True)
class Feature:

    id: UUID
    feature_type: str
    name: str
    locations: tuple[Location, ...]
    qualifiers: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.id, UUID):
            raise ModelValidationError("feature id must be a UUID")
        if not self.feature_type.strip():
            raise ModelValidationError("feature_type cannot be empty")
        if not self.name.strip():
            raise ModelValidationError("feature name cannot be empty")
        locations = tuple(self.locations)
        if not locations:
            raise ModelValidationError("a feature requires at least one location")
        if not all(isinstance(location, Location) for location in locations):
            raise ModelValidationError("feature locations must be Location values")
        object.__setattr__(self, "locations", locations)
        object.__setattr__(self, "qualifiers", _freeze_value(self.qualifiers))

    @classmethod
    def create(
        cls,
        feature_type: str,
        name: str,
        locations: tuple[Location, ...] | list[Location],
        qualifiers: Mapping[str, Any] | None = None,
    ) -> "Feature":
        return cls(
            id=uuid4(),
            feature_type=feature_type,
            name=name,
            locations=tuple(locations),
            qualifiers=qualifiers or {},
        )

    def validate(self, sequence_length: int, topology: Topology) -> None:
        for location in self.locations:
            location.validate(sequence_length, topology)

    def with_updates(self, **changes: Any) -> "Feature":
        return replace(self, **changes)


@dataclass(frozen=True, slots=True)
class SequenceDocument:

    id: UUID
    name: str
    molecule_type: MoleculeType
    topology: Topology
    sequence: str
    features: tuple[Feature, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=dict)
    revision: int = 1

    def __post_init__(self) -> None:
        if not isinstance(self.id, UUID):
            raise ModelValidationError("document id must be a UUID")
        if not self.name.strip():
            raise ModelValidationError("document name cannot be empty")
        try:
            molecule_type = MoleculeType(self.molecule_type)
            topology = Topology(self.topology)
        except ValueError as error:
            raise ModelValidationError(str(error)) from error
        object.__setattr__(self, "molecule_type", molecule_type)
        object.__setattr__(self, "topology", topology)
        object.__setattr__(self, "sequence", _normalise_sequence(self.sequence, molecule_type))
        features = tuple(self.features)
        if not all(isinstance(feature, Feature) for feature in features):
            raise ModelValidationError("document features must be Feature values")
        feature_ids = [feature.id for feature in features]
        if len(set(feature_ids)) != len(feature_ids):
            raise ModelValidationError("feature ids must be unique within a document")
        for feature in features:
            feature.validate(len(self.sequence), topology)
        object.__setattr__(self, "features", features)
        object.__setattr__(self, "metadata", _freeze_value(self.metadata))
        if isinstance(self.revision, bool) or not isinstance(self.revision, int) or self.revision < 1:
            raise ModelValidationError("revision must be a positive integer")

    @classmethod
    def create(
        cls,
        name: str,
        sequence: str,
        *,
        molecule_type: MoleculeType = MoleculeType.DNA,
        topology: Topology = Topology.CIRCULAR,
        features: tuple[Feature, ...] | list[Feature] = (),
        metadata: Mapping[str, Any] | None = None,
    ) -> "SequenceDocument":
        return cls(
            id=uuid4(),
            name=name,
            molecule_type=molecule_type,
            topology=topology,
            sequence=sequence,
            features=tuple(features),
            metadata=metadata or {},
        )

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def sequence_hash(self) -> str:
        return sha256(self.sequence.encode("ascii")).hexdigest()

    def feature(self, feature_id: UUID) -> Feature:
        for feature in self.features:
            if feature.id == feature_id:
                return feature
        raise KeyError(f"feature {feature_id} is not present in {self.id}")

    def revised(self, **changes: Any) -> "SequenceDocument":

        changes.setdefault("revision", self.revision + 1)
        return replace(self, **changes)

    def with_updates(self, **changes: Any) -> "SequenceDocument":

        return replace(self, **changes)
