from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from hashlib import sha256
import json
from types import MappingProxyType
from typing import Any


_IUPAC: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("GC"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
_COMPLEMENT = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N",
}


class CloningValidationError(ValueError):
    pass


class Topology(str, Enum):

    LINEAR = "linear"
    CIRCULAR = "circular"


class EndKind(str, Enum):

    BLUNT = "blunt"
    FIVE_PRIME = "5-prime"
    THREE_PRIME = "3-prime"


def normalize_dna(sequence: str, *, allow_ambiguous: bool = True) -> str:

    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    normalised = "".join(sequence.split()).upper().replace("U", "T")
    alphabet = _IUPAC if allow_ambiguous else {base: _IUPAC[base] for base in "ACGT"}
    invalid = sorted(set(normalised).difference(alphabet))
    if invalid:
        raise CloningValidationError(
            f"sequence contains unsupported DNA symbols: {''.join(invalid)}"
        )
    return normalised


def reverse_complement(sequence: str) -> str:

    return "".join(_COMPLEMENT[base] for base in reversed(normalize_dna(sequence)))


def sequence_hash(sequence: str) -> str:

    return sha256(normalize_dna(sequence).encode("ascii")).hexdigest()


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze(item) for key, item in value.items()})
    if isinstance(value, (tuple, list)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, set):
        return tuple(sorted((_freeze(item) for item in value), key=repr))
    return value


def _json_value(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    raise CloningValidationError(
        f"provenance parameter contains non-JSON value {type(value).__name__}"
    )


def _stable_identifier(*parts: str) -> str:
    payload = "\x1f".join(parts).encode("utf-8")
    return sha256(payload).hexdigest()[:24]


@dataclass(frozen=True, slots=True)
class CohesiveEnd:

    kind: EndKind | str
    sequence: str = ""

    def __post_init__(self) -> None:
        try:
            kind = EndKind(self.kind)
        except ValueError as error:
            raise CloningValidationError("unknown cohesive-end kind") from error
        sequence = normalize_dna(self.sequence, allow_ambiguous=False)
        if kind is EndKind.BLUNT and sequence:
            raise CloningValidationError("a blunt end cannot contain an overhang sequence")
        if kind is not EndKind.BLUNT and not sequence:
            raise CloningValidationError("a sticky end requires an overhang sequence")
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "sequence", sequence)

    @property
    def is_blunt(self) -> bool:
        return self.kind is EndKind.BLUNT

    def compatible_with(self, other: "CohesiveEnd | None") -> bool:
        if other is None or not isinstance(other, CohesiveEnd):
            return False
        if self.kind is EndKind.BLUNT or other.kind is EndKind.BLUNT:
            return self.kind is other.kind is EndKind.BLUNT
        return self.kind is other.kind and self.sequence == reverse_complement(other.sequence)


@dataclass(frozen=True, slots=True)
class DNAFragment:

    name: str
    sequence: str
    topology: Topology | str = Topology.LINEAR
    id: str = ""
    left_end: CohesiveEnd | None = None
    right_end: CohesiveEnd | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise CloningValidationError("fragment name cannot be empty")
        sequence = normalize_dna(self.sequence)
        if not sequence:
            raise CloningValidationError("fragment sequence cannot be empty")
        try:
            topology = Topology(self.topology)
        except ValueError as error:
            raise CloningValidationError("topology must be 'linear' or 'circular'") from error
        if self.left_end is not None and not isinstance(self.left_end, CohesiveEnd):
            raise CloningValidationError("left_end must be a CohesiveEnd or None")
        if self.right_end is not None and not isinstance(self.right_end, CohesiveEnd):
            raise CloningValidationError("right_end must be a CohesiveEnd or None")
        if topology is Topology.CIRCULAR and (self.left_end is not None or self.right_end is not None):
            raise CloningValidationError("a circular fragment cannot have free DNA ends")
        identifier = self.id.strip() if isinstance(self.id, str) else ""
        if not identifier:
            identifier = _stable_identifier(self.name.strip(), sequence, topology.value)
        object.__setattr__(self, "name", self.name.strip())
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "topology", topology)
        object.__setattr__(self, "id", identifier)
        object.__setattr__(self, "metadata", _freeze(self.metadata))

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def sha256(self) -> str:
        return sequence_hash(self.sequence)

    @classmethod
    def from_document(cls, document: object, *, name: str | None = None) -> "DNAFragment":

        if isinstance(document, DNAFragment):
            return document
        sequence = getattr(document, "sequence", None)
        if not isinstance(sequence, str):
            raise TypeError("document-like input must expose a string sequence")
        raw_molecule_type = getattr(document, "molecule_type", None)
        molecule_type = getattr(raw_molecule_type, "value", raw_molecule_type)
        if molecule_type is not None and str(molecule_type).upper() != "DNA":
            raise CloningValidationError("cloning inputs must be DNA documents")
        raw_topology = getattr(document, "topology", Topology.LINEAR)
        topology = getattr(raw_topology, "value", raw_topology)
        raw_id = getattr(document, "id", "")
        return cls(
            name=name or str(getattr(document, "name", "sequence")),
            sequence=sequence,
            topology=topology,
            id=str(raw_id) if raw_id else "",
        )


def as_fragment(value: DNAFragment | str | object, *, name: str = "sequence") -> DNAFragment:

    if isinstance(value, DNAFragment):
        return value
    if isinstance(value, str):
        return DNAFragment(name=name, sequence=value)
    return DNAFragment.from_document(value, name=name)


@dataclass(frozen=True, slots=True)
class EnzymeDefinition:

    name: str
    recognition_site: str
    cut_forward: int
    cut_reverse: int

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise CloningValidationError("enzyme name cannot be empty")
        if not isinstance(self.cut_forward, int) or not isinstance(self.cut_reverse, int):
            raise CloningValidationError("enzyme cut offsets must be integers")
        site = normalize_dna(self.recognition_site)
        if not site:
            raise CloningValidationError("enzyme recognition_site cannot be empty")
        object.__setattr__(self, "name", self.name.strip())
        object.__setattr__(self, "recognition_site", site)

    @property
    def reverse_recognition_site(self) -> str:
        return reverse_complement(self.recognition_site)

    @property
    def end_kind(self) -> EndKind:
        if self.cut_forward == self.cut_reverse:
            return EndKind.BLUNT
        return EndKind.FIVE_PRIME if self.cut_forward < self.cut_reverse else EndKind.THREE_PRIME

    @classmethod
    def from_any(cls, enzyme: "EnzymeDefinition | object") -> "EnzymeDefinition":

        if isinstance(enzyme, cls):
            return enzyme
        try:
            return cls(
                name=str(getattr(enzyme, "name")),
                recognition_site=str(getattr(enzyme, "recognition_site")),
                cut_forward=getattr(enzyme, "cut_forward"),
                cut_reverse=getattr(enzyme, "cut_reverse"),
            )
        except (AttributeError, TypeError, ValueError) as error:
            raise TypeError(
                "enzyme must expose name, recognition_site, cut_forward, and cut_reverse"
            ) from error


@dataclass(frozen=True, slots=True)
class RestrictionCut:

    enzyme: EnzymeDefinition
    start: int
    end: int
    strand: int
    cut_forward: int
    cut_reverse: int
    left_end: CohesiveEnd
    right_end: CohesiveEnd
    wraps_origin: bool = False

    @property
    def end_geometry(self) -> CohesiveEnd:

        return self.left_end


@dataclass(frozen=True, slots=True)
class OperationRecord:

    operation_id: str
    operation_type: str
    input_ids: tuple[str, ...]
    input_sequence_hashes: Mapping[str, str]
    parameters: Mapping[str, Any]
    engine_name: str = "opengene-portable"
    engine_version: str = "1"
    output_id: str | None = None
    output_sequence_hash: str | None = None
    warnings: tuple[str, ...] = ()
    reproducibility_status: str = "reproducible"

    def __post_init__(self) -> None:
        if not self.operation_type.strip():
            raise CloningValidationError("operation_type cannot be empty")
        if not self.engine_name.strip() or not self.engine_version.strip():
            raise CloningValidationError("engine name and version cannot be empty")
        inputs = tuple(str(identifier) for identifier in self.input_ids)
        hashes = {str(key): str(value) for key, value in self.input_sequence_hashes.items()}
        parameters = _freeze(self.parameters)
        object.__setattr__(self, "input_ids", inputs)
        object.__setattr__(self, "input_sequence_hashes", MappingProxyType(hashes))
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "warnings", tuple(str(warning) for warning in self.warnings))
        expected_id = _operation_identifier(
            self.operation_type,
            inputs,
            hashes,
            parameters,
            self.engine_name,
            self.engine_version,
            self.output_id,
            self.output_sequence_hash,
        )
        if self.operation_id and self.operation_id != expected_id:
            raise CloningValidationError("operation_id does not match deterministic operation content")
        object.__setattr__(self, "operation_id", expected_id)

    def to_dict(self) -> dict[str, Any]:
        return {
            "operation_id": self.operation_id,
            "operation_type": self.operation_type,
            "input_ids": list(self.input_ids),
            "input_sequence_hashes": dict(self.input_sequence_hashes),
            "parameters": _json_value(self.parameters),
            "engine_name": self.engine_name,
            "engine_version": self.engine_version,
            "output_id": self.output_id,
            "output_sequence_hash": self.output_sequence_hash,
            "warnings": list(self.warnings),
            "reproducibility_status": self.reproducibility_status,
        }


def _operation_identifier(
    operation_type: str,
    input_ids: tuple[str, ...],
    input_sequence_hashes: Mapping[str, str],
    parameters: Mapping[str, Any],
    engine_name: str,
    engine_version: str,
    output_id: str | None,
    output_sequence_hash: str | None,
) -> str:
    payload = {
        "operation_type": operation_type,
        "input_ids": list(input_ids),
        "input_sequence_hashes": _json_value(input_sequence_hashes),
        "parameters": _json_value(parameters),
        "engine_name": engine_name,
        "engine_version": engine_version,
        "output_id": output_id,
        "output_sequence_hash": output_sequence_hash,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return f"op-{sha256(encoded).hexdigest()[:32]}"


def make_operation_record(
    operation_type: str,
    inputs: tuple[DNAFragment, ...] | list[DNAFragment],
    *,
    parameters: Mapping[str, Any] | None = None,
    output: DNAFragment | None = None,
    warnings: tuple[str, ...] | list[str] = (),
    reproducibility_status: str = "reproducible",
    engine_name: str = "opengene-portable",
    engine_version: str = "1",
) -> OperationRecord:

    frozen_inputs = tuple(inputs)
    if not all(isinstance(item, DNAFragment) for item in frozen_inputs):
        raise TypeError("inputs must contain DNAFragment values")
    input_ids = tuple(item.id for item in frozen_inputs)
    input_hashes = {item.id: item.sha256 for item in frozen_inputs}
    parameters = parameters or {}
    return OperationRecord(
        operation_id="",
        operation_type=operation_type,
        input_ids=input_ids,
        input_sequence_hashes=input_hashes,
        parameters=parameters,
        engine_name=engine_name,
        engine_version=engine_version,
        output_id=output.id if output else None,
        output_sequence_hash=output.sha256 if output else None,
        warnings=tuple(warnings),
        reproducibility_status=reproducibility_status,
    )


@dataclass(frozen=True, slots=True)
class RestrictionDigestPlan:

    input_fragment: DNAFragment
    enzymes: tuple[EnzymeDefinition, ...]
    cuts: tuple[RestrictionCut, ...]
    fragments: tuple[DNAFragment, ...]
    warnings: tuple[str, ...]
    provenance: OperationRecord


@dataclass(frozen=True, slots=True)
class LigationJunction:

    left_fragment_id: str
    right_fragment_id: str
    left_end: CohesiveEnd | None
    right_end: CohesiveEnd | None
    compatible: bool


@dataclass(frozen=True, slots=True)
class RestrictionLigationPlan:

    fragments: tuple[DNAFragment, ...]
    junctions: tuple[LigationJunction, ...]
    circular_product: bool
    product: DNAFragment | None
    warnings: tuple[str, ...]
    provenance: OperationRecord

    @property
    def feasible(self) -> bool:
        return self.product is not None


@dataclass(frozen=True, slots=True)
class OverlapJunction:

    left_fragment_id: str
    right_fragment_id: str
    sequence: str
    length: int
    compatible: bool


@dataclass(frozen=True, slots=True)
class OverlapAssemblyPlan:

    fragments: tuple[DNAFragment, ...]
    minimum_overlap: int
    circular_product: bool
    junctions: tuple[OverlapJunction, ...]
    product: DNAFragment | None
    warnings: tuple[str, ...]
    provenance: OperationRecord

    @property
    def feasible(self) -> bool:
        return self.product is not None


@dataclass(frozen=True, slots=True)
class GoldenGatePart:

    fragment: DNAFragment
    left_overhang: str
    right_overhang: str

    def __post_init__(self) -> None:
        if self.fragment.topology is not Topology.LINEAR:
            raise CloningValidationError("Golden Gate parts must be linear payload fragments")
        left = normalize_dna(self.left_overhang, allow_ambiguous=False)
        right = normalize_dna(self.right_overhang, allow_ambiguous=False)
        if not left or not right:
            raise CloningValidationError("Golden Gate overhangs cannot be empty")
        object.__setattr__(self, "left_overhang", left)
        object.__setattr__(self, "right_overhang", right)


@dataclass(frozen=True, slots=True)
class GoldenGateJunction:

    left_part_id: str
    right_part_id: str
    left_overhang: str
    right_overhang: str
    compatible: bool


@dataclass(frozen=True, slots=True)
class GoldenGatePlan:

    enzyme_name: str
    parts: tuple[GoldenGatePart, ...]
    circular_product: bool
    junctions: tuple[GoldenGateJunction, ...]
    product: DNAFragment | None
    warnings: tuple[str, ...]
    provenance: OperationRecord

    @property
    def feasible(self) -> bool:
        return self.product is not None


@dataclass(frozen=True, slots=True)
class OligoAlignment:

    top_start: int
    bottom_start: int
    overlap_length: int
    mismatches: int


@dataclass(frozen=True, slots=True)
class OligoAnnealingPlan:

    top_oligo: str
    bottom_oligo: str
    alignment: OligoAlignment | None
    product: DNAFragment | None
    warnings: tuple[str, ...]
    provenance: OperationRecord

    @property
    def feasible(self) -> bool:
        return self.product is not None
