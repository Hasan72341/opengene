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
    "A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S",
    "W": "W", "K": "M", "M": "K", "B": "V", "D": "H", "H": "D", "V": "B", "N": "N",
}


class PCRValidationError(ValueError):
    pass


class Topology(str, Enum):

    LINEAR = "linear"
    CIRCULAR = "circular"


def normalize_dna(sequence: str, *, allow_ambiguous: bool = True) -> str:

    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    normalised = "".join(sequence.split()).upper().replace("U", "T")
    alphabet = _IUPAC if allow_ambiguous else {base: _IUPAC[base] for base in "ACGT"}
    invalid = sorted(set(normalised).difference(alphabet))
    if invalid:
        raise PCRValidationError(
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
    raise PCRValidationError(f"provenance contains non-JSON value {type(value).__name__}")


def _stable_identifier(*parts: str) -> str:
    return sha256("\x1f".join(parts).encode("utf-8")).hexdigest()[:24]


@dataclass(frozen=True, slots=True)
class PCRTemplate:

    name: str
    sequence: str
    topology: Topology | str = Topology.LINEAR
    id: str = ""
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise PCRValidationError("template name cannot be empty")
        sequence = normalize_dna(self.sequence)
        if not sequence:
            raise PCRValidationError("template sequence cannot be empty")
        try:
            topology = Topology(self.topology)
        except ValueError as error:
            raise PCRValidationError("topology must be 'linear' or 'circular'") from error
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
    def from_document(cls, document: object, *, name: str | None = None) -> "PCRTemplate":

        if isinstance(document, cls):
            return document
        sequence = getattr(document, "sequence", None)
        if not isinstance(sequence, str):
            raise TypeError("document-like input must expose a string sequence")
        raw_molecule_type = getattr(document, "molecule_type", None)
        molecule_type = getattr(raw_molecule_type, "value", raw_molecule_type)
        if molecule_type is not None and str(molecule_type).upper() != "DNA":
            raise PCRValidationError("PCR templates must be DNA documents")
        raw_topology = getattr(document, "topology", Topology.LINEAR)
        topology = getattr(raw_topology, "value", raw_topology)
        raw_id = getattr(document, "id", "")
        return cls(
            name=name or str(getattr(document, "name", "template")),
            sequence=sequence,
            topology=topology,
            id=str(raw_id) if raw_id else "",
        )


def as_template(value: PCRTemplate | str | object, *, name: str = "template") -> PCRTemplate:

    if isinstance(value, PCRTemplate):
        return value
    if isinstance(value, str):
        return PCRTemplate(name=name, sequence=value)
    return PCRTemplate.from_document(value, name=name)


@dataclass(frozen=True, slots=True)
class PCRPrimer:

    name: str
    annealing_sequence: str
    five_prime_extension: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise PCRValidationError("primer name cannot be empty")
        annealing = normalize_dna(self.annealing_sequence, allow_ambiguous=False)
        extension = normalize_dna(self.five_prime_extension, allow_ambiguous=False)
        if not annealing:
            raise PCRValidationError("primer annealing_sequence cannot be empty")
        object.__setattr__(self, "name", self.name.strip())
        object.__setattr__(self, "annealing_sequence", annealing)
        object.__setattr__(self, "five_prime_extension", extension)

    @property
    def full_sequence(self) -> str:
        return self.five_prime_extension + self.annealing_sequence

    @classmethod
    def from_any(cls, value: "PCRPrimer | str | object", *, default_name: str) -> "PCRPrimer":

        if isinstance(value, cls):
            return value
        if isinstance(value, str):
            return cls(default_name, value)
        try:
            return cls(
                name=str(getattr(value, "name")),
                annealing_sequence=str(getattr(value, "annealing_sequence")),
                five_prime_extension=str(getattr(value, "five_prime_extension", "")),
            )
        except (AttributeError, TypeError, ValueError) as error:
            raise TypeError(
                "primer must be DNA text or expose name, annealing_sequence, and five_prime_extension"
            ) from error


@dataclass(frozen=True, slots=True)
class PCRBindingSite:

    start: int
    end: int
    strand: int
    mismatches: int = 0
    wraps_origin: bool = False

    def __post_init__(self) -> None:
        if self.start < 0 or self.end < 0:
            raise PCRValidationError("binding coordinates cannot be negative")
        if self.strand not in (-1, 1):
            raise PCRValidationError("binding strand must be 1 or -1")
        if self.mismatches < 0:
            raise PCRValidationError("mismatches cannot be negative")


@dataclass(frozen=True, slots=True)
class PCRProduct:

    id: str
    name: str
    sequence: str
    forward_site: PCRBindingSite
    reverse_site: PCRBindingSite
    circular_template: bool = False

    def __post_init__(self) -> None:
        sequence = normalize_dna(self.sequence)
        if not sequence:
            raise PCRValidationError("PCR product sequence cannot be empty")
        if not self.id.strip() or not self.name.strip():
            raise PCRValidationError("PCR product id and name cannot be empty")
        object.__setattr__(self, "sequence", sequence)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def sha256(self) -> str:
        return sequence_hash(self.sequence)


@dataclass(frozen=True, slots=True)
class PCRProvenanceRecord:

    operation_id: str
    operation_type: str
    input_ids: tuple[str, ...]
    input_sequence_hashes: Mapping[str, str]
    parameters: Mapping[str, Any]
    output_sequence_hashes: tuple[str, ...] = ()
    warnings: tuple[str, ...] = ()
    engine_name: str = "opengene-pcr-planner"
    engine_version: str = "1"
    reproducibility_status: str = "reproducible"

    def __post_init__(self) -> None:
        if not self.operation_type.strip() or not self.engine_name.strip() or not self.engine_version.strip():
            raise PCRValidationError("provenance operation and engine values cannot be empty")
        input_ids = tuple(str(identifier) for identifier in self.input_ids)
        input_hashes = {str(key): str(value) for key, value in self.input_sequence_hashes.items()}
        parameters = _freeze(self.parameters)
        output_hashes = tuple(str(item) for item in self.output_sequence_hashes)
        object.__setattr__(self, "input_ids", input_ids)
        object.__setattr__(self, "input_sequence_hashes", MappingProxyType(input_hashes))
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "output_sequence_hashes", output_hashes)
        object.__setattr__(self, "warnings", tuple(str(item) for item in self.warnings))
        expected = _provenance_id(
            self.operation_type,
            input_ids,
            input_hashes,
            parameters,
            output_hashes,
            self.engine_name,
            self.engine_version,
        )
        if self.operation_id and self.operation_id != expected:
            raise PCRValidationError("operation_id does not match deterministic provenance content")
        object.__setattr__(self, "operation_id", expected)

    def to_dict(self) -> dict[str, Any]:
        return {
            "operation_id": self.operation_id,
            "operation_type": self.operation_type,
            "input_ids": list(self.input_ids),
            "input_sequence_hashes": dict(self.input_sequence_hashes),
            "parameters": _json_value(self.parameters),
            "output_sequence_hashes": list(self.output_sequence_hashes),
            "warnings": list(self.warnings),
            "engine_name": self.engine_name,
            "engine_version": self.engine_version,
            "reproducibility_status": self.reproducibility_status,
        }


def _provenance_id(
    operation_type: str,
    input_ids: tuple[str, ...],
    input_hashes: Mapping[str, str],
    parameters: Mapping[str, Any],
    output_hashes: tuple[str, ...],
    engine_name: str,
    engine_version: str,
) -> str:
    payload = {
        "operation_type": operation_type,
        "input_ids": list(input_ids),
        "input_hashes": _json_value(input_hashes),
        "parameters": _json_value(parameters),
        "output_hashes": list(output_hashes),
        "engine_name": engine_name,
        "engine_version": engine_version,
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return "op-" + sha256(encoded).hexdigest()[:32]


def make_pcr_provenance(
    operation_type: str,
    template: PCRTemplate,
    primers: tuple[PCRPrimer, ...] = (),
    *,
    parameters: Mapping[str, Any] | None = None,
    outputs: tuple[PCRProduct, ...] = (),
    warnings: tuple[str, ...] | list[str] = (),
    reproducibility_status: str = "reproducible",
) -> PCRProvenanceRecord:

    input_ids = (template.id, *(f"primer:{primer.name}" for primer in primers))
    input_hashes = {
        template.id: template.sha256,
        **{f"primer:{primer.name}": sequence_hash(primer.full_sequence) for primer in primers},
    }
    return PCRProvenanceRecord(
        operation_id="",
        operation_type=operation_type,
        input_ids=input_ids,
        input_sequence_hashes=input_hashes,
        parameters=parameters or {},
        output_sequence_hashes=tuple(product.sha256 for product in outputs),
        warnings=tuple(warnings),
        reproducibility_status=reproducibility_status,
    )


@dataclass(frozen=True, slots=True)
class PCRPlan:

    template: PCRTemplate
    forward_primer: PCRPrimer
    reverse_primer: PCRPrimer
    circular_template: bool
    max_mismatches: int
    products: tuple[PCRProduct, ...]
    warnings: tuple[str, ...]
    provenance: PCRProvenanceRecord

    @property
    def feasible(self) -> bool:
        return bool(self.products)


@dataclass(frozen=True, slots=True)
class SequenceMutation:

    start: int
    end: int
    replacement: str
    name: str = "mutation"

    def __post_init__(self) -> None:
        if isinstance(self.start, bool) or isinstance(self.end, bool) or not isinstance(self.start, int) or not isinstance(self.end, int):
            raise PCRValidationError("mutation coordinates must be integers")
        if self.start < 0 or self.end < self.start:
            raise PCRValidationError("mutation requires 0 <= start <= end")
        if not isinstance(self.name, str) or not self.name.strip():
            raise PCRValidationError("mutation name cannot be empty")
        object.__setattr__(self, "replacement", normalize_dna(self.replacement, allow_ambiguous=False))
        object.__setattr__(self, "name", self.name.strip())


@dataclass(frozen=True, slots=True)
class MutagenesisProduct:

    id: str
    name: str
    sequence: str
    topology: Topology
    mutations: tuple[SequenceMutation, ...]

    def __post_init__(self) -> None:
        sequence = normalize_dna(self.sequence)
        if not sequence:
            raise PCRValidationError("mutagenesis product sequence cannot be empty")
        if not self.id.strip() or not self.name.strip():
            raise PCRValidationError("mutagenesis product id and name cannot be empty")
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "mutations", tuple(self.mutations))

    @property
    def sha256(self) -> str:
        return sequence_hash(self.sequence)


@dataclass(frozen=True, slots=True)
class MutagenesisPlan:

    template: PCRTemplate
    mutations: tuple[SequenceMutation, ...]
    product: MutagenesisProduct | None
    pcr_plan: PCRPlan | None
    warnings: tuple[str, ...]
    provenance: PCRProvenanceRecord

    @property
    def feasible(self) -> bool:
        return self.product is not None
