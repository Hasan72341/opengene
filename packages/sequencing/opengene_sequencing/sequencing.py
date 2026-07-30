from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from importlib.util import find_spec
from math import isfinite
from pathlib import Path
from types import MappingProxyType
from typing import Any

from opengene_alignment import AlignmentResult, smith_waterman


class SequencingError(ValueError):
    pass


class Ab1UnavailableError(RuntimeError):
    pass


class Ab1ParseError(SequencingError):
    pass


_IUPAC_TO_BASES: dict[str, frozenset[str]] = {
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
_BASES_TO_IUPAC = {bases: symbol for symbol, bases in _IUPAC_TO_BASES.items()}
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
_CANONICAL_BASES = ("A", "C", "G", "T")


def normalize_dna(sequence: str, *, label: str = "sequence") -> str:

    if not isinstance(sequence, str):
        raise SequencingError(f"{label} must be a string")
    normalized = "".join(sequence.split()).upper().replace("U", "T")
    invalid = sorted(set(normalized).difference(_IUPAC_TO_BASES))
    if invalid:
        raise SequencingError(f"{label} contains unsupported DNA symbols: {''.join(invalid)}")
    return normalized


def reverse_complement(sequence: str) -> str:

    return "".join(_COMPLEMENT[base] for base in reversed(normalize_dna(sequence)))


def _freeze_metadata(metadata: Mapping[str, Any]) -> Mapping[str, Any]:
    return MappingProxyType(dict(metadata))


def decode_phred(quality_text: str, *, phred_offset: int = 33) -> tuple[int, ...]:

    if not isinstance(quality_text, str):
        raise SequencingError("quality_text must be a string")
    _validate_phred_offset(phred_offset)
    qualities = tuple(ord(character) - phred_offset for character in quality_text)
    _validate_qualities(qualities)
    return qualities


def encode_phred(qualities: Iterable[int], *, phred_offset: int = 33) -> str:

    _validate_phred_offset(phred_offset)
    normalized = tuple(qualities)
    _validate_qualities(normalized)
    return "".join(chr(quality + phred_offset) for quality in normalized)


def _validate_phred_offset(phred_offset: int) -> None:
    if isinstance(phred_offset, bool) or not isinstance(phred_offset, int):
        raise SequencingError("phred_offset must be an integer")
    if phred_offset < 0 or phred_offset > 127:
        raise SequencingError("phred_offset must be between 0 and 127")


def _validate_qualities(qualities: tuple[int, ...]) -> None:
    for quality in qualities:
        if isinstance(quality, bool) or not isinstance(quality, int) or not 0 <= quality <= 93:
            raise SequencingError("Phred qualities must be integers from 0 through 93")


@dataclass(frozen=True, slots=True)
class SequencingRead:

    name: str
    sequence: str
    qualities: tuple[int, ...] | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise SequencingError("read name cannot be empty")
        sequence = normalize_dna(self.sequence, label="read sequence")
        if not sequence:
            raise SequencingError("read sequence cannot be empty")
        object.__setattr__(self, "name", self.name.strip())
        object.__setattr__(self, "sequence", sequence)
        if self.qualities is not None:
            qualities = tuple(self.qualities)
            _validate_qualities(qualities)
            if len(qualities) != len(sequence):
                raise SequencingError("quality count must equal read sequence length")
            object.__setattr__(self, "qualities", qualities)
        if not isinstance(self.metadata, Mapping):
            raise SequencingError("read metadata must be a mapping")
        object.__setattr__(self, "metadata", _freeze_metadata(self.metadata))

    @classmethod
    def from_fastq(
        cls,
        name: str,
        sequence: str,
        quality_text: str,
        *,
        phred_offset: int = 33,
        metadata: Mapping[str, Any] | None = None,
    ) -> "SequencingRead":

        return cls(
            name=name,
            sequence=sequence,
            qualities=decode_phred(quality_text, phred_offset=phred_offset),
            metadata=metadata or {},
        )

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def mean_quality(self) -> float | None:

        return sum(self.qualities) / len(self.qualities) if self.qualities else None

    def fastq_quality(self, *, phred_offset: int = 33) -> str:

        if self.qualities is None:
            raise SequencingError(f"read {self.name!r} has no quality values")
        return encode_phred(self.qualities, phred_offset=phred_offset)

    def reverse_complement(self) -> "SequencingRead":

        return SequencingRead(
            name=self.name,
            sequence=reverse_complement(self.sequence),
            qualities=tuple(reversed(self.qualities)) if self.qualities is not None else None,
            metadata=self.metadata,
        )


SangerRead = SequencingRead


def parse_fastq(text: str, *, phred_offset: int = 33) -> tuple[SequencingRead, ...]:

    if not isinstance(text, str):
        raise SequencingError("FASTQ text must be a string")
    lines = text.splitlines()
    if len(lines) % 4:
        raise SequencingError("FASTQ input must contain complete four-line records")
    reads: list[SequencingRead] = []
    for line_index in range(0, len(lines), 4):
        header, sequence, separator, qualities = lines[line_index : line_index + 4]
        record_number = line_index // 4 + 1
        if not header.startswith("@") or len(header) == 1:
            raise SequencingError(f"FASTQ record {record_number} needs a non-empty @ header")
        if not separator.startswith("+"):
            raise SequencingError(f"FASTQ record {record_number} needs a + separator")
        repeated_name = separator[1:].strip()
        name = header[1:].strip()
        if repeated_name and repeated_name != name:
            raise SequencingError(f"FASTQ record {record_number} has mismatched + header")
        reads.append(
            SequencingRead.from_fastq(
                name,
                sequence,
                qualities,
                phred_offset=phred_offset,
                metadata={"source_format": "FASTQ"},
            )
        )
    return tuple(reads)


def write_fastq(reads: Iterable[SequencingRead], *, phred_offset: int = 33) -> str:

    rows: list[str] = []
    for read in reads:
        if not isinstance(read, SequencingRead):
            raise TypeError("reads must contain SequencingRead values")
        rows.extend((f"@{read.name}", read.sequence, "+", read.fastq_quality(phred_offset=phred_offset)))
    return "\n".join(rows) + ("\n" if rows else "")


@dataclass(frozen=True, slots=True)
class MappedRead:

    read: SequencingRead
    alignment: AlignmentResult
    strand: int = 1

    def __post_init__(self) -> None:
        if not isinstance(self.read, SequencingRead):
            raise SequencingError("mapped read requires a SequencingRead")
        if not isinstance(self.alignment, AlignmentResult):
            raise SequencingError("mapped read requires an AlignmentResult")
        if self.strand not in (-1, 1):
            raise SequencingError("mapped read strand must be 1 or -1")
        if self.alignment.query_end > self.read.length:
            raise SequencingError("alignment query coordinates exceed read length")

    @property
    def oriented_sequence(self) -> str:
        return self.read.sequence if self.strand == 1 else reverse_complement(self.read.sequence)

    @property
    def oriented_qualities(self) -> tuple[int, ...] | None:
        if self.read.qualities is None:
            return None
        return self.read.qualities if self.strand == 1 else tuple(reversed(self.read.qualities))


def align_read_to_reference(
    reference: str,
    read: SequencingRead,
    *,
    allow_reverse_complement: bool = True,
    match_score: float = 2,
    mismatch_score: float = -1,
    gap_score: float = -2,
) -> MappedRead:

    reference = normalize_dna(reference, label="reference")
    if not reference:
        raise SequencingError("reference cannot be empty")
    if not isinstance(read, SequencingRead):
        raise TypeError("read must be a SequencingRead")
    forward = smith_waterman(
        reference,
        read.sequence,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_score=gap_score,
        ambiguity_matches=True,
    )
    selected = forward
    strand = 1
    if allow_reverse_complement:
        reverse = smith_waterman(
            reference,
            reverse_complement(read.sequence),
            match_score=match_score,
            mismatch_score=mismatch_score,
            gap_score=gap_score,
            ambiguity_matches=True,
        )
        if (reverse.score, reverse.matches, reverse.alignment_length) > (
            forward.score,
            forward.matches,
            forward.alignment_length,
        ):
            selected = reverse
            strand = -1
    return MappedRead(read=read, alignment=selected, strand=strand)


map_read = align_read_to_reference


@dataclass(frozen=True, slots=True)
class ConsensusCall:

    position: int
    base: str
    coverage: int
    base_weights: tuple[tuple[str, float], ...]
    deletion_weight: float = 0.0

    def __post_init__(self) -> None:
        if isinstance(self.position, bool) or not isinstance(self.position, int) or self.position < 0:
            raise SequencingError("consensus position must be a non-negative integer")
        if self.base not in _IUPAC_TO_BASES:
            raise SequencingError("consensus base must be an IUPAC DNA symbol")
        if isinstance(self.coverage, bool) or not isinstance(self.coverage, int) or self.coverage < 0:
            raise SequencingError("coverage must be a non-negative integer")
        expected_bases = set(_CANONICAL_BASES)
        weights = tuple(self.base_weights)
        if {base for base, _ in weights} != expected_bases or len(weights) != len(expected_bases):
            raise SequencingError("base_weights must contain A, C, G, and T exactly once")
        for _, weight in weights:
            if not isinstance(weight, (int, float)) or not isfinite(weight) or weight < 0:
                raise SequencingError("base weights must be finite non-negative numbers")
        if (
            not isinstance(self.deletion_weight, (int, float))
            or not isfinite(self.deletion_weight)
            or self.deletion_weight < 0
        ):
            raise SequencingError("deletion_weight must be a finite non-negative number")

    def weight_for(self, base: str) -> float:

        normalized = base.upper()
        for observed_base, weight in self.base_weights:
            if observed_base == normalized:
                return weight
        raise KeyError(f"{base!r} is not a canonical DNA base")


@dataclass(frozen=True, slots=True)
class ConsensusInsertion:

    position: int
    sequence: str
    read_name: str
    strand: int

    def __post_init__(self) -> None:
        if isinstance(self.position, bool) or not isinstance(self.position, int) or self.position < 0:
            raise SequencingError("insertion position must be a non-negative integer")
        sequence = normalize_dna(self.sequence, label="insertion sequence")
        if not sequence:
            raise SequencingError("insertion sequence cannot be empty")
        if not self.read_name:
            raise SequencingError("insertion read_name cannot be empty")
        if self.strand not in (-1, 1):
            raise SequencingError("insertion strand must be 1 or -1")
        object.__setattr__(self, "sequence", sequence)


@dataclass(frozen=True, slots=True)
class ConsensusResult:

    reference: str
    sequence: str
    calls: tuple[ConsensusCall, ...]
    mapped_reads: tuple[MappedRead, ...]
    insertions: tuple[ConsensusInsertion, ...] = ()
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        reference = normalize_dna(self.reference, label="reference")
        sequence = normalize_dna(self.sequence, label="consensus sequence")
        if len(sequence) != len(reference):
            raise SequencingError("reference-anchored consensus must equal reference length")
        calls = tuple(self.calls)
        if len(calls) != len(reference):
            raise SequencingError("consensus requires one call per reference base")
        if any(call.position != index for index, call in enumerate(calls)):
            raise SequencingError("consensus calls must be ordered by reference position")
        if sequence != "".join(call.base for call in calls):
            raise SequencingError("consensus sequence does not match consensus calls")
        if not all(isinstance(read, MappedRead) for read in self.mapped_reads):
            raise SequencingError("mapped_reads must contain MappedRead values")
        if not all(isinstance(insertion, ConsensusInsertion) for insertion in self.insertions):
            raise SequencingError("insertions must contain ConsensusInsertion values")
        object.__setattr__(self, "reference", reference)
        object.__setattr__(self, "sequence", sequence)
        object.__setattr__(self, "calls", calls)
        object.__setattr__(self, "mapped_reads", tuple(self.mapped_reads))
        object.__setattr__(self, "insertions", tuple(self.insertions))
        object.__setattr__(self, "warnings", tuple(str(warning) for warning in self.warnings))

    @property
    def coverage(self) -> tuple[int, ...]:

        return tuple(call.coverage for call in self.calls)

    @property
    def ambiguous_positions(self) -> tuple[int, ...]:

        return tuple(
            call.position
            for call in self.calls
            if call.base not in _CANONICAL_BASES
        )


def _quality_weight(quality: int | None, *, quality_weighted: bool) -> float:
    if not quality_weighted or quality is None:
        return 1.0
    return max(0.001, 1.0 - 10 ** (-quality / 10))


def _choose_consensus_base(
    weights: Mapping[str, float],
    deletion_weight: float,
    *,
    minimum_base_fraction: float,
    ambiguity_fraction: float,
) -> str:
    total = sum(weights.values())
    if total <= 0:
        return "N"
    best_weight = max(weights.values())
    if deletion_weight > best_weight:
        return "N"
    best_base = next(base for base in _CANONICAL_BASES if weights[base] == best_weight)
    if best_weight / total >= minimum_base_fraction:
        return best_base
    supported = frozenset(
        base for base in _CANONICAL_BASES if weights[base] / total >= ambiguity_fraction
    )
    return _BASES_TO_IUPAC.get(supported, "N")


def consensus_from_mapped_reads(
    reference: str,
    mapped_reads: Iterable[MappedRead],
    *,
    minimum_coverage: int = 1,
    minimum_base_fraction: float = 0.6,
    ambiguity_fraction: float = 0.2,
    quality_weighted: bool = True,
) -> ConsensusResult:

    reference = normalize_dna(reference, label="reference")
    if not reference:
        raise SequencingError("reference cannot be empty")
    if isinstance(minimum_coverage, bool) or not isinstance(minimum_coverage, int) or minimum_coverage < 1:
        raise SequencingError("minimum_coverage must be a positive integer")
    for name, value in (
        ("minimum_base_fraction", minimum_base_fraction),
        ("ambiguity_fraction", ambiguity_fraction),
    ):
        if not isinstance(value, (int, float)) or not isfinite(value) or not 0 < value <= 1:
            raise SequencingError(f"{name} must be in the interval (0, 1]")
    if not isinstance(quality_weighted, bool):
        raise SequencingError("quality_weighted must be a boolean")

    reads = tuple(mapped_reads)
    if not all(isinstance(mapped_read, MappedRead) for mapped_read in reads):
        raise TypeError("mapped_reads must contain MappedRead values")
    reference_length = len(reference)
    weights = [{base: 0.0 for base in _CANONICAL_BASES} for _ in reference]
    coverage = [0 for _ in reference]
    deletion_weights = [0.0 for _ in reference]
    insertion_observations: list[ConsensusInsertion] = []
    deleted_positions = 0

    for mapped_read in reads:
        alignment = mapped_read.alignment
        if alignment.reference_end > reference_length:
            raise SequencingError("alignment reference coordinates exceed reference length")
        query_position = alignment.query_start
        reference_position = alignment.reference_start
        oriented_qualities = mapped_read.oriented_qualities
        insertion_start: int | None = None
        insertion_bases: list[str] = []

        def flush_insertion() -> None:
            nonlocal insertion_start, insertion_bases
            if insertion_start is not None and insertion_bases:
                insertion_observations.append(
                    ConsensusInsertion(
                        position=insertion_start,
                        sequence="".join(insertion_bases),
                        read_name=mapped_read.read.name,
                        strand=mapped_read.strand,
                    )
                )
            insertion_start = None
            insertion_bases = []

        for reference_base, query_base in zip(alignment.aligned_reference, alignment.aligned_query):
            if reference_base == "-":
                if query_base != "-":
                    if insertion_start is None:
                        insertion_start = reference_position
                    insertion_bases.append(query_base)
                    query_position += 1
                continue
            flush_insertion()
            if reference_position >= reference_length:
                raise SequencingError("alignment reference coordinate exceeds reference length")
            coverage[reference_position] += 1
            if query_base == "-":
                deletion_weights[reference_position] += 1.0
                deleted_positions += 1
            else:
                quality = oriented_qualities[query_position] if oriented_qualities is not None else None
                weight = _quality_weight(quality, quality_weighted=quality_weighted)
                possible_bases = _IUPAC_TO_BASES[query_base]
                distributed_weight = weight / len(possible_bases)
                for base in possible_bases:
                    weights[reference_position][base] += distributed_weight
                query_position += 1
            reference_position += 1
        flush_insertion()

    calls: list[ConsensusCall] = []
    for position in range(reference_length):
        call_weights = weights[position]
        if coverage[position] < minimum_coverage:
            base = "N"
        else:
            base = _choose_consensus_base(
                call_weights,
                deletion_weights[position],
                minimum_base_fraction=float(minimum_base_fraction),
                ambiguity_fraction=float(ambiguity_fraction),
            )
        calls.append(
            ConsensusCall(
                position=position,
                base=base,
                coverage=coverage[position],
                base_weights=tuple((canonical_base, call_weights[canonical_base]) for canonical_base in _CANONICAL_BASES),
                deletion_weight=deletion_weights[position],
            )
        )

    warnings: list[str] = []
    uncalled = sum(call.base == "N" for call in calls)
    if uncalled:
        warnings.append(f"{uncalled} reference position(s) have insufficient or conflicting evidence")
    if insertion_observations:
        warnings.append(
            f"{len(insertion_observations)} insertion observation(s) were retained separately from the reference-anchored consensus"
        )
    if deleted_positions:
        warnings.append(
            f"{deleted_positions} deletion observation(s) were not assembled into the reference-anchored consensus"
        )
    return ConsensusResult(
        reference=reference,
        sequence="".join(call.base for call in calls),
        calls=tuple(calls),
        mapped_reads=reads,
        insertions=tuple(insertion_observations),
        warnings=tuple(warnings),
    )


def build_consensus(
    reference: str,
    reads: Iterable[SequencingRead | MappedRead],
    *,
    minimum_coverage: int = 1,
    minimum_base_fraction: float = 0.6,
    ambiguity_fraction: float = 0.2,
    quality_weighted: bool = True,
    allow_reverse_complement: bool = True,
) -> ConsensusResult:

    reference = normalize_dna(reference, label="reference")
    collected: list[MappedRead] = []
    for read in reads:
        if isinstance(read, MappedRead):
            collected.append(read)
        elif isinstance(read, SequencingRead):
            collected.append(
                align_read_to_reference(
                    reference,
                    read,
                    allow_reverse_complement=allow_reverse_complement,
                )
            )
        else:
            raise TypeError("reads must contain SequencingRead or MappedRead values")
    return consensus_from_mapped_reads(
        reference,
        collected,
        minimum_coverage=minimum_coverage,
        minimum_base_fraction=minimum_base_fraction,
        ambiguity_fraction=ambiguity_fraction,
        quality_weighted=quality_weighted,
    )


call_consensus = build_consensus
consensus_from_reads = build_consensus


@dataclass(frozen=True, slots=True)
class Ab1ReaderStatus:

    available: bool
    provider: str | None
    message: str


def ab1_reader_status() -> Ab1ReaderStatus:

    try:
        available = find_spec("Bio.SeqIO") is not None
    except (ImportError, ModuleNotFoundError, ValueError):
        available = False
    if available:
        return Ab1ReaderStatus(
            available=True,
            provider="Biopython SeqIO (abi)",
            message="AB1 import is available through the optional Biopython reader.",
        )
    return Ab1ReaderStatus(
        available=False,
        provider=None,
        message="AB1 import requires the optional Biopython dependency (pip install biopython).",
    )


def require_ab1_reader() -> Ab1ReaderStatus:

    status = ab1_reader_status()
    if not status.available:
        raise Ab1UnavailableError(status.message)
    return status


def read_ab1(path: str | Path) -> SequencingRead:

    require_ab1_reader()
    try:
        from Bio import SeqIO  # type: ignore[import-not-found]

        record = SeqIO.read(str(path), "abi")
        qualities = record.letter_annotations.get("phred_quality")
        return SequencingRead(
            name=record.id or Path(path).stem,
            sequence=str(record.seq),
            qualities=tuple(int(quality) for quality in qualities) if qualities is not None else None,
            metadata={"source_format": "AB1", "source_name": Path(path).name},
        )
    except Ab1UnavailableError:
        raise
    except Exception as error:
        raise Ab1ParseError(f"could not read AB1 file {Path(path).name!r}: {error}") from error
