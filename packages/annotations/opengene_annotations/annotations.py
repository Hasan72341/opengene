from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from itertools import product
from types import MappingProxyType
from typing import Any, Literal


class AnnotationError(ValueError):
    pass


_DNA_IUPAC: Mapping[str, frozenset[str]] = MappingProxyType(
    {
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
)

_COMPLEMENT: Mapping[str, str] = MappingProxyType(
    {
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
)


def normalize_dna(sequence: str, *, allow_ambiguous: bool = True) -> str:

    if not isinstance(sequence, str):
        raise AnnotationError("sequence must be a string")
    normalised = "".join(sequence.split()).upper().replace("U", "T")
    permitted = set(_DNA_IUPAC) if allow_ambiguous else set("ACGT")
    invalid = "".join(sorted(set(normalised).difference(permitted)))
    if invalid:
        raise AnnotationError(f"sequence contains unsupported DNA symbols: {invalid}")
    return normalised


def reverse_complement(sequence: str) -> str:

    return "".join(_COMPLEMENT[base] for base in reversed(normalize_dna(sequence)))


STANDARD_GENETIC_CODE: Mapping[str, str] = MappingProxyType(
    {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
        "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
)


@dataclass(frozen=True, slots=True)
class CodonPreference:

    codon: str
    weight: float

    def __post_init__(self) -> None:
        codon = normalize_dna(self.codon, allow_ambiguous=False)
        if len(codon) != 3:
            raise AnnotationError("a codon preference requires exactly three bases")
        if isinstance(self.weight, bool) or not isinstance(self.weight, (int, float)):
            raise AnnotationError("codon preference weight must be numeric")
        if self.weight <= 0:
            raise AnnotationError("codon preference weight must be greater than zero")
        object.__setattr__(self, "codon", codon)
        object.__setattr__(self, "weight", float(self.weight))


_GENERIC_CODON_USAGE: dict[str, tuple[tuple[str, float], ...]] = {
    "A": (("GCT", 1.00), ("GCC", 0.95), ("GCA", 0.90), ("GCG", 0.85)),
    "R": (("CGT", 1.00), ("CGC", 0.95), ("CGA", 0.90), ("CGG", 0.85), ("AGA", 0.80), ("AGG", 0.75)),
    "N": (("AAT", 1.00), ("AAC", 0.95)),
    "D": (("GAT", 1.00), ("GAC", 0.95)),
    "C": (("TGT", 1.00), ("TGC", 0.95)),
    "Q": (("CAA", 1.00), ("CAG", 0.95)),
    "E": (("GAA", 1.00), ("GAG", 0.95)),
    "G": (("GGT", 1.00), ("GGC", 0.95), ("GGA", 0.90), ("GGG", 0.85)),
    "H": (("CAT", 1.00), ("CAC", 0.95)),
    "I": (("ATT", 1.00), ("ATC", 0.95), ("ATA", 0.90)),
    "L": (("CTG", 1.00), ("CTC", 0.95), ("CTT", 0.90), ("CTA", 0.85), ("TTG", 0.80), ("TTA", 0.75)),
    "K": (("AAA", 1.00), ("AAG", 0.95)),
    "M": (("ATG", 1.00),),
    "F": (("TTT", 1.00), ("TTC", 0.95)),
    "P": (("CCT", 1.00), ("CCC", 0.95), ("CCA", 0.90), ("CCG", 0.85)),
    "S": (("TCT", 1.00), ("TCC", 0.95), ("TCA", 0.90), ("TCG", 0.85), ("AGT", 0.80), ("AGC", 0.75)),
    "T": (("ACT", 1.00), ("ACC", 0.95), ("ACA", 0.90), ("ACG", 0.85)),
    "W": (("TGG", 1.00),),
    "Y": (("TAT", 1.00), ("TAC", 0.95)),
    "V": (("GTT", 1.00), ("GTC", 0.95), ("GTA", 0.90), ("GTG", 0.85)),
    "*": (("TAA", 1.00), ("TAG", 0.95), ("TGA", 0.90)),
}

GENERIC_CODON_USAGE: Mapping[str, tuple[CodonPreference, ...]] = MappingProxyType(
    {
        amino_acid: tuple(CodonPreference(codon, weight) for codon, weight in preferences)
        for amino_acid, preferences in _GENERIC_CODON_USAGE.items()
    }
)


@dataclass(frozen=True, slots=True)
class TranslationResult:

    protein_sequence: str
    codons: tuple[str, ...]
    frame: int
    strand: int
    trailing_bases: str = ""
    terminated: bool = False
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.frame not in (0, 1, 2):
            raise AnnotationError("translation frame must be 0, 1, or 2")
        if self.strand not in (-1, 1):
            raise AnnotationError("translation strand must be 1 or -1")
        codons = tuple(normalize_dna(codon, allow_ambiguous=True) for codon in self.codons)
        if any(len(codon) != 3 for codon in codons):
            raise AnnotationError("translation codons must each have length three")
        if any(len(amino_acid) != 1 or not (amino_acid.isalpha() or amino_acid == "*") for amino_acid in self.protein_sequence):
            raise AnnotationError("protein_sequence contains unsupported symbols")
        trailing = normalize_dna(self.trailing_bases, allow_ambiguous=True)
        object.__setattr__(self, "protein_sequence", self.protein_sequence.upper())
        object.__setattr__(self, "codons", codons)
        object.__setattr__(self, "trailing_bases", trailing)
        object.__setattr__(self, "warnings", tuple(str(warning) for warning in self.warnings))

    @property
    def translated_base_count(self) -> int:

        return 3 * len(self.codons)


@dataclass(frozen=True, slots=True)
class AnnotationLocation:

    start: int
    end: int
    strand: int = 1
    wraps_origin: bool = False

    def __post_init__(self) -> None:
        if isinstance(self.start, bool) or isinstance(self.end, bool):
            raise AnnotationError("annotation coordinates must be integers")
        if not isinstance(self.start, int) or not isinstance(self.end, int):
            raise AnnotationError("annotation coordinates must be integers")
        if self.start < 0 or self.end < 0:
            raise AnnotationError("annotation coordinates cannot be negative")
        if self.strand not in (-1, 1):
            raise AnnotationError("annotation strand must be 1 or -1")
        if self.wraps_origin and self.start <= self.end:
            raise AnnotationError("an origin-spanning annotation requires start > end")
        if not self.wraps_origin and self.start > self.end:
            raise AnnotationError("a non-wrapping annotation requires start <= end")

    def validate(self, sequence_length: int, *, circular: bool = False) -> None:
        if isinstance(sequence_length, bool) or not isinstance(sequence_length, int) or sequence_length < 0:
            raise AnnotationError("sequence_length must be a non-negative integer")
        if self.start > sequence_length or self.end > sequence_length:
            raise AnnotationError("annotation coordinates exceed sequence length")
        if self.wraps_origin:
            if not circular:
                raise AnnotationError("origin-spanning annotations require circular sequence topology")
            if self.start >= sequence_length or self.end == 0:
                raise AnnotationError("invalid origin-spanning annotation for sequence length")

    def segments(self, sequence_length: int, *, circular: bool = False) -> tuple[tuple[int, int], ...]:

        self.validate(sequence_length, circular=circular)
        if not self.wraps_origin:
            return ((self.start, self.end),)
        return ((self.start, sequence_length), (0, self.end))

    def length(self, sequence_length: int, *, circular: bool = False) -> int:

        return sum(end - start for start, end in self.segments(sequence_length, circular=circular))


@dataclass(frozen=True, slots=True)
class ORF:

    location: AnnotationLocation
    frame: int
    nucleotide_sequence: str
    protein_sequence: str
    start_codon: str
    stop_codon: str | None = None
    circular: bool = False

    def __post_init__(self) -> None:
        if self.frame not in (0, 1, 2):
            raise AnnotationError("ORF frame must be 0, 1, or 2")
        nucleotide = normalize_dna(self.nucleotide_sequence, allow_ambiguous=True)
        start = normalize_dna(self.start_codon, allow_ambiguous=False)
        if len(nucleotide) < 3 or len(nucleotide) % 3:
            raise AnnotationError("an ORF sequence must contain complete codons")
        if len(start) != 3 or nucleotide[:3] != start:
            raise AnnotationError("ORF start_codon must equal the first sequence codon")
        if self.stop_codon is not None:
            stop = normalize_dna(self.stop_codon, allow_ambiguous=False)
            if len(stop) != 3 or nucleotide[-3:] != stop:
                raise AnnotationError("ORF stop_codon must equal the final sequence codon")
            object.__setattr__(self, "stop_codon", stop)
        protein = "".join(self.protein_sequence.split()).upper()
        if any(len(amino_acid) != 1 or not amino_acid.isalpha() for amino_acid in protein):
            raise AnnotationError("ORF protein_sequence must contain one-letter amino acids only")
        object.__setattr__(self, "nucleotide_sequence", nucleotide)
        object.__setattr__(self, "start_codon", start)
        object.__setattr__(self, "protein_sequence", protein)

    @property
    def start(self) -> int:
        return self.location.start

    @property
    def end(self) -> int:
        return self.location.end

    @property
    def strand(self) -> int:
        return self.location.strand

    @property
    def wraps_origin(self) -> bool:
        return self.location.wraps_origin

    @property
    def reading_frame(self) -> int:

        return (self.frame + 1) if self.strand == 1 else -(self.frame + 1)

    @property
    def nucleotide_length(self) -> int:
        return len(self.nucleotide_sequence)

    @property
    def protein_length(self) -> int:
        return len(self.protein_sequence)

    @property
    def has_terminal_stop(self) -> bool:
        return self.stop_codon is not None


def _freeze_metadata(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze_metadata(item) for key, item in value.items()})
    if isinstance(value, (tuple, list)):
        return tuple(_freeze_metadata(item) for item in value)
    if isinstance(value, set):
        return tuple(sorted((_freeze_metadata(item) for item in value), key=repr))
    return value


@dataclass(frozen=True, slots=True)
class DetectedFeature:

    feature_type: str
    name: str
    location: AnnotationLocation
    qualifiers: Mapping[str, Any] = field(default_factory=dict)
    confidence: float = 1.0
    evidence: tuple[str, ...] = ()
    source: str = "opengene_annotations"

    def __post_init__(self) -> None:
        if not self.feature_type.strip():
            raise AnnotationError("detected feature_type cannot be empty")
        if not self.name.strip():
            raise AnnotationError("detected feature name cannot be empty")
        if isinstance(self.confidence, bool) or not isinstance(self.confidence, (int, float)):
            raise AnnotationError("feature confidence must be numeric")
        if not 0.0 <= float(self.confidence) <= 1.0:
            raise AnnotationError("feature confidence must be between zero and one")
        if not self.source.strip():
            raise AnnotationError("detected feature source cannot be empty")
        object.__setattr__(self, "feature_type", self.feature_type.strip())
        object.__setattr__(self, "name", self.name.strip())
        object.__setattr__(self, "qualifiers", _freeze_metadata(self.qualifiers))
        object.__setattr__(self, "confidence", float(self.confidence))
        object.__setattr__(self, "evidence", tuple(str(item) for item in self.evidence))
        object.__setattr__(self, "source", self.source.strip())


@dataclass(frozen=True, slots=True)
class FeatureDetectionResult:

    sequence_length: int
    features: tuple[DetectedFeature, ...]
    orfs: tuple[ORF, ...] = ()
    warnings: tuple[str, ...] = ()
    algorithm: str = "opengene-orf-scan-v1"

    def __post_init__(self) -> None:
        if isinstance(self.sequence_length, bool) or not isinstance(self.sequence_length, int) or self.sequence_length < 0:
            raise AnnotationError("sequence_length must be a non-negative integer")
        features = tuple(self.features)
        orfs = tuple(self.orfs)
        if not all(isinstance(feature, DetectedFeature) for feature in features):
            raise AnnotationError("features must be DetectedFeature values")
        if not all(isinstance(orf, ORF) for orf in orfs):
            raise AnnotationError("orfs must be ORF values")
        for feature in features:
            feature.location.validate(self.sequence_length, circular=feature.location.wraps_origin)
        for orf in orfs:
            orf.location.validate(self.sequence_length, circular=orf.circular)
        if not self.algorithm.strip():
            raise AnnotationError("feature detection algorithm cannot be empty")
        object.__setattr__(self, "features", features)
        object.__setattr__(self, "orfs", orfs)
        object.__setattr__(self, "warnings", tuple(str(warning) for warning in self.warnings))
        object.__setattr__(self, "algorithm", self.algorithm.strip())


@dataclass(frozen=True, slots=True)
class CodonOptimizationResult:

    protein_sequence: str
    dna_sequence: str
    codons: tuple[str, ...]
    usage_table_name: str
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        protein = _normalize_protein(self.protein_sequence)
        dna = normalize_dna(self.dna_sequence, allow_ambiguous=False)
        codons = tuple(normalize_dna(codon, allow_ambiguous=False) for codon in self.codons)
        if any(len(codon) != 3 for codon in codons):
            raise AnnotationError("optimisation codons must each have length three")
        if dna != "".join(codons):
            raise AnnotationError("optimisation DNA must equal the concatenated codons")
        if len(codons) != len(protein):
            raise AnnotationError("optimisation requires one codon per amino acid")
        if not self.usage_table_name.strip():
            raise AnnotationError("usage_table_name cannot be empty")
        object.__setattr__(self, "protein_sequence", protein)
        object.__setattr__(self, "dna_sequence", dna)
        object.__setattr__(self, "codons", codons)
        object.__setattr__(self, "usage_table_name", self.usage_table_name.strip())
        object.__setattr__(self, "warnings", tuple(str(warning) for warning in self.warnings))


def _normalise_genetic_code(genetic_code: Mapping[str, str] | None) -> Mapping[str, str]:
    if genetic_code is None or genetic_code is STANDARD_GENETIC_CODE:
        return STANDARD_GENETIC_CODE
    normalised: dict[str, str] = {}
    for raw_codon, raw_amino_acid in genetic_code.items():
        codon = normalize_dna(str(raw_codon), allow_ambiguous=False)
        amino_acid = str(raw_amino_acid).strip().upper()
        if len(codon) != 3 or len(amino_acid) != 1 or not (amino_acid.isalpha() or amino_acid == "*"):
            raise AnnotationError("genetic code entries must map DNA codons to one-letter amino acids")
        normalised[codon] = amino_acid
    if not normalised:
        raise AnnotationError("genetic_code cannot be empty")
    return MappingProxyType(normalised)


def _translate_codon(codon: str, genetic_code: Mapping[str, str], unknown_symbol: str) -> str:

    if codon in genetic_code:
        return genetic_code[codon]
    possibilities = ("".join(parts) for parts in product(*(_DNA_IUPAC[base] for base in codon)))
    translated = {genetic_code.get(possibility, unknown_symbol) for possibility in possibilities}
    return translated.pop() if len(translated) == 1 else unknown_symbol


def translate_dna(
    sequence: str,
    frame: int = 0,
    *,
    strand: int = 1,
    genetic_code: Mapping[str, str] | None = None,
    stop_at_stop: bool = False,
    include_stop: bool = True,
    unknown_symbol: str = "X",
) -> TranslationResult:

    if frame not in (0, 1, 2):
        raise AnnotationError("translation frame must be 0, 1, or 2")
    if strand not in (-1, 1):
        raise AnnotationError("translation strand must be 1 or -1")
    unknown = str(unknown_symbol).strip().upper()
    if len(unknown) != 1 or not (unknown.isalpha() or unknown == "*"):
        raise AnnotationError("unknown_symbol must be one amino-acid character")

    normalised = normalize_dna(sequence, allow_ambiguous=True)
    oriented = normalised if strand == 1 else reverse_complement(normalised)
    code = _normalise_genetic_code(genetic_code)
    codons: list[str] = []
    amino_acids: list[str] = []
    terminated = False

    offset = frame
    while offset + 3 <= len(oriented):
        codon = oriented[offset : offset + 3]
        amino_acid = _translate_codon(codon, code, unknown)
        codons.append(codon)
        if amino_acid == "*":
            terminated = True
            if include_stop:
                amino_acids.append(amino_acid)
            if stop_at_stop:
                offset += 3
                break
        else:
            amino_acids.append(amino_acid)
        offset += 3

    trailing = oriented[offset:]
    warnings: list[str] = []
    if trailing:
        if stop_at_stop and terminated:
            warnings.append(f"did not translate {len(trailing)} base(s) after a stop codon")
        else:
            warnings.append(f"ignored {len(trailing)} trailing base(s) outside complete codons")
    return TranslationResult(
        protein_sequence="".join(amino_acids),
        codons=tuple(codons),
        frame=frame,
        strand=strand,
        trailing_bases=trailing,
        terminated=terminated,
        warnings=tuple(warnings),
    )


def translate(sequence: str, frame: int = 0, **kwargs: Any) -> str:

    return translate_dna(sequence, frame, **kwargs).protein_sequence


def _normalise_triplets(values: Iterable[str], *, label: str) -> frozenset[str]:
    if isinstance(values, str):
        values = (values,)
    triplets = frozenset(normalize_dna(value, allow_ambiguous=False) for value in values)
    if not triplets or any(len(triplet) != 3 for triplet in triplets):
        raise AnnotationError(f"{label} must contain at least one three-base DNA codon")
    return triplets


def _resolve_circular(topology: Literal["linear", "circular"] | str, circular: bool | None) -> bool:
    raw_topology = getattr(topology, "value", topology)
    normalised_topology = str(raw_topology).strip().lower()
    if normalised_topology not in {"linear", "circular"}:
        raise AnnotationError("topology must be 'linear' or 'circular'")
    from_topology = normalised_topology == "circular"
    if circular is None:
        return from_topology
    if not isinstance(circular, bool):
        raise AnnotationError("circular must be a boolean when supplied")
    if normalised_topology != "linear" and circular != from_topology:
        raise AnnotationError("topology and circular specify conflicting sequence topologies")
    return circular


def _oriented_slice(sequence: str, start: int, length: int, *, circular: bool) -> str:
    if not circular or start + length <= len(sequence):
        return sequence[start : start + length]
    split = len(sequence) - start
    return sequence[start:] + sequence[: length - split]


def _location_from_oriented_span(
    oriented_start: int,
    span: int,
    *,
    sequence_length: int,
    strand: int,
    circular: bool,
) -> AnnotationLocation:

    if not circular:
        if strand == 1:
            return AnnotationLocation(oriented_start, oriented_start + span, strand=1)
        return AnnotationLocation(
            sequence_length - (oriented_start + span),
            sequence_length - oriented_start,
            strand=-1,
        )

    remains_linear = oriented_start + span <= sequence_length
    if strand == 1:
        if remains_linear:
            return AnnotationLocation(oriented_start, oriented_start + span, strand=1)
        return AnnotationLocation(
            oriented_start,
            (oriented_start + span) % sequence_length,
            strand=1,
            wraps_origin=True,
        )
    if remains_linear:
        return AnnotationLocation(
            sequence_length - (oriented_start + span),
            sequence_length - oriented_start,
            strand=-1,
        )
    return AnnotationLocation(
        (sequence_length - (oriented_start + span)) % sequence_length,
        sequence_length - oriented_start,
        strand=-1,
        wraps_origin=True,
    )


def _build_orf(
    oriented: str,
    start: int,
    span: int,
    *,
    frame: int,
    strand: int,
    circular: bool,
    stop_codon: str | None,
    genetic_code: Mapping[str, str],
) -> ORF:
    nucleotide = _oriented_slice(oriented, start, span, circular=circular)
    coding_sequence = nucleotide[:-3] if stop_codon is not None else nucleotide
    protein = translate_dna(coding_sequence, genetic_code=genetic_code, include_stop=False).protein_sequence
    return ORF(
        location=_location_from_oriented_span(
            start,
            span,
            sequence_length=len(oriented),
            strand=strand,
            circular=circular,
        ),
        frame=frame,
        nucleotide_sequence=nucleotide,
        protein_sequence=protein,
        start_codon=nucleotide[:3],
        stop_codon=stop_codon,
        circular=circular,
    )


def _scan_oriented_orfs(
    oriented: str,
    *,
    strand: int,
    circular: bool,
    start_codons: frozenset[str],
    stop_codons: frozenset[str],
    min_protein_length: int,
    include_partial: bool,
    genetic_code: Mapping[str, str],
) -> list[ORF]:

    sequence_length = len(oriented)
    results: list[ORF] = []
    for frame in range(3):
        for start in range(frame, sequence_length, 3):
            first_codon = _oriented_slice(oriented, start, 3, circular=circular)
            if len(first_codon) != 3 or first_codon not in start_codons:
                continue

            terminal_stop: str | None = None
            span: int | None = None
            candidate_span = 6
            while True:
                if circular:
                    if candidate_span >= sequence_length:
                        break
                elif start + candidate_span > sequence_length:
                    break
                codon = _oriented_slice(
                    oriented,
                    start + candidate_span - 3,
                    3,
                    circular=circular,
                )
                if codon in stop_codons:
                    terminal_stop = codon
                    span = candidate_span
                    break
                candidate_span += 3

            if span is None and include_partial:
                if circular:
                    span = ((sequence_length - 1) // 3) * 3
                else:
                    span = ((sequence_length - start) // 3) * 3
                if span < 3:
                    continue

            if span is None:
                continue
            orf = _build_orf(
                oriented,
                start,
                span,
                frame=frame,
                strand=strand,
                circular=circular,
                stop_codon=terminal_stop,
                genetic_code=genetic_code,
            )
            if orf.protein_length >= min_protein_length:
                results.append(orf)
    return results


def find_orfs(
    sequence: str,
    *,
    topology: Literal["linear", "circular"] | str = "linear",
    circular: bool | None = None,
    start_codons: Iterable[str] = ("ATG",),
    stop_codons: Iterable[str] | None = None,
    min_protein_length: int = 0,
    include_partial: bool = False,
    genetic_code: Mapping[str, str] | None = None,
) -> tuple[ORF, ...]:

    if isinstance(min_protein_length, bool) or not isinstance(min_protein_length, int) or min_protein_length < 0:
        raise AnnotationError("min_protein_length must be a non-negative integer")
    normalised = normalize_dna(sequence, allow_ambiguous=True)
    is_circular = _resolve_circular(topology, circular)
    code = _normalise_genetic_code(genetic_code)
    starts = _normalise_triplets(start_codons, label="start_codons")
    stops = (
        _normalise_triplets(stop_codons, label="stop_codons")
        if stop_codons is not None
        else frozenset(codon for codon, amino_acid in code.items() if amino_acid == "*")
    )
    if not stops:
        raise AnnotationError("genetic_code does not define stop codons; provide stop_codons")
    if len(normalised) < 3:
        return ()

    forward = _scan_oriented_orfs(
        normalised,
        strand=1,
        circular=is_circular,
        start_codons=starts,
        stop_codons=stops,
        min_protein_length=min_protein_length,
        include_partial=include_partial,
        genetic_code=code,
    )
    reverse = _scan_oriented_orfs(
        reverse_complement(normalised),
        strand=-1,
        circular=is_circular,
        start_codons=starts,
        stop_codons=stops,
        min_protein_length=min_protein_length,
        include_partial=include_partial,
        genetic_code=code,
    )
    return tuple(
        sorted(
            (*forward, *reverse),
            key=lambda orf: (
                0 if orf.strand == 1 else 1,
                orf.start,
                orf.end,
                orf.frame,
                orf.nucleotide_length,
            ),
        )
    )


def _orf_name(orf: ORF) -> str:
    if orf.wraps_origin:
        span = f"{orf.start + 1}..origin..{orf.end}"
    else:
        span = f"{orf.start + 1}..{orf.end}"
    return f"ORF {orf.reading_frame:+d} {span}"


def detect_orf_features(sequence: str, **kwargs: Any) -> FeatureDetectionResult:

    normalised = normalize_dna(sequence, allow_ambiguous=True)
    orfs = find_orfs(normalised, **kwargs)
    features: list[DetectedFeature] = []
    for orf in orfs:
        terminal = "terminal stop codon found" if orf.has_terminal_stop else "partial ORF without terminal stop"
        qualifiers: dict[str, Any] = {
            "translation": orf.protein_sequence,
            "reading_frame": orf.reading_frame,
            "start_codon": orf.start_codon,
            "terminal_stop_codon": orf.stop_codon,
            "partial": not orf.has_terminal_stop,
        }
        features.append(
            DetectedFeature(
                feature_type="CDS",
                name=_orf_name(orf),
                location=orf.location,
                qualifiers=qualifiers,
                confidence=1.0 if orf.has_terminal_stop else 0.55,
                evidence=("recognised in-frame start codon", terminal),
                source="opengene_annotations.orf_scan",
            )
        )
    return FeatureDetectionResult(
        sequence_length=len(normalised),
        features=tuple(features),
        orfs=orfs,
    )


def detect_features(sequence: str, **kwargs: Any) -> FeatureDetectionResult:

    return detect_orf_features(sequence, **kwargs)


def _normalize_protein(protein: str) -> str:
    if not isinstance(protein, str):
        raise AnnotationError("protein must be a string")
    normalised = "".join(protein.split()).upper()
    unsupported = "".join(sorted(set(normalised).difference(GENERIC_CODON_USAGE)))
    if unsupported:
        raise AnnotationError(
            "protein contains amino acids without deterministic DNA codons: " + unsupported
        )
    return normalised


CodonUsageInput = Mapping[
    str,
    Mapping[str, float] | Iterable[CodonPreference | tuple[str, float] | str],
]


def _normalise_codon_usage(codon_usage: CodonUsageInput | None) -> Mapping[str, tuple[CodonPreference, ...]]:
    if codon_usage is None or codon_usage is GENERIC_CODON_USAGE:
        return GENERIC_CODON_USAGE
    normalised: dict[str, tuple[CodonPreference, ...]] = {}
    for raw_amino_acid, raw_preferences in codon_usage.items():
        amino_acid = str(raw_amino_acid).strip().upper()
        if len(amino_acid) != 1 or amino_acid not in GENERIC_CODON_USAGE:
            raise AnnotationError("codon usage keys must be standard one-letter amino acids or '*'")
        if isinstance(raw_preferences, Mapping):
            preferences = tuple(
                CodonPreference(str(codon), weight)
                for codon, weight in raw_preferences.items()
            )
        else:
            collected: list[CodonPreference] = []
            for preference in raw_preferences:
                if isinstance(preference, CodonPreference):
                    collected.append(preference)
                elif isinstance(preference, str):
                    collected.append(CodonPreference(preference, 1.0))
                else:
                    try:
                        codon, weight = preference
                    except (TypeError, ValueError) as error:
                        raise AnnotationError("codon usage values must contain codons or (codon, weight) pairs") from error
                    collected.append(CodonPreference(str(codon), weight))
            preferences = tuple(collected)
        if not preferences:
            raise AnnotationError(f"codon usage for {amino_acid} cannot be empty")
        if len({preference.codon for preference in preferences}) != len(preferences):
            raise AnnotationError(f"codon usage for {amino_acid} contains duplicate codons")
        for preference in preferences:
            if STANDARD_GENETIC_CODE.get(preference.codon) != amino_acid:
                raise AnnotationError(
                    f"codon {preference.codon} does not encode amino acid {amino_acid} in the standard genetic code"
                )
        normalised[amino_acid] = tuple(
            sorted(preferences, key=lambda preference: (-preference.weight, preference.codon))
        )
    return MappingProxyType(normalised)


def _normalise_forbidden_motifs(forbidden_motifs: Iterable[str]) -> tuple[str, ...]:
    if isinstance(forbidden_motifs, str):
        forbidden_motifs = (forbidden_motifs,)
    motifs: list[str] = []
    for motif in forbidden_motifs:
        normalised = normalize_dna(motif, allow_ambiguous=False)
        if not normalised:
            raise AnnotationError("forbidden motifs cannot be empty")
        motifs.append(normalised)
    return tuple(dict.fromkeys(motifs))


def optimize_codon_usage(
    protein: str,
    *,
    codon_usage: CodonUsageInput | None = None,
    forbidden_motifs: Iterable[str] = (),
    usage_table_name: str = "generic-balanced-v1",
) -> CodonOptimizationResult:

    amino_acids = _normalize_protein(protein)
    table = _normalise_codon_usage(codon_usage)
    motifs = _normalise_forbidden_motifs(forbidden_motifs)
    chosen: list[str] = []
    warnings: list[str] = []
    current = ""

    for position, amino_acid in enumerate(amino_acids, start=1):
        preferences = table.get(amino_acid)
        if not preferences:
            raise AnnotationError(
                f"codon usage table has no candidate codon for amino acid {amino_acid} at position {position}"
            )
        selected: CodonPreference | None = None
        for preference in preferences:
            candidate = current + preference.codon
            if not any(candidate.endswith(motif) for motif in motifs):
                selected = preference
                break
        if selected is None:
            selected = preferences[0]
            violating = tuple(motif for motif in motifs if (current + selected.codon).endswith(motif))
            warnings.append(
                f"position {position} ({amino_acid}) could not avoid forbidden motif(s): {', '.join(violating)}"
            )
        chosen.append(selected.codon)
        current += selected.codon

    return CodonOptimizationResult(
        protein_sequence=amino_acids,
        dna_sequence=current,
        codons=tuple(chosen),
        usage_table_name=usage_table_name,
        warnings=tuple(warnings),
    )


def reverse_translate(
    protein: str,
    *,
    codon_usage: CodonUsageInput | None = None,
    forbidden_motifs: Iterable[str] = (),
) -> str:

    return optimize_codon_usage(
        protein,
        codon_usage=codon_usage,
        forbidden_motifs=forbidden_motifs,
    ).dna_sequence


find_open_reading_frames = find_orfs
codon_optimize = optimize_codon_usage
