from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Literal


class AlignmentError(ValueError):
    pass


class AlignmentMode(str, Enum):

    GLOBAL = "global"
    LOCAL = "local"


class VariantKind(str, Enum):

    SUBSTITUTION = "substitution"
    INSERTION = "insertion"
    DELETION = "deletion"


_DNA_IUPAC: dict[str, frozenset[str]] = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "U": frozenset("U"),
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


def _normalise_sequence(sequence: str, *, label: str) -> str:

    if not isinstance(sequence, str):
        raise AlignmentError(f"{label} must be a string")
    normalised = "".join(sequence.split()).upper()
    if "-" in normalised:
        raise AlignmentError(f"{label} must not contain gap characters")
    if any(base not in "ABCDEFGHIJKLMNOPQRSTUVWXYZ*" for base in normalised):
        raise AlignmentError(f"{label} contains unsupported characters")
    return normalised


def _bases_match(reference_base: str, query_base: str, *, ambiguity_matches: bool) -> bool:

    if reference_base == query_base:
        return True
    if not ambiguity_matches:
        return False
    reference_options = _DNA_IUPAC.get(reference_base)
    query_options = _DNA_IUPAC.get(query_base)
    return bool(reference_options and query_options and reference_options & query_options)


@dataclass(frozen=True, slots=True)
class AlignmentVariant:

    kind: VariantKind
    reference_start: int
    reference_end: int
    reference_bases: str
    alternate_bases: str
    query_start: int
    query_end: int

    def __post_init__(self) -> None:
        try:
            kind = VariantKind(self.kind)
        except ValueError as error:
            raise AlignmentError(f"unknown variant kind: {self.kind!r}") from error
        object.__setattr__(self, "kind", kind)
        coordinate_names = (
            "reference_start",
            "reference_end",
            "query_start",
            "query_end",
        )
        for name in coordinate_names:
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise AlignmentError(f"{name} must be a non-negative integer")
        if self.reference_end < self.reference_start or self.query_end < self.query_start:
            raise AlignmentError("variant end coordinates cannot precede start coordinates")
        if kind is VariantKind.SUBSTITUTION:
            if not self.reference_bases or not self.alternate_bases:
                raise AlignmentError("a substitution requires reference and alternate bases")
            if self.reference_end - self.reference_start != len(self.reference_bases):
                raise AlignmentError("substitution reference coordinates do not match bases")
            if self.query_end - self.query_start != len(self.alternate_bases):
                raise AlignmentError("substitution query coordinates do not match bases")
        elif kind is VariantKind.INSERTION:
            if self.reference_bases or not self.alternate_bases:
                raise AlignmentError("an insertion requires only alternate bases")
            if self.reference_start != self.reference_end:
                raise AlignmentError("an insertion must have an empty reference interval")
            if self.query_end - self.query_start != len(self.alternate_bases):
                raise AlignmentError("insertion query coordinates do not match bases")
        else:
            if not self.reference_bases or self.alternate_bases:
                raise AlignmentError("a deletion requires only reference bases")
            if self.reference_end - self.reference_start != len(self.reference_bases):
                raise AlignmentError("deletion reference coordinates do not match bases")
            if self.query_start != self.query_end:
                raise AlignmentError("a deletion must have an empty query interval")


@dataclass(frozen=True, slots=True)
class AlignmentResult:

    mode: AlignmentMode
    score: float
    aligned_reference: str
    aligned_query: str
    reference_start: int
    reference_end: int
    query_start: int
    query_end: int
    matches: int
    mismatches: int
    gap_columns: int
    ambiguity_matches: bool = False

    def __post_init__(self) -> None:
        try:
            mode = AlignmentMode(self.mode)
        except ValueError as error:
            raise AlignmentError(f"unknown alignment mode: {self.mode!r}") from error
        object.__setattr__(self, "mode", mode)
        if len(self.aligned_reference) != len(self.aligned_query):
            raise AlignmentError("aligned sequences must have equal lengths")
        if any(
            reference_base == "-" and query_base == "-"
            for reference_base, query_base in zip(self.aligned_reference, self.aligned_query)
        ):
            raise AlignmentError("an alignment cannot contain a double-gap column")
        coordinate_names = (
            "reference_start",
            "reference_end",
            "query_start",
            "query_end",
            "matches",
            "mismatches",
            "gap_columns",
        )
        for name in coordinate_names:
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise AlignmentError(f"{name} must be a non-negative integer")
        if self.reference_end < self.reference_start or self.query_end < self.query_start:
            raise AlignmentError("alignment end coordinates cannot precede start coordinates")
        reference_bases = sum(base != "-" for base in self.aligned_reference)
        query_bases = sum(base != "-" for base in self.aligned_query)
        if reference_bases != self.reference_end - self.reference_start:
            raise AlignmentError("reference coordinates do not match aligned reference bases")
        if query_bases != self.query_end - self.query_start:
            raise AlignmentError("query coordinates do not match aligned query bases")
        matches, mismatches, gap_columns = _alignment_statistics(
            self.aligned_reference,
            self.aligned_query,
            ambiguity_matches=self.ambiguity_matches,
        )
        if (self.matches, self.mismatches, self.gap_columns) != (matches, mismatches, gap_columns):
            raise AlignmentError("alignment statistics do not match aligned sequences")

    @property
    def alignment_length(self) -> int:

        return len(self.aligned_reference)

    @property
    def identity(self) -> float:

        paired = self.matches + self.mismatches
        return self.matches / paired if paired else 0.0

    @property
    def variants(self) -> tuple[AlignmentVariant, ...]:

        return variants_from_alignment(self)

    @property
    def cigar(self) -> str:

        if not self.aligned_reference:
            return ""
        operations: list[str] = []
        run_operation: str | None = None
        run_length = 0
        for reference_base, query_base in zip(self.aligned_reference, self.aligned_query):
            if reference_base == "-":
                operation = "I"
            elif query_base == "-":
                operation = "D"
            elif _bases_match(
                reference_base,
                query_base,
                ambiguity_matches=self.ambiguity_matches,
            ):
                operation = "="
            else:
                operation = "X"
            if operation == run_operation:
                run_length += 1
            else:
                if run_operation is not None:
                    operations.append(f"{run_length}{run_operation}")
                run_operation = operation
                run_length = 1
        if run_operation is not None:
            operations.append(f"{run_length}{run_operation}")
        return "".join(operations)


PairwiseAlignment = AlignmentResult


def _alignment_statistics(
    aligned_reference: str,
    aligned_query: str,
    *,
    ambiguity_matches: bool,
) -> tuple[int, int, int]:
    matches = 0
    mismatches = 0
    gap_columns = 0
    for reference_base, query_base in zip(aligned_reference, aligned_query):
        if reference_base == "-" or query_base == "-":
            gap_columns += 1
        elif _bases_match(reference_base, query_base, ambiguity_matches=ambiguity_matches):
            matches += 1
        else:
            mismatches += 1
    return matches, mismatches, gap_columns


def _validate_scores(match_score: float, mismatch_score: float, gap_score: float) -> None:
    for name, value in (
        ("match_score", match_score),
        ("mismatch_score", mismatch_score),
        ("gap_score", gap_score),
    ):
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise AlignmentError(f"{name} must be numeric")


def _result(
    *,
    mode: AlignmentMode,
    score: float,
    aligned_reference: str,
    aligned_query: str,
    reference_start: int,
    reference_end: int,
    query_start: int,
    query_end: int,
    ambiguity_matches: bool,
) -> AlignmentResult:
    matches, mismatches, gap_columns = _alignment_statistics(
        aligned_reference,
        aligned_query,
        ambiguity_matches=ambiguity_matches,
    )
    return AlignmentResult(
        mode=mode,
        score=score,
        aligned_reference=aligned_reference,
        aligned_query=aligned_query,
        reference_start=reference_start,
        reference_end=reference_end,
        query_start=query_start,
        query_end=query_end,
        matches=matches,
        mismatches=mismatches,
        gap_columns=gap_columns,
        ambiguity_matches=ambiguity_matches,
    )


def needleman_wunsch(
    reference: str,
    query: str,
    *,
    match_score: float = 2,
    mismatch_score: float = -1,
    gap_score: float = -2,
    ambiguity_matches: bool = False,
) -> AlignmentResult:

    reference = _normalise_sequence(reference, label="reference")
    query = _normalise_sequence(query, label="query")
    _validate_scores(match_score, mismatch_score, gap_score)
    if not isinstance(ambiguity_matches, bool):
        raise AlignmentError("ambiguity_matches must be a boolean")

    reference_length = len(reference)
    query_length = len(query)
    scores = [[0.0] * (query_length + 1) for _ in range(reference_length + 1)]
    trace: list[list[Literal["diag", "up", "left"] | None]] = [
        [None] * (query_length + 1) for _ in range(reference_length + 1)
    ]
    for reference_index in range(1, reference_length + 1):
        scores[reference_index][0] = scores[reference_index - 1][0] + gap_score
        trace[reference_index][0] = "up"
    for query_index in range(1, query_length + 1):
        scores[0][query_index] = scores[0][query_index - 1] + gap_score
        trace[0][query_index] = "left"

    for reference_index in range(1, reference_length + 1):
        reference_base = reference[reference_index - 1]
        for query_index in range(1, query_length + 1):
            query_base = query[query_index - 1]
            diagonal = scores[reference_index - 1][query_index - 1] + (
                match_score
                if _bases_match(reference_base, query_base, ambiguity_matches=ambiguity_matches)
                else mismatch_score
            )
            deletion = scores[reference_index - 1][query_index] + gap_score
            insertion = scores[reference_index][query_index - 1] + gap_score
            best = max(diagonal, deletion, insertion)
            scores[reference_index][query_index] = best
            if best == diagonal:
                trace[reference_index][query_index] = "diag"
            elif best == deletion:
                trace[reference_index][query_index] = "up"
            else:
                trace[reference_index][query_index] = "left"

    reference_index = reference_length
    query_index = query_length
    reverse_reference: list[str] = []
    reverse_query: list[str] = []
    while reference_index > 0 or query_index > 0:
        direction = trace[reference_index][query_index]
        if direction == "diag":
            reverse_reference.append(reference[reference_index - 1])
            reverse_query.append(query[query_index - 1])
            reference_index -= 1
            query_index -= 1
        elif direction == "up":
            reverse_reference.append(reference[reference_index - 1])
            reverse_query.append("-")
            reference_index -= 1
        elif direction == "left":
            reverse_reference.append("-")
            reverse_query.append(query[query_index - 1])
            query_index -= 1
        else:
            raise AlignmentError("incomplete global alignment traceback")
    aligned_reference = "".join(reversed(reverse_reference))
    aligned_query = "".join(reversed(reverse_query))
    return _result(
        mode=AlignmentMode.GLOBAL,
        score=scores[reference_length][query_length],
        aligned_reference=aligned_reference,
        aligned_query=aligned_query,
        reference_start=0,
        reference_end=reference_length,
        query_start=0,
        query_end=query_length,
        ambiguity_matches=ambiguity_matches,
    )


def smith_waterman(
    reference: str,
    query: str,
    *,
    match_score: float = 2,
    mismatch_score: float = -1,
    gap_score: float = -2,
    ambiguity_matches: bool = False,
) -> AlignmentResult:

    reference = _normalise_sequence(reference, label="reference")
    query = _normalise_sequence(query, label="query")
    _validate_scores(match_score, mismatch_score, gap_score)
    if gap_score > 0:
        raise AlignmentError("Smith-Waterman requires a non-positive gap_score")
    if not isinstance(ambiguity_matches, bool):
        raise AlignmentError("ambiguity_matches must be a boolean")

    reference_length = len(reference)
    query_length = len(query)
    scores = [[0.0] * (query_length + 1) for _ in range(reference_length + 1)]
    trace: list[list[Literal["diag", "up", "left"] | None]] = [
        [None] * (query_length + 1) for _ in range(reference_length + 1)
    ]
    best_score = 0.0
    best_reference_index = 0
    best_query_index = 0

    for reference_index in range(1, reference_length + 1):
        reference_base = reference[reference_index - 1]
        for query_index in range(1, query_length + 1):
            query_base = query[query_index - 1]
            diagonal = scores[reference_index - 1][query_index - 1] + (
                match_score
                if _bases_match(reference_base, query_base, ambiguity_matches=ambiguity_matches)
                else mismatch_score
            )
            deletion = scores[reference_index - 1][query_index] + gap_score
            insertion = scores[reference_index][query_index - 1] + gap_score
            best = max(0.0, diagonal, deletion, insertion)
            scores[reference_index][query_index] = best
            if best == 0.0:
                trace[reference_index][query_index] = None
            elif best == diagonal:
                trace[reference_index][query_index] = "diag"
            elif best == deletion:
                trace[reference_index][query_index] = "up"
            else:
                trace[reference_index][query_index] = "left"
            if best > best_score:
                best_score = best
                best_reference_index = reference_index
                best_query_index = query_index

    reference_end = best_reference_index
    query_end = best_query_index
    reference_index = best_reference_index
    query_index = best_query_index
    reverse_reference: list[str] = []
    reverse_query: list[str] = []
    while reference_index > 0 and query_index > 0 and scores[reference_index][query_index] > 0:
        direction = trace[reference_index][query_index]
        if direction == "diag":
            reverse_reference.append(reference[reference_index - 1])
            reverse_query.append(query[query_index - 1])
            reference_index -= 1
            query_index -= 1
        elif direction == "up":
            reverse_reference.append(reference[reference_index - 1])
            reverse_query.append("-")
            reference_index -= 1
        elif direction == "left":
            reverse_reference.append("-")
            reverse_query.append(query[query_index - 1])
            query_index -= 1
        else:
            break
    aligned_reference = "".join(reversed(reverse_reference))
    aligned_query = "".join(reversed(reverse_query))
    return _result(
        mode=AlignmentMode.LOCAL,
        score=best_score,
        aligned_reference=aligned_reference,
        aligned_query=aligned_query,
        reference_start=reference_index,
        reference_end=reference_end,
        query_start=query_index,
        query_end=query_end,
        ambiguity_matches=ambiguity_matches,
    )


def align_global(*args: object, **kwargs: object) -> AlignmentResult:

    return needleman_wunsch(*args, **kwargs)  # type: ignore[arg-type]


def align_local(*args: object, **kwargs: object) -> AlignmentResult:

    return smith_waterman(*args, **kwargs)  # type: ignore[arg-type]


def variants_from_alignment(alignment: AlignmentResult) -> tuple[AlignmentVariant, ...]:

    if not isinstance(alignment, AlignmentResult):
        raise TypeError("alignment must be an AlignmentResult")

    variants: list[AlignmentVariant] = []
    current_kind: VariantKind | None = None
    reference_start = 0
    query_start = 0
    reference_bases: list[str] = []
    alternate_bases: list[str] = []

    def flush() -> None:
        nonlocal current_kind, reference_bases, alternate_bases
        if current_kind is None:
            return
        variants.append(
            AlignmentVariant(
                kind=current_kind,
                reference_start=reference_start,
                reference_end=reference_start + len(reference_bases),
                reference_bases="".join(reference_bases),
                alternate_bases="".join(alternate_bases),
                query_start=query_start,
                query_end=query_start + len(alternate_bases),
            )
        )
        current_kind = None
        reference_bases = []
        alternate_bases = []

    reference_position = alignment.reference_start
    query_position = alignment.query_start
    for reference_base, query_base in zip(alignment.aligned_reference, alignment.aligned_query):
        if reference_base == "-":
            kind = VariantKind.INSERTION
        elif query_base == "-":
            kind = VariantKind.DELETION
        elif _bases_match(
            reference_base,
            query_base,
            ambiguity_matches=alignment.ambiguity_matches,
        ):
            flush()
            reference_position += 1
            query_position += 1
            continue
        else:
            kind = VariantKind.SUBSTITUTION

        can_extend = current_kind is kind
        if can_extend:
            expected_reference = reference_start + len(reference_bases)
            expected_query = query_start + len(alternate_bases)
            can_extend = reference_position == expected_reference and query_position == expected_query
        if not can_extend:
            flush()
            current_kind = kind
            reference_start = reference_position
            query_start = query_position
        if reference_base != "-":
            reference_bases.append(reference_base)
            reference_position += 1
        if query_base != "-":
            alternate_bases.append(query_base)
            query_position += 1
    flush()
    return tuple(variants)
