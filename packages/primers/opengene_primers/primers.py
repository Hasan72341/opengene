from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


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


def normalize_dna(sequence: str, *, allow_ambiguous: bool = True) -> str:

    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    normalised = "".join(sequence.split()).upper().replace("U", "T")
    permitted = _IUPAC if allow_ambiguous else {base: _IUPAC[base] for base in "ACGT"}
    invalid = sorted(set(normalised).difference(permitted))
    if invalid:
        raise ValueError(f"sequence contains unsupported DNA symbols: {''.join(invalid)}")
    return normalised


def reverse_complement(sequence: str) -> str:

    return "".join(_COMPLEMENT[base] for base in reversed(normalize_dna(sequence)))


def _expected_gc(base: str) -> float:

    bases = _IUPAC[base]
    return len(bases.intersection({"G", "C"})) / len(bases)


def gc_content(sequence: str) -> float:

    normalised = normalize_dna(sequence)
    if not normalised:
        raise ValueError("GC content is undefined for an empty sequence")
    return 100.0 * sum(_expected_gc(base) for base in normalised) / len(normalised)


def melting_temperature(sequence: str, *, method: str = "wallace") -> float:

    normalised = normalize_dna(sequence)
    if not normalised:
        raise ValueError("melting temperature is undefined for an empty sequence")
    expected_gc = sum(_expected_gc(base) for base in normalised)
    expected_at = len(normalised) - expected_gc
    normalized_method = method.strip().lower()
    if normalized_method == "wallace":
        return 2.0 * expected_at + 4.0 * expected_gc
    if normalized_method in {"basic", "long"}:
        if len(normalised) < 2:
            raise ValueError("basic Tm approximation requires at least two bases")
        return 64.9 + 41.0 * (expected_gc - 16.4) / len(normalised)
    raise ValueError("method must be 'wallace' or 'basic'")


def _has_long_homopolymer(sequence: str, limit: int = 4) -> bool:
    return any(base * limit in sequence for base in "ACGT")


def _is_self_complementary(sequence: str, min_overlap: int = 4) -> bool:

    reverse = reverse_complement(sequence)
    for length in range(min_overlap, len(sequence) + 1):
        for start in range(len(sequence) - length + 1):
            if sequence[start : start + length] in reverse:
                return True
    return False


@dataclass(frozen=True, slots=True)
class PrimerBindingSite:

    start: int
    end: int
    strand: int
    mismatches: int = 0
    wraps_origin: bool = False


@dataclass(frozen=True, slots=True)
class Primer:

    name: str
    annealing_sequence: str
    five_prime_extension: str = ""
    binding_site: PrimerBindingSite | None = None

    def __post_init__(self) -> None:
        if not self.name or not self.name.strip():
            raise ValueError("primer name cannot be empty")
        annealing = normalize_dna(self.annealing_sequence, allow_ambiguous=False)
        extension = normalize_dna(self.five_prime_extension, allow_ambiguous=False)
        if not annealing:
            raise ValueError("annealing_sequence cannot be empty")
        object.__setattr__(self, "name", self.name.strip())
        object.__setattr__(self, "annealing_sequence", annealing)
        object.__setattr__(self, "five_prime_extension", extension)

    @property
    def full_sequence(self) -> str:
        return self.five_prime_extension + self.annealing_sequence

    @property
    def gc_percent(self) -> float:
        return gc_content(self.annealing_sequence)

    @property
    def tm_celsius(self) -> float:
        return melting_temperature(self.annealing_sequence)


@dataclass(frozen=True, slots=True)
class PrimerPair:

    forward: Primer
    reverse: Primer
    product_size: int
    score: float


@dataclass(frozen=True, slots=True)
class PCRProduct:

    sequence: str
    forward_site: PrimerBindingSite
    reverse_site: PrimerBindingSite
    circular_template: bool = False

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True, slots=True)
class PCRSimulation:

    template_length: int
    products: tuple[PCRProduct, ...]
    warnings: tuple[str, ...] = ()


def _matches_base(template_base: str, primer_base: str) -> bool:
    return bool(_IUPAC[template_base] & _IUPAC[primer_base])


def _mismatch_count(window: str, primer: str) -> int:
    return sum(not _matches_base(template_base, primer_base) for template_base, primer_base in zip(window, primer))


def _window(template: str, start: int, length: int, circular: bool) -> str:
    if not circular:
        return template[start : start + length]
    return "".join(template[(start + offset) % len(template)] for offset in range(length))


def find_primer_binding_sites(
    template: str,
    primer: Primer | str,
    *,
    circular: bool = False,
    max_mismatches: int = 0,
) -> tuple[PrimerBindingSite, ...]:

    if max_mismatches < 0:
        raise ValueError("max_mismatches cannot be negative")
    normalised_template = normalize_dna(template)
    annealing = primer.annealing_sequence if isinstance(primer, Primer) else normalize_dna(primer, allow_ambiguous=False)
    if not normalised_template:
        return ()
    length = len(annealing)
    if not circular and length > len(normalised_template):
        return ()
    starts = range(len(normalised_template)) if circular else range(len(normalised_template) - length + 1)
    reverse = reverse_complement(annealing)
    sites: list[PrimerBindingSite] = []
    for start in starts:
        observed = _window(normalised_template, start, length, circular)
        raw_end = start + length
        end = raw_end if not circular or raw_end <= len(normalised_template) else raw_end % len(normalised_template)
        forward_mismatches = _mismatch_count(observed, annealing)
        if forward_mismatches <= max_mismatches:
            sites.append(
                PrimerBindingSite(
                    start=start,
                    end=end,
                    strand=1,
                    mismatches=forward_mismatches,
                    wraps_origin=circular and raw_end > len(normalised_template),
                )
            )
        reverse_mismatches = _mismatch_count(observed, reverse)
        if reverse_mismatches <= max_mismatches:
            sites.append(
                PrimerBindingSite(
                    start=start,
                    end=end,
                    strand=-1,
                    mismatches=reverse_mismatches,
                    wraps_origin=circular and raw_end > len(normalised_template),
                )
            )
    return tuple(sorted(sites, key=lambda site: (site.start, site.strand, site.mismatches)))


def _candidate_score(sequence: str, *, tm_target: float, gc_target: float) -> float:

    score = abs(melting_temperature(sequence) - tm_target)
    score += abs(gc_content(sequence) - gc_target) / 10.0
    if sequence[-1] not in "GC":
        score += 1.0
    if _has_long_homopolymer(sequence):
        score += 5.0
    if _is_self_complementary(sequence):
        score += 3.0
    return score


def _primer_candidates(
    template: str,
    *,
    starts: Iterable[int],
    lengths: range,
    reverse: bool,
    tm_range: tuple[float, float],
    gc_range: tuple[float, float],
    extension: str,
    name_prefix: str,
) -> list[tuple[Primer, int, int, float]]:
    candidates: list[tuple[Primer, int, int, float]] = []
    target_tm = (tm_range[0] + tm_range[1]) / 2.0
    target_gc = (gc_range[0] + gc_range[1]) / 2.0
    for anchor in starts:
        for length in lengths:
            start, end = (anchor - length, anchor) if reverse else (anchor, anchor + length)
            if start < 0 or end > len(template):
                continue
            binding_sequence = template[start:end]
            if set(binding_sequence).difference("ACGT"):
                continue
            annealing = reverse_complement(binding_sequence) if reverse else binding_sequence
            tm = melting_temperature(annealing)
            gc = gc_content(annealing)
            if not tm_range[0] <= tm <= tm_range[1] or not gc_range[0] <= gc <= gc_range[1]:
                continue
            strand = -1 if reverse else 1
            site = PrimerBindingSite(start, end, strand)
            primer = Primer(
                name=f"{name_prefix}_{'R' if reverse else 'F'}_{start + 1}_{length}",
                annealing_sequence=annealing,
                five_prime_extension=extension,
                binding_site=site,
            )
            candidates.append((primer, start, end, _candidate_score(annealing, tm_target=target_tm, gc_target=target_gc)))
    return candidates


def design_primers(
    template: str,
    *,
    target_start: int = 0,
    target_end: int | None = None,
    product_size_range: tuple[int, int] = (70, 10_000),
    primer_length_range: tuple[int, int] = (18, 25),
    tm_range: tuple[float, float] = (52.0, 68.0),
    gc_range: tuple[float, float] = (35.0, 65.0),
    search_window: int = 80,
    forward_extension: str = "",
    reverse_extension: str = "",
    max_pairs: int = 10,
    name_prefix: str = "PCR",
) -> tuple[PrimerPair, ...]:

    normalised = normalize_dna(template)
    if not normalised:
        raise ValueError("template cannot be empty")
    target_end = len(normalised) if target_end is None else target_end
    if not (0 <= target_start < target_end <= len(normalised)):
        raise ValueError("target interval must lie within the template")
    min_product, max_product = product_size_range
    min_length, max_length = primer_length_range
    if min_product <= 0 or min_product > max_product:
        raise ValueError("product_size_range must be ordered positive values")
    if min_length <= 0 or min_length > max_length:
        raise ValueError("primer_length_range must be ordered positive values")
    if search_window < 0 or max_pairs < 0:
        raise ValueError("search_window and max_pairs cannot be negative")
    if tm_range[0] > tm_range[1] or gc_range[0] > gc_range[1]:
        raise ValueError("Tm and GC ranges must be ordered")

    lengths = range(min_length, max_length + 1)
    forward_starts = range(max(0, target_start - search_window), target_start + 1)
    reverse_ends = range(target_end, min(len(normalised), target_end + search_window) + 1)
    forward_candidates = _primer_candidates(
        normalised,
        starts=forward_starts,
        lengths=lengths,
        reverse=False,
        tm_range=tm_range,
        gc_range=gc_range,
        extension=forward_extension,
        name_prefix=name_prefix,
    )
    reverse_candidates = _primer_candidates(
        normalised,
        starts=reverse_ends,
        lengths=lengths,
        reverse=True,
        tm_range=tm_range,
        gc_range=gc_range,
        extension=reverse_extension,
        name_prefix=name_prefix,
    )

    pairs: list[PrimerPair] = []
    for forward, f_start, _f_end, f_score in forward_candidates:
        for reverse, _r_start, r_end, r_score in reverse_candidates:
            product_size = r_end - f_start
            if not min_product <= product_size <= max_product:
                continue
            score = f_score + r_score + abs(forward.tm_celsius - reverse.tm_celsius) / 2.0
            pairs.append(PrimerPair(forward, reverse, product_size, score))
    pairs.sort(key=lambda pair: (pair.score, pair.product_size, pair.forward.name, pair.reverse.name))
    return tuple(pairs[:max_pairs])


def _positive_sites(template: str, primer: Primer, *, circular: bool, mismatches: int) -> tuple[PrimerBindingSite, ...]:
    return tuple(
        site
        for site in find_primer_binding_sites(template, primer, circular=circular, max_mismatches=mismatches)
        if site.strand == 1
    )


def _negative_sites(template: str, primer: Primer, *, circular: bool, mismatches: int) -> tuple[PrimerBindingSite, ...]:
    return tuple(
        site
        for site in find_primer_binding_sites(template, primer, circular=circular, max_mismatches=mismatches)
        if site.strand == -1
    )


def _circular_between(sequence: str, start: int, end: int) -> str:
    if start == end:
        return ""
    return sequence[start:end] if start < end else sequence[start:] + sequence[:end]


def simulate_pcr(
    template: str,
    forward_primer: Primer | str,
    reverse_primer: Primer | str,
    *,
    circular: bool = False,
    max_mismatches: int = 0,
) -> PCRSimulation:

    normalised = normalize_dna(template)
    if not normalised:
        raise ValueError("template cannot be empty")
    forward = forward_primer if isinstance(forward_primer, Primer) else Primer("Forward", forward_primer)
    reverse = reverse_primer if isinstance(reverse_primer, Primer) else Primer("Reverse", reverse_primer)
    forward_sites = _positive_sites(normalised, forward, circular=circular, mismatches=max_mismatches)
    reverse_sites = _negative_sites(normalised, reverse, circular=circular, mismatches=max_mismatches)
    products: list[PCRProduct] = []
    seen: set[tuple[int, int]] = set()
    for forward_site in forward_sites:
        for reverse_site in reverse_sites:
            if not circular:
                if forward_site.start >= reverse_site.end:
                    continue
                interior = normalised[forward_site.end : reverse_site.start]
            else:
                if forward_site.start == reverse_site.start and forward_site.end == reverse_site.end:
                    continue
                forward_end = (forward_site.start + len(forward.annealing_sequence)) % len(normalised)
                interior = _circular_between(normalised, forward_end, reverse_site.start)
            key = (forward_site.start, reverse_site.start)
            if key in seen:
                continue
            seen.add(key)
            sequence = forward.full_sequence + interior + reverse_complement(reverse.full_sequence)
            products.append(PCRProduct(sequence, forward_site, reverse_site, circular_template=circular))
    products.sort(key=lambda product: (product.length, product.forward_site.start, product.reverse_site.start))
    warnings: list[str] = []
    if not products:
        warnings.append("No oriented primer pair produced an amplicon")
    elif len(products) > 1:
        warnings.append(f"{len(products)} possible products found; inspect primer specificity")
    if max_mismatches:
        warnings.append("Products include bindings allowed by max_mismatches")
    return PCRSimulation(len(normalised), tuple(products), tuple(warnings))
