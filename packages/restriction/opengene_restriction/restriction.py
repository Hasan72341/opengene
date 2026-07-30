from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
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


class Topology(str, Enum):

    LINEAR = "linear"
    CIRCULAR = "circular"


def normalize_dna(sequence: str) -> str:

    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    normalised = "".join(sequence.split()).upper().replace("U", "T")
    invalid = sorted(set(normalised).difference(_IUPAC))
    if invalid:
        raise ValueError(f"sequence contains unsupported DNA symbols: {''.join(invalid)}")
    return normalised


def reverse_complement(sequence: str) -> str:

    return "".join(_COMPLEMENT[base] for base in reversed(normalize_dna(sequence)))


def _matches(window: str, motif: str) -> bool:

    return all(_IUPAC[observed] & _IUPAC[expected] for observed, expected in zip(window, motif))


def _normalise_topology(topology: Topology | str | bool) -> Topology:
    if isinstance(topology, bool):
        return Topology.CIRCULAR if topology else Topology.LINEAR
    try:
        return Topology(topology)
    except ValueError as exc:
        raise ValueError("topology must be 'linear' or 'circular'") from exc


@dataclass(frozen=True, slots=True)
class RestrictionEnzyme:

    name: str
    recognition_site: str
    cut_forward: int
    cut_reverse: int

    def __post_init__(self) -> None:
        if not self.name or not self.name.strip():
            raise ValueError("enzyme name cannot be empty")
        site = normalize_dna(self.recognition_site)
        if not site:
            raise ValueError("recognition_site cannot be empty")
        if not isinstance(self.cut_forward, int) or not isinstance(self.cut_reverse, int):
            raise TypeError("cut offsets must be integers")
        object.__setattr__(self, "recognition_site", site)
        object.__setattr__(self, "name", self.name.strip())

    @property
    def site_length(self) -> int:
        return len(self.recognition_site)

    @property
    def reverse_recognition_site(self) -> str:
        return reverse_complement(self.recognition_site)

    @property
    def is_palindromic(self) -> bool:
        return self.recognition_site == self.reverse_recognition_site

    @property
    def overhang_length(self) -> int:
        return abs(self.cut_forward - self.cut_reverse)

    @property
    def end_type(self) -> str:

        if self.cut_forward == self.cut_reverse:
            return "blunt"
        return "5-prime" if self.cut_forward < self.cut_reverse else "3-prime"


@dataclass(frozen=True, slots=True)
class RestrictionSite:

    enzyme: RestrictionEnzyme
    start: int
    end: int
    strand: int
    matched_sequence: str
    cut_forward: int
    cut_reverse: int
    wraps_origin: bool = False

    @property
    def cut_position(self) -> int:

        return self.cut_forward


@dataclass(frozen=True, slots=True)
class DigestFragment:

    start: int
    end: int
    sequence: str
    wraps_origin: bool = False
    is_circular: bool = False
    left_cut: int | None = None
    right_cut: int | None = None

    @property
    def length(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True, slots=True)
class DigestResult:

    sequence: str
    topology: Topology
    sites: tuple[RestrictionSite, ...]
    cut_positions: tuple[int, ...]
    fragments: tuple[DigestFragment, ...]
    warnings: tuple[str, ...] = ()


def _window(sequence: str, start: int, length: int, circular: bool) -> str:
    if not circular:
        return sequence[start : start + length]
    size = len(sequence)
    return "".join(sequence[(start + offset) % size] for offset in range(length))


def _cut_coordinate(value: int, length: int, circular: bool) -> int:
    return value % length if circular and length else value


def scan_restriction_sites(
    sequence: str,
    enzymes: Iterable[RestrictionEnzyme],
    *,
    topology: Topology | str | bool = Topology.LINEAR,
) -> tuple[RestrictionSite, ...]:

    sequence = normalize_dna(sequence)
    selected_topology = _normalise_topology(topology)
    circular = selected_topology is Topology.CIRCULAR
    enzyme_list = tuple(enzymes)
    for enzyme in enzyme_list:
        if not isinstance(enzyme, RestrictionEnzyme):
            raise TypeError("enzymes must contain RestrictionEnzyme objects")
    if not sequence:
        return ()

    sites: list[RestrictionSite] = []
    sequence_length = len(sequence)
    for enzyme in enzyme_list:
        oriented_motifs = [(enzyme.recognition_site, 1)]
        if not enzyme.is_palindromic:
            oriented_motifs.append((enzyme.reverse_recognition_site, -1))
        starts = range(sequence_length) if circular else range(max(0, sequence_length - enzyme.site_length + 1))
        for start in starts:
            observed = _window(sequence, start, enzyme.site_length, circular)
            for motif, strand in oriented_motifs:
                if not _matches(observed, motif):
                    continue
                if strand == 1:
                    forward_cut = start + enzyme.cut_forward
                    reverse_cut = start + enzyme.cut_reverse
                else:
                    forward_cut = start + enzyme.site_length - enzyme.cut_reverse
                    reverse_cut = start + enzyme.site_length - enzyme.cut_forward
                raw_end = start + enzyme.site_length
                sites.append(
                    RestrictionSite(
                        enzyme=enzyme,
                        start=start,
                        end=raw_end if not circular or raw_end <= sequence_length else raw_end % sequence_length,
                        strand=strand,
                        matched_sequence=observed,
                        cut_forward=_cut_coordinate(forward_cut, sequence_length, circular),
                        cut_reverse=_cut_coordinate(reverse_cut, sequence_length, circular),
                        wraps_origin=circular and raw_end > sequence_length,
                    )
                )
    return tuple(sorted(sites, key=lambda site: (site.start, site.enzyme.name.casefold(), site.strand)))


def _circular_slice(sequence: str, start: int, end: int) -> str:

    if start == end:
        return sequence[start:] + sequence[:start]
    if start < end:
        return sequence[start:end]
    return sequence[start:] + sequence[:end]


def digest(
    sequence: str,
    enzymes: Iterable[RestrictionEnzyme],
    *,
    topology: Topology | str | bool = Topology.LINEAR,
) -> DigestResult:

    normalised = normalize_dna(sequence)
    selected_topology = _normalise_topology(topology)
    circular = selected_topology is Topology.CIRCULAR
    sites = scan_restriction_sites(normalised, enzymes, topology=selected_topology)
    size = len(normalised)
    if not size:
        return DigestResult(normalised, selected_topology, sites, (), ())

    warnings: list[str] = []
    if circular:
        cuts = tuple(sorted({site.cut_forward % size for site in sites}))
    else:
        usable_sites: list[RestrictionSite] = []
        for site in sites:
            if 0 <= site.cut_forward <= size and 0 <= site.cut_reverse <= size:
                usable_sites.append(site)
            else:
                warnings.append(
                    f"{site.enzyme.name} at {site.start} cuts outside the linear sequence and was not fragmented"
                )
        cuts = tuple(sorted({site.cut_forward for site in usable_sites}))

    fragments: list[DigestFragment] = []
    if circular:
        if not cuts:
            fragments.append(DigestFragment(0, size, normalised, is_circular=True))
        else:
            for index, start in enumerate(cuts):
                end = cuts[(index + 1) % len(cuts)]
                fragments.append(
                    DigestFragment(
                        start=start,
                        end=end,
                        sequence=_circular_slice(normalised, start, end),
                        wraps_origin=end <= start,
                        left_cut=start,
                        right_cut=end,
                    )
                )
    else:
        interior_cuts = tuple(cut for cut in cuts if 0 < cut < size)
        boundaries = (0, *interior_cuts, size)
        for index, (start, end) in enumerate(zip(boundaries, boundaries[1:])):
            fragments.append(
                DigestFragment(
                    start=start,
                    end=end,
                    sequence=normalised[start:end],
                    left_cut=boundaries[index] if index else None,
                    right_cut=boundaries[index + 1] if index + 1 < len(boundaries) - 1 else None,
                )
            )
        if not fragments:
            fragments.append(DigestFragment(0, size, normalised))

    return DigestResult(normalised, selected_topology, sites, cuts, tuple(fragments), tuple(warnings))


def unique_cutters(
    sequence: str,
    enzymes: Iterable[RestrictionEnzyme],
    *,
    topology: Topology | str | bool = Topology.LINEAR,
) -> tuple[RestrictionEnzyme, ...]:

    selected = tuple(enzymes)
    return tuple(
        enzyme
        for enzyme in selected
        if len(digest(sequence, (enzyme,), topology=topology).cut_positions) == 1
    )


def find_unique_cutters(
    sequence: str,
    enzymes: Iterable[RestrictionEnzyme],
    *,
    topology: Topology | str | bool = Topology.LINEAR,
) -> tuple[RestrictionEnzyme, ...]:

    return unique_cutters(sequence, enzymes, topology=topology)


COMMON_ENZYMES: dict[str, RestrictionEnzyme] = {
    enzyme.name: enzyme
    for enzyme in (
        RestrictionEnzyme("EcoRI", "GAATTC", 1, 5),
        RestrictionEnzyme("BamHI", "GGATCC", 1, 5),
        RestrictionEnzyme("HindIII", "AAGCTT", 1, 5),
        RestrictionEnzyme("XhoI", "CTCGAG", 1, 5),
        RestrictionEnzyme("NotI", "GCGGCCGC", 2, 6),
        RestrictionEnzyme("PstI", "CTGCAG", 5, 1),
        RestrictionEnzyme("SmaI", "CCCGGG", 3, 3),
        RestrictionEnzyme("BsaI", "GGTCTC", 7, 11),
    )
}


def enzyme_by_name(name: str) -> RestrictionEnzyme:

    normalized = name.casefold()
    for enzyme_name, enzyme in COMMON_ENZYMES.items():
        if enzyme_name.casefold() == normalized:
            return enzyme
    raise KeyError(f"unknown built-in restriction enzyme: {name}")
