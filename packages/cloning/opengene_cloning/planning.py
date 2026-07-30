from __future__ import annotations

from collections.abc import Iterable
from dataclasses import replace
from typing import Any

from .models import (
    CloningValidationError,
    CohesiveEnd,
    DNAFragment,
    EndKind,
    EnzymeDefinition,
    GoldenGateJunction,
    GoldenGatePart,
    GoldenGatePlan,
    LigationJunction,
    OligoAlignment,
    OligoAnnealingPlan,
    OverlapAssemblyPlan,
    OverlapJunction,
    RestrictionCut,
    RestrictionDigestPlan,
    RestrictionLigationPlan,
    Topology,
    as_fragment,
    make_operation_record,
    normalize_dna,
    reverse_complement,
)


def _matches_iupac(observed: str, expected: str) -> bool:

    allowed = {
        "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"), "T": frozenset("T"),
        "R": frozenset("AG"), "Y": frozenset("CT"), "S": frozenset("GC"), "W": frozenset("AT"),
        "K": frozenset("GT"), "M": frozenset("AC"), "B": frozenset("CGT"), "D": frozenset("AGT"),
        "H": frozenset("ACT"), "V": frozenset("ACG"), "N": frozenset("ACGT"),
    }
    return all(allowed[left] & allowed[right] for left, right in zip(observed, expected))


def _window(sequence: str, start: int, length: int, *, circular: bool) -> str:
    if not circular:
        return sequence[start : start + length]
    return "".join(sequence[(start + offset) % len(sequence)] for offset in range(length))


def _clockwise_slice(sequence: str, start: int, end: int) -> str:

    if start == end:
        return sequence[start:] + sequence[:start]
    return sequence[start:end] if start < end else sequence[start:] + sequence[:end]


def _coerce_enzymes(enzymes: Iterable[EnzymeDefinition | object]) -> tuple[EnzymeDefinition, ...]:
    selected = tuple(EnzymeDefinition.from_any(enzyme) for enzyme in enzymes)
    if not selected:
        raise CloningValidationError("at least one enzyme is required")
    return selected


def _cut_ends(
    sequence: str,
    *,
    raw_forward_cut: int,
    raw_reverse_cut: int,
    circular: bool,
) -> tuple[CohesiveEnd, CohesiveEnd]:

    if raw_forward_cut == raw_reverse_cut:
        blunt = CohesiveEnd(EndKind.BLUNT)
        return blunt, blunt
    start = min(raw_forward_cut, raw_reverse_cut)
    length = abs(raw_forward_cut - raw_reverse_cut)
    if circular:
        exposed = _window(sequence, start, length, circular=True)
    else:
        exposed = sequence[start : start + length]
    if len(exposed) != length:
        raise CloningValidationError("cut overhang lies outside the supplied sequence")
    kind = EndKind.FIVE_PRIME if raw_forward_cut < raw_reverse_cut else EndKind.THREE_PRIME
    return CohesiveEnd(kind, exposed), CohesiveEnd(kind, reverse_complement(exposed))


def _scan_restriction_cuts(
    fragment: DNAFragment,
    enzymes: tuple[EnzymeDefinition, ...],
) -> tuple[tuple[RestrictionCut, ...], tuple[str, ...]]:
    sequence = fragment.sequence
    length = len(sequence)
    circular = fragment.topology is Topology.CIRCULAR
    cuts: list[RestrictionCut] = []
    warnings: list[str] = []

    for enzyme in enzymes:
        motifs = [(enzyme.recognition_site, 1)]
        if enzyme.reverse_recognition_site != enzyme.recognition_site:
            motifs.append((enzyme.reverse_recognition_site, -1))
        starts = range(length) if circular else range(max(0, length - len(enzyme.recognition_site) + 1))
        for start in starts:
            observed = _window(sequence, start, len(enzyme.recognition_site), circular=circular)
            for motif, strand in motifs:
                if not _matches_iupac(observed, motif):
                    continue
                if strand == 1:
                    raw_forward = start + enzyme.cut_forward
                    raw_reverse = start + enzyme.cut_reverse
                else:
                    raw_forward = start + len(enzyme.recognition_site) - enzyme.cut_reverse
                    raw_reverse = start + len(enzyme.recognition_site) - enzyme.cut_forward
                if not circular and not (0 <= raw_forward <= length and 0 <= raw_reverse <= length):
                    warnings.append(
                        f"{enzyme.name} site at {start} cuts outside linear fragment and was not used"
                    )
                    continue
                try:
                    left_end, right_end = _cut_ends(
                        sequence,
                        raw_forward_cut=raw_forward,
                        raw_reverse_cut=raw_reverse,
                        circular=circular,
                    )
                except CloningValidationError as error:
                    warnings.append(f"{enzyme.name} site at {start} was not used: {error}")
                    continue
                raw_end = start + len(enzyme.recognition_site)
                cuts.append(
                    RestrictionCut(
                        enzyme=enzyme,
                        start=start,
                        end=raw_end if not circular or raw_end <= length else raw_end % length,
                        strand=strand,
                        cut_forward=raw_forward % length if circular else raw_forward,
                        cut_reverse=raw_reverse % length if circular else raw_reverse,
                        left_end=left_end,
                        right_end=right_end,
                        wraps_origin=circular and raw_end > length,
                    )
                )

    by_boundary: dict[int, list[RestrictionCut]] = {}
    for cut in cuts:
        by_boundary.setdefault(cut.cut_forward, []).append(cut)
    selected: list[RestrictionCut] = []
    for boundary in sorted(by_boundary):
        candidates = sorted(
            by_boundary[boundary],
            key=lambda cut: (cut.enzyme.name.casefold(), cut.start, cut.strand),
        )
        first = candidates[0]
        if any(
            candidate.left_end != first.left_end or candidate.right_end != first.right_end
            for candidate in candidates[1:]
        ):
            names = ", ".join(candidate.enzyme.name for candidate in candidates)
            warnings.append(
                f"conflicting ends at cut coordinate {boundary} ({names}); using {first.enzyme.name}"
            )
        selected.append(first)
    return tuple(selected), tuple(warnings)


def plan_restriction_digest(
    input_fragment: DNAFragment | str | object,
    enzymes: Iterable[EnzymeDefinition | object],
    *,
    name: str | None = None,
) -> RestrictionDigestPlan:

    source = as_fragment(input_fragment, name=name or "digest input")
    selected_enzymes = _coerce_enzymes(enzymes)
    cuts, warning_tuple = _scan_restriction_cuts(source, selected_enzymes)
    warnings = list(warning_tuple)
    sequence = source.sequence
    size = len(sequence)
    cut_by_boundary = {cut.cut_forward: cut for cut in cuts}
    boundaries = tuple(sorted(cut_by_boundary))
    fragments: list[DNAFragment] = []

    if not cuts:
        warnings.append("No usable restriction cuts were found; the input remains intact")
        fragments.append(source)
    elif source.topology is Topology.CIRCULAR:
        for index, start in enumerate(boundaries):
            end = boundaries[(index + 1) % len(boundaries)]
            fragment_sequence = _clockwise_slice(sequence, start, end)
            left_cut = cut_by_boundary[start]
            right_cut = cut_by_boundary[end]
            fragments.append(
                DNAFragment(
                    name=f"{source.name} digest {index + 1}",
                    sequence=fragment_sequence,
                    topology=Topology.LINEAR,
                    id=f"{source.id}:digest:{start}:{end}",
                    left_end=left_cut.right_end,
                    right_end=right_cut.left_end,
                    metadata={
                        "source_fragment_id": source.id,
                        "start": start,
                        "end": end,
                        "wraps_origin": end <= start,
                    },
                )
            )
    else:
        all_boundaries = (0, *tuple(boundary for boundary in boundaries if 0 < boundary < size), size)
        for index, (start, end) in enumerate(zip(all_boundaries, all_boundaries[1:])):
            if start == end:
                warnings.append(f"zero-length digest fragment at coordinate {start} was omitted")
                continue
            left_end = cut_by_boundary[start].right_end if start in cut_by_boundary else source.left_end
            right_end = cut_by_boundary[end].left_end if end in cut_by_boundary else source.right_end
            fragments.append(
                DNAFragment(
                    name=f"{source.name} digest {index + 1}",
                    sequence=sequence[start:end],
                    topology=Topology.LINEAR,
                    id=f"{source.id}:digest:{start}:{end}",
                    left_end=left_end,
                    right_end=right_end,
                    metadata={"source_fragment_id": source.id, "start": start, "end": end},
                )
            )

    provenance = make_operation_record(
        "restriction_digest",
        (source,),
        parameters={
            "enzymes": [
                {
                    "name": enzyme.name,
                    "recognition_site": enzyme.recognition_site,
                    "cut_forward": enzyme.cut_forward,
                    "cut_reverse": enzyme.cut_reverse,
                }
                for enzyme in selected_enzymes
            ],
            "topology": source.topology.value,
        },
        warnings=warnings,
        engine_name="opengene-restriction-planner",
    )
    return RestrictionDigestPlan(
        input_fragment=source,
        enzymes=selected_enzymes,
        cuts=cuts,
        fragments=tuple(fragments),
        warnings=tuple(warnings),
        provenance=provenance,
    )


def plan_restriction_ligation(
    fragments: Iterable[DNAFragment | str | object],
    *,
    circular_product: bool = True,
    product_name: str = "ligation product",
) -> RestrictionLigationPlan:

    selected = tuple(as_fragment(item, name=f"fragment {index + 1}") for index, item in enumerate(fragments))
    if not selected:
        raise CloningValidationError("at least one fragment is required for ligation")
    warnings: list[str] = []
    if any(fragment.topology is Topology.CIRCULAR for fragment in selected):
        warnings.append("Circular inputs cannot be directly ligated; digest them to linear fragments first")

    pairs = list(zip(selected, selected[1:]))
    if circular_product:
        pairs.append((selected[-1], selected[0]))
    junctions = tuple(
        LigationJunction(
            left_fragment_id=left.id,
            right_fragment_id=right.id,
            left_end=left.right_end,
            right_end=right.left_end,
            compatible=left.right_end.compatible_with(right.left_end) if left.right_end else False,
        )
        for left, right in pairs
    )
    for junction in junctions:
        if not junction.compatible:
            warnings.append(
                f"incompatible ligation ends between {junction.left_fragment_id} and {junction.right_fragment_id}"
            )
    if len(selected) == 1 and not circular_product:
        warnings.append("A single linear fragment has no ligation junction")

    feasible = (
        not any(fragment.topology is Topology.CIRCULAR for fragment in selected)
        and (bool(junctions) or (len(selected) == 1 and not circular_product))
        and all(junction.compatible for junction in junctions)
    )
    product: DNAFragment | None = None
    if feasible:
        product = DNAFragment(
            name=product_name,
            sequence="".join(fragment.sequence for fragment in selected),
            topology=Topology.CIRCULAR if circular_product else Topology.LINEAR,
            left_end=None if circular_product else selected[0].left_end,
            right_end=None if circular_product else selected[-1].right_end,
            metadata={"assembled_from": tuple(fragment.id for fragment in selected)},
        )
    provenance = make_operation_record(
        "restriction_ligation",
        selected,
        parameters={
            "circular_product": circular_product,
            "ordered_fragment_ids": [fragment.id for fragment in selected],
            "junctions": [
                {
                    "left_fragment_id": junction.left_fragment_id,
                    "right_fragment_id": junction.right_fragment_id,
                    "compatible": junction.compatible,
                }
                for junction in junctions
            ],
        },
        output=product,
        warnings=warnings,
        engine_name="opengene-restriction-planner",
        reproducibility_status="reproducible" if feasible else "infeasible",
    )
    return RestrictionLigationPlan(
        fragments=selected,
        junctions=junctions,
        circular_product=circular_product,
        product=product,
        warnings=tuple(warnings),
        provenance=provenance,
    )


def _longest_suffix_prefix(
    left: str,
    right: str,
    *,
    maximum: int | None = None,
) -> str:
    limit = min(len(left), len(right)) if maximum is None else min(len(left), len(right), maximum)
    for length in range(limit, 0, -1):
        candidate = left[-length:]
        if candidate == right[:length]:
            return candidate
    return ""


def plan_overlap_assembly(
    fragments: Iterable[DNAFragment | str | object],
    *,
    minimum_overlap: int = 20,
    circular_product: bool = True,
    product_name: str = "overlap assembly product",
) -> OverlapAssemblyPlan:

    if minimum_overlap < 1:
        raise CloningValidationError("minimum_overlap must be at least one base")
    selected = tuple(as_fragment(item, name=f"fragment {index + 1}") for index, item in enumerate(fragments))
    if not selected:
        raise CloningValidationError("at least one fragment is required for overlap assembly")
    warnings: list[str] = []
    if any(fragment.topology is Topology.CIRCULAR for fragment in selected):
        warnings.append("Overlap assembly inputs must be linearized before assembly")

    pairs = list(zip(selected, selected[1:]))
    if circular_product and len(selected) > 1:
        pairs.append((selected[-1], selected[0]))
    elif circular_product and len(selected) == 1:
        pairs.append((selected[0], selected[0]))
    junctions: list[OverlapJunction] = []
    for left, right in pairs:
        maximum = left.length - 1 if left is right else None
        overlap = _longest_suffix_prefix(left.sequence, right.sequence, maximum=maximum)
        compatible = len(overlap) >= minimum_overlap
        junctions.append(
            OverlapJunction(left.id, right.id, overlap, len(overlap), compatible)
        )
        if not compatible:
            warnings.append(
                f"{left.name} → {right.name} has {len(overlap)} bp terminal homology; "
                f"requires at least {minimum_overlap} bp"
            )

    feasible = (
        not any(fragment.topology is Topology.CIRCULAR for fragment in selected)
        and all(junction.compatible for junction in junctions)
    )
    product: DNAFragment | None = None
    if feasible:
        assembled = selected[0].sequence
        internal_junctions = junctions[:-1] if circular_product and junctions else junctions
        for fragment, junction in zip(selected[1:], internal_junctions):
            assembled += fragment.sequence[junction.length :]
        if circular_product:
            closure = junctions[-1]
            if closure.length >= len(assembled):
                warnings.append("closing overlap consumes the complete assembled sequence")
                feasible = False
            else:
                assembled = assembled[: -closure.length]
        if feasible:
            product = DNAFragment(
                name=product_name,
                sequence=assembled,
                topology=Topology.CIRCULAR if circular_product else Topology.LINEAR,
                metadata={"assembled_from": tuple(fragment.id for fragment in selected)},
            )
    provenance = make_operation_record(
        "overlap_assembly",
        selected,
        parameters={
            "minimum_overlap": minimum_overlap,
            "circular_product": circular_product,
            "ordered_fragment_ids": [fragment.id for fragment in selected],
            "overlaps": [junction.sequence for junction in junctions],
        },
        output=product,
        warnings=warnings,
        engine_name="opengene-overlap-planner",
        reproducibility_status="reproducible" if product else "infeasible",
    )
    return OverlapAssemblyPlan(
        fragments=selected,
        minimum_overlap=minimum_overlap,
        circular_product=circular_product,
        junctions=tuple(junctions),
        product=product,
        warnings=tuple(warnings),
        provenance=provenance,
    )


def plan_gibson_assembly(
    fragments: Iterable[DNAFragment | str | object],
    *,
    minimum_overlap: int = 20,
    circular_product: bool = True,
    product_name: str = "Gibson assembly product",
) -> OverlapAssemblyPlan:

    plan = plan_overlap_assembly(
        fragments,
        minimum_overlap=minimum_overlap,
        circular_product=circular_product,
        product_name=product_name,
    )
    record = make_operation_record(
        "gibson_assembly",
        plan.fragments,
        parameters={
            "minimum_overlap": minimum_overlap,
            "circular_product": circular_product,
            "ordered_fragment_ids": [fragment.id for fragment in plan.fragments],
            "overlaps": [junction.sequence for junction in plan.junctions],
        },
        output=plan.product,
        warnings=plan.warnings,
        engine_name="opengene-overlap-planner",
        reproducibility_status="reproducible" if plan.product else "infeasible",
    )
    return replace(plan, provenance=record)


def _canonical_overhang(overhang: str) -> str:
    return min(overhang, reverse_complement(overhang))


def plan_golden_gate(
    parts: Iterable[GoldenGatePart],
    *,
    enzyme_name: str = "Type IIS enzyme",
    circular_product: bool = True,
    overhang_length: int = 4,
    product_name: str = "Golden Gate product",
) -> GoldenGatePlan:

    if not isinstance(enzyme_name, str) or not enzyme_name.strip():
        raise CloningValidationError("enzyme_name cannot be empty")
    if overhang_length < 1:
        raise CloningValidationError("overhang_length must be at least one")
    selected = tuple(parts)
    if not selected:
        raise CloningValidationError("at least one Golden Gate part is required")
    if not all(isinstance(part, GoldenGatePart) for part in selected):
        raise TypeError("parts must contain GoldenGatePart values")
    warnings: list[str] = []
    for part in selected:
        if len(part.left_overhang) != overhang_length or len(part.right_overhang) != overhang_length:
            warnings.append(
                f"{part.fragment.name} has overhang length(s) other than {overhang_length} bp"
            )
    pairs = list(zip(selected, selected[1:]))
    if circular_product:
        pairs.append((selected[-1], selected[0]))
    junctions = tuple(
        GoldenGateJunction(
            left_part_id=left.fragment.id,
            right_part_id=right.fragment.id,
            left_overhang=left.right_overhang,
            right_overhang=right.left_overhang,
            compatible=left.right_overhang == reverse_complement(right.left_overhang),
        )
        for left, right in pairs
    )
    for junction in junctions:
        if not junction.compatible:
            warnings.append(
                f"incompatible Golden Gate overhangs between {junction.left_part_id} and {junction.right_part_id}"
            )
    seen_junctions: dict[str, int] = {}
    for index, junction in enumerate(junctions, start=1):
        if not junction.compatible:
            continue
        key = _canonical_overhang(junction.left_overhang)
        if key in seen_junctions:
            warnings.append(
                f"junction {index} reuses overhang class from junction {seen_junctions[key]}; misassembly is possible"
            )
        else:
            seen_junctions[key] = index
    reused = len(seen_junctions) != sum(junction.compatible for junction in junctions)
    feasible = all(junction.compatible for junction in junctions) and not reused and not any(
        len(part.left_overhang) != overhang_length or len(part.right_overhang) != overhang_length
        for part in selected
    )
    product: DNAFragment | None = None
    if feasible:
        product = DNAFragment(
            name=product_name,
            sequence="".join(part.fragment.sequence for part in selected),
            topology=Topology.CIRCULAR if circular_product else Topology.LINEAR,
            left_end=None if circular_product else CohesiveEnd(EndKind.FIVE_PRIME, selected[0].left_overhang),
            right_end=None if circular_product else CohesiveEnd(EndKind.FIVE_PRIME, selected[-1].right_overhang),
            metadata={"assembled_from": tuple(part.fragment.id for part in selected)},
        )
    inputs = tuple(part.fragment for part in selected)
    provenance = make_operation_record(
        "golden_gate_assembly",
        inputs,
        parameters={
            "enzyme_name": enzyme_name.strip(),
            "circular_product": circular_product,
            "overhang_length": overhang_length,
            "parts": [
                {
                    "fragment_id": part.fragment.id,
                    "left_overhang": part.left_overhang,
                    "right_overhang": part.right_overhang,
                }
                for part in selected
            ],
        },
        output=product,
        warnings=warnings,
        engine_name="opengene-golden-gate-planner",
        reproducibility_status="reproducible" if product else "infeasible",
    )
    return GoldenGatePlan(
        enzyme_name=enzyme_name.strip(),
        parts=selected,
        circular_product=circular_product,
        junctions=junctions,
        product=product,
        warnings=tuple(warnings),
        provenance=provenance,
    )


def _best_oligo_alignment(
    top: str,
    bottom_on_top: str,
    *,
    minimum_overlap: int,
) -> OligoAlignment | None:
    candidates: list[OligoAlignment] = []
    for shift in range(-len(bottom_on_top) + minimum_overlap, len(top) - minimum_overlap + 1):
        start = max(0, shift)
        end = min(len(top), shift + len(bottom_on_top))
        overlap = end - start
        if overlap < minimum_overlap:
            continue
        mismatches = sum(
            top[position] != bottom_on_top[position - shift]
            for position in range(start, end)
        )
        candidates.append(OligoAlignment(0, shift, overlap, mismatches))
    if not candidates:
        return None
    return min(
        candidates,
        key=lambda alignment: (
            -alignment.overlap_length,
            alignment.mismatches,
            abs(alignment.bottom_start),
            alignment.bottom_start,
        ),
    )


def _duplex_forward_sequence(top: str, bottom_on_top: str, shift: int) -> str:
    start = min(0, shift)
    end = max(len(top), shift + len(bottom_on_top))
    bases: list[str] = []
    for position in range(start, end):
        top_base = top[position] if 0 <= position < len(top) else None
        bottom_base = bottom_on_top[position - shift] if 0 <= position - shift < len(bottom_on_top) else None
        bases.append(top_base or bottom_base or "N")
    return "".join(bases)


def plan_oligo_annealing(
    top_oligo: str,
    bottom_oligo: str,
    *,
    minimum_overlap: int = 12,
    maximum_mismatches: int = 0,
    product_name: str = "annealed oligos",
) -> OligoAnnealingPlan:

    if minimum_overlap < 1:
        raise CloningValidationError("minimum_overlap must be at least one")
    if maximum_mismatches < 0:
        raise CloningValidationError("maximum_mismatches cannot be negative")
    top = normalize_dna(top_oligo, allow_ambiguous=False)
    bottom = normalize_dna(bottom_oligo, allow_ambiguous=False)
    if not top or not bottom:
        raise CloningValidationError("oligos cannot be empty")
    top_input = DNAFragment("top oligo", top, id="oligo:top:" + top)
    bottom_input = DNAFragment("bottom oligo", bottom, id="oligo:bottom:" + bottom)
    bottom_on_top = reverse_complement(bottom)
    alignment = _best_oligo_alignment(top, bottom_on_top, minimum_overlap=minimum_overlap)
    warnings: list[str] = []
    feasible = alignment is not None and alignment.mismatches <= maximum_mismatches
    if alignment is None:
        warnings.append(f"No alignment meets the {minimum_overlap} bp minimum overlap")
    elif alignment.mismatches > maximum_mismatches:
        warnings.append(
            f"Best {alignment.overlap_length} bp alignment has {alignment.mismatches} mismatch(es); "
            f"maximum is {maximum_mismatches}"
        )
    if alignment and (alignment.bottom_start != 0 or len(bottom_on_top) != len(top)):
        warnings.append("Annealed product contains one or more single-stranded terminal extensions")
    product: DNAFragment | None = None
    if feasible and alignment:
        product = DNAFragment(
            name=product_name,
            sequence=_duplex_forward_sequence(top, bottom_on_top, alignment.bottom_start),
            topology=Topology.LINEAR,
            metadata={
                "top_oligo": top,
                "bottom_oligo": bottom,
                "overlap_length": alignment.overlap_length,
                "mismatches": alignment.mismatches,
            },
        )
    provenance = make_operation_record(
        "oligo_annealing",
        (top_input, bottom_input),
        parameters={
            "minimum_overlap": minimum_overlap,
            "maximum_mismatches": maximum_mismatches,
            "alignment": {
                "top_start": alignment.top_start,
                "bottom_start": alignment.bottom_start,
                "overlap_length": alignment.overlap_length,
                "mismatches": alignment.mismatches,
            }
            if alignment
            else None,
        },
        output=product,
        warnings=warnings,
        engine_name="opengene-oligo-planner",
        reproducibility_status="reproducible" if product else "infeasible",
    )
    return OligoAnnealingPlan(
        top_oligo=top,
        bottom_oligo=bottom,
        alignment=alignment,
        product=product,
        warnings=tuple(warnings),
        provenance=provenance,
    )
