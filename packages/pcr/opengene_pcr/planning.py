from __future__ import annotations

from hashlib import sha256
from typing import Iterable

from .models import (
    MutagenesisPlan,
    MutagenesisProduct,
    PCRBindingSite,
    PCRPlan,
    PCRPrimer,
    PCRProduct,
    PCRProvenanceRecord,
    PCRTemplate,
    PCRValidationError,
    SequenceMutation,
    Topology,
    as_template,
    make_pcr_provenance,
    normalize_dna,
    reverse_complement,
)


_IUPAC: dict[str, frozenset[str]] = {
    "A": frozenset("A"), "C": frozenset("C"), "G": frozenset("G"), "T": frozenset("T"),
    "R": frozenset("AG"), "Y": frozenset("CT"), "S": frozenset("GC"), "W": frozenset("AT"),
    "K": frozenset("GT"), "M": frozenset("AC"), "B": frozenset("CGT"), "D": frozenset("AGT"),
    "H": frozenset("ACT"), "V": frozenset("ACG"), "N": frozenset("ACGT"),
}


def _matches(window: str, primer: str) -> int:

    return sum(not (_IUPAC[template_base] & _IUPAC[primer_base]) for template_base, primer_base in zip(window, primer))


def _window(sequence: str, start: int, length: int, *, circular: bool) -> str:
    if not circular:
        return sequence[start : start + length]
    return "".join(sequence[(start + offset) % len(sequence)] for offset in range(length))


def _clockwise_between(sequence: str, start: int, end: int) -> str:
    if start == end:
        return ""
    return sequence[start:end] if start < end else sequence[start:] + sequence[:end]


def find_primer_binding_sites(
    template: PCRTemplate | str | object,
    primer: PCRPrimer | str | object,
    *,
    circular: bool | None = None,
    max_mismatches: int = 0,
) -> tuple[PCRBindingSite, ...]:

    if max_mismatches < 0:
        raise PCRValidationError("max_mismatches cannot be negative")
    source = as_template(template)
    selected_primer = PCRPrimer.from_any(primer, default_name="Primer")
    circular_template = source.topology is Topology.CIRCULAR if circular is None else bool(circular)
    sequence = source.sequence
    primer_sequence = selected_primer.annealing_sequence
    if len(primer_sequence) > len(sequence):
        return ()
    starts = range(len(sequence)) if circular_template else range(len(sequence) - len(primer_sequence) + 1)
    reverse = reverse_complement(primer_sequence)
    sites: list[PCRBindingSite] = []
    for start in starts:
        observed = _window(sequence, start, len(primer_sequence), circular=circular_template)
        raw_end = start + len(primer_sequence)
        end = raw_end if not circular_template or raw_end <= len(sequence) else raw_end % len(sequence)
        forward_mismatches = _matches(observed, primer_sequence)
        if forward_mismatches <= max_mismatches:
            sites.append(
                PCRBindingSite(
                    start=start,
                    end=end,
                    strand=1,
                    mismatches=forward_mismatches,
                    wraps_origin=circular_template and raw_end > len(sequence),
                )
            )
        reverse_mismatches = _matches(observed, reverse)
        if reverse_mismatches <= max_mismatches:
            sites.append(
                PCRBindingSite(
                    start=start,
                    end=end,
                    strand=-1,
                    mismatches=reverse_mismatches,
                    wraps_origin=circular_template and raw_end > len(sequence),
                )
            )
    return tuple(sorted(sites, key=lambda site: (site.start, site.strand, site.mismatches)))


def _product_identifier(
    template: PCRTemplate,
    forward: PCRPrimer,
    reverse: PCRPrimer,
    forward_site: PCRBindingSite,
    reverse_site: PCRBindingSite,
    sequence: str,
) -> str:
    payload = "|".join(
        (
            template.id,
            forward.full_sequence,
            reverse.full_sequence,
            str(forward_site.start),
            str(reverse_site.start),
            sequence,
        )
    )
    return "pcr-" + sha256(payload.encode("utf-8")).hexdigest()[:24]


def plan_pcr(
    template: PCRTemplate | str | object,
    forward_primer: PCRPrimer | str | object,
    reverse_primer: PCRPrimer | str | object,
    *,
    circular: bool | None = None,
    max_mismatches: int = 0,
    max_products: int | None = 100,
    product_name: str = "PCR product",
) -> PCRPlan:

    if max_mismatches < 0:
        raise PCRValidationError("max_mismatches cannot be negative")
    if max_products is not None and max_products < 1:
        raise PCRValidationError("max_products must be positive or None")
    source = as_template(template)
    forward = PCRPrimer.from_any(forward_primer, default_name="Forward")
    reverse = PCRPrimer.from_any(reverse_primer, default_name="Reverse")
    circular_template = source.topology is Topology.CIRCULAR if circular is None else bool(circular)
    forward_sites = tuple(
        site
        for site in find_primer_binding_sites(
            source, forward, circular=circular_template, max_mismatches=max_mismatches
        )
        if site.strand == 1
    )
    reverse_sites = tuple(
        site
        for site in find_primer_binding_sites(
            source, reverse, circular=circular_template, max_mismatches=max_mismatches
        )
        if site.strand == -1
    )
    candidates: list[PCRProduct] = []
    seen: set[tuple[int, int]] = set()
    sequence = source.sequence
    for forward_site in forward_sites:
        for reverse_site in reverse_sites:
            key = (forward_site.start, reverse_site.start)
            if key in seen:
                continue
            if not circular_template:
                if forward_site.end > reverse_site.start:
                    continue
                interior = sequence[forward_site.end : reverse_site.start]
            else:
                if forward_site.start == reverse_site.start and forward_site.end == reverse_site.end:
                    continue
                forward_end = (forward_site.start + len(forward.annealing_sequence)) % len(sequence)
                interior = _clockwise_between(sequence, forward_end, reverse_site.start)
            product_sequence = (
                forward.full_sequence + interior + reverse_complement(reverse.full_sequence)
            )
            seen.add(key)
            candidates.append(
                PCRProduct(
                    id=_product_identifier(
                        source, forward, reverse, forward_site, reverse_site, product_sequence
                    ),
                    name=product_name,
                    sequence=product_sequence,
                    forward_site=forward_site,
                    reverse_site=reverse_site,
                    circular_template=circular_template,
                )
            )
    candidates.sort(
        key=lambda product: (
            product.length,
            product.forward_site.start,
            product.reverse_site.start,
            product.forward_site.mismatches + product.reverse_site.mismatches,
            product.sequence,
        )
    )
    warnings: list[str] = []
    if not candidates:
        warnings.append("No correctly oriented primer pair produced an amplicon")
    elif len(candidates) > 1:
        warnings.append(f"{len(candidates)} possible amplicons found; inspect primer specificity")
    if max_mismatches:
        warnings.append("Amplicons include primer bindings permitted by max_mismatches")
    if max_products is not None and len(candidates) > max_products:
        warnings.append(f"Amplicon results were limited to the first {max_products} deterministic candidates")
        candidates = candidates[:max_products]
    named_products = tuple(
        PCRProduct(
            id=product.id,
            name=f"{product_name} {index}",
            sequence=product.sequence,
            forward_site=product.forward_site,
            reverse_site=product.reverse_site,
            circular_template=product.circular_template,
        )
        for index, product in enumerate(candidates, start=1)
    )
    provenance = make_pcr_provenance(
        "pcr",
        source,
        (forward, reverse),
        parameters={
            "circular_template": circular_template,
            "max_mismatches": max_mismatches,
            "max_products": max_products,
            "forward_primer": {
                "name": forward.name,
                "annealing_sequence": forward.annealing_sequence,
                "five_prime_extension": forward.five_prime_extension,
            },
            "reverse_primer": {
                "name": reverse.name,
                "annealing_sequence": reverse.annealing_sequence,
                "five_prime_extension": reverse.five_prime_extension,
            },
        },
        outputs=named_products,
        warnings=warnings,
        reproducibility_status="reproducible" if named_products else "infeasible",
    )
    return PCRPlan(
        template=source,
        forward_primer=forward,
        reverse_primer=reverse,
        circular_template=circular_template,
        max_mismatches=max_mismatches,
        products=named_products,
        warnings=tuple(warnings),
        provenance=provenance,
    )


def _validate_mutations(
    template: PCRTemplate,
    mutations: Iterable[SequenceMutation],
) -> tuple[SequenceMutation, ...]:
    selected = tuple(mutations)
    if not selected:
        raise PCRValidationError("at least one mutation is required")
    if not all(isinstance(mutation, SequenceMutation) for mutation in selected):
        raise TypeError("mutations must contain SequenceMutation values")
    ordered = tuple(sorted(selected, key=lambda mutation: (mutation.start, mutation.end, mutation.name)))
    previous_end = -1
    previous_start = -1
    for mutation in ordered:
        if mutation.end > template.length:
            raise PCRValidationError(
                f"mutation {mutation.name!r} lies outside template length {template.length}"
            )
        if mutation.start < previous_end or (
            mutation.start == previous_start and mutation.start == mutation.end
        ):
            raise PCRValidationError("mutations overlap or have ambiguous shared insertion coordinates")
        previous_end = mutation.end
        previous_start = mutation.start
    return ordered


def _apply_mutations(sequence: str, mutations: tuple[SequenceMutation, ...]) -> str:
    cursor = 0
    pieces: list[str] = []
    for mutation in mutations:
        pieces.append(sequence[cursor : mutation.start])
        pieces.append(mutation.replacement)
        cursor = mutation.end
    pieces.append(sequence[cursor:])
    return "".join(pieces)


def _mutation_primer_warnings(
    mutations: tuple[SequenceMutation, ...],
    primers: tuple[PCRPrimer, PCRPrimer] | None,
    pcr_plan: PCRPlan | None,
) -> list[str]:
    warnings: list[str] = []
    if primers is None:
        warnings.append("No mutagenic primers were supplied for PCR validation")
        return warnings
    primer_sequences = tuple(
        sequence
        for primer in primers
        for sequence in (primer.full_sequence, reverse_complement(primer.full_sequence))
    )
    for mutation in mutations:
        if mutation.replacement and not any(mutation.replacement in sequence for sequence in primer_sequences):
            warnings.append(
                f"Replacement for {mutation.name!r} was not found in either supplied primer sequence"
            )
    if pcr_plan is not None and not pcr_plan.products:
        warnings.append("Supplied mutagenic primers did not yield a simulated oriented amplicon")
    return warnings


def plan_site_directed_mutagenesis(
    template: PCRTemplate | str | object,
    mutations: Iterable[SequenceMutation],
    *,
    forward_primer: PCRPrimer | str | object | None = None,
    reverse_primer: PCRPrimer | str | object | None = None,
    circular: bool | None = None,
    max_mismatches: int = 0,
    product_name: str = "mutagenized product",
) -> MutagenesisPlan:

    if (forward_primer is None) != (reverse_primer is None):
        raise PCRValidationError("supply both mutagenic primers or neither")
    source = as_template(template)
    topology = source.topology if circular is None else (Topology.CIRCULAR if circular else Topology.LINEAR)
    selected_mutations = _validate_mutations(source, mutations)
    mutated_sequence = _apply_mutations(source.sequence, selected_mutations)
    if not mutated_sequence:
        raise PCRValidationError("mutagenesis cannot produce an empty sequence")
    pcr_plan: PCRPlan | None = None
    primers: tuple[PCRPrimer, PCRPrimer] | None = None
    if forward_primer is not None and reverse_primer is not None:
        forward = PCRPrimer.from_any(forward_primer, default_name="Mutagenic forward")
        reverse = PCRPrimer.from_any(reverse_primer, default_name="Mutagenic reverse")
        primers = (forward, reverse)
        pcr_plan = plan_pcr(
            source,
            forward,
            reverse,
            circular=topology is Topology.CIRCULAR,
            max_mismatches=max_mismatches,
            product_name="mutagenic PCR product",
        )
    warnings = _mutation_primer_warnings(selected_mutations, primers, pcr_plan)
    product_payload = "|".join(
        (source.id, mutated_sequence, topology.value, *(f"{item.start}:{item.end}:{item.replacement}" for item in selected_mutations))
    )
    product = MutagenesisProduct(
        id="mut-" + sha256(product_payload.encode("utf-8")).hexdigest()[:24],
        name=product_name,
        sequence=mutated_sequence,
        topology=topology,
        mutations=selected_mutations,
    )
    input_ids = (source.id, *(f"primer:{primer.name}" for primer in primers or ()))
    input_hashes = {
        source.id: source.sha256,
        **{
            f"primer:{primer.name}": sha256(primer.full_sequence.encode("ascii")).hexdigest()
            for primer in primers or ()
        },
    }
    parameters = {
        "topology": topology.value,
        "max_mismatches": max_mismatches,
        "mutations": [
            {
                "name": mutation.name,
                "start": mutation.start,
                "end": mutation.end,
                "replacement": mutation.replacement,
            }
            for mutation in selected_mutations
        ],
        "primer_validation_operation_id": pcr_plan.provenance.operation_id if pcr_plan else None,
    }
    provenance = PCRProvenanceRecord(
        operation_id="",
        operation_type="site_directed_mutagenesis",
        input_ids=input_ids,
        input_sequence_hashes=input_hashes,
        parameters=parameters,
        output_sequence_hashes=(product.sha256,),
        warnings=tuple(warnings),
        reproducibility_status="reproducible" if not warnings else "requires_review",
    )
    return MutagenesisPlan(
        template=source,
        mutations=selected_mutations,
        product=product,
        pcr_plan=pcr_plan,
        warnings=tuple(warnings),
        provenance=provenance,
    )


def plan_mutagenesis(
    template: PCRTemplate | str | object,
    mutations: Iterable[SequenceMutation],
    **kwargs: object,
) -> MutagenesisPlan:

    return plan_site_directed_mutagenesis(template, mutations, **kwargs)
