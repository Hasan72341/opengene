import json

from opengene_pcr import (
    PCRPrimer,
    PCRTemplate,
    PCRValidationError,
    SequenceMutation,
    Topology,
    find_primer_binding_sites,
    plan_pcr,
    plan_site_directed_mutagenesis,
)


def test_pcr_incorporates_extensions_and_retains_template() -> None:
    template = PCRTemplate("template", "AAATGCGTACGTTAGCGGTT")
    forward = PCRPrimer("F", "ATGCG", "GG")
    reverse = PCRPrimer("R", "AACCG", "TT")

    plan = plan_pcr(template, forward, reverse)

    assert template.sequence == "AAATGCGTACGTTAGCGGTT"
    assert plan.feasible
    assert len(plan.products) == 1
    assert plan.products[0].sequence == "GGATGCGTACGTTAGCGGTTAA"
    assert plan.products[0].length == 22
    assert plan.provenance.operation_id.startswith("op-")
    assert json.loads(json.dumps(plan.provenance.to_dict()))["operation_type"] == "pcr"


def test_binding_and_circular_pcr_cross_origin() -> None:
    template = PCRTemplate("circular", "TTAAACCCGGGATGCG", topology=Topology.CIRCULAR)
    forward = PCRPrimer("F", "ATGCG")
    reverse = PCRPrimer("R", "TTTAA")

    sites = find_primer_binding_sites(template, reverse)
    plan = plan_pcr(template, forward, reverse)

    assert any(site.strand == -1 and site.start == 0 for site in sites)
    assert plan.feasible
    assert plan.products[0].sequence == "ATGCGTTAAA"
    assert plan.products[0].circular_template is True


def test_pcr_reports_multiple_candidates_and_deterministically_caps_them() -> None:
    template = PCRTemplate("repeat", "ATGCGAAAAACGGTTATGCGAAAAACGGTT")
    forward = PCRPrimer("F", "ATGCG")
    reverse = PCRPrimer("R", "AACCG")

    uncapped = plan_pcr(template, forward, reverse)
    capped = plan_pcr(template, forward, reverse, max_products=1)

    assert len(uncapped.products) > 1
    assert any("possible amplicons" in warning for warning in uncapped.warnings)
    assert len(capped.products) == 1
    assert any("limited" in warning for warning in capped.warnings)


def test_site_directed_mutagenesis_applies_ordered_half_open_replacements() -> None:
    template = PCRTemplate("template", "AAATGCGTACGTTAGCGGTT", topology=Topology.CIRCULAR)
    mutations = (
        SequenceMutation(3, 4, "C", "substitution"),
        SequenceMutation(7, 7, "GG", "insertion"),
        SequenceMutation(14, 16, "", "deletion"),
    )

    plan = plan_site_directed_mutagenesis(template, mutations)

    assert plan.feasible and plan.product is not None
    assert template.sequence == "AAATGCGTACGTTAGCGGTT"
    assert plan.product.sequence == "AAACGCGGGTACGTTAGGTT"
    assert plan.product.topology is Topology.CIRCULAR
    assert plan.provenance.reproducibility_status == "requires_review"
    assert "No mutagenic primers" in plan.warnings[0]


def test_mutagenesis_can_retain_pcr_validation_and_rejects_overlaps() -> None:
    template = PCRTemplate("template", "AAATGCGTACGTTAGCGGTT")
    forward = PCRPrimer("F", "ATGCG")
    reverse = PCRPrimer("R", "AACCG")
    plan = plan_site_directed_mutagenesis(
        template,
        [SequenceMutation(3, 4, "C", "change")],
        forward_primer=forward,
        reverse_primer=reverse,
    )

    assert plan.pcr_plan is not None
    assert plan.pcr_plan.feasible
    assert plan.product is not None
    try:
        plan_site_directed_mutagenesis(
            template,
            [SequenceMutation(2, 5, "A"), SequenceMutation(4, 6, "T")],
        )
    except PCRValidationError:
        pass
    else:
        raise AssertionError("overlapping mutations should be rejected")
