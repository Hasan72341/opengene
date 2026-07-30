import json

from opengene_cloning import (
    DNAFragment,
    EnzymeDefinition,
    GoldenGatePart,
    Topology,
    plan_gibson_assembly,
    plan_golden_gate,
    plan_oligo_annealing,
    plan_overlap_assembly,
    plan_restriction_digest,
    plan_restriction_ligation,
)


def test_restriction_digest_exposes_complementary_ends_and_ligates_without_mutating_input() -> None:
    source = DNAFragment("vector", "AAAAGAATTCCCC")
    enzyme = EnzymeDefinition("EcoRI", "GAATTC", 1, 5)

    digest = plan_restriction_digest(source, [enzyme])

    assert source.sequence == "AAAAGAATTCCCC"
    assert [fragment.sequence for fragment in digest.fragments] == ["AAAAG", "AATTCCCC"]
    assert digest.fragments[0].right_end is not None
    assert digest.fragments[0].right_end.compatible_with(digest.fragments[1].left_end)
    ligation = plan_restriction_ligation(digest.fragments, circular_product=False)
    assert ligation.feasible
    assert ligation.product is not None
    assert ligation.product.sequence == source.sequence
    assert ligation.provenance.operation_id.startswith("op-")


def test_circular_digest_can_religate_and_provenance_is_deterministic() -> None:
    source = DNAFragment("plasmid", "AAAAGAATTCCCC", topology=Topology.CIRCULAR)
    enzyme = EnzymeDefinition("EcoRI", "GAATTC", 1, 5)

    first = plan_restriction_digest(source, [enzyme])
    second = plan_restriction_digest(source, [enzyme])
    ligation = plan_restriction_ligation(first.fragments, circular_product=True)

    assert first.provenance.operation_id == second.provenance.operation_id
    assert ligation.feasible
    assert ligation.product is not None
    assert ligation.product.topology is Topology.CIRCULAR
    assert ligation.product.sequence == source.sequence[5:] + source.sequence[:5]
    assert json.loads(json.dumps(ligation.provenance.to_dict()))["operation_type"] == "restriction_ligation"


def test_overlap_and_named_gibson_plans_trim_declared_homology() -> None:
    left = DNAFragment("left", "AAACCC")
    right = DNAFragment("right", "CCCGGGAAA")

    plan = plan_overlap_assembly([left, right], minimum_overlap=3, circular_product=True)
    gibson = plan_gibson_assembly([left, right], minimum_overlap=3, circular_product=True)

    assert plan.feasible and gibson.feasible
    assert plan.product is not None and gibson.product is not None
    assert plan.product.sequence == "AAACCCGGG"
    assert gibson.product.sequence == plan.product.sequence
    assert [junction.length for junction in plan.junctions] == [3, 3]
    assert gibson.provenance.operation_type == "gibson_assembly"


def test_overlap_plan_reports_insufficient_homology_without_product() -> None:
    plan = plan_overlap_assembly(
        [DNAFragment("left", "AAAAAC"), DNAFragment("right", "CGGGGG")],
        minimum_overlap=3,
        circular_product=False,
    )

    assert not plan.feasible
    assert plan.product is None
    assert "requires at least 3 bp" in plan.warnings[0]
    assert plan.provenance.reproducibility_status == "infeasible"


def test_golden_gate_requires_complementary_unique_overhangs() -> None:
    first = GoldenGatePart(DNAFragment("A", "AAAA"), "AAGC", "AATG")
    second = GoldenGatePart(DNAFragment("B", "CCCC"), "CATT", "GCTT")

    plan = plan_golden_gate([first, second], enzyme_name="BsaI", circular_product=True)

    assert plan.feasible
    assert plan.product is not None
    assert plan.product.sequence == "AAAACCCC"
    assert all(junction.compatible for junction in plan.junctions)

    reused = plan_golden_gate(
        [
            GoldenGatePart(DNAFragment("one", "AAAA"), "AATG", "AATG"),
            GoldenGatePart(DNAFragment("two", "CCCC"), "CATT", "CATT"),
        ],
        circular_product=True,
    )
    assert not reused.feasible
    assert any("reuses overhang" in warning for warning in reused.warnings)


def test_oligo_annealing_reports_extensions_and_mismatch_limits() -> None:
    exact = plan_oligo_annealing("ATGCG", "CGCAT", minimum_overlap=5)
    staggered = plan_oligo_annealing("GGATGCG", "CGCAT", minimum_overlap=5)
    mismatch = plan_oligo_annealing("ATGCG", "CGCAA", minimum_overlap=5)

    assert exact.feasible and exact.product is not None
    assert exact.product.sequence == "ATGCG"
    assert staggered.feasible and staggered.product is not None
    assert any("single-stranded" in warning for warning in staggered.warnings)
    assert not mismatch.feasible
    assert mismatch.alignment is not None and mismatch.alignment.mismatches == 1
