from __future__ import annotations

from collections.abc import Iterable
from typing import Any

from opengene_cloning import OverlapAssemblyPlan, plan_gibson_assembly
from opengene_core import SequenceDocument, Topology


class AssemblyInfeasibleError(ValueError):
    pass


def gibson_plan(
    fragments: Iterable[SequenceDocument | str | object],
    *,
    minimum_overlap: int = 20,
    circular_product: bool = True,
    product_name: str = "Gibson assembly product",
) -> OverlapAssemblyPlan:

    return plan_gibson_assembly(
        tuple(fragments),
        minimum_overlap=minimum_overlap,
        circular_product=circular_product,
        product_name=product_name,
    )


def gibson_assembly(
    fragments: Iterable[SequenceDocument | str | object],
    *,
    minimum_overlap: int = 20,
    circular_product: bool = True,
    product_name: str = "Gibson assembly product",
) -> SequenceDocument:

    source_values = tuple(fragments)
    plan = gibson_plan(
        source_values,
        minimum_overlap=minimum_overlap,
        circular_product=circular_product,
        product_name=product_name,
    )
    if not plan.product:
        warning = "; ".join(plan.warnings) or "assembly has no feasible product"
        raise AssemblyInfeasibleError(warning)
    input_ids = [
        str(getattr(fragment, "id"))
        for fragment in source_values
        if getattr(fragment, "id", None) is not None
    ]
    record: dict[str, Any] = plan.provenance.to_dict()
    recipe = {
        "operation_type": record["operation_type"],
        "input_document_ids": input_ids,
        "parameters": {"engine_record": record},
        "software_version": "OpenGene Studio 0.1.0",
        "engine_version": record.get("engine_version", "opengene-overlap-planner"),
        "warnings": list(plan.warnings),
        "reproducibility_status": record.get("reproducibility_status", "reproducible"),
    }
    return SequenceDocument.create(
        name=product_name,
        sequence=plan.product.sequence,
        topology=Topology(plan.product.topology.value),
        metadata={
            "assembly_method": "gibson",
            "source_fragment_hashes": {fragment.id: fragment.sha256 for fragment in plan.fragments},
            "opengene.provenance": recipe,
        },
    )


__all__ = ["AssemblyInfeasibleError", "gibson_assembly", "gibson_plan"]
