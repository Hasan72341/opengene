from __future__ import annotations

import argparse
from collections.abc import Mapping
from dataclasses import fields, is_dataclass
from enum import Enum
import json
import sys
from pathlib import Path
from uuid import UUID


def _bootstrap_local_packages() -> None:

    root = Path(__file__).resolve().parents[2]
    for package_dir in (
        root / "packages" / "core",
        root / "packages" / "io",
        root / "packages" / "restriction",
        root / "packages" / "primers",
        root / "packages" / "annotations",
        root / "packages" / "alignment",
        root / "packages" / "sequencing",
        root / "packages" / "cloning",
        root / "packages" / "pcr",
    ):
        value = str(package_dir)
        if value not in sys.path:
            sys.path.insert(0, value)


_bootstrap_local_packages()

from opengene_alignment import align_global, align_local  # noqa: E402
from opengene_annotations import find_orfs  # noqa: E402
from opengene_cloning import DNAFragment, plan_gibson_assembly  # noqa: E402
from opengene_core import Project  # noqa: E402
from opengene_io import (  # noqa: E402
    parse_fasta,
    parse_genbank,
    read_ogs,
    write_ogs,
)
from opengene_pcr import plan_pcr  # noqa: E402
from opengene_restriction import COMMON_ENZYMES, digest, scan_restriction_sites  # noqa: E402
from opengene_sequencing import build_consensus, parse_fastq  # noqa: E402


def _json_value(value: object) -> object:

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, (UUID, Path)):
        return str(value)
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value) and not isinstance(value, type):
        return {item.name: _json_value(getattr(value, item.name)) for item in fields(value)}
    if isinstance(value, Mapping):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list, set, frozenset)):
        return [_json_value(item) for item in value]
    raise TypeError(f"cannot serialize {type(value).__name__} as JSON")


def emit(payload: object) -> None:
    print(json.dumps(_json_value(payload), indent=2, sort_keys=True))


def read_document(path: str | Path):

    input_path = Path(path)
    if not input_path.is_file():
        raise ValueError(f"input file does not exist: {input_path}")
    text = input_path.read_text(encoding="utf-8")
    if input_path.suffix.lower() in {".gb", ".gbk", ".genbank"} or text.lstrip().startswith("LOCUS"):
        return parse_genbank(text)
    records = parse_fasta(text)
    if not records:
        raise ValueError("no FASTA records were found")
    if len(records) > 1:
        raise ValueError("this command accepts one FASTA/GenBank record per input file")
    return records[0]


def command_info(args: argparse.Namespace) -> int:
    document = read_document(args.input)
    emit({
        "name": document.name,
        "length": document.length,
        "topology": document.topology.value,
        "molecule_type": document.molecule_type.value,
        "feature_count": len(document.features),
        "sequence_sha256": document.sequence_hash,
    })
    return 0


def command_digest(args: argparse.Namespace) -> int:
    document = read_document(args.input)
    requested = {name.lower() for name in args.enzymes}
    enzymes = [enzyme for enzyme in COMMON_ENZYMES.values() if enzyme.name.lower() in requested]
    unknown = requested - {enzyme.name.lower() for enzyme in enzymes}
    if unknown:
        raise ValueError(f"unknown built-in enzyme(s): {', '.join(sorted(unknown))}")
    sites = scan_restriction_sites(document.sequence, enzymes, topology=document.topology)
    result = digest(document.sequence, enzymes, topology=document.topology)
    emit({
        "sequence": document.name,
        "topology": document.topology.value,
        "cuts": [
            {"enzyme": site.enzyme.name, "position": site.start, "cut": site.cut_position}
            for site in sites
        ],
        "fragments": [
            {"length": fragment.length, "start": fragment.start, "end": fragment.end}
            for fragment in result.fragments
        ],
    })
    return 0


def command_align(args: argparse.Namespace) -> int:
    reference = read_document(args.reference)
    query = read_document(args.query)
    alignment = (
        align_local(reference.sequence, query.sequence)
        if args.local
        else align_global(reference.sequence, query.sequence)
    )
    emit({
        "reference": reference.name,
        "query": query.name,
        "mode": alignment.mode.value,
        "score": alignment.score,
        "identity": alignment.identity,
        "cigar": alignment.cigar,
        "reference_span": [alignment.reference_start, alignment.reference_end],
        "query_span": [alignment.query_start, alignment.query_end],
        "matches": alignment.matches,
        "mismatches": alignment.mismatches,
        "gap_columns": alignment.gap_columns,
        "variants": alignment.variants,
    })
    return 0


def command_orfs(args: argparse.Namespace) -> int:
    document = read_document(args.input)
    orfs = find_orfs(
        document.sequence,
        topology=document.topology.value,
        min_protein_length=args.min_amino_acids,
        include_partial=args.include_partial,
    )
    emit({
        "sequence": document.name,
        "topology": document.topology.value,
        "minimum_amino_acids": args.min_amino_acids,
        "orfs": orfs,
    })
    return 0


def command_pcr(args: argparse.Namespace) -> int:
    document = read_document(args.input)
    plan = plan_pcr(
        document,
        args.forward_primer,
        args.reverse_primer,
        max_mismatches=args.max_mismatches,
        max_products=args.max_products,
        product_name=args.product_name,
    )
    emit({
        "template": document.name,
        "feasible": plan.feasible,
        "forward_primer": plan.forward_primer,
        "reverse_primer": plan.reverse_primer,
        "products": plan.products,
        "warnings": plan.warnings,
        "provenance": plan.provenance,
    })
    return 0 if plan.feasible else 3


def command_gibson(args: argparse.Namespace) -> int:
    fragments = [DNAFragment.from_document(read_document(path)) for path in args.inputs]
    plan = plan_gibson_assembly(
        fragments,
        minimum_overlap=args.minimum_overlap,
        circular_product=not args.linear_product,
        product_name=args.product_name,
    )
    emit({
        "feasible": plan.feasible,
        "fragments": plan.fragments,
        "junctions": plan.junctions,
        "product": plan.product,
        "warnings": plan.warnings,
        "provenance": plan.provenance,
    })
    return 0 if plan.feasible else 3


def command_verify(args: argparse.Namespace) -> int:
    reference = read_document(args.reference)
    reads_path = Path(args.reads)
    if not reads_path.is_file():
        raise ValueError(f"FASTQ file does not exist: {reads_path}")
    reads = parse_fastq(reads_path.read_text(encoding="utf-8"))
    consensus = build_consensus(
        reference.sequence,
        reads,
        minimum_coverage=args.minimum_coverage,
        minimum_base_fraction=args.minimum_base_fraction,
    )
    consensus_payload: dict[str, object] = {
        "sequence": consensus.sequence,
        "length": len(consensus.sequence),
        "coverage": consensus.coverage,
        "ambiguous_positions": consensus.ambiguous_positions,
        "insertions": consensus.insertions,
        "warnings": consensus.warnings,
    }
    if args.include_calls:
        consensus_payload["calls"] = consensus.calls
    emit({
        "reference": reference.name,
        "read_count": len(reads),
        "consensus": consensus_payload,
        "mapped_reads": [
            {
                "name": mapped.read.name,
                "strand": mapped.strand,
                "identity": mapped.alignment.identity,
                "cigar": mapped.alignment.cigar,
                "reference_span": [mapped.alignment.reference_start, mapped.alignment.reference_end],
                "variants": mapped.alignment.variants,
            }
            for mapped in consensus.mapped_reads
        ],
    })
    return 0


def command_project_create(args: argparse.Namespace) -> int:
    documents = tuple(read_document(path) for path in args.inputs)
    project = Project.create(args.name, documents, metadata={"created_by": "opengene CLI"})
    output = Path(args.output)
    write_ogs(output, project, application_version="0.1.0")
    emit({
        "output": str(output),
        "project_id": project.id,
        "name": project.name,
        "document_count": len(project.documents),
        "format": "OpenGene Studio .ogs",
    })
    return 0


def command_project_info(args: argparse.Namespace) -> int:
    project = read_ogs(args.input)
    emit({
        "project_id": project.id,
        "name": project.name,
        "schema_version": project.schema_version,
        "metadata": project.metadata,
        "documents": [
            {
                "id": document.id,
                "name": document.name,
                "length": document.length,
                "topology": document.topology.value,
                "sequence_sha256": document.sequence_hash,
                "feature_count": len(document.features),
            }
            for document in project.documents
        ],
        "provenance_count": len(project.provenance),
    })
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="opengene",
        description="Portable, clean-room molecular design utilities",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    info = subparsers.add_parser("info", help="report basic FASTA or GenBank metadata")
    info.add_argument("input", help="FASTA or GenBank input")
    info.set_defaults(handler=command_info)

    digest_parser = subparsers.add_parser("digest", help="run an independent virtual restriction digest")
    digest_parser.add_argument("input", help="FASTA or GenBank input")
    digest_parser.add_argument("enzymes", nargs="+", help="built-in enzyme names, e.g. EcoRI BamHI")
    digest_parser.set_defaults(handler=command_digest)

    align = subparsers.add_parser("align", help="globally or locally align two sequence documents")
    align.add_argument("reference", help="reference FASTA or GenBank")
    align.add_argument("query", help="query FASTA or GenBank")
    align.add_argument("--local", action="store_true", help="use Smith-Waterman local alignment")
    align.set_defaults(handler=command_align)

    orfs = subparsers.add_parser("orfs", help="scan all six reading frames for ORFs")
    orfs.add_argument("input", help="FASTA or GenBank input")
    orfs.add_argument("--min-amino-acids", type=int, default=30, help="minimum translated amino-acid length (default: 30)")
    orfs.add_argument("--include-partial", action="store_true", help="include ORFs without a terminal stop")
    orfs.set_defaults(handler=command_orfs)

    pcr = subparsers.add_parser("pcr", help="plan exact IUPAC-aware PCR products")
    pcr.add_argument("input", help="template FASTA or GenBank")
    pcr.add_argument("forward_primer", help="forward primer, 5′ to 3′")
    pcr.add_argument("reverse_primer", help="reverse primer, 5′ to 3′")
    pcr.add_argument("--max-mismatches", type=int, default=0)
    pcr.add_argument("--max-products", type=int, default=100)
    pcr.add_argument("--product-name", default="PCR product")
    pcr.set_defaults(handler=command_pcr)

    gibson = subparsers.add_parser("gibson", help="plan exact-overlap Gibson-style assembly")
    gibson.add_argument("inputs", nargs="+", help="ordered fragment FASTA/GenBank files")
    gibson.add_argument("--minimum-overlap", type=int, default=20)
    gibson.add_argument("--linear-product", action="store_true", help="do not require circular closure")
    gibson.add_argument("--product-name", default="Gibson assembly product")
    gibson.set_defaults(handler=command_gibson)

    verify = subparsers.add_parser("verify", help="map FASTQ reads and build a reference-anchored consensus")
    verify.add_argument("reference", help="reference FASTA or GenBank")
    verify.add_argument("reads", help="strict single-line FASTQ reads")
    verify.add_argument("--minimum-coverage", type=int, default=1)
    verify.add_argument("--minimum-base-fraction", type=float, default=0.6)
    verify.add_argument("--include-calls", action="store_true", help="include per-position evidence calls in JSON")
    verify.set_defaults(handler=command_verify)

    project = subparsers.add_parser("project", help="create or inspect portable .ogs projects")
    project_subparsers = project.add_subparsers(dest="project_command", required=True)
    project_create = project_subparsers.add_parser("create", help="create a verified .ogs ZIP + SQLite project")
    project_create.add_argument("output", help="target .ogs path")
    project_create.add_argument("name", help="project name")
    project_create.add_argument("inputs", nargs="+", help="FASTA or GenBank inputs")
    project_create.set_defaults(handler=command_project_create)
    project_info = project_subparsers.add_parser("info", help="inspect and verify a .ogs project")
    project_info.add_argument("input", help="input .ogs path")
    project_info.set_defaults(handler=command_project_info)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.handler(args)
    except (OSError, ValueError, TypeError) as error:
        print(f"opengene: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover - direct source-checkout use
    raise SystemExit(main())
