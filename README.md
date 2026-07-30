# OpenGene Studio

OpenGene Studio is a clean-room, local-first molecular-biology workstation. It is intentionally its own product: it does not copy SnapGene code, artwork, proprietary databases, or interaction design. It provides familiar scientific workflows through an independently designed, touch-first interface.

This repository ships a runnable, offline-friendly touch workstation, a portable Python domain core, a scriptable CLI, an implemented `.ogs` ZIP + SQLite project container, and an optional Qt Quick native host. The browser shell is deliberately dependency-free so it can be opened on a desktop, tablet, or lab touchscreen without a cloud account.

## Run it

```bash
npm run dev
# open http://localhost:4174
```

The starter project is loaded locally and changes are autosaved in browser storage. Use **Import** for FASTA or GenBank, **Export** to download standard interchange files, and **Save project** to retain an editable portable project snapshot. No sequence is sent anywhere.

For automation and the native project container:

```bash
python3 apps/cli/main.py --help
python3 apps/cli/main.py project create project.ogs "My project" construct.gb insert.fasta
python3 apps/cli/main.py project info project.ogs
```

The optional Qt Quick host can be checked without installing Qt:

```bash
python3 -m apps.desktop --check
```

## Run it in Docker

Three build targets share one `Dockerfile`: `web` (nginx serving the browser shell), `cli` (the installed `opengene` command), and `test` (both suites plus the comment check).

```bash
# browser workstation on http://localhost:4174
docker compose up -d web

# scriptable CLI; ./tests/fixtures is mounted read-only at /work/fixtures
docker compose --profile tools run --rm cli info fixtures/cli-demo.fasta
docker compose --profile tools run --rm cli digest fixtures/cli-demo.fasta EcoRI BamHI

# full suite inside the container
docker compose --profile ci run --rm test
```

`scripts/docker-verify.sh` builds every target, runs the suites, exercises the CLI, then serves the web image and probes each asset for a 200 and the right content type. The Qt Quick desktop host is not containerized — it needs a host GPU and display, so run it natively.

## Included workflows

- Adaptive map, sequence, features, tools, and provenance workspaces
- Circular plasmid and linear feature-map rendering with selection synchronization
- Sequence edits with undo/redo and feature-coordinate remapping
- FASTA, basic GenBank, FASTQ read, JSON snapshot, and validated `.ogs` project interchange boundaries
- Restriction-site scan, digest fragment table, virtual-gel visualization, and portable sticky/blunt-end planning
- Primer suggestions, exact IUPAC-aware PCR products, extensions, and site-directed-mutagenesis planning
- Restriction/ligation, Gibson/overlap, Golden Gate overhang, and oligo-annealing feasibility planners with deterministic provenance
- Six-frame ORF suggestions, translation, reverse translation, and transparent generic codon-design candidates
- Global/local alignment, read mapping, FASTQ/Phred consensus, substitutions/indels, and optional real-AB1 decoding boundary
- Non-mutating cloning operation records, candidate documents, graphical provenance, and verified project snapshots
- Touch-safe navigation, large controls, bottom sheets, keyboard parity, and reduced-motion support
- Portable Python core for immutable document models, commands, provenance DAGs, ZIP+SQLite projects, and scientific engines

## Project shape

```text
apps/desktop/         optional PySide6 / Qt Quick host
apps/cli/             automation and project-container CLI
packages/core/        stable document model, edits, history, projects
packages/io/          FASTA and GenBank interchange
packages/restriction/ independent restriction analysis
packages/primers/     simple primer/PCR utilities
packages/annotations/ translation, ORFs, feature suggestions, codon design
packages/alignment/   independently implemented global/local pairwise alignment
packages/sequencing/  FASTQ, consensus, optional AB1 boundary
packages/cloning/     assembly and ligation feasibility planners
packages/pcr/         PCR and mutagenesis planning
packages/api/         high-level `opengene` automation facade
src/                  responsive local browser workstation
tests/                CLI fixtures and browser wiring tests
```

## Scientific and interoperability notes

Internal coordinates are zero-based and half-open. Circular origin-spanning annotations are represented explicitly. The built-in enzyme list is deliberately small and independently maintained; production distribution should replace it with a separately audited, redistributable dataset. Primer scoring is an explainable heuristic when Primer3 is unavailable; it is not a substitute for experimental validation.

The `.ogs` project format is implemented as a verified ZIP + SQLite container with content-addressed sequence blobs and SHA-256 checksums. The dependency-free browser shell uses a browser-specific `.ogs.json` snapshot for portable local save/load; use the CLI or native Python layer for the binary `.ogs` container. Standard FASTA and GenBank remain the primary sequence interchange formats.

## Test it

```bash
# all portable packages, container, and CLI coverage
python3 tests/run_all.py

# dependency-free browser wiring check
npm run test:web

# everything, in containers
./scripts/docker-verify.sh
```

## Source style

The tree is deliberately comment-free: `scripts/strip_comments.py` removes `#` comments and docstrings from Python (keeping functional pragmas such as `# noqa` and `# type: ignore`), `//` and `/* */` from JavaScript, QML, and CSS, and `<!-- -->` from HTML. It is string- and regex-literal aware, and verifies every Python file's AST is unchanged before writing. Run `python3 scripts/strip_comments.py --check` to assert the tree is clean; the containerized test target runs it as a gate.

Use a qualified molecular biologist and independently licensed validation datasets before treating any plan, annotation, consensus, or compatibility path as release-grade laboratory evidence.

## Automation API

```python
from opengene import Project
from opengene.cloning import gibson_assembly

project = Project.open("example.ogs")
product = gibson_assembly(
    [project.sequence("backbone"), project.sequence("insert")],
    minimum_overlap=20,
)
project.add(product, name="assembled plasmid")
project.save()
```

The adapter creates a separate immutable output and turns its JSON-ready
assembly recipe into a provenance DAG node when it is added to the project.

## Clean-room commitment

This project is not affiliated with or endorsed by SnapGene. “SnapGene-compatible” should only ever describe independently implemented file interoperability after legal and compatibility review. Do not use proprietary assets, copy proprietary database content, or imply product affiliation.

See [feature-status.md](docs/feature-status.md) for the precise delivered
scope, integrations that remain optional, and the validation required before a
production laboratory release.

## Validation boundary

This is software for planning and documentation, not a guarantee of laboratory outcome. Validate primer designs, assembly plans, digests, and sequence verification with qualified scientific review and appropriate laboratory controls.
