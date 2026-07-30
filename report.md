# OpenGene Studio — Capability Report and SnapGene Gap Analysis

**Repository:** `opengene-studio` v0.1.0 · **Report date:** 2026-07-30
**Basis:** direct source audit of this working tree (7,712 LOC of library Python, 1,073 LOC of tests, 2,039 LOC of browser application, 61 test functions, all suites passing).

---

## 1. Executive summary

OpenGene Studio is a clean-room, local-first molecular biology workstation. It is a **credible alpha with an unusually strong architectural core and a deliberately thin scientific surface.**

The honest headline: **the repository implements roughly a quarter of SnapGene's functional surface, but it implements that quarter with more rigor than SnapGene exposes.** The immutable document model, the provenance DAG with cycle rejection, the checksum-validated `.ogs` container, and the deterministic command/undo system are genuinely well-built — closer to what a regulated lab or a reproducible-research pipeline needs than SnapGene's own file format offers.

What it is *not* close to is SnapGene's **data breadth and biophysical fidelity**. SnapGene's practical value comes substantially from ~700 commercial restriction enzymes, thousands of curated auto-annotation features, nearest-neighbor thermodynamics, real chromatogram rendering, and a decade of file-format compatibility. This repository ships **8 enzymes, 18 browser-side common features, a Wallace-rule Tm, no trace rendering, and 3 interchange formats.**

The gap is therefore **not primarily an engineering-difficulty gap — it is a data-licensing, curation, and validation gap.** That is good news for planning: most of the remaining work is well-understood, parallelizable, and bounded, but a meaningful portion of it is legal/curatorial rather than something that can be coded.

**Overall parity estimate: ~25–30% of SnapGene's functional surface; ~15% weighted by day-to-day bench usage** (because the missing enzyme and feature databases sit on the critical path of nearly every real workflow).

---

## 2. What is implemented

### 2.1 Architecture

Five layers, cleanly separated, with the domain logic entirely dependency-free:

| Layer | Location | Notes |
|---|---|---|
| Browser workstation | `index.html`, `src/app.js`, `src/styles.css` | 2,039 LOC, zero dependencies, offline-capable (service worker + manifest) |
| Portable domain packages | `packages/*` | 10 packages, pure standard library |
| Automation façade | `packages/api/opengene` | High-level `Project` API with automatic provenance capture |
| CLI | `apps/cli/main.py` | 10 commands, JSON output |
| Native host | `apps/desktop` | Qt Quick / PySide6 — **presentation shell only** |
| Containers | `Dockerfile`, `docker-compose.yml` | `web`, `cli`, `test` targets |

### 2.2 Core data model — `opengene_core`

Genuinely the strongest part of the codebase.

- Immutable `SequenceDocument`, `Feature`, `Location` with frozen/thawed value semantics.
- `MoleculeType`: DNA, RNA, **PROTEIN** — protein documents are modelled, though tooling around them is thin.
- `Topology`: linear and circular, with correct origin-spanning location handling (`wraps_origin`, `start > end`).
- Zero-based, half-open `[start, end)` coordinates throughout, with conversion isolated at the GenBank import/export boundary. This is a deliberate and correct decision that avoids an entire class of off-by-one bugs.
- 12 command types with full undo/redo via `CommandHistory`: `InsertBases`, `DeleteBases`, `ReplaceRegion`, `ReverseComplement`, `RotateOrigin`, `ChangeTopology`, `TranslateSequence`, `AddFeature`, `UpdateFeature`, `RemoveFeature`.
- Feature coordinate remapping survives edits.
- `Project` + `ProvenanceNode` form a **directed acyclic graph with active cycle rejection** (`_validate_provenance_dag`). Operations never overwrite inputs. Nodes carry input IDs, output ID, parameters, engine version, warnings, sequence hashes, and a `reproducibility_status` flag that is set to non-reproducible rather than silently guessing when plugin or database requirements are unresolved.

**This exceeds SnapGene.** SnapGene records a linear-ish history; it does not expose a validated provenance DAG with per-operation engine versioning and explicit reproducibility status.

### 2.3 File formats — `opengene_io`

- **FASTA** — read/write/parse.
- **GenBank** — read/write, with a real location-expression parser: `join()`, `complement()`, nested/mixed-strand operators, origin-wrap merging, qualifier continuation lines, LOCUS metadata extraction, 79-column output wrapping.
- **FASTQ** — parse/write with Phred encode/decode.
- **`.ogs`** — the native container: ZIP + SQLite, deterministic member names, content-addressed SHA-256 sequence blobs (identical sequences stored once), `manifest.json`, `checksums.json`, atomic replace-on-write. The reader rejects unknown format versions, unsafe ZIP member paths (zip-slip), checksum failures, malformed SQLite structure, and manifest/blob/row discrepancies.
- **`.ogs.json`** — a browser-side portable snapshot, because the browser shell has no ZIP/SQLite dependency.

**The `.ogs` container is more rigorous than SnapGene's `.dna`** in verifiability and tamper-evidence. It is, of course, understood by nothing else in the world.

### 2.4 Restriction analysis — `opengene_restriction`

- Site scanning with IUPAC ambiguity, circular-aware wrap detection.
- `digest()` producing fragments with correct circular/linear semantics.
- Unique-cutter identification (`find_unique_cutters`).
- Enzyme model carries recognition site, cut offset, and complementary cut offset — so **overhang polarity (5′ / 3′ / blunt) is correctly modelled**, not just cut position.
- **Built-in enzymes: 8** — EcoRI, BamHI, HindIII, XhoI, NotI, PstI, SmaI, BsaI.

### 2.5 Cloning — `opengene_cloning`

Six planners, all returning structured feasibility plans rather than asserted outcomes:

- `plan_restriction_digest`, `plan_restriction_ligation` — cohesive/blunt end compatibility, junction modelling.
- `plan_gibson_assembly`, `plan_overlap_assembly` — exact-overlap detection with junction records.
- `plan_golden_gate` — Type IIS, junction/fusion-site analysis.
- `plan_oligo_annealing` — oligo duplex planning with alignment records.

Each emits an `OperationRecord` suitable for provenance attachment.

### 2.6 PCR and mutagenesis — `opengene_pcr`

- `plan_pcr` — exact IUPAC-aware binding-site location, product generation, circular template support.
- `plan_mutagenesis`, `plan_site_directed_mutagenesis` — mutation application with explicit `SequenceMutation` records.
- `find_primer_binding_sites`.
- `PCRProvenanceRecord` for DAG attachment.

### 2.7 Primer design — `opengene_primers`

- `design_primers` with a real parameter surface: target window, product size range, primer length range, Tm range, GC range, search window, 5′ extensions (for cloning tails), max pairs.
- Heuristic screens including homopolymer-run rejection.
- `simulate_pcr`, `melting_temperature`, `gc_content`.

### 2.8 Annotation and translation — `opengene_annotations`

- Six-frame ORF detection (`find_open_reading_frames`) with circular awareness.
- Translation with the standard genetic code; deterministic reverse translation.
- **Codon optimization** (`optimize_codon_usage`) with pluggable codon-usage tables, forbidden-motif avoidance, and warning emission. Notably, this is a capability **SnapGene does not offer natively**.
- `detect_features` — currently an alias for `detect_orf_features`, i.e. ORF-based only.

### 2.9 Alignment — `opengene_alignment`

- Needleman–Wunsch (global) and Smith–Waterman (local), full dynamic programming with traceback.
- Configurable match/mismatch/gap scores; optional IUPAC ambiguity matching.
- CIGAR string generation.
- `variants_from_alignment` — typed variant extraction (substitution / insertion / deletion) with reference and query coordinate spans.

### 2.10 Sequencing — `opengene_sequencing`

- FASTQ parse/write, Phred encode/decode.
- Reverse-complement-aware read mapping with correctly oriented quality coordinates.
- Quality-weighted consensus calling with IUPAC output for conflicting high-quality observations, explicit uncovered-base marking, and explicit insertion records.
- **AB1 reading** via an optional Biopython delegation, with a clean capability-probe API (`ab1_reader_status`, `require_ab1_reader`, `Ab1UnavailableError`) — it degrades honestly instead of failing opaquely.

### 2.11 Browser workstation

~130 functions covering: circular map with interactive zoom/pan, linear map, sequence view with touch selection, features table and detail, tools panel (restriction, primer, cloning, gel, verification, annotation), project library, history/provenance view, undo/redo, FASTA/GenBank/FASTQ import, FASTA/GenBank export, local autosave to browser storage, codon-optimization modal, assembly modal with Golden Gate cut handling, read inspection, variant focus, toast notifications, and an icon sprite system.

### 2.12 CLI

`info`, `digest`, `align`, `orfs`, `pcr`, `gibson`, `verify`, `project create`, `project info`.

### 2.13 Infrastructure

- Comment-free source tree enforced by `scripts/strip_comments.py --check`, run as a CI gate.
- Three Docker targets; `scripts/docker-verify.sh` builds all, runs both suites, exercises CLI workflows, and probes every served asset for status and content type, gated on container health.

---

## 3. Comparison to SnapGene

Parity scale: **●●●** at or above parity · **●●○** partial · **●○○** token/demo · **○○○** absent.

> SnapGene figures below are approximate, drawn from public product documentation. They should be re-verified before being used in any external-facing document.

### 3.1 Where OpenGene meets or beats SnapGene

| Capability | Parity | Notes |
|---|---|---|
| Provenance / reproducibility | ●●● **beats** | Validated DAG, engine versioning, explicit reproducibility status. SnapGene has no equivalent. |
| File-format integrity | ●●● **beats** | Content-addressed blobs, SHA-256 checksums, zip-slip rejection, atomic writes. |
| Coordinate-model discipline | ●●● **beats** | Uniform zero-based half-open internally; conversion isolated to I/O. |
| Scriptability / automation | ●●● **beats** | First-class Python API + CLI with JSON output. SnapGene's automation is limited and largely GUI-bound. |
| Codon optimization | ●●● **beats** | Not a native SnapGene feature. |
| Local-first / privacy | ●●● **beats** | Zero network egress by construction; no account required. |
| Licensing | ●●● **beats** | MPL-2.0 vs commercial per-seat. |
| Undo/redo command model | ●●○ | 12 typed commands with coordinate remapping; SnapGene's is broader in UI reach. |

### 3.2 Where OpenGene is partial

| Capability | Parity | Gap |
|---|---|---|
| Circular/linear map view | ●●○ | Renders well; lacks SnapGene's label collision avoidance, print/high-res export, and layout polish. |
| Sequence editing | ●●○ | Core edits present; no find/replace across documents, no bulk operations, no ORF/translation track overlays in-editor. |
| Restriction analysis | ●●○ | Algorithms are sound; **8 enzymes vs ~700**. No enzyme-set management, no supplier/isoschizomer/methylation-sensitivity data. |
| Auto-annotation | ●●○ | Mechanism works; **18 browser features, ORF-only in the Python layer, vs SnapGene's curated thousands**. |
| Primer design | ●●○ | Reasonable parameter surface; **Wallace/basic Tm only — no nearest-neighbor thermodynamics, no salt/Mg²⁺/dNTP corrections, no hairpin/dimer ΔG**. |
| Cloning planning | ●●○ | Six methods with real end-compatibility logic; missing Gateway, TOPO, In-Fusion, TA/blunt, MultiSite. |
| Sanger verification | ●●○ | Mapping and consensus are solid and quality-aware; **no chromatogram rendering at all**. |
| Alignment | ●●○ | Correct DP with variant extraction, but **linear gap penalty only (no affine), O(n·m) memory, pairwise only**. |
| GenBank support | ●●○ | Parser is respectable; not validated against the long tail of real-world malformed files. |
| Virtual gel | ●○○ | Log-scale size-to-position preview, hardcoded to EcoRI+BamHI+HindIII. Not a migration model. |
| Native desktop app | ●○○ | **Presentation shell only** — the Qt controller exposes app name, version, core status, and a status message. No editing. |

### 3.3 Where OpenGene has nothing

| Capability | SnapGene | OpenGene |
|---|---|---|
| Enzyme database (~700 commercial, REBASE-derived) | ●●● | ○○○ |
| Curated common-features database (thousands) | ●●● | ○○○ |
| `.dna` native SnapGene format | ●●● | ○○○ (deliberate — clean-room boundary) |
| EMBL, GFF/GTF, BED, SBOL, ApE, Vector NTI, Clone Manager, DNA Strider, MacVector | ●●● | ○○○ |
| `.scf` chromatogram formats | ●●● | ○○○ |
| Chromatogram trace rendering | ●●● | ○○○ |
| CRISPR guide RNA design and off-target scanning | ●●● | ○○○ |
| NCBI / online database import | ●●● | ○○○ |
| Multiple sequence alignment (MAFFT/MUSCLE-class) | ●●● | ○○○ |
| BLAST / genome-wide off-target search | ●●● | ○○○ |
| Contig assembly | ●●● | ○○○ |
| High-resolution image / PDF / print export | ●●● | ○○○ |
| Shared plasmid libraries, team/cloud collaboration | ●●● | ○○○ |
| Primer database and ordering integrations | ●●● | ○○○ |
| Gateway / TOPO / In-Fusion / TA cloning | ●●● | ○○○ |
| Protein analysis tooling (beyond the data model) | ●●○ | ○○○ |
| Agarose gel simulation with realistic migration, ladders, dyes | ●●● | ●○○ |
| Enzyme set management, methylation sensitivity, star activity | ●●● | ○○○ |
| LIMS / ELN connectors, inventory | ●●○ | ○○○ |

---

## 4. Depth gaps inside features that *look* finished

This is the most important section for planning, because these are the items most likely to be mistaken for done.

1. **Melting temperature is not usable for real primer design.** `melting_temperature()` offers only the Wallace rule (`2AT + 4GC`) and the basic GC formula. Both diverge substantially from measured Tm for primers above ~20 nt and neither accounts for salt, Mg²⁺, or dNTP concentration. Every primer, PCR, mutagenesis, and Gibson plan that reports a Tm is reporting a number that should not drive a bench decision. **Nearest-neighbor thermodynamics (SantaLucia unified parameters) with salt correction is the single highest-value scientific fix in the repository.**

2. **No secondary-structure screening.** There is homopolymer rejection, but no hairpin, self-dimer, or hetero-dimer ΔG calculation. Primer pairs can be emitted that will not work.

3. **Alignment uses a linear gap penalty.** Real biological alignment needs affine gaps (Gotoh). Current output will over-fragment indels into many short gaps rather than one long one, which directly distorts `variants_from_alignment` results.

4. **Alignment is O(n·m) in memory.** Two 100 kb sequences would require ~10¹⁰ cells. There is no banding, no Hirschberg linear-space variant, no seed-and-extend. **Plasmid-scale works; genome-scale will exhaust memory.**

5. **`detect_features` is ORF detection wearing a general name.** `detect_features()` is a direct alias for `detect_orf_features()`. The name promises motif/database annotation that does not exist in the Python layer. The 18-entry feature database lives only in the browser JavaScript and is not shared with the Python packages — **the two layers will disagree about what a sequence contains.**

6. **The browser and Python layers duplicate science.** Tm, ORF finding, restriction scanning, and the genetic code are implemented twice, in two languages, with different data. This is a correctness hazard, not just duplication: results shown in the UI are not the results the CLI produces.

7. **The virtual gel is a visualization, not a model.** Band position is `log10`-interpolated; the gel-percentage selector does not affect migration; the enzyme set is hardcoded to three enzymes and ignores user selection.

8. **AB1 support is read-only metadata extraction.** `read_ab1` pulls sequence and Phred scores via Biopython but discards the trace channels entirely. Chromatogram review — a primary reason labs open SnapGene — is impossible.

9. **The Qt desktop host does nothing.** The controller exposes four read-only properties and one `announce` slot. It cannot open, edit, or display a sequence.

10. **Test coverage is thin relative to claims.** 61 test functions / 40 unittest cases across ~7,700 LOC of scientific code. There are no property-based tests, no round-trip fuzzing of the GenBank parser, no malformed-input corpus, and no numerical validation against reference implementations.

---

## 5. Roadmap to meaningful parity

Ordered by value-per-unit-effort.

### Phase 1 — Make the existing science trustworthy (highest value, low risk)
1. Nearest-neighbor Tm with salt/Mg²⁺ correction; deprecate Wallace to an explicit opt-in.
2. Hairpin and dimer ΔG screening in `design_primers`.
3. Affine gap penalties (Gotoh) in both alignment modes.
4. Unify the browser and Python science layers — compile the Python core to WebAssembly, or reduce the browser to a rendering layer over a single shared implementation. This eliminates an entire class of "the UI says X, the CLI says Y" bugs.
5. Numerical validation suite against published reference values.

### Phase 2 — Close the data gap (highest user-visible value, mostly non-engineering)
6. **REBASE-derived enzyme database.** REBASE is free for academic use; commercial redistribution requires review. This is a licensing task first and an ingestion task second. Target: 500+ enzymes with supplier, isoschizomer, and methylation-sensitivity data.
7. **Curated common-features database.** Must be built independently — the clean-room boundary in `docs/clean-room.md` explicitly forbids porting SnapGene's. Realistically 1,000+ curated entries; this is months of curation, not a sprint.
8. Sequence-motif and database-driven `detect_features` built on (7).
9. Enzyme set management in the UI.

### Phase 3 — Close the format gap
10. EMBL, GFF3, BED — straightforward, high compatibility payoff.
11. `.scf` and full AB1 trace channel extraction.
12. SBOL3 for synthetic-biology interoperability.
13. **`.dna` import remains a legal question, not a technical one.** Treat as out of scope until counsel advises otherwise.

### Phase 4 — Close the workflow gap
14. Chromatogram trace rendering (depends on 11) — this is the single biggest missing *UI* capability.
15. CRISPR guide design with on/off-target scoring.
16. Multiple sequence alignment.
17. Gateway, TOPO, In-Fusion, TA/blunt cloning planners.
18. High-resolution image and PDF export.
19. Contig assembly.

### Phase 5 — Product surface
20. Make the Qt host real, or formally deprecate it and commit to the browser shell. **Maintaining a non-functional second UI is a net negative today.**
21. NCBI import, shared libraries, primer databases, ordering integrations.

---

## 6. Metrics

| Metric | Value |
|---|---|
| Library Python LOC | 7,712 |
| Test LOC | 1,073 |
| Test functions / unittest cases | 61 / 40 |
| Browser application LOC | 2,039 (`app.js`) + 835 (CSS) + 366 (HTML) |
| Domain packages | 10 |
| CLI commands | 10 |
| Runtime dependencies | **0** (Biopython optional, AB1 only) |
| Built-in restriction enzymes | 8 |
| Browser common-features entries | 18 |
| Interchange formats read | 3 (FASTA, GenBank, FASTQ) + 2 native |
| Cloning methods planned | 6 |
| Edit command types | 12 |
| Docker targets | 3 |
| Comments/docstrings in source | 0 (enforced) |

---

## 7. Verdict

**As a product competing with SnapGene today: not close.** A bench scientist cannot switch, because 8 enzymes, 18 features, an approximate Tm, and no chromatogram view do not support daily work.

**As a foundation to build that product on: unusually good.** The data model, provenance system, container format, and coordinate discipline are the parts that are expensive to retrofit and painful to get wrong, and they are already right. The remaining work is largely additive — databases, formats, thermodynamics, and UI surface — rather than structural.

The critical strategic observation is that **the largest single blocker is curatorial and legal, not technical.** The enzyme and feature databases gate the majority of real workflows, cannot be copied from SnapGene under the project's own clean-room policy, and will take longer to build correctly than any code in the roadmap. **That work should start now, in parallel with Phase 1, rather than being sequenced after it.**

Finally, the repository's own `docs/feature-status.md` is commendably honest about these boundaries. That honesty is a genuine asset and should be preserved — nothing in this report contradicts it; this document quantifies it.

---

## 8. Method and caveats

- Findings are derived from direct source inspection of this working tree, not from documentation claims. Where `docs/feature-status.md` and the code disagreed, the code was taken as authoritative.
- All test suites were executed (natively and containerized) and pass at the time of writing.
- **SnapGene capability figures are approximate**, drawn from public product documentation and general familiarity. They are directionally reliable but should be independently verified before external use.
- Parity percentages are informed judgment, not measurement. They are intended for prioritization, not for marketing.
- No SnapGene source, database, or asset was consulted or referenced in producing this report; the comparison is behavioral and based on publicly documented capability.
