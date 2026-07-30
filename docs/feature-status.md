# Delivery status and validation boundary

OpenGene Studio is an independently designed, local-first **alpha**. The
items below distinguish runnable functionality from intentionally unclaimed
scientific or compatibility scope.

## Implemented in this repository

- Touch-first browser workstation with synchronized map, linear map, sequence,
  features, tools, project library, and provenance views.
- Local sequence edits, undo/redo, coordinate remapping, circular features,
  FASTA/basic GenBank import/export, FASTQ import, and browser-local snapshots.
- Portable immutable Python document and provenance-DAG model.
- Validated `.ogs` ZIP + SQLite container with checksum validation and
  content-addressed sequence blobs.
- Restriction scanning/digest primitives, virtual-gel size visualisation,
  primer heuristics, exact PCR/mutagenesis planning, and cloning feasibility
  plans for restriction/ligation, Gibson/overlap, Golden Gate and oligo
  annealing.
- Six-frame ORF detection, translation, deterministic reverse translation and
  generic codon optimisation.
- Pairwise global/local alignment, FASTQ/Phred handling, reverse-aware read
  mapping, reference-anchored consensus and explicit indel observations.
- CLI and `opengene` automation façade with provenance attachment.

## Deliberately optional or not yet a release claim

- Primer3 thermodynamics, minimap2/Parasail/MAFFT/MUSCLE integrations,
  genome-wide off-target search, AB1 trace rendering, and contig assembly
  require separately installed/validated engines and a dedicated UI binding.
- EMBL, GFF/GTF, BED, SBOL3, ApE, SnapGene `.dna`, and comprehensive
  chromatogram/alignment-format compatibility are not claimed by this build.
- The optional PySide6 host is a touch-ready native presentation shell; the
  runnable feature-complete touch workflow today is the browser workstation.
- Enzyme and feature datasets are intentionally small and independently
  maintained. Production packaging requires redistribution and data-quality
  review.
- Plugin isolation, collaboration service, LIMS/ELN connectors, inventory,
  ordering integrations and regulated-sequence workflow review remain product
  extensions, not hidden placeholders.

## Required before a production laboratory release

1. Domain-expert validation against independently licensed reference datasets
   and wet-lab controls.
2. Dependency, data, trademark and clean-room legal audit.
3. Compatibility testing for each advertised file format, including malformed
   and round-trip fixtures users are entitled to contribute.
4. Platform packaging, accessibility, backup/recovery, security and update
   testing for the intended deployment environments.

No in-silico plan, consensus, primer pair, digest, or assembly output should
be represented as a laboratory guarantee.
