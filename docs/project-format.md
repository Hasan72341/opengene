# OpenGene Project Format (`.ogs`)

## Status

Version 1 is implemented by the dependency-free `opengene_io` package and
the `opengene project create` / `opengene project info` CLI commands. The
runnable browser workstation exports a browser-specific JSON snapshot with the
extension `.ogs.json` because it deliberately has no bundled ZIP or SQLite
dependency. FASTA and GenBank remain the primary sequence interchange formats.

## Container

An `.ogs` file is a ZIP container with deterministic member names:

```text
project.ogs
├── manifest.json
├── project.sqlite
├── sequences/
│   └── sha256-<lowercase hex>.seq
└── checksums.json
```

Future versions may add `attachments/` and `previews/`; readers must not
assume those optional directories exist in version 1.

All textual members use UTF-8. Sequence blobs are normalized uppercase IUPAC
text without line wrapping. Hashes are SHA-256 of the exact stored bytes.

## Coordinates

- Coordinates in SQLite and JSON are zero-based, half-open `[start, end)`.
- Linear locations require `start <= end`.
- Circular origin-spanning locations require `wraps_origin: true` and
  `start > end`; they represent `[start, sequence_length)` plus `[0, end)`.
- Imports and exports alone convert GenBank’s one-based inclusive notation.

## `manifest.json`

```json
{
  "format": "OpenGene Studio project",
  "formatVersion": 1,
  "applicationVersion": "0.1.0",
  "createdAt": "2026-07-23T00:00:00Z",
  "hashAlgorithm": "sha256",
  "projectId": "uuid",
  "projectSchemaVersion": 1,
  "documentIds": ["uuid"],
  "sequenceMembers": {
    "uuid": "sequences/sha256-<lowercase hex>.seq"
  },
  "pluginRequirements": [],
  "databaseVersions": {}
}
```

Readers reject unknown required format versions, unsafe ZIP member names,
checksum failures, malformed SQLite structure, and discrepancies between the
manifest, sequence blobs, and SQLite project rows. The writer stores a
content-addressed sequence blob once even when multiple documents have
identical sequence content, and atomically replaces the target only after a
complete archive has been built.

## Provenance

Each operation is a directed-acryclic-graph node containing its input document
IDs, output document ID, parameters, engine/version identifiers, warnings, and
sequence hashes. Operations never overwrite an input document. An operation
with unresolved plugin or database requirements is marked non-reproducible,
not silently reinterpreted.

## Interoperability and clean room

The format is independently specified. It does not reuse an application’s
proprietary binary layout, visual assets, terminology, or biological database.
Any optional third-party file importer belongs at the interchange boundary and
must preserve unknown source data only when that data was legally supplied by
the user.
