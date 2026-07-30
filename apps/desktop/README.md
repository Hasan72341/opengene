# OpenGene Desktop Host

This folder contains an optional, native-ready PySide6 / Qt Quick host for
OpenGene Studio. It is deliberately a UI scaffold, not a second scientific
implementation: sequence, project, edit-command, restriction, primer, and
provenance logic belong in the portable Python packages.

## What is implemented

- An adaptive Qt Quick application shell with a desktop workspace rail and a
  touch-friendly bottom bar below 900 dp.
- 48–52 dp interactive targets, readable compact layouts, status feedback,
  modal help, and accessible navigation labels.
- A native map-preview component with pinch and explicit zoom controls.
- Clearly marked controller boundaries for map, sequence, annotations, tools,
  and construct history.
- A launcher that never crashes merely because the optional Qt runtime is not
  installed.

The sample map, sequence, features, and cards are presentation fixtures only.
They do **not** claim to run molecular-biology calculations or modify project
data.

## Run the scaffold

Install it as a separate optional package from the repository root:

```bash
python -m pip install -e ./apps/desktop
opengene-desktop
```

Or, when working from a source checkout whose dependencies are already
installed:

```bash
python -m apps.desktop
```

`PySide6` is declared by this folder's standalone `pyproject.toml`. If it is
not installed, the command exits clearly and points to the working browser
workstation instead:

```bash
npm run dev
# open http://localhost:4174
```

Useful non-GUI checks:

```bash
python -m apps.desktop --check
python -m apps.desktop --qml-path
```

## Architecture boundary

`main.py` keeps the PySide6 imports lazy, injects a minimal `AppController`
into `Main.qml`, and refuses to start if the QML entry point is absent. The
controller currently only communicates UI status. Future work should add an
explicit service layer that exposes immutable document snapshots, commands,
background job state, and file-dialog requests; it should not put algorithms
or direct mutation in QML.

The QML source stays as regular files rather than generated resources during
this scaffolding phase so it can be reviewed and styled independently. The
standalone package metadata includes those files for both wheels and source
distributions.

## Clean-room and validation boundary

This interface uses original layout, text, components, and map drawing. It is
not affiliated with SnapGene and makes no compatibility claim. Before a
scientific feature is made actionable, wire it to a validated portable engine,
add tests with independently licensed data, and preserve all parameters and
warnings in provenance.
