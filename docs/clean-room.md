# Clean-room implementation boundary

OpenGene Studio is an independently designed molecular-biology workstation.
The project may support common workflows and standard scientific interchange
formats, but it must not copy or port SnapGene source code, documentation,
icons, proprietary UI, proprietary enzyme/feature databases, or protected
file-format implementation details.

Compatibility work is isolated at import/export boundaries. A compatibility
test corpus must contain only files that contributors are entitled to use.
Unknown binary blocks should be treated as opaque user data; their presence is
not permission to infer, copy, or distribute proprietary semantics.

Scientific dependencies and bundled data each need a license review. An
open-source library does not automatically grant redistribution rights for its
datasets. Production releases should publish a dependency and data-attribution
manifest alongside the build.
