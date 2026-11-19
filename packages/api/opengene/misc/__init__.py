"""OpenGene Miscellaneous Submodules."""

from __future__ import annotations

import sys
from . import assgn_two

# Primary module aliases
assgn = assgn_two
assgn2 = assgn_two

# Register hyphenated and shorthand sys.modules aliases
sys.modules["opengene.misc.assgn-two"] = assgn_two
sys.modules["opengene.misc.assgn2"] = assgn_two
sys.modules["opengene.misc.assgn"] = assgn_two
sys.modules["opengene.assgn-two"] = assgn_two
sys.modules["opengene.assgn_two"] = assgn_two
sys.modules["opengene.assgn2"] = assgn_two
sys.modules["opengene.assgn"] = assgn_two

# Direct functions
target_gene = assgn_two.target_gene
AssignmentResult = assgn_two.AssignmentResult

__all__ = ["assgn_two", "assgn", "assgn2", "target_gene", "AssignmentResult"]
