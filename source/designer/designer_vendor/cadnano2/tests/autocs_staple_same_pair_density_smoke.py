"""Regression checks for Square same-helix-pair AutoCS density rules.

Load only the three pure functions under test from ``part.py``.  This keeps
the smoke test independent of a display server and of the Qt architecture.
"""

import ast
from collections import defaultdict
import os


PART_PATH = os.path.join(os.path.dirname(os.path.dirname(__file__)),
                         "model", "parts", "part.py")
with open(PART_PATH, "r", encoding="utf-8") as source:
    tree = ast.parse(source.read(), PART_PATH)
selected = [node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name in (
                "_circularIndexDistance",
                "_indexInsideExactPairInterval",
                "_squareStapleCrossoverIsOverdense")]
namespace = {"_SCAFFOLD_STAPLE_PAIR_EXCLUSION_DISTANCE": 10}
exec(compile(ast.Module(body=selected, type_ignores=[]), PART_PATH, "exec"),
     namespace)
_squareStapleCrossoverIsOverdense = namespace[
    "_squareStapleCrossoverIsOverdense"]


class _SquarePart(object):
    _step = 32


def blocked(index, pair_positions, pair=(0, 1), circular_size=None):
    positions = defaultdict(set)
    for helix_pair, values in pair_positions.items():
        positions[tuple(sorted(helix_pair))].update(values)
    return _squareStapleCrossoverIsOverdense(
        _SquarePart(), pair[0], pair[1], index, defaultdict(set),
        circular_size, defaultdict(set), positions)


# Exact 32-bp scaffold intervals suppress every internal staple site,
# regardless of whether either endpoint is an 8-multiple.
assert blocked(21, {(0, 1): (5, 36)})
assert blocked(16, {(0, 1): (0, 31)})

# At 64 bp, only the normal <10-bp endpoint exclusion applies.  The central
# legal site remains available so sparse 1/64 scaffold routing can recover
# the normal 1/32 total density.
assert blocked(8, {(0, 1): (0, 63)})
assert not blocked(31, {(0, 1): (0, 63)})
assert not blocked(32, {(0, 1): (0, 63)})

# A scaffold crossover belonging to another helix pair never affects A-B.
assert not blocked(16, {(0, 2): (0, 31)}, pair=(0, 1))

# Circular intervals obey the same rule, including the wraparound interval.
assert blocked(4, {(0, 1): (80, 15)}, circular_size=96)
assert not blocked(48, {(0, 1): (16, 79)}, circular_size=96)

# Exact coordinates from the reported 1.json: the 0--1 scaffold pair has
# two inclusive 32-bp modules and both midpoint AutoCS sites must be blocked.
assert blocked(15, {(0, 1): (0, 31, 32, 63)})
assert blocked(16, {(0, 1): (0, 31, 32, 63)})
assert blocked(47, {(0, 1): (0, 31, 32, 63)})
assert blocked(48, {(0, 1): (0, 31, 32, 63)})

print("PASS: Square AutoCS same-pair 32/64-bp density rules")
