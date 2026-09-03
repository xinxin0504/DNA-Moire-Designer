"""Shared auxiliary-channel routing for ultra-small two-layer SST spacing.

The caDNAno canvas cannot store two nucleotides at the same
``(virtual-helix, base)`` coordinate.  At Z2=0/8 the two otherwise valid SST
layers can overlap at a small number of boundary nodes.  This module keeps
the accepted polymer topology and absolute phase, but draws only those
colliding layer-2 nodes on h64--79 (internal and public numbering).  Direct
primary/auxiliary links make the real 5'->3' path explicit for sequence
assignment.  Geometry exporters may use ``ideal_helix`` to project those
nodes back onto h48--63.
"""

from __future__ import annotations

from copy import deepcopy
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


EMPTY = [-1, -1, -1, -1]
PRIMARY_INTERNAL = tuple(range(48, 64))
AUXILIARY_INTERNAL = tuple(range(64, 80))
AUXILIARY_SHIFT = 16


def _nonempty(record: Sequence[int]) -> bool:
    return list(record) != EMPTY


def occupied_nodes(layer: Mapping[str, Mapping[int, Sequence[Sequence[int]]]]
                   ) -> set[Tuple[int, int]]:
    return {
        (int(helix), int(base))
        for field in ("scaf", "stap")
        for helix, records in layer[field].items()
        for base, record in enumerate(records)
        if _nonempty(record)
    }


def _components(field_rows: Mapping[int, Sequence[Sequence[int]]]
                ) -> List[set[Tuple[int, int]]]:
    nodes = {
        (int(helix), int(base))
        for helix, records in field_rows.items()
        for base, record in enumerate(records)
        if _nonempty(record)
    }
    unseen = set(nodes)
    output = []
    while unseen:
        start = next(iter(unseen))
        component = set()
        stack = [start]
        while stack:
            node = stack.pop()
            if node in component:
                continue
            component.add(node)
            helix, base = node
            record = field_rows[helix][base]
            for slot in (0, 2):
                partner = tuple(map(int, record[slot:slot + 2]))
                if partner in nodes and partner not in component:
                    stack.append(partner)
        unseen -= component
        output.append(component)
    return output


def _compress_bases(values: Iterable[int]) -> List[List[int]]:
    values = sorted(set(map(int, values)))
    if not values:
        return []
    output = []
    low = high = values[0]
    for value in values[1:]:
        if value == high + 1:
            high = value
        else:
            output.append([low, high])
            low = high = value
    output.append([low, high])
    return output


def route_layer2_conflicts(layer1, layer2, array_length: int):
    """Return merged-primary fields, auxiliary fields and audit metadata."""
    occupied1 = occupied_nodes(layer1)
    destination = {}
    for field in ("scaf", "stap"):
        for helix, records in layer2[field].items():
            for base, record in enumerate(records):
                if not _nonempty(record):
                    continue
                actual = (int(helix) + AUXILIARY_SHIFT
                          if (int(helix), int(base)) in occupied1
                          else int(helix))
                destination[(field, int(helix), int(base))] = (
                    actual, int(base))

    primary2 = {
        field: {helix: [EMPTY[:] for unused in range(array_length)]
                for helix in PRIMARY_INTERNAL}
        for field in ("scaf", "stap")
    }
    auxiliary = {
        field: {helix: [EMPTY[:] for unused in range(array_length)]
                for helix in AUXILIARY_INTERNAL}
        for field in ("scaf", "stap")
    }
    detours = []
    for field in ("scaf", "stap"):
        for component_index, component in enumerate(
                _components(layer2[field])):
            collisions = sorted(component & occupied1)
            if not collisions:
                continue
            boundary_edges = set()
            for helix, base in component:
                here_aux = (helix, base) in occupied1
                record = layer2[field][helix][base]
                for slot in (0, 2):
                    partner = tuple(map(int, record[slot:slot + 2]))
                    if partner in component and \
                            here_aux != (partner in occupied1):
                        boundary_edges.add(tuple(sorted(
                            ((helix, base), partner))))
            detours.append({
                "field": field,
                "component_index": int(component_index),
                "component_length_nt": len(component),
                "detoured_node_count": len(collisions),
                "detoured_helices": sorted({h for h, unused in collisions}),
                "detoured_base_range": [
                    min(base for unused, base in collisions),
                    max(base for unused, base in collisions)],
                "primary_auxiliary_boundary_links": [
                    [list(first), list(second)]
                    for first, second in sorted(boundary_edges)],
            })

        for helix, records in layer2[field].items():
            for base, source_record in enumerate(records):
                if not _nonempty(source_record):
                    continue
                actual_helix, actual_base = destination[
                    (field, int(helix), int(base))]
                translated = list(source_record)
                for slot in (0, 2):
                    partner, partner_base = map(
                        int, translated[slot:slot + 2])
                    if partner < 0:
                        continue
                    target = destination.get(
                        (field, partner, partner_base))
                    if target is None:
                        raise ValueError(
                            "SST辅助绕行遇到缺失的layer-2互反端点。")
                    translated[slot:slot + 2] = list(target)
                target_rows = (auxiliary if actual_helix >= 64 else primary2)
                target_rows[field][actual_helix][actual_base] = translated

    merged = deepcopy(layer1)
    for field in ("scaf", "stap"):
        for helix, records in primary2[field].items():
            for base, record in enumerate(records):
                if not _nonempty(record):
                    continue
                if _nonempty(merged[field][helix][base]):
                    raise ValueError("SST辅助绕行后仍有未解决的主通道占位冲突。")
                merged[field][helix][base] = list(record)

    intervals = {}
    for logical in PRIMARY_INTERNAL:
        bases = {
            base for field in ("scaf", "stap")
            for base, record in enumerate(
                auxiliary[field][logical + AUXILIARY_SHIFT])
            if _nonempty(record)}
        if bases:
            intervals[str(logical)] = _compress_bases(bases)
    metadata = {
        "enabled": bool(intervals),
        "policy": (
            "only colliding layer-2 nodes detour to h64-79; direct links "
            "preserve one continuous caDNAno 5prime-to-3prime path"),
        "layer": 2,
        "primary_internal_helices": list(PRIMARY_INTERNAL),
        "auxiliary_internal_helices": list(AUXILIARY_INTERNAL),
        "auxiliary_to_ideal_internal_helix": {
            str(number): number - AUXILIARY_SHIFT
            for number in AUXILIARY_INTERNAL},
        "detoured_intervals_by_logical_helix": intervals,
        "detoured_components": detours,
        "geometry_export_policy": (
            "project auxiliary h64-79 back to ideal h48-63; retain actual "
            "caDNAno helix/base ids in sequence exports"),
    }
    return merged, auxiliary, destination, metadata


def actual_helix(layout: Mapping[str, object], layer: int, field: str,
                 logical_helix: int, base: int) -> int:
    """Resolve a logical primary SST node to its caDNAno display helix."""
    routing = layout.get("auxiliary_sst_routing", {})
    if int(layer) != 2 or not isinstance(routing, Mapping) or not \
            routing.get("enabled"):
        return int(logical_helix)
    intervals = routing.get(
        "detoured_intervals_by_logical_helix", {}).get(
            str(int(logical_helix)), [])
    if any(int(low) <= int(base) <= int(high) for low, high in intervals):
        return int(logical_helix) + AUXILIARY_SHIFT
    return int(logical_helix)


def ideal_helix(metadata: Mapping[str, object], actual: int) -> int:
    """Return the ideal 3D helix for one real caDNAno auxiliary helix."""
    routing = metadata.get("auxiliary_sst_routing", {})
    mapping = (routing.get("auxiliary_to_ideal_internal_helix", {})
               if isinstance(routing, Mapping) else {})
    return int(mapping.get(str(int(actual)), int(actual)))

