"""Process-isolated entry point for the bundled DNAxiS routing engine."""

import base64
import gzip
import hashlib
import json
import math
import os
import random
import sys


ROOT = os.path.abspath(os.path.dirname(__file__))
os.environ["DNAXIS_HEADLESS"] = "1"
os.chdir(ROOT)
sys.path.insert(0, ROOT)

from app import config  # noqa: E402
from app.driver import driver  # noqa: E402
from app.routing.helper import log, dnaconnector  # noqa: E402
from app.routing import sequence  # noqa: E402
from app.shapeutils.custom import CustomInput  # noqa: E402


OXDNA_LENGTH_NM = 0.8518
CURVED_SCAFFOLD_MAX_BASES = 25000


class _VirtualScaffoldSequence(object):
    """Length-only scaffold used while Curved Design skips sequence export."""
    def __init__(self, length):
        self._length = max(1, int(length))

    def __len__(self):
        return self._length

    def __iter__(self):
        # Curved Design always calls the driver with skip_sequence=True.  Keep
        # iteration well-defined as a defensive fallback without allocating a
        # potentially very large string.
        for unused_index in range(self._length):
            yield "A"


_ORIGINAL_GETSEQ = sequence.getseq


def _install_curved_scaffold_capacity(required_bases):
    """Make DNAxiS routing independent of bundled sequence file lengths.

    The headless Curved Design route does not assign a DNA sequence.  DNAxiS
    nevertheless asks its sequence module for a scaffold length to decide
    whether to split the route.  Supplying one virtual length-only scaffold
    preserves a single route up to Curved Design's 25,000-base limit while
    leaving real sequence assignment to cadnano after the JSON has opened.
    """
    required_bases = int(required_bases)
    if required_bases > CURVED_SCAFFOLD_MAX_BASES:
        raise ValueError(
            "Curved design requires %d scaffold bases; the supported maximum "
            "is %d bases." %
            (required_bases, CURVED_SCAFFOLD_MAX_BASES))
    virtual_name = "curvedvirtual"
    virtual_sequence = _VirtualScaffoldSequence(required_bases)

    def curved_getseq(name="m13mp18"):
        if str(name) == virtual_name:
            return virtual_sequence
        return _ORIGINAL_GETSEQ(name)

    sequence.getseq = curved_getseq
    config.AVAIL_SEQUENCES = [virtual_name]
    return virtual_name


def _write_json(path, value):
    with open(path, "w", encoding="utf-8", newline="\n") as output:
        json.dump(value, output, separators=(",", ":"))


def _balanced_nominal_bp(target, period):
    """Choose a low-edit lattice circumference, then minimize local strain.

    First minimize the number of chemically disruptive indel sites.  Among
    equally small edits, use the candidate with the lowest worst and average
    count per lattice crossover period.  This is especially important for
    honeycomb: switching a 94-bp ring from 84+10 insertions to 105-11
    deletions is both a larger edit and substantially slower to route.
    """
    target = int(target)
    period = int(period)
    lower = max(period, (target // period) * period)
    candidates = {lower}
    if lower != target:
        candidates.add(lower + period)

    def score(nominal):
        period_count = max(1, nominal // period)
        residual = target - nominal
        magnitude = abs(residual)
        peak = int(math.ceil(magnitude / float(period_count)))
        average = magnitude / float(period_count)
        # In a complete tie, retain the larger lattice circumference.  This
        # avoids changing the established deletion-side tie behaviour.
        return magnitude, peak, average, -nominal

    return min(candidates, key=score)


def _unit(values):
    length = math.sqrt(sum(float(value) ** 2 for value in values))
    if length <= 1e-12:
        raise ValueError("DNAxiS produced a zero-length nucleotide axis.")
    return [float(value) / length for value in values]


def _geometry(origami):
    frames = {}
    for module in origami.get_modules():
        for strand_type, nucleotides in (
                ("scaffold", module.helix.scaf),
                ("staple", module.helix.stap)):
            for nucleotide in nucleotides:
                key = "%s:%d:%d" % (
                    strand_type, int(nucleotide.cdna_num),
                    int(nucleotide.cdna_ind))
                frames[key] = {
                    "pos": [float(value) / OXDNA_LENGTH_NM
                            for value in nucleotide.position],
                    "a1": _unit(nucleotide.backbonevec),
                    "a3": _unit(nucleotide.normvec)}
    keys = sorted(frames)
    return {
        "coordinate_units": "oxDNA",
        "source": "DNAxiS ideal ring geometry",
        "unrelaxed": True,
        "base_count": len(frames),
        "base_set_fingerprint": hashlib.sha256(
            "\n".join(keys).encode("utf-8")).hexdigest(),
        "frames": frames}


def _encode(value):
    raw = json.dumps(value, separators=(",", ":"),
                     sort_keys=True).encode("utf-8")
    return base64.b64encode(gzip.compress(raw, 9)).decode("ascii")


def _strand_report(json_obj, strand_type):
    helices = {int(row["num"]): row for row in json_obj["vstrands"]}
    nodes = set()
    starts = []
    ends = []
    transitions = set()
    for helix_id, row in helices.items():
        for base_index, connection in enumerate(row[strand_type]):
            if connection == [-1, -1, -1, -1]:
                continue
            nodes.add((helix_id, base_index))
            if int(connection[0]) < 0:
                starts.append((helix_id, base_index))
            if int(connection[2]) < 0:
                ends.append((helix_id, base_index))
            elif int(connection[2]) != helix_id:
                transitions.add(tuple(sorted((helix_id,
                                              int(connection[2])))))

    lengths = []
    visited = set()
    for initial in starts + list(nodes):
        if initial in visited:
            continue
        current = initial
        length = 0
        while current in nodes and current not in visited:
            visited.add(current)
            skip = int(helices[current[0]].get("skip", [])[current[1]])
            loop = int(helices[current[0]].get("loop", [])[current[1]])
            length += max(0, 1 + skip + loop)
            connection = helices[current[0]][strand_type][current[1]]
            current = (int(connection[2]), int(connection[3]))
        lengths.append(length)
    return {
        "helices": helices, "nodes": nodes, "starts": starts,
        "ends": ends, "lengths": lengths, "transitions": transitions}


def _cadnano_array_errors(json_obj, strand_type):
    """Apply the legacy caDNAno decoder's local topology invariants."""
    rows = {int(row["num"]): row for row in json_obj["vstrands"]}
    errors = []
    for helix_id, row in rows.items():
        strand = row[strand_type]
        offset = 1 if strand_type == "scaf" else -1
        endpoints = []
        for base_index, connection in enumerate(strand):
            if not isinstance(connection, list) or len(connection) != 4:
                errors.append("helix %d base %d has a malformed entry" %
                              (helix_id, base_index))
                continue
            five_helix, five_index, three_helix, three_index = [
                int(value) for value in connection]
            if five_helix == -1 and three_helix == -1:
                continue
            is_endpoint = (
                (five_helix == helix_id and three_helix != helix_id) or
                (five_helix != helix_id and three_helix == helix_id) or
                (helix_id % 2 == 0 and five_helix == helix_id and
                 five_index != base_index - offset) or
                (helix_id % 2 == 0 and three_helix == helix_id and
                 three_index != base_index + offset) or
                (helix_id % 2 == 1 and five_helix == helix_id and
                 five_index != base_index + offset) or
                (helix_id % 2 == 1 and three_helix == helix_id and
                 three_index != base_index - offset) or
                (five_helix == -1 and three_helix != -1) or
                (five_helix != -1 and three_helix == -1))
            if is_endpoint:
                endpoints.append(base_index)
            if five_helix != helix_id and three_helix != helix_id:
                endpoints.append(base_index)

            for direction, target_helix, target_index in (
                    ("5'", five_helix, five_index),
                    ("3'", three_helix, three_index)):
                if target_helix < 0:
                    continue
                target_row = rows.get(target_helix)
                if (target_row is None or target_index < 0 or
                        target_index >= len(target_row[strand_type])):
                    errors.append(
                        "helix %d base %d has an invalid %s target %d[%d]" %
                        (helix_id, base_index, direction,
                         target_helix, target_index))
                    continue
                target = target_row[strand_type][target_index]
                reciprocal = ((int(target[2]), int(target[3]))
                              if direction == "5'" else
                              (int(target[0]), int(target[1])))
                if reciprocal != (helix_id, base_index):
                    errors.append(
                        "helix %d base %d has a non-reciprocal %s link" %
                        (helix_id, base_index, direction))
        if len(endpoints) % 2:
            errors.append(
                "helix %d has an odd number of %s segment endpoints (%d)" %
                (helix_id, strand_type, len(endpoints)))

        occupied = [entry != [-1, -1, -1, -1] for entry in strand]

        def locally_connected(first, second):
            if first < 0 or second >= len(strand):
                return False
            if not occupied[first] or not occupied[second]:
                return False
            first_entry, second_entry = strand[first], strand[second]
            return ((helix_id, second) in
                    ((first_entry[0], first_entry[1]),
                     (first_entry[2], first_entry[3])) and
                    (helix_id, first) in
                    ((second_entry[0], second_entry[1]),
                     (second_entry[2], second_entry[3])))

        run_start = None
        for index, is_occupied in enumerate(occupied):
            if is_occupied and not locally_connected(index - 1, index):
                run_start = index
            if (is_occupied and not locally_connected(index, index + 1) and
                    run_start == index):
                errors.append(
                    "helix %d has a one-base %s segment at base %d" %
                    (helix_id, strand_type, index))
    return errors


def _normalize_cadnano_orientation(json_obj, origami):
    """Make DNAxiS ring arrays obey caDNAno's parity orientation rules.

    DNAxiS numbers every circular module in its native forward order.  Legacy
    caDNAno instead derives strand direction from virtual-helix parity.  The
    two conventions agree for only one of scaffold/staple on each row, so an
    unnormalised file can be interpreted as hundreds of one-base segments.
    """
    rows = {int(row["num"]): row for row in json_obj["vstrands"]}
    reverse = {}
    for helix_id, row in rows.items():
        for strand_type in ("scaf", "stap"):
            votes = 0
            for base_index, connection in enumerate(row[strand_type]):
                three_helix, three_index = int(connection[2]), int(connection[3])
                if three_helix != helix_id:
                    continue
                delta = three_index - base_index
                if delta == 1:
                    votes += 1
                elif delta == -1:
                    votes -= 1
            if not votes:
                raise RuntimeError(
                    "Cannot determine caDNAno direction for helix %d %s." %
                    (helix_id, strand_type))
            actual = 1 if votes > 0 else -1
            expected = (1 if helix_id % 2 == 0 else -1)
            if strand_type == "stap":
                expected *= -1
            reverse[(helix_id, strand_type)] = actual != expected

    for helix_id, row in rows.items():
        for strand_type in ("scaf", "stap"):
            source = row[strand_type]
            size = len(source)
            target = [[-1, -1, -1, -1] for unused in range(size)]
            for old_index, connection in enumerate(source):
                new_index = (size - 1 - old_index
                             if reverse[(helix_id, strand_type)]
                             else old_index)
                mapped = list(connection)
                for helix_position, index_position in ((0, 1), (2, 3)):
                    other_helix = int(mapped[helix_position])
                    other_index = int(mapped[index_position])
                    if (other_helix >= 0 and
                            reverse.get((other_helix, strand_type), False)):
                        other_size = len(rows[other_helix][strand_type])
                        mapped[index_position] = other_size - 1 - other_index
                target[new_index] = mapped
            row[strand_type] = target
        colors = []
        for index, color in row.get("stap_colors", []):
            if reverse[(helix_id, "stap")]:
                index = len(row["stap"]) - 1 - int(index)
            colors.append([int(index), color])
        row["stap_colors"] = colors

    for module in origami.get_modules():
        helix_id = int(module.cdna_num)
        for strand_type, nucleotides in (
                ("scaf", module.helix.scaf),
                ("stap", module.helix.stap)):
            if reverse[(helix_id, strand_type)]:
                size = len(rows[helix_id][strand_type])
                for nucleotide in nucleotides:
                    nucleotide.cdna_ind = size - 1 - int(nucleotide.cdna_ind)


def _assign_curved_lattice_coordinates(json_obj, ring_rows,
                                       lattice="square"):
    """Project layer/slice neighbours onto a valid caDNAno grid.

    The DNAxiS exporter places every ring in one visual column.  That makes
    radial neighbours in a multilayer surface appear many lattice cells
    apart, which SNUPI rejects even though their embedded 3-D coordinates are
    correct.  A layer-by-slice grid preserves both neighbour directions.
    """
    rows = {int(row["num"]): row for row in json_obj["vstrands"]}
    if len(rows) != len(ring_rows):
        raise RuntimeError(
            "Curved lattice projection cannot match %d helices to %d rings." %
            (len(rows), len(ring_rows)))
    lattice = str(lattice).lower()
    if lattice not in ("square", "honeycomb"):
        raise ValueError("Curved lattice must be square or honeycomb.")
    occupied = set()
    for helix_id, ring in enumerate(ring_rows):
        if helix_id not in rows:
            raise RuntimeError(
                "Curved lattice projection is missing helix %d." % helix_id)
        if lattice == "square":
            row_index = 2 + int(ring["slice"])
            col_index = 24 + int(ring["layer"])
        else:
            # Native brick-wall Honeycomb coordinates: every layer is the
            # mirror of the previous one.  Vertical/radial rungs therefore
            # occur at alternating slice indices and no node exceeds degree 3.
            row_index = 2 + int(ring["layer"])
            col_index = 24 + int(ring["slice"])
        coordinate = (row_index, col_index)
        if coordinate in occupied:
            raise RuntimeError(
                "Curved lattice projection produced duplicate coordinate %s." %
                (coordinate,))
        if (row_index + col_index) % 2 != helix_id % 2:
            raise RuntimeError(
                "Curved lattice projection parity mismatch for helix %d." %
                helix_id)
        occupied.add(coordinate)
        rows[helix_id]["row"] = row_index
        rows[helix_id]["col"] = col_index
    json_obj["lattice"] = lattice
    json_obj["num_bases"] = len(json_obj["vstrands"][0]["scaf"])


def _ensure_single_scaffold_nick(json_obj):
    """Open one native scaffold bond after topology-only DNAxiS routing."""
    report = _strand_report(json_obj, "scaf")
    if len(report["starts"]) == 1 and len(report["ends"]) == 1:
        return report["starts"][0]
    if report["starts"] or report["ends"]:
        raise RuntimeError("DNAxiS scaffold has an incomplete nick set.")
    rows = {int(row["num"]): row for row in json_obj["vstrands"]}

    def is_local_neighbour(connection, helix_id, index, side):
        position = (0, 1) if side == "five" else (2, 3)
        other_helix = int(connection[position[0]])
        other_index = int(connection[position[1]])
        return (other_helix == helix_id and
                abs(other_index - index) == 1)

    # A nick immediately beside a crossover would leave its terminal base as
    # a one-base caDNAno segment.  This occurred on the last ring of reinforced
    # (two-layer) curved designs.  Only cut an internal local bond whose two
    # resulting ends each retain a normal same-helix neighbour.
    for helix_id in sorted(rows, reverse=True):
        strand = rows[helix_id]["scaf"]
        for index in range(len(strand) - 1, -1, -1):
            entry = strand[index]
            target_helix, target_index = int(entry[2]), int(entry[3])
            if (target_helix != helix_id or target_index < 0 or
                    abs(target_index - index) != 1):
                continue
            target = strand[target_index]
            if (int(target[0]), int(target[1])) != (helix_id, index):
                continue
            if (not is_local_neighbour(entry, helix_id, index, "five") or
                    not is_local_neighbour(
                        target, helix_id, target_index, "three")):
                continue
            entry[2:4] = [-1, -1]
            target[0:2] = [-1, -1]
            return helix_id, target_index
    raise RuntimeError("Cannot place the single scaffold nick.")


def _curvature_indel_mapping(json_obj, ring_rows, lattice="square"):
    """Encode arbitrary ring circumferences on lattice-periodic helices.

    Each DNAxiS ring is represented by the nearest complete lattice period.
    The residual circumference is encoded by evenly distributed paired-base
    insertions or deletions.  Indels are selected only on native same-helix
    duplex runs and therefore never coincide with scaffold/staple crossovers.

    Returns the old-to-new index mapping used to re-key saved 3-D frames.
    """
    period = 21 if str(lattice).lower() == "honeycomb" else 32
    domain_size = 7 if str(lattice).lower() == "honeycomb" else 8
    maximum_per_domain = 3
    rows = {int(row["num"]): row for row in json_obj["vstrands"]}
    if set(rows) != set(range(len(ring_rows))):
        raise RuntimeError("Curved indel mapping cannot match ring helices.")

    nominal = {helix_id: int(ring["nominal_bp"])
               for helix_id, ring in enumerate(ring_rows)}
    max_size = max(nominal.values())
    mappings = {}
    originals = {}
    blank = [-1, -1, -1, -1]

    for helix_id, row in rows.items():
        originals[helix_id] = {
            strand_type: [list(entry) for entry in row[strand_type]]
            for strand_type in ("scaf", "stap")}
        occupied = set()
        for strand_type in ("scaf", "stap"):
            occupied.update(index for index, entry in enumerate(
                row[strand_type]) if entry != blank)
        old_indices = sorted(occupied)
        if len(old_indices) != nominal[helix_id] or any(
                right != left + 1
                for left, right in zip(old_indices, old_indices[1:])):
            raise RuntimeError(
                "DNAxiS ring %d does not occupy one contiguous %d-base "
                "lattice range." % (helix_id, nominal[helix_id]))
        mappings[helix_id] = {
            old_index: ordinal for ordinal, old_index in enumerate(old_indices)}

    for helix_id, row in rows.items():
        mapping = mappings[helix_id]
        new_arrays = {}
        for strand_type in ("scaf", "stap"):
            result = [list(blank) for unused in range(max_size)]
            for old_index, new_index in mapping.items():
                entry = list(originals[helix_id][strand_type][old_index])
                for helix_position, index_position in ((0, 1), (2, 3)):
                    other_helix = int(entry[helix_position])
                    other_index = int(entry[index_position])
                    if other_helix >= 0:
                        entry[index_position] = mappings[other_helix][
                            other_index]
                result[new_index] = entry
            new_arrays[strand_type] = result
        row["scaf"] = new_arrays["scaf"]
        row["stap"] = new_arrays["stap"]
        row["loop"] = [0] * max_size
        row["skip"] = [0] * max_size
        row["stap_colors"] = [
            [mapping[int(index)], color]
            for index, color in row.get("stap_colors", [])
            if int(index) in mapping]

    records = []
    maximum_insertion_per_crossover = 0
    maximum_deletion_per_crossover = 0
    maximum_insertion_per_domain = 0
    maximum_deletion_per_domain = 0
    for helix_id, row in rows.items():
        target = int(ring_rows[helix_id]["bp"])
        nominal_size = nominal[helix_id]
        residual = target - nominal_size
        crossover_periods = max(1, nominal_size // period)

        def unavailable(index):
            for strand_type in ("scaf", "stap"):
                entry = row[strand_type][index]
                if entry == blank:
                    return True
                if (int(entry[0]) < 0 or int(entry[2]) < 0 or
                        int(entry[0]) != helix_id or
                        int(entry[2]) != helix_id):
                    return True
            return False

        domain_count = max(1, int(math.ceil(
            nominal_size / float(domain_size))))
        indel_count = abs(residual)
        if indel_count > domain_count * maximum_per_domain:
            raise RuntimeError(
                "Ring %d requires %d indels across %d %d-bp domains; "
                "the hard limit is +/-3 per domain." %
                (helix_id, indel_count, domain_count, domain_size))
        available_by_domain = {}
        for domain in range(domain_count):
            candidates = [
                index for index in range(
                    max(1, domain * domain_size),
                    min(nominal_size - 1, (domain + 1) * domain_size))
                if not unavailable(index)]
            available_by_domain[domain] = candidates
        capacities = [
            (maximum_per_domain if residual > 0 and
             available_by_domain[domain] else
             min(maximum_per_domain,
                 len(available_by_domain[domain])))
            for domain in range(domain_count)]
        if sum(capacities) < indel_count:
            raise RuntimeError(
                "Ring %d cannot distribute %d curvature indels away from "
                "crossovers and nicks while retaining +/-3 per domain." %
                (helix_id, indel_count))

        quotient, remainder = divmod(indel_count, domain_count)
        quotas = [quotient] * domain_count
        if remainder:
            phase = helix_id % domain_count
            chosen_domains = []
            for order in range(remainder):
                desired = int(math.floor(
                    (order + 0.5) * domain_count / float(remainder)))
                candidate = (desired + phase) % domain_count
                while candidate in chosen_domains:
                    candidate = (candidate + 1) % domain_count
                chosen_domains.append(candidate)
                quotas[candidate] += 1
        for domain in range(domain_count):
            while quotas[domain] > capacities[domain]:
                targets = [
                    other for other in range(domain_count)
                    if quotas[other] < capacities[other]]
                if not targets:
                    raise RuntimeError(
                        "Ring %d cannot place its domain indel budget." %
                        helix_id)
                target_domain = min(targets, key=lambda other: (
                    min((other - domain) % domain_count,
                        (domain - other) % domain_count), other))
                quotas[domain] -= 1
                quotas[target_domain] += 1

        selected = []
        for domain, quota in enumerate(quotas):
            candidates = available_by_domain[domain]
            for order in range(quota):
                rank = min(len(candidates) - 1, int(math.floor(
                    (order + 0.5) * len(candidates) / float(quota))))
                chosen = candidates[rank]
                if residual < 0 and chosen in selected:
                    alternatives = [value for value in candidates
                                    if value not in selected]
                    if not alternatives:
                        raise RuntimeError(
                            "Ring %d cannot place distinct deletions in "
                            "domain %d." % (helix_id, domain))
                    chosen = alternatives[0]
                selected.append(chosen)
        if residual > 0:
            insertions = sorted(selected)
            for index in insertions:
                row["loop"][index] += 1
            deletions = []
            maximum_insertion_per_domain = max(
                maximum_insertion_per_domain, max(quotas or [0]))
        elif residual < 0:
            deletions = sorted(selected)
            for index in deletions:
                row["skip"][index] = -1
            insertions = []
            maximum_deletion_per_domain = max(
                maximum_deletion_per_domain, max(quotas or [0]))
        else:
            insertions, deletions = [], []
        insertion_average = len(insertions) / float(crossover_periods)
        deletion_average = len(deletions) / float(crossover_periods)
        insertion_peak = int(math.ceil(insertion_average))
        deletion_peak = int(math.ceil(deletion_average))
        maximum_insertion_per_crossover = max(
            maximum_insertion_per_crossover, insertion_peak)
        maximum_deletion_per_crossover = max(
            maximum_deletion_per_crossover, deletion_peak)
        records.append({
            "helix": helix_id, "target_bases": target,
            "nominal_bases": nominal_size,
            "crossover_periods": crossover_periods,
            "domain_size_bp": domain_size,
            "domain_indel_quotas": quotas,
            "maximum_indel_in_one_domain": max(quotas or [0]),
            "insertions": insertions, "deletions": deletions,
            "insertion_per_crossover_average": insertion_average,
            "deletion_per_crossover_average": deletion_average,
            "maximum_insertion_in_one_crossover_period": insertion_peak,
            "maximum_deletion_in_one_crossover_period": deletion_peak})

    json_obj["num_bases"] = max_size
    json_obj["curvature_indels"] = {
        "version": 2, "mode": "lattice-periodic-indels",
        "lattice": str(lattice).lower(), "period": period,
        "rings": records,
        "maximum_insertion_per_crossover":
            maximum_insertion_per_crossover,
        "maximum_deletion_per_crossover":
            maximum_deletion_per_crossover,
        "domain_size_bp": domain_size,
        "maximum_indel_per_domain_allowed": maximum_per_domain,
        "maximum_insertion_per_domain": maximum_insertion_per_domain,
        "maximum_deletion_per_domain": maximum_deletion_per_domain,
        "effective_crossover_spacing_minimum":
            period - maximum_deletion_per_crossover,
        "effective_crossover_spacing_maximum":
            period + maximum_insertion_per_crossover}
    return mappings


def _rekey_geometry(geometry, mappings, ring_rows):
    """Move DNAxiS physical frames to their lattice-indel design indices."""
    remapped = {}
    for key, frame in geometry.get("frames", {}).items():
        strand_type, helix_text, index_text = key.split(":")
        helix_id, old_index = int(helix_text), int(index_text)
        if old_index not in mappings[helix_id]:
            # This physical nucleotide is represented as an insertion subbase
            # and is reconstructed between its two retained target frames.
            continue
        new_index = mappings[helix_id][old_index]
        remapped["%s:%d:%d" % (
            strand_type, helix_id, new_index)] = frame
    geometry["frames"] = remapped
    keys = sorted(remapped)
    geometry["base_set_fingerprint"] = hashlib.sha256(
        "\n".join(keys).encode("utf-8")).hexdigest()
    geometry["curvature_encoding"] = "lattice-periodic-indels"
    geometry["base_count"] = 2 * sum(
        int(ring["bp"]) for ring in ring_rows)
    return geometry


def _validate_design(json_obj, required_edges):
    """Reject incomplete fast routes before cadnano opens them as designs."""
    scaffold = _strand_report(json_obj, "scaf")
    staples = _strand_report(json_obj, "stap")
    errors = (_cadnano_array_errors(json_obj, "scaf") +
              _cadnano_array_errors(json_obj, "stap"))

    if len(scaffold["lengths"]) != 1:
        errors.append("scaffold is split into %d components" %
                      len(scaffold["lengths"]))
    if len(scaffold["starts"]) != 1 or len(scaffold["ends"]) != 1:
        errors.append("scaffold must contain exactly one nick")
    elif scaffold["starts"][0][0] != scaffold["ends"][0][0]:
        errors.append("scaffold nick ends are on different helices")
    else:
        helix_id = scaffold["starts"][0][0]
        occupied = sorted(index for hid, index in scaffold["nodes"]
                          if hid == helix_id)
        start_index = scaffold["starts"][0][1]
        end_index = scaffold["ends"][0][1]
        step = 1 if helix_id % 2 == 0 else -1
        adjacent = (end_index + step == start_index)
        if step > 0:
            adjacent = adjacent or (
                end_index == occupied[-1] and start_index == occupied[0])
        else:
            adjacent = adjacent or (
                end_index == occupied[0] and start_index == occupied[-1])
        if not adjacent:
            errors.append("scaffold path does not close at its single nick")

    missing_staple = sorted(
        edge for edge in required_edges
        if edge not in staples["transitions"])
    if missing_staple:
        errors.append("missing staple crossovers on neighbouring helices: %s" %
                      ", ".join("%d-%d" % edge
                                for edge in missing_staple))
    bad_lengths = sorted(length for length in staples["lengths"]
                         if length < 21 or length > 58)
    if bad_lengths:
        errors.append("staple lengths outside 21-58 nt: %s" % bad_lengths)
    if errors:
        raise RuntimeError("Curved routing validation failed: " +
                           "; ".join(errors))


def _rebalance_staple_lengths(origami, minimum=21, maximum=58):
    """Apply Curved Design's lattice-neutral Autobreak rules.

    Nick placement shares the ordinary Autobreak length, crossover-clearance,
    and continuous-run rules, but deliberately has no square/honeycomb phase
    preference because a curved ring does not have a fixed lattice period.
    """

    modules = list(origami.get_modules())
    local_index = {}
    module_size = {}
    for module in modules:
        module_size[module] = len(module.helix.stap)
        for index, nucleotide in enumerate(module.helix.stap):
            local_index[id(nucleotide)] = index

    def external_positions(module, nucleotides):
        positions = set()
        for index, nucleotide in enumerate(nucleotides):
            for neighbour in (nucleotide.toFive, nucleotide.toThree):
                if (not isinstance(neighbour, int) and
                        neighbour.get_top_module() is not module):
                    positions.add(index)
        return positions

    staple_xovers = {
        module: external_positions(module, module.helix.stap)
        for module in modules}
    scaffold_xovers = {
        module: external_positions(module, module.helix.scaf)
        for module in modules}

    def ring_distance(first, second, size):
        distance = abs(first - second)
        return min(distance, size - distance)

    def boundary_penalty(left, right, minimum_xover_distance):
        module = left.get_top_module()
        if right.get_top_module() is not module:
            return None
        left_index = local_index[id(left)]
        right_index = local_index[id(right)]
        size = module_size[module]
        if ring_distance(left_index, right_index, size) != 1:
            return None

        staple_positions = staple_xovers[module]
        staple_clearance_reward = 0.0
        if staple_positions:
            nearest_staple = min(
                ring_distance(index, xover, size)
                for index in (left_index, right_index)
                for xover in staple_positions)
            # First solve with the square-Autobreak clearance of >=7.  If a
            # component has no such solution, the fallback has no advertised
            # lower threshold: it simply maximizes clearance and the affected
            # final staples are marked red for manual correction.
            if nearest_staple < minimum_xover_distance:
                return None
            staple_clearance_reward = (
                -1000.0 * min(float(nearest_staple), 7.0))

        scaffold_positions = scaffold_xovers[module]
        if not scaffold_positions:
            return staple_clearance_reward
        nearest_scaffold = min(
            ring_distance(index, xover, size)
            for index in (left_index, right_index)
            for xover in scaffold_positions)
        if nearest_scaffold == 0:
            return None
        # As in Autobreak, >3 nt is preferred rather than mandatory.  Once
        # both sites are safe, greater clearance is a late tie-breaker.
        return (staple_clearance_reward +
                (50000.0 if nearest_scaffold <= 3 else 0.0) -
                min(float(nearest_scaffold), 30.0))

    def solve_linear(path, cyclic, minimum_xover_distance=7):
        size = len(path)
        if size < minimum:
            return None

        def creates_single_base_segment(ordered, boundary):
            if boundary == 0 or boundary == size:
                if not cyclic:
                    return False
                left_previous = ordered[-2]
                left = ordered[-1]
                right = ordered[0]
                right_next = ordered[1]
            else:
                left_previous = ordered[boundary - 2]
                left = ordered[boundary - 1]
                right = ordered[boundary]
                right_next = (ordered[boundary + 1]
                              if boundary + 1 < size else None)
            left_is_single = (left_previous.get_top_module() is not
                              left.get_top_module())
            right_is_single = (
                right_next is not None and
                right.get_top_module() is not right_next.get_top_module())
            return left_is_single or right_is_single

        def has_continuous_run(segment, required=14):
            longest = current = 0
            previous = None
            for nucleotide in segment:
                module = nucleotide.get_top_module()
                index = local_index[id(nucleotide)]
                if previous is not None:
                    previous_module, previous_index = previous
                    adjacent = (
                        module is previous_module and
                        ring_distance(previous_index, index,
                                      module_size[module]) == 1)
                    current = current + 1 if adjacent else 1
                else:
                    current = 1
                longest = max(longest, current)
                previous = (module, index)
            return longest >= required

        starts = range(size) if cyclic else (0,)
        for start in starts:
            ordered = path[start:] + path[:start] if cyclic else path
            if cyclic and creates_single_base_segment(ordered, 0):
                continue
            start_penalty = (boundary_penalty(
                                ordered[-1], ordered[0],
                                minimum_xover_distance)
                             if cyclic else 0.0)
            if start_penalty is None:
                continue
            best = {0: (start_penalty, None)}
            for position in range(size):
                if position not in best:
                    continue
                for length in range(minimum, maximum + 1):
                    target = position + length
                    if target > size:
                        break
                    if target < size:
                        # caDNAno's legacy strand model cannot round-trip a
                        # one-base local segment.  Do not nick immediately
                        # after entering a module by crossover, or immediately
                        # before leaving it by crossover.
                        if creates_single_base_segment(ordered, target):
                            continue
                    penalty = (boundary_penalty(
                        ordered[target - 1], ordered[target],
                        minimum_xover_distance)
                        if target < size else 0.0)
                    if penalty is None:
                        continue
                    segment = ordered[position:target]
                    continuous_penalty = (
                        0.0 if has_continuous_run(segment) else 1000000.0)
                    soft_minimum_penalty = (
                        0.0 if length >= 28 else 100000.0)
                    preferred_range_penalty = (
                        1000.0 * max(0, 30 - length, length - 50))
                    score = (best[position][0] +
                             continuous_penalty + soft_minimum_penalty +
                             preferred_range_penalty +
                             (length - 40) ** 2 + penalty)
                    if target not in best or score < best[target][0]:
                        best[target] = (score, position)
            if size not in best:
                continue
            boundaries = []
            position = size
            while position:
                previous = best[position][1]
                if previous is None:
                    break
                if previous:
                    boundaries.append(previous)
                position = previous
            return ordered, sorted(boundaries)
        return None

    all_bases = [nucleotide
                 for module in origami.get_modules()
                 for nucleotide in module.helix.stap]
    successor = {}
    predecessor = {}
    for nucleotide in all_bases:
        next_nucleotide = (nucleotide.toThree
                           if nucleotide.toThree != -1
                           else nucleotide.__strand3__)
        successor[id(nucleotide)] = next_nucleotide
        if not isinstance(next_nucleotide, int):
            predecessor[id(next_nucleotide)] = nucleotide
    visited = set()
    components = []
    for base in all_bases:
        if id(base) in visited:
            continue
        start = base
        seen_back = set()
        while (id(start) in predecessor and
               id(predecessor[id(start)]) not in seen_back):
            seen_back.add(id(start))
            start = predecessor[id(start)]
            if start is base:
                break
        component = []
        current = start
        while not isinstance(current, int) and id(current) not in visited:
            component.append(current)
            visited.add(id(current))
            current = successor[id(current)]
        cyclic = bool(component and current is component[0])
        components.append((component, cyclic))

    for component, cyclic in components:
        if len(component) < minimum:
            raise RuntimeError(
                "A %d nt staple topology component cannot be repaired "
                "without changing the existing crossover layout." %
                len(component))
        solved = solve_linear(component, cyclic, 7)
        if solved is None:
            solved = solve_linear(component, cyclic, 0)
        if solved is None:
            raise RuntimeError(
                "Cannot partition an absolute staple path of %d nt into "
                "serializable %d-%d nt strands without changing the existing "
                "crossover layout." %
                (len(component), minimum, maximum))
        ordered, boundaries = solved
        for nucleotide in ordered:
            next_nucleotide = successor[id(nucleotide)]
            if not isinstance(next_nucleotide, int):
                dnaconnector.nuclconnect(nucleotide, next_nucleotide)
        if cyclic:
            dnaconnector.nuclbreak(ordered[-1])
        for boundary in boundaries:
            dnaconnector.nuclbreak(ordered[boundary - 1])


def _mark_close_xover_staples(json_obj, minimum_distance=7):
    """Color staples red when either end is too close to a staple xover."""
    rows = {int(row["num"]): row for row in json_obj["vstrands"]}
    xovers = {}
    starts = []
    nodes = set()
    for helix_id, row in rows.items():
        positions = set()
        for index, connection in enumerate(row["stap"]):
            if connection == [-1, -1, -1, -1]:
                continue
            nodes.add((helix_id, index))
            if int(connection[0]) < 0:
                starts.append((helix_id, index))
            if ((int(connection[0]) >= 0 and
                 int(connection[0]) != helix_id) or
                    (int(connection[2]) >= 0 and
                     int(connection[2]) != helix_id)):
                positions.add(index)
        xovers[helix_id] = positions

    warning_starts = set()
    visited = set()
    for start in starts:
        current = start
        path = []
        while current in nodes and current not in visited:
            visited.add(current)
            path.append(current)
            connection = rows[current[0]]["stap"][current[1]]
            current = (int(connection[2]), int(connection[3]))
        if not path:
            continue
        unsafe = False
        for helix_id, index in (path[0], path[-1]):
            positions = xovers.get(helix_id, ())
            if positions and min(abs(index - value)
                                 for value in positions) < minimum_distance:
                unsafe = True
        if unsafe:
            warning_starts.add(start)

    red = 0xCC0000
    for helix_id, index in warning_starts:
        colors = rows[helix_id].setdefault("stap_colors", [])
        colors[:] = [record for record in colors
                     if int(record[0]) != index]
        colors.append([index, red])
        colors.sort(key=lambda record: int(record[0]))
    return len(warning_starts)


def run(spec_path, result_path):
    with open(spec_path, "r", encoding="utf-8") as source:
        spec = json.load(source)
    output_dir = os.path.abspath(spec["output_dir"])
    os.makedirs(output_dir, exist_ok=True)
    ring_rows = [dict(row) for row in spec["rings"]]
    lattice = str(spec.get("lattice", "square")).lower()
    period = 21 if lattice == "honeycomb" else 32
    for row in ring_rows:
        target = int(row["bp"])
        if "nominal_bp" not in row:
            row["nominal_bp"] = _balanced_nominal_bp(target, period)
        row["nominal_bp"] = int(row["nominal_bp"])
    rings = [[int(row["nominal_bp"]), float(row["height_nm"]),
              bool(row["direction"]),
              float(row.get("geometry_radius_nm", row["radius_nm"]))]
             for row in ring_rows]
    count = len(rings)
    if count < 2:
        raise ValueError("A curved design requires at least two DNA rings.")

    connection_map = {}
    pathway = []
    # The scaffold pathway follows the serpentine input order.
    for index in range(count - 1):
        first, second = [index, 0], [index + 1, 0]
        pathway.append({"source": first, "target": second})
    # Staple/scaffold crossover candidates include every vertical and radial
    # neighbour, not just the single scaffold pathway.  This is important for
    # reinforced two- and three-layer targets.
    lookup = {(int(row["layer"]), int(row["slice"])): index
              for index, row in enumerate(ring_rows)}
    edge_set = set()
    for (layer, slice_index), index in lookup.items():
        neighbours = [(layer, slice_index + 1)]
        if lattice == "square" or (layer + slice_index) % 2 == 1:
            neighbours.append((layer + 1, slice_index))
        for neighbour in neighbours:
            other = lookup.get(neighbour)
            if other is not None:
                edge_set.add(tuple(sorted((index, other))))
    for first_index, second_index in sorted(edge_set):
        first, second = [first_index, 0], [second_index, 0]
        connection_map.setdefault((first_index, 0), []).append(second)
        connection_map.setdefault((second_index, 0), []).append(first)
    connections = {"connections": [
        {"source": list(source), "targets": targets}
        for source, targets in sorted(connection_map.items())]}
    _write_json(os.path.join(
        output_dir, "connections_3.0_0.0.json"), connections)
    _write_json(os.path.join(output_dir, "pathway.json"),
                {"pathway": pathway})

    log.new("curved", output_dir, console=True, debug=False,
            developermode=False, log=True)
    config.SCAF_NICKING = "asymmetric_single"
    config.VALIDXOVERTHRESHBP = 3.0
    config.VALIDXOVERSPACING_SAME = 14
    config.VALIDXOVERSPACING_ADJ = 7
    random.seed(1)
    # Curved Design leaves sequence assignment for later.  DNAxiS still uses
    # sequence length to decide whether a route must be split, so install one
    # length-only virtual scaffold spanning the complete requested topology,
    # up to the shared 25,000-base Curved Design limit.  This removes the
    # bundled-scaffold capacity limit without allocating or exporting a
    # synthetic sequence.
    base_total = sum(int(row["bp"]) for row in ring_rows)
    _install_curved_scaffold_capacity(base_total)
    shape = CustomInput(rings, spec["output_name"])
    origami = driver(
        spec["output_name"], output_dir, shape, 4, (3.0, 0.0), 1,
        # Build DNAxiS' global geometry-legal initial crossover set, then use
        # the bounded empty-edge repair added in stxover.  Full simulated
        # annealing remains disabled because the initial-set bug and repeated
        # global rescoring made it unbounded for interactive use.
        skip_sequence=True, optimize_xovers=False, use_old_routing=True,
        # Native cadnano AutoCS creates all final staples in the parent
        # process.  Skipping DNAxiS staple routing and nick optimization is
        # the principal performance improvement of the lattice-indel path.
        skip_routing_stap=True, skip_nicks=True,
        force_reseeding=False, enable_validate=False,
        enable_stats=False, twist_normalized=False)
    # DNAxiS supplies the target-space geometry and neighbour topology.  The
    # main cadnano process replaces these provisional staples with its native
    # lattice AutoCS/Autobreak result after this child process returns.
    json_obj = origami.print_caDNAno(spec["output_name"])
    _assign_curved_lattice_coordinates(json_obj, ring_rows, lattice)
    _normalize_cadnano_orientation(json_obj, origami)
    _ensure_single_scaffold_nick(json_obj)
    mappings = _curvature_indel_mapping(json_obj, ring_rows, lattice)
    indel_summary = {
        key: json_obj["curvature_indels"][key]
        for key in (
            "maximum_insertion_per_crossover",
            "maximum_deletion_per_crossover",
            "effective_crossover_spacing_minimum",
            "effective_crossover_spacing_maximum")}
    # At this stage only scaffold topology is final; native staple AutoCS and
    # Autobreak are intentionally applied in the parent cadnano runtime.
    scaffold_errors = _cadnano_array_errors(json_obj, "scaf")
    if len(_strand_report(json_obj, "scaf")["lengths"]) != 1:
        scaffold_errors.append("scaffold is not a single component")
    if scaffold_errors:
        raise RuntimeError("Curved scaffold validation failed: " +
                           "; ".join(scaffold_errors))
    warning_staples = 0
    geometry = _rekey_geometry(_geometry(origami), mappings, ring_rows)
    metadata = dict(spec["metadata"])
    metadata.update({
        "metadata_version": 1,
        "lattice": lattice,
        "curvature_encoding": "lattice-periodic-indels",
        "native_staple_rules_pending": True,
        "red_staple_warning_count": warning_staples,
        "red_staple_warning": (
            "Red staples have at least one end fewer than 7 indices from "
            "a staple crossover and require manual review."),
        "geometry_encoding": "gzip+base64+json",
        "geometry_data": _encode(geometry)})
    metadata.update(indel_summary)
    json_obj["curved_metadata"] = metadata
    json_path = os.path.join(
        os.path.abspath(spec["project_root"]),
        spec["output_name"] + ".json")
    _write_json(json_path, json_obj)
    _write_json(result_path, {
        "json_path": json_path,
        "base_count": geometry["base_count"],
        "ring_count": count,
        "indel_summary": indel_summary})


if __name__ == "__main__":
    try:
        run(sys.argv[1], sys.argv[2])
    except Exception as error:
        import traceback
        traceback.print_exc()
        sys.exit(1)
