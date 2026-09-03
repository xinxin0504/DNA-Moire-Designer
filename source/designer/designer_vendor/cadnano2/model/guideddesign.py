"""Pure helpers for cadnano's beginner Guided Design workflow.

The guide treats the uploaded picture as a planar silhouette.  Bases run
along the picture's horizontal axis and sampled image rows select helices.
Additional layers affect only the cross-section lattice; they never mirror
the uploaded silhouette.
"""

from math import ceil


BASE_RISE_NM = 0.34
HELIX_SPACING_NM = 2.8

def guided_base_offset(lattice):
    """Use a phase-aligned margin that keeps the first xover past base 50."""
    return 42 if lattice == 'honeycomb' else 58


def base_count_for_width(width_nm):
    return max(2, int(round(float(width_nm) / BASE_RISE_NM)))


def profile_count_for_height(height_nm):
    return max(1, int(round(float(height_nm) / HELIX_SPACING_NM)))


def max_steps_for_bases(base_count, lattice):
    step = 21 if lattice == 'honeycomb' else 32
    # One additional period is reserved for legal scaffold-only perimeter
    # closure; it is not part of the image-derived duplex target.
    return max(4, int(ceil(float(
        base_count + guided_base_offset(lattice) + step) / step)))


def cross_section_coords(profile_count, layers, lattice,
                         honeycomb_direction='z'):
    """Return lattice coordinates and their logical profile/layer source.

    In the default honeycomb Z direction consecutive rows are A, mirror-A,
    A, mirror-A.  This follows from honeycomb row parity; the target profile
    itself is reused unchanged for every layer.
    """
    raw = []
    for layer in range(max(1, int(layers))):
        for profile in range(max(1, int(profile_count))):
            if lattice != 'honeycomb' or honeycomb_direction == 'z':
                row, col = layer, profile
            elif honeycomb_direction == 'up-right':
                row, col = profile, profile + layer
            else:  # down-right
                row, col = profile + layer, profile
            raw.append((row, col, profile, layer))
    min_row = min(record[0] for record in raw)
    min_col = min(record[1] for record in raw)
    return [(row - min_row + 1, col - min_col + 1, profile, layer)
            for row, col, profile, layer in raw]


def boolean_runs(values, target_bases):
    """Scale a boolean scanline into inclusive cadnano base intervals."""
    values = list(values)
    if not values:
        return []
    sampled = []
    for base in range(max(2, int(target_bases))):
        source = min(len(values) - 1,
                     int(float(base) * len(values) / target_bases))
        sampled.append(bool(values[source]))
    runs = []
    start = None
    for index, occupied in enumerate(sampled + [False]):
        if occupied and start is None:
            start = index
        elif not occupied and start is not None:
            if index - start >= 2:
                runs.append((start, index - 1))
            start = None
    return runs


def estimate_scaffold_length(profile_runs, layer_count):
    one_layer = sum(high - low + 1 for runs in profile_runs
                    for low, high in runs)
    return one_layer * max(1, int(layer_count))
