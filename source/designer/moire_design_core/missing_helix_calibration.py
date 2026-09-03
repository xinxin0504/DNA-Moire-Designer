"""Frozen Square missing-helix calibration used by design prediction.

The nodes originate from the finalized 2026-08-02 SNUPI matrix; the
S8-R4x4C reference node uses its converged nine-point 2026-08-13 extension.
Each profile is
``(x0, T0, k_minus, k_plus)`` and each descriptor contains missing fraction,
eccentricity, polar-moment ratio, boundary fraction, mean neighbour degree
and mean wall depth.  This copy keeps DNA Moiré Designer independent from a
particular cadnano installation.
"""

from __future__ import annotations

import math


# Extended S8-R4x4C SNUPI calibration (2026-08-13).  Only converged points
# participate in the calibration.  +8 and +10 indels/helix are retained as
# failed-boundary metadata and are never treated as observations.
S8_R4X4C_CALIBRATION_VERSION = "2026-08-14-S8-R4x4C-extended"
S8_R4X4C_TWIST_POINTS = (
    (-10.0, -0.16367829784857435),
    (-8.0, -0.11353790057673968),
    (-6.0, -0.05083101129187896),
    (-4.0, -0.017839876583065965),
    (-2.0, 0.03304955169447221),
    (0.0, 0.10302360948900995),
    (2.0, 0.15808502575673034),
    (4.0, 0.22000278208910148),
    (6.0, 0.24969345505563809),
)
S8_R4X4C_VALIDATED_INDEL_RANGE = (-10.0, 6.0)
S8_R4X4C_FAILED_INDEL_POINTS = (8.0, 10.0)
S8_R4X4C_EXTENDED_PROFILE = (
    0.08165743050280329,
    0.10302360948900995,
    1.035727869051505,
    0.9977630088056791,
)
S8_R4X4C_BRANCH_RMSE = {
    "deletion": 0.009110560088416491,
    "insertion": 0.008136450352592136,
}


NODES = (
    ('S4-R1x2C',4,'C',(.125,.235702260,.968537415,1.,2.428571429,1.),(.561007367,.477297178,.706329094,.751327784)),
    ('S4-R1x2E',4,'E',(.125,.707106781,.859693878,.857142857,2.714285714,1.142857143),(.576740809,.337425277,.864077486,.908157786)),
    ('S4-R1x2K',4,'K',(.125,.849836586,.805272109,.785714286,2.857142857,1.214285714),(.585371794,.462729718,.753310395,.751668885)),
    ('S4-R2x2C',4,'C',(.25,0.,.940476190,1.,2.,1.),(.564887684,.449354789,.634653214,1.094956113)),
    ('S4-R2x2E',4,'E',(.25,.471404521,.813492063,1.,2.333333333,1.),(.584031414,.664087596,1.165721935,.861220959)),
    ('S4-R2x2K',4,'K',(.25,.666666667,.686507937,.916666667,2.666666667,1.083333333),(.606430859,.417798383,1.035124266,.846163656)),
    ('S4-R2x2K_2',4,'K',(.25,0.,.559523810,.666666667,2.666666667,1.333333333),(.633424967,.504541394,.947602332,1.096323450)),
    ('S6-R2x2C',6,'C',(.111111111,0.,.988344988,.875,3.,1.125),(.158019720,.193332061,1.515268594,1.511554344)),
    ('S6-R2x2E',6,'E',(.111111111,.565685425,.904428904,.6875,3.125,1.3125),(.179865712,.225743994,1.861063933,1.847519078)),
    ('S6-R2x2E_2',6,'E',(.222222222,0.,.827505828,.857142857,2.857142857,1.142857143),(.201750959,.422124446,2.289586380,2.046250970)),
    ('S6-R2x2K',6,'K',(.111111111,.8,.820512821,.59375,3.25,1.5),(.203840472,.135060512,1.489215175,1.188170556)),
    ('S6-R2x2K_4',6,'K',(.444444444,0.,.356643357,.8,2.8,1.2),(.408983263,.424903801,1.373639477,.971571534)),
    ('S6-R3x3C',6,'C',(.25,.2,.910839161,.962962963,2.666666667,1.037037037),(.178126803,.188588430,1.177075608,1.943512564)),
    ('S6-R3x3K',6,'K',(.25,.6,.687062937,.703703704,3.111111111,1.296296296),(.247544128,.418436067,1.539729047,1.137400163)),
    ('S6-R4x4C',6,'C',(.444444444,0.,.804195804,1.,2.,1.),(.208786088,.226542264,.874960048,6.587030388)),
    ('S6-R4x4K',6,'K',(.444444444,.4,.535664336,.95,2.8,1.05),(.308830953,.490870291,1.499822692,1.112018931)),
    ('S8-R2x2C',8,'C',(.0625,0.,.996323529,.6,3.333333333,1.4),(.077606413,.041856563,1.168361281,1.158166519)),
    ('S8-R2x2E',8,'E',(.0625,.606091527,.939852941,.5,3.4,1.666666667),(.081541260,.096648093,1.976614340,1.807923208)),
    ('S8-R2x2K',8,'K',(.0625,.857142857,.883382353,.45,3.466666667,1.866666667),(.085720005,.084860134,1.356187341,1.342270350)),
    ('S8-R2x2K_4',8,'K',(.25,0.,.561764706,.5,3.333333333,1.833333333),(.116246984,.117873057,2.054917119,2.170422291)),
    ('S8-R3x3C',8,'C',(.140625,.142857143,.972997995,.727272727,3.2,1.272727273),(.079203993,.107998106,1.402842367,1.304674655)),
    ('S8-R3x3K',8,'K',(.140625,.714285714,.788185160,.490909091,3.418181818,1.709090909),(.093409497,.138568454,1.494108746,1.687341235)),
    ('S8-R4x4C',8,'C',(.25,0.,.938235294,.916666667,3.,1.083333333),S8_R4X4C_EXTENDED_PROFILE),
    ('S8-R4x4E',8,'E',(.25,.404061018,.812745098,.708333333,3.166666667,1.291666667),(.091340232,.180880459,2.308877319,2.780836880)),
    ('S8-R4x4K',8,'K',(.25,.571428571,.687254902,.5625,3.333333333,1.5),(.102650219,.169049651,1.766018257,1.870029124)),
    ('S8-R5x5C',8,'C',(.390625,.142857143,.818179676,.974358974,2.666666667,1.025641026),(.090890803,.088834685,1.011932824,1.372376015)),
    ('S8-R5x5K',8,'K',(.390625,.428571429,.576852376,.692307692,3.179487179,1.307692308),(.114459685,.268136369,2.591560728,2.311867627)),
    ('S8-R6x6C',8,'C',(.5625,0.,.684558824,1.,2.,1.),(.102915292,.137579998,1.220530038,14.506472760)),
    ('S8-R6x6K',8,'K',(.5625,.285714286,.442542017,.964285714,2.857142857,1.035714286),(.132333888,.398841487,3.379947717,2.711179989)),
)

VALIDATION_RMSE = .09194196915609779
REFERENCE_NODE = 'S8-R4x4C'


def default_seed_cells(size=8):
    """Return the calibrated S8-R4x4C 48-helix frame."""
    return tuple((row, col) for row in range(size) for col in range(size)
                 if not (2 <= row <= 5 and 2 <= col <= 5))


def _polar_moment(cells):
    points = [(float(row) * 2.0, float(col) * 2.0)
              for row, col in cells]
    if not points:
        return 0.0
    cx = sum(point[0] for point in points) / len(points)
    cy = sum(point[1] for point in points) / len(points)
    intrinsic = math.pi / 2.0
    area = math.pi
    return sum(intrinsic + area * ((x-cx)**2 + (y-cy)**2)
               for x, y in points)


def square_section_descriptor(cells, size=8):
    """Describe a user-selected Square-grid Seed cross-section."""
    size = int(size)
    if size not in (4, 6, 8):
        raise ValueError("Square Seed校准仅支持4×4、6×6或8×8外框。")
    present = {(int(row), int(col)) for row, col in cells
               if 0 <= int(row) < size and 0 <= int(col) < size}
    if len(present) < 4:
        raise ValueError("Seed截面至少需要4根helix。")
    expected = {(row, col) for row in range(size) for col in range(size)}
    missing = sorted(expected - present)
    center = (size - 1) / 2.0
    if missing:
        mr = sum(point[0] for point in missing) / len(missing)
        mc = sum(point[1] for point in missing) / len(missing)
        eccentricity = math.hypot(mr-center, mc-center) / max(
            math.sqrt(2.0) * center, 1e-12)
    else:
        eccentricity = 0.0
    degrees, depths = [], []
    boundary = 0
    for row, col in present:
        degree = sum((row+dr, col+dc) in present
                     for dr, dc in ((-1, 0), (1, 0), (0, -1), (0, 1)))
        degrees.append(degree)
        boundary += int(degree < 4)
        depths.append(min(abs(row-r) + abs(col-c)
                          for r in range(-1, size+1)
                          for c in range(-1, size+1)
                          if (r, c) not in present))
    descriptor = (
        len(missing) / float(size * size), eccentricity,
        _polar_moment(present) / max(_polar_moment(expected), 1e-12),
        boundary / float(len(present)),
        sum(degrees) / float(len(degrees)),
        sum(depths) / float(len(depths)),
    )
    return descriptor


def calibration_profile(cells, size=8):
    """Interpolate the frozen topology-aware calibration profile."""
    descriptor = square_section_descriptor(cells, size)
    candidates = [node for node in NODES if node[1] == int(size)]
    scales = [max(max(node[3][index] for node in candidates) -
                  min(node[3][index] for node in candidates), 1e-9)
              for index in range(len(descriptor))]
    ranked = []
    for node in candidates:
        distance = math.sqrt(sum(
            ((descriptor[index]-node[3][index])/scales[index])**2
            for index in range(len(descriptor))))
        ranked.append((distance, node))
    ranked.sort(key=lambda item: item[0])
    exact = ranked[0][0] < 1e-6
    excluded = {'S6-R3x3C', 'S6-R4x4C', 'S8-R5x5C', 'S8-R6x6C'}
    pool = ranked if exact else [item for item in ranked
                                 if item[1][0] not in excluded]
    chosen = [ranked[0]] if exact else pool[:min(4, len(pool))]
    if exact:
        weights = [1.0]
    else:
        weights = [(distance + .15) ** -2 for distance, unused in chosen]
        total = sum(weights)
        weights = [value / total for value in weights]
    profile = tuple(sum(weights[index] * chosen[index][1][4][field]
                        for index in range(len(chosen)))
                    for field in range(4))
    return {
        "profile": profile,
        "descriptor": descriptor,
        "exact": exact,
        "nearest_distance": chosen[0][0],
        "extrapolated": chosen[0][0] > 1.5,
        "neighbors": tuple((item[1][0], weights[index], item[0])
                           for index, item in enumerate(chosen)),
        "validation_rmse_deg_per_base": VALIDATION_RMSE,
    }


def reference_profile():
    return next(node[4] for node in NODES if node[0] == REFERENCE_NODE)
