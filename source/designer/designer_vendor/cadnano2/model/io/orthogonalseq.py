"""Offline orthogonal DNA sequence design and reporting.

Sequence length/count, global GC, homopolymer, and pairwise same/reverse-
complement substring limits are always active. Additional biochemical screens
are optional and disabled by default.
"""

import math
import random
from collections import Counter
from xml.sax.saxutils import escape
from zipfile import ZIP_DEFLATED, ZipFile


DNA_BASES = "ACGT"
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

MELTING_NA_MM = 50.0
MELTING_MG_MM = 10.0
MELTING_STRAND_NM = 100.0
MELTING_TARGET_NM = 100.0

# Allawi & SantaLucia (1997) DNA/DNA nearest-neighbor values.  Values are
# (delta H kcal/mol, delta S cal/mol/K), matching the widely used DNA_NN3
# model.  Salt is corrected in entropy using the SantaLucia 1998 method.
_DNA_NN = {
    "AA/TT": (-7.9, -22.2), "AT/TA": (-7.2, -20.4),
    "TA/AT": (-7.2, -21.3), "CA/GT": (-8.5, -22.7),
    "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2),
    "GC/CG": (-9.8, -24.4), "GG/CC": (-8.0, -19.9),
}

_DNA_RESIDUE_MASS = {
    "A": 313.21, "C": 289.18, "G": 329.21, "T": 304.20,
}

_DNA_DINUCLEOTIDE_EXTINCTION = {
    "AA": 27400, "AC": 21200, "AG": 25000, "AT": 22800,
    "CA": 21200, "CC": 14600, "CG": 18000, "CT": 15200,
    "GA": 25200, "GC": 17600, "GG": 21600, "GT": 20000,
    "TA": 23400, "TC": 16200, "TG": 19000, "TT": 16800,
}
_DNA_BASE_EXTINCTION = {
    "A": 15400, "C": 7400, "G": 11500, "T": 8700,
}

ISSUE_LABELS_CN = {
    "global_gc": "全局GC超出范围",
    "local_gc": "局部GC超出范围",
    "homopolymer": "连续相同碱基过长",
    "entropy": "序列复杂度过低",
    "self_complement": "自身互补过长",
    "hairpin": "发卡茎过长",
    "forbidden_motif": "包含禁用片段",
    "same_substring": "同向相同片段过长",
    "cross_complement": "链间互补片段过长",
    "hamming": "汉明距离不足",
}

SETTING_LABELS_CN = {
    "length": "新序列长度",
    "count": "请求生成数量",
    "gc_min": "全局GC下限",
    "gc_max": "全局GC上限",
    "max_homopolymer": "最大连续相同碱基",
    "max_same_substring": "最大同向相同片段",
    "max_cross_complement": "最大链间互补片段",
    "scaffold_max_same_substring": "与骨架最大同向相同片段",
    "scaffold_max_cross_complement": "与骨架最大链间互补片段",
    "use_local_gc": "启用局部GC规则",
    "local_gc_window": "局部GC窗口长度",
    "local_gc_min": "局部GC下限",
    "local_gc_max": "局部GC上限",
    "use_entropy": "启用序列熵规则",
    "min_entropy": "最低序列熵",
    "use_self_complement": "启用自身互补规则",
    "max_self_complement": "最大自身互补长度",
    "use_hairpin": "启用发卡规则",
    "max_hairpin_stem": "最大发卡茎长度",
    "use_hamming": "启用汉明距离规则",
    "min_hamming_fraction": "最小汉明距离比例",
    "use_forbidden_motifs": "启用禁用片段规则",
    "forbidden_motifs": "禁用片段",
    "candidate_pool": "每轮候选池大小",
    "attempts_per_sequence": "每条序列最大尝试次数",
}
DEFAULT_SETTINGS = {
    "length": 24,
    "count": 20,
    "gc_min": 0.40,
    "gc_max": 0.60,
    "max_homopolymer": 3,
    "use_local_gc": False,
    "local_gc_window": 8,
    "local_gc_min": 0.25,
    "local_gc_max": 0.75,
    "use_entropy": False,
    "min_entropy": 1.70,
    "use_self_complement": True,
    "max_self_complement": 5,
    "use_hairpin": True,
    "max_hairpin_stem": 4,
    "max_same_substring": 6,
    "max_cross_complement": 5,
    "scaffold_max_same_substring": 7,
    "scaffold_max_cross_complement": 7,
    "use_hamming": False,
    "min_hamming_fraction": 0.25,
    "use_forbidden_motifs": False,
    "forbidden_motifs": (),
    "candidate_pool": 16,
    "attempts_per_sequence": 20000,
}


class GenerationCancelled(Exception):
    """Raised when the UI asks the generator to stop."""


def normalized_settings(overrides=None):
    settings = dict(DEFAULT_SETTINGS)
    if overrides:
        settings.update(overrides)
    # Compatibility with reports or callers from versions that exposed a
    # configurable hairpin-loop threshold.  Hairpins are now screened only
    # by stem length.
    settings.pop("hairpin_min_loop", None)
    motifs = settings.get("forbidden_motifs", ())
    if isinstance(motifs, str):
        motifs = motifs.replace(";", ",").split(",")
    settings["forbidden_motifs"] = tuple(
        motif.strip().upper() for motif in motifs if motif.strip())
    if settings["length"] < 4:
        raise ValueError("Sequence length must be at least 4 nt")
    if settings["count"] < 1:
        raise ValueError("Sequence count must be at least 1")
    if not 0 <= settings["gc_min"] <= settings["gc_max"] <= 1:
        raise ValueError("GC limits must satisfy 0 <= minimum <= maximum <= 1")
    if settings["use_local_gc"] and not \
            0 <= settings["local_gc_min"] <= \
            settings["local_gc_max"] <= 1:
        raise ValueError(
            "Local GC limits must satisfy 0 <= minimum <= maximum <= 1")
    if settings["attempts_per_sequence"] < 1:
        raise ValueError("Attempts per sequence must be positive")
    invalid = [m for m in settings["forbidden_motifs"]
               if any(base not in DNA_BASES for base in m)]
    if settings["use_forbidden_motifs"] and invalid:
        raise ValueError("Forbidden motifs may contain only A, C, G and T: %s"
                         % ", ".join(invalid))
    return settings


def sanitize_sequence(sequence):
    """Return an uppercase A/C/G/T sequence, or ``None`` when invalid."""
    sequence = "".join(str(sequence or "").upper().split())
    if not sequence or any(base not in DNA_BASES for base in sequence):
        return None
    return sequence


def reverse_complement(sequence):
    return sequence.translate(_COMPLEMENT)[::-1]


def gc_fraction(sequence):
    return sum(base in "GC" for base in sequence) / float(len(sequence))


def max_homopolymer(sequence):
    best = current = 0
    previous = None
    for base in sequence:
        current = current + 1 if base == previous else 1
        previous = base
        best = max(best, current)
    return best


def shannon_entropy(sequence):
    counts = Counter(sequence)
    length = float(len(sequence))
    return -sum((count / length) * math.log(count / length, 2)
                for count in counts.values())


def melting_temperature(sequence, na_mm=MELTING_NA_MM,
                        mg_mm=MELTING_MG_MM,
                        strand_nm=MELTING_STRAND_NM,
                        target_nm=MELTING_TARGET_NM):
    """Return perfect-duplex DNA Tm in Celsius using nearest neighbors.

    ``strand_nm`` and ``target_nm`` are the concentrations of the two
    distinct, non-self-complementary strands.  Mg is converted to a sodium
    equivalent and the SantaLucia entropy salt correction is then applied.
    """
    sequence = sanitize_sequence(sequence)
    if sequence is None or len(sequence) < 2:
        return None
    delta_h = 0.0
    delta_s = 0.0

    # General initiation and terminal A/T or G/C contributions (DNA_NN3).
    delta_h += 0.0
    delta_s += 0.0
    for terminal in (sequence[0], sequence[-1]):
        if terminal in "AT":
            delta_h += 2.3
            delta_s += 4.1
        else:
            delta_h += 0.1
            delta_s -= 2.8

    complement = sequence.translate(_COMPLEMENT)  # template is 3' -> 5'
    for index in range(len(sequence) - 1):
        neighbors = (sequence[index:index + 2] + "/" +
                     complement[index:index + 2])
        values = _DNA_NN.get(neighbors)
        if values is None:
            values = _DNA_NN[neighbors[::-1]]
        delta_h += values[0]
        delta_s += values[1]

    # von Ahsen sodium-equivalent conversion, followed by SantaLucia's
    # entropy correction (dNTP is fixed at zero for this oligo application).
    monovalent_mm = float(na_mm)
    if mg_mm > 0:
        monovalent_mm += 120.0 * math.sqrt(float(mg_mm))
    monovalent_m = monovalent_mm * 1e-3
    delta_s += 0.368 * (len(sequence) - 1) * math.log(monovalent_m)

    concentration = (float(strand_nm) - float(target_nm) / 2.0) * 1e-9
    if concentration <= 0:
        raise ValueError("DNA strand concentrations do not define a valid Tm")
    gas_constant = 1.987  # cal/(K mol)
    return ((1000.0 * delta_h) /
            (delta_s + gas_constant * math.log(concentration)) - 273.15)


def oligo_molecular_weight(sequence):
    """Return average molar mass for unmodified 5'/3'-OH ssDNA."""
    sequence = sanitize_sequence(sequence)
    if sequence is None:
        return None
    return sum(_DNA_RESIDUE_MASS[base] for base in sequence) - 61.96


def oligo_extinction_coefficient(sequence):
    """Return ssDNA epsilon-260 using the nearest-neighbor method."""
    sequence = sanitize_sequence(sequence)
    if sequence is None:
        return None
    if len(sequence) == 1:
        return _DNA_BASE_EXTINCTION[sequence]
    doublets = sum(_DNA_DINUCLEOTIDE_EXTINCTION[sequence[index:index + 2]]
                   for index in range(len(sequence) - 1))
    internal_bases = sum(_DNA_BASE_EXTINCTION[base]
                         for base in sequence[1:-1])
    return doublets - internal_bases


def longest_common_substring(first, second):
    """Length of the longest exact contiguous substring, including edges."""
    if not first or not second:
        return 0
    shorter, longer = ((first, second) if len(first) <= len(second)
                       else (second, first))

    def has_common(size):
        if size == 0:
            return True
        words = {shorter[index:index + size]
                 for index in range(len(shorter) - size + 1)}
        return any(longer[index:index + size] in words
                   for index in range(len(longer) - size + 1))

    low, high = 0, len(shorter)
    while low < high:
        middle = (low + high + 1) // 2
        if has_common(middle):
            low = middle
        else:
            high = middle - 1
    return low


def max_hairpin_stem(sequence):
    """Return the longest non-overlapping perfect intramolecular stem."""
    best = 0
    length = len(sequence)
    complements = dict(zip(DNA_BASES, "TGCA"))
    for left in range(length):
        for right in range(left + 1, length):
            stem = 0
            while left + stem < right - stem and \
                    complements[sequence[left + stem]] == sequence[right - stem]:
                stem += 1
                best = max(best, stem)
    return best


def hamming_distance(first, second):
    if len(first) != len(second):
        return None
    return sum(left != right for left, right in zip(first, second))


def _kmers(sequence, size):
    if size <= 0 or size > len(sequence):
        return set()
    return {sequence[index:index + size]
            for index in range(len(sequence) - size + 1)}


def pair_metrics(first, second):
    hamming = hamming_distance(first, second)
    return {
        "same_substring": longest_common_substring(first, second),
        "cross_complement": longest_common_substring(
            first, reverse_complement(second)),
        "hamming": hamming,
        "hamming_fraction": (hamming / float(len(first))
                             if hamming is not None else None),
    }


def sequence_metrics(sequence, settings):
    return {
        "length": len(sequence),
        "gc": gc_fraction(sequence),
        "homopolymer": max_homopolymer(sequence),
        "entropy": shannon_entropy(sequence),
        "self_complement": longest_common_substring(
            sequence, reverse_complement(sequence)),
        "hairpin_stem": max_hairpin_stem(sequence),
    }


def _local_gc_is_valid(sequence, settings):
    window = min(int(settings["local_gc_window"]), len(sequence))
    if window <= 0:
        return True
    for index in range(len(sequence) - window + 1):
        value = gc_fraction(sequence[index:index + window])
        if value < settings["local_gc_min"] or \
                value > settings["local_gc_max"]:
            return False
    return True


def single_sequence_failures(sequence, settings):
    metrics = sequence_metrics(sequence, settings)
    failures = []
    if not settings["gc_min"] <= metrics["gc"] <= settings["gc_max"]:
        failures.append("global_gc")
    if settings["use_local_gc"] and not _local_gc_is_valid(
            sequence, settings):
        failures.append("local_gc")
    if metrics["homopolymer"] > settings["max_homopolymer"]:
        failures.append("homopolymer")
    if settings["use_entropy"] and \
            metrics["entropy"] < settings["min_entropy"]:
        failures.append("entropy")
    if settings["use_self_complement"] and \
            metrics["self_complement"] > settings["max_self_complement"]:
        failures.append("self_complement")
    if settings["use_hairpin"] and \
            metrics["hairpin_stem"] > settings["max_hairpin_stem"]:
        failures.append("hairpin")
    if settings["use_forbidden_motifs"]:
        for motif in settings["forbidden_motifs"]:
            if motif in sequence or reverse_complement(motif) in sequence:
                failures.append("forbidden_motif")
                break
    return failures, metrics


def pair_failures(first, second, settings):
    """Return hard pairwise failures using edge-complete k-mer checks."""
    failures = []
    same_size = settings["max_same_substring"] + 1
    if _kmers(first, same_size) & _kmers(second, same_size):
        failures.append("same_substring")
    complement_size = settings["max_cross_complement"] + 1
    if _kmers(first, complement_size) & \
            _kmers(reverse_complement(second), complement_size):
        failures.append("cross_complement")
    distance = hamming_distance(first, second)
    if settings["use_hamming"] and distance is not None and \
            distance / float(len(first)) < \
            settings["min_hamming_fraction"]:
        failures.append("hamming")
    return failures


def _scaffold_pair_settings(settings):
    scaffold_settings = dict(settings)
    scaffold_settings["max_same_substring"] = \
        settings["scaffold_max_same_substring"]
    scaffold_settings["max_cross_complement"] = \
        settings["scaffold_max_cross_complement"]
    return scaffold_settings


def _pair_profile(sequence, settings):
    same_size = settings["max_same_substring"] + 1
    complement_size = settings["max_cross_complement"] + 1
    return {
        "sequence": sequence,
        "settings": settings,
        "reverse_sequence": reverse_complement(sequence),
        "same_size": same_size,
        "same_kmers": _kmers(sequence, same_size),
        "complement_size": complement_size,
        "reverse_kmers": _kmers(reverse_complement(sequence),
                                 complement_size),
    }


def _profile_pair_failures(candidate, profile):
    """Fast candidate check using k-mers cached for long input/scaffolds."""
    settings = profile["settings"]
    failures = []
    if _kmers(candidate, profile["same_size"]) & profile["same_kmers"]:
        failures.append("same_substring")
    if _kmers(candidate, profile["complement_size"]) & \
            profile["reverse_kmers"]:
        failures.append("cross_complement")
    distance = hamming_distance(candidate, profile["sequence"])
    if settings["use_hamming"] and distance is not None and \
            distance / float(len(candidate)) < \
            settings["min_hamming_fraction"]:
        failures.append("hamming")
    return failures


def _orthogonality_priority(sequence, profiles):
    """Return a symmetric, non-additive priority for both pairwise risks.

    The same-substring and cross-complement objectives remain independent:
    first minimize whichever objective is closest to its own hard threshold,
    then minimize the other one.  Raw worst-case lengths break utilization
    ties.  Consequently, an excellent value in one objective can never
    compensate for a near-threshold value in the other.
    """
    worst_same = 0
    worst_cross = 0
    worst_same_utilization = 0.0
    worst_cross_utilization = 0.0
    for profile in profiles:
        same = longest_common_substring(sequence, profile["sequence"])
        cross = longest_common_substring(
            sequence, profile["reverse_sequence"])
        settings = profile["settings"]
        same_limit = settings["max_same_substring"]
        cross_limit = settings["max_cross_complement"]
        worst_same = max(worst_same, same)
        worst_cross = max(worst_cross, cross)
        worst_same_utilization = max(
            worst_same_utilization, same / float(max(1, same_limit)))
        worst_cross_utilization = max(
            worst_cross_utilization, cross / float(max(1, cross_limit)))

    # Sort each objective pair before comparison so neither same-substring nor
    # cross-complement receives a fixed tie-breaking advantage.
    utilization = sorted(
        (worst_same_utilization, worst_cross_utilization), reverse=True)
    raw_lengths = sorted((worst_same, worst_cross), reverse=True)
    return (-utilization[0], -utilization[1],
            -raw_lengths[0], -raw_lengths[1])


def _candidate_score(sequence, peers, metrics, settings, profiles=None):
    if profiles is None:
        profiles = [_pair_profile(peer, settings) for peer in peers]
    same_length_distances = [hamming_distance(sequence, peer)
                             for peer in peers if len(peer) == len(sequence)]
    min_distance = min(same_length_distances) if same_length_distances else \
        len(sequence)
    target_gc = (settings["gc_min"] + settings["gc_max"]) / 2.0
    score = list(_orthogonality_priority(sequence, profiles))
    if settings["use_hamming"]:
        score.append(min_distance)
    score.append(-abs(metrics["gc"] - target_gc))
    if settings["use_entropy"]:
        score.append(metrics["entropy"])
    return tuple(score)


def generate_sequences(overrides=None, background_sequences=(),
                       scaffold_sequences=(),
                       progress=None, cancelled=None):
    """Generate a diverse orthogonal set using bounded pool-based selection.

    The best candidate from a small valid pool is selected at every step,
    rather than permanently accepting the first random hit.  This retains the
    simplicity of rejection sampling while reducing order-dependent collapse.
    """
    settings = normalized_settings(overrides)
    input_sequences = []
    for raw_sequence in background_sequences:
        sequence = sanitize_sequence(raw_sequence)
        if sequence:
            input_sequences.append(sequence)
    # Duplicate input rows remain visible in the report, while one copy is
    # enough for candidate rejection.
    background = list(dict.fromkeys(input_sequences))
    scaffolds = []
    for index, item in enumerate(scaffold_sequences, 1):
        if isinstance(item, (tuple, list)) and len(item) == 2:
            name, raw_sequence = item
        else:
            name, raw_sequence = "骨架链-%03d" % index, item
        sequence = sanitize_sequence(raw_sequence)
        if sequence:
            scaffolds.append((str(name), sequence))
    regular_profiles = [_pair_profile(sequence, settings)
                        for sequence in background]
    scaffold_settings = _scaffold_pair_settings(settings)
    scaffold_profiles = [_pair_profile(sequence, scaffold_settings)
                         for unused_name, sequence in scaffolds]
    profiles = regular_profiles + scaffold_profiles
    # SystemRandom deliberately gives a fresh result on every run, even when
    # every visible parameter is identical.
    rng = random.SystemRandom()
    selected = []
    rejection_counts = Counter()
    total_attempts = 0
    complete = True

    for sequence_index in range(settings["count"]):
        pool = []
        for unused_attempt in range(settings["attempts_per_sequence"]):
            if cancelled and cancelled():
                raise GenerationCancelled()
            total_attempts += 1
            candidate = "".join(rng.choice(DNA_BASES)
                                for unused in range(settings["length"]))
            failures, metrics = single_sequence_failures(candidate, settings)
            if failures:
                rejection_counts.update(failures)
            else:
                pair_error = []
                for profile in profiles:
                    pair_error.extend(_profile_pair_failures(
                        candidate, profile))
                    if pair_error:
                        break
                if pair_error:
                    rejection_counts.update(pair_error)
                else:
                    pool.append((
                        _candidate_score(candidate, background + selected,
                                         metrics, settings, profiles),
                        candidate, metrics))
                    if len(pool) >= settings["candidate_pool"]:
                        break
            if progress and total_attempts % 100 == 0:
                progress(sequence_index, total_attempts)
        if not pool:
            complete = False
            break
        unused_score, sequence, metrics = max(pool, key=lambda item: item[0])
        selected.append(sequence)
        profiles.append(_pair_profile(sequence, settings))
        if progress:
            progress(len(selected), total_attempts)

    # Input and scaffold order are preserved.  Newly generated oligos are
    # numbered and exported from the lowest intended-duplex Tm to the highest.
    selected.sort(key=melting_temperature)
    rows = _sequence_report_rows(
        input_sequences, scaffolds, selected, settings)
    pairs = _pair_report_rows(
        input_sequences, scaffolds, selected, settings)
    return {
        "sequences": selected,
        "input_sequences": input_sequences,
        "scaffold_sequences": scaffolds,
        "rows": rows,
        "pairs": pairs,
        "settings": settings,
        "attempts": total_attempts,
        "rejections": dict(rejection_counts),
        "background_count": len(input_sequences),
        "complete": complete and len(selected) == settings["count"],
    }


SEQUENCE_HEADERS = (
    "来源", "名称", "序列（5′→3′）", "长度", "GC（%）",
    "互补序列（5′→3′）", "熔解温度（°C）", "分子量（g/mol）",
    "消光系数 ε260（L·mol⁻¹·cm⁻¹）",
    "最长连续相同碱基",
    "输入/新生成链间最差同向相同片段（nt）",
    "输入/新生成链间最差链间互补片段（nt）",
    "与骨架链最差同向相同片段（nt）",
    "与骨架链最差链间互补片段（nt）",
    "自身互补长度（nt）", "发卡茎长度（bp）",
    "最小汉明距离（%）", "序列熵（bits）", "状态", "问题说明",
)
PAIR_HEADERS = (
    "来源1", "序列1", "来源2", "序列2",
    "同向相同片段（nt）", "链间互补片段（nt）",
    "汉明距离", "汉明距离（%）",
    "状态", "问题说明",
)
SETTING_HEADERS = ("参数", "数值")


def _named_entries(input_sequences, scaffold_sequences, generated_sequences):
    entries = []
    entries.extend(("输入", "输入-%03d" % index, sequence)
                   for index, sequence in enumerate(input_sequences, 1))
    entries.extend(("骨架链", name, sequence)
                   for name, sequence in scaffold_sequences)
    entries.extend(("新生成", "新序列-%03d" % index, sequence)
                   for index, sequence in enumerate(generated_sequences, 1))
    return entries


def _entry_pair_settings(first_entry, second_entry, settings):
    if first_entry[0] == "骨架链" or second_entry[0] == "骨架链":
        return _scaffold_pair_settings(settings)
    return settings


def _sequence_report_rows(input_sequences, scaffold_sequences,
                          generated_sequences, settings):
    rows = []
    entries = _named_entries(
        input_sequences, scaffold_sequences, generated_sequences)
    for entry_index, (source, name, sequence) in enumerate(entries):
        if source == "骨架链":
            # A multi-kilobase scaffold is a fixed background, not a candidate.
            # Hairpin enumeration is both biologically unhelpful at this scale
            # and quadratic, so only composition and cross-sequence metrics are
            # reported for this row.
            metrics = {
                "gc": gc_fraction(sequence),
                "homopolymer": max_homopolymer(sequence),
                "entropy": shannon_entropy(sequence),
                "self_complement": "",
                "hairpin_stem": "",
            }
            own_failures = []
            complement_sequence = ""
            melt_temp = ""
            molecular_weight = ""
            extinction_coefficient = ""
        else:
            metrics = sequence_metrics(sequence, settings)
            own_failures, unused_metrics = single_sequence_failures(
                sequence, settings)
            complement_sequence = reverse_complement(sequence)
            calculated_tm = melting_temperature(sequence)
            melt_temp = (round(calculated_tm, 2)
                         if calculated_tm is not None else "")
            molecular_weight = round(oligo_molecular_weight(sequence), 2)
            extinction_coefficient = oligo_extinction_coefficient(sequence)
        oligo_entries = [
            entry for peer_index, entry in enumerate(entries)
            if peer_index != entry_index and entry[0] != "骨架链"]
        scaffold_entries = [entry for entry in entries
                            if entry[0] == "骨架链"]

        if source == "骨架链":
            # Scaffold is a fixed background.  Its own row identifies the
            # selected sequence, while pair diagnostics belong to each input
            # or newly generated oligo row.
            oligo_values = []
            scaffold_values = []
            failure_messages = []
            has_failures = False
        else:
            oligo_values = [pair_metrics(sequence, entry[2])
                            for entry in oligo_entries]
            scaffold_values = [pair_metrics(sequence, entry[2])
                               for entry in scaffold_entries]
            oligo_failures = []
            for peer_entry in oligo_entries:
                oligo_failures.extend(pair_failures(
                    sequence, peer_entry[2], settings))
            scaffold_failures = []
            scaffold_settings = _scaffold_pair_settings(settings)
            for scaffold_entry in scaffold_entries:
                scaffold_failures.extend(pair_failures(
                    sequence, scaffold_entry[2], scaffold_settings))
            failure_messages = [
                ISSUE_LABELS_CN.get(code, code) for code in own_failures]
            failure_messages.extend(
                "输入/新生成链间：%s" % ISSUE_LABELS_CN.get(code, code)
                for code in oligo_failures)
            failure_messages.extend(
                "与骨架链：%s" % ISSUE_LABELS_CN.get(code, code)
                for code in scaffold_failures)
            has_failures = bool(failure_messages)

        def worst(values, key, empty):
            return max([value[key] for value in values] or [empty])

        oligo_empty = "" if source == "骨架链" else 0
        oligo_same = worst(oligo_values, "same_substring", oligo_empty)
        oligo_cross = worst(oligo_values, "cross_complement", oligo_empty)
        scaffold_same = worst(scaffold_values, "same_substring", "")
        scaffold_cross = worst(scaffold_values, "cross_complement", "")
        distance_fractions = [
            value["hamming_fraction"] for value in oligo_values
            if value["hamming_fraction"] is not None]
        failure_messages = sorted(set(failure_messages))
        rows.append((
            source, name, sequence, len(sequence),
            round(metrics["gc"] * 100.0, 2),
            complement_sequence, melt_temp, molecular_weight,
            extinction_coefficient,
            metrics["homopolymer"],
            oligo_same, oligo_cross,
            scaffold_same, scaffold_cross,
            metrics["self_complement"], metrics["hairpin_stem"],
            round(min(distance_fractions) * 100.0, 2)
            if distance_fractions else "",
            round(metrics["entropy"], 3),
            "通过" if not has_failures else "需复核",
            "；".join(failure_messages),
        ))
    return rows


def _pair_report_rows(input_sequences, scaffold_sequences,
                      generated_sequences, settings):
    rows = []
    # The scaffold is a screening background, not an oligo that needs a
    # pair-by-pair purchasing report.  Keep enforcing scaffold orthogonality
    # during generation and in the per-sequence summary, but omit scaffold
    # combinations from the "两两分析" worksheet.
    entries = _named_entries(
        input_sequences, (), generated_sequences)
    for first_index, first_entry in enumerate(entries):
        for second_index in range(first_index + 1, len(entries)):
            second_entry = entries[second_index]
            metrics = pair_metrics(first_entry[2], second_entry[2])
            failures = pair_failures(
                first_entry[2], second_entry[2],
                _entry_pair_settings(first_entry, second_entry, settings))
            rows.append((
                first_entry[0], first_entry[1],
                second_entry[0], second_entry[1],
                metrics["same_substring"], metrics["cross_complement"],
                metrics["hamming"],
                round(metrics["hamming_fraction"] * 100.0, 2)
                if metrics["hamming_fraction"] is not None else "",
                "通过" if not failures else "需复核",
                "；".join(ISSUE_LABELS_CN.get(code, code)
                         for code in failures),
            ))
    return rows


def read_sequence_text(filename):
    """Read a one-sequence-per-line TXT file with detailed line errors."""
    sequences = []
    errors = []
    with open(filename, "r", encoding="utf-8-sig") as input_file:
        for line_number, raw_line in enumerate(input_file, 1):
            value = raw_line.strip().upper()
            if not value or value.startswith("#"):
                continue
            if any(base not in DNA_BASES for base in value):
                errors.append(
                    "第%d行包含A、C、G、T以外的字符：%s"
                    % (line_number, raw_line.rstrip()))
                continue
            sequences.append(value)
    return sequences, errors


def write_input_template(filename):
    """Write a comment-only TXT template; comments and blanks are ignored."""
    content = (
        "# 正交序列输入模板\n"
        "# 请在下方每行粘贴一条只包含A、C、G、T的已有序列。\n"
        "# 空行及以#开头的说明行会被忽略。\n"
        "# 不同输入序列可以具有不同长度。\n\n")
    with open(filename, "w", encoding="utf-8") as output_file:
        output_file.write(content)


def _column_name(index):
    name = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        name = chr(65 + remainder) + name
    return name


def _xlsx_cell(reference, value, style_id=0):
    style = ' s="%d"' % style_id if style_id else ""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return '<c r="%s"%s><v>%s</v></c>' % (reference, style, value)
    text = escape(str(value))
    return ('<c r="%s" t="inlineStr"%s><is><t>%s</t></is></c>' %
            (reference, style, text))


def _xlsx_sheet(headers, rows, widths, row_styles=None):
    all_rows = [headers] + list(rows)
    xml_rows = []
    for row_number, row in enumerate(all_rows, 1):
        style_id = (1 if row_number == 1 else
                    (row_styles or {}).get(row[0], 0))
        cells = [_xlsx_cell("%s%d" % (_column_name(column), row_number), value,
                            style_id)
                 for column, value in enumerate(row, 1)]
        xml_rows.append('<row r="%d">%s</row>' %
                        (row_number, "".join(cells)))
    columns = "".join(
        '<col min="%d" max="%d" width="%s" customWidth="1"/>' %
        (index, index, width) for index, width in enumerate(widths, 1))
    last_column = _column_name(len(headers))
    last_row = len(all_rows)
    return '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<dimension ref="A1:%s%d"/><sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>
<cols>%s</cols><sheetData>%s</sheetData><autoFilter ref="A1:%s%d"/>
</worksheet>''' % (last_column, last_row, columns, "".join(xml_rows),
                    last_column, last_row)


def write_orthogonal_workbook(filename, result):
    """Write sequences, pairwise diagnostics and exact settings to XLSX."""
    settings = result["settings"]
    setting_rows = [
        ("筛选模式", "基础规则＋可选高级规则"),
        ("请求生成数量", settings["count"]),
        ("实际生成数量", len(result["sequences"])),
        ("输入序列数量", len(result.get("input_sequences", ()))),
        ("输入文件", result.get("input_file", "") or "无"),
        ("骨架链数量", len(result.get("scaffold_sequences", ()))),
        ("骨架链名称", "，".join(
            name for name, unused_sequence in
            result.get("scaffold_sequences", ())) or "无"),
        ("熔解温度模型", "SantaLucia DNA最近邻模型"),
        ("熔解温度Na⁺浓度", "50 mM"),
        ("熔解温度Mg²⁺浓度", "10 mM"),
        ("熔解温度链浓度", "100 nM"),
        ("熔解温度互补目标链浓度", "100 nM"),
        ("是否完整生成", "是" if result["complete"] else "否"),
        ("随机方式", "每次运行使用新的系统随机数"),
        ("候选评价次数", result["attempts"]),
    ]
    for key in sorted(settings):
        value = settings[key]
        if isinstance(value, (tuple, list)):
            value = "，".join(value)
        elif isinstance(value, bool):
            value = "启用" if value else "未启用"
        setting_rows.append((SETTING_LABELS_CN.get(key, key), value))
    for key, value in sorted(result["rejections"].items()):
        setting_rows.append(("因“%s”淘汰的候选数" %
                             ISSUE_LABELS_CN.get(key, key), value))

    sequence_widths = (
        13, 15, 52, 9, 10, 52, 17, 19, 30, 16,
        31, 32, 27, 28,
        20, 18, 20, 15, 12, 46)
    sequence_headers = SEQUENCE_HEADERS
    sequence_rows = result["rows"]
    if not result.get("scaffold_sequences"):
        # Without a selected scaffold, omit the two meaningless scaffold
        # comparison columns entirely instead of exporting blank columns.
        scaffold_columns = {12, 13}
        keep = [index for index in range(len(SEQUENCE_HEADERS))
                if index not in scaffold_columns]
        sequence_headers = tuple(SEQUENCE_HEADERS[index] for index in keep)
        sequence_widths = tuple(sequence_widths[index] for index in keep)
        sequence_rows = [tuple(row[index] for index in keep)
                         for row in result["rows"]]

    source_row_styles = {"输入": 2, "骨架链": 3, "新生成": 4}
    sheets = (
        ("序列分析", sequence_headers, sequence_rows, sequence_widths,
         source_row_styles),
        ("两两分析", PAIR_HEADERS, result["pairs"],
         (13, 15, 13, 15, 20, 22, 18, 14, 12, 38), None),
        ("设置", SETTING_HEADERS, setting_rows, (38, 60), None),
    )
    content_overrides = "".join(
        '<Override PartName="/xl/worksheets/sheet%d.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>' % index
        for index in range(1, len(sheets) + 1))
    content_types = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>%s<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/></Types>''' % content_overrides
    package_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>'''
    workbook_sheets = "".join(
        '<sheet name="%s" sheetId="%d" r:id="rId%d"/>' %
        (escape(sheet[0]), index, index)
        for index, sheet in enumerate(sheets, 1))
    workbook = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets>%s</sheets></workbook>''' % workbook_sheets
    sheet_rels = "".join(
        '<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet%d.xml"/>' % (index, index)
        for index in range(1, len(sheets) + 1))
    workbook_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">%s<Relationship Id="rId%d" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>''' % (sheet_rels, len(sheets) + 1)
    styles = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<fonts count="2"><font/><font><b/><color rgb="FFFFFFFF"/></font></fonts>
<fills count="6">
<fill><patternFill patternType="none"/></fill>
<fill><patternFill patternType="gray125"/></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FF285C8E"/><bgColor indexed="64"/></patternFill></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FFEAF3FF"/><bgColor indexed="64"/></patternFill></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FFF3F1ED"/><bgColor indexed="64"/></patternFill></fill>
<fill><patternFill patternType="solid"><fgColor rgb="FFECF8EE"/><bgColor indexed="64"/></patternFill></fill>
</fills>
<borders count="1"><border/></borders>
<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>
<cellXfs count="5">
<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>
<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/>
<xf numFmtId="0" fontId="0" fillId="3" borderId="0" xfId="0" applyFill="1"/>
<xf numFmtId="0" fontId="0" fillId="4" borderId="0" xfId="0" applyFill="1"/>
<xf numFmtId="0" fontId="0" fillId="5" borderId="0" xfId="0" applyFill="1"/>
</cellXfs>
<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>
</styleSheet>'''
    with ZipFile(filename, "w", ZIP_DEFLATED) as workbook_file:
        workbook_file.writestr("[Content_Types].xml", content_types)
        workbook_file.writestr("_rels/.rels", package_rels)
        workbook_file.writestr("xl/workbook.xml", workbook)
        workbook_file.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        workbook_file.writestr("xl/styles.xml", styles)
        for index, unused_sheet in enumerate(sheets, 1):
            name, headers, rows, widths, row_styles = unused_sheet
            workbook_file.writestr(
                "xl/worksheets/sheet%d.xml" % index,
                _xlsx_sheet(headers, rows, widths, row_styles))
