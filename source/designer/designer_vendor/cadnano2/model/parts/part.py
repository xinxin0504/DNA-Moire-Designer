#!/usr/bin/env python
# encoding: utf-8

from heapq import heapify, heappush, heappop
from itertools import product
from collections import defaultdict
import random
import re

from cadnano2.model.enum import StrandType
from cadnano2.model.virtualhelix import VirtualHelix
from cadnano2.model.strand import Strand
from cadnano2.model.oligo import Oligo
from cadnano2.model.strandset import StrandSet
from cadnano2.views import styles

import cadnano2.util as util

util.qtWrapImport('QtCore', globals(), ['pyqtSignal', 'QObject'])
util.qtWrapImport('QtGui', globals(), ['QUndoCommand'])


def _sequenceRowSortKey(record):
    """Group exported rows by color, then sort their 5' start numerically."""
    color = str(record[4]).strip().lower()
    start = str(record[0]).strip()
    # Hybrid records add a lattice prefix (``H:`` or ``S:``).  Parsing the
    # complete bracketed integer prevents 100 from sorting before 20.
    match = re.match(r'^(?:[HS]\s*:\s*)?(-?\d+)\s*\[\s*(-?\d+)\s*\]$',
                     start, re.IGNORECASE)
    if match is not None:
        helixNumber = int(match.group(1))
        baseIndex = int(match.group(2))
    else:
        # sequenceEndpoints() always emits ``helix[base]``.  Keep malformed
        # third-party records deterministic instead of breaking an export.
        helixNumber = float('inf')
        baseIndex = float('inf')
    return (color, helixNumber, baseIndex, start)


_SCAFFOLD_STAPLE_PAIR_EXCLUSION_DISTANCE = 10


def _circularIndexDistance(first, second, circularSize):
    """Return the shortest index distance on a circular helix."""
    direct = abs(int(first) - int(second)) % int(circularSize)
    return min(direct, int(circularSize) - direct)


def _indexInsideExactPairInterval(index, positions, intervalBases,
                                  circularSize=None):
    """Return whether *index* lies inside an exact same-pair bp interval.

    ``positions`` must contain scaffold-crossover coordinates for one
    unordered helix pair only.  In a circular design, the final coordinate
    and the first coordinate are consecutive across the periodic boundary.
    cadnano spacing is inclusive: indices 0 and 31 delimit a 32-bp interval,
    while indices 0 and 63 delimit a 64-bp interval.
    Endpoints themselves are excluded because the ordinary scaffold
    clearance rule handles them.
    """
    intervalBases = int(intervalBases)
    indexDelta = intervalBases - 1
    if circularSize:
        circularSize = int(circularSize)
        ordered = sorted(set(int(position) % circularSize
                             for position in positions))
        if len(ordered) < 2:
            return False
        candidate = int(index) % circularSize
        for left, right in zip(ordered, ordered[1:] +
                               [ordered[0] + circularSize]):
            if right - left != indexDelta:
                continue
            offset = (candidate - left) % circularSize
            if 0 < offset < indexDelta:
                return True
        return False

    ordered = sorted(set(int(position) for position in positions))
    candidate = int(index)
    return any(right - left == indexDelta and left < candidate < right
               for left, right in zip(ordered, ordered[1:]))


def _squareStapleCrossoverIsOverdense(part, fromHelix, toHelix, index,
                                       scaffoldPositionsByPairSide,
                                       circularSize=None,
                                       scaffoldPositionsByHelix=None,
                                       scaffoldPositionsByPair=None):
    """Apply Square scaffold/staple conflicts without losing 1/32 sites.

    Ordinary AutoCS suppresses a native staple site less than 10 bases from
    a scaffold crossover on the same helix pair.  A stricter Square rule has
    priority: two consecutive scaffold crossovers on that pair exactly 32 bp
    apart already provide the maximum total crossover density, so no staple
    crossover may be placed between them.  At a 64-bp (or larger) scaffold
    interval, legal internal staple sites remain available after the normal
    10-bp endpoint-clearance filter.  Curved routes apply the same rules with
    circular coordinate distance.
    """
    # The density family is one *unordered neighbouring helix pair*.  A
    # scaffold crossover from either helix to some other neighbour must not
    # suppress this pair's legal staple site merely because it shares one
    # helix.  This is essential when sparse scaffold routing (for example
    # 1/64) is completed to the normal 1/32 staple register.
    pair = tuple(sorted((fromHelix, toHelix)))
    positions = sorted((scaffoldPositionsByPair or {}).get(pair, ()))
    if _indexInsideExactPairInterval(index, positions, 32, circularSize):
        return True
    if circularSize:
        return any(_circularIndexDistance(
                       index, position, circularSize) <
                   _SCAFFOLD_STAPLE_PAIR_EXCLUSION_DISTANCE
                   for position in positions)
    if any(abs(index - position) <
           _SCAFFOLD_STAPLE_PAIR_EXCLUSION_DISTANCE
           for position in positions):
        return True

    return False


def _honeycombStapleCrossoverIsOverdense(part, fromHelix, toHelix, index,
                                          scaffoldPositionsByPairSide,
                                          circularSize=None,
                                          scaffoldPositionsByHelix=None,
                                          scaffoldPositionsByPair=None):
    """Return whether scaffold xovers leave no legal 21-bp slot.

    The density family is the same unordered neighbouring helix pair.
    Scaffold crossovers to a different neighbour must not suppress normal
    AutoCS staple crossovers on this pair.
    """
    pair = tuple(sorted((fromHelix, toHelix)))
    positions = (scaffoldPositionsByPair or {}).get(pair, ())
    if circularSize:
        return any(_circularIndexDistance(
                       index, position, circularSize) <
                   _SCAFFOLD_STAPLE_PAIR_EXCLUSION_DISTANCE
                   for position in positions)
    return any(abs(index - position) <
               _SCAFFOLD_STAPLE_PAIR_EXCLUSION_DISTANCE
               for position in positions)


def _stapleCrossoverIsOverdense(part, fromHelix, toHelix, index,
                                scaffoldPositionsByPairSide,
                                circularSize=None,
                                scaffoldPositionsByHelix=None,
                                scaffoldPositionsByPair=None):
    """Apply the active lattice's local staple-crossover density rule."""
    if part._step == 32:
        return _squareStapleCrossoverIsOverdense(
            part, fromHelix, toHelix, index, scaffoldPositionsByPairSide,
            circularSize, scaffoldPositionsByHelix,
            scaffoldPositionsByPair)
    if part._step == 21:
        return _honeycombStapleCrossoverIsOverdense(
            part, fromHelix, toHelix, index, scaffoldPositionsByPairSide,
            circularSize, scaffoldPositionsByHelix,
            scaffoldPositionsByPair)
    return False


def _crossoverPositionsByHelix(part, strandType):
    """Return all crossover endpoint indices, grouped by virtual helix."""
    positions = defaultdict(set)
    for vh in part.getVirtualHelices():
        strandSet = vh.getStrandSetByType(strandType)
        for strand in strandSet:
            if strand.connection5p() is not None and \
                    strand.connection5p().part() is part:
                positions[vh.number()].add(strand.idx5Prime())
            if strand.connection3p() is not None and \
                    strand.connection3p().part() is part:
                positions[vh.number()].add(strand.idx3Prime())
    return positions


def _stapleOligoBaseRecords(oligo):
    """Return one record per design position in an oligo's 5'-to-3' order."""
    records = []
    for strand in oligo.strand5p().generator3pStrand():
        vh = strand.virtualHelix()
        step = 1 if strand.isDrawn5to3() else -1
        insertionMap = strand.part().insertions().get(vh.coord(), {})
        for index in range(strand.idx5Prime(), strand.idx3Prime() + step,
                           step):
            insertion = insertionMap.get(index)
            actualLength = 1 + (insertion.length() if insertion else 0)
            records.append((vh.number(), index, strand,
                            max(0, actualLength)))
    return records


def _hasContinuousPositionRun(records, minimum=16):
    """Test actual consecutive bases, including insertions and deletions.

    Insertions add their real nucleotide count to the continuous region;
    deletions contribute zero.  A helix change, coordinate gap, crossover,
    or nick still starts a new run.
    """
    longest = run = 0
    previous = None
    for helix, index, unused_strand, actual_length in records:
        if previous is not None and previous[0] == helix and \
                abs(previous[1] - index) == 1:
            run += max(0, int(actual_length))
        else:
            run = max(0, int(actual_length))
        longest = max(longest, run)
        previous = (helix, index)
    return longest >= minimum


def _legalStapleNickBoundaries(records, stapleXovers, scaffoldXovers,
                                minimumIndexDistance=7,
                                preferredPhase=8,
                                ignoreIndels=False):
    """Return legal native boundaries keyed by the right-hand list offset."""
    candidates = {}
    for offset in range(1, len(records)):
        left = records[offset - 1]
        right = records[offset]
        if left[2] is not right[2] or left[0] != right[0] or \
                abs(left[1] - right[1]) != 1:
            continue
        helix = left[0]
        lowIndex, highIndex = sorted((left[1], right[1]))
        strandSet = left[2].strandSet()
        # A nick is the boundary between lowIndex and highIndex.  Keep both
        # endpoint bases free of insertions and deletions/skips: placing an
        # indel on either side makes that structural adjustment coincide with
        # a staple terminus even though the crossover-only checks below pass.
        insertionMap = left[2].part().insertions().get(
            left[2].virtualHelix().coord(), {})
        if not ignoreIndels and (lowIndex in insertionMap or
                                 highIndex in insertionMap):
            continue
        splitIndex = lowIndex if strandSet.isDrawn5to3() else highIndex
        if not strandSet.strandCanBeSplit(left[2], splitIndex):
            continue
        if lowIndex in scaffoldXovers.get(helix, ()) or \
                highIndex in scaffoldXovers.get(helix, ()):
            continue
        if any(min(abs(lowIndex - xover), abs(highIndex - xover)) <
               minimumIndexDistance
               for xover in stapleXovers.get(helix, ())):
            continue
        scaffoldDistances = [
            min(abs(lowIndex - xover), abs(highIndex - xover))
            for xover in scaffoldXovers.get(helix, ())]
        # The actual nearest distance is retained as a late tie-breaker.  If
        # two sites both satisfy phase and >3-nt clearance, this maximizes the
        # smaller of their distances to the surrounding scaffold xovers.
        scaffoldClearance = (min(scaffoldDistances) if scaffoldDistances
                             else len(records) + 1)
        farFromScaffold = scaffoldClearance > 3
        candidates[offset] = (helix, highIndex,
                              highIndex % preferredPhase == 0,
                              farFromScaffold, scaffoldClearance)
    return candidates


def _occupiedEdgeIndices(strandSet):
    """Return low/high indices of every contiguous occupied staple run."""
    intervals = sorted(strand.idxs() for strand in strandSet)
    merged = []
    for low, high in intervals:
        if merged and low <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], high)
        else:
            merged.append([low, high])
    return set(index for low, high in merged for index in (low, high))


def _linearStapleBreakPlan(records, boundaries, closingBoundary=None,
                           continuousMinimum=16, softMinimum=30,
                           terminalMaximum=49, preferredMinimum=30,
                           preferredMaximum=50, targetLength=40,
                           preferDeletionDense=False,
                           requireDeletionDenseMinimum=False,
                           hardMaximum=57,
                           densePreferredMinimum=40,
                           densePreferredMaximum=60,
                           denseTargetLength=50,
                           denseTerminalMaximum=60):
    """Find the best 21--``hardMaximum`` nt partition.

    Normal Autobreak uses 57 nt.  Curved/Frame deletion-dense products use
    60 nt so that their 40--60 nt preferred range is fully available.  The
    caller may make one explicit 64-nt fallback pass only after the normal
    problem has no legal solution.
    """
    count = len(records)
    edgeIndicesByStrandSet = {}
    atPhysicalEdge = []
    for unused_helix, index, strand, unused_length in records:
        strandSet = strand.strandSet()
        if strandSet not in edgeIndicesByStrandSet:
            edgeIndicesByStrandSet[strandSet] = \
                                        _occupiedEdgeIndices(strandSet)
        atPhysicalEdge.append(index in edgeIndicesByStrandSet[strandSet])
    # states[position][segment_count] = (additive score, chosen boundaries)
    states = [dict() for unused_index in range(count + 1)]
    states[0][0] = ((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), [])
    for start in range(count):
        if not states[start]:
            continue
        actualLength = 0
        for end in range(start + 1, count + 1):
            actualLength += records[end - 1][3]
            if actualLength > hardMaximum:
                break
            if actualLength < 21:
                continue
            if end != count and end not in boundaries:
                continue
            segment = records[start:end]
            deletionDenseSegment = int(
                preferDeletionDense and
                sum(int(record[3]) == 0 for record in segment) >= 2)
            # Curved/Frame first exhaust a strict local problem in which
            # every final product containing at least two deletions is at
            # least 40 nt.  Only an explicit later relaxation may retain a
            # shorter dense product when immutable topology makes 40 nt
            # impossible.
            if deletionDenseSegment and requireDeletionDenseMinimum and \
                    actualLength < densePreferredMinimum:
                continue
            # In the normal pass, ordinary segments end at 57 nt while only
            # a segment which itself contains at least two deletions may use
            # the deletion-dense 60-nt ceiling.  A caller's explicit 64-nt
            # pass is the sole no-solution exception and therefore bypasses
            # these two normal ceilings.
            if hardMaximum <= densePreferredMaximum:
                segmentMaximum = (densePreferredMaximum
                                  if deletionDenseSegment else 57)
                if actualLength > segmentMaximum:
                    continue
            segmentSoftMinimum = (densePreferredMinimum
                                  if deletionDenseSegment else softMinimum)
            segmentPreferredMinimum = (densePreferredMinimum
                                       if deletionDenseSegment
                                       else preferredMinimum)
            segmentPreferredMaximum = (densePreferredMaximum
                                       if deletionDenseSegment
                                       else preferredMaximum)
            segmentTargetLength = (denseTargetLength
                                   if deletionDenseSegment
                                   else targetLength)
            segmentTerminalMaximum = (denseTerminalMaximum
                                      if deletionDenseSegment
                                      else terminalMaximum)
            # A final 5'/3' end at a physical occupied-run edge obeys the
            # local segment's terminal maximum.  Internal nick endpoints are
            # not occupied-run edges and therefore do not trigger this rule.
            if actualLength > segmentTerminalMaximum and \
                    (atPhysicalEdge[start] or atPhysicalEdge[end - 1]):
                continue

            hasRun = int(_hasContinuousPositionRun(
                                        segment, continuousMinimum))
            meetsSoftMinimum = int(actualLength >= segmentSoftMinimum)
            preferredLength = int(
                segmentPreferredMinimum <= actualLength <=
                segmentPreferredMaximum)
            outsideDistance = max(
                0, segmentPreferredMinimum - actualLength,
                actualLength - segmentPreferredMaximum)
            deviation = abs(actualLength - segmentTargetLength)
            deletionDensePreferred = int(
                deletionDenseSegment and
                densePreferredMinimum <= actualLength <=
                densePreferredMaximum)
            preferredCut = 0
            scaffoldSafeCut = 0
            scaffoldClearance = 0
            chosenBoundary = None
            if end != count:
                chosenBoundary = end
                preferredCut = int(boundaries[end][2])
                scaffoldSafeCut = int(boundaries[end][3])
                scaffoldClearance = boundaries[end][4]

            for segmentCount, (score, cuts) in list(states[start].items()):
                newCount = segmentCount + 1
                newScore = (score[0] + hasRun,
                            score[1] + scaffoldSafeCut + preferredCut,
                            score[2] + scaffoldSafeCut,
                            score[3] + preferredCut,
                            score[4] + meetsSoftMinimum,
                            score[5] + preferredLength,
                            score[6] + scaffoldClearance,
                            score[7] - outsideDistance,
                            score[8] - deviation,
                            score[9] + deletionDensePreferred,
                            score[10] + deletionDenseSegment)
                newCuts = cuts + ([chosenBoundary]
                                  if chosenBoundary is not None else [])
                old = states[end].get(newCount)
                if old is None or newScore > old[0]:
                    states[end][newCount] = (newScore, newCuts)

    best = None
    for segmentCount, (score, cuts) in states[count].items():
        cutCount = segmentCount - 1
        combinedPlacement = score[1]
        scaffoldSafeCuts = score[2]
        preferredCuts = score[3]
        softMinimumLengths = score[4]
        preferredLengths = score[5]
        totalScaffoldClearance = score[6]
        if closingBoundary is not None:
            cutCount += 1
            combinedPlacement += (int(closingBoundary[3]) +
                                  int(closingBoundary[2]))
            scaffoldSafeCuts += int(closingBoundary[3])
            preferredCuts += int(closingBoundary[2])
            totalScaffoldClearance += closingBoundary[4]
        # Compare proportions first, then avoid unnecessary fragmentation.
        denseKey = (float(score[9]) / score[10]
                    if score[10] else 1.0)
        key = (denseKey,
               float(score[0]) / segmentCount,
               float(preferredLengths) / segmentCount,
               float(softMinimumLengths) / segmentCount,
               float(combinedPlacement) / max(1, 2 * cutCount),
               float(scaffoldSafeCuts) / max(1, cutCount),
               float(preferredCuts) / max(1, cutCount),
               -segmentCount,
               float(totalScaffoldClearance) / max(1, cutCount),
               score[7], score[8])
        if best is None or key > best[0]:
            best = (key, cuts)
    return None if best is None else best[1]


def _bestStapleBreakPlan(oligo, stapleXovers, scaffoldXovers,
                         minimumIndexDistance=7, continuousMinimum=16,
                         preferredPhase=8, softMinimum=30,
                         ignoreIndels=False, preferredMinimum=30,
                         preferredMaximum=50, targetLength=40,
                         terminalMaximum=49,
                         preferDeletionDense=False,
                         requireDeletionDenseMinimum=False,
                         hardMaximum=57,
                         densePreferredMinimum=40,
                         densePreferredMaximum=60,
                         denseTargetLength=50,
                         denseTerminalMaximum=60):
    """Return ``(helix, upper_index)`` nick positions for one staple oligo."""
    records = _stapleOligoBaseRecords(oligo)
    if not records:
        return []
    candidates = _legalStapleNickBoundaries(
                            records, stapleXovers, scaffoldXovers,
                            minimumIndexDistance, preferredPhase,
                            ignoreIndels=ignoreIndels)
    if not oligo.isLoop():
        cuts = _linearStapleBreakPlan(
                    records, candidates,
                    continuousMinimum=continuousMinimum,
                    softMinimum=softMinimum,
                    preferredMinimum=preferredMinimum,
                    preferredMaximum=preferredMaximum,
                    targetLength=targetLength,
                    terminalMaximum=terminalMaximum,
                    preferDeletionDense=preferDeletionDense,
                    requireDeletionDenseMinimum=
                        requireDeletionDenseMinimum,
                    hardMaximum=hardMaximum,
                    densePreferredMinimum=densePreferredMinimum,
                    densePreferredMaximum=densePreferredMaximum,
                    denseTargetLength=denseTargetLength,
                    denseTerminalMaximum=denseTerminalMaximum)
        if cuts is None:
            return None
        return [(candidates[offset][0], candidates[offset][1])
                for offset in cuts]

    best = None
    # A circular oligo needs a closing nick as well as all internal cuts.
    # In every hard-valid partition, the first cut encountered after the
    # arbitrary record-0 origin is at most ``hardMaximum`` actual nucleotides
    # away.  It is therefore sufficient—and still exhaustive—to try only
    # those rotations.
    # Trying every legal boundary made a 2--3 kb Curved staple repeat the same
    # dynamic program hundreds of times.
    prefixActualLength = 0
    startCandidates = []
    for offset in range(1, len(records)):
        prefixActualLength += records[offset - 1][3]
        if prefixActualLength > hardMaximum:
            break
        if offset in candidates:
            startCandidates.append((offset, candidates[offset]))
    for startOffset, closing in startCandidates:
        rotated = records[startOffset:] + records[:startOffset]
        rotatedBoundaries = {}
        for offset, candidate in candidates.items():
            if offset == startOffset:
                continue
            rotatedOffset = (offset - startOffset) % len(records)
            if rotatedOffset:
                rotatedBoundaries[rotatedOffset] = candidate
        cuts = _linearStapleBreakPlan(
                    rotated, rotatedBoundaries, closing,
                    continuousMinimum=continuousMinimum,
                    softMinimum=softMinimum,
                    preferredMinimum=preferredMinimum,
                    preferredMaximum=preferredMaximum,
                    targetLength=targetLength,
                    terminalMaximum=terminalMaximum,
                    preferDeletionDense=preferDeletionDense,
                    requireDeletionDenseMinimum=
                        requireDeletionDenseMinimum,
                    hardMaximum=hardMaximum,
                    densePreferredMinimum=densePreferredMinimum,
                    densePreferredMaximum=densePreferredMaximum,
                    denseTargetLength=denseTargetLength,
                    denseTerminalMaximum=denseTerminalMaximum)
        if cuts is None:
            continue
        selectedOffsets = [startOffset] + [
                    (startOffset + offset) % len(records) for offset in cuts]
        selectedOffsets = sorted(selectedOffsets)
        selected = [candidates[offset] for offset in selectedOffsets]
        segmentCount = len(selectedOffsets)
        goodRuns = softMinimumLengths = preferredLengths = 0
        denseSegments = densePreferred = 0
        outsideDistance = deviation = 0
        for index, segmentStart in enumerate(selectedOffsets):
            segmentEnd = selectedOffsets[(index + 1) % segmentCount]
            if segmentEnd > segmentStart:
                segment = records[segmentStart:segmentEnd]
            else:
                segment = records[segmentStart:] + records[:segmentEnd]
            actualLength = sum(item[3] for item in segment)
            goodRuns += int(_hasContinuousPositionRun(
                                        segment, continuousMinimum))
            segmentDense = int(
                preferDeletionDense and
                sum(int(record[3]) == 0 for record in segment) >= 2)
            segmentSoftMinimum = (densePreferredMinimum
                                  if segmentDense else softMinimum)
            segmentPreferredMinimum = (densePreferredMinimum
                                       if segmentDense
                                       else preferredMinimum)
            segmentPreferredMaximum = (densePreferredMaximum
                                       if segmentDense
                                       else preferredMaximum)
            segmentTargetLength = (denseTargetLength
                                   if segmentDense else targetLength)
            softMinimumLengths += int(
                actualLength >= segmentSoftMinimum)
            preferredLengths += int(
                segmentPreferredMinimum <= actualLength <=
                segmentPreferredMaximum)
            denseSegments += segmentDense
            densePreferred += int(
                segmentDense and
                densePreferredMinimum <= actualLength <=
                densePreferredMaximum)
            outsideDistance += max(
                0, segmentPreferredMinimum - actualLength,
                actualLength - segmentPreferredMaximum)
            deviation += abs(actualLength - segmentTargetLength)
        preferredCuts = sum(int(item[2]) for item in selected)
        scaffoldSafeCuts = sum(int(item[3]) for item in selected)
        totalScaffoldClearance = sum(item[4] for item in selected)
        combinedPlacement = preferredCuts + scaffoldSafeCuts
        denseKey = (float(densePreferred) / denseSegments
                    if denseSegments else 1.0)
        key = (denseKey,
               float(goodRuns) / segmentCount,
               float(preferredLengths) / segmentCount,
               float(softMinimumLengths) / segmentCount,
               float(combinedPlacement) / (2 * segmentCount),
               float(scaffoldSafeCuts) / segmentCount,
               float(preferredCuts) / segmentCount,
               -segmentCount,
               float(totalScaffoldClearance) / segmentCount,
               -outsideDistance, -deviation,
               tuple((item[0], item[1]) for item in selected))
        if best is None or key > best[0]:
            best = (key, selected)
    if best is None:
        return None
    return [(item[0], item[1]) for item in best[1]]


def _occupiedRunBounds(strandSet, index):
    """Return the contiguous occupied interval containing ``index``."""
    intervals = sorted(strand.idxs() for strand in strandSet)
    merged = []
    for low, high in intervals:
        if merged and low <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], high)
        else:
            merged.append([low, high])
    for low, high in merged:
        if low <= index <= high:
            return low, high
    return None


def _existingStapleNickBoundaries(part):
    """Return native same-helix staple nicks present before Autobreak."""
    boundaries = set()
    for vh in part.getVirtualHelices():
        strands = sorted(vh.stapleStrandSet(),
                         key=lambda strand: strand.lowIdx())
        for left, right in zip(strands, strands[1:]):
            if left.highIdx() + 1 != right.lowIdx():
                continue
            if left.hasXoverAt(left.highIdx()) or \
                    right.hasXoverAt(right.lowIdx()):
                continue
            boundaries.add((vh.number(), right.lowIdx()))
    return boundaries


def _shortEdgeStapleXover(part, protectedNicks=(), excluded=()):
    """Find a staple xover 1--4 bases from a physical occupied edge."""
    records = []
    for vh in part.getVirtualHelices():
        for strand5p in vh.stapleStrandSet():
            strand3p = strand5p.connection3p()
            if strand3p is None or strand3p.part() is not part or \
                    not strand3p.strandSet().isStaple():
                continue
            records.append((strand5p.virtualHelix().number(),
                            strand5p.idx3Prime(),
                            strand3p.virtualHelix().number(),
                            strand3p.idx5Prime(), strand5p, strand3p))
    for unused_h5, unused_i5, unused_h3, unused_i3, strand5p, strand3p in \
            sorted(records, key=lambda item: item[:4]):
        if (id(strand5p), id(strand3p)) in excluded:
            continue
        sides = []
        qualifies = False
        for strand, index, isThreePrime in ((strand5p,
                                             strand5p.idx3Prime(), True),
                                            (strand3p,
                                             strand3p.idx5Prime(), False)):
            strandSet = strand.strandSet()
            bounds = _occupiedRunBounds(strandSet, index)
            if bounds is not None and \
                    min(index - bounds[0] + 1,
                        bounds[1] - index + 1) <= 4:
                qualifies = True
            if isThreePrime:
                delta = 1 if strandSet.isDrawn5to3() else -1
            else:
                delta = -1 if strandSet.isDrawn5to3() else 1
            nativeNeighbor = strandSet.getStrand(index + delta)
            if nativeNeighbor is None or nativeNeighbor is strand:
                sides = []
                break
            sides.append((strandSet, index, delta))
        sideBoundaries = set(
            (strandSet.virtualHelix().number(), max(index, index + delta))
            for strandSet, index, delta in sides)
        if sideBoundaries.intersection(protectedNicks):
            continue
        if qualifies and len(sides) == 2:
            return strand5p, strand3p, sides
    return None


def _removeShortEdgeStapleXovers(part, protectedNicks=()):
    """Remove generated staple xovers within 1--4 bases of an edge.

    Removing the crossover alone would leave two artificial nicks on the
    native helices.  Heal both same-helix joins as part of the same operation.
    The caller owns the surrounding undo macro.
    """
    removed = 0
    excluded = set()
    while True:
        edgeXover = _shortEdgeStapleXover(
                                part, protectedNicks, excluded)
        if edgeXover is None:
            break
        strand5p, strand3p, sides = edgeXover

        # AutoCS normally creates a reciprocal pair at one crossover site.
        # The native neighbor on each helix can therefore still be connected
        # to the other native neighbor.  Capture that partner before removing
        # the first crossover; otherwise strandsCanBeMerged() correctly says
        # the native joins are still occupied and the old cleanup stops after
        # only half-removing the site.
        nativeNeighbors = [
            strandSet.getStrand(index + delta)
            for strandSet, index, delta in sides]

        def boundaryConnection(strand, adjacent):
            if strand.highIdx() < adjacent.lowIdx():
                return strand.connectionHigh(), strand.highIdx()
            return strand.connectionLow(), strand.lowIdx()

        partner = None
        if len(nativeNeighbors) == 2 and all(nativeNeighbors):
            firstNeighbor, secondNeighbor = nativeNeighbors
            firstCurrent = sides[0][0].getStrand(sides[0][1])
            secondCurrent = sides[1][0].getStrand(sides[1][1])
            firstBoundary, firstBoundaryIndex = boundaryConnection(
                                        firstNeighbor, firstCurrent)
            secondBoundary, secondBoundaryIndex = boundaryConnection(
                                        secondNeighbor, secondCurrent)
            if firstBoundary is secondNeighbor and \
                    secondBoundary is firstNeighbor and \
                    firstNeighbor.idx3Prime() == firstBoundaryIndex and \
                    secondNeighbor.idx5Prime() == secondBoundaryIndex and \
                    firstNeighbor.connection3p() is secondNeighbor:
                partner = firstNeighbor, secondNeighbor
            elif firstBoundary is secondNeighbor and \
                    secondBoundary is firstNeighbor and \
                    secondNeighbor.idx3Prime() == secondBoundaryIndex and \
                    firstNeighbor.idx5Prime() == firstBoundaryIndex and \
                    secondNeighbor.connection3p() is firstNeighbor:
                partner = secondNeighbor, firstNeighbor

            # A native boundary occupied by some unrelated crossover cannot
            # be healed safely by this edge rule.  Skip it and continue
            # scanning instead of half-mutating it and aborting all cleanup.
            if partner is None and \
                    (firstBoundary is not None or
                     secondBoundary is not None):
                excluded.add((id(strand5p), id(strand3p)))
                continue

        part.removeXover(strand5p, strand3p)
        if partner is not None:
            part.removeXover(partner[0], partner[1])
        mergedBothSides = True
        for strandSet, index, delta in sides:
            first = strandSet.getStrand(index)
            second = strandSet.getStrand(index + delta)
            if first is None or second is None or first is second or \
                    not strandSet.strandsCanBeMerged(first, second):
                mergedBothSides = False
                break
            strandSet.mergeStrands(first, second)
        if not mergedBothSides:
            break
        removed += 1
    return removed


def _autoScaffoldCandidateFits(candidate, pairEvents, helixEvents,
                                densitySpacing, directionSpacing,
                                densityExemptPairs=()):
    """Return whether one scaffold crossover satisfies all hard spacings.

    ``directionSpacing`` is an index difference: 7 in Square represents an
    inclusive 8-bp segment, and 6 in Honeycomb represents 7 bp.  Both ends of
    every candidate are checked against all already selected directions.
    """
    fromHelix, toHelix, fromIndex, toIndex = candidate
    pairKey = (fromHelix, toHelix)
    if tuple(sorted((fromHelix, toHelix))) not in densityExemptPairs:
        for existingIndex in pairEvents.get(pairKey, ()):
            if abs(fromIndex - existingIndex) < densitySpacing:
                return False

    for helix, index, neighbor in ((fromHelix, fromIndex, toHelix),
                                   (toHelix, toIndex, fromHelix)):
        for existingIndex, existingNeighbor in helixEvents.get(helix, ()):
            if existingIndex == index:
                return False
            if existingNeighbor != neighbor and \
                    abs(index - existingIndex) < directionSpacing:
                return False
    return True


def _differentDirectionClearanceVector(helixEvents,
                                        noConflictDistance=0):
    """Return sorted nearest-opposite-direction distances for every event.

    Lexicographically maximizing this vector maximizes the worst clearance
    first, then the second-worst, and so on.  Events on helices that contain
    no other crossover direction receive ``noConflictDistance`` so they do
    not create an artificial zero-distance penalty.
    """
    clearances = []
    for unused_helix, events in helixEvents.items():
        for index, neighbor in events:
            distances = [
                abs(index - otherIndex)
                for otherIndex, otherNeighbor in events
                if otherNeighbor != neighbor]
            clearances.append(
                min(distances) if distances else noConflictDistance)
    return tuple(sorted(clearances))


def _autoScaffoldSnakePaths(helixRecords):
    """Return deterministic, non-branching paths through adjacent helices.

    ``helixRecords`` contains ``(number, row, column, neighbors)`` tuples.
    Rows are traversed in alternating directions, matching the conventional
    serpentine scaffold route.  Holes or disconnected helix groups start a
    new path rather than introducing a branch or closing a ring.
    """
    if not helixRecords:
        return []

    byNumber = dict((record[0], record) for record in helixRecords)
    rows = defaultdict(list)
    for number, row, column, unused_neighbors in helixRecords:
        rows[row].append((column, number))

    geometricOrder = []
    for rowIndex, row in enumerate(sorted(rows)):
        rowRecords = sorted(rows[row], reverse=bool(rowIndex % 2))
        geometricOrder.extend(number for unused_column, number in rowRecords)

    rank = dict((number, index) for index, number in
                enumerate(geometricOrder))

    adjacency = {}
    for number, unused_row, unused_column, neighbors in helixRecords:
        adjacency[number] = set(neighbor for neighbor in neighbors
                                if neighbor in byNumber)

    def extendPath(path, available):
        """Greedily extend one seeded route without branches."""
        while True:
            currentRank = rank[path[-1]]
            choices = [number for number in adjacency[path[-1]]
                       if number in available]
            if not choices:
                break

            # Prefer the next helix in serpentine order.  When a design has
            # holes, continue to the closest later helix before backtracking
            # in the ordering; this remains deterministic and branch-free.
            choices.sort(key=lambda number: (
                rank[number] <= currentRank,
                abs(rank[number] - currentRank),
                rank[number]))
            nextNumber = choices[0]
            path.append(nextNumber)
            available.remove(nextNumber)
        return path

    unvisited = set(byNumber)
    paths = []
    while unvisited:
        start = min(unvisited, key=lambda number: rank[number])
        path = [start]
        unvisited.remove(start)
        paths.append(extendPath(path, unvisited))
    return paths


def _connectAutoScaffoldSnakePaths(paths, helixRecords):
    """Join route fragments through physical neighbor edges.

    A stepped outline can have an articulation helix, so no Hamiltonian path
    through the helices exists without revisiting that articulation.  The old
    non-branching route consequently emitted one closed scaffold loop per
    fragment.  Represent each extra fragment as an out-and-back excursion
    from the first already-routed neighboring helix.  Repeated helix numbers
    are intentional: they reserve two crossover registers on the bridge and
    allow one global scaffold loop to traverse the branch.
    """
    connected = [list(path) for path in paths if path]
    if len(connected) < 2:
        return connected

    adjacency = dict(
        (number, set(neighbor for neighbor in neighbors))
        for number, unused_row, unused_column, neighbors in helixRecords)

    def closedExcursion(path, position):
        """Walk every edge of a path out and back to ``position``."""
        anchor = path[position]
        excursion = [anchor]
        excursion.extend(reversed(path[:position]))
        excursion.extend(path[1:position + 1])
        excursion.extend(path[position + 1:])
        excursion.extend(reversed(path[position:-1]))
        return excursion

    mergedRoutes = []
    while connected:
        route = connected.pop(0)
        while True:
            choices = []
            for fragmentIndex, fragment in enumerate(connected):
                for routeIndex, routeHelix in enumerate(route):
                    for fragmentPosition, fragmentHelix in \
                            enumerate(fragment):
                        if fragmentHelix in adjacency.get(routeHelix, ()):
                            choices.append((
                                routeIndex, fragmentIndex, fragmentPosition,
                                fragmentHelix))
            if not choices:
                break
            routeIndex, fragmentIndex, fragmentPosition, unused_helix = \
                min(choices)
            fragment = connected.pop(fragmentIndex)
            excursion = closedExcursion(fragment, fragmentPosition)
            route = (route[:routeIndex + 1] + excursion +
                     [route[routeIndex]] + route[routeIndex + 1:])
        mergedRoutes.append(route)
    return mergedRoutes


def _autoScaffoldStraightPaths(helixRecords):
    """Prefer long straight Square-lattice runs over zigzag traversal.

    Both row-wise and column-wise decompositions are evaluated.  The
    direction containing the longest continuous run wins; remaining runs are
    kept intact and later attached as branches.  This prevents a connected
    outline from repeatedly moving up/down and left/right merely to obtain a
    Hamiltonian helix order.
    """
    if not helixRecords:
        return []
    records = dict((number, (row, column, set(neighbors)))
                   for number, row, column, neighbors in helixRecords)

    def pathsForAxis(groupIndex, orderIndex):
        groups = defaultdict(list)
        for number, row, column, unused_neighbors in helixRecords:
            coordinates = (row, column)
            groups[coordinates[groupIndex]].append(
                (coordinates[orderIndex], number))
        paths = []
        for groupValue in sorted(groups):
            ordered = sorted(groups[groupValue])
            path = []
            for unused_position, number in ordered:
                if path and number not in records[path[-1]][2]:
                    paths.append(path)
                    path = []
                path.append(number)
            if path:
                paths.append(path)
        return paths

    rowPaths = pathsForAxis(0, 1)
    columnPaths = pathsForAxis(1, 0)
    choices = [rowPaths, columnPaths]

    def score(paths):
        lengths = sorted((len(path) for path in paths), reverse=True)
        return (lengths[0] if lengths else 0,
                sum(length * length for length in lengths),
                -len(paths))

    rowLengths = sorted((len(path) for path in rowPaths), reverse=True)
    hasWideTwoRowCap = (
        len(rowLengths) >= 2 and
        rowLengths[0] == rowLengths[1] and
        rowLengths[0] >= 4 and
        all(length * 2 <= rowLengths[0]
            for length in rowLengths[2:]))
    paths = rowPaths if hasWideTwoRowCap else max(choices, key=score)
    # The longest straight run is the main scaffold trunk.  Within a run,
    # increasing panel coordinates leave the later/rightmost helices at the
    # end, where the final pair can serve as the clean closing seam.
    return sorted(paths, key=lambda path: (
        -len(path),
        min(records[number][0] for number in path),
        min(records[number][1] for number in path),
        min(path)))


def _autoScaffoldRightPanelModulePaths(part, helixRecords):
    """Return explicit eight-helix Square modules in Path-view order.

    Large modular sheets can contain many geometric shortcuts between
    neighboring modules.  The generic longest-straight-run decomposition then
    mixes module members (Sp.json) and closes only the first resulting ring.
    When the user-maintained Path order cleanly consists of connected
    eight-helix blocks, preserve those blocks and connect them sequentially.
    Smaller/irregular designs keep the established geometric routing.
    """
    if part._step != 32 or not part._importedVHelixOrder:
        return None
    records = dict(
        (number, set(neighbors))
        for number, unused_row, unused_column, neighbors in helixRecords)
    ordered = []
    for coord in part._importedVHelixOrder:
        vh = part.virtualHelixAtCoord(coord)
        if vh is not None and vh.number() in records:
            ordered.append(vh.number())
    if len(ordered) < 32 or len(ordered) % 8:
        return None
    if len(set(ordered)) != len(records) or set(ordered) != set(records):
        return None
    paths = []
    for moduleIndex, index in enumerate(range(0, len(ordered), 8)):
        path = ordered[index:index + 8]
        phase = moduleIndex % 4
        if phase == 1:
            path = list(reversed(path))
        elif phase == 3:
            path = path[1:] + path[:1]
        paths.append(path)
    if any(
            second not in records.get(first, ())
            for path in paths
            for first, second in zip(path, path[1:])):
        return None
    if any(
            not any(second in records.get(first, ())
                    for first in firstPath for second in secondPath)
            for firstPath, secondPath in zip(paths, paths[1:])):
        return None
    return paths


def _autoScaffoldImportedOrder(part, helixRecords):
    """Return scaffold helix numbers in the user's current Path-view order."""
    available = set(record[0] for record in helixRecords)
    ordered = []
    for coord in part._importedVHelixOrder or ():
        vh = part.virtualHelixAtCoord(coord)
        if vh is not None and vh.number() in available and \
                vh.number() not in ordered:
            ordered.append(vh.number())
    ordered.extend(sorted(available.difference(ordered)))
    return ordered


def _autoScaffoldGeneralModulePaths(part, helixRecords):
    """Build general straight modules ordered by the Path panel.

    Unlike the established eight-helix Square special case, this candidate
    accepts arbitrary module sizes.  Geometry defines each uninterrupted
    straight module; the user's Path-view order determines module order and
    orientation.  It is returned as an alternate candidate rather than
    replacing the proven default route.
    """
    if not helixRecords:
        return None
    adjacency = dict(
        (number, set(neighbors))
        for number, unused_row, unused_column, neighbors in helixRecords)
    imported = _autoScaffoldImportedOrder(part, helixRecords)
    rank = dict((number, index) for index, number in enumerate(imported))
    geometricPaths = (
        _autoScaffoldStraightPaths(helixRecords)
        if part._step == 32 else
        _autoScaffoldSnakePaths(helixRecords))
    paths = []
    for path in geometricPaths:
        path = list(path)
        if len(path) > 1 and \
                rank.get(path[-1], len(rank)) < \
                rank.get(path[0], len(rank)):
            path.reverse()
        paths.append(path)
    paths.sort(key=lambda path: (
        min(rank.get(number, len(rank)) for number in path),
        -len(path), tuple(path)))
    if len(paths) < 2:
        return None
    if any(
            second not in adjacency.get(first, ())
            for path in paths
            for first, second in zip(path, path[1:])):
        return None
    # Every consecutive module must have at least one physical bridge.
    if any(
            not any(second in adjacency.get(first, ())
                    for first in firstPath for second in secondPath)
            for firstPath, secondPath in zip(paths, paths[1:])):
        return None
    return paths


def _autoScaffoldImportedRunPaths(part, helixRecords):
    """Split the Path-view order only where consecutive helices are not adjacent."""
    if not helixRecords:
        return None
    adjacency = dict(
        (number, set(neighbors))
        for number, unused_row, unused_column, neighbors in helixRecords)
    ordered = _autoScaffoldImportedOrder(part, helixRecords)
    paths = []
    for number in ordered:
        if not paths or number not in adjacency.get(paths[-1][-1], ()):
            paths.append([number])
        else:
            paths[-1].append(number)
    if len(paths) < 2 or any(len(path) < 2 for path in paths):
        return None
    if any(
            not any(second in adjacency.get(first, ())
                    for first in firstPath for second in secondPath)
            for firstPath, secondPath in zip(paths, paths[1:])):
        return None
    return paths


def _autoScaffoldModuleOrderPenalty(part, paths, helixRecords):
    """Return Path-order fragmentation/inversion penalties for route modules."""
    imported = _autoScaffoldImportedOrder(part, helixRecords)
    rank = dict((number, index) for index, number in enumerate(imported))
    spans = []
    fragmentation = 0
    orientation = 0
    for path in paths:
        ranks = [rank.get(number, len(rank)) for number in path]
        if not ranks:
            continue
        spans.append(min(ranks))
        fragmentation += max(ranks) - min(ranks) + 1 - len(set(ranks))
        orientation += sum(
            1 for left, right in zip(ranks, ranks[1:])
            if right < left)
    inversions = sum(
        1 for index, left in enumerate(spans)
        for right in spans[index + 1:] if right < left)
    return fragmentation, inversions, orientation


def _autoScaffoldSequentialModuleBridgePairs(paths, helixRecords):
    """Allow every physical bridge between consecutive ordered modules."""
    adjacency = dict(
        (number, set(neighbors))
        for number, unused_row, unused_column, neighbors in helixRecords)
    return set(
        tuple(sorted((first, second)))
        for firstPath, secondPath in zip(paths, paths[1:])
        for first in firstPath for second in secondPath
        if second in adjacency.get(first, ()))


def _autoScaffoldStraightBridgePairs(paths, helixRecords):
    """Return symmetric outer bridges between neighboring straight runs.

    The generic fragment connector represents every shorter run as an
    out-and-back excursion attached at the first available neighbor.  On a
    trapezoidal Square outline that biases every attachment to the same side
    and later requires irregular interior loop merges.  Straight runs instead
    have a natural pair of outer connections: use the first and last physical
    adjacency between each neighboring run.
    """
    coordinates = dict(
        (number, (row, column))
        for number, row, column, unused_neighbors in helixRecords)
    adjacency = dict(
        (number, set(neighbors))
        for number, unused_row, unused_column, neighbors in helixRecords)
    bridges = set()
    for firstIndex, firstPath in enumerate(paths):
        for secondPath in paths[firstIndex + 1:]:
            candidates = set(
                tuple(sorted((first, second)))
                for first in firstPath for second in secondPath
                if second in adjacency.get(first, ()))
            if not candidates:
                continue
            firstStart = coordinates[firstPath[0]]
            firstEnd = coordinates[firstPath[-1]]
            horizontal = (abs(firstEnd[1] - firstStart[1]) >=
                          abs(firstEnd[0] - firstStart[0]))
            axis = 1 if horizontal else 0
            ordered = sorted(candidates, key=lambda pair: (
                (coordinates[pair[0]][axis] +
                 coordinates[pair[1]][axis]),
                pair))
            bridges.add(ordered[0])
            bridges.add(ordered[-1])
    return bridges


def _autoScaffoldSingleStraightBridgePairs(paths, helixRecords):
    """Choose one early physical bridge between consecutive straight runs."""
    adjacency = dict(
        (number, set(neighbors))
        for number, unused_row, unused_column, neighbors in helixRecords)
    positions = dict(
        (number, (pathIndex, position))
        for pathIndex, path in enumerate(paths)
        for position, number in enumerate(path))
    bridges = set()
    for firstPath, secondPath in zip(paths, paths[1:]):
        candidates = [
            tuple(sorted((first, second)))
            for first in firstPath for second in secondPath
            if second in adjacency.get(first, ())]
        if not candidates:
            continue
        bridges.add(min(candidates, key=lambda pair: (
            positions[pair[0]][1] + positions[pair[1]][1],
            max(positions[pair[0]][1], positions[pair[1]][1]),
            pair)))
    return bridges


def _autoScaffoldPhaseFamilies(records, step):
    """Group reciprocal native crossover phases belonging to one register."""
    remaining = set(record[2] % step for record in records)
    families = []
    while remaining:
        phase = min(remaining)
        remaining.remove(phase)
        family = set([phase])
        pending = [phase]
        while pending:
            current = pending.pop()
            adjacent = [other for other in remaining
                        if min((current - other) % step,
                               (other - current) % step) <= 1]
            for other in adjacent:
                remaining.remove(other)
                family.add(other)
                pending.append(other)
        families.append(family)
    return families


def _autoScaffoldPhaseDistance(first, second, step):
    return min(min((left - right) % step, (right - left) % step)
               for left in first for right in second)


def _filterAutoScaffoldSparseCandidatesForPaths(candidates, paths, step,
                                           avoidMultiplesOfEight=False,
                                           minBase=None, maxBase=None,
                                           helixCoordinates=None,
                                           alignedPaths=(),
                                           seedUnevenFinalSeam=False,
                                           legacyEdgeTails=False,
                                           densitySpacing=None,
                                           preferBoundaryRegisters=False,
                                           registerBoundaryTrim=0,
                                           seamPairs=None,
                                           reversePhaseOrder=False,
                                           phasePreferenceOffset=None):
    """Keep one well-spaced crossover phase family per snake-path edge.

    At a physical edge, omit a crossover when it would close a short terminal
    section.  The final helix pair is the global routing seam, so it keeps one
    boundary block; the normal longitudinal backbone at later legal blocks
    joins the otherwise separate period-sized rings into one closed scaffold.
    Scaffold bases are never deleted.
    """
    recordsByPair = defaultdict(list)
    for record in candidates:
        pair = tuple(sorted((record[0], record[1])))
        recordsByPair[pair].append(record)

    selected = []
    for pathIndex, path in enumerate(paths):
        previousFamily = None
        previousRecords = None
        alignmentReferences = {}
        for edgeIndex, (first, second) in enumerate(zip(path, path[1:])):
            pair = tuple(sorted((first, second)))
            records = recordsByPair.get(pair, ())
            if not records:
                previousFamily = None
                previousRecords = None
                continue
            families = _autoScaffoldPhaseFamilies(records, step)

            familyPool = families
            if previousFamily is not None:
                hardPhaseSpacing = 6 if step == 21 else 7
                legalFamilies = [
                    family for family in families
                    if _autoScaffoldPhaseDistance(
                        family, previousFamily, step) >= hardPhaseSpacing]
                if legalFamilies:
                    familyPool = legalFamilies

            def familyKey(family):
                avoided = sum(1 for phase in family
                              if avoidMultiplesOfEight and phase % 8 == 0)
                if previousFamily is None:
                    # Honeycomb uses the second native scaffold register as
                    # its route seed.  Subsequent registers are then chosen
                    # by maximum separation on the shared helix.
                    if step == 21:
                        return (0, avoided, -min(family))
                    if phasePreferenceOffset is not None:
                        phaseDistance = min(
                            min((phase - phasePreferenceOffset) % step,
                                (phasePreferenceOffset - phase) % step)
                            for phase in family)
                        return (phaseDistance, avoided, min(family))
                    return (0, avoided,
                            -max(family) if reversePhaseOrder else
                            min(family))
                spacing = _autoScaffoldPhaseDistance(
                                            family, previousFamily, step)
                phaseTie = (min(
                    min((phase - phasePreferenceOffset) % step,
                        (phasePreferenceOffset - phase) % step)
                    for phase in family)
                    if phasePreferenceOffset is not None else 0)
                if step == 32 and avoidMultiplesOfEight:
                    # Hard inclusive 8-bp spacing was enforced in familyPool.
                    # Among those legal families, prefer non-x8 before
                    # maximizing any additional different-direction clearance.
                    return (avoided, -spacing, phaseTie,
                            -max(family) if reversePhaseOrder else
                            min(family))
                return (-spacing, avoided, phaseTie,
                        -max(family) if reversePhaseOrder else
                        min(family))

            chosenFamily = None
            chosenRecords = None
            if preferBoundaryRegisters and legacyEdgeTails and \
                    step == 21 and densitySpacing and \
                    densitySpacing > step:
                # A sparse honeycomb phase family contains several native
                # 21-bp blocks.  Choose its density register here, together
                # with the family, instead of leaving that decision to a
                # later record-order tie break.  This reproduces the compact
                # 1/42 route in ho-ai3 without rebuilding many alternatives.
                registerCount = max(1, densitySpacing // step)
                registerChoices = []
                allLow = min(record[2] for record in records)
                allHigh = max(record[2] for record in records)
                targetLow = min(allHigh,
                                allLow + int(registerBoundaryTrim))
                targetHigh = max(allLow,
                                 allHigh - int(registerBoundaryTrim))
                firstRow, firstColumn = helixCoordinates[first]
                secondRow, secondColumn = helixCoordinates[second]
                directionKey = (secondRow - firstRow,
                                secondColumn - firstColumn)
                alignmentKey = (edgeIndex % 2, directionKey,
                                allLow, allHigh)
                alignmentReference = alignmentReferences.get(alignmentKey)
                for family in families:
                    familyRecords = sorted(
                        [record for record in records
                         if record[2] % step in family],
                        key=lambda record: (record[2],
                                            record[0], record[1]))
                    blocks = []
                    for record in familyRecords:
                        if blocks and \
                                record[2] - blocks[-1][-1][2] <= 1:
                            blocks[-1].append(record)
                        else:
                            blocks.append([record])
                    for offset in range(min(registerCount, len(blocks))):
                        registerRecords = [
                            record for block in blocks[offset::registerCount]
                            for record in block]
                        if not registerRecords:
                            continue
                        low = min(record[2] for record in registerRecords)
                        high = max(record[2] for record in registerRecords)
                        boundaryPenalty = (abs(low - targetLow) +
                                           abs(targetHigh - high))
                        if previousRecords:
                            separation = min(
                                abs(record[2] - previousRecord[2])
                                for record in registerRecords
                                for previousRecord in previousRecords)
                        else:
                            separation = 0
                        registerChoices.append((
                            -separation, boundaryPenalty,
                            min(family), offset,
                            family, registerRecords))
                if registerChoices:
                    if alignmentReference:
                        referenceLow = min(
                            record[2] for record in alignmentReference)
                        referenceHigh = max(
                            record[2] for record in alignmentReference)
                        chosenChoice = min(
                            registerChoices,
                            key=lambda choice: (
                                abs(min(record[2]
                                        for record in choice[5]) -
                                    referenceLow) +
                                abs(max(record[2]
                                        for record in choice[5]) -
                                    referenceHigh),
                                choice[1], choice[0],
                                choice[2], choice[3]))
                    elif edgeIndex % 2 == 0:
                        chosenChoice = min(
                            registerChoices,
                            key=lambda choice: (
                                choice[1], choice[0],
                                min(record[2]
                                    for record in choice[5]),
                                choice[2], choice[3]))
                    else:
                        maximumSeparation = max(
                            -choice[0] for choice in registerChoices)
                        nearMaximum = [
                            choice for choice in registerChoices
                            if -choice[0] >= maximumSeparation - 3]
                        movesRight = (firstRow == secondRow and
                                      secondColumn > firstColumn)
                        movesLeft = (firstRow == secondRow and
                                     secondColumn < firstColumn)
                        if movesRight:
                            chosenChoice = min(
                                nearMaximum,
                                key=lambda choice: (
                                    min(record[2]
                                        for record in choice[5]),
                                    choice[1], choice[0],
                                    choice[2], choice[3]))
                        elif movesLeft:
                            chosenChoice = min(
                                nearMaximum,
                                key=lambda choice: (
                                    -min(record[2]
                                         for record in choice[5]),
                                    choice[1], choice[0],
                                    choice[2], choice[3]))
                        else:
                            chosenChoice = min(
                                nearMaximum,
                                key=lambda choice: choice[:4])
                    unused_separation, unused_boundary, unused_phase, \
                        unused_offset, chosenFamily, chosenRecords = \
                            chosenChoice
                    alignmentReferences.setdefault(
                                        alignmentKey, chosenRecords)
            if chosenRecords is None:
                chosenFamily = min(familyPool, key=familyKey)
                chosenRecords = [record for record in records
                                 if record[2] % step in chosenFamily]
            if step == 21 and minBase is not None and maxBase is not None \
                    and helixCoordinates is not None:
                firstRow, firstColumn = helixCoordinates[first]
                secondRow, secondColumn = helixCoordinates[second]
                sameRow = firstRow == secondRow
                movesRight = sameRow and secondColumn > firstColumn
                movesLeft = sameRow and secondColumn < firstColumn
                edgeMargin = 7
                chosenRecords = [record for record in chosenRecords
                                 if not (record[2] - minBase < edgeMargin and
                                         not movesLeft) and
                                    not (maxBase - record[2] < edgeMargin and
                                         not movesRight)]

            isClosingSeam = (
                tuple(sorted((first, second))) in seamPairs
                if seamPairs is not None else
                second == path[-1])
            if isClosingSeam and \
                    len(chosenRecords) > 2:
                orderedRecords = sorted(chosenRecords,
                                        key=lambda record: record[2])
                blocks = []
                for record in orderedRecords:
                    if blocks and record[2] - blocks[-1][-1][2] <= 1:
                        blocks[-1].append(record)
                    else:
                        blocks.append([record])
                if len(blocks) > 2:
                    if legacyEdgeTails:
                        # Original manual-design behavior: retain the two
                        # legal boundary blocks and leave unmatched scaffold
                        # beyond them as single-stranded edge tails.
                        chosenRecords = blocks[0] + blocks[-1]
                    elif tuple(path) in alignedPaths or \
                            not seedUnevenFinalSeam:
                        # When every helix already shares the same two legal
                        # perimeter endpoints, this pair is a pure seam and
                        # needs no internal crossover block (sq-ai2).
                        chosenRecords = []
                    else:
                        # Uneven outlines need one seed block so the endpoint
                        # stage can absorb their staggered edge fragments.
                        chosenRecords = blocks[0]
            selected.extend(chosenRecords)
            previousFamily = chosenFamily
            previousRecords = chosenRecords
    return selected


def _selectHoneycombSparseLoopRecords(candidates, paths, spacing,
                                      helixCoordinates,
                                      directionSpacing=6,
                                      drawn5to3=None):
    """Select a periodic locally closed scaffold route.

    Alternating path gaps play different topological roles.  Even gaps keep
    outward-facing single crossovers at the two perimeter registers and
    reciprocal pairs between them.  Odd gaps keep reciprocal internal blocks.
    The final even gap follows the same rule: dropping its internal reciprocal
    blocks opens one of the otherwise closed lattice-period units.
    """
    recordsByPair = defaultdict(list)
    for record in candidates:
        recordsByPair[tuple(sorted(record[:2]))].append(record)

    selected = []
    previousIndices = []
    for path in paths:
        for edgeIndex, (first, second) in enumerate(zip(path, path[1:])):
            records = recordsByPair.get(tuple(sorted((first, second))), ())
            blocks = []
            used = set()
            for forward in records:
                if forward in used:
                    continue
                partners = [
                    reverse for reverse in records
                    if reverse not in used and
                    reverse[0] == forward[1] and
                    reverse[1] == forward[0] and
                    abs(reverse[2] - forward[2]) == 1]
                if not partners:
                    continue
                reverse = min(partners,
                              key=lambda item: (item[2], item[0], item[1]))
                block = tuple(sorted((forward, reverse),
                                     key=lambda item: item[2]))
                blocks.append(block)
                used.update(block)
            blocks = sorted(set(blocks),
                            key=lambda block: block[0][2])
            if not blocks:
                previousIndices = []
                continue

            byLow = dict((block[0][2], block) for block in blocks)
            chains = []
            for block in blocks:
                low = block[0][2]
                chain = []
                while low in byLow:
                    chain.append(byLow[low])
                    low += spacing
                if chain:
                    chains.append(chain)
            maximumBlocks = max(len(chain) for chain in chains)
            chains = [chain for chain in chains
                      if len(chain) == maximumBlocks]

            def chainSeparation(chain):
                indices = [record[2] for block in chain
                           for record in block]
                return (
                    min(abs(index - previous)
                        for index in indices for previous in previousIndices)
                    if previousIndices else 0)
            maximumSeparation = max(
                chainSeparation(chain) for chain in chains)
            firstCoord = helixCoordinates.get(first)
            secondCoord = helixCoordinates.get(second)
            movesLeft = (
                firstCoord is not None and secondCoord is not None and
                firstCoord[0] == secondCoord[0] and
                secondCoord[1] < firstCoord[1])
            if edgeIndex % 2 == 1 and movesLeft:
                legal = [chain for chain in chains
                         if chainSeparation(chain) >= directionSpacing]
                chain = min(legal or chains,
                            key=lambda item: item[0][0][2])
            else:
                nearMaximum = [
                    chain for chain in chains
                    if chainSeparation(chain) >= maximumSeparation - 3]
                chain = min(nearMaximum,
                            key=lambda item: item[0][0][2])
            isEvenGap = edgeIndex % 2 == 0
            edgeRecords = []
            if isEvenGap:
                lowBlock = chain[0]
                highBlock = chain[-1]
                firstRunsForward = bool(
                    drawn5to3 and drawn5to3.get(first))
                lowDirection = (
                    (second, first) if firstRunsForward
                    else (first, second))
                highDirection = tuple(reversed(lowDirection))
                lowRecords = [
                    record for record in lowBlock
                    if record[:2] == lowDirection]
                highRecords = [
                    record for record in highBlock
                    if record[:2] == highDirection]
                if lowRecords:
                    edgeRecords.append(lowRecords[0])
                for block in chain[1:-1]:
                    edgeRecords.extend(block)
                if highRecords and highBlock is not lowBlock:
                    edgeRecords.append(highRecords[0])
            else:
                # The last legal block is reserved for the adjacent perimeter
                # gap, leaving a clean and aligned right edge.
                for block in chain[:-1]:
                    edgeRecords.extend(block)
            selected.extend(edgeRecords)
            previousIndices = [record[2] for record in edgeRecords]
    return sorted(set(selected),
                  key=lambda item: (item[2], item[0], item[1]))


def _filterAutoScaffoldCandidatesForPaths(candidates, paths, step,
                                           avoidMultiplesOfEight=False,
                                           minBase=None, maxBase=None,
                                           helixCoordinates=None,
                                           legacyEdgeTails=False,
                                           reversePhaseOrder=False,
                                           densitySpacing=None):
    """Keep one well-spaced crossover phase family per snake-path edge.

    At a physical edge, omit a crossover when it would close a short terminal
    section.  Every path gap, including the final helix pair, retains its
    periodic internal blocks so each module is composed of independent
    lattice-period scaffold loops; there is no special seam.  Scaffold bases
    are never resized or deleted.
    """
    recordsByPair = defaultdict(list)
    for record in candidates:
        pair = tuple(sorted((record[0], record[1])))
        recordsByPair[pair].append(record)

    # At reduced crossover densities the register repeats over the requested
    # density period, not over one native lattice period.  Comparing phases
    # modulo ``step`` loses that register information: for example, in a
    # 1/42 Honeycomb design the {11, 12} family in the following 21-bp block
    # can be farther from {15, 16} than the apparently farther native phase
    # {1, 2}.  Keep the full density-period phase while optimizing clearance.
    registerStep = max(step, int(densitySpacing or step))

    selected = []
    for path in paths:
        previousFamily = None
        for edgeIndex, (first, second) in enumerate(zip(path, path[1:])):
            pair = tuple(sorted((first, second)))
            records = recordsByPair.get(pair, ())
            if not records:
                previousFamily = None
                continue
            families = _autoScaffoldPhaseFamilies(records, registerStep)

            familyPool = families
            if previousFamily is not None:
                hardPhaseSpacing = 6 if step == 21 else 7
                legalFamilies = [
                    family for family in families
                    if _autoScaffoldPhaseDistance(
                        family, previousFamily,
                        registerStep) >= hardPhaseSpacing]
                if legalFamilies:
                    familyPool = legalFamilies

            def familyKey(family):
                avoidModulo = (7 if step == 21 else
                               8 if avoidMultiplesOfEight else None)
                avoided = sum(1 for phase in family
                              if avoidModulo is not None and
                              phase % avoidModulo == 0)
                if previousFamily is None:
                    # Honeycomb uses the second native scaffold register as
                    # its route seed.  Subsequent registers are then chosen
                    # by maximum separation on the shared helix.
                    if step == 21:
                        nativePhases = [phase % step for phase in family]
                        return (0, avoided, -min(nativePhases), min(family))
                    nativePhases = [phase % step for phase in family]
                    nativeOrder = (-max(nativePhases)
                                   if reversePhaseOrder else
                                   min(nativePhases))
                    registerOrder = (-max(family)
                                     if reversePhaseOrder else
                                     min(family))
                    return (0, avoided, nativeOrder, registerOrder)
                spacing = _autoScaffoldPhaseDistance(
                    family, previousFamily, registerStep)
                phaseOrder = (-max(family) if reversePhaseOrder else
                              min(family))
                if step == 32 and avoidMultiplesOfEight:
                    # Once the inclusive 8-bp hard spacing is satisfied,
                    # avoiding x8 phases outranks maximizing extra distance.
                    return (avoided, -spacing, phaseOrder)
                return (-spacing, avoided, phaseOrder)

            chosenFamily = min(familyPool, key=familyKey)
            chosenRecords = [record for record in records
                             if record[2] % registerStep in chosenFamily]
            if step == 21 and minBase is not None and maxBase is not None \
                    and helixCoordinates is not None:
                firstRow, firstColumn = helixCoordinates[first]
                secondRow, secondColumn = helixCoordinates[second]
                sameRow = firstRow == secondRow
                movesRight = sameRow and secondColumn > firstColumn
                movesLeft = sameRow and secondColumn < firstColumn
                edgeMargin = 7
                chosenRecords = [record for record in chosenRecords
                                 if not (record[2] - minBase < edgeMargin and
                                         not movesLeft) and
                                    not (maxBase - record[2] < edgeMargin and
                                         not movesRight)]

            selected.extend(chosenRecords)
            previousFamily = chosenFamily
    return selected


def _selectAutoScaffoldCrossoverRecords(candidates, existing,
                                         densitySpacing, directionSpacing,
                                         avoidMultiplesOfEight=False,
                                         variant=0,
                                         densityExemptPairs=(),
                                         helixOrder=None):
    """Choose a dense deterministic set under scaffold crossover rules.

    Several phase and traversal orders are evaluated.  Cardinality is the
    primary score, followed by avoiding square-lattice multiples of eight,
    covering more directed helix pairs, and keeping same-direction spacings
    close to the lattice period.
    """
    candidates = sorted(set(candidates))
    existing = list(existing)
    if not candidates:
        return []
    helixOrder = helixOrder or {}
    fallbackRank = len(helixOrder)

    def helixRank(number):
        return helixOrder.get(number, fallbackRank + number)

    def directedPairKey(item):
        return (helixRank(item[0]), helixRank(item[1]))

    def geometricRecordKey(item):
        return (item[2], helixRank(item[0]), helixRank(item[1]),
                item[3])

    def isAvoided(candidate):
        return bool(avoidMultiplesOfEight and candidate[2] % 8 == 0)

    def phaseDistance(index, phase):
        delta = (index - phase) % densitySpacing
        return min(delta, densitySpacing - delta)

    orders = [
        sorted(candidates,
               key=lambda item: (directedPairKey(item),
                                 isAvoided(item), item[2])),
        sorted(candidates,
               key=lambda item: (directedPairKey(item),
                                 isAvoided(item), -item[2])),
        sorted(candidates,
               key=lambda item: (isAvoided(item), item[2],
                                 directedPairKey(item))),
        sorted(candidates,
               key=lambda item: (isAvoided(item), -item[2],
                                 directedPairKey(item))),
    ]
    for phase in range(densitySpacing):
        orders.append(sorted(
            candidates,
            key=lambda item: (directedPairKey(item),
                              isAvoided(item),
                              phaseDistance(item[2], phase), item[2])))
        orders.append(sorted(
            candidates,
            key=lambda item: (isAvoided(item),
                              phaseDistance(item[2], phase), item[2],
                              directedPairKey(item))))
        orders.append(sorted(
            candidates,
            key=lambda item: (isAvoided(item),
                              phaseDistance(item[2], phase), -item[2],
                              directedPairKey(item))))

    alternatives = {}
    for orderedCandidates in orders:
        pairEvents = defaultdict(list)
        helixEvents = defaultdict(list)
        for fromHelix, toHelix, fromIndex, toIndex in existing:
            pairEvents[(fromHelix, toHelix)].append(fromIndex)
            helixEvents[fromHelix].append((fromIndex, toHelix))
            helixEvents[toHelix].append((toIndex, fromHelix))

        selected = []
        for candidate in orderedCandidates:
            if not _autoScaffoldCandidateFits(
                    candidate, pairEvents, helixEvents,
                    densitySpacing, directionSpacing,
                    densityExemptPairs):
                continue
            fromHelix, toHelix, fromIndex, toIndex = candidate
            selected.append(candidate)
            pairEvents[(fromHelix, toHelix)].append(fromIndex)
            helixEvents[fromHelix].append((fromIndex, toHelix))
            helixEvents[toHelix].append((toIndex, fromHelix))

        avoidedCount = sum(1 for item in selected if isAvoided(item))
        selectedPairs = set((item[0], item[1]) for item in selected)
        spacingPenalty = 0
        for pairKey, indices in pairEvents.items():
            if tuple(sorted(pairKey)) in densityExemptPairs:
                continue
            orderedIndices = sorted(indices)
            spacingPenalty += sum(
                max(0, right - left - densitySpacing)
                for left, right in zip(orderedIndices, orderedIndices[1:]))
        directionClearance = _differentDirectionClearanceVector(
            helixEvents, noConflictDistance=densitySpacing * 2)
        score = (len(selected), -avoidedCount, directionClearance,
                 len(selectedPairs), -spacingPenalty)
        normalized = tuple(sorted(selected, key=geometricRecordKey))
        previous = alternatives.get(normalized)
        if previous is None or score > previous:
            alternatives[normalized] = score

    ranked = sorted(
        alternatives.items(),
        key=lambda item: (
            item[1],
            tuple(geometricRecordKey(record) for record in item[0])),
        reverse=True)
    if not ranked:
        return []
    variant = min(max(0, int(variant)), len(ranked) - 1)
    return list(ranked[variant][0])


def _selectMinimumDensityScaffoldRecords(candidates, paths, originalSpans,
                                         step, drawn5to3, existing=()):
    """Select one connection group per consecutive helix gap.

    Even-numbered gaps use one crossover at each edge.  Odd-numbered gaps use
    one reciprocal adjacent pair near the overlap center.  This alternating
    pattern is the explicit minimum-density scaffold route.
    """
    byPair = defaultdict(list)
    for record in candidates:
        byPair[tuple(sorted(record[:2]))].append(record)
    occupiedPairs = set(tuple(sorted(record[:2])) for record in existing)
    avoidModulo = 7 if step == 21 else 8

    def avoided(record):
        return int(record[2] % avoidModulo == 0)

    selected = []
    for path in paths:
        for edgeIndex, (first, second) in enumerate(
                zip(path, path[1:])):
            pair = tuple(sorted((first, second)))
            if pair in occupiedPairs:
                continue
            records = byPair.get(pair, ())
            forward = [record for record in records
                       if record[0] == first and record[1] == second]
            reverse = [record for record in records
                       if record[0] == second and record[1] == first]
            firstSpan = originalSpans.get(first)
            secondSpan = originalSpans.get(second)
            if not forward or not reverse or \
                    firstSpan is None or secondSpan is None:
                continue
            overlapLow = max(firstSpan[0], secondSpan[0])
            overlapHigh = min(firstSpan[1], secondSpan[1])
            if overlapLow > overlapHigh:
                continue
            if edgeIndex % 2 == 0:
                # Direction is structural, not a distance tie-breaker:
                # it follows the 5'→3' direction of the first helix.
                forwardAtHigh = drawn5to3.get(first, True)
                highRecords = forward if forwardAtHigh else reverse
                lowRecords = reverse if forwardAtHigh else forward
                high = min(
                    highRecords,
                    key=lambda record:
                    (avoided(record),
                     abs(overlapHigh - record[2]), -record[2]))
                low = min(
                    lowRecords,
                    key=lambda record:
                    (avoided(record),
                     abs(record[2] - overlapLow), record[2]))
                if high[2] != low[2]:
                    selected.extend((low, high))
            else:
                reciprocalPairs = [
                    (forwardRecord, reverseRecord)
                    for forwardRecord in forward
                    for reverseRecord in reverse
                    if abs(forwardRecord[2] - reverseRecord[2]) == 1 and
                    ((drawn5to3.get(first, True) and
                      forwardRecord[2] < reverseRecord[2]) or
                     (not drawn5to3.get(first, True) and
                      forwardRecord[2] > reverseRecord[2]))]
                if reciprocalPairs:
                    center = (overlapLow + overlapHigh) / 2.0
                    chosen = min(
                        reciprocalPairs,
                        key=lambda pairRecords:
                        (sum(1 for record in pairRecords
                             if record[2] % avoidModulo == 0),
                         abs(((pairRecords[0][2] +
                               pairRecords[1][2]) / 2.0) - center),
                         min(pairRecords[0][2], pairRecords[1][2])))
                    selected.extend(chosen)
    return sorted(set(selected),
                  key=lambda record: (record[2], record[0], record[1]))


def _densifyClosedScaffoldLoop(part, candidates, spacing,
                                directionSpacing,
                                avoidMultiplesOfEight=False,
                                densityExemptPairs=()):
    """Add reciprocal crossover pairs without changing a closed topology.

    Once a single scaffold loop exists, either member of a reciprocal
    crossover pair temporarily splits that loop.  The second member joins it
    again, so stopping immediately when the first loop is found leaves an
    unnecessarily sparse design.  Select additions against the *final*
    perimeter crossovers, then install only complete adjacent reciprocal
    pairs.  This also prevents a dense pair from being placed too close to an
    endpoint-closing crossover.
    """
    existing = part._existingScaffoldCrossoverRecords()
    additions = _selectAutoScaffoldCrossoverRecords(
        candidates, existing, spacing, directionSpacing,
        avoidMultiplesOfEight, 0, densityExemptPairs)
    additions = set(additions)
    paired = []
    used = set()
    for record in sorted(additions, key=lambda item:
                         (item[2], item[0], item[1])):
        if record in used:
            continue
        reverse = [
            other for other in additions
            if other not in used and
            other[0] == record[1] and other[1] == record[0] and
            abs(other[2] - record[2]) == 1]
        if not reverse:
            continue
        partner = min(reverse, key=lambda item:
                      (abs(item[2] - record[2]), item[2]))
        paired.append(tuple(sorted((record, partner),
                                   key=lambda item: item[2])))
        used.add(record)
        used.add(partner)

    helixByNumber = dict((vh.number(), vh)
                         for vh in part.getVirtualHelices())
    loopOligos = [oligo for oligo in _scaffoldOligos(part)
                  if oligo.isLoop()]
    if not loopOligos:
        return 0
    targetLoop = max(loopOligos, key=lambda oligo: oligo.length())
    targetLength = targetLoop.length()
    eligiblePairs = []
    for pair in paired:
        pairStrands = []
        for fromNumber, toNumber, idx5p, idx3p in pair:
            fromHelix = helixByNumber.get(fromNumber)
            toHelix = helixByNumber.get(toNumber)
            strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                        if fromHelix is not None else None)
            strand3p = (toHelix.scaffoldStrandSet().getStrand(idx3p)
                        if toHelix is not None else None)
            pairStrands.append((strand5p, strand3p))
        # Regular designs may retain ordinary open scaffold tails.  Never
        # densify through one of those tails: doing so can cut a short ring
        # out of the otherwise valid global loop.
        if any(strand5p is None or strand3p is None or
               strand5p.oligo() is not targetLoop or
               strand3p.oligo() is not targetLoop
               for strand5p, strand3p in pairStrands):
            continue
        eligiblePairs.append(pair)

    def installPairs(pairs):
        installed = []
        records = [record for pair in pairs for record in pair]
        for fromNumber, toNumber, idx5p, idx3p in records:
            fromHelix = helixByNumber.get(fromNumber)
            toHelix = helixByNumber.get(toNumber)
            strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                        if fromHelix is not None else None)
            strand3p = (toHelix.scaffoldStrandSet().getStrand(idx3p)
                        if toHelix is not None else None)
            if not part._canCreateScaffoldXover(
                    strand5p, strand3p, idx5p):
                break
            part.createXover(strand5p, idx5p, strand3p, idx3p)
            currentStrand5p = \
                fromHelix.scaffoldStrandSet().getStrand(idx5p)
            currentStrand3p = (currentStrand5p.connection3p()
                               if currentStrand5p is not None else None)
            if currentStrand3p is None or \
                    currentStrand3p.virtualHelix().number() != toNumber or \
                    currentStrand3p.idx5Prime() != idx3p:
                break
            installed.append((fromNumber, toNumber, idx5p, idx3p))
        accepted = (len(installed) == len(records) and
                    any(oligo.isLoop() and
                        oligo.length() == targetLength
                        for oligo in _scaffoldOligos(part)))
        if accepted:
            return len(installed)
        removedEndpoints = set()
        for fromNumber, toNumber, idx5p, idx3p in reversed(installed):
            fromHelix = helixByNumber.get(fromNumber)
            strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                        if fromHelix is not None else None)
            strand3p = (strand5p.connection3p()
                        if strand5p is not None else None)
            if strand3p is not None and \
                    strand3p.virtualHelix().number() == toNumber and \
                    strand3p.idx5Prime() == idx3p:
                part.removeXover(strand5p, strand3p)
                removedEndpoints.update(
                    ((fromNumber, idx5p), (toNumber, idx3p)))
        removedBoundaries = set()
        for vh in part.getVirtualHelices():
            strands = sorted(list(vh.scaffoldStrandSet()),
                             key=lambda strand: strand.lowIdx())
            for left, right in zip(strands, strands[1:]):
                if left.highIdx() + 1 == right.lowIdx() and \
                        ((vh.number(), left.highIdx()) in
                         removedEndpoints or
                         (vh.number(), right.lowIdx()) in
                         removedEndpoints):
                    removedBoundaries.add(
                        (vh.number(), left.highIdx(), right.lowIdx()))
        _mergeRemovedScaffoldXoverBoundaries(
            part, removedBoundaries)
        return 0

    # Some valid dense layouts require several reciprocal pairs to be
    # installed together: an individual pair can temporarily split the
    # loop, while the complete periodic set restores it (7.json).
    created = installPairs(eligiblePairs) if eligiblePairs else 0
    if created:
        return created
    if 1 < len(eligiblePairs) <= 12:
        # A boundary candidate can be the one pair that changes a valid
        # periodic loop into two components.  Search largest batches first;
        # the 7.json honeycomb target selects eight of nine eligible pairs.
        pairCount = len(eligiblePairs)
        for subsetSize in range(pairCount - 1, 1, -1):
            for mask in range(1, 1 << pairCount):
                if bin(mask).count('1') != subsetSize:
                    continue
                subset = [
                    eligiblePairs[index]
                    for index in range(pairCount)
                    if mask & (1 << index)]
                created = installPairs(subset)
                if created:
                    return created
    for pair in eligiblePairs:
        created += installPairs([pair])
    return created


def _singleSeamScaffoldRecordSet(part, paths, targetLength):
    """Return a dense multi-module route with one global seam.

    Closing every module before splicing the rings together leaves one seam
    per module.  Fill the non-global seams and greedily move complete
    reciprocal crossover blocks until the original main-loop length is
    restored with only one sparse route pair.  The topology search is shared
    by Honeycomb and Square and keeps their respective hard spacings.
    """
    if part._step not in (21, 32) or \
            len(paths) < 2 or targetLength <= 0:
        return None
    directionSpacing = 6 if part._step == 21 else 7
    minimum = part.minBaseIdx()
    maximum = part.maxBaseIdx()
    length = maximum - minimum + 1
    pathNumbers = set(number for path in paths for number in path)
    virtualHelices = [
        vh for vh in part.getVirtualHelices()
        if vh.number() in pathNumbers]
    for vh in virtualHelices:
        strands = list(vh.scaffoldStrandSet())
        if not strands or min(strand.lowIdx() for strand in strands) != \
                minimum or max(strand.highIdx() for strand in strands) != \
                maximum:
            return None

    directions = dict(
        (vh.number(), vh.scaffoldStrandSet().isDrawn5to3())
        for vh in virtualHelices)
    existing = set(part._existingScaffoldCrossoverRecords())
    potentials = set()
    for vh in virtualHelices:
        strandSet = vh.scaffoldStrandSet()
        drawn5to3 = strandSet.isDrawn5to3()
        for neighbor, index, strandType, isLowIdx in \
                part.potentialCrossoverList(vh):
            if strandType != StrandType.Scaffold:
                continue
            fromHelixIs5p = ((isLowIdx and drawn5to3) or
                             (not isLowIdx and not drawn5to3))
            if fromHelixIs5p:
                potentials.add(
                    (vh.number(), neighbor.number(), index, index))

    seamPairs = [
        tuple(sorted(path[-2:])) for path in paths if len(path) >= 2]
    if len(seamPairs) < 2:
        return None
    nonGlobalSeams = set(seamPairs[:-1])
    denseSeamRecords = set()
    for pair in nonGlobalSeams:
        current = sorted(
            record for record in existing
            if tuple(sorted(record[:2])) == pair)
        if len(current) < 2:
            return None
        low = min(record[2] for record in current)
        high = max(record[2] for record in current)
        denseSeamRecords.update(current)
        denseSeamRecords.update(
            record for record in potentials
            if tuple(sorted(record[:2])) == pair and
            low <= record[2] <= high and
            min((record[2] - low) % part._step,
                (low - record[2]) % part._step) <= 1)

    helixNumbers = sorted(directions)
    rank = dict((number, index)
                for index, number in enumerate(helixNumbers))
    nodeCount = len(helixNumbers) * length

    def topology(records):
        outgoing = [-1] * nodeCount
        for number, drawn5to3 in directions.items():
            offset = rank[number] * length
            for relative in range(length):
                nextRelative = relative + (1 if drawn5to3 else -1)
                outgoing[offset + relative] = (
                    offset + nextRelative
                    if 0 <= nextRelative < length else -1)
        for unused_from, toNumber, unused_fromIndex, toIndex in records:
            relative = toIndex - minimum
            predecessor = relative - (
                1 if directions[toNumber] else -1)
            if 0 <= predecessor < length:
                outgoing[rank[toNumber] * length + predecessor] = -1
        for fromNumber, toNumber, fromIndex, toIndex in records:
            outgoing[
                rank[fromNumber] * length + fromIndex - minimum] = \
                rank[toNumber] * length + toIndex - minimum

        seen = bytearray(nodeCount)
        cycles = []
        openCount = 0
        for start in range(nodeCount):
            if seen[start]:
                continue
            order = {}
            route = []
            current = start
            while current >= 0 and current not in order and \
                    not seen[current]:
                order[current] = len(route)
                route.append(current)
                current = outgoing[current]
            for node in route:
                seen[node] = 1
            if current in order:
                cycles.append(len(route) - order[current])
            else:
                openCount += 1
        return sorted(cycles), openCount

    def recordsAreValid(records):
        sources = set()
        targets = set()
        directedEvents = defaultdict(list)
        helixEvents = defaultdict(list)
        for record in records:
            source = (record[0], record[2])
            target = (record[1], record[3])
            if source in sources or target in targets:
                return False
            sources.add(source)
            targets.add(target)
            directedKey = (record[0], record[1])
            if any(abs(record[2] - index) < part._step
                   for index in directedEvents[directedKey]):
                return False
            directedEvents[directedKey].append(record[2])
            for helix, index, neighbor in (
                    (record[0], record[2], record[1]),
                    (record[1], record[3], record[0])):
                if any(otherNeighbor != neighbor and
                       abs(index - otherIndex) < directionSpacing
                       for otherIndex, otherNeighbor in
                       helixEvents[helix]):
                    return False
                helixEvents[helix].append((index, neighbor))
        return True

    localPairs = set(
        tuple(sorted(pair)) for path in paths
        for pair in zip(path, path[1:]))

    def seamCount(records):
        counts = dict((pair, 0) for pair in localPairs)
        for record in records:
            pair = tuple(sorted(record[:2]))
            if pair in counts:
                counts[pair] += 1
        return sum(1 for count in counts.values() if count <= 2)

    reciprocalBlocks = set()
    pool = list(potentials | existing)
    for first in pool:
        for second in pool:
            if first[0] == second[1] and \
                    first[1] == second[0] and \
                    abs(first[2] - second[2]) == 1:
                reciprocalBlocks.add(tuple(sorted(
                    (first, second), key=lambda item: item[2])))
    reciprocalBlocks = sorted(reciprocalBlocks)
    forced = frozenset(denseSeamRecords)
    state = frozenset(existing | denseSeamRecords)

    seenStates = set([state])
    for unused_depth in range(16):
        cycles, unused_open = topology(state)
        if cycles == [targetLength] and seamCount(state) == 1:
            return set(state)
        choices = []
        stateSet = set(state)
        for block in reciprocalBlocks:
            blockSet = set(block)
            if blockSet.issubset(stateSet):
                if blockSet.intersection(forced):
                    continue
                candidate = stateSet.difference(blockSet)
            elif blockSet.isdisjoint(stateSet):
                candidate = stateSet.union(blockSet)
                if not recordsAreValid(candidate):
                    continue
            else:
                continue
            frozen = frozenset(candidate)
            if frozen in seenStates:
                continue
            candidateCycles, openCount = topology(candidate)
            if sum(candidateCycles) != targetLength:
                continue
            score = (-len(candidateCycles), -seamCount(candidate),
                     candidateCycles[-1], -openCount, -len(candidate))
            choices.append((score, frozen, block))
        if not choices:
            return None
        unused_score, state, unused_chosenBlock = max(
            choices, key=lambda item: item[0])
        seenStates.add(state)
    return None


def _applyScaffoldCrossoverRecordSet(part, targetRecords):
    """Replace scaffold crossovers while preserving native strand bases."""
    targetRecords = set(targetRecords)
    existing = set(part._existingScaffoldCrossoverRecords())

    # Install complete reciprocal additions first.  Removing the displaced
    # blocks first opens the perimeter into many oligos; subsequent splits
    # can then strand the short interval between adjacent reciprocal sites.
    # With the additions in place, removing the old blocks is a pure
    # topology rebalance and every removed longitudinal boundary can be
    # healed deterministically.
    helixByNumber = dict(
        (vh.number(), vh) for vh in part.getVirtualHelices())
    for fromNumber, toNumber, idx5p, idx3p in sorted(
            targetRecords.difference(existing),
            key=lambda record: (record[2], record[0], record[1])):
        fromHelix = helixByNumber.get(fromNumber)
        toHelix = helixByNumber.get(toNumber)
        strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                    if fromHelix is not None else None)
        strand3p = (toHelix.scaffoldStrandSet().getStrand(idx3p)
                    if toHelix is not None else None)
        if not part._canCreateScaffoldXover(
                strand5p, strand3p, idx5p):
            return False
        part.createXover(strand5p, idx5p, strand3p, idx3p)

    removable = existing.difference(targetRecords)
    removedEndpoints = set()
    removableEndpoints = set(
        (number, index)
        for fromNumber, toNumber, idx5p, idx3p in removable
        for number, index in ((fromNumber, idx5p), (toNumber, idx3p)))
    removedBoundaries = set()
    for vh in part.getVirtualHelices():
        strands = sorted(list(vh.scaffoldStrandSet()),
                         key=lambda strand: strand.lowIdx())
        for left, right in zip(strands, strands[1:]):
            if left.highIdx() + 1 == right.lowIdx() and \
                    ((vh.number(), left.highIdx()) in
                     removableEndpoints or
                     (vh.number(), right.lowIdx()) in
                     removableEndpoints):
                removedBoundaries.add(
                    (vh.number(), left.highIdx(), right.lowIdx()))

    for fromNumber, toNumber, idx5p, unused_idx3p in sorted(
            removable):
        fromHelix = part.virtualHelix(fromNumber)
        strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                    if fromHelix is not None else None)
        strand3p = (strand5p.connection3p()
                    if strand5p is not None else None)
        if strand3p is None or \
                strand3p.virtualHelix().number() != toNumber:
            continue
        part.removeXover(strand5p, strand3p)
        removedEndpoints.update(
            ((fromNumber, idx5p), (toNumber, strand3p.idx5Prime())))

    # A neighboring removal can expose an additional matching split only
    # after its predecessor has been removed, so supplement the boundaries
    # captured above with the post-removal strand layout.
    for vh in part.getVirtualHelices():
        strands = sorted(list(vh.scaffoldStrandSet()),
                         key=lambda strand: strand.lowIdx())
        for left, right in zip(strands, strands[1:]):
            if left.highIdx() + 1 == right.lowIdx() and \
                    ((vh.number(), left.highIdx()) in removedEndpoints or
                     (vh.number(), right.lowIdx()) in removedEndpoints):
                removedBoundaries.add(
                    (vh.number(), left.highIdx(), right.lowIdx()))
    _mergeRemovedScaffoldXoverBoundaries(part, removedBoundaries)
    return set(part._existingScaffoldCrossoverRecords()) == targetRecords


def _regularizeScaffoldBridgePairs(part, modulePaths, spacing,
                                    directionSpacing,
                                    avoidMultiplesOfEight=False):
    """Align and maximize periodic xovers on selected module bridge pairs.

    Bridge minimization concerns distinct helix-pair types, not the number of
    reciprocal blocks on a chosen pair.  Keep the selected pair identities
    fixed, then fill their best legal periodic register.
    """
    if len(modulePaths) < 2 or spacing <= 0:
        return 0
    groupByHelix = dict(
        (number, groupIndex)
        for groupIndex, path in enumerate(modulePaths)
        for number in path)
    existing = set(part._existingScaffoldCrossoverRecords())
    bridgePairs = sorted(set(
        tuple(sorted(record[:2]))
        for record in existing
        if groupByHelix.get(record[0]) !=
           groupByHelix.get(record[1])))
    if not bridgePairs:
        return 0

    loopLengthsBefore = sorted(
        oligo.length() for oligo in _scaffoldOligos(part)
        if oligo.isLoop())
    if not loopLengthsBefore:
        return 0

    potentials = set()
    for vh in part.getVirtualHelices():
        strandSet = vh.scaffoldStrandSet()
        drawn5to3 = strandSet.isDrawn5to3()
        for neighbor, index, strandType, isLowIdx in \
                part.potentialCrossoverList(vh):
            if strandType != StrandType.Scaffold:
                continue
            fromHelixIs5p = (
                (isLowIdx and drawn5to3) or
                (not isLowIdx and not drawn5to3))
            if fromHelixIs5p:
                potentials.add(
                    (vh.number(), neighbor.number(), index, index))

    def respectsDirectionSpacing(records):
        events = defaultdict(list)
        for fromNumber, toNumber, fromIndex, toIndex in records:
            for number, index, neighbor in (
                    (fromNumber, fromIndex, toNumber),
                    (toNumber, toIndex, fromNumber)):
                if any(
                        index == otherIndex or
                        (neighbor != otherNeighbor and
                         abs(index - otherIndex) < directionSpacing)
                        for otherIndex, otherNeighbor in events[number]):
                    return False
                events[number].append((index, neighbor))
        return True

    target = set(existing)
    for bridgePair in bridgePairs:
        pairExisting = set(
            record for record in existing
            if tuple(sorted(record[:2])) == bridgePair)
        pool = set(pairExisting)
        pool.update(
            record for record in potentials
            if tuple(sorted(record[:2])) == bridgePair)

        reciprocalBlocks = set()
        for first in pool:
            for second in pool:
                if first[0] == second[1] and \
                        first[1] == second[0] and \
                        abs(first[2] - second[2]) == 1:
                    reciprocalBlocks.add(tuple(sorted(
                        (first, second), key=lambda item: item[2])))
        if not reciprocalBlocks:
            continue

        otherRecords = set(
            record for record in target
            if tuple(sorted(record[:2])) != bridgePair)
        registers = defaultdict(list)
        for block in reciprocalBlocks:
            lowerRecord = min(block, key=lambda record: record[2])
            # A bridge belongs between completed periodic horizontal blocks.
            # Do not seed it in the incomplete first edge period, where an
            # otherwise legal reciprocal block can split the global loop.
            if lowerRecord[2] < part.minBaseIdx() + part._step:
                continue
            registers[(
                lowerRecord[2] % spacing,
                lowerRecord[0], lowerRecord[1])].append(block)

        choices = []
        for register, blocks in registers.items():
            phase = register[0]
            blocks = sorted(set(blocks),
                            key=lambda block: min(
                                record[2] for record in block))
            accepted = []
            for block in blocks:
                trial = set(otherRecords)
                trial.update(
                    record for acceptedBlock in accepted
                    for record in acceptedBlock)
                trial.update(block)
                if respectsDirectionSpacing(trial):
                    accepted.append(block)
            if not accepted:
                continue
            # The full register can change loop parity at a physical edge.
            # Evaluate its maximal contiguous subsets as well; this keeps as
            # many exactly periodic bridge blocks as topology permits.
            for firstIndex in range(len(accepted)):
                for lastIndex in range(
                        len(accepted), firstIndex, -1):
                    subset = accepted[firstIndex:lastIndex]
                    acceptedRecords = set(
                        record for block in subset for record in block)
                    byDirection = defaultdict(list)
                    for record in acceptedRecords:
                        byDirection[record[:2]].append(record[2])
                    if any(
                            right - left < spacing
                            for indices in byDirection.values()
                            for left, right in
                            zip(sorted(indices), sorted(indices)[1:])):
                        continue
                    preservedBlocks = sum(
                        1 for block in subset
                        if set(block).issubset(pairExisting))
                    avoided = sum(
                        1 for record in acceptedRecords
                        if avoidMultiplesOfEight and record[2] % 8 == 0)
                    indices = sorted(
                        min(record[2] for record in block)
                        for block in subset)
                    edgeWaste = (
                        indices[0] - part.minBaseIdx() +
                        part.maxBaseIdx() - indices[-1])
                    choices.append((
                        -len(subset), -preservedBlocks, avoided,
                        edgeWaste, phase, -firstIndex,
                        acceptedRecords))
        if not choices:
            continue
        currentTarget = set(target)
        for choice in sorted(choices, key=lambda item: item[:-1]):
            chosenRecords = choice[-1]
            trialTarget = set(currentTarget)
            trialTarget.difference_update(pairExisting)
            trialTarget.update(chosenRecords)
            if trialTarget == currentTarget or \
                    not respectsDirectionSpacing(trialTarget):
                continue
            if not _applyScaffoldCrossoverRecordSet(
                    part, trialTarget):
                _applyScaffoldCrossoverRecordSet(
                    part, currentTarget)
                continue
            loopLengthsAfter = sorted(
                oligo.length() for oligo in _scaffoldOligos(part)
                if oligo.isLoop())
            if loopLengthsAfter == loopLengthsBefore and \
                    _scaffoldDensityIsValid(part, part._step):
                target = trialTarget
                break
            _applyScaffoldCrossoverRecordSet(
                part, currentTarget)

    return len(target.symmetric_difference(existing))


def _autoScaffoldDenseRegularFast(part, minimumIndex, directionSpacing,
                                   avoidMultiplesOfEight,
                                   rebuildExisting):
    """Run the original deterministic dense route for regular full panels.

    This is the stable 1/32 Square and 1/21 Honeycomb behavior: route one
    serpentine path, keep both legal boundary blocks on its final seam, and
    install the complete selected phase.  Later boundary work replaced this
    with a multi-variant global search even when no explicit boundary exists,
    which was both much slower and could fall back to minimum density.
    Return ``None`` when the simple route genuinely does not cover the
    connected design, allowing the more general fallback to handle it.
    """
    undoStack = part.undoStack()
    beforeIndex = undoStack.index()
    initialRecords = part._existingScaffoldCrossoverRecords()
    removed = 0
    util.beginSuperMacro(part, desc="AutoCS_scaffolds dense regular")
    try:
        if rebuildExisting:
            removedBoundaries = set()
            removableEndpoints = set(
                (number, index)
                for fromNumber, toNumber, idx5p, idx3p in initialRecords
                for number, index in (
                    (fromNumber, idx5p), (toNumber, idx3p)))
            for vh in part.getVirtualHelices():
                strands = sorted(list(vh.scaffoldStrandSet()),
                                 key=lambda strand: strand.lowIdx())
                for left, right in zip(strands, strands[1:]):
                    if left.highIdx() + 1 == right.lowIdx() and \
                            ((vh.number(), left.highIdx()) in
                             removableEndpoints or
                             (vh.number(), right.lowIdx()) in
                             removableEndpoints):
                        removedBoundaries.add(
                            (vh.number(), left.highIdx(), right.lowIdx()))
            for fromNumber, toNumber, idx5p, unused_idx3p in initialRecords:
                fromHelix = part.virtualHelix(fromNumber)
                strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                            if fromHelix is not None else None)
                strand3p = (strand5p.connection3p()
                            if strand5p is not None else None)
                if strand3p is None or \
                        strand3p.virtualHelix().number() != toNumber:
                    continue
                part.removeXover(strand5p, strand3p)
                removed += 1
            _mergeRemovedScaffoldXoverBoundaries(
                part, removedBoundaries)

        existing = part._existingScaffoldCrossoverRecords()
        virtualHelices = sorted(
            part.getVirtualHelices(), key=lambda item: item.number())
        helixRecords = []
        helixCoordinates = {}
        for vh in virtualHelices:
            if not any(True for unused_strand in
                       vh.scaffoldStrandSet()):
                continue
            row, column = vh.coord()
            helixCoordinates[vh.number()] = (row, column)
            neighbors = [neighbor.number() for neighbor in
                         part.getVirtualHelixNeighbors(vh)
                         if neighbor is not None]
            helixRecords.append(
                (vh.number(), row, column, neighbors))
        paths = _autoScaffoldSnakePaths(helixRecords)
        routePairs = set(
            tuple(sorted(pair))
            for path in paths for pair in zip(path, path[1:]))
        candidates = []
        seen = set()
        for vh in virtualHelices:
            strandSet = vh.scaffoldStrandSet()
            is5to3 = strandSet.isDrawn5to3()
            for neighbor, idx, strandType, isLowIdx in \
                    part.potentialCrossoverList(vh):
                if strandType != StrandType.Scaffold or \
                        (minimumIndex is not None and idx < minimumIndex):
                    continue
                fromHelixIs5p = ((isLowIdx and is5to3) or
                                 (not isLowIdx and not is5to3))
                record = (vh.number(), neighbor.number(), idx, idx)
                if not fromHelixIs5p or record in seen or \
                        tuple(sorted(record[:2])) not in routePairs:
                    continue
                strand5p = strandSet.getStrand(idx)
                strand3p = neighbor.scaffoldStrandSet().getStrand(idx)
                if not part._canCreateScaffoldXover(
                        strand5p, strand3p, idx):
                    continue
                seen.add(record)
                candidates.append(record)
        candidates = _filterAutoScaffoldCandidatesForPaths(
            candidates, paths, part._step, avoidMultiplesOfEight,
            part.minBaseIdx(), part.maxBaseIdx(), helixCoordinates,
            legacyEdgeTails=True,
            reversePhaseOrder=bool(
                paths and not
                part.virtualHelix(paths[0][0]).scaffoldStrandSet(
                ).isDrawn5to3()))
        selected = _selectAutoScaffoldCrossoverRecords(
            candidates, existing, part._step, directionSpacing,
            avoidMultiplesOfEight)
        helixByNumber = dict(
            (vh.number(), vh) for vh in virtualHelices)
        for fromNumber, toNumber, idx5p, idx3p in selected:
            fromHelix = helixByNumber.get(fromNumber)
            toHelix = helixByNumber.get(toNumber)
            strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                        if fromHelix is not None else None)
            strand3p = (toHelix.scaffoldStrandSet().getStrand(idx3p)
                        if toHelix is not None else None)
            if not part._canCreateScaffoldXover(
                    strand5p, strand3p, idx5p):
                continue
            part.createXover(strand5p, idx5p, strand3p, idx3p)
        # A generated crossover is useful only when it belongs to the
        # retained main scaffold loop.  Open tails and secondary components
        # remain ordinary longitudinal scaffold without crossovers.
        _removeAutoCrossoversOutsideScaffoldLoops(
            part, protectedRecords=existing)
    finally:
        util.endSuperMacro(part)

    oligos = _scaffoldOligos(part)
    loops = [oligo for oligo in oligos if oligo.isLoop()]
    scaffoldHelices = set(record[0] for record in helixRecords)
    totalLength = sum(oligo.length() for oligo in oligos)
    minimumCoverage = max(
        0, totalLength - part._step * len(scaffoldHelices))
    coveringLoops = [
        oligo for oligo in loops
        if scaffoldHelices.issubset(set(
            strand.virtualHelix().number()
            for strand in oligo.strand5p().generator3pStrand())) and
        oligo.length() >= minimumCoverage]
    densityValid = _scaffoldDensityIsValid(part, part._step)
    if not coveringLoops or not densityValid:
        while undoStack.index() > beforeIndex:
            undoStack.undo()
        return None

    finalRecords = part._existingScaffoldCrossoverRecords()
    created = max(0, len(finalRecords) -
                  (0 if rebuildExisting else len(initialRecords)))
    mainLoopLength = max(oligo.length() for oligo in coveringLoops)
    return {
        'success': True,
        'created': created,
        'removed': removed,
        'components': len(oligos),
        'is_loop': len(oligos) == 1 and oligos[0].isLoop(),
        'spacing': part._step,
        'requested_spacing': part._step,
        'main_loop_length': mainLoopLength,
        'minimum_spacing': _minimumDirectedScaffoldSpacing(part),
        'density_exceptions': 0,
        'density_deficit': 0,
        'extension_bases': 0,
        'message': (
            '已用快速确定性路线生成 1/%d bp scaffold 布局：'
            '%d 个 crossover，主闭环 %d nt。' %
            (part._step, created, mainLoopLength))}


def _mergeClosedScaffoldLoops(part, spacing, directionSpacing,
                               minimumIndex=None,
                               avoidMultiplesOfEight=False,
                               allowedPairs=None):
    """Splice separate scaffold loops with one reciprocal crossover pair."""
    createdTotal = 0
    while True:
        loopOligos = [oligo for oligo in _scaffoldOligos(part)
                      if oligo.isLoop()]
        if len(loopOligos) < 2:
            return createdTotal

        existing = part._existingScaffoldCrossoverRecords()
        pairEvents = defaultdict(list)
        helixEvents = defaultdict(list)
        for fromHelix, toHelix, fromIndex, toIndex in existing:
            pairEvents[(fromHelix, toHelix)].append(fromIndex)
            helixEvents[fromHelix].append((fromIndex, toHelix))
            helixEvents[toHelix].append((toIndex, fromHelix))

        candidates = []
        seen = set()
        for vh in part.getVirtualHelices():
            strandSet = vh.scaffoldStrandSet()
            is5to3 = strandSet.isDrawn5to3()
            for neighbor, index, strandType, isLowIdx in \
                    part.potentialCrossoverList(vh):
                if strandType != StrandType.Scaffold or \
                        (minimumIndex is not None and index < minimumIndex):
                    continue
                fromHelixIs5p = ((isLowIdx and is5to3) or
                                 (not isLowIdx and not is5to3))
                record = (vh.number(), neighbor.number(), index, index)
                pair = tuple(sorted(record[:2]))
                if not fromHelixIs5p or record in seen:
                    continue
                if allowedPairs is not None and pair not in allowedPairs:
                    continue
                strand5p = strandSet.getStrand(index)
                strand3p = neighbor.scaffoldStrandSet().getStrand(index)
                if not part._canCreateScaffoldXover(
                                            strand5p, strand3p, index):
                    continue
                if strand5p.oligo() is strand3p.oligo() or \
                        not strand5p.oligo().isLoop() or \
                        not strand3p.oligo().isLoop():
                    continue
                if not _autoScaffoldCandidateFits(
                        record, pairEvents, helixEvents,
                        spacing, directionSpacing):
                    continue
                seen.add(record)
                candidates.append(record)

        candidateSet = set(candidates)
        reciprocalPairs = []
        seenPairs = set()
        for record in sorted(candidates,
                             key=lambda item: (item[2],
                                               item[0], item[1])):
            reverse = [
                other for other in candidateSet
                if other[0] == record[1] and other[1] == record[0] and
                abs(other[2] - record[2]) == 1]
            for partner in reverse:
                pair = tuple(sorted((record, partner),
                                    key=lambda item: item[2]))
                if pair in seenPairs:
                    continue
                seenPairs.add(pair)
                lowerRecord = pair[0]
                lowerHelix = part.virtualHelix(lowerRecord[0])
                if lowerHelix is None or \
                        not lowerHelix.scaffoldStrandSet().isDrawn5to3():
                    continue
                isolation = min(
                    [abs(item[2] - existingIndex)
                     for item in pair
                     for helix in (item[0], item[1])
                     for existingIndex, unused_neighbor in
                     helixEvents.get(helix, ())]
                    or [part.maxBaseIdx() + 1])
                reciprocalPairs.append((
                    abs(record[0] - record[1]),
                    sum(1 for item in pair
                        if avoidMultiplesOfEight and item[2] % 8 == 0),
                    -isolation,
                    min(item[2] for item in pair),
                    pair))

        if not reciprocalPairs:
            return createdTotal

        merged = False
        beforeLoopCount = len(loopOligos)
        for unused_numberGap, unused_avoided, unused_isolation, \
                unused_index, pair in sorted(reciprocalPairs,
                                             key=lambda item: item[:4]):
            installed = []
            firstRecord = pair[0]
            firstFromHelix = part.virtualHelix(firstRecord[0])
            firstToHelix = part.virtualHelix(firstRecord[1])
            firstStrand5p = firstFromHelix.scaffoldStrandSet().getStrand(
                                                        firstRecord[2])
            firstStrand3p = firstToHelix.scaffoldStrandSet().getStrand(
                                                        firstRecord[3])
            expectedMergedLength = (firstStrand5p.oligo().length() +
                                    firstStrand3p.oligo().length())
            for fromNumber, toNumber, idx5p, idx3p in pair:
                fromHelix = part.virtualHelix(fromNumber)
                toHelix = part.virtualHelix(toNumber)
                strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                            if fromHelix is not None else None)
                strand3p = (toHelix.scaffoldStrandSet().getStrand(idx3p)
                            if toHelix is not None else None)
                if not part._canCreateScaffoldXover(
                                            strand5p, strand3p, idx5p):
                    break
                part.createXover(strand5p, idx5p,
                                 strand3p, idx3p)
                currentStrand5p = \
                    fromHelix.scaffoldStrandSet().getStrand(idx5p)
                currentStrand3p = (currentStrand5p.connection3p()
                                   if currentStrand5p is not None else None)
                if currentStrand3p is None or \
                        currentStrand3p.virtualHelix().number() != \
                        toNumber or \
                        currentStrand3p.idx5Prime() != idx3p:
                    break
                installed.append((fromNumber, toNumber,
                                  idx5p, idx3p))
            loopsAfter = [oligo for oligo in _scaffoldOligos(part)
                          if oligo.isLoop()]
            if len(installed) == 2 and \
                    len(loopsAfter) == beforeLoopCount - 1 and \
                    any(oligo.length() == expectedMergedLength
                        for oligo in loopsAfter):
                createdTotal += 2
                merged = True
                break
            for fromNumber, toNumber, idx5p, idx3p in reversed(installed):
                fromHelix = part.virtualHelix(fromNumber)
                strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                            if fromHelix is not None else None)
                strand3p = (strand5p.connection3p()
                            if strand5p is not None else None)
                if strand3p is not None and \
                        strand3p.virtualHelix().number() == toNumber and \
                        strand3p.idx5Prime() == idx3p:
                    part.removeXover(strand5p, strand3p)
        if not merged:
            return createdTotal


def _sparsifyClosedScaffoldLoop(part, protectedRecords=()):
    """Remove redundant reciprocal pairs while preserving one global loop."""
    protected = set(protectedRecords)

    def reciprocalPairs():
        records = set(part._existingScaffoldCrossoverRecords())
        pairs = []
        used = set()
        for record in sorted(records, key=lambda item:
                             (item[2], item[0], item[1])):
            if record in used or record in protected:
                continue
            reverse = [
                other for other in records
                if other not in used and other not in protected and
                other[0] == record[1] and other[1] == record[0] and
                abs(other[2] - record[2]) == 1]
            if not reverse:
                continue
            partner = min(reverse, key=lambda item:
                          (abs(item[2] - record[2]), item[2]))
            pairs.append(tuple(sorted((record, partner),
                                      key=lambda item: item[2])))
            used.add(record)
            used.add(partner)
        return pairs

    def tryRemove(pairs):
        removedRecords = []
        for pair in pairs:
            for fromNumber, toNumber, idx5p, idx3p in pair:
                fromHelix = part.virtualHelix(fromNumber)
                strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                            if fromHelix is not None else None)
                strand3p = (strand5p.connection3p()
                            if strand5p is not None else None)
                if strand3p is None or \
                        strand3p.virtualHelix().number() != toNumber or \
                        strand3p.idx5Prime() != idx3p:
                    break
                part.removeXover(strand5p, strand3p)
                removedRecords.append(
                    (fromNumber, toNumber, idx5p, idx3p))
            else:
                continue
            break
        expected = 2 * len(pairs)
        oligos = _scaffoldOligos(part)
        if len(removedRecords) == expected and \
                len(oligos) == 1 and oligos[0].isLoop():
            return expected
        for fromNumber, toNumber, idx5p, idx3p in removedRecords:
            fromHelix = part.virtualHelix(fromNumber)
            toHelix = part.virtualHelix(toNumber)
            strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                        if fromHelix is not None else None)
            strand3p = (toHelix.scaffoldStrandSet().getStrand(idx3p)
                        if toHelix is not None else None)
            if part._canCreateScaffoldXover(
                    strand5p, strand3p, idx5p):
                part.createXover(strand5p, idx5p, strand3p, idx3p)
        return 0

    pairs = reciprocalPairs()
    # Removing one reciprocal pair can temporarily change loop parity even
    # when two such pairs are jointly redundant.  Try the maximal simultaneous
    # removal first; this is the common "two edge crossovers only" topology.
    removed = tryRemove(pairs) if pairs else 0
    if removed:
        return removed

    # For irregular outlines, retain any indispensable pairs and remove the
    # remaining pairs greedily when each deletion independently preserves the
    # single loop.
    removed = 0
    while True:
        changed = False
        for pair in reciprocalPairs():
            count = tryRemove([pair])
            if count:
                removed += count
                changed = True
                break
        if not changed:
            return removed


def _thinClosedScaffoldLoopToSpacing(
        part, preferredSpacing, minimumSpacing,
        protectedRecords=()):
    """Thin a proven dense loop toward a sparse preferred spacing.

    Reciprocal crossover pairs are removed only when the same largest loop
    and its length survive.  The score first eliminates gaps below the
    preferred spacing, then avoids overshooting it.  This preserves the
    reliable dense route while implementing the agreed "preferred sparse
    period, native period hard floor" rule.
    """
    protected = set(protectedRecords)
    loops = [oligo for oligo in _scaffoldOligos(part) if oligo.isLoop()]
    if not loops:
        return 0
    targetLength = max(oligo.length() for oligo in loops)

    # Test large removal configurations in a compact base-level topology
    # model before touching Qt objects.  The former implementation tried up
    # to 4096 configurations by repeatedly removing and recreating real
    # crossovers, which made sparse fallback unacceptably slow.
    virtualHelices = [
        vh for vh in part.getVirtualHelices()
        if any(True for unused_strand in vh.scaffoldStrandSet())]
    minimum = part.minBaseIdx()
    maximum = part.maxBaseIdx()
    spanLength = maximum - minimum + 1
    simulationAvailable = bool(virtualHelices)
    directions = {}
    for vh in virtualHelices:
        strands = list(vh.scaffoldStrandSet())
        if not strands or min(strand.lowIdx() for strand in strands) != \
                minimum or max(strand.highIdx() for strand in strands) != \
                maximum:
            simulationAvailable = False
            break
        directions[vh.number()] = \
            vh.scaffoldStrandSet().isDrawn5to3()
    helixNumbers = sorted(directions)
    rank = dict((number, index)
                for index, number in enumerate(helixNumbers))

    def simulatedCycles(records):
        if not simulationAvailable:
            return None
        nodeCount = len(helixNumbers) * spanLength
        outgoing = [-1] * nodeCount
        for number, drawn5to3 in directions.items():
            offset = rank[number] * spanLength
            for relative in range(spanLength):
                nextRelative = relative + (1 if drawn5to3 else -1)
                outgoing[offset + relative] = (
                    offset + nextRelative
                    if 0 <= nextRelative < spanLength else -1)
        for unused_from, toNumber, unused_fromIndex, toIndex in records:
            if toNumber not in rank:
                return []
            relative = toIndex - minimum
            predecessor = relative - (
                1 if directions[toNumber] else -1)
            if 0 <= predecessor < spanLength:
                outgoing[rank[toNumber] * spanLength + predecessor] = -1
        for fromNumber, toNumber, fromIndex, toIndex in records:
            if fromNumber not in rank or toNumber not in rank:
                return []
            outgoing[
                rank[fromNumber] * spanLength +
                fromIndex - minimum] = \
                rank[toNumber] * spanLength + toIndex - minimum

        seen = bytearray(nodeCount)
        cycles = []
        for start in range(nodeCount):
            if seen[start]:
                continue
            order = {}
            route = []
            current = start
            while current >= 0 and current not in order and \
                    not seen[current]:
                order[current] = len(route)
                route.append(current)
                current = outgoing[current]
            for node in route:
                seen[node] = 1
            if current in order:
                cycles.append(len(route) - order[current])
        return sorted(cycles)

    def reciprocalPairs():
        records = set(part._existingScaffoldCrossoverRecords())
        pairs = []
        used = set()
        for record in sorted(records, key=lambda item:
                             (item[2], item[0], item[1])):
            if record in used or record in protected:
                continue
            reverse = [
                other for other in records
                if other not in used and other not in protected and
                other[0] == record[1] and other[1] == record[0] and
                abs(other[2] - record[2]) == 1]
            if not reverse:
                continue
            partner = min(reverse, key=lambda item:
                          (abs(item[2] - record[2]), item[2]))
            pair = tuple(sorted((record, partner),
                                key=lambda item: item[2]))
            fromHelix = part.virtualHelix(pair[0][0])
            strand = (fromHelix.scaffoldStrandSet().getStrand(pair[0][2])
                      if fromHelix is not None else None)
            if strand is None or not strand.oligo().isLoop() or \
                    strand.oligo().length() != targetLength:
                continue
            pairs.append(pair)
            used.update(pair)
        return pairs

    def densityScore(records):
        byPair = defaultdict(list)
        for fromHelix, toHelix, fromIndex, unused_toIndex in records:
            byPair[(fromHelix, toHelix)].append(fromIndex)
        deficit = 0
        overshoot = 0
        for indices in byPair.values():
            for left, right in zip(sorted(indices), sorted(indices)[1:]):
                gap = right - left
                if gap < minimumSpacing:
                    return None
                if gap < preferredSpacing:
                    deficit += preferredSpacing - gap
                else:
                    overshoot += gap - preferredSpacing
        return (deficit, overshoot, len(records))

    def removePairs(pairs):
        flatRecords = [
            record for pair in pairs for record in pair]
        removed = []
        boundaries = set()
        endpoints = set(
            (number, index)
            for fromNumber, toNumber, idx5p, idx3p in flatRecords
            for number, index in ((fromNumber, idx5p),
                                  (toNumber, idx3p)))
        for vh in part.getVirtualHelices():
            strands = sorted(list(vh.scaffoldStrandSet()),
                             key=lambda strand: strand.lowIdx())
            for left, right in zip(strands, strands[1:]):
                if left.highIdx() + 1 == right.lowIdx() and \
                        ((vh.number(), left.highIdx()) in endpoints or
                         (vh.number(), right.lowIdx()) in endpoints):
                    boundaries.add(
                        (vh.number(), left.highIdx(), right.lowIdx()))
        for fromNumber, toNumber, idx5p, idx3p in flatRecords:
            fromHelix = part.virtualHelix(fromNumber)
            strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                        if fromHelix is not None else None)
            strand3p = (strand5p.connection3p()
                        if strand5p is not None else None)
            if strand3p is None or \
                    strand3p.virtualHelix().number() != toNumber or \
                    strand3p.idx5Prime() != idx3p:
                return False
            part.removeXover(strand5p, strand3p)
            removed.append((fromNumber, toNumber, idx5p, idx3p))
        _mergeRemovedScaffoldXoverBoundaries(part, boundaries)
        loopsAfter = [oligo for oligo in _scaffoldOligos(part)
                      if oligo.isLoop()]
        accepted = any(oligo.length() == targetLength
                       for oligo in loopsAfter)
        if accepted:
            return True
        for fromNumber, toNumber, idx5p, idx3p in removed:
            fromHelix = part.virtualHelix(fromNumber)
            toHelix = part.virtualHelix(toNumber)
            strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                        if fromHelix is not None else None)
            strand3p = (toHelix.scaffoldStrandSet().getStrand(idx3p)
                        if toHelix is not None else None)
            if part._canCreateScaffoldXover(
                    strand5p, strand3p, idx5p):
                part.createXover(strand5p, idx5p,
                                 strand3p, idx3p)
        return False

    removedTotal = 0
    baseRecords = part._existingScaffoldCrossoverRecords()
    groupedPairs = defaultdict(list)
    for pair in reciprocalPairs():
        groupedPairs[tuple(sorted(pair[0][:2]))].append(pair)
    ratio = max(1, int(round(
        float(preferredSpacing) / float(minimumSpacing))))
    optionGroups = []
    for unused_helixPair, pairs in sorted(groupedPairs.items()):
        ordered = sorted(pairs, key=lambda pair:
                         min(record[2] for record in pair))
        if len(ordered) < 2:
            continue
        options = []
        for keptOffset in range(min(ratio, len(ordered))):
            removals = tuple(
                pair for index, pair in enumerate(ordered)
                if index % ratio != keptOffset)
            if removals:
                options.append(removals)
        if options:
            optionGroups.append(options)

    configurations = [()]
    for options in optionGroups:
        if len(configurations) * len(options) > 4096:
            break
        configurations = [
            previous + option
            for previous in configurations for option in options]
    rankedConfigurations = []
    for configuration in configurations:
        removedRecords = set(
            record for pair in configuration for record in pair)
        score = densityScore([
            record for record in baseRecords
            if record not in removedRecords])
        if score is not None:
            rankedConfigurations.append((
                score, -len(removedRecords), configuration))
    for unused_score, unused_removedCount, configuration in sorted(
            rankedConfigurations,
            key=lambda item: (item[0], item[1])):
        remainingRecords = [
            record for record in baseRecords
            if all(record not in pair for pair in configuration)]
        if simulationAvailable and \
                simulatedCycles(remainingRecords) != [targetLength]:
            continue
        if configuration and removePairs(configuration):
            return 2 * len(configuration)

    while True:
        currentRecords = part._existingScaffoldCrossoverRecords()
        currentScore = densityScore(currentRecords)
        choices = []
        for pair in reciprocalPairs():
            remaining = [
                record for record in currentRecords
                if record not in pair]
            score = densityScore(remaining)
            if score is None or score >= currentScore:
                continue
            choices.append((score, pair))
        changed = False
        for unused_score, pair in sorted(
                choices, key=lambda item: (
                    item[0], min(record[2] for record in item[1]),
                    item[1])):
            if removePairs([pair]):
                removedTotal += 2
                changed = True
                break
        if not changed:
            return removedTotal


def _scaffoldOligos(part):
    return [oligo for oligo in part.oligos()
            if oligo.strand5p() is not None and not oligo.isStaple()]


def _removeAutoCrossoversOutsideScaffoldLoops(part, protectedRecords=()):
    """Keep generated crossovers only on the largest scaffold loop.

    A regular AutoCS result may intentionally preserve unmatched scaffold
    outside its closed route.  Periodic crossovers on open tails or secondary
    local rings add no connectivity to the unique main loop and make one edge
    appear unnecessarily busy.  Existing user crossovers are protected when
    AutoCS was requested without rebuilding them.
    """
    protected = set(protectedRecords)
    loops = [oligo for oligo in _scaffoldOligos(part) if oligo.isLoop()]
    targetLoop = max(loops, key=lambda oligo: oligo.length()) \
        if loops else None
    removable = []
    for record in part._existingScaffoldCrossoverRecords():
        if record in protected:
            continue
        fromHelix = part.virtualHelix(record[0])
        strand5p = (fromHelix.scaffoldStrandSet().getStrand(record[2])
                    if fromHelix is not None else None)
        if strand5p is not None and strand5p.oligo() is not targetLoop:
            removable.append(record)

    removableEndpoints = set(
        (number, index)
        for fromNumber, toNumber, idx5p, idx3p in removable
        for number, index in ((fromNumber, idx5p), (toNumber, idx3p)))
    removedBoundaries = set()
    for vh in part.getVirtualHelices():
        strands = sorted(list(vh.scaffoldStrandSet()),
                         key=lambda strand: strand.lowIdx())
        for left, right in zip(strands, strands[1:]):
            if left.highIdx() + 1 != right.lowIdx():
                continue
            if (vh.number(), left.highIdx()) in removableEndpoints or \
                    (vh.number(), right.lowIdx()) in removableEndpoints:
                removedBoundaries.add(
                    (vh.number(), left.highIdx(), right.lowIdx()))

    removed = 0
    for fromNumber, toNumber, idx5p, idx3p in removable:
        fromHelix = part.virtualHelix(fromNumber)
        strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                    if fromHelix is not None else None)
        strand3p = (strand5p.connection3p()
                    if strand5p is not None else None)
        if strand3p is None or \
                strand3p.virtualHelix().number() != toNumber or \
                strand3p.idx5Prime() != idx3p:
            continue
        # Removing a crossover from an open oligo cannot alter a separate
        # closed oligo, but resolve every record afresh because strand objects
        # can be merged by each removal.
        if strand5p.oligo() is targetLoop:
            continue
        part.removeXover(strand5p, strand3p)
        removed += 1
    _mergeRemovedScaffoldXoverBoundaries(part, removedBoundaries)
    return removed


def _removeScaffoldCrossoversOutsideClosedLoops(part,
                                                protectedRecords=()):
    """Remove generated crossovers from every non-loop scaffold component."""
    protected = set(protectedRecords)
    removable = []
    for record in part._existingScaffoldCrossoverRecords():
        if record in protected:
            continue
        fromHelix = part.virtualHelix(record[0])
        strand5p = (fromHelix.scaffoldStrandSet().getStrand(record[2])
                    if fromHelix is not None else None)
        if strand5p is not None and not strand5p.oligo().isLoop():
            removable.append(record)

    removableEndpoints = set(
        (number, index)
        for fromNumber, toNumber, idx5p, idx3p in removable
        for number, index in ((fromNumber, idx5p), (toNumber, idx3p)))
    removedBoundaries = set()
    for vh in part.getVirtualHelices():
        strands = sorted(list(vh.scaffoldStrandSet()),
                         key=lambda strand: strand.lowIdx())
        for left, right in zip(strands, strands[1:]):
            if left.highIdx() + 1 == right.lowIdx() and (
                    (vh.number(), left.highIdx()) in removableEndpoints or
                    (vh.number(), right.lowIdx()) in removableEndpoints):
                removedBoundaries.add(
                    (vh.number(), left.highIdx(), right.lowIdx()))

    removed = 0
    for fromNumber, toNumber, idx5p, idx3p in removable:
        fromHelix = part.virtualHelix(fromNumber)
        strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                    if fromHelix is not None else None)
        strand3p = (strand5p.connection3p()
                    if strand5p is not None else None)
        if strand3p is None or \
                strand3p.virtualHelix().number() != toNumber or \
                strand3p.idx5Prime() != idx3p:
            continue
        if strand5p.oligo().isLoop():
            continue
        part.removeXover(strand5p, strand3p)
        removed += 1
    _mergeRemovedScaffoldXoverBoundaries(part, removedBoundaries)
    return removed


def _minimumDirectedScaffoldSpacing(part):
    """Return the smallest repeated same-direction scaffold xover gap."""
    byPair = defaultdict(list)
    for fromHelix, toHelix, fromIndex, unused_toIndex in \
            part._existingScaffoldCrossoverRecords():
        byPair[(fromHelix, toHelix)].append(fromIndex)
    gaps = [right - left
            for indices in byPair.values()
            for left, right in zip(sorted(indices), sorted(indices)[1:])]
    return min(gaps) if gaps else None


def _scaffoldDensityIsValid(part, spacing, exemptPairs=()):
    """Check the final route, including endpoint-closing crossovers."""
    byPair = defaultdict(list)
    for fromHelix, toHelix, fromIndex, unused_toIndex in \
            part._existingScaffoldCrossoverRecords():
        byPair[(fromHelix, toHelix)].append(fromIndex)
    for pair, indices in byPair.items():
        if tuple(sorted(pair)) in exemptPairs:
            continue
        ordered = sorted(indices)
        if any(right - left < spacing
               for left, right in zip(ordered, ordered[1:])):
            return False
    return True


def _hasLegalScaffoldCrossoverCandidate(part, minimumIndex=None):
    """Return whether the occupied scaffold exposes any creatable xover."""
    for vh in part.getVirtualHelices():
        strandSet = vh.scaffoldStrandSet()
        drawn5to3 = strandSet.isDrawn5to3()
        for neighbor, index, strandType, isLowIdx in \
                part.potentialCrossoverList(vh):
            if strandType != StrandType.Scaffold or \
                    (minimumIndex is not None and index < minimumIndex):
                continue
            fromHelixIs5p = (
                (isLowIdx and drawn5to3) or
                (not isLowIdx and not drawn5to3))
            if not fromHelixIs5p:
                continue
            strand5p = strandSet.getStrand(index)
            strand3p = neighbor.scaffoldStrandSet().getStrand(index)
            if part._canCreateScaffoldXover(
                    strand5p, strand3p, index):
                return True
    return False


def _autoScaffoldGlobalLoopQuality(part, preferredSpacing, paths=()):
    """Measure whether AutoCS retained the intended global scaffold loop.

    A small loop that merely visits every helix is not a successful route:
    it can leave most scaffold bases in open fragments.  Treat the largest
    physically connected scaffold region as the design target, require one
    loop to cover every helix in that region, and allow at most one preferred
    density period of open edge scaffold per target helix.  Disconnected
    regions are deliberately excluded so AutoCS does not add crossovers to
    scaffold that cannot belong to the selected global loop.
    """
    scaffoldHelices = set()
    adjacency = {}
    strandLengthByHelix = defaultdict(int)
    for vh in part.getVirtualHelices():
        strands = list(vh.scaffoldStrandSet())
        if not strands:
            continue
        number = vh.number()
        scaffoldHelices.add(number)
        strandLengthByHelix[number] = sum(
            strand.totalLength() for strand in strands)
        adjacency[number] = set(
            neighbor.number()
            for neighbor in part.getVirtualHelixNeighbors(vh)
            if neighbor is not None)

    regions = []
    remaining = set(scaffoldHelices)
    while remaining:
        seed = min(remaining)
        region = set([seed])
        pending = [seed]
        remaining.remove(seed)
        while pending:
            number = pending.pop()
            neighbors = adjacency.get(number, set()).intersection(remaining)
            region.update(neighbors)
            pending.extend(neighbors)
            remaining.difference_update(neighbors)
        regions.append(region)

    pathRank = dict(
        (number, rank)
        for rank, number in enumerate(
            number for path in paths for number in path))
    targetRegion = (
        max(
            regions,
            key=lambda region: (
                len(region),
                sum(strandLengthByHelix[number] for number in region),
                -min(pathRank.get(number, number) for number in region)))
        if regions else set())

    oligos = _scaffoldOligos(part)
    loops = [oligo for oligo in oligos if oligo.isLoop()]
    coveringLoops = []
    targetLength = sum(
        strandLengthByHelix[number] for number in targetRegion)
    minimumCoverage = max(
        0, targetLength - int(preferredSpacing) * len(targetRegion))
    for oligo in loops:
        loopHelices = set(
            strand.virtualHelix().number()
            for strand in oligo.strand5p().generator3pStrand())
        if targetRegion.issubset(loopHelices) and \
                oligo.length() >= minimumCoverage:
            coveringLoops.append(oligo)
    mainLoopLength = max(
        [oligo.length() for oligo in coveringLoops] or [0])

    recordsByPair = defaultdict(list)
    for record in part._existingScaffoldCrossoverRecords():
        recordsByPair[tuple(sorted(record[:2]))].append(record)
    orderedRoutePairs = []
    for path in paths:
        for pair in zip(path, path[1:]):
            pair = tuple(sorted(pair))
            if pair not in orderedRoutePairs:
                orderedRoutePairs.append(pair)
    # A module seam is a sparse final pair relative to the normal pairs in
    # that same module.  Boundary pairs at the beginning of a path and sparse
    # inter-module bridges are not seams.  Counting every module final pair
    # here exposes the old "one seam per local ring" failure instead of
    # hiding all but the designated global seam.
    seamPairs = []
    for path in paths:
        pathPairs = [tuple(sorted(pair))
                     for pair in zip(path, path[1:])]
        pathCounts = [
            len(recordsByPair.get(pair, ())) for pair in pathPairs]
        if not pathPairs or not pathCounts:
            continue
        finalCount = pathCounts[-1]
        if finalCount and finalCount < max(pathCounts):
            seamPairs.append(pathPairs[-1])
    # If two route descriptions share the same final edge, it is still one
    # physical seam.
    seamPairs = list(dict.fromkeys(seamPairs))
    seamPositions = [
        orderedRoutePairs.index(pair) for pair in seamPairs]

    hardDensityValid = _scaffoldDensityIsValid(part, part._step)
    return {
        'valid': bool(coveringLoops) and hardDensityValid,
        'target_helices': sorted(targetRegion),
        'covered_length': mainLoopLength,
        'minimum_covered_length': minimumCoverage,
        'hard_density_valid': hardDensityValid,
        'seam_count': len(seamPairs),
        'seam_pairs': seamPairs,
        'seam_position': max(seamPositions) if seamPositions else None,
    }


def _autoScaffoldUnifiedCandidateScore(part, details, paths, helixRecords,
                                        preferredSpacing):
    """Score one completed AutoCS candidate in the agreed priority order."""
    referenceModules = (
        _autoScaffoldGeneralModulePaths(part, helixRecords)
        if part._step == 32 else None) or paths
    groupByHelix = dict(
        (number, groupIndex)
        for groupIndex, path in enumerate(referenceModules)
        for number in path)
    bridgePairs = set(
        tuple(sorted(record[:2]))
        for record in part._existingScaffoldCrossoverRecords()
        if groupByHelix.get(record[0]) != groupByHelix.get(record[1]))
    bridgeRecords = defaultdict(list)
    for record in part._existingScaffoldCrossoverRecords():
        pair = tuple(sorted(record[:2]))
        if pair in bridgePairs:
            bridgeRecords[pair].append(record)
    bridgeBlockIndices = {}
    for pair, records in bridgeRecords.items():
        blocks = set()
        for first in records:
            for second in records:
                if first[0] == second[1] and \
                        first[1] == second[0] and \
                        abs(first[2] - second[2]) == 1:
                    blocks.add(min(first[2], second[2]))
        bridgeBlockIndices[pair] = sorted(blocks)
    bridgeCounts = [
        len(indices) for indices in bridgeBlockIndices.values()]
    bridgeCountImbalance = (
        max(bridgeCounts) - min(bridgeCounts)
        if bridgeCounts else 0)
    bridgePeriodError = sum(
        abs((right - left) - preferredSpacing)
        for indices in bridgeBlockIndices.values()
        for left, right in zip(indices, indices[1:]))
    bridgeBlockCount = sum(bridgeCounts)
    fragmentation, inversions, orientation = \
        _autoScaffoldModuleOrderPenalty(part, paths, helixRecords)
    densityExceptions = _scaffoldDensityExceptionCount(
        part, preferredSpacing)
    densityDeficit = _scaffoldDensityDeficit(
        part, preferredSpacing)
    seamCount = int(details.get('seam_count', 0))
    seamPosition = details.get('seam_position')
    seamPenalty = abs(seamCount - 1)
    score = (
        0 if details.get('success') else 1,
        0 if details.get('hard_density_valid', False) else 1,
        densityExceptions,
        densityDeficit,
        fragmentation,
        inversions,
        len(bridgePairs),
        bridgePeriodError,
        bridgeCountImbalance,
        -bridgeBlockCount,
        seamPenalty,
        -(seamPosition if seamPosition is not None else -1),
        orientation)
    metrics = {
        'candidate_score': score,
        'module_fragmentation': fragmentation,
        'module_inversions': inversions,
        'module_orientation_reversals': orientation,
        'longitudinal_bridge_pairs': sorted(bridgePairs),
        'longitudinal_bridge_count': len(bridgePairs),
        'longitudinal_bridge_blocks': [
            [pair[0], pair[1], list(indices)]
            for pair, indices in sorted(bridgeBlockIndices.items())
        ],
        'longitudinal_bridge_block_count': bridgeBlockCount,
        'longitudinal_bridge_count_imbalance':
            bridgeCountImbalance,
        'longitudinal_bridge_period_error': bridgePeriodError,
        'preferred_density_exceptions': densityExceptions,
        'preferred_density_deficit': densityDeficit,
    }
    return score, metrics


def _scaffoldDensityExceptionCount(part, preferredSpacing,
                                     exemptPairs=()):
    """Count repeated directed gaps that miss the preferred spacing."""
    byPair = defaultdict(list)
    for fromHelix, toHelix, fromIndex, unused_toIndex in \
            part._existingScaffoldCrossoverRecords():
        byPair[(fromHelix, toHelix)].append(fromIndex)
    return sum(
        1
        for pair, indices in byPair.items()
        if tuple(sorted(pair)) not in exemptPairs
        for left, right in zip(sorted(indices), sorted(indices)[1:])
        if right - left < preferredSpacing)


def _scaffoldDensityDeficit(part, preferredSpacing, exemptPairs=()):
    """Return total bp by which repeated gaps miss the preferred spacing."""
    byPair = defaultdict(list)
    for fromHelix, toHelix, fromIndex, unused_toIndex in \
            part._existingScaffoldCrossoverRecords():
        byPair[(fromHelix, toHelix)].append(fromIndex)
    return sum(
        preferredSpacing - (right - left)
        for pair, indices in byPair.items()
        if tuple(sorted(pair)) not in exemptPairs
        for left, right in zip(sorted(indices), sorted(indices)[1:])
        if right - left < preferredSpacing)


def _extendScaffoldPathEdges(part, paths, minimumIndex=None,
                             legalOnly=False):
    """Extend near-aligned path edges to the final pair's legal phases.

    The last adjacent helix pair supplies the two perimeter closure sites.
    Only components whose physical ends differ by at most one lattice period
    are normalized; strongly angled or unrelated regions are left untouched.
    Returned records identify the added scaffold-only intervals.
    """
    helixByNumber = dict((vh.number(), vh)
                         for vh in part.getVirtualHelices())
    extensionRecords = defaultdict(list)
    addedBases = 0
    for path in paths:
        if len(path) < 2:
            continue
        pathHelices = [helixByNumber[number] for number in path
                       if number in helixByNumber]
        boundary = []
        for vh in pathHelices:
            strands = list(vh.scaffoldStrandSet())
            if not strands:
                continue
            boundary.append((vh,
                             min(strands, key=lambda strand: strand.lowIdx()),
                             max(strands, key=lambda strand: strand.highIdx())))
        if len(boundary) != len(pathHelices):
            continue
        lows = [item[1].lowIdx() for item in boundary]
        highs = [item[2].highIdx() for item in boundary]
        edgesAligned = len(set(lows)) == 1 and len(set(highs)) == 1

        lastFirst = helixByNumber[path[-2]]
        lastSecond = helixByNumber[path[-1]]
        legalPositions = set()
        for source, target in ((lastFirst, lastSecond),
                               (lastSecond, lastFirst)):
            for neighbor, index, strandType, unused_isLow in \
                    part.potentialCrossoverList(source):
                if neighbor is target and strandType == StrandType.Scaffold \
                        and (minimumIndex is None or index >= minimumIndex):
                    legalPositions.add(index)
        lowCandidates = [index for index in legalPositions
                         if part.minBaseIdx() <= index and
                         (index <= min(lows) if edgesAligned else
                          index < min(lows))]
        highCandidates = [index for index in legalPositions
                          if index <= part.maxBaseIdx() and
                          (index >= max(highs) if edgesAligned else
                           index > max(highs))]
        if part._step == 32:
            lowPhaseExtension, highPhaseExtension = 9, 8
        else:
            lowPhaseExtension = highPhaseExtension = 7
        fallbackLow = max(part.minBaseIdx(),
                          min(lows) - lowPhaseExtension)
        fallbackHigh = min(part.maxBaseIdx(),
                           max(highs) + highPhaseExtension)
        if legalOnly:
            targetLow = max(lowCandidates) if lowCandidates else min(lows)
            targetHigh = (min(highCandidates) if highCandidates else
                          max(highs))
        else:
            targetLow = (max(lowCandidates) if lowCandidates else
                         (fallbackLow if minimumIndex is None or
                          fallbackLow >= minimumIndex else min(lows)))
            targetHigh = (min(highCandidates) if highCandidates else
                          fallbackHigh)
        if min(lows) - targetLow > part._step or \
                targetHigh - max(highs) > part._step:
            continue

        for vh, lowStrand, highStrand in boundary:
            oldLow = lowStrand.lowIdx()
            oldHigh = highStrand.highIdx()
            if lowStrand is highStrand:
                newBounds = (min(oldLow, targetLow),
                             max(oldHigh, targetHigh))
                if newBounds != lowStrand.idxs():
                    lowStrand.resize(newBounds)
                    if targetLow < oldLow:
                        extensionRecords[vh.coord()].append(
                                                    (targetLow, oldLow - 1))
                        addedBases += oldLow - targetLow
                    if targetHigh > oldHigh:
                        extensionRecords[vh.coord()].append(
                                                    (oldHigh + 1, targetHigh))
                        addedBases += targetHigh - oldHigh
            else:
                if targetLow < oldLow:
                    lowStrand.resize((targetLow, lowStrand.highIdx()))
                    extensionRecords[vh.coord()].append(
                                                (targetLow, oldLow - 1))
                    addedBases += oldLow - targetLow
                if targetHigh > oldHigh:
                    highStrand.resize((highStrand.lowIdx(), targetHigh))
                    extensionRecords[vh.coord()].append(
                                                (oldHigh + 1, targetHigh))
                    addedBases += targetHigh - oldHigh
    return addedBases, dict(extensionRecords)


def _bridgeGuidedScaffoldEdgeGaps(part):
    """Join a manually extended Guided edge back to its duplex scaffold.

    A user may draw ordinary scaffold beyond the image-derived duplex mask.
    If that extension is separated from the mask-bearing strand by at most one
    lattice period, fill only the intervening gap.  The user's extension
    remains ordinary scaffold; only these automatically inserted bases are
    returned as scaffold-only metadata.
    """
    targets = part._guidedDuplexRegions
    if not targets:
        return 0, {}
    addedBases = 0
    added = defaultdict(list)
    for vh in part.getVirtualHelices():
        targetIntervals = targets.get(vh.coord(), ())
        if not targetIntervals:
            continue
        strandSet = vh.scaffoldStrandSet()
        while True:
            strands = sorted(list(strandSet), key=lambda strand:
                             strand.lowIdx())
            bridged = False
            for left, right in zip(strands, strands[1:]):
                gapLow = left.highIdx() + 1
                gapHigh = right.lowIdx() - 1
                if gapLow > gapHigh or gapHigh - gapLow + 1 > part._step:
                    continue
                leftInTarget = any(left.lowIdx() <= high and
                                   left.highIdx() >= low
                                   for low, high in targetIntervals)
                rightInTarget = any(right.lowIdx() <= high and
                                    right.highIdx() >= low
                                    for low, high in targetIntervals)
                # Exactly one side must be the image-derived body.  This
                # avoids silently healing intentional nicks inside the design.
                if leftInTarget == rightInTarget:
                    continue
                if strandSet.createStrand(gapLow, gapHigh) < 0:
                    continue
                bridge = strandSet.getStrand(gapLow)
                priority = left if leftInTarget else right
                strandSet.mergeStrands(priority, bridge)
                merged = strandSet.getStrand(gapLow)
                other = right if leftInTarget else left
                strandSet.mergeStrands(merged, other)
                added[vh.coord()].append((gapLow, gapHigh))
                addedBases += gapHigh - gapLow + 1
                bridged = True
                break
            if not bridged:
                break
    return addedBases, dict(added)


def _mergeRemovedScaffoldXoverBoundaries(part, boundaries):
    """Heal only same-helix splits that belonged to removed crossovers.

    Removing a crossover leaves its two longitudinal neighbors as separate
    Strand objects.  Rebuilding AutoCS without healing those boundaries turns
    one scaffold run into dozens of artificial oligos and encourages the
    endpoint router to make geometric shortcuts.  ``boundaries`` is captured
    before removal, so intentional user nicks are not included.
    """
    helixByNumber = dict((vh.number(), vh)
                         for vh in part.getVirtualHelices())
    for number, leftIndex, rightIndex in sorted(boundaries):
        vh = helixByNumber.get(number)
        if vh is None:
            continue
        strandSet = vh.scaffoldStrandSet()
        left = strandSet.getStrand(leftIndex)
        right = strandSet.getStrand(rightIndex)
        if left is None or right is None or left is right or \
                left.highIdx() + 1 != right.lowIdx():
            continue
        strandSet.mergeStrands(left, right)


def _pruneCrossoversBetweenStraightScaffoldRuns(
        part, paths, protectedRecords=()):
    """Keep one physical helix-pair bridge between straight scaffold runs.

    Periodic crossovers belong inside each long straight run.  If every
    physical neighbor between two runs is also treated as a density edge, a
    two-row sheet repeatedly jumps between rows (0-5, 1-6, 2-7, 3-8) instead
    of using one clean group bridge.  Build a deterministic spanning tree of
    the straight-run graph and retain only one earliest neighbor pair for
    each tree edge.  Longitudinal registers on that chosen helix pair remain
    periodic.
    """
    paths = [list(path) for path in paths if path]
    if len(paths) < 2:
        return 0
    groupByHelix = {}
    positionByHelix = {}
    for groupIndex, path in enumerate(paths):
        for position, number in enumerate(path):
            groupByHelix[number] = groupIndex
            positionByHelix[number] = position

    physicalEdges = set()
    for vh in part.getVirtualHelices():
        first = vh.number()
        firstGroup = groupByHelix.get(first)
        if firstGroup is None:
            continue
        for neighbor in part.getVirtualHelixNeighbors(vh):
            if neighbor is None:
                continue
            second = neighbor.number()
            secondGroup = groupByHelix.get(second)
            if secondGroup is None or secondGroup == firstGroup:
                continue
            physicalEdges.add(tuple(sorted((first, second))))
    if not physicalEdges:
        return 0

    edgeChoices = []
    for first, second in physicalEdges:
        firstGroup = groupByHelix[first]
        secondGroup = groupByHelix[second]
        groupPair = tuple(sorted((firstGroup, secondGroup)))
        edgeChoices.append((
            groupPair,
            positionByHelix[first] + positionByHelix[second],
            max(positionByHelix[first], positionByHelix[second]),
            first, second))

    chosenPairs = set()
    connectedGroups = set([0])
    remainingGroups = set(range(1, len(paths)))
    while remainingGroups:
        choices = [
            item for item in edgeChoices
            if ((item[0][0] in connectedGroups and
                 item[0][1] in remainingGroups) or
                (item[0][1] in connectedGroups and
                 item[0][0] in remainingGroups))]
        if not choices:
            break
        chosen = min(choices, key=lambda item:
                     (item[1], item[2], item[0], item[3], item[4]))
        chosenPairs.add(tuple(sorted((chosen[3], chosen[4]))))
        connectedGroups.update(chosen[0])
        remainingGroups.difference_update(chosen[0])

    crossRunPairs = physicalEdges
    protected = set(protectedRecords)
    removable = [
        record for record in part._existingScaffoldCrossoverRecords()
        if tuple(sorted(record[:2])) in crossRunPairs and
        tuple(sorted(record[:2])) not in chosenPairs and
        record not in protected]
    if not removable:
        return 0

    removableEndpoints = set(
        (number, index)
        for fromNumber, toNumber, idx5p, idx3p in removable
        for number, index in ((fromNumber, idx5p), (toNumber, idx3p)))
    removedBoundaries = set()
    for vh in part.getVirtualHelices():
        strands = sorted(list(vh.scaffoldStrandSet()),
                         key=lambda strand: strand.lowIdx())
        for left, right in zip(strands, strands[1:]):
            if left.highIdx() + 1 == right.lowIdx() and \
                    ((vh.number(), left.highIdx()) in
                     removableEndpoints or
                     (vh.number(), right.lowIdx()) in
                     removableEndpoints):
                removedBoundaries.add(
                    (vh.number(), left.highIdx(), right.lowIdx()))

    removed = 0
    for fromNumber, toNumber, idx5p, unused_idx3p in removable:
        fromHelix = part.virtualHelix(fromNumber)
        strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                    if fromHelix is not None else None)
        strand3p = (strand5p.connection3p()
                    if strand5p is not None else None)
        if strand3p is None or \
                strand3p.virtualHelix().number() != toNumber:
            continue
        part.removeXover(strand5p, strand3p)
        removed += 1
    _mergeRemovedScaffoldXoverBoundaries(part, removedBoundaries)
    return removed


def _removeScaffoldCrossoversOnPairs(part, pairs, protectedRecords=()):
    """Remove scaffold crossovers on selected helix pairs and heal nicks."""
    pairs = set(tuple(sorted(pair)) for pair in pairs)
    protected = set(protectedRecords)
    removable = [
        record for record in part._existingScaffoldCrossoverRecords()
        if tuple(sorted(record[:2])) in pairs and
        record not in protected]
    if not removable:
        return 0

    removableEndpoints = set(
        (number, index)
        for fromNumber, toNumber, idx5p, idx3p in removable
        for number, index in ((fromNumber, idx5p), (toNumber, idx3p)))
    removedBoundaries = set()
    for vh in part.getVirtualHelices():
        strands = sorted(list(vh.scaffoldStrandSet()),
                         key=lambda strand: strand.lowIdx())
        for left, right in zip(strands, strands[1:]):
            if left.highIdx() + 1 == right.lowIdx() and \
                    ((vh.number(), left.highIdx()) in
                     removableEndpoints or
                     (vh.number(), right.lowIdx()) in
                     removableEndpoints):
                removedBoundaries.add(
                    (vh.number(), left.highIdx(), right.lowIdx()))

    removed = 0
    for fromNumber, toNumber, idx5p, unused_idx3p in removable:
        fromHelix = part.virtualHelix(fromNumber)
        strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                    if fromHelix is not None else None)
        strand3p = (strand5p.connection3p()
                    if strand5p is not None else None)
        if strand3p is None or \
                strand3p.virtualHelix().number() != toNumber:
            continue
        part.removeXover(strand5p, strand3p)
        removed += 1
    _mergeRemovedScaffoldXoverBoundaries(part, removedBoundaries)
    return removed


def _alternatingSquareRowBridgePairSets(part, paths):
    """Return all and selected bridges for alternating Square rows.

    Translating the same outline by one checkerboard parity reverses every
    scaffold direction.  A valid dense phase can then include both outer
    bridges between each pair of horizontal rows.  The reference s1 layout
    instead snakes cleanly: top-to-middle uses the left bridge,
    middle-to-bottom the right bridge, alternating thereafter.
    """
    coordinates = dict(
        (vh.number(), vh.coord()) for vh in part.getVirtualHelices())
    horizontal = []
    for path in paths:
        if not path:
            continue
        rows = [coordinates[number][0] for number in path]
        if len(set(rows)) != 1:
            return set(), set()
        horizontal.append(list(path))
    if len(horizontal) < 2:
        return set(), set()
    horizontal.sort(key=lambda path: (
        coordinates[path[0]][0],
        min(coordinates[number][1] for number in path)))

    groupByHelix = dict(
        (number, groupIndex)
        for groupIndex, path in enumerate(horizontal)
        for number in path)
    crossPairs = set()
    candidatesByGroupPair = defaultdict(set)
    for vh in part.getVirtualHelices():
        first = vh.number()
        firstGroup = groupByHelix.get(first)
        if firstGroup is None:
            continue
        for neighbor in part.getVirtualHelixNeighbors(vh):
            if neighbor is None:
                continue
            second = neighbor.number()
            secondGroup = groupByHelix.get(second)
            if secondGroup is None or secondGroup == firstGroup:
                continue
            pair = tuple(sorted((first, second)))
            crossPairs.add(pair)
            candidatesByGroupPair[
                tuple(sorted((firstGroup, secondGroup)))].add(pair)

    keepPairs = set()
    for upperIndex in range(len(horizontal) - 1):
        groupPair = (upperIndex, upperIndex + 1)
        candidates = candidatesByGroupPair.get(groupPair, ())
        if not candidates:
            continue
        ordered = sorted(candidates, key=lambda pair: (
            (coordinates[pair[0]][1] + coordinates[pair[1]][1]),
            pair))
        keepPairs.add(
            ordered[0] if upperIndex % 2 == 0 else ordered[-1])
    return crossPairs, keepPairs


def _pruneAlternatingSquareRowBridges(
        part, paths, protectedRecords=()):
    """Keep only the selected alternating-side Square row bridges."""
    crossPairs, keepPairs = \
        _alternatingSquareRowBridgePairSets(part, paths)
    removablePairs = crossPairs.difference(keepPairs)
    if not removablePairs:
        return False
    _removeScaffoldCrossoversOnPairs(
        part, removablePairs, protectedRecords)
    return True


def _normalizeHoneycombScaffoldPathGroups(
        part, paths, spacing, directionSpacing,
        minimumIndex=None, protectedRecords=()):
    """Build each honeycomb snake group locally before joining the groups.

    A global loop merge can consume the only legal register needed by one
    snake group.  In the 5.json regression, the early 7-12 group bridge
    blocked the local 4-7/6-7 seam and forced an unrelated 3-6 shortcut.
    Remove those global shortcuts, close every snake path only on its own
    planned edges, then splice the completed path loops through one physical
    inter-group helix pair.
    """
    paths = [list(path) for path in paths if path]
    if len(paths) < 2:
        return False

    groupByHelix = {}
    for groupIndex, path in enumerate(paths):
        for number in path:
            groupByHelix[number] = groupIndex
    routePairs = set(
        tuple(sorted(pair))
        for path in paths for pair in zip(path, path[1:]))
    crossGroupPairs = set()
    for vh in part.getVirtualHelices():
        first = vh.number()
        firstGroup = groupByHelix.get(first)
        if firstGroup is None:
            continue
        for neighbor in part.getVirtualHelixNeighbors(vh):
            if neighbor is None:
                continue
            second = neighbor.number()
            secondGroup = groupByHelix.get(second)
            if secondGroup is None or secondGroup == firstGroup:
                continue
            crossGroupPairs.add(tuple(sorted((first, second))))
    if not crossGroupPairs:
        return False

    offRouteInternalPairs = set()
    hasCrossGroupRecord = False
    for record in part._existingScaffoldCrossoverRecords():
        first, second = record[:2]
        firstGroup = groupByHelix.get(first)
        secondGroup = groupByHelix.get(second)
        pair = tuple(sorted((first, second)))
        if firstGroup is None or secondGroup is None:
            continue
        if firstGroup != secondGroup:
            hasCrossGroupRecord = True
        elif pair not in routePairs:
            offRouteInternalPairs.add(pair)
    if not offRouteInternalPairs or not hasCrossGroupRecord:
        return False

    _removeScaffoldCrossoversOnPairs(
        part, offRouteInternalPairs.union(crossGroupPairs),
        protectedRecords)
    _mergeClosedScaffoldLoops(
        part, spacing, directionSpacing, minimumIndex,
        allowedPairs=routePairs)
    _mergeClosedScaffoldLoops(
        part, spacing, directionSpacing, minimumIndex,
        allowedPairs=crossGroupPairs)

    # The chosen inter-group bridge becomes a normal periodic 1/21 edge.
    # Earlier groups no longer need their local closing seam; only the final
    # snake path keeps its last pair sparse as the one global seam.
    chosenBridgePairs = set(
        tuple(sorted(record[:2]))
        for record in part._existingScaffoldCrossoverRecords()
        if tuple(sorted(record[:2])) in crossGroupPairs)
    finalSeamPairs = set(
        [tuple(sorted(paths[-1][-2:]))] if len(paths[-1]) > 1 else [])
    densePairs = routePairs.difference(finalSeamPairs).union(
                                                chosenBridgePairs)
    densityCandidates = []
    seen = set()
    for vh in part.getVirtualHelices():
        strandSet = vh.scaffoldStrandSet()
        is5to3 = strandSet.isDrawn5to3()
        for neighbor, index, strandType, isLowIdx in \
                part.potentialCrossoverList(vh):
            if strandType != StrandType.Scaffold or \
                    (minimumIndex is not None and index < minimumIndex):
                continue
            fromHelixIs5p = ((isLowIdx and is5to3) or
                             (not isLowIdx and not is5to3))
            record = (vh.number(), neighbor.number(), index, index)
            if not fromHelixIs5p or record in seen or \
                    tuple(sorted(record[:2])) not in densePairs:
                continue
            seen.add(record)
            densityCandidates.append(record)
    _densifyClosedScaffoldLoop(
        part, densityCandidates, spacing, directionSpacing)
    return True


def _endpointCanExtend(strand, index, prime):
    """Return whether one free endpoint can grow outward to ``index``."""
    endpoint = strand.idx3Prime() if prime == '3p' else strand.idx5Prime()
    if index == endpoint:
        return True
    low, high = strand.idxs()
    if endpoint == low:
        if index >= low:
            return False
        lowNeighbor, unused_high = strand.strandSet().getNeighbors(strand)
        return (index >= 0 and
                (lowNeighbor is None or index > lowNeighbor.highIdx()))
    if endpoint == high:
        if index <= high:
            return False
        unused_low, highNeighbor = strand.strandSet().getNeighbors(strand)
        return (index <= strand.part().maxBaseIdx() and
                (highNeighbor is None or index < highNeighbor.lowIdx()))
    return False


def _endpointCrossoverCandidates(part, minimumIndex=None,
                                  maximumExtension=None, allowedPairs=None):
    """Return legal directed links between open scaffold components."""
    oligos = [oligo for oligo in _scaffoldOligos(part)
              if not oligo.isLoop()]
    if not oligos:
        return oligos, {}
    if maximumExtension is None:
        maximumExtension = part._step * 4
    endpoints = {}
    for oligo in oligos:
        strand5p = oligo.strand5p()
        strand3p = list(strand5p.generator3pStrand())[-1]
        endpoints[oligo] = (strand3p, strand5p)

    candidates = defaultdict(list)
    for source in oligos:
        strand3p = endpoints[source][0]
        sourceVh = strand3p.virtualHelix()
        sourceSS = strand3p.strandSet()
        neighbors = part.getVirtualHelixNeighbors(sourceVh)
        for target in oligos:
            if target is source and len(oligos) > 1:
                continue
            strand5p = endpoints[target][1]
            targetVh = strand5p.virtualHelix()
            pair = tuple(sorted((sourceVh.number(), targetVh.number())))
            if allowedPairs is not None and pair not in allowedPairs:
                continue
            try:
                neighborIndex = neighbors.index(targetVh)
            except ValueError:
                continue
            phases = (part._scafL[neighborIndex] if sourceSS.isDrawn5to3()
                      else part._scafH[neighborIndex])
            best = []
            for block in range(0, part.maxBaseIdx() + 1, part._step):
                for phase in phases:
                    index = block + phase
                    if index > part.maxBaseIdx() or \
                            (minimumIndex is not None and
                             index < minimumIndex):
                        continue
                    if not _endpointCanExtend(strand3p, index, '3p') or \
                            not _endpointCanExtend(strand5p, index, '5p'):
                        continue
                    sourceExtension = abs(index - strand3p.idx3Prime())
                    targetExtension = abs(index - strand5p.idx5Prime())
                    if max(sourceExtension, targetExtension) > maximumExtension:
                        continue
                    score = (sourceExtension + targetExtension,
                             max(sourceExtension, targetExtension), index)
                    best.append((score, source, target, index,
                                 strand3p, strand5p))
            # Connectivity depends on the target component, not on which of
            # several equivalent phases was used. Keep only the shortest to
            # prevent duplicate branches from making cycle search explode.
            candidates[source].extend(sorted(best)[:1])
        candidates[source].sort(key=lambda item: item[0])
    return oligos, candidates


def _findEndpointCrossoverCycle(oligos, candidates):
    """Find one directed Hamiltonian cycle through open scaffold oligos."""
    if not oligos:
        return None
    if len(oligos) == 1:
        source = oligos[0]
        selfLinks = [record for record in candidates.get(source, ())
                     if record[2] is source]
        return selfLinks[:1] or None

    def stableKey(oligo):
        strand = oligo.strand5p()
        return (strand.virtualHelix().number(), strand.idx5Prime())

    start = min(oligos, key=stableKey)
    path = [start]
    links = []
    visited = set([start])
    searchCount = [0]

    def visit(current):
        searchCount[0] += 1
        if searchCount[0] > 20000:
            return False
        if len(path) == len(oligos):
            closing = [record for record in candidates.get(current, ())
                       if record[2] is start]
            if closing:
                links.append(closing[0])
                return True
            return False
        options = [record for record in candidates.get(current, ())
                   if record[2] not in visited]
        options.sort(key=lambda record: (
            len([item for item in candidates.get(record[2], ())
                 if item[2] not in visited]), record[0]))
        for record in options:
            target = record[2]
            visited.add(target)
            path.append(target)
            links.append(record)
            if visit(target):
                return True
            links.pop()
            path.pop()
            visited.remove(target)
        return False

    return list(links) if visit(start) else None


def _closeLargestRegularScaffoldComponent(
        part, directionSpacing, minimumIndex=None, allowedPairs=None):
    """Close the longest regular scaffold component without growing bases.

    Regular designs may contain edge scaffold that cannot be part of the main
    closed route.  A legal crossover between two positions of the same open
    oligo closes that main route and leaves the excluded edge pieces as
    ordinary open scaffold tails.
    """
    existing = part._existingScaffoldCrossoverRecords()
    pairEvents = defaultdict(list)
    helixEvents = defaultdict(list)
    for fromHelix, toHelix, fromIndex, toIndex in existing:
        pairEvents[(fromHelix, toHelix)].append(fromIndex)
        helixEvents[fromHelix].append((fromIndex, toHelix))
        helixEvents[toHelix].append((toIndex, fromHelix))

    candidates = []
    orderByOligo = {}
    for vh in part.getVirtualHelices():
        strandSet = vh.scaffoldStrandSet()
        is5to3 = strandSet.isDrawn5to3()
        for neighbor, index, strandType, isLowIdx in \
                part.potentialCrossoverList(vh):
            if strandType != StrandType.Scaffold or \
                    (minimumIndex is not None and index < minimumIndex):
                continue
            pair = tuple(sorted((vh.number(), neighbor.number())))
            if allowedPairs is not None and pair not in allowedPairs:
                continue
            fromHelixIs5p = ((isLowIdx and is5to3) or
                             (not isLowIdx and not is5to3))
            if not fromHelixIs5p:
                continue
            strand5p = strandSet.getStrand(index)
            strand3p = neighbor.scaffoldStrandSet().getStrand(index)
            if not part._canCreateScaffoldXover(
                    strand5p, strand3p, index):
                continue
            oligo = strand5p.oligo()
            if oligo is not strand3p.oligo() or oligo.isLoop():
                continue
            if oligo not in orderByOligo:
                positionOrder = {}
                offset = 0
                for oligoStrand in oligo.strand5p().generator3pStrand():
                    step = 1 if oligoStrand.isDrawn5to3() else -1
                    for baseIndex in range(
                            oligoStrand.idx5Prime(),
                            oligoStrand.idx3Prime() + step, step):
                        positionOrder[(
                            oligoStrand.virtualHelix().number(),
                            baseIndex)] = offset
                        offset += 1
                orderByOligo[oligo] = positionOrder
            positionOrder = orderByOligo[oligo]
            sourceOrder = positionOrder.get((vh.number(), index))
            targetOrder = positionOrder.get((neighbor.number(), index))
            # Only the backward-directed connection closes a loop.  The
            # reciprocal forward connection cuts out an additional open
            # segment instead.
            if sourceOrder is None or targetOrder is None or \
                    targetOrder >= sourceOrder:
                continue
            record = (vh.number(), neighbor.number(), index, index)
            if not _autoScaffoldCandidateFits(
                    record, pairEvents, helixEvents,
                    part._step, directionSpacing):
                continue
            edgeDistance = min(index - part.minBaseIdx(),
                               part.maxBaseIdx() - index)
            loopSpan = sourceOrder - targetOrder + 1
            candidates.append((
                -loopSpan, -oligo.length(), edgeDistance, index,
                vh.number(), neighbor.number(), strand5p, strand3p))
    if not candidates:
        return 0
    unused_span, unused_length, unused_edge, index, unused_from, unused_to, \
        strand5p, strand3p = min(candidates)
    part.createXover(strand5p, index, strand3p, index)
    return int(strand5p.connection3p() is strand3p)


def _applyEndpointCrossoverCycle(part, links):
    """Resize selected edges and close all scaffold paths into one loop."""
    desiredBounds = {}
    originalBounds = {}
    for unused_score, unused_source, unused_target, index, strand3p, strand5p \
            in links:
        for strand in (strand3p, strand5p):
            low, high = desiredBounds.get(strand, strand.idxs())
            originalBounds.setdefault(strand, strand.idxs())
            desiredBounds[strand] = (min(low, index), max(high, index))
    extensionBases = 0
    extensionRecords = defaultdict(list)
    for strand, bounds in desiredBounds.items():
        oldLow, oldHigh = originalBounds[strand]
        extensionBases += max(0, oldLow - bounds[0]) + \
                          max(0, bounds[1] - oldHigh)
        if bounds[0] < oldLow:
            extensionRecords[strand.virtualHelix().coord()].append(
                (bounds[0], oldLow - 1))
        if bounds[1] > oldHigh:
            extensionRecords[strand.virtualHelix().coord()].append(
                (oldHigh + 1, bounds[1]))
        if bounds != strand.idxs():
            strand.resize(bounds)
    created = 0
    for unused_score, unused_source, unused_target, index, strand3p, strand5p \
            in links:
        if strand3p.idx3Prime() != index or strand5p.idx5Prime() != index:
            return created, extensionBases, False, dict(extensionRecords)
        part.createXover(strand3p, index, strand5p, index)
        if strand3p.connection3p() is strand5p:
            created += 1
    oligos = _scaffoldOligos(part)
    success = len(oligos) == 1 and oligos[0].isLoop()
    return created, extensionBases, success, dict(extensionRecords)


class _ScaffoldOnlyRegionsCommand(QUndoCommand):
    """Make scaffold-only metadata participate in the design undo stack."""
    def __init__(self, part, newRecords):
        super(_ScaffoldOnlyRegionsCommand, self).__init__()
        self._part = part
        self._oldRecords = dict(
            (coord, list(intervals))
            for coord, intervals in part._scaffoldOnlyRegions.items())
        self._newRecords = dict(
            (coord, list(intervals))
            for coord, intervals in newRecords.items())

    def redo(self):
        self._part.setScaffoldOnlyRegions(self._newRecords)

    def undo(self):
        self._part.setScaffoldOnlyRegions(self._oldRecords)


class _AutobreakAppliedCommand(QUndoCommand):
    """Keep Autobreak's once-only state in the same undo macro as its nicks."""
    def __init__(self, part, applied):
        super(_AutobreakAppliedCommand, self).__init__()
        self._part = part
        self._oldApplied = part._autobreakStaplesApplied
        self._newApplied = bool(applied)

    def redo(self):
        self._part._autobreakStaplesApplied = self._newApplied

    def undo(self):
        self._part._autobreakStaplesApplied = self._oldApplied


def _autoScaffoldAdjacentOnly(part, minimumIndex=None, densityMultiple=1,
                              minimumDensity=False,
                              rebuildExisting=False):
    """Add scaffold crossovers only between adjacent Path-view helices.

    This intentionally has no route, loop, seam, module, bridge, extension,
    fallback-density, or helix-reordering policy.  Geometry and consecutive
    Path-view order define the eligible helix pairs; the remaining decisions
    are only crossover density and the native direction-spacing constraints.
    """
    initialRecords = part._existingScaffoldCrossoverRecords()
    virtualHelices = [
        vh for vh in part.getVirtualHelices()
        if any(True for unused_strand in vh.scaffoldStrandSet())]
    helixByNumber = dict((vh.number(), vh) for vh in virtualHelices)
    available = set(helixByNumber)

    ordered = []
    for coord in part._importedVHelixOrder or ():
        vh = part.virtualHelixAtCoord(coord)
        if vh is not None and vh.number() in available and \
                vh.number() not in ordered:
            ordered.append(vh.number())
    ordered.extend(sorted(available.difference(ordered)))
    rank = dict((number, index) for index, number in enumerate(ordered))

    physicalNeighbors = dict(
        (vh.number(), set(
            neighbor.number() for neighbor in
            part.getVirtualHelixNeighbors(vh)
            if neighbor is not None and neighbor.number() in available))
        for vh in virtualHelices)
    helixCoordinates = dict(
        (vh.number(), vh.coord()) for vh in virtualHelices)
    eligiblePairs = set(
        tuple(sorted((first, second)))
        for first, second in zip(ordered, ordered[1:])
        if second in physicalNeighbors.get(first, ()))
    if not eligiblePairs:
        return {
            'success': False, 'created': 0, 'removed': 0,
            'spacing': None, 'eligible_pairs': [],
            'message': (
                '左侧几何视图与右侧 Path 顺序中没有同时相邻的 '
                'scaffold helix；未修改设计。')}

    paths = []
    for number in ordered:
        if not paths or tuple(sorted(
                (paths[-1][-1], number))) not in eligiblePairs:
            paths.append([number])
        else:
            paths[-1].append(number)
    paths = [path for path in paths if len(path) > 1]

    removed = 0
    created = 0
    densitySpacing = max(1, part._step * max(1, int(densityMultiple)))
    directionSpacing = 6 if part._step == 21 else 7
    avoidMultiplesOfEight = part._step == 32
    util.beginSuperMacro(part, desc="AutoCS_scaffolds adjacent only")
    try:
        if rebuildExisting and initialRecords:
            if _applyScaffoldCrossoverRecordSet(part, set()):
                removed = len(initialRecords)

        existing = part._existingScaffoldCrossoverRecords()
        drawn5to3 = dict(
            (vh.number(), vh.scaffoldStrandSet().isDrawn5to3())
            for vh in virtualHelices)
        candidates = []
        seen = set()
        for vh in virtualHelices:
            strandSet = vh.scaffoldStrandSet()
            isDrawn5to3 = strandSet.isDrawn5to3()
            for neighbor, index, strandType, isLowIdx in \
                    part.potentialCrossoverList(vh):
                if strandType != StrandType.Scaffold or \
                        neighbor.number() not in available or \
                        tuple(sorted((vh.number(), neighbor.number()))) \
                        not in eligiblePairs or \
                        (minimumIndex is not None and index < minimumIndex):
                    continue
                fromHelixIs5p = (
                    (isLowIdx and isDrawn5to3) or
                    (not isLowIdx and not isDrawn5to3))
                record = (vh.number(), neighbor.number(), index, index)
                if fromHelixIs5p and record not in seen:
                    seen.add(record)
                    candidates.append(record)

        if minimumDensity:
            originalSpans = {}
            for vh in virtualHelices:
                strands = list(vh.scaffoldStrandSet())
                if not strands:
                    continue
                originalSpans[vh.number()] = (
                    min(strand.lowIdx() for strand in strands),
                    max(strand.highIdx() for strand in strands))
            selected = _selectMinimumDensityScaffoldRecords(
                candidates, paths, originalSpans, part._step,
                drawn5to3, existing)
            selectedIndices = [record[2] for record in selected]
            densitySpacing = (
                max(selectedIndices) - min(selectedIndices)
                if len(selectedIndices) > 1 else 1)
        else:
            # Both lattices use the same seam-free periodic construction:
            # 21-bp Honeycomb or 32-bp Square units close independently.
            candidates = _filterAutoScaffoldCandidatesForPaths(
                candidates, paths, part._step, avoidMultiplesOfEight,
                part.minBaseIdx(), part.maxBaseIdx(), helixCoordinates,
                legacyEdgeTails=True,
                reversePhaseOrder=bool(
                    paths and not
                    part.virtualHelix(paths[0][0]).scaffoldStrandSet(
                    ).isDrawn5to3()),
                densitySpacing=densitySpacing)
            selected = _selectAutoScaffoldCrossoverRecords(
                candidates, existing, densitySpacing, directionSpacing,
                avoidMultiplesOfEight, helixOrder=rank)
        for fromNumber, toNumber, idx5p, idx3p in sorted(
                selected,
                key=lambda record: (
                    record[2], rank.get(record[0], len(rank)),
                    rank.get(record[1], len(rank)))):
            fromHelix = helixByNumber.get(fromNumber)
            toHelix = helixByNumber.get(toNumber)
            strand5p = (fromHelix.scaffoldStrandSet().getStrand(idx5p)
                        if fromHelix is not None else None)
            strand3p = (toHelix.scaffoldStrandSet().getStrand(idx3p)
                        if toHelix is not None else None)
            if not part._canCreateScaffoldXover(
                    strand5p, strand3p, idx5p):
                continue
            part.createXover(strand5p, idx5p, strand3p, idx3p)
            created += 1
        _removeScaffoldCrossoversOutsideClosedLoops(
            part, protectedRecords=(
                initialRecords if not rebuildExisting else ()))
    finally:
        util.endSuperMacro(part)

    finalRecords = part._existingScaffoldCrossoverRecords()
    created = (len(finalRecords) if rebuildExisting else
               len(set(finalRecords).difference(initialRecords)))
    removed = len(set(initialRecords).difference(finalRecords))
    oligos = _scaffoldOligos(part)
    loops = [oligo for oligo in oligos if oligo.isLoop()]
    modeText = ('最低密度（最大可用间距）'
                if minimumDensity else '1/%d bp' % densitySpacing)
    return {
        'success': bool(created or removed or finalRecords),
        'created': created,
        'removed': removed,
        'components': len(oligos),
        'is_loop': len(oligos) == 1 and bool(loops),
        'spacing': densitySpacing,
        'requested_spacing': densitySpacing,
        'minimum_spacing': _minimumDirectedScaffoldSpacing(part),
        'hard_density_valid':
            _scaffoldDensityIsValid(part, part._step),
        'eligible_pairs': [list(pair) for pair in sorted(eligiblePairs)],
        'message': (
            '仅在左侧和右侧均相邻的 %d 对 helix 之间，按%s添加 '
            '%d 个 scaffold crossover；删除 %d 个原 crossover。'
            '两端未使用的 scaffold 保持不变。' %
            (len(eligiblePairs), modeText, created, removed))}


class Part(QObject):
    """
    A Part is a group of VirtualHelix items that are on the same lattice.
    Parts are the model component that most directly corresponds to a
    DNA origami design.

    Parts are always parented to the document.
    Parts know about their oligos, and the internal geometry of a part
    Copying a part recursively copies all elements in a part:
        VirtualHelices, Strands, etc

    PartInstances are parented to either the document or an assembly
    PartInstances know global position of the part
    Copying a PartInstance only creates a new PartInstance with the same
    Part(), with a mutable parent and position field.
    """

    _step = 21  # this is the period (in bases) of the part lattice
    _radius = 1.125  # nanometers
    _turnsPerStep = 2
    _helicalPitch = _step / _turnsPerStep
    _twistPerBase = 360 / _helicalPitch  # degrees

    def __init__(self, *args, **kwargs):
        """
        Sets the parent document, sets bounds for part dimensions, and sets up
        bookkeeping for partInstances, Oligos, VirtualHelix's, and helix ID
        number assignment.
        """
        if self.__class__ == Part:
            e = "This class is abstract. Perhaps you want HoneycombPart."
            raise NotImplementedError(e)
        self._document = kwargs.get('document', None)
        super(Part, self).__init__(parent=self._document)
        # Data structure
        self._insertions = defaultdict(dict)  # dict of insertions per virtualhelix
        self._oligos = set()
        self._coordToVirtualHelix = {}
        self._numberToVirtualHelix = {}
        # Dimensions
        self._maxRow = 50  # subclass overrides based on prefs
        self._maxCol = 50
        self._minBase = 0
        self._maxBase = int(2 * self._step - 1)
        # ID assignment
        self.oddRecycleBin, self.evenRecycleBin = [], []
        self.reserveBin = set()
        self._highestUsedOdd = -1  # Used in _reserveHelixIDNumber
        self._highestUsedEven = -2  # same
        self._importedVHelixOrder = None
        # ``None`` means a normal design: staples follow every scaffold base.
        # Guided Design stores its original silhouette so edge extensions can
        # remain unpaired scaffold ssDNA.
        self._guidedDuplexRegions = None
        # AutoCS edge extensions are scaffold-only in both regular and Guided
        # Design workflows.
        self._scaffoldOnlyRegions = {}
        # Autobreak is once per generated staple layout. AutoCS/Add staples
        # reset it, and the state participates in undo with the topology.
        self._autobreakStaplesApplied = False
        # Runtime state
        self._activeBaseIndex = self._step
        self._activeVirtualHelix = None
        self._activeVirtualHelixIdx = None

    # end def

    def __repr__(self):
        clsName = self.__class__.__name__
        return "<%s %s>" % (clsName, str(id(self))[-4:])

    ### SIGNALS ###
    partActiveSliceIndexSignal = pyqtSignal(QObject, int)  # self, index
    partActiveSliceResizeSignal = pyqtSignal(QObject)      # self
    partDimensionsChangedSignal = pyqtSignal(QObject)      # self
    partInstanceAddedSignal = pyqtSignal(QObject)          # self
    partParentChangedSignal = pyqtSignal(QObject)          # self
    partPreDecoratorSelectedSignal = pyqtSignal(object, int, int, int)  # row,col,idx
    partRemovedSignal = pyqtSignal(QObject)                # self
    partStrandChangedSignal = pyqtSignal(object, QObject)          # self, virtualHelix
    partVirtualHelixAddedSignal = pyqtSignal(object, QObject)      # self, virtualhelix
    partVirtualHelixRenumberedSignal = pyqtSignal(object, tuple)   # self, coord
    partVirtualHelixResizedSignal = pyqtSignal(object, tuple)      # self, coord
    partVirtualHelicesReorderedSignal = pyqtSignal(object, list)   # self, list of coords
    partHideSignal = pyqtSignal(QObject)
    partActiveVirtualHelixChangedSignal = pyqtSignal(QObject, QObject)

    ### SLOTS ###

    ### ACCESSORS ###
    def document(self):
        return self._document
    # end def

    def oligos(self):
        return self._oligos
    # end def

    def setDocument(self, document):
        self._document = document
    # end def

    def stepSize(self):
        return self._step
    # end def

    def subStepSize(self):
        """Note: _subStepSize is defined in subclasses."""
        return self._subStepSize
    # end def

    def undoStack(self):
        return self._document.undoStack()
    # end def

    ### PUBLIC METHODS FOR QUERYING THE MODEL ###
    def virtualHelix(self, vhref, returnNoneIfAbsent=True):
        # vhrefs are the shiny new way to talk to part about its constituent
        # virtualhelices. Wherever you see f(...,vhref,...) you can
        # f(...,27,...)         use the virtualhelix's id number
        # f(...,vh,...)         use an actual virtualhelix
        # f(...,(1,42),...)     use the coordinate representation of its position
        """A vhref is the number of a virtual helix, the (row, col) of a virtual helix,
        or the virtual helix itself. For conveniece, CRUD should now work with any of them."""
        vh = None
        if type(vhref) in (int, int):
            vh = self._numberToVirtualHelix.get(vhref, None)
        elif type(vhref) in (tuple, list):
            vh = self._coordToVirtualHelix.get(vhref, None)
        else:
            vh = vhref
        if not isinstance(vh, VirtualHelix):
            if returnNoneIfAbsent:
                return None
            else:
                err = "Couldn't find the virtual helix in part %s "+\
                      "referenced by index %s" % (self, vhref)
                raise IndexError(err)
        return vh

    def activeBaseIndex(self):
        return self._activeBaseIndex
    # end def

    def activeVirtualHelix(self):
        return self._activeVirtualHelix
     # end def

    def activeVirtualHelixIdx(self):
        return self._activeVirtualHelixIdx
     # end def

    def dimensions(self):
        """Returns a tuple of the max X and maxY coordinates of the lattice."""
        return self.latticeCoordToPositionXY(self._maxRow, self._maxCol)
    # end def

    def getStapleSequences(self, includeHeader=True, includeLattice=True):
        """Return the unmodified legacy cadnano staple CSV output."""
        s = ("Start,End,Sequence,Length,Color\n"
             if includeHeader else "")
        for oligo in self._oligos:
            if oligo.strand5p().strandSet().isStaple():
                s = s + oligo.sequenceExport(
                    includeLattice=includeLattice)
        return s

    def getInputSequenceRows(self):
        """Return applied scaffold sequences in the staple export schema."""
        records = []
        for oligo in self._oligos:
            strand5p = oligo.strand5p()
            if strand5p is None or not strand5p.strandSet().isScaffold():
                continue
            sequence = oligo.sequence()
            if sequence and sequence.strip():
                records.append(oligo.sequenceRecord())
        return sorted(records, key=_sequenceRowSortKey)

    def getScaffoldSequenceTemplateRows(self):
        """Return every importable scaffold with a blank Sequence cell."""
        records = []
        for oligo in self._oligos:
            strand5p = oligo.strand5p()
            if (strand5p is None or oligo.isLoop() or
                    not strand5p.strandSet().isScaffold()):
                continue
            start, end = oligo.sequenceEndpoints()
            records.append((start, end, '', oligo.actualLength(),
                            oligo.color()))
        return sorted(records, key=lambda record: record[0])

    def getOutputSequenceRows(self):
        """Return staple sequences in the common spreadsheet schema."""
        records = []
        for oligo in self._oligos:
            strand5p = oligo.strand5p()
            if strand5p is not None and strand5p.strandSet().isStaple():
                records.append(oligo.sequenceRecord())
        return sorted(records, key=_sequenceRowSortKey)

    def stapleContinuousRunStatistics(self, minimum=16):
        """Return fraction with a continuous original-coordinate run."""
        total = qualified = 0
        for oligo in self._oligos:
            strand5p = oligo.strand5p()
            if strand5p is None or not strand5p.strandSet().isStaple():
                continue
            if oligo.isHybrid():
                continue
            total += 1
            if _hasContinuousPositionRun(
                    _stapleOligoBaseRecords(oligo), minimum):
                qualified += 1
        fraction = float(qualified) / total if total else 0.0
        return qualified, total, fraction

    def stapleLengthDistribution(self, binSize=10, firstBin=10):
        """Return real-time staple-length bins and proportions.

        Bins start at 10--19 by design.  Percentages use every staple as the
        denominator, so an unusual staple shorter than 10 bases is not hidden
        by changing the denominator.
        """
        lengths = [oligo.actualLength() for oligo in self._oligos
                   if oligo.strand5p() is not None and oligo.isStaple()
                   and not oligo.isHybrid()]
        total = len(lengths)
        if not lengths or max(lengths) < firstBin:
            return total, []
        lastBin = (max(lengths) // binSize) * binSize
        records = []
        for low in range(firstBin, lastBin + 1, binSize):
            high = low + binSize - 1
            count = sum(1 for length in lengths if low <= length <= high)
            fraction = float(count) / total if total else 0.0
            records.append((low, high, count, fraction))
        return total, records

    def getScaffoldSequenceRecords(self):
        """Return applied scaffold sequences and their 5' start positions.

        The records are intentionally made only from model state so they can
        be used by both the JSON encoder and sequence workbook exporter.
        """
        records = []
        for oligo in self._oligos:
            strand5p = oligo.strand5p()
            if strand5p is None or not strand5p.strandSet().isScaffold():
                continue
            sequence = oligo.sequence()
            if not sequence or not sequence.strip():
                continue
            records.append({
                'start_vh': strand5p.virtualHelix().number(),
                'start_idx': strand5p.idx5Prime(),
                'sequence': sequence
            })
        return sorted(records,
                      key=lambda record: (record['start_vh'],
                                          record['start_idx']))

    def getScaffoldColorRecords(self):
        """Return the persistent color and model location of each scaffold.

        Scaffold colors are design appearance data, not sequence data, so
        callers must save these records regardless of whether sequences are
        included in the file.
        """
        records = []
        for oligo in self._oligos:
            strand5p = oligo.strand5p()
            if strand5p is None or not strand5p.strandSet().isScaffold():
                continue
            records.append({
                'start_vh': strand5p.virtualHelix().number(),
                'start_idx': strand5p.idx5Prime(),
                'color': str(oligo.color())
            })
        return sorted(records,
                      key=lambda record: (record['start_vh'],
                                          record['start_idx']))

    def getVirtualHelices(self):
        """yield an iterator to the virtualHelix references in the part"""
        return list(self._coordToVirtualHelix.values())
    # end def

    def indexOfRightmostNonemptyBase(self):
        """
        During reduction of the number of bases in a part, the first click
        removes empty bases from the right hand side of the part (red
        left-facing arrow). This method returns the new numBases that will
        effect that reduction.
        """
        ret = self._step - 1
        for vh in self.getVirtualHelices():
            ret = max(ret, vh.indexOfRightmostNonemptyBase())
        return ret
    # end def

    def insertions(self):
        """Return dictionary of insertions."""
        return self._insertions
    # end def

    def setGuidedDuplexRegions(self, records):
        """Set/clear the image-derived staple mask for Guided Design."""
        if records is None:
            self._guidedDuplexRegions = None
            self._refreshScaffoldOnlyAppearance()
            return
        normalized = defaultdict(list)
        for coord, intervals in records.items():
            coord = tuple(coord)
            for low, high in intervals:
                low, high = int(low), int(high)
                if high >= low:
                    normalized[coord].append((low, high))
        self._guidedDuplexRegions = dict(
            (coord, sorted(intervals))
            for coord, intervals in normalized.items())
        self._refreshScaffoldOnlyAppearance()

    def guidedDuplexRegionRecords(self):
        if self._guidedDuplexRegions is None:
            return None
        records = []
        for coord, intervals in sorted(self._guidedDuplexRegions.items()):
            vh = self.virtualHelixAtCoord(coord)
            if vh is None:
                continue
            records.append({'vh': vh.number(), 'row': coord[0],
                            'col': coord[1],
                            'intervals': [[low, high]
                                          for low, high in intervals]})
        return records

    def setScaffoldOnlyRegions(self, records):
        # scaffold-only is a Guided Design boundary annotation.  Older
        # regular files may contain records created by previous AutoCS
        # versions; discard those records so their scaffold is again ordinary
        # staple-eligible scaffold.
        if self._guidedDuplexRegions is None:
            records = {}
        normalized = defaultdict(list)
        for coord, intervals in (records or {}).items():
            for low, high in intervals:
                low, high = int(low), int(high)
                if high >= low:
                    normalized[tuple(coord)].append((low, high))
        merged = {}
        for coord, intervals in normalized.items():
            compact = []
            for low, high in sorted(intervals):
                if compact and low <= compact[-1][1] + 1:
                    compact[-1] = (compact[-1][0],
                                   max(compact[-1][1], high))
                else:
                    compact.append((low, high))
            if compact:
                merged[coord] = compact
        self._scaffoldOnlyRegions = merged
        self._refreshScaffoldOnlyAppearance()

    def addScaffoldOnlyRegions(self, records, useUndoStack=False):
        if self._guidedDuplexRegions is None:
            if self._scaffoldOnlyRegions:
                self.setScaffoldOnlyRegions({})
            return
        combined = defaultdict(list)
        for coord, intervals in self._scaffoldOnlyRegions.items():
            combined[coord].extend(intervals)
        for coord, intervals in (records or {}).items():
            combined[tuple(coord)].extend(intervals)
        merged = {}
        for coord, intervals in combined.items():
            normalized = []
            for low, high in sorted(intervals):
                if normalized and low <= normalized[-1][1] + 1:
                    normalized[-1] = (normalized[-1][0],
                                      max(normalized[-1][1], high))
                else:
                    normalized.append((low, high))
            merged[coord] = normalized
        if useUndoStack:
            command = _ScaffoldOnlyRegionsCommand(self, merged)
            util.execCommandList(self, [command], desc=None,
                                 useUndoStack=True)
        else:
            self.setScaffoldOnlyRegions(merged)

    def takeScaffoldOnlyIntervals(self, coord, low, high):
        """Remove and return explicit scaffold-only metadata in one range."""
        coord = tuple(coord)
        removed = []
        remaining = []
        for intervalLow, intervalHigh in self._scaffoldOnlyRegions.get(
                                                    coord, ()):
            overlapLow = max(low, intervalLow)
            overlapHigh = min(high, intervalHigh)
            if overlapLow > overlapHigh:
                remaining.append((intervalLow, intervalHigh))
                continue
            removed.append((overlapLow, overlapHigh))
            if intervalLow < overlapLow:
                remaining.append((intervalLow, overlapLow - 1))
            if overlapHigh < intervalHigh:
                remaining.append((overlapHigh + 1, intervalHigh))
        if removed:
            updated = dict(
                (itemCoord, list(intervals))
                for itemCoord, intervals in
                self._scaffoldOnlyRegions.items())
            if remaining:
                updated[coord] = remaining
            else:
                updated.pop(coord, None)
            self.setScaffoldOnlyRegions(updated)
        return removed

    def _refreshScaffoldOnlyAppearance(self):
        """Refresh Path-view strands after scaffold-only metadata changes."""
        for vh in self.getVirtualHelices():
            for strand in vh.scaffoldStrandSet():
                strand.strandUpdateSignal.emit(strand)

    def scaffoldOnlyIntervals(self, virtualHelix):
        """Return only intervals explicitly created by an automatic action."""
        if self._guidedDuplexRegions is None:
            return ()
        return tuple(self._scaffoldOnlyRegions.get(
                                            virtualHelix.coord(), ()))

    def scaffoldOnlyRegionRecords(self):
        if self._guidedDuplexRegions is None:
            return []
        records = []
        for coord, intervals in sorted(self._scaffoldOnlyRegions.items()):
            vh = self.virtualHelixAtCoord(coord)
            if vh is not None:
                records.append({'vh': vh.number(), 'row': coord[0],
                                'col': coord[1],
                                'intervals': [[low, high]
                                              for low, high in intervals]})
        return records

    def stapleTargetSegments(self, virtualHelix):
        """Return staple-eligible scaffold intervals on one helix."""
        scaffoldSegments = []
        for strand in virtualHelix.scaffoldStrandSet():
            low, high = strand.idxs()
            if scaffoldSegments and scaffoldSegments[-1][1] == low - 1:
                scaffoldSegments[-1] = (scaffoldSegments[-1][0], high)
            else:
                scaffoldSegments.append((low, high))
        if self._guidedDuplexRegions is None:
            segments = scaffoldSegments
        else:
            targets = self._guidedDuplexRegions.get(
                                                virtualHelix.coord(), ())
            segments = [(max(scafLow, targetLow),
                         min(scafHigh, targetHigh))
                        for scafLow, scafHigh in scaffoldSegments
                        for targetLow, targetHigh in targets
                        if max(scafLow, targetLow) <=
                           min(scafHigh, targetHigh)]
        exclusions = (self._scaffoldOnlyRegions.get(
                            virtualHelix.coord(), ())
                      if self._guidedDuplexRegions is not None else ())
        for excludeLow, excludeHigh in exclusions:
            remaining = []
            for low, high in segments:
                if excludeHigh < low or excludeLow > high:
                    remaining.append((low, high))
                    continue
                if low < excludeLow:
                    remaining.append((low, excludeLow - 1))
                if excludeHigh < high:
                    remaining.append((excludeHigh + 1, high))
            segments = remaining
        return segments

    def isEvenParity(self, row, column):
        """Should be overridden when subclassing."""
        raise NotImplementedError
    # end def

    def getStapleLoopOligos(self):
        """
        Returns staple oligos with no 5'/3' ends. Used by
        actionExportStaplesSlot in documentcontroller to validate before
        exporting staple sequences.
        """
        stapLoopOlgs = []
        for o in list(self.oligos()):
            if o.isStaple() and o.isLoop():
                stapLoopOlgs.append(o)
        return stapLoopOlgs

    def hasVirtualHelixAtCoord(self, coord):
        return coord in self._coordToVirtualHelix
    # end def

    def maxBaseIdx(self):
        return self._maxBase
    # end def

    def minBaseIdx(self):
        return self._minBase
    # end def

    def numberOfVirtualHelices(self):
        return len(self._coordToVirtualHelix)
    # end def

    def radius(self):
        return self._radius
    # end def

    def helicalPitch(self):
        return self._helicalPitch
    # end def

    def twistPerBase(self):
        return self._twistPerBase
    # end def

    def virtualHelixAtCoord(self, coord):
        """
        Looks for a virtualHelix at the coordinate, coord = (row, colum)
        if it exists it is returned, else None is returned
        """
        try:
            return self._coordToVirtualHelix[coord]
        except:
            return None
    # end def

    ### PUBLIC METHODS FOR EDITING THE MODEL ###
    def autoFillWithoutCrossovers(part, strandType):
        """Add unconnected strands only where the requested type is empty.

        Existing strands, crossovers, colors, sequences, and nicks are left
        untouched. Scaffold targets every complete virtual helix; staple
        targets only bases occupied by scaffold. New strands fill the empty
        runs inside those target intervals without merging adjacent strands.
        """
        virtualHelices = part.getVirtualHelices()
        if not virtualHelices:
            return 0

        if strandType == StrandType.Scaffold:
            description = "Add scaffolds"
        elif strandType == StrandType.Staple:
            description = "Add staples"
        else:
            raise ValueError("Unsupported strand type")

        emptyRuns = []
        lowIdx = part.minBaseIdx()
        highIdx = part.maxBaseIdx()
        for vh in virtualHelices:
            strandSet = vh.getStrandSetByType(strandType)
            if strandType == StrandType.Scaffold:
                segments = [(lowIdx, highIdx)]
            else:
                segments = part.stapleTargetSegments(vh)
            for lo, hi in segments:
                idx = lo
                while idx <= hi:
                    existing = strandSet.getStrand(idx)
                    if existing is not None:
                        idx = existing.highIdx() + 1
                        continue
                    runLow = idx
                    idx += 1
                    while idx <= hi and strandSet.getStrand(idx) is None:
                        idx += 1
                    emptyRuns.append((strandSet, runLow, idx - 1))

        if not emptyRuns:
            return 0

        created = 0
        util.beginSuperMacro(part, desc=description)
        try:
            for strandSet, lo, hi in emptyRuns:
                canInsert, strandSetIdx = strandSet.getIndexToInsert(lo, hi)
                if not canInsert:
                    continue
                command = StrandSet.CreateStrandCommand(
                    strandSet, lo, hi, strandSetIdx,
                    color=(styles.AUTOMATIC_STAP_COLOR
                           if strandType == StrandType.Staple else None))
                util.execCommandList(part, [command],
                                     desc="Fill empty strand region")
                created += 1
            if created and strandType == StrandType.Staple:
                util.execCommandList(
                    part, [_AutobreakAppliedCommand(part, False)],
                    desc="Reset Autobreak state")
        finally:
            util.endSuperMacro(part)
        return created

    def _existingScaffoldCrossoverRecords(part):
        records = []
        for vh in part.getVirtualHelices():
            for strand in vh.scaffoldStrandSet():
                connected = strand.connection3p()
                if connected is None or \
                        connected.part() is not part or \
                        not connected.strandSet().isScaffold():
                    continue
                records.append((vh.number(),
                                connected.virtualHelix().number(),
                                strand.idx3Prime(),
                                connected.idx5Prime()))
        return records

    def _canCreateScaffoldXover(part, strand5p, strand3p, idx):
        """Return whether a crossover can be created without resizing DNA."""
        if strand5p is None or strand3p is None:
            return False
        if strand5p.hasXoverAt(idx) or strand3p.hasXoverAt(idx):
            return False

        strandSet5p = strand5p.strandSet()
        strandSet3p = strand3p.strandSet()
        if strand5p.idx3Prime() != idx and \
                not strandSet5p.strandCanBeSplit(strand5p, idx):
            return False

        if strand3p.idx5Prime() != idx:
            offset3p = -1 if strandSet3p.isDrawn5to3() else 1
            if not strandSet3p.strandCanBeSplit(strand3p,
                                                idx + offset3p):
                return False
        return True

    def autoScaffoldPathOrder(part):
        """Return scaffold-bearing helix coordinates in geometric route order."""
        virtualHelices = sorted(part.getVirtualHelices(),
                                key=lambda item: item.number())
        helixRecords = []
        helixByNumber = {}
        for vh in virtualHelices:
            if not any(True for unused_strand in vh.scaffoldStrandSet()):
                continue
            row, column = vh.coord()
            helixByNumber[vh.number()] = vh
            neighbors = [neighbor.number() for neighbor in
                         part.getVirtualHelixNeighbors(vh)
                         if neighbor is not None]
            helixRecords.append((vh.number(), row, column, neighbors))
        return [helixByNumber[number].coord()
                for path in _autoScaffoldSnakePaths(helixRecords)
                for number in path]

    def _autoScaffoldCrossoversPreGuided(
            part, densitySpacing=None, minimumIndex=None,
            minimumDensity=False, explicitPaths=None):
        """Route scaffold along dense, branch-free serpentine helix paths."""
        densitySpacing = densitySpacing or part._step
        if part._step == 21:
            # An inclusive seven-base segment has endpoint index difference 6.
            directionSpacing = 6
            avoidMultiplesOfEight = False
        elif part._step == 32:
            # An inclusive eight-base segment has endpoint index difference 7.
            directionSpacing = 7
            avoidMultiplesOfEight = True
        else:
            return 0

        existing = part._existingScaffoldCrossoverRecords()
        candidates = []
        seen = set()
        virtualHelices = sorted(part.getVirtualHelices(),
                                key=lambda item: item.number())
        helixRecords = []
        helixCoordinates = {}
        for vh in virtualHelices:
            if not any(True for unused_strand in vh.scaffoldStrandSet()):
                continue
            row, column = vh.coord()
            helixCoordinates[vh.number()] = (row, column)
            neighbors = [neighbor.number() for neighbor in
                         part.getVirtualHelixNeighbors(vh)
                         if neighbor is not None]
            helixRecords.append((vh.number(), row, column, neighbors))
        if explicitPaths is not None:
            paths = [list(path) for path in explicitPaths if path]
        elif part._step == 32:
            paths = (
                _autoScaffoldRightPanelModulePaths(
                    part, helixRecords) or
                _autoScaffoldStraightPaths(helixRecords))
        else:
            paths = _autoScaffoldSnakePaths(helixRecords)
        filterPaths = paths
        routePairs = set(tuple(sorted(pair)) for path in paths
                         for pair in zip(path, path[1:]))
        originalSpans = {}
        drawn5to3 = {}
        for vh in virtualHelices:
            strands = list(vh.scaffoldStrandSet())
            if not strands:
                continue
            originalSpans[vh.number()] = (
                min(strand.lowIdx() for strand in strands),
                max(strand.highIdx() for strand in strands))
            drawn5to3[vh.number()] = \
                vh.scaffoldStrandSet().isDrawn5to3()

        for vh in virtualHelices:
            strandSet = vh.scaffoldStrandSet()
            is5to3 = strandSet.isDrawn5to3()
            for neighbor, idx, strandType, isLowIdx in \
                    part.potentialCrossoverList(vh):
                if strandType != StrandType.Scaffold or \
                        (minimumIndex is not None and idx < minimumIndex):
                    continue
                fromHelixIs5p = ((isLowIdx and is5to3) or
                                 (not isLowIdx and not is5to3))
                if not fromHelixIs5p:
                    continue
                record = (vh.number(), neighbor.number(), idx, idx)
                if tuple(sorted(record[:2])) not in routePairs:
                    continue
                if record in seen:
                    continue
                strand5p = strandSet.getStrand(idx)
                strand3p = neighbor.scaffoldStrandSet().getStrand(idx)
                if not part._canCreateScaffoldXover(
                                                strand5p, strand3p, idx):
                    continue
                seen.add(record)
                candidates.append(record)

        sparseHoneycombSelection = None
        minimumDensitySelection = None
        if minimumDensity:
            # Lowest density is a geometric optimization, not an artificial
            # large lattice-period multiple.  Use the legal native sites
            # nearest the two overlap boundaries (or the reciprocal pair
            # nearest the overlap center) so each route gap is connected with
            # the theoretically widest available spacing.
            minimumDensitySelection = \
                _selectMinimumDensityScaffoldRecords(
                    candidates, filterPaths, originalSpans,
                    part._step, drawn5to3, existing)
        elif densitySpacing > part._step and part._step == 21:
            sparseHoneycombSelection = \
                _selectHoneycombSparseLoopRecords(
                    candidates, filterPaths, densitySpacing,
                    helixCoordinates)
        elif densitySpacing > part._step:
            candidates = _filterAutoScaffoldSparseCandidatesForPaths(
                candidates, filterPaths, part._step,
                avoidMultiplesOfEight,
                part.minBaseIdx(), part.maxBaseIdx(), helixCoordinates,
                legacyEdgeTails=True,
                densitySpacing=densitySpacing,
                preferBoundaryRegisters=True,
                registerBoundaryTrim=21)
        else:
            candidates = _filterAutoScaffoldCandidatesForPaths(
                candidates, filterPaths, part._step,
                avoidMultiplesOfEight,
                part.minBaseIdx(), part.maxBaseIdx(), helixCoordinates)
        selected = (
            minimumDensitySelection
            if minimumDensitySelection is not None else
            sparseHoneycombSelection
            if sparseHoneycombSelection is not None else
            _selectAutoScaffoldCrossoverRecords(
                candidates, existing, densitySpacing, directionSpacing,
                avoidMultiplesOfEight,
                helixOrder=dict(
                    (number, rank)
                    for rank, number in enumerate(
                        number for path in paths for number in path))))
        if not selected:
            return 0

        helixByNumber = dict((vh.number(), vh) for vh in virtualHelices)
        created = 0
        util.beginSuperMacro(part, desc="AutoCS_scaffolds")
        try:
            for fromNumber, toNumber, idx5p, idx3p in selected:
                currentOligos = _scaffoldOligos(part)
                if len(currentOligos) == 1 and \
                        currentOligos[0].isLoop():
                    break
                fromHelix = helixByNumber.get(fromNumber)
                toHelix = helixByNumber.get(toNumber)
                if fromHelix is None or toHelix is None:
                    continue
                strand5p = fromHelix.scaffoldStrandSet().getStrand(idx5p)
                strand3p = toHelix.scaffoldStrandSet().getStrand(idx3p)
                if not part._canCreateScaffoldXover(
                                            strand5p, strand3p, idx5p):
                    continue
                part.createXover(strand5p, idx5p, strand3p, idx3p)
                created += 1
        finally:
            util.endSuperMacro(part)
        return created

    def autoScaffoldCrossovers(part, minimumIndex=None, densityMultiple=1,
                               routeOnly=False, rebuildExisting=False,
                               returnDetails=False, **unused_options):
        """Add scaffold crossovers using only the three simple user rules."""
        if part._step not in (21, 32):
            details = {
                'success': False, 'created': 0, 'removed': 0,
                'components': 0,
                'message': '当前点阵不支持 AutoCS_scaffolds。'}
            return details if returnDetails else 0
        details = _autoScaffoldAdjacentOnly(
            part, minimumIndex=minimumIndex,
            densityMultiple=densityMultiple,
            minimumDensity=routeOnly,
            rebuildExisting=rebuildExisting)
        return details if returnDetails else details.get('created', 0)

        # Legacy global-loop implementation retained below only as inert
        # source history.  The return above is the sole AutoCS scaffold path.
        directionSpacing = 6 if part._step == 21 else 7
        avoidMultiplesOfEight = part._step == 32

        densityMultiple = max(1, int(densityMultiple))
        requestedSpacing = part._step * densityMultiple
        unifiedInternal = bool(unused_options.get('_unifiedInternal'))
        densitySearchDisabled = bool(
            unused_options.get('_densitySearchDisabled'))
        explicitPaths = unused_options.get('_explicitPaths')
        explicitModules = bool(unused_options.get('_explicitModules'))
        if routeOnly and part._step == 21:
            panelLength = part.maxBaseIdx() - part.minBaseIdx() + 1
            # A minimum Honeycomb route still needs two separated legal
            # registers on every path gap: one alone cannot close all
            # scaffold-bearing helices into the same perimeter loop.
            densitySpacing = part._step * max(
                1, panelLength // part._step - 2)
        else:
            densitySpacing = (
                max(requestedSpacing,
                    part.maxBaseIdx() - part.minBaseIdx() + part._step)
                if routeOnly else requestedSpacing)
        initialRecords = part._existingScaffoldCrossoverRecords()
        if not initialRecords and not \
                _hasLegalScaffoldCrossoverCandidate(
                    part, minimumIndex):
            details = {
                'success': False,
                'created': 0,
                'removed': 0,
                'components': len(_scaffoldOligos(part)),
                'is_loop': False,
                'spacing': densitySpacing,
                'requested_spacing': requestedSpacing,
                'main_loop_length': 0,
                'minimum_spacing': None,
                'density_exceptions': 0,
                'density_deficit': 0,
                'hard_density_valid': True,
                'seam_count': 0,
                'seam_pairs': [],
                'evaluated_candidates': 0,
                'message': (
                    '当前 scaffold 区域没有任何合法 crossover 位点；'
                    '未修改设计。')}
            return details if returnDetails else 0
        removed = 0
        beforeAutoIndex = part.undoStack().index()
        sparseFallbackUsed = False

        util.beginSuperMacro(part, desc="AutoCS_scaffolds")
        try:
            if rebuildExisting:
                removedBoundaries = set()
                removableEndpoints = set(
                    (number, index)
                    for fromNumber, toNumber, idx5p, idx3p in initialRecords
                    for number, index in (
                        (fromNumber, idx5p), (toNumber, idx3p)))
                for vh in part.getVirtualHelices():
                    strands = sorted(list(vh.scaffoldStrandSet()),
                                     key=lambda strand: strand.lowIdx())
                    for left, right in zip(strands, strands[1:]):
                        if left.highIdx() + 1 == right.lowIdx() and \
                                ((vh.number(), left.highIdx()) in
                                 removableEndpoints or
                                 (vh.number(), right.lowIdx()) in
                                 removableEndpoints):
                            removedBoundaries.add(
                                (vh.number(), left.highIdx(),
                                 right.lowIdx()))
                for fromNumber, toNumber, idx5p, unused_idx3p in \
                        initialRecords:
                    fromHelix = part.virtualHelix(fromNumber)
                    strand5p = (
                        fromHelix.scaffoldStrandSet().getStrand(idx5p)
                        if fromHelix is not None else None)
                    strand3p = (
                        strand5p.connection3p()
                        if strand5p is not None else None)
                    if strand3p is None or \
                            strand3p.virtualHelix().number() != toNumber:
                        continue
                    part.removeXover(strand5p, strand3p)
                    removed += 1
                _mergeRemovedScaffoldXoverBoundaries(
                    part, removedBoundaries)

            helixRecords = []
            for vh in part.getVirtualHelices():
                if not any(True for unused_strand in
                           vh.scaffoldStrandSet()):
                    continue
                row, column = vh.coord()
                neighbors = [
                    neighbor.number()
                    for neighbor in part.getVirtualHelixNeighbors(vh)
                    if neighbor is not None]
                helixRecords.append(
                    (vh.number(), row, column, neighbors))
            orderedSquareModules = None
            if explicitPaths is not None:
                paths = [list(path) for path in explicitPaths if path]
                if explicitModules:
                    orderedSquareModules = paths
            elif part._step == 32:
                orderedSquareModules = \
                    _autoScaffoldRightPanelModulePaths(
                        part, helixRecords)
                paths = (
                    orderedSquareModules or
                    _autoScaffoldStraightPaths(helixRecords))
            else:
                paths = _autoScaffoldSnakePaths(helixRecords)
            if part._step == 32 and not orderedSquareModules:
                seamPath = paths[0] if paths else []
            else:
                seamPath = paths[-1] if paths else []
            globalSeamPair = (
                tuple(sorted(seamPath[-2:]))
                if len(seamPath) >= 2 else None)

            if orderedSquareModules:
                moduleCreated = [
                    part._autoScaffoldCrossoversPreGuided(
                        densitySpacing=densitySpacing,
                        minimumIndex=minimumIndex,
                        minimumDensity=routeOnly,
                        explicitPaths=[modulePath])
                    for modulePath in orderedSquareModules]
                created = sum(moduleCreated)
            else:
                created = part._autoScaffoldCrossoversPreGuided(
                    densitySpacing=densitySpacing,
                    minimumIndex=minimumIndex,
                    minimumDensity=routeOnly)
            bridgePairs = (
                (_autoScaffoldSequentialModuleBridgePairs(
                    paths, helixRecords)
                 if orderedSquareModules else
                 _autoScaffoldSingleStraightBridgePairs(
                    paths, helixRecords))
                if part._step == 32 else set())
            finalizeSingleLoop = True

            if finalizeSingleLoop and bridgePairs:
                _mergeClosedScaffoldLoops(
                    part, densitySpacing, directionSpacing,
                    minimumIndex=minimumIndex,
                    avoidMultiplesOfEight=avoidMultiplesOfEight,
                    allowedPairs=bridgePairs)
            if finalizeSingleLoop and part._step == 21:
                _mergeClosedScaffoldLoops(
                    part, densitySpacing, directionSpacing,
                    minimumIndex=minimumIndex,
                    avoidMultiplesOfEight=avoidMultiplesOfEight)
            if finalizeSingleLoop:
                _removeAutoCrossoversOutsideScaffoldLoops(part)

            if finalizeSingleLoop and not routeOnly:
                groupByHelix = dict(
                    (number, groupIndex)
                    for groupIndex, path in enumerate(paths)
                    for number in path)
                actualBridgePairs = set(
                    tuple(sorted(record[:2]))
                    for record in
                    part._existingScaffoldCrossoverRecords()
                    if groupByHelix.get(record[0]) !=
                       groupByHelix.get(record[1]))
                routePairs = set(
                    tuple(sorted(pair))
                    for path in paths for pair in zip(path, path[1:]))
                globalSeam = (
                    set([globalSeamPair])
                    if globalSeamPair is not None else set())
                densePairs = routePairs.difference(globalSeam)
                densePairs.update(actualBridgePairs)

                densityCandidates = []
                seen = set()
                for vh in part.getVirtualHelices():
                    strandSet = vh.scaffoldStrandSet()
                    is5to3 = strandSet.isDrawn5to3()
                    for neighbor, index, strandType, isLowIdx in \
                            part.potentialCrossoverList(vh):
                        if strandType != StrandType.Scaffold or \
                                (minimumIndex is not None and
                                 index < minimumIndex):
                            continue
                        fromHelixIs5p = (
                            (isLowIdx and is5to3) or
                            (not isLowIdx and not is5to3))
                        record = (
                            vh.number(), neighbor.number(),
                            index, index)
                        if not fromHelixIs5p or record in seen or \
                                tuple(sorted(record[:2])) not in \
                                densePairs:
                            continue
                        seen.add(record)
                        densityCandidates.append(record)
                _densifyClosedScaffoldLoop(
                    part, densityCandidates, densitySpacing,
                    directionSpacing, avoidMultiplesOfEight)
                _removeAutoCrossoversOutsideScaffoldLoops(part)

                # A modular route can close each subpath independently before
                # joining them, leaving one sparse seam per module.  At the
                # native density, rebalance complete reciprocal blocks so the
                # joined global loop retains only its final seam.
                # Two-path Honeycomb sheets use the established paired-row
                # layout (1.json -> 2.json).  Rebalancing is needed only for
                # the three-or-more-path case that otherwise accumulates a
                # seam on every additional band.
                if densityMultiple == 1 and (
                        (part._step == 21 and len(paths) > 2) or
                        (part._step == 32 and len(paths) > 1)):
                    closedLoops = [
                        oligo for oligo in _scaffoldOligos(part)
                        if oligo.isLoop()]
                    targetLength = max(
                        [oligo.length() for oligo in closedLoops] or [0])
                    targetRecords = _singleSeamScaffoldRecordSet(
                        part, paths, targetLength)
                    if targetRecords is not None:
                        _applyScaffoldCrossoverRecordSet(
                            part, targetRecords)
                        _removeAutoCrossoversOutsideScaffoldLoops(part)
                bridgeModules = (
                    _autoScaffoldGeneralModulePaths(
                        part, helixRecords) or paths)
                _regularizeScaffoldBridgePairs(
                    part, bridgeModules, densitySpacing,
                    directionSpacing, avoidMultiplesOfEight)
                _removeAutoCrossoversOutsideScaffoldLoops(part)
        finally:
            util.endSuperMacro(part)

        oligos = _scaffoldOligos(part)
        loops = [oligo for oligo in oligos if oligo.isLoop()]

        def sparseRouteCoversAllHelices(unused_closedLoops):
            return _autoScaffoldGlobalLoopQuality(
                part, requestedSpacing, paths).get('valid', False)

        if not routeOnly and densityMultiple > 1 and \
                not sparseRouteCoversAllHelices(loops):
            # The requested sparse route can occasionally have no legal
            # closing phase at exactly N lattice periods.  Restore the input,
            # build the proven native-period loop, then remove reciprocal
            # blocks only when the same loop survives.  This makes the
            # requested spacing a preference while the native period remains
            # the hard maximum-density limit.
            while part.undoStack().index() > beforeAutoIndex and \
                    part.undoStack().canUndo():
                part.undoStack().undo()
            util.beginSuperMacro(
                part, desc="AutoCS_scaffolds sparse fallback")
            try:
                denseDetails = part.autoScaffoldCrossovers(
                    minimumIndex=minimumIndex, densityMultiple=1,
                    routeOnly=False, rebuildExisting=rebuildExisting,
                    returnDetails=True, _unifiedInternal=True,
                    _explicitPaths=paths,
                    _explicitModules=bool(orderedSquareModules))
                denseLoops = [
                    oligo for oligo in _scaffoldOligos(part)
                    if oligo.isLoop()]
                if sparseRouteCoversAllHelices(denseLoops):
                    _thinClosedScaffoldLoopToSpacing(
                        part, requestedSpacing, part._step,
                        protectedRecords=(
                            initialRecords if not rebuildExisting else ()))
                    _removeAutoCrossoversOutsideScaffoldLoops(
                        part, protectedRecords=(
                            initialRecords if not rebuildExisting else ()))
                    sparseFallbackUsed = True
                    removed = denseDetails.get('removed', removed)
            finally:
                util.endSuperMacro(part)
            oligos = _scaffoldOligos(part)
            loops = [oligo for oligo in oligos if oligo.isLoop()]

        globalQuality = _autoScaffoldGlobalLoopQuality(
            part, requestedSpacing, paths)
        disableGlobalRecovery = bool(
            unused_options.get('_disableGlobalRecovery'))
        if not routeOnly and densityMultiple == 1 and \
                not globalQuality.get('valid') and \
                not disableGlobalRecovery:
            # Restore the pre-rollback deterministic global candidate only as
            # a rescue path.  A valid current module layout is never replaced.
            while part.undoStack().index() > beforeAutoIndex and \
                    part.undoStack().canUndo():
                part.undoStack().undo()
            recovered = _autoScaffoldDenseRegularFast(
                part, minimumIndex, directionSpacing,
                avoidMultiplesOfEight, rebuildExisting)
            globalRecoveryUsed = False
            if recovered is not None:
                recoveredPaths = _autoScaffoldSnakePaths(
                    helixRecords)
                recoveredSeamPair = (
                    tuple(sorted(recoveredPaths[-1][-2:]))
                    if recoveredPaths and
                    len(recoveredPaths[-1]) >= 2 else None)
                recoveredQuality = _autoScaffoldGlobalLoopQuality(
                    part, part._step, recoveredPaths)
                if recoveredQuality.get('valid'):
                    globalRecoveryUsed = True
                    paths = recoveredPaths
                    globalSeamPair = recoveredSeamPair
                    globalQuality = recoveredQuality
                    oligos = _scaffoldOligos(part)
                    loops = [
                        oligo for oligo in oligos if oligo.isLoop()]
                    removed = recovered.get('removed', removed)
            # Keep behavior deterministic even when neither candidate can
            # satisfy the hard global test: rebuild the current route once
            # without re-entering this recovery branch.
            if not globalRecoveryUsed:
                while part.undoStack().index() > beforeAutoIndex and \
                        part.undoStack().canUndo():
                    part.undoStack().undo()
                return part.autoScaffoldCrossovers(
                    minimumIndex=minimumIndex,
                    densityMultiple=densityMultiple,
                    routeOnly=routeOnly,
                    rebuildExisting=rebuildExisting,
                    returnDetails=returnDetails,
                    _disableGlobalRecovery=True,
                    _unifiedInternal=unifiedInternal,
                    _densitySearchDisabled=densitySearchDisabled,
                    _explicitPaths=paths,
                    _explicitModules=bool(orderedSquareModules))

        finalRecords = part._existingScaffoldCrossoverRecords()
        created = (len(finalRecords) if rebuildExisting else
                   max(0, len(finalRecords) - len(initialRecords)))
        minimumSpacing = _minimumDirectedScaffoldSpacing(part)
        densityDeficit = (
            max(0, requestedSpacing - minimumSpacing)
            if sparseFallbackUsed and minimumSpacing is not None else 0)
        globalQuality = _autoScaffoldGlobalLoopQuality(
            part, requestedSpacing, paths)
        details = {
            'success': (True if routeOnly else
                        globalQuality.get('valid', False)),
            'created': created,
            'removed': removed,
            'components': len(oligos),
            'is_loop': len(oligos) == 1 and bool(loops),
            'spacing': densitySpacing,
            'requested_spacing': requestedSpacing,
            'main_loop_length':
                max([oligo.length() for oligo in loops] or [0]),
            'minimum_spacing': minimumSpacing,
            'density_exceptions': 0,
            'density_deficit': densityDeficit,
            'extension_bases': 0,
            'target_helices': globalQuality.get('target_helices', []),
            'covered_length': globalQuality.get('covered_length', 0),
            'minimum_covered_length':
                globalQuality.get('minimum_covered_length', 0),
            'hard_density_valid':
                globalQuality.get('hard_density_valid', False),
            'seam_count': globalQuality.get('seam_count', 0),
            'seam_pairs': globalQuality.get('seam_pairs', []),
            'seam_position': globalQuality.get('seam_position'),
            'message': (
                '已按确定性路线并执行单闭环收尾，%s生成 '
                '%d 个 scaffold crossover（当前共 %d 个）。' %
                ('采用最低密度' if routeOnly else
                 ('优先采用 1/%d bp；精确周期不能闭环时已在 '
                  '1/%d bp 硬上限内回退' %
                  (requestedSpacing, part._step))
                 if sparseFallbackUsed else
                 '以 1/%d bp 密度上限' % requestedSpacing,
                 created, len(finalRecords)))}

        baselineScore, baselineMetrics = \
            _autoScaffoldUnifiedCandidateScore(
                part, details, paths, helixRecords, requestedSpacing)
        details.update(baselineMetrics)
        canonicalModules = (
            _autoScaffoldGeneralModulePaths(part, helixRecords) or paths)
        details['route_modules'] = [
            list(path) for path in canonicalModules]
        details['route_paths'] = [list(path) for path in paths]

        def tryAlternativeDensity():
            if unifiedInternal or densitySearchDisabled:
                return None
            hasScaffoldCandidate = \
                _hasLegalScaffoldCrossoverCandidate(
                    part, minimumIndex)
            if not hasScaffoldCandidate:
                return None
            panelPeriods = max(
                1, (part.maxBaseIdx() - part.minBaseIdx() + 1) //
                part._step)
            maximumMultiple = min(
                max(3, densityMultiple + 2), panelPeriods)
            alternatives = [
                multiplier for multiplier in
                range(1, maximumMultiple + 1)
                if multiplier != densityMultiple]
            alternatives.sort(key=lambda multiplier: (
                abs(multiplier - densityMultiple), multiplier))
            evaluatedDensity = []
            for alternative in alternatives:
                while part.undoStack().index() > beforeAutoIndex and \
                        part.undoStack().canUndo():
                    part.undoStack().undo()
                alternativeDetails = part.autoScaffoldCrossovers(
                    minimumIndex=minimumIndex,
                    densityMultiple=alternative,
                    routeOnly=False,
                    rebuildExisting=rebuildExisting,
                    returnDetails=True,
                    _densitySearchDisabled=True)
                if alternativeDetails.get('success'):
                    densityExceptions = \
                        _scaffoldDensityExceptionCount(
                            part, requestedSpacing)
                    densityDeficit = _scaffoldDensityDeficit(
                        part, requestedSpacing)
                    structuralTail = tuple(
                        alternativeDetails.get(
                            'candidate_score', (0,) * 10))[4:]
                    evaluatedDensity.append((
                        (densityExceptions, densityDeficit,
                         abs(alternative - densityMultiple),
                         alternative) + structuralTail,
                        alternative))
            while part.undoStack().index() > beforeAutoIndex and \
                    part.undoStack().canUndo():
                part.undoStack().undo()
            if not evaluatedDensity:
                return None
            winningMultiple = min(evaluatedDensity)[1]
            alternativeDetails = part.autoScaffoldCrossovers(
                minimumIndex=minimumIndex,
                densityMultiple=winningMultiple,
                routeOnly=False,
                rebuildExisting=rebuildExisting,
                returnDetails=True,
                _densitySearchDisabled=True)
            alternativeDetails['requested_spacing'] = requestedSpacing
            alternativeDetails['density_fallback_multiple'] = \
                winningMultiple
            alternativeDetails['message'] = (
                '首选 1/%d bp 路线无法通过唯一闭环硬准则；'
                '已选择距离最近且合法的 1/%d bp 闭环。' %
                (requestedSpacing, part._step * winningMultiple) +
                alternativeDetails.get('message', ''))
            return alternativeDetails

        if not unifiedInternal and not routeOnly:
            candidatePaths = []
            baselineStructural = (
                baselineMetrics.get('module_fragmentation', 0),
                baselineMetrics.get('module_inversions', 0),
                baselineMetrics.get('longitudinal_bridge_count', 0),
                abs(int(details.get('seam_count', 0)) - 1),
                -(details.get('seam_position')
                  if details.get('seam_position') is not None else -1),
                baselineMetrics.get(
                    'module_orientation_reversals', 0))
            for candidate in (
                    _autoScaffoldGeneralModulePaths(part, helixRecords),
                    _autoScaffoldImportedRunPaths(part, helixRecords),
                    _autoScaffoldSnakePaths(helixRecords)):
                if not candidate or candidate == paths or \
                        candidate in candidatePaths:
                    continue
                fragmentation, inversions, orientation = \
                    _autoScaffoldModuleOrderPenalty(
                        part, candidate, helixRecords)
                minimumBridgeCount = max(0, len(candidate) - 1)
                candidateSeamPosition = max(
                    -1, sum(max(0, len(path) - 1)
                            for path in candidate) - 1)
                candidateStructural = (
                    fragmentation, inversions, minimumBridgeCount, 0,
                    -candidateSeamPosition, orientation)
                # If the current valid route is already at least as good on
                # every remaining structural objective, evaluating another
                # full DNA topology only adds latency and cannot win the
                # agreed lexicographic ranking.
                if details.get('success') and \
                        int(details.get('seam_count', 0)) <= 1 and \
                        candidateStructural >= baselineStructural:
                    continue
                candidatePaths.append(candidate)

            evaluated = [(
                baselineScore + (0,), None, False, baselineMetrics)]
            for candidateIndex, candidate in enumerate(candidatePaths, 1):
                while part.undoStack().index() > beforeAutoIndex and \
                        part.undoStack().canUndo():
                    part.undoStack().undo()
                candidateDetails = part.autoScaffoldCrossovers(
                    minimumIndex=minimumIndex,
                    densityMultiple=densityMultiple,
                    routeOnly=False,
                    rebuildExisting=rebuildExisting,
                    returnDetails=True,
                    _unifiedInternal=True,
                    _explicitPaths=candidate,
                    _explicitModules=True)
                candidateScore, candidateMetrics = \
                    _autoScaffoldUnifiedCandidateScore(
                        part, candidateDetails, candidate,
                        helixRecords, requestedSpacing)
                evaluated.append((
                    candidateScore + (candidateIndex,),
                    candidate, True, candidateMetrics))

            if candidatePaths:
                winning = min(evaluated, key=lambda item: item[0])
                while part.undoStack().index() > beforeAutoIndex and \
                        part.undoStack().canUndo():
                    part.undoStack().undo()
                winningPaths = winning[1]
                finalDetails = part.autoScaffoldCrossovers(
                    minimumIndex=minimumIndex,
                    densityMultiple=densityMultiple,
                    routeOnly=False,
                    rebuildExisting=rebuildExisting,
                    returnDetails=True,
                    _unifiedInternal=True,
                    _explicitPaths=winningPaths,
                    _explicitModules=winning[2])
                finalScore, finalMetrics = \
                    _autoScaffoldUnifiedCandidateScore(
                        part, finalDetails,
                        winningPaths if winningPaths is not None else paths,
                        helixRecords, requestedSpacing)
                finalDetails.update(finalMetrics)
                finalDetails['candidate_score'] = finalScore
                finalPaths = (
                    winningPaths if winningPaths is not None else paths)
                finalModules = (
                    _autoScaffoldGeneralModulePaths(
                        part, helixRecords) or finalPaths)
                finalDetails['route_modules'] = [
                    list(path) for path in finalModules]
                finalDetails['route_paths'] = [
                    list(path) for path in finalPaths]
                finalDetails['evaluated_candidates'] = len(evaluated)
                finalDetails['message'] = (
                    '已统一按“闭环、密度、Path 顺序/模块、纵向桥、'
                    'seam”排序 %d 个候选。' % len(evaluated) +
                    finalDetails.get('message', ''))
                if not finalDetails.get('success'):
                    densityAlternative = tryAlternativeDensity()
                    if densityAlternative is not None:
                        return (densityAlternative if returnDetails else
                                densityAlternative.get('created', 0))
                    while part.undoStack().index() > beforeAutoIndex and \
                            part.undoStack().canUndo():
                        part.undoStack().undo()
                    finalDetails.update({
                        'created': 0,
                        'removed': 0,
                        'components': len(_scaffoldOligos(part)),
                        'main_loop_length': 0,
                    })
                    finalDetails['message'] = (
                        '所有候选均未通过唯一主闭环硬准则，'
                        '已恢复运行前设计；未保留局部环或临时 '
                        'crossover。')
                return finalDetails if returnDetails else \
                    finalDetails.get('created', 0)

        details['evaluated_candidates'] = 1
        if not unifiedInternal and not details.get('success'):
            densityAlternative = tryAlternativeDensity()
            if densityAlternative is not None:
                return (densityAlternative if returnDetails else
                        densityAlternative.get('created', 0))
            while part.undoStack().index() > beforeAutoIndex and \
                    part.undoStack().canUndo():
                part.undoStack().undo()
            details.update({
                'created': 0,
                'removed': 0,
                'components': len(_scaffoldOligos(part)),
                'main_loop_length': 0,
            })
            details['message'] = (
                '当前路线未通过唯一主闭环硬准则，已恢复运行前设计；'
                '未保留局部环或临时 crossover。')
        return details if returnDetails else details.get('created', created)

    def autoStaple(part, preservePeriodicCrossovers=False,
                   allowedCrossoverPairs=None):
        """Autostaple does the following:
        1. Clear existing staple strands by iterating over each strand
        and calling RemoveStrandCommand on each. The next strand to remove
        is always at index 0.
        2. Create temporary strands that span regions where scaffold is present.
        3. Determine where actual strands will go based on strand overlap with
        prexovers.
        4. Delete temporary strands and create new strands.
        """
        # A cross-lattice staple is intentionally outside either lattice's
        # local AutoCS rule set. Preserve it instead of allowing AutoStaple's
        # clear-and-rebuild phase to delete strands from both panels.
        if any(strand.oligo().isHybrid()
               for vh in part.getVirtualHelices()
               for strand in vh.stapleStrandSet()):
            return False
        epDict = {}  # keyed on StrandSet
        cmds = []
        allowedCrossoverPairs = (
            set(tuple(sorted(pair)) for pair in allowedCrossoverPairs)
            if allowedCrossoverPairs is not None else None)

        def hasIndelAt(virtualHelix, index):
            return index in part.insertions().get(
                virtualHelix.coord(), {})

        scaffoldPositionsByPairSide = defaultdict(set)
        scaffoldPositionsByHelix = defaultdict(set)
        scaffoldPositionsByPair = defaultdict(set)
        circularDensitySize = (
            part.maxBaseIdx() + 1 if preservePeriodicCrossovers else None)
        if part._step in (21, 32):
            for fromHelix, toHelix, fromIndex, toIndex in \
                    part._existingScaffoldCrossoverRecords():
                scaffoldPositionsByPairSide[
                    (fromHelix, "outgoing", toHelix)].add(fromIndex)
                scaffoldPositionsByPairSide[
                    (toHelix, "incoming", fromHelix)].add(toIndex)
                scaffoldPositionsByHelix[fromHelix].add(fromIndex)
                scaffoldPositionsByHelix[toHelix].add(toIndex)
                pair = tuple(sorted((fromHelix, toHelix)))
                scaffoldPositionsByPair[pair].update(
                    (fromIndex, toIndex))

        # clear existing staple strands
        # part.verifyOligos()

        for o in list(part.oligos()):
            if not o.isStaple():
                continue
            c = Oligo.RemoveOligoCommand(o)
            cmds.append(c)
        # end for
        util.execCommandList(part, cmds, desc="Clear staples")
        cmds = []

        # Normal designs follow all scaffold bases. Guided Design uses its
        # saved image silhouette and omits scaffold-only edge extensions.
        for vh in part.getVirtualHelices():
            segments = part.stapleTargetSegments(vh)
            stapSS = vh.stapleStrandSet()
            epDict[stapSS] = []
            for i in range(len(segments)):
                lo, hi = segments[i]
                epDict[stapSS].extend(segments[i])
                c = StrandSet.CreateStrandCommand(
                    stapSS, lo, hi, i,
                    color=styles.AUTOMATIC_STAP_COLOR)
                cmds.append(c)
        util.execCommandList(part, cmds, desc="Add tmp strands", useUndoStack=False)
        cmds = []

        # determine where xovers should be installed
        for vh in part.getVirtualHelices():
            stapSS = vh.stapleStrandSet()
            scafSS = vh.scaffoldStrandSet()
            is5to3 = stapSS.isDrawn5to3()
            potentialXovers = part.potentialCrossoverList(vh)
            for neighborVh, idx, strandType, isLowIdx in potentialXovers:
                if strandType != StrandType.Staple:
                    continue
                if allowedCrossoverPairs is not None and tuple(sorted(
                        (vh.number(), neighborVh.number()))) not in \
                        allowedCrossoverPairs:
                    continue
                if isLowIdx and is5to3:
                    # Curved periodic topology has priority over provisional
                    # indels.  Its final rebalance pass moves every indel away
                    # from the completed scaffold/staple crossover set.
                    if not preservePeriodicCrossovers and \
                            (hasIndelAt(vh, idx) or
                             hasIndelAt(neighborVh, idx)):
                        continue
                    if _stapleCrossoverIsOverdense(
                            part, vh.number(), neighborVh.number(), idx,
                            scaffoldPositionsByPairSide,
                            circularDensitySize,
                            scaffoldPositionsByHelix,
                            scaffoldPositionsByPair):
                        continue
                    strand = stapSS.getStrand(idx)
                    neighborSS = neighborVh.stapleStrandSet()
                    nStrand = neighborSS.getStrand(idx)
                    if strand == None or nStrand == None:
                        continue
                    # check for bases on both strands at [idx-1:idx+3]
                    if not (strand.lowIdx() < idx and strand.highIdx() > idx + 1):
                        continue
                    if not (nStrand.lowIdx() < idx and nStrand.highIdx() > idx + 1):
                        continue

                    # disable edge xovers
                    scafStrandL1 = scafSS.getStrand(idx-1)
                    scafStrandM = scafSS.getStrand(idx)
                    scafStrandH1 = scafSS.getStrand(idx+1)
                    if scafStrandL1:
                        if scafStrandL1.hasXoverAt(idx-1) and not vh.hasStrandAtIdx(idx-2):
                            continue
                        if scafStrandL1.hasXoverAt(idx-2) and not vh.hasStrandAtIdx(idx-3):
                            continue
                    if scafStrandM:
                        if scafStrandM.hasXoverAt(idx-1) and not vh.hasStrandAtIdx(idx-2):
                            continue
                        if scafStrandM.hasXoverAt(idx+1) and not vh.hasStrandAtIdx(idx+2):
                            continue
                    if scafStrandH1:
                        if scafStrandH1.hasXoverAt(idx+1) and not vh.hasStrandAtIdx(idx+2):
                            continue
                        if scafStrandH1.hasXoverAt(idx+2) and not vh.hasStrandAtIdx(idx+3):
                            continue

                    # Finally, add the xovers to install
                    epDict[stapSS].extend([idx, idx+1])
                    epDict[neighborSS].extend([idx, idx+1])

        # clear temporary staple strands
        for vh in part.getVirtualHelices():
            stapSS = vh.stapleStrandSet()
            for strand in stapSS:
                c = StrandSet.RemoveStrandCommand(stapSS, strand, 0)
                cmds.append(c)
        util.execCommandList(part, cmds, desc="Rm tmp strands", useUndoStack=False)
        cmds = []

        util.beginSuperMacro(part, desc="AutoCS_staples")

        for stapSS, epList in epDict.items():
            assert (len(epList) % 2 == 0)
            epList = sorted(epList)
            ssIdx = 0
            for i in range(0, len(epList),2):
                lo, hi = epList[i:i+2]
                c = StrandSet.CreateStrandCommand(
                    stapSS, lo, hi, ssIdx,
                    color=styles.AUTOMATIC_STAP_COLOR)
                cmds.append(c)
                ssIdx += 1
        util.execCommandList(part, cmds, desc="Create strands")
        cmds = []

        # create crossovers wherever possible (from strand5p only)
        for vh in part.getVirtualHelices():
            stapSS = vh.stapleStrandSet()
            is5to3 = stapSS.isDrawn5to3()
            potentialXovers = part.potentialCrossoverList(vh)
            for neighborVh, idx, strandType, isLowIdx in potentialXovers:
                if strandType != StrandType.Staple:
                    continue
                if allowedCrossoverPairs is not None and tuple(sorted(
                        (vh.number(), neighborVh.number()))) not in \
                        allowedCrossoverPairs:
                    continue
                if (isLowIdx and is5to3) or (not isLowIdx and not is5to3):
                    if not preservePeriodicCrossovers and \
                            (hasIndelAt(vh, idx) or
                             hasIndelAt(neighborVh, idx)):
                        continue
                    if _stapleCrossoverIsOverdense(
                            part, vh.number(), neighborVh.number(), idx,
                            scaffoldPositionsByPairSide,
                            circularDensitySize,
                            scaffoldPositionsByHelix,
                            scaffoldPositionsByPair):
                        continue
                    strand = stapSS.getStrand(idx)
                    neighborSS = neighborVh.stapleStrandSet()
                    nStrand = neighborSS.getStrand(idx)
                    if strand == None or nStrand == None:
                        continue
                    if idx in strand.idxs() and idx in nStrand.idxs():
                        # only install xovers on pre-split strands
                        part.createXover(strand, idx, nStrand, idx, updateOligo=False)

        c = Part.RefreshOligosCommand(part)
        cmds.append(c)
        util.execCommandList(part, cmds, desc="Assign oligos")

        for oligo in list(part.oligos()):
            if oligo.isStaple() and not oligo.isHybrid():
                oligo.applyColor(styles.AUTOMATIC_STAP_COLOR)
        util.execCommandList(
            part, [_AutobreakAppliedCommand(part, False)],
            desc="Reset Autobreak state")

        cmds = []
        util.endSuperMacro(part)
        return True

    # end def

    def autoBreakStaples(part, preserveCrossovers=False,
                         markUnbreakable=False,
                         preferDeletionDense=False):
        """Break staples using the active lattice's spacing and phase rules.

        Curved Design sets ``preserveCrossovers`` because its crossover set
        has already been selected from native lattice sites; an unbreakable
        staple is reported in red instead of deleting or moving topology.
        """
        if part._step == 32:
            latticeName = 'square'
            minimumIndexDistance = 7
            continuousMinimum = 16
            preferredPhase = 8
            softMinimum = 32
        elif part._step == 21:
            latticeName = 'honeycomb'
            minimumIndexDistance = 6
            continuousMinimum = 14
            preferredPhase = 7
            softMinimum = 28
        else:
            return {'supported': False, 'square': False,
                    'lattice': 'unknown', 'nicks': 0,
                    'removed_xovers': 0, 'skipped': 0,
                    'protected_nicks': 0, 'already_applied': False}

        if part._autobreakStaplesApplied:
            return {'supported': True, 'square': latticeName == 'square',
                    'lattice': latticeName, 'nicks': 0,
                    'removed_xovers': 0, 'skipped': 0,
                    'protected_nicks': len(
                        _existingStapleNickBoundaries(part)),
                    'already_applied': True}

        createdNicks = 0
        removedXovers = 0
        skipped = 0
        toleratedLong = 0
        relaxedDenseMinimum = 0
        unbreakableBases = set()
        hasStaples = False
        protectedNicks = _existingStapleNickBoundaries(part)
        util.beginSuperMacro(part, desc="Autobreak staples")
        try:
            # For both lattices, remove a crossover only 1--4 nt from a
            # physical staple edge and restore both native same-helix joins.
            # A nick that existed before Autobreak is a fixed user boundary:
            # never merge across it while cleaning up an edge crossover.
            if not preserveCrossovers:
                removedXovers += _removeShortEdgeStapleXovers(
                                                part, protectedNicks)

            stapleXovers = _crossoverPositionsByHelix(
                                                part, StrandType.Staple)
            scaffoldXovers = _crossoverPositionsByHelix(
                                                part, StrandType.Scaffold)
            plans = []
            stapleOligos = sorted(
                (oligo for oligo in part.oligos()
                 if oligo.isStaple() and not oligo.isHybrid()),
                key=lambda oligo: (oligo.strand5p().virtualHelix().number(),
                                   oligo.strand5p().idx5Prime(),
                                   oligo.actualLength()))
            hasStaples = bool(stapleOligos)
            for oligo in stapleOligos:
                deletionCount = 0
                if preferDeletionDense:
                    deletionCount = sum(
                        int(record[3]) == 0
                        for record in _stapleOligoBaseRecords(oligo))
                # Normal products are 21--57 nt and prefer 30--50 nt.  A
                # deletion-dense Curved/Frame component has already lost
                # local contiguous hybridization, so its first pass instead
                # prefers 40--60 nt (target 50 nt).  Only if this complete
                # first-pass problem has no legal solution may a second pass
                # use the exceptional 64-nt ceiling.  Crossover clearance and
                # all indel/nick exclusions remain unchanged in both passes.
                deletionDense = deletionCount >= 2
                # Preferences remain segment-local.  Enabling the dense mode
                # only lets candidate final staples which themselves contain
                # at least two deletions use 40--60 nt; ordinary segments from
                # the same source oligo retain 30--50 nt and the 57-nt normal
                # ceiling.
                planArguments = dict(
                    preferredMinimum=30,
                    preferredMaximum=50,
                    targetLength=40,
                    terminalMaximum=49,
                    preferDeletionDense=deletionDense,
                    requireDeletionDenseMinimum=deletionDense)
                plan = _bestStapleBreakPlan(
                    oligo, stapleXovers, scaffoldXovers,
                    minimumIndexDistance, continuousMinimum, preferredPhase,
                    softMinimum,
                    hardMaximum=(60 if deletionDense else 57),
                    **planArguments)
                if plan is None:
                    plan = _bestStapleBreakPlan(
                        oligo, stapleXovers, scaffoldXovers,
                        minimumIndexDistance, continuousMinimum,
                        preferredPhase,
                        softMinimum,
                        hardMaximum=64, **planArguments)
                    if plan is not None:
                        # Count oligos which genuinely required the exceptional
                        # 64-nt ceiling; the ordinary/deletion-dense pass was
                        # proven unsatisfiable immediately above.
                        toleratedLong += 1
                # If immutable physical edges/crossovers make the strict
                # 40-nt dense minimum impossible, retain the best relaxed
                # native partition rather than changing topology.  Curved
                # and Frame then heal neighbouring same-helix nicks and rerun
                # this strict problem locally before accepting the exception.
                if plan is None and deletionDense:
                    relaxedArguments = dict(planArguments)
                    relaxedArguments["requireDeletionDenseMinimum"] = False
                    plan = _bestStapleBreakPlan(
                        oligo, stapleXovers, scaffoldXovers,
                        minimumIndexDistance, continuousMinimum,
                        preferredPhase, softMinimum,
                        hardMaximum=60, **relaxedArguments)
                    if plan is None:
                        plan = _bestStapleBreakPlan(
                            oligo, stapleXovers, scaffoldXovers,
                            minimumIndexDistance, continuousMinimum,
                            preferredPhase, softMinimum,
                            hardMaximum=64, **relaxedArguments)
                        if plan is not None:
                            toleratedLong += 1
                    if plan is not None:
                        relaxedDenseMinimum += 1
                if plan is None:
                    # An already-open 58--64 nt terminal product may have no
                    # legal internal nick because of fixed edge/crossover or
                    # indel exclusions.  Preserve it only after both complete
                    # partition passes failed.  A closed product still needs
                    # one best legal opening nick to remain synthesizable.
                    if 58 <= oligo.actualLength() <= 64:
                        if oligo.isLoop():
                            records = _stapleOligoBaseRecords(oligo)
                            toleratedCandidates = \
                                _legalStapleNickBoundaries(
                                    records, stapleXovers, scaffoldXovers,
                                    minimumIndexDistance, preferredPhase)
                            if toleratedCandidates:
                                unused_offset, candidate = max(
                                    toleratedCandidates.items(),
                                    key=lambda item: (
                                        int(item[1][3]), int(item[1][2]),
                                        item[1][4], -item[1][1],
                                        -item[1][0]))
                                plans.append((candidate[0], candidate[1]))
                                toleratedLong += 1
                                continue
                        else:
                            toleratedLong += 1
                            continue
                    skipped += 1
                    warningRecords = _stapleOligoBaseRecords(oligo)
                    unbreakableBases.update(
                        (helix, index) for helix, index,
                        unused_strand, unused_length in warningRecords)
                    # A closed warning oligo has no serializable 5' color
                    # anchor.  Open it at the best otherwise-legal nick so
                    # the user gets one visible red chain for manual repair;
                    # do not move or remove any crossover.
                    if markUnbreakable and oligo.isLoop():
                        warningCandidates = _legalStapleNickBoundaries(
                            warningRecords, stapleXovers, scaffoldXovers,
                            minimumIndexDistance, preferredPhase)
                        if warningCandidates:
                            unused_offset, candidate = max(
                                warningCandidates.items(),
                                key=lambda item: (
                                    int(item[1][3]), int(item[1][2]),
                                    item[1][4], -item[1][1],
                                    -item[1][0]))
                            plans.append((candidate[0], candidate[1]))
                if plan is not None:
                    plans.extend(plan)

            # Coordinates remain stable while strands/oligos are split, so
            # plans can be applied deterministically after all optimization.
            for helixNumber, upperIndex in sorted(set(plans)):
                vh = part.virtualHelix(helixNumber)
                if vh is None:
                    continue
                strandSet = vh.stapleStrandSet()
                lowerIndex = upperIndex - 1
                strand = strandSet.getStrand(lowerIndex)
                if strand is None or strand is not \
                        strandSet.getStrand(upperIndex):
                    continue
                splitIndex = (lowerIndex if strandSet.isDrawn5to3()
                              else upperIndex)
                if strandSet.splitStrand(strand, splitIndex):
                    createdNicks += 1
            # Keep all Autobreak products visually identical to Add staples
            # and AutoCS staples, even when the input was manually colored.
            if hasStaples:
                for oligo in list(part.oligos()):
                    if oligo.isStaple() and not oligo.isHybrid():
                        oligo.applyColor(styles.AUTOMATIC_STAP_COLOR)
                if markUnbreakable:
                    # Splitting other staple oligos can rebuild the oligo
                    # registry, so the pre-plan Oligo objects are not stable
                    # identifiers.  Re-identify warnings by occupied bases.
                    for oligo in list(part.oligos()):
                        if not oligo.isStaple() or oligo.isHybrid():
                            continue
                        if any((helix, index) in unbreakableBases
                               for helix, index, unused_strand,
                               unused_length
                               in _stapleOligoBaseRecords(oligo)):
                            oligo.applyColor('#cc0000')
                util.execCommandList(
                    part, [_AutobreakAppliedCommand(part, True)],
                    desc="Mark Autobreak applied")
        finally:
            util.endSuperMacro(part)
        # QUndoStack.push() normally runs redo synchronously.  Keep the
        # runtime guard explicit as well so a second button click can never
        # enter the optimizer again before the view/event queue catches up.
        if hasStaples:
            part._autobreakStaplesApplied = True
        return {'supported': True, 'square': latticeName == 'square',
                'lattice': latticeName, 'nicks': createdNicks,
                'removed_xovers': removedXovers, 'skipped': skipped,
                'tolerated_long_staples': toleratedLong,
                'relaxed_deletion_dense_minimum': relaxedDenseMinimum,
                'protected_nicks': len(protectedNicks),
                'already_applied': False}

    # end def

    def verifyOligoStrandCounts(self):
        total_stap_strands = 0
        stapOligos = set()
        total_stap_oligos = 0

        for vh in self.getVirtualHelices():
            stapSS = vh.stapleStrandSet()
            total_stap_strands += len(stapSS._strandList)
            for strand in stapSS:
                stapOligos.add(strand.oligo())
        # print "# stap oligos:", len(stapOligos), "# stap strands:", total_stap_strands


    def verifyOligos(self):
        total_errors = 0
        total_passed = 0

        for o in list(self.oligos()):
            oL = o.length()
            a = 0
            gen = o.strand5p().generator3pStrand()

            for s in gen:
                a += s.totalLength()
            # end for
            if oL != a:
                total_errors += 1
                # print "wtf", total_errors, "oligoL", oL, "strandsL", a, "isStaple?", o.isStaple()
                o.applyColor('#ff0000')
            else:
                total_passed += 1
        # end for
        # print "Total Passed: ", total_passed, "/", total_passed+total_errors
    # end def

    def removeVirtualHelices(self, useUndoStack=True):
        vhs = [vh for vh in self._coordToVirtualHelix.values()]
        for vh in vhs:
            vh.remove(useUndoStack)
        # end for
    # end def

    # def remove(self, useUndoStack=True):
    #     """
    #     This method uses the slow method of removing each element one at a time
    #     while maintaining state while the command is executed
    #     """
    #     self.partHideSignal.emit(self)
    #     self._activeVirtualHelix = None
    #     if useUndoStack:
    #         self.undoStack().beginMacro("Delete Part")
    #     self.removeVirtualHelices(useUndoStack)
    #     c = Part.RemovePartCommand(self)
    #     if useUndoStack:
    #         self.undoStack().push(c)
    #         self.undoStack().endMacro()
    #     else:
    #         c.redo()
    # # end def

    def remove(self, useUndoStack=True):
        """
        This method assumes all strands are and all VirtualHelices are
        going away, so it does not maintain a valid model state while
        the command is being executed.
        Everything just gets pushed onto the undostack more or less as is.
        Except that strandSets are actually cleared then restored, but this
        is neglible performance wise.  Also, decorators/insertions are assumed
        to be parented to strands in the view so their removal Signal is
        not emitted.  This causes problems with undo and redo down the road
        but works as of now.
        """
        self.partHideSignal.emit(self)
        self._activeVirtualHelix = None
        if useUndoStack:
            self.undoStack().beginMacro("Delete Part")
        # remove strands and oligos
        self.removeAllOligos(useUndoStack)
        # remove VHs
        vhs = list(self._coordToVirtualHelix.values())
        for vh in vhs:
            d = VirtualHelix.RemoveVirtualHelixCommand(self, vh)
            if useUndoStack:
                self.undoStack().push(d)
            else:
                d.redo()
        # end for
        # remove the part
        e = Part.RemovePartCommand(self)
        if useUndoStack:
            self.undoStack().push(e)
            self.undoStack().endMacro()
        else:
            e.redo()
    # end def
    
    def removeAllOligos(self, useUndoStack=True):
        # clear existing oligos
        cmds = []
        for o in list(self.oligos()):
            cmds.append(Oligo.RemoveOligoCommand(o))
        # end for
        util.execCommandList(self, cmds, desc="Clear oligos", useUndoStack=useUndoStack)
    # end def

    def addOligo(self, oligo):
        self._oligos.add(oligo)

    # end def

    def createVirtualHelix(self, row, col, useUndoStack=True):
        c = Part.CreateVirtualHelixCommand(self, row, col)
        util.execCommandList(self, [c], desc="Add VirtualHelix", \
                                                useUndoStack=useUndoStack)
    # end def

    def createXover(self, strand5p, idx5p, strand3p, idx3p, updateOligo=True, useUndoStack=True):
        # prexoveritem needs to store left or right, and determine
        # locally whether it is from or to
        # pass that info in here in and then do the breaks
        ss5p = strand5p.strandSet()
        ss3p = strand3p.strandSet()
        if ss5p.strandType() != ss3p.strandType():
            return
        if useUndoStack:
            self.undoStack().beginMacro("Create Xover")
        if ss5p.isScaffold() and useUndoStack:  # ignore on import
            strand5p.oligo().applySequence(None)
            strand3p.oligo().applySequence(None)
        if strand5p == strand3p:
            """
            This is a complicated case basically we need a truth table.
            1 strand becomes 1, 2 or 3 strands depending on where the xover is
            to.  1 and 2 strands happen when the xover is to 1 or more existing
            endpoints.  Since SplitCommand depends on a StrandSet index, we need
            to adjust this strandset index depending which direction the crossover is
            going in.

            Below describes the 3 strand process
            1) Lookup the strands strandset index (ssIdx)
            1) Split attempted on the 3 prime strand, AKA 5prime endpoint of
            one of the new strands.  We have now created 2 strands, and the ssIdx
            is either the same as the first lookup, or one more than it depending
            on which way the the strand is drawn (isDrawn5to3).  If a split occured
            the 5prime strand is definitely part of the 3prime strand created in this step
            2) Split is attempted on the resulting 2 strands.  There is
            now 3 strands, and the final 3 prime strand may be one of the two new strands
            created in this step. Check it.
            3) Create the Xover
            """
            c = None
            # lookup the initial strandset index
            found, overlap, ssIdx3p = ss3p._findIndexOfRangeFor(strand3p)
            if strand3p.idx5Prime() == idx3p:  # yes, idx already matches
                temp5 = xoStrand3 = strand3p
            else:
                offset3p = -1 if ss3p.isDrawn5to3() else 1
                if ss3p.strandCanBeSplit(strand3p, idx3p + offset3p):
                    c = ss3p.SplitCommand(strand3p, idx3p + offset3p, ssIdx3p)
                    # cmds.append(c)
                    xoStrand3 = c._strandHigh if ss3p.isDrawn5to3() else c._strandLow
                    # adjust the target 5prime strand, always necessary if a split happens here
                    if idx5p > idx3p and ss3p.isDrawn5to3():
                        temp5 = xoStrand3
                    elif idx5p < idx3p and not ss3p.isDrawn5to3():
                        temp5 = xoStrand3
                    else:
                        temp5 = c._strandLow if ss3p.isDrawn5to3() else c._strandHigh
                    if useUndoStack:
                        self.undoStack().push(c)
                    else:
                        c.redo()
                else:
                    if useUndoStack:
                        self.undoStack().endMacro()
                        # unclear the applied sequence
                        if self.undoStack().canUndo() and ss5p.isScaffold():
                            self.undoStack().undo()
                    return
                # end if
            if xoStrand3.idx3Prime() == idx5p:
                xoStrand5 = temp5
            else:
                ssIdx5p = ssIdx3p
                # if the strand was split for the strand3p, then we need to adjust the strandset index
                if c:
                    # the insertion index into the set is increases
                    if ss3p.isDrawn5to3():
                        ssIdx5p = ssIdx3p + 1 if idx5p > idx3p else ssIdx3p
                    else:
                        ssIdx5p = ssIdx3p + 1 if idx5p > idx3p else ssIdx3p
                if ss5p.strandCanBeSplit(temp5, idx5p):
                    d = ss5p.SplitCommand(temp5, idx5p, ssIdx5p)
                    # cmds.append(d)
                    xoStrand5 = d._strandLow if ss5p.isDrawn5to3() else d._strandHigh
                    if useUndoStack:
                        self.undoStack().push(d)
                    else:
                        d.redo()
                    # adjust the target 3prime strand, IF necessary
                    if idx5p > idx3p and ss3p.isDrawn5to3():
                        xoStrand3 = xoStrand5
                    elif idx5p < idx3p and not ss3p.isDrawn5to3():
                        xoStrand3 = xoStrand5
                else:
                    if useUndoStack:
                        self.undoStack().endMacro()
                        # unclear the applied sequence
                        if self.undoStack().canUndo() and ss5p.isScaffold():
                            self.undoStack().undo()
                    return
        # end if
        else:  # Do the following if it is in fact a different strand
            # is the 5' end ready for xover installation?
            if strand3p.idx5Prime() == idx3p:  # yes, idx already matches
                xoStrand3 = strand3p
            else:  # no, let's try to split
                offset3p = -1 if ss3p.isDrawn5to3() else 1
                if ss3p.strandCanBeSplit(strand3p, idx3p + offset3p):
                    found, overlap, ssIdx = ss3p._findIndexOfRangeFor(strand3p)
                    if found:
                        c = ss3p.SplitCommand(strand3p, idx3p + offset3p, ssIdx)
                        # cmds.append(c)
                        xoStrand3 = c._strandHigh if ss3p.isDrawn5to3() else c._strandLow
                        if useUndoStack:
                            self.undoStack().push(c)
                        else:
                            c.redo()
                else:  # can't split... abort
                    if useUndoStack:
                        self.undoStack().endMacro()
                        # unclear the applied sequence
                        if self.undoStack().canUndo() and ss5p.isScaffold():
                            self.undoStack().undo()
                    return

            # is the 3' end ready for xover installation?
            if strand5p.idx3Prime() == idx5p:  # yes, idx already matches
                xoStrand5 = strand5p
            else:
                if ss5p.strandCanBeSplit(strand5p, idx5p):
                    found, overlap, ssIdx = ss5p._findIndexOfRangeFor(strand5p)
                    if found:
                        d = ss5p.SplitCommand(strand5p, idx5p, ssIdx)
                        # cmds.append(d)
                        xoStrand5 = d._strandLow if ss5p.isDrawn5to3() else d._strandHigh
                        if useUndoStack:
                            self.undoStack().push(d)
                        else:
                            d.redo()
                else:  # can't split... abort
                    if useUndoStack:
                        self.undoStack().endMacro()
                        # unclear the applied sequence
                        if self.undoStack().canUndo() and ss5p.isScaffold():
                            self.undoStack().undo()
                    return
        # end else

        e = Part.CreateXoverCommand(self, xoStrand5, idx5p, xoStrand3, idx3p, updateOligo=updateOligo)
        if useUndoStack:
            self.undoStack().push(e)
            self.undoStack().endMacro()
        else:
            e.redo()

    # end def

    def removeXover(self, strand5p, strand3p, useUndoStack=True):
        cmds = []
        if strand5p.connection3p() == strand3p:
            c = Part.RemoveXoverCommand(self, strand5p, strand3p)
            cmds.append(c)
            util.execCommandList(self, cmds, desc="Remove Xover", \
                                                    useUndoStack=useUndoStack)
    # end def

    def destroy(self):
        self.setParent(None)
        self.deleteLater()  # QObject also emits a destroyed() Signal
    # end def

    def generatorFullLattice(self):
        """
        Returns a generator that yields the row, column lattice points to draw
        relative to the part origin.
        """
        return product(list(range(self._maxRow)), list(range(self._maxCol)))
    # end def

    def generatorSpatialLattice(self, scaleFactor=1.0):
        """
        Returns a generator that yields the XY spatial lattice points to draw
        relative to the part origin.
        """
        # nested for loop in one line
        latticeCoordToPositionXY = self.latticeCoordToPositionXY
        for latticeCoord in product(list(range(self._maxRow)), list(range(self._maxCol))):
            row, col = latticeCoord
            x, y = latticeCoordToPositionXY(row, col, scaleFactor)
            yield x, y, row, col
    # end def

    def getPreXoversHigh(self, strandType, neighborType, minIdx=0, maxIdx=None):
        """
        Returns all prexover positions for neighborType that are below
        maxIdx. Used in emptyhelixitem.py.
        """
        preXO = self._scafH if strandType == StrandType.Scaffold else self._stapH
        if maxIdx == None:
            maxIdx = self._maxBase
        steps = (self._maxBase // self._step) + 1
        ret = [i * self._step + j for i in range(steps) for j in preXO[neighborType]]
        return [x for x in ret if x >= minIdx and x <= maxIdx]

    def getPreXoversLow(self, strandType, neighborType, minIdx=0, maxIdx=None):
        """
        Returns all prexover positions for neighborType that are above
        minIdx. Used in emptyhelixitem.py.
        """
        preXO = self._scafL if strandType == StrandType.Scaffold \
                                else self._stapL
        if maxIdx == None:
            maxIdx = self._maxBase
        steps = (self._maxBase // self._step) + 1
        ret = [i * self._step + j for i in range(steps) for j in preXO[neighborType]]
        return [x for x in ret if x >= minIdx and x <= maxIdx]

    def latticeCoordToPositionXY(self, row, col, scaleFactor=1.0):
        """
        Returns a tuple of the (x,y) position for a given lattice row and
        column.

        Note: The x,y position is the upperLeftCorner for the given
        coordinate, and relative to the part instance.
        """
        raise NotImplementedError  # To be implemented by Part subclass
    # end def

    def positionToCoord(self, x, y, scaleFactor=1.0):
        """
        Returns a tuple (row, column) lattice coordinate for a given
        x and y position that is within +/- 0.5 of a true valid lattice
        position.

        Note: mapping should account for int-to-float rounding errors.
        x,y is relative to the Part Instance Position.
        """
        raise NotImplementedError  # To be implemented by Part subclass
    # end def

    def newPart(self):
        return Part(self._document)
    # end def

    def removeOligo(self, oligo):
        # Not a designated method
        # (there exist methods that also directly
        # remove parts from self._oligos)
        try:
            self._oligos.remove(oligo)
        except KeyError:
            print(util.trace(5))
            # print "error removing oligo", oligo
    # end def

    def renumber(self, coordList, useUndoStack=True):
        if useUndoStack:
            self.undoStack().beginMacro("Renumber VirtualHelices")
        c = Part.RenumberVirtualHelicesCommand(self, coordList)
        if useUndoStack:
            self.undoStack().push(c)
            self.undoStack().endMacro()
        else:
            c.redo()
    # end def
    
    class RenumberVirtualHelicesCommand(QUndoCommand):
        """
        """
        def __init__(self, part, coordList):
            super(Part.RenumberVirtualHelicesCommand, self).__init__()
            self._part = part
            self._vhs = [part.virtualHelixAtCoord(coord) for coord in coordList]
            self._oldNumbers = [vh.number() for vh in self._vhs]
        # end def
            
        def redo(self):
            even = 0
            odd = 1
            for vh in self._vhs:
                if vh.isEvenParity():
                    vh.setNumber(even)
                    even += 2
                else:
                    vh.setNumber(odd)
                    odd += 2
            # end for
            part = self._part
            aVH =  part.activeVirtualHelix()
            if aVH:
                part.partStrandChangedSignal.emit(part, aVH)
            for oligo in part._oligos:
                for strand in oligo.strand5p().generator3pStrand():
                    strand.strandUpdateSignal.emit(strand)
        # end def
            
        def undo(self):
            for vh, num in zip(self._vhs, self._oldNumbers):
                vh.setNumber(num)
            # end for
            part = self._part
            aVH =  part.activeVirtualHelix()
            if aVH:
                part.partStrandChangedSignal.emit(part, aVH)
            for oligo in part._oligos:
                for strand in oligo.strand5p().generator3pStrand():
                    strand.strandUpdateSignal.emit(strand)
        # end def
    # end def

    def resizeLattice(self):
        """docstring for resizeLattice"""
        pass
    # end def

    def resizeVirtualHelices(self, minDelta, maxDelta, useUndoStack=True):
        """docstring for resizeVirtualHelices"""
        c = Part.ResizePartCommand(self, minDelta, maxDelta)
        util.execCommandList(self, [c], desc="Resize part", \
                                                    useUndoStack=useUndoStack)
    # end def

    def setActiveBaseIndex(self, idx):
        self._activeBaseIndex = idx
        self.partActiveSliceIndexSignal.emit(self, idx)
    # end def

    def setActiveVirtualHelix(self, virtualHelix, idx=None):
        self._activeVirtualHelix = virtualHelix
        self._activeVirtualHelixIdx = idx
        self.partStrandChangedSignal.emit(self, virtualHelix)
    # end def

    def selectPreDecorator(self, selectionList):
        """
        Handles view notifications that a predecorator has been selected.
        """
        if (len(selectionList) == 0):
            return
            # print "all PreDecorators were unselected"
            # partPreDecoratorUnSelectedSignal.emit()
        sel = selectionList[0]
        (row, col, baseIdx) = (sel[0], sel[1], sel[2])
        self.partPreDecoratorSelectedSignal.emit(self, row, col, baseIdx)

    def xoverSnapTo(self, strand, idx, delta):
        """
        Returns the nearest xover position to allow snap-to behavior in
        resizing strands via dragging selected xovers.
        """
        strandType = strand.strandType()
        if delta > 0:
            minIdx, maxIdx = idx - delta, idx + delta
        else:
            minIdx, maxIdx = idx + delta, idx - delta

        # determine neighbor strand and bind the appropriate prexover method
        lo, hi = strand.idxs()
        if idx == lo:
            connectedStrand = strand.connectionLow()
            preXovers = self.getPreXoversHigh
        else:
            connectedStrand = strand.connectionHigh()
            preXovers = self.getPreXoversLow
        connectedVh = connectedStrand.virtualHelix()

        # determine neighbor position, if any
        neighbors = self.getVirtualHelixNeighbors(strand.virtualHelix())
        if connectedVh in neighbors:
            neighborIdx = neighbors.index(connectedVh)
            try:
                newIdx = util.nearest(idx + delta,
                                    preXovers(strandType,
                                                neighborIdx,
                                                minIdx=minIdx,
                                                maxIdx=maxIdx)
                                    )
                return newIdx
            except ValueError:
                return None  # nearest not found in the expanded list
        else:  # no neighbor (forced xover?)... don't snap, just return
            return idx + delta

    ### PRIVATE SUPPORT METHODS ###
    def _addVirtualHelix(self, virtualHelix):
        """
        private method for adding a virtualHelix to the Parts data structure
        of virtualHelix references
        """
        self._coordToVirtualHelix[virtualHelix.coord()] = virtualHelix
    # end def

    def _removeVirtualHelix(self, virtualHelix):
        """
        private method for adding a virtualHelix to the Parts data structure
        of virtualHelix references
        """
        del self._coordToVirtualHelix[virtualHelix.coord()]
    # end def

    def _reserveHelixIDNumber(self, parityEven=True, requestedIDnum=None):
        """
        Reserves and returns a unique numerical label appropriate for a
        virtualhelix of a given parity. If a specific index is preferable
        (say, for undo/redo) it can be requested in num.
        """
        num = requestedIDnum
        if num != None:  # We are handling a request for a particular number
            assert num >= 0, int(num) == num
            # assert not num in self._numberToVirtualHelix
            if num in self.oddRecycleBin:
                self.oddRecycleBin.remove(num)
                heapify(self.oddRecycleBin)
                return num
            if num in self.evenRecycleBin:
                self.evenRecycleBin.remove(num)
                heapify(self.evenRecycleBin)
                return num
            self.reserveBin.add(num)
            return num
        # end if
        else:
            # Hybrid parts share one visible helix-number namespace. Keep the
            # native even/odd parity rule, but skip numbers already used by
            # the other lattice part.
            document = self.document()
            otherNumbers = set()
            if document is not None and document.isHybrid():
                otherNumbers = set(
                    vh.number()
                    for otherPart in document.parts()
                    if otherPart is not self
                    for vh in otherPart.getVirtualHelices())

            def recycledNumber(heap):
                held = []
                available = None
                while heap:
                    candidate = heappop(heap)
                    if candidate not in otherNumbers:
                        available = candidate
                        break
                    held.append(candidate)
                for candidate in held:
                    heappush(heap, candidate)
                return available

            # Just find any valid index (subject to parity constraints)
            if parityEven:
                if len(self.evenRecycleBin):
                    available = recycledNumber(self.evenRecycleBin)
                    if available is not None:
                        return available
                while self._highestUsedEven + 2 in self.reserveBin or \
                        self._highestUsedEven + 2 in otherNumbers:
                    self._highestUsedEven += 2
                self._highestUsedEven += 2
                return self._highestUsedEven
            else:
                if len(self.oddRecycleBin):
                    available = recycledNumber(self.oddRecycleBin)
                    if available is not None:
                        return available
                # use self._highestUsedOdd iff the recycle bin is empty
                # and highestUsedOdd+2 is not reserved or used cross-lattice
                while self._highestUsedOdd + 2 in self.reserveBin or \
                        self._highestUsedOdd + 2 in otherNumbers:
                    self._highestUsedOdd += 2
                self._highestUsedOdd += 2
                return self._highestUsedOdd
        # end else
    # end def

    def _recycleHelixIDNumber(self, n):
        """
        The caller's contract is to ensure that n is not used in *any* helix
        at the time of the calling of this function (or afterwards, unless
        reserveLabelForHelix returns the label again).
        """
        if n % 2 == 0:
            heappush(self.evenRecycleBin, n)
        else:
            heappush(self.oddRecycleBin, n)
    # end def

    def _splitBeforeAutoXovers(self, vh5p, vh3p, idx, useUndoStack=True):
        # prexoveritem needs to store left or right, and determine
        # locally whether it is from or to
        # pass that info in here in and then do the breaks
        ss5p = strand5p.strandSet()
        ss3p = strand3p.strandSet()
        cmds = []

        # is the 5' end ready for xover installation?
        if strand3p.idx5Prime() == idx5p:  # yes, idx already matches
            xoStrand3 = strand3p
        else:  # no, let's try to split
            offset3p = -1 if ss3p.isDrawn5to3() else 1
            if ss3p.strandCanBeSplit(strand3p, idx3p + offset3p):
                found, overlap, ssIdx = ss3p._findIndexOfRangeFor(strand3p)
                if found:
                    c = ss3p.SplitCommand(strand3p, idx3p + offset3p, ssIdx)
                    cmds.append(c)
                    xoStrand3 = c._strandHigh if ss3p.isDrawn5to3() else c._strandLow
            else:  # can't split... abort
                return

        # is the 3' end ready for xover installation?
        if strand5p.idx3Prime() == idx5p:  # yes, idx already matches
            xoStrand5 = strand5p
        else:
            if ss5p.strandCanBeSplit(strand5p, idx5p):
                found, overlap, ssIdx = ss5p._findIndexOfRangeFor(strand5p)
                if found:
                    d = ss5p.SplitCommand(strand5p, idx5p, ssIdx)
                    cmds.append(d)
                    xoStrand5 = d._strandLow if ss5p.isDrawn5to3() \
                                                else d._strandHigh
            else:  # can't split... abort
                return
        c = Part.CreateXoverCommand(self, xoStrand5, idx5p, xoStrand3, idx3p)
        cmds.append(c)
        util.execCommandList(self, cmds, desc="Create Xover", \
                                                useUndoStack=useUndoStack)
    # end def

    ### PUBLIC SUPPORT METHODS ###
    def shallowCopy(self):
        part = self.newPart()
        part._virtualHelices = dict(self._virtualHelices)
        part._oligos = set(self._oligos)
        part._maxBase = self._maxBase
        return part
    # end def

    def deepCopy(self):
        """
        1) Create a new part
        2) copy the VirtualHelices
        3) Now you need to map the ORIGINALs Oligos onto the COPY's Oligos
        To do this you can for each Oligo in the ORIGINAL
            a) get the strand5p() of the ORIGINAL
            b) get the corresponding strand5p() in the COPY based on
                i) lookup the hash idNum of the ORIGINAL strand5p() VirtualHelix
                ii) get the StrandSet() that you created in Step 2 for the
                StrandType of the original using the hash idNum
        """
        # 1) new part
        part = self.newPart()
        for key, vhelix in self._virtualHelices:
            # 2) Copy VirtualHelix
            part._virtualHelices[key] = vhelix.deepCopy(part)
        # end for
        # 3) Copy oligos
        for oligo, val in self._oligos:
            strandGenerator = oligo.strand5p().generator3pStrand()
            strandType = oligo.strand5p().strandType()
            newOligo = oligo.deepCopy(part)
            lastStrand = None
            for strand in strandGenerator:
                idNum = strand.virtualHelix().number()
                newVHelix = part._virtualHelices[idNum]
                newStrandSet = newVHelix().getStrandSetByType(strandType)
                newStrand = strand.deepCopy(newStrandSet, newOligo)
                if lastStrand:
                    lastStrand.setConnection3p(newStrand)
                else:
                    # set the first condition
                    newOligo.setStrand5p(newStrand)
                newStrand.setConnection5p(lastStrand)
                newStrandSet.addStrand(newStrand)
                lastStrand = newStrand
            # end for
            # check loop condition
            if oligo.isLoop():
                s5p = newOligo.strand5p()
                lastStrand.set3pconnection(s5p)
                s5p.set5pconnection(lastStrand)
            # add to part
            oligo.add()
        # end for
        return part
    # end def

    def areSameOrNeighbors(self, virtualHelixA, virtualHelixB):
        """
        returns True or False
        """
        return virtualHelixB in self.getVirtualHelixNeighbors(virtualHelixA) or \
            virtualHelixA == virtualHelixB
    # end def

    def potentialCrossoverList(self, virtualHelix, idx=None):
        """
        Returns a list of tuples
            (neighborVirtualHelix, index, strandType, isLowIdx)

        where:

        neighborVirtualHelix is a virtualHelix neighbor of the arg virtualHelix
        index is the index where a potential Xover might occur
        strandType is from the enum (StrandType.Scaffold, StrandType.Staple)
        isLowIdx is whether or not it's the at the low index (left in the Path
        view) of a potential Xover site
        """
        vh = virtualHelix
        ret = []  # LUT = Look Up Table
        part = self
        # these are the list of crossover points simplified
        # they depend on whether the strandType is scaffold or staple
        # create a list of crossover points for each neighbor of the form
        # [(_scafL[i], _scafH[i], _stapL[i], _stapH[i]), ...]
        lutsNeighbor = list(
                            zip(
                                part._scafL,
                                part._scafH,
                                part._stapL,
                                part._stapH
                                )
                            )

        sTs = (StrandType.Scaffold, StrandType.Staple)
        numBases = part.maxBaseIdx()

        # create a range for the helical length dimension of the Part,
        # incrementing by the lattice step size.
        baseRange = list(range(0, numBases, part._step))

        if idx != None:
            baseRange = [x for x in baseRange if x >= idx - 3 * part._step and \
                                        x <= idx + 2 * part._step]

        if vh is None:
            return

        fromStrandSets = vh.getStrandSets()
        neighbors = self.getVirtualHelixNeighbors(vh)

        # print neighbors, lutsNeighbor
        for neighbor, lut in zip(neighbors, lutsNeighbor):
            if not neighbor:
                continue

            # now arrange again for iteration
            # (_scafL[i], _scafH[i]), (_stapL[i], _stapH[i]) )
            # so we can pair by StrandType
            lutScaf = lut[0:2]
            lutStap = lut[2:4]
            lut = (lutScaf, lutStap)

            toStrandSets = neighbor.getStrandSets()
            for fromSS, toSS, pts, st in zip(fromStrandSets, toStrandSets, lut, sTs):
                # test each period of each lattice for each StrandType
                for pt, isLowIdx in zip(pts, (True, False)):
                    for i, j in product(baseRange, pt):
                        index = i + j
                        # ``numBases`` is historically named but is actually
                        # maxBaseIdx().  Include that final physical base.
                        # Excluding it drops the reciprocal half of a native
                        # crossover block when a circular Curved helix length
                        # is an exact lattice-period multiple (335/0 for a
                        # 336-base Honeycomb ring, and likewise 31/0 in
                        # Square).
                        if index <= numBases:
                            if fromSS.hasNoStrandAtOrNoXover(index) and \
                                    toSS.hasNoStrandAtOrNoXover(index):
                                ret.append((neighbor, index, st, isLowIdx))
                            # end if
                        # end if
                    # end for
                # end for
            # end for
        # end for
        return ret
    # end def

    def supplementalStapleCrossoverList(self, virtualHelix, idx=None):
        """Return denser, UI-only staple crossover candidates.

        These candidates supplement the native staple lookup tables so the
        total manual staple density matches scaffold density. They are kept
        out of ``potentialCrossoverList`` intentionally, which means the
        existing AutoCS_staples algorithm continues using native sites
        only.
        """
        vh = virtualHelix
        if vh is None:
            return []

        if self._step == 21:       # Honeycomb: 1 native + 1 supplemental
            phaseOffsets = (10,)
        elif self._step == 32:     # Square: 1 native + 2 supplemental
            phaseOffsets = (11, 21)
        else:
            return []

        ret = []
        numBases = self.maxBaseIdx()
        baseRange = list(range(0, numBases, self._step))
        if idx is not None:
            baseRange = [base for base in baseRange
                         if base >= idx - 3 * self._step and
                         base <= idx + 2 * self._step]

        fromSS = vh.stapleStrandSet()
        neighbors = self.getVirtualHelixNeighbors(vh)
        for neighbor, nativeLow, nativeHigh in zip(
                            neighbors, self._stapL, self._stapH):
            if neighbor is None:
                continue
            toSS = neighbor.stapleStrandSet()
            for nativePoints, isLowIdx in ((nativeLow, True),
                                           (nativeHigh, False)):
                nativeSet = set(nativePoints)
                supplementalPoints = set()
                for nativePoint in nativePoints:
                    for phaseOffset in phaseOffsets:
                        point = (nativePoint + phaseOffset) % self._step
                        if point not in nativeSet:
                            supplementalPoints.add(point)
                for base, point in product(baseRange,
                                           sorted(supplementalPoints)):
                    index = base + point
                    if index < numBases and \
                            fromSS.hasNoStrandAtOrNoXover(index) and \
                            toSS.hasNoStrandAtOrNoXover(index):
                        ret.append((neighbor, index, StrandType.Staple,
                                    isLowIdx))
        return ret
    # end def

    def possibleXoverAt(self, fromVirtualHelix, toVirtualHelix, strandType, idx):
        fromSS = fromVirtualHelix.getStrandSetByType(strandType)
        toSS = toVirtualHelix.getStrandSetByType(strandType)
        return fromSS.hasStrandAtAndNoXover(idx) and \
                toSS.hasStrandAtAndNoXover(idx)
    # end def

    def setImportedVHelixOrder(self, orderedCoordList, emitSignal=True):
        """Used on file import to store the order of the virtual helices."""
        self._importedVHelixOrder = list(orderedCoordList)
        if emitSignal:
            self.partVirtualHelicesReorderedSignal.emit(
                self, self._importedVHelixOrder)

    ### COMMANDS ###
    class CreateVirtualHelixCommand(QUndoCommand):
        def __init__(self, part, row, col):
            super(Part.CreateVirtualHelixCommand, self).__init__()
            self._part = part
            self._parityEven = part.isEvenParity(row, col)
            idNum = part._reserveHelixIDNumber(self._parityEven,
                                                requestedIDnum=None)
            self._vhelix = VirtualHelix(part, row, col, idNum)
            self._idNum = idNum
        # end def

        def redo(self):
            vh = self._vhelix
            part = self._part
            idNum = self._idNum
            vh.setPart(part)
            part._addVirtualHelix(vh)
            vh.setNumber(idNum)
            if not vh.number():
                part._reserveHelixIDNumber(self._parityEven,
                                            requestedIDnum=idNum)
            # end if
            part.partVirtualHelixAddedSignal.emit(part, vh)
            part.partActiveSliceResizeSignal.emit(part)
        # end def

        def undo(self):
            vh = self._vhelix
            part = self._part
            idNum = self._idNum
            part._removeVirtualHelix(vh)
            part._recycleHelixIDNumber(idNum)
            # clear out part references
            vh.setNumber(None)  # must come before setPart(None)
            vh.setPart(None)
            vh.virtualHelixRemovedSignal.emit(vh)
            part.partActiveSliceResizeSignal.emit(part)
        # end def
    # end class

    class CreateXoverCommand(QUndoCommand):
        """
        Creates a Xover from the 3' end of strand5p to the 5' end of strand3p
        this needs to
        1. preserve the old oligo of strand3p
        2. install the crossover
        3. apply the strand5p oligo to the strand3p
        """
        def __init__(self, part, strand5p, strand5pIdx, strand3p, strand3pIdx, updateOligo=True):
            super(Part.CreateXoverCommand, self).__init__()
            self._part = part
            self._strand5p = strand5p
            self._strand5pIdx = strand5pIdx
            self._strand3p = strand3p
            self._strand3pIdx = strand3pIdx
            self._oldOligo3p = strand3p.oligo()
            self._oldOligo3pPart = strand3p.oligo().part()
            self._updateOligo = updateOligo
        # end def

        def redo(self):
            part = self._part
            strand5p = self._strand5p
            strand5pIdx = self._strand5pIdx
            strand3p = self._strand3p
            strand3pIdx = self._strand3pIdx
            olg5p = strand5p.oligo()
            oldOlg3p = self._oldOligo3p

            # 0. Deselect the involved strands
            doc = strand5p.document()
            doc.removeStrandFromSelection(strand5p)
            doc.removeStrandFromSelection(strand3p)

            if self._updateOligo:
                # Test for Loopiness
                if olg5p == strand3p.oligo():
                    olg5p.setLoop(True)
                else:
                    # 1. update preserved oligo length
                    olg5p.incrementLength(oldOlg3p.length())
                    # 2. Remove the old oligo and apply the 5' oligo to the 3' strand
                    oldOlg3p.removeFromPart()
                    for strand in strand3p.generator3pStrand():
                        # emits strandHasNewOligoSignal
                        Strand.setOligo(strand, olg5p)

            # 3. install the Xover
            strand5p.setConnection3p(strand3p)
            strand3p.setConnection5p(strand5p)

            ss5 = strand5p.strandSet()
            vh5p = ss5.virtualHelix()
            st5p = ss5.strandType()
            ss3 = strand3p.strandSet()
            vh3p = ss3.virtualHelix()
            st3p = ss3.strandType()

            part.partActiveVirtualHelixChangedSignal.emit(part, vh5p)
            # strand5p.strandXover5pChangedSignal.emit(strand5p, strand3p)
            if self._updateOligo:
                strand5p.strandUpdateSignal.emit(strand5p)
                strand3p.strandUpdateSignal.emit(strand3p)
        # end def

        def undo(self):
            part = self._part
            strand5p = self._strand5p
            strand5pIdx = self._strand5pIdx
            strand3p = self._strand3p
            strand3pIdx = self._strand3pIdx
            oldOlg3p = self._oldOligo3p
            olg5p = strand5p.oligo()

            # 0. Deselect the involved strands
            doc = strand5p.document()
            doc.removeStrandFromSelection(strand5p)
            doc.removeStrandFromSelection(strand3p)

            # 1. uninstall the Xover
            strand5p.setConnection3p(None)
            strand3p.setConnection5p(None)

            if self._updateOligo:
                # Test Loopiness
                if oldOlg3p.isLoop():
                    oldOlg3p.setLoop(False)
                else:
                    # 2. restore the modified oligo length
                    olg5p.decrementLength(oldOlg3p.length())
                    # 3. apply the old oligo to strand3p
                    oldOlg3p.addToPart(self._oldOligo3pPart)
                    for strand in strand3p.generator3pStrand():
                        # emits strandHasNewOligoSignal
                        Strand.setOligo(strand, oldOlg3p)

            ss5 = strand5p.strandSet()
            vh5p = ss5.virtualHelix()
            st5p = ss5.strandType()
            ss3 = strand3p.strandSet()
            vh3p = ss3.virtualHelix()
            st3p = ss3.strandType()

            part.partActiveVirtualHelixChangedSignal.emit(part, vh5p)
            # strand5p.strandXover5pChangedSignal.emit(strand5p, strand3p)
            if self._updateOligo:
                strand5p.strandUpdateSignal.emit(strand5p)
                strand3p.strandUpdateSignal.emit(strand3p)
        # end def
    # end class

    class RefreshOligosCommand(QUndoCommand):
        """
        RefreshOligosCommand is a post-processing step for AutoStaple.

        Normally when an xover is created, all strands in the 3' direction are
        assigned the oligo of the 5' strand. This becomes very expensive
        during autoStaple, because the Nth xover requires updating up to N-1
        strands.

        Hence, we disable oligo assignment during the xover creation step,
        and then do it all in one pass at the end with this command.
        """
        def __init__(self, part):
            super(Part.RefreshOligosCommand, self).__init__()
            self._part = part
        # end def

        def redo(self):
            visited = {}
            for vh in self._part.getVirtualHelices():
                stapSS = vh.stapleStrandSet()
                for strand in stapSS:
                    visited[strand] = False

            for strand in list(visited.keys()):
                if visited[strand]:
                    continue
                visited[strand] = True
                startOligo = strand.oligo()
                strand5gen = strand.generator5pStrand()
                # this gets the oligo and burns a strand in the generator
                strand5 = next(strand5gen)
                for strand5 in strand5gen:
                    oligo5 = strand5.oligo()
                    if oligo5 != startOligo:
                        oligo5.removeFromPart()
                        Strand.setOligo(strand5, startOligo)  # emits strandHasNewOligoSignal
                    visited[strand5] = True
                # end for
                startOligo.setStrand5p(strand5)
                # is it a loop?
                if strand.connection3p() == strand5:
                    startOligo.setLoop(True)
                else:
                    strand3gen = strand.generator3pStrand()
                    strand3 = next(strand3gen)   # burn one
                    for strand3 in strand3gen:
                        oligo3 = strand3.oligo()
                        if oligo3 != startOligo:
                            oligo3.removeFromPart()
                            Strand.setOligo(strand3, startOligo)  # emits strandHasNewOligoSignal
                        visited[strand3] = True
                    # end for
                startOligo.refreshLength()
            # end for

            oligoSet = set()
            for strand in list(visited.keys()):
                oligoSet.add(strand.oligo())
                strand.strandUpdateSignal.emit(strand)
        # end def

        def undo(self):
            """Doesn't reassign """
            pass
        # end def
    # end class

    class RemoveXoverCommand(QUndoCommand):
        """
        Removes a Xover from the 3' end of strand5p to the 5' end of strand3p
        this needs to
        1. preserve the old oligo of strand3p
        2. install the crossover
        3. update the oligo length
        4. apply the new strand3p oligo to the strand3p
        """
        def __init__(self, part, strand5p, strand3p):
            super(Part.RemoveXoverCommand, self).__init__()
            self._part = part
            self._strand5p = strand5p
            self._strand5pIdx = strand5p.idx3Prime()
            self._strand3p = strand3p
            self._strand3pIdx = strand3p.idx5Prime()
            nO3p = self._newOligo3p = strand3p.oligo().shallowCopy()
            self._newOligo3pPart = strand3p.part()
            # Both staple and scaffold products retain the color of the
            # original oligo when a crossover is removed.
            nO3p.setColor(strand3p.oligo().color())
            nO3p.setLength(0)
            for strand in strand3p.generator3pStrand():
                nO3p.incrementLength(strand.totalLength())
            # end def
            nO3p.setStrand5p(strand3p)
            
            self._isLoop = strand3p.oligo().isLoop()
        # end def

        def redo(self):
            part = self._part
            strand5p = self._strand5p
            strand5pIdx = self._strand5pIdx
            strand3p = self._strand3p
            strand3pIdx = self._strand3pIdx
            newOlg3p = self._newOligo3p
            olg5p = self._strand5p.oligo()

            # 0. Deselect the involved strands
            doc = strand5p.document()
            doc.removeStrandFromSelection(strand5p)
            doc.removeStrandFromSelection(strand3p)

            # 1. uninstall the Xover
            strand5p.setConnection3p(None)
            strand3p.setConnection5p(None)

            if self._isLoop:
                olg5p.setLoop(False)
                olg5p.setStrand5p(strand3p)
            else:
                # 2. restore the modified oligo length
                olg5p.decrementLength(newOlg3p.length())
                # 3. apply the old oligo to strand3p
                newOlg3p.addToPart(self._newOligo3pPart)
                for strand in strand3p.generator3pStrand():
                    # emits strandHasNewOligoSignal
                    Strand.setOligo(strand, newOlg3p)

            ss5 = strand5p.strandSet()
            vh5p = ss5.virtualHelix()
            st5p = ss5.strandType()
            ss3 = strand3p.strandSet()
            vh3p = ss3.virtualHelix()
            st3p = ss3.strandType()

            part.partActiveVirtualHelixChangedSignal.emit(part, vh5p)
            # strand5p.strandXover5pChangedSignal.emit(strand5p, strand3p)
            strand5p.strandUpdateSignal.emit(strand5p)
            strand3p.strandUpdateSignal.emit(strand3p)
        # end def

        def undo(self):
            part = self._part
            strand5p = self._strand5p
            strand5pIdx = self._strand5pIdx
            strand3p = self._strand3p
            strand3pIdx = self._strand3pIdx
            olg5p = strand5p.oligo()
            newOlg3p = self._newOligo3p

            # 0. Deselect the involved strands
            doc = strand5p.document()
            doc.removeStrandFromSelection(strand5p)
            doc.removeStrandFromSelection(strand3p)

            if self._isLoop:
                olg5p.setLoop(True)
                # No need to restore whatever the old Oligo._strand5p was
            else:
                # 1. update preserved oligo length
                olg5p.incrementLength(newOlg3p.length())
                # 2. Remove the old oligo and apply the 5' oligo to the 3' strand
                newOlg3p.removeFromPart()
                for strand in strand3p.generator3pStrand():
                    # emits strandHasNewOligoSignal
                    Strand.setOligo(strand, olg5p)
            # end else

            # 3. install the Xover
            strand5p.setConnection3p(strand3p)
            strand3p.setConnection5p(strand5p)

            ss5 = strand5p.strandSet()
            vh5p = ss5.virtualHelix()
            st5p = ss5.strandType()
            ss3 = strand3p.strandSet()
            vh3p = ss3.virtualHelix()
            st3p = ss3.strandType()

            part.partActiveVirtualHelixChangedSignal.emit(part, vh5p)
            # strand5p.strandXover5pChangedSignal.emit(strand5p, strand3p)
            strand5p.strandUpdateSignal.emit(strand5p)
            strand3p.strandUpdateSignal.emit(strand3p)
        # end def
    # end class

    class RemovePartCommand(QUndoCommand):
        """
        RemovePartCommand deletes a part. Emits partRemovedSignal.
        """
        def __init__(self, part):
            super(Part.RemovePartCommand, self).__init__()
            self._part = part
            self._doc = part.document()
        # end def

        def redo(self):
            # Remove the strand
            part = self._part
            doc = self._doc
            doc.removePart(part)
            part.setDocument(None)
            part.partRemovedSignal.emit(part)
        # end def

        def undo(self):
            part = self._part
            doc = self._doc
            doc._addPart(part)
            part.setDocument(doc)
            doc.documentPartAddedSignal.emit(doc, part)
        # end def
    # end class

    class RemoveAllStrandsCommand(QUndoCommand):
        """
        1. Remove all strands. Emits strandRemovedSignal for each.
        2. Remove all oligos. 
        """
        def __init__(self, part):
            super(Part.RemoveAllStrandsCommand, self).__init__()
            self._part = part
            self._vhs = vhs = part.getVirtualHelices()
            self._strandSets = []
            for vh in self._vhs:
                x = vh.getStrandSets()
                self._strandSets.append(x[0])
                self._strandSets.append(x[1])
            self._strandSetListCopies = \
                        [[y for y in x._strandList] for x in self._strandSets]
            self._oligos = set(part.oligos())
        # end def

        def redo(self):
            part = self._part
            # Remove the strand
            for sSet in self._strandSets:
                sList = sSet._strandList
                for strand in sList:
                    sSet.removeStrand(strand)
                # end for
                sSet._strandList = []
            #end for
            for vh in self._vhs:
                # for updating the Slice View displayed helices
                part.partStrandChangedSignal.emit(part, vh)
            # end for
            self._oligos.clear()
        # end def

        def undo(self):
            part = self._part
            # Remove the strand
            sListCopyIterator = iter(self._strandSetListCopies)
            for sSet in self._strandSets:
                sList = next(sListCopyIterator)
                for strand in sList:
                    sSet.strandsetStrandAddedSignal.emit(sSet, strand)
                # end for
                sSet._strandList = sList
            #end for
            for vh in self._vhs:
                # for updating the Slice View displayed helices
                part.partStrandChangedSignal.emit(part, vh)
            # end for
            for olg in self._oligos:
                part.addOligo(olg)
        # end def
    # end class

    class ResizePartCommand(QUndoCommand):
        """
        set the maximum and mininum base index in the helical direction

        need to adjust all subelements in the event of a change in the
        minimum index
        """
        def __init__(self, part, minHelixDelta, maxHelixDelta):
            super(Part.ResizePartCommand, self).__init__()
            self._part = part
            self._minDelta = minHelixDelta
            self._maxDelta = maxHelixDelta
            self._oldActiveIdx = part.activeBaseIndex()
        # end def

        def redo(self):
            part = self._part
            part._minBase += self._minDelta
            part._maxBase += self._maxDelta
            if self._minDelta != 0:
                self.deltaMinDimension(part, self._minDelta)
            for vh in part._coordToVirtualHelix.values():
                part.partVirtualHelixResizedSignal.emit(part, vh.coord())
            if self._oldActiveIdx > part._maxBase:
                part.setActiveBaseIndex(part._maxBase)
            part.partDimensionsChangedSignal.emit(part)
        # end def

        def undo(self):
            part = self._part
            part._minBase -= self._minDelta
            part._maxBase -= self._maxDelta
            if self._minDelta != 0:
                self.deltaMinDimension(part, self._minDelta)
            for vh in part._coordToVirtualHelix.values():
                part.partVirtualHelixResizedSignal.emit(part, vh.coord())
            if self._oldActiveIdx != part.activeBaseIndex():
                part.setActiveBaseIndex(self._oldActiveIdx)
            part.partDimensionsChangedSignal.emit(part)
        # end def

        def deltaMinDimension(self, part, minDimensionDelta):
            """
            Need to update:
            strands
            insertions
            """
            for vhDict in part._insertions.values():
                for insertion in vhDict:
                    insertion.updateIdx(minDimensionDelta)
                # end for
            # end for
            for vh in part._coordToVirtualHelix.values():
                for strand in vh.scaffoldStrand().generatorStrand():
                    strand.updateIdxs(minDimensionDelta)
                for strand in vh.stapleStrand().generatorStrand():
                    strand.updateIdxs(minDimensionDelta)
            # end for
        # end def
    # end class
# end class
