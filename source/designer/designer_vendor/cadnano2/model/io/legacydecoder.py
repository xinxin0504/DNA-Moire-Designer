from collections import defaultdict

from cadnano2.model.document import Document
from cadnano2.model.enum import LatticeType, StrandType
from cadnano2.model.parts.honeycombpart import HoneycombPart
from cadnano2.model.parts.squarepart import SquarePart
from cadnano2.model.virtualhelix import VirtualHelix
from cadnano2.views import styles
import cadnano2.util as util
import cadnano2.cadnano as cadnano
# import Qt stuff into the module namespace with PySide, PyQt4 independence
util.qtWrapImport('QtGui', globals(),  ['QColor'])
if cadnano.app().isGui():
    from cadnano2.ui.dialogs.ui_latticetype import Ui_LatticeType
    util.qtWrapImport('QtWidgets', globals(),  ['QDialog', 'QDialogButtonBox'])

NODETAG = "node"
NAME = "name"
OBJ_ID = "objectid"
INST_ID = "instanceid"
DONE = "done"
CHECKED = "check"
LOCKED = "locked"

VHELIX = "vhelix"
NUM = "num"
COL = "col"
ROW = "row"
SCAFFOLD = "scaffold"
STAPLE = "staple"
INSERTION = "insertion"
DELETION = "deletion"


def import_legacy_dict(document, obj, latticeType=LatticeType.Honeycomb,
                       forceLatticeType=False):
    """
    Parses a dictionary (obj) created from reading a json file and uses it
    to populate the given document with model data.
    """
    vstrands = obj.get('vstrands', [])
    defaultBases = (21 if latticeType == LatticeType.Honeycomb else 32)
    numBases = (len(vstrands[0]['scaf']) if vstrands else
                int(obj.get('num_bases', defaultBases)))
    if cadnano.app().isGui() and not forceLatticeType:
        # from ui.dialogs.ui_latticetype import Ui_LatticeType
        # util.qtWrapImport('QtGui', globals(),  ['QDialog', 'QDialogButtonBox'])
        dialog = QDialog()
        dialogLT = Ui_LatticeType()
        dialogLT.setupUi(dialog)
        # DETERMINE LATTICE TYPE
        if numBases % 21 == 0 and numBases % 32 == 0:
            if dialog.exec() == 1:
                latticeType = LatticeType.Square
            else:
                latticeType = LatticeType.Honeycomb
        elif numBases % 32 == 0:
            latticeType = LatticeType.Square
        elif numBases % 21 == 0:
            latticeType = LatticeType.Honeycomb
        else:
            if dialog.exec() == 1:
                latticeType = LatticeType.Square
            else:
                latticeType = LatticeType.Honeycomb
    else:  # Headless, assume the latticeType arg was meaningful
        pass

    # DETERMINE MAX ROW,COL
    maxRowJson = maxColJson = 0
    for helix in vstrands:
        maxRowJson = max(maxRowJson, int(helix['row'])+1)
        maxColJson = max(maxColJson, int(helix['col'])+1)

    # CREATE PART ACCORDING TO LATTICE TYPE
    if latticeType == LatticeType.Honeycomb:
        steps = numBases/21
        nRows = max(30, maxRowJson, cadnano.app().prefs.honeycombRows)
        nCols = max(32, maxColJson, cadnano.app().prefs.honeycombCols)
        part = HoneycombPart(document=document, maxRow=nRows, maxCol=nCols, maxSteps=steps)
    elif latticeType == LatticeType.Square:
        isSQ100 = bool(vstrands) and not forceLatticeType
        for helix in vstrands:
            if helix['col'] != 0:
                isSQ100 = False
                break
        if isSQ100 and cadnano.app().isGui():
            dialogLT.label.setText("Is this a SQ100 file?")
            if dialog.exec() == 1:
                nRows, nCols = 100, 1
            else:
                nRows, nCols = 40, 30
        else:
            nRows, nCols = 40, 30
        steps = numBases/32
        nRows = max(30, maxRowJson, cadnano.app().prefs.squareRows)
        nCols = max(32, maxColJson, cadnano.app().prefs.squareCols)
        part = SquarePart(document=document, maxRow=nRows, maxCol=nCols, maxSteps=steps)
    else:
        raise TypeError("Lattice type not recognized")
    document._addPart(part, useUndoStack=False)

    # POPULATE VIRTUAL HELICES
    orderedCoordList = []
    vhNumToCoord = {}
    for helix in vstrands:
        vhNum = helix['num']
        row = helix['row']
        col = helix['col']
        scaf= helix['scaf']
        coord = (row, col)
        vhNumToCoord[vhNum] = coord
        orderedCoordList.append(coord)
    # make sure we retain the original order
    for vhNum in sorted(vhNumToCoord.keys()):
        row, col = vhNumToCoord[vhNum]
        part.createVirtualHelix(row, col, useUndoStack=False)
    part.setImportedVHelixOrder(orderedCoordList)

    # INSTALL STRANDS AND COLLECT XOVER LOCATIONS
    numHelixes = len(vstrands)-1
    scaf_seg = defaultdict(list)
    scaf_xo = defaultdict(list)
    stap_seg = defaultdict(list)
    stap_xo = defaultdict(list)
    try:
        for helix in vstrands:
            vhNum = helix['num']
            row = helix['row']
            col = helix['col']
            scaf = helix['scaf']
            stap = helix['stap']
            insertions = helix['loop']
            skips = helix['skip']
            vh = part.virtualHelixAtCoord((row, col))
            scafStrandSet = vh.scaffoldStrandSet()
            stapStrandSet = vh.stapleStrandSet()
            assert(len(scaf)==len(stap) and len(stap)==part.maxBaseIdx()+1 and\
                   len(scaf)==len(insertions) and len(insertions)==len(skips))
            # Build model strands from actual local adjacency.  The historic
            # endpoint heuristic confuses a base that has an external
            # crossover on one side and a nick (-1) on the other with a
            # double crossover, causing the nick to be joined on re-save.
            scaf_seg[vhNum].extend(segmentBounds(scaf, vhNum))
            # read scaffold xovers
            for i in range(len(scaf)):
                fiveVH, fiveIdx, threeVH, threeIdx = scaf[i]
                if fiveVH == -1 and threeVH == -1:
                    continue  # null base
                if is3primeXover(StrandType.Scaffold, vhNum, i, threeVH, threeIdx):
                    scaf_xo[vhNum].append((i, threeVH, threeIdx))
            assert (len(scaf_seg[vhNum]) % 2 == 0)
            # install scaffold segments
            for i in range(0, len(scaf_seg[vhNum]), 2):
                lowIdx = scaf_seg[vhNum][i]
                highIdx = scaf_seg[vhNum][i+1]
                scafStrandSet.createStrand(lowIdx, highIdx, useUndoStack=False)
            stap_seg[vhNum].extend(segmentBounds(stap, vhNum))
            # read staple xovers
            for i in range(len(stap)):
                fiveVH, fiveIdx, threeVH, threeIdx = stap[i]
                if fiveVH == -1 and threeVH == -1:
                    continue  # null base
                if is3primeXover(StrandType.Staple, vhNum, i, threeVH, threeIdx):
                    stap_xo[vhNum].append((i, threeVH, threeIdx))
            assert (len(stap_seg[vhNum]) % 2 == 0)
            # install staple segments
            for i in range(0, len(stap_seg[vhNum]), 2):
                lowIdx = stap_seg[vhNum][i]
                highIdx = stap_seg[vhNum][i+1]
                stapStrandSet.createStrand(lowIdx, highIdx, useUndoStack=False)
    except AssertionError:
        if not cadnano.app().isGui():
            print("Unrecognized file format.")
        else:
            dialogLT.label.setText("Unrecognized file format.")
            standardButton = getattr(QDialogButtonBox, 'StandardButton',
                                     QDialogButtonBox)
            dialogLT.buttonBox.setStandardButtons(standardButton.Ok)
            dialog.exec()
        return

    # INSTALL XOVERS
    for helix in vstrands:
        vhNum = helix['num']
        row = helix['row']
        col = helix['col']
        scaf = helix['scaf']
        stap = helix['stap']
        insertions = helix['loop']
        skips = helix['skip']
        fromVh = part.virtualHelixAtCoord((row, col))
        scafStrandSet = fromVh.scaffoldStrandSet()
        stapStrandSet = fromVh.stapleStrandSet()
        # install scaffold xovers
        for (idx5p, toVhNum, idx3p) in scaf_xo[vhNum]:
            # idx3p is 3' end of strand5p, idx5p is 5' end of strand3p
            strand5p = scafStrandSet.getStrand(idx5p)
            toVh = part.virtualHelixAtCoord(vhNumToCoord[toVhNum])
            strand3p = toVh.scaffoldStrandSet().getStrand(idx3p)
            part.createXover(strand5p, idx5p, strand3p, idx3p, useUndoStack=False)
        # install staple xovers
        for (idx5p, toVhNum, idx3p) in stap_xo[vhNum]:
            # idx3p is 3' end of strand5p, idx5p is 5' end of strand3p
            strand5p = stapStrandSet.getStrand(idx5p)
            toVh = part.virtualHelixAtCoord(vhNumToCoord[toVhNum])
            strand3p = toVh.stapleStrandSet().getStrand(idx3p)
            part.createXover(strand5p, idx5p, strand3p, idx3p, useUndoStack=False)

    # SET DEFAULT COLOR
    for oligo in part.oligos():
        if oligo.isStaple():
            defaultColor = styles.DEFAULT_STAP_COLOR
        else:
            defaultColor = styles.DEFAULT_SCAF_COLOR
        oligo.applyColor(defaultColor, useUndoStack=False)

    # COLORS, INSERTIONS, SKIPS
    for helix in vstrands:
        vhNum = helix['num']
        row = helix['row']
        col = helix['col']
        scaf = helix['scaf']
        stap = helix['stap']
        insertions = helix['loop']
        skips = helix['skip']
        vh = part.virtualHelixAtCoord((row, col))
        scafStrandSet = vh.scaffoldStrandSet()
        stapStrandSet = vh.stapleStrandSet()
        # install insertions and skips
        for baseIdx in range(len(stap)):
            sumOfInsertSkip = insertions[baseIdx] + skips[baseIdx]
            if sumOfInsertSkip != 0:
                scaf_strand = scafStrandSet.getStrand(baseIdx)
                stap_strand = stapStrandSet.getStrand(baseIdx)
                if scaf_strand:
                    scaf_strand.addInsertion(baseIdx, sumOfInsertSkip, useUndoStack=False)
                elif stap_strand:
                    stap_strand.addInsertion(baseIdx, sumOfInsertSkip, useUndoStack=False)
        # end for
        # populate colors
        for baseIdx, colorNumber in helix['stap_colors']:
            color = QColor((colorNumber>>16)&0xFF, (colorNumber>>8)&0xFF, colorNumber&0xFF).name()
            strand = stapStrandSet.getStrand(baseIdx)
            strand.oligo().applyColor(color, useUndoStack=False)

    # Restore scaffold colors after all scaffold crossovers have established
    # the final oligos. Files created by older cadnano versions omit this
    # extension field and retain the default scaffold color above.
    for record in obj.get('scaffold_colors', []):
        if not isinstance(record, dict):
            continue
        try:
            vhNum = int(record['start_vh'])
            startIdx = int(record['start_idx'])
            color = record['color']
        except (KeyError, TypeError, ValueError):
            continue
        if not isinstance(color, str) or len(color) != 7 or \
                not color.startswith('#'):
            continue
        try:
            int(color[1:], 16)
        except ValueError:
            continue
        vh = part.virtualHelix(vhNum)
        if vh is None:
            continue
        strand = vh.scaffoldStrandSet().getStrand(startIdx)
        if strand is None:
            continue
        strand.oligo().applyColor(color.lower(), useUndoStack=False)

    # Restore applied scaffold sequences only after the complete strand and
    # crossover topology exists. Unknown or stale locations are ignored so
    # older and externally edited files remain loadable.
    for record in obj.get('scaffold_sequences', []):
        if not isinstance(record, dict):
            continue
        try:
            vhNum = int(record['start_vh'])
            startIdx = int(record['start_idx'])
            sequence = record['sequence']
        except (KeyError, TypeError, ValueError):
            continue
        if not isinstance(sequence, str) or not sequence.strip():
            continue
        vh = part.virtualHelix(vhNum)
        if vh is None:
            continue
        strand = vh.scaffoldStrandSet().getStrand(startIdx)
        if strand is None or strand.idx5Prime() != startIdx:
            continue
        oligo = strand.oligo()
        if oligo.strand5p() is not strand:
            continue
        oligo.applySequence(sequence, useUndoStack=False)
    guidedRegions = obj.get('guided_duplex_regions', None)
    if isinstance(guidedRegions, list):
        restored = defaultdict(list)
        for record in guidedRegions:
            if not isinstance(record, dict):
                continue
            try:
                coord = (int(record['row']), int(record['col']))
                intervals = record['intervals']
            except (KeyError, TypeError, ValueError):
                continue
            for interval in intervals if isinstance(intervals, list) else ():
                try:
                    low, high = int(interval[0]), int(interval[1])
                except (TypeError, ValueError, IndexError):
                    continue
                restored[coord].append((low, high))
        part.setGuidedDuplexRegions(restored)
    scaffoldOnlyRegions = obj.get('scaffold_only_regions', None)
    if isinstance(scaffoldOnlyRegions, list):
        restored = defaultdict(list)
        for record in scaffoldOnlyRegions:
            if not isinstance(record, dict):
                continue
            try:
                coord = (int(record['row']), int(record['col']))
                intervals = record['intervals']
            except (KeyError, TypeError, ValueError):
                continue
            for interval in intervals if isinstance(intervals, list) else ():
                try:
                    low, high = int(interval[0]), int(interval[1])
                except (TypeError, ValueError, IndexError):
                    continue
                restored[coord].append((low, high))
        part.setScaffoldOnlyRegions(restored)
    part._autobreakStaplesApplied = bool(
        obj.get('autobreak_staples_applied', False))
    return part

def segmentBounds(strand, vhNum):
    """Return low/high bounds of maximal locally connected base runs."""
    occupied = [connection != [-1, -1, -1, -1]
                for connection in strand]

    def connected(first, second):
        if first < 0 or second >= len(strand):
            return False
        if not occupied[first] or not occupied[second]:
            return False
        firstConnection = strand[first]
        secondConnection = strand[second]
        firstTargets = ((firstConnection[0], firstConnection[1]),
                        (firstConnection[2], firstConnection[3]))
        secondTargets = ((secondConnection[0], secondConnection[1]),
                         (secondConnection[2], secondConnection[3]))
        return ((vhNum, second) in firstTargets and
                (vhNum, first) in secondTargets)

    bounds = []
    for index, isOccupied in enumerate(occupied):
        if not isOccupied:
            continue
        if not connected(index - 1, index):
            bounds.append(index)
        if not connected(index, index + 1):
            bounds.append(index)
    return bounds


def isSegmentStartOrEnd(strandType, vhNum, baseIdx, fiveVH, fiveIdx, threeVH, threeIdx):
    """Returns True if the base is a breakpoint or crossover."""
    if strandType == StrandType.Scaffold:
        offset = 1
    else:
        offset = -1
    if (fiveVH == vhNum and threeVH != vhNum):
        return True
    if (fiveVH != vhNum and threeVH == vhNum):
        return True
    if (vhNum % 2 == 0 and fiveVH == vhNum and fiveIdx != baseIdx-offset):
        return True
    if (vhNum % 2 == 0 and threeVH == vhNum and threeIdx != baseIdx+offset):
        return True
    if (vhNum % 2 == 1 and fiveVH == vhNum and fiveIdx != baseIdx+offset):
        return True
    if (vhNum % 2 == 1 and threeVH == vhNum and threeIdx != baseIdx-offset):
        return True
    if (fiveVH == -1 and threeVH != -1):
        return True
    if (fiveVH != -1 and threeVH == -1):
        return True
    return False

def is3primeXover(strandType, vhNum, baseIdx, threeVH, threeIdx):
    """Returns True of the threeVH doesn't match vhNum, or threeIdx
    is not a natural neighbor of baseIdx."""
    if threeVH == -1:
        return False
    if vhNum != threeVH:
        return True
    if strandType == StrandType.Scaffold:
        offset = 1
    else:
        offset = -1
    if (vhNum % 2 == 0 and threeVH == vhNum and threeIdx != baseIdx+offset):
        return True
    if (vhNum % 2 == 1 and threeVH == vhNum and threeIdx != baseIdx-offset):
        return True
    return False
