#!/usr/bin/env python
# encoding: utf-8

from .parts.honeycombpart import HoneycombPart
from .parts.squarepart import SquarePart
from .parts.part import Part
from .strand import Strand
from operator import itemgetter
import cadnano2.util as util
import cadnano2.cadnano as cadnano
util.qtWrapImport('QtCore', globals(), ['pyqtSignal', 'QObject'])
util.qtWrapImport('QtGui', globals(), ['QUndoCommand', 'QUndoStack'])


class Document(QObject):
    """
    The Document class is the root of the model. It has two main purposes:
    1. Serve as the parent all Part objects within the model.
    2. Track all sub-model actions on its undoStack.
    """
    def __init__(self):
        super(Document, self).__init__()
        self._undoStack = QUndoStack()
        self._parts = []
        self._assemblies = []
        self._controller = None
        self._selectedPart = None
        self._hybridMode = False
        self._athenaMetadata = None
        self._curvedMetadata = None
        self._twistBendMetadata = None
        # the dictionary maintains what is selected
        self._selectionDict = {}
        # the added list is what was recently selected or deselected
        self._selectedChangedDict = {}
        cadnano.app().documentWasCreatedSignal.emit(self)

    ### SIGNALS ###
    documentPartAddedSignal = pyqtSignal(object, QObject)  # doc, part

    # dict of tuples of objects using the reference as the key,
    # and the value is a tuple with meta data
    # in the case of strands the metadata would be which endpoints of selected
    # e.g. { objectRef: (value0, value1),  ...}
    documentSelectedChangedSignal = pyqtSignal(dict)  # tuples of items + data
    documentSelectionFilterChangedSignal = pyqtSignal(list) # doc, filterlist

    documentViewResetSignal = pyqtSignal(QObject)
    documentClearSelectionsSignal = pyqtSignal(QObject)
    documentHybridConnectionChangedSignal = pyqtSignal(QObject)


    ### SLOTS ###

    ### ACCESSORS ###
    def undoStack(self):
        """
        This is the actual undoStack to use for all commands. Any children
        needing to perform commands should just ask their parent for the
        undoStack, and eventually the request will get here.
        """
        return self._undoStack

    def parts(self):
        """Returns a list of parts associated with the document."""
        return self._parts

    def isHybrid(self):
        return self._hybridMode

    def setHybrid(self, enabled):
        self._hybridMode = bool(enabled)

    def athenaMetadata(self):
        """Return the self-contained ATHENA target-geometry metadata."""
        return self._athenaMetadata

    def setAthenaMetadata(self, metadata):
        self._athenaMetadata = (dict(metadata)
                                if isinstance(metadata, dict) else None)

    def curvedMetadata(self):
        """Return the self-contained DNAxiS curved-geometry metadata."""
        return self._curvedMetadata

    def setCurvedMetadata(self, metadata):
        self._curvedMetadata = (dict(metadata)
                                if isinstance(metadata, dict) else None)

    def twistBendMetadata(self):
        """Return saved Twist and Bend regions and preview parameters."""
        return self._twistBendMetadata

    def setTwistBendMetadata(self, metadata):
        self._twistBendMetadata = (dict(metadata)
                                   if isinstance(metadata, dict) else None)

    def partsByLattice(self):
        """Return the document parts keyed by their stable lattice label."""
        parts = {}
        for part in self._parts:
            if part._step == 21:
                parts['honeycomb'] = part
            elif part._step == 32:
                parts['square'] = part
        return parts

    def oligos(self):
        """Return every physical oligo once, including cross-part oligos."""
        oligos = set()
        for part in self._parts:
            oligos.update(part.oligos())
        return oligos

    def hybridConnectionRecords(self):
        """Describe every 3'-to-5' connection crossing lattice parts."""
        records = []
        for part in self._parts:
            for vh in part.getVirtualHelices():
                for strandSet in (vh.scaffoldStrandSet(),
                                  vh.stapleStrandSet()):
                    for strand in strandSet:
                        connected = strand.connection3p()
                        if connected is None or connected.part() is part:
                            continue
                        records.append({
                            'strand_type': ('scaffold' if strand.isScaffold()
                                            else 'staple'),
                            'three_prime': self._hybridEndpointRecord(strand),
                            'five_prime': self._hybridEndpointRecord(connected)
                        })
        return records

    def _hybridEndpointRecord(self, strand):
        part = strand.part()
        lattice = ('honeycomb' if part._step == 21 else 'square')
        index = (strand.idx3Prime() if strand.connection3p() is not None and
                 strand.connection3p().part() is not part
                 else strand.idx5Prime())
        return {'lattice': lattice,
                'helix': strand.virtualHelix().number(),
                'index': index}

    def createHybridConnection(self, firstStrand, firstEnd,
                               secondStrand, secondEnd,
                               useUndoStack=True):
        """Connect complementary free endpoints in different lattice parts."""
        if firstStrand is None or secondStrand is None:
            return False, 'Both endpoints must contain a strand.'
        if firstStrand.document() is not self or \
                secondStrand.document() is not self:
            return False, 'Both endpoints must belong to this design.'
        if firstStrand.part() is secondStrand.part():
            return False, 'Choose endpoints from different lattice panels.'
        if firstStrand.part()._step == secondStrand.part()._step:
            return False, 'Hybrid connections require different lattices.'
        if firstStrand.strandType() != secondStrand.strandType():
            return False, 'Scaffold can only connect to scaffold, and staple to staple.'
        if set((firstEnd, secondEnd)) != set(('3p', '5p')):
            return False, "Choose one free 3' endpoint and one free 5' endpoint."
        strand5p = firstStrand if firstEnd == '3p' else secondStrand
        strand3p = secondStrand if secondEnd == '5p' else firstStrand
        if strand5p.connection3p() is not None:
            return False, "The selected 3' endpoint is already connected."
        if strand3p.connection5p() is not None:
            return False, "The selected 5' endpoint is already connected."
        strand5p.part().createXover(
            strand5p, strand5p.idx3Prime(),
            strand3p, strand3p.idx5Prime(),
            useUndoStack=useUndoStack)
        created = strand5p.connection3p() is strand3p
        if created:
            strand5p.part().partStrandChangedSignal.emit(
                strand5p.part(), strand5p.virtualHelix())
            strand3p.part().partStrandChangedSignal.emit(
                strand3p.part(), strand3p.virtualHelix())
            self.documentHybridConnectionChangedSignal.emit(self)
            return True, ''
        return False, 'The endpoints could not be connected.'

    def removeHybridConnection(self, strand, end, useUndoStack=True):
        """Remove the cross-part connection attached to one endpoint."""
        connected = (strand.connection3p() if end == '3p'
                     else strand.connection5p())
        if connected is None or connected.part() is strand.part():
            return False
        strand5p = strand if end == '3p' else connected
        strand3p = connected if end == '3p' else strand
        strand5p.part().removeXover(
            strand5p, strand3p, useUndoStack=useUndoStack)
        removed = strand5p.connection3p() is not strand3p
        if removed:
            strand5p.part().partStrandChangedSignal.emit(
                strand5p.part(), strand5p.virtualHelix())
            strand3p.part().partStrandChangedSignal.emit(
                strand3p.part(), strand3p.virtualHelix())
            self.documentHybridConnectionChangedSignal.emit(self)
        return removed

    def getInputSequenceRows(self):
        return self._sequenceRows(False)

    def getOutputSequenceRows(self):
        return self._sequenceRows(True)

    def _sequenceRows(self, staples):
        categoryOrder = {'honeycomb': 0, 'square': 1, 'hybrid': 2}
        rows = []
        for oligo in self.oligos():
            strand5p = oligo.strand5p()
            if strand5p is None or oligo.isLoop() or \
                    oligo.isStaple() != staples:
                continue
            if not staples:
                sequence = oligo.sequence()
                if not sequence or not sequence.strip():
                    continue
            record = oligo.sequenceRecord()
            rows.append((categoryOrder[oligo.latticeCategory()], record))
        from .parts.part import _sequenceRowSortKey
        rows.sort(key=lambda item: (item[0],) + _sequenceRowSortKey(item[1]))
        return [record for unused_category, record in rows]

    def getScaffoldSequenceTemplateRows(self):
        categoryOrder = {'honeycomb': 0, 'square': 1, 'hybrid': 2}
        rows = []
        for oligo in self.oligos():
            strand5p = oligo.strand5p()
            if strand5p is None or oligo.isLoop() or oligo.isStaple():
                continue
            start, end = oligo.sequenceEndpoints()
            record = (start, end, '', oligo.actualLength(), oligo.color())
            rows.append((categoryOrder[oligo.latticeCategory()], record))
        from .parts.part import _sequencePositionSortKey
        rows.sort(key=lambda item: ((item[0],) +
                                    _sequencePositionSortKey(item[1])))
        return [record for unused_category, record in rows]

    def assemblies(self):
        """Returns a list of assemblies associated with the document."""
        return self._assemblies

    ### PUBLIC METHODS FOR QUERYING THE MODEL ###
    def selectedPart(self):
        return self._selectedPart

    def addToSelection(self, obj, value):
        self._selectionDict[obj] = value
        self._selectedChangedDict[obj] = value
    # end def

    def removeFromSelection(self, obj):
        if obj in self._selectionDict:
            del self._selectionDict[obj]
            self._selectedChangedDict[obj] = (False, False)
            return True
        else:
            return False
    # end def

    def clearSelections(self):
        """
        Only clear the dictionary
        """
        self._selectionDict = {}
    # end def

    def addStrandToSelection(self, strand, value):
        sS = strand.strandSet()
        if sS in self._selectionDict:
            self._selectionDict[sS][strand] = value
        else:
            self._selectionDict[sS] = {strand: value}
        self._selectedChangedDict[strand] = value
    # end def

    def removeStrandFromSelection(self, strand):
        sS = strand.strandSet()
        if sS in self._selectionDict:
            temp = self._selectionDict[sS]
            if strand in temp:
                del temp[strand]
                if len(temp) == 0:
                    del self._selectionDict[sS]
                self._selectedChangedDict[strand] = (False, False)
                return True
            else:
                return False
        else:
            return False
    # end def

    def selectionDict(self):
        return self._selectionDict
    # end def
    
    def selectedOligos(self):
        """
        as long as one endpoint of a strand is in the selection, then the oligo
        is considered selected
        """
        sDict = self._selectionDict
        selectedOs = set()
        for sS in sDict.keys():
            for strand in sS:
                 selectedOs.add(strand.oligo())
            # end for
        # end for
        return selectedOs if len(selectedOs) > 0 else None
    #end def
    
    def clearAllSelected(self):
        self._selectionDict = {}
        # the added list is what was recently selected or deselected
        self._selectedChangedDict = {}
        self.documentClearSelectionsSignal.emit(self)
    # end def

    def isModelSelected(self, obj):
        return obj in self._selectionDict
    # end def

    def isModelStrandSelected(self, strand):
        sS = strand.strandSet()
        if sS in self._selectionDict:
            if strand in self._selectionDict[sS]:
                return True
            else:
                return False
        else:
            return False
    # end def

    def getSelectedValue(self, obj):
        """
        obj is an objects to look up
        it is prevetted to be in the dictionary
        """
        return self._selectionDict[obj]

    def getSelectedStrandValue(self, strand):
        """
        strand is an objects to look up
        it is prevetted to be in the dictionary
        """
        return self._selectionDict[strand.strandSet()][strand]
    # end def

    def sortedSelectedStrands(self, strandSet):
        # outList = self._selectionDict[strandSet].keys()
        # outList.sort(key=Strand.lowIdx)
        outList = list(self._selectionDict[strandSet].items())
        getLowIdx = lambda x: Strand.lowIdx(itemgetter(0)(x))
        outList.sort(key=getLowIdx)
        return outList
    # end def

    def determineStrandSetBounds(self, selectedStrandList, strandSet):
        minLowDelta = strandSet.partMaxBaseIdx()
        minHighDelta = strandSet.partMaxBaseIdx()  # init the return values
        sSDict = self._selectionDict[strandSet]
        # get the StrandSet index of the first item in the list
        sSIdx = strandSet._findIndexOfRangeFor(selectedStrandList[0][0])[2]
        sSList = strandSet._strandList
        lenSSList = len(sSList)
        maxSSIdx = lenSSList - 1
        i = 0
        for strand, value in selectedStrandList:
            while strand != sSList[sSIdx]:
                # incase there are gaps due to double xovers
                sSIdx += 1
            # end while
            idxL, idxH = strand.idxs()
            if value[0]:    # the end is selected
                if sSIdx > 0:
                    lowNeighbor = sSList[sSIdx - 1]
                    if lowNeighbor in sSDict:
                        valueN = sSDict[lowNeighbor]
                        # we only care if the low neighbor is not selected
                        temp = minLowDelta if valueN[1] \
                                        else idxL - lowNeighbor.highIdx() - 1
                    # end if
                    else:  # not selected
                        temp = idxL - lowNeighbor.highIdx() - 1
                    # end else
                else:
                    temp = idxL - 0
                # end else
                if temp < minLowDelta:
                    minLowDelta = temp
                # end if
                # check the other end of the strand
                if not value[1]:
                    temp = idxH - idxL - 1
                    if temp < minHighDelta:
                        minHighDelta = temp
            # end if
            if value[1]:
                if sSIdx < maxSSIdx:
                    highNeighbor = sSList[sSIdx + 1]
                    if highNeighbor in sSDict:
                        valueN = sSDict[highNeighbor]
                        # we only care if the low neighbor is not selected
                        temp = minHighDelta if valueN[0] \
                                        else highNeighbor.lowIdx() - idxH - 1
                    # end if
                    else:  # not selected
                        temp = highNeighbor.lowIdx() - idxH - 1
                    # end else
                else:
                    temp = strandSet.partMaxBaseIdx() - idxH
                # end else
                if temp < minHighDelta:
                    minHighDelta = temp
                # end if
                # check the other end of the strand
                if not value[0]:
                    temp = idxH - idxL - 1
                    if temp < minLowDelta:
                        minLowDelta = temp
            # end if
            # increment counter
            sSIdx += 1
        # end for
        return (minLowDelta, minHighDelta)
    # end def

    def getSelectionBounds(self):
        minLowDelta = -1
        minHighDelta = -1
        for strandSet in self._selectionDict.keys():
            selectedList = self.sortedSelectedStrands(strandSet)
            tempLow, tempHigh = self.determineStrandSetBounds(
                                                    selectedList, strandSet)
            if tempLow < minLowDelta or minLowDelta < 0:
                minLowDelta = tempLow
            if tempHigh < minHighDelta or minHighDelta < 0:
                minHighDelta = tempHigh
        # end for Mark train bus to metro
        return (minLowDelta, minHighDelta)
    # end def

    # def operateOnStrandSelection(self, method, arg, both=False):
    #     pass
    # # end def

    def deleteSelection(self, useUndoStack=True):
        """
        Delete selected strands. First iterates through all selected strands
        and extracts refs to xovers and strands. Next, calls removeXover
        on xoverlist as part of its own macroed command for isoluation
        purposes. Finally, calls removeStrand on all strands that were 
        fully selected (low and high), or had at least one non-xover
        endpoint selected.
        """
        xoList = []
        strandDict = {}
        for strandSetDict in list(self._selectionDict.values()):
            for strand, selected in list(strandSetDict.items()):
                part = strand.virtualHelix().part()
                idxL, idxH = strand.idxs()
                strand5p = strand.connection5p()
                strand3p = strand.connection3p()
                # both ends are selected
                strandDict[strand] = selected[0] and selected[1]

                # only look at 3' ends to handle xover deletion
                sel3p = selected[0] if idxL == strand.idx3Prime() else selected[1]
                if sel3p:  # is idx3p selected?
                    if strand3p:  # is there an xover
                        xoList.append((part, strand, strand3p, useUndoStack))
                    else:  # idx3p is a selected endpoint
                        strandDict[strand] = True
                else:
                    if not strand5p:  # idx5p is a selected endpoint
                        strandDict[strand] = True


        if useUndoStack and xoList:
            self.undoStack().beginMacro("Delete xovers")
        for part, strand, strand3p, useUndo in xoList:
            Part.removeXover(part, strand, strand3p, useUndo)
            self.removeStrandFromSelection(strand)
            self.removeStrandFromSelection(strand3p)
        self._selectionDict = {}
        self.documentClearSelectionsSignal.emit(self)
        if useUndoStack:
            if xoList: # end xover macro if it was started
                self.undoStack().endMacro()
            if True in list(strandDict.values()):
                self.undoStack().beginMacro("Delete selection")
            else:
                return  # nothing left to do
        for strand, delete in list(strandDict.items()):
            if delete:
                strand.strandSet().removeStrand(strand)
        if useUndoStack:
            self.undoStack().endMacro()

    def paintSelection(self, scafColor, stapColor, useUndoStack=True):
        """Delete xovers if present. Otherwise delete everything."""
        scafOligos = {}
        stapOligos = {}
        for strandSetDict in list(self._selectionDict.values()):
            for strand, value in list(strandSetDict.items()):
                if strand.isScaffold():
                    scafOligos[strand.oligo()] = True
                else:
                    stapOligos[strand.oligo()] = True

        if useUndoStack:
            self.undoStack().beginMacro("Paint strands")
        for olg in list(scafOligos.keys()):
            olg.applyColor(scafColor)
        for olg in list(stapOligos.keys()):
            olg.applyColor(stapColor)
        if useUndoStack:
            self.undoStack().endMacro()

    def resizeSelection(self, delta, useUndoStack=True):
        """
        Moves the selected idxs by delta by first iterating over all strands
        to calculate new idxs (method will return if snap-to behavior would
        create illegal state), then applying a resize command to each strand.
        """
        resizeList = []

        # calculate new idxs
        for strandSetDict in self._selectionDict.values():
            for strand, selected in strandSetDict.items():
                part = strand.virtualHelix().part()
                idxL, idxH = strand.idxs()
                newL, newH = strand.idxs()
                deltaL = deltaH = delta

                # process xovers to get revised delta
                if selected[0] and strand.connectionLow():
                    newL = part.xoverSnapTo(strand, idxL, delta)
                    if newL == None:
                        return
                    deltaH = newL-idxL
                if selected[1] and strand.connectionHigh():
                    newH = part.xoverSnapTo(strand, idxH, delta)
                    if newH == None:
                        return
                    deltaL = newH-idxH

                # process endpoints
                if selected[0] and not strand.connectionLow():
                    newL = idxL + deltaL
                if selected[1] and not strand.connectionHigh():
                    newH = idxH + deltaH

                if newL > newH:  # check for illegal state
                    return

                resizeList.append((strand, newL, newH))
            # end for
        # end for

        # execute the resize commands
        if useUndoStack:
            self.undoStack().beginMacro("Resize Selection")

        for strand, idxL, idxH in resizeList:
            Strand.resize(strand, (idxL, idxH), useUndoStack)

        if useUndoStack:
            self.undoStack().endMacro()
    # end def

    def updateSelection(self):
        """
        do it this way in the future when we have
        a better signaling architecture between views
        """
        # self.documentSelectedChangedSignal.emit(self._selectedChangedDict)
        """
        For now, individual objects need to emit signals
        """
        for obj, value in self._selectedChangedDict.items():
            obj.selectedChangedSignal.emit(obj, value)
        # end for
        self._selectedChangedDict = {}
        # for sS in self._selectionDict:
        #     print self.sortedSelectedStrands(sS)
    # end def

    def resetViews(self):
        # This is a fast way to clear selections and the views.
        # We could manually deselect each item from the Dict, but we'll just
        # let them be garbage collect
        # the dictionary maintains what is selected
        self._selectionDict = {}
        # the added list is what was recently selected or deselected
        self._selectedChangedDict = {}
        self.documentViewResetSignal.emit(self)
    # end def

    ### PUBLIC METHODS FOR EDITING THE MODEL ###
    def addHoneycombPart(self):
        """
        Create and store a new DNAPart and instance, and return the instance.
        """
        dnapart = None
        if len(self._parts) == 0:
            dnapart = HoneycombPart(document=self)
            self._addPart(dnapart)
        return dnapart

    def addSquarePart(self):
        """
        Create and store a new DNAPart and instance, and return the instance.
        """
        dnapart = None
        if len(self._parts) == 0:
            dnapart = SquarePart(document=self)
            self._addPart(dnapart)
        return dnapart

    def addGuidedPart(self, lattice, maxRow, maxCol, maxSteps):
        """Create a correctly sized single-lattice part for Guided Design."""
        if self._parts:
            return None
        partClass = HoneycombPart if lattice == 'honeycomb' else SquarePart
        dnapart = partClass(document=self, maxRow=maxRow, maxCol=maxCol,
                            maxSteps=maxSteps)
        self._addPart(dnapart)
        return dnapart

    def addHybridParts(self):
        """Create one honeycomb and one square part in the same document."""
        if self._parts:
            return None
        self._hybridMode = True
        honeycomb = HoneycombPart(document=self)
        square = SquarePart(document=self)
        self._addPart(honeycomb)
        self._addPart(square)
        self.setSelectedPart(honeycomb)
        return honeycomb, square

    def removeAllParts(self):
        """Used to reset the document. Not undoable."""
        self.documentClearSelectionsSignal.emit(self)
        for part in list(self._parts):
            part.remove(useUndoStack=False)
        self._hybridMode = False
    # end def

    def removePart(self, part):
        self.documentClearSelectionsSignal.emit(self)
        self._parts.remove(part)
        if self._selectedPart is part:
            self._selectedPart = self._parts[-1] if self._parts else None
        

    ### PUBLIC SUPPORT METHODS ###
    def controller(self):
        return self._controller

    def setController(self, controller):
        """Called by DocumentController setDocument method."""
        self._controller = controller
    # end def

    def setSelectedPart(self, newPart):
        if self._selectedPart == newPart:
            return
        self._selectedPart = newPart
    # end def

    ### PRIVATE SUPPORT METHODS ###
    def _addPart(self, part, useUndoStack=True):
        """Add part to the document via AddPartCommand."""
        c = self.AddPartCommand(self, part)
        util.execCommandList(
                        self, [c], desc="Add part", useUndoStack=useUndoStack)
        return c.part()
    # end def

    ### COMMANDS ###
    class AddPartCommand(QUndoCommand):
        """
        Undo ready command for deleting a part.
        """
        def __init__(self, document, part):
            QUndoCommand.__init__(self)
            self._doc = document
            self._part = part
        # end def

        def part(self):
            return self._part
        # end def

        def redo(self):
            if self._part not in self._doc._parts:
                self._doc._parts.append(self._part)
                self._part.setDocument(self._doc)
                self._doc.setSelectedPart(self._part)
                self._doc.documentPartAddedSignal.emit(self._doc, self._part)
        # end def

        def undo(self):
            self._doc.removePart(self._part)
            self._part.setDocument(None)
            self._doc.setSelectedPart(
                self._doc._parts[-1] if self._doc._parts else None)
            self._part.partRemovedSignal.emit(self._part)
            # self._doc.documentPartAddedSignal.emit(self._doc, self._part)
        # end def
    # end class
# end class
