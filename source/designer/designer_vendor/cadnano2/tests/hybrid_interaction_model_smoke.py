"""Hybrid numbering, cross-lattice toggle, and endpoint-format smoke test."""

from PyQt6.QtWidgets import QApplication

import cadnano2.cadnano as cadnano


qtApp = QApplication.instance() or QApplication([])
app = cadnano.initAppWithoutGui([])
app.prefs.honeycombRows = 50
app.prefs.honeycombCols = 50
app.prefs.honeycombSteps = 2
app.prefs.squareRows = 50
app.prefs.squareCols = 50
app.prefs.squareSteps = 2


class DummySignal:
    def emit(self, *unused_args):
        pass


app.documentWasCreatedSignal = DummySignal()

from cadnano2.model.document import Document


document = Document()
honeycomb, square = document.addHybridParts()
for column in range(6):
    square.createVirtualHelix(0, column, useUndoStack=False)
for column in range(4):
    honeycomb.createVirtualHelix(0, column, useUndoStack=False)

squareNumbers = sorted(vh.number() for vh in square.getVirtualHelices())
honeycombNumbers = sorted(vh.number() for vh in honeycomb.getVirtualHelices())
assert squareNumbers == list(range(6)), squareNumbers
assert honeycombNumbers == list(range(6, 10)), honeycombNumbers
assert not set(squareNumbers).intersection(honeycombNumbers)

squareStrandSet = square.virtualHelixAtCoord((0, 0)).stapleStrandSet()
honeyStrandSet = honeycomb.virtualHelixAtCoord((0, 0)).stapleStrandSet()
squareStrandSet.createStrand(0, 10, useUndoStack=False)
honeyStrandSet.createStrand(0, 10, useUndoStack=False)
squareStrand = squareStrandSet.getStrand(5)
honeyStrand = honeyStrandSet.getStrand(5)

created, error = document.createHybridConnection(
    squareStrand, '3p', honeyStrand, '5p', useUndoStack=False)
assert created, error
assert squareStrand.connection3p() is honeyStrand
assert document.removeHybridConnection(
    squareStrand, '3p', useUndoStack=True)
assert squareStrand.connection3p() is None
assert honeyStrand.connection5p() is None
document.undoStack().undo()
assert squareStrand.connection3p() is honeyStrand
document.undoStack().redo()
assert squareStrand.connection3p() is None

startWithLattice, unused_end = squareStrand.oligo().sequenceEndpoints()
startRaw, unused_end = squareStrand.oligo().sequenceEndpoints(
    includeLattice=False)
assert startWithLattice.startswith('S:')
assert not startRaw.startswith(('H:', 'S:'))

print("hybrid interaction model smoke: OK")
