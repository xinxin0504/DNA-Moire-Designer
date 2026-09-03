"""Regression checks for scaffold-only display metadata and undo behavior."""

import json
import sys

def main():
    from PyQt6.QtWidgets import QApplication
    qtApp = QApplication.instance() or QApplication([])
    import cadnano2.cadnano as cadnano
    app = cadnano.initAppWithoutGui([])
    app.prefs.squareRows = 50
    app.prefs.squareCols = 50
    app.prefs.squareSteps = 4

    class DummySignal:
        def emit(self, *unused_args):
            pass

    app.documentWasCreatedSignal = DummySignal()
    from cadnano2.model.document import Document

    document = Document()
    part = document.addGuidedPart('square', 5, 5, 4)
    part.createVirtualHelix(1, 1, useUndoStack=False)
    vh = part.virtualHelixAtCoord((1, 1))
    strandSet = vh.scaffoldStrandSet()
    strandSet.createStrand(10, 30, useUndoStack=False)
    strand = strandSet.getStrand(20)
    part.addScaffoldOnlyRegions({(1, 1): [(10, 15)]})

    strand.resize((16, 30))
    afterResize = part.scaffoldOnlyRegionRecords()
    part.undoStack().undo()
    afterResizeUndo = part.scaffoldOnlyRegionRecords()
    part.undoStack().redo()
    afterResizeRedo = part.scaffoldOnlyRegionRecords()
    part.undoStack().undo()

    strand = strandSet.getStrand(20)
    strandSet.removeStrand(strand)
    afterRemove = part.scaffoldOnlyRegionRecords()
    part.undoStack().undo()
    afterRemoveUndo = part.scaffoldOnlyRegionRecords()
    part.undoStack().redo()
    afterRemoveRedo = part.scaffoldOnlyRegionRecords()

    print(json.dumps({
        'after_resize': afterResize,
        'after_resize_undo': afterResizeUndo,
        'after_resize_redo': afterResizeRedo,
        'after_remove': afterRemove,
        'after_remove_undo': afterRemoveUndo,
        'after_remove_redo': afterRemoveRedo,
    }, indent=2))


if __name__ == '__main__':
    main()
