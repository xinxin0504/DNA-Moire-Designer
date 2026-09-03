"""Smoke test: Autobreak must preserve every pre-existing staple nick."""

import json

import cadnano2.cadnano as cadnano


def main():
    from PyQt6.QtWidgets import QApplication
    qtApp = QApplication.instance() or QApplication([])
    app = cadnano.initAppWithoutGui([])
    app.prefs.squareRows = 50
    app.prefs.squareCols = 50
    app.prefs.squareSteps = 4

    class DummySignal:
        def emit(self, *unused_args):
            pass

    app.documentWasCreatedSignal = DummySignal()
    from cadnano2.model.document import Document
    from cadnano2.model.parts.part import _existingStapleNickBoundaries

    document = Document()
    part = document.addGuidedPart('square', 6, 14, 4)
    part.createVirtualHelix(1, 2, useUndoStack=False)
    vh = part.virtualHelixAtCoord((1, 2))
    strandSet = vh.stapleStrandSet()
    strandSet.createStrand(0, 100, useUndoStack=False)
    strand = strandSet.getStrand(40)
    strandSet.splitStrand(strand, 40, useUndoStack=False)

    before = _existingStapleNickBoundaries(part)
    result = part.autoBreakStaples()
    after = _existingStapleNickBoundaries(part)
    secondResult = part.autoBreakStaples()
    afterSecond = _existingStapleNickBoundaries(part)
    assert secondResult['already_applied']
    assert afterSecond == after
    print(json.dumps({
        'result': result,
        'second_result': secondResult,
        'before': sorted(before),
        'after': sorted(after),
        'after_second': sorted(afterSecond),
        'preserved': before.issubset(after),
        'added': sorted(after - before),
    }, indent=2))


if __name__ == '__main__':
    main()
