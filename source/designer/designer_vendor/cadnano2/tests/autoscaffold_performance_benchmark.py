"""Core AutoCS benchmark for 64-helix, 200-base lattice panels."""

import time

import cadnano2.cadnano as cadnano


def build_part(document, lattice):
    if lattice == 'square':
        part = document.addSquarePart()
    else:
        part = document.addHoneycombPart()
    coords = []
    for rowOffset in range(8):
        columns = range(8) if rowOffset % 2 == 0 else range(7, -1, -1)
        for columnOffset in columns:
            coord = (20 + rowOffset, 20 + columnOffset)
            part.createVirtualHelix(*coord, useUndoStack=False)
            part.virtualHelixAtCoord(coord).scaffoldStrandSet().createStrand(
                0, 199, useUndoStack=False)
            coords.append(coord)
    part.setImportedVHelixOrder(coords, emitSignal=False)
    document.setSelectedPart(part)
    return part


def main():
    from PyQt6.QtWidgets import QApplication
    qtApp = QApplication.instance() or QApplication([])
    app = cadnano.initAppWithoutGui([])
    app.prefs.squareRows = 80
    app.prefs.squareCols = 80
    app.prefs.squareSteps = 7
    app.prefs.honeycombRows = 80
    app.prefs.honeycombCols = 80
    app.prefs.honeycombSteps = 10

    class DummySignal:
        def emit(self, *unused_args):
            pass

    app.documentWasCreatedSignal = DummySignal()
    from cadnano2.model.document import Document

    for lattice in ('square', 'honeycomb'):
        document = Document()
        part = build_part(document, lattice)
        start = time.perf_counter()
        details = part.autoScaffoldCrossovers(
            densityMultiple=1, rebuildExisting=True, returnDetails=True)
        elapsed = time.perf_counter() - start
        loops = sum(
            1 for oligo in part.oligos()
            if not oligo.isStaple() and oligo.isLoop())
        print('%s %.6f seconds, %d xovers, %d loops' % (
            lattice, elapsed, len(part._existingScaffoldCrossoverRecords()),
            loops))
        assert details['hard_density_valid']


if __name__ == '__main__':
    main()
