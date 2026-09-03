"""Regression: Add scaffolds/staples uses the current strand command API."""

from PyQt6.QtWidgets import QApplication

import cadnano2.cadnano as cadnano


qtApp = QApplication.instance() or QApplication([])
app = cadnano.initAppWithoutGui([])
app.prefs.squareRows = 50
app.prefs.squareCols = 50
app.prefs.squareSteps = 2


class DummySignal:
    def emit(self, *unused_args):
        pass


app.documentWasCreatedSignal = DummySignal()

from cadnano2.model.document import Document
from cadnano2.model.enum import StrandType
from cadnano2.views import styles


document = Document()
part = document.addSquarePart()
part.createVirtualHelix(0, 0, useUndoStack=False)
assert part.autoFillWithoutCrossovers(StrandType.Scaffold) == 1
assert part.autoFillWithoutCrossovers(StrandType.Staple) == 1
staples = [oligo for oligo in part.oligos() if oligo.isStaple()]
assert len(staples) == 1
assert staples[0].color().lower() == styles.AUTOMATIC_STAP_COLOR
print("automatic fill/color smoke: OK")
