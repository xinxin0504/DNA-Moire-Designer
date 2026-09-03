"""Regression checks for the compact PDB/mmCIF + oxDNA export."""

import os
import tempfile

from PyQt6.QtWidgets import QApplication

import cadnano2.cadnano as cadnano


qtApp = QApplication.instance() or QApplication([])
app = cadnano.initAppWithoutGui([])
app.prefs.honeycombRows = 50
app.prefs.honeycombCols = 50
app.prefs.honeycombSteps = 2


class DummySignal:
    def emit(self, *unused_args):
        pass


app.documentWasCreatedSignal = DummySignal()

from cadnano2.model.document import Document
from cadnano2.model.io.oxdnaexport import export_structure_bundle

document = Document()
part = document.addGuidedPart('honeycomb', 3, 3, 2)
part.createVirtualHelix(1, 1, useUndoStack=False)
part.virtualHelixAtCoord((1, 1)).scaffoldStrandSet().createStrand(
    0, 5, useUndoStack=False)

with tempfile.TemporaryDirectory() as directory:
    regular = os.path.join(directory, 'regular')
    summary = export_structure_bundle(
        document, regular, 'design', include_oxdna=True)
    assert summary['structure_format'] == 'PDB'
    assert sorted(os.listdir(regular)) == [
        'design.dat', 'design.pdb', 'design.top']

    guided = os.path.join(directory, 'guided')
    summary = export_structure_bundle(
        document, guided, 'design', include_oxdna=False, pdb_atom_limit=1)
    assert summary['structure_format'] == 'mmCIF'
    assert os.listdir(guided) == ['design.cif']

print('compact structure export smoke: OK')
