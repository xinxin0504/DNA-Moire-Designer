"""Regression: Hybrid structure export must keep lattice bodies separated."""

import json
import sys

from PyQt6.QtWidgets import QApplication

import cadnano2.cadnano as cadnano


path = sys.argv[1]
qtApp = QApplication.instance() or QApplication([])
app = cadnano.initAppWithoutGui([])
app.prefs.honeycombRows = 50
app.prefs.honeycombCols = 50
app.prefs.honeycombSteps = 20
app.prefs.squareRows = 50
app.prefs.squareCols = 50
app.prefs.squareSteps = 20


class DummySignal:
    def emit(self, *unused_args):
        pass


app.documentWasCreatedSignal = DummySignal()

from cadnano2.model.document import Document
from cadnano2.model.io.decoder import decode
from cadnano2.model.io.oxdnaexport import (
    _collect, _dat_text, _number_records, OXDNA_LENGTH_NM)


document = Document()
with open(path, encoding="utf-8") as source:
    decode(document, source.read())
records, strands, unused_assigned, unused_residual = _collect(
    document, 2.8)
bounds = {}
for lattice in ("honeycomb", "square"):
    latticeRecords = [
        record for record in records
        if ("honeycomb" if record["part"]._step == 21 else "square") ==
        lattice]
    bounds[lattice] = (
        min(record["pos"][0] for record in latticeRecords),
        max(record["pos"][0] for record in latticeRecords))

gap = bounds["square"][0] - bounds["honeycomb"][1]
assert gap > 2.8 / OXDNA_LENGTH_NM, (bounds, gap)
ordered = _number_records(strands)
dat = _dat_text(ordered)
box = tuple(float(value) for value in dat.splitlines()[1].split()[2:])
by_index = dict((record["global_index"], record) for record in ordered)
for record in ordered:
    neighbor = by_index.get(record["five_neighbor"])
    if neighbor is None or neighbor["part"] is record["part"]:
        continue
    for axis in range(3):
        component = abs(neighbor["pos"][axis] - record["pos"][axis])
        assert component < box[axis] / 4.0, (
            axis, component, box[axis], record["global_index"])
print(json.dumps({"bounds": bounds, "gap_oxdna": gap,
                  "gap_nm": gap * OXDNA_LENGTH_NM,
                  "periodic_box": box}, indent=2))
