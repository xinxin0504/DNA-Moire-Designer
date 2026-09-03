"""Print model/export chain composition for one Hybrid design."""

import json
import math
import sys

from PyQt6.QtWidgets import QApplication

import cadnano2.cadnano as cadnano


qtApp = QApplication.instance() or QApplication([])
app = cadnano.initAppWithoutGui([])
app.prefs.honeycombRows = app.prefs.honeycombCols = 50
app.prefs.squareRows = app.prefs.squareCols = 50
app.prefs.honeycombSteps = app.prefs.squareSteps = 20


class DummySignal:
    def emit(self, *unused_args):
        pass


app.documentWasCreatedSignal = DummySignal()
from cadnano2.model.document import Document
from cadnano2.model.io.decoder import decode
from cadnano2.model.io.oxdnaexport import _collect, _number_records


document = Document()
with open(sys.argv[1], encoding="utf-8") as source:
    decode(document, source.read())
records, strands, unused_assigned, unused_residual = _collect(document, 2.8)
ordered = _number_records(strands)
summary = []
for strand_id, info in enumerate(strands, 1):
    runs = []
    for rec in info["records"]:
        lattice = "H" if rec["part"]._step == 21 else "S"
        marker = (lattice, rec["strand"].virtualHelix().number())
        if not runs or tuple(runs[-1][:2]) != marker:
            runs.append([marker[0], marker[1], 1])
        else:
            runs[-1][2] += 1
    summary.append({
        "id": strand_id,
        "type": "staple" if info["oligo"].isStaple() else "scaffold",
        "loop": info["loop"], "length": len(info["records"]),
        "runs": runs})
print(json.dumps([item for item in summary
                  if item["type"] == "scaffold" or
                  any(run[0] == "S" for run in item["runs"])],
                 indent=2))

long_bonds = []
for first, second in zip(ordered, ordered[1:]):
    if first["five_neighbor"] != second["global_index"]:
        continue
    distance = math.sqrt(sum(
        (first["pos"][axis] - second["pos"][axis]) ** 2
        for axis in range(3)))
    if distance > 1.0:
        long_bonds.append({
            "distance": round(distance, 4),
            "delta": [round(second["pos"][axis] - first["pos"][axis], 4)
                      for axis in range(3)],
            "from": ["H" if first["part"]._step == 21 else "S",
                     first["strand"].virtualHelix().number(), first["idx"]],
            "to": ["H" if second["part"]._step == 21 else "S",
                   second["strand"].virtualHelix().number(), second["idx"]]})
print("longest bonds")
print(json.dumps(sorted(long_bonds, key=lambda item: -item["distance"])[:30],
                 indent=2))

if len(sys.argv) > 2:
    with open(sys.argv[2], encoding="utf-8") as source:
        top_lines = source.read().splitlines()
    count, strand_count = map(int, top_lines[0].split())
    assert count == len(ordered)
    assert strand_count == len(strands)
    for index, (line, rec) in enumerate(zip(top_lines[1:], ordered)):
        sid, base, n3, n5 = line.split()
        assert int(sid) == rec["strand_id"], (index, line, rec)
        assert int(n3) == rec["three_neighbor"], (index, line, rec)
        assert int(n5) == rec["five_neighbor"], (index, line, rec)
    print("topology matches model export records")
