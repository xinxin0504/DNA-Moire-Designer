from os.path import basename
from copy import deepcopy
from cadnano2.model.enum import StrandType


def _scaffold_records(part):
    """Return persistent scaffold colours and applied 5'-anchored sequences."""
    colors = []
    sequences = []
    for oligo in part.oligos():
        strand = oligo.strand5p()
        if strand is None or oligo.isStaple():
            continue
        start_vh = int(strand.virtualHelix().number())
        start_idx = int(strand.idx5Prime())
        colors.append({
            'start_vh': start_vh,
            'start_idx': start_idx,
            'color': str(oligo.color()),
        })
        sequence = oligo.sequence()
        if sequence and sequence.strip():
            sequences.append({
                'start_vh': start_vh,
                'start_idx': start_idx,
                'sequence': sequence,
            })
    key = lambda record: (record['start_vh'], record['start_idx'])
    return sorted(colors, key=key), sorted(sequences, key=key)


def legacy_dict_from_doc(document, fname, helixOrderList):
    part = document.selectedPart()
    numBases = part.maxBaseIdx()+1

    # iterate through virtualhelix list
    vhList = []
    for row, col in helixOrderList:
        vh = part.virtualHelixAtCoord((row, col))
        # insertions and skips
        insertionDict = part.insertions()[(row, col)]
        insts = [0 for i in range(numBases)]
        skips = [0 for i in range(numBases)]
        for idx, insertion in insertionDict.items():
            if insertion.isSkip():
                skips[idx] = insertion.length()
            else:
                insts[idx] = insertion.length()
        # colors
        stapColors = []
        for strand in vh.stapleStrandSet():
            if strand.connection5p() == None:
                c = str(strand.oligo().color())[1:]  # drop the hash
                stapColors.append([strand.idx5Prime(), int(c, 16)])
        scafColors = set()
        for strand in vh.scaffoldStrandSet():
            if strand.connection5p() == None or \
               (strand == strand.oligo().strand5p() and strand.oligo().isLoop()):
                c = str(strand.oligo().color())[1:]  # drop the hash
                scafColors.add((strand.idx5Prime(), int(c, 16)))

        vhDict = {"row": row,
                  "col": col,
                  "num": vh.number(),
                  "scaf": vh.getLegacyStrandSetArray(StrandType.Scaffold),
                  "stap": vh.getLegacyStrandSetArray(StrandType.Staple),
                  "loop": insts,
                  "skip": skips,
                  "scafLoop": [],
                  "stapLoop": [],
                  "stap_colors": stapColors,
                  "scaf_colors": list(scafColors)}
        vhList.append(vhDict)
    bname = basename(str(fname))
    extensions = deepcopy(getattr(document, '_moireJsonExtensions', {}))
    obj = extensions
    obj.update({"name": bname, "vstrands": vhList,
                "num_bases": numBases,
                "lattice": ("honeycomb" if part.stepSize() == 21
                            else "square")})
    scaffold_colors, scaffold_sequences = _scaffold_records(part)
    if scaffold_colors:
        obj["scaffold_colors"] = scaffold_colors
    if scaffold_sequences:
        # Ordinary Save/Save As intentionally preserves applied sequence data.
        # Standard caDNAno readers safely ignore this extension field.
        obj["scaffold_sequences"] = scaffold_sequences
    return obj
