from os.path import basename
from cadnano2.model.enum import StrandType


def legacy_dict_from_doc(document, fname, helixOrderList,
                         includeSequences=False):
    part = document.selectedPart()
    return legacy_dict_from_part(part, fname, helixOrderList,
                                 includeSequences=includeSequences)


def legacy_dict_from_part(part, fname, helixOrderList=None,
                          includeSequences=False):
    """Encode one lattice part; cross-part connections remain endpoints."""
    numBases = part.maxBaseIdx()+1

    if helixOrderList is None:
        helixOrderList = part._importedVHelixOrder
    if helixOrderList is None:
        helixOrderList = [vh.coord() for vh in sorted(
            part.getVirtualHelices(),
            key=lambda virtualHelix: virtualHelix.number())]

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
        stapStrandSet = vh.stapleStrandSet()
        for strand in stapStrandSet:
            if strand.connection5p() == None:
                c = str(strand.oligo().color())[1:]  # drop the hash
                stapColors.append([strand.idx5Prime(), int(c, 16)])

        vhDict = {"row": row,
                  "col": col,
                  "num": vh.number(),
                  "scaf": vh.getLegacyStrandSetArray(StrandType.Scaffold),
                  "stap": vh.getLegacyStrandSetArray(StrandType.Staple),
                  "loop": insts,
                  "skip": skips,
                  "scafLoop": [],
                  "stapLoop": [],
                  "stap_colors": stapColors}
        vhList.append(vhDict)
    bname = basename(str(fname))
    obj = {"name": bname, "vstrands": vhList,
           "num_bases": numBases,
           "lattice": "honeycomb" if part._step == 21 else "square"}
    scaffoldColors = part.getScaffoldColorRecords()
    if scaffoldColors:
        # Extension field: scaffold colors are always part of the design,
        # independently of the optional applied-sequence data below.
        obj["scaffold_colors"] = scaffoldColors
    scaffoldSequences = part.getScaffoldSequenceRecords()
    if includeSequences and scaffoldSequences:
        # Extension field: legacy caDNAno readers ignore unknown top-level
        # keys, while this version uses it to restore applied sequences.
        obj["scaffold_sequences"] = scaffoldSequences
    guidedRegions = part.guidedDuplexRegionRecords()
    if guidedRegions is not None:
        obj["guided_duplex_regions"] = guidedRegions
    scaffoldOnlyRegions = part.scaffoldOnlyRegionRecords()
    if scaffoldOnlyRegions:
        obj["scaffold_only_regions"] = scaffoldOnlyRegions
    if part._autobreakStaplesApplied:
        obj["autobreak_staples_applied"] = True
    return obj
