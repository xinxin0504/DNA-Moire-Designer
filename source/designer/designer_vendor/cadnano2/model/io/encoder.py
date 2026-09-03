from json import dumps
from .legacyencoder import legacy_dict_from_doc, legacy_dict_from_part


HYBRID_FORMAT = 'cadnano2-hybrid-v1'


def encode(document, helixOrderList, io, includeSequences=False):
    if document.isHybrid():
        parts = document.partsByLattice()
        obj = {'.format': HYBRID_FORMAT,
               'name': getattr(io, 'name', 'hybrid.json'),
               'parts': {},
               'hybrid_connections': document.hybridConnectionRecords()}
        for lattice in ('honeycomb', 'square'):
            part = parts.get(lattice)
            if part is not None:
                obj['parts'][lattice] = legacy_dict_from_part(
                    part, io.name, includeSequences=includeSequences)
        if includeSequences:
            hybridSequences = []
            for oligo in document.oligos():
                if not oligo.isHybrid() or oligo.isStaple() or \
                        oligo.isLoop():
                    continue
                sequence = oligo.sequence()
                if sequence and sequence.strip():
                    hybridSequences.append({
                        'start': oligo.sequenceEndpoints()[0],
                        'sequence': sequence})
            if hybridSequences:
                obj['hybrid_sequences'] = hybridSequences
    else:
        obj = legacy_dict_from_doc(document, io.name, helixOrderList,
                                   includeSequences=includeSequences)
    athenaMetadata = document.athenaMetadata()
    if athenaMetadata:
        # ATHENA geometry is design metadata, not applied sequence data, so
        # both ordinary Save and Save As with Sequences retain it.
        obj['athena_metadata'] = athenaMetadata
    curvedMetadata = document.curvedMetadata()
    if curvedMetadata:
        # DNAxiS target-space geometry is independent of whether sequence
        # information is included, so every Save mode retains it.
        obj['curved_metadata'] = curvedMetadata
    twistBendMetadata = document.twistBendMetadata()
    if twistBendMetadata:
        # Parameters are retained by every save mode so the interactive
        # preview and task list can be restored when the design is reopened.
        obj['twist_bend_metadata'] = twistBendMetadata
    json_string = dumps(obj, separators=(',',':'))  # compact encoding
    io.write(json_string)
