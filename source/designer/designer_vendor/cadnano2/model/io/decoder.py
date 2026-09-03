import json
from .legacydecoder import import_legacy_dict
from .encoder import HYBRID_FORMAT
from cadnano2.model.enum import LatticeType, StrandType
import cadnano2.util as util
import cadnano2.cadnano as cadnano
if cadnano.app().isGui():  # headless:
    from cadnano2.ui.dialogs.ui_latticetype import Ui_LatticeType
    util.qtWrapImport('QtWidgets', globals(),  ['QDialog', 'QDialogButtonBox'])


def decode(document, string):
    if cadnano.app().isGui():
        # from ui.dialogs.ui_latticetype import Ui_LatticeType
        # util.qtWrapImport('QtGui', globals(),  ['QDialog', 'QDialogButtonBox'])
        dialog = QDialog()
        dialogLT = Ui_LatticeType()  # reusing this dialog, should rename
        dialogLT.setupUi(dialog)

    # try:  # try to do it fast
    #     try:
    #         import cjson
    #         packageObject = cjson.decode(string)
    #     except:  # fall back to if cjson not available or on decode error
    #         packageObject = json.loads(string)
    # except ValueError:
    #     dialogLT.label.setText("Error decoding JSON object.")
    #     dialogLT.buttonBox.setStandardButtons(QDialogButtonBox.Ok)
    #     dialog.exec()
    #     return
    packageObject = json.loads(string)
    document.setAthenaMetadata(packageObject.get('athena_metadata'))
    document.setCurvedMetadata(packageObject.get('curved_metadata'))
    document.setTwistBendMetadata(packageObject.get('twist_bend_metadata'))

    if packageObject.get('.format', None) == HYBRID_FORMAT:
        document.setHybrid(True)
        parts = {}
        encodedParts = packageObject.get('parts', {})
        for lattice, latticeType in (
                ('honeycomb', LatticeType.Honeycomb),
                ('square', LatticeType.Square)):
            encoded = encodedParts.get(lattice)
            if encoded is not None:
                parts[lattice] = import_legacy_dict(
                    document, encoded, latticeType,
                    forceLatticeType=True)
        for record in packageObject.get('hybrid_connections', []):
            try:
                three = record['three_prime']
                five = record['five_prime']
                strandType = (StrandType.Scaffold if
                              record['strand_type'] == 'scaffold' else
                              StrandType.Staple)
                part3 = parts[three['lattice']]
                part5 = parts[five['lattice']]
                vh3 = part3.virtualHelix(int(three['helix']))
                vh5 = part5.virtualHelix(int(five['helix']))
                strand3 = vh3.getStrandSetByType(strandType).getStrand(
                    int(three['index']))
                strand5 = vh5.getStrandSetByType(strandType).getStrand(
                    int(five['index']))
            except (KeyError, TypeError, ValueError, AttributeError):
                continue
            document.createHybridConnection(strand3, '3p', strand5, '5p',
                                            useUndoStack=False)
        targets = dict((oligo.sequenceEndpoints()[0], oligo)
                       for oligo in document.oligos()
                       if not oligo.isLoop() and not oligo.isStaple())
        for record in packageObject.get('hybrid_sequences', []):
            oligo = targets.get(record.get('start'))
            sequence = record.get('sequence')
            if oligo is not None and isinstance(sequence, str):
                oligo.applySequence(sequence, useUndoStack=False)
        if parts:
            document.setSelectedPart(parts.get('honeycomb') or
                                     next(iter(parts.values())))
        return
    if packageObject.get('.format', None) != 'caDNAno2':
        import_legacy_dict(document, packageObject)
