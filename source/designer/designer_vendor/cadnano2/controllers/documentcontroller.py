import io
import os
import traceback
from copy import deepcopy
from math import ceil
from ..data.dnasequences import sequences as standardScaffoldSequences
from ..cadnano import app
from ..model.document import Document
from ..model.enum import StrandType
from ..model.io.decoder import decode
from ..model.io.encoder import encode
from ..model.io.sequenceimport import validate_sequence_import
from ..model.io.sequencexlsx import (read_sequence_template,
                                     write_sequence_template,
                                     write_sequence_workbook)
from ..model.io.orthogonalseq import (
    DEFAULT_SETTINGS as ORTHOGONAL_SEQUENCE_DEFAULTS,
    GenerationCancelled, generate_sequences, normalized_settings,
    read_sequence_text, write_orthogonal_workbook)
from ..model.io.oxdnaexport import (export_structure_bundle,
                                    structure_bundle_paths,
                                    athena_structure_bundle_paths,
                                    export_athena_structure_bundle,
                                    curved_structure_bundle_paths,
                                    export_curved_structure_bundle)
from ..model.io.athena import create_project, wireframe_output_name
from ..model.io.curved import (build_rings, create_curved_project,
                               curved_indel_plan, curved_output_name)
from ..model.io.frame import create_frame_project
from ..model.io.curvedreport import start_curved_report_export
from ..model.io.indelanalysis import (
    single_helix_distribution_data, write_pair_curvature_csv,
    write_pair_curvature_svg, write_single_helix_distribution_csv)
from ..views.documentwindow import DocumentWindow
from ..views.curvedreport import CurvedReportDialog
from ..views.twistbendreport import (TwistBendReportDialog,
                                     create_twistbend_report_image)
from ..views.orthogonalsequences import OrthogonalSequenceDialog
from ..views.primer3analysis import Primer3AnalysisDialog
from ..views import styles
from ..views.pathview.hybridoverlay import (drawHybridConnection,
                                             hybridConnections,
                                             hybridSceneEndpoint)
import cadnano2.util as util
util.qtWrapImport('QtCore', globals(), ['QDir', 'QFileInfo', 'QPointF', 'QRect',
                                        'QRectF',
                                        'QSettings',
                                        'QSize', 'Qt', 'QTimer'])
util.qtWrapImport('QtGui', globals(), [
                                       'QFont',
                                       'QIcon',
                                       'QImage',
                                       'QFontMetricsF',
                                       'QKeySequence',
                                       'QPainter',
                                       'QUndoCommand',
                                       'QTransform'
                                       ])
util.qtWrapImport('QtWidgets', globals(), ['QApplication',
                                           'QCheckBox',
                                           'QComboBox',
                                           'QDialog',
                                           'QDialogButtonBox',
                                           'QDockWidget',
                                           'QFileDialog',
                                           'QFormLayout',
                                           'QGraphicsItem',
                                           'QGraphicsSimpleTextItem',
                                           'QInputDialog',
                                           'QLabel',
                                           'QMainWindow',
                                           'QMessageBox',
                                           'QProgressDialog',
                                           'QSpinBox',
                                           'QStyleOptionGraphicsItem'])
util.qtWrapImport('QtSvg', globals(), ['QSvgGenerator'])


class _SetTwistBendMetadataCommand(QUndoCommand):
    """Keep saved editor parameters in the same undo step as its indels."""

    def __init__(self, document, metadata):
        super(_SetTwistBendMetadataCommand, self).__init__()
        self._document = document
        self._old = deepcopy(document.twistBendMetadata())
        self._new = deepcopy(metadata)

    def redo(self):
        self._document.setTwistBendMetadata(deepcopy(self._new))

    def undo(self):
        self._document.setTwistBendMetadata(deepcopy(self._old))


class DocumentController():
    """
    Connects UI buttons to their corresponding actions in the model.
    """
    ### INIT METHODS ###
    def __init__(self):
        """docstring for __init__"""
        # initialize variables
        self._document = Document()
        self._document.setController(self)
        self._activePart = None
        self._filename = None
        self._fileOpenPath = None  # will be set in _readSettings
        self._hasNoAssociatedFile = True
        self._pathViewInstance = None
        self._sliceViewInstance = None
        self._undoStack = None
        self.win = None
        self.fileopendialog = None
        self.filesavedialog = None
        self.saveCsvDialog = None
        self._saveSequencesOnDialog = False
        self._continuousStapleLabel = None
        self._continuousStapleStatsPending = False
        self._lastPrimer3Entries = []
        self._curvedReportWatchers = []
        self._curvedReportDialogs = []

        self.settings = QSettings()
        self._readSettings()

        QDir.addSearchPath('icons', 'ui/mainwindow/images/')

        # call other init methods
        self._initWindow()
        self._updateAthenaActions()
        self._continuousStapleLabel = QLabel(self.win)
        self._continuousStapleLabel.setMinimumWidth(720)
        self._continuousStapleLabel.setWordWrap(False)
        statisticsFont = self._continuousStapleLabel.font()
        statisticsFont.setPointSize(9)
        self._continuousStapleLabel.setFont(statisticsFont)
        self._continuousStapleLabel.setToolTip(
            "Percentage of staple strands containing at least one "
            "14/16-base continuous actual-nucleotide region (including "
            "insertions and deletions), plus the real-time "
            "staple-length distribution in 10-base intervals.")
        self.win.statusBar().addPermanentWidget(
                                        self._continuousStapleLabel)
        # The undo-stack signal covers all normal edits, including crossover
        # creation/removal and Undo/Redo.  The part signal also covers the few
        # bulk/import operations that deliberately bypass the undo stack.
        self._document.undoStack().indexChanged.connect(
                                self._scheduleContinuousStapleStatistics)
        self._updateContinuousStapleStatistics()
        if app().isInMaya():
            self._initMaya()
        app().documentControllers.add(self)

    def _initWindow(self):
        """docstring for initWindow"""
        self.win = DocumentWindow(docCtrlr=self)
        self.win.setWindowIcon(QIcon('icons:cadnano2-app-icon.png'))
        app().documentWindowWasCreatedSignal.emit(self._document, self.win)
        self._connectWindowSignalsToSelf()
        self.win.show()

    def _initMaya(self):
        """
        Initialize Maya-related state. Delete Maya nodes if there
        is an old document left over from the same session. Set up
        the Maya window.
        """
        # There will only be one document
        if (app().activeDocument and app().activeDocument.win and
                                not app().activeDocument.win.close()):
            return
        del app().activeDocument
        app().activeDocument = self

        import maya.OpenMayaUI as OpenMayaUI
        import sip
        ptr = OpenMayaUI.MQtUtil.mainWindow()
        mayaWin = sip.wrapinstance(int(ptr), QMainWindow)
        self.windock = QDockWidget("cadnano")
        self.windock.setFeatures(QDockWidget.DockWidgetMovable
                                 | QDockWidget.DockWidgetFloatable)
        self.windock.setAllowedAreas(Qt.LeftDockWidgetArea
                                     | Qt.RightDockWidgetArea)
        self.windock.setWidget(self.win)
        mayaWin.addDockWidget(Qt.DockWidgetArea(Qt.LeftDockWidgetArea),
                                self.windock)
        self.windock.setVisible(True)

    def destroyDC(self):
        self.disconnectSignalsToSelf()
        if self.win is not None:
            self.win.destroyWin()
            self.win = None
    # end def

    def disconnectSignalsToSelf(self):
        win = self.win
        if win is not None:
            win.actionNew.triggered.disconnect(self.actionNewSlot)
            win.actionOpen.triggered.disconnect(self.actionOpenSlot)
            win.actionDrop.triggered.disconnect(self.actionDropSlot)
            win.actionClose.triggered.disconnect(self.actionCloseSlot)
            win.actionSave.triggered.disconnect(self.actionSaveSlot)
            win.actionSaveWithSequences.triggered.disconnect(
                                        self.actionSaveWithSequencesSlot)
            win.actionExportSequenceTemplate.triggered.disconnect(
                                        self.actionExportSequenceTemplateSlot)
            win.actionImportSequenceXlsx.triggered.disconnect(
                                        self.actionImportSequenceXlsxSlot)
            win.actionOrthogonalSequences.triggered.disconnect(
                                        self.actionOrthogonalSequencesSlot)
            win.actionX3D.toggled.disconnect(self.actionToggleThreeDSlot)
            win.actionSave_As.triggered.disconnect(self.actionSaveAsSlot)
            win.actionSVG.triggered.disconnect(self.actionSVGSlot)
            win.actionIllustrator.triggered.disconnect(
                                        self.actionIllustratorSlot)
            win.actionAutoScaffoldWithoutCS.triggered.disconnect(
                                self.actionAutoScaffoldWithoutCSSlot)
            win.actionAutoScaffoldCrossovers.triggered.disconnect(
                                self.actionAutoScaffoldCrossoversSlot)
            win.actionAutoStapleWithoutCS.triggered.disconnect(
                                self.actionAutoStapleWithoutCSSlot)
            win.actionAutoStaple.triggered.disconnect(self.actionAutostapleSlot)
            win.actionAutoBreakStaples.triggered.disconnect(
                                        self.actionAutoBreakStaplesSlot)
            win.actionExportStaples.triggered.disconnect(self.actionExportStaplesSlot)
            win.actionExportCsv.triggered.disconnect(self.actionExportCsvSlot)
            win.actionExportPdb.triggered.disconnect(self.actionExportPdbSlot)
            win.actionAthenaDesign.triggered.disconnect(
                                        self.actionAthenaDesignSlot)
            win.actionAthenaExport.triggered.disconnect(
                                        self.actionAthenaExportSlot)
            win.actionCurvedDesign.triggered.disconnect(
                                        self.actionCurvedDesignSlot)
            win.actionFrameDesign.triggered.disconnect(
                                        self.actionFrameDesignSlot)
            win.actionTwistBend.triggered.disconnect(
                                        self.actionTwistBendSlot)
            win.actionCurvedExport.triggered.disconnect(
                                        self.actionCurvedExportSlot)
            win.actionPreferences.triggered.disconnect(self.actionPrefsSlot)
            win.actionModify.triggered.disconnect(self.actionModifySlot)
            win.actionNewHoneycombPart.triggered.disconnect(self.actionAddHoneycombPartSlot)
            win.actionNewSquarePart.triggered.disconnect(self.actionAddSquarePartSlot)
            win.actionNewHybridPart.triggered.disconnect(
                                        self.actionAddHybridPartSlot)
            win.closeEvent = self.windowCloseEventHandler
            win.actionAbout.triggered.disconnect(self.actionAboutSlot)
            win.actionCadnanoWebsite.triggered.disconnect(self.actionCadnanoWebsiteSlot)
            win.actionFeedback.triggered.disconnect(self.actionFeedbackSlot)
            win.actionFilterHandle.triggered.disconnect(self.actionFilterHandleSlot)
            win.actionFilterEndpoint.triggered.disconnect(self.actionFilterEndpointSlot)
            win.actionFilterStrand.triggered.disconnect(self.actionFilterStrandSlot)
            win.actionFilterXover.triggered.disconnect(self.actionFilterXoverSlot)
            win.actionFilterScaf.triggered.disconnect(self.actionFilterScafSlot)
            win.actionFilterStap.triggered.disconnect(self.actionFilterStapSlot)
            win.actionRenumber.triggered.disconnect(self.actionRenumberSlot)
    # end def

    def _connectWindowSignalsToSelf(self):
        """This method serves to group all the signal & slot connections
        made by DocumentController"""
        self.win.actionNew.triggered.connect(self.actionNewSlot)
        self.win.actionOpen.triggered.connect(self.actionOpenSlot)
        self.win.actionDrop.triggered.connect(self.actionDropSlot)
        self.win.actionClose.triggered.connect(self.actionCloseSlot)
        self.win.actionSave.triggered.connect(self.actionSaveSlot)
        self.win.actionSaveWithSequences.triggered.connect(
                                        self.actionSaveWithSequencesSlot)
        self.win.actionExportSequenceTemplate.triggered.connect(
                                        self.actionExportSequenceTemplateSlot)
        self.win.actionImportSequenceXlsx.triggered.connect(
                                        self.actionImportSequenceXlsxSlot)
        self.win.actionOrthogonalSequences.triggered.connect(
                                        self.actionOrthogonalSequencesSlot)
        self.win.actionX3D.toggled.connect(self.actionToggleThreeDSlot)
        self.win.actionSave_As.triggered.connect(self.actionSaveAsSlot)
        self.win.actionSVG.triggered.connect(self.actionSVGSlot)
        self.win.actionIllustrator.triggered.connect(
                                        self.actionIllustratorSlot)
        self.win.actionAutoScaffoldWithoutCS.triggered.connect(
                                self.actionAutoScaffoldWithoutCSSlot)
        self.win.actionAutoScaffoldCrossovers.triggered.connect(
                                self.actionAutoScaffoldCrossoversSlot)
        self.win.actionAutoStapleWithoutCS.triggered.connect(
                                self.actionAutoStapleWithoutCSSlot)
        self.win.actionAutoStaple.triggered.connect(self.actionAutostapleSlot)
        self.win.actionAutoBreakStaples.triggered.connect(
                                        self.actionAutoBreakStaplesSlot)
        self.win.actionExportStaples.triggered.connect(self.actionExportStaplesSlot)
        self.win.actionExportCsv.triggered.connect(self.actionExportCsvSlot)
        self.win.actionExportPdb.triggered.connect(self.actionExportPdbSlot)
        self.win.actionAthenaDesign.triggered.connect(
                                        self.actionAthenaDesignSlot)
        self.win.actionAthenaExport.triggered.connect(
                                        self.actionAthenaExportSlot)
        self.win.actionCurvedDesign.triggered.connect(
                                        self.actionCurvedDesignSlot)
        self.win.actionFrameDesign.triggered.connect(
                                        self.actionFrameDesignSlot)
        self.win.actionTwistBend.triggered.connect(
                                        self.actionTwistBendSlot)
        self.win.actionCurvedExport.triggered.connect(
                                        self.actionCurvedExportSlot)
        self.win.actionPreferences.triggered.connect(self.actionPrefsSlot)
        self.win.actionModify.triggered.connect(self.actionModifySlot)
        self.win.actionNewHoneycombPart.triggered.connect(self.actionAddHoneycombPartSlot)
        self.win.actionNewSquarePart.triggered.connect(self.actionAddSquarePartSlot)
        self.win.actionNewHybridPart.triggered.connect(
                                        self.actionAddHybridPartSlot)
        self.win.closeEvent = self.windowCloseEventHandler
        self.win.actionAbout.triggered.connect(self.actionAboutSlot)
        self.win.actionCadnanoWebsite.triggered.connect(self.actionCadnanoWebsiteSlot)
        self.win.actionFeedback.triggered.connect(self.actionFeedbackSlot)
        self.win.actionFilterHandle.triggered.connect(self.actionFilterHandleSlot)
        self.win.actionFilterEndpoint.triggered.connect(self.actionFilterEndpointSlot)
        self.win.actionFilterStrand.triggered.connect(self.actionFilterStrandSlot)
        self.win.actionFilterXover.triggered.connect(self.actionFilterXoverSlot)
        self.win.actionFilterScaf.triggered.connect(self.actionFilterScafSlot)
        self.win.actionFilterStap.triggered.connect(self.actionFilterStapSlot)
        self.win.actionRenumber.triggered.connect(self.actionRenumberSlot)


    ### SLOTS ###
    def undoStackCleanChangedSlot(self):
        """The title changes to include [*] on modification."""
        self.win.setWindowModified(not self.undoStack().isClean())
        self.win.setWindowTitle(self.documentTitle())

    def actionAboutSlot(self):
        """Displays the about cadnano dialog."""
        from cadnano2.ui.dialogs.ui_about import Ui_About
        dialog = QDialog()
        dialogAbout = Ui_About()  # reusing this dialog, should rename
        dialog.setStyleSheet("QDialog { background-image: url(ui/dialogs/images/cadnano2-about.png); background-repeat: none; }")
        dialogAbout.setupUi(dialog)
        dialog.exec()

    filterList = ["strand", "endpoint", "xover", "virtualHelix"]
    def actionFilterHandleSlot(self):
        """Disables all other selection filters when active."""
        fH = self.win.actionFilterHandle
        fE = self.win.actionFilterEndpoint
        fS = self.win.actionFilterStrand
        fX = self.win.actionFilterXover
        fH.setChecked(True)
        if fE.isChecked():
            fE.setChecked(False)
        if fS.isChecked():
            fS.setChecked(False)
        if fX.isChecked():
            fX.setChecked(False)
        self._document.documentSelectionFilterChangedSignal.emit(["virtualHelix"])

    def actionFilterEndpointSlot(self):
        """
        Disables handle filters when activated.
        Remains checked if no other item-type filter is active.
        """
        fH = self.win.actionFilterHandle
        fE = self.win.actionFilterEndpoint
        fS = self.win.actionFilterStrand
        fX = self.win.actionFilterXover
        if fH.isChecked():
            fH.setChecked(False)
        if not fS.isChecked() and not fX.isChecked():
            fE.setChecked(True)
        self._strandFilterUpdate()
    # end def

    def actionFilterStrandSlot(self):
        """
        Disables handle filters when activated.
        Remains checked if no other item-type filter is active.
        """
        fH = self.win.actionFilterHandle
        fE = self.win.actionFilterEndpoint
        fS = self.win.actionFilterStrand
        fX = self.win.actionFilterXover
        if fH.isChecked():
            fH.setChecked(False)
        if not fE.isChecked() and not fX.isChecked():
            fS.setChecked(True)
        self._strandFilterUpdate()
    # end def

    def actionFilterXoverSlot(self):
        """
        Disables handle filters when activated.
        Remains checked if no other item-type filter is active.
        """
        fH = self.win.actionFilterHandle
        fE = self.win.actionFilterEndpoint
        fS = self.win.actionFilterStrand
        fX = self.win.actionFilterXover
        if fH.isChecked():
            fH.setChecked(False)
        if not fE.isChecked() and not fS.isChecked():
            fX.setChecked(True)
        self._strandFilterUpdate()
    # end def

    def actionFilterScafSlot(self):
        """Remains checked if no other strand-type filter is active."""
        fSc = self.win.actionFilterScaf
        fSt = self.win.actionFilterStap
        if not fSc.isChecked() and not fSt.isChecked():
            fSc.setChecked(True)
        self._strandFilterUpdate()

    def actionFilterStapSlot(self):
        """Remains checked if no other strand-type filter is active."""
        fSc = self.win.actionFilterScaf
        fSt = self.win.actionFilterStap
        if not fSc.isChecked() and not fSt.isChecked():
            fSt.setChecked(True)
        self._strandFilterUpdate()
    # end def

    def _strandFilterUpdate(self):
        win = self.win
        filterList = []
        if win.actionFilterEndpoint.isChecked():
            filterList.append("endpoint")
        if win.actionFilterStrand.isChecked():
            filterList.append("strand")
        if win.actionFilterXover.isChecked():
            filterList.append("xover")
        if win.actionFilterScaf.isChecked():
            filterList.append("scaffold")
        if win.actionFilterStap.isChecked():
            filterList.append("staple")
        self._document.documentSelectionFilterChangedSignal.emit(filterList)
    # end def

    def actionNewSlot(self):
        """
        1. If document is has no parts, do nothing.
        2. If document is dirty, call maybeSave and continue if it succeeds.
        3. Create a new document and swap it into the existing ctrlr/window.
        """
        # clear/reset the view!

        if len(self._document.parts()) == 0:
            return  # no parts
        if self.maybeSave() == False:
            return  # user canceled in maybe save
        else:  # user did not cancel
            if self.filesavedialog != None:
                self.filesavedialog.finished.connect(self.newClickedCallback)
            else:  # user did not save
                self.newClickedCallback()  # finalize new

    def actionOpenSlot(self):
        """
        1. If document is untouched, proceed to open dialog.
        2. If document is dirty, call maybesave and continue if it succeeds.
        Downstream, the file is selected in openAfterMaybeSave,
        and the selected file is actually opened in openAfterMaybeSaveCallback.
        """
        if self.maybeSave() == False:
            return  # user canceled in maybe save
        else:  # user did not cancel
            if hasattr(self, "filesavedialog"): # user did save
                if self.filesavedialog != None:
                    self.filesavedialog.finished.connect(self.openAfterMaybeSave)
                else:
                    self.openAfterMaybeSave()  # windows
            else:  # user did not save
                self.openAfterMaybeSave()  # finalize new

    def actionDropSlot(self):
        """
        1. If document is untouched, proceed to open dialog.
        2. If document is dirty, call maybesave and continue if it succeeds.
        Equivalent to actionOpenSlot() for Drag and Drop Event.
        """
        if self.maybeSave() is False:
            return  # user canceled in maybe save
        else:  # user did not cancel
            if hasattr(self, "filesavedialog"):  # user did save
                if self.filesavedialog is not None:
                    self.filesavedialog.finished.connect(self.openDropAfterMaybeSave)
                else:
                    self.openDropAfterMaybeSave()  # windows
            else:  # user did not save
                self.openDropAfterMaybeSave()  # finalize new

    def actionCloseSlot(self):
        """This will trigger a Window closeEvent."""
        if util.isWindows():
            #print "close win"
            if self.win is not None:
                self.win.close()
            if not app().isInMaya():
                #print "exit app"
                import sys
                sys.exit(1)

    def actionSaveSlot(self):
        """Save design structure without sequence information."""
        if self._hasNoAssociatedFile:
            self.saveFileDialog(includeSequences=False)
            return
        self.writeDocumentToFile(includeSequences=False)

    def actionSaveWithSequencesSlot(self):
        """Save a new design file with sequences and 5' positions."""
        self.saveFileDialog(includeSequences=True)

    def actionSaveAsSlot(self):
        """Save design structure to a chosen file, without sequences."""
        self.saveFileDialog(includeSequences=False)

    def actionToggleThreeDSlot(self, isVisible):
        """Create/destroy the optional 3D dock on demand."""
        self.win.setThreeDVisible(isVisible)

    def actionExportSequenceTemplateSlot(self):
        """Export every importable scaffold chain to an XLSX template."""
        part = self.activePart()
        if part is None:
            self._showSequenceMessage(
                QMessageBox.Icon.Warning, "Sequence Template",
                "The current document does not contain a design.")
            return
        rows = (self._document.getScaffoldSequenceTemplateRows()
                if self._document.isHybrid() else
                part.getScaffoldSequenceTemplateRows())
        if not rows:
            self._showSequenceMessage(
                QMessageBox.Icon.Warning, "Sequence Template",
                "The design does not contain any non-circular scaffold "
                "chains that can accept a sequence.")
            return
        directory = ("." if self.filename() is None else
                     QFileInfo(self.filename()).path())
        fname = QFileDialog.getSaveFileName(
                    self.win,
                    "%s - Export Sequence Template" %
                    QApplication.applicationName(),
                    directory,
                    "Sequence Import Template (*.xlsx)")
        if isinstance(fname, (list, tuple)):
            fname = fname[0]
        if not fname or os.path.isdir(fname):
            return
        if not fname.lower().endswith(".xlsx"):
            fname += ".xlsx"
        try:
            write_sequence_template(fname, rows)
        except (IOError, OSError) as error:
            self._showSequenceMessage(
                QMessageBox.Icon.Critical, "Sequence Template",
                "Could not write the XLSX template.", str(error))
            return
        self.win.statusBar().showMessage(
            "Exported %d scaffold chains to %s" % (len(rows), fname),
            8000)

    def actionImportSequenceXlsxSlot(self):
        """Strictly validate and atomically apply sequences from XLSX."""
        part = self.activePart()
        if part is None:
            self._showSequenceMessage(
                QMessageBox.Icon.Warning, "Sequence Import",
                "The current document does not contain a design.")
            return
        directory = ("." if self.filename() is None else
                     QFileInfo(self.filename()).path())
        fname = QFileDialog.getOpenFileName(
                    self.win,
                    "%s - Import Sequences" %
                    QApplication.applicationName(),
                    directory,
                    "Sequence Import Workbook (*.xlsx)")
        if isinstance(fname, (list, tuple)):
            fname = fname[0]
        if not fname or os.path.isdir(fname):
            return
        try:
            headers, workbookRows = read_sequence_template(fname)
            operations, errors = validate_sequence_import(
                self._document if self._document.isHybrid() else part,
                headers, workbookRows)
        except (IOError, OSError, ValueError) as error:
            self._showSequenceMessage(
                QMessageBox.Icon.Critical, "Sequence Import Failed",
                "No sequences were imported because the XLSX file could "
                "not be read.", str(error))
            return
        if errors:
            self._showSequenceMessage(
                QMessageBox.Icon.Critical, "Sequence Import Failed",
                "No sequences were imported. %d validation error(s) were "
                "found; expand Details to see the chain, field, and exact "
                "location of each error." % len(errors),
                '\n'.join(errors))
            return
        if not operations:
            self._showSequenceMessage(
                QMessageBox.Icon.Information, "Sequence Import",
                "No sequences were imported because every Sequence cell "
                "was blank.")
            return

        undoStack = self.undoStack()
        undoStack.beginMacro("Import Scaffold Sequences")
        try:
            for oligo, sequence, start, rowNumber in operations:
                oligo.applySequence(sequence)
        finally:
            undoStack.endMacro()
        self._showSequenceMessage(
            QMessageBox.Icon.Information, "Sequence Import Complete",
            "Successfully imported %d scaffold sequence(s). The entire "
            "import can be reverted with one Undo." % len(operations))

    def actionOrthogonalSequencesSlot(self):
        """Design and export a biochemically screened orthogonal DNA set."""
        directory = (os.path.expanduser('~/Desktop') if self.filename() is None
                     else QFileInfo(self.filename()).path())
        dialog = OrthogonalSequenceDialog(
            ORTHOGONAL_SEQUENCE_DEFAULTS, self.win,
            primer3_entries=self._lastPrimer3Entries,
            suggested_directory=directory)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        requested = dialog.settings()
        inputFilename = requested.pop('input_file', '')
        scaffoldName = requested.pop('scaffold_name', '')
        try:
            settings = normalized_settings(requested)
        except ValueError as error:
            self._showSequenceMessage(
                QMessageBox.Icon.Warning, '正交序列设计', str(error))
            return

        inputSequences = []
        inputErrors = []
        if inputFilename:
            try:
                inputSequences, inputErrors = read_sequence_text(inputFilename)
            except (IOError, OSError, UnicodeError) as error:
                self._showSequenceMessage(
                    QMessageBox.Icon.Critical, '正交序列输入失败',
                    '无法读取所选 TXT 文件。', str(error))
                return
        if inputErrors:
            self._showSequenceMessage(
                QMessageBox.Icon.Critical, '正交序列输入失败',
                'TXT 中发现 %d 个格式错误；没有开始生成。' %
                len(inputErrors), '\n'.join(inputErrors))
            return

        scaffoldBackground = []
        if scaffoldName:
            scaffoldSequence = standardScaffoldSequences.get(scaffoldName)
            if not scaffoldSequence:
                self._showSequenceMessage(
                    QMessageBox.Icon.Critical, '正交序列设计',
                    '无法找到所选骨架链：%s' % scaffoldName)
                return
            scaffoldBackground.append((scaffoldName, scaffoldSequence))

        designName = ('cadnano' if self.filename() is None else
                      os.path.splitext(os.path.basename(
                          str(self.filename())))[0])
        suggested = os.path.join(
            directory, designName + '_orthogonal_sequences.xlsx')
        filename = QFileDialog.getSaveFileName(
            self.win, '正交序列设计 — 保存结果', suggested,
            'Orthogonal Sequence Workbook (*.xlsx)')
        if isinstance(filename, (tuple, list)):
            filename = filename[0]
        if not filename or os.path.isdir(filename):
            return
        if not filename.lower().endswith('.xlsx'):
            filename += '.xlsx'

        progressDialog = QProgressDialog(
            '正在生成并评价候选序列…', '取消', 0, 0, self.win)
        progressDialog.setWindowTitle('正交序列设计')
        progressDialog.setWindowModality(Qt.WindowModality.WindowModal)
        progressDialog.setMinimumDuration(0)

        def updateProgress(accepted, attempts):
            progressDialog.setLabelText(
                '已生成 %d / %d 条；已评价 %d 个候选…' %
                (accepted, settings['count'], attempts))
            QApplication.processEvents()

        try:
            result = generate_sequences(
                settings, background_sequences=inputSequences,
                scaffold_sequences=scaffoldBackground,
                progress=updateProgress,
                cancelled=progressDialog.wasCanceled)
            result['input_file'] = inputFilename
        except GenerationCancelled:
            progressDialog.close()
            self.win.statusBar().showMessage('正交序列生成已取消。', 5000)
            return
        except Exception as error:
            progressDialog.close()
            traceback.print_exc()
            self._showSequenceMessage(
                QMessageBox.Icon.Critical, '正交序列设计失败',
                '生成过程中发生错误。', str(error))
            return
        finally:
            progressDialog.close()

        if not result['sequences']:
            self._showSequenceMessage(
                QMessageBox.Icon.Warning, '正交序列设计',
                '在当前约束和最大尝试次数内没有找到合格序列。请放宽参数，'
                '尤其是 GC、同向片段或链间互补长度；若启用了高级规则，也可'
                '逐项关闭或放宽。',
                self._orthogonalRejectionDetails(result))
            return
        if not result['complete']:
            choice = QMessageBox.question(
                self.win, '正交序列设计未完全完成',
                '请求 %d 条，但在有限尝试次数内只找到 %d 条。是否保存这些'
                '已经通过全部检查的序列？' %
                (settings['count'], len(result['sequences'])),
                QMessageBox.StandardButton.Yes |
                QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.Yes)
            if choice != QMessageBox.StandardButton.Yes:
                return
        try:
            write_orthogonal_workbook(filename, result)
        except (IOError, OSError, ValueError) as error:
            self._showSequenceMessage(
                QMessageBox.Icon.Critical, '正交序列导出失败',
                '无法写入 XLSX 文件。', str(error))
            return
        self._lastPrimer3Entries = [
            ('输入', '输入-%03d' % index, sequence)
            for index, sequence in enumerate(
                result.get('input_sequences', ()), 1)]
        self._lastPrimer3Entries.extend(
            ('新生成', '新序列-%03d' % index, sequence)
            for index, sequence in enumerate(result.get('sequences', ()), 1))
        completion = QMessageBox(
            QMessageBox.Icon.Information, '正交序列设计完成',
            '已读取 %d 条输入序列、采用 %d 条骨架链、生成 %d 条新序列，'
            '并保存到：\n%s\n\n'
            '“序列分析”会区分输入、骨架链和新生成序列；骨架链仅参与筛选，'
            '不会列入“两两分析”。' %
            (len(result.get('input_sequences', ())),
             len(result.get('scaffold_sequences', ())),
             len(result['sequences']), filename),
            QMessageBox.StandardButton.Ok, self.win)
        completion.setDetailedText(self._orthogonalRejectionDetails(result))
        primer3Button = completion.addButton(
            'Primer3 热力学分析…', QMessageBox.ButtonRole.ActionRole)
        completion.exec()
        if completion.clickedButton() is primer3Button:
            Primer3AnalysisDialog(
                self._lastPrimer3Entries, directory, self.win).exec()

    @staticmethod
    def _orthogonalRejectionDetails(result):
        labels = {
            'global_gc': '全局 GC', 'local_gc': '局部 GC',
            'homopolymer': '均聚碱基',
            'entropy': '低复杂度', 'self_complement': '自身互补',
            'hairpin': '发卡', 'forbidden_motif': '禁用 motif',
            'same_substring': '同向相同片段',
            'cross_complement': '链间互补',
            'hamming': '汉明距离',
        }
        lines = [
            '候选评价次数：%d' % result.get('attempts', 0),
            '输入序列：%d' % result.get('background_count', 0),
            '', '各规则淘汰次数：']
        for key, count in sorted(
                result.get('rejections', {}).items(),
                key=lambda item: (-item[1], item[0])):
            lines.append('%s：%d' % (labels.get(key, key), count))
        return '\n'.join(lines)

    def _showSequenceMessage(self, icon, title, text, details=None):
        message = QMessageBox(icon, title, text,
                              QMessageBox.StandardButton.Ok,
                              self.win)
        if details:
            message.setDetailedText(details)
        message.exec()

    def actionSVGSlot(self):
        """docstring for actionSVGSlot"""
        fname = os.path.basename(str(self.filename()))
        if fname == None:
            directory = "."
        else:
            directory = QFileInfo(fname).path()

        fdialog = QFileDialog(
                    self.win,
                    "%s - Save As" % QApplication.applicationName(),
                    directory,
                    "%s (*.svg)" % QApplication.applicationName())
        fdialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
        fdialog.setWindowFlags(Qt.WindowType.Sheet)
        fdialog.setWindowModality(Qt.WindowModality.WindowModal)
        self.svgsavedialog = fdialog
        self.svgsavedialog.filesSelected.connect(self.saveSVGDialogCallback)
        fdialog.open()

    class DummyChild(QGraphicsItem):
        def boundingRect(self):
            return QRect(200, 200) # self.parentObject().boundingRect()
        def paint(self, painter, option, widget=None):
            pass

    def saveSVGDialogCallback(self, selected):
        if isinstance(selected, (list, tuple)):
            fname = selected[0]
        else:
            fname = selected
        if not fname or fname is None or os.path.isdir(fname):
            return False
        fname = str(fname)
        if not fname.lower().endswith(".svg"):
            fname += ".svg"
        if self.svgsavedialog != None:
            self.svgsavedialog.filesSelected.disconnect(self.saveSVGDialogCallback)
            del self.svgsavedialog  # prevents hang
            self.svgsavedialog = None

        generator = QSvgGenerator()
        generator.setFileName(fname)
        generator.setSize(QSize(200, 200))
        generator.setViewBox(QRect(0, 0, 2000, 2000))
        painter = QPainter()

        # Render through scene
        # painter.begin(generator)
        # self.win.pathscene.render(painter)
        # painter.end()

        # Render item-by-item
        painter = QPainter()
        styleOption = QStyleOptionGraphicsItem()
        q = [self.win.pathroot]
        painter.begin(generator)
        while q:
            graphicsItem = q.pop()
            transform = graphicsItem.itemTransform(self.win.sliceroot)[0]
            painter.setTransform(transform)
            if graphicsItem.isVisible():
                graphicsItem.paint(painter, styleOption, None)
                q.extend(graphicsItem.childItems())
        painter.end()

    def actionIllustratorSlot(self):
        """Export both occupied Slice View and Path View as editable SVG."""
        if self.activePart() is None:
            QMessageBox.warning(self.win, "Illustrator Export",
                                "The current document does not contain a design.")
            return

        currentName = self.filename()
        if currentName:
            info = QFileInfo(currentName)
            directory = info.absolutePath()
            baseName = info.completeBaseName()
        else:
            directory = "."
            baseName = "cadnano-design"
        defaultPath = os.path.join(directory, baseName + "-illustrator.svg")
        fname = QFileDialog.getSaveFileName(
            self.win,
            "%s - Export for Illustrator" % QApplication.applicationName(),
            defaultPath,
            "Illustrator-compatible SVG (*.svg)")
        if isinstance(fname, (list, tuple)):
            fname = fname[0]
        if not fname or os.path.isdir(fname):
            return
        if not str(fname).lower().endswith(".svg"):
            fname = str(fname) + ".svg"
        try:
            self._writeIllustratorSVG(str(fname))
        except Exception as error:
            QMessageBox.critical(
                self.win, "Illustrator Export",
                "The vector drawing could not be exported.\n\n%s" % error)
            return
        self.win.statusBar().showMessage(
            "Exported Slice + Path Illustrator SVG with editable Monaco "
            "sequence text to %s" % fname, 8000)

    def _illustratorContentRect(self):
        """Return the real Path View design bounds in scene coordinates."""
        root = self.win.pathroot
        contentRect = QRectF()
        for partItem in root.partItems():
            localRect = getattr(partItem, '_vHRect', None)
            if localRect is None or localRect.isNull():
                localRect = partItem.childrenBoundingRect()
            sceneRect = partItem.mapRectToScene(localRect)
            contentRect = (sceneRect if contentRect.isNull()
                           else contentRect.united(sceneRect))
        if contentRect.isNull():
            contentRect = root.mapRectToScene(root.childrenBoundingRect())
        # Hybrid endpoint labels extend horizontally beyond their anchored
        # bases; reserve extra room so neither text nor the rounded 3' accent
        # is clipped by the Illustrator artboard.
        margin = float(styles.PATH_BASE_WIDTH) * (
                                    5.0 if self._document.isHybrid() else 1.0)
        return contentRect.adjusted(-margin, -margin, margin, margin)

    def _sequenceTextItems(self):
        """Find the main strand sequence labels (insertions are paths already)."""
        labels = []
        for item in self.win.pathscene.items():
            if not isinstance(item, QGraphicsSimpleTextItem):
                continue
            parent = item.parentItem()
            if (parent is not None and
                    getattr(parent, '_seqLabel', None) is item and
                    item.isVisible() and item.text()):
                labels.append(item)
        return labels

    def _illustratorTransientItems(self):
        """Return visible Path View controls that are not design artwork."""
        items = []
        for partItem in self.win.pathroot.partItems():
            for attrName in ('_activeSliceItem', '_addBasesButton',
                             '_removeBasesButton', '_modRect'):
                item = getattr(partItem, attrName, None)
                if item is not None and item.isVisible():
                    items.append(item)
            for item in getattr(partItem, '_preXoverItems', []):
                if item is not None and item.isVisible():
                    items.append(item)
        # Hide only the selection rectangles. Selected strands may be
        # temporarily parented to these groups and must remain exportable.
        root = self.win.pathroot
        for groupName in ('_vhiHSelectionGroup', '_strandItemSelectionGroup'):
            group = getattr(root, groupName, None)
            box = getattr(group, 'selectionbox', None)
            if box is not None and box.isVisible():
                items.append(box)
        return items

    def _illustratorSliceData(self):
        """Return occupied Slice View bounds and items to hide for export."""
        root = self.win.sliceroot
        activePart = self.activePart()
        partItem = None
        for candidate in getattr(root, '_instanceItems', {}).keys():
            if candidate.part() is activePart:
                partItem = candidate
                break
        if partItem is None:
            return root, QRectF(), []

        virtualHelixItems = [item for item in
                             partItem._virtualHelixHash.values()
                             if item is not None and item.isVisible()]
        contentRect = QRectF()
        for item in virtualHelixItems:
            itemRect = item.sceneBoundingRect()
            contentRect = (itemRect if contentRect.isNull()
                           else contentRect.united(itemRect))
        if not contentRect.isNull():
            margin = float(styles.SLICE_HELIX_RADIUS) * 0.5
            contentRect = contentRect.adjusted(-margin, -margin,
                                               margin, margin)

        # Empty lattice sites and editing overlays are deliberately omitted.
        hiddenItems = []
        for emptyItem in partItem._emptyhelixhash.values():
            if (emptyItem.virtualHelixItem() is None and
                    emptyItem.isVisible()):
                hiddenItems.append(emptyItem)
        # Scaffold and staple orientation markers are design artwork and stay
        # visible in the Illustrator cross-section export.
        modifier = getattr(partItem, '_modCirc', None)
        if modifier is not None and modifier.isVisible():
            hiddenItems.append(modifier)
        return root, contentRect, hiddenItems

    @staticmethod
    def _otherVisibleTopLevelItems(scene, root):
        return [item for item in scene.items()
                if item.parentItem() is None and item is not root and
                item.isVisible()]

    @staticmethod
    def _drawEditableSequenceLabels(painter, sequenceLabels,
                                    sceneToExport):
        """Draw one positioned Monaco text object per base for Illustrator."""
        painter.resetTransform()
        for label in sequenceLabels:
            itemTransform = label.sceneTransform()

            def mapLocal(x, y):
                return sceneToExport.map(itemTransform.map(QPointF(x, y)))

            origin = mapLocal(0.0, 0.0)
            unitX = mapLocal(1.0, 0.0)
            unitY = mapLocal(0.0, 1.0)
            localToExport = QTransform(
                unitX.x() - origin.x(), unitX.y() - origin.y(),
                unitY.x() - origin.x(), unitY.y() - origin.y(),
                origin.x(), origin.y())
            painter.setWorldTransform(localToExport)
            originalFont = label.font()
            originalMetrics = QFontMetricsF(originalFont)
            advance = originalMetrics.horizontalAdvance('A')
            editableFont = QFont(originalFont)
            editableFont.setLetterSpacing(
                QFont.SpacingType.AbsoluteSpacing, 0.0)
            editableMetrics = QFontMetricsF(editableFont)
            painter.setFont(editableFont)
            painter.setPen(label.brush().color())
            baseline = editableMetrics.ascent()
            for index, character in enumerate(label.text()):
                if not character.isspace():
                    painter.drawText(QPointF(index * advance, baseline),
                                     character)
            painter.resetTransform()

    def _writeIllustratorSVG(self, fname):
        """Render occupied Slice View left of Path View in one vector SVG."""
        pathScene = self.win.pathscene
        pathRoot = self.win.pathroot
        pathRect = self._illustratorContentRect()
        if pathRect.isNull() or pathRect.width() <= 0 or pathRect.height() <= 0:
            raise ValueError("The Path View does not contain exportable artwork.")

        sliceRoot, sliceRect, hiddenSliceItems = self._illustratorSliceData()
        hasSlice = (not sliceRect.isNull() and sliceRect.width() > 0 and
                    sliceRect.height() > 0)
        gap = float(styles.PATH_BASE_WIDTH) * 2.0 if hasSlice else 0.0
        # The two on-screen views fit independently. Match their content
        # heights so the combined export preserves that left/right layout.
        sliceFit = (pathRect.height() / sliceRect.height()
                    if hasSlice else 1.0)
        sliceWidth = sliceRect.width() * sliceFit if hasSlice else 0.0
        layoutWidth = sliceWidth + gap + pathRect.width()
        layoutHeight = pathRect.height()

        # Keep large designs inside Illustrator's practical artboard range.
        longestSide = max(layoutWidth, layoutHeight)
        exportScale = min(1.0, 15000.0 / longestSide)
        width = max(1, int(ceil(layoutWidth * exportScale)))
        height = max(1, int(ceil(layoutHeight * exportScale)))

        generator = QSvgGenerator()
        generator.setFileName(fname)
        generator.setSize(QSize(width, height))
        generator.setViewBox(QRect(0, 0, width, height))
        generator.setResolution(96)
        generator.setTitle("cadnano Slice + Path Views")
        generator.setDescription(
            "Illustrator-compatible vector export. The left cross-section "
            "contains occupied virtual helices only; the right panel is the "
            "cadnano Path View with editable Monaco sequence letters.")

        hiddenTopLevelItems = self._otherVisibleTopLevelItems(pathScene,
                                                               pathRoot)
        if hasSlice:
            hiddenTopLevelItems.extend(self._otherVisibleTopLevelItems(
                                            self.win.slicescene, sliceRoot))
        for item in hiddenTopLevelItems:
            item.hide()

        hiddenTransientItems = self._illustratorTransientItems()
        for item in hiddenTransientItems:
            item.hide()
        for item in hiddenSliceItems:
            item.hide()

        # Always replace Qt's single spaced text run. SVG/Illustrator ignores
        # Qt's absolute letter-spacing metadata, which makes sequences bunch
        # together even when Monaco is installed.
        sequenceLabels = self._sequenceTextItems()
        for label in sequenceLabels:
            label.hide()

        painter = QPainter()
        try:
            if not painter.begin(generator):
                raise IOError("Qt could not open the SVG output file.")
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
            if hasSlice:
                sliceTarget = QRectF(
                    0.0, 0.0,
                    sliceWidth * exportScale,
                    layoutHeight * exportScale)
                self.win.slicescene.render(
                    painter, sliceTarget, sliceRect,
                    Qt.AspectRatioMode.IgnoreAspectRatio)

            pathOffsetX = (sliceWidth + gap) * exportScale
            pathTarget = QRectF(pathOffsetX, 0.0,
                                pathRect.width() * exportScale,
                                pathRect.height() * exportScale)
            pathScene.render(painter, pathTarget, pathRect,
                             Qt.AspectRatioMode.IgnoreAspectRatio)

            pathSceneToExport = QTransform(
                exportScale, 0.0, 0.0, exportScale,
                pathOffsetX - pathRect.left() * exportScale,
                -pathRect.top() * exportScale)
            self._drawEditableSequenceLabels(
                painter, sequenceLabels, pathSceneToExport)
            if self._document.isHybrid():
                # The two on-screen Hybrid panels use a shared widget
                # overlay.  Recreate the same base-anchored vector paths in
                # export coordinates so Illustrator receives real Beziers,
                # not a rasterized overlay or clipped per-panel stubs.
                painter.setWorldTransform(pathSceneToExport)
                for strand3Prime, strand5Prime in hybridConnections(
                                                    self._document):
                    point3 = hybridSceneEndpoint(pathRoot, strand3Prime,
                                                 '3p')
                    point5 = hybridSceneEndpoint(pathRoot, strand5Prime,
                                                 '5p')
                    if point3 is None or point5 is None:
                        continue
                    drawHybridConnection(
                        painter, point3, point5, strand3Prime,
                        selected=False, drawLabels=True)
                painter.resetTransform()
        finally:
            if painter.isActive():
                painter.end()
            for label in sequenceLabels:
                label.show()
            for item in hiddenTransientItems:
                item.show()
            for item in hiddenSliceItems:
                item.show()
            for item in hiddenTopLevelItems:
                item.show()

    def actionExportStaplesSlot(self):
        """
        Export input scaffold and output staple sequences to an XLSX workbook.
        """
        # Validate that no staple oligos are loops.
        part = self.activePart()
        if part is None:
            return
        if self._document.isHybrid():
            stapLoopOlgs = [oligo for oligo in self._document.oligos()
                            if oligo.isStaple() and oligo.isLoop()]
        else:
            stapLoopOlgs = part.getStapleLoopOligos()
        if stapLoopOlgs:
            from cadnano2.ui.dialogs.ui_warning import Ui_Warning
            dialog = QDialog()
            dialogWarning = Ui_Warning()
            dialog.setStyleSheet("QDialog { background-image: url(ui/dialogs/images/cadnano2-about.png); background-repeat: none; }")
            dialogWarning.setupUi(dialog, name="Warning-Circular")

            locs = ", ".join([o.locString() for o in stapLoopOlgs])
            msg = "Part contains staple loop(s) at %s.\n\nUse the break tool to introduce 5' & 3' ends before exporting. Loops have been colored red; use undo to revert." % locs
            dialogWarning.title.setText("Staple validation failed")
            dialogWarning.message.setText(msg)
            for o in stapLoopOlgs:
                o.applyColor(styles.stapColors[0].name())
            dialog.exec()
            return

        # Proceed with staple export.
        fname = self.filename()
        if fname == None:
            directory = "."
        else:
            directory = QFileInfo(fname).path()
        if util.isWindows():  # required for native looking file window
            fname = QFileDialog.getSaveFileName(
                            self.win,
                            "%s - Export As" % QApplication.applicationName(),
                            directory,
                            "Excel Workbook (*.xlsx)")
            self.saveStaplesDialog = None
            self.exportStaplesCallback(fname)
        else:  # access through non-blocking callback
            fdialog = QFileDialog(
                            self.win,
                            "%s - Export As" % QApplication.applicationName(),
                            directory,
                            "Excel Workbook (*.xlsx)")
            fdialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
            fdialog.setWindowFlags(Qt.WindowType.Sheet)
            fdialog.setWindowModality(Qt.WindowModality.WindowModal)
            self.saveStaplesDialog = fdialog
            self.saveStaplesDialog.filesSelected.connect(self.exportStaplesCallback)
            fdialog.open()
    # end def

    def actionExportCsvSlot(self):
        """Run cadnano's original staple-only CSV export workflow."""
        part = self.activePart()
        if part is None:
            return
        if self._document.isHybrid():
            stapLoopOlgs = [oligo for oligo in self._document.oligos()
                            if oligo.isStaple() and oligo.isLoop()]
        else:
            stapLoopOlgs = part.getStapleLoopOligos()
        if stapLoopOlgs:
            from cadnano2.ui.dialogs.ui_warning import Ui_Warning
            dialog = QDialog()
            dialogWarning = Ui_Warning()
            dialog.setStyleSheet(
                "QDialog { background-image: "
                "url(ui/dialogs/images/cadnano2-about.png); "
                "background-repeat: none; }")
            dialogWarning.setupUi(dialog, name="Warning-Circular")
            locs = ", ".join([oligo.locString()
                              for oligo in stapLoopOlgs])
            msg = ("Part contains staple loop(s) at %s.\n\nUse the break "
                   "tool to introduce 5' & 3' ends before exporting. "
                   "Loops have been colored red; use undo to revert." % locs)
            dialogWarning.title.setText("Staple validation failed")
            dialogWarning.message.setText(msg)
            for oligo in stapLoopOlgs:
                oligo.applyColor(styles.stapColors[0].name())
            dialog.exec()
            return

        fname = self.filename()
        directory = "." if fname is None else QFileInfo(fname).path()
        if util.isWindows():
            selected = QFileDialog.getSaveFileName(
                self.win,
                "%s - Export As" % QApplication.applicationName(),
                directory, "CSV (*.csv)")
            self.saveCsvDialog = None
            self.exportCsvCallback(selected)
        else:
            fdialog = QFileDialog(
                self.win,
                "%s - Export As" % QApplication.applicationName(),
                directory, "CSV (*.csv)")
            fdialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
            fdialog.setWindowFlags(Qt.WindowType.Sheet)
            fdialog.setWindowModality(Qt.WindowModality.WindowModal)
            self.saveCsvDialog = fdialog
            self.saveCsvDialog.filesSelected.connect(self.exportCsvCallback)
            fdialog.open()
    # end def

    def actionExportPdbSlot(self, checked=False, guided=False):
        """Export all-atom PDB/mmCIF and an optional oxDNA pair."""
        if not self._document.oligos():
            QMessageBox.warning(
                self.win, "结构导出" if guided else "PDB / oxDNA Export",
                ("当前设计中没有 DNA 链。" if guided else
                 "The current document does not contain any DNA strands."))
            return

        spacing, accepted = QInputDialog.getDouble(
            self.win, "结构导出" if guided else "PDB / oxDNA Export",
            ("螺旋中心间距（nm）：" if guided else
             "Helix center spacing (nm):"), 2.80, 1.00, 10.00, 2)
        if not accepted:
            return

        includeOxdna = True
        if guided:
            choice = QMessageBox.question(
                self.win, "可选的 oxDNA 文件",
                "是否同时导出用于 oxView 和 oxDNA 模拟的 .top 与 .dat 文件？\n\n"
                "选择“否”时只导出全原子 PDB；结构过大时自动改为 mmCIF。",
                QMessageBox.StandardButton.Yes |
                QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            includeOxdna = choice == QMessageBox.StandardButton.Yes

        currentName = self.filename()
        if currentName:
            info = QFileInfo(currentName)
            designName = info.completeBaseName()
        else:
            designName = "cadnano-design"
        desktop = os.path.join(QDir.homePath(), "Desktop")
        selected = QFileDialog.getExistingDirectory(
            self.win,
            (("请选择 %s 的导出文件夹" % designName) if guided else
             "%s - Choose Parent Folder for %s" %
             (QApplication.applicationName(), designName)), desktop)
        if isinstance(selected, (list, tuple)):
            selected = selected[0]
        if not selected:
            return
        selected = os.path.abspath(str(selected))
        if os.path.basename(selected) == designName:
            outputRoot = selected
        else:
            outputRoot = os.path.join(selected, designName)
        planned = structure_bundle_paths(
            outputRoot, designName, include_oxdna=includeOxdna)
        # Either PDB or mmCIF can be selected automatically after atom count
        # is known, so both possible structure paths participate in the
        # overwrite check.
        outputPaths = list(planned["all_atom"].values())
        if includeOxdna:
            outputPaths.extend(planned["oxdna"].values())
        existing = [path for path in outputPaths if os.path.exists(path)]
        if existing:
            box = QMessageBox(self.win)
            box.setWindowTitle("结构导出" if guided else
                               "PDB / oxDNA Export")
            box.setIcon(QMessageBox.Icon.Warning)
            box.setText("是否替换已有导出文件？" if guided else
                        "Replace the existing export file(s)?")
            box.setInformativeText("\n".join(existing))
            replaceButton = box.addButton(
                "替换" if guided else "Replace",
                QMessageBox.ButtonRole.AcceptRole)
            box.addButton("取消" if guided else "Cancel",
                          QMessageBox.ButtonRole.RejectRole)
            box.exec()
            if box.clickedButton() is not replaceButton:
                return
            del box

        try:
            summary = export_structure_bundle(
                self._document, outputRoot, designName,
                spacing_nm=spacing, include_oxdna=includeOxdna)
        except Exception as error:
            QMessageBox.critical(
                self.win, "结构导出" if guided else "PDB / oxDNA Export",
                (("无法导出当前结构。\n\n%s" % error) if guided else
                 "The model could not be exported.\n\n%s" % error))
            return

        if guided:
            message = ("已导出 %d 个核苷酸、%d 条链；螺旋中心间距 %.2f nm。"
                       "\n\n全原子 %s：%d 个原子\n%s\n\n%s" %
                       (summary["nucleotides"], summary["strands"], spacing,
                        summary["structure_format"],
                        summary["all_atom_count"],
                        ("oxDNA：已导出 .top 和 .dat" if includeOxdna else
                         "oxDNA：未导出"), summary["paths"]["root"]))
        else:
            message = ("Exported %d nucleotides in %d strands at %.2f nm "
                       "helix spacing.\n\nAll-atom %s: %d atoms\n%s\n\n%s" %
                       (summary["nucleotides"], summary["strands"], spacing,
                        summary["structure_format"],
                        summary["all_atom_count"], "oxDNA: .top + .dat",
                        summary["paths"]["root"]))
        if summary["assigned_bases"]:
            message += (("\n\n%d 个未指定碱基已按固定规则补全。" %
                         summary["assigned_bases"]) if guided else
                        ("\n\n%d unspecified bases were assigned "
                         "deterministically." % summary["assigned_bases"]))
        if summary["hybrid_residual"] > 1.0:
            message += ("\n\n跨点阵结构按 3D 视图间距放置；oxDNA "
                        "周期盒已扩展以保持两个点阵完整显示，"
                        "后续仍需要进行结构松弛。"
                        if guided else
                        "\n\nHybrid lattices use the 3D-view spacing; the "
                        "oxDNA periodic box was expanded to keep each lattice "
                        "in one image. Relaxation is still required.")
        QMessageBox.information(
            self.win, "结构导出" if guided else "PDB / oxDNA Export",
            message)
        self.win.statusBar().showMessage(
            "Exported all-atom %s%s to %s" %
            (summary["structure_format"],
             " + oxDNA" if includeOxdna else "", outputRoot), 10000)
    # end def

    def _updateAthenaActions(self):
        if self.win is not None:
            self.win.actionAthenaExport.setEnabled(
                bool(self._document.athenaMetadata()))
            self.win.actionCurvedExport.setEnabled(
                bool(self._document.curvedMetadata()))

    def actionAthenaDesignSlot(self):
        """Run an official ATHENA backend and open its saved cadnano JSON."""
        from ..views.athenadesign import AthenaDesignDialog
        dialog = AthenaDesignDialog(self.win)
        if not dialog.exec():
            return
        spec = dialog.spec()
        del dialog

        projectRoot = spec["project_root"]
        if os.path.isdir(projectRoot) and os.listdir(projectRoot):
            box = QMessageBox(self.win)
            box.setWindowTitle("Wireframe Design")
            box.setIcon(QMessageBox.Icon.Warning)
            box.setText("The selected project folder is not empty.")
            box.setInformativeText(
                "Known ATHENA output files with the same names will be "
                "replaced. Other files will be kept.\n\n%s" % projectRoot)
            continueButton = box.addButton(
                "Use Folder", QMessageBox.ButtonRole.AcceptRole)
            box.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
            box.exec()
            if box.clickedButton() is not continueButton:
                return
            del box

        if not self.maybeSave():
            return

        progress = QProgressDialog(
            "Preparing wireframe design…", "Cancel", 0, 0, self.win)
        progress.setWindowTitle("Wireframe Design")
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setMinimumDuration(0)
        progress.show()

        latestMessage = [""]

        def updateProgress(line):
            text = str(line).strip()
            if text and (text.startswith("+") or
                         text.lstrip().startswith(("*", "|"))):
                latestMessage[0] = text.strip("| ")
                progress.setLabelText(
                    latestMessage[0][:180] or "Running ATHENA…")
            QApplication.processEvents()

        try:
            result = create_project(
                spec, progress=updateProgress,
                cancelled=progress.wasCanceled)
        except Exception as error:
            progress.close()
            QMessageBox.critical(
                self.win, "Wireframe Design",
                "ATHENA could not complete the design.\n\n%s" % error)
            return
        progress.close()

        # The current document was already dealt with by maybeSave(), so the
        # generated design can be opened synchronously without another prompt.
        self.openAfterMaybeSaveCallback(result["json_path"])
        self._updateAthenaActions()
        QMessageBox.information(
            self.win, "Wireframe Design",
            "ATHENA design completed and opened in cadnano.\n\n"
            "JSON:\n%s\n\n"
            "Cylindrical, routing and pseudo-atomic BILD models were saved "
            "in:\n%s" %
            (result["json_path"],
             os.path.join(result["project_root"], "models")))
        self.win.statusBar().showMessage(
            "Wireframe design opened: %s" % result["json_path"], 12000)

    def _athenaSequenceReady(self):
        scaffoldOligos = [
            oligo for oligo in self._document.oligos()
            if not oligo.isStaple()]
        if not scaffoldOligos:
            return False
        for oligo in scaffoldOligos:
            sequence = (oligo.sequence() or "").upper()
            if not sequence or any(char not in "ACGT" for char in sequence):
                return False
        return True

    def actionAthenaExportSlot(self):
        """Export current sequence/topology in the saved ATHENA geometry."""
        metadata = self._document.athenaMetadata()
        if not metadata:
            QMessageBox.warning(
                self.win, "Wireframe 3D Export",
                "The current JSON does not contain ATHENA 3D mapping data.")
            return
        if not self._athenaSequenceReady():
            QMessageBox.warning(
                self.win, "Wireframe 3D Export",
                "Assign a complete scaffold sequence with Seq before "
                "exporting the sequence-accurate ATHENA model.")
            return

        currentName = self.filename()
        designName = metadata.get("output_name") or wireframe_output_name(
            metadata.get("name") or
            (QFileInfo(currentName).completeBaseName()
             if currentName else "wireframe-design"),
            metadata.get("edge_type", "DX"))
        projectRoot = metadata.get("project_root")
        if not projectRoot or not os.path.isdir(projectRoot):
            projectRoot = (QFileInfo(currentName).path()
                           if currentName else
                           os.path.join(QDir.homePath(), "Desktop"))
        defaultFolder = os.path.join(
            str(projectRoot), designName + "_3D_Export")
        os.makedirs(defaultFolder, exist_ok=True)
        selected = QFileDialog.getExistingDirectory(
            self.win, "Wireframe 3D Export — Choose output folder",
            defaultFolder)
        if isinstance(selected, (tuple, list)):
            selected = selected[0]
        if not selected:
            return
        selected = os.path.abspath(str(selected))

        exportName = designName
        planned = athena_structure_bundle_paths(selected, exportName)
        existing = [path for key, path in planned.items()
                    if key != "root" and os.path.exists(path)]
        if existing:
            box = QMessageBox(self.win)
            box.setWindowTitle("Wireframe 3D Export")
            box.setIcon(QMessageBox.Icon.Warning)
            box.setText("Files with this design name already exist.")
            box.setInformativeText("\n".join(existing))
            replaceButton = box.addButton(
                "Replace", QMessageBox.ButtonRole.AcceptRole)
            versionButton = box.addButton(
                "Create New Version", QMessageBox.ButtonRole.ActionRole)
            box.addButton("Cancel", QMessageBox.ButtonRole.RejectRole)
            box.exec()
            clicked = box.clickedButton()
            if clicked is versionButton:
                version = 2
                while True:
                    candidate = "%s_v%d" % (designName, version)
                    paths = athena_structure_bundle_paths(
                        selected, candidate)
                    if not any(os.path.exists(path)
                               for key, path in paths.items()
                               if key != "root"):
                        exportName = candidate
                        break
                    version += 1
            elif clicked is not replaceButton:
                return
            del box

        try:
            summary = export_athena_structure_bundle(
                self._document, metadata, selected, exportName)
        except Exception as error:
            QMessageBox.critical(
                self.win, "Wireframe 3D Export",
                "The ATHENA 3D model could not be exported.\n\n%s" % error)
            return
        QMessageBox.information(
            self.win, "Wireframe 3D Export",
            "Exported %d nucleotides in %d strands.\n\n"
            "%s:\n%s\n\noxDNA:\n%s\n%s\n\n"
            "This is an unrelaxed ATHENA target-geometry configuration." %
            (summary["nucleotides"], summary["strands"],
             summary["structure_format"], summary["structure_path"],
             summary["paths"]["top"], summary["paths"]["dat"]))
        self.win.statusBar().showMessage(
            "Wireframe 3D Export saved to %s" % selected, 12000)

    def actionCurvedDesignSlot(self):
        """Run DNAxiS in a child process and open its editable JSON."""
        from ..views.curveddesign import CurvedDesignDialog
        dialog = CurvedDesignDialog(self.win)
        if not dialog.exec():
            return
        spec = dialog.spec()
        del dialog

        if not self.maybeSave():
            return
        progress = QProgressDialog(
            "正在准备DNAxiS曲面设计…", "取消", 0, 0, self.win)
        progress.setWindowTitle("Curved Design")
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setMinimumDuration(0)
        progress.show()

        def updateProgress(line):
            text = str(line).strip()
            if text:
                progress.setLabelText(text[-180:])
            QApplication.processEvents()

        try:
            result = create_curved_project(
                spec, progress=updateProgress,
                cancelled=progress.wasCanceled)
        except Exception as error:
            progress.close()
            QMessageBox.critical(
                self.win, "Curved Design",
                "DNAxiS无法完成该设计。\n\n%s" % error)
            return
        progress.close()
        try:
            self.openAfterMaybeSaveCallback(result["json_path"])
        except Exception as error:
            QMessageBox.critical(
                self.win, "Curved Design",
                "设计文件已经生成，但cadnano无法安全打开它。\n\n%s" %
                error)
            return
        self._updateAthenaActions()
        indelSummary = result.get("indel_summary", {})
        maximumInsertion = int(indelSummary.get(
            "maximum_insertion_per_domain", 0))
        maximumDeletion = int(indelSummary.get(
            "maximum_deletion_per_domain", 0))
        curvedMetadata = result.get("metadata", {})
        domainSize = int(indelSummary.get("domain_size_bp", (
            7 if curvedMetadata.get("lattice") == "honeycomb" else 8)))
        actualHeight = float(curvedMetadata.get(
            "actual_outer_height_nm", 0.0))
        actualMinimumDiameter = float(curvedMetadata.get(
            "actual_minimum_outer_diameter_nm", 0.0))
        actualMaximumDiameter = float(curvedMetadata.get(
            "actual_maximum_outer_diameter_nm", 0.0))
        requestedHeight = float(curvedMetadata.get("height_nm", 0.0))
        requestedMinimumDiameter = float(curvedMetadata.get(
            "minimum_diameter_nm", 0.0))
        requestedMaximumDiameter = float(curvedMetadata.get(
            "maximum_diameter_nm", 0.0))
        minimumRingBp = int(curvedMetadata.get("minimum_ring_bp", 0))
        maximumRingBp = int(curvedMetadata.get("maximum_ring_bp", 0))
        scaffoldDensityMode = str(curvedMetadata.get(
            "scaffold_crossover_density_mode", "maximum"))
        scaffoldRequestedSpacing = int(curvedMetadata.get(
            "scaffold_crossover_requested_spacing_bp", 0))
        if scaffoldDensityMode != "minimum" and not scaffoldRequestedSpacing:
            scaffoldRequestedSpacing = (
                21 if curvedMetadata.get("lattice") == "honeycomb" else 32)
        scaffoldNativeMinimum = int(curvedMetadata.get(
            "scaffold_crossover_native_spacing_minimum_bp", 0))
        scaffoldNativeMaximum = int(curvedMetadata.get(
            "scaffold_crossover_native_spacing_maximum_bp", 0))
        scaffoldDensityText = (
            "最低路由密度" if scaffoldDensityMode == "minimum" else
            "目标1/%d bp" % scaffoldRequestedSpacing)
        scaffoldDensityText += "；最终原生坐标间距%d–%d bp" % (
            scaffoldNativeMinimum, scaffoldNativeMaximum)
        plannedIndels = dict(curvedMetadata.get(
            "planned_indel_summary") or {})
        parentReference = int(plannedIndels.get(
            "optimized_parent_reference_bp", 0))
        commonNominal = int(plannedIndels.get("common_nominal_bp", 0))
        summaryLines = [
            "DNA圆环：%d" % result["ring_count"],
            "目标外高：%.2f nm" % requestedHeight,
            "目标外径：%.2f–%.2f nm" % (
                requestedMinimumDiameter, requestedMaximumDiameter),
            "生成后实际外高：%.2f nm" % actualHeight,
            "生成后实际全部helix外径：%.2f–%.2f nm" % (
                actualMinimumDiameter, actualMaximumDiameter),
            "单环长度：%d–%d bp" % (minimumRingBp, maximumRingBp),
            "Scaffold crossover密度：%s" % scaffoldDensityText,
            "Indel domain：%d bp；硬上限：±3/domain" % domainSize,
            "最大insertion/domain：%d" % maximumInsertion,
            "最大deletion/domain：%d" % maximumDeletion,
            ("公共parent：参考%d bp → 整周期%d bp；始终最小化跨helix"
             "最大绝对indel负载" % (parentReference, commonNominal)),
            "AutoCS后crossover：固定不变（不再次增加或删除）"]
        outputName = str(curvedMetadata.get("output_name") or
                         os.path.splitext(os.path.basename(
                             result["json_path"]))[0])
        reportPath = os.path.join(
            result["project_root"], outputName + "_design_report.png")
        try:
            reportProcess = start_curved_report_export(
                result["json_path"], summaryLines, reportPath)
            self._watchCurvedReport(reportPath, reportProcess)
        except Exception as error:
            QMessageBox.warning(
                self.win, "Curved Design",
                "曲面设计已完成并在cadnano中打开。\n\n%s\n\n"
                "无法启动后台报告导出：%s" % (
                    "\n".join(summaryLines), error))
        self.win.statusBar().showMessage(
            "Curved design opened; report exporting in background: %s" %
            reportPath, 12000)

    def actionFrameDesignSlot(self):
        """Create a straight-edge frame with vertex-local indels."""
        from ..views.framedesign import FrameDesignDialog
        dialog = FrameDesignDialog(self.win)
        if not dialog.exec():
            return
        spec = dialog.spec()
        del dialog
        if not self.maybeSave():
            return
        progress = QProgressDialog(
            "正在生成Frame闭环与顶点局部曲率…", "取消", 0, 0, self.win)
        progress.setWindowTitle("Frame Design")
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setMinimumDuration(0)
        progress.show()

        def updateProgress(line):
            text = str(line).strip()
            if text:
                progress.setLabelText(text[-180:])
            QApplication.processEvents()

        try:
            result = create_frame_project(
                spec, progress=updateProgress,
                cancelled=progress.wasCanceled)
        except Exception as error:
            progress.close()
            QMessageBox.critical(
                self.win, "Frame Design",
                "Frame Design无法完成该设计。AutoCS结果未被修改。\n\n%s" %
                error)
            return
        progress.close()
        try:
            self.openAfterMaybeSaveCallback(result["json_path"])
        except Exception as error:
            QMessageBox.critical(
                self.win, "Frame Design",
                "Frame设计已生成，但cadnano无法安全打开。\n\n%s" % error)
            return
        self._updateAthenaActions()
        plan = result["frame_plan"]
        summaryLines = [
            "Frame边数：%d" % len(plan["vertices_nm"]),
            "原生闭环：%d bp（%.2f nm）" % (
                plan["nominal_perimeter_bp"],
                plan["nominal_perimeter_nm"]),
            "顶点转角：%s" % ", ".join(
                "%.1f°" % value for value in
                plan["turn_angles_degrees"]),
            "弯曲窗口：%s bp" % ", ".join(
                str(value) for value in plan["bend_length_bp"]),
            "最大insertion/domain：%d" % int(result[
                "indel_summary"]["maximum_insertion_per_domain"]),
            "最大deletion/domain：%d" % int(result[
                "indel_summary"]["maximum_deletion_per_domain"]),
            "AutoCS/crossover：完全沿用Curved Design并保持不变",
            "Indel：仅在顶点弯曲窗口内重新定位；直边为0曲率"]
        metadata = result.get("metadata", {})
        outputName = str(metadata.get("output_name") or
                         os.path.splitext(os.path.basename(
                             result["json_path"]))[0])
        reportPath = os.path.join(
            result["project_root"], outputName + "_design_report.png")
        try:
            reportProcess = start_curved_report_export(
                result["json_path"], summaryLines, reportPath)
            self._watchCurvedReport(reportPath, reportProcess)
        except Exception as error:
            QMessageBox.warning(
                self.win, "Frame Design",
                "Frame设计已完成。\n\n%s\n\n报告导出失败：%s" %
                ("\n".join(summaryLines), error))
        self.win.statusBar().showMessage(
            "Frame design opened; report exporting in background: %s" %
            reportPath, 12000)

    def actionTwistBendSlot(self):
        """Edit selected single-lattice regions with distributed indels."""
        if self._document.isHybrid() or len(self._document.parts()) != 1:
            QMessageBox.information(
                self.win, "Twist and Bend",
                "第一版仅支持普通单点阵 Square 或 Honeycomb 设计。\n"
                "Hybrid 暂不支持；Curved Design JSON 可作为普通点阵打开后使用。")
            return
        part = self._activePart or self._document.selectedPart()
        if part is None or not part.getVirtualHelices():
            QMessageBox.information(
                self.win, "Twist and Bend", "请先打开或创建一个含 helix 的设计。")
            return
        from ..views.twistbend import TwistBendDialog
        dialog = TwistBendDialog(self._document, part, self.win)
        if not dialog.exec():
            del dialog
            return
        try:
            result = dialog.resultData()
        except Exception as error:
            del dialog
            QMessageBox.critical(self.win, "Twist and Bend", str(error))
            return
        del dialog

        # Resolve every target before changing the model.  This guarantees
        # that a bad/stale region cannot leave a partially applied edit.
        targets = []
        for edit in result['edits']:
            vh = part.virtualHelix(int(edit['helix']))
            idx = int(edit['idx'])
            strand = None
            if vh is not None:
                strand = vh.scaffoldStrandSet().getStrand(idx)
                if strand is None:
                    strand = vh.stapleStrandSet().getStrand(idx)
            if strand is None:
                QMessageBox.critical(
                    self.win, "Twist and Bend",
                    "helix %s[%d] 不再包含可编辑链；设计可能在对话框打开后发生了变化。" %
                    (edit['helix'], idx))
                return
            operation = edit.get('operation', 'add')
            if operation == 'remove_existing':
                if not strand.hasInsertionAt(idx):
                    QMessageBox.critical(
                        self.win, "Twist and Bend",
                        "helix %s[%d] 的 insertion/deletion 已不存在；设计可能在"
                        "对话框打开后发生了变化。" % (edit['helix'], idx))
                    return
            elif strand.hasInsertionAt(idx):
                QMessageBox.critical(
                    self.win, "Twist and Bend",
                    "helix %s[%d] 已有 insertion/deletion；未覆盖现有修改。" %
                    (edit['helix'], idx))
                return
            targets.append((strand, idx, int(edit['length']), operation))

        stack = self._document.undoStack()
        # Adding every indel through Strand.addInsertion() would clear the
        # complete oligo sequence once per edit and create a nested undo macro
        # each time.  Bend/Add Twist commonly produce hundreds of edits, so
        # that path becomes extremely slow and may make Qt appear to hang.
        # Clear each affected oligo once, then push the primitive insertion
        # commands into one atomic user-visible macro.
        from ..model.strand import Strand
        affected_oligos = {}
        for strand, idx, length, operation in targets:
            affected_oligos[id(strand.oligo())] = strand.oligo()
            for complement in strand.getComplementStrands():
                affected_oligos[id(complement.oligo())] = complement.oligo()
        stack.beginMacro("Twist and Bend")
        applied = False
        try:
            stack.push(_SetTwistBendMetadataCommand(
                self._document, result['metadata']))
            for oligo in affected_oligos.values():
                stack.push(oligo.applySequenceCMD(None))
            for strand, idx, length, operation in targets:
                if operation == 'remove_existing':
                    stack.push(Strand.RemoveInsertionCommand(strand, idx))
                else:
                    stack.push(Strand.AddInsertionCommand(strand, idx, length))
            applied = True
        except Exception as error:
            apply_error = error
        finally:
            stack.endMacro()
        if not applied:
            # endMacro() commits the commands pushed before the exception;
            # immediately undo that one macro so the model is never left in a
            # partially modified state.
            if stack.canUndo():
                stack.undo()
            QMessageBox.critical(
                self.win, "Twist and Bend",
                "应用 indel 失败，已撤回本次全部修改。\n\n%s" % apply_error)
            return
        insertions = sum(1 for unused_s, unused_i, length, operation in targets
                         if operation != 'remove_existing' and length > 0)
        deletions = sum(1 for unused_s, unused_i, length, operation in targets
                        if operation != 'remove_existing' and length < 0)
        removed = sum(1 for unused_s, unused_i, unused_length, operation in
                      targets if operation == 'remove_existing')
        self.win.statusBar().showMessage(
            "Twist and Bend applied: %d insertion, %d deletion, %d existing "
            "indel removed; one Undo step." %
            (insertions, deletions, removed), 12000)
        # Reports are deliberately final-only: changing a slider or pressing
        # Apply Parameters must remain interactive and must not write files.
        # The accepted plans already contain the final pair-aware refinement.
        try:
            plans = list(result.get('metadata', {}).get('last_plans', ()))
            single_rows = single_helix_distribution_data(plans)
            current_name = self.filename()
            if current_name:
                report_directory = os.path.dirname(os.path.abspath(
                    str(current_name)))
                report_stem = os.path.splitext(os.path.basename(
                    str(current_name)))[0]
            else:
                report_directory = os.path.expanduser('~/Desktop')
                report_stem = 'untitled'
            report_base = os.path.join(
                report_directory, report_stem + '_twist_bend_report')
            report_path = report_base + '.png'
            write_single_helix_distribution_csv(
                single_rows, report_base + '_single_helix_distribution.csv')
            bend_index = 0
            for plan in plans:
                if plan.get('kind') != 'bend':
                    continue
                bend_index += 1
                suffix = '_bend_%02d_pair_curvature' % bend_index
                write_pair_curvature_csv(
                    plan.get('pair_curvature_rows', ()),
                    report_base + suffix + '.csv')
                write_pair_curvature_svg(
                    plan.get('pair_curvature_rows', ()),
                    plan.get('pair_curvature_summary', {}),
                    report_base + suffix + '.svg')
            image = create_twistbend_report_image(
                plans, single_rows, report_path)
            report_dialog = TwistBendReportDialog(image, self.win)
            report_dialog.setModal(False)
            report_dialog.setAttribute(
                Qt.WidgetAttribute.WA_DeleteOnClose, True)
            self._curvedReportDialogs.append(report_dialog)

            def twistBendReportClosed(unused_object=None):
                if report_dialog in self._curvedReportDialogs:
                    self._curvedReportDialogs.remove(report_dialog)

            report_dialog.destroyed.connect(twistBendReportClosed)
            report_dialog.show()
            report_dialog.raise_()
            self.win.statusBar().showMessage(
                'Twist and Bend applied; final report saved: %s' %
                report_path, 12000)
        except Exception as report_error:
            QMessageBox.warning(
                self.win, 'Twist and Bend',
                'Twist/Bend 已应用，但最终分布报告导出失败：\n%s' %
                report_error)

    def _watchCurvedReport(self, reportPath, reportProcess=None):
        """Show a finished background report without blocking design edits."""
        # ``openAfterMaybeSaveCallback`` can replace the document window
        # immediately before a Frame/Curved report finishes.  Parenting the
        # watcher and dialog to that old window silently destroyed the popup
        # even though the PNG was written successfully.  Keep the watcher on
        # the application and open a real top-level report window instead.
        application = QApplication.instance()
        timer = QTimer(application)
        timer.setInterval(600)
        self._curvedReportWatchers.append(timer)

        def reportReady():
            if not os.path.isfile(reportPath):
                if (reportProcess is not None and
                        reportProcess.poll() is not None):
                    timer.stop()
                    if timer in self._curvedReportWatchers:
                        self._curvedReportWatchers.remove(timer)
                    returnCode = reportProcess.returncode
                    QMessageBox.warning(
                        self.win, "Curved / Frame Design",
                        "设计已经生成，但后台报告导出失败（返回码 %s）。\n\n"
                        "预期报告：%s" % (returnCode, reportPath))
                    self.win.statusBar().showMessage(
                        "Curved / Frame Design report export failed.",
                        12000)
                return
            image = QImage(reportPath)
            if image.isNull():
                return
            timer.stop()
            if timer in self._curvedReportWatchers:
                self._curvedReportWatchers.remove(timer)
            dialog = CurvedReportDialog(image, None)
            # The report is generated in the background so the design opens
            # first.  Once ready it is intentionally window-modal: on macOS a
            # non-modal child completed by a detached process can remain
            # behind the main cadnano window even after raise_/activateWindow.
            # ``open`` below is asynchronous and therefore does not block the
            # event loop or prevent editing after the user closes the report.
            dialog.setWindowModality(Qt.WindowModality.ApplicationModal)
            dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
            self._curvedReportDialogs.append(dialog)

            def dialogClosed(unused_object=None):
                if dialog in self._curvedReportDialogs:
                    self._curvedReportDialogs.remove(dialog)

            dialog.destroyed.connect(dialogClosed)
            dialog.open()
            dialog.raise_()
            dialog.activateWindow()
            QApplication.alert(dialog, 0)
            self.win.statusBar().showMessage(
                "Curved / Frame Design report saved: %s" % reportPath,
                12000)

        timer.timeout.connect(reportReady)
        timer.start()
        reportReady()

    def actionCurvedExportSlot(self):
        """Export final sequence/topology in embedded DNAxiS coordinates."""
        metadata = self._document.curvedMetadata()
        if not metadata:
            QMessageBox.warning(
                self.win, "Curved / Frame 3D Export",
                "当前JSON不包含Curved或Frame真实三维坐标。")
            return
        if not self._athenaSequenceReady():
            QMessageBox.warning(
                self.win, "Curved / Frame 3D Export",
                "请先使用Import sequence为全部scaffold赋予完整序列。")
            return
        currentName = self.filename()
        designName = metadata.get("output_name") or curved_output_name(
            metadata.get("name") or
            (QFileInfo(currentName).completeBaseName()
             if currentName else "curved-design"),
            metadata.get("shape", "curved"),
            metadata.get("lattice"), metadata.get("layers"),
            metadata.get("height_nm"),
            metadata.get("maximum_diameter_nm"),
            metadata.get("minimum_diameter_nm"))
        projectRoot = metadata.get("project_root")
        if not projectRoot or not os.path.isdir(projectRoot):
            projectRoot = (QFileInfo(currentName).path()
                           if currentName else
                           os.path.join(QDir.homePath(), "Desktop",
                                        designName))
        structureDir = os.path.join(str(projectRoot), "structure")
        os.makedirs(structureDir, exist_ok=True)
        exportName = designName
        planned = curved_structure_bundle_paths(structureDir, exportName)
        existing = [path for key, path in planned.items()
                    if key != "root" and os.path.exists(path)]
        if existing:
            version = 1
            while True:
                candidate = "%s_%d" % (designName, version)
                paths = curved_structure_bundle_paths(
                    structureDir, candidate)
                if not any(os.path.exists(path)
                           for key, path in paths.items()
                           if key != "root"):
                    exportName = candidate
                    break
                version += 1
        try:
            summary = export_curved_structure_bundle(
                self._document, metadata, structureDir, exportName)
        except Exception as error:
            QMessageBox.critical(
                self.win, "Curved / Frame 3D Export",
                "无法导出Curved或Frame真实三维结构。\n\n%s" % error)
            return
        QMessageBox.information(
            self.win, "Curved / Frame 3D Export",
            "已导出%d个核苷酸、%d条链。\n\n"
            "%s：\n%s\n\noxDNA：\n%s\n%s\n\n"
            "BILD模型：\n%s\n%s\n\n"
            "PDB/mmCIF和DAT保留Curved或Frame真实三维坐标，不是平行helix构型。" %
            (summary["nucleotides"], summary["strands"],
             summary["structure_format"], summary["structure_path"],
             summary["paths"]["top"], summary["paths"]["dat"],
             summary["paths"]["cylindrical_model"],
             summary["paths"]["routing_model_multi"]))
        self.win.statusBar().showMessage(
            "Curved / Frame 3D Export saved to %s" % structureDir, 12000)

    def actionPrefsSlot(self):
        app().prefsClicked()

    def _autoFillWithoutCrossovers(self, strandType, label):
        parts = self._automaticTargetParts()
        if not parts:
            return
        self.win.pathGraphicsView.setViewportUpdateOn(False)
        try:
            try:
                count = sum(part.autoFillWithoutCrossovers(strandType)
                            for part in parts)
            except Exception as error:
                # PyQt aborts the complete process when an exception escapes
                # a QAction slot.  Keep an automatic-tool failure inside the
                # window and leave a full traceback in the launcher log.
                traceback.print_exc()
                QMessageBox.critical(
                    self.win, label,
                    "%s could not complete:\n%s" % (label, error))
                return
        finally:
            self.win.pathGraphicsView.setViewportUpdateOn(True)
        self.win.statusBar().showMessage(
            "%s created %d strand(s)." % (label, count), 8000)

    def actionAutoScaffoldWithoutCSSlot(self):
        self._autoFillWithoutCrossovers(
                    StrandType.Scaffold, "Add scaffolds")

    def actionAutoScaffoldCrossoversSlot(self):
        parts = self._automaticTargetParts()
        if not parts:
            return
        # A newly drawn design can have the correct visual Path order before
        # that order has ever been serialized into the model.  Snapshot the
        # rows directly from the live right-hand view so AutoCS never falls
        # back to numeric helix order and misses pairs such as 0/3.
        pathRoot = self.win.pathGraphicsView.sceneRootItem
        for part in parts:
            try:
                orderedCoords = pathRoot.partItemForPart(
                    part).getOrderedVirtualHelixList()
            except (AttributeError, KeyError):
                continue
            part.setImportedVHelixOrder(
                orderedCoords, emitSignal=False)
        settings = self._autoScaffoldSettingsDialog(parts)
        if settings is None:
            return
        densityMultiple, routeOnly, rebuildExisting = settings
        self.win.pathGraphicsView.setViewportUpdateOn(False)
        try:
            try:
                results = [part.autoScaffoldCrossovers(
                    densityMultiple=densityMultiple, routeOnly=routeOnly,
                    rebuildExisting=rebuildExisting, returnDetails=True)
                    for part in parts]
            except Exception as error:
                traceback.print_exc()
                QMessageBox.critical(
                    self.win, 'AutoCS_scaffolds',
                    "AutoCS_scaffolds could not complete:\n%s" % error)
                return
        finally:
            self.win.pathGraphicsView.setViewportUpdateOn(True)
        failures = [result for result in results if not result['success']]
        if failures:
            message = '\n\n'.join(result['message'] for result in failures)
            if not rebuildExisting:
                message += ('\n\n现有 scaffold crossover 已保留；'
                            '需要全部按新规则重建时，请勾选'
                            '“允许重新安排现有 crossover”。')
            QMessageBox.warning(self.win, 'AutoCS_scaffolds', message)
        else:
            message = ' | '.join(result['message'] for result in results)
            self.win.statusBar().showMessage(message, 15000)

    def _autoScaffoldSettingsDialog(self, parts):
        """Ask for the regular-design scaffold crossover density cap."""
        dialog = QDialog(self.win)
        dialog.setWindowTitle('AutoCS_scaffolds 设置')
        layout = QFormLayout(dialog)
        latticeText = ' + '.join(
            'Honeycomb（21 bp）' if part._step == 21 else 'Square（32 bp）'
            for part in parts)
        layout.addRow('当前点阵：', QLabel(latticeText, dialog))
        density = QComboBox(dialog)
        base = parts[0]._step if len(parts) == 1 else None
        for multiplier in range(1, 3):
            if base is None:
                text = '点阵周期 × %d' % multiplier
            else:
                text = '最大密度 1/%d bp' % (base * multiplier)
            density.addItem(text, multiplier)
        density.addItem('自定义点阵周期倍数', -1)
        density.addItem(
            '最低密度（仅添加全局连接所需的 crossover）', 0)
        custom = QSpinBox(dialog)
        custom.setRange(1, 100)
        custom.setValue(3)
        custom.setSuffix(' × 点阵周期')
        custom.setEnabled(False)
        density.currentIndexChanged.connect(
            lambda unused_index: custom.setEnabled(
                                        density.currentData() == -1))
        existingCount = sum(len(part._existingScaffoldCrossoverRecords())
                            for part in parts)
        rebuild = QCheckBox(
            '允许 AutoCS 重新安排现有 scaffold crossover', dialog)
        rebuild.setChecked(existingCount > 0)
        rebuild.setEnabled(existingCount > 0)
        rebuild.setToolTip(
            '不勾选时，现有 crossover 被视为固定约束。')
        note = QLabel(
            'AutoCS 仅处理在左侧几何视图和右侧 Path 顺序中都相邻的 '
            'helix 对，不建立闭环、不重排 helix，也不延长或删除 scaffold。'
            '所选数值控制 crossover 密度；“最低密度”采用当前设计中'
            '理论上最大的可用 crossover 间距：每个双重相邻连续段的'
            '第一对 helix 先用左右两侧最远的两个 crossover 封闭，'
            '之后沿 Path 顺序在中间连接与两侧连接之间交替。'
            '同一 helix 上连接不同邻居的 crossover，Square 保持至少 '
            '8 bases、Honeycomb 保持至少 7 bases。Square 合法性允许时'
            '优先避开右侧坐标为 8 的倍数的位置。'
            '勾选重新安排时会先删除现有 scaffold crossover。'
            '现有 crossover：%d 个。' %
            existingCount, dialog)
        note.setWordWrap(True)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok |
            QDialogButtonBox.StandardButton.Cancel, parent=dialog)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addRow('最大 crossover 密度：', density)
        layout.addRow('自定义倍数：', custom)
        layout.addRow('', rebuild)
        layout.addRow('', note)
        layout.addRow(buttons)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return None
        selected = density.currentData()
        routeOnly = selected == 0
        multiplier = custom.value() if selected == -1 else max(1, selected)
        return multiplier, routeOnly, rebuild.isChecked()

    def actionAutoStapleWithoutCSSlot(self):
        self._autoFillWithoutCrossovers(
                    StrandType.Staple, "Add staples")

    def actionAutostapleSlot(self):
        parts = self._automaticTargetParts()
        if parts:
            self.win.pathGraphicsView.setViewportUpdateOn(False)
            try:
                skippedHybrid = []
                for part in parts:
                    if part.autoStaple() is False:
                        skippedHybrid.append(
                            'Honeycomb' if part._step == 21 else 'Square')
            finally:
                self.win.pathGraphicsView.setViewportUpdateOn(True)
            if skippedHybrid:
                self.win.statusBar().showMessage(
                    "AutoCS_staples preserved existing cross-lattice "
                    "staples and skipped: %s." % ', '.join(skippedHybrid),
                    10000)

    def actionAutoBreakStaplesSlot(self):
        parts = self._automaticTargetParts()
        if not parts:
            return
        self.win.pathGraphicsView.setViewportUpdateOn(False)
        try:
            results = [part.autoBreakStaples() for part in parts]
        finally:
            self.win.pathGraphicsView.setViewportUpdateOn(True)
        if len(results) > 1:
            message = " | ".join(
                ("%s: %d new nick(s), %d manual nick(s) preserved, "
                 "%d edge crossover(s) removed, %d skipped") %
                (result['lattice'].capitalize(), result['nicks'],
                 result['protected_nicks'], result['removed_xovers'],
                 result['skipped'])
                for result in results)
        else:
            result = results[0]
            if not result['supported']:
                message = "Autobreak staples does not support this lattice."
            elif result.get('already_applied'):
                message = ("Autobreak staples already ran for this staple "
                           "layout; no additional nicks were created.")
            elif result['square']:
                message = ("Autobreak staples (square) created %(nicks)d nick(s), "
                           "preserved %(protected_nicks)d manual nick(s), "
                           "removed %(removed_xovers)d short-edge crossover(s)"
                           " and skipped %(skipped)d unbreakable oligo(s)." % result)
            else:
                message = ("Autobreak staples (honeycomb) created %(nicks)d "
                           "nick(s), preserved %(protected_nicks)d manual "
                           "nick(s), removed %(removed_xovers)d short-edge "
                           "crossover(s)"
                           " and skipped %(skipped)d unbreakable oligo(s)." % result)
        self._updateContinuousStapleStatistics()
        self.win.statusBar().showMessage(message, 10000)

    def _automaticTargetParts(self):
        """Use both local rule engines in Hybrid mode, otherwise active part."""
        if self._document.isHybrid():
            byLattice = self._document.partsByLattice()
            return [byLattice[name] for name in ('honeycomb', 'square')
                    if name in byLattice]
        part = self.activePart()
        return [part] if part is not None else []

    def actionModifySlot(self):
        """
        Notifies that part root items that parts should respond to modifier
        selection signals.
        """
        # uncomment for debugging
        # isChecked = self.win.actionModify.isChecked()
        # self.win.pathroot.setModifyState(isChecked)
        # self.win.sliceroot.setModifyState(isChecked)
        if app().isInMaya():
            isChecked = self.win.actionModify.isChecked()
            self.win.pathroot.setModifyState(isChecked)
            self.win.sliceroot.setModifyState(isChecked)
            self.win.solidroot.setModifyState(isChecked)

    def actionAddHoneycombPartSlot(self):
        """docstring for actionAddHoneycombPartSlot"""
        part = self._document.addHoneycombPart()
        self.setActivePart(part)

    def actionAddSquarePartSlot(self):
        """docstring for actionAddSquarePartSlot"""
        part = self._document.addSquarePart()
        self.setActivePart(part)

    def actionAddHybridPartSlot(self):
        """Create both lattice parts and reveal the four-panel workspace."""
        parts = self._document.addHybridParts()
        if parts is None:
            self.win.statusBar().showMessage(
                "Start a new document before choosing Hybrid mode.", 8000)
            return
        honeycomb, unused_square = parts
        self.win.setHybridMode(True)
        self.setActivePart(honeycomb)

    def actionRenumberSlot(self):
        coordList = self.win.pathroot.getSelectedPartOrderedVHList()
        part = self.activePart()
        part.renumber(coordList)
    # end def

    ### ACCESSORS ###
    def document(self):
        return self._document

    def window(self):
        return self.win

    def setDocument(self, doc):
        """
        Sets the controller's document, and informs the document that
        this is its controller.
        """
        self._document = doc
        doc.setController(self)

    def activePart(self):
        if self._activePart == None:
            selectedPart = self._document.selectedPart()
            if selectedPart is not None:
                self.setActivePart(selectedPart)
        return self._activePart

    def setActivePart(self, part):
        previousPart = self._activePart
        if previousPart is not None and previousPart is not part:
            try:
                previousPart.partStrandChangedSignal.disconnect(
                                self._scheduleContinuousStapleStatistics)
            except (TypeError, RuntimeError):
                pass
        self._activePart = part
        if part is not None and previousPart is not part:
            part.partStrandChangedSignal.connect(
                                self._scheduleContinuousStapleStatistics)
        self._updateContinuousStapleStatistics()

    def _scheduleContinuousStapleStatistics(self, *unused_args):
        """Coalesce rapid model edits into one inexpensive UI refresh."""
        if self._continuousStapleStatsPending:
            return
        self._continuousStapleStatsPending = True
        QTimer.singleShot(0, self._updateContinuousStapleStatistics)

    def _updateContinuousStapleStatistics(self):
        self._continuousStapleStatsPending = False
        label = self._continuousStapleLabel
        if label is None:
            return
        if self._document.isHybrid():
            parts = self._document.partsByLattice()
            texts = []
            for lattice, minimum, labelText in (
                    ('honeycomb', 14, 'Honeycomb staples'),
                    ('square', 16, 'Square staples')):
                part = parts.get(lattice)
                if part is None:
                    continue
                qualified, total, fraction = \
                    part.stapleContinuousRunStatistics(minimum)
                unused_total, bins = part.stapleLengthDistribution(10, 10)
                binTexts = ["%d–%d: %.1f%%" %
                            (low, high, binFraction * 100.0)
                            for low, high, count, binFraction in bins
                            if count > 0]
                texts.append("%s ≥%d consecutive: %.1f%% (%d/%d); %s" %
                             (labelText, minimum, fraction * 100.0,
                              qualified, total,
                              " | ".join(binTexts) if binTexts else
                              "no staples"))
            label.setText("        ||        ".join(texts) +
                          "        |        Hybrid staples excluded")
            return
        part = self._activePart
        if part is None:
            label.setText(
                "Staples with ≥14 consecutive bases: 0.0% (0/0)"
                "        |        "
                "Staple lengths (% of all staples): no staples")
            return
        minimum = 14 if part._step == 21 else 16
        qualified, total, fraction = \
            part.stapleContinuousRunStatistics(minimum)
        unused_total, bins = part.stapleLengthDistribution(10, 10)
        binTexts = ["%d–%d: %.1f%% (%d)" %
                    (low, high, binFraction * 100.0, count)
                    for low, high, count, binFraction in bins
                    if count > 0]
        if binTexts:
            distributionText = " | ".join(binTexts)
        else:
            distributionText = "no staples ≥10 bases"
        label.setText(
            "Staples with ≥%d consecutive bases: %.1f%% (%d/%d)"
            "        |        "
            "Staple lengths (%% of all staples): %s" %
            (minimum, fraction * 100.0, qualified, total,
             distributionText))

    def undoStack(self):
        return self._document.undoStack()

    ### PRIVATE SUPPORT METHODS ###
    def newDocument(self, doc=None, fname=None):
        """Creates a new Document, reusing the DocumentController."""
        self.setActivePart(None)
        self._document.resetViews()
        self._document.removeAllParts()  # clear out old parts
        self._document.setAthenaMetadata(None)
        self._document.setTwistBendMetadata(None)
        self.win.setHybridMode(False)
        self._updateAthenaActions()
        self._document.undoStack().clear()  # reset undostack
        self._filename = fname if fname else "untitled.json"
        self._hasNoAssociatedFile = fname == None
        self.win.setWindowTitle(self.documentTitle() + '[*]')

    def saveFileDialog(self, includeSequences=False):
        if self.filesavedialog is not None:
            # A cancelled non-blocking sheet must not remain as a stale
            # target that swallows the next Save action.
            try:
                if self.filesavedialog.isVisible():
                    self.filesavedialog.raise_()
                    self.filesavedialog.activateWindow()
                    return
            except RuntimeError:
                pass
            staleDialog = self.filesavedialog
            self.filesavedialog = None
            try:
                staleDialog.filesSelected.disconnect(
                                                self.saveFileDialogCallback)
            except (TypeError, RuntimeError):
                pass
            try:
                staleDialog.finished.disconnect(
                                                self.saveFileDialogFinished)
            except (TypeError, RuntimeError):
                pass
            try:
                staleDialog.deleteLater()
            except RuntimeError:
                pass
        self._saveSequencesOnDialog = includeSequences
        fname = self.filename()
        if fname == None:
            directory = "."
        else:
            directory = QFileInfo(fname).path()
        if util.isWindows():  # required for native looking file window
            fname = QFileDialog.getSaveFileName(
                            self.win,
                            "%s - Save As" % QApplication.applicationName(),
                            directory,
                            "%s (*.json)" % QApplication.applicationName())
            self.writeDocumentToFile(fname,
                                     includeSequences=includeSequences)
        else:  # access through non-blocking callback
            fdialog = QFileDialog(
                            self.win,
                            "%s - Save As" % QApplication.applicationName(),
                            directory,
                            "%s (*.json)" % QApplication.applicationName())
            fdialog.setAcceptMode(QFileDialog.AcceptMode.AcceptSave)
            fdialog.setWindowFlags(Qt.WindowType.Sheet)
            fdialog.setWindowModality(Qt.WindowModality.WindowModal)
            self.filesavedialog = fdialog
            self.filesavedialog.filesSelected.connect(
                                                self.saveFileDialogCallback)
            self.filesavedialog.finished.connect(
                                                self.saveFileDialogFinished)
            fdialog.open()
    # end def

    def saveFileDialogFinished(self, unused_result):
        """Release Save As state after both acceptance and cancellation."""
        dialog = self.filesavedialog
        if dialog is None:
            return
        try:
            dialog.filesSelected.disconnect(self.saveFileDialogCallback)
        except (TypeError, RuntimeError):
            pass
        try:
            dialog.finished.disconnect(self.saveFileDialogFinished)
        except (TypeError, RuntimeError):
            pass
        self.filesavedialog = None
        self._saveSequencesOnDialog = False

    def _readSettings(self):
        self.settings.beginGroup("FileSystem")
        self._fileOpenPath = self.settings.value("openpath", QDir().homePath())
        self.settings.endGroup()

    def _writeFileOpenPath(self, path):
        """docstring for _writePath"""
        self._fileOpenPath = path
        self.settings.beginGroup("FileSystem")
        self.settings.setValue("openpath", path)
        self.settings.endGroup()

    ### SLOT CALLBACKS ###
    def actionNewSlotCallback(self):
        """
        Gets called on completion of filesavedialog after newClicked's
        maybeSave. Removes the dialog if necessary, but it was probably
        already removed by saveFileDialogCallback.
        """
        if self.filesavedialog != None:
            self.filesavedialog.finished.disconnect(self.actionNewSlotCallback)
            del self.filesavedialog  # prevents hang (?)
            self.filesavedialog = None
        self.newDocument()

    def exportStaplesCallback(self, selected):
        """Export input and output sequences to a two-sheet XLSX workbook.

        Args:
            selected (Tuple, List or str): if a List or Tuple, the filename
            should be the first element
        """
        if isinstance(selected, (list, tuple)):
            fname = selected[0]
        else:
            fname = selected
        # Return if fname is '', None, or a directory path
        if not fname or fname is None or os.path.isdir(fname):
            return False
        if not fname.lower().endswith(".xlsx"):
            fname += ".xlsx"
        if self.saveStaplesDialog is not None:
            self.saveStaplesDialog.filesSelected.disconnect(self.exportStaplesCallback)
            # manual garbage collection to prevent hang (in osx)
            del self.saveStaplesDialog
            self.saveStaplesDialog = None
        # write the file
        ap = self.activePart()
        if ap is not None:
            if self._document.isHybrid():
                inputRows = self._document.getInputSequenceRows()
                outputRows = self._document.getOutputSequenceRows()
            else:
                inputRows = ap.getInputSequenceRows()
                outputRows = ap.getOutputSequenceRows()
            write_sequence_workbook(fname, inputRows, outputRows)
    # end def

    def exportCsvCallback(self, selected):
        """Write the exact original cadnano staple CSV text."""
        if isinstance(selected, (list, tuple)):
            fname = selected[0]
        else:
            fname = selected
        if not fname or fname is None or os.path.isdir(fname):
            return False
        fname = str(fname)
        if not fname.lower().endswith(".csv"):
            fname += ".csv"
        if self.saveCsvDialog is not None:
            self.saveCsvDialog.filesSelected.disconnect(
                                                    self.exportCsvCallback)
            del self.saveCsvDialog
            self.saveCsvDialog = None
        part = self.activePart()
        if part is not None:
            if self._document.isHybrid():
                # The lower-right raw CSV intentionally keeps cadnano's
                # original endpoint syntax.  Lattice classification remains
                # exclusive to the XLSX exporter.
                output = "Start,End,Sequence,Length,Color\n"
                stapleOligos = [
                    oligo for oligo in self._document.oligos()
                    if oligo.isStaple()]
                stapleOligos.sort(
                    key=lambda oligo: (
                        oligo.strand5p().virtualHelix().number(),
                        oligo.strand5p().idx5Prime()))
                for oligo in stapleOligos:
                    output += oligo.sequenceExport(includeLattice=False)
            else:
                output = part.getStapleSequences()
            with open(fname, 'w') as csvFile:
                csvFile.write(output)
            self.win.statusBar().showMessage(
                "Exported original cadnano staple CSV to %s" % fname, 8000)
            return True
        return False
    # end def

    def newClickedCallback(self):
        """
        Gets called on completion of filesavedialog after newClicked's
        maybeSave. Removes the dialog if necessary, but it was probably
        already removed by saveFileDialogCallback.
        """

        if self.filesavedialog != None:
            self.filesavedialog.finished.disconnect(self.newClickedCallback)
            del self.filesavedialog  # prevents hang (?)
            self.filesavedialog = None
        self.newDocument()

    def openAfterMaybeSaveCallback(self, selected):
        """
        Receives file selection info from the dialog created by
        openAfterMaybeSave, following user input.

        Extracts the file name and passes it to the decode method, which
        returns a new document doc, which is then set as the open document
        by newDocument. Calls finalizeImport and disconnects dialog signaling.
        """
        if isinstance(selected, (list, tuple)):
            fname = selected[0]
        else:
            fname = selected
        if fname is None or fname == '' or os.path.isdir(fname):
            return False
        if not os.path.exists(fname):
            return False
        self._writeFileOpenPath(os.path.dirname(fname))
        self.newDocument(fname=fname)

        with io.open(fname, 'r', encoding='utf-8') as fd:
            decode(self._document, fd.read())
        self.win.setHybridMode(self._document.isHybrid())
        self.setActivePart(self._document.selectedPart())
        self._updateAthenaActions()

        if hasattr(self, "filesavedialog"):  # user did save
            if self.fileopendialog is not None:
                self.fileopendialog.filesSelected.disconnect(self.openAfterMaybeSaveCallback)
            # manual garbage collection to prevent hang (in osx)
            del self.fileopendialog
            self.fileopendialog = None

    def saveFileDialogCallback(self, selected):
        """If the user chose to save, write to that file."""
        if isinstance(selected, (list, tuple)):
            fname = selected[0]
        else:
            fname = selected
        if not fname or os.path.isdir(fname):
            self._saveSequencesOnDialog = False
            return False
        if not fname.lower().endswith(".json"):
            fname += ".json"
        includeSequences = self._saveSequencesOnDialog
        self.writeDocumentToFile(fname,
                                 includeSequences=includeSequences)
        self._writeFileOpenPath(os.path.dirname(fname))

    ### EVENT HANDLERS ###
    def windowCloseEventHandler(self, event):
        """Intercept close events when user attempts to close the window."""
        if self.maybeSave():
            event.accept()
            if app().isInMaya():
                self.windock.setVisible(False)
                del self.windock
                self.windock = None
            the_app = app()
            self.destroyDC()
            if the_app.documentControllers:
                the_app.destroyApp()
        else:
            event.ignore()
        self.actionCloseSlot()

    ### FILE INPUT ##
    def documentTitle(self):
        fname = os.path.basename(str(self.filename()))
        if not self.undoStack().isClean():
            fname += '[*]'
        return fname

    def filename(self):
        return self._filename

    def setFilename(self, proposedFName):
        if self._filename == proposedFName:
            return True
        self._filename = proposedFName
        self._hasNoAssociatedFile = False
        self.win.setWindowTitle(self.documentTitle())
        return True

    def openDropAfterMaybeSave(self):
        """
        This is the method that initiates file opening after a drag and Drop event.
        It is called by actionDropSlot.
        """
        fname = self.win._dropped_file
        if fname.endswith(".json"):
            self.openAfterMaybeSaveCallback(fname)
        else:
            print(f"Ignoring dropped file {fname}. Use a cadnano.json file.")

    def openAfterMaybeSave(self):
        """
        This is the method that initiates file opening. It is called by
        actionOpenSlot to spawn a QFileDialog and connect it to a callback
        method.
        """
        path = self._fileOpenPath
        if util.isWindows():  # required for native looking file window#"/",
            fname = QFileDialog.getOpenFileName(
                        None,
                        "Open Document", path,
                        "cadnano1 / cadnano2 Files (*.nno *.json *.cadnano)")
            self.filesavedialog = None
            self.openAfterMaybeSaveCallback(fname)
        else:  # access through non-blocking callback
            fdialog = QFileDialog(
                        self.win,
                        "Open Document",
                        path,
                        "cadnano1 / cadnano2 Files (*.nno *.json *.cadnano)")
            fdialog.setAcceptMode(QFileDialog.AcceptMode.AcceptOpen)
            fdialog.setWindowFlags(Qt.WindowType.Sheet)
            fdialog.setWindowModality(Qt.WindowModality.WindowModal)
            self.fileopendialog = fdialog
            self.fileopendialog.filesSelected.connect(self.openAfterMaybeSaveCallback)
            fdialog.open()
    # end def

    ### FILE OUTPUT ###
    def maybeSave(self):
        """Save on quit, check if document changes have occured."""
        if app().dontAskAndJustDiscardUnsavedChanges:
            return True
        if not self.undoStack().isClean():    # document dirty?
            part = self.activePart()
            sequenceCount = 0
            if self._document.isHybrid():
                sequenceCount = len(self._document.getInputSequenceRows())
            elif part is not None:
                sequenceCount = len(part.getScaffoldSequenceRecords())
            if sequenceCount:
                message = ("The document has been modified and contains "
                           "%d applied scaffold sequence(s).\n"
                           "Choose whether sequence information should be "
                           "included in the saved design." %
                           sequenceCount)
            else:
                message = ("The document has been modified.\n"
                           "Do you want to save your changes?")
            savebox = QMessageBox(QMessageBox.Icon.Question,   "Application",
                message,
                QMessageBox.StandardButton.NoButton,
                self.win,
                Qt.WindowType.Dialog | Qt.WindowType.MSWindowsFixedSizeDialogHint | Qt.WindowType.Sheet)
            savebox.setWindowModality(Qt.WindowModality.WindowModal)
            saveWithSequences = None
            if sequenceCount:
                saveWithSequences = savebox.addButton(
                    "Save with Sequences", QMessageBox.ButtonRole.AcceptRole)
                saveWithSequences.setShortcut("Ctrl+Alt+S")
                save = savebox.addButton(
                    "Save without Sequences", QMessageBox.ButtonRole.AcceptRole)
            else:
                save = savebox.addButton(
                    "Save", QMessageBox.ButtonRole.AcceptRole)
            discard = savebox.addButton(
                "Discard", QMessageBox.ButtonRole.DestructiveRole)
            cancel = savebox.addButton(
                "Cancel", QMessageBox.ButtonRole.RejectRole)
            save.setShortcut("Ctrl+S")
            discard.setShortcut(QKeySequence("D,Ctrl+D"))
            cancel.setShortcut(QKeySequence("C,Ctrl+C,.,Ctrl+."))
            savebox.exec()
            clickedButton = savebox.clickedButton()
            del savebox  # manual garbage collection to prevent hang (in osx)
            if clickedButton is saveWithSequences:
                return self._saveFromMaybeSave(includeSequences=True)
            elif clickedButton is save:
                return self._saveFromMaybeSave(includeSequences=False)
            elif clickedButton is cancel or clickedButton is None:
                return False
        return True

    def _saveFromMaybeSave(self, includeSequences=False):
        """Synchronously save during close/new/open confirmation."""
        if self._hasNoAssociatedFile:
            fname = QFileDialog.getSaveFileName(
                        self.win,
                        "%s - Save As" % QApplication.applicationName(),
                        QFileInfo(self.filename()).path(),
                        "%s (*.json)" % QApplication.applicationName())
            if isinstance(fname, (list, tuple)):
                fname = fname[0]
            if not fname or os.path.isdir(fname):
                return False
            if not fname.lower().endswith(".json"):
                fname += ".json"
            return self.writeDocumentToFile(
                        fname, includeSequences=includeSequences)
        return self.writeDocumentToFile(
                        includeSequences=includeSequences)

    def writeDocumentToFile(self, filename=None, includeSequences=False):
        helixOrderList = self.win.pathroot.getSelectedPartOrderedVHList()

        hasHybridContent = (self._document.isHybrid() and any(
            part.numberOfVirtualHelices() > 0
            for part in self._document.parts()))
        if helixOrderList == None and not hasHybridContent:
            print("Cannot save empty document.")
            return False
        if self._document.isHybrid():
            helixOrderList = []

        if filename == None:
            assert(not self._hasNoAssociatedFile)
            filename = self.filename()
        try:
            if util.isWindows() and isinstance(filename, (list,tuple)):
                filename = filename[0]
            with open(filename, 'w') as f:
                encode(self._document, helixOrderList, f,
                       includeSequences=includeSequences)
        except (IOError, OSError, TypeError):
            flags = Qt.WindowType.Dialog | Qt.MSWindowsFixedSizeDialogHint | Qt.WindowType.Sheet
            errorbox = QMessageBox(QMessageBox.Critical,
                                   "cadnano",
                                   "Could not write to '%s'." % filename,
                                   QMessageBox.Ok,
                                   self.win,
                                   flags)
            errorbox.setWindowModality(Qt.WindowModality.WindowModal)
            errorbox.open()
            return False
        self.undoStack().setClean()
        self.setFilename(filename)
        self.win.statusBar().showMessage(
            "Saved design to %s" % filename, 6000)
        return True

    def actionCadnanoWebsiteSlot(self):
        import webbrowser
        webbrowser.open("http://cadnano.org/")

    def actionFeedbackSlot(self):
        import webbrowser
        webbrowser.open("http://cadnano.org/feedback")
