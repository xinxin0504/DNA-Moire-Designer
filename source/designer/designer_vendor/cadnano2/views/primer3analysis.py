"""Independent Primer3 thermodynamic-analysis window."""

import os

from cadnano2.model.io.primer3analysis import (
    Primer3Cancelled, Primer3Unavailable, normalized_entries,
    read_primer3_sequences, run_primer3_analysis, write_primer3_workbook)
from cadnano2.views.primer3structure import Primer3StructureView
import cadnano2.util as util

util.qtWrapImport('QtCore', globals(), ['Qt'])
util.qtWrapImport('QtGui', globals(), ['QBrush', 'QColor', 'QFont'])
util.qtWrapImport('QtWidgets', globals(), [
    'QAbstractItemView', 'QApplication', 'QCheckBox', 'QDialog',
    'QDialogButtonBox', 'QFileDialog', 'QGroupBox',
    'QGridLayout', 'QHeaderView', 'QLabel', 'QMessageBox', 'QPlainTextEdit',
    'QProgressDialog', 'QPushButton', 'QSplitter', 'QTabWidget', 'QTableWidget',
    'QTableWidgetItem', 'QVBoxLayout', 'QWidget'])


class Primer3AnalysisDialog(QDialog):
    """Select short strands and inspect Primer3 thermodynamic structures."""

    def __init__(self, entries=(), suggested_directory='', parent=None):
        super(Primer3AnalysisDialog, self).__init__(parent)
        self.setWindowTitle('Primer3 Thermodynamic Analysis')
        self.resize(1440, 840)
        self._entries = normalized_entries(entries)
        self._results = []
        self._suggestedDirectory = suggested_directory or \
            os.path.expanduser('~/Desktop')

        outer = QVBoxLayout(self)
        self.mainSplitter = QSplitter(Qt.Orientation.Horizontal, self)
        self.mainSplitter.setChildrenCollapsible(False)
        outer.addWidget(self.mainSplitter, 1)

        sequencePanel = QWidget(self.mainSplitter)
        sequenceLayout = QVBoxLayout(sequencePanel)
        setupTitle = QLabel('<b>Analysis Setup</b>', sequencePanel)
        sequenceLayout.addWidget(setupTitle)

        self.setupSplitter = QSplitter(
            Qt.Orientation.Vertical, sequencePanel)
        self.setupSplitter.setChildrenCollapsible(False)
        sequenceLayout.addWidget(self.setupSplitter, 1)

        selectionPanel = QGroupBox('Sequence Selection', self.setupSplitter)
        selectionLayout = QVBoxLayout(selectionPanel)
        description = QLabel(
            'Analyze orthogonal short sequences independently of scaffold '
            'sequences and the current cadnano design. Select the sequences '
            'to include in the analysis.', selectionPanel)
        description.setWordWrap(True)
        selectionLayout.addWidget(description)

        sourceButtons = QGridLayout()
        self.importButton = QPushButton(
            'Import TXT / Orthogonal Sequence XLSX…', selectionPanel)
        self.selectAllButton = QPushButton('Select All', selectionPanel)
        self.selectNoneButton = QPushButton('Select None', selectionPanel)
        self.removeButton = QPushButton(
            'Remove Selected Rows', selectionPanel)
        sourceButtons.addWidget(self.importButton, 0, 0, 1, 2)
        sourceButtons.addWidget(self.selectAllButton, 1, 0)
        sourceButtons.addWidget(self.selectNoneButton, 1, 1)
        sourceButtons.addWidget(self.removeButton, 2, 0, 1, 2)
        selectionLayout.addLayout(sourceButtons)

        self.sequenceTable = QTableWidget(selectionPanel)
        self.sequenceTable.setColumnCount(5)
        self.sequenceTable.setHorizontalHeaderLabels(
            ('Analyze', 'Source', 'Name', 'Sequence (5′→3′)', 'Length (nt)'))
        self.sequenceTable.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.sequenceTable.setSelectionMode(
            QAbstractItemView.SelectionMode.ExtendedSelection)
        self.sequenceTable.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        header = self.sequenceTable.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setMinimumSectionSize(56)
        for column, width in enumerate((76, 125, 145, 310, 88)):
            self.sequenceTable.setColumnWidth(column, width)
        self.sequenceTable.verticalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Interactive)
        selectionLayout.addWidget(self.sequenceTable, 1)
        self.setupSplitter.addWidget(selectionPanel)

        configurationPanel = QGroupBox(
            'Analysis Configuration', self.setupSplitter)
        configurationLayout = QVBoxLayout(configurationPanel)
        modeGroup = QGroupBox('Analyses', configurationPanel)
        modeLayout = QVBoxLayout(modeGroup)
        self.hairpin = QCheckBox('Hairpin', modeGroup)
        self.homodimer = QCheckBox('Homodimer', modeGroup)
        self.heterodimer = QCheckBox(
            'Heterodimers (all pairwise combinations of selected sequences)',
            modeGroup)
        for checkbox in (self.hairpin, self.homodimer, self.heterodimer):
            checkbox.setChecked(True)
            modeLayout.addWidget(checkbox)
        configurationLayout.addWidget(modeGroup)

        fixedSettings = QLabel(
            'Fixed conditions: monovalent cations, 50 mM; divalent cations, '
            '10 mM; dNTPs, 0 mM; DNA, 100 nM; ΔG temperature, 37 °C. '
            'Primer3 internal structure-search parameters are not displayed '
            'as screening criteria.',
            configurationPanel)
        fixedSettings.setWordWrap(True)
        configurationLayout.addWidget(fixedSettings)

        self.runButton = QPushButton(
            'Run Primer3 Analysis', configurationPanel)
        self.exportButton = QPushButton(
            'Export Analysis XLSX…', configurationPanel)
        self.exportButton.setEnabled(False)
        self.summary = QLabel('', configurationPanel)
        self.summary.setWordWrap(True)
        configurationLayout.addWidget(self.runButton)
        configurationLayout.addWidget(self.exportButton)
        configurationLayout.addWidget(self.summary)
        self.setupSplitter.addWidget(configurationPanel)
        self.setupSplitter.setSizes((460, 300))
        self.mainSplitter.addWidget(sequencePanel)

        resultPanel = QWidget(self.mainSplitter)
        resultLayout = QVBoxLayout(resultPanel)
        resultTitle = QLabel('<b>Analysis Results</b>', resultPanel)
        resultLayout.addWidget(resultTitle)
        resultDescription = QLabel(
            'Results are sorted by ascending ΔG. More-negative values '
            'indicate more stable predicted structures and therefore higher '
            'interaction risk.', resultPanel)
        resultDescription.setWordWrap(True)
        resultLayout.addWidget(resultDescription)

        self.resultSplitter = QSplitter(Qt.Orientation.Vertical, resultPanel)
        self.resultSplitter.setChildrenCollapsible(False)
        resultLayout.addWidget(self.resultSplitter, 1)

        resultTablePanel = QGroupBox('Thermodynamic Results',
                                     self.resultSplitter)
        resultTableLayout = QVBoxLayout(resultTablePanel)
        self.resultTable = QTableWidget(resultTablePanel)
        self.resultTable.setColumnCount(9)
        self.resultTable.setHorizontalHeaderLabels((
            'Analysis Type', 'Sequence 1', 'Sequence 2', 'Tm (°C)',
            'ΔG (kcal/mol)', 'ΔH (kcal/mol)', 'ΔS (cal/(K·mol))',
            'Structure', 'Notes'))
        self.resultTable.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self.resultTable.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.resultTable.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection)
        resultHeader = self.resultTable.horizontalHeader()
        resultHeader.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        resultHeader.setMinimumSectionSize(68)
        for column, width in enumerate(
                (135, 145, 145, 90, 125, 125, 145, 105, 220)):
            self.resultTable.setColumnWidth(column, width)
        self.resultTable.verticalHeader().setSectionResizeMode(
            QHeaderView.ResizeMode.Interactive)
        resultTableLayout.addWidget(self.resultTable)
        self.resultSplitter.addWidget(resultTablePanel)

        structurePanel = QGroupBox('Predicted Pairing Structure',
                                   self.resultSplitter)
        structureLayout = QVBoxLayout(structurePanel)
        structureLabel = QLabel(
            'Select a result above to inspect its structure. Use the mouse '
            'wheel to zoom and drag to pan.', structurePanel)
        structureLabel.setWordWrap(True)
        structureLayout.addWidget(structureLabel)
        structureTabs = QTabWidget(structurePanel)
        self.structureGraphic = Primer3StructureView(structureTabs)
        structureTabs.addTab(self.structureGraphic, '2D Pairing Diagram')
        self.structure = QPlainTextEdit(structurePanel)
        self.structure.setReadOnly(True)
        structureFont = QFont('Monaco')
        structureFont.setStyleHint(QFont.StyleHint.Monospace)
        structureFont.setPointSize(11)
        self.structure.setFont(structureFont)
        structureTabs.addTab(self.structure, 'Primer3 Raw Structure')
        structureLayout.addWidget(structureTabs)
        self.resultSplitter.addWidget(structurePanel)
        self.resultSplitter.setSizes((340, 400))
        self.mainSplitter.addWidget(resultPanel)
        self.mainSplitter.setStretchFactor(0, 2)
        self.mainSplitter.setStretchFactor(1, 3)
        self.mainSplitter.setSizes((560, 840))

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close,
                                   parent=self)
        buttons.rejected.connect(self.reject)
        outer.addWidget(buttons)

        self.importButton.clicked.connect(self._importFiles)
        self.selectAllButton.clicked.connect(lambda: self._checkAll(True))
        self.selectNoneButton.clicked.connect(lambda: self._checkAll(False))
        self.removeButton.clicked.connect(self._removeSelectedRows)
        self.runButton.clicked.connect(self._run)
        self.exportButton.clicked.connect(self._export)
        self.resultTable.currentCellChanged.connect(self._showStructure)
        self._populateSequences()

    def _populateSequences(self):
        self.sequenceTable.setRowCount(len(self._entries))
        colors = {
            'Input': QColor('#eaf3ff'),
            'Newly generated': QColor('#ecf8ee'),
            'Imported': QColor('#fff8e6'),
            # Backward-compatible colors for projects saved before the
            # English-only interface was introduced.
            '输入': QColor('#eaf3ff'), '新生成': QColor('#ecf8ee'),
            '导入': QColor('#fff8e6')}
        for row, entry in enumerate(self._entries):
            check = QTableWidgetItem('')
            check.setFlags(check.flags() |
                           Qt.ItemFlag.ItemIsUserCheckable)
            check.setCheckState(Qt.CheckState.Checked)
            values = (entry['source'], entry['name'], entry['sequence'],
                      len(entry['sequence']))
            self.sequenceTable.setItem(row, 0, check)
            background = colors.get(entry['source'], QColor('#fff8e6'))
            for column, value in enumerate(values, 1):
                item = QTableWidgetItem(str(value))
                item.setBackground(QBrush(background))
                self.sequenceTable.setItem(row, column, item)
        self.summary.setText('%d short sequences loaded.' % len(self._entries))

    def _checkAll(self, checked):
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        for row in range(self.sequenceTable.rowCount()):
            self.sequenceTable.item(row, 0).setCheckState(state)

    def _removeSelectedRows(self):
        selected = sorted(set(index.row()
                              for index in self.sequenceTable.selectedIndexes()),
                          reverse=True)
        for row in selected:
            del self._entries[row]
        if selected:
            self._populateSequences()

    def _importFiles(self):
        filenames = QFileDialog.getOpenFileNames(
            self, 'Import Sequences for Primer3 Analysis',
            self._suggestedDirectory,
            'Sequence files (*.txt *.xlsx);;Text files (*.txt);;'
            'Orthogonal sequence workbooks (*.xlsx)')
        if isinstance(filenames, tuple):
            filenames = filenames[0]
        if not filenames:
            return
        imported = []
        messages = []
        for filename in filenames:
            try:
                entries, errors = read_primer3_sequences(filename)
                imported.extend(entries)
                messages.extend('%s: %s' % (os.path.basename(filename), error)
                                for error in errors)
            except (IOError, OSError, UnicodeError, ValueError) as error:
                messages.append('%s: %s' % (os.path.basename(filename), error))
        before = len(self._entries)
        self._entries = normalized_entries(self._entries + imported)
        self._populateSequences()
        self._suggestedDirectory = os.path.dirname(filenames[0])
        if messages:
            message = QMessageBox(
                QMessageBox.Icon.Warning, 'Primer3 Sequence Import',
                '%d sequences were added, but some entries could not be '
                'imported.' %
                (len(self._entries) - before),
                QMessageBox.StandardButton.Ok, self)
            message.setDetailedText('\n'.join(messages))
            message.exec()

    def _selectedEntries(self):
        return [entry for row, entry in enumerate(self._entries)
                if self.sequenceTable.item(row, 0).checkState() ==
                Qt.CheckState.Checked]

    def _selectedModes(self):
        modes = []
        if self.hairpin.isChecked():
            modes.append('hairpin')
        if self.homodimer.isChecked():
            modes.append('homodimer')
        if self.heterodimer.isChecked():
            modes.append('heterodimer')
        return modes

    def _run(self):
        entries = self._selectedEntries()
        modes = self._selectedModes()
        if not entries:
            QMessageBox.warning(
                self, 'Primer3 Analysis', 'Select at least one sequence.')
            return
        if not modes:
            QMessageBox.warning(
                self, 'Primer3 Analysis', 'Select at least one analysis.')
            return
        if modes == ['heterodimer'] and len(entries) < 2:
            QMessageBox.warning(
                self, 'Primer3 Analysis',
                'Heterodimer analysis requires at least two sequences.')
            return
        total = ((len(entries) if 'hairpin' in modes else 0) +
                 (len(entries) if 'homodimer' in modes else 0) +
                 (len(entries) * (len(entries) - 1) // 2
                  if 'heterodimer' in modes else 0))
        progress = QProgressDialog(
            'Initializing Primer3…', 'Cancel', 0, total, self)
        progress.setWindowTitle('Primer3 Thermodynamic Analysis')
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setMinimumDuration(0)

        def update(completed, unused_total, label):
            progress.setValue(completed)
            progress.setLabelText(
                'Completed %d of %d\n%s' % (completed, total, label))
            QApplication.processEvents()

        try:
            self._results = run_primer3_analysis(
                entries, modes=modes, progress=update,
                cancelled=progress.wasCanceled)
        except Primer3Cancelled:
            self.summary.setText('Primer3 analysis was canceled.')
            return
        except Primer3Unavailable as error:
            QMessageBox.critical(self, 'Primer3 Unavailable', str(error))
            return
        except Exception as error:
            QMessageBox.critical(self, 'Primer3 Analysis Failed', str(error))
            return
        finally:
            progress.close()
        self._populateResults()

    @staticmethod
    def _displayNumber(value):
        return '' if value == '' else ('%.3f' % float(value)).rstrip('0').rstrip('.')

    def _populateResults(self):
        self.resultTable.setRowCount(len(self._results))
        danger = QColor('#ffe9e9')
        for row, result in enumerate(self._results):
            second = result['second']
            values = (
                result['kind_label'], result['first']['name'],
                second['name'] if second else '',
                self._displayNumber(result['tm']),
                self._displayNumber(result['dg']),
                self._displayNumber(result['dh']),
                self._displayNumber(result['ds']),
                'Found' if result['structure_found'] else 'Not found',
                result['error'])
            for column, value in enumerate(values):
                item = QTableWidgetItem(str(value))
                item.setData(Qt.ItemDataRole.UserRole, row)
                if row < min(10, len(self._results)) and result['dg'] != '':
                    item.setBackground(QBrush(danger))
                self.resultTable.setItem(row, column, item)
        failed = sum(bool(result['error']) for result in self._results)
        self.summary.setText(
            '%d analyses completed; %d could not be completed because of '
            'sequence-length or calculation errors. Results are sorted by '
            'ascending ΔG (highest predicted risk first).' %
            (len(self._results), failed))
        self.exportButton.setEnabled(bool(self._results))
        if self._results:
            self.resultTable.selectRow(0)
            self._showStructure(0, 0, -1, -1)

    def _showStructure(self, row, unused_column, unused_old_row,
                       unused_old_column):
        if row < 0 or row >= len(self._results):
            self.structure.clear()
            self.structureGraphic.showResult(None)
            return
        result = self._results[row]
        self.structureGraphic.showResult(result)
        text = result['structure']
        if not text:
            text = result['error'] or \
                'Primer3 did not return a displayable pairing structure.'
        self.structure.setPlainText(text)

    def _export(self):
        if not self._results:
            return
        filename = QFileDialog.getSaveFileName(
            self, 'Export Primer3 Analysis Results',
            os.path.join(self._suggestedDirectory, 'primer3_analysis.xlsx'),
            'Excel workbook (*.xlsx)')
        if isinstance(filename, (tuple, list)):
            filename = filename[0]
        if not filename:
            return
        if not filename.lower().endswith('.xlsx'):
            filename += '.xlsx'
        try:
            write_primer3_workbook(filename, self._results)
        except (IOError, OSError, ValueError) as error:
            QMessageBox.critical(self, 'Primer3 Export Failed', str(error))
            return
        self._suggestedDirectory = os.path.dirname(filename)
        QMessageBox.information(
            self, 'Primer3 Analysis', 'Results saved to:\n%s' % filename)
