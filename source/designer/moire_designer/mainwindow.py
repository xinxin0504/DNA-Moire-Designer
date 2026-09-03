"""Standalone and embeddable DNA Moiré bilayer designer window."""

from __future__ import annotations

import math
import os
import csv
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path

from PyQt6.QtCore import (QPointF, QRect, QRectF, QSize, QSignalBlocker,
                          QTimer, Qt, pyqtSignal)
from PyQt6.QtGui import (
    QAction, QColor, QCursor, QFont, QIcon, QImage, QPainter,
    QPen, QPixmap, QPolygonF)
from PyQt6.QtWidgets import (
    QAbstractSpinBox,
    QAbstractItemView,
    QApplication,
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QLayout,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QProgressDialog,
    QScrollArea,
    QSizePolicy,
    QSlider,
    QSpinBox,
    QSplitter,
    QStackedWidget,
    QStatusBar,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextBrowser,
    QToolBar,
    QVBoxLayout,
    QWidget,
)
from PyQt6.QtSvg import QSvgGenerator
from cadnano2.model.io.orthogonalseq import (
    longest_common_substring, reverse_complement)
from cadnano2.model.io.sequencexlsx import write_sequence_template

from moire_design_core import (
    MAX_SEED_DELETIONS_PER_DOMAIN,
    MAX_SEED_INSERTION_PER_HELIX,
    MoireProject,
    SquareBilayerSettings,
    analyze_sequence_design,
    assign_standard_scaffold_sequence,
    build_sequenced_design,
    compatible_growth_values,
    compatible_z2_values,
    export_final_package,
    export_sequence_variants,
    export_sst_input_template,
    estimate_scaffold_capacity,
    fixed_seed_overlap_layout,
    extract_scaffold_sequence,
    finalize_structure,
    generate_scaffold_review,
    import_sst_input_template,
    list_standard_scaffolds,
    load_project,
    maximum_seed_insertion_per_helix,
    minimum_seed_deletion_per_helix,
    save_project,
    solve_square_bilayer,
    structure_layout,
    validate_sst,
    validate_structure,
    write_shifted_sst,
)
from moire_design_core.project import add_measurement, export_capture_map
from moire_design_core.template import export_reference_seed, reference_seed_path
from moire_runtime import cadnano_executable, worker_command

from .preview import (BilayerPreview, MoireTopViewPreview,
                      SeedCrossSectionPicker,
                      StructureDesignPreview)
from .analysis_bulk import AnalysisBulkMixin
from .orthogonal_sequence_tool import (
    generate_orthogonal_sequences_automatic,
    run_orthogonal_sequence_designer,
)
from .project_session import ProjectSetupDialog
from .i18n import (UiLocalizer, current_language, install_dialog_hooks,
                   install_painter_hook, localize_svg, localize_xlsx, set_language,
                   translate)


APP_ROOT = Path(__file__).resolve().parent
CADNANO_EXECUTABLE = cadnano_executable()


class LocalizedStatusBar(QStatusBar):
    """Translate dynamic status messages at the point they are displayed."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._source_message = ""
        self._source_timeout = 0
        self.setSizeGripEnabled(False)
        # A blank QStatusBar still reserves roughly 22 px below every page.
        # Collapse it whenever no message is present so the page's own 8-px
        # lower inset is the only visible whitespace. Timed status messages
        # temporarily reopen the bar and collapse it again when they expire.
        self.messageChanged.connect(self._sync_visibility)
        self.hide()

    def _sync_visibility(self, message):
        self.setVisible(bool(str(message).strip()))

    def showMessage(self, message, timeout=0):
        self._source_message = str(message)
        self._source_timeout = int(timeout)
        if self._source_message.strip():
            self.show()
        super().showMessage(translate(self._source_message), timeout)

    def retranslate(self):
        if self.currentMessage() and self._source_message:
            super().showMessage(
                translate(self._source_message), self._source_timeout)


class ResultCard(QFrame):
    def __init__(self, title, parent=None):
        super().__init__(parent)
        self.setObjectName("resultCard")
        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 10, 14, 10)
        layout.setSpacing(2)
        caption = QLabel(title)
        caption.setObjectName("cardCaption")
        self.value = QLabel("—")
        self.value.setObjectName("cardValue")
        self.value.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(caption)
        layout.addWidget(self.value)


def _preview_parameter_html(twist_text="—", period_text="—"):
    """Compact white metrics shown together under the Top-view heading."""
    return (
        '<table width="100%%" cellspacing="0" cellpadding="0">'
        '<tr><td width="58%%" align="center">'
        '<span style="color:#f3f6f9;font-size:11px;font-weight:700;">'
        'Twist angle: </span>'
        '<span style="color:#ffffff;font-size:11px;font-weight:700;">%s</span>'
        '</td><td width="42%%" align="center">'
        '<span style="color:#f3f6f9;font-size:11px;font-weight:700;">'
        'Moiré period: </span>'
        '<span style="color:#ffffff;font-size:11px;font-weight:700;">%s</span>'
        '</td></tr></table>' % (twist_text, period_text))


def _side_preview_parameter_html(
        sst_z1="—", spacing="—", sst_z3="—",
        seed_z1="—", seed_z2="—", seed_z3="—"):
    """Two compact dimensional rows shown under the Side-view heading."""
    return (
        '<table width="90%%" align="center" cellspacing="0" cellpadding="0" '
        'style="color:#f3f6f9;font-size:11px;font-weight:650;">'
        '<tr>'
        '<td width="20%%" align="left">SST sublattice</td>'
        '<td width="26%%" align="right" style="color:#2a78d1;">'
        '1st layer:&nbsp;%s</td>'
        '<td width="27%%" align="right" style="color:#8a61bb;">'
        'Spacing:&nbsp;%s</td>'
        '<td width="27%%" align="right" style="color:#d65b74;">'
        '2nd layer:&nbsp;%s</td>'
        '</tr>'
        '<tr><td height="11" colspan="4"></td></tr>'
        '<tr>'
        '<td width="20%%" align="left">Seed</td>'
        '<td width="26%%" align="right" style="color:#2a78d1;">'
        'Z1:&nbsp;%s</td>'
        '<td width="27%%" align="right" style="color:#8a61bb;">'
        'Z2:&nbsp;%s</td>'
        '<td width="27%%" align="right" style="color:#d65b74;">'
        'Z3:&nbsp;%s</td>'
        '</tr></table>' % (
            sst_z1, spacing, sst_z3, seed_z1, seed_z2, seed_z3))


def _format_bp(value):
    """Show physical bp lengths without hiding fractional mean indels."""
    if value is None:
        return "—"
    number = float(value)
    if abs(number - round(number)) < 1e-9:
        return str(int(round(number)))
    return "%.1f" % number


def _twist_handedness(angle_deg):
    angle_deg = float(angle_deg)
    if angle_deg > 0.0:
        return "right-handed"
    if angle_deg < 0.0:
        return "left-handed"
    return "no handedness"


def _sst_lattice_for_symmetry(symmetry):
    """Map a bilayer UI selection to its stage-2 SST generator."""
    symmetry = str(symmetry)
    if symmetry == "square_kagome":
        return "square_kagome"
    if symmetry == "kagome_kagome":
        return "kagome"
    return "square"


class SignedAngleSpinBox(QDoubleSpinBox):
    """Display an explicit sign while retaining normal numeric editing."""

    def textFromValue(self, value):
        precision = self.decimals()
        if abs(float(value)) < 0.5 * (10.0 ** -precision):
            value = 0.0
        return ("%%+.%df" % precision) % float(value)


class PhaseLengthCombo(QComboBox):
    """A phase-filtered list that can become the next phase anchor."""

    aboutToOpen = pyqtSignal()

    def showPopup(self):
        self.aboutToOpen.emit()
        super().showPopup()


class MoireDesignerWindow(AnalysisBulkMixin, QMainWindow):
    """Design/prediction UI with a calibrated selectable Square Seed."""

    def __init__(self, parent=None, cadnano_controller=None):
        super().__init__(parent)
        # The distributable application is English-only.
        set_language("en")
        self.setStatusBar(LocalizedStatusBar(self))
        install_painter_hook()
        install_dialog_hooks()
        self.cadnano_controller = cadnano_controller
        self.project = None
        self.project_path = None
        self._updating = False
        self._phase_sync = False
        self._target_driver = "indel"
        self._orthogonal_primer3_entries = []
        self._sequence_analysis = None
        self._sequence_assignments = {}
        self._auto_input_design_running = False
        self.structure_root = None
        self._app_mode = None
        self._startup_complete = False
        self._last_design_step = 0
        self._last_analysis_module = 0
        self._history = []
        self._history_index = -1
        self._restoring_history = False
        self._history_restore_token = 0
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.setWindowTitle("DNA Moiré Designer")
        screen = QApplication.primaryScreen()
        if screen is not None:
            available = screen.availableGeometry()
            width = min(1420, max(900, available.width()-80))
            height = min(880, max(600, available.height()-70))
            self.resize(width, height)
            self.setMinimumSize(min(900, width), min(600, height))
        else:
            self.resize(1280, 720)
        icon = APP_ROOT / "assets" / "moire-design.svg"
        if icon.is_file():
            self.setWindowIcon(QIcon(str(icon)))
        self._build_ui()
        self._localizer = UiLocalizer(self)
        self._apply_style()
        self._connect_signals()
        self.apply_paper_preset()
        # Startup is intentionally non-modal: the application opens directly
        # in Design mode.  Analysis is an in-window workspace selected from
        # the workflow bar and never requires a project file.
        self._apply_app_mode("design")
        # Apply the complete English catalog before the first frame is shown.
        # The timer remains active for labels/reports created after clicks.
        self._localizer.retranslate()
        self._startup_complete = True
        self._autosave_timer = QTimer(self)
        self._autosave_timer.setInterval(30000)
        self._autosave_timer.timeout.connect(self._autosave_project)
        self._autosave_timer.start()

    def _reset_downstream_design_ui(self):
        """Immediately show Steps 2 and 3 as pending after a design edit."""
        if hasattr(self, "generate_simple_design_button"):
            self.generate_simple_design_button.setEnabled(False)
            self.inspect_final_design_button.setEnabled(False)
            self.structure_next_button.setEnabled(False)
            self._set_acceptance_button(
                self.accept_structure_button, False,
                "Accept current DNA design")
            self.accept_structure_button.setEnabled(False)
            self.accepted_parameters_summary.setText(
                "No accepted Moiré parameters")
            self.accepted_design_summary.setText("No accepted DNA design")
            self.design_generation_action_status.hide()
            self.cadnano_edit_action_status.hide()
            self.structure_accept_action_status.hide()
        self._reset_sequence_assignment_ui()

    def _reset_sequence_assignment_ui(self):
        """Restore sequence assignment controls to their initial state."""
        if not hasattr(self, "accept_added_scaffold_button"):
            return
        self._set_acceptance_button(
            self.accept_added_scaffold_button, False,
            "Accept assigned scaffold sequences")
        self.accept_added_scaffold_button.setEnabled(False)
        self._set_acceptance_button(
            self.accept_added_sst_button, False,
            "Accept assigned SST sublattice input sequences")
        self.accept_added_sst_button.setEnabled(False)
        self.detect_scaffold_sequences_button.setEnabled(False)
        self.detect_sst_inputs_button.setEnabled(False)
        self.auto_design_sst_inputs_button.setEnabled(False)
        self.export_sst_input_template_button.setEnabled(False)
        self.import_sst_input_template_button.setEnabled(False)
        self.final_sequence_export_button.setEnabled(False)
        self._clear_sequence_cards(self.scaffold_cards_layout)
        self.scaffold_cards_layout.addWidget(QLabel(
            "The structure has not been read."))
        self._clear_sequence_cards(self.sst_cards_layout)
        self.sst_cards_layout.addWidget(QLabel(
            "SST sublattice input positions have not been detected."))
        self.scaffold_detection_action_status.hide()
        self.scaffold_sequence_status.hide()
        self.sst_detection_status.setText(
            "Accept the assigned scaffold sequences first.")
        for label in (self.sst_auto_import_status,
                      self.sst_expert_import_status,
                      self.sst_acceptance_status):
            label.clear()
            label.hide()
        self.sequence_export_status.setText(
            "Sequence assignment is not complete.")
        self.sequence_preview.clear()
        self.sequence_preview_status.setText(
            "Detect an accepted design to display its sequence routes.")

    def _invalidate_final_design_acceptance(self):
        """Require final-design acceptance again and reset Step 3."""
        workflow = self._workflow()
        self._drop_keys(workflow, (
            "structure_accepted_at", "sequence_analysis",
            "sequence_assignments", "sequence_design_json",
            "sequence_scaffold_accepted",
            "sequence_scaffold_accepted_at", "sequence_sst_accepted",
            "sequence_sst_accepted_at", "sequence_sst_detected",
            "sequence_sst_detection_status",
            "sequence_sst_import_method", "sequence_sst_import_status",
            "sequence_sst_imported_at", "sequence_sst_import_source",
            "sequence_sst_acceptance_status", "sequence_source",
            "sequence_exports"))
        workflow["structure_accepted"] = False
        self._sequence_analysis = None
        self._sequence_assignments = {}
        self._reset_sequence_assignment_ui()
        self.structure_next_button.setEnabled(False)
        self._set_acceptance_button(
            self.accept_structure_button, False,
            "Accept current DNA design")
        self.accept_structure_button.setEnabled(bool(
            workflow.get("structure_complete")))
        self.structure_accept_action_status.hide()

    def _build_ui(self):
        self.new_action = QAction("新建项目…", self)
        self.open_action = QAction("打开 .moire.json", self)
        self.save_action = QAction("保存", self)
        self.save_as_action = QAction("另存为…", self)
        self.export_action = QAction("导出原型项目", self)
        self.cadnano_action = QAction("在 cadnano 中打开当前结构", self)

        menu_bar = self.menuBar()
        menu_bar.setNativeMenuBar(True)
        file_menu = menu_bar.addMenu("文件")
        file_menu.addActions([self.new_action, self.open_action,
                              self.save_action, self.save_as_action,
                              self.export_action])
        self.paper_preset_menu_action = QAction("应用论文参数预设", self)
        self.view_design_action = QAction("设计与预测", self)
        self.view_capture_action = QAction("Automated DNA Design", self)
        self.view_sequence_action = QAction("序列与导出", self)
        self.analysis_crystal_action = QAction(
            "Moiré analysis", self)

        root = QWidget(self)
        root_layout = QVBoxLayout(root)
        # Each workflow page owns its 8-px lower inset. Do not stack another
        # root-level lower margin beneath it.
        root_layout.setContentsMargins(10, 8, 10, 0)
        root_layout.setSpacing(6)
        workflow = QFrame()
        workflow.setObjectName("workflowBar")
        workflow_layout = QHBoxLayout(workflow)
        workflow_layout.setContentsMargins(5, 4, 5, 4)
        workflow_layout.setSpacing(4)
        self.workflow_group = QButtonGroup(self)
        self.workflow_group.setExclusive(True)
        self.workflow_buttons = []
        self.workflow_design_separators = []
        labels = (
            "1 Moiré 参数输入",
            "1.2 层长度与层间距",
            "2 Automated DNA Design", "3  序列导出",
            "Moiré analysis")
        object_names = (
            "workflowButton", "workflowButton", "workflowButton",
            "workflowButton", "workflowButton")
        for index, label in enumerate(labels):
            button = QPushButton(label)
            button.setObjectName(object_names[index])
            button.setCheckable(True)
            if index < 4:
                # The three visible Design workflow stages use one identical
                # frame.  The retired hidden 1.2 compatibility slot keeps the
                # same geometry but never appears.
                button.setFixedSize(190, 32)
            else:
                button.setFixedHeight(self.ANALYSIS_BUTTON_HEIGHT)
            button.setSizePolicy(QSizePolicy.Policy.Fixed,
                                 QSizePolicy.Policy.Fixed)
            if index < 4:
                button.clicked.connect(
                    lambda checked=False, step=index: self._go_to_step(step))
            else:
                button.clicked.connect(
                    lambda checked=False, module=index-4:
                    self._open_analysis_module(module))
            self.workflow_group.addButton(button, index)
            self.workflow_buttons.append(button)
            workflow_layout.addWidget(button)
            if index <= 3:
                arrow = QLabel("│" if index == 3 else "›")
                arrow.setObjectName("workflowArrow")
                workflow_layout.addWidget(arrow)
                self.workflow_design_separators.append(arrow)
            if index == 3:
                # Keep the history implementation available internally for
                # compatibility with saved sessions, but remove the visible
                # Undo/Redo controls from the workflow bar.
                self.history_back_button = QPushButton("←", workflow)
                self.history_back_button.setObjectName("historyButton")
                self.history_back_button.setToolTip("撤回上一个设计操作")
                self.history_back_button.setEnabled(False)
                self.history_back_button.hide()
                self.history_forward_button = QPushButton("→", workflow)
                self.history_forward_button.setObjectName("historyButton")
                self.history_forward_button.setToolTip("恢复下一个设计操作")
                self.history_forward_button.setEnabled(False)
                self.history_forward_button.hide()
        workflow_layout.addStretch(1)
        self.mode_switch_button = QPushButton("Switch to Analysis Mode")
        self.mode_switch_button.setObjectName("primaryButton")
        self.mode_switch_button.setFixedHeight(self.ANALYSIS_BUTTON_HEIGHT)
        self.mode_switch_button.setMinimumWidth(190)
        self.mode_switch_button.setSizePolicy(
            QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)
        workflow_layout.addWidget(self.mode_switch_button)
        # Step 1 is now a single combined parameter-and-preview page.  Keep
        # the former 1.2 button alive only as a compatibility slot for old
        # saved history indices; it is never shown or navigated to.
        self.workflow_buttons[1].hide()
        self.workflow_design_separators[0].hide()
        self.workflow_buttons[0].setChecked(True)
        root_layout.addWidget(workflow)

        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_design_tab(), "设计与预测")
        self.tabs.addTab(self._build_capture_tab(), "Automated DNA Design")
        self.tabs.addTab(self._build_sequence_tab(), "序列与导出")
        self.tabs.addTab(self._build_analysis_tab(), "Analysis modules")
        self.tabs.tabBar().hide()
        self.tabs.currentChanged.connect(self._tab_changed)
        root_layout.addWidget(self.tabs, 1)
        self.setCentralWidget(root)
        self.statusBar().showMessage(
            "Square Seed截面校准节点已载入", 5000)

    def _scroll_panel(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        content = QWidget()
        layout = QVBoxLayout(content)
        layout.setContentsMargins(2, 2, 8, 8)
        layout.setSpacing(12)
        scroll.setWidget(content)
        return scroll, layout

    @staticmethod
    def _action_feedback_label():
        """Create a compact, hidden success/next-step message."""
        label = QLabel()
        label.setWordWrap(True)
        label.setObjectName("successStatus")
        policy = QSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        policy.setHeightForWidth(True)
        label.setSizePolicy(policy)
        label.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        label.hide()
        return label

    @staticmethod
    def _show_action_feedback(label, message):
        label.setObjectName("successStatus")
        # Some feedback targets predate the dedicated helper. Normalize their
        # policy here as well: Maximum prevents a QLabel from absorbing spare
        # vertical space, while height-for-width lets wrapped text request all
        # lines it actually needs. Avoid fixed heights so nothing is clipped.
        policy = QSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        policy.setHeightForWidth(True)
        label.setSizePolicy(policy)
        label.setMinimumHeight(0)
        label.setMaximumHeight(16777215)
        label.setWordWrap(True)
        label.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        label.setText(translate(str(message)))
        label.style().unpolish(label)
        label.style().polish(label)
        label.show()
        label.updateGeometry()
        # QLabel's wrapped height can be underestimated until the containing
        # scroll panel has completed its layout.  Reserve the measured text
        # height immediately and once more on the next event-loop pass so a
        # multi-line result card is never vertically clipped.
        def fit_height():
            if not label.isVisible():
                return
            width = max(160, label.contentsRect().width())
            required = label.heightForWidth(width)
            if required > 0:
                label.setMinimumHeight(required)
                label.updateGeometry()
        fit_height()
        QTimer.singleShot(0, fit_height)

    def _build_design_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        # Keep the preview frame equally inset from the top and bottom.
        layout.setContentsMargins(0, 8, 0, 8)
        self.design_stack = QStackedWidget()
        selection_page = self._build_design_selection_page()
        parameter_page = self._build_design_parameter_page()
        self._merge_design_parameter_pages()
        self.design_stack.addWidget(selection_page)
        # Retain an inert compatibility page so a history snapshot written by
        # an older release can still be restored safely. Navigation redirects
        # its old step number back to the combined page.
        self.design_stack.addWidget(parameter_page)
        self.design_stack.setCurrentIndex(0)
        layout.addWidget(self.design_stack, 1)
        return tab

    def _build_design_selection_page(self):
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(0, 0, 0, 0)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        layout.addWidget(splitter, 1)

        left_panel = QWidget()
        left = QVBoxLayout(left_panel)
        self.design_basis_left_layout = left
        # The scroll area, rather than a forced minimum-size layout, owns the
        # viewport height.  SetMinimumSize caused the macOS scroll thumb and
        # splitter drag to jump while accordion sections were opened/closed.
        left.setSizeConstraint(QLayout.SizeConstraint.SetDefaultConstraint)
        left.setContentsMargins(4, 0, 10, 6)
        left.setSpacing(9)
        left_panel.setMinimumWidth(455)
        left_panel.setMaximumWidth(535)

        self.design_basis_section_group = QButtonGroup(self)
        self.design_basis_section_group.setExclusive(True)
        self.symmetry_step_button = QPushButton("1.1 选择双层对称性")
        self.symmetry_step_button.setObjectName("parameterStepButton")
        self.symmetry_step_button.setCheckable(True)
        self.symmetry_step_button.setMinimumHeight(34)
        self.design_basis_section_group.addButton(
            self.symmetry_step_button, 0)
        left.addWidget(self.symmetry_step_button)

        symmetry_box = QGroupBox("双层对称性")
        self.symmetry_box = symmetry_box
        symmetry_box.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        symmetry_layout = QVBoxLayout(symmetry_box)
        symmetry_layout.setContentsMargins(9, 9, 9, 8)
        self.bilayer_symmetry_selector = QComboBox()
        self.bilayer_symmetry_selector.addItem(
            "Square–Square", "square_square_c4")
        self.bilayer_symmetry_selector.addItem(
            "Kagome–Kagome", "kagome_kagome")
        self.bilayer_symmetry_selector.addItem(
            "Square–Kagome", "square_kagome")
        symmetry_layout.addWidget(self.bilayer_symmetry_selector)
        self.bilayer_symmetry_note = QLabel(symmetry_box)
        self.bilayer_symmetry_note.hide()
        left.addWidget(symmetry_box)

        self.twist_period_step_button = QPushButton(
            "1.2 输入 Twist 或 Moiré period")
        self.twist_period_step_button.setObjectName("parameterStepButton")
        self.twist_period_step_button.setCheckable(True)
        self.twist_period_step_button.setMinimumHeight(34)
        self.design_basis_section_group.addButton(
            self.twist_period_step_button, 1)
        left.addWidget(self.twist_period_step_button)

        target_box = QGroupBox("Twist / Moiré period")
        self.target_box = target_box
        target_box.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        target_layout = QGridLayout(target_box)
        target_layout.setContentsMargins(9, 9, 9, 8)
        target_layout.setHorizontalSpacing(7)
        target_layout.setVerticalSpacing(5)
        target_layout.setColumnMinimumWidth(0, 118)
        target_layout.setColumnStretch(1, 1)
        self.project_name = QLineEdit("square_moire_bilayer", page)
        self.project_name.hide()
        self.target_definition = QComboBox()
        self.target_definition.addItem(
            "局部生长面角度", "local_surface")
        self.target_definition.addItem(
            "实验总角度（自动计入capture）", "experimental_total")
        self.angle = SignedAngleSpinBox()
        # Twist is intentionally not constrained by an arbitrary UI range.
        # Feasibility is checked later from the required insertion/deletion
        # load, where the user receives the structure-specific warning.
        self.angle.setRange(-1.0e9, 1.0e9)
        self.angle.setDecimals(1)
        self.angle.setSuffix("° (no handedness)")
        self.period = QDoubleSpinBox()
        self.period.setRange(0.0, 5000.0)
        self.period.setDecimals(1)
        self.period.setSuffix(" nm")
        self.period.setSpecialValueText("∞")
        self.period_label = QLabel("Moiré period")
        self.lattice_symmetry = QLabel("Square–Square")
        self.lattice_symmetry.setObjectName("fixedValue")
        self.lattice_context = QComboBox()
        self.lattice_context.addItem(
            "溶液 / cryo-EM", ("solution_cryo", 2.8, 5.4))
        self.lattice_context.addItem(
            "干燥 TEM", ("dried_tem", 2.2, 4.4))
        self.lattice_constant = QDoubleSpinBox()
        self.lattice_constant.setRange(.1, 20.0)
        self.lattice_constant.setDecimals(1)
        self.lattice_constant.setSuffix(" nm")
        self.lattice_constant_2 = QDoubleSpinBox()
        self.lattice_constant_2.setRange(.1, 20.0)
        self.lattice_constant_2.setDecimals(1)
        self.lattice_constant_2.setSuffix(" nm")
        self.lattice_constant_label = QLabel("1st layer a (square)")
        self.lattice_constant_2_label = QLabel("2nd layer a (square)")
        self.lattice_constant_fixed = QLabel()
        self.lattice_constant_fixed.setObjectName("fixedValue")
        self.lattice_constant_2_fixed = QLabel()
        self.lattice_constant_2_fixed.setObjectName("fixedValue")
        twist_label = QLabel("Twist")
        twist_label.setStyleSheet("font-weight:700")
        self.period_label.setStyleSheet("font-weight:700")
        target_layout.addWidget(twist_label, 0, 0)
        target_layout.addWidget(self.angle, 0, 1, 1, 2)
        target_layout.addWidget(self.period_label, 1, 0)
        target_layout.addWidget(self.period, 1, 1, 1, 2)
        target_layout.addWidget(QLabel("测量环境"), 2, 0)
        target_layout.addWidget(self.lattice_context, 2, 1, 1, 2)
        target_layout.addWidget(self.lattice_constant_label, 3, 0)
        target_layout.addWidget(self.lattice_constant_fixed, 3, 1, 1, 2)
        target_layout.addWidget(self.lattice_constant_2_label, 4, 0)
        target_layout.addWidget(self.lattice_constant_2_fixed, 4, 1, 1, 2)
        self.target_definition.hide()
        self.lattice_constant.hide()
        self.lattice_constant_2.hide()
        left.addWidget(target_box)

        seed_box = QGroupBox("Seed 截面")
        self.seed_cross_section_box = seed_box
        seed_box.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        seed_layout = QVBoxLayout(seed_box)
        seed_layout.setContentsMargins(9, 8, 9, 8)
        seed_layout.setSpacing(5)
        self.seed_cross_section_preset = QComboBox()
        self.seed_cross_section_preset.setObjectName("editableParameter")
        self.seed_cross_section_preset.hide()
        self.seed_cross_section_preset_display = QLabel(
            "8×8 + 4×4 pore")
        self.seed_cross_section_preset_display.setObjectName("fixedValue")
        self.seed_cross_section_preset_display.setAlignment(
            Qt.AlignmentFlag.AlignCenter)
        self.seed_cross_section_preset_display.setMinimumHeight(28)
        self.seed_cross_section_preset_display.setToolTip(
            "Seed 截面固定为 8×8 + 4×4 pore，不可选择或编辑。")
        seed_layout.addWidget(self.seed_cross_section_preset_display)
        self.seed_cross_section_picker = SeedCrossSectionPicker()
        self.seed_cross_section_picker.set_interactive(False)
        self.seed_cross_section_picker.setMinimumHeight(210)
        self.seed_cross_section_picker.setMaximumHeight(210)
        self.seed_cross_section_picker.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        seed_layout.addWidget(self.seed_cross_section_picker)
        self.seed_cross_section_status = QLabel(seed_box)
        self.seed_cross_section_status.hide()
        self.reset_seed_cross_section_button = QPushButton("恢复当前预设")
        self.reset_seed_cross_section_button.hide()
        left.addWidget(seed_box)

        # The Seed cross-section remains visible while the two compact
        # parameter sections above behave like a one-at-a-time accordion.
        self.symmetry_step_button.toggled.connect(
            symmetry_box.setVisible)
        self.twist_period_step_button.toggled.connect(
            target_box.setVisible)
        self.symmetry_step_button.setChecked(True)
        symmetry_box.show()
        target_box.hide()

        design_buttons = QVBoxLayout()
        design_buttons.setSpacing(7)
        self.confirm_design_basis_button = QPushButton(
            "4. 接受当前对称性与 Twist")
        self.confirm_design_basis_button.setObjectName("primaryButton")
        self.confirm_design_basis_button.setMinimumHeight(34)
        self.design_basis_next_button = QPushButton(
            "下一步：层长度与层间距")
        self.design_basis_next_button.setObjectName("primaryButton")
        self.design_basis_next_button.setEnabled(False)
        design_buttons.addWidget(self.confirm_design_basis_button)
        design_buttons.addWidget(self.design_basis_next_button)
        self.design_basis_legacy_buttons_layout = design_buttons
        left.addLayout(design_buttons)
        self.design_basis_action_status = self._action_feedback_label()
        left.addWidget(self.design_basis_action_status)
        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        left_scroll.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        left_scroll.setMinimumWidth(455)
        left_scroll.setMaximumWidth(535)
        left_scroll.setWidget(left_panel)
        splitter.addWidget(left_scroll)

        self.setup_preview_box = QGroupBox("二维点阵与 Moiré 预览")
        self.setup_preview_box.setObjectName("previewGroupBox")
        right = self.setup_preview_box
        right_layout = QVBoxLayout(right)
        self.setup_preview_layout = right_layout
        right_layout.setContentsMargins(10, 14, 8, 8)
        right_layout.setSpacing(0)
        self.setup_preview_parameters = QLabel(_preview_parameter_html())
        self.setup_preview_parameters.setObjectName("previewParameterBanner")
        self.setup_preview_parameters.setTextFormat(Qt.TextFormat.RichText)
        self.setup_preview_parameters.setAlignment(
            Qt.AlignmentFlag.AlignCenter)
        self.setup_preview_parameters.setMinimumHeight(82)
        right_layout.addWidget(self.setup_preview_parameters)
        self.setup_preview = MoireTopViewPreview()
        right_layout.addWidget(self.setup_preview, 1)
        splitter.addWidget(right)
        self.design_selection_splitter = splitter
        splitter.setSizes([460, 1000])
        return page

    def _build_design_parameter_page(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(0, 0, 0, 0)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        layout.addWidget(splitter)

        # Keep the whole design form visible without a scroll gesture.
        left_panel = QWidget()
        left = QVBoxLayout(left_panel)
        self.design_parameter_left_layout = left
        left.setContentsMargins(2, 0, 8, 4)
        left.setSpacing(7)
        left_panel.setMinimumWidth(455)
        left_panel.setMaximumWidth(535)
        splitter.addWidget(left_panel)

        sst_box = QGroupBox("层长度与层间距")
        self.sst_box = sst_box
        sst_box.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        sst_form = QFormLayout(sst_box)
        self._compact_form(sst_form)
        sst_form.setVerticalSpacing(8)
        self.sst_z1 = PhaseLengthCombo()
        self.sst_spacing = PhaseLengthCombo()
        self.sst_z3 = PhaseLengthCombo()
        self.layers_identical = QComboBox()
        self.layers_identical.addItem("一致", True)
        self.layers_identical.addItem("不一致", False)
        # Retain the computed phase text as internal state for project and
        # validation updates, but do not consume a visible form row.
        self.phase_hint = QTextBrowser(sst_box)
        self.phase_hint.hide()
        self.phase_hint.setObjectName("phaseHint")
        self.phase_hint.setFrameShape(QFrame.Shape.NoFrame)
        self.phase_hint.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.phase_hint.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        # Keep the explanatory row compact so the SST rows do not look
        # looser than the equally sized Seed parameter rows below.
        self.phase_hint.setMinimumHeight(48)
        self.phase_hint.setMaximumHeight(56)
        self.phase_hint.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.phase_hint.document().setDocumentMargin(0)
        sst_form.addRow("双层设计和序列", self.layers_identical)
        self.sst_z1_label = QLabel("1st layer length")
        self.sst_z3_label = QLabel("2nd layer length")
        self.sst_z1_label.setStyleSheet(
            "color:#2a78d1;font-weight:700")
        self.sst_z3_label.setStyleSheet(
            "color:#d65b74;font-weight:700")
        self.sst_spacing_label = QLabel("Layer spacing")
        self.sst_spacing_label.setStyleSheet(
            "color:#8a61bb;font-weight:700")
        sst_form.addRow(self.sst_z3_label, self.sst_z3)
        sst_form.addRow(self.sst_spacing_label, self.sst_spacing)
        sst_form.addRow(self.sst_z1_label, self.sst_z1)
        left.addWidget(sst_box)

        seed_box = QGroupBox("Seed S(F) 分区")
        self.seed_box = seed_box
        seed_box.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
        # Use an explicit grid here instead of QFormLayout.WrapLongRows so
        # every Seed label and its concise overlap result remain on one row
        # as the left pane is resized.
        seed_form = QGridLayout(seed_box)
        seed_form.setContentsMargins(9, 12, 9, 10)
        seed_form.setHorizontalSpacing(14)
        seed_form.setVerticalSpacing(9)
        seed_form.setColumnStretch(0, 0)
        seed_form.setColumnStretch(1, 1)
        # Compatibility-only hidden controls: legacy project files may still
        # contain Seed Z1/Z3 values, but the accepted two-layer Seed is now
        # immutable and those values never affect generation.
        self.seed_z1 = QSpinBox()
        self.seed_z1.setRange(128, 128)
        self.seed_z1.setValue(128)
        self.seed_z1.hide()
        self.seed_z2 = PhaseLengthCombo()
        self.seed_z2.setToolTip(
            "与SST superlattice spacing是同一个参数，"
            "任意一处修改都会同步。")
        self.seed_z2.setEnabled(False)
        self.seed_z2_readout = QLabel()
        self.seed_z2_readout.setObjectName("fixedValue")
        self.seed_z3 = QSpinBox()
        self.seed_z3.setRange(128, 128)
        self.seed_z3.setValue(128)
        self.seed_z3.hide()
        self.seed_z1_overlap_readout = QLabel(
            "128 bp · 8 capture columns (minimum 4)")
        self.seed_z1_overlap_readout.setObjectName("fixedValue")
        self.seed_z3_overlap_readout = QLabel(
            "128 bp · 8 capture columns (minimum 4)")
        self.seed_z3_overlap_readout.setObjectName("fixedValue")
        for readout in (self.seed_z1_overlap_readout,
                        self.seed_z2_readout,
                        self.seed_z3_overlap_readout):
            readout.setMinimumHeight(28)
            readout.setAlignment(
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            readout.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.mean_indel = QDoubleSpinBox()
        # This is a read-only result, not an input constraint.  The former
        # +/-12 range silently clipped the calibrated value while Twist is
        # intentionally unrestricted (for example, a 32-bp Z2 needs about
        # +17.1/helix at +20 degrees).  Keep the full solved value visible so
        # both this readout and the physical-Z2 readout update continuously.
        self.mean_indel.setRange(-8192.0, 8192.0)
        self.mean_indel.setDecimals(1)
        self.mean_indel.setSingleStep(0.1)
        self.mean_indel.setSuffix(
            " / helix (minimum -12, maximum +10)")
        self.mean_indel.setReadOnly(True)
        self.mean_indel.setButtonSymbols(
            QAbstractSpinBox.ButtonSymbols.NoButtons)
        self.actual_z2_spacing = QDoubleSpinBox()
        self.actual_z2_spacing.setRange(0.0, 8192.0)
        self.actual_z2_spacing.setDecimals(1)
        self.actual_z2_spacing.setSuffix(" bp")
        self.actual_z2_spacing.setReadOnly(True)
        self.actual_z2_spacing.setButtonSymbols(
            QAbstractSpinBox.ButtonSymbols.NoButtons)
        self.actual_z2_spacing.setToolTip(
            "名义spacing加上平均insertion/deletion；"
            "SST superlattice spacing与Seed Z2"
            "始终共享该实际长度。小数表示不同helices分配整数增删后的平均值。")
        self.seed_z1_label = QLabel("1st support · Z1")
        self.seed_z1_label.setStyleSheet(
            "color:#2a78d1;font-weight:700")
        self.seed_z2_label = QLabel("Seed Z2")
        self.seed_z2_label.setStyleSheet(
            "color:#8a61bb;font-weight:700")
        self.seed_z3_label = QLabel("2nd support · Z3")
        self.seed_z3_label.setStyleSheet(
            "color:#d65b74;font-weight:700")
        seed_form.addWidget(self.seed_z3_label, 0, 0)
        seed_form.addWidget(self.seed_z3_overlap_readout, 0, 1)
        seed_form.addWidget(self.seed_z2_label, 1, 0)
        seed_form.addWidget(self.seed_z2_readout, 1, 1)
        self.seed_z2.hide()
        seed_form.addWidget(self.seed_z1_label, 2, 0)
        seed_form.addWidget(self.seed_z1_overlap_readout, 2, 1)
        mean_indel_label = QLabel("Mean insertion/deletion")
        actual_spacing_label = QLabel("Actual Z2 / spacing")
        for label in (self.seed_z1_label, self.seed_z2_label,
                      self.seed_z3_label, mean_indel_label,
                      actual_spacing_label):
            label.setMinimumHeight(26)
            label.setAlignment(
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.mean_indel.setMinimumHeight(30)
        self.actual_z2_spacing.setMinimumHeight(30)
        seed_form.addWidget(mean_indel_label, 3, 0)
        seed_form.addWidget(self.mean_indel, 3, 1)
        seed_form.addWidget(actual_spacing_label, 4, 0)
        seed_form.addWidget(self.actual_z2_spacing, 4, 1)
        self.seed_scaffold_capacity_status = QLabel(
            "接受参数时按合法routing精确核算")
        self.seed_scaffold_capacity_status.setWordWrap(True)
        self.seed_scaffold_capacity_status.setObjectName("structureNote")
        self.seed_scaffold_capacity_status.hide()
        left.addWidget(seed_box)

        button_row = QVBoxLayout()
        button_row.setSpacing(7)
        self.change_design_basis_button = QPushButton("重新选择点阵 / Seed 截面")
        self.accept_parameters_button = QPushButton("接受当前 Moiré 参数")
        self.accept_parameters_button.setObjectName("primaryButton")
        self.parameters_next_button = QPushButton(
            "Next: Automated DNA Design")
        self.parameters_next_button.setObjectName("primaryButton")
        self.parameters_next_button.setEnabled(False)
        button_row.addWidget(self.change_design_basis_button)
        left.addLayout(button_row)
        parameter_navigation = QHBoxLayout()
        parameter_navigation.setSpacing(7)
        self.accept_parameters_button.setMinimumHeight(34)
        parameter_navigation.addWidget(self.accept_parameters_button)
        parameter_navigation.addWidget(self.parameters_next_button)
        left.addLayout(parameter_navigation)
        self.parameters_action_status = self._action_feedback_label()
        left.addWidget(self.parameters_action_status)

        self.side_preview_box = QGroupBox("三维 Seed/SST 预览")
        self.side_preview_box.setObjectName("previewGroupBox")
        right = self.side_preview_box
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(10, 14, 8, 8)
        right_layout.setSpacing(0)
        # Legacy value holders remain available to recalculate(), but must be
        # children of the 1.2 panel. Parentless QFrames are native top-level
        # windows on macOS and used to reappear when the main window closed.
        self.angle_card = ResultCard("当前 Twist angle", right)
        self.period_card = ResultCard("当前 Moiré period", right)
        self.angle_card.hide()
        self.period_card.hide()
        self.side_preview_parameters = QLabel(_side_preview_parameter_html())
        self.side_preview_parameters.setObjectName("previewParameterBanner")
        self.side_preview_parameters.setTextFormat(Qt.TextFormat.RichText)
        self.side_preview_parameters.setAlignment(
            Qt.AlignmentFlag.AlignCenter)
        self.side_preview_parameters.setMinimumHeight(82)
        self.preview = BilayerPreview()
        self.moire_preview = self.setup_preview
        right_layout.addWidget(self.side_preview_parameters)
        right_layout.addWidget(self.preview, 1)
        validation_header = QFrame()
        validation_header.setObjectName("validationHeader")
        validation_header_layout = QHBoxLayout(validation_header)
        validation_header_layout.setContentsMargins(9, 3, 5, 3)
        validation_header_layout.setSpacing(8)
        validation_title = QLabel("设计验证说明")
        validation_title.setObjectName("validationTitle")
        validation_header_layout.addWidget(validation_title)
        validation_header_layout.addStretch(1)
        self.validation_toggle = QPushButton("展开验证说明  ▴")
        self.validation_toggle.setObjectName("validationToggle")
        self.validation_toggle.setCheckable(True)
        self.validation_toggle.setChecked(True)
        self.validation_toggle.setToolTip("显示设计验证说明")
        validation_header_layout.addWidget(self.validation_toggle)
        validation_header.hide()
        self.validation = QTextBrowser()
        self.validation.setObjectName("validationPanel")
        self.validation.setMinimumHeight(0)
        self.validation.setMaximumHeight(0)
        self.validation.hide()
        splitter.addWidget(right)
        self.design_parameter_splitter = splitter
        splitter.setSizes([460, 1000])
        return tab

    def _merge_design_parameter_pages(self):
        """Compose the former 1.1/1.2 controls into one guided page."""
        self.sst_parameter_step_button = QPushButton(
            "1.3 输入 SST superlattice 参数")
        self.sst_parameter_step_button.setObjectName("parameterStepButton")
        self.sst_parameter_step_button.setCheckable(True)
        self.sst_parameter_step_button.setMinimumHeight(34)
        self.design_basis_section_group.addButton(
            self.sst_parameter_step_button, 2)

        self.sst_seed_parameter_container = QWidget()
        parameter_layout = QVBoxLayout(self.sst_seed_parameter_container)
        parameter_layout.setContentsMargins(0, 0, 0, 0)
        parameter_layout.setSpacing(7)
        parameter_layout.addWidget(self.sst_box)
        parameter_layout.addWidget(self.seed_box)

        # Arrange the three parameter steps like the compact 3.1/3.2/3.3
        # workflow on the sequence page.  Each expanded panel belongs directly
        # below its own button; the outer spacing stays identical between step
        # sections, while the smaller inner spacing keeps a button visually
        # attached to the content it controls.
        self.design_basis_section_group.setExclusive(False)
        self.design_parameter_sections = QWidget()
        self.design_parameter_sections.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        sections_layout = QVBoxLayout(self.design_parameter_sections)
        sections_layout.setContentsMargins(0, 0, 0, 0)
        sections_layout.setSpacing(12)
        sections_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        section_pairs = (
            (self.symmetry_step_button, self.symmetry_box),
            (self.twist_period_step_button, self.target_box),
            (self.sst_parameter_step_button,
             self.sst_seed_parameter_container),
        )
        self.design_parameter_section_widgets = []
        for button, content in section_pairs:
            # Identical frame height to the page-2 workflow buttons and the
            # Accept/Next controls below.
            button.setFixedHeight(34)
            button.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            section = QWidget()
            section.setSizePolicy(
                QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
            section_layout = QVBoxLayout(section)
            section_layout.setContentsMargins(0, 0, 0, 0)
            section_layout.setSpacing(7)
            section_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
            section_layout.addWidget(button)
            section_layout.addWidget(content)
            sections_layout.addWidget(section)
            self.design_parameter_section_widgets.append(section)
        self.design_basis_left_layout.insertWidget(
            0, self.design_parameter_sections)
        self.project_action_bar = QWidget()
        project_action_layout = QHBoxLayout(self.project_action_bar)
        project_action_layout.setContentsMargins(0, 0, 0, 0)
        project_action_layout.setSpacing(7)
        self.new_project_button = QPushButton("New Project")
        self.open_project_button = QPushButton("Open Project")
        for button in (self.new_project_button, self.open_project_button):
            button.setObjectName("projectActionButton")
            button.setFixedHeight(34)
            button.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            project_action_layout.addWidget(button, 1)
        self.design_basis_left_layout.insertWidget(
            0, self.project_action_bar)
        self.design_basis_left_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        # Match the 12-px outer rhythm used on page 2.  The 7-px navigation
        # spacing below is reserved for controls in the same action group.
        self.design_basis_left_layout.setSpacing(12)
        # The former split 1A/1B navigation layout contains only retired,
        # hidden buttons.  A childless QBoxLayout still advertises vertical
        # expansion on some Qt/macOS builds, producing the large gaps seen
        # above 1.1 and before Accept.  Remove the layout itself as well.
        self.design_basis_left_layout.removeItem(
            self.design_basis_legacy_buttons_layout)
        self.sst_parameter_step_button.toggled.connect(
            self.sst_seed_parameter_container.setVisible)
        for button in (self.symmetry_step_button,
                       self.twist_period_step_button,
                       self.sst_parameter_step_button):
            button.setChecked(False)
        self.symmetry_box.hide()
        self.target_box.hide()
        self.sst_seed_parameter_container.hide()

        self.confirm_design_basis_button.hide()
        self.design_basis_next_button.hide()
        self.design_basis_action_status.hide()
        self.change_design_basis_button.hide()

        self.accept_parameters_button.setText(
            "接受当前 Moiré 参数")
        self.accept_parameters_button.setFixedHeight(34)
        self.parameters_next_button.setText(
            "Next: Automated DNA Design")
        self.parameters_next_button.setFixedHeight(34)
        combined_navigation = QVBoxLayout()
        combined_navigation.setSpacing(7)
        combined_navigation.addWidget(self.accept_parameters_button)
        # Keep one acceptance explanation only, directly after Accept.  The
        # legacy design-basis message belongs to the retired split-page flow
        # and must never reappear on the combined parameter page.
        combined_navigation.addWidget(self.parameters_action_status)
        combined_navigation.addWidget(self.parameters_next_button)
        self.design_parameter_navigation_layout = combined_navigation
        self.design_basis_left_layout.addLayout(combined_navigation)
        self.design_basis_left_layout.addStretch(1)

        # Build one combined preview area.  The 2D and 3D canvases are shown
        # side by side by default, with a single shared title and a single
        # Twist/Moiré-period banner above them.
        # Restore the original outer preview frame/title while keeping the
        # complete 2D+3D field visually continuous and black inside it.
        combined_preview = QGroupBox("二维点阵与 Moiré 预览")
        combined_preview.setObjectName("previewGroupBox")
        self.design_combined_preview_box = combined_preview
        combined_preview_layout = QVBoxLayout(combined_preview)
        combined_preview_layout.setContentsMargins(10, 10, 10, 8)
        # The title, parameter banner, and two-view matrix form one continuous
        # black preview surface, without white bands between rows.
        combined_preview_layout.setSpacing(0)
        self.design_preview_surface = QFrame()
        self.design_preview_surface.setObjectName("previewSurface")
        preview_surface_layout = QVBoxLayout(self.design_preview_surface)
        preview_surface_layout.setContentsMargins(0, 0, 0, 0)
        preview_surface_layout.setSpacing(0)
        combined_preview_layout.addWidget(self.design_preview_surface, 1)
        self.seed_cross_section_box.hide()
        self.design_preview_title = QLabel(
            "Square–Square Bilayer DNA Moiré Superlattice")
        self.design_preview_title.setObjectName("previewCommonTitle")
        self.design_preview_title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.design_preview_title.setMinimumHeight(44)
        self.design_preview_title.setContentsMargins(0, 4, 0, 0)
        preview_surface_layout.addWidget(self.design_preview_title)

        # Retain compatibility value holders for recalculate(), but remove the
        # visible dimension-summary row requested by the user.
        self.design_preview_summary = QFrame()
        self.design_preview_summary.setObjectName("previewSummary")
        summary_layout = QHBoxLayout(self.design_preview_summary)
        summary_layout.setContentsMargins(12, 5, 12, 5)
        summary_layout.setSpacing(18)
        self.design_preview_geometry_summary = QLabel()
        self.design_preview_geometry_summary.setAlignment(
            Qt.AlignmentFlag.AlignCenter)
        self.design_preview_length_summary = QLabel()
        self.design_preview_length_summary.setAlignment(
            Qt.AlignmentFlag.AlignCenter)
        summary_layout.addWidget(self.design_preview_geometry_summary, 1)
        summary_layout.addWidget(self.design_preview_length_summary, 1)
        self.design_preview_summary.hide()

        def view_heading(text, pane):
            heading = QFrame()
            heading.setObjectName("previewViewHeading")
            heading.setProperty("previewPane", pane)
            heading_layout = QHBoxLayout(heading)
            heading_layout.setContentsMargins(10, 4, 10, 4)
            heading_layout.setSpacing(0)
            subtitle = QLabel(text)
            subtitle.setObjectName("previewViewSubtitle")
            subtitle.setProperty("previewPane", pane)
            subtitle.setAlignment(Qt.AlignmentFlag.AlignCenter)
            heading_layout.addWidget(subtitle, 1)
            return heading

        def legend_panel(items, swatch_shape, pane):
            legend = QFrame()
            legend.setObjectName("previewLegend")
            legend.setProperty("previewPane", pane)
            legend_layout = QHBoxLayout(legend)
            legend_layout.setContentsMargins(8, 6, 8, 6)
            legend_layout.setSpacing(7)
            legend_layout.addStretch(1)
            for label, color in items:
                swatch = QFrame()
                if swatch_shape == "circle":
                    swatch.setFixedSize(10, 10)
                    swatch.setStyleSheet(
                        "background:%s;border:0;border-radius:5px" % color)
                elif swatch_shape == "vertical":
                    swatch.setFixedSize(4, 20)
                    swatch.setStyleSheet(
                        "background:%s;border:0;border-radius:2px" % color)
                else:
                    swatch.setFixedSize(22, 4)
                    swatch.setStyleSheet(
                        "background:%s;border:0;border-radius:2px" % color)
                legend_layout.addWidget(swatch)
                legend_layout.addWidget(QLabel(label))
                legend_layout.addSpacing(7)
            legend_layout.addStretch(1)
            return legend

        self.design_preview_matrix = QFrame()
        self.design_preview_matrix.setObjectName("previewMatrix")
        matrix = QGridLayout(self.design_preview_matrix)
        self.design_preview_grid_layout = matrix
        # Leave a small quiet band under the common title. Both parameter
        # banners begin on this same baseline.
        matrix.setContentsMargins(0, 4, 0, 0)
        matrix.setHorizontalSpacing(0)
        matrix.setVerticalSpacing(0)
        matrix.setColumnStretch(0, 1)
        matrix.setColumnStretch(2, 1)
        matrix.setRowStretch(0, 1)

        self.design_preview_top_heading = view_heading("Top view", "top")
        self.design_preview_side_heading = view_heading("Side view", "side")
        # Keep the compatibility objects for older integrations, but the
        # canvases no longer repeat the redundant Top/Side view captions.
        self.design_preview_top_heading.hide()
        self.design_preview_side_heading.hide()
        # Values are now annotated directly on their respective canvases.
        # Keep these labels as hidden compatibility holders because
        # recalculate() still refreshes them for older integrations.
        self.setup_preview_parameters.hide()
        self.side_preview_parameters.hide()
        self.design_preview_vertical_separator = QFrame()
        self.design_preview_vertical_separator.setObjectName(
            "previewVerticalSeparator")
        self.design_preview_vertical_separator.setFixedWidth(1)
        self.design_preview_vertical_separator.hide()

        # The two preview widgets keep their existing rendering code, but the
        # inner white QGroupBox chrome and subtitles are replaced by the shared
        # plain heading row above.
        self.setup_preview_box.setTitle("")
        self.side_preview_box.setTitle("")
        self.setup_preview_box.setObjectName("previewCanvasBox")
        self.side_preview_box.setObjectName("previewCanvasBox")
        self.setup_preview_box.setProperty("previewPane", "top")
        self.side_preview_box.setProperty("previewPane", "side")
        self.setup_preview.set_preview_background("#070b10")
        self.preview.set_preview_background("#070b10")
        self.setup_preview_box.layout().setContentsMargins(0, 0, 0, 0)
        self.side_preview_box.layout().setContentsMargins(0, 0, 0, 0)
        # Both panes occupy one grid row and use identical vertical geometry.
        # Explicit expanding policies prevent either legacy page's size hint
        # from making the 2D canvas open shorter than the 3D canvas.
        for preview_box, preview_canvas in (
                (self.setup_preview_box, self.setup_preview),
                (self.side_preview_box, self.preview)):
            preview_box.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
            preview_canvas.setMinimumHeight(260)
            preview_canvas.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        matrix.addWidget(self.setup_preview_box, 0, 0)
        matrix.addWidget(self.side_preview_box, 0, 2)

        self.design_preview_2d_legend = legend_panel((
            ("SST sublattice 1st layer", "#2a78d1"),
            ("Seed", "#ffffff"),
            ("SST sublattice 2nd layer", "#d65b74"),
        ), "circle", "top")
        self.design_preview_3d_legend = legend_panel((
            ("Seed Z1", "#2a78d1"),
            ("Seed Z2", "#d9dee3"),
            ("Seed Z3", "#d65b74"),
        ), "vertical", "side")
        matrix.addWidget(self.design_preview_2d_legend, 1, 0)
        matrix.addWidget(self.design_preview_3d_legend, 1, 2)
        preview_surface_layout.addWidget(self.design_preview_matrix, 1)
        # Moving setup_preview_box into the vertical splitter removes the old
        # second widget from design_selection_splitter automatically.
        self.design_selection_splitter.addWidget(combined_preview)
        self.design_selection_splitter.setSizes([460, 1000])

    @staticmethod
    def _compact_form(form):
        form.setContentsMargins(9, 9, 9, 8)
        form.setHorizontalSpacing(8)
        form.setVerticalSpacing(5)
        form.setFormAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
        form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)

    def _build_capture_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        # Match the lower page inset to the upper page inset.
        layout.setContentsMargins(0, 8, 0, 8)
        self.capture_splitter = QSplitter(Qt.Orientation.Horizontal)
        layout.addWidget(self.capture_splitter, 1)

        scroll, right = self._scroll_panel()
        scroll.setMinimumWidth(455)
        scroll.setMaximumWidth(535)
        scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.capture_splitter.addWidget(scroll)
        self.capture_preview_box = QGroupBox("Embedded caDNAno design view")
        self.capture_preview_box.setObjectName("previewGroupBox")
        preview_panel = self.capture_preview_box
        preview_panel.setMinimumWidth(360)
        preview_layout = QVBoxLayout(preview_panel)
        # Match page 1: the preview canvas begins immediately below the group
        # title with no extra channel-selection row or blank header band.
        preview_layout.setContentsMargins(10, 10, 8, 8)
        preview_layout.setSpacing(4)
        # Keep the generated-stage selector as hidden internal state so the
        # newest available design is still selected automatically.
        self.structure_preview_channel = QComboBox(preview_panel)
        self.structure_preview_channel.setMinimumWidth(255)
        self.structure_preview_channel.setToolTip(
            "只列出当前已生成的结构层级")
        self.structure_preview_channel.hide()
        # Keep the design report in a full-width, vertically resizable row
        # below the three preview panels.  This matches the SST input report
        # typography and leaves the Path viewport at its full width.
        self.capture_results_splitter = QSplitter(
            Qt.Orientation.Vertical, preview_panel)
        self.capture_results_splitter.setChildrenCollapsible(True)
        self.capture_results_splitter.setHandleWidth(10)
        self.capture_results_splitter.setOpaqueResize(True)
        self.capture_preview = StructureDesignPreview()
        self.capture_preview.setMinimumHeight(0)
        self.capture_preview.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)
        self.capture_results_splitter.addWidget(self.capture_preview)
        self.structure_preview_status = QLabel("", preview_panel)
        self.structure_preview_status.setWordWrap(True)
        self.structure_preview_status.setObjectName("structureNote")
        self.structure_preview_status.setMinimumHeight(0)
        self.structure_preview_status.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)
        self.structure_preview_status.hide()
        self.capture_results_splitter.addWidget(
            self.structure_preview_status)
        self.capture_results_splitter.setCollapsible(0, False)
        self.capture_results_splitter.setCollapsible(1, True)
        self.capture_results_splitter.setStretchFactor(0, 5)
        self.capture_results_splitter.setStretchFactor(1, 1)
        self.capture_results_splitter.setSizes([760, 130])
        preview_layout.addWidget(self.capture_results_splitter, 1)
        self.capture_splitter.addWidget(preview_panel)
        self.capture_splitter.setCollapsible(0, True)
        self.capture_splitter.setCollapsible(1, False)
        self.capture_splitter.setStretchFactor(0, 2)
        self.capture_splitter.setStretchFactor(1, 5)
        self.capture_splitter.setSizes([460, 1000])

        intro = QLabel(
            "本模块可独立打开。导入功能仅接受 DNA Moiré Designer "
            "工程（.moire.json），用于恢复设计参数与流程。\n"
            "1. 生成 Scaffold routing 需要已接受的设计参数\n"
            "2. 生成 Staple / Capture 需要已接受或已导入的 Scaffold routing\n"
            "3. Z1/Z2/Z3 按8-bp domain和原Square相位规则生成；"
            "长度增加时会自动增加容量安全的scaffold分段")
        intro.setWordWrap(True)
        intro.setObjectName("structureIntro")
        intro.hide()
        self.import_capture_project_button = QPushButton(
            "导入 Moiré 工程 (.moire.json)")
        self.import_capture_project_button.setObjectName("primaryButton")
        self.import_capture_project_button.hide()
        self.capture_import_status = QLabel(
            "仅可导入本软件保存的 .moire.json 工程文件。")
        self.capture_import_status.setWordWrap(True)
        self.capture_import_status.hide()
        self.accepted_parameters_summary = QLabel("暂无设计参数")
        self.accepted_parameters_summary.setWordWrap(True)
        self.accepted_parameters_summary.setObjectName("structureIntro")
        right.addWidget(self.accepted_parameters_summary)
        self.generate_simple_design_button = QPushButton(
            "生成并导出全部 3 个设计文件")
        self.generate_simple_design_button.setObjectName("primaryButton")
        self.generate_simple_design_button.setEnabled(False)
        self.generate_simple_design_button.setToolTip(
            "Generate the SST sublattice-only, SST sublattice + scaffold, "
            "and final SST sublattice + scaffold + staple + capture JSON "
            "files.")
        right.addWidget(self.generate_simple_design_button)
        self.design_generation_action_status = self._action_feedback_label()
        right.addWidget(self.design_generation_action_status)
        self.structure_expert_button = QPushButton(
            "Optional · expert mode")
        self.structure_expert_button.setCheckable(True)
        self.structure_expert_button.setObjectName("optionalButton")
        self.structure_expert_button.hide()
        self.back_to_parameters_button = QPushButton()
        self.back_to_parameters_button.hide()

        scaffold_box = QGroupBox("2.1 生成并审核 Scaffold routing")
        self.scaffold_expert_box = scaffold_box
        scaffold_layout = QVBoxLayout(scaffold_box)
        self.generate_scaffold_button = QPushButton(
            "生成固定 SST superlattice + Scaffold routing")
        self.generate_scaffold_button.setObjectName("primaryButton")
        self.generate_scaffold_button.setEnabled(True)
        self.open_scaffold_button = QPushButton(
            "在 cadnano 内专家编辑完成 Scaffold routing")
        self.open_scaffold_button.setEnabled(False)
        self.open_scaffold_button.setToolTip(
            "cadnano 会打开当前待审核文件；请直接保存该文件，返回后点击接受。")
        self.load_scaffold_button = QPushButton(
            "载入 cadnano 专家编辑后的 JSON")
        self.load_scaffold_button.hide()
        self.accept_scaffold_button = QPushButton(
            "接受 cadnano 当前保存的 Scaffold routing")
        self.accept_scaffold_button.setObjectName("primaryButton")
        self.accept_scaffold_button.setEnabled(False)
        self.accept_scaffold_button.setToolTip(
            "自动重新读取、验证并接受刚才在 cadnano 中保存的同一个文件。")
        self.scaffold_status = QLabel("等待生成。")
        self.scaffold_status.setWordWrap(True)
        scaffold_layout.addWidget(self.generate_scaffold_button)
        scaffold_layout.addWidget(self.open_scaffold_button)
        scaffold_layout.addWidget(self.accept_scaffold_button)
        scaffold_layout.addWidget(self.scaffold_status)
        right.addWidget(scaffold_box)
        scaffold_box.hide()

        staple_box = QGroupBox("2.2 自动生成 Staple / Capture 结构")
        self.staple_expert_box = staple_box
        staple_layout = QVBoxLayout(staple_box)
        self.back_to_scaffold_button = QPushButton()
        self.back_to_scaffold_button.hide()
        self.generate_structure_button = QPushButton(
            "生成 Staple / Capture 设计")
        self.generate_structure_button.setEnabled(True)
        self.open_structure_button = QPushButton(
            "在 cadnano 内专家编辑完成结构")
        self.open_structure_button.setEnabled(False)
        self.accept_structure_button = QPushButton(
            "接受当前设计图")
        self.accept_structure_button.setObjectName("primaryButton")
        self.accept_structure_button.setEnabled(False)
        self.structure_status = QLabel(
            "接受scaffold后才可生成。capture每个pair使用同一种颜色。")
        self.structure_status.setWordWrap(True)
        staple_layout.addWidget(self.generate_structure_button)
        staple_layout.addWidget(self.open_structure_button)
        staple_layout.addWidget(self.structure_status)
        right.addWidget(staple_box)
        staple_box.hide()

        self.inspect_final_design_button = QPushButton(
            "Optional：打开 cadnano")
        self.inspect_final_design_button.setObjectName("optionalButton")
        # cadnano expert editing is meaningful only after the final
        # SST + Scaffold + Staple + Capture design has been generated.
        self.inspect_final_design_button.setEnabled(False)
        self.inspect_final_design_button.setToolTip(
            "直接打开 cadnano；请在 cadnano 中打开并修改当前项目 cadnano design "
            "文件夹内的最终 JSON。点击接受时，软件会自动选择该文件夹中"
            "修改时间最新的合法 SST + Scaffold + Staple + Capture 文件。")
        right.addWidget(self.inspect_final_design_button)
        self.cadnano_edit_action_status = self._action_feedback_label()
        right.addWidget(self.cadnano_edit_action_status)

        self.accept_structure_button.setMinimumHeight(34)
        right.addWidget(self.accept_structure_button)
        self.structure_accept_action_status = self._action_feedback_label()
        right.addWidget(self.structure_accept_action_status)
        self.structure_next_button = QPushButton("下一步：序列导出")
        self.structure_next_button.setObjectName("primaryButton")
        self.structure_next_button.setMinimumHeight(34)
        self.structure_next_button.setEnabled(False)
        right.addWidget(self.structure_next_button)

        note = QLabel(
            "固定SST superlattice在生成Scaffold routing时自动"
            "加入，无需单独审核。"
            "序列与导出模块可载入带序列 JSON，并分别导出"
            "Capture，以及只有完整32-nt SST的 JSON、XLSX 和 SVG。")
        note.setWordWrap(True)
        note.setObjectName("structureNote")
        note.hide()
        right.addStretch(1)
        return tab

    def _build_sequence_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        # Match the lower page inset to the upper page inset.
        layout.setContentsMargins(0, 8, 0, 8)
        self.back_to_structure_button = QPushButton()
        self.back_to_structure_button.hide()
        self.sequence_splitter = QSplitter(Qt.Orientation.Horizontal)
        layout.addWidget(self.sequence_splitter, 1)

        left_scroll, left = self._scroll_panel()
        left_scroll.setMinimumWidth(455)
        left_scroll.setMaximumWidth(535)
        left_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left_scroll.widget().setMinimumWidth(410)
        self.sequence_splitter.addWidget(left_scroll)

        self.accepted_design_summary = QLabel("尚未接受设计图")
        self.accepted_design_summary.setWordWrap(True)
        self.accepted_design_summary.setObjectName("structureIntro")
        left.addWidget(self.accepted_design_summary)

        self.sequence_section_group = QButtonGroup(self)
        self.sequence_section_group.setExclusive(False)

        def sequence_step_button(text):
            button = QPushButton(text)
            button.setObjectName("parameterStepButton")
            button.setCheckable(True)
            button.setFixedHeight(34)
            self.sequence_section_group.addButton(button)
            return button

        # Keep the hidden legacy container alive until the sequence tab itself
        # is destroyed.  Without this ownership, Qt deletes its child button
        # at the end of this method and _connect_signals() receives a stale
        # Python wrapper during application startup.
        import_box = QGroupBox("导入 Moiré 工程", tab)
        self.sequence_import_box = import_box
        import_layout = QVBoxLayout(import_box)
        self.import_sequence_project_button = QPushButton(
            "导入 Moiré 工程 (.moire.json)")
        self.import_sequence_project_button.setObjectName("primaryButton")
        import_layout.addWidget(self.import_sequence_project_button)
        self.sequence_design_status = QLabel(
            "仅可导入本软件保存的 .moire.json，恢复参数与工作流。")
        self.sequence_design_status.setWordWrap(True)
        import_layout.addWidget(self.sequence_design_status)
        import_box.hide()

        scaffold_box = QGroupBox()
        self.sequence_scaffold_section_content = scaffold_box
        scaffold_layout = QVBoxLayout(scaffold_box)
        scaffold_note = QLabel(
            "Detect the position and length of each scaffold route in the "
            "accepted design. Then assign one built-in caDNAno scaffold "
            "sequence to each route.")
        scaffold_note.setWordWrap(True)
        self.detect_scaffold_sequences_button = QPushButton(
            "自动读取 Scaffold 位置和长度")
        self.detect_scaffold_sequences_button.setObjectName("primaryButton")
        self.detect_scaffold_sequences_button.setEnabled(False)
        self.accept_added_scaffold_button = QPushButton(
            "Accept assigned scaffold sequences")
        self.accept_added_scaffold_button.setObjectName("primaryButton")
        self.accept_added_scaffold_button.setEnabled(False)
        self.scaffold_sequence_status = QLabel("等待读取结构。")
        self.scaffold_sequence_status.setWordWrap(True)
        self.scaffold_sequence_status.setObjectName("successStatus")
        self.scaffold_sequence_status.hide()
        self.scaffold_detection_action_status = \
            self._action_feedback_label()
        scaffold_note.hide()
        scaffold_layout.addWidget(self.detect_scaffold_sequences_button)
        scaffold_layout.addWidget(self.scaffold_detection_action_status)
        self.scaffold_cards_widget = QWidget()
        self.scaffold_cards_layout = QVBoxLayout(
            self.scaffold_cards_widget)
        self.scaffold_cards_layout.setContentsMargins(0, 0, 0, 0)
        self.scaffold_cards_layout.addWidget(QLabel("尚未读取结构。"))
        scaffold_layout.addWidget(self.scaffold_cards_widget)
        scaffold_layout.addWidget(self.accept_added_scaffold_button)
        scaffold_layout.addWidget(self.scaffold_sequence_status)

        sst_box = QGroupBox()
        self.sequence_sst_section_content = sst_box
        sst_layout = QVBoxLayout(sst_box)
        self.sst_sequence_layout = sst_layout
        sst_note = QLabel(
            "Detect the position, count, and length of the SST sublattice "
            "inputs in each layer. Assign sequences automatically, or use "
            "expert mode to export a template and import specified "
            "sequences. Identical layers require one input set, which is "
            "mapped automatically to the corresponding positions in the "
            "other layer.")
        sst_note.setWordWrap(True)
        self.detect_sst_inputs_button = QPushButton(
            "Detect SST sublattice input positions and lengths")
        self.detect_sst_inputs_button.setObjectName("primaryButton")
        self.auto_design_sst_inputs_button = QPushButton(
            "Design and assign SST sublattice inputs automatically")
        self.auto_design_sst_inputs_button.setObjectName("primaryButton")
        self.auto_design_sst_inputs_button.setEnabled(False)
        self.sequence_expert_button = QPushButton(
            "Optional: expert mode")
        self.sequence_expert_button.setObjectName("optionalButton")
        self.sequence_expert_button.setCheckable(True)
        # Orthogonal Sequence Design is independent of the structure and
        # remains reachable before scaffold or SST-input preparation.
        self.sequence_expert_button.setEnabled(True)
        self.export_sst_input_template_button = QPushButton(
            "Export input template")
        self.export_sst_input_template_button.setObjectName("primaryButton")
        self.import_sst_input_template_button = QPushButton(
            "Import and assign input sequences")
        self.import_sst_input_template_button.setObjectName("primaryButton")
        self.accept_added_sst_button = QPushButton(
            "Accept assigned SST sublattice input sequences")
        self.accept_added_sst_button.setObjectName("primaryButton")
        self.accept_added_sst_button.setEnabled(False)
        self.sst_detection_status = QLabel(
            "Accept the assigned scaffold sequences first.")
        self.sst_detection_status.setWordWrap(True)
        self.sst_detection_status.setObjectName("cardCaption")
        # Compatibility alias for older call sites and saved UI integrations.
        self.sst_sequence_status = self.sst_detection_status
        self.sst_auto_import_status = QLabel()
        self.sst_auto_import_status.setWordWrap(True)
        self.sst_auto_import_status.setObjectName("successStatus")
        self.sst_auto_import_status.hide()
        self.sst_expert_import_status = QLabel()
        self.sst_expert_import_status.setWordWrap(True)
        self.sst_expert_import_status.setObjectName("successStatus")
        self.sst_expert_import_status.hide()
        self.sst_acceptance_status = QLabel()
        self.sst_acceptance_status.setWordWrap(True)
        self.sst_acceptance_status.setObjectName("successStatus")
        self.sst_acceptance_status.hide()
        for widget in (self.detect_sst_inputs_button,
                       self.export_sst_input_template_button,
                       self.import_sst_input_template_button):
            widget.setEnabled(False)
        sst_note.hide()
        sst_layout.addWidget(self.detect_sst_inputs_button)
        sst_layout.addWidget(self.sst_detection_status)
        self.sst_cards_widget = QWidget()
        self.sst_cards_layout = QVBoxLayout(self.sst_cards_widget)
        self.sst_cards_layout.setContentsMargins(0, 0, 0, 0)
        self.sst_cards_layout.addWidget(
            QLabel("SST sublattice inputs have not been detected."))
        sst_layout.addWidget(self.sst_cards_widget)
        sst_layout.addWidget(self.auto_design_sst_inputs_button)
        sst_layout.addWidget(self.sst_auto_import_status)
        sst_layout.addWidget(self.sequence_expert_button)
        sst_layout.addWidget(self.sst_expert_import_status)
        sst_layout.addWidget(self.accept_added_sst_button)
        sst_layout.addWidget(self.sst_acceptance_status)
        self.open_orthogonal_sequences_button = QPushButton(
            "Orthogonal sequence design")
        self.open_orthogonal_sequences_button.setObjectName("primaryButton")

        export_box = QGroupBox()
        self.sequence_export_section_content = export_box
        export_layout = QVBoxLayout(export_box)
        export_note = QLabel(
            "Export the final caDNAno design, oligonucleotide sequence "
            "tables, input parameters, and 3D structure files.")
        export_note.setWordWrap(True)
        self.final_sequence_export_button = QPushButton(
            "Export final package")
        self.final_sequence_export_button.setObjectName("primaryButton")
        self.final_sequence_export_button.setEnabled(False)
        self.sequence_export_status = QLabel(
            "Sequence assignment is not complete.")
        self.sequence_export_status.setWordWrap(True)
        export_note.hide()
        export_layout.addWidget(self.final_sequence_export_button)
        export_layout.addWidget(self.sequence_export_status)

        self.sequence_scaffold_step_button = sequence_step_button(
            "3.1 Assign scaffold sequences")
        self.sequence_sst_step_button = sequence_step_button(
            "3.2 Assign SST sublattice input sequences")
        self.sequence_export_step_button = sequence_step_button(
            "3.3 Final export")
        self.sequence_section_widgets = []
        for button, content in (
                (self.sequence_scaffold_step_button, scaffold_box),
                (self.sequence_sst_step_button, sst_box),
                (self.sequence_export_step_button, export_box)):
            section = QWidget()
            section.setSizePolicy(
                QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
            section_layout = QVBoxLayout(section)
            section_layout.setContentsMargins(0, 0, 0, 0)
            section_layout.setSpacing(7)
            section_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
            content.setSizePolicy(
                QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
            section_layout.addWidget(button)
            section_layout.addWidget(content)
            button.toggled.connect(content.setVisible)
            button.setChecked(False)
            content.hide()
            left.addWidget(section)
            self.sequence_section_widgets.append(section)
        left.addStretch(1)

        right_panel = QWidget()
        right_panel.setMinimumWidth(360)
        right = QVBoxLayout(right_panel)
        # The tab itself owns the symmetric top/bottom page inset; do not
        # stack another lower margin on the preview column.
        right.setContentsMargins(0, 0, 0, 0)
        right.setSpacing(6)
        self.sequence_results_splitter = QSplitter(Qt.Orientation.Vertical)
        self.sequence_results_splitter.setChildrenCollapsible(True)
        self.sequence_results_splitter.setHandleWidth(12)
        self.sequence_results_splitter.setOpaqueResize(True)
        right.addWidget(self.sequence_results_splitter, 1)

        self.sequence_preview_box = QGroupBox(
            "Sequence position and structure preview")
        self.sequence_preview_box.setObjectName("previewGroupBox")
        self.sequence_preview_box.setMinimumHeight(0)
        self.sequence_preview_box.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Ignored)
        preview_box_layout = QVBoxLayout(self.sequence_preview_box)
        # Match pages 1 and 2: no explanatory banner above the canvas.
        preview_box_layout.setContentsMargins(10, 10, 8, 8)
        preview_box_layout.setSpacing(4)
        self.sequence_preview_status = QLabel(
            "After detection, the scaffold-only or SST sublattice "
            "input-only design is displayed here.")
        self.sequence_preview_status.setWordWrap(True)
        self.sequence_preview_status.setObjectName("structureNote")
        self.sequence_preview_status.hide()
        self.sequence_preview = StructureDesignPreview()
        self.sequence_preview.setMinimumHeight(0)
        self.sequence_preview.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)
        preview_box_layout.addWidget(self.sequence_preview, 1)
        self.sequence_results_splitter.addWidget(self.sequence_preview_box)
        self.sst_sequence_table_box = QGroupBox(
            "SST sublattice input sequence analysis")
        table_layout = QVBoxLayout(self.sst_sequence_table_box)
        table_layout.setContentsMargins(10, 14, 8, 8)
        table_layout.setSpacing(6)
        self.sst_sequence_report = QLabel(
            "No SST sublattice input sequences have been imported.")
        self.sst_sequence_report.setWordWrap(True)
        self.sst_sequence_report.setObjectName("structureNote")
        self.sst_sequence_report.setMinimumHeight(0)
        self.sst_sequence_report.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)
        self.sst_analysis_splitter = QSplitter(Qt.Orientation.Vertical)
        self.sst_analysis_splitter.setChildrenCollapsible(True)
        self.sst_analysis_splitter.setHandleWidth(10)
        self.sst_analysis_splitter.setOpaqueResize(True)
        self.sst_analysis_splitter.addWidget(self.sst_sequence_report)
        self.sst_sequence_table = QTableWidget(0, 6)
        self.sst_sequence_table.setHorizontalHeaderLabels([
            "Position (5′→3′)",
            "Sequence",
            "Length (nt)",
            "GC content (%)",
            "Maximum same-orientation exact match (nt)",
            "Maximum interstrand complementarity (nt)",
        ])
        self.sst_sequence_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self.sst_sequence_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.sst_sequence_table.setAlternatingRowColors(True)
        self.sst_sequence_table.verticalHeader().setVisible(False)
        header = self.sst_sequence_table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setMinimumSectionSize(70)
        header.setStretchLastSection(False)
        for column, width in enumerate((210, 320, 120, 160, 360, 360)):
            header.resizeSection(column, width)
        self.sst_sequence_table.setMinimumHeight(0)
        self.sst_sequence_table.setMaximumHeight(16777215)
        self.sst_sequence_table.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Ignored)
        self.sst_analysis_splitter.addWidget(self.sst_sequence_table)
        self.sst_analysis_splitter.setCollapsible(0, True)
        self.sst_analysis_splitter.setCollapsible(1, True)
        self.sst_analysis_splitter.setStretchFactor(0, 1)
        self.sst_analysis_splitter.setStretchFactor(1, 3)
        self.sst_analysis_splitter.setSizes([90, 270])
        table_layout.addWidget(self.sst_analysis_splitter, 1)
        self.sst_sequence_table_box.setMinimumHeight(0)
        self.sst_sequence_table_box.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Ignored)
        self.sst_sequence_table_box.hide()
        self.sequence_results_splitter.addWidget(
            self.sst_sequence_table_box)
        self.sst_template_actions = QGroupBox(
            "SST sublattice input template")
        template_layout = QVBoxLayout(self.sst_template_actions)
        template_layout.setContentsMargins(10, 14, 8, 10)
        template_layout.setSpacing(9)
        self.sst_template_note = QLabel(
            "专家流程：先导出模板，也可用正交序列设计生成序列；"
            "填入 Sequence 列后再导入。行顺序按完整 base 数值优先、"
            "helix 数值其次排列。")
        self.sst_template_note.setWordWrap(True)
        self.sst_template_note.setMinimumHeight(68)
        self.sst_template_note.setContentsMargins(0, 5, 0, 8)
        self.sst_template_note.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        note_policy = QSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        note_policy.setHeightForWidth(True)
        self.sst_template_note.setSizePolicy(note_policy)
        template_layout.addWidget(self.sst_template_note)
        self.export_template_action_status = self._action_feedback_label()
        self.orthogonal_action_status = self._action_feedback_label()
        self.import_template_action_status = self._action_feedback_label()
        for button, status in (
                (self.export_sst_input_template_button,
                 self.export_template_action_status),
                (self.open_orthogonal_sequences_button,
                 self.orthogonal_action_status),
                (self.import_sst_input_template_button,
                 self.import_template_action_status)):
            button.setMinimumHeight(40)
            button.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            template_layout.addWidget(button)
            template_layout.addWidget(status)
        self.sst_template_actions.setMinimumHeight(0)
        self.sst_template_actions.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Ignored)
        self.sst_template_actions.hide()
        self.sequence_orthogonal_box = self.sst_template_actions
        self.sequence_results_splitter.addWidget(self.sst_template_actions)
        for index in range(3):
            self.sequence_results_splitter.setCollapsible(index, True)
        self.sequence_results_splitter.setStretchFactor(0, 5)
        self.sequence_results_splitter.setStretchFactor(1, 3)
        self.sequence_results_splitter.setStretchFactor(2, 3)
        # Give the expert instructions enough initial height for wrapped
        # academic English; users can still resize or collapse every panel.
        self.sequence_results_splitter.setSizes([430, 260, 340])
        self.sequence_splitter.addWidget(right_panel)
        self.sequence_splitter.setCollapsible(0, True)
        self.sequence_splitter.setCollapsible(1, False)
        self.sequence_splitter.setStretchFactor(0, 2)
        self.sequence_splitter.setStretchFactor(1, 5)
        self.sequence_splitter.setSizes([460, 1000])

        # Old attribute names remain aliases so older saved UI state and
        # integration code do not crash while the staged flow is adopted.
        self.load_sequence_design_button = self.detect_scaffold_sequences_button
        self.export_sequence_variants_button = self.final_sequence_export_button
        # Compatibility-only action target. Data Analysis is entered from the
        # independent main navigation drop-down, not as a continuation of step 3.
        self.go_to_analysis_button = QPushButton()
        self.go_to_analysis_button.hide()
        return tab

    def _legacy_build_analysis_tab(self):
        tab = QWidget()
        outer = QVBoxLayout(tab)
        outer.setContentsMargins(8, 10, 8, 8)
        navigation = QHBoxLayout()
        self.back_to_sequence_button = QPushButton("返回第3步：序列与导出")
        navigation.addWidget(self.back_to_sequence_button)
        navigation.addStretch(1)
        outer.addLayout(navigation)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(splitter, 1)
        scroll, layout = self._scroll_panel()
        scroll.setMinimumWidth(440)
        splitter.addWidget(scroll)
        self.analysis_intro = QLabel(
            "从一张原始TEM图自动计算真实FFT。单层模式只分析点阵常数、"
            "取向和可能的孪晶域；双层模式同时分析两层点阵、TEM实空间Moiré"
            "周期与相对twist。FFT仅用于交叉验证，最终双层twist以TEM周期为准。")
        self.analysis_intro.setWordWrap(True)
        layout.addWidget(self.analysis_intro)

        input_box = QGroupBox("图像输入与分析设置")
        input_layout = QVBoxLayout(input_box)
        self.bulk_analysis_checkbox = QCheckBox("Bulk analysis")
        self.bulk_analysis_checkbox.setChecked(False)
        self.bulk_analysis_checkbox.setToolTip(
            "开启后可一次选择多个TEM样品，并导出CSV汇总统计。")
        input_layout.addWidget(self.bulk_analysis_checkbox)
        self.image_analysis_mode = QComboBox()
        self.image_analysis_mode.addItem("双层分析", "bilayer")
        self.image_analysis_mode.addItem("单层分析", "single")
        input_layout.addWidget(self.image_analysis_mode)
        tem_row = QHBoxLayout()
        self.select_tem_button = QPushButton("上传TEM图")
        self.tem_path_label = QLabel("尚未选择")
        self.tem_path_label.setWordWrap(True)
        tem_row.addWidget(self.select_tem_button)
        tem_row.addWidget(self.tem_path_label, 1)
        input_layout.addLayout(tem_row)
        self.analysis_selector_widget = QWidget()
        selector_row = QHBoxLayout(self.analysis_selector_widget)
        selector_row.setContentsMargins(0, 0, 0, 0)
        selector_row.addWidget(QLabel("当前预览"))
        self.analysis_file_selector = QComboBox()
        self.analysis_file_selector.setEnabled(False)
        selector_row.addWidget(self.analysis_file_selector, 1)
        input_layout.addWidget(self.analysis_selector_widget)
        self.analysis_selector_widget.setVisible(False)
        overlay_row = QHBoxLayout()
        self.analysis_overlay_label = QLabel("Overlay 显示与导出比例")
        overlay_row.addWidget(self.analysis_overlay_label)
        self.ifft_view_mode = QComboBox()
        self.ifft_view_mode.addItem("Overlay", "overlay")
        self.ifft_view_mode.addItem("Pure inverse FFT", "pure")
        overlay_row.addWidget(self.ifft_view_mode)
        self.ifft_strength = QSlider(Qt.Orientation.Horizontal)
        self.ifft_strength.setRange(0, 100)
        self.ifft_strength.setValue(65)
        self.ifft_strength.setSingleStep(1)
        self.ifft_strength_value = QLabel("65%")
        overlay_row.addWidget(self.ifft_strength, 1)
        overlay_row.addWidget(self.ifft_strength_value)
        input_layout.addLayout(overlay_row)
        self.run_image_analysis_button = QPushButton("3. 自动识别并分析")
        self.run_image_analysis_button.setObjectName("primaryButton")
        input_layout.addWidget(self.run_image_analysis_button)
        layout.addWidget(input_box)

        calibration_box = QGroupBox("当前图像结果与标尺校正")
        calibration = QFormLayout(calibration_box)
        self.scale_bar_pixels = QDoubleSpinBox()
        self.scale_bar_pixels.setRange(0, 100000)
        self.scale_bar_pixels.setDecimals(1)
        self.scale_bar_pixels.setSuffix(" px")
        self.scale_bar_nm = QDoubleSpinBox()
        self.scale_bar_nm.setRange(0, 100000)
        self.scale_bar_nm.setDecimals(2)
        self.scale_bar_nm.setSuffix(" nm")
        self.measured_lattice = QDoubleSpinBox()
        self.measured_lattice.setRange(0, 100000)
        self.measured_lattice.setDecimals(2)
        self.measured_lattice.setSuffix(" nm")
        self.measured_moire = QDoubleSpinBox()
        self.measured_moire.setRange(0, 100000)
        self.measured_moire.setDecimals(2)
        self.measured_moire.setSuffix(" nm")
        self.tem_twist_result = QLabel("—")
        self.fft_twist_result = QLabel("—")
        self.final_twist_result = QLabel("—")
        self.final_twist_result.setObjectName("cardValue")
        calibration.addRow("Scale bar像素长度", self.scale_bar_pixels)
        calibration.addRow("Scale bar标注值", self.scale_bar_nm)
        calibration.addRow("平均 lattice constant", self.measured_lattice)
        calibration.addRow("TEM平均 Moiré period", self.measured_moire)
        calibration.addRow("TEM period计算Twist", self.tem_twist_result)
        calibration.addRow("FFT一阶峰拟合Twist", self.fft_twist_result)
        calibration.addRow("最终采用Twist", self.final_twist_result)
        self.image_analysis_status = QLabel("等待TEM图像。")
        self.image_analysis_status.setWordWrap(True)
        calibration.addRow("识别状态", self.image_analysis_status)
        self.save_image_analysis_button = QPushButton(
            "手动导出当前分析结果（PNG、SVG、JSON）")
        self.save_image_analysis_button.setEnabled(False)
        calibration.addRow("", self.save_image_analysis_button)
        layout.addWidget(calibration_box)
        layout.addStretch(1)

        # Internal compatibility sink. The former manual “添加测量” section
        # is intentionally not shown in the data-analysis workflow.
        self.measurements_table = QTableWidget(0, 5)
        self.measurements_table.setHorizontalHeaderLabels([
            "来源", "Twist angle", "Moiré period", "预测角度", "角度误差"])

        preview = QWidget()
        preview_layout = QVBoxLayout(preview)
        preview_title = QLabel("自动分析结果（所有图像保持原始宽高比）")
        preview_title.setObjectName("title")
        preview_layout.addWidget(preview_title)
        preview_splitter = QSplitter(Qt.Orientation.Horizontal)
        preview_layout.addWidget(preview_splitter, 1)

        def preview_panel(title, placeholder):
            panel = QWidget()
            panel_layout = QVBoxLayout(panel)
            panel_layout.setContentsMargins(2, 2, 2, 2)
            heading = QLabel(title)
            heading.setAlignment(Qt.AlignmentFlag.AlignCenter)
            heading.setObjectName("analysisFigureTitle")
            image_label = QLabel(placeholder)
            image_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
            image_label.setMinimumSize(300, 260)
            image_scroll = QScrollArea()
            image_scroll.setWidgetResizable(True)
            image_scroll.setWidget(image_label)
            panel_layout.addWidget(heading)
            panel_layout.addWidget(image_scroll, 1)
            preview_splitter.addWidget(panel)
            return image_label

        self.tem_analysis_image = preview_panel(
            "Original TEM", "Original TEM将在这里显示")
        self.reconstructed_analysis_image = preview_panel(
            "Reconstructed bilayer", "理论双层重构将在这里显示")
        self.ifft_analysis_image = preview_panel(
            "Selected-spot inverse FFT · phase-matched reconstruction",
            "Selected-spot inverse FFT将在这里显示")
        preview_splitter.setSizes([1, 1, 1])
        self.analysis_preview_note = QLabel(
            "Original TEM右上角使用该TEM的真实FFT；Square用虚线四边形、"
            "Honeycomb/Kagome用六边形标记一阶峰。Selected-spot孔径按真实"
            "峰形拟合并保持复相位。分析完成后不会自动写出SVG；点击左侧"
            "手动导出按钮后，才会按原始像素尺寸和当前Overlay比例生成文件。")
        self.analysis_preview_note.setWordWrap(True)
        preview_layout.addWidget(self.analysis_preview_note)
        splitter.addWidget(preview)
        splitter.setSizes([470, 960])

        self._tem_image_path = None
        self._tem_image_paths = []
        self._fft_image_path = None
        self._tem_analysis = None
        self._fft_analysis = None
        self._analysis_records = []
        self._analysis_output_paths = {}
        self._final_analysis_angle = None
        self._tem_analysis_angle = None
        self._fft_analysis_angle = None
        return tab

    def _apply_style(self):
        self.setStyleSheet("""
        QMainWindow, QWidget { background: #f5f7fa; color: #1c2936; font-family: "Arial", "PingFang SC"; font-size: 13px; }
        QMenuBar { background: white; border-bottom: 1px solid #dce2e9; }
        QMenuBar::item { padding: 3px 8px; }
        QLabel#title { font-size: 24px; font-weight: 700; }
        QGroupBox#previewGroupBox { margin-top: 12px; }
        QGroupBox#previewGroupBox::title {
            subcontrol-origin: margin;
            left: 12px;
            top: 2px;
            padding: 0 5px;
            font-size: 17px;
            font-weight: 750;
        }
        QGroupBox#previewCanvasBox {
            background: #070b10;
            border: 0;
            border-radius: 0;
            margin-top: 0;
            padding: 0;
        }
        QGroupBox#previewCanvasBox[previewPane="top"], QGroupBox#previewCanvasBox[previewPane="side"] { background: #070b10; }
        QLabel#analysisFigureTitle { color: #203142; font-size: 14px; font-weight: 750; padding: 4px; }
        QLabel#analysisColumnHeader { color: #315f88; font-size: 13px; font-weight: 750; }
        QFrame#analysisFigurePanel { background: #0b1014; border: 1px solid #c8d1d9; }
        QFrame#analysisControls { background: #f7f9fa; border: 1px solid #c8d1d9; border-radius: 6px; }
        QLabel#subtitle { color: #66778a; font-size: 13px; }
        QLabel#badge { color: #335f9c; background: #e5effc; border: 1px solid #bad1ee; border-radius: 9px; padding: 5px 10px; font-size: 12px; font-weight: 700; }
        QFrame#workflowBar { background: #e6edf3; border: 1px solid #c2d0dc; border-radius: 6px; }
        QPushButton#workflowButton { color: #23445f; background: #d4e3ef; border: 1px solid #aec3d4; border-radius: 4px; padding: 4px 10px; font-size: 12px; font-weight: 700; text-align: center; }
        QPushButton#workflowButton:hover { color: #17324f; background: #c4dbea; border-color: #8eabc1; }
        QPushButton#workflowButton:checked { color: white; background: #356fab; border-color: #2c629a; font-weight: 750; }
        QPushButton#projectActionButton { color: #234d70; background: #e8f0f7; border: 1px solid #a9c2d8; border-radius: 6px; padding: 5px 12px; font-size: 13px; font-weight: 700; }
        QPushButton#projectActionButton:hover { background: #d9e8f4; border-color: #82a8c5; }
        QPushButton#modeSwitchButton { color: white; background: #6d4c8d; border: 1px solid #5a3e76; border-radius: 4px; padding: 4px 12px; font-size: 12px; font-weight: 750; }
        QPushButton#modeSwitchButton:hover { background: #7d5aa0; border-color: #493260; }
        QPushButton#historyButton { min-width: 32px; max-width: 32px; padding: 3px; background: #f8fafc; color: #23445f; }
        QPushButton#historyButton[historySuccess="true"] { background: #dff4e8; color: #08713f; border-color: #65b88a; font-weight: 800; }
        QPushButton#yieldAnalysisButton { color: #245d43; background: #dcefe5; border: 1px solid #8ebba2; border-radius: 4px; padding: 4px 10px; font-size: 12px; font-weight: 700; }
        QPushButton#yieldAnalysisButton:hover { background: #c9e6d5; border-color: #6ca384; }
        QPushButton#yieldAnalysisButton:checked { color: white; background: #2f855a; border-color: #256d49; font-weight: 750; }
        QPushButton#particleAnalysisButton { color: #744616; background: #f8ead8; border: 1px solid #d8aa73; border-radius: 4px; padding: 4px 10px; font-size: 12px; font-weight: 700; }
        QPushButton#particleAnalysisButton:hover { background: #f3ddc1; border-color: #c28b4c; }
        QPushButton#particleAnalysisButton:checked { color: white; background: #c56a1a; border-color: #a95712; font-weight: 750; }
        QPushButton#twistAnalysisButton { color: #5d3f7c; background: #ebe3f5; border: 1px solid #b8a0d2; border-radius: 4px; padding: 4px 10px; font-size: 12px; font-weight: 700; }
        QPushButton#twistAnalysisButton:hover { background: #dfd1ef; border-color: #9d7dbc; }
        QPushButton#twistAnalysisButton:checked { color: white; background: #7b55a6; border-color: #65458a; font-weight: 750; }
        QLabel#workflowArrow { color: #7890a5; background: transparent; font-size: 14px; font-weight: 700; }
        QLabel#currentProject { color: #214b6d; background: #f4f8fb; border: 1px solid #bfd0de; border-radius: 5px; padding: 4px 9px; font-size: 12px; font-weight: 700; }
        QGroupBox { background: white; border: 1px solid #dce2e9; border-radius: 8px; margin-top: 8px; padding: 8px 6px 5px 6px; font-size: 13px; font-weight: 650; }
        QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 5px; }
        QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox { background: white; border: 1px solid #cfd8e2; border-radius: 5px; padding: 3px 5px; min-height: 20px; font-size: 13px; }
        QPushButton { background: #e9edf2; border: 1px solid #cad3dd; border-radius: 6px; padding: 6px 10px; font-size: 13px; font-weight: 600; }
        QPushButton:hover { background: #dde7f1; }
        QPushButton:disabled { background: #d8dde3; color: #8b949e; border-color: #c8cdd3; }
        QPushButton#primaryButton { background: #356fab; color: white; border-color: #2c629a; }
        QPushButton#primaryButton:disabled { background: #d8dde3; color: #8b949e; border-color: #c8cdd3; }
        QPushButton#parameterStepButton { background: #e8f0f7; color: #234d70; border-color: #a9c2d8; text-align: left; padding: 7px 12px; font-weight: 700; }
        QPushButton#parameterStepButton:hover { background: #d9e8f4; border-color: #82a8c5; }
        QPushButton#parameterStepButton:checked { background: #356fab; color: white; border-color: #2c629a; }
        QPushButton#optionalButton { background: #ece4f6; color: #65458a; border-color: #b8a0d2; font-weight: 700; }
        QPushButton#optionalButton:hover { background: #dfd1ef; border-color: #9d7dbc; }
        QPushButton#optionalButton:checked { background: #7b55a6; color: white; border-color: #65458a; }
        QPushButton#optionalButton:disabled { background: #d8dde3; color: #8b949e; border-color: #c8cdd3; }
        QPushButton#acceptedButton { background: #2f855a; color: white; border-color: #256d49; font-weight: 750; }
        QComboBox#editableParameter, QLabel#editableParameter { color: #126a88; font-weight: 700; }
        QFrame#resultCard { background: white; border: 1px solid #dce2e9; border-radius: 8px; }
        QFrame#validationHeader { background: #edf2f7; border: 1px solid #d4dde7; border-radius: 6px; }
        QLabel#validationTitle { background: transparent; color: #36536f; font-size: 12px; font-weight: 650; }
        QPushButton#validationToggle { background: transparent; color: #315f88; border: 0; padding: 3px 7px; font-size: 12px; }
        QPushButton#validationToggle:hover { background: #dce8f2; }
        QLabel#cardCaption { color: #738396; font-size: 12px; }
        QLabel#cardValue { color: #153b64; font-size: 18px; font-weight: 700; }
        QLabel#fixedValue { color: #315f88; font-size: 12px; font-weight: 650; padding: 2px 0; }
        QLabel#previewParameterBanner { color: #f3f6f9; background: #070b10; border: 1px dashed #52616f; border-radius: 0; padding: 2px 6px; margin: 0; }
        QLabel#previewParameterBanner[previewPane="top"] { background: #070b10; border-right: 0; }
        QLabel#previewParameterBanner[previewPane="side"] { background: #070b10; border-left: 0; }
        QLabel#previewCommonTitle { color: #f3f6f9; background: #070b10; border: 0; padding: 5px 12px; margin: 0; font-size: 19px; font-weight: 750; }
        QFrame#previewSummary { background: #070b10; border: 0; }
        QFrame#previewSummary QLabel { color: #c9d2dc; background: transparent; font-size: 10px; }
        QFrame#previewDashedSeparator { background: transparent; border: 0; border-top: 1px dashed #52616f; }
        QFrame#previewSurface, QFrame#previewMatrix { background: #070b10; border: 0; }
        QFrame#previewViewHeading { border: 0; }
        QFrame#previewViewHeading[previewPane="top"], QFrame#previewViewHeading[previewPane="side"] { background: #070b10; }
        QFrame#previewVerticalSeparator { background: transparent; border: 0; border-left: 1px dashed #52616f; }
        QLabel#previewViewSubtitle { color: #f3f6f9; background: transparent; border: 0; font-size: 11px; font-weight: 750; }
        QFrame#previewLegend { border: 0; border-radius: 0; }
        QFrame#previewLegend[previewPane="top"], QFrame#previewLegend[previewPane="side"] { background: #070b10; }
        QFrame#previewLegend QLabel { color: #f3f6f9; background: transparent; font-size: 11px; }
        QLabel#structureIntro { color: #264f73; background: #eaf3fb; border: 1px solid #c7dceb; border-radius: 7px; padding: 9px; }
        QLabel#structureNote { color: #6d5b2d; background: #fff8dd; border: 1px solid #ead89b; border-radius: 7px; padding: 9px; }
        QLabel#successStatus { color: #17623b; background: #e9f7ef; border: 1px solid #8bc7a6; border-radius: 6px; padding: 7px 9px; font-weight: 650; }
        QTabWidget::pane { border: 0; }
        QTabBar::tab { background: #e7ebf0; padding: 8px 18px; margin-right: 2px; border-radius: 5px; font-size: 13px; }
        QTabBar::tab:selected { background: #356fab; color: white; }
        QTextBrowser, QTableWidget { background: white; border: 1px solid #dce2e9; border-radius: 7px; font-size: 13px; }
        QTextBrowser#phaseHint { background: transparent; border: 0; color: #365b7e; font-size: 12px; }
        """)

    def _set_current_step(self, step):
        if step == 1:
            step = 0
        if 0 <= step < len(self.workflow_buttons):
            self.workflow_buttons[step].setChecked(True)

    def _current_workflow_step(self):
        if self.tabs.currentIndex() == 0:
            return 0
        if self.tabs.currentIndex() in (1, 2):
            return self.tabs.currentIndex() + 1
        return 4 + getattr(self, "analysis_module_stack", self.tabs).currentIndex()

    def _record_history(self, label):
        """Store an in-project undo point; exported files are never deleted."""
        if self._restoring_history or self.project is None or \
                self._app_mode != "design":
            return
        snapshot = {
            "label": str(label),
            "project": self.project.to_dict(),
            "project_path": self.project_path,
            "step": self._current_workflow_step(),
        }
        encoded = json.dumps(snapshot["project"], sort_keys=True)
        if self._history_index >= 0:
            previous = json.dumps(
                self._history[self._history_index]["project"], sort_keys=True)
            if encoded == previous and snapshot["step"] == \
                    self._history[self._history_index]["step"]:
                return
        self._history = self._history[:self._history_index + 1]
        self._history.append(snapshot)
        self._history = self._history[-80:]
        self._history_index = len(self._history) - 1
        self._update_history_buttons()

    def _update_history_buttons(self):
        self.history_back_button.setEnabled(self._history_index > 0)
        self.history_forward_button.setEnabled(
            0 <= self._history_index < len(self._history) - 1)

    def _restore_history_snapshot(self, index, action=None):
        if not (0 <= index < len(self._history)):
            return False
        snapshot = self._history[index]
        self._history_restore_token += 1
        restore_token = self._history_restore_token
        self._restoring_history = True
        try:
            self._sequence_analysis = None
            self.project = MoireProject.from_dict(snapshot["project"])
            self.project_path = snapshot.get("project_path")
            self._load_settings(self.project.settings)
            restored = MoireProject.from_dict(snapshot["project"])
            self.project.seed_plan = restored.seed_plan
            self.project.capture_plan = restored.capture_plan
            self.project.measurements = restored.measurements
            self._restore_structure_workflow()
            self._update_current_project_display()
            self._go_to_step(snapshot.get("step", 0))
            self._history_index = index
            self._update_history_buttons()
            label = snapshot.get("label", "设计状态")
            if action == "undo":
                message = "撤销成功：已恢复到“%s”。" % label
            elif action == "redo":
                message = "重做成功：已恢复到“%s”。" % label
            else:
                message = "已恢复：%s（已导出的磁盘文件不会被删除）" % label
            self.statusBar().showMessage(message, 5000)
        finally:
            # Some combo-box/model updates post signals to the Qt event queue.
            # Releasing the guard immediately lets those delayed signals call
            # _record_history and truncate the redo branch created by Undo.
            # Keep the guard through the current event-loop turn; the token
            # prevents an older callback from ending a newer restore.
            QTimer.singleShot(
                0, lambda token=restore_token:
                self._finish_history_restore(token))
        return True

    def _finish_history_restore(self, token):
        if token == self._history_restore_token:
            self._restoring_history = False

    def _flash_history_success(self, button, default_text):
        """Make a successful Undo/Redo click visible next to the cursor."""
        button.setProperty("historySuccess", True)
        button.setText("✓")
        button.style().unpolish(button)
        button.style().polish(button)

        def restore_button():
            button.setProperty("historySuccess", False)
            button.setText(default_text)
            button.style().unpolish(button)
            button.style().polish(button)

        QTimer.singleShot(900, restore_button)

    def _history_back(self):
        if self._restore_history_snapshot(
                self._history_index - 1, action="undo"):
            self._flash_history_success(self.history_back_button, "←")

    def _history_forward(self):
        if self._restore_history_snapshot(
                self._history_index + 1, action="redo"):
            self._flash_history_success(self.history_forward_button, "→")

    def _toggle_validation(self, collapsed):
        geometry = self.geometry()
        state = self.windowState()
        if collapsed:
            self.validation.setMinimumHeight(0)
            self.validation.setMaximumHeight(0)
        else:
            self.validation.setMinimumHeight(82)
            self.validation.setMaximumHeight(118)
        self.validation_toggle.setText(
            "展开验证说明  ▴" if collapsed else "收起验证说明  ▾")
        self.validation_toggle.setToolTip(
            "显示设计验证说明" if collapsed else "最小化设计验证说明")
        # Toggling only redistributes space inside the existing window.  It
        # must never resize, un-maximize, or extend the page beyond the screen.
        if state & Qt.WindowState.WindowMaximized:
            self.setWindowState(state)
        else:
            self.setGeometry(geometry)

    def _go_to_step(self, step):
        step = int(step)
        # History from the former split 1.1/1.2 workflow may still reference
        # step 1.  It now resolves to the unified Design and Prediction page.
        if step == 1:
            step = 0
        if self._app_mode == "analysis" and step < 4:
            return
        if self._app_mode == "design" and step >= 4:
            return
        if step == 0:
            self.tabs.setCurrentIndex(0)
            self.design_stack.setCurrentIndex(0)
        elif step == 2:
            self.tabs.setCurrentIndex(1)
        elif step == 3:
            self.tabs.setCurrentIndex(2)
        elif 4 <= step <= 5:
            self.tabs.setCurrentIndex(3)
            if hasattr(self, "analysis_module_stack"):
                self.analysis_module_stack.setCurrentIndex(step - 4)
        else:
            return
        self._set_current_step(step)

    def _tab_changed(self, index):
        if index == 3 and hasattr(self, "analysis_module_stack"):
            self._set_current_step(
                4 + self.analysis_module_stack.currentIndex())
        elif index == 0:
            if self.design_stack.currentIndex() != 0:
                self.design_stack.setCurrentIndex(0)
            self._set_current_step(0)
        elif index in (1, 2):
            self._set_current_step(index + 1)
        else:
            self._set_current_step(index)

    def _connect_signals(self):
        self.new_action.triggered.connect(self.new_project)
        self.open_action.triggered.connect(self.open_project)
        self.new_project_button.clicked.connect(self.new_project)
        self.open_project_button.clicked.connect(self.open_project)
        self.save_action.triggered.connect(self.save_project)
        self.save_as_action.triggered.connect(self.save_project_as)
        self.export_action.triggered.connect(self.export_project)
        self.cadnano_action.triggered.connect(self.open_in_cadnano)
        self.paper_preset_menu_action.triggered.connect(
            self.apply_paper_preset)
        self.view_design_action.triggered.connect(lambda: self._go_to_step(0))
        self.view_capture_action.triggered.connect(lambda: self._go_to_step(2))
        self.view_sequence_action.triggered.connect(lambda: self._go_to_step(3))
        self.analysis_crystal_action.triggered.connect(
            lambda: self._open_analysis_module(0))
        self.bilayer_symmetry_selector.currentIndexChanged.connect(
            lambda unused_index: self._apply_symmetry_ui())
        self.seed_cross_section_picker.selectionChanged.connect(
            self._basis_selection_changed)
        self.seed_cross_section_preset.currentIndexChanged.connect(
            self._cross_section_preset_changed)
        self.reset_seed_cross_section_button.clicked.connect(
            self.seed_cross_section_picker.reset_default)
        self.confirm_design_basis_button.clicked.connect(
            self._accept_design_basis)
        self.design_basis_next_button.clicked.connect(
            lambda: self._go_to_step(1))
        self.change_design_basis_button.clicked.connect(
            lambda: self.design_stack.setCurrentIndex(0))
        self.accept_parameters_button.clicked.connect(self.accept_parameters)
        self.parameters_next_button.clicked.connect(
            lambda: self._go_to_step(2))
        self.history_back_button.clicked.connect(self._history_back)
        self.history_forward_button.clicked.connect(self._history_forward)
        self.mode_switch_button.clicked.connect(self._switch_app_mode)
        self.go_to_analysis_button.clicked.connect(
            lambda: self._open_analysis_module(
                self.analysis_module_stack.currentIndex()))
        self.select_tem_button.clicked.connect(self.select_tem_image)
        self.run_image_analysis_button.clicked.connect(self.run_image_analysis)
        self.save_image_analysis_button.clicked.connect(
            self.save_image_analysis)
        self.image_analysis_mode.currentIndexChanged.connect(
            self._analysis_mode_changed)
        self.bulk_analysis_checkbox.toggled.connect(
            self._bulk_analysis_changed)
        self.bulk_scale_mode.currentIndexChanged.connect(
            self._bulk_scale_mode_changed)
        self.analysis_file_selector.currentIndexChanged.connect(
            self._activate_analysis_record)
        self.ifft_strength.valueChanged.connect(
            self._analysis_overlay_changed)
        self.ifft_view_mode.currentIndexChanged.connect(
            self._analysis_overlay_changed)
        self.bulk_ifft_strength.valueChanged.connect(
            self._analysis_overlay_changed)
        self.bulk_ifft_view_mode.currentIndexChanged.connect(
            self._analysis_overlay_changed)
        self.scale_bar_pixels.valueChanged.connect(
            self._scale_calibration_changed)
        self.scale_bar_nm.valueChanged.connect(
            self._scale_calibration_changed)
        self.measured_lattice.valueChanged.connect(
            self._analysis_values_changed)
        self.measured_moire.valueChanged.connect(
            self._analysis_values_changed)
        self.back_to_parameters_button.clicked.connect(
            self.return_to_parameters)
        self.back_to_scaffold_button.clicked.connect(
            self.return_to_scaffold)
        self.back_to_structure_button.clicked.connect(
            self.return_to_structure)
        self.import_capture_project_button.clicked.connect(
            self.import_capture_moire_project)
        self.import_sequence_project_button.clicked.connect(
            self.import_sequence_moire_project)
        self.structure_preview_channel.currentIndexChanged.connect(
            self._structure_preview_channel_changed)
        self.open_orthogonal_sequences_button.clicked.connect(
            self.open_orthogonal_sequence_designer)
        self.generate_scaffold_button.clicked.connect(
            self.generate_scaffold_design)
        self.generate_simple_design_button.clicked.connect(
            self.generate_simple_structure_design)
        self.structure_expert_button.toggled.connect(
            self._toggle_structure_expert)
        self.open_scaffold_button.clicked.connect(
            lambda: self._open_structure_file("scaffold_review"))
        self.load_scaffold_button.clicked.connect(self.load_expert_scaffold)
        self.accept_scaffold_button.clicked.connect(self.accept_scaffold)
        self.generate_structure_button.clicked.connect(
            self.generate_complete_structure)
        self.open_structure_button.clicked.connect(
            lambda: self._open_structure_file("structure_complete"))
        self.inspect_final_design_button.clicked.connect(
            self.open_cadnano_for_optional_editing)
        self.accept_structure_button.clicked.connect(
            self.accept_complete_structure)
        self.structure_next_button.clicked.connect(
            lambda: self._go_to_step(3))
        self.load_sequence_design_button.clicked.connect(
            self.detect_sequence_scaffolds)
        self.accept_added_scaffold_button.clicked.connect(
            self.accept_added_scaffolds)
        self.detect_sst_inputs_button.clicked.connect(
            self.detect_sequence_sst_inputs)
        self.auto_design_sst_inputs_button.clicked.connect(
            self.auto_design_and_add_sst_inputs)
        self.sequence_expert_button.toggled.connect(
            self._toggle_sequence_expert)
        self.export_sst_input_template_button.clicked.connect(
            self.export_sequence_sst_template)
        self.import_sst_input_template_button.clicked.connect(
            self.import_sequence_sst_template)
        self.accept_added_sst_button.clicked.connect(
            self.accept_added_sst_inputs)
        self.export_sequence_variants_button.clicked.connect(
            self.export_sequence_final_package)
        self.lattice_context.currentIndexChanged.connect(self._context_changed)
        self.layers_identical.currentIndexChanged.connect(
            self._layers_identical_changed)
        self.target_definition.currentIndexChanged.connect(
            lambda: self._control_changed("indel", 0))
        self.angle.valueChanged.connect(
            lambda: self._control_changed("angle", 0))
        self.angle.valueChanged.connect(
            lambda unused_value: self._basis_selection_changed())
        self.period.valueChanged.connect(
            lambda: self._control_changed("period", 0))
        self.lattice_constant.valueChanged.connect(
            lambda: self._control_changed("angle", 0))
        self.lattice_constant_2.valueChanged.connect(
            lambda: self._control_changed("angle", 0))
        for name, widget in (
                ("sst_z1", self.sst_z1),
                ("sst_spacing", self.sst_spacing),
                ("sst_z3", self.sst_z3),
                ("seed_z2", self.seed_z2)):
            widget.aboutToOpen.connect(
                lambda n=name: self._prepare_phase_anchor(n))
            widget.currentIndexChanged.connect(
                lambda unused, n=name: self._phase_changed(n))
        self._populate_cross_section_presets()

        self._basis_selection_changed()

    def _selected_symmetry(self):
        return str(self.bilayer_symmetry_selector.currentData() or
                   "square_square_c4")

    @staticmethod
    def _seed_preset_cells(key):
        if key == "s6_r2x2":
            return [[row, col] for row in range(1, 7) for col in range(1, 7)
                    if not (3 <= row <= 4 and 3 <= col <= 4)]
        if key == "s6x5_r2x1":
            return [[row, col] for row in range(1, 7) for col in range(1, 6)
                    if not (3 <= row <= 4 and col == 3)]
        return [[row, col] for row in range(8) for col in range(8)
                if not (2 <= row <= 5 and 2 <= col <= 5)]

    @classmethod
    def _preset_key_for_cells(cls, cells):
        selected = {tuple(map(int, cell)) for cell in cells}
        for key in ("s8_r4x4", "s6_r2x2", "s6x5_r2x1"):
            if selected == {tuple(item) for item in cls._seed_preset_cells(key)}:
                return key
        return None

    def _populate_cross_section_presets(self):
        if not hasattr(self, "seed_cross_section_preset"):
            return
        symmetry = self._selected_symmetry()
        current = self.seed_cross_section_preset.currentData()
        options = [("8×8 + 4×4 pore", "s8_r4x4")]
        with QSignalBlocker(self.seed_cross_section_preset):
            self.seed_cross_section_preset.clear()
            option_colors = (QColor("#2a78d1"),)
            for option_index, (label, key) in enumerate(options):
                self.seed_cross_section_preset.addItem(label, key)
                self.seed_cross_section_preset.setItemData(
                    option_index, option_colors[option_index],
                    Qt.ItemDataRole.ForegroundRole)
            index = self.seed_cross_section_preset.findData(current)
            self.seed_cross_section_preset.setCurrentIndex(
                index if index >= 0 else 0)
        if hasattr(self, "seed_cross_section_preset_display"):
            self.seed_cross_section_preset_display.setText(
                self.seed_cross_section_preset.currentText())
        self._cross_section_preset_changed()

    def _cross_section_preset_changed(self):
        if self._updating or not hasattr(self, "seed_cross_section_preset"):
            return
        key = str(self.seed_cross_section_preset.currentData() or
                  "s8_r4x4")
        self.seed_cross_section_picker.set_cells(
            self._seed_preset_cells(key))
        self._basis_selection_changed()
        if self.project is not None:
            self.recalculate()
            self._record_history("更改 Seed 截面预设")

    @staticmethod
    def _set_acceptance_button(button, accepted, pending_text,
                               accepted_text="已接受"):
        button.setObjectName("acceptedButton" if accepted else "primaryButton")
        button.setText(accepted_text if accepted else pending_text)
        button.style().unpolish(button)
        button.style().polish(button)

    def _invalidate_design_basis_acceptance(self):
        if self._updating or self.project is None:
            return
        workflow = self._workflow()
        self._invalidate_downstream_design_state(workflow)
        workflow["design_basis_accepted"] = False
        workflow["parameters_accepted"] = False
        workflow.pop("design_basis_accepted_at", None)
        workflow.pop("parameters_accepted_at", None)
        workflow["parameters_editing"] = True
        self.design_basis_next_button.setEnabled(False)
        self.parameters_next_button.setEnabled(False)
        self._set_acceptance_button(
            self.confirm_design_basis_button, False,
            "4. 接受当前对称性与 Twist")
        self._set_acceptance_button(
            self.accept_parameters_button, False,
            "接受当前 Moiré 参数")
        if hasattr(self, "design_basis_action_status"):
            self.design_basis_action_status.hide()
        if hasattr(self, "parameters_action_status"):
            self.parameters_action_status.hide()

    def _invalidate_parameter_acceptance(self):
        if self._updating or self.project is None:
            return
        workflow = self._workflow()
        self._invalidate_downstream_design_state(workflow)
        workflow["parameters_accepted"] = False
        workflow.pop("parameters_accepted_at", None)
        workflow["parameters_editing"] = True
        self.parameters_next_button.setEnabled(False)
        self._set_acceptance_button(
            self.accept_parameters_button, False,
            "接受当前 Moiré 参数")
        workflow.pop("scaffold_capacity_precheck", None)
        if hasattr(self, "seed_scaffold_capacity_status"):
            self.seed_scaffold_capacity_status.setText(
                "参数已改变；接受时重新精确核算")
        if hasattr(self, "parameters_action_status"):
            self.parameters_action_status.hide()

    @staticmethod
    def _symmetry_label(symmetry):
        return {
            "square_square_c4": "Square–Square",
            "kagome_kagome": "Kagome–Kagome",
            "square_kagome": "Square–Kagome",
        }.get(str(symmetry), str(symmetry))

    def _context_lattice_constants(self, symmetry=None):
        symmetry = symmetry or self._selected_symmetry()
        data = self.lattice_context.currentData()
        if not data:
            context, square_a, kagome_a = "solution_cryo", 2.8, 5.4
        else:
            context, square_a, kagome_a = data
        if context == "custom":
            first = self.lattice_constant.value()
            return ((first, self.lattice_constant_2.value())
                    if symmetry == "square_kagome" else (first, first))
        if symmetry == "square_square_c4":
            return square_a, square_a
        if symmetry == "kagome_kagome":
            return kagome_a, kagome_a
        return square_a, kagome_a

    def _basis_selection_changed(self):
        self._invalidate_design_basis_acceptance()
        symmetry = self._selected_symmetry()
        cells = self.seed_cross_section_picker.cells()
        count = len(cells)
        notes = {
            "square_square_c4": (
                "两层均为Square；可由Twist与Square a计算Moiré period。"),
            "kagome_kagome": (
                "两层均为Kagome；helix间距仍为2.8 nm，a使用"
                "cryo-EM 5.4 nm或干燥TEM 4.4 nm。"),
            "square_kagome": (
                "1st layer为Square，2nd layer为Kagome；可设置Twist，"
                "但不同点阵之间不定义Moiré period。"),
        }
        self.bilayer_symmetry_note.setText(notes[symmetry])
        preset = self.seed_cross_section_preset.currentText()
        self.seed_cross_section_status.setText(
            "当前预设：%s，%d根helix。%s" % (
                preset, count,
                ("可进入1.2 Seed/SST参数。" if count >= 4 else
                 "至少需要4根helix。")))
        self.confirm_design_basis_button.setEnabled(count >= 4)
        self.accept_parameters_button.setEnabled(count >= 4)
        self.setup_preview.set_configuration(
            symmetry, cells,
            self._context_lattice_constants(symmetry),
            self.angle.value(),
            None if symmetry == "square_kagome" else self.period.value())

    def _apply_symmetry_ui(self):
        symmetry = self._selected_symmetry()
        label = self._symmetry_label(symmetry)
        self.lattice_symmetry.setText(label)
        if hasattr(self, "design_preview_title"):
            self.design_preview_title.setText(
                "%s Bilayer DNA Moiré Superlattice" % label)
        mixed = symmetry == "square_kagome"
        kagome = symmetry == "kagome_kagome"
        zero_spacing = (
            hasattr(self, "sst_spacing") and
            int(self.sst_spacing.currentData() or 0) == 0)
        self.angle.setEnabled(not zero_spacing)
        self.period.setEnabled(not mixed and not zero_spacing)
        self.period.setVisible(not mixed)
        self.period_label.setVisible(not mixed)
        self.period_card.hide()
        self.period.setToolTip(
            "Square–Kagome由两种不同点阵组成，不定义单一Moiré period。"
            if mixed else "由Twist和当前点阵a计算。")
        # Lattice constants are measured-context constants rather than user
        # parameters.  Keep the former spin boxes hidden for compatibility
        # and expose only the fixed readouts (without adjustment arrows).
        self.lattice_constant.setVisible(False)
        self.lattice_constant_2.setVisible(False)
        # Always show the two physical lattice constants separately.  Even
        # identical bilayers retain one readout per layer so the presentation
        # is consistent with the mixed-lattice mode.
        self.lattice_constant_2_label.setVisible(True)
        self.lattice_constant_2_fixed.setVisible(True)
        self._populate_cross_section_presets()
        if mixed:
            with QSignalBlocker(self.layers_identical):
                different_index = self.layers_identical.findData(False)
                self.layers_identical.setCurrentIndex(
                    max(0, different_index))
            self.layers_identical.setEnabled(False)
            self.lattice_constant_label.setText("1st layer a (square)")
            self.lattice_constant_2_label.setText("2nd layer a (Kagome)")
            self.sst_z1_label.setText("1st layer length（Square）")
            self.sst_z3_label.setText("2nd layer length（Kagome）")
            self._target_driver = "angle"
        else:
            self.layers_identical.setEnabled(True)
            lattice_name = "Kagome" if kagome else "square"
            self.lattice_constant_label.setText(
                "1st layer a (%s)" % lattice_name)
            self.lattice_constant_2_label.setText(
                "2nd layer a (%s)" % lattice_name)
            self.sst_z1_label.setText("1st layer length")
            self.sst_z3_label.setText("2nd layer length")
        self._context_changed()
        self.accept_parameters_button.setEnabled(
            len(self.seed_cross_section_picker.cells()) >= 4)

    def _accept_design_basis(self):
        if len(self.seed_cross_section_picker.cells()) < 4:
            QMessageBox.warning(
                self, "Seed截面不足", "请至少选择4根Square网格helix。")
            return
        self._apply_symmetry_ui()
        self.recalculate()
        workflow = self._workflow()
        workflow["design_basis_accepted"] = True
        workflow["design_basis_accepted_at"] = datetime.now().isoformat()
        workflow["parameters_accepted"] = False
        self.design_basis_next_button.setEnabled(True)
        self.accept_parameters_button.setEnabled(True)
        self._set_acceptance_button(
            self.confirm_design_basis_button, True,
            "4. 接受当前对称性与 Twist",
            "✓ 1.1 已接受")
        self._set_acceptance_button(
            self.accept_parameters_button, False,
            "接受当前 Moiré 参数")
        self.design_basis_action_status.hide()
        self._record_history("接受 1.1 对称性与 Twist 参数")

    def _confirm_design_basis(self):
        """Compatibility alias retained for older signal wiring."""
        self._accept_design_basis()

    @staticmethod
    def _all_growth_values():
        return list(range(64, 401, 8))

    @staticmethod
    def _all_z2_values():
        return list(range(0, 161, 8))

    @staticmethod
    def _nearest(values, current):
        return min(values, key=lambda value: (abs(value-current), value))

    def _fill_length_combo(self, combo, values, current):
        values = list(values)
        if not values:
            return
        selected = self._nearest(values, int(current))
        with QSignalBlocker(combo):
            combo.clear()
            for value in values:
                combo.addItem("%d bp" % value, value)
            combo.setCurrentIndex(combo.findData(selected))

    def _prepare_phase_anchor(self, name):
        """Let the opened control choose any 8-bp domain phase."""
        if self._phase_sync:
            return
        combo = getattr(self, name)
        current = combo.currentData()
        values = (self._all_z2_values()
                  if name in ("sst_spacing", "seed_z2") else
                  self._all_growth_values())
        self._fill_length_combo(combo, values, current)

    def _phase_changed(self, source):
        if self._phase_sync or self._updating:
            return
        combo = getattr(self, source)
        value = combo.currentData()
        if value is None:
            return
        self._phase_sync = True
        try:
            if source in ("sst_spacing", "seed_z2"):
                self._fill_length_combo(
                    self.sst_spacing, self._all_z2_values(), value)
                self._fill_length_combo(
                    self.seed_z2, self._all_z2_values(), value)
            if not bool(self.layers_identical.currentData()):
                # Independent layers preserve every selected length.  The
                # combos still expose only integral 8-bp steps.
                pass
            elif source in ("sst_spacing", "seed_z2"):
                compatible = compatible_growth_values(value, maximum=400)
                self._fill_length_combo(
                    self.sst_z1, compatible,
                    self.sst_z1.currentData() or 128)
                self._fill_length_combo(
                    self.sst_z3, compatible,
                    self.sst_z3.currentData() or 128)
            else:
                compatible_z2 = [item for item in
                                 compatible_z2_values(value, maximum=160)
                                 if item <= 160]
                self._fill_length_combo(
                    self.sst_spacing, compatible_z2,
                    self.sst_spacing.currentData() or 0)
                self._fill_length_combo(
                    self.seed_z2, compatible_z2,
                    self.sst_spacing.currentData() or 0)
                other = (self.sst_z3 if source == "sst_z1"
                         else self.sst_z1)
                self._fill_length_combo(
                    other, compatible_growth_values(
                        self.sst_spacing.currentData(), maximum=400),
                    value)
        finally:
            self._phase_sync = False
        self.seed_z2_readout.setText(
            "%d bp（由 Layer spacing 强制联动）" %
            int(self.seed_z2.currentData() or 0))
        self._update_phase_hint()
        # Length is independent of the target twist.  Recalculate only the
        # indel needed to preserve the currently displayed angle.
        self._control_changed("angle", 1)

    def _layers_identical_changed(self):
        if self._phase_sync or self._updating:
            return
        self._phase_sync = True
        try:
            z1 = int(self.sst_z1.currentData() or 128)
            z2 = int(self.sst_spacing.currentData() or 0)
            z3 = int(self.sst_z3.currentData() or 128)
            self._fill_length_combo(
                self.sst_z1, self._all_growth_values(), z1)
            if bool(self.layers_identical.currentData()):
                self._fill_length_combo(
                    self.sst_spacing,
                    [item for item in compatible_z2_values(
                        z1, maximum=160) if item <= 160], z2)
                self._fill_length_combo(
                    self.seed_z2,
                    [item for item in compatible_z2_values(
                        z1, maximum=160) if item <= 160],
                    self.sst_spacing.currentData())
                self._fill_length_combo(
                    self.sst_z3,
                    compatible_growth_values(int(
                        self.sst_spacing.currentData() or 0), maximum=400),
                    z1)
            else:
                self._fill_length_combo(
                    self.sst_spacing, self._all_z2_values(), z2)
                self._fill_length_combo(
                    self.seed_z2, self._all_z2_values(), z2)
                self._fill_length_combo(
                    self.sst_z3, self._all_growth_values(), z3)
        finally:
            self._phase_sync = False
        self.seed_z2_readout.setText(
            "%d bp（由 Layer spacing 强制联动）" %
            int(self.seed_z2.currentData() or 0))
        self._update_phase_hint()
        self._control_changed("angle", 1)

    def _update_phase_hint(self):
        if not bool(self.layers_identical.currentData()):
            self.phase_hint.setText(
                "两层长度与spacing独立，步长均为8 bp。\n"
                "不强制32 bp相位联动。")
            return
        z1 = int(self.sst_z1.currentData() or 0)
        z2_values = [item for item in compatible_z2_values(
            z1, maximum=160) if item <= 160]
        growth_residue = z1 % 32
        self.phase_hint.setText(
            "SST superlattice两层长度相同：32n%+d；"
            "Layer spacing/Seed spacing：32n%+d\n"
            "spacing：%d–%d bp，步长32 bp" %
            (growth_residue, (-growth_residue) % 32,
             z2_values[0], z2_values[-1]))

    def apply_paper_preset(self):
        self._updating = True
        blockers = [QSignalBlocker(widget) for widget in (
            self.target_definition, self.angle, self.period,
            self.lattice_context, self.lattice_constant,
            self.lattice_constant_2, self.mean_indel,
            self.layers_identical)]
        try:
            self.target_definition.setCurrentIndex(0)
            self.angle.setValue(3.2967555036483183)
            self.period.setValue(48.669158335514105)
            self.lattice_context.setCurrentIndex(0)
            self.lattice_constant.setValue(2.8)
            self.lattice_constant_2.setValue(2.8)
            self.mean_indel.setValue(0.0)
            self.layers_identical.setCurrentIndex(0)
            self._fill_length_combo(
                self.sst_z1, self._all_growth_values(), 128)
            self._fill_length_combo(
                self.sst_spacing,
                [item for item in compatible_z2_values(
                    128, maximum=160) if item <= 160], 32)
            self._fill_length_combo(
                self.seed_z2,
                [item for item in compatible_z2_values(
                    128, maximum=160) if item <= 160], 32)
            self._fill_length_combo(
                self.sst_z3,
                compatible_growth_values(32, maximum=400), 128)
            self.seed_z1.setValue(128)
            self.seed_z3.setValue(128)
        finally:
            del blockers
            self._updating = False
        self._target_driver = "indel"
        self._update_phase_hint()
        self._apply_symmetry_ui()
        self._go_to_step(0)
        self.recalculate()

    def _project_directory(self):
        if not self.project_path:
            return None
        return Path(self.project_path).expanduser().resolve().parent

    def _project_output_dir(self, category=None):
        # Analysis exports are deliberately independent from a Design task.
        # Returning None makes each analysis module ask for its own export
        # destination instead of writing beside a loaded .moire.json file.
        if self._app_mode == "analysis":
            return None
        root = self._project_directory()
        if root is None:
            return None
        output = root if not category else root/Path(category)
        output.mkdir(parents=True, exist_ok=True)
        return output

    def _update_current_project_display(self):
        if self._app_mode == "analysis":
            self.setWindowTitle("DNA Moiré Designer — Analysis Mode")
            return
        name = self.project_name.text().strip() or "Untitled"
        self.setWindowTitle(
            "DNA Moiré Designer — Current Project: %s" % name)

    def _project_setup_selection(self, title, accept_text,
                                 choose_language=False):
        project_directory = self._project_directory()
        directory = (project_directory.parent if project_directory else
                     (Path.home()/"Desktop"))
        dialog = ProjectSetupDialog(
            self, title=title,
            project_name=self.project_name.text() or "moire_project",
            directory=directory, accept_text=accept_text,
            show_language=choose_language)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return None
        return dialog.selection()

    def _prompt_startup(self):
        """Compatibility entry point: startup now always opens Design mode."""
        self._apply_app_mode("design")
        self._startup_complete = True

    def _prompt_design_entry(self):
        dialog = QMessageBox(self)
        dialog.setWindowTitle(translate("Design · 设计"))
        dialog.setText(translate("新建项目或打开已有 Moiré 项目。"))
        new_button = dialog.addButton(
            translate("新建项目"), QMessageBox.ButtonRole.AcceptRole)
        open_button = dialog.addButton(
            translate("打开项目"), QMessageBox.ButtonRole.ActionRole)
        dialog.addButton(QMessageBox.StandardButton.Cancel)
        dialog.exec()
        if dialog.clickedButton() is new_button:
            if not self.new_project():
                self.close()
                return
        elif dialog.clickedButton() is open_button:
            if not self.open_project():
                self.close()
                return
        else:
            self.close()
            return
        self._startup_complete = True
        self._record_history("打开设计项目")

    def _apply_app_mode(self, mode):
        self._app_mode = str(mode)
        design = self._app_mode == "design"
        for button in self.workflow_buttons[:4]:
            button.setVisible(design)
        for button in self.workflow_buttons[4:]:
            button.setVisible(not design)
        for separator in self.workflow_design_separators:
            separator.setVisible(design)
        # The old split 1.2 slot remains compatibility-only in Design mode.
        self.workflow_buttons[1].setVisible(False)
        self.workflow_design_separators[0].setVisible(False)
        self.history_back_button.hide()
        self.history_forward_button.hide()
        self.new_action.setVisible(design)
        self.open_action.setVisible(design)
        self.view_design_action.setVisible(design)
        self.view_capture_action.setVisible(design)
        self.view_sequence_action.setVisible(design)
        self.analysis_crystal_action.setVisible(not design)
        self.save_action.setVisible(design)
        self.save_as_action.setVisible(design)
        self.export_action.setVisible(design)
        self.cadnano_action.setVisible(design)
        self.paper_preset_menu_action.setVisible(design)
        if design:
            self.mode_switch_button.setText("Switch to Analysis Mode")
            self._go_to_step(self._last_design_step)
            self._update_current_project_display()
        else:
            self.mode_switch_button.setText("Switch to Design Mode")
            self.tabs.setCurrentIndex(3)
            self.analysis_module_stack.setCurrentIndex(
                self._last_analysis_module)
            self._set_current_step(4 + self._last_analysis_module)
            self._update_current_project_display()

    def _switch_app_mode(self):
        """Switch workspaces without replacing or saving Design state."""
        if self._app_mode == "analysis":
            self._last_analysis_module = \
                self.analysis_module_stack.currentIndex()
            self._apply_app_mode("design")
            self.statusBar().showMessage(
                "Design Mode restored. Continue the existing design or "
                "open a project.", 5000)
            return

        current_step = self._current_workflow_step()
        if current_step < 4:
            self._last_design_step = current_step
        self._apply_app_mode("analysis")
        self.statusBar().showMessage(
            "Analysis Mode. No project file is required or saved.", 5000)

    def _open_analysis_module(self, index):
        """Open Moiré analysis after the four-step design workflow."""
        index = 0
        self._last_analysis_module = index
        self.analysis_module_stack.setCurrentIndex(index)
        self._go_to_step(4 + index)

    def new_project(self, unused_checked=False):
        selection = self._project_setup_selection(
            "新建 DNA Moiré 项目", "创建项目", choose_language=False)
        if not selection:
            return False
        name, filename, unused_language = selection
        self.project = None
        self.project_path = filename
        self.structure_root = None
        self._sequence_analysis = None
        self._sequence_assignments = {}
        self.project_name.setText(name)
        self.bilayer_symmetry_selector.setCurrentIndex(0)
        self.seed_cross_section_picker.reset_default()
        self.apply_paper_preset()
        self.design_stack.setCurrentIndex(0)
        self._refresh_structure_preview()
        self.recalculate()
        self._save_current_project(silent=True)
        self._update_current_project_display()
        self._record_history("新建项目")
        return True

    def _ensure_project_for_parameter_acceptance(self):
        """Create a save target without resetting the edited parameters."""
        if self.project_path:
            return True
        selection = self._project_setup_selection(
            "Create DNA Moiré Project", "Create Project",
            choose_language=False)
        if not selection:
            return False
        name, filename, unused_language = selection
        self.project_name.setText(name)
        self.project_path = filename
        self.structure_root = None
        self.recalculate()
        if self.project is None or not self._save_current_project(silent=True):
            self.project_path = None
            QMessageBox.critical(
                self, "Project could not be created",
                "The selected project file could not be saved. Please "
                "choose another name or location.")
            return False
        self._update_current_project_display()
        self.statusBar().showMessage(
            "Project created: %s" % self.project_path, 6000)
        return True

    def accept_parameters(self):
        if not self._ensure_project_for_parameter_acceptance():
            return
        if self.project is None:
            self.recalculate()
        workflow = self._workflow()
        # The unified page has one acceptance action.  Preserve the former
        # design_basis fields in saved projects for backwards compatibility,
        # but commit them atomically with the complete Moiré parameter set.
        if not workflow.get("design_basis_accepted"):
            if len(self.seed_cross_section_picker.cells()) < 4:
                QMessageBox.warning(
                    self, "Seed截面不足",
                    "请至少选择4根Square网格helix。")
                return
            self._apply_symmetry_ui()
            self.recalculate()
            workflow = self._workflow()
            workflow["design_basis_accepted"] = True
            workflow["design_basis_accepted_at"] = \
                datetime.now().isoformat()
        if self.project.prediction.get("seed_deletion_limit_exceeded"):
            required = self.project.settings.mean_indel_per_helix
            minimum = self.project.prediction[
                "minimum_seed_deletion_per_helix"]
            spacing = self.project.settings.spacer_bp_z2
            QMessageBox.warning(
                self, "Seed deletion exceeds the spacing-dependent limit",
                "The current design requires %.1f bases/helix, below the "
                "%.1f limit for %d-bp spacing. Each 8-bp domain permits at "
                "most %d evenly distributed deletions. Reduce the Twist "
                "magnitude or increase the spacing." % (
                    required, minimum, spacing,
                    MAX_SEED_DELETIONS_PER_DOMAIN))
            return
        if self.project.prediction.get("seed_insertion_limit_exceeded"):
            required = self.project.settings.mean_indel_per_helix
            maximum = self.project.prediction[
                "maximum_seed_insertion_per_helix"]
            spacing = self.project.settings.spacer_bp_z2
            QMessageBox.warning(
                self, "Seed insertion exceeds the spacing-dependent limit",
                "The current design requires %.1f bases/helix, above the "
                "+%.1f limit for %d-bp spacing. Each 8-bp domain permits "
                "at most %d evenly distributed insertions, with a global "
                "+%.1f cap. Reduce the Twist magnitude or increase the "
                "spacing." % (
                    required, maximum, spacing,
                    MAX_SEED_DELETIONS_PER_DOMAIN,
                    MAX_SEED_INSERTION_PER_HELIX))
            return
        self._busy(True)
        try:
            cfg = self.project.settings
            with tempfile.TemporaryDirectory(
                    prefix="moire_capacity_") as folder:
                preview_sst = Path(folder) / "capacity_input.json"
                seed_preset = str(
                    self.seed_cross_section_preset.currentData() or
                    "s8_r4x4")
                sst_lattice = _sst_lattice_for_symmetry(
                    cfg.lattice_symmetry)
                write_shifted_sst(
                    str(preview_sst), cfg.sst_growth_bp_z1,
                    cfg.spacer_bp_z2, cfg.sst_growth_bp_z3,
                    cfg.growth_bp_z1, cfg.growth_bp_z3,
                    32 if sst_lattice in ("kagome", "square_kagome")
                    else 16,
                    sst_lattice, seed_preset,
                    layers_design_sequence_identical=bool(
                        cfg.layers_design_sequence_identical),
                    mean_indel_per_helix=cfg.mean_indel_per_helix)
                capacity_report = estimate_scaffold_capacity(
                    str(preview_sst))
                requested_insertions = int(round(max(
                    0.0, float(cfg.mean_indel_per_helix)) * 48.0))
                insertion_headroom = int(capacity_report.get(
                    "seed_insertion_headroom_nt", 0))
                if requested_insertions > insertion_headroom:
                    raise RuntimeError(
                        "The selected Twist requires %d Seed insertions, "
                        "but the two fixed scaffold routes have room for "
                        "only %d nt while retaining the 7557-nt limit. "
                        "Reduce the positive Twist slightly or increase "
                        "the spacing." % (
                            requested_insertions, insertion_headroom))
        except Exception as error:
            self.seed_scaffold_capacity_status.setText(
                "容量不足或无法合法分配：%s" % error)
            match = re.search(r"总长度\D*(\d+)\s*nt", str(error))
            if match:
                QMessageBox.warning(
                    self, "Seed scaffold容量不足",
                    "骨架链实际长度：%d nt\n"
                    "3条骨架链总长度限制：22671 nt\n"
                    "请减小 Seed 长度。" % int(match.group(1)))
            else:
                QMessageBox.critical(
                    self, "设计参数验证失败", str(error))
            return
        finally:
            self._busy(False)
        signature = self._structure_signature(self.project.settings)
        previous = workflow.get("settings_signature")
        has_downstream = any(workflow.get(key) for key in (
            "sst_accepted", "scaffold_accepted", "structure_complete"))
        if previous and previous != signature and has_downstream:
            answer = QMessageBox.question(
                self, "参数修改会影响后续设计",
                "检测到影响SST/Seed结构的参数发生变化。继续后必须从SST设计"
                "重新生成；以前的JSON文件会保留，但不再作为当前接受版本。\n\n"
                "是否接受新参数并重新开始结构设计？",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes:
                return
            self._drop_keys(workflow, (
                "sst_review", "sst_two_layer", "sst_accepted",
                "sst_dependency_fingerprint", "scaffold_review",
                "scaffold_accepted", "scaffold_dependency_fingerprint",
                "structure_complete", "structure_accepted",
                "sequence_analysis", "sequence_assignments",
                "sequence_scaffold_accepted", "sequence_sst_accepted",
                "sequence_sst_accepted_at", "sequence_sst_detected",
                "sequence_sst_detection_status",
                "sequence_sst_import_method", "sequence_sst_import_status",
                "sequence_sst_imported_at", "sequence_sst_import_source",
                "sequence_sst_acceptance_status",
                "sequence_source", "sequence_exports"))
        unchanged = previous == signature and has_downstream
        workflow["settings_signature"] = signature
        workflow["parameters_accepted"] = True
        workflow["parameters_accepted_at"] = datetime.now().isoformat()
        workflow["scaffold_capacity_precheck"] = capacity_report
        self.seed_scaffold_capacity_status.setText("")
        workflow.pop("stale_notice", None)
        workflow.pop("parameters_editing", None)
        if unchanged and workflow.get("structure_accepted"):
            self.statusBar().showMessage(
                "参数没有影响后续设计的变化，已继续使用原接受结构。", 8000)
        else:
            self.statusBar().showMessage(
                "设计参数已接受，可直接生成Scaffold routing。", 8000)
            self._restore_structure_workflow()
        self.parameters_next_button.setEnabled(True)
        self._set_acceptance_button(
            self.accept_parameters_button, True,
            "接受当前 Moiré 参数",
            "✓ 当前 Moiré 参数已接受")
        self._show_action_feedback(
            self.parameters_action_status,
            "Moiré 参数已接受。下一步：打开 Automated DNA Design "
            "并生成三个设计文件。")
        self._record_history("接受当前 Moiré 参数")

    @staticmethod
    def _drop_keys(workflow, keys):
        for key in keys:
            workflow.pop(key, None)

    def _invalidate_downstream_design_state(self, workflow=None):
        """Reset Steps 2 and 3 after any accepted parameter changes.

        Generated files remain on disk, but no previously generated design,
        sequence assignment, or export is allowed to remain active in the
        current project workflow.  This is deliberately the same downstream
        state as a newly created project.
        """
        if workflow is None:
            workflow = self._workflow()
        self._drop_keys(workflow, (
            "scaffold_capacity_precheck",
            "sst_review", "sst_two_layer", "sst_accepted",
            "sst_dependency_fingerprint",
            "scaffold_review", "scaffold_accepted",
            "scaffold_dependency_fingerprint", "scaffold_editing",
            "structure_complete", "structure_accepted",
            "structure_accepted_at", "automatic_design_exports",
            "cadnano_inspection_files",
            "sequence_analysis", "sequence_assignments",
            "sequence_design_json", "sequence_scaffold_accepted",
            "sequence_scaffold_accepted_at", "sequence_sst_accepted",
            "sequence_sst_accepted_at", "sequence_sst_detected",
            "sequence_sst_detection_status",
            "sequence_sst_import_method", "sequence_sst_import_status",
            "sequence_sst_imported_at", "sequence_sst_import_source",
            "sequence_sst_acceptance_status", "sequence_source",
            "sequence_exports", "stale_notice"))
        self._sequence_analysis = None
        self._sequence_assignments = {}
        self._reset_downstream_design_ui()

    def _invalidate_after_sst(self):
        self._drop_keys(self._workflow(), (
            "sst_accepted", "scaffold_review", "scaffold_accepted",
            "structure_complete", "structure_accepted",
            "structure_accepted_at",
            "sequence_analysis", "sequence_assignments",
            "sequence_scaffold_accepted",
            "sequence_scaffold_accepted_at", "sequence_sst_accepted",
            "sequence_sst_accepted_at", "sequence_sst_detected",
            "sequence_sst_detection_status",
            "sequence_sst_import_method", "sequence_sst_import_status",
            "sequence_sst_imported_at", "sequence_sst_import_source",
            "sequence_sst_acceptance_status",
            "sequence_source", "sequence_exports"))
        self._sequence_analysis = None
        self._sequence_assignments = {}
        self._reset_sequence_assignment_ui()

    def _invalidate_after_scaffold(self):
        self._drop_keys(self._workflow(), (
            "scaffold_accepted", "structure_complete",
            "structure_accepted", "structure_accepted_at",
            "sequence_analysis",
            "sequence_assignments", "sequence_scaffold_accepted",
            "sequence_scaffold_accepted_at",
            "sequence_sst_accepted", "sequence_sst_accepted_at",
            "sequence_sst_detected", "sequence_sst_detection_status",
            "sequence_sst_import_method", "sequence_sst_import_status",
            "sequence_sst_imported_at", "sequence_sst_import_source",
            "sequence_sst_acceptance_status", "sequence_source",
            "sequence_exports"))
        self._sequence_analysis = None
        self._sequence_assignments = {}
        self._reset_sequence_assignment_ui()

    def _invalidate_after_sequence_scaffold_change(self):
        """Revoke sequence-stage acceptance after a scaffold reassignment."""
        workflow = self._workflow()
        self._drop_keys(workflow, (
            "sequence_scaffold_accepted",
            "sequence_scaffold_accepted_at",
            "sequence_sst_accepted", "sequence_sst_accepted_at",
            "sequence_sst_detected", "sequence_sst_detection_status",
            "sequence_sst_import_method", "sequence_sst_import_status",
            "sequence_sst_imported_at", "sequence_sst_import_source",
            "sequence_sst_acceptance_status", "sequence_source",
            "sequence_exports"))
        self._sequence_assignments = {
            target_id: assignment
            for target_id, assignment in self._sequence_assignments.items()
            if assignment.get("category") == "seed_scaffold"
        }
        self._store_sequence_assignments()
        self._reset_sequence_assignment_ui()
        self.detect_scaffold_sequences_button.setEnabled(bool(
            workflow.get("structure_accepted")))
        if self._sequence_analysis:
            targets = self._sequence_analysis.get(
                "targets", {}).get("seed_scaffold", [])
            self._render_scaffold_cards()
            complete = bool(targets) and all(
                target["id"] in self._sequence_assignments
                for target in targets)
            self.accept_added_scaffold_button.setEnabled(complete)
            design = self._sequence_design_path()
            if design:
                self.sequence_preview.set_source(
                    design, "sequence_scaffold", "Scaffold routing only")
                self.sequence_preview.set_sequence_scaffold_targets(targets)
                self.sequence_preview.set_sequence_scaffold_assignments(
                    self._sequence_assignments.values())
                self.sequence_preview_status.setText(
                    "Scaffold-only routing with assigned sequences.")
        self.sequence_export_status.setText(
            "Scaffold assignments changed. Accept them again before "
            "continuing to SST sublattice input sequences.")
        self._update_sequence_expert_actions()
        if self.project_path:
            self._save_current_project(silent=True)

    def _accepted_version(self, source, label):
        root = Path(self._workflow()["root"])
        base = self.project.settings.project_name+"_"+label+"_accepted_v"
        version = 1
        while (root/(base+"%03d.json" % version)).exists():
            version += 1
        target = root/(base+"%03d.json" % version)
        shutil.copy2(source, target)
        return str(target.resolve())

    @staticmethod
    def _dependency_fingerprint(filename, stage):
        """Hash only fields that can alter a later generated design."""
        payload = json.loads(Path(filename).read_text(encoding="utf-8"))
        rows = []
        for row in sorted(payload.get("vstrands", []),
                          key=lambda item: int(item["num"])):
            number = int(row["num"])
            item = {
                "num": number,
                "row": int(row["row"]),
                "col": int(row["col"]),
                "scaf": row.get("scaf", []),
                "loop": row.get("loop", []),
                "skip": row.get("skip", []),
            }
            if stage == "sst" or number >= 48:
                item["stap"] = row.get("stap", [])
            rows.append(item)
        encoded = json.dumps(rows, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(encoded.encode("utf-8")).hexdigest()

    @staticmethod
    def _propagate_visual_fields(source, targets, helix_filter=None):
        """Carry color-only expert edits forward without rerouting."""
        source_payload = json.loads(Path(source).read_text(encoding="utf-8"))
        source_rows = {int(row["num"]): row
                       for row in source_payload.get("vstrands", [])}
        for target in targets:
            if not target or not Path(target).is_file():
                continue
            payload = json.loads(Path(target).read_text(encoding="utf-8"))
            changed = False
            for row in payload.get("vstrands", []):
                number = int(row["num"])
                if helix_filter is not None and number not in helix_filter:
                    continue
                source_row = source_rows.get(number)
                if source_row is None:
                    continue
                colors = source_row.get("stap_colors", [])
                if row.get("stap_colors", []) != colors:
                    row["stap_colors"] = colors
                    changed = True
            if changed:
                Path(target).write_text(
                    json.dumps(payload, ensure_ascii=False, indent=2),
                    encoding="utf-8")

    def return_to_parameters(self):
        workflow = self._workflow()
        workflow["parameters_editing"] = True
        workflow["parameters_accepted"] = False
        self._restore_structure_workflow()
        self.tabs.setCurrentIndex(0)
        self.design_stack.setCurrentIndex(1)
        self._set_current_step(0)

    def return_to_scaffold(self):
        workflow = self._workflow()
        accepted = workflow.get("scaffold_accepted")
        workflow["scaffold_editing"] = True
        if accepted:
            workflow["scaffold_review"] = accepted
        self._restore_structure_workflow()
        self._open_structure_file("scaffold_review")

    def return_to_structure(self):
        self._invalidate_final_design_acceptance()
        self.tabs.setCurrentIndex(1)
        self._set_current_step(1)
        self._restore_structure_workflow()

    @staticmethod
    def _structure_signature(settings):
        return {
            "project_name": settings.project_name,
            "sst_z1": int(settings.sst_growth_bp_z1),
            "spacing_seed_z2": int(settings.spacer_bp_z2),
            "sst_z3": int(settings.sst_growth_bp_z3),
            "seed_z1": int(settings.growth_bp_z1),
            "seed_z3": int(settings.growth_bp_z3),
            "identical": bool(settings.layers_design_sequence_identical),
            "angle": round(float(settings.target_angle_deg), 6),
            "mean_indel": round(float(settings.mean_indel_per_helix), 6),
            "lattice_symmetry": settings.lattice_symmetry,
            "layer_lattice_constants": [
                round(float(settings.layer1_lattice_constant_nm), 6),
                round(float(settings.layer2_lattice_constant_nm), 6)],
            "seed_cross_section_cells": sorted(
                [list(map(int, cell))
                 for cell in settings.seed_cross_section_cells]),
        }

    def _settings(self):
        symmetry = self._selected_symmetry()
        layer1_a, layer2_a = self._context_lattice_constants(symmetry)
        seed_cells = self.seed_cross_section_picker.cells()
        return SquareBilayerSettings(
            project_name=self.project_name.text().strip() or "square_moire_bilayer",
            interface_language=current_language(),
            target_mode=self._target_driver,
            target_definition="local_surface",
            target_angle_deg=self.angle.value(),
            target_period_nm=self.period.value(),
            lattice_constant_nm=layer1_a,
            lattice_context=self.lattice_context.currentData()[0],
            lattice_symmetry=symmetry,
            layer1_lattice_constant_nm=layer1_a,
            layer2_lattice_constant_nm=layer2_a,
            seed_cross_section_size=8,
            seed_cross_section_cells=seed_cells,
            seed_template=(self.seed_cross_section_preset.currentText() +
                           " (%d helices)" % len(seed_cells)),
            growth_bp_z1=128,
            spacer_bp_z2=int(self.sst_spacing.currentData()),
            growth_bp_z3=128,
            sst_growth_bp_z1=int(self.sst_z1.currentData()),
            sst_growth_bp_z3=int(self.sst_z3.currentData()),
            mean_indel_per_helix=self.mean_indel.value(),
            layers_design_sequence_identical=bool(
                self.layers_identical.currentData()),
            auto_solve_spacer=False,
        )

    def _context_changed(self):
        context, square_a, kagome_a = self.lattice_context.currentData()
        symmetry = self._selected_symmetry()
        if symmetry == "square_square_c4":
            values = square_a, square_a
        elif symmetry == "kagome_kagome":
            values = kagome_a, kagome_a
        else:
            values = square_a, kagome_a
        with QSignalBlocker(self.lattice_constant), \
                QSignalBlocker(self.lattice_constant_2):
            self.lattice_constant.setValue(values[0])
            self.lattice_constant_2.setValue(values[1])
        self.lattice_constant.setEnabled(False)
        self.lattice_constant_2.setEnabled(False)
        if hasattr(self, "lattice_constant_fixed"):
            self.lattice_constant_fixed.setText("%.1f nm" % values[0])
            self.lattice_constant_2_fixed.setText("%.1f nm" % values[1])
        self._basis_selection_changed()
        self._control_changed("angle", 0)

    def _control_changed(self, driver, step=None):
        if self._updating:
            return
        if step is not None:
            if int(step) == 0:
                self._invalidate_design_basis_acceptance()
            elif int(step) == 1:
                self._invalidate_parameter_acceptance()
            self._set_current_step(step)
        self._target_driver = driver
        self.recalculate()
        self._record_history("更改设计参数")

    def recalculate(self):
        previous_workflow = None
        previous_measurements = []
        if self.project is not None:
            previous_workflow = self.project.seed_plan.get(
                "structure_workflow")
            previous_measurements = list(self.project.measurements)
        try:
            project = solve_square_bilayer(self._settings())
        except Exception as error:
            QMessageBox.warning(self, "参数错误", str(error))
            return
        if previous_workflow:
            project.seed_plan["structure_workflow"] = previous_workflow
        project.measurements = previous_measurements
        self.project = project
        cfg, prediction = project.settings, project.prediction
        self._updating = True
        try:
            with QSignalBlocker(self.angle), QSignalBlocker(self.period), \
                    QSignalBlocker(self.mean_indel):
                self.angle.setValue(cfg.target_angle_deg)
                if prediction.get("period_available"):
                    self.period.setValue(
                        cfg.target_period_nm if
                        math.isfinite(cfg.target_period_nm) else 0.0)
                self.angle.setSuffix(
                    "° (%s)" % _twist_handedness(cfg.target_angle_deg))
                self.mean_indel.setValue(cfg.mean_indel_per_helix)
                minimum_deletion = prediction.get(
                    "minimum_seed_deletion_per_helix",
                    minimum_seed_deletion_per_helix(cfg.spacer_bp_z2))
                maximum_insertion = prediction.get(
                    "maximum_seed_insertion_per_helix",
                    maximum_seed_insertion_per_helix(cfg.spacer_bp_z2))
                self.mean_indel.setSuffix(
                    " / helix (minimum %+.0f, maximum +%.0f)" % (
                        minimum_deletion,
                        maximum_insertion))
                self.mean_indel.setStyleSheet(
                    "color:#c62828;font-weight:700;"
                    if prediction.get("seed_indel_limit_exceeded") else
                    "")
                self.actual_z2_spacing.setValue(
                    prediction["actual_z2_spacing_bp"])
                self.seed_z2_readout.setText(
                    "%d bp（由 Layer spacing 强制联动）" %
                    int(cfg.spacer_bp_z2))
                try:
                    # The support readout is the actual Seed/SST duplex
                    # intersection on each side of Z2, not the requested SST
                    # length.  The shared preview partition already places
                    # the fixed Seed boundaries around the current Z2.
                    partition = prediction["preview_seed_partition"]
                    seed_partitions = partition["seed_partition_ranges"]
                    seed_support_ranges = [
                        seed_partitions[0], seed_partitions[2]]
                    canvas_shift = int(partition.get(
                        "coordinate_shift_bp", 0))
                    fixed_capture_grid = [
                        list(range(56 + canvas_shift,
                                   329 + canvas_shift, 16)),
                        list(range(56 + canvas_shift,
                                   329 + canvas_shift, 16)),
                    ]
                    overlap = fixed_seed_overlap_layout(
                        partition["sst_layer_ranges"],
                        lattice_type=cfg.lattice_symmetry,
                        seed_layer_ranges=seed_support_ranges,
                        seed_capture_positions_by_layer=fixed_capture_grid)
                    overlap_bp = [
                        int(partition["sst_overlap_z1_bp"]),
                        int(partition["sst_overlap_z3_bp"]),
                    ]
                    columns = overlap.get(
                        "capture_columns_by_layer",
                        [len(values) for values in overlap.get(
                            "capture_positions_by_layer", [[], []])])
                    self.seed_z1_overlap_readout.setText(
                        "%d bp · %d capture columns (minimum 4)" %
                        (overlap_bp[0], columns[0]))
                    self.seed_z3_overlap_readout.setText(
                        "%d bp · %d capture columns (minimum 4)" %
                        (overlap_bp[1], columns[1]))
                except (ValueError, IndexError):
                    self.seed_z1_overlap_readout.setText("—")
                    self.seed_z3_overlap_readout.setText("—")
        finally:
            self._updating = False
        zero_spacing = int(cfg.spacer_bp_z2) == 0
        mixed_lattice = cfg.lattice_symmetry == "square_kagome"
        self.angle.setEnabled(not zero_spacing)
        self.period.setEnabled(not mixed_lattice and not zero_spacing)
        self.angle.setToolTip(
            "0-bp spacing contains no 8-bp domain, so Twist is fixed at 0°."
            if zero_spacing else "")
        self.period.setToolTip(
            "0-bp spacing fixes Twist at 0° and the Moiré period at infinity."
            if zero_spacing else
            "Square–Kagome consists of different lattices and has no single "
            "Moiré period." if mixed_lattice else
            "Calculated from Twist and the current lattice constant.")
        reported_angle = prediction["reported_angle_deg"]
        handedness = _twist_handedness(reported_angle)
        self.angle_card.value.setText(
            "%+.1f° (%s)" % (reported_angle, handedness))
        period = prediction["predicted_moire_period_nm"]
        self.period_card.value.setText(
            ("—（不同点阵）" if period is None else
             "∞" if not math.isfinite(period) else "%.1f nm" % period))
        preview_parameters = _preview_parameter_html(
            "%+.1f° (%s)" % (reported_angle, handedness),
            "Not applicable" if period is None else
            "∞" if not math.isfinite(period) else "%.1f nm" % period)
        self.setup_preview_parameters.setText(preview_parameters)
        preview_partition = prediction["preview_seed_partition"]
        actual_z2_text = "%s bp" % _format_bp(
            prediction["actual_z2_spacing_bp"])
        self.side_preview_parameters.setText(_side_preview_parameter_html(
            "%d bp" % cfg.sst_growth_bp_z1,
            actual_z2_text,
            "%d bp" % cfg.sst_growth_bp_z3,
            "%d bp" % int(preview_partition["sst_overlap_z1_bp"]),
            actual_z2_text,
            "%d bp" % int(preview_partition["sst_overlap_z3_bp"])))
        self.design_preview_geometry_summary.setText(
            "dsDNA helix Ø 2.0 nm · Seed spacing 2.8 nm · "
            "lattice a1/a2 %.1f/%.1f nm" % (
                cfg.layer1_lattice_constant_nm,
                cfg.layer2_lattice_constant_nm))
        self.design_preview_length_summary.setText(
            "SST sublattice %d/%d bp · Seed Z1/Z2/Z3 %d/%s/%d bp "
            "(physical total %s bp)" % (
                cfg.sst_growth_bp_z1, cfg.sst_growth_bp_z3,
                int(prediction["preview_seed_partition"][
                    "sst_overlap_z1_bp"]),
                _format_bp(prediction["actual_z2_spacing_bp"]),
                int(prediction["preview_seed_partition"][
                    "sst_overlap_z3_bp"]),
                _format_bp(
                    int(prediction["preview_seed_partition"][
                        "sst_overlap_z1_bp"]) +
                    prediction["actual_z2_spacing_bp"] +
                    int(prediction["preview_seed_partition"][
                        "sst_overlap_z3_bp"]))))
        self.preview.set_design(project)
        self.moire_preview.set_design(project)
        self._update_validation()
        self._update_measurements()
        self._restore_structure_workflow()
        self.statusBar().showMessage(
            "实时更新：名义Z2=%d bp，实际Z2/spacing=%.1f bp，Twist %+.1f° (%s)，period %s" %
            (cfg.spacer_bp_z2, prediction["actual_z2_spacing_bp"],
             reported_angle, handedness,
             ("不适用" if period is None else
              "∞" if not math.isfinite(period) else "%.1f nm" % period)),
            4000)

    def _update_validation(self):
        colors = {"pass": "#217a4b", "warning": "#a86500",
                  "error": "#b43a3a", "info": "#356fab"}
        icons = {"pass": "✓", "warning": "!", "error": "×", "info": "i"}
        rows = []
        for item in self.project.validation:
            level = item["level"]
            rows.append(
                "<p style='margin:5px 2px'><b style='color:%s'>%s %s</b>　%s</p>" %
                (colors[level], icons[level], item["title"], item["detail"]))
        self.validation.setHtml("".join(rows))

    def _workflow(self):
        if self.project is None:
            self.recalculate()
        return self.project.seed_plan.setdefault("structure_workflow", {})

    def _ensure_structure_root(self):
        workflow = self._workflow()
        existing = workflow.get("root")
        project_root = self._project_directory()
        if existing and Path(existing).is_dir() and project_root and \
                project_root in Path(existing).resolve().parents:
            self.structure_root = Path(existing)
            return self.structure_root
        root = self._project_output_dir("cadnano design")
        if root is None:
            QMessageBox.warning(
                self, "尚未创建项目", "请先新建或打开 .moire.json 项目。")
            return None
        workflow["root"] = str(root)
        workflow["settings_signature"] = self._structure_signature(
            self.project.settings)
        self.structure_root = root
        return root

    def _structure_template_supported(self):
        cfg = self.project.settings
        selected_cells = {tuple(map(int, cell))
                          for cell in cfg.seed_cross_section_cells}
        preset = self._preset_key_for_cells(selected_cells)
        supported = (
            cfg.lattice_symmetry in (
                "square_square_c4", "kagome_kagome", "square_kagome")
            and preset == "s8_r4x4")
        if not supported:
            QMessageBox.information(
                self, "当前组合仅支持设计预测",
                "当前已开放 Square–Square、Kagome–Kagome 和 "
                "Square–Kagome "
                "S8–R4×4C 的后续设计。")
            return False
        sst_values = tuple(map(int, (
            cfg.sst_growth_bp_z1, cfg.spacer_bp_z2,
            cfg.sst_growth_bp_z3)))
        if sst_values[0] < 64 or sst_values[2] < 64:
            QMessageBox.warning(
                self, "Z长度不足",
                "SST superlattice 1st layer和2nd layer至少需要64 bp。")
            return False
        if any(value % 8 for value in sst_values):
            QMessageBox.warning(
                self, "SST superlattice长度不合法",
                "SST superlattice 1st layer、spacing和2nd layer必须是"
                "8 bp整数倍。")
            return False
        return True

    @staticmethod
    def _validation_text(report):
        parts = []
        if report.get("seed_scaffold_lengths"):
            parts.append("Seed scaffold：%s nt" % ", ".join(
                str(value) for value in report["seed_scaffold_lengths"]))
        if report.get("sst_ranges"):
            parts.append("两层Capture SST superlattice：%s" % "、".join(
                "%d–%d" % tuple(item["range"] if isinstance(item, dict)
                                    else item)
                for item in report["sst_ranges"]))
        if report.get("seed_requested_lengths"):
            parts.append("Seed设置Z1/Z3：%s bp；实际routing：%s bp" % (
                "/".join(map(str, report["seed_requested_lengths"])),
                "/".join(map(str, report["seed_routing_lengths"]))))
        if report.get("capture_columns_by_layer"):
            parts.append("Capture列：1st %d、2nd %d（两列为一组，末列可单独保留）" %
                         tuple(report["capture_columns_by_layer"]))
        if report.get("maximum_edge_stagger_bp") is not None:
            parts.append("Seed helix边界最大错位：%d bp" %
                         report["maximum_edge_stagger_bp"])
        if report.get("protected_capture_gap_endpoints"):
            parts.append("已保护capture缺口端点：%d个" %
                         report["protected_capture_gap_endpoints"])
        if report.get("capture_bridge_component_count"):
            parts.append("Seed–SST superlattice capture桥：%d条；"
                         "颜色组：%d" % (
                report["capture_bridge_component_count"],
                report.get("capture_color_count", 0)))
        if report.get("minimum_staple_length"):
            parts.append(
                "Seed scaffold staple覆盖缺失：%d bp；最短staple：%d nt" % (
                    report.get("seed_staple_missing_base_count", 0),
                    report["minimum_staple_length"]))
        if report.get("warnings"):
            parts.extend("注意："+item for item in report["warnings"])
        if report.get("errors"):
            parts.extend("错误："+item for item in report["errors"])
        return "\n".join(parts)

    def _busy(self, active):
        if active:
            QApplication.setOverrideCursor(QCursor(Qt.CursorShape.WaitCursor))
            if not hasattr(self, "_task_progress_dialog") or \
                    self._task_progress_dialog is None:
                self._task_progress_dialog = QProgressDialog(
                    "Processing, please wait…", "", 0, 0, self)
                self._task_progress_dialog.setWindowTitle(
                    "Processing")
                self._task_progress_dialog.setWindowModality(
                    Qt.WindowModality.NonModal)
                self._task_progress_dialog.setCancelButton(None)
                self._task_progress_dialog.setMinimumDuration(0)
                self._task_progress_dialog.setAutoClose(False)
                self._task_progress_dialog.show()
        else:
            QApplication.restoreOverrideCursor()
            dialog = getattr(self, "_task_progress_dialog", None)
            if dialog is not None:
                dialog.close()
                dialog.deleteLater()
                self._task_progress_dialog = None
        QApplication.processEvents()

    @staticmethod
    def _staple_analysis_text(report):
        # Normal staples are summarized by the validator's hard extrema. The
        # detailed histogram is shown when available in newer payloads.
        histogram = report.get("normal_staple_length_histogram", {})
        histogram_text = "  ".join(
            "%s:%s%%" % (key, value) for key, value in histogram.items()
            if float(value) > 0.0)
        parts = [
            translate("Normal staple：%s–%s nt" % (
                report.get("minimum_normal_staple_length", "—"),
                report.get("maximum_normal_staple_length", "—")))]
        if histogram_text:
            parts.append(translate("长度分布：") + histogram_text)
        if report.get("continuous_16_base_percentage") is not None:
            parts.append(translate(
                "具有连续16-base区域：%.1f%%" %
                float(report["continuous_16_base_percentage"])))
        return ("；" if current_language() == "zh_CN" else "; ").join(parts)

    def _toggle_structure_expert(self, enabled):
        self.scaffold_expert_box.setVisible(bool(enabled))
        self.staple_expert_box.setVisible(bool(enabled))
        self.structure_expert_button.setText(
            "Close expert mode" if enabled else
            "Optional · expert mode")

    def generate_simple_structure_design(self):
        """Generate all three design files; only the final file is accepted."""
        workflow = self._workflow()
        if not workflow.get("parameters_accepted"):
            QMessageBox.warning(
                self, "设计参数缺失", "请先接受 1.1 与 1.2 参数。")
            return
        if not self._structure_template_supported():
            return
        self.generate_scaffold_design()
        if not workflow.get("scaffold_review"):
            return
        self.generate_complete_structure(
            scaffold_source=workflow.get("scaffold_review"))
        if workflow.get("structure_complete"):
            workflow["automatic_design_exports"] = {
                "sst": workflow.get("sst_two_layer"),
                "sst_scaffold_routing": workflow.get("scaffold_review"),
                "sst_scaffold_routing_staple_capture": workflow.get(
                    "structure_complete"),
            }
            self._refresh_structure_preview("complete")
            self.inspect_final_design_button.setEnabled(True)
            self.accept_structure_button.setEnabled(True)
            self.statusBar().showMessage(
                "All three design files were exported. Only the final "
                "Staple/Capture design can be accepted.",
                8000)
            self._show_action_feedback(
                self.design_generation_action_status,
                "Generated the SST sublattice-only, SST sublattice + "
                "scaffold, and final staple/capture design files. Next: "
                "optionally inspect the final JSON in caDNAno, then accept "
                "it.")
            self._record_history("Generate and export all three design files")

    def generate_scaffold_design(self):
        workflow = self._workflow()
        if not workflow.get("parameters_accepted"):
            QMessageBox.warning(
                self, "设计参数缺失",
                "无法设计 Scaffold routing：尚未导入或接受设计参数。")
            return
        if not self._structure_template_supported():
            return
        if self.project.prediction.get("seed_indel_limit_exceeded"):
            minimum = self.project.prediction[
                "minimum_seed_deletion_per_helix"]
            maximum = self.project.prediction[
                "maximum_seed_insertion_per_helix"]
            required = self.project.settings.mean_indel_per_helix
            QMessageBox.warning(
                self, "Twist and spacing are incompatible",
                "The current Twist requires %+.1f bases/helix, outside the "
                "allowed %+.1f to %+.1f interval for the selected spacing. "
                "Return to the Moiré parameters and choose a compatible "
                "Twist or spacing before generating the design." % (
                    required, minimum, maximum))
            return
        if any(workflow.get(key) for key in (
                "scaffold_accepted", "structure_complete",
                "structure_accepted", "sequence_source")):
            answer = QMessageBox.question(
                self, "重新生成会使下游设计失效",
                "当前工程已有接受的Scaffold、Staple/Capture或序列结果。"
                "重新生成Scaffold后，这些下游结果必须重新生成；"
                "旧文件会保留。\n\n是否继续？",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes:
                return
        root = self._ensure_structure_root()
        if root is None:
            return
        workflow = self._workflow()
        sst_path = root/(self.project.settings.project_name+"_sst.json")
        path = root/(self.project.settings.project_name+"_sst_scaffold.json")
        self._busy(True)
        try:
            settings = self.project.settings
            seed_preset = str(
                self.seed_cross_section_preset.currentData() or "s8_r4x4")
            sst_lattice = _sst_lattice_for_symmetry(
                settings.lattice_symmetry)
            extension_nt = (32 if sst_lattice in
                            ("kagome", "square_kagome") else 16)
            write_shifted_sst(
                str(sst_path), settings.sst_growth_bp_z1,
                settings.spacer_bp_z2, settings.sst_growth_bp_z3,
                settings.growth_bp_z1, settings.growth_bp_z3,
                extension_nt, sst_lattice, seed_preset,
                layers_design_sequence_identical=bool(
                    settings.layers_design_sequence_identical),
                mean_indel_per_helix=settings.mean_indel_per_helix)
            sst_report = validate_sst(str(sst_path))
            if not sst_report["valid"]:
                raise RuntimeError("固定SST superlattice内部验证失败：%s" %
                                   "；".join(sst_report["errors"]))
            generate_scaffold_review(str(path), str(sst_path))
            report = validate_structure(str(path))
        except Exception as error:
            QMessageBox.critical(self, "Scaffold生成失败", str(error))
            return
        finally:
            self._busy(False)
        self._invalidate_after_sst()
        workflow["sst_two_layer"] = str(sst_path)
        workflow["sst_review"] = str(sst_path)
        workflow["sst_dependency_fingerprint"] = \
            self._dependency_fingerprint(str(sst_path), "sst")
        workflow["scaffold_review"] = str(path)
        workflow.pop("scaffold_editing", None)
        self.scaffold_status.setText(
            "固定两层SST superlattice已在后台生成并验证。\n"
            "已生成待审核Scaffold文件：%s\n"
            "在cadnano中编辑后请直接保存；接受时会自动读取该文件。\n%s" %
            (path.name, self._validation_text(report)))
        self.open_scaffold_button.setEnabled(True)
        self.accept_scaffold_button.setEnabled(report["valid"])
        self.generate_structure_button.setEnabled(True)
        self.open_structure_button.setEnabled(False)
        self.accept_structure_button.setEnabled(False)
        if hasattr(self, "sequence_export_status"):
            self.sequence_export_status.setText(
                "The structure changed. Assign the scaffold and SST "
                "sublattice input sequences again.")
            self.final_sequence_export_button.setEnabled(False)
        self._refresh_structure_preview("scaffold")
        self.statusBar().showMessage("Scaffold routing已生成，请先检查并接受。", 8000)

    def load_expert_scaffold(self):
        start = (str(self.structure_root) if self.structure_root else
                 str(Path.home()/"Desktop"))
        filename, unused = QFileDialog.getOpenFileName(
            self, "载入专家编辑后的 Scaffold JSON", start,
            "caDNAno design (*.json)")
        if not filename:
            return
        try:
            report = validate_structure(filename)
        except Exception as error:
            QMessageBox.critical(self, "JSON验证失败", str(error))
            return
        self.scaffold_status.setText(
            "已载入专家文件：%s\n%s" %
            (Path(filename).name, self._validation_text(report)))
        if not report["valid"]:
            QMessageBox.warning(
                self, "Scaffold不合格", "\n".join(report["errors"]))
            return
        self._workflow()["scaffold_review"] = str(Path(filename).resolve())
        self._workflow()["scaffold_editing"] = True
        self.open_scaffold_button.setEnabled(True)
        self.accept_scaffold_button.setEnabled(True)
        self._refresh_structure_preview("scaffold")

    def accept_scaffold(self):
        workflow = self._workflow()
        filename = workflow.get("scaffold_review")
        if not filename:
            return
        self._busy(True)
        try:
            report = validate_structure(filename)
        except Exception as error:
            QMessageBox.critical(self, "Scaffold 验证失败", str(error))
            return
        finally:
            self._busy(False)
        if not report["valid"]:
            QMessageBox.warning(
                self, "不能接受", "\n".join(report["errors"]))
            return
        fingerprint = self._dependency_fingerprint(filename, "scaffold")
        previous_fingerprint = workflow.get("scaffold_dependency_fingerprint")
        has_downstream = bool(workflow.get("structure_complete"))
        changed = bool(previous_fingerprint and
                       previous_fingerprint != fingerprint)
        if changed and has_downstream:
            answer = QMessageBox.question(
                self, "Scaffold修改需要重新生成Staple",
                "检测到Scaffold routing、crossover、indel或有效区间发生变化。"
                "现有Staple/Capture结构必须重新生成。\n\n"
                "是否接受修改并重新生成？旧JSON文件仍会保留。",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes:
                return
            self._invalidate_after_scaffold()
        elif previous_fingerprint == fingerprint and has_downstream:
            self._propagate_visual_fields(
                filename, (workflow.get("structure_complete"),))
        accepted = self._accepted_version(filename, "scaffold")
        workflow["scaffold_accepted"] = accepted
        workflow["scaffold_review"] = accepted
        workflow["scaffold_dependency_fingerprint"] = fingerprint
        workflow.pop("scaffold_editing", None)
        self._restore_structure_workflow()
        self.generate_structure_button.setEnabled(True)
        detail = self._validation_text(report)
        self.scaffold_status.setText(
            "已接受 Scaffold routing：%s\n%s" %
            (Path(accepted).name, detail))
        if report["warnings"]:
            QMessageBox.information(
                self, "Scaffold已接受",
                "结构硬检查已通过。\n"+"\n".join(report["warnings"]))
        if previous_fingerprint == fingerprint and has_downstream:
            self.statusBar().showMessage(
                "Scaffold没有影响后续设计的变化，已继续沿用现有结构。", 8000)
            if workflow.get("structure_accepted"):
                self._go_to_step(3)
        else:
            self.statusBar().showMessage(
                "Scaffold routing已锁定，可生成staple/capture。", 8000)

    def generate_complete_structure(self, scaffold_source=None):
        workflow = self._workflow()
        scaffold = scaffold_source or workflow.get("scaffold_accepted")
        if not scaffold or not Path(scaffold).is_file():
            QMessageBox.warning(
                self, "Scaffold routing 缺失",
                "无法设计 Staple / Capture：请先导入或接受 Scaffold routing JSON。")
            return
        if scaffold_source is None and workflow.get("structure_complete"):
            answer = QMessageBox.question(
                self, "重新生成 Staple / Capture",
                "已存在Staple / Capture结构。继续将生成新版本，"
                "并使已接受结构和序列导出状态失效；旧JSON会保留。\n\n"
                "是否继续？",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if answer != QMessageBox.StandardButton.Yes:
                return
            self._drop_keys(workflow, (
                "structure_accepted", "sequence_analysis",
                "sequence_assignments", "sequence_scaffold_accepted",
                "sequence_sst_accepted", "sequence_sst_accepted_at",
                "sequence_sst_detected", "sequence_sst_detection_status",
                "sequence_sst_import_method", "sequence_sst_import_status",
                "sequence_sst_imported_at", "sequence_sst_import_source",
                "sequence_sst_acceptance_status", "sequence_source",
                "sequence_exports"))
        root = Path(workflow["root"])
        path = root/(
            self.project.settings.project_name+
            "_sst_scaffold_staple_capture.json")
        self._busy(True)
        try:
            finalize_structure(scaffold, str(path))
            report = validate_structure(str(path), require_staples=True)
        except Exception as error:
            QMessageBox.critical(self, "结构生成失败", str(error))
            return
        finally:
            self._busy(False)
        workflow["structure_complete"] = str(path)
        workflow["structure_accepted"] = False
        workflow["automatic_design_exports"] = {
            "sst": workflow.get("sst_two_layer"),
            "sst_scaffold_routing": workflow.get("scaffold_review"),
            "sst_scaffold_routing_staple_capture": str(path),
        }
        self._invalidate_final_design_acceptance()
        self.structure_status.setText(
            "已导出全部设计文件：\n"
            "1. %s\n2. %s\n3. %s（最终可接受文件）\n%s" %
            (Path(workflow.get("sst_two_layer", "—")).name,
             Path(workflow.get("scaffold_review", "—")).name,
             path.name, self._validation_text(report)))
        self.open_structure_button.setEnabled(True)
        self.inspect_final_design_button.setEnabled(True)
        self.accept_structure_button.setEnabled(report["valid"])
        self._refresh_structure_preview("complete")
        self.statusBar().showMessage("Staple / Capture结构已生成，请检查。", 8000)

    def accept_complete_structure(self):
        workflow = self._workflow()
        filename = self._latest_complete_structure_file()
        if not filename:
            QMessageBox.warning(
                self, "没有可接受的最终设计",
                "当前最终文件和选择过的文件中，没有通过验证的 "
                "SST + Scaffold + Staple + Capture JSON。")
            return
        self._busy(True)
        try:
            report = validate_structure(filename, require_staples=True)
        except Exception as error:
            QMessageBox.critical(self, "结构验证失败", str(error))
            return
        finally:
            self._busy(False)
        if not report["valid"]:
            QMessageBox.warning(self, "不能接受", "\n".join(report["errors"]))
            return
        accepted = self._accepted_version(filename, "staple_capture")
        self._drop_keys(workflow, (
            "sequence_analysis", "sequence_assignments",
            "sequence_scaffold_accepted", "sequence_sst_accepted",
            "sequence_sst_accepted_at", "sequence_sst_detected",
            "sequence_sst_detection_status",
            "sequence_sst_import_method", "sequence_sst_import_status",
            "sequence_sst_imported_at", "sequence_sst_import_source",
            "sequence_sst_acceptance_status",
            "sequence_source", "sequence_exports"))
        self._sequence_analysis = None
        self._sequence_assignments = {}
        workflow["structure_complete"] = accepted
        workflow["structure_accepted"] = True
        workflow["structure_accepted_at"] = datetime.now().isoformat(
            timespec="seconds")
        self.final_sequence_export_button.setEnabled(False)
        self.scaffold_sequence_status.clear()
        self.scaffold_sequence_status.hide()
        self.scaffold_detection_action_status.clear()
        self.scaffold_detection_action_status.hide()
        self.sst_detection_status.setText(
            "Please accept the Scaffold sequences first.")
        for label in (self.sst_auto_import_status,
                      self.sst_expert_import_status,
                      self.sst_acceptance_status):
            label.clear()
            label.hide()
        self.sequence_export_status.setText("尚未完成序列导入。")
        self.structure_status.setText(
            "最终结构设计已接受：%s。\nSST 与 SST + Scaffold 是过程导出文件，"
            "不作为接受版本。第3步将只使用此最终JSON进行"
            "SST superlattice/capture序列设计。" %
            Path(accepted).name)
        # Make 3.1 immediately actionable. Previously this state changed only
        # after reopening/restoring the project, which made the button appear
        # to do nothing in the current session.
        self.detect_scaffold_sequences_button.setEnabled(True)
        self.accepted_design_summary.setText(
            "已接受设计图：%s\n接受时间：%s" %
            (Path(accepted).name, workflow["structure_accepted_at"]))
        self._refresh_structure_preview("complete")
        self._set_acceptance_button(
            self.accept_structure_button, True, "接受当前设计图",
            "✓ 当前设计图已接受")
        self.structure_next_button.setEnabled(True)
        self._show_action_feedback(
            self.structure_accept_action_status,
            "Final Staple/Capture design accepted. Next: open Sequence "
            "export and assign the scaffold and SST sublattice input "
            "sequences.")
        self._record_history("接受结构设计")

    def _activate_imported_project(self, filename):
        """Load a DNA Moiré Designer project and restore its workflow."""
        project = load_project(filename)
        self._sequence_analysis = None
        self._sequence_assignments = {}
        self.project = project
        self.project_path = str(Path(filename).resolve())
        self._load_settings(project.settings)
        self.recalculate()
        self.project.measurements = project.measurements
        self._update_measurements()
        self._restore_structure_workflow()

    def _import_moire_project(self, target):
        title = ("导入 Moiré 工程到 Scaffold / Capture 模块"
                 if target == "capture" else
                 "导入 Moiré 工程到序列模块")
        filename, unused = QFileDialog.getOpenFileName(
            self, title, str(Path.home()/"Desktop"),
            "Moiré project (*.moire.json)")
        if not filename:
            return
        try:
            payload = json.loads(Path(filename).read_text(encoding="utf-8"))
        except Exception as error:
            QMessageBox.critical(self, "JSON 无法读取", str(error))
            return
        if "seed_plan" not in payload or "settings" not in payload:
            QMessageBox.warning(
                self, "不是 Moiré 工程",
                "所选文件是 cadnano/普通 JSON，不含 DNA Moiré Designer 的参数与工作流。\n"
                "此入口只接受本软件保存的 .moire.json 工程。")
            return
        try:
            self._activate_imported_project(filename)
        except Exception as error:
            QMessageBox.critical(self, "Moiré 工程无法读取", str(error))
            return
        message = "已导入 Moiré 工程并恢复参数/流程：%s" % Path(filename).name
        if target == "capture":
            self.capture_import_status.setText(message)
        else:
            self.sequence_design_status.setText(message)
        self.statusBar().showMessage(message, 8000)

    def import_capture_moire_project(self):
        self._import_moire_project("capture")

    def import_sequence_moire_project(self):
        self._import_moire_project("sequence")

    def open_orthogonal_sequence_designer(self):
        """Open the designer inside this window, without cadnano startup."""
        sequence_root = self._project_output_dir(
            Path("SST sublattice input"))
        try:
            result = run_orthogonal_sequence_designer(
                self,
                project_filename=self.project_path,
                primer3_entries=self._orthogonal_primer3_entries,
                suggested_directory=sequence_root,
            )
        except Exception as error:
            QMessageBox.critical(self, "正交序列设计无法打开", str(error))
            return
        if result:
            self._orthogonal_primer3_entries = result["primer3_entries"]
            self._show_action_feedback(
                self.orthogonal_action_status,
                "Orthogonal input sequences exported. Next: copy them "
                "into the input template and import the completed workbook.")
            self.statusBar().showMessage(
                "正交序列已导出：%s" % result["filename"], 8000)

    @staticmethod
    def _clear_sequence_cards(layout):
        while layout.count():
            item = layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
            elif item.layout() is not None:
                MoireDesignerWindow._clear_sequence_cards(item.layout())

    def _sequence_design_path(self):
        workflow = self._workflow()
        design = workflow.get("structure_complete")
        if not design or not Path(design).is_file():
            QMessageBox.warning(
                self, "设计图缺失",
                "请先完成或导入 Scaffold + Staple + Capture + SST 设计。")
            return None
        return str(Path(design).resolve())

    def _store_sequence_assignments(self):
        self._workflow()["sequence_assignments"] = list(
            self._sequence_assignments.values())

    def _render_scaffold_cards(self):
        self._clear_sequence_cards(self.scaffold_cards_layout)
        self.sst_template_actions.hide()
        targets = (self._sequence_analysis or {}).get(
            "targets", {}).get("seed_scaffold", [])
        if not targets:
            self.scaffold_cards_layout.addWidget(QLabel("没有识别到 scaffold。"))
            return
        for index, target in enumerate(targets, 1):
            card = QGroupBox("Scaffold %d" % index)
            body = QVBoxLayout(card)
            color = target.get("color", "#1769aa")
            label = QLabel(
                "颜色：<span style='color:%s'>■</span> %s<br>"
                "位置：%s → %s<br>长度：<b>%d nt</b>" %
                (color, color, target["start"], target["end"],
                 target["length"]))
            label.setWordWrap(True)
            assignment = self._sequence_assignments.get(target["id"])
            assignment_name = ((assignment or {}).get("scaffold_name") or
                               (assignment or {}).get("source"))
            status = QLabel(
                ("Scaffold sequence assigned: %s. Next: assign sequences "
                 "to any remaining routes, or accept the assigned scaffold "
                 "sequences." % assignment_name
                 if assignment else
                 "No scaffold sequence assigned. Next: click Assign "
                 "scaffold sequence."))
            status.setWordWrap(True)
            button = QPushButton(
                "Change scaffold sequence" if assignment else
                "Assign scaffold sequence")
            button.setObjectName("primaryButton")
            button.clicked.connect(
                lambda checked=False, item=target:
                self._add_standard_sequence_scaffold(item))
            body.addWidget(label)
            body.addWidget(button)
            body.addWidget(status)
            self.scaffold_cards_layout.addWidget(card)
        complete = all(target["id"] in self._sequence_assignments
                       for target in targets)
        self.accept_added_scaffold_button.setEnabled(complete)

    def detect_sequence_scaffolds(self):
        design = self._sequence_design_path()
        if not design:
            return
        QApplication.processEvents()
        self._busy(True)
        try:
            self._sequence_analysis = analyze_sequence_design(design)
        except Exception as error:
            QMessageBox.critical(self, "Scaffold 读取失败", str(error))
            return
        finally:
            self._busy(False)
        targets = self._sequence_analysis["targets"]["seed_scaffold"]
        self.sequence_preview.set_source(
            design, "sequence_scaffold", "Scaffold routing only")
        self.sequence_preview.set_sequence_scaffold_targets(targets)
        self.sequence_preview.set_sequence_scaffold_assignments(
            self._sequence_assignments.values())
        self.sequence_preview_status.setText(
            "Scaffold-only routing；不同 scaffold 使用不同颜色。")
        self._render_scaffold_cards()
        self._show_action_feedback(
            self.scaffold_detection_action_status,
            "Detected %d scaffold routes. Next: assign a scaffold sequence "
            "to every route below." % len(targets))
        self._workflow()["sequence_analysis"] = {
            "path": design, "summary": self._sequence_analysis["summary"]}

    def _add_standard_sequence_scaffold(self, target):
        targets = (self._sequence_analysis or {}).get(
            "targets", {}).get("seed_scaffold", [])
        multiple = len(targets) > 1
        used_names = [
            assignment.get("scaffold_name")
            for target_id, assignment in self._sequence_assignments.items()
            if target_id != target["id"] and assignment.get("scaffold_name")]
        try:
            report = list_standard_scaffolds(
                target["length"], multiple, used_names)
        except Exception as error:
            QMessageBox.critical(self, "Scaffold 列表读取失败", str(error))
            return
        options = report.get("scaffolds", [])
        if not options:
            QMessageBox.warning(
                self, "没有可用 Scaffold",
                ("当前routing需要%d nt。多scaffold模式只允许 CS3L、"
                 "CS4、P7560，且同一名称不能重复使用。" %
                 target["length"]))
            return
        labels = ["%s · %s nt" %
                  (item["name"], format(item["length"], ","))
                  for item in options]
        current_name = (self._sequence_assignments.get(
            target["id"], {}).get("scaffold_name"))
        current_index = next((index for index, item in enumerate(options)
                              if item["name"] == current_name), 0)
        selected_label, accepted = QInputDialog.getItem(
            self, "Assign scaffold sequence",
            ("Select a built-in caDNAno scaffold for the %d-nt route:" %
             target["length"]), labels, current_index, False)
        if not accepted:
            return
        selected = options[labels.index(selected_label)]
        try:
            assignment = assign_standard_scaffold_sequence(
                target, selected["name"], multiple, used_names)
        except Exception as error:
            QMessageBox.critical(
                self, "Unable to assign scaffold sequence", str(error))
            return
        self._sequence_assignments[target["id"]] = assignment
        self._invalidate_after_sequence_scaffold_change()
        self.sequence_preview.set_sequence_scaffold_assignments(
            self._sequence_assignments.values())
        self._render_scaffold_cards()
        self._record_history("Assign built-in scaffold sequence")

    def _add_sequence_scaffold(self, target):
        design = self._sequence_design_path()
        if not design:
            return
        start = str(Path(design).parent)
        filename, unused = QFileDialog.getOpenFileName(
            self, "选择 cadnano 保存的带 Scaffold 序列 JSON", start,
            "caDNAno design (*.json)")
        if not filename:
            return
        self._busy(True)
        try:
            assignment = extract_scaffold_sequence(filename, target)
        except Exception as error:
            QMessageBox.critical(
                self, "Unable to assign scaffold sequence", str(error))
            return
        finally:
            self._busy(False)
        assignment.update({
            "target_id": target["id"],
            "start_vh": target["start_vh"],
            "start_idx": target["start_idx"],
            "length": target["length"],
            "category": "seed_scaffold",
            "layer": None,
        })
        self._sequence_assignments[target["id"]] = assignment
        self._invalidate_after_sequence_scaffold_change()
        self._render_scaffold_cards()
        self._record_history("Assign scaffold sequence")

    def accept_added_scaffolds(self):
        if not self._sequence_analysis:
            self.detect_sequence_scaffolds()
            if not self._sequence_analysis:
                return
        targets = self._sequence_analysis["targets"]["seed_scaffold"]
        missing = [item for item in targets
                   if item["id"] not in self._sequence_assignments]
        if missing:
            QMessageBox.warning(
                self, "Scaffold 尚未完整",
                "仍有 %d 条 scaffold 未添加序列。" % len(missing))
            return
        workflow = self._workflow()
        workflow["sequence_scaffold_accepted"] = True
        workflow["sequence_scaffold_accepted_at"] = datetime.now().isoformat(
            timespec="seconds")
        self._store_sequence_assignments()
        self.detect_sst_inputs_button.setEnabled(True)
        # Template rows do not exist until SST input positions are detected.
        # Orthogonal Sequence Design remains independently available.
        self.export_sst_input_template_button.setEnabled(False)
        self.import_sst_input_template_button.setEnabled(False)
        self._show_action_feedback(
            self.scaffold_sequence_status,
            "Accepted %d scaffold sequences. Next: detect the SST "
            "sublattice input positions and lengths." % len(targets))
        self._set_acceptance_button(
            self.accept_added_scaffold_button, True,
            "Accept assigned scaffold sequences",
            "✓ Scaffold sequences accepted")
        self._record_history("Accept assigned scaffold sequences")

    def _sst_layers_identical(self):
        if not self.project:
            return False
        # Linked lengths permit one shared input set only when both layers
        # also use the same SST topology. Square--Kagome has two different
        # strand graphs and therefore always keeps independent assignments.
        return bool(
            self.project.settings.layers_design_sequence_identical and
            self.project.settings.lattice_symmetry != "square_kagome")

    @staticmethod
    def _show_sequence_status(label, text):
        label.setText(str(text or ""))
        label.setVisible(bool(text))

    def _record_sst_import_status(self, method, message, source=None):
        workflow = self._workflow()
        self._drop_keys(workflow, (
            "sequence_sst_accepted", "sequence_sst_accepted_at",
            "sequence_sst_acceptance_status", "sequence_source",
            "sequence_exports"))
        workflow["sequence_sst_import_method"] = str(method)
        workflow["sequence_sst_import_status"] = str(message)
        workflow["sequence_sst_imported_at"] = datetime.now().isoformat(
            timespec="seconds")
        if source:
            workflow["sequence_sst_import_source"] = str(source)
        if method == "expert":
            self._show_sequence_status(
                self.sst_expert_import_status, message)
            self.sst_auto_import_status.hide()
        else:
            self._show_sequence_status(
                self.sst_auto_import_status, message)
            self.sst_expert_import_status.hide()
        self.sst_acceptance_status.hide()
        self.final_sequence_export_button.setEnabled(False)
        self._set_acceptance_button(
            self.accept_added_sst_button, False,
            "Accept assigned SST sublattice input sequences")
        if self.project_path:
            self._save_current_project(silent=True)

    def _sst_target_map(self):
        if not self._sequence_analysis:
            return {}
        return {
            item["id"]: item
            for layer in (1, 2)
            for item in self._sequence_analysis["targets"].get(
                "sst_input_layer_%d" % layer, [])
        }

    def _enrich_sst_assignments(self, assignments):
        """Attach stable structure positions omitted by the import worker."""
        target_map = self._sst_target_map()
        enriched = []
        for assignment in assignments:
            item = dict(assignment)
            target = target_map.get(item.get("target_id"), {})
            for key in ("start", "end", "color"):
                if key in target:
                    item[key] = target[key]
            enriched.append(item)
        return enriched

    def _sst_assignments_complete(self):
        target_map = self._sst_target_map()
        return bool(target_map) and all(
            target_id in self._sequence_assignments
            for target_id in target_map)

    @staticmethod
    def _sst_pairwise_quality(assignments):
        """Return independent worst-case pair metrics for each unique input."""
        sequences = list(dict.fromkeys(
            str(item.get("sequence", "")).upper()
            for item in assignments if item.get("sequence")))
        metrics = {sequence: {"same": 0, "complement": 0}
                   for sequence in sequences}
        for left_index, left in enumerate(sequences):
            for right in sequences[left_index + 1:]:
                same = longest_common_substring(left, right)
                complement = longest_common_substring(
                    left, reverse_complement(right))
                metrics[left]["same"] = max(
                    metrics[left]["same"], same)
                metrics[right]["same"] = max(
                    metrics[right]["same"], same)
                metrics[left]["complement"] = max(
                    metrics[left]["complement"], complement)
                metrics[right]["complement"] = max(
                    metrics[right]["complement"], complement)
        return metrics

    def _populate_sst_sequence_table(self, assignments):
        rows = self._enrich_sst_assignments(assignments)
        visible_layers = (1,) if self._sst_layers_identical() else (1, 2)
        visible_categories = {
            "sst_input_layer_%d" % layer for layer in visible_layers}
        rows = [item for item in rows if str(item.get("category", ""))
                in visible_categories and item.get("sequence")]
        rows.sort(key=lambda item: (
            int(item.get("start_idx", -1)),
            int(item.get("start_vh", -1)),
            int(item.get("length", 0))))
        self.sst_sequence_table.setRowCount(len(rows))
        if not rows:
            self.sst_sequence_report.setText(
                "No SST sublattice input sequences have been imported.")
            self.sst_sequence_table_box.hide()
            return
        quality = self._sst_pairwise_quality(rows)
        gc_values = []
        for row_index, item in enumerate(rows):
            sequence = str(item["sequence"]).upper()
            start = item.get(
                "start",
                "helix %d, base %d" % (
                    int(item.get("start_vh", -1)),
                    int(item.get("start_idx", -1))))
            end = item.get("end", "")
            position = "%s → %s" % (start, end) if end else str(start)
            gc_content = (100.0 * sum(
                base in "GC" for base in sequence) / len(sequence))
            gc_values.append(gc_content)
            values = (
                position,
                sequence,
                str(len(sequence)),
                "%.1f" % gc_content,
                str(quality[sequence]["same"]),
                str(quality[sequence]["complement"]),
            )
            for column, value in enumerate(values):
                cell = QTableWidgetItem(value)
                if column in (2, 3, 4, 5):
                    cell.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                self.sst_sequence_table.setItem(row_index, column, cell)
        maximum_same = max(
            (item["same"] for item in quality.values()), default=0)
        maximum_complement = max(
            (item["complement"] for item in quality.values()), default=0)
        self.sst_sequence_report.setText(
            "SST sublattice input sequences imported successfully. "
            "Inputs displayed: %d; total length: %d nt; GC content: "
            "%.1f–%.1f%%; maximum same-orientation exact match: %d nt; "
            "maximum interstrand complementarity: %d nt." % (
                len(rows), sum(len(str(item["sequence"])) for item in rows),
                min(gc_values), max(gc_values), maximum_same,
                maximum_complement))
        self.sst_sequence_table_box.show()
        QTimer.singleShot(0, self._balance_sequence_result_panels)
        self.sequence_preview.set_sst_input_assignments(rows)
        self.sequence_preview_status.setText(
            "SST sublattice input-only view with the assigned nucleotide "
            "sequences. Zoom the path view to inspect individual bases.")

    def _render_sst_cards(self):
        self._clear_sequence_cards(self.sst_cards_layout)
        self._update_sequence_expert_actions()
        expert_visible = bool(self.sequence_expert_button.isChecked())
        self.sst_template_actions.setVisible(expert_visible)
        if expert_visible:
            QTimer.singleShot(0, self._balance_sequence_result_panels)
        if not self._sequence_analysis:
            self.sst_cards_layout.addWidget(
                QLabel("SST sublattice inputs have not been detected."))
            return
        layers = (1,) if self._sst_layers_identical() else (1, 2)
        for layer in layers:
            key = "sst_input_layer_%d" % layer
            summary = self._sequence_analysis["summary"][key]
            lengths = "，".join(
                "%s nt × %s" % (length, count)
                for length, count in summary["lengths"].items()) or "—"
            title = ("SST sublattice input" if len(layers) == 1 else
                     "SST sublattice %s layer input" %
                     ("1st" if layer == 1 else "2nd"))
            card = QGroupBox(title)
            card.setSizePolicy(
                QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed)
            body = QVBoxLayout(card)
            body.setSizeConstraint(QLayout.SizeConstraint.SetMinimumSize)
            if len(layers) == 1:
                identical_note = QLabel(
                    "两层相同：只需导入一次，另一层自动映射。")
                identical_note.setWordWrap(True)
                identical_note.setObjectName("structureNote")
                body.addWidget(identical_note)
            metrics = QLabel(
                "数量：<b>%d</b>　总长度：<b>%d nt</b><br>%s" %
                (summary["count"], summary["total_nt"], lengths))
            metrics.setWordWrap(True)
            metrics.setMinimumHeight(42)
            body.addWidget(metrics)
            target_ids = [item["id"] for item in
                          self._sequence_analysis["targets"][key]]
            added = sum(item in self._sequence_assignments
                        for item in target_ids)
            status = QLabel("Assigned inputs: %d / %d" %
                            (added, len(target_ids)))
            status.setWordWrap(True)
            status.setMinimumHeight(24)
            body.addWidget(status)
            self.sst_cards_layout.addWidget(card)
        expected = [item["id"] for layer in (1, 2) for item in
                    self._sequence_analysis["targets"][
                        "sst_input_layer_%d" % layer]]
        self.accept_added_sst_button.setEnabled(
            bool(expected) and all(item in self._sequence_assignments
                                   for item in expected))
        self._populate_sst_sequence_table(
            self._sequence_assignments.values())

    def detect_sequence_sst_inputs(self):
        if not self._workflow().get("sequence_scaffold_accepted"):
            QMessageBox.warning(
                self, "Scaffold sequences not accepted",
                "Assign and accept the scaffold sequences first.")
            return
        design = self._sequence_design_path()
        if not design:
            return
        if not self._sequence_analysis or \
                self._sequence_analysis.get("path") != design:
            self._busy(True)
            try:
                self._sequence_analysis = analyze_sequence_design(design)
            except Exception as error:
                QMessageBox.critical(
                    self, "Unable to detect SST sublattice inputs",
                    str(error))
                return
            finally:
                self._busy(False)
        self.sequence_preview.set_source(
            design, "sst_input", "SST sublattice input only")
        layers = (1,) if self._sst_layers_identical() else (1, 2)
        targets = [
            item for layer in layers
            for item in self._sequence_analysis["targets"][
                "sst_input_layer_%d" % layer]
        ]
        self.sequence_preview.set_sst_input_targets(targets)
        self.sequence_preview_status.setText(
            "SST sublattice input-only view. Identical bilayers use one "
            "assignment group, but the two layers are exported separately.")
        self._render_sst_cards()
        detection_message = (
            "SST sublattice input positions and lengths detected. "
            "Next: design and assign the inputs automatically, or open "
            "expert mode.")
        self._show_action_feedback(
            self.sst_detection_status, detection_message)
        workflow = self._workflow()
        workflow["sequence_sst_detected"] = True
        workflow["sequence_sst_detection_status"] = detection_message
        workflow["sequence_analysis"] = {
            "path": design, "summary": self._sequence_analysis["summary"]}
        self.auto_design_sst_inputs_button.setEnabled(bool(targets))
        self.sequence_expert_button.setEnabled(True)
        self._update_sequence_expert_actions()
        if self.project_path:
            self._save_current_project(silent=True)

    def _sequence_expert_template_ready(self):
        workflow = self._workflow()
        return bool(
            workflow.get("sequence_scaffold_accepted") and
            workflow.get("sequence_sst_detected") and
            self._sequence_analysis)

    def _update_sequence_expert_actions(self):
        """Enable templates only after detection; orthogonal stays free."""
        ready = self._sequence_expert_template_ready()
        self.export_sst_input_template_button.setEnabled(ready)
        self.import_sst_input_template_button.setEnabled(ready)
        self.open_orthogonal_sequences_button.setEnabled(True)
        return ready

    def _toggle_sequence_expert(self, enabled):
        visible = bool(enabled)
        template_ready = self._update_sequence_expert_actions()
        self.sst_template_actions.setVisible(visible)
        self.sequence_expert_button.setText(
            "Close expert mode" if enabled else
            "Optional: expert mode")
        if visible:
            if not self.sst_expert_import_status.text().strip():
                self._show_sequence_status(
                    self.sst_expert_import_status,
                    ("Expert workflow opened. Export or import an input "
                     "template, or design orthogonal sequences."
                     if template_ready else
                     "Expert mode opened. Orthogonal sequence design is "
                     "available now. Accept the scaffold sequences and "
                     "detect the SST sublattice input positions to enable "
                     "template "
                     "export and import."))
            else:
                self.sst_expert_import_status.show()
            QTimer.singleShot(0, self._balance_sequence_result_panels)
        else:
            # Preserve an expert-import result for reopening/project restore,
            # but never leave its status card visible below a closed panel.
            self.sst_expert_import_status.hide()

    def _balance_sequence_result_panels(self):
        """Initialize a newly shown pane without resetting user resizing."""
        splitter = getattr(self, "sequence_results_splitter", None)
        if splitter is None:
            return
        widgets = (
            self.sequence_preview_box,
            self.sst_sequence_table_box,
            self.sst_template_actions)
        visible = [index for index, widget in enumerate(widgets)
                   if widget.isVisible()]
        if not visible:
            return
        visibility_state = tuple(widget.isVisible() for widget in widgets)
        if visibility_state == getattr(
                self, "_sequence_result_visibility_state", None):
            return
        self._sequence_result_visibility_state = visibility_state
        sizes = splitter.sizes()
        # Preserve every user-selected ratio once all visible panes have a
        # non-zero extent. This method is called repeatedly as sequence data
        # and Expert Mode are refreshed, so unconditional setSizes() made the
        # handles appear draggable while immediately restoring old sizes.
        if len(sizes) == len(widgets) and all(sizes[index] > 0
                                              for index in visible):
            return
        available = max(splitter.height(), sum(sizes), 360)
        preferred_weights = (5, 3, 3)
        weight_total = sum(preferred_weights[index] for index in visible)
        initialized = [0, 0, 0]
        for index in visible:
            initialized[index] = max(
                1, int(available * preferred_weights[index] / weight_total))
        splitter.setSizes(initialized)

    def auto_design_and_add_sst_inputs(self):
        """Qt-safe entry point: a slot exception must never abort the app."""
        if self._auto_input_design_running:
            return
        self._auto_input_design_running = True
        self.auto_design_sst_inputs_button.setEnabled(False)
        try:
            self._auto_design_and_add_sst_inputs_impl()
        except Exception as error:
            # PyQt6 aborts the process when an exception escapes a button
            # callback. Keep the entire workflow inside one explicit boundary.
            if self._sst_assignments_complete():
                self._render_sst_cards()
                self._record_sst_import_status(
                    "automatic", "SST sublattice input sequences assigned "
                    "successfully. "
                    "Next: review the analysis on the right and accept the "
                    "assigned inputs.")
                QMessageBox.information(
                    self, "Done",
                    "SST sublattice input sequences were assigned "
                    "successfully. "
                    "The detailed report is displayed on the right.")
            else:
                self._show_sequence_status(
                    self.sst_auto_import_status,
                    "Automatic SST sublattice input-sequence design and "
                    "assignment did not "
                    "complete. See the error dialog for details.")
                dialog = QMessageBox(
                    QMessageBox.Icon.Critical,
                    "Automatic SST sublattice input-sequence design",
                    "The SST sublattice input sequences could not be "
                    "designed and assigned.",
                    QMessageBox.StandardButton.Ok, self)
                dialog.setDetailedText(str(error))
                dialog.exec()
        finally:
            self._auto_input_design_running = False
            has_targets = bool(
                self._sequence_analysis and any(
                    self._sequence_analysis.get("targets", {}).get(
                        "sst_input_layer_%d" % layer, [])
                    for layer in (1, 2)))
            self.auto_design_sst_inputs_button.setEnabled(has_targets)

    def _auto_design_and_add_sst_inputs_impl(self):
        if not self._sequence_analysis:
            self.detect_sequence_sst_inputs()
            if not self._sequence_analysis:
                return
        layers = (1,) if self._sst_layers_identical() else (1, 2)
        targets = [
            item for layer in layers
            for item in self._sequence_analysis["targets"][
                "sst_input_layer_%d" % layer]
        ]
        targets.sort(key=lambda item: (
            int(item["start_idx"]), int(item["start_vh"]),
            int(item["length"])))
        if not targets:
            QMessageBox.warning(
                self, "No SST sublattice inputs",
                "No SST sublattice input positions were detected.")
            return
        length_counts = {}
        for target in targets:
            length_counts[int(target["length"])] = \
                length_counts.get(int(target["length"]), 0) + 1
        output_root = self._project_output_dir(Path("SST sublattice input"))
        if output_root is None:
            return
        automatic_root = output_root / "automatic_orthogonal_input"
        self._busy(True)
        try:
            generated = generate_orthogonal_sequences_automatic(
                length_counts, automatic_root, self)
            pools = {int(length): list(values) for length, values in
                     generated["by_length"].items()}
            rows = []
            for target in targets:
                length = int(target["length"])
                if not pools.get(length):
                    raise RuntimeError("%d-nt 序列分配数量不足。" % length)
                sequence = pools[length].pop(0)
                rows.append((target["start"], target["end"], sequence,
                             length, target["color"]))
            template = automatic_root / (
                self.project.settings.project_name +
                "_sst_input_filled_template.xlsx")
            write_sequence_template(str(template), rows)
            localize_xlsx(template)
            report = import_sst_input_template(
                self._sequence_design_path(), str(template),
                self._sst_layers_identical())
        except Exception as error:
            QMessageBox.critical(self, "自动 Input 序列设计失败", str(error))
            return
        finally:
            self._busy(False)
        enriched_assignments = self._enrich_sst_assignments(
            report["assignments"])
        for assignment in enriched_assignments:
            self._sequence_assignments[assignment["target_id"]] = assignment
        self._store_sequence_assignments()
        self._render_sst_cards()
        success_message = (
            "SST sublattice input sequences were designed and assigned "
            "successfully. Next: review the analysis on the right and "
            "accept the assigned inputs.")
        self._record_sst_import_status(
            "automatic", success_message, str(template))
        QMessageBox.information(
            self, "Done",
            "SST sublattice input sequences were designed and assigned "
            "successfully. "
            "The detailed report is displayed on the right.")
        self._record_history(
            "Design and assign SST sublattice inputs automatically")

    def export_sequence_sst_template(self):
        design = self._sequence_design_path()
        if not design:
            return
        output_root = self._project_output_dir(Path("SST sublattice input"))
        if output_root is None:
            return
        filename = str(output_root /
                       (Path(design).stem + "_sst_input_template.xlsx"))
        try:
            report = export_sst_input_template(
                design, filename, self._sst_layers_identical())
        except Exception as error:
            QMessageBox.critical(self, "Template 导出失败", str(error))
            return
        self._show_sequence_status(
            self.sst_expert_import_status,
            "Exported a template for %d SST sublattice input sequences: %s. "
            "Enter the sequences in the Sequence column and then import "
            "the workbook; sequences can also be generated with the "
            "orthogonal-sequence designer." %
            (report["row_count"], Path(filename).name))
        self._show_action_feedback(
            self.export_template_action_status,
            "Input template exported. Next: fill its Sequence column, then "
            "use Import and assign input sequences.")

    def import_sequence_sst_template(self):
        design = self._sequence_design_path()
        if not design:
            return
        QMessageBox.information(
            self, "Import and assign SST sublattice inputs",
            "Enter the sequences in the Sequence column of the exported "
            "input template. Input sequences can also be generated with "
            "the independent orthogonal sequence design tool below.")
        filename, unused = QFileDialog.getOpenFileName(
            self, "Import and assign SST sublattice inputs",
            str(Path(design).parent),
            "Excel workbook (*.xlsx)")
        if not filename:
            return
        self._busy(True)
        try:
            report = import_sst_input_template(
                design, filename, self._sst_layers_identical())
        except Exception as error:
            QMessageBox.critical(
                self, "Unable to assign SST sublattice inputs", str(error))
            return
        finally:
            self._busy(False)
        enriched_assignments = self._enrich_sst_assignments(
            report["assignments"])
        for assignment in enriched_assignments:
            self._sequence_assignments[assignment["target_id"]] = assignment
        self._store_sequence_assignments()
        self._render_sst_cards()
        success_message = (
            "SST sublattice input sequences assigned successfully.")
        self._record_sst_import_status(
            "expert", success_message, str(Path(filename).resolve()))
        self._show_action_feedback(
            self.import_template_action_status,
            "SST sublattice input sequences assigned successfully. Next: "
            "review the sequence analysis, then accept the assigned inputs.")
        self._record_history("Import and assign SST sublattice inputs")

    def accept_added_sst_inputs(self):
        if not self._sequence_analysis:
            QMessageBox.warning(
                self, "SST sublattice inputs not detected",
                "Detect the SST sublattice input positions first.")
            return
        all_targets = [item for group in self._sequence_analysis["targets"].values()
                       for item in group]
        missing = [item for item in all_targets
                   if item["id"] not in self._sequence_assignments]
        if missing:
            QMessageBox.warning(
                self, "序列尚未完整",
                "%d scaffold or SST sublattice input sequences have not "
                "been assigned." % len(missing))
            return
        design = self._sequence_design_path()
        if not design:
            return
        root = Path(self._workflow().get("root") or Path(design).parent)
        output = root / (self.project.settings.project_name +
                         "_with_sequences.json")
        self._busy(True)
        try:
            report = build_sequenced_design(
                design, str(output), self._sequence_assignments.values())
        except Exception as error:
            QMessageBox.critical(self, "序列应用失败", str(error))
            return
        finally:
            self._busy(False)
        if report.get("unresolved_output_bases"):
            QMessageBox.warning(
                self, "仍有未定输出序列",
                "结构中仍有 %d 个无法由 input 互补得到的碱基。" %
                report["unresolved_output_bases"])
            return
        workflow = self._workflow()
        workflow["sequence_sst_accepted"] = True
        workflow["sequence_sst_accepted_at"] = datetime.now().isoformat(
            timespec="seconds")
        workflow["sequence_source"] = report["path"]
        self.final_sequence_export_button.setEnabled(True)
        acceptance_message = (
            "Assigned SST sublattice inputs accepted. The nucleotide "
            "sequences have been written to both layers. Next: export the "
            "final package.")
        workflow["sequence_sst_acceptance_status"] = acceptance_message
        self._show_sequence_status(
            self.sst_acceptance_status, acceptance_message)
        self.sequence_export_status.setText(
            "Sequence assignment is complete. The final package is ready "
            "to export.")
        self._set_acceptance_button(
            self.accept_added_sst_button, True,
            "Accept assigned SST sublattice input sequences",
            "✓ SST sublattice input sequences accepted")
        if self.project_path:
            self._save_current_project(silent=True)
        self._record_history(
            "Accept assigned SST sublattice input sequences")

    @staticmethod
    def _render_export_widget(widget, png_path, svg_path):
        widget.ensurePolished()
        size = widget.size()
        if size.width() < 2 or size.height() < 2:
            widget.resize(max(900, widget.minimumWidth()),
                          max(600, widget.minimumHeight()))
            size = widget.size()
        pixmap = QPixmap(size.width() * 2, size.height() * 2)
        pixmap.setDevicePixelRatio(2.0)
        pixmap.fill(Qt.GlobalColor.white)
        painter = QPainter(pixmap)
        widget.render(painter)
        painter.end()
        localize_svg(svg_path)
        pixmap.save(str(png_path), "PNG")
        generator = QSvgGenerator()
        generator.setFileName(str(svg_path))
        generator.setSize(size)
        generator.setViewBox(QRect(0, 0, size.width(), size.height()))
        generator.setTitle(svg_path.stem)
        painter = QPainter(generator)
        painter.fillRect(QRect(0, 0, size.width(), size.height()),
                         Qt.GlobalColor.white)
        widget.render(painter)
        painter.end()

    def _write_final_input_figures(self, input_root):
        input_root = Path(input_root)
        input_root.mkdir(parents=True, exist_ok=True)
        save_project(
            self.project,
            str(input_root /
                (self.project.settings.project_name + ".moire.json")))
        summary = QWidget()
        summary.resize(1200, 900)
        layout = QVBoxLayout(summary)
        heading = QLabel("DNA Moiré Designer — Design Parameters")
        heading.setObjectName("sectionTitle")
        layout.addWidget(heading)
        text = QTextBrowser()
        text.setPlainText(json.dumps(
            self.project.to_dict(), ensure_ascii=False, indent=2))
        text.document().setTextWidth(1140)
        text.setFixedHeight(max(
            700, int(text.document().size().height()) + 28))
        layout.addWidget(text)
        summary.resize(1200, text.height() + 90)
        self._render_export_widget(
            summary, input_root / "design_parameters.png",
            input_root / "design_parameters.svg")
        for widget, stem in (
                (self.preview, "design_3d_side"),
                (self.moire_preview, "design_2d_top"),
                (self.seed_cross_section_picker, "seed_cross_section")):
            self._render_export_widget(
                widget, input_root / (stem + ".png"),
                input_root / (stem + ".svg"))

    def export_sequence_final_package(self):
        workflow = self._workflow()
        source = workflow.get("sequence_source")
        if not workflow.get("sequence_sst_accepted") or not source or \
                not Path(source).is_file():
            QMessageBox.warning(
                self, "Sequence assignment incomplete",
                "Accept the assigned scaffold sequences and SST sublattice "
                "inputs first.")
            return
        parent = self._project_output_dir("final_export")
        if parent is None:
            return
        self._busy(True)
        try:
            report = export_final_package(
                self.project_path or "", source, str(parent), workflow)
            self._write_final_input_figures(
                Path(report["root"]) / "Input parameters")
        except Exception as error:
            QMessageBox.critical(self, "Final export failed", str(error))
            return
        finally:
            self._busy(False)
        workflow["sequence_exports"] = report
        export_message = (
            "Final export completed successfully.\n"
            "Output folder: %s\n"
            "Sequences: Oligonucleotide sequences\n"
            "Structures: PDB/oxView files\n"
            "Next: open the output folder and review the exported files." %
            report["root"])
        self._show_action_feedback(
            self.sequence_export_status, export_message)
        QMessageBox.information(
            self, "Final export complete", report["root"])

    def load_sequence_design(self):
        workflow = self._workflow()
        design = workflow.get("structure_complete")
        if (not workflow.get("structure_accepted") or not design or
                not Path(design).is_file()):
            QMessageBox.warning(
                self, "设计图缺失",
                "无法自动读取序列：请先导入或接受完整的 Staple / Capture 设计 JSON。")
            return
        start = (str(Path(workflow.get("structure_complete", "")).parent)
                 if workflow.get("structure_complete") else
                 str(Path.home()/"Desktop"))
        suggested = workflow.get("sequence_design_json")
        if suggested and Path(suggested).is_file():
            filename = suggested
        else:
            filename, unused = QFileDialog.getOpenFileName(
                self, "载入使用 Save as with Sequences 保存的设计", start,
                "caDNAno design (*.json)")
        if not filename:
            return
        try:
            payload = json.loads(Path(filename).read_text(encoding="utf-8"))
            report = validate_structure(filename, require_staples=True)
        except Exception as error:
            QMessageBox.critical(self, "无法载入", str(error))
            return
        if not report["valid"]:
            QMessageBox.warning(
                self, "结构不符合当前设计", "\n".join(report["errors"]))
            return
        sequences = payload.get("scaffold_sequences", [])
        if not sequences:
            QMessageBox.warning(
                self, "没有保存序列",
                "该JSON不含scaffold序列。请在cadnano中导入序列后使用"
                " Save as with Sequences 保存。")
            return
        workflow["sequence_source"] = str(Path(filename).resolve())
        self.sequence_export_status.setText(
            "已载入：%s\n检测到 %d 条已保存scaffold序列；可后台导出两种状态。" %
            (Path(filename).name, len(sequences)))
        self.export_sequence_variants_button.setEnabled(True)

    def export_capture_output_sequences(self):
        workflow = self._workflow()
        design = workflow.get("structure_complete")
        if (not workflow.get("structure_accepted") or not design or
                not Path(design).is_file()):
            QMessageBox.warning(
                self, "设计图缺失",
                "无法导出序列：请先导入或接受完整的 Staple / Capture 设计 JSON。")
            return
        source = workflow.get("sequence_source")
        if not source or not Path(source).is_file():
            QMessageBox.warning(self, "尚未载入", "请先载入带序列的JSON。")
            return
        parent = self._project_output_dir(Path("SST sublattice input"))
        if parent is None:
            return
        base = self.project.settings.project_name if self.project else Path(source).stem
        self._busy(True)
        try:
            report = export_sequence_variants(source, str(parent), base)
        except Exception as error:
            QMessageBox.critical(self, "导出失败", str(error))
            return
        finally:
            self._busy(False)
        capture = report["variants"]["capture"]
        complete_sst = report["variants"]["complete_sst"]
        workflow["sequence_exports"] = report
        self.sequence_export_status.setText(
            "导出完成，当前设计未改变。\n"
            "Capture：%d条序列 + 对应SVG\n"
            "完整SST-only：%d条32-nt SST序列 + 对应JSON/SVG"
            "（未定单链碱基：%d）\n%s" %
            (capture["sequence_count"], complete_sst["sequence_count"],
             complete_sst.get("unresolved_base_count", 0), parent))
        QMessageBox.information(
            self, "导出完成",
            "已分别导出 capture/完整SST-only 的 JSON、XLSX 与同步设计SVG：\n%s" %
            parent)

    def _open_structure_file(self, key):
        filename = self._workflow().get(key)
        if not filename:
            return
        if not CADNANO_EXECUTABLE.is_file():
            QMessageBox.warning(self, "caDNAno", "未找到当前 cadnano 启动程序。")
            return
        subprocess.Popen(
            [str(CADNANO_EXECUTABLE), str(filename)],
            start_new_session=True, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL)

    def open_cadnano_for_optional_editing(self):
        """Open cadnano directly; users choose the design inside cadnano."""
        workflow = self._workflow()
        latest = workflow.get("structure_complete")
        if not latest or not Path(latest).is_file():
            QMessageBox.warning(
                self, "Design file not generated",
                "Generate all three design files before opening caDNAno "
                "expert mode.")
            self.inspect_final_design_button.setEnabled(False)
            return
        if not CADNANO_EXECUTABLE.is_file():
            QMessageBox.warning(self, "caDNAno", "未找到当前 cadnano 启动程序。")
            return
        self._invalidate_final_design_acceptance()
        start = (Path(latest).parent if latest else
                 Path(workflow.get("root")) if workflow.get("root") else
                 self.structure_root if self.structure_root else
                 Path.home() / "Desktop")
        start = Path(start).expanduser()
        kwargs = {
            "start_new_session": True,
            "stdout": subprocess.DEVNULL,
            "stderr": subprocess.DEVNULL,
        }
        if start.is_dir():
            kwargs["cwd"] = str(start)
        subprocess.Popen([str(CADNANO_EXECUTABLE)], **kwargs)
        self._show_action_feedback(
            self.cadnano_edit_action_status,
            "caDNAno opened. Modify and save the final JSON there, then "
            "return and click Accept Current Design.")
        self.statusBar().showMessage(
            "cadnano 已打开。请在其中选择并修改最终 JSON；保存后点击接受，"
            "软件会使用当前设计文件夹中修改时间最新的合法最终文件。", 8000)

    def _latest_complete_structure_file(self):
        """Return the newest valid final design among generated/edited files."""
        workflow = self._workflow()
        automatic = workflow.get("automatic_design_exports", {})
        candidates = [
            workflow.get("structure_complete"),
            automatic.get("sst_scaffold_routing_staple_capture"),
        ]
        candidates.extend(workflow.get("cadnano_inspection_files", []))
        # cadnano is opened without a filename so its own Open command remains
        # in control.  Include every JSON saved in the current design folder;
        # validation below filters out SST-only and Scaffold-only intermediates.
        folders = {
            Path(item).expanduser() for item in (
                workflow.get("root"), self.structure_root)
            if item
        }
        for filename in tuple(candidates):
            if filename:
                folders.add(Path(filename).expanduser().parent)
        for folder in folders:
            try:
                candidates.extend(folder.glob("*.json"))
            except OSError:
                continue
        valid = []
        seen = set()
        for filename in candidates:
            if not filename:
                continue
            path = Path(filename).expanduser().resolve()
            if path in seen or not path.is_file():
                continue
            seen.add(path)
            try:
                report = validate_structure(str(path), require_staples=True)
            except Exception:
                continue
            if report.get("capture_bridge_component_count", 0) > 0:
                valid.append((path.stat().st_mtime_ns, str(path)))
        return max(valid)[1] if valid else None

    def _structure_preview_sources(self):
        if self.project is None:
            return []
        workflow = self._workflow()
        sources = []

        def append(channel, label, filename, accepted=False,
                   generated_only=False):
            if filename and Path(filename).is_file():
                state = ("已接受" if accepted else
                         "已生成" if generated_only else "待审核")
                sources.append({
                    "channel": channel,
                    "label": translate(label) + " · " + translate(state),
                    "filename": str(Path(filename).resolve()),
                    "accepted": bool(accepted),
                })

        sst = (workflow.get("sst_accepted") or
               workflow.get("sst_review") or
               workflow.get("sst_two_layer"))
        append("sst", "单独 SST", sst,
               bool(workflow.get("sst_accepted")), generated_only=True)
        scaffold = (workflow.get("scaffold_review") or
                    workflow.get("scaffold_accepted"))
        append("scaffold", "Scaffold + SST sublattice", scaffold,
               bool(scaffold == workflow.get("scaffold_accepted") and
                    not workflow.get("scaffold_editing")),
               generated_only=not bool(workflow.get("scaffold_accepted")))
        complete = workflow.get("structure_complete")
        append("complete", "Scaffold + Staple + Capture + SST sublattice",
               complete, bool(workflow.get("structure_accepted")))
        return sources

    def _set_structure_preview_report(self, text):
        text = str(text or "").strip()
        self.structure_preview_status.setText(text)
        self.structure_preview_status.setVisible(bool(text))
        # The report now occupies a full-width bottom row; clear the legacy
        # report painted inside the Path viewport.
        self.capture_preview.set_path_report("")
        if text:
            height = max(1, self.capture_results_splitter.height())
            report_height = max(90, min(180, height // 5))
            self.capture_results_splitter.setSizes(
                [max(1, height - report_height), report_height])

    def _refresh_structure_preview(self, preferred=None):
        if not hasattr(self, "structure_preview_channel"):
            return
        previous = self.structure_preview_channel.currentData()
        previous_channel = (previous.get("channel")
                            if isinstance(previous, dict) else None)
        sources = self._structure_preview_sources()
        blocker = QSignalBlocker(self.structure_preview_channel)
        self.structure_preview_channel.clear()
        if not sources:
            self.structure_preview_channel.addItem(
                "尚无可显示通道", None)
            self.structure_preview_channel.setEnabled(False)
            self.capture_preview.clear()
            self.capture_preview.setMinimumHeight(0)
            # Before a design exists, show no yellow placeholder report.
            self._set_structure_preview_report("")
            del blocker
            return
        self.structure_preview_channel.setEnabled(True)
        for source in sources:
            self.structure_preview_channel.addItem(
                source["label"], source)
        available = [item["channel"] for item in sources]
        selected = (preferred if preferred in available else
                    previous_channel if previous_channel in available else
                    available[-1])
        self.structure_preview_channel.setCurrentIndex(
            available.index(selected))
        del blocker
        self._structure_preview_channel_changed()

    def _structure_preview_channel_changed(self, unused_index=None):
        if not hasattr(self, "structure_preview_channel"):
            return
        source = self.structure_preview_channel.currentData()
        if not isinstance(source, dict):
            self.capture_preview.clear()
            self.capture_preview.setMinimumHeight(0)
            self._set_structure_preview_report("")
            return
        try:
            self.capture_preview.set_source(
                source["filename"], source["channel"], source["label"])
            self.capture_preview.setMinimumHeight(0)
        except Exception as error:
            self.capture_preview.clear("设计预览无法读取。")
            self.capture_preview.setMinimumHeight(0)
            self._set_structure_preview_report(
                "预览读取失败：%s" % error)
            return
        detail = ""
        if source["channel"] == "complete":
            try:
                report = validate_structure(
                    source["filename"], require_staples=True)
                payload = json.loads(Path(source["filename"]).read_text(
                    encoding="utf-8"))
                layout = payload.get("moire_structure_metadata", {}).get(
                    "variable_length_layout", {})
                z1 = report.get("seed_z1_overlap_bp")
                z2 = report.get("seed_z2_actual_bp")
                z3 = report.get("seed_z3_overlap_bp")
                lengths = report.get("seed_scaffold_lengths", [])
                result = translate(
                    "实际 Z1 / Z2 / Z3：%s / %s / %s bp\n"
                    "Scaffold：%d条 · %s nt\n%s" % (
                        _format_bp(z1), _format_bp(z2), _format_bp(z3),
                        len(lengths), " / ".join(map(str, lengths)) or "—",
                        self._staple_analysis_text(report)))
                detail = result
            except Exception as error:
                detail = translate("统计读取失败：%s" % error)
        report_text = Path(source["filename"]).name
        if detail:
            report_text += "\n" + detail
        self._set_structure_preview_report(report_text)

    def _restore_structure_workflow(self):
        if not hasattr(self, "scaffold_status") or self.project is None:
            return
        workflow = self.project.seed_plan.get("structure_workflow", {})
        current_signature = self._structure_signature(self.project.settings)
        if (workflow.get("settings_signature") and
                workflow.get("settings_signature") != current_signature):
            had_downstream = any(workflow.get(key) for key in (
                "sst_review", "scaffold_review", "scaffold_accepted",
                "structure_complete", "structure_accepted",
                "sequence_source", "sequence_assignments"))
            self._invalidate_downstream_design_state(workflow)
            workflow["parameters_accepted"] = False
            workflow.pop("parameters_accepted_at", None)
            workflow["parameters_editing"] = True
            notice_key = json.dumps(current_signature, sort_keys=True)
            if had_downstream:
                workflow["stale_notice"] = notice_key
                QMessageBox.warning(
                    self, "Design parameters changed — regeneration required",
                    "The updated lengths, insertion/deletion settings, or "
                    "other structural parameters no longer match the "
                    "previous DNA design. Existing files remain on disk, "
                    "but Steps 2 and 3 have been reset. Accept the updated "
                    "parameters, then regenerate and accept each stage in "
                    "order.")
        basis_ok = bool(workflow.get("design_basis_accepted"))
        parameters_ok = bool(workflow.get("parameters_accepted"))
        self.design_basis_next_button.setEnabled(basis_ok)
        self.accept_parameters_button.setEnabled(
            len(self.seed_cross_section_picker.cells()) >= 4)
        self.parameters_next_button.setEnabled(parameters_ok)
        self._set_acceptance_button(
            self.confirm_design_basis_button, basis_ok,
            "4. 接受当前对称性与 Twist",
            "✓ 1.1 已接受")
        self._set_acceptance_button(
            self.accept_parameters_button, parameters_ok,
            "接受当前 Moiré 参数",
            "✓ 当前 Moiré 参数已接受")
        # The former 1.1/1.2 split workflow displayed a second green status
        # card.  The combined page has only the final Moiré acceptance card.
        self.design_basis_action_status.hide()
        if parameters_ok:
            self._show_action_feedback(
                self.parameters_action_status,
                "Moiré 参数已接受。下一步：打开 Automated DNA Design "
                "并生成三个设计文件。")
        else:
            self.parameters_action_status.hide()
        self.scaffold_status.setText(
            "等待生成；固定两层SST会在后台同时加入。")
        self.structure_status.setText(
            "接受scaffold后才可生成。capture每个pair使用同一种颜色。")
        self.generate_scaffold_button.setEnabled(parameters_ok)
        self.generate_simple_design_button.setEnabled(parameters_ok)
        self.open_scaffold_button.setEnabled(False)
        self.accept_scaffold_button.setEnabled(False)
        self.generate_structure_button.setEnabled(False)
        self.open_structure_button.setEnabled(False)
        self.accept_structure_button.setEnabled(False)
        self.design_generation_action_status.hide()
        self.cadnano_edit_action_status.hide()
        self.structure_accept_action_status.hide()
        self.structure_root = (Path(workflow["root"])
                               if workflow.get("root") else None)
        sst_accepted = workflow.get("sst_accepted")
        scaffold = workflow.get("scaffold_review")
        accepted = workflow.get("scaffold_accepted")
        complete = workflow.get("structure_complete")
        if parameters_ok:
            cfg = self.project.settings
            partition = self.project.prediction.get(
                "preview_seed_partition", {})
            seed_z1 = int(partition.get(
                "sst_overlap_z1_bp", cfg.growth_bp_z1))
            seed_z3 = int(partition.get(
                "sst_overlap_z3_bp", cfg.growth_bp_z3))
            period = self.project.prediction.get("predicted_moire_period_nm")
            accepted_angle = self.project.prediction.get(
                "reported_angle_deg", 0.0)
            self.accepted_parameters_summary.setText(
                translate(
                    "Accepted parameters · twist %+.1f° (%s) · moiré "
                    "period %s · "
                    "SST sublattice %d / %d / %d bp · "
                    "seed %d / %d / %d bp") % (
                    accepted_angle, _twist_handedness(accepted_angle),
                    ("—" if period is None else
                     "∞" if not math.isfinite(period) else "%.1f nm" % period),
                    cfg.sst_growth_bp_z1, cfg.spacer_bp_z2,
                    cfg.sst_growth_bp_z3, seed_z1,
                    cfg.spacer_bp_z2, seed_z3))
            self.seed_scaffold_capacity_status.setText("")
        else:
            self.accepted_parameters_summary.setText(
                translate("暂无设计参数"))
            if not workflow.get("scaffold_capacity_precheck"):
                self.seed_scaffold_capacity_status.setText(
                    translate("接受参数时按合法routing精确核算"))
        if sst_accepted and not workflow.get("sst_dependency_fingerprint"):
            try:
                workflow["sst_dependency_fingerprint"] = \
                    self._dependency_fingerprint(sst_accepted, "sst")
            except Exception:
                pass
        if accepted and not workflow.get("scaffold_dependency_fingerprint"):
            try:
                workflow["scaffold_dependency_fingerprint"] = \
                    self._dependency_fingerprint(accepted, "scaffold")
            except Exception:
                pass
        if scaffold:
            self.scaffold_status.setText(
                ("重新审核：" if workflow.get("scaffold_editing") else
                 ("已接受：" if accepted else "待审核："))+
                Path(scaffold).name)
        if complete:
            self.structure_status.setText(
                ("已接受：" if workflow.get("structure_accepted") else
                 "待审核：")+Path(complete).name)
            self._show_action_feedback(
                self.design_generation_action_status,
                "The three design files are available. Next: optionally "
                "inspect the final JSON in caDNAno, then accept it.")
        self.open_scaffold_button.setEnabled(bool(scaffold))
        self.accept_scaffold_button.setEnabled(bool(
            scaffold and (not accepted or workflow.get("scaffold_editing"))))
        self.generate_structure_button.setEnabled(bool(accepted))
        self.open_structure_button.setEnabled(bool(complete))
        self.inspect_final_design_button.setEnabled(bool(
            complete and Path(complete).is_file()))
        self.accept_structure_button.setEnabled(bool(
            complete and not workflow.get("structure_accepted")))
        structure_ok = bool(workflow.get("structure_accepted"))
        self.structure_next_button.setEnabled(structure_ok)
        self.detect_scaffold_sequences_button.setEnabled(structure_ok)
        self._set_acceptance_button(
            self.accept_structure_button, structure_ok, "接受当前设计图",
            "✓ 当前设计图已接受")
        if structure_ok and complete:
            accepted_at = workflow.get("structure_accepted_at", "—")
            self.accepted_design_summary.setText(
                "已接受设计图：%s\n接受时间：%s" %
                (Path(complete).name, accepted_at))
            self._show_action_feedback(
                self.structure_accept_action_status,
                "Final design accepted. Next: assign sequences in "
                "Sequence Export.")
        else:
            self.accepted_design_summary.setText("尚未接受设计图")
        if hasattr(self, "sequence_export_status"):
            # Always rebuild the sequence UI from the persisted workflow.
            # Otherwise cards left in the Qt layouts from a previous final
            # export can survive even after their workflow data is invalid.
            self._sequence_analysis = None
            self._clear_sequence_cards(self.scaffold_cards_layout)
            self.scaffold_cards_layout.addWidget(QLabel(
                "The structure has not been read."))
            self._clear_sequence_cards(self.sst_cards_layout)
            self.sst_cards_layout.addWidget(QLabel(
                "SST sublattice input positions have not been detected."))
            self.sequence_preview.clear()
            self.sequence_preview_status.setText(
                "Detect an accepted design to display its sequence routes.")
            scaffold_name_aliases = {"CS3": "CS3L", "CS4-L": "CS4"}
            self._sequence_assignments = {}
            for saved_item in workflow.get("sequence_assignments", []):
                if not saved_item.get("target_id"):
                    continue
                item = dict(saved_item)
                old_name = item.get("scaffold_name")
                if old_name in scaffold_name_aliases:
                    item["scaffold_name"] = scaffold_name_aliases[old_name]
                self._sequence_assignments[item["target_id"]] = item
            self.sst_auto_import_status.hide()
            self.sst_expert_import_status.hide()
            self.sst_acceptance_status.hide()
            scaffold_accepted = bool(
                workflow.get("sequence_scaffold_accepted"))
            self.scaffold_detection_action_status.hide()
            self.scaffold_sequence_status.hide()
            if scaffold_accepted:
                self._show_action_feedback(
                    self.scaffold_sequence_status,
                    "Scaffold sequences accepted. Next: detect the SST "
                    "sublattice input positions and lengths.")
            self.detect_sst_inputs_button.setEnabled(scaffold_accepted)
            self.export_sst_input_template_button.setEnabled(False)
            self.import_sst_input_template_button.setEnabled(False)
            self.accept_added_scaffold_button.setEnabled(False)
            self._set_acceptance_button(
                self.accept_added_scaffold_button, scaffold_accepted,
                "Accept assigned scaffold sequences",
                "✓ Scaffold sequences accepted")
            self._set_acceptance_button(
                self.accept_added_sst_button, False,
                "Accept assigned SST sublattice input sequences",
                "✓ SST sublattice input sequences accepted")
            self.auto_design_sst_inputs_button.setEnabled(False)
            self.sequence_expert_button.setEnabled(True)
            self.sst_detection_status.setText(
                "SST sublattice inputs have not been detected."
                if scaffold_accepted else
                "Please accept the scaffold sequences first.")

            design = (str(Path(complete).resolve())
                      if structure_ok and complete and
                      Path(complete).is_file() else None)
            if scaffold_accepted and design:
                try:
                    if (not self._sequence_analysis or
                            self._sequence_analysis.get("path") != design):
                        self._sequence_analysis = analyze_sequence_design(
                            design)
                except Exception as error:
                    self.sst_detection_status.setText(
                        "The saved SST sublattice input analysis could not "
                        "be restored: "
                        "%s" % error)
                if self._sequence_analysis:
                    self._render_scaffold_cards()
                    has_sst_assignments = any(
                        str(item.get("category", "")).startswith(
                            "sst_input_layer_")
                        for item in self._sequence_assignments.values())
                    restore_sst = bool(
                        workflow.get("sequence_sst_detected") or
                        has_sst_assignments or
                        workflow.get("sequence_sst_accepted"))
                    if restore_sst:
                        self.sequence_preview.set_source(
                            design, "sst_input",
                            "SST sublattice input only")
                        layers = ((1,) if self._sst_layers_identical()
                                  else (1, 2))
                        targets = [
                            item for layer in layers for item in
                            self._sequence_analysis["targets"].get(
                                "sst_input_layer_%d" % layer, [])]
                        self.sequence_preview.set_sst_input_targets(targets)
                        self._render_sst_cards()
                        self._show_action_feedback(
                            self.sst_detection_status,
                            workflow.get(
                                "sequence_sst_detection_status") or
                            "SST sublattice input positions and lengths "
                            "restored from the saved project. Next: design "
                            "and assign the input sequences.")
                        self.auto_design_sst_inputs_button.setEnabled(
                            bool(targets))
                        self.sequence_expert_button.setEnabled(True)
                        self._update_sequence_expert_actions()

            import_status = workflow.get("sequence_sst_import_status")
            import_method = workflow.get("sequence_sst_import_method")
            if not import_status and self._sequence_analysis and \
                    self._sst_assignments_complete():
                import_status = (
                    "SST sublattice input sequences assigned successfully.")
                assignment_source = next((
                    str(item.get("source", "")) for item in
                    self._sequence_assignments.values()
                    if str(item.get("category", "")).startswith(
                        "sst_input_layer_") and item.get("source")), "")
                import_method = (
                    "automatic" if "automatic_orthogonal_input" in
                    assignment_source else "expert")
            if import_status:
                target_label = (
                    self.sst_expert_import_status
                    if import_method == "expert" else
                    self.sst_auto_import_status)
                self._show_sequence_status(target_label, import_status)
                if import_method == "expert" and not \
                        self.sequence_expert_button.isChecked():
                    self.sst_expert_import_status.hide()

            sequence_source = workflow.get("sequence_source")
            ready = bool(workflow.get("sequence_sst_accepted") and
                         sequence_source and Path(sequence_source).is_file())
            self.final_sequence_export_button.setEnabled(ready)
            accepted_input = bool(workflow.get("sequence_sst_accepted"))
            if accepted_input:
                self._show_action_feedback(
                    self.scaffold_sequence_status,
                    "Accepted scaffold sequences restored. Next: review "
                    "or accept the restored SST sublattice input sequences.")
                self._show_sequence_status(
                    self.sst_acceptance_status,
                    workflow.get("sequence_sst_acceptance_status") or
                    "Assigned SST sublattice inputs accepted. The nucleotide "
                    "sequences have been written to both layers.")
                if ready:
                    self.sequence_export_status.setText(
                        "Sequence assignment is complete: %s. The final "
                        "package can be exported again." %
                        Path(sequence_source).name)
                self._set_acceptance_button(
                    self.accept_added_sst_button, True,
                    "Accept assigned SST sublattice input sequences",
                    "✓ SST sublattice input sequences accepted")
            self._update_sequence_expert_actions()
        self._refresh_structure_preview()

    def _save_current_project(self, silent=False):
        if self._app_mode != "design":
            return False
        if self.project is None:
            self.recalculate()
        if self.project is None or not self.project_path:
            return False
        try:
            self.project_path = str(save_project(
                self.project, self.project_path))
        except Exception as error:
            if not silent:
                QMessageBox.critical(self, "项目保存失败", str(error))
            return False
        self._update_current_project_display()
        if not silent:
            self.statusBar().showMessage(
                "项目已保存：%s" % self.project_path, 8000)
        return True

    def _autosave_project(self):
        if self._app_mode != "design":
            return
        if self.project_path and self.project is not None:
            if not self._save_current_project(silent=True):
                self.statusBar().showMessage("项目自动保存失败。", 6000)

    def save_project(self):
        if not self.project_path:
            return self.save_project_as()
        return self._save_current_project(silent=False)

    def save_project_as(self):
        selection = self._project_setup_selection(
            "Moiré 项目另存为", "另存为")
        if not selection:
            return False
        name, filename, unused_language = selection
        self.project_name.setText(name)
        self.project_path = filename
        self.structure_root = None
        self.recalculate()
        workflow = self._workflow()
        workflow["root"] = str(
            self._project_directory()/"cadnano design")
        saved = self._save_current_project(silent=False)
        if saved:
            self.statusBar().showMessage(
                "当前项目已另存为：%s" % self.project_path, 8000)
        return saved

    def open_project(self):
        filename, unused = QFileDialog.getOpenFileName(
            self, "打开 Moiré 工程", str(Path.home()/"Desktop"),
            "Moiré project (*.moire.json)")
        if not filename:
            return False
        try:
            project = load_project(filename)
        except Exception as error:
            QMessageBox.critical(self, "无法打开工程", str(error))
            return False
        self._sequence_analysis = None
        self._sequence_assignments = {}
        self.project = project
        self.project_path = str(Path(filename).resolve())
        project.settings.interface_language = current_language()
        self._load_settings(project.settings)
        self.recalculate()
        self.project.measurements = project.measurements
        self._update_measurements()
        self._update_current_project_display()
        self._restore_structure_workflow()
        self._record_history("打开项目")
        return True

    def _load_settings(self, cfg):
        self._updating = True
        blockers = [QSignalBlocker(widget) for widget in (
            self.target_definition, self.angle, self.period,
            self.lattice_context, self.lattice_constant,
            self.lattice_constant_2, self.mean_indel,
            self.layers_identical, self.bilayer_symmetry_selector,
            self.seed_cross_section_picker)]
        try:
            self.project_name.setText(cfg.project_name)
            self.bilayer_symmetry_selector.setCurrentIndex(max(
                0, self.bilayer_symmetry_selector.findData(
                    getattr(cfg, "lattice_symmetry", "square_square_c4"))))
            # Populate the choices for the loaded lattice before mapping the
            # serialized cells back to their named preset.  Otherwise a
            # Kagome rectangular preset could be replaced by the 8×8 default.
            self._populate_cross_section_presets()
            self.seed_cross_section_picker.set_cells(getattr(
                cfg, "seed_cross_section_cells",
                [[row, col] for row in range(8) for col in range(8)
                 if not (2 <= row <= 5 and 2 <= col <= 5)]))
            preset_key = self._preset_key_for_cells(
                self.seed_cross_section_picker.cells())
            if preset_key is not None:
                preset_index = self.seed_cross_section_preset.findData(
                    preset_key)
                if preset_index >= 0:
                    self.seed_cross_section_preset.setCurrentIndex(
                        preset_index)
            # The staged designer exposes only the local growth-surface
            # angle.  Legacy projects using the former experimental-total
            # selector are normalized on load so a hidden setting can never
            # silently alter the prediction.
            self.target_definition.setCurrentIndex(0)
            angle_value = float(cfg.target_angle_deg)
            if getattr(cfg, "target_definition", "local_surface") != \
                    "local_surface" and self.project is not None:
                angle_value = float(self.project.prediction.get(
                    "predicted_local_surface_angle_deg", angle_value))
            self.angle.setValue(angle_value)
            if getattr(cfg, "lattice_symmetry", "square_square_c4") != \
                    "square_kagome" and math.isfinite(cfg.target_period_nm):
                self.period.setValue(cfg.target_period_nm)
            self.lattice_constant.setValue(getattr(
                cfg, "layer1_lattice_constant_nm", cfg.lattice_constant_nm))
            self.lattice_constant_2.setValue(getattr(
                cfg, "layer2_lattice_constant_nm", cfg.lattice_constant_nm))
            context = getattr(cfg, "lattice_context", "solution_cryo")
            context_index = next((
                index for index in range(self.lattice_context.count())
                if self.lattice_context.itemData(index)[0] == context), 0)
            self.lattice_context.setCurrentIndex(context_index)
            self.mean_indel.setValue(getattr(cfg, "mean_indel_per_helix", 0.0))
            identical = getattr(
                cfg, "layers_design_sequence_identical", True)
            self.layers_identical.setCurrentIndex(0 if identical else 1)
            sst_z1 = getattr(cfg, "sst_growth_bp_z1", cfg.growth_bp_z1)
            sst_z3 = getattr(cfg, "sst_growth_bp_z3", cfg.growth_bp_z3)
            self._fill_length_combo(
                self.sst_z1, self._all_growth_values(), sst_z1)
            z2_values = ([item for item in compatible_z2_values(
                int(self.sst_z1.currentData()), maximum=160)
                          if item <= 160]
                         if identical else self._all_z2_values())
            self._fill_length_combo(
                self.sst_spacing, z2_values, cfg.spacer_bp_z2)
            self._fill_length_combo(
                self.seed_z2, z2_values, cfg.spacer_bp_z2)
            growth_values = (compatible_growth_values(
                int(self.sst_spacing.currentData()), maximum=400)
                             if identical else self._all_growth_values())
            self._fill_length_combo(
                self.sst_z3, growth_values,
                sst_z1 if identical else sst_z3)
            self.seed_z1.setValue(128)
            self.seed_z3.setValue(128)
        finally:
            del blockers
            self._updating = False
        self._target_driver = "angle"
        self._update_phase_hint()
        self._apply_symmetry_ui()
        self.design_stack.setCurrentIndex(1)

    def export_project(self):
        if self.project is None:
            self.recalculate()
        parent = QFileDialog.getExistingDirectory(
            self, "选择原型项目保存位置", str(Path.home()/"Desktop"))
        if not parent:
            return
        root = Path(parent)/(self.project.settings.project_name)
        if root.exists():
            suffix = 1
            while Path(str(root)+"_%d" % suffix).exists():
                suffix += 1
            root = Path(str(root)+"_%d" % suffix)
        design = root/"design"
        sequences = root/"SST sublattice input"
        design.mkdir(parents=True)
        sequences.mkdir(parents=True)
        project_file = save_project(
            self.project, str(root/(self.project.settings.project_name+".moire.json")))
        export_capture_map(self.project, str(sequences/"capture_map.csv"))
        seed_message = ""
        workflow = self.project.seed_plan.get("structure_workflow", {})
        copied = []
        for key in ("scaffold_accepted", "structure_complete"):
            source = workflow.get(key)
            if source and Path(source).is_file():
                target = design/Path(source).name
                shutil.copy2(source, target)
                copied.append(target.name)
        if not copied:
            try:
                export_reference_seed(
                    self.project,
                    str(design/(self.project.settings.project_name+"_seed.json")))
            except Exception as error:
                seed_message = "\n参考Seed JSON未导出：%s" % error
        readme = root/"README.txt"
        readme.write_text(
            translate("DNA Moiré Designer prototype\n\n"
            "This project uses the calibrated S8-R4x4C Square-Square bilayer model.\n"
            "Structure workflow: the fixed two-layer SST sublattice is "
            "generated "
            "internally together with reviewed scaffold routing, followed "
            "by staple/capture generation and dual-state sequence export.\n",
            getattr(self.project.settings, "interface_language", "en")),
            encoding="utf-8")
        QMessageBox.information(
            self, "导出完成", "原型项目已保存到：\n%s%s" % (root, seed_message))

    def open_in_cadnano(self):
        if self.project is None:
            self.recalculate()
        if not CADNANO_EXECUTABLE.is_file():
            QMessageBox.warning(self, "caDNAno", "未找到当前 cadnano 启动程序。")
            return
        workflow = self.project.seed_plan.get("structure_workflow", {})
        preferred = (workflow.get("structure_complete") or
                     workflow.get("scaffold_accepted") or
                     workflow.get("scaffold_review"))
        if preferred and Path(preferred).is_file():
            subprocess.Popen([str(CADNANO_EXECUTABLE), str(preferred)],
                             start_new_session=True,
                             stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL)
            self.statusBar().showMessage("已将当前结构发送到 cadnano", 8000)
            return
        temp_root = Path.home()/"Desktop"/(self.project.settings.project_name+"_cadnano")
        suffix = 1
        candidate = temp_root
        while candidate.exists():
            candidate = Path(str(temp_root)+"_%d" % suffix); suffix += 1
        candidate.mkdir(parents=True)
        try:
            seed = export_reference_seed(
                self.project, str(candidate/(self.project.settings.project_name+"_seed.json")))
            save_project(self.project, str(candidate/(self.project.settings.project_name+".moire.json")))
            subprocess.Popen([str(CADNANO_EXECUTABLE), str(seed)],
                             start_new_session=True,
                             stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL)
        except Exception as error:
            QMessageBox.critical(self, "无法打开 cadnano", str(error))
            return
        self.statusBar().showMessage("已将参考 Seed S 发送到 cadnano", 8000)

    def _legacy_analysis_mode_changed(self):
        mode = self.image_analysis_mode.currentData()
        self.select_tem_button.setEnabled(mode in ("tem", "combined"))
        self.select_fft_button.setEnabled(mode in ("fft", "combined"))
        if hasattr(self, "_tem_analysis"):
            self._analysis_values_changed()

    def _legacy_select_tem_image(self):
        filename, unused = QFileDialog.getOpenFileName(
            self, "选择TEM图", str(Path.home()/"Desktop"),
            "Image (*.png *.jpg *.jpeg *.tif *.tiff *.bmp)")
        if filename:
            self._tem_image_path = str(Path(filename).resolve())
            self.tem_path_label.setText(Path(filename).name)

    def select_fft_image(self):
        filename, unused = QFileDialog.getOpenFileName(
            self, "选择FFT图", str(Path.home()/"Desktop"),
            "Image (*.png *.jpg *.jpeg *.tif *.tiff *.bmp)")
        if filename:
            self._fft_image_path = str(Path(filename).resolve())
            self.fft_path_label.setText(Path(filename).name)

    @staticmethod
    def _image_to_pgm(source, target):
        image = QImage(source)
        if image.isNull():
            raise ValueError("Qt无法读取图像：%s" % source)
        gray = image.convertToFormat(QImage.Format.Format_Grayscale8)
        pointer = gray.constBits()
        pointer.setsize(gray.sizeInBytes())
        raw = bytes(pointer)
        rows = b"".join(
            raw[y * gray.bytesPerLine():y * gray.bytesPerLine() + gray.width()]
            for y in range(gray.height()))
        Path(target).write_bytes(
            ("P5\n%d %d\n255\n" % (gray.width(), gray.height())).encode("ascii") +
            rows)

    @staticmethod
    def _legacy_run_image_worker(mode, source):
        with tempfile.TemporaryDirectory(prefix="moire-image-analysis-") as folder:
            pgm = Path(folder)/"image.pgm"
            MoireDesignerWindow._image_to_pgm(source, pgm)
            command = worker_command(
                "image-analysis", mode, str(pgm),
                "--original", str(source))
            result = subprocess.run(
                command, check=False, text=True, capture_output=True,
                timeout=180)
            if result.returncode:
                raise RuntimeError(result.stderr.strip() or
                                   "图像分析进程失败。")
            try:
                return json.loads(result.stdout)
            except Exception as error:
                raise RuntimeError("图像分析没有返回有效结果。") from error

    @staticmethod
    def _unique_analysis_path(source, suffix):
        source = Path(source)
        target = source.with_name(source.stem + suffix)
        counter = 1
        while target.exists():
            target = source.with_name(
                source.stem + suffix.rsplit(".", 1)[0] +
                "_%d." % counter + suffix.rsplit(".", 1)[1])
            counter += 1
        return target

    def _legacy_run_image_analysis(self):
        mode = self.image_analysis_mode.currentData()
        if mode in ("tem", "combined") and not self._tem_image_path:
            QMessageBox.information(self, "缺少TEM图", "请先上传TEM图。")
            return
        if mode in ("fft", "combined") and not self._fft_image_path:
            QMessageBox.information(self, "缺少FFT图", "请先上传FFT图。")
            return
        self._busy(True)
        messages = []
        try:
            self._tem_analysis = (self._run_image_worker(
                "tem", self._tem_image_path)
                if mode in ("tem", "combined") else None)
            self._fft_analysis = (self._run_image_worker(
                "fft", self._fft_image_path)
                if mode in ("fft", "combined") else None)
        except Exception as error:
            QMessageBox.critical(self, "图像分析失败", str(error))
            return
        finally:
            self._busy(False)

        blockers = [QSignalBlocker(widget) for widget in (
            self.scale_bar_pixels, self.scale_bar_nm,
            self.measured_lattice, self.measured_moire)]
        try:
            if self._tem_analysis:
                bar = self._tem_analysis.get("scale_bar") or {}
                self.scale_bar_pixels.setValue(float(
                    bar.get("pixel_length") or 0.0))
                self.scale_bar_nm.setValue(float(
                    self._tem_analysis.get("scale_value_nm") or 0.0))
                self.measured_lattice.setValue(float(
                    self._tem_analysis.get("lattice_constant_nm") or 0.0))
                self.measured_moire.setValue(float(
                    self._tem_analysis.get("moire_period_nm") or 0.0))
                if not bar:
                    messages.append("未可靠识别scale bar横线，请手动输入像素长度。")
                if not self._tem_analysis.get("scale_value_nm"):
                    messages.append("OCR未读出scale bar数值，请手动输入nm值。")
                if not self._tem_analysis.get("moire_period_px"):
                    messages.append("TEM频谱没有稳定识别到moiré周期峰。")
                lattice_fft = self._tem_analysis.get("lattice_fft") or {}
                real_space = self._tem_analysis.get("moire_real_space") or {}
                if not lattice_fft.get("valid"):
                    messages.append(lattice_fft.get(
                        "error", "TEM的FFT没有稳定识别两套一阶Square峰。"))
                if not real_space.get("valid"):
                    messages.append(real_space.get(
                        "error", "TEM实空间没有识别到足够的moiré单元。"))
                difference = self._tem_analysis.get(
                    "moire_consistency_percent")
                if difference is not None:
                    if abs(float(difference)) <= 10.0:
                        messages.append(
                            "TEM中心距与FFT预测period匹配（差异%+.1f%%）。" %
                            float(difference))
                    else:
                        messages.append(
                            "TEM中心距与FFT预测period不匹配（差异%+.1f%%）；"
                            "请检查一阶峰选择或手动校正scale bar。" %
                            float(difference))
            if self._fft_analysis:
                angle = self._fft_analysis.get("twist_angle_deg")
                if angle is None:
                    messages.append(self._fft_analysis.get(
                        "error", "FFT没有稳定识别到两套点阵。"))
        finally:
            del blockers
        self._analysis_output_paths = {}
        if self._tem_image_path and self._tem_analysis:
            self._analysis_output_paths["tem"] = str(
                self._unique_analysis_path(
                    self._tem_image_path, "_moire_analysis.png"))
        if self._fft_image_path and self._fft_analysis:
            self._analysis_output_paths["fft"] = str(
                self._unique_analysis_path(
                    self._fft_image_path, "_moire_analysis.png"))
        self.image_analysis_status.setText(
            "自动分析完成。" + ("\n" + "\n".join(messages) if messages else
                           " 结果可人工校正。"))
        self._analysis_values_changed()
        self.save_image_analysis_button.setEnabled(
            self._final_analysis_angle is not None)

    def _legacy_scale_calibration_changed(self):
        if not self._tem_analysis:
            return
        pixels = self.scale_bar_pixels.value()
        value_nm = self.scale_bar_nm.value()
        if pixels <= 0 or value_nm <= 0:
            return
        pixel_size = value_nm / pixels
        blockers = [QSignalBlocker(self.measured_lattice),
                    QSignalBlocker(self.measured_moire)]
        try:
            lattice_px = self._tem_analysis.get("lattice_constant_px")
            moire_px = self._tem_analysis.get("moire_period_px")
            if lattice_px:
                self.measured_lattice.setValue(lattice_px * pixel_size)
            if moire_px:
                self.measured_moire.setValue(moire_px * pixel_size)
        finally:
            del blockers
        self._analysis_values_changed()

    def _legacy_analysis_values_changed(self):
        lattice = self.measured_lattice.value()
        period = self.measured_moire.value()
        tem_angle = None
        if lattice > 0 and period >= lattice / 2.0:
            tem_angle = math.degrees(2.0 * math.asin(
                min(1.0, lattice / (2.0 * period))))
        fft_angle = (float(self._fft_analysis.get("twist_angle_deg"))
                     if self._fft_analysis and
                     self._fft_analysis.get("twist_angle_deg") is not None
                     else (float(self._tem_analysis.get("fft_twist_angle_deg"))
                           if self._tem_analysis and
                           self._tem_analysis.get("fft_twist_angle_deg")
                           is not None else None))
        self.tem_twist_result.setText(
            "—" if tem_angle is None else "%.3f°" % tem_angle)
        self.fft_twist_result.setText(
            "—" if fft_angle is None else "%.3f°" % fft_angle)
        mode = self.image_analysis_mode.currentData()
        # Combined analysis deliberately uses the TEM real-space period.  FFT
        # remains an independent orientation check because peak fitting can be
        # biased by windowing, strain, and finite field of view.
        if mode in ("tem", "combined") and tem_angle is not None:
            final = tem_angle
            authority = "TEM period"
        else:
            final = fft_angle
            authority = "FFT orientation" if fft_angle is not None else None
        self._final_analysis_angle = final
        self._tem_analysis_angle = tem_angle
        self._fft_analysis_angle = fft_angle
        self.final_twist_result.setText(
            "—" if final is None else "%.3f°（%s）" % (final, authority))
        self._refresh_analysis_previews()

    def _annotated_analysis_image(self, source, mode, result, target):
        image = QImage(source).convertToFormat(QImage.Format.Format_RGB32)
        if image.isNull():
            return None
        painter = QPainter(image)
        scale = max(1.0, min(image.width(), image.height()) / 900.0)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setFont(QFont("Arial", max(11, int(15 * scale)), QFont.Weight.Bold))
        lines = []
        if mode == "tem":
            bar = result.get("scale_bar") or {}
            if bar:
                painter.setPen(QPen(QColor("#ffd21f"), max(2, int(3 * scale))))
                painter.drawRect(int(bar["x0"]), int(bar["y0"]),
                                 max(1, int(bar["x1"] - bar["x0"])),
                                 max(2, int(bar["y1"] - bar["y0"] + 1)))
            if self.measured_lattice.value() > 0:
                lines.append("Lattice constant: %.3f nm" %
                             self.measured_lattice.value())
            if self.measured_moire.value() > 0:
                lines.append("TEM-derived moiré period: %.3f nm" %
                             self.measured_moire.value())
            if self._final_analysis_angle is not None:
                lines.append("Twist: %.3f deg" % self._final_analysis_angle)
            lattice_fft = result.get("lattice_fft") or {}
            real_space = result.get("moire_real_space") or {}
            pixel_size = (self.scale_bar_nm.value() /
                          self.scale_bar_pixels.value()
                          if self.scale_bar_nm.value() > 0 and
                          self.scale_bar_pixels.value() > 0 else None)
            if lattice_fft.get("valid"):
                lines.append("FFT-derived twist: %.3f deg" %
                             float(lattice_fft["twist_angle_deg"]))
                if pixel_size:
                    lines.append("FFT-derived moiré period: %.3f nm" %
                                 (float(lattice_fft[
                                     "predicted_moire_period_px"]) *
                                  pixel_size))
            difference = result.get("moire_consistency_percent")
            if difference is not None:
                lines.append("TEM−FFT period deviation: %+.2f%%" %
                             float(difference))
            if real_space.get("valid"):
                basis = real_space.get("basis_vectors_px") or []
                if len(basis) == 2:
                    first, second = basis
                    painter.setPen(QPen(QColor(77, 255, 126, 205),
                                        max(1, int(2 * scale))))
                    for center in real_space.get("centers", []):
                        cx, cy = float(center["x"]), float(center["y"])
                        polygon = QPolygonF([
                            QPointF(cx - first[0] / 2 - second[0] / 2,
                                    cy - first[1] / 2 - second[1] / 2),
                            QPointF(cx + first[0] / 2 - second[0] / 2,
                                    cy + first[1] / 2 - second[1] / 2),
                            QPointF(cx + first[0] / 2 + second[0] / 2,
                                    cy + first[1] / 2 + second[1] / 2),
                            QPointF(cx - first[0] / 2 + second[0] / 2,
                                    cy - first[1] / 2 + second[1] / 2),
                        ])
                        painter.drawPolygon(polygon)
                pair = real_space.get("representative_pair") or {}
                if pair:
                    first_point = pair.get("first") or {}
                    second_point = pair.get("second") or {}
                    x1, y1 = float(first_point.get("x", 0)), float(
                        first_point.get("y", 0))
                    x2, y2 = float(second_point.get("x", 0)), float(
                        second_point.get("y", 0))
                    painter.setPen(QPen(QColor("#ffcf33"),
                                        max(2, int(4 * scale))))
                    painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
                    radius = max(3, int(5 * scale))
                    painter.drawEllipse(QPointF(x1, y1), radius, radius)
                    painter.drawEllipse(QPointF(x2, y2), radius, radius)
                    if pixel_size:
                        label = "%.2f nm" % (float(
                            pair.get("distance_px", 0)) * pixel_size)
                        painter.setPen(QColor("#ffcf33"))
                        painter.drawText(QPointF((x1 + x2) / 2 + 8 * scale,
                                                 (y1 + y2) / 2 - 8 * scale),
                                         label)
        else:
            center = result.get("center") or [image.width()/2, image.height()/2]
            colors = (QColor("#ff5252"), QColor("#22d3ee"))
            for point in result.get("peaks", []):
                color = colors[int(point.get("layer", 0)) % 2]
                painter.setPen(QPen(color, max(2, int(3 * scale))))
                painter.drawEllipse(int(point["x"] - 6 * scale),
                                    int(point["y"] - 6 * scale),
                                    int(12 * scale), int(12 * scale))
                painter.drawLine(int(center[0]), int(center[1]),
                                 int(point["x"]), int(point["y"]))
            angle = result.get("twist_angle_deg")
            if angle is not None:
                lines.append("FFT lattice twist: %.3f deg" % angle)
        if lines:
            metrics = painter.fontMetrics()
            margin = int(14 * scale)
            line_height = metrics.height() + int(5 * scale)
            width = max(metrics.horizontalAdvance(line) for line in lines) + margin * 2
            height = line_height * len(lines) + margin
            painter.fillRect(8, 8, width, height, QColor(0, 0, 0, 175))
            painter.setPen(QColor("white"))
            for index, line in enumerate(lines):
                painter.drawText(8 + margin, 8 + margin +
                                 (index + 1) * line_height - 4, line)
        painter.end()
        image.save(str(target), "PNG")
        return QPixmap.fromImage(image)

    def _legacy_refresh_analysis_previews(self):
        if not hasattr(self, "_analysis_output_paths"):
            return
        for mode, source, result, label in (
                ("tem", self._tem_image_path, self._tem_analysis,
                 self.tem_analysis_image),
                ("fft", self._fft_image_path, self._fft_analysis,
                 self.fft_analysis_image)):
            target = self._analysis_output_paths.get(mode)
            if not source or not result or not target:
                continue
            pixmap = self._annotated_analysis_image(
                source, mode, result, target)
            if pixmap is not None:
                label.setPixmap(pixmap.scaled(
                    1000, 760, Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation))
                label.adjustSize()

    def _legacy_save_image_analysis(self):
        if self._final_analysis_angle is None:
            QMessageBox.information(self, "没有结果", "请先完成图像分析。")
            return
        self._refresh_analysis_previews()
        primary = self._tem_image_path or self._fft_image_path
        report_path = self._unique_analysis_path(
            primary, "_moire_analysis.json")
        mode = self.image_analysis_mode.currentData()
        report = {
            "mode": mode,
            "tem_source": self._tem_image_path,
            "fft_source": self._fft_image_path,
            "scale_bar_pixels": self.scale_bar_pixels.value(),
            "scale_bar_nm": self.scale_bar_nm.value(),
            "lattice_constant_nm": self.measured_lattice.value() or None,
            "moire_period_nm": self.measured_moire.value() or None,
            "tem_twist_deg": self._tem_analysis_angle,
            "fft_twist_deg": self._fft_analysis_angle,
            "final_twist_deg": self._final_analysis_angle,
            "final_authority": ("TEM period" if self._tem_analysis_angle
                                is not None and mode in ("tem", "combined") else
                                "FFT orientation"),
            "annotated_images": self._analysis_output_paths,
            "tem_automatic": self._tem_analysis,
            "fft_automatic": self._fft_analysis,
        }
        Path(report_path).write_text(
            json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        if self.project is None:
            self.recalculate()
        source = {"tem": "TEM automatic", "fft": "FFT automatic",
                  "combined": "TEM + FFT combined"}[mode]
        add_measurement(
            self.project, self._final_analysis_angle,
            self.measured_moire.value() or None, source)
        self.project.measurements[-1].update({
            "lattice_constant_nm": self.measured_lattice.value() or None,
            "analysis_report": str(report_path),
            "annotated_images": dict(self._analysis_output_paths),
        })
        self._update_measurements()
        QMessageBox.information(
            self, "分析已保存",
            "分析JSON：\n%s\n\n标注图：\n%s" %
            (report_path, "\n".join(self._analysis_output_paths.values())))

    def add_measurement(self):
        if self.project is None:
            self.recalculate()
        period = self.measure_period.value()
        add_measurement(
            self.project, self.measure_angle.value(),
            None if period <= 0 else period,
            self.measure_source.currentText())
        self._update_measurements()

    def _update_measurements(self):
        self.measurements_table.setRowCount(0)
        if self.project is None:
            return
        predicted = float(self.project.prediction["reported_angle_deg"])
        for measurement in self.project.measurements:
            row = self.measurements_table.rowCount()
            self.measurements_table.insertRow(row)
            values = [
                measurement["source"],
                "%.1f°" % measurement["angle_deg"],
                ("—" if measurement.get("period_nm") is None else
                 "%.1f nm" % measurement["period_nm"]),
                "%.1f°" % predicted,
                "%+.1f°" % measurement["prediction_error_deg"],
            ]
            for column, value in enumerate(values):
                self.measurements_table.setItem(row, column, QTableWidgetItem(value))


def create_window(parent=None, cadnano_controller=None):
    return MoireDesignerWindow(parent, cadnano_controller)
