"""Project-name and destination chooser for DNA Moiré Designer."""

from __future__ import annotations

from pathlib import Path

from PyQt6.QtWidgets import (
    QDialog, QDialogButtonBox, QFileDialog, QFormLayout,
    QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPushButton, QVBoxLayout)

from .i18n import UiLocalizer


def normalized_project_name(value):
    """Return a project name suitable for both settings and file names."""
    name = str(value or "").strip()
    if name.lower().endswith(".moire.json"):
        name = name[:-len(".moire.json")].strip()
    return name


class ProjectSetupDialog(QDialog):
    """Collect a project name and the folder that owns all project outputs."""

    def __init__(self, parent=None, title="新建 Moiré 项目",
                 project_name="moire_project", directory=None,
                 accept_text="创建项目", show_language=False):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setModal(True)
        self.setMinimumWidth(560)
        root = QVBoxLayout(self)
        intro = QLabel(
            "软件将在所选地址内创建“项目名”文件夹，并生成同名的 "
            ".moire.json 项目文件。结构、序列、最终导出和分析文件夹仅在"
            "首次产生对应文件时创建；项目状态每30秒自动保存。")
        intro.setWordWrap(True)
        root.addWidget(intro)

        form = QFormLayout()
        self.name_edit = QLineEdit(normalized_project_name(project_name))
        default_directory = Path(directory or (Path.home()/"Desktop"))
        if not default_directory.exists():
            default_directory = Path.home()
        self.directory_edit = QLineEdit(str(default_directory))
        browse_row = QHBoxLayout()
        browse_row.addWidget(self.directory_edit, 1)
        browse = QPushButton("选择文件夹…")
        browse.clicked.connect(self._browse)
        browse_row.addWidget(browse)
        form.addRow("项目名", self.name_edit)
        form.addRow("保存地址（项目文件夹的上一级）", browse_row)
        root.addLayout(form)

        self.preview = QLabel()
        self.preview.setWordWrap(True)
        root.addWidget(self.preview)
        self.name_edit.textChanged.connect(self._update_preview)
        self.directory_edit.textChanged.connect(self._update_preview)
        self._update_preview()

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save |
            QDialogButtonBox.StandardButton.Cancel)
        buttons.button(QDialogButtonBox.StandardButton.Save).setText(
            accept_text)
        buttons.accepted.connect(self._accept_if_valid)
        buttons.rejected.connect(self.reject)
        root.addWidget(buttons)
        self._localizer = UiLocalizer(self, interval_ms=100)

    def _browse(self):
        selected = QFileDialog.getExistingDirectory(
            self, "选择项目文件夹的保存地址", self.directory_edit.text())
        if selected:
            self.directory_edit.setText(selected)

    def _target_path(self):
        name = normalized_project_name(self.name_edit.text())
        parent = Path(self.directory_edit.text()).expanduser()
        return parent/name/(name+".moire.json") if name else parent

    def _update_preview(self):
        self.preview.setText("项目文件：%s" % self._target_path())

    def _accept_if_valid(self):
        name = normalized_project_name(self.name_edit.text())
        if not name:
            QMessageBox.warning(self, "项目名缺失", "请输入项目名。")
            return
        if any(character in name for character in ("/", "\\", ":")):
            QMessageBox.warning(
                self, "项目名不合法", "项目名不能包含 /、\\ 或 :。")
            return
        parent = Path(self.directory_edit.text()).expanduser()
        folder = parent/name
        try:
            folder.mkdir(parents=True, exist_ok=True)
        except Exception as error:
            QMessageBox.critical(self, "无法创建保存地址", str(error))
            return
        target = folder/(name+".moire.json")
        if target.exists():
            choice = QMessageBox.question(
                self, "项目已存在",
                "该项目文件已存在，是否覆盖？\n%s" % target,
                QMessageBox.StandardButton.Yes |
                QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No)
            if choice != QMessageBox.StandardButton.Yes:
                return
        self.name_edit.setText(name)
        self.accept()

    def selection(self):
        return (normalized_project_name(self.name_edit.text()),
                str(self._target_path().expanduser().resolve()),
                "en")
