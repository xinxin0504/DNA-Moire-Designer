"""Parameter dialog for the orthogonal-sequence designer."""

import os

from cadnano2.data.dnasequences import sequences as scaffold_sequences
from cadnano2.model.io.orthogonalseq import write_input_template
from cadnano2.views.primer3analysis import Primer3AnalysisDialog
import cadnano2.util as util

util.qtWrapImport('QtCore', globals(), ['Qt'])
util.qtWrapImport('QtGui', globals(), ['QIcon'])
util.qtWrapImport('QtWidgets', globals(), [
    'QCheckBox', 'QComboBox', 'QDialog', 'QDialogButtonBox', 'QDoubleSpinBox',
    'QFileDialog', 'QFormLayout', 'QHBoxLayout', 'QLabel', 'QLineEdit',
    'QMessageBox', 'QPushButton', 'QScrollArea', 'QSpinBox', 'QVBoxLayout',
    'QWidget'])


class OrthogonalSequenceDialog(QDialog):
    """Collect core rules and independently optional advanced rules."""

    def __init__(self, defaults, parent=None, primer3_entries=(),
                 suggested_directory=''):
        super(OrthogonalSequenceDialog, self).__init__(parent)
        self.setWindowTitle('正交序列设计')
        # The English labels and the three-part input-file control need a
        # wider initial viewport than the earlier compact Chinese dialog.
        # Keep a practical minimum, while allowing the form to wrap cleanly
        # on smaller displays instead of introducing horizontal scrolling.
        self.resize(1050, 760)
        self.setMinimumSize(900, 650)
        self._primer3Entries = tuple(primer3_entries)
        self._suggestedDirectory = suggested_directory or \
            os.path.expanduser('~/Desktop')

        outer = QVBoxLayout(self)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        content = QWidget(scroll)
        form = QFormLayout(content)
        form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)

        self.inputPath = QLineEdit(content)
        self.inputPath.setReadOnly(True)
        browse = QPushButton('选择 TXT…', content)
        template = QPushButton('导出输入模板…', content)
        inputWidget = QWidget(content)
        inputLayout = QHBoxLayout(inputWidget)
        inputLayout.setContentsMargins(0, 0, 0, 0)
        inputLayout.addWidget(self.inputPath, 1)
        inputLayout.addWidget(browse)
        inputLayout.addWidget(template)
        browse.clicked.connect(self._browseInput)
        template.clicked.connect(self._exportTemplate)

        self.length = self._spin(8, 300, defaults['length'])
        self.count = self._spin(1, 1000, defaults['count'])
        self.gcMin = self._double(0, 100, defaults['gc_min'] * 100, 1, ' %')
        self.gcMax = self._double(0, 100, defaults['gc_max'] * 100, 1, ' %')
        self.maxRun = self._spin(1, 12, defaults['max_homopolymer'])
        self.useSelfComplement, selfWidget = self._optionalWidget(
            defaults['use_self_complement'])
        self.selfComplement = self._spin(
            0, 30, defaults['max_self_complement'])
        self._finishOptional(
            selfWidget, self.useSelfComplement, (self.selfComplement,))

        self.useHairpin, hairpinWidget = self._optionalWidget(
            defaults['use_hairpin'])
        self.hairpin = self._spin(0, 30, defaults['max_hairpin_stem'])
        self._finishOptional(hairpinWidget, self.useHairpin, (
            QLabel('茎≤', content), self.hairpin))

        self.sameSubstring = self._spin(
            0, 30, defaults['max_same_substring'])
        self.crossComplement = self._spin(
            0, 30, defaults['max_cross_complement'])
        self.scaffold = QComboBox(content)
        self.scaffold.addItem('不使用骨架链', '')
        for name in sorted(scaffold_sequences):
            self.scaffold.addItem(
                '%s（%d nt）' % (name, len(scaffold_sequences[name])), name)
        self.scaffoldSame = self._spin(
            0, 30, defaults['scaffold_max_same_substring'])
        self.scaffoldCross = self._spin(
            0, 30, defaults['scaffold_max_cross_complement'])
        self.scaffoldSame.setEnabled(False)
        self.scaffoldCross.setEnabled(False)
        self.scaffold.currentIndexChanged.connect(
            lambda unused_index: self._setScaffoldControlsEnabled())
        self.pool = self._spin(1, 100, defaults['candidate_pool'])
        self.attempts = self._spin(
            100, 1000000, defaults['attempts_per_sequence'])

        form.addRow('已有序列输入文件：', inputWidget)
        form.addRow('新序列长度：', self.length)
        form.addRow('新生成数量：', self.count)
        form.addRow('全局 GC 下限：', self.gcMin)
        form.addRow('全局 GC 上限：', self.gcMax)
        form.addRow('最大连续相同碱基：', self.maxRun)
        form.addRow('最大自身互补长度：', selfWidget)
        form.addRow('发卡结构：', hairpinWidget)
        form.addRow('最大同向相同片段：', self.sameSubstring)
        form.addRow('最大链间互补片段：', self.crossComplement)
        form.addRow('骨架链：', self.scaffold)
        form.addRow('与骨架最大同向相同片段：', self.scaffoldSame)
        form.addRow('与骨架最大链间互补片段：', self.scaffoldCross)

        advanced = QLabel('高级规则（左侧勾选后参与筛选）', content)
        advanced.setStyleSheet('font-weight: 600; margin-top: 10px;')
        form.addRow(advanced)

        self.useLocalGc, localWidget = self._optionalWidget(
            defaults['use_local_gc'])
        self.localWindow = self._spin(4, 50, defaults['local_gc_window'])
        self.localMin = self._double(
            0, 100, defaults['local_gc_min'] * 100, 1, ' %')
        self.localMax = self._double(
            0, 100, defaults['local_gc_max'] * 100, 1, ' %')
        self._finishOptional(localWidget, self.useLocalGc, (
            QLabel('窗口', content), self.localWindow,
            QLabel('范围', content), self.localMin,
            QLabel('—', content), self.localMax))

        self.useEntropy, entropyWidget = self._optionalWidget(
            defaults['use_entropy'])
        self.entropy = self._double(
            0, 2, defaults['min_entropy'], 2, ' bits')
        self._finishOptional(
            entropyWidget, self.useEntropy, (self.entropy,))

        self.useHamming, hammingWidget = self._optionalWidget(
            defaults['use_hamming'])
        self.hamming = self._double(
            0, 100, defaults['min_hamming_fraction'] * 100, 1, ' %')
        self._finishOptional(hammingWidget, self.useHamming, (self.hamming,))

        self.useMotifs, motifWidget = self._optionalWidget(
            defaults['use_forbidden_motifs'])
        self.motifs = QLineEdit(
            ', '.join(defaults['forbidden_motifs']), content)
        self.motifs.setPlaceholderText('例如：GAATTC, GGTCTC')
        self._finishOptional(motifWidget, self.useMotifs, (self.motifs,))

        form.addRow('局部 GC：', localWidget)
        form.addRow('最低序列熵：', entropyWidget)
        form.addRow('最小汉明距离：', hammingWidget)
        form.addRow('禁用 motif：', motifWidget)

        runtime = QLabel('运行设置（防止卡死并改善候选选择）', content)
        runtime.setStyleSheet('font-weight: 600; margin-top: 10px;')
        form.addRow(runtime)
        form.addRow('每轮候选池大小：', self.pool)
        form.addRow('每条序列最大尝试次数：', self.attempts)

        tooltips = (
            ((self.inputPath, inputWidget, browse),
             '可选。TXT中每行一条已有序列。新序列必须满足与这些输入序列'
             '之间的同向相同片段和链间互补片段限制；输入序列也会保留在'
             '最终分析中。'),
            ((template,),
             '保存一个TXT模板。模板说明每行填写一条只含A、C、G、T的序列。'),
            ((self.length,),
             '每条新生成序列的碱基数。只影响新序列，输入TXT中的已有序列'
             '可以具有不同长度。'),
            ((self.count,), '本次需要新增的正交序列数量，不包含TXT输入序列。'),
            ((self.gcMin,), '新序列整条链中G和C所占比例的允许下限。'),
            ((self.gcMax,), '新序列整条链中G和C所占比例的允许上限。'),
            ((self.maxRun,),
             '允许连续出现的相同碱基最大数量。设置为3时，AAA允许，AAAA及'
             '更长不允许。'),
            ((self.sameSubstring,),
             '任意两条序列最多允许多少nt连续完全相同。输入5表示允许到5 nt，'
             '出现6 nt或更长即淘汰新序列。逐条和两两分析都会报告实际值。'),
            ((self.crossComplement,),
             '任意两条序列最多允许多少nt连续反向互补，用于限制非特异性杂交。'
             '输入5表示允许到5 nt，出现6 nt或更长即淘汰新序列。'),
            ((self.scaffold,),
             '可选一条cadnano内置骨架链。新延伸链必须同时避开这条骨架、TXT'
             '已有延伸链和本次已生成的其他延伸链。骨架链会在报告中单独标记。'),
            ((self.scaffoldSame,),
             '新延伸链与所选长骨架之间允许的最长同向相同片段。长骨架组合空间'
             '很大，推荐允许到7 nt；出现8 nt或更长时淘汰候选。'),
            ((self.scaffoldCross,),
             '新延伸链与所选骨架之间允许的最长反向互补片段，这是防止延伸链'
             '错误结合骨架的主要规则。推荐允许到7 nt；出现8 nt或更长时淘汰。'),
            ((self.useLocalGc, localWidget, self.localWindow,
              self.localMin, self.localMax),
             '可选规则。沿序列滑动指定长度的窗口，每个窗口的GC比例都必须在'
             '所设范围内，避免局部过度富含AT或GC。'),
            ((self.useEntropy, entropyWidget, self.entropy),
             '可选规则。序列熵反映A/T/C/G分布复杂度，理论最大值为2 bits。'
             '推荐最低1.70；越高越不容易出现低复杂度组成。'),
            ((self.useSelfComplement, selfWidget, self.selfComplement),
             '可选规则。限制序列与自身反向互补序列的最长连续匹配，降低自二聚体'
             '风险。一般推荐最大5 nt；更严格可设4 nt。'),
            ((self.useHairpin, hairpinWidget, self.hairpin),
             '检测由完全互补碱基形成的发卡茎，不再单独限制loop长度。一般推荐'
             '茎最大4 bp；更严格可设为3 bp。该长度模型不计算自由能。'),
            ((self.useHamming, hammingWidget, self.hamming),
             '可选规则。仅对等长序列逐位比较，要求至少有指定比例的位置不同。'
             '一般可从25%开始。'),
            ((self.useMotifs, motifWidget, self.motifs),
             '可选规则。禁止新序列包含所列片段或其反向互补片段；可用于排除'
             '酶切位点等实验相关序列，多个片段用逗号分隔。'),
            ((self.pool,),
             '每加入一条序列前先收集多少条合格候选，再选择与现有序列更分散的'
             '候选。数值越大通常结果越好，但运行更慢。'),
            ((self.attempts,),
             '为每条新序列允许评价的候选上限。达到上限仍找不到时停止，防止'
             '约束过严导致程序无限运行。'),
        )
        for widgets, tooltip in tooltips:
            for widget in widgets:
                widget.setToolTip(tooltip)

        scroll.setWidget(content)
        outer.addWidget(scroll, 1)
        note = QLabel(
            '每次运行都使用新的系统随机数，因此相同参数也会得到不同结果。'
            '高级规则即使没有勾选，相关数值仍会写入报告供查看，但不会淘汰序列。',
            self)
        note.setWordWrap(True)
        outer.addWidget(note)
        primer3Layout = QHBoxLayout()
        self.primer3Button = QPushButton('Primer3 热力学分析…', self)
        self.primer3Button.setIcon(QIcon('icons:primer3-analysis.svg'))
        self.primer3Button.setToolTip(
            '分析输入序列和新生成序列的发卡、自身二聚体及链间二聚体；'
            '不会读取骨架链或当前cadnano设计。')
        self.primer3Button.clicked.connect(self._openPrimer3)
        primer3Layout.addWidget(self.primer3Button)
        primer3Layout.addStretch(1)
        outer.addLayout(primer3Layout)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok |
            QDialogButtonBox.StandardButton.Cancel, parent=self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        outer.addWidget(buttons)

    def _openPrimer3(self):
        dialog = Primer3AnalysisDialog(
            self._primer3Entries, self._suggestedDirectory, self)
        dialog.exec()

    def _spin(self, minimum, maximum, value):
        widget = QSpinBox(self)
        widget.setRange(minimum, maximum)
        widget.setValue(int(value))
        return widget

    def _double(self, minimum, maximum, value, decimals, suffix=''):
        widget = QDoubleSpinBox(self)
        widget.setRange(minimum, maximum)
        widget.setDecimals(decimals)
        widget.setValue(float(value))
        widget.setSuffix(suffix)
        return widget

    def _optionalWidget(self, checked):
        checkbox = QCheckBox(self)
        checkbox.setChecked(bool(checked))
        widget = QWidget(self)
        layout = QHBoxLayout(widget)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(checkbox)
        widget._optionLayout = layout
        return checkbox, widget

    @staticmethod
    def _finishOptional(widget, checkbox, controls):
        for control in controls:
            widget._optionLayout.addWidget(control)
            control.setEnabled(checkbox.isChecked())
            checkbox.toggled.connect(control.setEnabled)
        widget._optionLayout.addStretch(1)

    def _browseInput(self):
        filename = QFileDialog.getOpenFileName(
            self, '选择已有正交序列 TXT', '', 'Text files (*.txt);;All files (*)')
        if isinstance(filename, (tuple, list)):
            filename = filename[0]
        if filename:
            self.inputPath.setText(filename)

    def _setScaffoldControlsEnabled(self):
        enabled = bool(self.scaffold.currentData())
        self.scaffoldSame.setEnabled(enabled)
        self.scaffoldCross.setEnabled(enabled)

    def _exportTemplate(self):
        filename = QFileDialog.getSaveFileName(
            self, '导出正交序列输入模板',
            os.path.expanduser('~/Desktop/orthogonal_sequence_input.txt'),
            'Text files (*.txt)')
        if isinstance(filename, (tuple, list)):
            filename = filename[0]
        if not filename:
            return
        if not filename.lower().endswith('.txt'):
            filename += '.txt'
        try:
            write_input_template(filename)
        except (IOError, OSError) as error:
            QMessageBox.critical(self, '输入模板', '无法写入模板：%s' % error)
            return
        QMessageBox.information(self, '输入模板', '模板已保存到：\n%s' % filename)

    def settings(self):
        return {
            'input_file': self.inputPath.text().strip(),
            'length': self.length.value(),
            'count': self.count.value(),
            'gc_min': self.gcMin.value() / 100.0,
            'gc_max': self.gcMax.value() / 100.0,
            'max_homopolymer': self.maxRun.value(),
            'max_same_substring': self.sameSubstring.value(),
            'max_cross_complement': self.crossComplement.value(),
            'scaffold_name': str(self.scaffold.currentData() or ''),
            'scaffold_max_same_substring': self.scaffoldSame.value(),
            'scaffold_max_cross_complement': self.scaffoldCross.value(),
            'use_local_gc': self.useLocalGc.isChecked(),
            'local_gc_window': self.localWindow.value(),
            'local_gc_min': self.localMin.value() / 100.0,
            'local_gc_max': self.localMax.value() / 100.0,
            'use_entropy': self.useEntropy.isChecked(),
            'min_entropy': self.entropy.value(),
            'use_self_complement': self.useSelfComplement.isChecked(),
            'max_self_complement': self.selfComplement.value(),
            'use_hairpin': self.useHairpin.isChecked(),
            'max_hairpin_stem': self.hairpin.value(),
            'use_hamming': self.useHamming.isChecked(),
            'min_hamming_fraction': self.hamming.value() / 100.0,
            'use_forbidden_motifs': self.useMotifs.isChecked(),
            'forbidden_motifs': self.motifs.text(),
            'candidate_pool': self.pool.value(),
            'attempts_per_sequence': self.attempts.value(),
        }
