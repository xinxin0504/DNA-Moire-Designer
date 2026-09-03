"""Final-only Twist and Bend indel-distribution report."""

import os

import cadnano2.util as util

from .curvedreport import (_draw_pair_curvature_chart,
                           _font, _pair_chart_height, _text)

util.qtWrapImport('QtCore', globals(), ['QPointF', 'QRectF', 'Qt'])
util.qtWrapImport('QtGui', globals(), [
    'QBrush', 'QColor', 'QImage', 'QPainter', 'QPen', 'QPixmap'])
util.qtWrapImport('QtWidgets', globals(), [
    'QDialog', 'QDialogButtonBox', 'QLabel', 'QScrollArea', 'QVBoxLayout'])


REPORT_WIDTH = 1400
REPORT_SCALE = 2.0
MARGIN = 38


def _single_chart_height(row_count):
    return 105 + max(1, int(row_count))*42


def _draw_single_helix_chart(painter, rows, top, foreground, muted, border):
    _text(painter, MARGIN, top+29,
          'Add/Remove Twist 与 Bend：单 helix 等分分布',
          foreground, 15, True)
    _text(painter, MARGIN, top+51,
          '圆点为最终新增 indel 坐标；CV 仅为末级均匀性诊断，'
          '不优先于曲率方向和理论区间。', muted, 9)
    if not rows:
        _text(painter, MARGIN, top+87, '没有新增 indel。', muted, 10)
        return top + _single_chart_height(0)
    minimum = min(min(row['positions']) for row in rows if row['positions'])
    maximum = max(max(row['positions']) for row in rows if row['positions'])
    span = max(1, maximum-minimum)
    plot_left = MARGIN+185
    plot_width = REPORT_WIDTH-plot_left-MARGIN-180
    chart_top = top+78
    row_height = 42
    colors = {
        'add_twist': QColor('#2A9D8F'),
        'remove_twist': QColor('#457B9D'),
        'bend': QColor('#8E5AB5')}
    for row_index, row in enumerate(rows):
        y = chart_top+row_index*row_height+row_height/2.0
        color = colors.get(row['kind'], QColor('#555555'))
        _text(painter, MARGIN, y+4,
              'P%d %s · H%d' % (
                  int(row['plan']), row['kind'], int(row['helix'])),
              foreground, 9, True)
        painter.setPen(QPen(border, 1.0))
        painter.drawLine(QPointF(plot_left, y),
                         QPointF(plot_left+plot_width, y))
        for position in row['positions']:
            x = plot_left+(float(position)-minimum)/span*plot_width
            painter.setPen(QPen(QColor('#FFFFFF'), .8))
            painter.setBrush(QBrush(color))
            painter.drawEllipse(QRectF(x-4.5, y-4.5, 9, 9))
        _text(painter, plot_left+plot_width+15, y+4,
              'n=%d · mean %.1f bp · CV %.1f%%' % (
                  int(row['count']), float(row['mean_spacing_bp']),
                  100.0*float(row['spacing_cv'])), muted, 8)
    return top + _single_chart_height(len(rows))


def create_twistbend_report_image(plans, single_rows, output_path):
    bend_plans = [plan for plan in plans if plan.get('kind') == 'bend']
    summary_height = 118
    single_height = _single_chart_height(len(single_rows))
    bend_height = sum(_pair_chart_height(len(plan.get(
        'pair_curvature_rows', ())))+24 for plan in bend_plans)
    image_height = summary_height+single_height+bend_height+45
    image = QImage(
        int(REPORT_WIDTH*REPORT_SCALE),
        int(image_height*REPORT_SCALE),
        QImage.Format.Format_ARGB32)
    image.fill(QColor(250, 251, 253))
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    painter.scale(REPORT_SCALE, REPORT_SCALE)
    foreground = QColor(36, 39, 44)
    muted = QColor(100, 108, 118)
    border = QColor(196, 201, 209)
    _text(painter, MARGIN, 45, 'Twist and Bend 最终设计报告',
          foreground, 18, True)
    _text(painter, MARGIN, 76,
          '任务 %d；Add Twist %d；Remove Twist %d；Add Bending %d' % (
              len(plans),
              sum(plan.get('kind') == 'add_twist' for plan in plans),
              sum(plan.get('kind') == 'remove_twist' for plan in plans),
              len(bend_plans)), foreground, 11, True)
    _text(painter, MARGIN, 101,
          '所有图均来自用户最终点击 OK 时的已接受方案；实时预览不导出报告。',
          muted, 9)
    top = _draw_single_helix_chart(
        painter, single_rows, summary_height, foreground, muted, border)
    for plan_index, plan in enumerate(bend_plans, 1):
        _text(painter, MARGIN, top+22,
              'Add Bending pair analysis %d' % plan_index,
              foreground, 12, True)
        top = _draw_pair_curvature_chart(
            painter, {
                'pair_curvature_rows': plan.get(
                    'pair_curvature_rows', ()),
                'pair_curvature_summary': plan.get(
                    'pair_curvature_summary', {})},
            top+20, foreground, muted, border)
    painter.end()
    directory = os.path.dirname(os.path.abspath(output_path))
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    if not image.save(output_path, 'PNG'):
        raise IOError('无法保存 Twist and Bend 报告：%s' % output_path)
    return image


class TwistBendReportDialog(QDialog):
    def __init__(self, image, parent=None):
        super(TwistBendReportDialog, self).__init__(parent)
        self.setWindowTitle('Twist and Bend 报告')
        self.resize(1180, 820)
        layout = QVBoxLayout(self)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(False)
        label = QLabel(scroll)
        pixmap = QPixmap.fromImage(image)
        pixmap.setDevicePixelRatio(REPORT_SCALE)
        label.setPixmap(pixmap)
        label.resize(int(image.width()/REPORT_SCALE),
                     int(image.height()/REPORT_SCALE))
        scroll.setWidget(label)
        layout.addWidget(scroll, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok, self)
        buttons.accepted.connect(self.accept)
        layout.addWidget(buttons)
