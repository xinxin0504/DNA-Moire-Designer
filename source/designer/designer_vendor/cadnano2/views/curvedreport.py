"""Scrollable image report shown after a Curved Design completes."""

import os

import cadnano2.util as util

util.qtWrapImport('QtCore', globals(), [
    'QPointF', 'QRectF', 'Qt'])
util.qtWrapImport('QtGui', globals(), [
    'QBrush', 'QColor', 'QFont', 'QImage', 'QPainter', 'QPen', 'QPixmap'])
util.qtWrapImport('QtWidgets', globals(), [
    'QDialog', 'QDialogButtonBox', 'QLabel', 'QScrollArea', 'QVBoxLayout',
    'QWidget'])


REPORT_WIDTH = 1400
REPORT_SCALE = 2.0
MARGIN = 38


def _font(point_size, bold=False):
    font = QFont()
    font.setPointSize(point_size)
    font.setBold(bool(bold))
    return font


def _text(painter, x, baseline, value, color, size=11, bold=False):
    painter.setPen(QPen(color))
    painter.setFont(_font(size, bold))
    painter.drawText(QPointF(float(x), float(baseline)), str(value))


def _centered_text(painter, rectangle, value, color, size=10,
                   bold=False):
    painter.setPen(QPen(color))
    painter.setFont(_font(size, bold))
    painter.drawText(
        QRectF(rectangle),
        int(Qt.AlignmentFlag.AlignCenter), str(value))


def _domain_chart_height(helix_count):
    return 112 + max(1, int(helix_count)) * 34


def _staple_chart_height(helix_count):
    return 122 + max(1, int(helix_count)) * 58


def _pair_chart_height(pair_count):
    return 126 + max(1, int(pair_count)) * 31


def _single_helix_chart_height(row_count):
    return 112 + max(1, int(row_count)) * 29


def _draw_summary(painter, lines, top, foreground, unused_muted,
                  report_mode="curved"):
    title = ("Frame Design 设计报告" if report_mode == "frame" else
             "Curved Design 设计报告")
    _text(painter, MARGIN, top + 30, title,
          foreground, 18, True)
    y = top + 67
    for line in lines:
        _text(painter, MARGIN, y, line, foreground, 12, True)
        y += 29
    return y + 10


def _draw_frame_straight_twist(painter, data, top, foreground, muted,
                               border):
    audit = dict(data.get(
        "frame_straight_common_mode_remove_twist", {}) or {})
    if not audit:
        return top
    painter.setPen(QPen(border, 1.0))
    painter.setBrush(QBrush(QColor(244, 247, 250)))
    painter.drawRoundedRect(
        QRectF(MARGIN, top, REPORT_WIDTH-2*MARGIN, 122), 7, 7)
    _text(painter, MARGIN+18, top+27,
          "非弯曲直边：共模 Remove Twist", foreground, 14, True)
    baseline = dict(audit.get("baseline_prediction", {}) or {})
    final = dict(audit.get("final_prediction", {}) or {})
    _text(
        painter, MARGIN+18, top+52,
        ("直边 %d bp；%s 平均 %.3f indel/helix（范围 %d–%d；总计 %d）；"
         "twist %.5f → %.5f °/base") % (
             int(audit.get("straight_native_bp", 0)),
             {"insertion": "insertion", "deletion": "deletion"}.get(
                 audit.get("correction"), "无需新增"),
             float(audit.get("common_indels_per_helix", 0.0)),
             int((audit.get("per_helix_indel_count_range") or [0, 0])[0]),
             int((audit.get("per_helix_indel_count_range") or [0, 0])[-1]),
             int(audit.get("total_new_indels", 0)),
             float(baseline.get("twist_per_base_deg", 0.0)),
             float(final.get("twist_per_base_deg", 0.0))),
        muted, 10)
    _text(
        painter, MARGIN+18, top+76,
        ("物理直边平均配额 %s；各边整数配额范围 %s；"
         "分数余数按截面一阶矩最小化并在等长边间轮换") % (
             "/".join("%.3f" % float(value) for value in audit.get(
                 "edge_common_quota",
                 audit.get("interval_common_quota", ()))),
             "/".join("%d–%d" % tuple(map(int, values))
                      for values in audit.get(
                          "edge_helix_quota_range", ()))),
        muted, 10)
    _text(
        painter, MARGIN+18, top+100,
        ("正值保护：%s；跨 base 0 的首尾坐标段按同一物理直边处理") %
        ("已阻止过校正" if audit.get("handedness_guard_applied") else
         "未触发"),
        muted, 10)
    return top+122


def _legend_box(painter, x, y, color, label, foreground, alpha=95,
                diamond=False):
    fill = QColor(color)
    fill.setAlpha(alpha)
    painter.setPen(QPen(color, 1.2))
    painter.setBrush(QBrush(fill))
    if diamond:
        painter.save()
        painter.translate(float(x + 7), float(y + 7))
        painter.rotate(45.0)
        painter.drawRect(QRectF(-5, -5, 10, 10))
        painter.restore()
    else:
        painter.drawRoundedRect(QRectF(x, y, 18, 13), 3, 3)
    _text(painter, x + 26, y + 12, label, foreground, 9)


def _draw_domain_chart(painter, data, top, foreground, muted, border):
    rows = list(data.get("domain_indels", []))
    domain_size = int(rows[0].get("domain_size_bp", 0)) if rows else (
        7 if str(data.get("lattice", "honeycomb")).lower() ==
        "honeycomb" else 8)
    colors = {
        -3: QColor("#2166AC"),
        -2: QColor("#67A9CF"),
        -1: QColor("#D1E5F0"),
        0: QColor("#111111"),
        1: QColor("#FDDBC7"),
        2: QColor("#EF8A62"),
        3: QColor("#B2182B")}
    _text(painter, MARGIN, top + 29,
          "每条 helix 的完整 domain 实际增删", foreground, 15, True)
    scope = ("仅显示顶点弯曲窗口；直边区域不计" if
             data.get("report_mode") == "frame" else "统计完整闭环")
    _text(painter, MARGIN, top + 50,
          "每根短竖线代表一个完整 %d-bp domain；%s；不完整尾端不计" %
          (domain_size, scope), muted, 9)
    legend_y = top + 61
    legend_x = MARGIN
    for value in range(-3, 4):
        label = "%+d" % value if value else "0"
        painter.setPen(QPen(colors[value], 2.2))
        painter.drawLine(QPointF(legend_x + 7, legend_y),
                         QPointF(legend_x + 7, legend_y + 16))
        _text(painter, legend_x + 17, legend_y + 13,
              label, foreground, 9)
        legend_x += 73

    plot_left = MARGIN + 67
    plot_width = REPORT_WIDTH - plot_left - MARGIN
    chart_top = top + 101
    row_height = 34
    maximum_domains = max(
        [len(row.get("values", [])) for row in rows] or [0])
    if not rows or maximum_domains <= 0:
        _text(painter, MARGIN, chart_top + 23,
              "没有可报告的完整 domain。", muted, 10)
        return top + _domain_chart_height(0)

    for row_index, row in enumerate(rows):
        center_y = chart_top + row_index * row_height + row_height / 2.0
        helix = int(row.get("helix", row_index))
        _text(painter, MARGIN, center_y + 4,
              "H%d" % helix, foreground, 9, True)
        painter.setPen(QPen(QColor(230, 233, 238), 0.8))
        painter.drawLine(QPointF(plot_left, center_y),
                         QPointF(plot_left + plot_width, center_y))
        values = list(row.get("values", []))
        for domain, value in enumerate(values):
            x = (plot_left +
                 (float(domain) + 0.5) * plot_width / maximum_domains)
            color = colors.get(int(value), QColor("#7A0177"))
            painter.setPen(QPen(color, 2.2))
            painter.drawLine(QPointF(x, center_y - 8),
                             QPointF(x, center_y + 8))
    return top + _domain_chart_height(len(rows))


def _draw_pair_curvature_chart(painter, data, top, foreground, muted,
                               border):
    rows = list(data.get("pair_curvature_rows", []))
    summary = dict(data.get("pair_curvature_summary", {}) or {})
    _text(painter, MARGIN, top + 29,
          "相邻 helix 对的局部曲率分布", foreground, 15, True)
    scope_note = ("；仅统计顶点弯曲窗口" if
                  data.get("report_mode") == "frame" else "；统计完整闭环")
    _text(
        painter, MARGIN, top + 51,
        ("曲率 helix 对：%d；crossover 区间：%d；反向：%d；"
         "超出理论区间：%d；%d-bp 轴向曲率 CV：%.2f%%%s") % (
            int(summary.get("physical_adjacent_curvature_pairs", 0)),
            int(summary.get("crossover_intervals", 0)),
            int(summary.get("reverse_curvature_intervals", 0)),
            int(summary.get("outside_floor_ceiling_intervals", 0)),
            int(summary.get(
                "normalized_axial_curvature_bin_width_bp", 42)),
            100.0 * float(summary.get(
                "normalized_axial_curvature_42bp_bin_cv", 0.0)),
            scope_note),
        muted, 9)
    legend_y = top + 65
    valid = QColor("#2A9D8F")
    reverse = QColor("#D73027")
    below = QColor("#F4A261")
    above = QColor("#7B2CBF")
    _legend_box(painter, MARGIN, legend_y, valid,
                "理论 floor/ceiling 内", foreground, 180)
    _legend_box(painter, MARGIN + 205, legend_y, reverse,
                "反向弯曲", foreground, 180)
    _legend_box(painter, MARGIN + 340, legend_y, below,
                "低于理论区间", foreground, 180)
    _legend_box(painter, MARGIN + 500, legend_y, above,
                "高于理论区间", foreground, 180)

    chart_top = top + 105
    label_width = 112
    plot_left = MARGIN + label_width
    plot_width = REPORT_WIDTH - plot_left - MARGIN
    row_height = 31
    if not rows:
        _text(painter, MARGIN, chart_top + 23,
              "没有检测到具有设计差异应变的物理相邻 helix 对。",
              muted, 10)
        return top + _pair_chart_height(0)
    for row_index, row in enumerate(rows):
        y = chart_top + row_index * row_height
        pair = tuple(row.get("pair", (0, 0)))
        floor_value = int(row.get("ideal_floor", 0))
        ceiling = int(row.get("ideal_ceiling", floor_value))
        _text(painter, MARGIN, y + 18,
              "%s  [%d–%d]" % (
                  row.get("label", "H%d–H%d" % pair),
                  floor_value, ceiling), foreground, 9, True)
        values = list(row.get("differences", []))
        cell_width = plot_width / float(max(1, len(values)))
        for index, difference in enumerate(values):
            color = (reverse if difference < 0 else
                     below if difference < floor_value else
                     above if difference > ceiling else valid)
            x = plot_left + index * cell_width
            fill = QColor(color)
            fill.setAlpha(190)
            rectangle = QRectF(x, y, max(1.0, cell_width), 23)
            painter.setPen(QPen(QColor(250, 251, 253), 0.7))
            painter.setBrush(QBrush(fill))
            painter.drawRoundedRect(rectangle, 2, 2)
            if cell_width >= 14:
                _centered_text(
                    painter, rectangle, difference, QColor("#FFFFFF"),
                    7, True)
    return top + _pair_chart_height(len(rows))


def _draw_single_helix_chart(painter, data, top, foreground, muted,
                             border):
    rows = list(data.get("single_helix_distribution", []))
    frame_mode = data.get("report_mode") == "frame"
    _text(painter, MARGIN, top + 29,
          "单 helix indel 轴向均匀性", foreground, 15, True)
    scope = ("仅统计每个顶点的弯曲窗口；直边区域排除" if frame_mode
             else "统计完整闭环")
    mean_cv = (sum(float(row.get("spacing_cv", 0.0)) for row in rows) /
               float(len(rows))) if rows else 0.0
    _text(painter, MARGIN, top + 51,
          "%s；平均间距 CV：%.2f%%" % (scope, 100.0*mean_cv),
          muted, 9)
    chart_top = top + 81
    row_height = 29
    label_width = 92
    stats_width = 245
    plot_left = MARGIN + label_width
    plot_width = REPORT_WIDTH - plot_left - MARGIN - stats_width
    if not rows:
        _text(painter, MARGIN, chart_top + 23,
              "没有可报告的 indel 位点。", muted, 10)
        return top + _single_helix_chart_height(0)
    marker = QColor("#2A9D8F")
    for row_index, row in enumerate(rows):
        center_y = chart_top + row_index*row_height + row_height/2.0
        vertex = int(row.get("vertex", -1))
        label = ("H%d · V%d" % (int(row["helix"]), vertex+1)
                 if vertex >= 0 else "H%d" % int(row["helix"]))
        _text(painter, MARGIN, center_y+4, label, foreground, 9, True)
        painter.setPen(QPen(border, 1.0))
        painter.drawLine(QPointF(plot_left, center_y),
                         QPointF(plot_left+plot_width, center_y))
        start = float(row.get("window_start", 0.0))
        end = float(row.get("window_end", data.get("base_count", 1)))
        span = max(1.0, end-start)
        painter.setPen(QPen(marker, 1.0))
        painter.setBrush(QBrush(marker))
        for position in row.get("positions", ()):
            x = plot_left + (float(position)-start)/span*plot_width
            x = min(plot_left+plot_width, max(plot_left, x))
            painter.drawEllipse(QRectF(x-3.0, center_y-3.0, 6.0, 6.0))
        _text(painter, plot_left+plot_width+16, center_y+4,
              "n=%d  mean=%.1f bp  CV=%.1f%%  min/max=%d/%d" % (
                  int(row.get("count", 0)),
                  float(row.get("mean_spacing_bp", 0.0)),
                  100.0*float(row.get("spacing_cv", 0.0)),
                  int(row.get("minimum_spacing_bp", 0)),
                  int(row.get("maximum_spacing_bp", 0))),
              muted, 8)
    return top + _single_helix_chart_height(len(rows))


def _draw_staple_chart(painter, data, top, foreground, muted, border):
    helices = list(data.get("helices", []))
    staples = list(data.get("staples", []))
    base_count = max(1, int(data.get("base_count", 1)))
    _text(painter, MARGIN, top + 29,
          "Staple 中心位置与实际长度", foreground, 15, True)
    legend_y = top + 46
    orange = QColor(242, 142, 43)
    green = QColor(89, 161, 79)
    light_blue = QColor(126, 190, 230)
    _legend_box(painter, MARGIN, legend_y, orange,
                "21–29 nt", foreground, 145)
    _legend_box(painter, MARGIN + 138, legend_y, green,
                "30–50 nt", foreground, 145)
    _legend_box(painter, MARGIN + 282, legend_y, light_blue,
                ">50 nt", foreground, 110)
    _text(painter, MARGIN + 420, legend_y + 12,
          "标记位置＝实际序列中心；标记内数字＝实际总长度（nt）",
          muted, 9)

    audit = dict(data.get("staple_nick_optimization", {}) or {})
    if audit:
        optimization_name = ("Frame" if
                             data.get("report_mode") == "frame" else
                             "Curved")
        dense_total = int(audit.get("deletion_dense_staples", 0))
        dense_good = int(audit.get("deletion_dense_in_40_60", 0))
        dense_outside = list(
            audit.get("deletion_dense_outside_40_60_nt", ()) or ())
        hard_invalid = int(audit.get("unbreakable_staples", 0))
        _text(
            painter, MARGIN, top + 84,
            ("%s 终态 nick 局部重切分：全部 staple %d–%d nt；"
             "deletion 密集区 40–60 nt：%d/%d；"
             "未能进入目标区间：%d（拓扑保留例外 %s）") % (
                 optimization_name,
                 int(audit.get("minimum_staple_nt", 0)),
                 int(audit.get("maximum_staple_nt", 0)),
                 dense_good, dense_total,
                 hard_invalid + len(dense_outside),
                 ",".join(map(str, dense_outside)) if dense_outside else
                 "无"),
            muted, 9)

    plot_left = MARGIN + 62
    plot_width = REPORT_WIDTH - plot_left - MARGIN
    chart_top = top + 115
    row_height = 58
    lattice = str(data.get("lattice", "honeycomb")).lower()
    tick_step = 32 if lattice == "square" else 21

    def x_for_base(base):
        return (plot_left +
                ((float(base) + 0.5) / float(base_count)) * plot_width)

    ticks = list(range(0, base_count, tick_step))
    if not ticks or ticks[-1] != base_count - 1:
        ticks.append(base_count - 1)
    chart_bottom = chart_top + max(1, len(helices)) * row_height
    for base in ticks:
        x = x_for_base(base)
        painter.setPen(QPen(QColor(215, 219, 225), 0.8))
        painter.drawLine(QPointF(x, chart_top - 9),
                         QPointF(x, chart_bottom - 7))
        _centered_text(
            painter, QRectF(x - 25, chart_top - 30, 50, 18),
            base, muted, 8)

    by_helix = {helix: [] for helix in helices}
    for item in staples:
        by_helix.setdefault(int(item["helix"]), []).append(item)
    for row_index, helix in enumerate(helices):
        row_top = chart_top + row_index * row_height
        axis_y = row_top + row_height / 2.0
        _text(painter, MARGIN, axis_y + 4,
              "H%d" % helix, foreground, 9, True)
        painter.setPen(QPen(border, 1.0))
        painter.drawLine(QPointF(plot_left, axis_y),
                         QPointF(plot_left + plot_width, axis_y))
        items = sorted(by_helix.get(helix, []),
                       key=lambda item: item["base"])
        previous_x = -10000.0
        previous_lane = 1
        for item in items:
            actual_x = x_for_base(item["base"])
            lane = 1 - previous_lane if actual_x - previous_x < 44 else 0
            previous_x = actual_x
            previous_lane = lane
            marker_y = row_top + (13 if lane == 0 else 45)
            marker_x = min(plot_left + plot_width - 19,
                           max(plot_left + 19, actual_x))
            length = int(item["length"])
            color = (orange if length < 30 else
                     green if length <= 50 else light_blue)
            fill = QColor(color)
            fill.setAlpha(150 if length <= 50 else 105)
            painter.setPen(QPen(border, 0.8))
            painter.drawLine(QPointF(actual_x, axis_y),
                             QPointF(marker_x, marker_y))
            rectangle = QRectF(marker_x - 18, marker_y - 11, 36, 22)
            painter.setPen(QPen(color, 1.3))
            painter.setBrush(QBrush(fill))
            painter.drawRoundedRect(rectangle, 11, 11)
            _centered_text(painter, rectangle, length, foreground, 8, True)
    return top + _staple_chart_height(len(helices))


def create_curved_report_image(summary_lines, report_data, output_path):
    """Render and save the complete textual and graphical PNG report."""
    summary_height = 115 + 29 * len(summary_lines)
    domain_height = _domain_chart_height(
        len(report_data.get("domain_indels", [])))
    pair_height = _pair_chart_height(
        len(report_data.get("pair_curvature_rows", [])))
    single_height = _single_helix_chart_height(
        len(report_data.get("single_helix_distribution", [])))
    staple_height = _staple_chart_height(
        len(report_data.get("helices", [])))
    twist_height = (140 if report_data.get(
        "frame_straight_common_mode_remove_twist") else 0)
    image_height = (summary_height + domain_height + single_height +
                    pair_height + twist_height +
                    staple_height + 75)
    image = QImage(
        int(REPORT_WIDTH * REPORT_SCALE),
        int(image_height * REPORT_SCALE),
        QImage.Format.Format_ARGB32)
    image.fill(QColor(250, 251, 253))
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
    painter.scale(REPORT_SCALE, REPORT_SCALE)
    foreground = QColor(36, 39, 44)
    muted = QColor(100, 108, 118)
    border = QColor(196, 201, 209)
    top = _draw_summary(
        painter, summary_lines, 18, foreground, muted,
        report_data.get("report_mode", "curved"))
    top = _draw_frame_straight_twist(
        painter, report_data, top+12, foreground, muted, border)
    top = _draw_domain_chart(
        painter, report_data, top + 18, foreground, muted, border)
    top = _draw_single_helix_chart(
        painter, report_data, top + 24, foreground, muted, border)
    top = _draw_pair_curvature_chart(
        painter, report_data, top + 24, foreground, muted, border)
    _draw_staple_chart(
        painter, report_data, top + 24, foreground, muted, border)
    painter.end()
    directory = os.path.dirname(os.path.abspath(output_path))
    if directory and not os.path.isdir(directory):
        os.makedirs(directory)
    if not image.save(output_path, "PNG"):
        raise IOError("无法保存Curved Design报告图片：%s" % output_path)
    return image


class CurvedReportDialog(QDialog):
    def __init__(self, image, parent=None):
        super(CurvedReportDialog, self).__init__(parent)
        self.setWindowTitle("Curved Design 报告")
        self.resize(1180, 820)
        layout = QVBoxLayout(self)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(False)
        image_label = QLabel(scroll)
        pixmap = QPixmap.fromImage(image)
        pixmap.setDevicePixelRatio(REPORT_SCALE)
        image_label.setPixmap(pixmap)
        image_label.resize(
            int(image.width() / REPORT_SCALE),
            int(image.height() / REPORT_SCALE))
        scroll.setWidget(image_label)
        layout.addWidget(scroll, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok, self)
        buttons.accepted.connect(self.accept)
        layout.addWidget(buttons)
