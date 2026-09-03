"""Central, offline localization for the Moiré Designer UI and exports.

The application historically contained literal Chinese/English UI strings.
This module deliberately translates at the presentation boundary so existing
design algorithms and serialized scientific identifiers stay unchanged.
"""

from __future__ import annotations

import csv
import html
import json
import os
from pathlib import Path
import re
import tempfile
import warnings
from zipfile import ZipFile


LANGUAGES = (
    ("en", "English", "English"),
)
LANGUAGE_CODES = {item[0] for item in LANGUAGES}
_language = os.environ.get("MOIRE_INTERFACE_LANGUAGE", "en")
if _language not in LANGUAGE_CODES:
    _language = "en"

_catalog_path = Path(__file__).with_name("translations.json")
try:
    _payload = json.loads(_catalog_path.read_text(encoding="utf-8"))
    _catalogs = dict(_payload.get("languages", {}))
except Exception:
    _catalogs = {}

# The project chooser must be usable before the full generated catalog is
# loaded.  These entries also make localization tests deterministic.
_BUILTIN = {
    "en": {
        "新建 Moiré 项目": "New Moiré Project",
        "新建 DNA Moiré 项目": "New DNA Moiré Project",
        "创建项目": "Create Project",
        "项目名": "Project name",
        "保存地址（项目文件夹的上一级）": "Save location (parent of project folder)",
        "选择文件夹…": "Choose Folder…",
        "项目文件：%s": "Project file: %s",
        "文件": "File", "编辑": "Edit", "视图": "View", "工具": "Tools",
        "新建项目…": "New Project…", "保存": "Save", "另存为…": "Save As…",
        "设计与预测": "Design and Prediction", "序列与导出": "Sequences and Export",
        "当前项目：未创建": "Current project: not created",
    },
}

# Human-reviewed overrides for long mixed-language strings and domain terms
# that generic translation engines commonly misinterpret.  The terms here
# follow cadnano/DNA-origami usage rather than everyday-language synonyms.
_ENGLISH_REVIEWED = {
    "已接受参数": "Accepted parameters",
    "撤销成功：已恢复到“%s”。": "Undo successful: restored “%s”.",
    "重做成功：已恢复到“%s”。": "Redo successful: restored “%s”.",
    "Seed insertion超过上限": "Seed insertion limit exceeded",
    "Seed insertion上限": "Seed insertion limit",
    "当前需要增加 %.1f bases/helix，超过上限 %.1f。请降低 Twist 或增大 spacing。":
        "The design requires %.1f added bases per helix, exceeding the %.1f limit. Reduce Twist or increase spacing.",
    "SST模板缺少helix：": "The SST template is missing helices: ",
    "没有结果": "No results",
    "序列已完成：%s；可重新导出 Final Export。":
        "Sequences are complete: %s; Final Export can be generated again.",
    "无法为Moiré scaffold找到稳定的单一nick位置。":
        "No stable single-nick position was found for the Moiré scaffold.",
    "两层均为Kagome；helix间距仍为2.8 nm，a使用cryo-EM 5.4 nm或干燥TEM 4.4 nm。":
        "Both layers are Kagome. The helix center spacing remains 2.8 nm; use a = 5.4 nm for cryo-EM or a = 4.4 nm for dried TEM.",
    "先 Export 模板，将序列填写到 Sequence 列，再 Import。行顺序按完整 base 数值优先、helix 数值其次排列。":
        "Export the template first, enter sequences in the Sequence column, and then import it. Rows are sorted first by the full base index and then by helix number.",
    "请先把序列复制到 Export Input Template 的 Sequence 列。Input 序列也可由下方独立的正交序列设计生成。":
        "First copy the sequences into the Sequence column of Export Input Template. Input sequences can also be generated with the independent Orthogonal Sequence Design tool below.",
    "已自动识别并添加颗粒 %s 到 Group %s – %s。继续点击可新增；按住 Option/Alt 可移动或旋转已有框，松开后恢复新增。":
        "Particle %s was detected automatically and added to Group %s – %s. Continue clicking to add particles. Hold Option/Alt to move or rotate existing boxes; release it to resume adding.",
    "颜色：<span style='color:%s'>■</span> %s<br>位置：%s → %s<br>长度：<b>%d nt</b>":
        "Color: <span style='color:%s'>■</span> %s<br>Position: %s → %s<br>Length: <b>%d nt</b>",
    "1. 上传并调整图像\n2. 右侧直接拖出选框\n3. 右键赋值为 Lane / Target / Staples\n4. 移动或拉伸选框后点击分析":
        "1. Upload and adjust the image\n2. Draw selection boxes directly on the right\n3. Right-click to assign Lane / Target / Staples\n4. Move or resize the boxes, then click Analyze",
    "固定两层SST已在后台生成并验证。\n已生成待审核Scaffold文件：%s\n在cadnano中编辑后请直接保存；接受时会自动读取该文件。\n%s":
        "The fixed two-layer SST was generated and validated in the background.\nScaffold file ready for review: %s\nAfter editing in cadnano, save the file directly. It will be read automatically when you accept it.\n%s",
    "普通状态：左键拖动选框；点击选框后拖动圆形手柄旋转。添加模式点击后自动识别并放置；按住 Option/Alt 时暂停新增，只移动或旋转已有框。":
        "Normal mode: drag a selection box with the left mouse button; click a box and drag its circular handle to rotate it. In Add mode, clicking detects and places an object automatically. Hold Option/Alt to pause adding and only move or rotate existing boxes.",
    "Seed目标Z1/Z3=%d/%d bp；最近合法routing=%d/%d bp。Seed短于SST时允许悬伸并只在重叠区生成capture。":
        "Seed target Z1/Z3 = %d/%d bp; nearest legal routing = %d/%d bp. If the Seed is shorter than the SST, overhangs are allowed and capture is generated only in the overlap region.",
    "Group %s 线长度统计（中心路径，mean ± STDEV.S，n=%d）：%s nm\n线型对象只统计长度，不参与颗粒产率、尺寸或面积密度。":
        "Group %s line-length statistics (center path, mean ± STDEV.S, n=%d): %s nm\nLine objects contribute only to length statistics, not to particle yield, size, or area density.",
    "SST 1st = 2nd（双向联动）：32n%+d；spacing/Seed Z2：32n%+d\nspacing：%d–%d bp，步长32 bp":
        "SST 1st = 2nd (bidirectionally linked): 32n%+d; spacing/Seed Z2: 32n%+d\nSpacing: %d–%d bp in 32-bp steps",
    "先自动读取当前已接受结构中的每条 scaffold 位置与长度。然后在右侧为每条 scaffold 分别添加 cadnano 保存的带序列 JSON。":
        "First read the position and length of every scaffold in the accepted structure automatically. Then add a sequence-containing JSON saved by cadnano for each scaffold on the right.",
    "双层设计与序列一致；1st=%d bp、spacing/Seed Z2=%d bp、2nd=%d bp；满足8 bp domain与32 bp重复相位。":
        "The two layers use identical designs and sequences: 1st = %d bp, spacing/Seed Z2 = %d bp, and 2nd = %d bp; the 8-bp domains and 32-bp repeat phase are satisfied.",
    "所有样品使用真实FFT和相位保持的Selected-spot孔径，并自动生成原始分辨率PNG、保持比例的SVG、JSON和TEM/FFT分列的CSV统计。":
        "All samples use the measured FFT and phase-preserving selected-spot apertures. The software automatically generates full-resolution PNG files, aspect-preserving SVG files, JSON data, and CSV statistics with separate TEM and FFT columns.",
    "当前 Bulk 模式要求每个文件以‘<数值> <单位>_’开头。\n例：5 nm_sample.tif、0.2 µm_test.tif\n\n请重命名后重新选择：":
        "In Bulk mode, every file name must begin with '<value> <unit>_'.\nExamples: 5 nm_sample.tif and 0.2 µm_test.tif\n\nRename and select these files again:",
    "自动识别每层 SST input 的位置、数量和长度。可先导出模板，将序列填入后再导入；模板操作显示在右侧识别结果下方。相同双层只填写一次并按对应位置复制。":
        "The position, count, and length of SST inputs in each layer are detected automatically. Export a template, enter the sequences, and import it; template controls appear below the detected results on the right. For identical layers, enter each sequence once and copy it to the corresponding position.",
    "1. 上传 TEM 图；可以‘<数值> <单位>_’开头命名（如 20 nm_xx）\n2. 选择单层或双层分析\n3. 点击分析，自动识别 scale bar、TEM 与 FFT\n4. 核对尺度与结果；无 OCR/文件名尺度时才手动输入":
        "1. Upload a TEM image; its name may begin with '<value> <unit>_' (for example, 20 nm_xx)\n2. Select single-layer or bilayer analysis\n3. Click Analyze to detect the scale bar and analyze the TEM image and FFT automatically\n4. Verify the scale and results; enter the scale manually only when OCR and the file name provide no scale",
    "<b>这些参数通常保持默认即可。</b><br>• 分类后端：内置模式无需外部软件；EMAN2 仅在安装后可选。<br>• 主页面直接调整信号方向、阈值、背景展平和降噪强度。<br>• 自动分组数量：0 为软件建议；指定数字会强制分成相应组数。<br>• 模板匹配门槛：越高越严格、误检更少；越低召回更多。<br>• 密度定义：材料面积使用真实掩膜；外包络会填充内部空洞。":
        "<b>These parameters can usually remain at their defaults.</b><br>• Classification backend: Built-in mode requires no external software; EMAN2 is available only when installed.<br>• Signal polarity, threshold, background flattening, and denoising strength are adjusted directly on the main page.<br>• Automatic group count: 0 uses the software recommendation; another value forces that number of groups.<br>• Template-matching threshold: higher values are stricter and reduce false positives; lower values increase recall.<br>• Density definition: material area uses the actual mask; outer-envelope area fills internal holes.",
    "无法找到所选骨架链：%s":
        "The selected scaffold could not be found: %s",
    "已读取 %d 条输入序列、采用 %d 条骨架链、生成 %d 条新序列，并保存到：\n%s\n\n“序列分析”会区分输入、骨架链和新生成序列；骨架链仅参与筛选，不会列入“两两分析”。":
        "Read %d input sequences, used %d scaffolds, generated %d new sequences, and saved them to:\n%s\n\nSequence Analysis distinguishes input sequences, scaffolds, and newly generated sequences. Scaffolds are used only for screening and are excluded from pairwise analysis.",
    "Capture协同性": "Capture cooperativity",
    "非-staple 比例": "Non-staple fraction",
    "Scaffold已接受": "Scaffold accepted",
    "Scaffold不合格": "Scaffold validation failed",
    "双层点阵对称性": "Bilayer lattice symmetry",
    "1. SST 双层点阵对称性": "1. SST bilayer lattice symmetry",
    "请至少选择4根Square网格helix。":
        "Select at least four helices on the Square lattice.",
    "%d×%d Square网格，%d根helix；%s。":
        "%d×%d Square lattice, %d helices; %s.",
    "当前选择：8×8 Square网格，%d根helix。%s":
        "Current selection: 8×8 Square lattice, %d helices. %s",
    # Orthogonal-sequence and Primer3 terminology.  These strings are also
    # used in exported workbooks, so prefer standard oligonucleotide terms.
    "发卡": "Hairpin",
    "发卡结构：": "Hairpin structure:",
    "启用发卡规则": "Enable hairpin rule",
    "汉明距离": "Hamming distance",
    "骨架链": "Scaffold",
    "骨架链：": "Scaffold:",
    "骨架链名称": "Scaffold name",
    "骨架链数量": "Number of scaffolds",
    "不使用骨架链": "Do not use a scaffold",
    "骨架链-%03d": "Scaffold-%03d",
    "禁用片段": "Forbidden sequence",
    "包含禁用片段": "Contains a forbidden sequence",
    "启用禁用片段规则": "Enable forbidden-sequence rule",
    "禁用 motif：": "Forbidden motifs:",
    "链间互补": "Interstrand complementarity",
    "同向相同片段": "Identical same-orientation segment",
    "最大同向相同片段":
        "Maximum identical same-orientation segment",
    "最大同向相同片段：":
        "Maximum identical same-orientation segment:",
    "同向相同片段（nt）":
        "Identical same-orientation segment (nt)",
    "同向相同片段过长":
        "Identical same-orientation segment is too long",
    "最大链间互补片段":
        "Maximum interstrand complementary segment",
    "最大链间互补片段：":
        "Maximum interstrand complementary segment:",
    "链间互补片段（nt）":
        "Interstrand complementary segment (nt)",
    "链间互补片段过长":
        "Interstrand complementary segment is too long",
    "自身互补": "Self-complementarity",
    "最大自身互补长度": "Maximum self-complementarity length",
    "最大自身互补长度：": "Maximum self-complementarity length:",
    "自身互补长度（nt）": "Self-complementarity length (nt)",
    "启用自身互补规则": "Enable self-complementarity rule",
    "新生成数量：": "Number to generate:",
    "输入序列数量": "Number of input sequences",
    "序列（5′→3′）": "Sequence (5′→3′)",
    "输入-%03d": "Input-%03d",
    "导入": "Import",
    "发现": "Found",
    "启用": "Enabled",
    "未发现": "Not found",
    "未启用": "Disabled",
    "全不选": "Select none",
    "需复核": "Review required",
    "二维配对图": "2D pairing diagram",
    "输入文件": "Input file",
    "输入模板": "Input template",
    "已有序列输入文件：": "Existing-sequence input file:",
    "请求生成数量": "Requested number of sequences",
    "实际生成数量": "Number of sequences generated",
    "是否完整生成": "Generation completed",
    "与骨架链：%s": "Against scaffold: %s",
    "输入/新生成链间：%s":
        "Between input/generated sequences: %s",
    "熔解温度链浓度":
        "Strand concentration for melting-temperature calculation",
    "熔解温度互补目标链浓度":
        "Complementary-target concentration for melting-temperature calculation",
    "熔解温度Na⁺浓度":
        "Na⁺ concentration for melting-temperature calculation",
    "熔解温度Mg²⁺浓度":
        "Mg²⁺ concentration for melting-temperature calculation",
    "正交序列设计": "Orthogonal Sequence Design",
    "序列分析": "Sequence Analysis",
    "两两分析": "Pairwise Analysis",
    "链间二聚体（所选序列的全部两两组合）":
        "Interstrand dimers (all pairwise combinations of selected sequences)",
    "链间二聚体分析至少需要两条序列。":
        "Interstrand-dimer analysis requires at least two sequences.",
    # Current TEM/domain-analysis UI added by the parallel analysis update.
    "撤回": "Undo",
    "重新分析": "Reanalyze",
    "完成当前编辑": "Finish Editing",
    "清空所有选框": "Clear All Selection Boxes",
    "统计尺寸依据": "Dimension Measurement Basis",
    "选框几何尺寸": "Selection-box Geometry",
    "框内有效颗粒": "Valid Particle Region Within the Box",
    # Particle-analysis table terminology.  These labels are also rebuilt
    # dynamically when the selected particle geometry changes.
    "Group归属": "Group ID",
    "颗粒归属": "Particle Classification",
    "长/外径 nm": "Length/Outer Diameter (nm)",
    "宽/内径 nm": "Width/Inner Diameter (nm)",
    "面积 nm²": "Area (nm²)",
    "圆度": "Circularity",
    "边缘": "Edge Contact",
    "主轴长度 / 次轴宽度": "Major-axis Length / Minor-axis Width",
    "外径 / 内径": "Outer Diameter / Inner Diameter",
    "外径 / 环宽": "Outer Diameter / Ring Width",
    "最长边 / 最短边": "Longest Side / Shortest Side",
    "平均边长 / 周长": "Mean Side Length / Perimeter",
    "等效外径 / 周长": "Equivalent Diameter / Perimeter",
    "最大 Feret 外径 / 最小 Feret 外径":
        "Maximum Feret Diameter / Minimum Feret Diameter",
    "中心线长度": "Centerline Length",
    "中心路径长度": "Centerline Path Length",
    "长 / 宽": "Length / Width",
    "主尺寸": "Primary Dimension",
    "次尺寸": "Secondary Dimension",
    "1. 上传并调整图像；在右侧拖出选框并右键赋值\n2. 移动或拉伸选框后点击分析":
        "1. Upload and adjust the image; draw boxes on the right and right-click to assign them\n2. Move or resize the boxes, then click Analyze",
    "右侧框选说明": "Selection-box Instructions",
    "2. 分析": "2. Analyze",
    "2. 开始分析": "2. Start analysis",
    "3. 自动识别并分析": "3. Automatically identify and analyze",
    "3. 批量自动识别并分析":
        "3. Run batch detection and analysis",
    "当前为双层分析": "Current analysis mode: Bilayer",
    "单层分析不计算Moiré周期。":
        "Single-layer analysis does not calculate the Moiré period.",
    "TEM内可用Moiré单元过少，无法可靠拟合周期。":
        "Too few Moiré unit cells are visible in the TEM image to fit the period reliably.",
    "Bulk 分析进行中…\n准备逐张分析并导出。":
        "Bulk analysis in progress…\nPreparing to analyze and export each image.",
    "Bulk 分析进行中…\n%d/%d：%s":
        "Bulk analysis in progress…\n%d/%d: %s",
    "已撤回上一次 domain / background 编辑。":
        "The last domain/background edit has been undone.",
    "增加：手动勾勒最终边界；删除：悬停识别连续区域后点击。":
        "Add: manually outline the final boundary. Delete: hover to identify a contiguous region, then click it.",
    "正在进行局部FFT、边界细化并更新全部分析结果…":
        "Running local FFT analysis, refining boundaries, and updating all results…",
    "重新分析中…\n正在更新点阵、面积、边界和 FFT 结果。":
        "Reanalyzing…\nUpdating lattice, area, boundary, and FFT results.",
    "编辑完成：面积、点阵、a、取向、boundary、FFT摘要及导出结果已同步更新。":
        "Editing complete: area, lattice, a, orientation, boundaries, FFT summary, and exported results have been updated consistently.",
    "右键赋值 Lane / Target / Staples · 拖动或拉伸调整 · 靠近边界自动吸附":
        "Right-click to assign Lane / Target / Staples · Drag or resize to adjust · Snap automatically near a boundary",
    "Original image（在这里拖框）":
        "Original image (draw selection boxes here)",
    "左键创建代表框；选中后直接拖动框体、角点和圆形旋转手柄":
        "Left-click to create a representative box. After selecting it, drag the box, corner handles, or circular rotation handle directly.",
    "已形成%d个手动区域（最新边界%d个采样点）。边界已冻结且不会拟合；可继续勾勒，点击完成后才重新分析。":
        "%d manual regions created (the latest boundary has %d sampled points). The boundary is frozen and will not be refitted; continue outlining regions or click Finish to reanalyze.",
    # Values that may arrive dynamically from cadnano palettes or analysis
    # records rather than as literal UI strings.
    "红色": "Red", "蓝色": "Blue", "绿色": "Green",
    "黄色": "Yellow", "橙色": "Orange", "紫色": "Purple",
    "粉色": "Pink", "黑色": "Black", "灰色": "Gray",
    "白色": "White", "青色": "Cyan",
    "目标颗粒": "Target particle", "副产物": "Byproduct",
    "聚集体": "Aggregate",
    "错误识别（排除统计）": "Misidentified (excluded from statistics)",
    "样品边界": "Sample boundary", "晶畴边界": "Domain boundary",
    "请至少选择一种分析内容。":
        "Select at least one analysis type.",
    "高级规则（左侧勾选后参与筛选）":
        "Advanced rules (enable a checkbox on the left to apply a rule during screening)",
    "完成 %d 项；%d 项因长度或计算错误未完成；已按 ΔG 最危险优先排序":
        "%d analyses completed; %d could not be completed because of sequence length or calculation errors; sorted by the most negative (highest-risk) ΔG first",
    "所选 Seed 截面与 SST 双层点阵预览":
        "Selected Seed cross-section and SST bilayer lattice preview",
    "SST input-only 视图；相同双层只显示一个导入组，最终仍分别导出两层。":
        "SST input-only view. For identical bilayers, only one import group is shown; the two layers are still exported separately.",
    "由Twist和当前点阵a计算。":
        "Calculated from twist and the current lattice constant a.",
    "FFT没有稳定识别到两套点阵。":
        "FFT did not reliably resolve two lattice sets.",
    "# 正交序列输入模板\n# 请在下方每行粘贴一条只包含A、C、G、T的已有序列。\n# 空行及以#开头的说明行会被忽略。\n# 不同输入序列可以具有不同长度。":
        "# Orthogonal-sequence input template\n# Paste one existing sequence containing only A, C, G, and T on each line below.\n# Blank lines and comment lines beginning with # are ignored.\n# Input sequences may have different lengths.",
}

# Human-reviewed text for the staged Design workflow.  These entries cover
# both widgets and exported/status prose; no machine translation is used at
# runtime, so switching to English remains deterministic and scientifically
# unambiguous.
_ENGLISH_REVIEWED.update({
    "三维 Seed/SST 预览": "3D Seed/SST Preview",
    "内嵌 cadnano 设计视野": "Embedded cadnano design vision",
    "序列位置与结构预览": "Sequence position and structure preview",
    "1.1 对称性与 Twist": "1.1 Symmetry and Twist",
    "1.2 Seed/SST 参数": "1.2 Seed/SST Parameters",
    "下一步：1.2 Seed/SST 参数": "Next: 1.2 Seed/SST Parameters",
    "接受 1.1 参数": "Accept 1.1 Parameters",
    "接受 1.2 参数": "Accept 1.2 Parameters",
    "✓ 1.1 已接受": "✓ 1.1 Accepted",
    "✓ 1.2 已接受": "✓ 1.2 Accepted",
    "可进入1.2 Seed/SST参数。": "Ready for 1.2 Seed/SST Parameters.",
    "接受 1.1 对称性与 Twist 参数":
        "Accept 1.1 Symmetry and Twist Parameters",
    "请先接受 1.1 对称性与 Twist":
        "Accept 1.1 Symmetry and Twist first",
    "接受 1.2 Seed/SST 参数": "Accept 1.2 Seed/SST Parameters",
    "请先接受 1.1 与 1.2 参数。":
        "Accept the 1.1 and 1.2 parameters first.",
    "1A-设计-top view": "1A-Design-Top View",
    "1B-设计-side view": "1B-Design-Side View",
    "下一步：1B-设计-side view": "Next: 1B-Design-Side View",
    "接受 1A 参数": "Accept 1A Parameters",
    "接受 1B 参数": "Accept 1B Parameters",
    "1.1 设计 · Top View": "1.1 Design · Top View",
    "1 Moiré 参数输入": "1 Moiré Parameter Input",
    "1.2 层长度与层间距":
        "1.2 Layer Length and Interlayer Spacing",
    "1.1 选择双层对称性": "1.1 Select Bilayer Symmetry",
    "双层对称性": "Bilayer Symmetry",
    "1.2 输入 Twist 或 Moiré period":
        "1.2 Enter Twist or Moiré Period",
    "1.3 输入 SST superlattice 参数":
        "1.3 Enter SST Superlattice Parameters",
    "Twist / Moiré period": "Twist / Moiré Period",
    "3. Seed 截面": "3. Seed Cross-section",
    "4. 接受当前对称性与 Twist":
        "4. Accept Current Symmetry and Twist",
    "下一步：层长度与层间距":
        "Next: Layer Length and Interlayer Spacing",
    "层长度与层间距": "Layer Length and Interlayer Spacing",
    "接受当前 Moiré 参数": "Accept Current Moiré Parameters",
    "✓ 当前 Moiré 参数已接受":
        "✓ Current Moiré Parameters Accepted",
    "Moiré 参数已接受。下一步：打开 Automated DNA Design 并生成三个设计文件。":
        "Moiré parameters accepted. Next: open Automated DNA Design "
        "and generate the three design files.",
    "预览通道": "Preview Channel",
    "同时显示 2D + 3D": "Show 2D + 3D",
    "仅显示 2D": "2D Only",
    "仅显示 3D": "3D Only",
    "2D Moiré 预览": "2D Moiré Preview",
    "Top view": "Top view",
    "Side view": "Side view",
    "Seed helix": "Seed helix",
    "1st SST helix": "1st SST helix",
    "2nd SST helix": "2nd SST helix",
    "SST superlattice 第二层 / Seed Z3":
        "2nd SST Superlattice / Seed Z3",
    "SST superlattice 第一层 / Seed Z1":
        "1st SST Superlattice / Seed Z1",
    "✓ 1A 已接受": "✓ 1A Accepted",
    "✓ 1B 已接受": "✓ 1B Accepted",
    "1st support · Z1": "1st support · Z1",
    "Seed Z2": "Seed Z2",
    "2nd support · Z3": "2nd support · Z3",
    "8×8 + 4×4 pore": "8×8 + 4×4 pore",
    "Optional · Expert mode": "Optional · Expert Mode",
    "收起 Optional Expert mode": "Collapse Optional Expert Mode",
    "在 cadnano 内检查（Optional）": "Inspect in cadnano (Optional)",
    "Optional：选择 JSON 在 cadnano 内检查":
        "Optional: Choose a JSON to Inspect in cadnano",
    "Optional：打开 cadnano": "Optional: Open cadnano",
    "直接打开 cadnano；请在 cadnano 中打开并修改当前项目 cadnano design 文件夹内的最终 JSON。点击接受时，软件会自动选择该文件夹中修改时间最新的合法 SST + Scaffold + Staple + Capture 文件。":
        "Open cadnano directly, then open and edit the final JSON in the current project's cadnano design folder. When you accept, the software automatically selects the most recently modified valid SST + Scaffold + Staple + Capture file in that folder.",
    "cadnano 已打开。请在其中选择并修改最终 JSON；保存后点击接受，软件会使用当前设计文件夹中修改时间最新的合法最终文件。":
        "cadnano is open. Select and edit the final JSON there; after saving, click Accept to use the most recently modified valid final file in the current design folder.",
    "正在处理，请稍候…": "Processing, please wait…",
    "任务进行中": "Task in Progress",
    "正在读取 Scaffold 位置和长度…":
        "Reading Scaffold positions and lengths…",
    "骨架链实际长度：%d nt\n3条骨架链总长度限制：22671 nt\n请减小 Seed 长度。":
        "Actual Scaffold length: %d nt\nTotal limit for three Scaffolds: 22,671 nt\nReduce the Seed length.",
    "实际 Z1 / Z2 / Z3：%s / %s / %s bp\nScaffold：%d条 · %s nt\n%s":
        "Actual Z1 / Z2 / Z3: %s / %s / %s bp\nScaffolds: %d · %s nt\n%s",
    "Normal staple：%s–%s nt": "Normal staples: %s–%s nt",
    "长度分布：": "Length distribution: ",
    "Staple-capture：%d–%d nt": "Staple-capture: %d–%d nt",
    "具有连续16-base区域：%.1f%%":
        "Staples with a continuous 16-base region: %.1f%%",
    "右侧可上下滚动：cadnano 左视野（Slice）、Seed 侧视 capture、cadnano 右视野（Path）；三个视野均可独立缩放和平移，Path 单击可定位具体 base。":
        "The right panel contains vertically scrollable cadnano Slice, Seed capture side-view, and cadnano Path views. Each view supports independent zoom and pan; click the Path view to locate a specific base.",
    "中间区：两个 Seed capture face。":
        "The center panel shows the two Seed capture faces.",
    "Seed capture face 仅在完整结构通道显示。":
        "Seed capture faces are shown only in the complete-structure channel.",
    "统计读取失败：%s": "Failed to read statistics: %s",
    "1.1  设计预测–2D": "1.1  Design Prediction–2D",
    "1.2  设计预测–3D": "1.2  Design Prediction–3D",
    "撤回上一个设计操作": "Undo the previous design operation",
    "恢复下一个设计操作": "Redo the next design operation",
    "选择双层点阵、Seed截面、Twist与成像环境。右侧二维图会随Twist、Moiré period及截面预设实时更新。":
        "Select the bilayer lattice, Seed cross-section, twist, and imaging environment. The 2D view updates in real time with twist, Moiré period, and the selected cross-section preset.",
    "点阵对称性": "Lattice symmetry",
    "Twist 与点阵测量环境": "Twist and lattice measurement environment",
    "一致": "Consistent",
    "不一致": "Inconsistent",
    "Seed 截面预设": "Seed cross-section preset",
    "恢复当前预设": "Restore Current Preset",
    "接受 1.1 参数": "Accept 1.1 Parameters",
    "下一步：1.2 设计预测–3D": "Next: 1.2 Design Prediction–3D",
    "二维点阵与 Moiré 预览": "2D Lattice and Moiré Preview",
    "SST superlattice 参数": "SST Superlattice Parameters",
    "Seed S(F) 分区": "Seed S(F) Segments",
    "Seed spacing（与SST联动）": "Seed spacing (linked to SST)",
    "接受 1.2 参数": "Accept 1.2 Parameters",
    "下一步：Scaffold / Capture": "Next: Scaffold / Capture",
    "暂无设计参数": "No accepted Moiré parameters",
    "生成 Seed + SST superlattice 设计图": "Generate Seed + SST Superlattice Design",
    "生成并导出全部 3 个设计文件": "Generate and Export All 3 Design Files",
    "重新选择点阵 / Seed 截面": "Reselect Lattice / Seed Cross-section",
    "展开验证说明  ▴": "Expand Validation Notes  ▴",
    "导入 Moiré 工程 (.moire.json)": "Import Moiré Project (.moire.json)",
    "在 cadnano 内专家编辑完成 Scaffold routing":
        "Edit Scaffold Routing in cadnano",
    "载入 cadnano 专家编辑后的 JSON":
        "Load Expert-Edited cadnano JSON",
    "接受 cadnano 当前保存的 Scaffold routing":
        "Accept Current Scaffold Routing Saved in cadnano",
    "生成 Staple / Capture 设计": "Generate Staple/Capture Design",
    "在 cadnano 内专家编辑完成结构": "Edit Structure in cadnano",
    "接受当前 Added Scaffold": "Accept Current Added Scaffold",
    "一次生成 SST、SST + Scaffold，以及最终的 SST + Scaffold + Staple + Capture JSON。":
        "Generate the SST, SST + Scaffold, and final SST + Scaffold + Staple + Capture JSON files in one run.",
    "专家模式": "Expert Mode",
    "收起专家模式": "Collapse Expert Mode",
    "接受当前设计图": "Accept Current Design",
    "下一步：序列导出": "Next: Sequence Export",
    "尚未接受设计图": "No accepted DNA design",
    "3.1  Add Scaffold 序列": "3.1  Add Scaffold Sequences",
    "下一步：SST superlattice Input": "Next: SST Superlattice Input",
    "3.2  Add SST superlattice Input 序列": "3.2  Add SST Superlattice Input Sequences",
    "自动识别 SST superlattice Input 位置和长度": "Detect SST Superlattice Input Positions and Lengths",
    "自动设计并 Add Input 序列": "Automatically Design and Add Input Sequences",
    "Import Input Sequences": "Import Input Sequences",
    "正交序列设计": "Orthogonal Sequence Design",
    "接受当前 Input": "Accept Current Input",
    "✓ 当前 Input 已接受": "✓ Current Input Accepted",
    "接受当前 Added SST superlattice Input": "Accept the Added SST Superlattice Input",
    "导出订购序列和单独的 Seed、SST superlattice Layer 1、Layer 2 三维结构。":
        "Export ordering sequences and separate 3D structures for the Seed, SST superlattice Layer 1, and SST superlattice Layer 2.",
    "已恢复：%s（已导出的磁盘文件不会被删除）": "Restored: %s (files already exported to disk are not deleted)",
    "设计状态": "design state",
    "已接受": "Accepted",
    "已生成": "Generated",
    "待审核": "Pending review",
    "单独 SST": "SST only",
    "当前预设：%s，%d根helix。%s": "Current preset: %s, %d helices. %s",
    "可进入1.2设计。": "Ready for step 1.2.",
    "至少需要4根helix。": "At least four helices are required.",
    "✓ 1.1 已接受": "✓ 1.1 Accepted",
    "✓ 1.2 已接受": "✓ 1.2 Accepted",
    "接受 1.1 设计参数": "Accept 1.1 Design Parameters",
    "%d bp（由 Layer spacing 强制联动）": "%d bp (locked to Layer spacing)",
    "两层长度与spacing独立，步长均为8 bp。\n不强制32 bp相位联动。": "The two layer lengths and spacing are independent, all in 8-bp increments.\nNo 32-bp phase linkage is enforced.",
    "SST superlattice两层长度相同：32n%+d；Layer spacing/Seed spacing：32n%+d\nspacing：%d–%d bp，步长32 bp": "The two SST superlattice layers have equal lengths: 32n%+d; Layer spacing/Seed spacing: 32n%+d\nSpacing: %d–%d bp in 32-bp increments",
    "选择工作模式": "Select Work Mode",
    "请选择设计或分析模式。": "Select Design or Analysis mode.",
    "Design · 设计": "Design",
    "Analysis · 分析": "Analysis",
    "新建项目或打开已有 Moiré 项目。": "Create a new project or open an existing Moiré project.",
    "新建项目": "New Project",
    "打开项目": "Open Project",
    "打开设计项目": "Open design project",
    "请先接受 1.1": "Accept 1.1 First",
    "请先接受点阵、Seed截面、Twist与测量环境。": "First accept the lattice, Seed cross-section, twist, and measurement environment.",
    "接受 1.2 设计参数": "Accept 1.2 Design Parameters",
    "请先接受 1.1 与 1.2 参数。": "Accept the 1.1 and 1.2 parameters first.",
    "完整 Seed + SST superlattice 设计已生成，请检查并接受。": "The complete Seed + SST superlattice design has been generated. Review and accept it.",
    "生成 Seed + SST superlattice 设计": "Generate Seed + SST Superlattice Design",
    "✓ 当前设计图已接受": "✓ Current Design Accepted",
    "已接受设计图：%s\n接受时间：%s": "Accepted design: %s\nAccepted at: %s",
    "接受结构设计": "Accept Structure Design",
    "✓ Added Scaffold 已接受": "✓ Added Scaffold Accepted",
    "接受 Added Scaffold": "Accept Added Scaffold",
    "已识别 SST superlattice Input；可自动设计并Add，或打开专家模式。": "SST superlattice inputs detected. Design and add them automatically or open Expert Mode.",
    "没有 Input": "No Inputs",
    "没有识别到SST input位置。": "No SST input positions were detected.",
    "%d-nt 序列分配数量不足。": "Insufficient %d-nt sequences were available for assignment.",
    "自动 Input 序列设计失败": "Automatic Input Sequence Design Failed",
    "自动设计并Add完成：%d条；GC 40–60%%；最大同向相同片段7；最大链间互补片段7。": "Automatic design and addition complete: %d sequences; GC 40–60%%; maximum identical same-orientation segment 7; maximum interstrand complementary segment 7.",
    "SST superlattice Input 已完成": "SST Superlattice Input Complete",
    "已Add %d条Input序列。\n分析工作簿：%s\n填充模板：%s": "Added %d input sequences.\nAnalysis workbooks: %s\nPopulated template: %s",
    "参数：GC %.0f–%.0f%%；最大同向相同片段%d；最大链间互补片段%d\n\n%s": "Parameters: GC %.0f–%.0f%%; maximum identical same-orientation segment %d; maximum interstrand complementary segment %d\n\n%s",
    "自动设计并 Add SST Input 序列": "Automatically Design and Add SST Input Sequences",
    "✓ Added SST superlattice Input 已接受": "✓ Added SST Superlattice Input Accepted",
    "接受 Added SST superlattice Input": "Accept Added SST Superlattice Input",
    "Final Export 完成：%s\n全部序列结果直接保存在 All Sequences 文件夹；PDB CIF oxView Files 文件夹包含单独的 Seed、SST superlattice Layer 1、Layer 2 全原子 PDB/mmCIF、oxDNA TOP/DAT 及纯柱 BILD 模型。\n导出已完成；下一步请打开目标文件夹检查交付文件。": "Final Export complete: %s\nAll sequence results are stored directly in the All Sequences folder. The PDB CIF oxView Files folder contains separate all-atom PDB/mmCIF, oxDNA TOP/DAT, and pure-cylinder BILD models for the Seed, SST superlattice Layer 1, and SST superlattice Layer 2.\nExport is complete. Next: open the destination folder and review the deliverables.",
    "capture延伸长度目前只支持16 nt或32 nt。": "The capture extension currently supports 16 nt or 32 nt.",
    "capture延伸核心staple短于%d nt：": "The capture-extension staple core is shorter than %d nt:",
    "staple+capture短于%d nt（capture为%d nt，其单独staple核心不得短于32 nt）：": "The staple+capture product is shorter than %d nt (capture: %d nt; the staple core must be at least 32 nt):",
    "staple+capture超过64 nt：": "The staple+capture product exceeds 64 nt:",
    "普通staple半crossover回退会删除整对crossover，已拒绝该候选。": "The ordinary-staple half-crossover fallback would remove the full crossover pair; the candidate was rejected.",
    "capture staple半crossover回退会删除整对crossover，已拒绝该候选。": "The capture-staple half-crossover fallback would remove the full crossover pair; the candidate was rejected.",
    "单条staple-capture的32–48-nt硬范围不能满足。": "The mandatory 32–48-nt range for an individual staple-capture cannot be satisfied.",
    "staple crossover回退位点不互反。": "The fallback staple-crossover site is not reciprocal.",
    "已接受参数 · Twist %.1f° · Moiré period %s · SST superlattice %d / %d / %d bp · Seed %d / %d / %d bp": "Accepted parameters · Twist %.1f° · Moiré period %s · SST superlattice %d / %d / %d bp · Seed %d / %d / %d bp",
    "已接受参数 · Twist %+.1f° (%s) · Moiré period %s · SST superlattice %d / %d / %d bp · Seed %d / %d / %d bp": "Accepted parameters · Twist %+.1f° (%s) · Moiré period %s · SST superlattice %d / %d / %d bp · Seed %d / %d / %d bp",
    "自动设计 SST superlattice Input：%d / %d；当前长度 %d nt，已评价 %d 个候选": "Automatically designing SST superlattice inputs: %d / %d; current length %d nt; %d candidates evaluated",
    "%d-nt 序列只生成 %d / %d 条；请使用专家模式调整参数。": "%d-nt sequences: only %d / %d were generated; adjust the parameters in Expert Mode.",
    "SST superlattice Layer 1 / Seed 第一支撑区": "SST superlattice Layer 1 / Seed first support",
    "SST superlattice Layer 2 / Seed 第二支撑区": "SST superlattice Layer 2 / Seed second support",
    "SST superlattice Layer 1 /\nSeed 第一支撑区":
        "SST superlattice Layer 1 /\nSeed first support",
    "SST superlattice Layer 2 /\nSeed 第二支撑区":
        "SST superlattice Layer 2 /\nSeed second support",
    "SST superlattice 第一层 /\nSeed 第一支撑区 · Z1":
        "SST superlattice first layer /\nSeed first support · Z1",
    "SST superlattice 第二层 /\nSeed 第二支撑区 · Z3":
        "SST superlattice second layer /\nSeed second support · Z3",
    "Layer spacing /\nSeed Z2": "Layer spacing /\nSeed Z2",
    "更改 Seed 截面预设": "Change Seed cross-section preset",
    "更改设计参数": "Change design parameters",
    "Add 内置 Scaffold 序列": "Add a Built-in Scaffold Sequence",
    "Add Scaffold 序列": "Add a Scaffold Sequence",
    "Import SST Input 序列": "Import SST Input Sequences",
    "当前截面由预设选择；请使用上方截面选项切换。":
        "The cross-section is selected by a preset; use the cross-section options above to switch presets.",
    "存在%d条超过%d nt的普通staple；最长为%d nt。":
        "%d normal staples exceed %d nt; the longest is %d nt.",
    "存在%d条短于21 nt的普通staple；最短为%d nt。":
        "%d normal staples are shorter than 21 nt; the shortest is %d nt.",
    "识别到 %d 条 scaffold；请在上方逐条 Add scaffold。":
        "Detected %d scaffolds; add each scaffold above.",
    "与SST superlattice spacing是同一个参数，任意一处修改都会同步。":
        "This is the same parameter as SST superlattice spacing; changing either control updates both.",
    "名义spacing加上平均insertion/deletion；SST superlattice spacing与Seed Z2始终共享该实际长度。小数表示不同helices分配整数增删后的平均值。":
        "The actual length is the nominal spacing plus the mean insertion/deletion. SST superlattice spacing and Seed Z2 always share this actual length. A decimal is the mean after integer indels are distributed among helices.",
    "Seed spacing（与SST superlattice联动）":
        "Seed spacing (linked to SST superlattice)",
    "生成固定 SST superlattice + Scaffold routing":
        "Generate Fixed SST Superlattice + Scaffold Routing",
    "固定SST superlattice在生成Scaffold routing时自动加入，无需单独审核。序列与导出模块可载入带序列 JSON，并分别导出Capture，以及只有完整32-nt SST的 JSON、XLSX 和 SVG。":
        "The fixed SST superlattice is added automatically during scaffold-routing generation and does not require separate review. The sequence/export module can load a sequenced JSON and separately export Capture products and an intact 32-nt-SST-only JSON, XLSX, and SVG.",
    "读取后显示 scaffold-only 或 SST superlattice input-only 结构。":
        "After detection, show the scaffold-only or SST superlattice input-only structure.",
    "SST superlattice input-only 视图；相同双层只显示一个导入组，最终仍分别导出两层。":
        "SST superlattice input-only view. Identical bilayers show one import group; both layers are still exported separately.",
    "已接受当前 Added SST superlattice Input；两层真实序列均已写入结构。":
        "The added SST superlattice inputs have been accepted; the real sequences for both layers were written to the structure.",
    "尚未读取 SST superlattice input。":
        "SST superlattice inputs have not been detected.",
    "本轮只开放Square–Square + S8–R4×4C的后续SST superlattice/Scaffold/Capture生成。Kagome或其他Square Seed截面的后续设计将在相应routing规则加入后开放。":
        "Downstream SST superlattice/Scaffold/Capture generation currently supports Square–Square + S8–R4×4C. Kagome and the other Square Seed presets will be enabled after their routing rules are added.",
    "SST superlattice 1st layer和2nd layer至少需要64 bp。":
        "The first and second SST superlattice layers must each be at least 64 bp.",
    "SST superlattice长度不合法": "Invalid SST Superlattice Length",
    "SST superlattice 1st layer、spacing和2nd layer必须是8 bp整数倍。":
        "The first-layer length, spacing, and second-layer length of the SST superlattice must be multiples of 8 bp.",
    "固定两层SST superlattice已在后台生成并验证。\n已生成待审核Scaffold文件：%s\n在cadnano中编辑后请直接保存；接受时会自动读取该文件。\n%s":
        "The fixed two-layer SST superlattice was generated and validated in the background.\nScaffold file ready for review: %s\nAfter editing in cadnano, save that file directly; it will be reloaded when accepted.\n%s",
    "结构已改变，请重新完成 Scaffold 与 SST superlattice input 序列导入。":
        "The structure changed. Add the Scaffold and SST superlattice input sequences again.",
    "结构设计已接受：%s。\n第3步将使用此JSON进行SST superlattice/capture序列设计。":
        "Structure design accepted: %s.\nStep 3 will use this JSON for SST superlattice/capture sequence design.",
    "SST superlattice Input（两层相同，导入一次并自动映射）":
        "SST Superlattice Input (identical layers; import once and map automatically)",
    "SST superlattice Input": "SST Superlattice Input",
    "两层相同：只需导入一次，另一层自动映射。":
        "Identical layers: import once; the other layer is mapped automatically.",
    "请先读取 SST superlattice input。":
        "Detect the SST superlattice inputs first.",
    "两层Capture SST superlattice：%s":
        "Two-layer Capture SST superlattice: %s",
    "Seed–SST superlattice capture桥：%d条；颜色组：%d":
        "Seed–SST superlattice capture bridges: %d; color groups: %d",
    "仍有 %d 条 Scaffold/SST superlattice input 未添加。":
        "%d Scaffold/SST superlattice inputs are still missing.",
    "固定SST superlattice内部验证失败：%s":
        "Internal validation of the fixed SST superlattice failed: %s",
    "SST superlattice Input 读取失败":
        "SST Superlattice Input Detection Failed",
    "普通staple切分位点不互反。":
        "The normal-staple cut sites are not reciprocal.",
    "无法在不移动crossover的前提下断开过长环形staple。":
        "An overlong circular staple cannot be nicked without moving a crossover.",
    "环形staple的nick位点不互反。":
        "The circular-staple nick sites are not reciprocal.",
    "无法在保留全部crossover且两段均为21–58 nt的条件下切分%d-nt普通staple。":
        "The %d-nt normal staple cannot be split into two 21–58-nt segments while preserving every crossover.",
    "Scaffold容量": "Scaffold capacity",
    "接受参数时按合法routing精确核算":
        "Exact capacity will be calculated from the legal routing when the parameters are accepted.",
    "容量不足或无法合法分配：%s":
        "Insufficient capacity or no legal partition: %s",
    "Seed scaffold容量不足": "Insufficient Seed Scaffold Capacity",
    "%s\n\n最多允许3条scaffold（3 × 7557 = 22671 nt）。请减小Seed长度后重新接受参数。":
        "%s\n\nAt most three scaffolds are allowed (3 × 7557 = 22671 nt). Reduce the Seed length and accept the parameters again.",
    "总长度 %d nt → %d条scaffold；计划约为 %s nt":
        "Total length %d nt → %d scaffolds; planned lengths approximately %s nt",
    "Seed同侧边缘错位已最小化，但超过优选的10/11 bp；当前左侧%d bp、右侧%d bp，仍在21-bp硬上限内。":
        "The same-side Seed edge stagger was minimized but exceeds the preferred 10/11 bp: left %d bp, right %d bp; both remain within the 21-bp hard limit.",
})

# English-only Windows release coverage for late-created reports, validation
# errors, worker-process diagnostics, and analysis annotations.  These strings
# are not limited to startup widgets; they can surface after a design is
# generated, while a sequence package is exported, or when validation fails.
_ENGLISH_REVIEWED.update({
    "FFT中的两组同阶twist峰太近或证据不足，无法可靠区分。":
        "The two same-order twist-peak sets in the FFT are too close or have insufficient evidence for reliable separation.",
    "不报告：无法唯一归属于某一层对":
        "Not reported: cannot be assigned uniquely to one layer pair",
    "无法识别：可靠Moiré单元少于2个":
        "Not detected: fewer than two reliable Moiré units",
    "无法识别：两组同阶峰太近或证据不足":
        "Not detected: the two same-order peak sets are too close or have insufficient evidence",
    "仅显示晶格常数 a": "Lattice constant a only",
    "未报告：一阶峰不可可靠分离":
        "Not reported: first-order peaks cannot be separated reliably",
    "Seed 截面": "Seed cross-section",
    "当前已开放 Square–Square、Kagome–Kagome 和 Square–Kagome S8–R4×4C 的后续设计。":
        "Downstream design is available for Square–Square, Kagome–Kagome, and Square–Kagome S8–R4×4C.",
    "已导出全部设计文件：\n1. %s\n2. %s\n3. %s（最终可接受文件）\n%s":
        "All design files were exported:\n1. %s\n2. %s\n3. %s (final file eligible for acceptance)\n%s",
    "没有可接受的最终设计": "No acceptable final design",
    "当前最终文件和选择过的文件中，没有通过验证的 SST + Scaffold + Staple + Capture JSON。":
        "Neither the current final file nor the selected files contain a validated SST + Scaffold + Staple + Capture JSON.",
    "最终结构设计已接受：%s。\nSST 与 SST + Scaffold 是过程导出文件，不作为接受版本。第3步将只使用此最终JSON进行SST superlattice/capture序列设计。":
        "Final structure design accepted: %s.\nThe SST and SST + Scaffold files are intermediate exports and are not accepted versions. Step 3 will use only this final JSON for SST/capture sequence design.",
    "Scaffold 验证失败": "Scaffold validation failed",
    "结构验证失败": "Structure validation failed",

    "Twist预测的名义选区长度必须大于0 bp。":
        "The nominal region length used for twist prediction must be greater than 0 bp.",
    "Twist反算的名义选区长度必须大于0 bp。":
        "The nominal region length used for inverse twist calculation must be greater than 0 bp.",
    "当前Seed截面的弹性Twist响应接近零，无法反算indel。":
        "The elastic twist response of the current Seed cross-section is near zero; the indel cannot be calculated inversely.",
    "当前预览只支持固定288-bp Seed。":
        "The current preview supports only the fixed 288-bp Seed.",
    "Seed使用固定2L范本；scaffold、staple、capture core、nick和seam均不修改，只计算与当前SST的实际重叠。":
        "The Seed uses the fixed 2L template. Scaffold, staples, capture cores, nicks, and seams remain unchanged; only the actual overlap with the current SST is calculated.",
    "当前indel密度折算为96-bp校准区后为%.2f indel/helix，超出已收敛模拟的−10到+6范围；当前弹性模型后的SNUPI校正属于校准域外推，+8和+10 insertion模拟未收敛且未参与拟合。":
        "After conversion to the 96-bp calibration region, the current indel density is %.2f indels per helix, outside the converged simulation range of −10 to +6. The SNUPI correction after the current elastic model is therefore an extrapolation; the +8 and +10 insertion simulations did not converge and were excluded from the fit.",

    "Kagome SST长度和spacing必须位于8-bp网格。":
        "Kagome SST lengths and spacing must lie on the 8-bp grid.",
    "Kagome SST spacing不能小于0 bp。":
        "Kagome SST spacing cannot be less than 0 bp.",
    "Kagome SST全局平移必须是32 bp整数倍。":
        "The global Kagome SST translation must be an integer multiple of 32 bp.",
    "Kagome线型SST活动区间必须是至少32 nt的16-nt整数倍。":
        "An active linear Kagome SST interval must be at least 32 nt and an integer multiple of 16 nt.",
    "Kagome线型SST边界无法形成合法32/48-nt组件。":
        "The linear Kagome SST boundary cannot form a valid 32/48-nt component.",
    "找不到kagome_resource_128.json。":
        "kagome_resource_128.json could not be found.",
    "Kagome SST长度至少64 bp、spacing为0-160 bp，且均须为8 bp整数倍。":
        "Kagome SST lengths must be at least 64 bp and spacing must be 0–160 bp; all values must be multiples of 8 bp.",
    "Kagome SST无法在范本nick相位上形成合法边界组件。":
        "The Kagome SST cannot form a valid boundary component at the template nick phase.",
    "请先将SST-first payload转换为内部helix编号。":
        "Convert the SST-first payload to internal helix numbers first.",
    "Kagome线型SST活动区间必须是至少32 nt的16-nt整数倍。":
        "The active interval of a linear Kagome SST must be at least 32 nt and a multiple of 16 nt.",
    "Kagome SST缺少staple活动区间元数据。":
        "The Kagome SST is missing staple active-interval metadata.",
    "Kagome SST缺少与staple层对应的duplex区间。":
        "The Kagome SST is missing the duplex interval corresponding to the staple layer.",
    "Kagome capture重叠必须按两层提供SST和Seed坐标。":
        "Kagome capture overlap requires SST and Seed coordinates for both layers.",
    "Kagome Seed支撑区必须恰好包含两层。":
        "The Kagome Seed support region must contain exactly two layers.",
    "Kagome SST待切连接不是唯一互反连接。":
        "The Kagome SST connection to be cut is not a unique reciprocal connection.",
    "Kagome SST缺少lattice_type=kagome元数据。":
        "The Kagome SST is missing lattice_type=kagome metadata.",
    "Kagome SST必须包含两层duplex range。":
        "The Kagome SST must contain duplex ranges for both layers.",
    "Kagome SST组件必须恰有一个5′端。":
        "A Kagome SST component must have exactly one 5′ end.",
    "Kagome边缘状态不唯一：%s -> %s。":
        "The Kagome edge state is not unique: %s -> %s.",
    "Kagome capture延伸产生非16/32/48-nt组件：%s。":
        "The Kagome capture extension produced a component that is not 16, 32, or 48 nt: %s.",
    "Kagome SST helix集合错误：应为%s。":
        "The Kagome SST helix set is incorrect; expected %s.",
    "Kagome SST起始双链位置必须大于等于32。":
        "The starting duplex position of the Kagome SST must be at least 32.",
    "Kagome SST层间距与Z2不一致。":
        "The Kagome SST layer spacing does not match Z2.",
    "Kagome %s绝对相位非法：%d。":
        "The absolute phase of Kagome %s is invalid: %d.",
    "Kagome SST模板包含跨出SST截面的连接。":
        "The Kagome SST template contains a connection outside the SST cross-section.",
    "%s逻辑活动helix错误：%s。":
        "The logical active helices for %s are incorrect: %s.",
    "%s含非法SST组件长度：%s。":
        "%s contains invalid SST component lengths: %s.",
    "Kagome SST duplex长度与参数不一致。":
        "The Kagome SST duplex length does not match the parameters.",
    "Kagome SST非法capture候选 h%d:%d slot%d=%s。":
        "Invalid Kagome SST capture candidate h%d:%d slot%d=%s.",
    "Kagome辅助capture helix不存在：h%d。":
        "The auxiliary Kagome capture helix does not exist: h%d.",
    "Kagome空位helix %d不得含有strand。":
        "Vacant Kagome helix %d must not contain a strand.",
    "%s越界连接 h%d:%d。":
        "%s has an out-of-bounds connection at h%d:%d.",
    "%s非互反连接 h%d:%d。":
        "%s has a nonreciprocal connection at h%d:%d.",

    "Square SST两层长度均须至少为64 bp。":
        "Both Square SST layer lengths must be at least 64 bp.",
    "Square SST spacing必须位于0-160 bp。":
        "Square SST spacing must be between 0 and 160 bp.",
    "Square SST长度与spacing必须采用8 bp步长。":
        "Square SST lengths and spacing must use 8-bp increments.",
    "Square SST边界必须位于8 bp网格。":
        "Square SST boundaries must lie on the 8-bp grid.",
    "当前SST位置无法在固定288-bp Seed内定义完整Z1/Z2/Z3分区。":
        "The current SST position cannot define complete Z1/Z2/Z3 regions within the fixed 288-bp Seed.",
    "Square SST公共画布平移必须保持32-bp范本相位。":
        "The shared-canvas translation of the Square SST must preserve the 32-bp template phase.",
    "Square SST边界与完整32-nt U型范本相位不一致。":
        "The Square SST boundary is inconsistent with the phase of the complete 32-nt U-shaped template.",
    "Square capture绝对相位原点发生偏移。":
        "The absolute-phase origin of the Square capture sites has shifted.",
    "固定288-bp Seed内无法放置当前SST spacing。":
        "The current SST spacing cannot be placed within the fixed 288-bp Seed.",
    "Square capture列未保持固定Seed范本相位。":
        "The Square capture columns do not preserve the fixed Seed-template phase.",
    "Square capture列不在固定16-bp接触相位上。":
        "The Square capture columns are not on the fixed 16-bp contact phase.",
    "Square SST边界相位不可解析。":
        "The Square SST boundary phase cannot be resolved.",
    "SST辅助绕行后仍有未解决的主通道占位冲突。":
        "The auxiliary SST detour left an unresolved occupancy conflict in the main channel.",
    "SST辅助绕行遇到缺失的layer-2互反端点。":
        "The auxiliary SST detour encountered a missing reciprocal endpoint in layer 2.",

    "Seed scaffold总长度必须大于0 nt。":
        "The total Seed scaffold length must be greater than 0 nt.",
    "Seed scaffold总长度%d nt超过3条正交scaffold的总容量%d nt；请减小Seed长度。":
        "The total Seed scaffold length of %d nt exceeds the %d-nt capacity of three orthogonal scaffolds; reduce the Seed length.",
    "SST必须恰好包含两层实际双链范围。":
        "The SST must contain exactly two actual duplex ranges.",
    "capture合法位点必须按两层提供。":
        "Valid capture sites must be provided for both layers.",
    "Seed支撑区与capture网格必须恰好包含两层。":
        "The Seed support regions and capture grid must contain exactly two layers.",
    "Square–Kagome第二层无法在合法Kagome相位形成32/48-nt组件。":
        "The second Square–Kagome layer cannot form a 32/48-nt component at a valid Kagome phase.",
    "双层设计与序列一致时，Z1/Z2/Z3必须满足Square 32-bp相位联动规则。":
        "When the two layer designs and sequences are identical, Z1/Z2/Z3 must satisfy the Square 32-bp phase-linkage rule.",
    "helix 48–63及辅助helix 64–79":
        "helices 48–63 and auxiliary helices 64–79",
    "SST-only中的scaffold必须是本点阵合法、互反的32/48-nt单元。":
        "Scaffold components in SST-only must be valid reciprocal 32/48-nt units for the selected lattice.",
    "SST-only中的SST必须是本点阵合法、互反的32/48-nt单元；Kagome 16-nt仅允许范本相位的64-nt边界48+16特例。":
        "SST components in SST-only must be valid reciprocal 32/48-nt units for the selected lattice. A 16-nt Kagome component is allowed only as the template-phase 48+16 boundary case of a 64-nt unit.",
    "范本中有%d条超过%d nt的普通staple；最长为%d nt。":
        "The template contains %d normal staples longer than %d nt; the longest is %d nt.",
    "每层至少需要2列capture pair；当前为%s。固定2L Seed不会缩进、增长或改变routing，请减小SST spacing。":
        "Each layer requires at least two capture-pair columns; the current value is %s. The fixed 2L Seed will not contract, extend, or change routing; reduce the SST spacing.",
    "SST活动区间超出caDNAno坐标范围。":
        "The SST active interval exceeds the caDNAno coordinate range.",
    "Square SST组件必须恰有一个5′端。":
        "A Square SST component must have exactly one 5′ end.",
    "Kagome双层设计与序列一致时，Z1/Z2/Z3必须满足32-bp相位联动规则。":
        "When the Kagome layer designs and sequences are identical, Z1/Z2/Z3 must satisfy the 32-bp phase-linkage rule.",
    "未知SST点阵类型：%s。": "Unknown SST lattice type: %s.",
    "Kagome SST缺少variable_length_layout元数据。":
        "The Kagome SST is missing variable_length_layout metadata.",
    "Kagome SST缺少理论capture候选目录。":
        "The Kagome SST is missing the theoretical capture-candidate catalog.",
    "SST文件必须只包含%s，且不能有重复编号。":
        "The SST file must contain only %s and must not contain duplicate helix numbers.",
    "helix编号必须为0–63%s。":
        "Helix numbers must be 0–63%s.",
    "Seed scaffold存在%d个非法crossover（非相邻%d、端点索引不一致%d、Square相位非法%d、非互反%d、同helix不同向间距小于8 bp %d）。":
        "The Seed scaffold contains %d invalid crossovers (nonadjacent: %d; mismatched endpoint indices: %d; invalid Square phase: %d; nonreciprocal: %d; opposite-direction spacing below 8 bp on one helix: %d).",
    "Seed的最大实际Z1/Z3长度必须由capture helix定义；检测到普通helix伸得更远。":
        "The maximum actual Z1/Z3 lengths of the Seed must be defined by capture helices; a normal helix extends farther.",
    "Seed Z1/Z3超过11 bp但没有记录8-bp输入参数经合法scaffold相位量化产生的边缘增长例外。":
        "Seed Z1/Z3 exceeds 11 bp without a recorded edge-growth exception produced by quantizing the 8-bp input parameters to a valid scaffold phase.",
    "Seed scaffold元数据分段数%d与总长度%d nt规定的%d条不一致。":
        "The %d Seed scaffold segments recorded in metadata do not match the %d scaffolds required by the total length of %d nt.",
    "scaffold超过当前分段允许的%d nt：%s":
        "A scaffold exceeds the %d-nt limit for the current segmentation: %s",
    "Seed scaffold总长度%d nt必须使用%d条；当前为%d条。":
        "A total Seed scaffold length of %d nt requires %d scaffolds; the current design has %d.",
    "保留%d条无法安全并入的Seed边缘短staple（相邻nick用于capture延伸）；最短为%d nt。":
        "Retained %d short Seed-edge staples that could not be merged safely because the adjacent nick is used for capture extension; the shortest is %d nt.",
    "capture base不得与Seed内部staple crossover共位；发现%s。":
        "A capture base must not coincide with an internal Seed staple crossover; found %s.",
    "Square SST capture端点没有唯一的U型crossover边。":
        "The Square SST capture endpoint does not have a unique U-shaped crossover edge.",
    "Square SST capture待切边不是唯一互反边。":
        "The Square SST capture edge to be cut is not a unique reciprocal edge.",
    "Square–Kagome capture待切边不是唯一互反边。":
        "The Square–Kagome capture edge to be cut is not a unique reciprocal edge.",
    "helix %d 的SST scaffold必须恰好覆盖本层活动区间：%s。":
        "The SST scaffold on helix %d must cover exactly the active interval of this layer: %s.",
    "，并在启用辅助绕行时包含64–79":
        ", and include 64–79 when the auxiliary detour is enabled",
    "Seed同侧边缘错位必须不超过21 bp；当前左侧%d bp、右侧%d bp。":
        "Same-side Seed edge offsets must not exceed 21 bp; current left offset: %d bp, right offset: %d bp.",
    " 固定2L Seed不重新break。":
        " The fixed 2L Seed is not broken again.",
    "%s未覆盖应与%s SST直接连接的%d根Seed helix；缺少%s。":
        "%s does not cover the %d Seed helices that should connect directly to the %s SST; missing %s.",
    "Seed已能在11-bp优选范围内完成，禁止保留无收益的21-bp放宽状态。":
        "The Seed can already be completed within the preferred 11-bp range; an unhelpful 21-bp relaxation must not be retained.",
    "冻结2L范本边缘错位为%d bp；保留范本routing，仍满足21-bp硬上限。":
        "The frozen 2L template has a %d-bp edge offset. The template routing is retained and still satisfies the 21-bp hard limit.",
    "Seed只有在21-bp候选能减少scaffold数量时才允许同侧错位超过11 bp；当前左侧%d bp、右侧%d bp，但没有记录有效的减链收益。":
        "A same-side Seed offset above 11 bp is allowed only when a 21-bp candidate reduces the scaffold count. Current left offset: %d bp; right offset: %d bp; no valid scaffold-count reduction was recorded.",
    "Seed为减少scaffold数量而使用21-bp合法放宽；当前左侧%d bp、右侧%d bp。":
        "The Seed uses the valid 21-bp relaxation to reduce the scaffold count; current left offset: %d bp, right offset: %d bp.",

    "第2步只支持8×8 Seed减4×4 pore；不支持其他Seed截面。":
        "Step 2 supports only an 8×8 Seed with a 4×4 pore; other Seed cross-sections are not supported.",
    "Seed scaffold总长度%d nt超过3条正交scaffold的容量%d nt；请减小Seed长度。":
        "The total Seed scaffold length of %d nt exceeds the %d-nt capacity of three orthogonal scaffolds; reduce the Seed length.",
    "capture base与Seed内部staple crossover共位；这表示SST/capture的32-bp整体平移或相位错误。不得删除或移动AutoCS crossover。":
        "A capture base coincides with an internal Seed staple crossover, indicating an incorrect 32-bp translation or phase of the SST/capture pattern. AutoCS crossovers must not be deleted or moved.",
    "没有%s %d-bp的冻结SST范本。":
        "No frozen %s SST template is available at %d bp.",

    "冻结2L Seed不再满足两条scaffold容量。":
        "The frozen 2L Seed no longer fits within the capacity of two scaffolds.",
    "Kagome固定capture投影意外改动Seed普通staple。":
        "The fixed Kagome capture projection unexpectedly modified normal Seed staples.",
    "不支持的SST点阵类型：%s。": "Unsupported SST lattice type: %s.",
    "固定Seed只能按非负32-bp完整周期整体平移；当前为%d bp。":
        "The fixed Seed can be translated only by a nonnegative whole number of 32-bp periods; current translation: %d bp.",
    "固定Seed Z2缺少安全indel位点：%s。":
        "The fixed Seed Z2 lacks safe indel sites: %s.",
    "Z2 insertion需要%d nt，但两条固定scaffold仅剩%d nt容量；每条scaffold必须≤7557 nt。":
        "The Z2 insertion requires %d nt, but the two fixed scaffolds have only %d nt of remaining capacity; each scaffold must be ≤7557 nt.",
    "Z2安全deletion位点不足：需要%d，只有%d。":
        "Insufficient safe Z2 deletion sites: %d required, %d available.",
    "固定Seed scaffold已超过7557 nt，不能继续加入insertion。":
        "A fixed Seed scaffold already exceeds 7557 nt; no further insertion can be added.",
    "Z2 indel后scaffold超过7557 nt：%s。":
        "A scaffold exceeds 7557 nt after the Z2 indel: %s.",
    "冻结2L Seed routing已被改动：scaffold lengths=%r。":
        "The frozen 2L Seed routing was modified: scaffold lengths=%r.",
    "冻结2L Seed边缘范围异常：low=%d high=%d。":
        "The frozen 2L Seed edge range is invalid: low=%d, high=%d.",
    "每层至少需要2列capture pair；当前为%s。请减小SST spacing。":
        "Each layer requires at least two capture-pair columns; the current value is %s. Reduce the SST spacing.",
    "已接受Seed scaffold与冻结2L范本不一致；请先重新生成或在专家模式确认，不得静默覆盖。":
        "The accepted Seed scaffold is inconsistent with the frozen 2L template. Regenerate it or confirm it in Expert Mode; it must not be overwritten silently.",
    "已接受Seed scaffold与冻结2L范本不一致；禁止为Kagome SST静默覆盖Seed。":
        "The accepted Seed scaffold is inconsistent with the frozen 2L template. The Seed must not be overwritten silently for a Kagome SST.",
    "固定Seed必须含1或2条scaffold，当前检测到%d条。":
        "The fixed Seed must contain one or two scaffolds; %d were detected.",
    "冻结Seed范本的每条staple必须且只能有一个颜色标记；当前组件有%d个。":
        "Each staple in the frozen Seed template must have exactly one color marker; the current component has %d.",
    "冻结Seed capture核心无法唯一映射到真实列：%s。":
        "A frozen Seed capture core cannot be mapped uniquely to an actual column: %s.",
    "一条Seed capture核心跨越了多个真实列：%s。":
        "A Seed capture core spans multiple actual columns: %s.",
    "Seed capture列base %d没有独立显示颜色。":
        "Base %d of a Seed capture column has no distinct display color.",
    "生成的Seed staple无法唯一映射到冻结范本颜色；检测到%d种候选颜色。":
        "A generated Seed staple cannot be mapped uniquely to a frozen-template color; %d candidate colors were detected.",
    "冻结Seed范本包含无法映射到capture列的颜色#%06x。":
        "The frozen Seed template contains color #%06x, which cannot be mapped to a capture column.",
    "Kagome capture位点的Seed端没有可用nick：h%d:%d。":
        "The Seed end of Kagome capture site h%d:%d has no available nick.",
    "Kagome SST capture端点未打开：h%d:%d。":
        "Kagome SST capture endpoint h%d:%d is not open.",
    "Kagome capture位点的Seed端没有合法开放槽或纵向连接：h%d:%d。":
        "The Seed end of Kagome capture site h%d:%d has neither a valid open slot nor a longitudinal connection.",
    "capture候选h%d[%d]落在AutoCS staple crossover碱基上；必须调整capture相位/位点，禁止在该base连接。":
        "Capture candidate h%d[%d] lies on an AutoCS staple-crossover base. Adjust the capture phase/site; connection at this base is prohibited.",
    "capture位点h%d[%d]不是合法纵向端点。":
        "Capture site h%d[%d] is not a valid longitudinal endpoint.",

    "最低密度（最大可用间距）": "Minimum density (maximum available spacing)",
    "已用快速确定性路线生成 1/%d bp scaffold 布局：%d 个 crossover，主闭环 %d nt。":
        "Generated a 1/%d-bp scaffold layout using the fast deterministic route: %d crossovers; main closed loop, %d nt.",
    "左侧几何视图与右侧 Path 顺序中没有同时相邻的 scaffold helix；未修改设计。":
        "No scaffold helices are adjacent in both the left geometry view and the right Path order; the design was not modified.",
    "仅在左侧和右侧均相邻的 %d 对 helix 之间，按%s添加 %d 个 scaffold crossover；删除 %d 个原 crossover。两端未使用的 scaffold 保持不变。":
        "Added %d scaffold crossovers using %s only between the %d helix pairs adjacent in both views; removed %d existing crossovers. Unused scaffold at both ends remains unchanged.",
    "当前路线未通过唯一主闭环硬准则，已恢复运行前设计；未保留局部环或临时 crossover。":
        "The current route failed the unique-main-loop hard criterion. The pre-run design was restored; no local loops or temporary crossovers were retained.",
    "当前点阵不支持 AutoCS_scaffolds。":
        "The current lattice does not support AutoCS_scaffolds.",
    "当前 scaffold 区域没有任何合法 crossover 位点；未修改设计。":
        "The current scaffold region has no valid crossover sites; the design was not modified.",
    "已按确定性路线并执行单闭环收尾，%s生成 %d 个 scaffold crossover（当前共 %d 个）。":
        "Using a deterministic route followed by single-loop closure, %s generated %d scaffold crossovers (%d total).",
    "首选 1/%d bp 路线无法通过唯一闭环硬准则；已选择距离最近且合法的 1/%d bp 闭环。":
        "The preferred 1/%d-bp route failed the unique-loop hard criterion; the nearest valid 1/%d-bp closed loop was selected.",
    "所有候选均未通过唯一主闭环硬准则，已恢复运行前设计；未保留局部环或临时 crossover。":
        "All candidates failed the unique-main-loop hard criterion. The pre-run design was restored; no local loops or temporary crossovers were retained.",
    "采用最低密度": "using minimum density",
    "已统一按“闭环、密度、Path 顺序/模块、纵向桥、seam”排序 %d 个候选。":
        "Ranked %d candidates consistently by loop closure, density, Path order/module, longitudinal bridge, and seam.",
    "优先采用 1/%d bp；精确周期不能闭环时已在 1/%d bp 硬上限内回退":
        "prefer 1/%d bp; when the exact period cannot close, fall back within the 1/%d-bp hard limit",
    "以 1/%d bp 密度上限": "with a 1/%d-bp density limit",
})

# Application-wide English copy review.  Keep this separate from the staged
# structure-design wording above: it covers startup, Design step 1, sequence
# assignment, Gel Analysis, and Crystal/Particle Analysis.  Design step 2 is
# intentionally left to its parallel implementation task.
_ENGLISH_REVIEWED.update({
    # Menus, startup, projects, and Design step 1.
    "语言": "Language",
    "界面语言：中文": "Interface language: Chinese",
    "文件": "File",
    "编辑": "Edit",
    "视图": "View",
    "工具": "Tools",
    "创建项目": "Create Project",
    "项目名": "Project name",
    "项目名不合法": "Invalid Project Name",
    "项目已存在": "Project Already Exists",
    "文件名": "File name",
    "选择工作模式": "Select a Work Mode",
    "设计与预测": "Design and Prediction",
    "Design · 设计": "Design",
    "Analysis · 分析": "Analysis",
    "Seed 截面固定为 8×8 + 4×4 pore，不可选择或编辑。":
        "The Seed cross-section is fixed at 8×8 with a 4×4 pore and cannot be edited.",
    "参数已改变；接受时重新精确核算":
        "Parameters changed; exact values will be recalculated when accepted",
    "实时更新：名义Z2=%d bp，实际Z2/spacing=%.1f bp，Twist %+.1f° (%s)，period %s":
        "Live update: nominal Z2=%d bp, actual Z2/spacing=%.1f bp, "
        "Twist %+.1f° (%s), period %s",
    "可进入1B设计。": "Ready for Design 1B.",
    "接受 1A 设计参数": "Accept Design 1A Parameters",
    "请先接受 1A": "Accept Design 1A first",
    "设计参数验证失败": "Design Parameter Validation Failed",
    "接受 1B 设计参数": "Accept Design 1B Parameters",
    "两层均为Square；可由Twist与Square a计算Moiré period。":
        "Both layers use the Square lattice; the Moiré period is calculated from twist and the Square lattice constant a.",
    "1st layer为Square，2nd layer为Kagome；可设置Twist，但不同点阵之间不定义Moiré period。":
        "Layer 1 uses the Square lattice and Layer 2 uses the Kagome lattice. Twist can be specified, but a single Moiré period is not defined for unlike lattices.",
    "Square–Kagome使用两个不同点阵常数；不定义单一Moiré period。":
        "Square–Kagome uses two different lattice constants; a single Moiré period is not defined.",
    "—（不同点阵）": "— (unlike lattices)",
    "100×100根helix俯视图；左键拖动旋转；右键拖动平移；滚轮缩放；双击复位":
        "100 × 100-helix top view · Left-drag to rotate · Right-drag to pan · Scroll to zoom · Double-click to reset",
    "左拖旋转 · 右拖平移 · 滚轮缩放 · 双击复位":
        "Left-drag to rotate · Right-drag to pan · Scroll to zoom · Double-click to reset",
    "当前组合仅支持设计预测": "This lattice combination currently supports prediction only",
    "理论双层重构将在这里显示": "The theoretical bilayer reconstruction will appear here",
    "干燥 TEM": "Dried TEM",
    "溶液 / cryo-EM": "Solution / cryo-EM",
    "平均 lattice constant": "Mean lattice constant",
    "平均 insertion/deletion": "Mean insertion/deletion",
    "最终采用Twist": "Final twist",
    "预测角度": "Predicted angle",
    "目标": "Target",
    "路径": "Path",
    "适合窗口": "Fit to Window",
    "↶ 撤回": "↶ Undo",

    # Sequence assignment and orthogonal-sequence design.
    "专家流程：先导出模板，也可用正交序列设计生成序列；填入 Sequence 列后再导入。行顺序按完整 base 数值优先、helix 数值其次排列。":
        "Expert workflow: export the template, enter sequences in the Sequence column, and import it. Sequences may also be generated with Orthogonal Sequence Design. Rows are sorted first by the full base index and then by helix number.",
    "已识别 SST superlattice Input；可自动设计并 Add，或打开专家模式，依次 Export template、正交序列设计、Import sequences，最后在右侧接受当前 Input。":
        "SST superlattice inputs detected. Design and add them automatically, or open Expert Mode and use Export Input Template, Orthogonal Sequence Design, and Import Input Sequences in order. Finally, accept the current inputs on the right.",
    "已识别 SST superlattice Input。可自动设计并 Add，或使用专家模式依次导出模板、设计正交序列并导入。完成后在下方接受当前 Input。":
        "SST superlattice inputs detected. Design and add them automatically, or use Expert Mode to export the template, design orthogonal sequences, and import them. Then accept the current inputs below.",
    "自动 Input 序列设计未完成：%s":
        "Automatic input-sequence design did not complete: %s",
    "已接受当前 Input；两层真实序列均已写入结构。":
        "The current inputs have been accepted; the actual sequences for both layers were written to the structure.",
    "自动读取 Scaffold 位置和长度": "Detect Scaffold Positions and Lengths",
    "正在读取 Scaffold 位置和长度…": "Detecting scaffold positions and lengths…",
    "自动识别 SST superlattice Input 位置和长度":
        "Detect SST Superlattice Input Positions and Lengths",
    "自动设计并 Add Input 序列": "Design and Add Input Sequences Automatically",
    "自动设计并 Add SST Input 序列": "Design and Add SST Input Sequences Automatically",
    "自动设计 SST superlattice Input：%d / %d；当前长度 %d nt，已评价 %d 个候选":
        "Designing SST superlattice inputs: %d/%d · Current length: %d nt · Candidates evaluated: %d",
    "自动设计并Add完成：%d条；GC 40–60%%；最大同向相同片段7；最大链间互补片段7。":
        "Automatic design complete: %d sequences added · GC: 40–60%% · Maximum identical same-orientation segment: 7 nt · Maximum interstrand complementary segment: 7 nt.",
    "正交序列设计完成": "Orthogonal Sequence Design Complete",
    "正交序列设计未完全完成": "Orthogonal Sequence Design Incomplete",
    "正交序列设计无法打开": "Cannot Open Orthogonal Sequence Design",
    "序列导入完成，可导出 Final Export。":
        "Sequence assignment is complete. The final package is ready to export.",
    "完成序列导入并导出 Final Export":
        "Export Final Package",
    "序列工作流没有完成。": "The sequence workflow did not complete.",
    "序列工作流没有返回有效报告。": "The sequence workflow returned no valid report.",
    "没有 Input": "No inputs",
    "尚未 Add scaffold": "No scaffold sequence assigned",
    "没有识别到 scaffold。": "No scaffolds were detected.",
    "没有识别到SST input位置。": "No SST input positions were detected.",
    "SST input未完整导入：需要%d条，当前%d条。":
        "SST input import is incomplete: %d strands required, %d currently assigned.",
    "第%d行包含非法碱基：%s。": "Row %d contains invalid bases: %s.",
    "第%d行需要%d nt，当前为%d nt。": "Row %d requires %d nt; the imported sequence has %d nt.",
    "第%d行End与当前结构不一致。": "The End value in row %d does not match the current structure.",
    "第%d行Start不属于当前SST input：%s":
        "The Start value in row %d does not belong to a current SST input: %s",
    "第%d行": "Row %d",
    "全局 GC": "Overall GC Content",
    "局部 GC": "Local GC Content",

    # Gel Analysis.
    "等待 gel 图片。": "Waiting for a gel image.",
    "选择 gel 图片": "Select a Gel Image",
    "缺少 gel 图片": "Gel Image Required",
    "请先上传 gel 图片。": "Upload a gel image first.",
    "显示黑白": "Show in Grayscale",
    "不展平背景": "No background flattening",
    "处理预览": "Processing Preview",
    "Lane 空白处单击选中 Lane、拖动新建 band · 点击选框后拖动/拉伸 · 右键赋值 · 滚轮缩放 · 中键平移":
        "Click an empty area of a lane to select it; drag to create a band box · Drag a selected box to move or resize it · Right-click to assign a role · Scroll to zoom · Middle-drag to pan",
    "滚轮缩放 · 中键平移": "Scroll to zoom · Middle-drag to pan",
    "右键赋值 Lane / Target / Staples · 拖动或拉伸调整 · 靠近边界自动吸附":
        "Right-click to assign Lane / Target / Staples · Drag to move or resize · Snap to nearby boundaries",
    "设为 Lane": "Assign as Lane",
    "设为 Target band": "Assign as Target Band",
    "设为 Staples band": "Assign as Staples Band",
    "设为普通 band": "Assign as Unclassified Band",
    "清除赋值": "Clear Assignment",
    "删除选框": "Delete Selection Box",
    "删除选中 Delete": "Delete Selected",
    "清空右侧所有选框": "Clear All Selection Boxes",
    "缺少 Lane": "Lane Required",
    "请先右键把至少一个选框赋值为 Lane。":
        "Right-click and assign at least one selection box as a lane.",
    "自动识别已完成 · 请在右侧人工校订":
        "Automatic detection complete · Review and correct the annotations on the right",
    "正在分析 lanes 与 bands…": "Analyzing lanes and bands…",
    "带 lane / band 标注的 gel 图": "Annotated Gel Image",
    "查看单 lane 曲线": "View Single-Lane Profile",
    "逐 lane band 面积、比例与 target yield":
        "Band Areas, Fractions, and Target Yield by Lane",
    "自动寻峰灵敏度": "Peak-detection sensitivity",
    "自动最小峰宽": "Minimum automatic peak width",
    "背景扣除半径": "Background-subtraction radius",
    "曲线平滑": "Profile smoothing",
    "非-staple 比例": "Non-staple fraction",
    "分析已保存": "Analysis Saved",
    "已导出全部结果（%d 个文件）：%s": "Exported all results (%d files): %s",

    # Crystal / Particle Analysis setup and interaction.
    "1. 上传 TEM 图；可以‘<数值> <单位>_’开头命名（如 20 nm_xx）\n2. 选择单层或双层分析\n3. 点击分析，自动识别 scale bar、TEM 与 FFT\n4. 核对尺度与结果；无 OCR/文件名尺度时才手动输入":
        "1. Upload a TEM image; the file name may begin with '<value> <unit>_' (for example, 20 nm_sample)\n2. Select single-layer or bilayer analysis\n3. Click Analyze to detect the scale bar and analyze the TEM image and FFT automatically\n4. Review the scale and results; enter the scale manually only when neither OCR nor the file name provides it",
    "2. 选择分析层数": "2. Select Single-Layer or Bilayer Analysis",
    "1. 上传图像": "1. Upload Image",
    "2. 图像处理与自动尺度": "2. Image Processing and Automatic Scaling",
    "3. 形状与尺寸": "3. Shape and Size",
    "4. 框选代表颗粒": "4. Draw a Box Around a Representative Particle",
    "5. 自动选择颗粒并分析": "5. Automatic Particle Selection and Analysis",
    "5. 自动选择颗粒并分析（可选）":
        "5. Automatic Particle Selection and Analysis (Optional)",
    "6. 手动选择颗粒并分析": "6. Manual Particle Selection and Analysis",
    "开始自动识别与统计": "Run Automatic Detection and Quantification",
    "1. 上传 TEM 图像": "1. Upload TEM Image",
    "2. 选择单层或双层分析": "2. Select Single-Layer or Bilayer Analysis",
    "3. 自动识别并分析": "3. Automatically identify and analyze",
    "3. 批量自动识别并分析": "3. Run batch detection and analysis",
    "4. 结果摘要": "4. Result Summary",
    "4. Scale bar 与结果摘要": "4. Scale Bar and Result Summary",
    "4. 批量运行状态": "4. Batch Run Status",
    "开始分析": "Run Analysis",
    "开始批量分析": "Run Batch Analysis",
    "3. 可选：框一个代表颗粒": "3. Optional: Draw a Box Around a Representative Particle",
    "4. 开始自动识别与统计": "4. Run Automatic Detection and Quantification",
    "单层分析": "Single-Layer Analysis",
    "双层分析": "Bilayer Analysis",
    "单层分析不计算Moiré周期。": "Single-layer analysis does not calculate the Moiré period.",
    "等待TEM图像。": "Waiting for a TEM image.",
    "等待表征图。": "Waiting for a microscopy image.",
    "必需：选择一张 TEM、AFM 或 SEM 图像": "Required: select one TEM, AFM, or SEM image",
    "选择 TEM / AFM / SEM 图": "Select a TEM, AFM, or SEM Image",
    "原始 TEM / AFM / SEM 图": "Original TEM / AFM / SEM Image",
    "当前图像结果与标尺校正": "Current Image Results and Scale Calibration",
    "scale bar": "Scale bar",
    "标尺横线": "Scale-bar line",
    "手动标尺横线 + 用户标称长度": "Manual scale-bar line + user-specified length",
    "自动标尺横线 + 用户标称长度": "Automatically detected scale-bar line + user-specified length",
    "无法识别scale bar横线。请输入该图scale bar的像素长度：":
        "The scale-bar line could not be detected. Enter its length in pixels:",
    "无法读取scale bar标注数值。请输入该图对应的实际nm值：":
        "The scale-bar label could not be read. Enter the represented length in nm:",
    "OCR未读出 scale bar 标注。请输入该 scale bar 代表的实际长度（nm）：":
        "OCR could not read the scale-bar label. Enter the represented length in nm:",
    "OCR=%g nm，%s=%g nm；已优先使用 OCR，请核对。":
        "OCR: %g nm; %s: %g nm. OCR was used; verify the value.",
    "已找到图中的标尺横线（约 %.1f px），但未读出标注。\n请输入这条 scale bar 代表的实际长度：":
        "A scale-bar line was detected (approximately %.1f px), but its label could not be read.\nEnter the represented length:",
    "文件中未读到物理像素尺度。请输入每像素对应的 nm：":
        "No physical pixel scale was found in the file metadata. Enter the scale in nm per pixel:",
    "文件名缺少尺度": "Scale Missing from File Name",
    "无法确定尺度": "Cannot Determine Scale",
    "需要标尺信息": "Scale Information Required",
    "确认 scale bar": "Confirm Scale Bar",
    "用于含物理像素尺度元数据的原始 TIFF。不读取图内 scale bar，直接使用文件内置尺度。":
        "For raw TIFF files containing physical pixel-scale metadata. The embedded scale is used directly; an image scale bar is not read.",
    "所有图片的 scale bar 数值相同": "All images share the same scale-bar value",
    "输入所有图共用的 scale bar 标注值，再批量选图":
        "Enter the scale-bar value shared by all images, then select the batch",
    "按‘<数值> <单位>_文件名’为每张图标注尺度":
        "Encode each image scale as '<value> <unit>_filename'",
    "Bulk scale 来源": "Bulk Scale Source",
    "选择 Bulk scale 模式": "Select a Bulk Scale Mode",
    "Bulk analysis在左侧逐张运行并自动导出；右侧不加载批量图片。":
        "Bulk analysis processes and exports images sequentially from the left panel; batch images are not loaded into the right preview.",
    "Bulk完成：%d个成功，%d个失败。\n自动导出：%s":
        "Bulk analysis complete: %d succeeded, %d failed.\nExports: %s",
    "已逐张分析并自动导出。\n成功：%d\n失败：%d\n\n%s":
        "Batch analysis and export complete.\nSucceeded: %d\nFailed: %d\n\n%s",
    "选择Bulk analysis自动导出文件夹": "Select the Bulk Analysis Export Folder",
    "批量自动识别并分析": "Run Batch Detection and Analysis",
    "批量分析已导出": "Batch Analysis Exported",
    "批量分析失败": "Batch Analysis Failed",
    "开启后可一次选择多个TEM样品，并导出CSV汇总统计。":
        "Enable this option to select multiple TEM samples and export summary statistics as CSV.",
    "开启后可选择多个TEM样品；运行前选择文件夹，随后逐图自动导出并生成CSV汇总。":
        "Enable this option to select multiple TEM samples. Choose an output folder before running; each image is exported automatically and a summary CSV is generated.",
    "自动分组数量": "Number of Automatic Groups",
    "自动建议": "Auto (recommended)",
    "请输入": "Enter a value",
    "代表模板匹配门槛": "Representative-Template Matching Threshold",
    "统计哪个尺寸": "Dimension to Report",
    "面积密度定义": "Area-Density Definition",
    "材料": "Material Area",
    "外包络": "Outer Envelope",
    "不规则外框（手绘）": "Irregular Outline (Freehand)",
    "旋转矩形": "Rotated Rectangle",
    "棒 / 长条": "Rod / Strip",
    "棒/长条": "Rod/Strip",
    "菱形": "Diamond",
    "自动紧密贴合": "Automatic Tight Fit",
    "代表颗粒贴合完成": "Representative-Particle Fit Complete",
    "框内未找到可靠颗粒": "No Reliable Particle Found in the Box",
    "无法自动贴合": "Automatic Fit Failed",
    "处理结果预览（阈值应用后的同一张图；代表框与后期新增都使用这里保留的像素）":
        "Processing Preview (after thresholding; representative boxes and subsequent additions use the pixels retained here)",
    "底图是本次分析冻结的处理通道；半透明色块是实际统计区域；选择 Group 后只显示该组":
        "The base image is the processing channel frozen for this analysis. Translucent overlays show the regions used for quantification; select a Group to display that Group only.",
    "当前显示全部 Group：外框给出形状，半透明色块是实际统计掩膜。":
        "Showing all Groups: outlines indicate shape; translucent overlays are the masks used for quantification.",
    "全部纳入统计颗粒": "All Particles Included in Statistics",
    "错误识别 Group": "Misidentified-Particle Group",
    "错误识别（排除统计）": "Misidentified (excluded from statistics)",
    "剔除模式：连续左键点击错误识别颗粒；删除后立即从显示和全部统计中排除。按住 Option/Alt 可调整已有选框。":
        "Exclusion mode: left-click misidentified particles to remove them immediately from the display and all statistics. Hold Option/Alt to adjust existing boxes.",
    "开启后可在下图连续左键剔除误识别颗粒":
        "Enable this option to remove misidentified particles by left-clicking them in the image below.",
    "添加模式：黄色框只是搜索范围；点击后使用处理图完整识别中心和旋转。按住 Option/Alt 时只编辑已有框。":
        "Add mode: the yellow box defines only the search region. Clicking detects the complete particle in the processed image and estimates its center and orientation. Hold Option/Alt to edit existing boxes only.",
    "自由选框模式：在下方处理图拖画选框；松开后自动识别框内颗粒并校正中心和旋转。":
        "Free-box mode: draw a box in the processed image. On release, the particle inside is detected and its center and orientation are refined automatically.",
    "规则形状用拖框创建；不规则外框用鼠标沿颗粒边缘自由勾勒":
        "Create regular shapes by dragging a box; trace irregular outlines freehand along the particle boundary.",
    "线测量添加模式：线段可拖动或依次点两点；任意线可按住绘制或连续点击，Esc 完成连续折线。":
        "Line-measurement mode: drag a segment or click its two endpoints. For a free path, drag continuously or click successive vertices; press Esc to finish.",
    "线型对象固定按中心路径长度统计": "Line objects are quantified by centerline path length",
    "分析JSON...": "Reading analysis JSON…",

    # FFT, TEM, and exported scientific reports.
    "FFT一阶峰不能稳定分成两套Square点阵。":
        "The first-order FFT peaks could not be reliably separated into two Square lattices.",
    "FFT一阶峰方向聚类失败。": "Clustering of first-order FFT peak orientations failed.",
    "FFT一阶峰拟合Twist": "Twist Fitted from First-Order FFT Peaks",
    "FFT两套Square点阵夹角不可靠。":
        "The relative angle between the two Square lattices in the FFT is unreliable.",
    "FFT中没有一阶点阵峰。": "No first-order lattice peaks were detected in the FFT.",
    "FFT候选峰不能组成两套正交Square一阶峰。":
        "The candidate FFT peaks do not form two orthogonal sets of first-order Square-lattice peaks.",
    "FFT点阵方向无法稳定分成两组。":
        "The FFT lattice orientations could not be reliably separated into two groups.",
    "FFT预测period超出TEM视野。": "The FFT-predicted Moiré period exceeds the TEM field of view.",
    "没有FFT预测period。": "No FFT-predicted Moiré period is available.",
    "TEM period计算Twist": "Twist Calculated from the TEM Moiré Period",
    "TEM实空间没有识别到足够的moiré单元。":
        "Too few Moiré unit cells were detected in the real-space TEM image.",
    "TEM有效FFT分块不足。": "Too few valid TEM tiles were available for local FFT analysis.",
    "TEM的FFT没有稳定识别两套一阶Square峰。":
        "The TEM FFT did not reliably resolve two sets of first-order Square-lattice peaks.",
    "TEM频谱没有稳定的Moiré周期峰。":
        "No stable Moiré-period peak was detected in the TEM spectrum.",
    "自动分析结果（所有图像保持原始宽高比）":
        "Automatic Analysis Results (all images retain their original aspect ratios)",
    "所有样品使用真实FFT和相位保持的Selected-spot孔径，并自动生成原始分辨率PNG、保持比例的SVG、JSON和TEM/FFT分列的CSV统计。":
        "All samples use the measured FFT and phase-preserving selected-spot apertures. Exports include full-resolution PNG, aspect-preserving SVG, JSON, and CSV statistics with separate TEM and FFT columns.",
    "已导出%d个样品。\n\n%s\n\nCSV包含每个样品TEM与FFT各自的a、twist和period；所有PNG保持原始像素，SVG保持相同比例并保留可编辑矢量标注。":
        "Exported %d samples.\n\n%s\n\nThe CSV contains separate TEM- and FFT-derived a, twist, and period values for each sample. PNG files retain the original pixel dimensions; SVG files retain the aspect ratio and editable vector annotations.",
    "手动导出当前分析结果（PNG、SVG、JSON）":
        "Export Current Analysis Results (PNG, SVG, JSON)",
    "已按当前Overlay比例手动导出PNG、SVG和JSON。\n\n%s\n\n分析阶段没有自动写出任何SVG。":
        "Exported PNG, SVG, and JSON using the current overlay scale.\n\n%s\n\nNo SVG files were written automatically during analysis.",
    "当前分析结果已导出": "Current Analysis Results Exported",
    "统计摘要": "Statistical Summary",
    "置信度": "Confidence",
    "角度误差": "Angular error",
    "无数据": "No data",
    "未定义": "Not defined",
    "任务进行中": "Task in Progress",
    "1. %s\n2. 选择单层或双层分析\n3. 开始 Bulk 分析，逐图识别 TEM / FFT\n4. 自动导出每图结果与 CSV 汇总":
        "1. %s\n2. Select single-layer or bilayer analysis\n3. Run Bulk Analysis to process each TEM image and its FFT\n4. Export per-image results and a summary CSV automatically",
    "4. 开始自动识别与统计（可选）": "4. Run Automatic Detection and Quantification (Optional)",
    "JSON 无法读取": "Cannot Read JSON",
    "Moiré 工程无法读取": "Cannot Read Moiré Project",
    "无法读取": "Cannot Read File",
    "无法载入": "Cannot Load File",
    "无法读取所选 TXT 文件。": "The selected TXT file could not be read.",
    "Seed helix–helix间距固定为2.8 nm。Kagome的a为cryo-EM 5.4 nm、干燥TEM 4.4 nm；Square–Kagome只预测Twist，不定义period。":
        "The Seed helix-center spacing is fixed at 2.8 nm. For Kagome, a = 5.4 nm in cryo-EM and 4.4 nm in dried TEM. Square–Kagome predicts twist but does not define a single Moiré period.",
    "● 普通调整模式 · 左键拖动/旋转选框":
        "● Standard editing mode · Left-drag to move or rotate a selection box",
    "⚠ 等待手动标尺长度": "⚠ Waiting for a manual scale-bar length",
    "✓ 元数据尺度识别成功：%.5g nm/px": "✓ Pixel scale read from metadata: %.5g nm/px",
    "✓ 坐标轴尺度识别成功：%.5g nm/px": "✓ Pixel scale read from the axis: %.5g nm/px",
    "仅改变尺寸统计；颗粒数量、Group归属和材料密度仍按实际识别区域计算":
        "Changes only the reported dimensions. Particle count, Group assignment, and material-area density remain based on the detected particle regions.",
    "先输入这批图共用的 scale bar 标注值，再选择图片。横线像素长度仍逐图自动识别。":
        "Enter the scale-bar value shared by the batch, then select the images. The scale-bar length in pixels is detected separately for each image.",
    "各图数值不同（文件名开头标注）": "Values differ by image (encoded at the start of each file name)",
    "图中 scale bar + OCR": "Image Scale Bar + OCR",
    "尺寸：选框几何尺寸": "Dimensions: Selection-Box Geometry",
    "尺度：上传后自动读取原始元数据或图中 scale bar":
        "Scale: read automatically from image metadata or an embedded scale bar after upload",
    "左键拖动：新建任意选框\n右键选框：设为 Lane 1/2…；Target 与 Staples 均为可选\n选中选框：拖动位置；拖动白色控制点改变四边界\n靠近其他边界：显示粉色对齐线并自动吸附":
        "Left-drag: create a selection box\nRight-click a box: assign Lane 1/2…; Target and Staples are optional\nSelected box: drag to move; drag a white handle to adjust an edge\nNear another edge: a magenta alignment guide appears and the edge snaps into alignment",
    "应用论文参数预设": "Apply Publication Parameter Preset",
    "所选文件是 cadnano/普通 JSON，不含 DNA Moiré Designer 的参数与工作流。\n此入口只接受本软件保存的 .moire.json 工程。":
        "The selected file is a cadnano or generic JSON file and does not contain DNA Moiré Designer parameters or workflow state.\nThis command accepts only .moire.json projects saved by DNA Moiré Designer.",
    "手动画出的标尺长度约为 %.1f px。\n请输入它代表的实际长度：":
        "The manually drawn scale-bar line is approximately %.1f px long.\nEnter the represented physical length:",
    "数量：<b>%d</b>　总长度：<b>%d nt</b><br>%s":
        "Count: <b>%d</b>　Total length: <b>%d nt</b><br>%s",
    "等待生成。": "Waiting for generation.",
    "设计预览无法读取。": "The design preview could not be read.",
    "识别状态": "Detection Status",
    "目标 %d · 副产物 %d · 排除区域 %d · 手动标尺 %d":
        "Targets: %d · Byproducts: %d · Exclusion regions: %d · Manual scale bars: %d",
    "有效识别 %d · 错误识别（已排除） %d · 目标 %d · 副产物 %d · 聚集体 %d · 未确定 %d · 单颗粒数量产率 %s · %s面积密度 %.2f%% · 目标密度 %.2f%% · 聚集体密度 %.2f%%\n%s尺寸（mean ± STDEV.S，n=%d）：%s %s nm；%s %s nm":
        "Valid detections: %d · Misidentified (excluded): %d · Targets: %d · Byproducts: %d · Aggregates: %d · Unclassified: %d · Single-particle number yield: %s · %s area density: %.2f%% · Target density: %.2f%% · Aggregate density: %.2f%%\n%s dimensions (mean ± STDEV.S, n = %d): %s %s nm; %s %s nm",
    "该组现已从数量产率、尺寸、密度和统计导出中排除。":
        "This Group is now excluded from number yield, dimensions, density, and statistical exports.",
    "请先输入这批图共用的 scale bar 标注值（nm），再选择图片。":
        "Enter the scale-bar value shared by this batch (nm), then select the images.",
    "请在右侧第一张 Original image 上，沿 scale bar 横线拖一个窄框，然后再次点击“开始自动识别与统计”。软件届时只会询问该标尺代表多少 nm。":
        "On the first Original Image at right, draw a narrow box along the scale-bar line, then click Run Automatic Detection and Quantification again. You will be asked only for the length represented by that scale bar in nm.",
    "这个手动 Group 还没有预设选框。\n\n选择“是”：请在上方图像框一个代表颗粒，之后使用该框的形状和尺寸连续添加。\n选择“否”：每次在下方处理图拖画不同大小的选框。":
        "This manual Group has no reference box.\n\nChoose Yes to draw a box around a representative particle in the upper image and reuse its shape and dimensions for subsequent additions.\nChoose No to draw a new box of any size in the processed image for each particle.",
    "选框调整未保存：框内重新识别失败（%s）。":
        "Selection-box adjustment was not saved because detection within the box failed (%s).",
    "颗粒 %s 已移动到 Group %s – %s；选框和实际掩膜保持不变。":
        "Particle %s was moved to Group %s – %s; its selection box and segmentation mask are unchanged.",
    "颗粒归属继承所在 Group；如需更改，请在图中右键移动到其他 Group。选框的位置和角度可在图中直接左键调整。错误识别 Group 仍可查看，但不参与任何统计。尺寸、面积和密度均由图中的半透明实际分割掩膜计算；外框仅用于提示形状，不参与数值统计。":
        "Each particle inherits its Group assignment. To change it, right-click the particle in the image and move it to another Group. Left-drag to adjust the selection-box position or angle. The Misidentified Group remains visible but is excluded from all statistics. Dimensions, area, and density are calculated from the translucent segmentation mask; the outline indicates shape only and is not used for quantification.",
})

# Publication-facing terminology used by the orthogonal-sequence analysis
# workbook.  The generator deliberately keeps stable Chinese source tokens;
# localization happens only when the report is written for the user.
_ENGLISH_REVIEWED.update({
    "来源": "Source",
    "名称": "Name",
    "序列（5′→3′）": "Sequence (5′→3′)",
    "长度": "Length",
    "GC（%）": "GC content (%)",
    "互补序列（5′→3′）": "Reverse complement (5′→3′)",
    "熔解温度（°C）": "Melting temperature (°C)",
    "分子量（g/mol）": "Molecular weight (g mol⁻¹)",
    "消光系数 ε260（L·mol⁻¹·cm⁻¹）":
        "Extinction coefficient ε260 (L mol⁻¹ cm⁻¹)",
    "最长连续相同碱基": "Longest homopolymeric run",
    "输入/新生成链间最差同向相同片段（nt）":
        "Worst same-orientation identity among input/generated sequences (nt)",
    "输入/新生成链间最差链间互补片段（nt）":
        "Worst interstrand complementarity among input/generated sequences (nt)",
    "与骨架链最差同向相同片段（nt）":
        "Worst same-orientation identity with scaffold (nt)",
    "与骨架链最差链间互补片段（nt）":
        "Worst interstrand complementarity with scaffold (nt)",
    "自身互补长度（nt）": "Self-complementarity length (nt)",
    "发卡茎长度（bp）": "Hairpin stem length (bp)",
    "最小汉明距离（%）": "Minimum Hamming distance (%)",
    "序列熵（bits）": "Sequence entropy (bits)",
    "状态": "Status",
    "问题说明": "Review notes",
    "来源1": "Source 1",
    "序列1": "Sequence 1",
    "来源2": "Source 2",
    "序列2": "Sequence 2",
    "同向相同片段（nt）": "Same-orientation identity (nt)",
    "链间互补片段（nt）": "Interstrand complementarity (nt)",
    "汉明距离": "Hamming distance",
    "汉明距离（%）": "Hamming distance (%)",
    "参数": "Parameter",
    "数值": "Value",
    "序列分析": "Sequence Analysis",
    "两两分析": "Pairwise Analysis",
    "设置": "Settings",
    "输入": "Input",
    "输入-%03d": "Input-%03d",
    "骨架链": "Scaffold",
    "新生成": "Newly generated",
    "新序列-%03d": "New sequence-%03d",
    "通过": "Pass",
    "需复核": "Review required",
    "输入/新生成链间：%s": "Input/generated-sequence pair: %s",
    "与骨架链：%s": "With scaffold: %s",
    "筛选模式": "Screening mode",
    "基础规则＋可选高级规则": "Core rules + optional advanced rules",
    "请求生成数量": "Requested sequence count",
    "实际生成数量": "Generated sequence count",
    "输入序列数量": "Input sequence count",
    "输入文件": "Input file",
    "无": "None",
    "骨架链数量": "Scaffold count",
    "骨架链名称": "Scaffold names",
    "熔解温度模型": "Melting-temperature model",
    "SantaLucia DNA最近邻模型": "SantaLucia nearest-neighbor DNA model",
    "熔解温度Na⁺浓度": "Na⁺ concentration for Tm calculation",
    "熔解温度Mg²⁺浓度": "Mg²⁺ concentration for Tm calculation",
    "熔解温度链浓度": "Oligonucleotide concentration for Tm calculation",
    "熔解温度互补目标链浓度":
        "Complementary-strand concentration for Tm calculation",
    "是否完整生成": "Generation complete",
    "是": "Yes",
    "否": "No",
    "随机方式": "Randomization method",
    "每次运行使用新的系统随机数":
        "New system-generated random seed for each run",
    "候选评价次数": "Candidates evaluated",
    "新序列长度": "New-sequence length",
    "全局GC下限": "Minimum global GC content",
    "全局GC上限": "Maximum global GC content",
    "最大连续相同碱基": "Maximum homopolymeric-run length",
    "最大同向相同片段": "Maximum same-orientation identity",
    "最大链间互补片段": "Maximum interstrand complementarity",
    "与骨架最大同向相同片段":
        "Maximum same-orientation identity with scaffold",
    "与骨架最大链间互补片段":
        "Maximum interstrand complementarity with scaffold",
    "启用局部GC规则": "Local-GC rule",
    "局部GC窗口长度": "Local-GC window length",
    "局部GC下限": "Minimum local GC content",
    "局部GC上限": "Maximum local GC content",
    "启用序列熵规则": "Sequence-entropy rule",
    "最低序列熵": "Minimum sequence entropy",
    "启用自身互补规则": "Self-complementarity rule",
    "最大自身互补长度": "Maximum self-complementarity length",
    "启用发卡规则": "Hairpin rule",
    "最大发卡茎长度": "Maximum hairpin-stem length",
    "启用汉明距离规则": "Hamming-distance rule",
    "最小汉明距离比例": "Minimum Hamming-distance fraction",
    "启用禁用片段规则": "Forbidden-motif rule",
    "禁用片段": "Forbidden motifs",
    "每轮候选池大小": "Candidate-pool size per round",
    "每条序列最大尝试次数": "Maximum attempts per sequence",
    "启用": "Enabled",
    "未启用": "Disabled",
    "全局GC超出范围": "global GC content outside the allowed range",
    "局部GC超出范围": "local GC content outside the allowed range",
    "连续相同碱基过长": "homopolymeric run exceeds the limit",
    "序列复杂度过低": "sequence complexity below the limit",
    "自身互补过长": "self-complementarity exceeds the limit",
    "发卡茎过长": "hairpin stem exceeds the limit",
    "包含禁用片段": "forbidden motif present",
    "同向相同片段过长": "same-orientation identity exceeds the limit",
    "链间互补片段过长": "interstrand complementarity exceeds the limit",
    "汉明距离不足": "Hamming distance below the limit",
    "因“%s”淘汰的候选数": "Candidates rejected: %s",
})
for _code, _items in _BUILTIN.items():
    _catalogs.setdefault(_code, {}).update({
        key: value for key, value in _items.items()
        if key not in _catalogs.get(_code, {})})
_ENGLISH_REVIEWED.update({'接受当前 Added Scaffold': 'Accept assigned scaffold sequences', '完成序列导入并导出 Final Export': 'Export final package', 'SST superlattice 第二层 / Seed Z3': 'SST sublattice 2nd layer / seed Z3', 'SST superlattice 第一层 / Seed Z1': 'SST sublattice 1st layer / seed Z1', 'Optional · Expert mode': 'Optional · expert mode', '收起 Optional Expert mode': 'Close optional expert mode', '专家模式': 'Expert mode', '收起专家模式': 'Close expert mode', '3.1  Add Scaffold 序列': '3.1 Assign scaffold sequences', '3.2  Add SST superlattice Input 序列': '3.2 Assign SST sublattice input sequences', '自动识别 SST superlattice Input 位置和长度': 'Detect SST sublattice input positions and lengths', '自动设计并 Add Input 序列': 'Design and assign SST sublattice inputs automatically', 'Import Input Sequences': 'Import and assign input sequences', '接受当前 Added SST superlattice Input': 'Accept assigned SST sublattice inputs', '✓ Added Scaffold 已接受': '✓ Scaffold sequences accepted', '接受 Added Scaffold': 'Accept assigned scaffold sequences', 'SST superlattice Input 已完成': 'SST sublattice inputs assigned', '自动设计并 Add SST Input 序列': 'Design and assign SST sublattice inputs automatically', '✓ Added SST superlattice Input 已接受': '✓ SST sublattice inputs accepted', '接受 Added SST superlattice Input': 'Accept assigned SST sublattice inputs', 'Final Export 完成：%s\n全部序列结果直接保存在 All Sequences 文件夹；PDB CIF oxView Files 文件夹包含单独的 Seed、SST superlattice Layer 1、Layer 2 全原子 PDB/mmCIF、oxDNA TOP/DAT 及纯柱 BILD 模型。\n导出已完成；下一步请打开目标文件夹检查交付文件。': 'Final export complete: %s\nAll sequence results are stored directly in the Oligonucleotide sequences folder. The PDB/oxView files folder contains separate all-atom PDB/mmCIF, oxDNA TOP/DAT, and BILD models for the seed and the 1st and 2nd SST sublattice layers.\nExport is complete. Open the destination folder and review the deliverables.', 'SST superlattice Layer 1 / Seed 第一支撑区': 'SST sublattice 1st layer / seed 1st support', 'SST superlattice Layer 2 / Seed 第二支撑区': 'SST sublattice 2nd layer / seed 2nd support', 'SST superlattice Layer 1 /\nSeed 第一支撑区': 'SST sublattice 1st layer /\nseed 1st support', 'SST superlattice Layer 2 /\nSeed 第二支撑区': 'SST sublattice 2nd layer /\nseed 2nd support', 'Add 内置 Scaffold 序列': 'Assign a built-in scaffold sequence', 'Add Scaffold 序列': 'Assign scaffold sequence', 'Import SST Input 序列': 'Import and assign input sequences', '识别到 %d 条 scaffold；请在上方逐条 Add scaffold。': 'Detected %d scaffold routes. Assign a scaffold sequence to each route above.', '已接受当前 Added SST superlattice Input；两层真实序列均已写入结构。': 'The assigned SST sublattice inputs were accepted, and the nucleotide sequences were written to both layers.', 'SST superlattice 1st layer和2nd layer至少需要64 bp。': 'The 1st and 2nd SST sublattice layers must each be at least 64 bp.', 'SST superlattice长度不合法': 'Invalid SST sublattice length', 'SST superlattice 1st layer、spacing和2nd layer必须是8 bp整数倍。': 'The 1st-layer length, spacing, and 2nd-layer length of the SST sublattice must each be a multiple of 8 bp.', '已识别 SST superlattice Input；可自动设计并 Add，或打开专家模式，依次 Export template、正交序列设计、Import sequences，最后在右侧接受当前 Input。': 'SST sublattice inputs detected. Design and assign them automatically, or open expert mode to export a template, design orthogonal sequences, and import and assign the sequences. Then accept the assigned inputs.', '已识别 SST superlattice Input。可自动设计并 Add，或使用专家模式依次导出模板、设计正交序列并导入。完成后在下方接受当前 Input。': 'SST sublattice inputs detected. Design and assign them automatically, or use expert mode to export a template, design orthogonal sequences, and import and assign the sequences. Then accept the assigned inputs below.', '无法识别：可靠Moiré单元少于2个': 'Not detected: fewer than two reliable moiré units', '无法识别：两组同阶峰太近或证据不足': 'Not detected: the two same-order peak sets are too close or the evidence is insufficient', '仅显示晶格常数 a': 'Only lattice constant a is shown', '当前已开放 Square–Square、Kagome–Kagome 和 Square–Kagome S8–R4×4C 的后续设计。': 'Downstream design is available for Square–square, Kagome–kagome, and Square–kagome S8–R4×4C.', '没有可接受的最终设计': 'No final design is available for acceptance', '当前最终文件和选择过的文件中，没有通过验证的 SST + Scaffold + Staple + Capture JSON。': 'Neither the current final file nor the selected files contain a validated SST sublattice + scaffold + staple + capture JSON.', '最终结构设计已接受：%s。\nSST 与 SST + Scaffold 是过程导出文件，不作为接受版本。第3步将只使用此最终JSON进行SST superlattice/capture序列设计。': 'Final structure design accepted: %s.\nThe SST sublattice-only and SST sublattice + scaffold files are intermediate exports and are not accepted versions. Step 3 uses only this final JSON for SST sublattice and capture-sequence design.', '3.3 Final Export': '3.3 Final export', 'FFT同阶峰不能可靠分成两套同对称性点阵。': 'The same-order FFT peaks cannot be separated reliably into two lattices with the same symmetry.', 'Kagome SST 1st layer和2nd layer至少需要64 bp。': 'The first and second Kagome SST layer lengths must each be at least 64 bp.', 'Kagome SST模板缺少helix：%s。': 'The Kagome SST template is missing helices: %s.', 'TEM实空间没有形成可独立验证的重复Moiré单元。': 'The TEM real-space image does not contain repeating moiré units that can be validated independently.', 'TEM视野包含的可靠Moiré单元少于2个，无法识别TEM twist。': 'The TEM field of view contains fewer than two reliable moiré units, so the TEM-derived twist cannot be determined.', 'Z2安全位点容量不足，无法完成截面平衡整数化。': 'The safe-site capacity in Z2 is insufficient to complete balanced integer allocation across the cross-section.', '检测到混合多层，TEM周期无法唯一归属于某一层对；仅报告FFT中可靠同对称性层对的twist。': 'Mixed multilayers were detected. The TEM-derived period cannot be assigned uniquely to one layer pair; only the twist of a reliable same-symmetry layer pair in the FFT is reported.', '混合多层图像的TEM周期无法唯一归属于某一层对。': 'The TEM-derived period of a mixed-multilayer image cannot be assigned uniquely to one layer pair.'})
_ENGLISH_REVIEWED.update({
    "Optional · Expert mode": "Optional: expert mode",
    "接受当前 Added SST superlattice Input":
        "Accept assigned SST sublattice input sequences",
    "接受 Added SST superlattice Input":
        "Accept assigned SST sublattice input sequences",
    "✓ Added SST superlattice Input 已接受":
        "✓ SST sublattice input sequences accepted",
})
_catalogs.setdefault("en", {}).update(_ENGLISH_REVIEWED)
_catalogs.setdefault("zh_CN", {})

# Repair recurring scientific terms inside longer automatically translated
# sentences.  Exact reviewed entries above always take precedence.
_ENGLISH_TERM_REPLACEMENTS = (
    ("CADNANO", "caDNAno"),
    ("Cadnano", "caDNAno"),
    ("cadnano", "caDNAno"),
    ("Card issuance", "Hairpin"),
    ("card issuance", "hairpin"),
    ("Skeleton chains", "Scaffolds"),
    ("skeleton chains", "scaffolds"),
    ("Skeleton chain", "Scaffold"),
    ("skeleton chain", "scaffold"),
    ("backbone chains", "scaffolds"),
    ("backbone chain", "scaffold"),
    ("backbone strands", "scaffolds"),
    ("backbone strand", "scaffold"),
    ("inter-chain complementation", "interstrand complementarity"),
    ("Inter-chain complementation", "Interstrand complementarity"),
    ("inter-chain complementary fragment", "interstrand complementary segment"),
    ("inter-strand complementary fragment", "interstrand complementary segment"),
    ("complementary fragment", "complementary segment"),
    ("same fragment in the same direction", "identical same-orientation segment"),
    ("same segment in the same direction", "identical same-orientation segment"),
    ("identical fragments in the same direction", "identical same-orientation segment"),
    ("Disabled motif", "Forbidden motif"),
    ("disabled fragments", "forbidden sequences"),
    ("disable fragment", "forbidden-sequence"),
)
for _source, _target in tuple(_catalogs.get("en", {}).items()):
    for _old, _new in _ENGLISH_TERM_REPLACEMENTS:
        _target = _target.replace(_old, _new)
    # Product terminology uses "SST sublattice" consistently.  Keep
    # scientific identifiers such as SST-a unchanged while normalizing UI,
    # report, dialog, and export prose that historically used bare SST or
    # SST superlattice.
    _target = _target.replace("SST Superlattice", "SST sublattice")
    _target = _target.replace("SST superlattice", "SST sublattice")
    _target = re.sub(
        r"\bSST\b(?![\w-]|\s+sublattice\b)", "SST sublattice", _target)
    _catalogs["en"][_source] = _target

_PLACEHOLDER = re.compile(
    r"%(?:\([^)]+\))?[-+#0-9.]*[diouxXeEfFgGcrsa]|\{[^{}]*\}")
_template_cache = {}
_language_runtime_cache = {}
_translation_result_cache = {}
_phrase_runtime_cache = {}


def current_language():
    return _language


def set_language(language):
    global _language
    # The distributable product is English-only.  In particular, an older
    # project carrying interface_language=zh_CN must not switch the installed
    # application, worker reports, or exports back to Chinese.
    language = "en"
    _language = language
    os.environ["MOIRE_INTERFACE_LANGUAGE"] = language
    return language


def _runtime_catalog(language):
    prepared = _language_runtime_cache.get(language)
    if prepared is not None:
        return prepared
    catalog = _catalogs.get(language, {})
    reverse = {}
    for items in _catalogs.values():
        for source, translated in items.items():
            # Identity entries from the Chinese source catalog must not make
            # an uncataloged rich-text block return early before its visible
            # child phrases are translated.
            if translated != source:
                reverse.setdefault(translated, source)
    templates = [
        (source, target, _template_pattern(source))
        for source, target in catalog.items()
        if _PLACEHOLDER.search(source) and source != target]
    # A catalog may contain generic entries such as ``%s``.  Match the most
    # specific sentence first so a no-context placeholder cannot swallow a
    # complete dynamic result/status message before its reviewed template.
    templates.sort(
        key=lambda item: len(_PLACEHOLDER.sub("", item[0])), reverse=True)
    prepared = catalog, reverse, templates
    _language_runtime_cache[language] = prepared
    return prepared


def _template_pattern(source):
    pieces, end = [], 0
    for match in _PLACEHOLDER.finditer(source):
        # Python %-format templates spell one literal percent as ``%%``;
        # widgets receive the already-formatted text containing only ``%``.
        pieces.append(re.escape(
            source[end:match.start()].replace("%%", "%")))
        pieces.append("(.+?)")
        end = match.end()
    pieces.append(re.escape(source[end:].replace("%%", "%")))
    return re.compile("^" + "".join(pieces) + "$", re.DOTALL)


def _translated_template(source, target, text):
    key = (source, text)
    pattern = _template_cache.get(key)
    if pattern is None:
        pattern = _template_pattern(source)
        _template_cache[key] = pattern
    match = pattern.match(text)
    if match is None:
        return None
    values = iter(match.groups())
    return _PLACEHOLDER.sub(
        lambda unused: translate(next(values, unused.group())),
        target).replace("%%", "%")


def _translate_catalog_fragments(value, language, catalog):
    """Translate reviewed phrases embedded in a larger legacy message.

    Some historical UI and worker messages were assembled from adjacent
    string literals.  Their complete value therefore does not always match a
    generated catalog key even though all of the visible phrases are already
    reviewed.  Longest-first replacement closes that presentation-boundary
    gap without altering scientific identifiers or numeric values.
    """
    phrases = _phrase_runtime_cache.get(language)
    if phrases is None:
        phrases = tuple(sorted(
            ((source, target) for source, target in catalog.items()
             if source != target and re.search(r"[\u3400-\u9fff]", source)
             and not re.search(r"[\u3400-\u9fff]", target)),
            key=lambda item: len(item[0]), reverse=True))
        _phrase_runtime_cache[language] = phrases
    translated = value
    for source, target in phrases:
        if source in translated:
            translated = translated.replace(source, target)
    return translated


def translate(text, language=None):
    """Translate UI prose while preserving values substituted at runtime."""
    if text is None:
        return text
    language = language or _language
    value = str(text)
    cache_key = language, value
    if cache_key in _translation_result_cache:
        return _translation_result_cache[cache_key]
    catalog, reverse, templates = _runtime_catalog(language)
    if value in catalog:
        result = catalog[value]
        _translation_result_cache[cache_key] = result
        return result
    leading = value[:len(value)-len(value.lstrip())]
    trailing = value[len(value.rstrip()):]
    core = value.strip()
    if core in catalog:
        result = leading + catalog[core] + trailing
        _translation_result_cache[cache_key] = result
        return result
    # Resolve strings that were already translated before a language switch.
    source = reverse.get(core)
    if source is not None:
        result = leading + catalog.get(source, source) + trailing
        _translation_result_cache[cache_key] = result
        return result
    # Match formatted status messages against source templates.
    for source, target, pattern in templates:
        match = pattern.match(core)
        if match is None:
            continue
        values = iter(match.groups())
        translated = _PLACEHOLDER.sub(
            lambda unused: translate(
                next(values, unused.group()), language), target).replace(
                    "%%", "%")
        result = leading + translated + trailing
        _translation_result_cache[cache_key] = result
        return result
    # Rich-text QLabel headings often combine several independently cataloged
    # phrases inside <b>/<span>/<br> markup.  Translate only visible text
    # nodes so colors, fonts, alignment and other HTML remain untouched.
    if "<" in core and ">" in core:
        def replace_markup_text(match):
            visible = html.unescape(match.group(2))
            localized = translate(visible, language)
            return (match.group(1) + html.escape(localized, quote=False) +
                    match.group(3))

        localized_markup = re.sub(
            r"(>)([^<>]+)(<)", replace_markup_text, core)
        if localized_markup != core:
            result = leading + localized_markup + trailing
            _translation_result_cache[cache_key] = result
            return result
    # Preserve line layout and translate any exact line fragments.
    if "\n" in core:
        localized_lines = "\n".join(
            catalog.get(line, line) for line in core.splitlines())
        localized_lines = _translate_catalog_fragments(
            localized_lines, language, catalog)
        if localized_lines != core:
            result = leading + localized_lines + trailing
            _translation_result_cache[cache_key] = result
            return result
    fragmented = _translate_catalog_fragments(core, language, catalog)
    if fragmented != core:
        result = leading + fragmented + trailing
        _translation_result_cache[cache_key] = result
        return result
    _translation_result_cache[cache_key] = value
    return value


t = translate


def localize_xlsx(filename, language=None):
    """Translate workbook prose without changing sequences or numbers."""
    language = "en"
    path = Path(filename)
    temporary_path = None
    try:
        with tempfile.NamedTemporaryFile(
                dir=str(path.parent), prefix=path.stem + "-localized-",
                suffix=path.suffix, delete=False) as temporary:
            temporary_path = Path(temporary.name)
        text_node_pattern = re.compile(
            r"(<t(?:\s[^>]*)?>)(.*?)(</t>)", re.DOTALL)
        sheet_name_pattern = re.compile(
            r'(<sheet\b[^>]*\bname=")([^"]*)(")', re.DOTALL)

        def localize_text_node(match):
            source = html.unescape(match.group(2))
            localized = translate(source, language)
            return (match.group(1) + html.escape(localized, quote=False) +
                    match.group(3))

        def localize_sheet_name(match):
            source = html.unescape(match.group(2))
            localized = translate(source, language)[:31]
            return (match.group(1) + html.escape(localized, quote=True) +
                    match.group(3))

        def ensure_column_width(xml, column, minimum):
            pattern = re.compile(
                r'(<col\b(?=[^>]*\bmin="%d")(?=[^>]*\bmax="%d")'
                r'[^>]*\bwidth=")([^"]+)(")' % (column, column))

            def widen(match):
                try:
                    width = max(float(match.group(2)), float(minimum))
                except ValueError:
                    width = float(minimum)
                return match.group(1) + ("%g" % width) + match.group(3)

            return pattern.sub(widen, xml)

        def widen_orthogonal_report(xml):
            widths = None
            if "<t>Screening mode</t>" in xml:
                widths = {1: 62, 2: 72}
            elif "<t>Source 1</t>" in xml:
                widths = {
                    1: 16, 2: 18, 3: 16, 4: 18, 5: 32,
                    6: 36, 7: 20, 8: 22, 9: 16, 10: 48,
                }
            elif ("<t>Source</t>" in xml and
                  "<t>Reverse complement (5′→3′)</t>" in xml):
                widths = {
                    1: 18, 2: 20, 3: 52, 4: 10, 5: 16, 6: 52,
                    7: 28, 8: 26, 9: 42, 10: 28, 11: 52, 12: 54,
                    13: 34, 14: 28, 15: 30, 16: 22, 17: 18, 18: 48,
                    19: 34, 20: 48,
                }
            for column, minimum in (widths or {}).items():
                xml = ensure_column_width(xml, column, minimum)
            return xml

        with ZipFile(path, "r") as source_archive, ZipFile(
                temporary_path, "w") as target_archive:
            for member in source_archive.infolist():
                payload = source_archive.read(member.filename)
                if member.filename.endswith(".xml"):
                    xml = payload.decode("utf-8")
                    xml = text_node_pattern.sub(localize_text_node, xml)
                    if member.filename == "xl/workbook.xml":
                        xml = sheet_name_pattern.sub(
                            localize_sheet_name, xml)
                    elif member.filename.startswith("xl/worksheets/"):
                        xml = widen_orthogonal_report(xml)
                    payload = xml.encode("utf-8")
                target_archive.writestr(member, payload)
        temporary_path.replace(path)
    except Exception as error:
        if temporary_path is not None:
            try:
                temporary_path.unlink()
            except OSError:
                pass
        warnings.warn(
            "Could not localize XLSX report %s: %s" % (path, error),
            RuntimeWarning)
    return filename


def localize_csv(filename, language=None):
    language = "en"
    path = Path(filename)
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.reader(handle))
        with path.open("w", encoding="utf-8-sig", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerows([
                [translate(cell, language) for cell in row] for row in rows])
    except Exception:
        pass
    return filename


def localize_svg(filename, language=None):
    language = "en"
    path = Path(filename)
    try:
        source = path.read_text(encoding="utf-8")

        def replace_text(match):
            value = html.unescape(match.group(2))
            localized = translate(value, language)
            if localized == value:
                return match.group(0)
            return match.group(1) + html.escape(
                localized, quote=False) + match.group(3)

        # Replace text-node contents only.  This keeps paths, namespaces,
        # styles, Illustrator metadata and editable text geometry byte-for-
        # byte unchanged outside the actual prose.
        source = re.sub(r"(>)([^<>]+)(<)", replace_text, source)
        path.write_text(source, encoding="utf-8")
    except Exception:
        pass
    return filename


def install_painter_hook():
    """Translate QPainter text, including SVG/PNG analysis annotations."""
    try:
        import sys
        import PyQt6.QtGui as qt_gui
        QPainter = qt_gui.QPainter
        if getattr(QPainter, "_moire_i18n_installed", False):
            return
        original = QPainter.drawText

        def draw_text(painter, *arguments):
            localized = list(arguments)
            for index in range(len(localized)-1, -1, -1):
                if isinstance(localized[index], str):
                    localized[index] = translate(localized[index])
                    break
            return original(painter, *localized)

        try:
            QPainter.drawText = draw_text
            QPainter._moire_i18n_installed = True
        except (AttributeError, TypeError):
            # PyQt builds commonly expose QPainter as an immutable SIP type.
            # A Python subclass is therefore installed into the already
            # imported Moiré modules and into QtGui for future imports.  Their
            # functions resolve QPainter through the module global at call
            # time, so vector/raster rendering remains otherwise unchanged.
            class LocalizedPainter(QPainter):
                _moire_i18n_installed = True

                def drawText(self, *arguments):
                    localized = list(arguments)
                    for index in range(len(localized)-1, -1, -1):
                        if isinstance(localized[index], str):
                            localized[index] = translate(localized[index])
                            break
                    return original(self, *localized)

            qt_gui.QPainter = LocalizedPainter
            for name, module in tuple(sys.modules.items()):
                if (name == "moire_designer" or
                        name.startswith("moire_designer.")) and \
                        getattr(module, "QPainter", None) is QPainter:
                    setattr(module, "QPainter", LocalizedPainter)
    except Exception:
        # Widgets still localize normally; explicit export finalizers remain
        # active even if a vendor-specific frozen Qt build rejects both forms.
        return


def install_dialog_hooks():
    """Localize modal-dialog arguments before a native dialog is created."""
    try:
        from PyQt6.QtWidgets import QFileDialog, QInputDialog, QMessageBox
        if getattr(QFileDialog, "_moire_i18n_installed", False):
            return

        def translate_filter(value):
            parts = []
            for item in str(value).split(";;"):
                match = re.match(r"^(.*?)(\s*\([^()]*\))$", item)
                if match:
                    parts.append(translate(match.group(1)) + match.group(2))
                else:
                    parts.append(translate(item))
            return ";;".join(parts)

        for name in ("information", "warning", "critical", "question"):
            original = getattr(QMessageBox, name)

            def message(parent, title, text, *args, _original=original):
                return _original(parent, translate(title), translate(text),
                                 *args)

            setattr(QMessageBox, name, staticmethod(message))

        for name in ("getOpenFileName", "getOpenFileNames",
                     "getSaveFileName"):
            original = getattr(QFileDialog, name)

            def file_dialog(parent=None, caption="", directory="",
                            file_filter="", *args, _original=original):
                return _original(parent, translate(caption), directory,
                                 translate_filter(file_filter), *args)

            setattr(QFileDialog, name, staticmethod(file_dialog))
        original_directory = QFileDialog.getExistingDirectory

        def directory_dialog(parent=None, caption="", directory="", *args):
            return original_directory(parent, translate(caption), directory,
                                      *args)

        QFileDialog.getExistingDirectory = staticmethod(directory_dialog)

        for name in ("getText", "getMultiLineText", "getItem", "getInt",
                     "getDouble"):
            original = getattr(QInputDialog, name)

            def input_dialog(parent, title, label, *args,
                             _original=original):
                return _original(parent, translate(title), translate(label),
                                 *args)

            setattr(QInputDialog, name, staticmethod(input_dialog))
        QFileDialog._moire_i18n_installed = True
    except Exception:
        return


class UiLocalizer:
    """Continuously localize new and dynamically updated Qt widgets."""

    def __init__(self, root, interval_ms=250):
        from PyQt6.QtCore import QTimer
        self.root = root
        self._item_sources = {}
        self._item_lasts = {}
        self.timer = QTimer(root)
        self.timer.setInterval(interval_ms)
        self.timer.timeout.connect(self.retranslate)
        self.timer.start()

    @staticmethod
    def _value(owner, key, getter, setter):
        current = getter()
        source_key = "_moire_i18n_source_" + key
        last_key = "_moire_i18n_last_" + key
        source = getattr(owner, source_key, None)
        last = getattr(owner, last_key, None)
        if source is None or current != last:
            source = current
            setattr(owner, source_key, source)
        localized = translate(source)
        if current != localized:
            setter(localized)
        setattr(owner, last_key, localized)

    def retranslate(self):
        try:
            from PyQt6.QtGui import QAction
            from PyQt6.QtWidgets import (
                QAbstractButton, QComboBox, QDockWidget, QGroupBox, QLabel,
                QLineEdit, QListWidget, QMainWindow, QTabWidget,
                QDoubleSpinBox, QSpinBox, QTableWidget, QTreeWidget,
                QGraphicsScene, QGraphicsSimpleTextItem, QGraphicsTextItem,
                QWidget)
            objects = [self.root] + self.root.findChildren(QWidget) + \
                self.root.findChildren(QAction)
            for item in objects:
                if getattr(item, "_moire_i18n_skip", False):
                    continue
                if isinstance(item, QWidget) and item.isWindow():
                    self._value(item, "windowTitle", item.windowTitle,
                                item.setWindowTitle)
                if isinstance(item, (QLabel, QAbstractButton)):
                    self._value(item, "text", item.text, item.setText)
                if isinstance(item, QGroupBox):
                    self._value(item, "title", item.title, item.setTitle)
                if isinstance(item, QDockWidget):
                    self._value(item, "windowTitle", item.windowTitle,
                                item.setWindowTitle)
                if isinstance(item, QLineEdit):
                    self._value(item, "placeholder", item.placeholderText,
                                item.setPlaceholderText)
                if isinstance(item, (QSpinBox, QDoubleSpinBox)):
                    self._value(
                        item, "specialValueText", item.specialValueText,
                        item.setSpecialValueText)
                if isinstance(item, QAction):
                    self._value(item, "text", item.text, item.setText)
                if hasattr(item, "toolTip") and hasattr(item, "setToolTip"):
                    self._value(item, "toolTip", item.toolTip,
                                item.setToolTip)
                if isinstance(item, QComboBox):
                    sources = getattr(item, "_moire_i18n_items", {})
                    lasts = getattr(item, "_moire_i18n_item_lasts", {})
                    for index in range(item.count()):
                        current = item.itemText(index)
                        if index not in sources or current != lasts.get(index):
                            sources[index] = current
                        localized = translate(sources[index])
                        if current != localized:
                            item.setItemText(index, localized)
                        lasts[index] = localized
                    item._moire_i18n_items = sources
                    item._moire_i18n_item_lasts = lasts
                if isinstance(item, QTabWidget):
                    sources = getattr(item, "_moire_i18n_tabs", {})
                    lasts = getattr(item, "_moire_i18n_tab_lasts", {})
                    for index in range(item.count()):
                        current = item.tabText(index)
                        if index not in sources or current != lasts.get(index):
                            sources[index] = current
                        localized = translate(sources[index])
                        if current != localized:
                            item.setTabText(index, localized)
                        lasts[index] = localized
                    item._moire_i18n_tabs = sources
                    item._moire_i18n_tab_lasts = lasts
                if isinstance(item, (QTableWidget, QTreeWidget)):
                    header = item.headerItem() if isinstance(
                        item, QTreeWidget) else None
                    if isinstance(item, QTableWidget):
                        for column in range(item.columnCount()):
                            cell = item.horizontalHeaderItem(column)
                            if cell is not None:
                                self._translate_item_text(cell, 0)
                        for row in range(item.rowCount()):
                            header_cell = item.verticalHeaderItem(row)
                            if header_cell is not None:
                                self._translate_item_text(header_cell, 0)
                            for column in range(item.columnCount()):
                                cell = item.item(row, column)
                                if cell is not None:
                                    self._translate_item_text(cell, 0)
                    elif header is not None:
                        for column in range(header.columnCount()):
                            self._translate_item_text(header, column)
                        for index in range(item.topLevelItemCount()):
                            self._translate_tree_item(
                                item.topLevelItem(index))
                elif isinstance(item, QListWidget):
                    for index in range(item.count()):
                        self._translate_item_text(item.item(index), 0)
            for scene in self.root.findChildren(QGraphicsScene):
                for graphics_item in scene.items():
                    if isinstance(graphics_item, QGraphicsSimpleTextItem):
                        current = graphics_item.text()
                        source = getattr(
                            graphics_item, "_moire_i18n_source", current)
                        last = getattr(
                            graphics_item, "_moire_i18n_last", None)
                        if current != last:
                            source = current
                        localized = translate(source)
                        if current != localized:
                            graphics_item.setText(localized)
                        graphics_item._moire_i18n_source = source
                        graphics_item._moire_i18n_last = localized
                    elif isinstance(graphics_item, QGraphicsTextItem):
                        current = graphics_item.toPlainText()
                        source = getattr(
                            graphics_item, "_moire_i18n_source", current)
                        last = getattr(
                            graphics_item, "_moire_i18n_last", None)
                        if current != last:
                            source = current
                        localized = translate(source)
                        if current != localized:
                            graphics_item.setPlainText(localized)
                        graphics_item._moire_i18n_source = source
                        graphics_item._moire_i18n_last = localized
        except RuntimeError:
            return

    def _translate_item_text(self, item, column=0):
        key = id(item), int(column)
        current = item.text(column) if hasattr(item, "columnCount") \
            else item.text()
        source = self._item_sources.get(key)
        last = self._item_lasts.get(key)
        if source is None or current != last:
            source = current
            self._item_sources[key] = source
        localized = translate(source)
        if current != localized:
            if hasattr(item, "columnCount"):
                item.setText(column, localized)
            else:
                item.setText(localized)
        self._item_lasts[key] = localized

    def _translate_tree_item(self, item):
        for column in range(item.columnCount()):
            self._translate_item_text(item, column)
        for index in range(item.childCount()):
            self._translate_tree_item(item.child(index))
