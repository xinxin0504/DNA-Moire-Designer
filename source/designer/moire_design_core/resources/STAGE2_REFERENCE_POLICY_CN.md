# DNA Moiré Designer 第二步范本与规则边界

## 范本用途

- `Square_Seed_2L_newtemplate.json`：唯一 Seed 基线。其 scaffold、普通 staple、capture 核心及颜色均为不可重建的权威数据；Square 和 Kagome SST 都使用这一 2L Seed。
- `Kagome_Seed_Ka-seed-pore_3L.json`：只用于提取 Kagome capture 延伸的物理连接拓扑和位点；不提供 Seed scaffold routing，不参与 scaffold 分子数、seam、边界或 staple routing 决策。
- `square_sst_*bp_fixture.json` 与 `kagome_sst_*bp_fixture.json`：已审核的 SST superlattice 冻结范本，开放 96/104/112/120/128 bp。不使用动态生成器重算它们。

## 当前开放范围

- Seed 截面仅开放 8×8 外形、4×4 pore。
- Seed 完整 routing 当前仅开放已审核的 128/32/128 2L 基线。
- Seed 的物理长度、scaffold、普通 staple、capture 核心、nick、seam 与颜色始终来自固定范本，不再接受 Z1/Z3 长度输入。
- SST 长度和 Z2 spacing 只改变 SST 的冻结范本选择、32-bp 相位平移、与固定 Seed 的实际重叠及 capture 列；不得选择任何 Seed crop/growth/edge-seam/AutoCS 分支。
- 8×7 / seed87 / 4×3 pore 分支、测试与例外规则均已移除。

## Seed scaffold

- 冻结基线含两条 scaffold：7300 nt 与 7336 nt，均不超过 7557 nt。
- Seed 不执行“小于8 bp”crossover间距检查；范本拓扑按原样接受。
- Moiré Designer 不允许缩进、增长、边缘seam、重新AutoCS、重新Autobreak、普通staple补洞/合并或capture核心长度修复。
- 唯一允许修改Seed内容的是用户Twist/period产生的indel：只能落在Z2的连续双链安全位点，避开所有nick和scaffold/staple crossover，均匀分布，并硬性保证两条scaffold各自不超过7557 nt。
- 较长SST导致左侧画布空间不足时，允许整体平移整个Seed坐标系；平移必须保持Seed内部全部拓扑、颜色和相对坐标不变。
- 固定范本的两条 scaffold 长度为 7300 nt 与 7336 nt；Z2 insertion 后仍分别硬性不超过 7557 nt。Moiré Designer 不再重新选择 1/2/3 条 scaffold，也不重新分 seam。
- Seed crossover、边缘错位与 scaffold-only 区域均按固定范本原样接受；不再用通用 AutoCS 或旧边缘规则重新判定或改写它们。

## Staple 与 capture

- Seed 普通 staple 与 capture 核心长度不再参与生成判定；它们按新范本原样保留，不重新 break、合并、补洞或按长度拒绝。
- caDNAno 自身的普通 staple 长度与 Autobreak 规则仍保留在 caDNAno 主程序中，但 Moiré Designer 固定 Seed 不再次执行它们。
- capture 候选 base 不应与原生 AutoCS staple crossover 冲突；出现冲突视为 SST 32-bp 相位/平移错误，不通过删除 AutoCS crossover 来遮掩。
- 新范本共144条非黑色Seed核心staple：64条物理Square capture、64条32-bp平移capture、16条Z2潜在capture。序列导出时全部归入`staple-capture`；只有实际连接SST的核心附加capture延伸。
- 每层至少 2 列 capture pair 是 SST spacing/重叠的硬性可行性检查；失败时要求减小 spacing，绝不通过改变 Seed 长度补偿。
- Kagome capture pair 使用 Kagome 实际活动 SST mask 和 3L capture 范本的物理拓扑；不复用 Square 的位点计数方式。

## 保护边界

- 上述逻辑仅属于 DNA Moiré Designer 第二步。
- caDNAno 主程序自带的 AutoCS scaffold、AutoCS staple、Autobreak 及其他功能不在本次清理范围内，未删除、未替换。
