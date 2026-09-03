# Kagome 8×8−4×4 pore Seed / SST Capture 规则

## 1. 参考文件与范围

- Kagome reference design: `Ka-seed-pore_3L.json`
- Square 对照：`S-seed-pore_3L.json`
- 本文总结 8×8 外框、4×4 pore Seed 与 Kagome SST 之间的 capture 结构和序列映射。
- 参考文件只画出 SST 单元位于 Seed capture face 左半区时的物理连接；右半区必须在序列导出时将 SST 单元平移后，按实际接触 helix 和 base 相位再次求解。

SHA-256：

```text
Kagome  5cf88ca67ea6c22db1544bd9da615832eba5957299f39e42903d5ff3669bad26
Square  be762eb2975f91a41d77b41d08cc9694ded4cbfc4e5046459e0e78b3c78d7d72
```

## 2. Seed 本体：Kagome 与 Square 相同

- Seed 使用 48 条 helix，编号 `0–47`，构成 8×8 外框减去中心 4×4 pore：
  - 外框范围：row `16–23`、col `18–25`；
  - pore：row `18–21`、col `20–23`；
  - 上下两行完整，左右各保留两列。
- Kagome 与 Square 的 `0–47` scaffold 记录逐 base 完全一致。
- 除了因 Kagome 空位而取消 capture 的 Seed staple `h0/h2` 外，其余 Seed staple 路由也逐 base 一致。
- 因此未来 Kagome Moiré 不需要重新发明 Seed scaffold/pore routing；可以复用 Square Seed 生成器。必须独立实现的是 Kagome SST 和 Seed–SST capture 接口。

## 3. Kagome SST 截面差异

Square SST 单元使用 `48–63` 全部 16 条 helix。Kagome 只使用：

```text
48,49,50,51,52,54,56,57,58,59,60,62
```

Kagome 空位：

```text
53,55,61,63
```

局部 4×4 排布：

```text
row24: 48 49 50 51
row25: 55 54 53 52   (55、53为空)
row26: 56 57 58 59
row27: 63 62 61 60   (63、61为空)
```

空位的影响有两类：

1. 原本连接到空位 `61/63` 的 capture 必须取消，Seed 端恢复为未延伸状态。
2. 原本以空位为 U 型搭档的活动 helix 变为线型 SST；其 capture 延伸长度必须从实际路由读取，不能统一假定为 16 nt。

本文件中的线型集合与前述 Kagome SST-only 模板一致，只是 helix 编号整体加 32：

- `scaf` 线型 helix：`52,54,60,62`
- `stap` 线型 helix：`48,50,56,58`

## 4. 参考文件中的三个重复区段

`Ka-seed-pore_3L.json` 中存在三个相同相位的 128-bp SST/capture 区段：

| 区段 | scaf 活动区 | stap 活动区 | capture positions |
|---|---|---|---|
| 1 | 48–175 | 40–183 | 56,72,88,104,120,136,152,168 |
| 2 | 208–335 | 200–343 | 216,232,248,264,280,296,312,328 |
| 3 | 368–495 | 360–503 | 376,392,408,424,440,456,472,488 |

- 每个区段 8 个 capture columns，间隔 16 bp。
- 相邻两个 16-bp column 构成一个完整 32-bp A/B 周期。
- 三个区段之间的规则相同；未来软件应按实际 SST/Seed 重叠区生成，不应把“必须有三个区段”写死。

### 理论候选与实际 capture 必须分开

- 上表的完整周期只描述 Kagome SST 自身的**理论合法候选**，不能由
  `128 bp` 或其他 SST 长度推导固定的实际 capture 数量。
- 物理 capture 必须逐层计算：固定 Seed 范本 capture 位点与当前平移后
  SST 双链区的真实交集，再检查该 base 的 Kagome 活动 mask、方向和路由。
- 所以实际第一列既可能是 4 个 U 型端点，也可能是 2 个线形端点；类型由
  实际重叠到的范本位点决定，不由 SST 左边界、长度或局部列编号决定。
- 没有与 Seed 双链区重叠的理论候选可在审阅文件中另色显示，但不得计入
  真实 bridge 数、capture pair 数、序列导出或最少列数检查。
- Square 同样遵守这一交集原则；Square 的 128-bp SST 也没有固定的实际
  capture 数量。

## 5. 已画出的左半 SST 单元：实际 Capture 映射

### Seed capture faces

- Face 1（上边）：Seed `h0–7`。
- Face 2（下边）：Seed `h24–31`。
- 当前物理结构只覆盖每个 face 的左四列：
  - Face 1 origin half：`h0,1,2,3`
  - Face 2 origin half：`h28,29,30,31`
- 尚未画出的右四列：
  - Face 1 translated half：`h4,5,6,7`
  - Face 2 translated half：`h24,25,26,27`

### 32-bp 周期

以第一个区段为例；后续每加 32 bp 重复：

| 相位 | base | Kagome 中实际存在的桥 | Square 额外存在、Kagome 因空位取消的桥 |
|---|---:|---|---|
| B0 | 56 (`mod32=24`) | `49→30`, `51→28`, `60→3`, `62→1` | 无 |
| A0 | 72 (`mod32=8`) | `48→31`, `50→29` | `61→2`, `63→0` |

- 箭头为 `SST helix → Seed helix`，两端使用同一 base 坐标。
- 所有 Seed–SST capture 都位于 `stap` 路由；`scaf` 没有 Seed–SST 跨组连接。
- 在**该范本完整重叠区**内，B0 每周期保留 4 条桥；A0 每周期只保留
  Face 2 的 2 条桥。因此范本中的完整 Kagome 32-bp 周期有 6 条桥，
  对照 Square 完整周期有 8 条桥。
- 上述 6/8 只是完整范本周期的理论/局部计数，不是任意 128-bp SST 的
  固定实际 capture 数量。实际数量始终由当前 Seed 范本位点与平移后
  SST 双链重叠区取交集后逐端点统计。
- 三层参考文件的三个完整区段共 12 个周期：Kagome 为 72 条，Square 为
  96 条；这只是该参考几何的审计结果，不得作为可变长度设计的固定目标。

重要区分：

- `h0/h2` 属于已经覆盖的 origin half，但因接触位置对应 Kagome 空位 `63/61` 而没有 capture。
- `h4–7` 与 `h24–27` 才是尚未覆盖、需要 SST 平移后第二次生成序列的另一半。

## 6. 线型 SST 造成的 Capture 长度差异

Square capture 在 SST 侧均为 16 nt。Kagome 不能使用这个固定值：

- `h49,h51,h60,h62` 等保留 U 型相位的物理桥，SST 侧为 16 nt。
- A0 相位的 `h48/h50` 是由 Kagome 空位形成的线型 SST：
  - 大多数内部位置沿实际线型路由延伸 32 nt；
  - 每个 128-bp 区段的末端位置受边界 nick 影响，为 16 nt。
- 参考文件中 `h48/h50` 各有 9 个 32-nt 延伸和 3 个 16-nt 延伸。
- 因此 capture extension 必须由“桥接 base 到该 SST 线型组件终点”的实际路径长度计算，不能仅按 helix 编号或固定 phase 返回 16。

与 Seed staple core 合并后，参考内部 capture 产品主要为：

- `32+16=48 nt`
- `40+16=56 nt`
- `32+32=64 nt`

这与现有长度策略一致：32-nt SST 延伸时 Seed core 必须控制为 32 nt，完整链达到 64 nt 硬上限；16-nt 延伸可使用 32–40 nt 的 Seed core。文件最右端少量 21-nt Seed core 属于路由终端边界，不应当作周期模板推广。

## 7. Kagome SST 原始连接位点的三种类型

这里的类型专指：以未加入 Seed capture 的 `kagome_resource_128.json` 为基准，Kagome SST 在待连接 base 上原本是什么拓扑。它不是 Seed 位点分类，也不能从加入 capture 后的文件反推。

| Kagome SST 原始类型 | capture 端点数 | 原始位置 | 正确改链动作 |
|---|---:|---|---|
| SST–SST crossover | 48 | `h49,h51,h60,h62`；base `56/88/120/152`、`216/248/280/312`、`376/408/440/472` | 原 crossover 的两个方向都要互反拆开，再将形成的两个 SST 端点分别接到各自 Seed 接触点 |
| 原有 nick | 18 | `h48,h50`；base `72/104/136`、`232/264/296`、`392/424/456` | 不再切链，直接把该 nick 上方向正确的自由端延伸到 Seed；另一侧原有连接保持不变 |
| 线型连续链 | 6 | `h48,h50`；base `168/328/488` | 先互反切开同一 helix 上的连续连接，形成新端点，再把方向正确的端点接到 Seed |

计数口径是“将接到 Seed 的 Kagome SST 端点”。其中 48 个 crossover 端点对应 24 对原始 SST–SST crossover；每拆开一对 crossover，会形成两个分别接向 Seed 的 capture 端点。

三类位点必须针对**不可变的原始 Kagome SST 路由**预先判断，然后再统一应用修改：

1. 先为全部候选点记录 `helix/base/strand/side(prev|next)/original_partner/type`。
2. 校验所有原始连接均为互反记录。
3. 对 crossover 成对拆除；对连续链成对切开；原有 nick 不做额外切割。
4. 最后建立 Seed–SST 互反连接，并重新遍历组件计算实际 capture 延伸长度。

不能边扫描边改：例如先处理相邻 capture 后，新产生的端点会让后面的线型连续位点看起来像“已有 nick”，从而改变切链位置与序列长度。translated-half 也必须在平移后的**未加 capture SST 基线**上重新判定这三类，不能直接复制 origin-half 的类型。

### Kagome capture pair 不能使用 Square 计数

Kagome 的一个完整 capture pair 仍覆盖一个 32-bp 重复周期，但其有效位点集合由当前保留的 **8×8 Seed截面 × Kagome SST活动mask × 相对平移/相位** 联合决定，而不是 Square 的固定 `4 bridges/column`：

- 对当前8×8−4×4 pore Seed接触：crossover column有4个端点，linear/nick column有2个端点，因此一个完整pair为`4+2=6`个物理端点；
- 同一周期中的另外2个接触落在Kagome空位，标记为`intrinsic_absent`，不能记录成被删除或不完整；
- translated half 同样由平移后的活动 mask 和相位重算，不能用 Square 的8端点目标补齐。

pair 是颜色、周期和统计分组，不是删除单位。若 Seed–SST 重叠边缘只容纳其中一列，必须保留该列全部实际合法端点；不能为了得到“完整pair”而把合法单列一并删除。

Capture 保留/删除优先级：

1. 根据当前8×8 Seed截面、Kagome SST活动mask、相对平移和base相位生成候选目录；排除实际命中的 `intrinsic_absent` 空位。
2. 保留所有同时满足实际接触、双链重叠和方向/相位合法的候选。
3. Seed routing 优先用合法10/11-bp端点微调、同侧≤21-bp硬错位范围和纵向 nick 求解来容纳边缘 capture。
4. 不得因沿用 Square 的列数、桥数、pair完整性或固定16-nt延伸假设而删除 Kagome capture。
5. 只有候选实际落在SST空位、已超出真实双链接触/重叠区，或穷尽上述合法路由仍违反不可放宽的结构硬约束时，才可把该**具体最外侧候选**标记为 `forced_omission`。
6. 每次删除必须输出原因、layer、pair、column、helix/base和尝试过的保留方案；不得静默删除，也不得连带删除同pair内仍合法的其他候选。

## 8. 平移后另一半 Capture 的生成规则

不能把 origin half 的序列按编号机械复制。正确流程是：

1. 将完整 Kagome SST 单元按合法点阵平移到 Seed face 的另一半；本参考截面对应从 col `18–21` 平移到 col `22–25`。
2. 同时应用该点阵平移要求的 base 相位偏移；如果平移导致 A/B 互换，必须以平移后的实际 base 相位为准。
3. 对每个 translated SST 活动 helix，使用几何接触关系寻找对应的 Seed face helix：
   - Face 1：origin `0–3` 映射到 translated `4–7`；
   - Face 2：origin `28–31` 映射到 translated `24–27`。
4. 查询平移后 SST 在该 base 的真实路由：
   - 接触落在 Kagome 空位时不生成 capture；
   - U 型 SST 取完整 U 单元对应的 16-nt 接触臂；
   - 线型 SST 沿实际组件追踪到终点，得到 16 或 32 nt。
5. 只在序列导出中增加 translated-half capture，不在结构 JSON 中同时画出两套互相重叠的 SST 单元。
6. 使用平移同源规则转移 SST 序列坐标，并重新验证总链长、方向、互补关系与 5′/3′ 端。

如果平移不改变 A/B phase，则由本参考可预期：

```text
B1: 49→26, 51→24, 60→7, 62→5
A1: 48→27, 50→25
```

`61→6`、`63→4` 仍落在 Kagome 空位，不能生成。但软件实现仍应走几何+实际相位求解，而不是把上表写死；不同 SST 平移方向或 base phase 可能改变接触集合。

## 9. 与 Square 可复用和必须分开的部分

### 可复用

- 8×8−4×4 pore Seed 截面、helix 编号与 scaffold routing。
- 两个 capture face 的定义。
- origin half / translated half 的工作流。
- 32-bp A/B 重复周期、同 base 桥接方式。
- capture 产品 64 nt 硬上限及 translated-half 仅在导出时加入的策略。

### Kagome 必须独立处理

- 12/16 活动 SST helix mask 与四个空位。
- U 型 SST 与线型 SST 的区分。
- 空位导致的 capture 缺失，不能为“每条 Seed face helix”强行补桥。
- 线型 SST 的 16/32-nt 实际 extension 长度。
- SST 点阵平移后的真实空间接触与 A/B phase 求解。
- 生成后按实际组件重新核验 capture 数量、序列长度及方向。

## 10. 软件回归要求

- Kagome 与 Square 的 Seed `0–47 scaf` 必须逐 base 相同。
- Kagome SST 活动 helix 必须恰为 12 条，`53/55/61/63` 完全为空。
- Origin reference 中必须有 72 条 Seed–SST staple bridges，且是 Square 96 条桥的严格子集。
- 缺少的 24 条必须全部对应 A0 phase 的 `61→2`、`63→0`。
- 当前结构不得提前包含 translated Seed helices `4–7,24–27` 的桥。
- translated-half 导出必须依据平移后的活动 mask、接触几何和 base phase重新计算。
- 所有 `prev/next` 互反；capture extension 只能来自实际 SST 组件路径。
- 对 32-nt extension，Seed core 必须为 32 nt；完整产品不得超过 64 nt。
- 以未修改的 Kagome SST 为基线，三类 capture 端点数必须为：crossover `48`、原有 nick `18`、线型连续 `6`。
- crossover 必须成对拆除且无悬空旧搭档；连续链必须双向互反切开；不能把修改过程中产生的新 nick 当成原始 nick。
