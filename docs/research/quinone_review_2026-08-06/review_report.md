# iYali26 Quinone Biosynthesis 多代理文献审查报告

日期：2026-08-06  
状态：只读证据审查；尚未授权模型、curation、代码或 Obsidian 修改  
冻结模型 SHA-256：`0f3a6c2b151e945b3461d3fa85f04575f8e8570ba817ed2879013aec91f62415`

## 结论摘要

这次审查把原先混在一起的三个问题拆开后，得到三个不同强度的结论：

1. **链长/终产物问题：对当前 CoQ6 表示构成强挑战。** 在已测试的
   *Yarrowia lipolytica* isolates 中，Q9 是主要检测产物；保留的 Yl Coq1 catalytic
   core 在异源宿主中也生成 Q9。当前 `R763` 生成 C30 hexaprenyl diphosphate、
   主路终止于 ubiquinol-6，因此当前 CoQ6 化学与这些证据冲突。目标模型所对应的
   W29/PO1f/CLIB89 背景仍缺少条件匹配的原生全细胞确认。
2. **净需求问题：方向可辩护，系数不可辩护。** 细胞中存在有限的 CoQ9 pool，
   PO1f SD-Leu 筛选得到的位点 fitness 与 `COQ1/COQ2` 候选注释的 pathway
   relevance 一致，但不能单独证明 pathway-specific growth requirement。若模型
   要让新生物量维持固定 CoQ9 浓度，growth-coupled dilution 是合理概念；但
   当前没有 W29/PO1f、SD-Leu 条件下的 total CoQ9 mmol/gDCW，因此不能给出
   正式数值系数。
3. **GPR 问题：数据库候选身份可解析，原生功能大多仍靠同源推断。** 当前八个位点
   可映射到数据库中的 COQ2–COQ9 candidate annotations。保守机制支持
   step-specific 直接催化分工和动态 CoQ synthome；相同七基因 `AND` 不能被解释为
   每一步都有七个直接催化酶。它是否恰好描述了原生 *Yarrowia* 中维持各步通量所需
   的完整 accessory dependency，仍未验证。

这些结论不等价于已经确定了一份可直接写入模型的补丁。尤其是原生定位、早期
中间体、`R39/R808` 区室与运输、COQ8/COQ9 的强制依赖条件以及 demand 系数仍未
解决。

## 审查设计

三条证据线分别独立检索：

- native CoQ9 化学、链长和定位；
- pool abundance、dilution、turnover、enzyme coupling 与条件匹配的定量数据；
- COQ2–COQ9 身份、step-specific 反应归属和 repeated-AND provenance。

所有证据线都要求优先使用直接 *Yarrowia* 实验，分别记录菌株、培养条件、样品
类型和方法。数据库只用于身份桥接；其他酵母用于机制边界；多个继承同一规则的
GEM 不计作多次独立验证。独立审计员随后逐条打开高影响来源，检查 exact locator、
条件、直接性和反证。代理一致性本身不计作证据。

## 1. CoQ9 化学与定位

### 直接或接近直接的 *Yarrowia* 证据

- Yamada et al. 1976 从完整细胞提取 quinone，以反相纸层析并在研究中结合质谱，
  将 *Endomycopsis lipolytica* `CBS 6124` 归为 Q9 system。
- Dröse et al. 2002 从 PIPO 的线粒体 complex I 制备物中以 RP-HPLC 检出 Q9，
  两批蛋白水解后的回收量为 `0.2` 与 `0.4 mol Q9/mol complex I`。
- Saeed et al. 2024 将 *Yarrowia* GeneID `2909963` 的 COQ1 catalytic region 放入
  *S. cerevisiae* COQ1 locus；TLC/HPLC 显示其生成 CoQ9 而非 CoQ6。

其中最重要的证据边界是：Saeed 构建物删除了 *Yarrowia* 原生 N 端，并使用
*S. cerevisiae* mitochondrial import signal。因此它直接支持 COQ1 catalytic core
决定九个 isoprene 单元，但不能证明原生 *Yarrowia* COQ1 的定位序列。

### 对当前模型的含义

当前模型上游为：

\[
R763:\ IPP_{mi}+pentaprenyl\text{-}PP_{mi}
\rightarrow PP_i+hexaprenyl\text{-}PP_{mi}
\]

`R763` 的 GPR 是 `YALI1C26017g`。这条反应只产生 C30/六异戊二烯 donor；保留主路
最终产生 ubiquinol-6。由此，**模型化学与 Q9 产品证据发生直接链长冲突**。

可以确认 Q9 在线粒体呼吸中被使用，但不能从 complex-I-associated Q9 推出每个
合成步骤都在线粒体同一侧完成。比较性 *S. cerevisiae* 证据将 Coq2 与 downstream
ring modification 放在线粒体内膜/基质侧，因此当前胞质 `R39` 加 `R808` 回运与
保守机制不一致；尚无 *Yarrowia* 原生定位实验足以把这种不一致提升为物种内直接
反证。

## 2. 净 CoQ9 demand

### 直接支持什么

- 全细胞/提取物、纯化 complex I 和 cryo-EM/MS 研究共同证明 CoQ9 pool 非零。
- PO1f、SD-Leu、2% glucose 的筛选中：
  - `YALI1C26017g — COQ1 — CoQ9 side-chain synthase candidate`
    （异源功能支持；原生定位未验证）获得强但平台不完全一致的 essential evidence；
  - `YALI1F08349g — COQ2 — 4-hydroxybenzoate polyprenyltransferase candidate`
    （同源/家族支持；原生位点未验证）在多屏共识中得到条件性支持；
  - `YALI1B20835g — COQ3 — ubiquinone O-methyltransferase candidate`
    （同源支持；原生位点未验证）在两种 CRISPR 中均未被判 essential。

这些结果与“CoQ pathway 可能影响生长”一致，但扰动的是候选位点；在多数位点缺少
原生 locus 功能验证时，不能只靠 fitness 结果完成 pathway-specific 归因。fitness
score 也不是 pool size，不能反推出 mmol/gDCW。

### 为什么不能直接用现有数值

| 数据 | 原始分母 | 可回答的问题 | 不能回答的问题 |
|---|---|---|---|
| Kogan 1985：0.026% yield | 原摘要未说明 wet/dry/input denominator | 旧名菌株可提取 CoQ9 | PO1f 的 gDCW coefficient |
| DuPont：油中 0.2–0.3% w/w | extracted oil | ATCC 20362 油相含 CoQ9 | 每克干细胞的 CoQ9 |
| Dröse：0.2–0.4 Q9/CI | purified complex I | 制备物有 substoichiometric Q9 | whole-cell pool |
| Parey：1.9 Q9/CI | purified complex I | 另一制备物的 Q9 occupancy | 固定生理酶占有量 |

两组 complex-I occupancy 来自不同菌株/构建、detergent、纯化和分析方法，数值不可
直接比较，也不能移植成固定系数。CoQ 是 complex I、NDH2、
complex III 和 AOX 共享的移动电子载体，不适合无证据地建成“每个酶固定结合一份
CoQ”的 prosthetic-group coupling。

### 允许与不允许的表示

若条件匹配的全细胞 pool 为 \(c_Q\) mmol/gDCW，生长率为 \(\mu\)，仅考虑生长
稀释时：

\[
v_Q\geq \mu c_Q
\]

若另有一阶 turnover 常数 \(k_{deg}\)：

\[
v_Q\geq(\mu+k_{deg})c_Q
\]

第一式是有限池质量守恒的合理表示，但当前缺少合格的 \(c_Q\)。第二式还额外缺少
*Yarrowia* 的半衰期或 \(k_{deg}\)。因此：

- growth-coupled dilution：**定性可辩护；正式数值尚不可辩护**；
- turnover maintenance：**缺少参数证据**；
- 固定 enzyme-pool stoichiometry：**受共享载体机制和 occupancy 差异挑战**；
- 无 explicit demand：**只能表示“定量耦联未解析”，不能解释为生物学零需求**；
- `0…1000` sink：作为**建模治理规则（非文献主张）**，只可用于 pathway
  reachability counterfactual，不能作为生物学修复。

## 3. COQ synthome 与 step-specific GPR

### 当前位点的数据库候选身份

| 当前基因 | 名称/符号 | 候选蛋白功能 | 证据状态 | direct-catalyst candidate only |
|---|---|---|---|---|
| `YALI1F08349g` | `COQ2` | 4-HB polyprenyltransferase；UbiA family | TrEMBL/同源；无原生 locus 实验 | `R407` |
| `YALI1F34625g` | `COQ4` | CoQ ring oxidative decarboxylase/scaffold | reviewed entry 但 Yarrowia 功能仍为同源推断 | `R40` candidate |
| `YALI1B20527g` | `COQ8` | ADCK-family regulatory ATPase/kinase | TrEMBL/保守机制；无原生实验 | 非通用原子转移催化步骤 |
| `YALI1A08781g` | `COQ6` | FAD monooxygenase；两个非连续羟化 | TrEMBL/保守机制 | `R39`, `R19` candidates |
| `YALI1F34675g` | `COQ9` | synthome-associated accessory；影响 Coq6/Coq7 功能 | TrEMBL/保守机制 | 可能影响 `R695`，强制性未定 |
| `YALI1C25352g` | `COQ5` | C-methyltransferase | TrEMBL/保守机制 | `R18` candidate |
| `YALI1B20835g` | `COQ3` | 两次非连续 O-methylation | TrEMBL/保守机制 | `R715`, `R385` candidates |
| `YALI1E18269g` | `COQ7` | DMQ hydroxylase；兼有复合体稳定作用 | TrEMBL/保守机制 | `R695` candidate |

这里的名称都是 **candidate labels**。标识桥接总体稳定，但当前 UniProt live entry
并未把 `Q6CEH1` 明确命名为 COQ8，也未把 `Q6C071` 明确命名为 COQ9；这两项尤其
依赖额外的 orthology/family 推断。除 `Q6C074` 为 reviewed entry 外，其余多数为
TrEMBL/自动规则注释，不能按数据库数量当作多份独立功能验证。

### 为什么 repeated `AND` 有问题

当前 `R715/R19/R18/R695/R385` 都使用完全相同的七基因规则：

`COQ4 and COQ8 and COQ6 and COQ9 and COQ5 and COQ3 and COQ7`

这把两个不同层面的事实混在一起：

1. 多个 Coq 蛋白可组成动态、底物依赖的 synthome；
2. 每条具体化学反应由哪个蛋白直接催化。

保守系统、化学 bypass、突变体中间体和体外重建显示：直接催化是 step-specific，
复合体组成可变；COQ8/COQ9 是调节或辅助成员，不是每步通用的 atom-transfer
enzyme。因而，相同七基因 `AND` **不能被解释为七种 direct catalysts 的身份规则**。
但 GPR 也可能编码维持反应通量所需的稳定、递送和电子供体依赖。作为**比较性机制
提醒，而非原生 Yarrowia locus claim**，其他系统中 Coq6 存在电子供体 partner，
Coq9 也会影响 Coq6/Coq7 功能（`CLM-GR-005`, `CLM-GR-009`）。因此单个
direct-catalyst candidate 不自动等于完整 GPR。

另一方面，原生 *Yarrowia* synthome 的组成、化学计量和条件依赖性没有被直接
测量。因此审查不能进一步断言 COQ8/COQ9 在所有生长条件下都“不必要”。正确的
证据标签是：**当前 repeated-AND 的 direct-catalyst interpretation 不受支持；
native reaction-flux/accessory dependency unresolved**。

### 模型继承不是独立证据

相同七位点 conjunction 可从 iYali4 追踪到 iYali5、iYli21 和 iYali26。iYali5 明确
以 iYali4 为源，iYli21 又以 iYali4 为模板做标识转换。因此四代模型的一致性是一条
复制的 provenance lineage，不是四次实验验证。

## 4. 冲突与不确定性

1. **COQ1 功能实验 vs KEGG 标签：** 实验生成 CoQ9；KEGG 仍写 putative
   hexaprenyl pyrophosphate synthase。直接功能证据优先。
2. **Q9 必须存在 vs 无净需求：** 两者并不矛盾。FBA 中 redox carrier 可无限循环，
   因而 biosynthesis 可为零；缺的是 pool dilution representation。
3. **synthome interdependence vs step-specific catalysis：** 两者可以同时为真；
   step-specific catalysis 不能单独证明或否定完整的 reaction-flux dependency。
4. **0.2–0.4 vs 1.9 Q9/complex I：** 两者来自不同菌株/构建、制备和分析方法，
   不能据此选择一个固定 holoenzyme coefficient。
5. **PO1f candidate-locus phenotypes heterogeneous：** COQ1/2 candidates 进入
   多屏共识，而 COQ3 candidate 在两种 CRISPR 中均未判 essential。平台效应、
   位点身份、残余活性、冗余或真实生物学差异等解释尚未区分。COQ1/2 不属于当前
   repeated seven-gene `AND`，所以这个表型模式既不验证也不反驳该 GPR。

## 5. 明确的证据缺口

本次有界检索未找到：

- W29/PO1f、SD-Leu 条件下 total CoQ9 mmol/gDCW；
- CoQ9 pool 随生长率的稀释、分配或同位素追踪数据；
- *Yarrowia* CoQ9 turnover/half-life；
- W29/PO1f 的原生 COQ1/COQ2 亚细胞定位；
- *Yarrowia* 中 HHB、DMQ 等逐步中间体和 `R39/R808` transport 的直接检测；
- 原生 *Yarrowia* Coq synthome 组成、stoichiometry 和条件依赖；
- COQ8/COQ9 对某个具体 ring-modification reaction 构成绝对依赖的原生实验。

上述均是“本次没有找到”，不是“生物学上不存在”。

## 6. 独立审计状态

独立来源审计已依据 `audit_protocol.md` 封卷：

- material-claim coverage：`45/45 = 100%`；身份 crosswalk 已逐 locus 拆分，
  身份/定位、链长/氧化还原比例、稀释概念/数值参数以及本地模型事实均为独立原子主张；
- high-impact-source coverage：`23/26 = 88.5%`；
- 若排除两件专利灰色来源，peer-reviewed primary-study coverage 为
  `21/24 = 87.5%`；
- verdict 分布：`supported 22`、`partially_supported 7`、`unsupported 5`、
  `contradicted 6`、`unverified 5`。

未完整打开的三项被明确降级处理：Kogan 1985 只核到书目信息，因此 `0.026%`
原始分母仍不可验证；Schwartz 2019 与 Patterson 2018 未由独立审计员打开，不计作
独立重复。PO1f 的相关 fitness 数值和培养条件仅使用审计员直接打开的 Ramesh 2023
正文及其 Supplementary Data 4/7 复核。

独立反方审查进一步要求并已落实以下收窄：

- Q9 结论限定于已测试 isolates 和异源 catalytic-core 功能，不冒充目标菌株全细胞
  定量；
- 数据库 crosswalk 与 native locus function 分开表述；
- step-specific direct catalyst 与完整 reaction-flux/accessory dependency 分开；
- Nicoll 2024 明确标为 ancestral-sequence-reconstructed tetrapod proteins 的人工
  体外体系；
- complex-I occupancy 差异不归因于单一原因；
- `US20090142322A1/US8815567B2` 与 `WO2008073367A1` 的专利族关系已纠正。

因此本报告可作为**只读、带证据等级的 evidence map**，但不构成任何模型修改、
参数选择或 GPR patch 的批准。

## 永久人类闸门

本轮没有修改 curated tables、pipeline code、反应化学、GPR、bounds、biomass、
demand、`model.xml` 或 Obsidian。即使独立审计通过，也必须由用户另行明确批准，
才能进入“设计 counterfactual/curation patch”的下一阶段。
