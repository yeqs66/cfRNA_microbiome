# v9.0 论文与代码核对清单

基准：`manu_v9.0.docx`，题目为 *Longitudinal Multi-omics Profiling Reveals Systematic Host–Microbiome Metabolic Dysregulation in Gestational Diabetes Mellitus*。

本清单记录待作者核对的问题，不替作者修改统计方法，也不判定论文结果错误。只读取论文及授权范围内的 R 脚本，没有读取其引用的真实数据文件。按作者 2026-11-03 的要求，不执行分析、模型训练或真实数值复现。**语法可解析不代表统计正确、数据可得或结果已复现。**

`reference/` 的编号对应 `script_inventory.json` 中的 `source`，源代码行号保持原有顺序；绝对路径替换为配置标记，注释中私人路径脱敏。下列链接均指向可审阅的代码副本。状态为“待作者核对”的条目在确认前不自动修正。

## 1. 缺失和版本完整性

| 编号 | 论文位置 | 代码证据/位置 | 差异与影响 | 核对建议 | 状态 |
|---|---|---|---|---|---|
| F01 | 全文及 Methods | [全部脚本清单](script_inventory.md) | 盘点 139 个 R 脚本；其中 54 个存在空字节或解码问题。部分 `script/` 副本仅说明头可读，不能作为可复现原稿。 | 从有效原稿或历史版本恢复；保留损坏文件的哈希和来源记录，不假定同名版本相同。 | 已排除损坏文件；待作者核对版本 |
| F02 | Methods: cfRNA differential expression；Fig.2 | [miRNA](../reference/cfrna/01.R)、[lncRNA](../reference/cfrna/02.R)、[敏感性脚本](../reference/cfrna/03.R) | 找到 miRNA/lncRNA 主分析和 BMI 调整脚本，但没有确认一份完整、独立、与 v9 所有条件一致的 mRNA 主分析及 design 1–4 全套入口。读取旧拟合对象不等于提供建模代码。 | 提供生成 40 个差异 mRNA 及三套敏感性结果的准确版本与调用参数。 | 待补原实现/确认 |
| F03 | Methods: preprocessing、Deconformer、GTDB-Tk、TissueEnrich | [细胞来源绘图](../reference/cell_origin/01.R)、[TissueEnrich 绘图](../reference/wgcna/04.R) | R 脚本主要读取相应结果；没有找到可确认的 fastp/bowtie2/MetaPhlAn/HUMAnN、Deconformer 推断、GTDB-Tk 和 TissueEnrich 在线提交的完整原始执行流程。R 脚本范围不包含可能存在的 shell/Python 文件。 | 补充实际命令、软件/数据库版本、参考文件来源或在线操作说明；未经授权不查找其他类型的本地文件。 | 待补外部步骤 |
| F04 | 论文数据和图表 | [输入契约](input_contracts.md) | 很多阶段依赖外部 CSV、Excel、RData、模型 HDF5 及先前交互式环境；这些文件未授权读取，也未随公开代码提交。 | 通过本地配置提供输入；核对对象结构及生成步骤。 | 输入缺口已登记 |

## 2. cfRNA 差异、富集与组织分析

| 编号 | 论文位置 | 代码证据/位置 | 差异与影响 | 核对建议 | 状态 |
|---|---|---|---|---|---|
| L02 | Methods: Sample and gene filtering | [miRNA 13 行](../reference/cfrna/01.R#L13)、[lncRNA 14 行](../reference/cfrna/02.R#L14) | 稿件为至少 1 count、至少 10% 样本；代码为 `sum(x>1) > floor(ncol/10)`，即 count 和样本比例边界都更严格。样本级 6000/1900/200 检出条件的原始生成步骤尚未完整确认。 | 确认稿件描述或代码条件哪一项对应最终结果；不直接替换。 | 原条件保留，真实运行需明确确认 |
| L03 | Methods: 各孕期 GDM 对比 | [miRNA 23–24 行](../reference/cfrna/01.R#L23)、[lncRNA 24–25 行](../reference/cfrna/02.R#L24) | 固定 7 项对比向量依赖设计矩阵列顺序和因子基准。第三孕期向量含两个交互项，而 `~trimester*disease+batch` 的处理编码下通常只应加对应第三孕期交互项。还需确认 disease 基准是否确为 Control。 | 用原 metadata 打印列名和 contrasts，逐项核对。当前接口保留向量、限制为 7 列满秩矩阵，不自动改成另一个对比。 | 高优先级待作者核对 |
| L04 | Methods: 两轮 voom / participant block | [miRNA 32–39 行](../reference/cfrna/01.R#L32) | 第二次 `voomWithQualityWeights` 传 correlation 未显式传 block；`lmFit` 中传入 ID block。与稿件“第二轮 block”描述是否一致待确认。 | 提供最终调用记录及 limma 版本。 | 保留原调用，待作者核对 |
| L05 | Methods / Results: DEG FDR、log2FC 0.15 | [miRNA 44–55 行](../reference/cfrna/01.R#L44)、[富集](../reference/cfrna/04.R) | `topTable` 导出全部基因；需区分导出表、差异基因过滤、富集排序、绘图阈值和 FDR 检验范围。不能用绘图筛选自动代替正式 DEG 定义。 | 确认每类 RNA、每个孕期的最终过滤与排序规则。 | 待作者核对 |
| W01 | Methods: WGCNA | [构网 21–43 行](../reference/wgcna/01.R#L21)、[矩阵准备](../reference/metadata/02.R) | β=4 可定位；TOMType/networkType/cor_method 来自位置参数。稿件叙述 mRNA、miRNA、lncRNA 合并，但也提“两矩阵”；1359 样本/12541 基因与模块编号需要原输入和调用确认。 | 提供实际命令、最终输入矩阵哈希和模块编号映射。不能通过基因数排序猜测哪个模块是论文 module 7/13。 | 待作者核对 |
| W02 | Fig.3、Supplementary Tables 7–9 | [构网 80–119 行](../reference/wgcna/01.R#L80)、[热图](../reference/wgcna/03.R) | 存在 `corPvalueStudent` 及后续绘图筛选。须区分原始 Pearson P 和 BH FDR，确认校正是在每个性状内还是整张表内进行，以及重复孕期观测如何说明。 | 确认用于星号与正文数值的是哪个表；保留两类值，避免标签混用。 | 待作者核对 |

## 3. 微生物、代谢物与 BCAA 网络

| 编号 | 论文位置 | 代码证据/位置 | 差异与影响 | 核对建议 | 状态 |
|---|---|---|---|---|---|
| K01 | Methods: taxa MaAsLin2；KO ZINB | [物种](../reference/microbiome/01.R)、[KO](../reference/microbiome/02.R) | 名为 `masslin2for_cfRNA.R` 的脚本实际需按输入确认用途；不能凭文件名分类。KO 原实现含拟合失败重试/不同分布分支，需说明最终采用哪一模型及失败记录。 | 核对分组方向、ID 随机效应、零膨胀式、协变量、RPK 是否被当作计数及总体对比。 | 待作者核对 |
| K02 | Methods: species prevalence | 微生物过滤脚本与 Methods | 稿件“zero prevalence >5%”含义不清，可能与“非零检出率至少5%”混淆；不能擅自翻转筛选方向。 | 用最终代码/物种清单确认阈值和分母定义。 | 待作者核对 |
| M01 | Results BCAA / Methods metabolome LMM | [LMM 54–71 行](../reference/metabolomics/01.R#L54)、[网络输入 80–97 行](../reference/bcaa_network/01.R#L80) | LMM 从外部名单选择 BCAA+SCFA；网络输入是外部 BCAA 名单与 9 个指定 V 编号的并集，不保证等于论文列出的 9 种代谢物。名单未授权读取。 | 确认最终九种 Compound/V 编号，以及差异分析 FDR 是否针对九种还是更广的代谢物集合。 | 待作者核对；不缩小校正集合 |
| M02 | Methods: adjusted LMM；Supplementary Table 14 | [LMM 83–128 行](../reference/metabolomics/01.R#L83) | 同时输出 primary/adjusted；正文应对应哪一模型需确认。FDR 按 model×trimester、对全部输入代谢物分别计算。`all` 是三个孕期等权边际对比，不是简单合并 Spearman 或原始浓度总体比较。 | 核对论文 P/FDR 与模型标签、对比权重。 | 原函数已抽取，未运行 |
| M03 | 代谢物效应量 | [LMM 76–81、94–115 行](../reference/metabolomics/01.R#L76) | 模型在原始浓度上拟合，另计算带 +1 的均值比 log2FC；log2FC 不是 LMM 回归系数。`n` 是整模型样本量，不是每孕期对比样本数；原函数抑制部分拟合 warning。 | 报告尺度与 n 含义；检查奇异拟合、收敛及协变量缺失造成的样本减少。 | 待作者核对；原统计表达保留 |
| N01 | Methods Network、Supplementary Fig.9 | [网络输入](../reference/bcaa_network/01.R)、[Spearman](../reference/bcaa_network/02.R) | 输入脚本含 PCA/历史交互式代码，不能直接假定最终网络使用了哪个对象；Spearman 函数本身与 v9 的 min_n=20、缺失≤50%、KO检出≥10%、BH、|rho|≥0.10 等基本一致。 | 确认最终输入是 log 转换后的单特征矩阵，准确的 KO/cfRNA/代谢物列表及样本来源。 | 独立 network 入口使用显式 feature_sets，不自动猜测 |
| N02 | Methods Network / Results | [网络函数 108–125 行](../reference/bcaa_network/02.R#L108)、[FDR 190–205 行](../reference/bcaa_network/02.R#L190) | 1562 为 cfRNA 样本池，1387 为 cfRNA/代谢组重叠样本；每条边的有效样本数不必相同。三孕期联合 BH，all 单独 BH；all 未建模同一孕妇重复观测。 | 以实际配对边 n、unique ID 和三组学重叠样本计数解释；all 的 P/FDR 标为探索性名义结果。 | 原函数已保留计数与说明 |
| N03 | 网络可视化 | [网络函数 233–238 行](../reference/bcaa_network/02.R#L233) | 前80条仅为展示子集，不是全部显著边。Cytoscape 3.10.1 的布局和样式文件未找到已确认来源。 | 提供最终 Cytoscape session/style，或说明布局需人工操作；不得将展示边数描述为总关联数。 | GraphML/边表导出接口已整理；未执行 |

## 4. MOFA、靶基因与预测

| 编号 | 论文位置 | 代码证据/位置 | 差异与影响 | 核对建议 | 状态 |
| O01 | Methods MOFA input | [模型 39–56、148–150 行](../reference/mofa/01.R#L39)、[后处理 25–84 行](../reference/mofa/02.R#L25) | 部分 KO 输入来自 rpk 命名文件并做 log2，稿件称 RPM；还有将缺失样本补0的历史逻辑。需核对最终四视图、高变特征筛选和1476样本交集。 | 明确 TPM/RPK/RPM 单位、变换、每视图特征数、缺失处理和最终模型文件。 | 待作者核对 |
| O02 | Methods/Results: factors clinical BH FDR | [154 行](../reference/mofa/02.R#L154) | `p.adjust(cor_mofa_p[[i]])` 未给 method，R 默认是 Holm，不是 BH；循环按列校正，不是自动整张表联合校正。 | 确认论文所用 FDR 表的生成版本；不要直接改方法后沿用稿件现有数值。 | 高优先级待作者核对 |
| O03 | Results: top10% absolute loadings；Methods top20 miRNA | [478–514 行](../reference/mofa/02.R#L478) | 存在 abs 权重、按十分之一截取及按全部行截取的不同分支；需核对每视图/孕期/因子取前10%还是合并排序，以及最终靶基因数据库支持数。 | 确认用来生成 Supplementary Tables 19–22 的具体分支、数据库版本与运行时间。 | 待作者核对 |
| P01 | Methods integrated GDM prediction；Fig.6B | [现有 Keras 代码](../reference/prediction/01.R) | 全部本地 R 脚本及远端既有预测文件中未找到可确认的 nnet 5×3 nested-CV、176/90特征、8/12隐层、skip配置、5次置换和seed74701全套原实现；现有代码主要是 Keras。 | 补交生成 AUC0.811 和经验 P0.167 的确切脚本。不可用旧 Keras 结果替代，也不可按稿件新写一个模型并冒称原实现。 | 原实现待补；无自动重写 |
| P02 | Methods offspring BMI；Fig.6C–E | [48、97–170 行](../reference/prediction/03.R#L48) | 按样本行生成 folds；需要确认每孕期是否每位孕妇仅一行。参数网格在同一套折上比较，并非可确认的独立外层评估。 | 核对受试者划分、特征选择是否在训练折内，以及调参/最终评价是否分离；保留历史实现，不自行改CV方案。 | 待作者核对 |
| P03 | Results R²0.596 | [170 行](../reference/prediction/03.R#L170) | 代码为 Pearson `cor(y_test,pred)^2`，不同于 `1-SSE/SST` 的预测决定系数。 | 明确图和正文 R² 定义、MSE/RMSE 对应关系及 pooled out-of-fold 计算方式。 | 高优先级待作者核对 |
| P04 | Methods external validation；Fig.6A | [预测 01.R](../reference/prediction/01.R) | 已找到外部数据处理片段，但重叠基因、缺失特征补0、训练均值标准化及外部合并评价是否对应最终版本需原数据结构确认；未读取外部数据。 | 提供特征清单及固定预处理参数来源，确认未在外部集重新拟合。 | 待作者核对 |

## 5. 文稿内部及交付范围

- Supplementary Table 14 标题仍写 Wilcoxon，但 Methods/Results 已写 LMM；Supplementary Table 15 标题留空。
- Supplemental PDF 列表写 S1–S13，但方法中引用至 Supplementary Figure 14；需统一编号。
- 旧仓库 README 声称全部脚本为模拟数据演示，实际也包含分析代码；新 README 明确区分历史代码、来源副本与整理入口。旧提交完整保留在备份分支。
- 当前独立封装：network、metabolomics_lmm、limma_mirna、limma_lncrna；仅完成代码整理与静态检查。其余阶段通过来源副本和需作者确认的 reference recipe 提供，不能宣称全流程已经可运行。
- 本仓库不发布真实数据或论文全文；输入文件路径和作者身份相关的本机路径被移除。`script_inventory.json` 保留原脚本相对路径、哈希和关键代码行，供逐项核对。
- **验收限制：没有执行本次新增入口的统计计算、模拟模型或真实队列复现，没有验证论文中的任何 AUC、R²、P/FDR 或样本数。** 缺失实现和未确认参数保留待核对状态。

## 作者核对记录（留待填写）

| 问题编号 | 确认的最终脚本/版本 | 确认的参数/定义 | 是否更新稿件 | 核对人/日期 |
|---|---|---|---|---|
| 待填写 | | | | |
