# RNA 终极 Checklist 复盘（结果导向）

日期：2026-04-13  
仓库：`/Users/pluviophile/graph/protian-entity`（`main`）

判定口径：
- **本地结果**：`data/output/**` 与 `pipelines/**/reports/*.validation.json` 是否存在且 PASS
- **仓库实现**：对应 `pipelines/**` 脚本 + contract 是否已入库（git tracked）
- 说明：按仓库策略，`data/output/**` 与 `pipelines/**/reports/**` 默认不入库（`.gitignore`）

---

## A. 总体结论

- RNA checklist 对应的数据集，当前已形成可交付结果链路。  
- 本地 RNA 相关 validation（非 sample）共 **15/15 PASS**。  
- “缺口”层面：
  - **结构性特征项**（ENST↔URS、Rfam全覆盖、PDB稀疏）= 已有结果但覆盖受限，非当前必解
  - **可持续优化项** = `tRNA anticodon`（v2 已上线，较 v1 小幅提升）

---

## B. 终极 Checklist 逐项复盘

| 组件 | 属性字段 | 本地结果产物（存在性） | 本地结果指标 | 仓库对应实现（tracked） | 结论 |
|---|---|---|---|---|---|
| 序列数据 | RNA序列 | `data/output/rna_master_v1.tsv` ✅ | 292,274 行；`sequence_non_empty_rate=1.0`；字符集规则 1.0 | `pipelines/rna/`（run+scripts+contract）✅ | ✅ 完成 |
| 序列数据 | RNAcentral URS | 主表+映射层 ✅ | 主表 URS 比例 0.9005%（miRNA 主）；ENST↔URS v2=10 行；canonical 非空率 100% | `rna_xref_enst_urs` + `rna_id_canonical` + `rna_external_xref` ✅ | 🟡 结果可用，覆盖受结构性限制 |
| 序列数据 | RNA类型标签 | 主表+分表 ✅ | 主表：`mrna=289,642`、`mirna=2,632`；分表补齐 `lnc/tRNA/rRNA` | `rna` + `rna_type_features` ✅ | ✅ 完成（分层实现） |
| 序列数据 | 物种来源 | 主表及分表 ✅ | `taxon_id=9606` 规则通过（主表 1.0） | 各 contract 均有 taxon 规则 ✅ | ✅ 完成 |
| 结构数据 | 二级结构 | `data/output/rna_structure_rfam_v1.tsv` ✅ | 1,096 行；validation PASS | `rna_rfam_structure` ✅ | 🟡 已交付，覆盖受限 |
| 结构数据 | Rfam家族ID | 同上 ✅ | Rfam 结果与二级结构同表；validation PASS | `rna_rfam_structure` ✅ | 🟡 已交付，覆盖受限 |
| 结构数据 | 协方差模型 | `data/output/rna_covariance_models_index_v1.tsv` ✅ | 4,178 行；validation PASS | `rna_structure_models` ✅ | ✅ 完成 |
| 结构数据 | 三维结构 | `data/output/rna_pdb_structures_v1.tsv` ✅ | 1 行；validation PASS | `rna_pdb` ✅ | 🟡 已交付，数据稀疏 |
| 结构数据 | 预测结构 | `data/output/rna_predicted_structures_v1.tsv` ✅ | 3,617 行；validation PASS | `rna_structure_models` ✅ | ✅ 完成 |
| 功能注释 | miRNA条目 | 主表 + xref ✅ | 主表 miRNA=2,632；外部映射含 MIRBASE 记录 | `rna` + `rna_external_xref` ✅ | ✅ 完成 |
| 功能注释 | lncRNA条目 | `data/output/rna_lnc_entries_v1.tsv` ✅ | 189,152 行；validation PASS | `rna_type_features` ✅ | ✅ 完成 |
| 功能注释 | 表达与证据 | `rna_expression_evidence_v1.tsv` / `rna_rbp_sites_v1.tsv` ✅ | 289,614 / 31,084 行；join_rate=1.0；validation PASS | `rna_expression_rbp` ✅ | ✅ 完成 |
| 类型特征 | tRNA反密码子 | `rna_trna_features_v1.tsv` + `rna_trna_features_v2.tsv` ✅ | v1=0.113781；v2=0.119176（+22 行）；v2 validation PASS | `rna_type_features`（含 v2 脚本/contract/test）✅ | 🔶 持续优化项（已上线 v2） |
| 类型特征 | rRNA位点 | `data/output/rna_rrna_loci_v1.tsv` ✅ | 69,005 行；validation PASS | `rna_type_features` ✅ | ✅ 完成 |
| 互作锚点 | RBP结合位点 | `data/output/rna_rbp_sites_v1.tsv` ✅ | 31,084 行；validation PASS | `rna_expression_rbp` ✅ | ✅ 完成 |
| 标识映射 | 外部ID映射 | `data/output/rna_external_xref_v1.tsv` ✅ | 297,813 行；主表回链率 1.0；validation PASS | `rna_external_xref` ✅ | ✅ 完成 |

---

## C. 本地结果文件清单（核心）

| 文件 | 行数 | 本地存在 | validation |
|---|---:|---|---|
| `data/output/rna_master_v1.tsv` | 292,274 | ✅ | （主表 contract 在 `pipelines/rna/contracts`） |
| `data/output/rna_xref_mrna_enst_urs_v2.tsv` | 10 | ✅ | PASS |
| `data/output/rna_id_canonical_map_v1.tsv` | 292,274 | ✅ | PASS |
| `data/output/rna_external_xref_v1.tsv` | 297,813 | ✅ | PASS |
| `data/output/rna_structure_rfam_v1.tsv` | 1,096 | ✅ | PASS |
| `data/output/rna_covariance_models_index_v1.tsv` | 4,178 | ✅ | PASS |
| `data/output/rna_predicted_structures_v1.tsv` | 3,617 | ✅ | PASS |
| `data/output/rna_pdb_structures_v1.tsv` | 1 | ✅ | PASS |
| `data/output/rna_expression_evidence_v1.tsv` | 289,614 | ✅ | PASS |
| `data/output/rna_rbp_sites_v1.tsv` | 31,084 | ✅ | PASS |
| `data/output/rna_lnc_entries_v1.tsv` | 189,152 | ✅ | PASS |
| `data/output/rna_trna_features_v1.tsv` | 4,078 | ✅ | PASS |
| `data/output/rna_trna_features_v2.tsv` | 4,078 | ✅ | PASS |
| `data/output/rna_rrna_loci_v1.tsv` | 69,005 | ✅ | PASS |

---

## D. 仓库资产状态（你问的“仓库是否有对应产物”）

- **有**：对应 pipeline 代码、run.sh、contracts、README 都已入库（tracked）。
- **无（按策略）**：`data/output/**` 与 `pipelines/**/reports/**` 结果文件不入库（本地产物）。

结论：
- 从“工程可复现能力”看：**仓库完备**。  
- 从“数据集结果可交付”看：**本地结果完备且验证通过**。

---

## E. 最终判定（结果导向）

1. 你给的终极 checklist，当前已具备结果型交付。  
2. 真正还需要持续迭代的只剩：`tRNA anticodon`（已从 v1 推进到 v2）。  
3. 其他低覆盖项（ENST↔URS、Rfam 全量、PDB）属于结构性边界，不建议当作当前阻塞项。
