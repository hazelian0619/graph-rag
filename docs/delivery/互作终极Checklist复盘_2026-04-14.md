# 互作终极 Checklist 复盘（结果导向）

日期：2026-04-14  
范围：PPI / PSI / RPI + Cross-validation + Ontology（互作线）  
本地结果源：`/Users/pluviophile/graph/1218`  
仓库实现源：`/Users/pluviophile/graph/protian-entity`（`main`）

判定口径：
- **本地结果**：`data/output/evidence/**` 是否存在，且对应 `pipelines/**/reports/*.validation.json` 为 PASS
- **仓库实现**：对应 pipeline 的 `README/run.sh/contracts/scripts/reports` 是否已入库（git tracked）
- **发布可交付**：GitHub Release 是否可下载、是否带 manifest 与 QA 证据

---

## A. 总体结论

- 互作线已形成“**可提审 + 可下载**”闭环：
  1) PR #8 + PR #10 已合并到 `main`（互作主线 + PSI 条件增强入库）  
  2) Release `interaction-l2-v1` 已发布（含数据资产 + manifest + QA 证据）  
- A/B/C/D/E 核心门禁全部 PASS，且 PSI 条件增强 v3 硬门禁已 PASS。
- 当前剩余的非阻塞缺口：
  - **RPI 位点坐标字段为空**（`chromosome/start/end` 覆盖率 0）
  - **PSI 条件字段已显著改善但仍未达目标门禁**（`condition_context=0.5969`，`pH=0.2970`，`temperature=0.2944`；context<0.65, pH<0.30）

---

## B. 终极 Checklist 逐项复盘（逐字段严谨对齐）

| 组件 | 属性字段 | 本地结果产物（存在性） | 本地结果指标（实测） | 仓库实现（tracked） | 结论 |
|---|---|---|---|---|---|
| PPI 实验证据 | 检测方法/通量 | `ppi_method_context_v2.tsv` ✅ | `rows=884,555`；`method=1.0`；`throughput=1.0` | `ppi_semantic_enrichment` ✅ | ✅ 完成 |
| PPI 实验证据 | PMID/DOI | 同上 ✅ | `pmid=1.0`；`doi=1.0` | 同上 ✅ | ✅ 完成 |
| PPI 置信度 | experimental/text-mining 评分 | 同上 ✅ | `experimental_score=1.0`；`text_mining_score=1.0` | 同上 ✅ | ✅ 完成 |
| PPI 功能上下文 | GO/Reactome/KEGG 共现 | `ppi_function_context_v2.tsv` ✅ | `any_go_rate=0.8124`；`any_pathway_rate=0.5639` | 同上 ✅ | ✅ 完成（覆盖合理） |
| PSI 活性证据 | 活性类型（IC50/Ki/Kd/EC50） | `psi_activity_context_v2.tsv` ✅ | `standard_type=1.0`；分布：IC50 1,456,252 / Ki 560,388 / Kd 139,548 / EC50 203,835 | `psi_activity_structure_enrichment` ✅ | ✅ 完成 |
| PSI 活性证据 | 关系符（= < > >= ~） | 同上 ✅ | `standard_relation=1.0`；`=`占比最高 | 同上 ✅ | ✅ 完成 |
| PSI assay | assay 类型/ID | 同上 ✅ | `assay_type=1.0`；`assay_id=1.0` | 同上 ✅ | ✅ 完成 |
| PSI 条件上下文 | pH/温度/condition | `psi_activity_context_v3.tsv` ✅ | v2→v3：`condition_context 0.4266→0.5969`；`pH 0.1522→0.2970`；`temperature 0.0984→0.2944` | `psi_condition_enrichment`（PR #10）✅ | ✅ 硬门禁通过；目标门禁部分未达 |
| PSI 结构证据 | PDB 结构证据 | `psi_structure_evidence_v2.tsv` ✅ | `rows=1,808,985`；`pdb_id=1.0`；占 PSI 活性边覆盖 `0.7665` | 同上 ✅ | ✅ 完成（达门禁） |
| RPI 方法 | 方法类型（RIP/CLIP） | `rpi_site_context_v2.tsv` ✅ | `method_type=1.0`；分布：CLIP=100,255 | `rpi_site_domain_enrichment` ✅ | ✅ 完成 |
| RPI 位点 | 基因组/转录本坐标 | 同上 ✅ | `chromosome/start/end/transcript_id/gene_id = 0.0` | 同上 ✅ | 🟡 表已交付，但位点坐标未实填 |
| RPI 结构域 | RRM/KH/Pfam | `rpi_domain_context_v2.tsv` ✅ | `domain_class=1.0`；`PFAM 21,344 / RRM 7,999 / KH 1,576` | 同上 ✅ | ✅ 完成 |
| RPI 功能关系 | 转录/翻译/剪切调控 | `rpi_function_context_v2.tsv` ✅ | `function_relation=1.0`；三类关系均有覆盖 | 同上 ✅ | ✅ 完成 |
| 跨库一致性 | consistent_across_n | `interaction_cross_validation_v2.tsv` ✅ | 分布：`1:575,485 / 2:2,524,166 / 3:119,400 / 4:74,138` | `interaction_cross_validation` ✅ | ✅ 完成（非全1） |
| 冲突标记 | conflict_flag | 同上 ✅ | `true=21,396`；`false=3,271,793` | 同上 ✅ | ✅ 完成 |
| 聚合评分 | 0-1 评分分布 | `interaction_aggregate_score_v2.tsv` ✅ | `min=0.0` `p10=0.4268` `mean=0.5351` `p90=0.6324` `max=0.7197` | 同上 ✅ | ✅ 完成（非常数） |
| 本体标准化 | PSI-MI / GO / 谓词统一 | `interaction_ontology_mapping_v2.tsv` ✅ | `rows=1,905,589`；`ontology_uri=1.0`；谓词含 `detected_by / supported_by / participates_in / standardized_as` | `interaction_ontology_mapping` ✅ | ✅ 完成 |

---

## C. 本地核心结果文件清单（逐表核对）

| 文件 | 行数 | 本地存在 | validation |
|---|---:|---|---|
| `data/output/evidence/ppi_method_context_v2.tsv` | 884,555 | ✅ | PASS |
| `data/output/evidence/ppi_function_context_v2.tsv` | 884,555 | ✅ | PASS |
| `data/output/evidence/psi_activity_context_v2.tsv` | 2,360,023 | ✅ | PASS |
| `data/output/evidence/psi_activity_context_v3.tsv` | 2,360,023 | ✅ | PASS |
| `data/output/evidence/psi_condition_parse_audit_v3.tsv` | 199,742 | ✅ | PASS |
| `data/output/evidence/psi_structure_evidence_v2.tsv` | 1,808,985 | ✅ | PASS |
| `data/output/evidence/rpi_site_context_v2.tsv` | 100,255 | ✅ | PASS |
| `data/output/evidence/rpi_domain_context_v2.tsv` | 30,919 | ✅ | PASS |
| `data/output/evidence/rpi_function_context_v2.tsv` | 39,254 | ✅ | PASS |
| `data/output/evidence/interaction_cross_validation_v2.tsv` | 3,293,189 | ✅ | PASS |
| `data/output/evidence/interaction_aggregate_score_v2.tsv` | 3,293,189 | ✅ | PASS |
| `data/output/evidence/interaction_ontology_mapping_v2.tsv` | 1,905,589 | ✅ | PASS |

---

## D. 仓库与发布状态（你关心的“是否真的可交付”）

### D.1 仓库实现（`main`）

- PR：`#8` 已合并（MERGED，互作主线）  
  - URL: `https://github.com/hazelian0619/protian-entity/pull/8`
  - merge commit: `c9a4665dc3db680147904b5e57c7ac89f7019e45`
- PR：`#10` 已合并（MERGED，PSI 条件增强）  
  - URL: `https://github.com/hazelian0619/protian-entity/pull/10`
  - merge commit: `e9c5acdebcbd8bb39383a0ba7487046ccd040773`
- 以下 7 条互作 pipeline 均已 tracked：
  - `pipelines/ppi_semantic_enrichment`
  - `pipelines/psi_activity_structure_enrichment`
  - `pipelines/psi_condition_enrichment`
  - `pipelines/rpi_site_domain_enrichment`
  - `pipelines/interaction_cross_validation`
  - `pipelines/interaction_ontology_mapping`
  - `pipelines/interaction_release_local`

### D.2 数据发布（Release）

- Release tag：`interaction-l2-v1`
- URL：`https://github.com/hazelian0619/protian-entity/releases/tag/interaction-l2-v1`
- 发布状态：已发布（非 draft）
- 口径说明：当前 release 仍以 v2 主数据资产为主，`psi_activity_context_v3` 尚未单独发布新 tag（建议后续发 `interaction-l2-v1.1` 或 `interaction-l2-v2`）
- 资产结构：
  - 常规资产：manifest、release notes、QA/gates JSON、中等体量 `.tsv.zst`
  - 大体量资产采用分片：共 17 个 part 文件  
    - `psi_activity_context_v2.tsv.zst.part.000..009`  
    - `psi_structure_evidence_v2.tsv.zst.part.000..002`  
    - `interaction_aggregate_score_v2.tsv.zst.part.000..003`  

重组示例：

```bash
cat psi_activity_context_v2.tsv.zst.part.* > psi_activity_context_v2.tsv.zst
cat psi_structure_evidence_v2.tsv.zst.part.* > psi_structure_evidence_v2.tsv.zst
cat interaction_aggregate_score_v2.tsv.zst.part.* > interaction_aggregate_score_v2.tsv.zst
```

---

## E. 最终判定（结果导向）

1. 互作线目前已经达到“**工程可复现 + 数据可下载 + 门禁可审计**”的交付状态。  
2. A/B/C/D/E 验收口径全部 PASS，且已完成 main 合并与 release 发布。  
3. 当前不构成阻塞但需后续优化的点：
   - RPI 位点坐标信息（site 字段实填率）
   - PSI 条件字段目标门禁（context≥0.65、pH≥0.30）未完全达成

> 结论：互作这条线可视为“本阶段可提审、可交付、可复现”的完成态。
