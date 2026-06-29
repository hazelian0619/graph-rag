# Disease-gene / disease-target 边管线（`disease_gene_target_v1`）

本管线构建 disease ↔ 分子主干的核心边：disease-gene（疾病关联 / 遗传机制入口）与 disease-target（治疗 / 干预语义）。按 D-GATE-3，实现层收敛为统一 edge set，产品语义通过四个保留字段保留（DECISIONS 007 / 014）。

## 语义保留字段（不可丢失）

| field | 取值 |
|-------|------|
| `relation_semantics` | `disease_gene` / `disease_target` |
| `target_type` | `gene` / `protein_target` |
| `source` | `DisGeNET` / `Open_Targets` / ... |
| `evidence_class` | `curated_association` / `inferred_mapping` / `therapeutic_evidence` |

## 数据源

- DisGeNET → disease-gene curated associations
- Open Targets → disease-target therapeutic evidence

## 输入 / 输出

- 输入：`data/processed/disease_gene_target_v1.tsv`（上游 ingest 后的统一边表）
- 输出：
  - `disease/pipelines/disease_gene_target/reports/disease_gene_target_v1.validation.json`
  - `disease/pipelines/disease_gene_target/reports/disease_gene_target_v1.manifest.json`

## 状态

**scaffold-only**：源数据（DisGeNET / Open Targets）尚未本地可用。`run.sh` 在输入缺失时 `exit 2`，不伪造表。

## 字段 / 契约

见 `contracts/disease_gene_target_v1.json` 与 `specs/disease-layer-v1/contracts/disease-gene-target.md`。

## 运行

```bash
bash disease/pipelines/disease_gene_target/scripts/run.sh
```
