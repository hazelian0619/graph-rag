# Disease-drug 管线（`disease_drug_v1`）

构建 disease ↔ drug 边：**indication / treatment relation ONLY**（DECISIONS 006）。不纳入 contraindication / adverse event / safety。drug 必须映射到现有 molecule 主干（FR-010）。

## 语义

| field | 取值 |
|-------|------|
| `relation_type` | `indication` / `treatment`（**仅此两类**） |
| `source` | `DrugCentral` / `DrugBank` |
| `evidence_class` | `curated_indication` / `inferred_mapping` |
| `molecule_backbone_mapped` | `True` / `False`（drug 是否接入 molecule 主干） |

## 数据源

- DrugCentral → indication 字段
- DrugBank → indication 字段（受许可约束）
- 抽取范围受 indication-only 约束（DECISIONS 006 Implication），不抽取 AE/contraindication/safety 字段

## 输入 / 输出

- 输入：`data/processed/disease_drug_v1.tsv`（上游 ingest 的 indication 边）
- 输出：`disease/pipelines/disease_drug/reports/disease_drug_v1.{validation,manifest}.json`

## 状态

**scaffold-only**：源数据（DrugCentral / DrugBank）尚未本地可用。`run.sh` 在输入缺失时 `exit 2`，不伪造表。

## 运行

```bash
bash disease/pipelines/disease_drug/scripts/run.sh
```

## 契约

见 `contracts/disease_drug_v1.json` 与 `specs/disease-layer-v1/contracts/disease-drug.md`。T022 反向校验拒绝任何 contraindication/AE/safety 边。
