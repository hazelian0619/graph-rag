# Biomarker first-pass 管线（`biomarker_edges_v1`）

Biomarker v1 以 **relation layer** 形态存在（D-GATE-2，default-applied），非独立主表。源从 disease↔gene/protein 高价值映射派生（D-GATE-4，default-applied），标注 `bridge_only`。

## 形态

关系层：`biomarker_edges_v1`，桥接 disease ↔ gene/protein/rna。保留 `biomarker_id` / `biomarker_type` 以便未来升级为独立节点层。

## 字段

| field | 取值 |
|-------|------|
| `biomarker_type` | `gene_biomarker` / `protein_biomarker` / `rna_biomarker` |
| `evidence_class` | `bridge_only` / `evidence_backed`（v1 默认 bridge_only） |
| `source` | `DisGeNET_derived` / `Open_Targets_derived` / ... |

## 输入 / 输出

- 输入：`data/processed/biomarker_edges_v1.tsv`（从 disease-gene/protein 映射派生）
- 输出：`disease/pipelines/biomarker/reports/biomarker_edges_v1.{validation,manifest}.json`

## 状态

**scaffold-only**：派生源数据（disease-gene/protein 映射）尚未本地可用。`run.sh` 在输入缺失时 `exit 2`，不伪造表。

## 运行

```bash
bash disease/pipelines/biomarker/scripts/run.sh
```

## 契约 / 校验

见 `contracts/biomarker_edges_v1.json` 与 `specs/disease-layer-v1/contracts/biomarker.md`。结构 + backbone 校验见 `tools/biomarker_firstpass_check.py`。
