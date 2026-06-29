# Disease 节点 canonicalization 管线（`disease_node_v1`）

本管线负责 disease 主表 canonicalization：MONDO-first 主键、别名 / 多源 id 去重收敛（FR-001 / FR-013 / FR-014）。不重建 MONDO 本身，只做 contract validation + manifest。

## 输入 / 输出

- 输入：`data/processed/disease_node_v1.tsv`（上游已 ingest 的 MONDO disease 表；原始下载不入库，见 REPO_STANDARDS）
- 输出：
  - `disease/pipelines/disease/reports/disease_node_v1.validation.json`
  - `disease/pipelines/disease/reports/disease_node_v1.manifest.json`

## 状态

**scaffold-only**：源数据（MONDO）尚未本地可用。`run.sh` 在输入缺失时 `exit 2`（与 protein run.sh 一致），不伪造表。

## 字段 / 契约

见 `contracts/disease_node_v1.json` 与 `specs/disease-layer-v1/contracts/disease-node.md`。Primary key = MONDO id（D-GATE-1，default-applied）。

## 运行

```bash
bash disease/pipelines/disease/scripts/run.sh
```
