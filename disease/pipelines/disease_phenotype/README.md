# Disease-phenotype 管线（`phenotype_node_v1` + `disease_phenotype_v1`）

本管线在一个目录内产出两层（T015 节点 + T016 边，按 tasks.md「同管线内或子管线」）：

- **phenotype_node_v1**：HPO-only canonical phenotype 节点（FR-002）
- **disease_phenotype_v1**：disease (MONDO) ↔ phenotype (HPO) 边（FR-005 / FR-009）

## 策略（DECISIONS 004 / constitution）

Phenotype 层 **HPO-only**，不引入 symptom / sign / diagnosis 多体系。T017 校验会拒绝任何不符合 `^HP:\d{7}$` 的 phenotype。

## 输入 / 输出

- 输入：`data/processed/phenotype_node_v1.tsv`、`data/processed/disease_phenotype_v1.tsv`（上游 ingest 的 HPO + HPO annotations）
- 输出：`disease/pipelines/disease_phenotype/reports/{phenotype_node_v1,disease_phenotype_v1}.{validation,manifest}.json`

## 状态

**scaffold-only**：源数据（HPO / HPO annotations）尚未本地可用。`run.sh` 在输入缺失时 `exit 2`，不伪造表。

## 运行

```bash
bash disease/pipelines/disease_phenotype/scripts/run.sh
```

## 契约

见 `contracts/` 下两份 JSON 与 `specs/disease-layer-v1/contracts/disease-phenotype.md`。
