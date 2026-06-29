# Disease Layer v1 Public Package

Disease-centered semantic layer 的对外发布包（feature: disease-layer-v1）。

## 当前状态

**v1 partial release** — 3 张真实表已发布，3 张待源数据。

| 表 | 状态 | 行数 | 来源 |
|---|---|---|---|
| disease/disease_node_v1.tsv | ✅ RELEASED | 29,868 | MONDO v2026-06-02 |
| phenotype/phenotype_node_v1.tsv | ✅ RELEASED | 19,836 | HPO v2026-06-23 |
| disease_phenotype/disease_phenotype_v1.tsv | ✅ RELEASED | 269,518 | HPO annotations + MONDO xref bridge |
| disease_gene_target/disease_gene_target_v1.tsv | ⏳ PENDING | — | DisGeNET / Open Targets |
| disease_drug/disease_drug_v1.tsv | ⏳ PENDING | — | DrugCentral / DrugBank |
| biomarker/biomarker_edges_v1.tsv | ⏳ PENDING | — | relation-layer (Decision 015) |

## Reports

- `reports/disease_v1_build.json` — validation summary (3 tables, 319,222 rows, status=PASS)
- `reports/disease_v1_query_paths.json` — query path structural validation (SC-008)

## Manifests

- `manifests/disease_v1_manifest.json` — SHA256 + size for released tables

## 版本边界

本包归属 **Disease Layer v1**。版本指针：`disease/products/disease/current.json` → `disease-v1`。
