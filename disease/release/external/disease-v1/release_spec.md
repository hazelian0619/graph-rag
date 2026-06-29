# Disease Layer v1 Release Specification & QA Gates

Cross-cutting quality + productization spec for the disease-v1 public package (Phase 7: T031–T039). Enforces constitution Principles III/IV/V and FR-013–FR-022.

## Data-quality gates (T031–T034)

| Gate | Rule | Enforced by |
|------|------|-------------|
| Canonicalization (FR-013) | 所有核心节点主键稳定，非裸文本；MONDO/HP 正则 | `disease_node_v1.json` / `phenotype_node_v1.json` pk 规则 |
| Dedup (FR-014) | 同义词/多源 id 受控去重 | `pk_unique` 规则 + 上游 canonicalization |
| Source traceability (FR-015) | 核心节点/边可说明来源 | 每个 contract 的 `source_non_empty` 规则 |
| Semantic clarity (FR-016) | 每类边语义标注 | `relation_semantics`/`relation_type`/`evidence_class` allowed-value 规则 |
| Controlled scope (FR-017) | 不为覆盖率引入不稳定源 | 仅接 D008 源族；PrimeKG/Hetionet reference-only |

## Productization gates (T035–T039)

- **T035 (release assembly)**: 待数据就绪后组装 `tables/{disease,disease_gene_target,phenotype,disease_phenotype,disease_drug,biomarker}/*.tsv` + `manifests/*` + `reports/*`。当前为目录骨架 + release README（已就位）。
- **T036 (field definitions)**: 每个 contract JSON 定义字段口径（已就位）。
- **T037 (version boundary)**: `products/disease/{current,product}.json` 指向 `disease-v1`（已就位）。
- **T038 (reusability)**: 结构可被 downstream 直接复用（schema + contracts 已就位）。
- **T039 (style alignment)**: 组织方式与 molecule/protein/rna release 一致（`<domain>/release/external/<domain>-v1/{tables,manifests,reports}`）— 已就位。

## Run-all gate script

见 `tools/disease_release_gate.py`：逐表跑 contract validation + scope/semantic 专项 check，输出汇总 report。数据未就位时 SKIP（exit 0），不伪造。

## Query readiness (T040)

disease-centered 多跳路径（disease → phenotype/target/drug/biomarker → backbone）结构性支持由表 + 边 schema 保证。T040 在数据就位后验证连通性（复用 `disease_backbone_connectivity.py`）。

## Status

scaffold-complete for contracts/validators/layout; **data-blocked** for actual TSV/manifest/report assembly.
