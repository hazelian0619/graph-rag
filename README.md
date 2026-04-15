# Graph-RAG Bio-KG Foundation (Public Result Repo)

面向对外交付的生物知识图谱结果仓库，聚焦 4 条主线：

- Protein Entity
- RNA Entity
- Small Molecule Entity
- Cross-entity Interactions (PPI / PSI / RPI)

本仓库是**结果导向**：优先提供可下载数据、可验证清单、可追溯 manifest；不保留大量过程性协作文档。

## 快速入口

1. 统一产品索引：`release/index.json`
2. 产品当前版本指针：`products/*/current.json`
3. 下载工具：`scripts/download_dataset.py`

## 对照《组件细节》架构

- 资源层：`products/`, `release/`
- 数据层：`data/`, `pipelines/`
- 算法层：`pipelines/ppi_*`, `pipelines/psi_*`, `pipelines/rpi_*`, `pipelines/interaction_*`
- 应用层（对外接口入口准备）：`scripts/` + release metadata

详细映射见：`docs/architecture/doc_component_mapping.md`

## 数据分发策略

- **Protein**：仓库快照（`data/processed/protein_master_v6_clean.tsv`）
- **RNA / Molecule / Interaction**：以 Release 为主，仓库内保留 manifest / QA / checksum
- 大文件不直接入 git（>100MB 一律走 Release 资产）


> 过渡说明：当前 RNA/Molecule/Interaction 的 `release_url` 暂指向既有稳定发布源 `hazelian0619/protian-entity`，后续会在 `graph-rag` 建立同版本 tag 后切换。

## 基本校验

```bash
python3 scripts/validate_release_index.py \
  --index release/index.json \
  --schema release/schema/index.schema.json \
  --repo-root .

python3 scripts/check_release_consistency.py \
  --index release/index.json \
  --out release/consistency_report.json
```
