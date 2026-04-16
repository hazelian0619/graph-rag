# Graph-RAG Bio-KG Foundation

面向对外的数据结果仓库，聚焦 4 条主线：Protein / RNA / Molecule / Interaction。

## 快速入口
1. `release/index.json`（统一产品索引）
2. `products/*/current.json`（各产品当前版本）
3. `scripts/download_dataset.py`（下载工具）

## 仓库定位
- 结果导向：优先提供可下载数据与可验证元数据
- 公开专业：不收录过程性协作文档
- 稳定入口：四产品统一版本契约

详细标准见：`docs/REPO_STANDARDS.md`

## 数据分发策略
- Protein：仓库快照为主（当前含核心结果表）
- RNA / Molecule / Interaction：Release 资产为主，仓库保留 manifest/checksum/QA

> 过渡说明：当前 RNA/Molecule/Interaction 的 `release_url` 暂指向既有稳定发布源 `hazelian0619/protian-entity`，后续会迁移为 `graph-rag` 自有 tag。

## 基本校验
```bash
python3 scripts/validate_release_index.py \
  --index release/index.json \
  --schema release/schema/index.schema.json \
  --repo-root .

python3 scripts/check_release_consistency.py \
  --index release/index.json \
  --out release/consistency_report.json

pytest -q tests/release
```
