# Graph-RAG Protein Public Repository

这是一个**Protein-only** 的公开结果仓库。

## 快速入口
1. `products/protein/current.json`
2. `release/index.json`
3. `release/external/protein-v6/`

## 目录主线
- `release/external/protein-v6/tables/`：按 doc 分区的数据表
- `release/external/protein-v6/reports/`：验证与QA报告
- `release/external/protein-v6/manifests/`：发布清单
- `pipelines/protein*`：可复现构建与校验管线

## 校验
```bash
python3 scripts/build_release_index.py --products-root products --out release/index.json
python3 scripts/validate_release_index.py --index release/index.json --schema release/schema/index.schema.json --repo-root .
python3 scripts/check_release_consistency.py --index release/index.json --out release/consistency_report.json
pytest -q tests/release
```
