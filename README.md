# Graph-RAG Protein Public Repository

这是一个 **Protein-only** 的公开结果仓库。

## 快速入口
1. `protein/products/protein/current.json`
2. `protein/release/index.json`
3. `protein/release/external/protein-v6/`

## 目录主线
- `protein/release/external/protein-v6/tables/`：按 doc 分区的数据表
- `protein/release/external/protein-v6/reports/`：公开验证与 QA 报告
- `protein/release/external/protein-v6/manifests/`：发布清单
- `protein/pipelines/protein*`：可复现构建与校验管线

## 校验
```bash
python3 scripts/build_release_index.py --products-root protein/products --out protein/release/index.json
python3 scripts/validate_release_index.py --index protein/release/index.json --schema protein/release/schema/index.schema.json --repo-root .
python3 scripts/check_release_consistency.py --index protein/release/index.json --out protein/release/consistency_report.json
pytest -q tests/release
```
