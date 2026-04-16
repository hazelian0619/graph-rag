# Graph-RAG

该仓库采用**按类目分区**的公开结构。顶层四个类目：

- `protein/`（已实装，当前主交付）
- `rna/`（骨架已建，待填充）
- `molecule/`（骨架已建，含已整理的小分子输出）
- `interaction/`（骨架已建，待填充）

## 统一目录标准（每个类目）

- `data/`：原始/处理中/输出数据
- `pipelines/`：构建与校验流水线
- `products/`：版本指针与产品元数据
- `release/`：对外发布包（external）与内部过程物（internal）

## 当前可直接使用（Protein）

1. `protein/products/protein/current.json`
2. `protein/release/index.json`
3. `protein/release/external/protein-v6/`

## Protein 校验命令

```bash
python3 scripts/build_release_index.py --products-root protein/products --out protein/release/index.json
python3 scripts/validate_release_index.py --index protein/release/index.json --schema protein/release/schema/index.schema.json --repo-root .
python3 scripts/check_release_consistency.py --index protein/release/index.json --out protein/release/consistency_report.json
pytest -q tests/release
```
