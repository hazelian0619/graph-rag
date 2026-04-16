# Protein v6 Public Package

这是对外公开消费的 Protein 数据包（doc 对齐版）。

## 分区
- `tables/sequence/`：序列与主实体
- `tables/structure/`：PDB / AlphaFold
- `tables/function/`：KEGG / 结构域 / 通路 / PTM
- `tables/physicochemical/`：理化属性
- `tables/isoform/`：异构体层
- `tables/provenance/`：来源版本追溯
- `tables/network/`：蛋白互作基础边
- `tables/reference/`：参考映射（HGNC / UniProt seed）

## 质量与清单
- `reports/`：仅保留对外核心 validation/qa 报告
- `manifests/protein_public_manifest_v1.json`：公开包资产清单

## 说明
- 过程性构建产物（build/smoke/sample/metrics/audit）已移到：
  `release/internal/protein-v6-process/`
- 对外下载优先使用 `products/protein/current.json` 指向的路径。
