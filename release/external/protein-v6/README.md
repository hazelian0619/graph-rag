# Protein v6 Public Package

该目录是 Protein 产品的**对外交付包**，按《组件细节》进行分区：

- `tables/sequence/`：序列与主实体
- `tables/structure/`：PDB/AlphaFold 结构数据
- `tables/function/`：GO/KEGG/结构域/通路/PTM 相关
- `tables/physicochemical/`：理化属性
- `tables/isoform/`：异构体层
- `tables/provenance/`：版本追溯
- `tables/network/`：蛋白互作网络基础边
- `tables/reference/`：参考映射（HGNC/UniProt seed）
- `reports/`：各子层 validation/qa/build 指标
- `manifests/`：对外交付 manifest

说明：
- 运行管线仍以 `data/processed/**` 与 `pipelines/**` 为工程路径。
- 对外消费优先使用本目录的 `tables/** + reports/** + manifests/**`。
