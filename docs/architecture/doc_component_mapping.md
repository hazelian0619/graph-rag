# 《组件细节》到仓库资产映射（主线版）

## 1) 三类实体

### Protein Entity
- 序列/ID/基因/物种：`pipelines/protein/`, `pipelines/protein_provenance/`
- 结构（PDB/AlphaFold）：`pipelines/protein_pdb/`, `data/processed/alphafold_quality.tsv`
- 功能注释（GO/KEGG/域）：`pipelines/protein_kegg/`, `pipelines/protein_domains/`
- 理化属性：`pipelines/protein_physchem/`, `pipelines/protein_physchem_extended/`
- 主产物：`data/processed/protein_master_v6_clean.tsv`

### RNA Entity
- 主表/类型/映射：`pipelines/rna/`, `pipelines/rna_external_xref/`, `pipelines/rna_id_canonical/`
- 结构（Rfam/PDB/预测/协方差）：`pipelines/rna_rfam_structure/`, `pipelines/rna_pdb/`, `pipelines/rna_structure_models/`
- 功能与证据：`pipelines/rna_expression_rbp/`, `pipelines/rna_type_features/`
- 发布资产：`release/external/rna-l1l2-v2/`

### Small Molecule Entity
- 结构标识：`pipelines/molecule_structure_identifiers/`
- 理化属性：`pipelines/molecule_physchem_descriptors/`
- 生物活性：`pipelines/molecule_activity_fusion/`
- PK/Tox：`pipelines/molecule_pk_tox_v2/`
- 三维结构：`pipelines/molecule_3d_registry/`, `pipelines/molecule_3d_experimental_linker/`
- 发布资产：`release/external/molecule-l1l2-v2/`

## 2) 分子间相互作用

- PPI：`pipelines/edges_ppi/`, `pipelines/ppi_semantic_enrichment/`
- PSI：`pipelines/psi_activity_structure_enrichment/`, `pipelines/psi_condition_enrichment/`
- RPI：`pipelines/rna_rpi/`, `pipelines/rpi_site_domain_enrichment/`
- 跨库一致性与本体：`pipelines/interaction_cross_validation/`, `pipelines/interaction_ontology_mapping/`

## 3) 交付与质量门禁

- 产品指针：`products/*/current.json`
- 统一索引：`release/index.json`
- QA/validation：各 pipeline `reports/*.json`
