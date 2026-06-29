# Disease Domain

Disease 目录是 Graph-RAG 在现有 protein / rna / molecule / interaction 主干之上新增的 **disease-centered semantic layer**（feature: `disease-layer-v1`）。

本域不是平行建设的第二套图谱，而是建立在分子主干之上的语义扩展层：disease、phenotype、biomarker 与核心 disease relations 成为正式结构层，并映射回现有主干。

## 目录结构

- `data/`：本地源数据占位（原始下载数据不收录，见 REPO_STANDARDS）
- `pipelines/`：disease 相关构建与校验流水线（disease 节点、disease-phenotype、disease-gene/target、disease-drug、biomarker）
- `products/disease/current.json`：当前版本指针
- `release/external/disease-v1/`：对外发布包（结果性产物：tables / manifests / reports）
- `release/internal/`：内部产物
- `release/schema/`：disease 域 schema（disease 节点、phenotype、biomarker）

## 范围（disease-layer-v1）

In scope：disease 节点（MONDO-first）、phenotype（HPO-only）、biomarker first-pass、核心 disease relations（disease-gene/target、disease-phenotype、disease-drug indication/treatment）、与现有 backbone 的映射。

Out of scope（本域不构建）：anatomy / tissue / organ 主模块、contraindication / adverse event / safety 语义、全量临床术语 / symptom / diagnosis 桥接、临床 NLP 抽取、前端 / 问答产品、PrimeKG / Hetionet 整库接入（仅作 reference）。

## 版本边界

当前版本：`disease-v1`。本域所有正式产物归属 Disease Layer v1，并通过 `products/disease/current.json` 指向 `release/external/disease-v1/`。

## 与现有主干的关系

所有新增节点与边须能映射回 protein / rna / molecule / interaction 主干，不形成大规模孤立结构（constitution Backbone-First）。

## 决策对齐

本域实现继承 DECISIONS 001–012 与 constitution v1.0.0：MONDO-first、HPO-only、biomarker first-pass、disease-drug 仅 indication/treatment、PrimeKG/Hetionet reference-only、disease-gene vs disease-target 语义须区分。
