# Graph-RAG Public Bootstrap Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 从本地既有成果中提炼 4 条主线（protein/rna/molecule/interaction）并落地到 graph-rag 对外仓库。

**Architecture:** 复用已验证 pipeline 与 release metadata，采用“仓库中小体量结果 + Release 大资产”双轨分发；以 `release/index.json` 作为统一入口。

**Tech Stack:** Git + Python scripts + JSON manifest/QA reports

---

- [x] 建立空仓并导入主线代码与结果资产（从 protian-entity 选择性复制）
- [x] 对齐《组件细节》层次映射（docs/architecture/doc_component_mapping.md）
- [x] 增加发布边界规范（docs/governance/public_boundary.md）
- [x] 对 RNA/Molecule 使用 external release metadata 承接大文件
- [x] 回写 products/release 索引路径并执行一致性校验
