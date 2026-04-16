# Interaction L2 v1 Public Package

结果性公开包，按 PPI / PSI / RPI / Integrated 分区：

- `tables/ppi`：蛋白互作主边与上下文证据
- `tables/psi`：蛋白-小分子互作（含 activity/structure 分块文件）
- `tables/rpi`：RNA-蛋白互作与上下文
- `tables/integrated`：本体映射与跨关系整合（含分块）
- `reports`：构建、QA、gates、validation 等 JSON 报告
- `manifests`：发布清单与 doc 对照自查

## 分块重组
```bash
cat interaction_aggregate_score_v2.tsv.zst.part.* > interaction_aggregate_score_v2.tsv.zst
cat interaction_cross_validation_v2.tsv.zst.part.* > interaction_cross_validation_v2.tsv.zst
cat psi_activity_context_v2.tsv.zst.part.* > psi_activity_context_v2.tsv.zst
cat psi_structure_evidence_v2.tsv.zst.part.* > psi_structure_evidence_v2.tsv.zst
```
