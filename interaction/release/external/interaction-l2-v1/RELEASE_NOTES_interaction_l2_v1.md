# Interaction L2 v1

Interaction-only dataset release aligned with PR #8.

## Included datasets
Direct `.tsv.zst` assets:
- ppi_method_context_v2.tsv.zst
- ppi_function_context_v2.tsv.zst
- rpi_site_context_v2.tsv.zst
- rpi_domain_context_v2.tsv.zst
- rpi_function_context_v2.tsv.zst
- interaction_cross_validation_v2.tsv.zst
- interaction_ontology_mapping_v2.tsv.zst

Chunked assets (due upload channel stability):
- psi_structure_evidence_v2.tsv.zst.part.000..002
- interaction_aggregate_score_v2.tsv.zst.part.000..003
- psi_activity_context_v2.tsv.zst.part.000..009

Reassemble examples:
```bash
cat psi_structure_evidence_v2.tsv.zst.part.* > psi_structure_evidence_v2.tsv.zst
cat interaction_aggregate_score_v2.tsv.zst.part.* > interaction_aggregate_score_v2.tsv.zst
cat psi_activity_context_v2.tsv.zst.part.* > psi_activity_context_v2.tsv.zst
```

## QA status (all PASS)
- A PPI semantic enrichment QA
- B PSI activity/structure enrichment QA
- C RPI site/domain/function gates
- D Cross-validation + aggregate score QA
- E Ontology mapping gates
- interaction_release_local package

See `manifest_interaction_l2_v1.json` for checksums.
