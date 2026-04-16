#!/usr/bin/env bash
set -euo pipefail

# Run from interaction/release/external/interaction-l2-v1 or pass target dir as arg
ROOT="${1:-.}"
cd "$ROOT"

# integrated
cat tables/integrated/chunks/interaction_aggregate_score_v2.tsv.zst.part.* \
  > tables/integrated/interaction_aggregate_score_v2.tsv.zst

cat tables/integrated/chunks/interaction_cross_validation_v2.tsv.zst.part.* \
  > tables/integrated/interaction_cross_validation_v2.tsv.zst

# psi
cat tables/psi/chunks/psi_activity_context_v2.tsv.zst.part.* \
  > tables/psi/psi_activity_context_v2.tsv.zst

cat tables/psi/chunks/psi_structure_evidence_v2.tsv.zst.part.* \
  > tables/psi/psi_structure_evidence_v2.tsv.zst

echo "[OK] reconstructed chunked assets under tables/"
