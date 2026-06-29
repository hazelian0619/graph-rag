#!/usr/bin/env bash
set -euo pipefail

# Disease canonical node pipeline (disease-layer-v1) — scaffold.
# Canonicalizes disease nodes: MONDO-first PK (D-GATE-1), alias/multi-source dedup (FR-013/FR-014).
# Input is expected to be a processed MONDO disease table; source data is ingested upstream
# (not part of this repo per REPO_STANDARDS). This script performs contract validation + manifest.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

INPUT_TABLE="data/processed/disease_node_v1.tsv"
CONTRACT="disease/pipelines/disease/contracts/disease_node_v1.json"
REPORT_VALIDATION="disease/pipelines/disease/reports/disease_node_v1.validation.json"
REPORT_MANIFEST="disease/pipelines/disease/reports/disease_node_v1.manifest.json"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"

if [[ ! -f "$INPUT_TABLE" ]]; then
  echo "[ERROR] missing input table: $INPUT_TABLE" >&2
  echo "        Source data (MONDO) not yet available locally. Scaffold-only." >&2
  exit 2
fi

python3 tools/kg_validate_table.py \
  --contract "$CONTRACT" \
  --table "$INPUT_TABLE" \
  --out "$REPORT_VALIDATION"

python3 tools/kg_make_manifest.py \
  --data-version "$DATA_VERSION" \
  --out "$REPORT_MANIFEST" \
  "$INPUT_TABLE"

echo "[DONE] Disease node v1 QA + manifest finished."
