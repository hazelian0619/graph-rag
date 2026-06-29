#!/usr/bin/env bash
set -euo pipefail

# Disease-phenotype pipeline (disease-layer-v1) — scaffold.
# Two outputs in one pipeline per task spec (T015 node + T016 edge):
#   1. phenotype_node_v1   (HPO-only canonical phenotype node; FR-002)
#   2. disease_phenotype_v1 (disease-MONDO <-> phenotype-HPO edge; FR-005/FR-009)
# HPO-only: no symptom/sign/diagnosis multi-system (DECISIONS 004 / constitution).
# Source data ingested upstream; this script does contract validation + manifest only.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"

run_table () {
  local NAME="$1" INPUT_TABLE="$2" CONTRACT="$3"
  if [[ ! -f "$INPUT_TABLE" ]]; then
    echo "[ERROR] missing input table: $INPUT_TABLE" >&2
    echo "        Source data (HPO / HPO annotations) not yet available locally. Scaffold-only." >&2
    return 2
  fi
  python3 tools/kg_validate_table.py \
    --contract "$CONTRACT" \
    --table "$INPUT_TABLE" \
    --out "disease/pipelines/disease_phenotype/reports/${NAME}.validation.json"
  python3 tools/kg_make_manifest.py \
    --data-version "$DATA_VERSION" \
    --out "disease/pipelines/disease_phenotype/reports/${NAME}.manifest.json" \
    "$INPUT_TABLE"
}

# T015 — phenotype node
run_table "phenotype_node_v1" \
  "data/processed/phenotype_node_v1.tsv" \
  "disease/pipelines/disease_phenotype/contracts/phenotype_node_v1.json"

# T016 — disease-phenotype edge
run_table "disease_phenotype_v1" \
  "data/processed/disease_phenotype_v1.tsv" \
  "disease/pipelines/disease_phenotype/contracts/disease_phenotype_v1.json"

echo "[DONE] Disease-phenotype v1 (node + edge) QA + manifest finished."
