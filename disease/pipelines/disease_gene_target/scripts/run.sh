#!/usr/bin/env bash
set -euo pipefail

# Disease-gene / disease-target edge pipeline (disease-layer-v1) — scaffold.
# Builds the unified disease_gene_target_v1 edge set (D-GATE-3: implementation converges,
# product semantics preserved via relation_semantics / target_type / source / evidence_class).
# Sources: DisGeNET (disease-gene curated), Open Targets (disease-target therapeutic).
# Performs contract validation + manifest only; ingestion of source data is upstream.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

INPUT_TABLE="data/processed/disease_gene_target_v1.tsv"
CONTRACT="disease/pipelines/disease_gene_target/contracts/disease_gene_target_v1.json"
REPORT_VALIDATION="disease/pipelines/disease_gene_target/reports/disease_gene_target_v1.validation.json"
REPORT_MANIFEST="disease/pipelines/disease_gene_target/reports/disease_gene_target_v1.manifest.json"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"

if [[ ! -f "$INPUT_TABLE" ]]; then
  echo "[ERROR] missing input table: $INPUT_TABLE" >&2
  echo "        Source data (DisGeNET / Open Targets) not yet available locally. Scaffold-only." >&2
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

echo "[DONE] Disease-gene/target v1 QA + manifest finished."
