#!/usr/bin/env bash
set -euo pipefail

# Biomarker first-pass pipeline (disease-layer-v1) — scaffold.
# Form: relation layer (D-GATE-2). Source: derived from disease-gene/protein high-value
# mappings (D-GATE-4), tagged bridge_only. Validation + manifest only.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

INPUT_TABLE="data/processed/biomarker_edges_v1.tsv"
CONTRACT="disease/pipelines/biomarker/contracts/biomarker_edges_v1.json"
REPORT_VALIDATION="disease/pipelines/biomarker/reports/biomarker_edges_v1.validation.json"
REPORT_MANIFEST="disease/pipelines/biomarker/reports/biomarker_edges_v1.manifest.json"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"

if [[ ! -f "$INPUT_TABLE" ]]; then
  echo "[ERROR] missing input table: $INPUT_TABLE" >&2
  echo "        Derived source (disease-gene/protein mappings) not yet available locally. Scaffold-only." >&2
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

echo "[DONE] Biomarker first-pass v1 QA + manifest finished."
