#!/usr/bin/env bash
set -euo pipefail

# Disease-drug pipeline (disease-layer-v1) — scaffold.
# Builds disease_drug_v1 edge table: indication/treatment relation ONLY (DECISIONS 006).
# No contraindication/AE/safety. Drug must map to existing molecule backbone (FR-010).
# Sources: DrugCentral / DrugBank (indication fields only). Validation + manifest only.

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
export KG_ROOT="$REPO_ROOT"
cd "$REPO_ROOT"

INPUT_TABLE="data/processed/disease_drug_v1.tsv"
CONTRACT="disease/pipelines/disease_drug/contracts/disease_drug_v1.json"
REPORT_VALIDATION="disease/pipelines/disease_drug/reports/disease_drug_v1.validation.json"
REPORT_MANIFEST="disease/pipelines/disease_drug/reports/disease_drug_v1.manifest.json"

DATA_VERSION="${DATA_VERSION:-kg-data-local}"

if [[ ! -f "$INPUT_TABLE" ]]; then
  echo "[ERROR] missing input table: $INPUT_TABLE" >&2
  echo "        Source data (DrugCentral / DrugBank indication) not yet available locally. Scaffold-only." >&2
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

echo "[DONE] Disease-drug v1 QA + manifest finished."
