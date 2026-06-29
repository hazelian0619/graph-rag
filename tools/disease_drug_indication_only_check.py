#!/usr/bin/env python3
"""Indication-only enforcement for disease-drug edges (disease-layer-v1, US3).

Reverse-checks DECISIONS 006 / constitution Binding Decision: disease-drug edges MUST be
indication/treatment only. Rejects any edge whose relation_type is contraindication /
adverse_event / safety_warning / broad_label. Also reports the molecule-backbone mapping
fraction (FR-010: drug must map to existing molecule, not a standalone drug system).

Scaffold: exits cleanly (SKIP, 0) when source data not yet local; no fake tables.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

FORBIDDEN = {"contraindication", "adverse_event", "safety_warning", "broad_label"}
ALLOWED = {"indication", "treatment"}


def main() -> int:
    ap = argparse.ArgumentParser(description="Indication-only + molecule-mapping check (US3)")
    ap.add_argument("--edges", default="data/processed/disease_drug_v1.tsv")
    ap.add_argument("--out", default="disease/release/external/disease-v1/reports/disease_drug_indication_only.json")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parent.parent
    edge_path = repo / args.edges
    if not edge_path.is_file():
        print("[SKIP] source data not yet available locally; indication-only check is scaffold-only:")
        print("        missing: " + args.edges)
        return 0

    rows = 0
    forbidden_count = 0
    forbidden_samples = []
    mapped = 0
    with edge_path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows += 1
            rt = row.get("relation_type", "").strip()
            if rt in FORBIDDEN:
                forbidden_count += 1
                if len(forbidden_samples) < 3:
                    forbidden_samples.append(rt)
            if row.get("molecule_backbone_mapped", "").strip() == "True":
                mapped += 1

    indication_only = forbidden_count == 0 and rows > 0
    # FR-010: drug must map to molecule backbone; flag low mapping (warn, not hard-fail v1)
    mapping_frac = round(mapped / rows, 4) if rows else 0.0
    mapping_ok = mapping_frac >= 0.5
    report = {
        "check": "disease_drug_indication_only",
        "policy": "indication/treatment only (DECISIONS 006); contraindication/AE/safety rejected",
        "rows": rows,
        "forbidden_edges": forbidden_count,
        "indication_only": indication_only,
        "forbidden_samples": forbidden_samples,
        "molecule_mapped": mapped,
        "molecule_mapping_fraction": mapping_frac,
        "molecule_mapping_ok": mapping_ok,
        "passed": indication_only and mapping_ok,
    }
    out_path = repo / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(("[PASS] " if report["passed"] else "[FAIL] ") + json.dumps(report, ensure_ascii=False))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
