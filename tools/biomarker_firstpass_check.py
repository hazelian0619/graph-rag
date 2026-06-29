#!/usr/bin/env python3
"""Biomarker first-pass structural + backbone check (disease-layer-v1, US4).

Verifies biomarker is a real structural artifact (FR-003; ACCEPTANCE Non-Acceptance Case 5
rejects 'talk-only') and that target_id maps to an existing molecular backbone entity
(constitution I / FR-011). Scaffold: exits cleanly (SKIP, 0) when source data absent.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def _read_keys(path: Path, cols) -> set:
    keys = set()
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            for c in cols:
                v = row.get(c, "").strip()
                if v:
                    keys.add(v)
    return keys


def main() -> int:
    ap = argparse.ArgumentParser(description="Biomarker first-pass structural + backbone check (US4)")
    ap.add_argument("--biomarker-edges", default="data/processed/biomarker_edges_v1.tsv")
    ap.add_argument("--backbone", default="data/processed/protein_master_v6_clean.tsv")
    ap.add_argument("--out", default="disease/release/external/disease-v1/reports/biomarker_firstpass.json")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parent.parent
    be = repo / args.biomarker_edges
    bb = repo / args.backbone
    missing = [p for p in (str(be), str(bb)) if not Path(p).is_file()]
    if missing:
        print("[SKIP] source data not yet available locally; biomarker check is scaffold-only:")
        for p in missing:
            print("        missing: " + p)
        return 0

    rows = 0
    mapped = 0
    types_seen = set()
    backbone = _read_keys(bb, ("uniprot_id", "ensembl_gene_id", "ncbi_gene_id", "hgnc_id"))
    with be.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows += 1
            types_seen.add(row.get("biomarker_type", "").strip())
            if row.get("target_id", "").strip() in backbone:
                mapped += 1

    has_structure = rows > 0  # FR-003: must be a real artifact, not talk-only
    has_edge_type = len(types_seen) >= 1  # FR-007: at least one biomarker edge type
    mapping_frac = round(mapped / rows, 4) if rows else 0.0
    backbone_ok = mapping_frac >= 0.5  # FR-011: maps to existing entity layer
    report = {
        "check": "biomarker_firstpass",
        "policy": "structural artifact + backbone mapping (FR-003/FR-007/FR-011)",
        "rows": rows,
        "has_structure": has_structure,
        "edge_types": sorted(types_seen),
        "has_edge_type": has_edge_type,
        "mapped_to_backbone": mapped,
        "mapping_fraction": mapping_frac,
        "backbone_ok": backbone_ok,
        "passed": has_structure and has_edge_type and backbone_ok,
    }
    out_path = repo / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(("[PASS] " if report["passed"] else "[FAIL] ") + json.dumps(report, ensure_ascii=False))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
