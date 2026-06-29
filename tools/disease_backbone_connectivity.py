#!/usr/bin/env python3
"""Disease -> backbone connectivity check (disease-layer-v1, US1 validation skeleton).

Verifies that disease-layer nodes are not isolated from the molecular backbone
(constitution I: Backbone-First; FR-008 / FR-012; SC-006). This is a SKELETON:
it enforces the connectivity contract once release tables exist. When source data
is not yet available locally it exits cleanly with a clear message (no fake tables).

Checks (when data present):
  1. Every disease_id referenced in disease_gene_target edges resolves to disease_node.
  2. Edges' target_id can be mapped to the existing protein / gene backbone.
  3. The fraction of backbone-unmappable disease nodes stays below the isolation threshold.

Inputs (all optional at scaffold stage):
  --disease-node data/processed/disease_node_v1.tsv
  --edges       data/processed/disease_gene_target_v1.tsv
  --backbone    data/processed/protein_master_v6_clean.tsv (protein/gene backbone key set)
  --threshold   max isolated-disease fraction (default 0.05)
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path


def _read_keys(path: Path, key_col: str) -> set:
    keys = set()
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            v = row.get(key_col, "").strip()
            if v:
                keys.add(v)
    return keys


def main() -> int:
    ap = argparse.ArgumentParser(description="Disease->backbone connectivity check (US1)")
    ap.add_argument("--disease-node", default="data/processed/disease_node_v1.tsv")
    ap.add_argument("--edges", default="data/processed/disease_gene_target_v1.tsv")
    ap.add_argument("--backbone", default="data/processed/protein_master_v6_clean.tsv")
    ap.add_argument("--threshold", type=float, default=0.05,
                    help="max isolated-disease fraction")
    ap.add_argument("--out", default="disease/release/external/disease-v1/reports/disease_backbone_connectivity.json")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parent.parent

    missing = [p for p in (args.disease_node, args.edges, args.backbone) if not (repo / p).is_file()]
    if missing:
        print("[SKIP] source data not yet available locally; connectivity check is scaffold-only:")
        for p in missing:
            print("        missing: " + p)
        return 0

    node_path = repo / args.disease_node
    edge_path = repo / args.edges
    backbone_path = repo / args.backbone

    diseases = _read_keys(node_path, "disease_id")

    # backbone keys: union of uniprot + gene ids present in the protein/gene master
    backbone = set()
    with backbone_path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            for col in ("uniprot_id", "ensembl_gene_id", "ncbi_gene_id", "hgnc_id"):
                v = row.get(col, "").strip()
                if v:
                    backbone.add(v)

    mapped_diseases = set()
    fk_ok = True
    edge_count = 0
    with edge_path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            edge_count += 1
            did = row.get("disease_id", "").strip()
            tid = row.get("target_id", "").strip()
            if did in diseases:
                if tid and tid in backbone:
                    mapped_diseases.add(did)
            else:
                fk_ok = False

    total = len(diseases) or 1
    isolated = diseases - mapped_diseases
    isolated_frac = round(len(isolated) / total, 4)
    passed = fk_ok and isolated_frac <= args.threshold and edge_count > 0

    report = {
        "check": "disease_backbone_connectivity",
        "disease_nodes": len(diseases),
        "edges": edge_count,
        "mapped_diseases": len(mapped_diseases),
        "isolated_diseases": len(isolated),
        "isolated_fraction": isolated_frac,
        "threshold": args.threshold,
        "fk_integrity": fk_ok,
        "passed": passed,
    }
    out_path = repo / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(("[PASS] " if passed else "[FAIL] ") + json.dumps(report, ensure_ascii=False))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
