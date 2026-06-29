#!/usr/bin/env python3
"""Ingest phenotype.hpoa -> disease_phenotype_v1.tsv (disease-layer-v1, US2 edges).

Bridges hpoa disease keys (OMIM:/ORPHA:) to MONDO via MONDO exactMatch xrefs.
HPO-only phenotype (DECISIONS 004). Evidence class derived from hpoa evidence column.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
MONDO = REPO / "data/raw/disease/mondo-base.json"
HPOA = REPO / "data/raw/disease/phenotype.hpoa"
OUT = REPO / "data/processed/disease_phenotype_v1.tsv"

# hpoa evidence -> our evidence_class
EVIDENCE_MAP = {
    "IEA": "inferred_mapping",
    "PCS": "curated_association",
    "TAS": "curated_association",
}


def build_xref_map() -> dict:
    """MONDO id -> {OMIM:..., ORPHA:...} and reverse lookup."""
    with MONDO.open() as f:
        d = json.load(f)
    rev = {}  # "OMIM:12345" -> "MONDO:..."
    for n in d["graphs"][0]["nodes"]:
        nid = n.get("id", "")
        if "MONDO_" not in nid:
            continue
        if n.get("meta", {}).get("deprecated"):
            continue
        mondo_id = "MONDO:" + nid.split("MONDO_")[1]
        for bp in n.get("meta", {}).get("basicPropertyValues", []):
            pred = bp.get("pred", "")
            if "exactMatch" not in pred:
                continue
            val = bp.get("val", "")
            if "omim.org" in val and "/" in val:
                key = "OMIM:" + val.rstrip("/").split("/")[-1]
                rev.setdefault(key, mondo_id)
            elif "orpha.net" in val:
                key = "ORPHA:" + val.split("Orphanet_")[-1] if "Orphanet_" in val else ""
                if key:
                    rev.setdefault(key, mondo_id)
    return rev


def main() -> int:
    rev = build_xref_map()
    print("[info] xref map: " + str(len(rev)) + " OMIM/ORPHA -> MONDO")

    edges = {}
    unmapped = 0
    total = 0
    with HPOA.open(encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            break
        # line now holds header; rewind handled by DictReader below
    with HPOA.open(encoding="utf-8") as f:
        reader = csv.DictReader(
            (ln for ln in f if not ln.startswith("#")), delimiter="\t"
        )
        for row in reader:
            db_id = row.get("database_id", "").strip()
            hpo_id = row.get("hpo_id", "").strip()
            if not db_id or not hpo_id:
                continue
            total += 1
            mondo = rev.get(db_id)
            if not mondo:
                unmapped += 1
                continue
            evidence = row.get("evidence", "").strip()
            ev_class = EVIDENCE_MAP.get(evidence, "curated_association")
            key = (mondo, hpo_id, ev_class)
            if key not in edges:
                edges[key] = {
                    "disease_id": mondo,
                    "hpo_id": hpo_id,
                    "source": "HPO_annotations",
                    "evidence_class": ev_class,
                }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["disease_id", "hpo_id", "source", "evidence_class"]
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in edges.values():
            w.writerow(r)

    mapped = total - unmapped
    print("[OK] disease_phenotype_v1: " + str(len(edges)) + " edges (" + str(mapped) + "/" + str(total) + " mapped; " + str(unmapped) + " unmapped) -> " + str(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
