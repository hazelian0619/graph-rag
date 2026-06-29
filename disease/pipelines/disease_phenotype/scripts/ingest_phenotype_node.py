#!/usr/bin/env python3
"""Ingest HPO hp.json -> phenotype_node_v1.tsv (disease-layer-v1, US2 node layer).

HPO-only (DECISIONS 004). Primary key = HP id. Filters obsolete nodes.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
RAW = REPO / "data/raw/disease/hp.json"
OUT = REPO / "data/processed/phenotype_node_v1.tsv"
FETCH_DATE = "2026-06-29"


def main() -> int:
    with RAW.open() as f:
        d = json.load(f)
    g = d["graphs"][0]
    nodes = g["nodes"]

    rows = {}
    for n in nodes:
        nid = n.get("id", "")
        if "HP_" not in nid:
            continue
        meta = n.get("meta", {})
        if meta.get("deprecated"):
            continue
        lbl = n.get("lbl")
        if not lbl:
            continue
        hp_id = "HP:" + nid.split("HP_")[1]

        synonyms = []
        for syn in meta.get("synonyms", []):
            if syn.get("pred") in ("hasExactSynonym", "hasRelatedSynonym"):
                v = syn.get("val", "").strip()
                if v:
                    synonyms.append(v)

        hierarchy_ref = ""
        for bp in meta.get("basicPropertyValues", []):
            pred = bp.get("pred", "")
            val = bp.get("val", "")
            if ("subClassOf" in pred or "is_a" in pred) and "HP_" in val:
                hierarchy_ref = "HP:" + val.split("HP_")[1]
                break

        rows[hp_id] = {
            "hpo_id": hp_id,
            "term": lbl,
            "synonyms": "|".join(sorted(set(synonyms))),
            "source": "HPO",
            "hierarchy_reference": hierarchy_ref,
            "fetch_date": FETCH_DATE,
        }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["hpo_id", "term", "synonyms", "source", "hierarchy_reference", "fetch_date"]
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sorted(rows.values(), key=lambda x: x["hpo_id"]):
            w.writerow(r)

    print("[OK] phenotype_node_v1: " + str(len(rows)) + " HPO terms -> " + str(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
