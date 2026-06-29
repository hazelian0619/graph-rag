#!/usr/bin/env python3
"""Ingest MONDO base.json -> disease_node_v1.tsv (disease-layer-v1, US1 node layer).

Real ingest (not scaffold). D-GATE-1: primary key = MONDO id directly.
Filters out deprecated/obsolete nodes. Extracts canonical_name (lbl), synonyms,
DOID xref (supplementary, DECISIONS 003). Source = MONDO.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
RAW = REPO / "data/raw/disease/mondo-base.json"
OUT = REPO / "data/processed/disease_node_v1.tsv"
FETCH_DATE = "2026-06-29"


def main() -> int:
    with RAW.open() as f:
        d = json.load(f)
    g = d["graphs"][0]
    nodes = g["nodes"]

    rows = []
    for n in nodes:
        nid = n.get("id", "")
        if "MONDO_" not in nid:
            continue
        meta = n.get("meta", {})
        if meta.get("deprecated"):
            continue
        lbl = n.get("lbl")
        if not lbl:
            continue
        mondo_id = "MONDO:" + nid.split("MONDO_")[1]

        synonyms = []
        for syn in meta.get("synonyms", []):
            if syn.get("pred") == "hasExactSynonym":
                v = syn.get("val", "").strip()
                if v:
                    synonyms.append(v)

        doid = ""
        external_ids = {}
        for bp in meta.get("basicPropertyValues", []):
            pred = bp.get("pred", "")
            val = bp.get("val", "")
            if val.startswith("http://purl.obolibrary.org/obo/DOID_"):
                doid = "DOID:" + val.split("DOID_")[1]
            elif "hasDbXref" in pred or "xref" in pred.lower():
                if val.startswith("http://purl.obolibrary.org/obo/"):
                    continue
                ext_src = val.split(":")[0] if ":" in val else ""
                if ext_src and ext_src not in ("http", "https"):
                    external_ids.setdefault(ext_src, val)

        rows.append({
            "disease_id": mondo_id,
            "canonical_name": lbl,
            "aliases": "|".join(sorted(set(synonyms))),
            "source": "MONDO",
            "xref_doid": doid,
            "external_ids": json.dumps(external_ids, ensure_ascii=False) if external_ids else "",
            "fetch_date": FETCH_DATE,
        })

    seen = {}
    for r in rows:
        seen.setdefault(r["disease_id"], r)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["disease_id", "canonical_name", "aliases", "source", "xref_doid", "external_ids", "fetch_date"]
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in sorted(seen.values(), key=lambda x: x["disease_id"]):
            w.writerow(r)

    print("[OK] disease_node_v1: " + str(len(seen)) + " disease nodes -> " + str(OUT))
    return 0


if __name__ == "__main__":
    sys.exit(main())
