#!/usr/bin/env python3
"""HPO-only phenotype enforcement (disease-layer-v1, US2 validation skeleton).

Reverse-checks DECISIONS 004 / constitution Binding Decision: the phenotype layer
MUST contain only HPO terms. Rejects any phenotype node or disease-phenotype edge whose
hpo_id does not match ^HP:\\d{7}$, i.e. no symptom / sign / diagnosis multi-system.

Scaffold: exits cleanly (SKIP, 0) when source data not yet local; no fake tables.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

HPO_RE = re.compile(r"^HP:\d{7}$")


def _check_file(path: Path, key_col: str) -> dict:
    rows = 0
    bad = 0
    bad_samples = []
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            rows += 1
            v = row.get(key_col, "").strip()
            if not HPO_RE.match(v):
                bad += 1
                if len(bad_samples) < 3:
                    bad_samples.append(v)
    return {
        "file": str(path),
        "rows": rows,
        "non_hpo": bad,
        "hpo_only": bad == 0 and rows > 0,
        "samples": bad_samples,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description="HPO-only phenotype enforcement (US2)")
    ap.add_argument("--phenotype-node", default="data/processed/phenotype_node_v1.tsv")
    ap.add_argument("--edges", default="data/processed/disease_phenotype_v1.tsv")
    ap.add_argument("--out", default="disease/release/external/disease-v1/reports/disease_phenotype_hpo_only.json")
    args = ap.parse_args()

    repo = Path(__file__).resolve().parent.parent
    targets = [
        (repo / args.phenotype_node, "hpo_id"),
        (repo / args.edges, "hpo_id"),
    ]
    missing = [str(p) for p, _ in targets if not p.is_file()]
    if missing:
        print("[SKIP] source data not yet available locally; HPO-only check is scaffold-only:")
        for p in missing:
            print("        missing: " + p)
        return 0

    results = [_check_file(p, col) for p, col in targets]
    passed = all(r["hpo_only"] for r in results)
    report = {
        "check": "disease_phenotype_hpo_only",
        "policy": "HPO-only (DECISIONS 004); symptom/sign/diagnosis multi-system rejected",
        "files": results,
        "passed": passed,
    }
    out_path = repo / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    print(("[PASS] " if passed else "[FAIL] ") + json.dumps(report, ensure_ascii=False))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
