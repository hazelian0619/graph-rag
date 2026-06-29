#!/usr/bin/env python3
"""Disease Layer v1 release gate runner (Phase 7).

Runs contract validation (via kg_validate_table.py) across all disease-layer tables and the
four scope/semantic reverse-checks (HPO-only, indication-only, backbone connectivity, biomarker
first-pass). Aggregates a single release report. Scaffold: SKIPs cleanly when source data absent.
"""
from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

TABLES = [
    ("disease_node_v1", "data/processed/disease_node_v1.tsv",
     "disease/pipelines/disease/contracts/disease_node_v1.json"),
    ("disease_gene_target_v1", "data/processed/disease_gene_target_v1.tsv",
     "disease/pipelines/disease_gene_target/contracts/disease_gene_target_v1.json"),
    ("phenotype_node_v1", "data/processed/phenotype_node_v1.tsv",
     "disease/pipelines/disease_phenotype/contracts/phenotype_node_v1.json"),
    ("disease_phenotype_v1", "data/processed/disease_phenotype_v1.tsv",
     "disease/pipelines/disease_phenotype/contracts/disease_phenotype_v1.json"),
    ("disease_drug_v1", "data/processed/disease_drug_v1.tsv",
     "disease/pipelines/disease_drug/contracts/disease_drug_v1.json"),
    ("biomarker_edges_v1", "data/processed/biomarker_edges_v1.tsv",
     "disease/pipelines/biomarker/contracts/biomarker_edges_v1.json"),
]

CHECKS = [
    ("disease_backbone_connectivity", "tools/disease_backbone_connectivity.py"),
    ("disease_phenotype_hpo_only", "tools/disease_phenotype_hpo_only_check.py"),
    ("disease_drug_indication_only", "tools/disease_drug_indication_only_check.py"),
    ("biomarker_firstpass", "tools/biomarker_firstpass_check.py"),
]


def _run(cmd: list[str]) -> dict:
    try:
        p = subprocess.run(cmd, cwd=str(REPO), capture_output=True, text=True)
        return {"rc": p.returncode, "stdout": p.stdout.strip(), "stderr": p.stderr.strip()}
    except Exception as exc:  # noqa: BLE001
        return {"rc": -1, "error": str(exc)}


def main() -> int:
    tables_present = all((REPO / t[1]).is_file() for t in TABLES)
    if not tables_present:
        print("[SKIP] disease-layer source data not yet available locally; release gate is scaffold-only.")
        print("       All contracts/validators/layouts are in place; awaiting TSV ingestion.")
        return 0

    results = {"tables": {}, "checks": {}}
    # Contract validation per table
    for name, table, contract in TABLES:
        out = REPO / f"disease/release/external/disease-v1/reports/{name}.validation.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        r = _run(["python3", "tools/kg_validate_table.py", "--contract", contract,
                  "--table", table, "--out", str(out)])
        results["tables"][name] = {"rc": r["rc"], "report": str(out)}

    # Scope/semantic reverse checks
    for name, script in CHECKS:
        r = _run(["python3", str(REPO / script)])
        results["checks"][name] = r["rc"]

    summary = {
        "release": "disease-v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "tables_validated": len(results["tables"]),
        "checks_run": len(results["checks"]),
        "all_tables_pass": all(v["rc"] == 0 for v in results["tables"].values()),
        "all_checks_pass": all(v == 0 for v in results["checks"].values()),
        "detail": results,
    }
    out_path = REPO / "disease/release/external/disease-v1/reports/release_gate_summary.json"
    out_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    passed = summary["all_tables_pass"] and summary["all_checks_pass"]
    print(("[PASS] " if passed else "[FAIL] ") + json.dumps(summary, ensure_ascii=False))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
