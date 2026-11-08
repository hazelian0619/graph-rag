#!/usr/bin/env python3
"""Verify a starter package and run three offline examples (standard library only)."""

import argparse
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import sys


FIELDS = {
    "proteins": "uniprot_id entry_name protein_name taxon_id source fetch_date".split(),
    "molecules": "inchikey chembl_id source source_version fetch_date".split(),
    "activity_evidence": (
        "evidence_id compound_inchikey compound_chembl_id target_uniprot_accession "
        "target_chembl_id source_db source_record_id assay_type standard_type "
        "standard_relation standard_value standard_unit normalized_nM reference_doi "
        "reference_pubmed_id source_raw source_version fetch_date"
    ).split(),
}
IDS = {"proteins": "uniprot_id", "molecules": "inchikey", "activity_evidence": "evidence_id"}
MOLECULE_SOURCES = {
    "chembl_compound_structures.canonical_smiles",
    "chembl_compound_structures.canonical_smiles;chembl_compound_structures.standard_inchi",
}


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_package(data_dir):
    """Fail closed before producing any example output."""
    if data_dir.name != "data" or data_dir.is_symlink():
        raise ValueError("--data-dir must be the package's regular data/ directory")
    root = data_dir.parent
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    expected_paths = {f"data/{name}.tsv" for name in FIELDS} | {
        "NOTICE.md", "README.md", "examples/run_examples.py"
    }
    assets = manifest.get("assets")
    if not isinstance(assets, list) or len(assets) != len(expected_paths):
        raise ValueError("manifest must list all six package payloads")
    by_path = {}
    for asset in assets:
        path = asset["path"]
        if path not in expected_paths or path in by_path:
            raise ValueError(f"unexpected or duplicate manifest path: {path}")
        file_path = root / path
        if file_path.is_symlink() or not file_path.is_file():
            raise ValueError(f"missing or non-regular payload: {path}")
        if sha256(file_path) != asset["sha256"]:
            raise ValueError(f"SHA256 mismatch: {path}")
        by_path[path] = asset
    tables = {}
    for name, fields in FIELDS.items():
        path = data_dir / f"{name}.tsv"
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != fields:
                raise ValueError(f"invalid schema: {name}")
            rows = list(reader)
        asset = by_path[f"data/{name}.tsv"]
        if asset.get("fields") != fields or asset.get("rows") != len(rows):
            raise ValueError(f"manifest schema/row count mismatch: {name}")
        indexed = {}
        for row in rows:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"malformed TSV row: {name}")
            key = row[IDS[name]]
            if not key.strip() or key in indexed:
                raise ValueError(f"empty or duplicate ID in {name}: {key}")
            indexed[key] = row
        if not indexed:
            raise ValueError(f"empty table: {name}")
        tables[name] = indexed
    proteins, molecules, evidence = (tables[name] for name in FIELDS)
    if any(row["source"] != "UniProtKB/Swiss-Prot" for row in proteins.values()):
        raise ValueError("unapproved protein source")
    if any(row["source"] not in MOLECULE_SOURCES or
           row["source_version"] != "ChEMBL:chembl_36.db" or not row["chembl_id"].strip()
           for row in molecules.values()):
        raise ValueError("unapproved molecule source/version or missing ChEMBL ID")
    for row in evidence.values():
        protein = row["target_uniprot_accession"]
        molecule = row["compound_inchikey"]
        if protein not in proteins or molecule not in molecules:
            raise ValueError(f"endpoint closure failed: {row['evidence_id']}")
        if row["compound_chembl_id"] != molecules[molecule]["chembl_id"]:
            raise ValueError(f"compound ChEMBL ID mismatch: {row['evidence_id']}")
        if row["source_db"] != "chembl" or row["source_version"] != "chembl_36":
            raise ValueError("unapproved evidence source/version")
        if not (row["reference_doi"].strip() or row["reference_pubmed_id"].strip()):
            raise ValueError(f"missing DOI/PMID: {row['evidence_id']}")
    if {row["target_uniprot_accession"] for row in evidence.values()} != set(proteins):
        raise ValueError("unreferenced protein endpoints")
    if {row["compound_inchikey"] for row in evidence.values()} != set(molecules):
        raise ValueError("unreferenced molecule endpoints")
    return tables, sha256(root / "manifest.json")


def comparison_records(evidence):
    groups = {}
    for row in sorted(evidence.values(), key=lambda item: item["evidence_id"]):
        key = tuple(row[field] for field in (
            "compound_inchikey", "target_uniprot_accession", "assay_type",
            "standard_type", "standard_relation", "standard_unit"
        ))
        if any(not part.strip() for part in key):
            continue
        try:
            value = Decimal(row["normalized_nM"])
        except InvalidOperation:
            continue
        if not value.is_finite():
            continue
        groups.setdefault(key, []).append((value, row))
    for key in sorted(groups):
        first_value, first = groups[key][0]
        for value, row in groups[key][1:]:
            if value != first_value:
                return [first, row]
    raise ValueError("no pair with same endpoints/type/comparator/unit and distinct numeric values")


def make_examples(tables, manifest_hash, protein_id=None, molecule_id=None):
    proteins, molecules, evidence = (tables[name] for name in FIELDS)
    pair = comparison_records(evidence)
    protein_id = protein_id if protein_id is not None else pair[0]["target_uniprot_accession"]
    molecule_id = molecule_id if molecule_id is not None else pair[0]["compound_inchikey"]
    ordered = sorted(evidence.values(), key=lambda row: row["evidence_id"])
    protein = {
        "status": "recorded" if protein_id in proteins else "unrecorded",
        "query": {"uniprot_id": protein_id},
        "protein": proteins.get(protein_id),
        "evidence": [row for row in ordered if row["target_uniprot_accession"] == protein_id],
        "scope": "本包筛选子集；未收录不表示不存在该蛋白或作用。",
    }
    molecule = {
        "status": "recorded" if molecule_id in molecules else "unrecorded",
        "query": {"inchikey": molecule_id},
        "molecule": molecules.get(molecule_id),
        "evidence": [row for row in ordered if row["compound_inchikey"] == molecule_id],
        "scope": "本包中的 ChEMBL 实验记录；不据此推断临床疗效。",
    }
    pmids = {row["reference_pubmed_id"] for row in pair if row["reference_pubmed_id"].strip()}
    dois = {row["reference_doi"] for row in pair if row["reference_doi"].strip()}
    comparison = {
        "status": "recorded",
        "warning": "实验可比性未知",
        "match_fields": {field: pair[0][field] for field in (
            "compound_inchikey", "target_uniprot_accession", "assay_type",
            "standard_type", "standard_relation", "standard_unit"
        )},
        "records": pair,
        "same_pubmed_id": len(pmids) == 1 and all(row["reference_pubmed_id"].strip() for row in pair),
        "same_doi": len(dois) == 1 and all(row["reference_doi"].strip() for row in pair),
        "interpretation": "仅展示原记录的数值差异，不判定科学冲突或临床疗效。",
        "replication_note": "同一 PMID 不算独立复现；不同引用本身也不证明实验独立或可比。",
    }
    outputs = dict(zip(("01_protein.json", "02_molecule.json", "03_record_comparison.json"),
                       (protein, molecule, comparison)))
    for result in outputs.values():
        result["manifest_sha256"] = manifest_hash
    return outputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parents[1] / "data")
    parser.add_argument("--out", type=Path, default=Path("example-results"))
    parser.add_argument("--protein-id", help="exact UniProt accession; absent IDs are reported as unrecorded")
    parser.add_argument("--molecule-id", help="exact InChIKey; absent IDs are reported as unrecorded")
    args = parser.parse_args(argv)
    try:
        tables, manifest_hash = load_package(args.data_dir)
        outputs = make_examples(tables, manifest_hash, args.protein_id, args.molecule_id)
        protected = {args.data_dir.parent / "manifest.json"} | {
            args.data_dir.parent / path for path in (
                "README.md", "NOTICE.md", "examples/run_examples.py",
                *(f"data/{name}.tsv" for name in FIELDS)
            )
        }
        if any((args.out / name).resolve() in {p.resolve() for p in protected} for name in outputs):
            raise ValueError("output would overwrite package payload")
        args.out.mkdir(parents=True, exist_ok=True)
        for name, result in outputs.items():
            (args.out / name).write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({"status": "verified", "files": list(outputs)}, ensure_ascii=False))
    except (OSError, ValueError, KeyError, TypeError, csv.Error) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
