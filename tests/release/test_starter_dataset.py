"""Small end-to-end fixtures for the independent starter, without other release tools."""

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tarfile

import pytest


ROOT = Path(__file__).resolve().parents[2]
BUILDER = ROOT / "scripts/build_starter_dataset.py"
RUNNER = ROOT / "examples/run_examples.py"
INPUTS = {
    "proteins": "protein/release/external/protein-v6/tables/sequence/protein_master_v6_clean.tsv",
    "molecules": "molecule/release/external/molecule-v1/tables/identifiers/molecule_structure_identifiers_v1.tsv",
    "activity_evidence": "molecule/release/external/molecule-v1/tables/activity/molecule_activity_evidence_v2.sample.tsv",
}
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
RESULTS = ["01_protein.json", "02_molecule.json", "03_record_comparison.json"]


def invoke(script, *args, cwd=None):
    return subprocess.run([sys.executable, "-B", str(script), *map(str, args)], cwd=cwd,
                          text=True, capture_output=True)


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


@pytest.fixture
def inputs(tmp_path):
    protein = dict(zip(FIELDS["proteins"], ["P1", "PROTEIN_1", "Protein α", "9606",
                                          "UniProtKB/Swiss-Prot", "2025-10-26"]))
    molecule = dict(zip(FIELDS["molecules"], ["KEY1", "CHEMBL1",
        "chembl_compound_structures.canonical_smiles;chembl_compound_structures.standard_inchi",
        "ChEMBL:chembl_36.db", "2026-04-13"]))
    evidence = dict(zip(FIELDS["activity_evidence"], [
        "e1", "KEY1", "CHEMBL1", "P1", "CHEMBLT1", "chembl", "32706", "F", "EC50",
        "=", "3000", "nM", "3000", "10.1021/jm020046y", "12213051", "chembl_36",
        "chembl_36", "2026-04-12"
    ]))
    tables = {
        "proteins": [protein, {**protein, "uniprot_id": "UNREFERENCED"},
                     {**protein, "uniprot_id": "P_UNKNOWN", "source": "unknown"}],
        "molecules": [molecule, {**molecule, "inchikey": "UNREFERENCED"},
                      {**molecule, "inchikey": "KEY_UNKNOWN", "source": "pubchem_pug_rest.inchi"}],
        "activity_evidence": [
            {**evidence, "evidence_id": "e2", "source_record_id": "32707",
             "standard_value": "14900", "normalized_nM": "14900"},
            evidence,
            {**evidence, "evidence_id": "excluded-reference", "reference_doi": "", "reference_pubmed_id": ""},
            {**evidence, "evidence_id": "excluded-protein", "target_uniprot_accession": "P_UNKNOWN"},
            {**evidence, "evidence_id": "excluded-molecule", "compound_inchikey": "KEY_UNKNOWN"},
            {**evidence, "evidence_id": "excluded-match", "compound_chembl_id": "OTHER"},
            {**evidence, "evidence_id": "excluded-source", "source_db": "unknown"},
            {**evidence, "evidence_id": "excluded-version", "source_version": "chembl_99"},
        ],
    }
    root = tmp_path / "source"
    for name, rows in tables.items():
        # Excluded fields must never leak into the distributed data.
        write_tsv(root / INPUTS[name], FIELDS[name] + ["not_approved"],
                  [{**row, "not_approved": "private-structure-or-score"} for row in rows])
    return root


def build_fixture(inputs, tmp_path, name="starter-v1"):
    out = tmp_path / "build" / name
    result = invoke(BUILDER, "--repo-root", inputs, "--out", out)
    assert result.returncode == 0, result.stderr
    return out, json.loads(result.stdout)


def refresh_asset(package, name):
    path = package / "manifest.json"
    manifest = json.loads(path.read_text())
    for asset in manifest["assets"]:
        if asset["path"] == f"data/{name}.tsv":
            asset["sha256"] = hashlib.sha256((package / asset["path"]).read_bytes()).hexdigest()
            asset["rows"] = len(read_tsv(package / asset["path"])[1])
    path.write_text(json.dumps(manifest))


def test_build_extract_run_reproducible_and_attributed(inputs, tmp_path):
    package, report = build_fixture(inputs, tmp_path)
    archive = package.with_name("starter-v1.tar.gz")
    first_bytes = archive.read_bytes()
    assert 0 < report["bytes"] <= 20 * 1024 * 1024
    assert report["rows"] == {"proteins": 1, "molecules": 1, "activity_evidence": 2}
    assert (archive.parent / "SHA256SUMS.txt").read_text() == f"{hashlib.sha256(first_bytes).hexdigest()}  starter-v1.tar.gz\n"
    for name in FIELDS:
        fields, rows = read_tsv(package / "data" / f"{name}.tsv")
        assert fields == FIELDS[name]
        assert [r[fields[0]] for r in rows] == sorted(r[fields[0]] for r in rows)
        assert "private-structure-or-score" not in (package / "data" / f"{name}.tsv").read_text()
    manifest_text = (package / "manifest.json").read_text()
    assert str(inputs) not in manifest_text and str(tmp_path) not in manifest_text
    manifest = json.loads(manifest_text)
    for asset in manifest["assets"]:
        assert hashlib.sha256((package / asset["path"]).read_bytes()).hexdigest() == asset["sha256"]
        assert asset["name"] == Path(asset["path"]).name
        if asset["path"].startswith("data/"):
            name = Path(asset["path"]).stem
            assert asset["input"]["path"] == INPUTS[name]
            assert asset["input"]["sha256"] == hashlib.sha256((inputs / INPUTS[name]).read_bytes()).hexdigest()
            assert asset["license"]["id"] == ("CC-BY-4.0" if name == "proteins" else "CC-BY-SA-3.0")
            assert asset["modifications"] and asset["attribution"] and asset["source_url"]
    notice = (package / "NOTICE.md").read_text()
    assert (ROOT / "LICENSE").read_text() in notice
    assert "upstream UniProt release not recorded" in notice
    # A different parent directory cannot change package bytes or leak local paths.
    second_inputs = tmp_path / "copied-source"
    for name in FIELDS:
        fields, rows = read_tsv(inputs / INPUTS[name])
        write_tsv(second_inputs / INPUTS[name], fields, rows)
    again = invoke(BUILDER, "--repo-root", second_inputs, "--out", tmp_path / "other" / "starter-v1")
    assert again.returncode == 0, again.stderr
    assert first_bytes == (tmp_path / "other" / "starter-v1.tar.gz").read_bytes()
    # Rebuilding the same output is supported too.
    repeat = invoke(BUILDER, "--repo-root", second_inputs, "--out", package)
    assert repeat.returncode == 0, repeat.stderr
    assert first_bytes == archive.read_bytes()
    clean = tmp_path / "clean"
    clean.mkdir()
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            assert member.uid == member.gid == member.mtime == 0
            assert member.isfile() and member.mode == 0o644
        tar.extractall(clean, filter="data")
    extracted = clean / "starter-v1"
    result = invoke(extracted / "examples/run_examples.py", "--data-dir", extracted / "data",
                    "--out", clean / "results", cwd=clean)
    assert result.returncode == 0, result.stderr
    assert sorted(p.name for p in (clean / "results").iterdir()) == RESULTS
    protein, molecule, comparison = [json.loads((clean / "results" / name).read_text()) for name in RESULTS]
    assert protein["protein"]["protein_name"] == "Protein α"
    assert molecule["molecule"]["chembl_id"] == "CHEMBL1"
    assert comparison["warning"] == "实验可比性未知" and comparison["same_pubmed_id"] is True
    assert [r["source_record_id"] for r in comparison["records"]] == ["32706", "32707"]
    assert [r["normalized_nM"] for r in comparison["records"]] == ["3000", "14900"]
    assert all(r["reference_pubmed_id"] == "12213051" for r in comparison["records"])
    assert "独立复现" in comparison["replication_note"]


def test_custom_output_name_and_unknown_queries(inputs, tmp_path):
    package, report = build_fixture(inputs, tmp_path, "custom-name")
    assert report["archive"] == "custom-name.tar.gz"
    assert (package.parent / "SHA256SUMS.txt").read_text().endswith("  custom-name.tar.gz\n")
    result = invoke(package / "examples/run_examples.py", "--out", tmp_path / "results",
                    "--protein-id", "ABSENT", "--molecule-id", "ABSENT", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    for name, entity in zip(RESULTS[:2], ["protein", "molecule"]):
        row = json.loads((tmp_path / "results" / name).read_text())
        assert row["status"] == "unrecorded" and row[entity] is None and row["evidence"] == []


@pytest.mark.parametrize("name", list(FIELDS))
def test_duplicate_source_ids_rejected_even_with_conflicting_provenance(inputs, tmp_path, name):
    path = inputs / INPUTS[name]
    fields, rows = read_tsv(path)
    duplicate = {**rows[0], "source_db" if name == "activity_evidence" else "source": "unknown"}
    write_tsv(path, fields, rows + [duplicate])
    result = invoke(BUILDER, "--repo-root", inputs, "--out", tmp_path / "build" / "bad")
    assert result.returncode != 0 and "duplicate ID" in result.stderr
    assert not (tmp_path / "build" / "bad.tar.gz").exists()


@pytest.mark.parametrize("damage,expected", [
    ("hash", "SHA256 mismatch"), ("missing", "missing"), ("schema", "invalid schema"),
    ("duplicate", "duplicate ID"), ("endpoint", "endpoint closure"),
    ("chembl_id", "ChEMBL ID mismatch"), ("source", "unapproved evidence"),
    ("reference", "missing DOI/PMID"), ("equal_values", "no pair"),
])
def test_runner_fails_before_output_for_damaged_package(inputs, tmp_path, damage, expected):
    package, _ = build_fixture(inputs, tmp_path)
    path = package / "data/activity_evidence.tsv"
    fields, rows = read_tsv(path)
    if damage == "hash":
        path.write_text(path.read_text() + "\n")
    elif damage == "missing":
        path.unlink()
    else:
        if damage == "schema":
            fields[0] = "unexpected"
            rows = [{("unexpected" if k == "evidence_id" else k): v for k, v in r.items()} for r in rows]
        elif damage == "duplicate":
            rows[1]["evidence_id"] = rows[0]["evidence_id"]
        elif damage == "endpoint":
            rows[0]["target_uniprot_accession"] = "ABSENT"
        elif damage == "chembl_id":
            rows[0]["compound_chembl_id"] = "OTHER"
        elif damage == "source":
            rows[0]["source_db"] = "unknown"
        elif damage == "reference":
            rows[0]["reference_doi"] = rows[0]["reference_pubmed_id"] = ""
        elif damage == "equal_values":
            rows[1]["normalized_nM"] = "3e3"
        write_tsv(path, fields, rows)
        refresh_asset(package, "activity_evidence")
    out = tmp_path / "results"
    result = invoke(RUNNER, "--data-dir", package / "data", "--out", out)
    assert result.returncode != 0 and expected in result.stderr
    assert not out.exists()
