#!/usr/bin/env python3
"""Build a deterministic, offline biomedical starter from existing local TSVs."""

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys
import tarfile


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
IDS = {"proteins": "uniprot_id", "molecules": "inchikey", "activity_evidence": "evidence_id"}
MOLECULE_SOURCES = {
    "chembl_compound_structures.canonical_smiles",
    "chembl_compound_structures.canonical_smiles;chembl_compound_structures.standard_inchi",
}
CODE_ROOT = Path(__file__).resolve().parents[1]
MAX_ARCHIVE_BYTES = 20 * 1024 * 1024


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_input(root, name):
    path = root / INPUTS[name]
    digest = sha256(path)
    rows = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = reader.fieldnames
        if not header or len(header) != len(set(header)) or not set(FIELDS[name]) <= set(header):
            raise ValueError(f"missing fields or duplicate column names: {INPUTS[name]}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"malformed row: {INPUTS[name]} line {reader.line_num}")
            key = row[IDS[name]]
            if not key.strip() or key in rows:
                raise ValueError(f"empty or duplicate ID in {name}: {key}; no duplicate is selected")
            rows[key] = {field: row[field] for field in FIELDS[name]}
    if sha256(path) != digest:
        raise ValueError(f"input changed while reading: {INPUTS[name]}")
    return rows, {"path": INPUTS[name], "sha256": digest, "rows": len(rows)}


def select_tables(tables):
    proteins, molecules, evidence = (tables[name] for name in FIELDS)
    proteins = {key: row for key, row in proteins.items() if row["source"] == "UniProtKB/Swiss-Prot"}
    molecules = {key: row for key, row in molecules.items() if row["source"] in MOLECULE_SOURCES
                 and row["source_version"] == "ChEMBL:chembl_36.db" and row["chembl_id"].strip()}
    excluded = Counter()
    kept = {}
    for key, row in evidence.items():
        molecule = molecules.get(row["compound_inchikey"])
        if row["source_db"] != "chembl" or row["source_version"] != "chembl_36":
            excluded["unapproved_evidence_source_or_version"] += 1
        elif row["target_uniprot_accession"] not in proteins or molecule is None:
            excluded["missing_or_unapproved_endpoint"] += 1
        elif row["compound_chembl_id"] != molecule["chembl_id"]:
            excluded["compound_chembl_id_mismatch"] += 1
        elif not (row["reference_doi"].strip() or row["reference_pubmed_id"].strip()):
            excluded["missing_doi_and_pmid"] += 1
        else:
            kept[key] = row
    if not kept:
        raise ValueError("no evidence satisfies the approved source, endpoint and reference filters")
    protein_ids = {row["target_uniprot_accession"] for row in kept.values()}
    molecule_ids = {row["compound_inchikey"] for row in kept.values()}
    selected = {
        "proteins": {key: proteins[key] for key in sorted(protein_ids)},
        "molecules": {key: molecules[key] for key in sorted(molecule_ids)},
        "activity_evidence": {key: kept[key] for key in sorted(kept)},
    }
    return selected, {
        "protein_source": "UniProtKB/Swiss-Prot",
        "molecule_sources": sorted(MOLECULE_SOURCES),
        "molecule_source_version": "ChEMBL:chembl_36.db",
        "evidence_source_db": "chembl", "evidence_source_version": "chembl_36",
        "endpoint_join": "exact target_uniprot_accession=uniprot_id and compound_inchikey=inchikey",
        "compound_id_join": "exact compound_chembl_id=chembl_id",
        "reference": "nonblank reference_doi OR reference_pubmed_id",
        "endpoints": "only endpoints referenced by retained evidence",
        "duplicate_policy": "reject all duplicate or empty primary IDs, including identical duplicates",
        "values": "exported strings unchanged; stable ascending primary-ID order",
        "excluded_protein_sources": len(tables["proteins"]) - len(proteins),
        "excluded_molecule_sources_or_versions_or_missing_chembl_id": len(tables["molecules"]) - len(molecules),
        "excluded_evidence_first_reason": dict(sorted(excluded.items())),
    }


def table_metadata(name, rows):
    if name == "proteins":
        attribution = "UniProt Consortium, UniProtKB/Swiss-Prot"
        source_url = "https://rest.uniprot.org/help/license"
        license_id, license_url = "CC-BY-4.0", "https://creativecommons.org/licenses/by/4.0/"
        versions = ["protein-v6 project snapshot; upstream UniProt release not recorded"]
        modifications = "Approved fields and referenced endpoints selected from the existing cleaned table."
    else:
        attribution = "ChEMBL team, EMBL-EBI"
        source_url = "https://www.ebi.ac.uk/chembldb/"
        license_id, license_url = "CC-BY-SA-3.0", "https://creativecommons.org/licenses/by-sa/3.0/"
        versions = sorted({row["source_version"] for row in rows.values()})
        modifications = "Approved fields and referenced endpoints selected; no structure text exported."
        if name == "activity_evidence":
            modifications = ("Subset selected from the existing activity sample. Project-generated evidence_id "
                             "and existing normalized_nM derivation retained; no values recomputed or rewritten.")
    return {
        "source": attribution, "source_url": source_url, "source_versions": versions,
        "fetch_dates": sorted({row["fetch_date"] for row in rows.values()}),
        "license": {"id": license_id, "url": license_url, "official_basis": source_url},
        "attribution": attribution, "modifications": modifications,
    }


def build(repo_root, out):
    if not out.name or out.name in {".", ".."}:
        raise ValueError("--out must name a package directory, for example build/starter-v1")
    tables, inputs = {}, {}
    for name in FIELDS:
        tables[name], inputs[name] = read_input(repo_root, name)
    selected, filters = select_tables(tables)
    payloads = {f"data/{name}.tsv" for name in FIELDS} | {"README.md", "NOTICE.md", "examples/run_examples.py"}
    # Refuse to package stale/unrelated files, and never overwrite source inputs.
    expected = payloads | {"manifest.json"}
    if out.exists():
        for path in out.rglob("*"):
            if path.is_symlink() or (path.is_file() and path.relative_to(out).as_posix() not in expected):
                raise ValueError(f"unexpected existing output entry: {path}")
    if out.is_symlink():
        raise ValueError("output directory must not be a symlink")
    (out / "data").mkdir(parents=True, exist_ok=True)
    (out / "examples").mkdir(exist_ok=True)
    assets = []
    for name, rows in selected.items():
        path = out / "data" / f"{name}.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=FIELDS[name], delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows.values())
        assets.append({"name": path.name, "path": f"data/{name}.tsv", "sha256": sha256(path),
                       "rows": len(rows), "fields": FIELDS[name], "input": inputs[name],
                       **table_metadata(name, rows)})
    for source, target in (("examples/run_examples.py", "examples/run_examples.py"),
                           ("examples/README.md", "README.md")):
        (out / target).write_bytes((CODE_ROOT / source).read_bytes())
    notice = ["# Attribution and licenses\n", "Data licenses are separate from the MIT code license.\n"]
    for asset in assets:
        notice.extend([
            f"## {asset['path']}\n",
            f"Attribution: {asset['attribution']}. Source: {asset['source_url']}.\n",
            f"Version: {', '.join(asset['source_versions'])}. Fetch dates: {', '.join(asset['fetch_dates'])}.\n",
            f"License: {asset['license']['id']} ({asset['license']['url']}).\n",
            f"Exported fields: {', '.join(asset['fields'])}.\n",
            f"Modifications: {asset['modifications']} Exported raw values are unchanged.\n",
        ])
    notice.extend([
        "Retain attribution, source and license links, and modification notices on redistribution. "
        "ChEMBL-derived adaptations are distributed under CC BY-SA 3.0. "
        "The source organizations do not endorse this project.\n",
        "This is a single-source activity subset, not independent replication, a scientific conflict "
        "assessment, new experimental research or evidence of clinical efficacy.\n",
        "## Code: MIT\n\n" + (CODE_ROOT / "LICENSE").read_text(encoding="utf-8"),
    ])
    (out / "NOTICE.md").write_text("\n".join(notice), encoding="utf-8")
    for path in sorted(payloads - {asset["path"] for asset in assets}):
        assets.append({"name": Path(path).name, "path": path, "sha256": sha256(out / path)})
    manifest = {"name": out.name, "assets": sorted(assets, key=lambda item: item["path"]),
                "filters": filters, "code_license": "MIT",
                "scope": "filtered single-source ChEMBL activity sample with referenced UniProt proteins"}
    (out / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    archive = out.with_name(out.name + ".tar.gz")
    with archive.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as zipped:
            with tarfile.open(fileobj=zipped, mode="w", format=tarfile.USTAR_FORMAT) as tar:
                for path in sorted(expected):
                    file_path = out / path
                    info = tarfile.TarInfo(f"{out.name}/{path}")
                    info.size = file_path.stat().st_size
                    info.mode = 0o644
                    info.uid = info.gid = info.mtime = 0
                    info.uname = info.gname = ""
                    with file_path.open("rb") as handle:
                        tar.addfile(info, handle)
    if archive.stat().st_size > MAX_ARCHIVE_BYTES:
        archive.unlink()
        raise ValueError("archive exceeds 20 MiB")
    digest = sha256(archive)
    (out.parent / "SHA256SUMS.txt").write_text(f"{digest}  {archive.name}\n", encoding="utf-8")
    return {"archive": archive.name, "bytes": archive.stat().st_size, "sha256": digest,
            "rows": {name: len(rows) for name, rows in selected.items()}}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--out", type=Path, default=Path("build/starter-v1"))
    args = parser.parse_args(argv)
    try:
        print(json.dumps(build(args.repo_root, args.out), sort_keys=True))
    except (OSError, ValueError, KeyError, csv.Error) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
