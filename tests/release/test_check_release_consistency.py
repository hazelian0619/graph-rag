from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "check_release_consistency.py"


def _write_tsv(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("c1\tc2\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def test_strict_local_pass(tmp_path: Path) -> None:
    table = tmp_path / "data" / "table.tsv"
    _write_tsv(table, [["a", "1"], ["b", "2"]])

    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps({"tables": [{"name": "table", "source_path": "data/table.tsv", "rows": 2}]}),
        encoding="utf-8",
    )

    val_report = tmp_path / "pipelines" / "x" / "reports" / "table.validation.json"
    val_report.parent.mkdir(parents=True, exist_ok=True)
    val_report.write_text(
        json.dumps({"table": "data/table.tsv", "row_count": 2, "passed": True}),
        encoding="utf-8",
    )

    index = tmp_path / "release" / "index.json"
    index.parent.mkdir(parents=True, exist_ok=True)
    index.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "generated_at_utc": "2026-01-01T00:00:00Z",
                "products": [
                    {
                        "product": "alpha",
                        "latest": {
                            "version": "v1",
                            "distribution_mode": "repository_snapshot",
                            "manifest_path": "manifest.json",
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    out = tmp_path / "release" / "consistency_report.json"
    cmd = [sys.executable, str(SCRIPT), "--index", "release/index.json", "--out", "release/consistency_report.json"]
    subprocess.run(cmd, cwd=tmp_path, check=True)

    report = json.loads(out.read_text(encoding="utf-8"))
    assert report["summary"]["status"] == "PASS"
    assert report["summary"]["issues"] == 0


def test_release_assets_fail_on_missing_sha(tmp_path: Path) -> None:
    manifest = tmp_path / "interaction_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "items": [
                    {
                        "table": "ppi_edges",
                        "release_available": True,
                        "release_assets": [{"name": "ppi_edges.tsv.zst"}],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    index = tmp_path / "release" / "index.json"
    index.parent.mkdir(parents=True, exist_ok=True)
    index.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "generated_at_utc": "2026-01-01T00:00:00Z",
                "products": [
                    {
                        "product": "interaction",
                        "latest": {
                            "version": "v1",
                            "distribution_mode": "github_release",
                            "consistency_mode": "release_assets",
                            "manifest_path": "interaction_manifest.json",
                        },
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    cmd = [sys.executable, str(SCRIPT), "--index", "release/index.json", "--out", "release/out.json"]
    proc = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)

    assert proc.returncode == 1
    report = json.loads((tmp_path / "release" / "out.json").read_text(encoding="utf-8"))
    assert report["summary"]["status"] == "FAIL"
    assert report["summary"]["issues"] >= 1


def _run_case(tmp_path, products, manifest=None, allow=False):
    if manifest is not None:
        (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    (tmp_path / "index.json").write_text(json.dumps({"products": products}))
    cmd = [sys.executable, str(SCRIPT), "--index", "index.json", "--out", "out.json"]
    if allow:
        cmd.append("--allow-issues")
    proc = subprocess.run(cmd, cwd=tmp_path, capture_output=True, text=True)
    return proc, json.loads((tmp_path / "out.json").read_text())


def test_empty_missing_and_unknown_manifests_fail(tmp_path):
    for products, manifest, issue in [
        ([], None, "no_products"),
        ([{"product": "x", "latest": {}}], None, "manifest_path_missing"),
        ([{"product": "x", "latest": {"manifest_path": "manifest.json"}}], {"other": []}, "no_tables_checked"),
        ([{"product": "x", "latest": {"manifest_path": "manifest.json", "consistency_mode": "release_assets"}}], {"other": []}, "no_assets_checked"),
        ([{"product": "x", "latest": {"manifest_path": "manifest.json", "consistency_mode": "unknown"}}], {}, "unknown_consistency_mode"),
    ]:
        proc, report = _run_case(tmp_path, products, manifest)
        assert proc.returncode == 1
        assert report["summary"]["status"] == "FAIL"
        assert any(i["type"] == issue for i in report["issues"])
    proc, report = _run_case(tmp_path, [], allow=True)
    assert proc.returncode == 0 and report["summary"]["status"] == "FAIL"


def test_assets_path_zero_rows_and_sha(tmp_path):
    import hashlib
    table = tmp_path / "empty.tsv"
    _write_tsv(table, [])
    product = [{"product": "x", "latest": {"manifest_path": "manifest.json"}}]
    manifest = {"assets": [{"path": "empty.tsv", "rows": 0,
                            "sha256": hashlib.sha256(table.read_bytes()).hexdigest()}]}
    proc, report = _run_case(tmp_path, product, manifest)
    assert proc.returncode == 0 and report["checks"][0]["tables_checked"] == 1
    _write_tsv(table, [["a", "1"]])
    proc, report = _run_case(tmp_path, product, manifest)
    assert proc.returncode == 1
    assert {i["type"] for i in report["issues"]} == {"checksum_mismatch", "manifest_row_mismatch"}


def test_release_assets_path_is_metadata_only(tmp_path):
    product = [{"product": "x", "latest": {"manifest_path": "manifest.json", "consistency_mode": "release_assets"}}]
    proc, report = _run_case(tmp_path, product, {"assets": [{"path": "elsewhere/file.tsv", "sha256": "abc"}]})
    assert proc.returncode == 0
    assert report["checks"][0]["verification_scope"] == "metadata_only"
    assert report["checks"][0]["assets_checked"] == 1


def test_duplicate_path_preserves_and_checks_sha(tmp_path):
    (tmp_path / "x.tsv").write_text("header\n")
    product = [{"product": "x", "latest": {"manifest_path": "manifest.json"}}]
    manifest = {"tables": [{"path": "x.tsv", "sha256": "bad"}], "assets": [{"path": "x.tsv"}]}
    proc, report = _run_case(tmp_path, product, manifest)
    assert proc.returncode == 1
    assert any(i["type"] == "checksum_mismatch" for i in report["issues"])


def test_existing_table_without_expectations_is_not_verified(tmp_path):
    (tmp_path / "x.tsv").write_text("header\na\n")
    product = [{"product": "x", "latest": {"manifest_path": "manifest.json"}}]
    proc, report = _run_case(tmp_path, product, {"assets": [{"path": "x.tsv"}]})
    assert proc.returncode == 1
    assert report["summary"]["status"] == "FAIL"
    assert report["checks"][0]["comparisons_checked"] == 0
    assert any(i["type"] == "table_expectations_missing" for i in report["issues"])


def test_duplicate_table_sha_ignores_hex_case_but_rejects_conflict(tmp_path):
    import hashlib
    (tmp_path / "x.tsv").write_text("header\n")
    digest = hashlib.sha256(b"header\n").hexdigest()
    product = [{"product": "x", "latest": {"manifest_path": "manifest.json"}}]
    manifest = {"tables": [{"path": "x.tsv", "rows": 0, "sha256": digest}],
                "artifacts": [{"path": "x.tsv", "rows": 0, "sha256": digest.upper()}]}
    proc, report = _run_case(tmp_path, product, manifest)
    assert proc.returncode == 0 and report["summary"]["status"] == "PASS"
    manifest["artifacts"][0]["sha256"] = "0" * 64
    proc, report = _run_case(tmp_path, product, manifest)
    assert proc.returncode == 1
    assert any(i["type"] == "manifest_conflict" for i in report["issues"])
