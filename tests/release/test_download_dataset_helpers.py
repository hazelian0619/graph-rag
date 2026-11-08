from __future__ import annotations

import importlib.util
import json
import hashlib
import sys

import pytest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "scripts" / "download_dataset.py"


spec = importlib.util.spec_from_file_location("download_dataset", MODULE_PATH)
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_parse_sha256sums(tmp_path: Path) -> None:
    p = tmp_path / "SHA256SUMS.txt"
    p.write_text(
        "aaabbbccc111  file_a.tsv.zst\n"
        "ddd222 *file_b.tsv.gz\n",
        encoding="utf-8",
    )
    got = mod.parse_sha256sums(p)
    assert got["file_a.tsv.zst"] == "aaabbbccc111"
    assert got["file_b.tsv.gz"] == "ddd222"


def test_merge_chunk_files(tmp_path: Path) -> None:
    d = tmp_path / "dl"
    d.mkdir()
    (d / "x.tsv.zst.part.000").write_bytes(b"abc")
    (d / "x.tsv.zst.part.001").write_bytes(b"def")
    (d / "x.tsv.zst.part.002").write_bytes(b"ghi")

    # A known complete release list and whole-file hash are necessary to detect a missing tail.
    names = [f"x.tsv.zst.part.{i:03d}" for i in range(3)]
    merged = mod.merge_chunk_files(d, names, {"x.tsv.zst": hashlib.sha256(b"abcdefghi").hexdigest()})
    assert len(merged) == 1
    out = d / "x.tsv.zst"
    assert out.exists()
    assert out.read_bytes() == b"abcdefghi"


def test_manifest_checksum_map(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "tables": [{"name": "a.tsv", "sha256": "111"}],
                "assets": [{"name": "b.tsv.zst", "sha256": "222"}],
                "items": [
                    {
                        "release_assets": [
                            {"name": "c.tsv.zst", "sha256": "333"},
                        ]
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    got = mod.manifest_checksum_map(manifest)
    assert got["a.tsv"] == "111"
    assert got["b.tsv.zst"] == "222"
    assert got["c.tsv.zst"] == "333"


def test_path_only_and_ambiguous_manifest(tmp_path: Path) -> None:
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"reports": [{"path": "dir/a.tsv", "sha256": "abc"}]}))
    assert mod.manifest_checksum_map(manifest) == {"a.tsv": "abc"}
    manifest.write_text(json.dumps({"assets": [{"path": "one/a.tsv", "sha256": "abc"},
                                                  {"path": "two/a.tsv", "sha256": "abc"}]}))
    with pytest.raises(ValueError, match="ambiguous"):
        mod.manifest_checksum_map(manifest)


@pytest.mark.parametrize("actual,expected,whole", [
    ([0, 2], [0, 1, 2], True),
    ([0, 1], [0, 1, 2], True),
    ([0, 1, 2], [0, 1], True),
    ([0, 2], [0, 2], True),
    ([0, 1], [0, 1], False),
])
def test_merge_rejects_unproven_or_wrong_chunk_set(tmp_path, actual, expected, whole):
    base = "data.tsv"
    (tmp_path / base).write_bytes(b"existing")
    for i in actual:
        (tmp_path / f"{base}.part.{i:03d}").write_bytes(bytes([i]))
    names = [f"{base}.part.{i:03d}" for i in expected]
    hashes = {base: hashlib.sha256(bytes(expected)).hexdigest()} if whole else {}
    with pytest.raises(ValueError):
        mod.merge_chunk_files(tmp_path, names, hashes)
    assert (tmp_path / base).read_bytes() == b"existing"


def test_merge_bad_whole_hash_preserves_output(tmp_path):
    (tmp_path / "x.part.000").write_bytes(b"a")
    (tmp_path / "x").write_bytes(b"old")
    with pytest.raises(ValueError, match="checksum mismatch"):
        mod.merge_chunk_files(tmp_path, ["x.part.000"], {"x": "bad"})
    assert (tmp_path / "x").read_bytes() == b"old"
    assert not list(tmp_path.glob(".merge-*"))


@pytest.mark.parametrize("assets,checksums,skip,result", [
    ([], {}, False, 1),
    (["a"], {}, False, 2),
    (["a", "b"], {"a": hashlib.sha256(b"a").hexdigest()}, False, 2),
    (["a"], {"a": "bad"}, False, 2),
    (["SHA256SUMS.txt"], {}, False, 2),
    (["a"], {}, True, 0),
])
def test_main_requires_verified_payloads(tmp_path, monkeypatch, capsys, assets, checksums, skip, result):
    index = tmp_path / "index.json"
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"assets": [{"name": k, "sha256": v} for k, v in checksums.items()]}))
    index.write_text(json.dumps({"products": [{"product": "protein", "latest": {
        "version": "v1", "distribution_mode": "github_release", "release_tag": "v1",
        "manifest_path": str(manifest)}}]}))
    dest = tmp_path / "dest"
    stale = dest / "protein" / "v1" / "SHA256SUMS.txt"
    stale.parent.mkdir(parents=True)
    stale.write_text(f"{hashlib.sha256(b'a').hexdigest()}  a\n")
    monkeypatch.setattr(mod, "github_release_assets", lambda *_: [
        {"name": name, "browser_download_url": name} for name in assets])
    monkeypatch.setattr(mod, "download_file", lambda url, path: path.write_bytes(url.encode()))
    monkeypatch.setattr(sys, "argv", ["download", "--product", "protein", "--index", str(index),
                                       "--dest", str(dest)] + (["--skip-verify"] if skip else []))
    if result == 1:
        with pytest.raises(SystemExit):
            mod.main()
    else:
        assert mod.main() == result
    output = capsys.readouterr().out
    if skip:
        assert "UNVERIFIED" in output and "verification passed" not in output


def test_download_uses_fresh_sums_for_payload(tmp_path, monkeypatch, capsys):
    index = tmp_path / "index.json"
    index.write_text(json.dumps({"products": [{"product": "protein", "latest": {
        "version": "v1", "distribution_mode": "github_release", "release_tag": "v1"}}]}))
    names = ["data.tsv", "SHA256SUMS.txt"]
    monkeypatch.setattr(mod, "github_release_assets", lambda *_: [
        {"name": n, "browser_download_url": n} for n in names])
    def download(url, path):
        content = b"payload" if url == "data.tsv" else (hashlib.sha256(b"payload").hexdigest() + "  data.tsv\n").encode()
        path.write_bytes(content)
    monkeypatch.setattr(mod, "download_file", download)
    monkeypatch.setattr(sys, "argv", ["download", "--product", "protein", "--index", str(index), "--dest", str(tmp_path / "out")])
    assert mod.main() == 0
    assert "verification passed for all selected payloads" in capsys.readouterr().out


def test_manifest_ambiguous_identity_without_checksum(tmp_path):
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"assets": [{"path": "one/a", "sha256": "abc"}, {"path": "two/a"}]}))
    with pytest.raises(ValueError, match="ambiguous"):
        mod.manifest_checksum_map(manifest)


def test_sums_reject_ambiguous_basename(tmp_path):
    sums = tmp_path / "SHA256SUMS.txt"
    sums.write_text("aaa one/a\naaa two/a\n")
    with pytest.raises(ValueError, match="ambiguous"):
        mod.parse_sha256sums(sums)


@pytest.mark.parametrize("mode", ["github_release", "repository_snapshot"])
@pytest.mark.parametrize("names", [("same.tsv", "same.tsv"), ("Data.tsv", "data.tsv"), ("\u00e9.tsv", "e\u0301.tsv")])
def test_duplicate_download_destinations_fail_before_writing(tmp_path, monkeypatch, mode, names):
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"assets": [{"name": name, "sha256": hashlib.sha256(b"second").hexdigest()}
                                              for name in dict.fromkeys(names)]}))
    artifacts = [{"name": name, "download_url": url, "browser_download_url": url}
                 for name, url in zip(names, ("first", "second"))]
    latest = {"version": "v1", "distribution_mode": mode, "release_tag": "v1",
              "manifest_path": str(manifest), "artifacts": artifacts}
    index = tmp_path / "index.json"
    index.write_text(json.dumps({"products": [{"product": "protein", "latest": latest}]}))
    downloads = []
    monkeypatch.setattr(mod, "github_release_assets", lambda *_: artifacts)
    monkeypatch.setattr(mod, "download_file", lambda url, target: downloads.append(url))
    monkeypatch.setattr(sys, "argv", ["download", "--product", "protein", "--index", str(index),
                                      "--dest", str(tmp_path / "out")])
    with pytest.raises(SystemExit, match="duplicate asset destination"):
        mod.main()
    assert downloads == []
