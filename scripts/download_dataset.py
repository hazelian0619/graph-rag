#!/usr/bin/env python3
"""Download datasets by product/version using protein/release/index.json.

Supports:
- GitHub Release assets download
- Repository snapshot direct URLs (for non-release products)
- Checksum verification from manifest and/or SHA256SUMS
- Interaction chunk merge (*.part.000 ...)
- Optional decompression (.gz, .zst)
"""

from __future__ import annotations

import argparse
import fnmatch
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
import tempfile
import unicodedata
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


GITHUB_API = "https://api.github.com"


def load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download_file(url: str, out_path: Path, timeout: int = 60) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "bio-entity-kg-foundation-downloader/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp, out_path.open("wb") as f:
        shutil.copyfileobj(resp, f)


def github_release_assets(repo: str, tag: str) -> List[Dict[str, Any]]:
    url = f"{GITHUB_API}/repos/{repo}/releases/tags/{tag}"
    req = urllib.request.Request(url, headers={"User-Agent": "bio-entity-kg-foundation-downloader/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    return payload.get("assets", [])


def manifest_checksum_map(manifest_path: Optional[Path]) -> Dict[str, str]:
    checksums: Dict[str, str] = {}
    if not manifest_path or not manifest_path.exists():
        return checksums

    data = load_json(manifest_path)

    identities: Dict[str, str] = {}

    def add(item: Any) -> None:
        if not isinstance(item, dict):
            return
        identity = str(item.get("path") or item.get("name") or "")
        name = str(item.get("name") or Path(identity).name)
        if not identity:
            return
        if not name or Path(name).name != name:
            raise ValueError(f"invalid manifest asset identity: {item}")
        checksum = str(item.get("sha256") or "")
        if name in identities and (identities[name] != identity or
                                   (checksum and name in checksums and checksums[name].lower() != checksum.lower())):
            raise ValueError(f"ambiguous manifest asset basename: {name}")
        identities[name] = identity
        if checksum:
            checksums[name] = checksum

    for key in ("tables", "reports", "assets"):
        for item in data.get(key, []) if isinstance(data.get(key), list) else []:
            add(item)
    for item in data.get("items", []) if isinstance(data.get("items"), list) else []:
        if isinstance(item, dict):
            for asset in item.get("release_assets", []):
                add(asset)

    return checksums


def parse_sha256sums(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    identities: Dict[str, str] = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) >= 2:
            identity = parts[1].lstrip("*")
            name = Path(identity).name
            if name in out and (identities[name] != identity or out[name].lower() != parts[0].lower()):
                raise ValueError(f"ambiguous checksum basename: {name}")
            identities[name] = identity
            out[name] = parts[0]
    return out


def merge_chunk_files(
    directory: Path, expected_names: Optional[List[str]] = None,
    checksums: Optional[Dict[str, str]] = None,
) -> List[Path]:
    pattern = re.compile(r"^(?P<base>.+)\.part\.(?P<num>\d+)$")
    groups: Dict[str, List[Tuple[int, Path]]] = {}
    for p in directory.iterdir():
        m = pattern.match(p.name) if p.is_file() else None
        if m:
            groups.setdefault(m.group("base"), []).append((int(m.group("num")), p))

    merged: List[Path] = []
    for base, parts in groups.items():
        names = [name for name in (expected_names or []) if (m := pattern.match(name)) and m.group("base") == base]
        numbers = sorted(int(pattern.match(name).group("num")) for name in names)
        actual = {p.name for _, p in parts}
        expected_hash = (checksums or {}).get(base)
        if not names or numbers != list(range(len(numbers))) or len(set(names)) != len(names):
            raise ValueError(f"missing or invalid expected chunk set for {base}")
        if actual != set(names):
            raise ValueError(f"chunk set mismatch for {base}: expected {sorted(names)}, found {sorted(actual)}")
        if not expected_hash:
            raise ValueError(f"missing whole-file checksum for {base}")
        out_path = directory / base
        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile(dir=directory, prefix=".merge-", delete=False) as out_f:
                tmp_path = Path(out_f.name)
                for _, part_path in sorted(parts):
                    with part_path.open("rb") as in_f:
                        shutil.copyfileobj(in_f, out_f)
            if sha256_file(tmp_path).lower() != expected_hash.lower():
                raise ValueError(f"whole-file checksum mismatch for {base}")
            tmp_path.replace(out_path)
            merged.append(out_path)
        finally:
            if tmp_path is not None:
                tmp_path.unlink(missing_ok=True)
    return merged


def decompress_files(directory: Path) -> List[Path]:
    created: List[Path] = []
    zstd_exists = shutil.which("zstd") is not None

    for p in sorted(directory.iterdir()):
        if not p.is_file():
            continue

        if p.name.endswith(".gz"):
            out = directory / p.name[:-3]
            with gzip.open(p, "rb") as src, out.open("wb") as dst:
                shutil.copyfileobj(src, dst)
            created.append(out)

        elif p.name.endswith(".zst"):
            out = directory / p.name[:-4]
            if zstd_exists:
                subprocess.run(["zstd", "-d", "-f", str(p), "-o", str(out)], check=True)
                created.append(out)
            else:
                print(f"[WARN] zstd not found; skip decompress: {p.name}")

    return created


def filter_assets(assets: List[Dict[str, Any]], patterns: List[str]) -> List[Dict[str, Any]]:
    if not patterns:
        return assets
    out = []
    for a in assets:
        name = a.get("name", "")
        if any(fnmatch.fnmatch(name, pat) for pat in patterns):
            out.append(a)
    return out


def validate_asset_names(assets: List[Dict[str, Any]]) -> None:
    """Require distinct NFC/casefold filenames before writing portable payloads."""
    names = set()
    for asset in assets:
        name = asset.get("name")
        if not isinstance(name, str) or not name or name in {".", ".."} or Path(name).name != name:
            raise SystemExit(f"[ERROR] invalid asset filename: {name!r}")
        identity = unicodedata.normalize("NFC", name).casefold()
        if identity in names:
            raise SystemExit(f"[ERROR] duplicate asset destination: {name}")
        names.add(identity)


def load_index_product(index_path: Path, product: str) -> Dict[str, Any]:
    idx = load_json(index_path)
    for p in idx.get("products", []):
        if p.get("product") == product:
            return p
    raise SystemExit(f"[ERROR] product not found in index: {product}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Download bio-entity KG dataset assets")
    parser.add_argument("--product", required=True, choices=["protein"])
    parser.add_argument("--version", default="latest", help="Version to fetch (currently supports latest)")
    parser.add_argument("--repo", default="hazelian0619/graph-rag", help="GitHub repo owner/name")
    parser.add_argument("--index", default="protein/release/index.json", help="Unified release index")
    parser.add_argument("--dest", default="downloads", help="Destination directory")
    parser.add_argument("--asset-glob", action="append", default=[], help="Filter assets by glob pattern")
    parser.add_argument("--skip-verify", action="store_true", help="Skip sha256 verification")
    parser.add_argument("--merge-chunks", action="store_true", help="Merge chunked assets (*.part.000...) after download")
    parser.add_argument("--decompress", action="store_true", help="Decompress .gz/.zst after download")
    args = parser.parse_args()

    if args.version != "latest":
        raise SystemExit("[ERROR] only --version latest is currently supported")

    product_meta = load_index_product(Path(args.index), args.product)
    latest = product_meta.get("latest", {})
    version = latest.get("version", "latest")

    out_dir = Path(args.dest) / args.product / version
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = Path(latest.get("manifest_path")) if latest.get("manifest_path") else None
    manifest_sha = manifest_checksum_map(manifest_path)

    downloaded: List[Path] = []
    release_names: List[str] = []

    if latest.get("distribution_mode") == "github_release":
        tag = latest.get("release_tag")
        if not tag:
            raise SystemExit(f"[ERROR] product {args.product} missing release_tag")

        assets = github_release_assets(args.repo, tag)
        release_names = [a["name"] for a in assets]
        assets = filter_assets(assets, args.asset_glob)
        if not assets:
            raise SystemExit(f"[ERROR] no assets matched for {args.product}:{tag}")
        validate_asset_names(assets)

        print(f"[INFO] downloading {len(assets)} assets from release {tag}")
        for asset in assets:
            name = asset["name"]
            url = asset["browser_download_url"]
            target = out_dir / name
            print(f"[DL] {name}")
            download_file(url, target)
            downloaded.append(target)

    else:
        artifacts = latest.get("artifacts", [])
        release_names = [a["name"] for a in artifacts]
        if args.asset_glob:
            artifacts = [
                a for a in artifacts
                if any(fnmatch.fnmatch(str(a.get("name", "")), pat) for pat in args.asset_glob)
            ]
        if not artifacts:
            raise SystemExit(f"[ERROR] product {args.product} has no artifacts in current.json")
        validate_asset_names(artifacts)

        print(f"[INFO] downloading {len(artifacts)} repository artifacts")
        for a in artifacts:
            name = a["name"]
            url = a["download_url"]
            target = out_dir / name
            print(f"[DL] {name}")
            download_file(url, target)
            downloaded.append(target)

    expected = dict(manifest_sha)
    # Only a checksum file downloaded during this invocation is trusted.
    sums = next((p for p in downloaded if p.name == "SHA256SUMS.txt"), None)
    if sums:
        for name, checksum in parse_sha256sums(sums).items():
            if name in expected and expected[name].lower() != checksum.lower():
                print(f"[FAIL] conflicting checksums for {name}")
                return 2
            expected[name] = checksum

    if args.skip_verify:
        print("[UNVERIFIED] checksum verification explicitly skipped")
    else:
        payloads = [p for p in downloaded if p.name != "SHA256SUMS.txt"]
        if not payloads:
            print("[FAIL] no payloads selected for verification")
            return 2
        for f in downloaded:
            exp = expected.get(f.name)
            if not exp and f.name == "SHA256SUMS.txt":
                continue
            if not exp:
                print(f"[FAIL] missing checksum for {f.name}")
                return 2
            got = sha256_file(f)
            if got.lower() != exp.lower():
                print(f"[FAIL] checksum mismatch for {f.name}: expected={exp} got={got}")
                return 2
        print("[OK] checksum verification passed for all selected payloads")

    if args.merge_chunks:
        try:
            merged = merge_chunk_files(out_dir, release_names, expected)
        except ValueError as exc:
            print(f"[FAIL] {exc}")
            return 2
        if merged:
            print(f"[OK] merged {len(merged)} chunked files")
            for p in merged:
                print("  -", p.name)

    if args.decompress:
        created = decompress_files(out_dir)
        if created:
            print(f"[OK] decompressed {len(created)} files")
            for p in created:
                print("  -", p.name)

    print(f"[DONE] downloaded assets to {out_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except urllib.error.HTTPError as e:
        print(f"[ERROR] HTTP {e.code}: {e.reason}", file=sys.stderr)
        raise
