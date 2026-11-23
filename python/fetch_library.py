#!/usr/bin/env python3
"""
SPECTRA — Free open-spectral-library fetcher (stdlib only, zero pip installs).

Subcommands
-----------
  nist <query> [--out-dir library/ftir]
      Search NIST Chemistry WebBook and download the experimental IR JCAMP-DX.
  nmr <smiles> --nucleus 1h|13c [--name LABEL] [--out-dir library/nmr]
      Query nmrshiftdb2 and save a predicted/experimental JCAMP-DX spectrum.
  nmr-search <name|cas|inchikey> [--lib-root library]
      Extract a compound's real 1H/13C NMR offline from the cached
      nmrshiftdb2withsignals.sd dump (build the cache once with
      `spectra fetch pack --nmr`).
  import <source_dir_or_file> [--out-dir library/ftir]
      Recursively copy FLOPP / FLOPP-e / Open Specy / Zenodo FTIR references.
  uvvis <cas> [--out-dir library/uvvis]
      Download an experimental UV-Vis spectrum from NIST WebBook (Type=UVVis).
  list
      List every local .jdx/.csv reference in library/ftir, library/nmr, and library/uvvis.
  pack [--nist|--zenodo|--nmr|--uvvis|--all]
      Delegate to download_polymer_pack.py / download_uvvis_pack.py for the curated bulk download.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
from pathlib import Path

USER_AGENT = "SPECTRA/0.3.0 (offline-first spectroscopy toolkit)"
NIST_BASE = "https://webbook.nist.gov/cgi/cbook.cgi"
NMR_BASE = (
    "https://nmrshiftdb.nmr.uni-koeln.de/NmrshiftdbServlet/"
    "nmrshiftdbaction/searchorpredict"
)
SPECTRUM_EXTS = {".csv", ".jdx", ".dx", ".jcamp", ".dpt", ".txt", ".tsv", ".asc", ".xy"}


def slugify(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_").lower()
    return slug or "spectrum"


def fetch_bytes(url: str, timeout: int = 45) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def save_spectrum(url: str, out_path: Path) -> Path:
    data = fetch_bytes(url)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_bytes(data)
    print(f"saved {len(data)} bytes -> {out_path}")
    return out_path


def find_jcamp_link(html: str) -> str | None:
    # Prefer an explicit IR download query parameter.
    candidates = re.findall(r'href="([^"]*JCAMP=[^"]*)"', html, flags=re.IGNORECASE)
    for link in candidates:
        if "type=ir" in link.lower() or "units=si" in link.lower():
            return link
    return candidates[0] if candidates else None


def cmd_nist(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    query = args.query.strip()
    if re.fullmatch(r"\d+-\d+-\d+", query):
        url = f"{NIST_BASE}?ID=C{query}&Units=SI"
        slug = f"cas_{query}"
    else:
        url = f"{NIST_BASE}?Name={urllib.parse.quote(query)}&Units=SI"
        slug = slugify(query)

    print(f"searching NIST: {url}")
    html = fetch_bytes(url).decode("utf-8", errors="replace")
    link = find_jcamp_link(html)
    if not link:
        print("ERROR: no IR JCAMP-DX link found on the NIST page.", file=sys.stderr)
        return 1
    if link.startswith("/"):
        link = "https://webbook.nist.gov" + link
    elif not link.startswith("http"):
        link = urllib.parse.urljoin(url, link)

    save_spectrum(link, out_dir / f"{slug}_nist.jdx")
    return 0


def cmd_uvvis(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    cas = args.cas.strip().replace("-", "")
    # NIST WebBook UV-Vis endpoint: JCAMP=C{cas}&Index=0&Type=UVVis
    url = f"{NIST_BASE}?JCAMP=C{cas}&Index=0&Type=UVVis"
    slug = f"cas_{cas}"

    print(f"fetching NIST UV-Vis: {url}")
    try:
        data = fetch_bytes(url)
        text = data.decode("latin-1", errors="ignore")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    # Check if we got valid JCAMP data
    if "##TITLE=" not in text or ("##XYDATA=" not in text and "##XYPOINTS=" not in text):
        print("ERROR: no UV-Vis JCAMP-DX data found on the NIST page.", file=sys.stderr)
        return 1

    out_dir.mkdir(parents=True, exist_ok=True)
    # Save as .jdx (JCAMP-DX format)
    out_path = out_dir / f"{slug}_nist_uvvis.jdx"
    out_path.write_bytes(data)
    print(f"saved {len(data)} bytes -> {out_path}")

    # Also save a clean CSV for easier use in SPECTRA
    csv_path = out_dir / f"{slug}_nist_uvvis.csv"
    # We'll just save the raw JCAMP; the UI can parse it
    print("Note: For direct plotting, use `spectra fetch pack --uvvis` to get processed CSVs.")
    return 0


def cmd_nmr(args: argparse.Namespace) -> int:
    out_dir = Path(args.out_dir)
    nucleus = args.nucleus.lower()
    spectrum_type = "13C" if nucleus in {"13c", "c", "carbon"} else "1H"
    label = args.name or slugify(args.smiles)
    encoded = urllib.parse.quote(args.smiles, safe="")
    url = f"{NMR_BASE}/smiles/{encoded}/spectrumtype/{spectrum_type}/format/jcamp"

    print(f"querying nmrshiftdb2: {url}")
    try:
        save_spectrum(url, out_dir / f"{slugify(label)}_{nucleus}.jdx")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


def cmd_import(args: argparse.Namespace) -> int:
    src = Path(args.source)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if not src.exists():
        print(f"ERROR: {src} does not exist", file=sys.stderr)
        return 1

    copied = 0
    if src.is_file():
        files = [src]
    else:
        files = sorted(p for p in src.rglob("*") if p.is_file())
    for path in files:
        if path.suffix.lower() not in SPECTRUM_EXTS:
            continue
        target = out_dir / path.name
        if target.exists():
            stem, suffix = target.stem, target.suffix
            n = 2
            while target.exists():
                target = out_dir / f"{stem}_{n}{suffix}"
                n += 1
        shutil.copy2(path, target)
        copied += 1
        print(f"imported {path} -> {target}")
    print(f"done: {copied} reference spectra indexed in {out_dir}")
    return 0


def cmd_list(_args: argparse.Namespace) -> int:
    total = 0
    # Resolve library root relative to this script so it works from any CWD
    lib_root = Path(__file__).resolve().parent.parent / "library"
    for sub in ("ftir", "nmr", "uvvis"):
        directory = lib_root / sub
        print(f"[{directory}]")
        if not directory.is_dir():
            print("  (missing)")
            continue
        files = sorted(p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in SPECTRUM_EXTS)
        if not files:
            print("  (empty)")
        for p in files:
            size = p.stat().st_size
            print(f"  {p.name}  ({size} bytes)")
            total += 1
    print(f"total: {total} reference spectra")
    return 0


def cmd_pack(extra: list[str]) -> int:
    """Delegate `spectra fetch pack [...]` to download_polymer_pack.py or download_uvvis_pack.py."""
    # Check for --uvvis flag
    if "--uvvis" in extra or any(e.startswith("--uvvis") for e in extra):
        script = Path(__file__).resolve().parent / "download_uvvis_pack.py"
        if not script.is_file():
            print(f"ERROR: {script} not found", file=sys.stderr)
            return 1
        # Filter out --uvvis from extra args for download_uvvis_pack.py
        uvvis_args = [e for e in extra if e != "--uvvis"]
        if not uvvis_args:
            uvvis_args = ["--all"]
        cmd = [sys.executable, str(script), *uvvis_args]
        print("running:", " ".join(cmd))
        return subprocess.call(cmd)

    # Default to polymer pack
    script = Path(__file__).resolve().parent / "download_polymer_pack.py"
    if not script.is_file():
        print(f"ERROR: {script} not found", file=sys.stderr)
        return 1
    cmd = [sys.executable, str(script), *extra]
    print("running:", " ".join(cmd))
    return subprocess.call(cmd)


def cmd_nmr_search(args: argparse.Namespace) -> int:
    """Extract one compound's real 1H/13C NMR from the cached nmrshiftdb2 SDF."""
    script = Path(__file__).resolve().parent / "download_polymer_pack.py"
    if not script.is_file():
        print(f"ERROR: {script} not found", file=sys.stderr)
        return 1
    cmd = [sys.executable, str(script), "--search-nmr", args.query]
    if args.lib_root:
        cmd += ["--lib-root", args.lib_root]
    print("running:", " ".join(cmd))
    return subprocess.call(cmd)


def main() -> int:
    if len(sys.argv) >= 2 and sys.argv[1] == "pack":
        return cmd_pack(sys.argv[2:])

    parser = argparse.ArgumentParser(description="SPECTRA open spectral-library fetcher")
    sub = parser.add_subparsers(dest="command", required=True)

    # Resolve library root relative to this script so defaults work from any CWD
    _lib_root = Path(__file__).resolve().parent.parent / "library"
    _default_ftir = str(_lib_root / "ftir")
    _default_nmr = str(_lib_root / "nmr")
    _default_uvvis = str(_lib_root / "uvvis")

    p_nist = sub.add_parser("nist", help="download an experimental IR spectrum from NIST WebBook")
    p_nist.add_argument("query", help="compound name or CAS registry number")
    p_nist.add_argument("--out-dir", default=_default_ftir)
    p_nist.set_defaults(func=cmd_nist)

    p_nmr = sub.add_parser("nmr", help="download a JCAMP-DX spectrum from nmrshiftdb2")
    p_nmr.add_argument("smiles")
    p_nmr.add_argument("--nucleus", default="1h", choices=["1h", "13c"])
    p_nmr.add_argument("--name", default=None)
    p_nmr.add_argument("--out-dir", default=_default_nmr)
    p_nmr.set_defaults(func=cmd_nmr)

    p_import = sub.add_parser("import", help="index local FLOPP / Open Specy reference folders")
    p_import.add_argument("source")
    p_import.add_argument("--out-dir", default=_default_ftir)
    p_import.set_defaults(func=cmd_import)

    p_uvvis = sub.add_parser("uvvis", help="download an experimental UV-Vis spectrum from NIST WebBook")
    p_uvvis.add_argument("cas", help="CAS registry number (e.g., 71-43-2)")
    p_uvvis.add_argument("--out-dir", default=_default_uvvis)
    p_uvvis.set_defaults(func=cmd_uvvis)

    p_list = sub.add_parser("list", help="list local reference spectra")
    p_list.set_defaults(func=cmd_list)

    p_nmr_search = sub.add_parser(
        "nmr-search",
        help="extract any compound's real 1H/13C NMR offline from the cached nmrshiftdb2 SDF",
    )
    p_nmr_search.add_argument("query", help="compound name, CAS number or InChIKey")
    p_nmr_search.add_argument("--lib-root", default=None)
    p_nmr_search.set_defaults(func=cmd_nmr_search)

    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
