#!/usr/bin/env python3
"""
SPECTRA — Mega Experimental FTIR Library Downloader (Organic + Inorganic + Minerals)
Downloads & normalizes ~2,000+ real experimental FTIR spectra straight into library/ftir/:
  1. RRUFF Infrared (~14 MB): real mineral & inorganic FTIR spectra (University of Arizona)
  2. INFRA-ART (Zenodo 21993883): 1,085 ATR-FTIR spectra (inorganic/organic pigments, resins, waxes)
  3. MNHN Paris (Zenodo 10658337, ~48 MB): Bruker Vertex 70 ATR-FTIR (pure minerals + organics)
  4. NIST Chemistry WebBook: ~145 essential organic compounds (alkanes incl. n-butane, ketones,
     acids, esters, amines, aromatics, heterocycles, pharmaceuticals)
"""

from __future__ import annotations

import argparse
import io
import json
import re
import time
import urllib.request
import zipfile
from pathlib import Path

import numpy as np

UA = "SPECTRA-Chemistry-Toolkit/0.3 (academic-noncommercial)"

# ~145 Classic Organic Compounds for NIST WebBook IR (slug, CAS)
NIST_ORGANIC_MEGA: list[tuple[str, str]] = [
    # Alkanes & cycloalkanes
    ("n_butane", "106-97-8"),
    ("isobutane", "75-28-5"),
    ("n_pentane", "109-66-0"),
    ("n_hexane", "110-54-3"),
    ("n_heptane", "142-82-5"),
    ("n_octane", "111-65-9"),
    ("n_decane", "124-18-5"),
    ("cyclopentane", "287-92-3"),
    ("methylcyclohexane", "108-87-2"),
    ("isooctane", "540-84-1"),
    # Alkenes & alkynes
    ("1_hexene", "592-41-6"),
    ("cyclohexene", "110-83-8"),
    ("1_octene", "111-66-0"),
    ("isoprene", "78-79-5"),
    ("1_3_butadiene", "106-99-0"),
    ("alpha_methylstyrene", "98-83-9"),
    ("phenylacetylene", "536-74-3"),
    ("1_hexyne", "693-02-7"),
    # Aromatics
    ("ethylbenzene", "100-41-4"),
    ("o_xylene", "95-47-6"),
    ("m_xylene", "108-38-3"),
    ("p_xylene", "106-42-3"),
    ("cumene", "98-82-8"),
    ("mesitylene", "108-67-8"),
    ("naphthalene", "91-20-3"),
    ("anthracene", "120-12-7"),
    ("biphenyl", "92-52-4"),
    ("chlorobenzene", "108-90-7"),
    ("bromobenzene", "108-86-1"),
    ("fluorobenzene", "462-06-6"),
    ("nitrobenzene", "98-95-3"),
    ("benzonitrile", "100-47-0"),
    ("aniline", "62-53-3"),
    ("benzaldehyde", "100-52-7"),
    ("acetophenone", "98-86-2"),
    ("benzophenone", "119-61-9"),
    ("benzoic_acid", "65-85-0"),
    ("salicylic_acid", "69-72-7"),
    ("methyl_benzoate", "93-58-3"),
    ("ethyl_benzoate", "93-89-0"),
    ("benzyl_alcohol", "100-51-6"),
    ("p_cresol", "106-44-5"),
    ("hydroquinone", "123-31-9"),
    ("catechol", "120-80-9"),
    ("resorcinol", "108-46-3"),
    ("mehq_4_methoxyphenol", "150-76-5"),
    ("bht_antioxidant", "128-37-0"),
    ("diphenyl_ether", "101-84-8"),
    # Alcohols, diols & ethers
    ("1_propanol", "71-23-8"),
    ("1_butanol", "71-36-3"),
    ("tert_butanol", "75-65-0"),
    ("1_hexanol", "111-27-3"),
    ("1_octanol", "111-87-5"),
    ("cyclohexanol", "108-93-0"),
    ("glycerol", "56-81-5"),
    ("1_4_butanediol", "110-63-4"),
    ("propylene_glycol", "57-55-6"),
    ("diethylene_glycol", "111-46-6"),
    ("diethyl_ether", "60-29-7"),
    ("mtbe", "1634-04-4"),
    ("dimethoxyethane", "110-71-4"),
    ("dibutyl_ether", "142-96-1"),
    # Aldehydes & ketones
    ("acetaldehyde", "75-07-0"),
    ("propanal", "123-38-6"),
    ("butanal", "123-72-8"),
    ("hexanal", "66-25-1"),
    ("furfural", "98-01-1"),
    ("cinnamaldehyde", "104-55-2"),
    ("vanillin", "121-33-5"),
    ("mek_2_butanone", "78-93-3"),
    ("mibk", "108-10-1"),
    ("3_pentanone", "96-22-0"),
    ("cyclopentanone", "120-92-3"),
    ("cyclohexanone", "108-94-1"),
    ("acetylacetone", "123-54-6"),
    ("camphor", "76-22-2"),
    # Acids, esters, anhydrides, lactones
    ("formic_acid", "64-18-6"),
    ("propionic_acid", "79-09-4"),
    ("butyric_acid", "107-92-6"),
    ("hexanoic_acid", "142-62-1"),
    ("stearic_acid", "57-11-4"),
    ("oleic_acid", "112-80-1"),
    ("oxalic_acid", "144-62-7"),
    ("succinic_acid", "110-15-6"),
    ("citric_acid", "77-92-9"),
    ("lactic_acid", "50-21-5"),
    ("methyl_formate", "107-31-3"),
    ("methyl_acetate", "79-20-9"),
    ("butyl_acetate", "123-86-4"),
    ("methyl_acrylate", "96-33-3"),
    ("diethyl_phthalate", "84-66-2"),
    ("dibutyl_phthalate", "84-74-2"),
    ("dimethyl_carbonate", "616-38-6"),
    ("gamma_butyrolactone", "96-48-0"),
    ("epsilon_caprolactone", "502-44-3"),
    ("acetic_anhydride", "108-24-7"),
    ("succinic_anhydride", "108-30-5"),
    ("benzoyl_chloride", "98-88-4"),
    # Amines, amides, nitriles, nitro, heterocycles
    ("ethylamine", "75-04-7"),
    ("diethylamine", "109-89-7"),
    ("n_butylamine", "109-73-9"),
    ("cyclohexylamine", "108-91-8"),
    ("ethylenediamine", "107-15-3"),
    ("ethanolamine", "141-43-5"),
    ("formamide", "75-12-7"),
    ("acetamide", "60-35-5"),
    ("acrylamide", "79-06-1"),
    ("urea", "57-13-6"),
    ("propionitrile", "107-12-0"),
    ("butyronitrile", "109-74-0"),
    ("nitromethane", "75-52-5"),
    ("nitroethane", "79-24-3"),
    ("furan", "110-00-9"),
    ("thiophene", "110-02-1"),
    ("pyrrole", "109-97-7"),
    ("imidazole", "288-32-4"),
    ("indole", "120-72-9"),
    ("quinoline", "91-22-5"),
    ("morpholine", "110-91-8"),
    ("piperidine", "110-89-4"),
    # Halogenated, sulfur, silicon
    ("carbon_tetrachloride", "56-23-5"),
    ("1_2_dichloroethane", "107-06-2"),
    ("trichloroethylene", "79-01-6"),
    ("tetrachloroethylene", "127-18-4"),
    ("1_bromobutane", "109-65-9"),
    ("1_chlorobutane", "109-69-3"),
    ("carbon_disulfide", "75-15-0"),
    ("ethanethiol", "75-08-1"),
    ("sulfolane", "126-33-0"),
    ("teos_silicate", "78-10-4"),
    ("hmds", "999-97-3"),
    # Drugs & bio standards
    ("aspirin", "50-78-2"),
    ("caffeine", "58-08-2"),
    ("ibuprofen", "15687-27-1"),
    ("acetaminophen", "103-90-2"),
    ("glucose", "50-99-7"),
    ("sucrose", "57-50-1"),
    ("glycine", "56-40-6"),
    ("cholesterol", "57-88-5"),
    ("menthol", "89-78-1"),
    ("eugenol", "97-53-0"),
    ("limonene", "138-86-3"),
]


def slugify(text: str, max_len: int = 48) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_").lower()
    return (s[:max_len].strip("_")) or "spectrum"


def http_get_retry(url: str, attempts: int = 4, timeout: int = 45) -> bytes:
    delay = 4.0
    last: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:
            last = exc
            if attempt == attempts:
                break
            time.sleep(delay)
            delay *= 2.0
    assert last is not None
    raise last


def download_cached(url: str, cache_path: Path, timeout: int = 240) -> Path:
    if cache_path.is_file() and cache_path.stat().st_size > 100_000:
        print(f"  Reusing cached {cache_path.name} ({cache_path.stat().st_size / 1e6:.1f} MB)")
        return cache_path
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    part = cache_path.with_name(cache_path.name + ".part")
    print(f"  Downloading {cache_path.name} ...")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp, part.open("wb") as fh:
        done = 0
        while True:
            chunk = resp.read(1 << 16)
            if not chunk:
                break
            fh.write(chunk)
            done += len(chunk)
    part.replace(cache_path)
    print(f"  Cached -> {cache_path.name} ({done / 1e6:.1f} MB)")
    return cache_path


def parse_xy_text(text: str) -> tuple[np.ndarray, np.ndarray] | None:
    xs, ys = [], []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("//"):
            continue
        parts = re.split(r"[;,\t\s]+", line)
        if len(parts) < 2:
            continue
        try:
            x = float(parts[0])
            y = float(parts[1])
        except ValueError:
            continue
        if np.isfinite(x) and np.isfinite(y):
            xs.append(x)
            ys.append(y)
    if len(xs) < 120:
        return None
    xa = np.asarray(xs, dtype=np.float64)
    ya = np.asarray(ys, dtype=np.float64)
    order = np.argsort(xa)
    xa, ya = xa[order], ya[order]
    # Keep mid-IR window 400..4000 cm-1
    mask = (xa >= 380.0) & (xa <= 4050.0)
    xa, ya = xa[mask], ya[mask]
    if xa.size < 120 or (xa[-1] - xa[0]) < 800.0:
        return None
    # Resample gently to <=900 points for blazing-fast Rust matching
    if xa.size > 950:
        grid = np.linspace(xa[0], xa[-1], 900)
        ya = np.interp(grid, xa, ya)
        xa = grid
    return xa, ya


def write_clean_csv(dest: Path, xa: np.ndarray, ya: np.ndarray, source_note: str, force_abs: bool = False) -> None:
    y_max = float(np.max(ya))
    y_min = float(np.min(ya))
    is_trans = (not force_abs) and (y_max > 10.0 and y_max <= 115.0 and y_min >= -5.0)
    y_col = "Transmittance (%T)" if is_trans else "Absorbance"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as fh:
        fh.write(f"# {source_note}\n")
        fh.write(f"Wavenumber (cm-1),{y_col}\n")
        for x, y in zip(xa[::-1], ya[::-1]):
            fh.write(f"{x:.2f},{y:.5f}\n")


def download_rruff_inorganic(out_dir: Path, cache_dir: Path, max_per_mineral: int = 2) -> int:
    print("\n=== [1/4] RRUFF Project — Real Inorganic & Mineral FTIR (~14 MB) ===")
    try:
        zpath = download_cached(
            "https://www.rruff.net/zipped_data_files/infrared/RAW.zip",
            cache_dir / "rruff_infrared_RAW.zip",
            timeout=180,
        )
    except Exception as exc:
        print(f"  ✗ RRUFF download failed: {exc}")
        return 0

    written = 0
    counts: dict[str, int] = {}
    with zipfile.ZipFile(zpath, "r") as zf:
        members = [m for m in zf.namelist() if m.lower().endswith(".txt")]
        print(f"  Scanning {len(members)} RRUFF infrared files...")
        for m in sorted(members):
            base = Path(m).name
            mineral = base.split("__")[0] if "__" in base else Path(m).stem
            mslug = slugify(mineral, 36)
            if counts.get(mslug, 0) >= max_per_mineral:
                continue
            try:
                text = zf.read(m).decode("latin-1", errors="ignore")
            except Exception:
                continue
            parsed = parse_xy_text(text)
            if parsed is None:
                continue
            xa, ya = parsed
            counts[mslug] = counts.get(mslug, 0) + 1
            idx = counts[mslug]
            suffix = "" if idx == 1 else f"_{idx}"
            dest = out_dir / f"inorg_rruff_{mslug}{suffix}.csv"
            if not dest.is_file():
                write_clean_csv(dest, xa, ya, f"Source: RRUFF Infrared | Mineral: {mineral}", force_abs=False)
                written += 1
    print(f"  ✓ RRUFF inorganic/mineral spectra ready: {written} new ({len(counts)} unique minerals)")
    return written


def download_zenodo_zip(record_id: str, cache_dir: Path, max_mb: float = 220.0) -> Path | None:
    meta = json.loads(http_get_retry(f"https://zenodo.org/api/records/{record_id}", timeout=30).decode("utf-8"))
    zips = []
    for f in meta.get("files", []):
        key = str(f.get("key", ""))
        size_mb = float(f.get("size", 0)) / (1024 * 1024)
        url = f.get("links", {}).get("self")
        if url and key.lower().endswith(".zip") and size_mb <= max_mb:
            zips.append((size_mb, key, url))
    zips.sort()
    if not zips:
        return None
    size_mb, key, url = zips[0]
    print(f"  Zenodo {record_id}: {key} ({size_mb:.1f} MB)")
    return download_cached(url, cache_dir / f"zenodo_{record_id}_{slugify(key)}.zip", timeout=300)


def download_infra_art(out_dir: Path, cache_dir: Path) -> int:
    print("\n=== [2/4] INFRA-ART (Zenodo 21993883) — 1,085 Inorganic/Organic Pigments, Resins & Waxes ===")
    try:
        zpath = download_zenodo_zip("21993883", cache_dir, max_mb=250.0)
    except Exception as exc:
        print(f"  ✗ INFRA-ART failed: {exc}")
        return 0
    if not zpath:
        print("  - No ZIP archive found on Zenodo 21993883.")
        return 0

    written = 0
    seen_slugs: set[str] = set()
    with zipfile.ZipFile(zpath, "r") as zf:
        asc_files = [m for m in zf.namelist() if m.lower().endswith((".asc", ".txt", ".csv")) and "metadata" not in m.lower() and "readme" not in m.lower() and "license" not in m.lower()]
        dx_files = [m for m in zf.namelist() if m.lower().endswith((".dx", ".jdx"))]
        print(f"  Archive contains {len(asc_files)} ASCII and {len(dx_files)} JCAMP-DX files...")

        # Prefer clean baseline-corrected ASCII (.asc) normalized to fast 900-pt CSV
        for m in sorted(asc_files):
            stem = slugify(Path(m).stem, 42)
            if stem in seen_slugs:
                continue
            try:
                text = zf.read(m).decode("latin-1", errors="ignore")
            except Exception:
                continue
            parsed = parse_xy_text(text)
            if parsed is None:
                continue
            xa, ya = parsed
            seen_slugs.add(stem)
            dest = out_dir / f"mat_infra_{stem}.csv"
            if not dest.is_file():
                write_clean_csv(dest, xa, ya, f"Source: INFRA-ART (Zenodo 21993883) | {Path(m).stem}", force_abs=False)
                written += 1

        # Fallback: copy .dx files if ASCII list was empty
        if written == 0 and dx_files:
            for m in sorted(dx_files):
                stem = slugify(Path(m).stem, 42)
                dest = out_dir / f"mat_infra_{stem}.dx"
                if not dest.is_file():
                    raw = zf.read(m)
                    if b"##XYDATA" in raw or b"##XYPOINTS" in raw:
                        dest.write_bytes(raw)
                        written += 1

    print(f"  ✓ INFRA-ART spectra added: {written}")
    return written


def download_mnhn_paris(out_dir: Path, cache_dir: Path) -> int:
    print("\n=== [3/4] MNHN Paris (Zenodo 10658337, 48 MB) — Pure Minerals + Organics ATR-FTIR ===")
    try:
        zpath = download_zenodo_zip("10658337", cache_dir, max_mb=150.0)
    except Exception as exc:
        print(f"  ✗ MNHN Paris failed: {exc}")
        return 0
    if not zpath:
        return 0

    written = 0
    with zipfile.ZipFile(zpath, "r") as zf:
        jdx_files = [m for m in zf.namelist() if m.lower().endswith((".dx", ".jdx"))]
        print(f"  Found {len(jdx_files)} JCAMP-DX spectra in MNHN archive...")
        for m in sorted(jdx_files):
            up = m.upper()
            prefix = "inorg_mnhn" if "/MIN/" in up or "MIN_" in up else "org_mnhn" if "/ORG/" in up else "mat_mnhn"
            stem = slugify(Path(m).stem, 40)
            dest = out_dir / f"{prefix}_{stem}.jdx"
            if dest.is_file():
                continue
            try:
                raw = zf.read(m)
                head = raw[:4000].decode("latin-1", errors="ignore")
                if "##XYDATA" in head or "##XYPOINTS" in head or "##TITLE" in head:
                    dest.write_bytes(raw)
                    written += 1
            except Exception:
                continue
    print(f"  ✓ MNHN Paris spectra added: {written}")
    return written


def download_nist_organic_mega(out_dir: Path) -> int:
    print(f"\n=== [4/4] NIST Chemistry WebBook — {len(NIST_ORGANIC_MEGA)} Real Organic IR Spectra (.jdx) ===")
    ok, skip, fail = 0, 0, 0
    for slug, cas in NIST_ORGANIC_MEGA:
        dest = out_dir / f"org_nist_{slug}.jdx"
        legacy = out_dir / f"{slug}_nist.jdx"
        if (dest.is_file() and dest.stat().st_size > 400) or (legacy.is_file() and legacy.stat().st_size > 400):
            skip += 1
            continue
        cas_clean = cas.replace("-", "")
        url = f"https://webbook.nist.gov/cgi/cbook.cgi?JCAMP=C{cas_clean}&Index=0&Type=IR"
        try:
            data = http_get_retry(url, attempts=3, timeout=20)
            text = data.decode("latin-1", errors="ignore")
            if "##TITLE=" in text and ("##XYDATA=" in text or "##XYPOINTS=" in text):
                dest.write_bytes(data)
                ok += 1
                if ok % 15 == 0 or slug in ("n_butane", "aspirin", "caffeine", "ibuprofen"):
                    print(f"  ✓ [{ok}] {slug} ({len(data)//1024} KB)")
            else:
                fail += 1
        except Exception:
            fail += 1
        time.sleep(0.8)
    print(f"  ✓ NIST Organic Mega: {ok} new, {skip} already cached, {fail} unavailable.")
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description="Download Mega Organic + Inorganic FTIR library for SPECTRA")
    parser.add_argument("--rruff", action="store_true", help="RRUFF Inorganic/Minerals (~14 MB)")
    parser.add_argument("--infra", action="store_true", help="INFRA-ART Pigments/Resins/Inorganics (Zenodo)")
    parser.add_argument("--mnhn", action="store_true", help="MNHN Paris Minerals + Organics (Zenodo 48 MB)")
    parser.add_argument("--nist", action="store_true", help="NIST ~145 Organic compounds (~4 MB)")
    parser.add_argument("--all", action="store_true", help="Download all 4 sources (~2,000+ spectra)")
    parser.add_argument("--clean-cache", action="store_true", help="Delete downloaded ZIP caches after extracting")
    args = parser.parse_args()

    if not (args.all or args.rruff or args.infra or args.mnhn or args.nist):
        args.all = True

    root = Path(__file__).resolve().parent.parent / "library"
    out_dir = root / "ftir"
    cache_dir = root / ".cache"
    out_dir.mkdir(parents=True, exist_ok=True)

    added = 0
    if args.all or args.rruff:
        added += download_rruff_inorganic(out_dir, cache_dir)
    if args.all or args.infra:
        added += download_infra_art(out_dir, cache_dir)
    if args.all or args.mnhn:
        added += download_mnhn_paris(out_dir, cache_dir)
    if args.all or args.nist:
        added += download_nist_organic_mega(out_dir)

    if args.clean_cache:
        for z in cache_dir.glob("*.zip"):
            z.unlink(missing_ok=True)

    total = len([p for p in out_dir.iterdir() if p.is_file()])
    size_mb = sum(p.stat().st_size for p in out_dir.iterdir() if p.is_file()) / (1024 * 1024)
    print(f"\n✦ DONE! Added {added} new spectra -> library/ftir/ now has {total} real FTIR spectra ({size_mb:.1f} MB)!")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
