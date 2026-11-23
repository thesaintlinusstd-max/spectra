#!/usr/bin/env python3
"""
SPECTRA — UV-Vis Reference Library Downloader (stdlib + numpy)
Downloads and processes real experimental UV-Vis spectra from:
  1. OMLC PhotochemCAD: https://omlc.org/spectra/PhotochemCAD/
  2. NIST Chemistry WebBook (Type=UVVis): https://webbook.nist.gov/

Outputs clean CSVs to library/uvvis/ with header: "Wavelength (nm),Absorbance"
Normalized to 0..1.0 Absorbance, resampled to ≤800 points.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

UA = "SPECTRA-Chemistry-Toolkit/0.5 (academic-noncommercial)"


def http_get(url: str, timeout: int = 30) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def http_get_retry(url: str, attempts: int = 4, timeout: int = 30) -> bytes:
    delay = 6.0
    last: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            return http_get(url, timeout=timeout)
        except Exception as exc:
            last = exc
            if attempt == attempts:
                break
            time.sleep(delay)
            delay *= 2.0
    assert last is not None
    raise last


def slugify(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_").lower()
    return slug or "spectrum"


def resample_spectrum(x: np.ndarray, y: np.ndarray, max_points: int = 800) -> tuple[np.ndarray, np.ndarray]:
    """Resample spectrum to at most max_points using linear interpolation."""
    if len(x) <= max_points:
        return x, y
    x_new = np.linspace(x[0], x[-1], max_points)
    y_new = np.interp(x_new, x, y)
    return x_new, y_new


def normalize_absorbance(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Normalize absorbance to 0..1.0 range."""
    y_min = float(np.nanmin(y))
    y_max = float(np.nanmax(y))
    if y_max <= y_min:
        return np.zeros_like(y)
    return (y - y_min) / (y_max - y_min)


def write_uvvis_csv(dest: Path, x: np.ndarray, y: np.ndarray, source: str, lambda_max: float | None = None, eps_max: float | None = None) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as fh:
        header = f"# Source: {source}"
        if lambda_max is not None:
            header += f" | lambda_max={lambda_max:.1f} nm"
        if eps_max is not None:
            header += f" | eps_max={eps_max:.0f} M-1cm-1"
        fh.write(header + "\n")
        fh.write("Wavelength (nm),Absorbance\n")
        for xi, yi in zip(x, y):
            fh.write(f"{xi:.2f},{yi:.6f}\n")


# ---------------------------------------------------------------------------
# Source 1: OMLC PhotochemCAD
# ---------------------------------------------------------------------------
PHOTOCHEMCAD_INDEX = "https://omlc.org/spectra/PhotochemCAD/index.html"
PHOTOCHEMCAD_DATA_INDEX = "https://omlc.org/spectra/PhotochemCAD/data/"


def download_photchemcad(out_dir: Path) -> int:
    """Download all valid *-abs.txt files from OMLC PhotochemCAD.
    
    First fetches index.html to get compound names, then downloads data files
    from the data/ directory and saves with compound names like 
    pcad_041_beta_carotene.csv
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== [1/2] Downloading OMLC PhotochemCAD UV-Vis spectra -> {out_dir} ===")

    # Step 1: Fetch index.html to get compound name mapping (NNN -> name)
    name_map = {}
    try:
        index_html = http_get(PHOTOCHEMCAD_INDEX).decode("utf-8", errors="ignore")
        # Find html/NNN.html links with compound names
        matches = re.findall(r'href="html/(\d+)\.html"[^>]*>([^<]+)<', index_html)
        for num, name in matches:
            name_map[num] = name.strip()
        print(f"  Found {len(name_map)} compound names from index")
    except Exception as exc:
        print(f"  ⚠ Failed to fetch compound names from index: {exc}")

    # Step 2: Fetch the data directory index (Apache directory listing)
    try:
        data_index_html = http_get(PHOTOCHEMCAD_DATA_INDEX).decode("utf-8", errors="ignore")
    except Exception as exc:
        print(f"  ✗ Failed to fetch data index: {exc}")
        return 0

    # Find all NNN-abs.txt links in the directory listing
    abs_files = re.findall(r'href="(\d+-abs\.txt)"', data_index_html)
    if not abs_files:
        print("  - No NNN-abs.txt files found in data directory index")
        return 0

    print(f"  Found {len(abs_files)} absorption spectra")

    ok = 0
    skip = 0
    fail = 0

    # Remove old numeric-only files (pcad_NNN.csv) 
    for old_file in out_dir.glob("pcad_*.csv"):
        stem = old_file.stem  # e.g., "pcad_041"
        if stem.startswith("pcad_") and stem[5:].isdigit():
            print(f"  🗑 Removing old numeric-only file: {old_file.name}")
            old_file.unlink()

    for fname in abs_files:
        num = fname.replace("-abs.txt", "")
        compound_name = name_map.get(num, f"compound_{num}")
        slug = slugify(compound_name)
        dest = out_dir / f"pcad_{num}_{slug}.csv"

        if dest.is_file() and dest.stat().st_size > 200:
            skip += 1
            continue

        url = PHOTOCHEMCAD_DATA_INDEX + fname
        try:
            data = http_get_retry(url, attempts=3, timeout=20)
            text = data.decode("utf-8", errors="ignore").strip()

            # Parse the data: typically two columns (wavelength, molar extinction)
            # Tab-separated, with header lines starting with #
            lines = text.splitlines()
            xs = []
            ys = []
            for line in lines:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # Split by whitespace (tabs or spaces)
                parts = re.split(r"\s+", line)
                if len(parts) >= 2:
                    try:
                        x = float(parts[0])
                        y = float(parts[1])
                        if 180.0 <= x <= 1100.0 and math.isfinite(y):
                            xs.append(x)
                            ys.append(y)
                    except ValueError:
                        continue

            if len(xs) < 30:
                print(f"  - {num} ({compound_name[:30]}) (only {len(xs)} valid points, skipping)")
                fail += 1
                time.sleep(0.3)
                continue

            xs_arr = np.array(xs, dtype=np.float64)
            ys_arr = np.array(ys, dtype=np.float64)

            # Sort by wavelength
            order = np.argsort(xs_arr)
            xs_arr = xs_arr[order]
            ys_arr = ys_arr[order]

            # Find lambda_max
            idx_max = int(np.argmax(ys_arr))
            lambda_max = float(xs_arr[idx_max])
            eps_max = float(ys_arr[idx_max])  # raw molar extinction at lambda_max

            # Normalize to 0..1.0
            ys_norm = normalize_absorbance(xs_arr, ys_arr)

            # Resample to <= 800 points
            xs_resamp, ys_resamp = resample_spectrum(xs_arr, ys_norm, max_points=800)

            # Validate: >= 30 finite points in 180..1100 nm
            mask = np.isfinite(ys_resamp) & (xs_resamp >= 180.0) & (xs_resamp <= 1100.0)
            if int(np.sum(mask)) < 30:
                print(f"  - {num} ({compound_name[:30]}) (insufficient valid points after processing, skipping)")
                fail += 1
                time.sleep(0.3)
                continue

            source_desc = f"OMLC PhotochemCAD ({compound_name}, {fname})"
            write_uvvis_csv(dest, xs_resamp, ys_resamp, source_desc, lambda_max, eps_max)
            print(f"  ✓ {num} ({compound_name[:30]}) ({len(xs_resamp)} pts, λmax={lambda_max:.1f} nm)")
            ok += 1

        except Exception as exc:
            print(f"  ✗ {num} ({compound_name[:30]}) ({exc})")
            fail += 1

        time.sleep(0.3)  # Be polite

    print(f"PhotochemCAD summary: {ok} downloaded, {skip} cached, {fail} failed.")
    return ok


# ---------------------------------------------------------------------------
# Source 2: NIST WebBook UV-Vis
# ---------------------------------------------------------------------------
NIST_UVVIS_TARGETS: list[tuple[str, str]] = [
    ("benzene", "71-43-2"),
    ("toluene", "108-88-3"),
    ("styrene", "100-42-5"),
    ("naphthalene", "91-20-3"),
    ("anthracene", "120-12-7"),
    ("pyrene", "129-00-0"),
    ("phenanthrene", "85-01-8"),
    ("fluorene", "86-73-7"),
    ("biphenyl", "92-52-4"),
    ("phenol", "108-95-2"),
    ("aniline", "62-53-3"),
    ("pyridine", "110-86-1"),
    ("indole", "120-72-9"),
    ("quinoline", "91-22-5"),
    ("carbazole", "86-74-8"),
    ("acetone", "67-64-1"),
    ("acetophenone", "98-86-2"),
    ("benzophenone", "119-61-9"),
    ("benzaldehyde", "100-52-7"),
    ("nitrobenzene", "98-95-3"),
    ("benzoic_acid", "65-85-0"),
    ("salicylic_acid", "69-72-7"),
    ("aspirin", "50-78-2"),
    ("caffeine", "58-08-2"),
    ("uracil", "66-22-8"),
    ("thymine", "65-71-4"),
    ("azobenzene", "103-33-3"),
    ("coumarin", "91-64-5"),
    ("vanillin", "121-33-5"),
    ("furfural", "98-01-1"),
    ("bisphenol_a", "80-05-7"),
    ("trans_stilbene", "103-30-0"),
    ("phthalic_anhydride", "85-44-9"),
]


def parse_jcamp_dx(text: str) -> tuple[np.ndarray, np.ndarray] | None:
    """Parse JCAMP-DX format UV-Vis data from NIST.
    
    Handles both ##XYDATA and ##XYPOINTS formats.
    NIST UV-Vis typically uses ##XYPOINTS=(XY..XY) with comma-separated pairs.
    """
    lines = text.splitlines()

    x_units = None
    y_units = None
    x_unit_scale = 1.0       # unit conversion: microns -> nm
    x_is_wavenumber = False  # flag: X values are 1/cm
    x_factor = 1.0           # ##XFACTOR scale multiplier
    y_factor = 1.0           # ##YFACTOR scale multiplier
    y_is_log_epsilon = False

    for line in lines:
        line = line.strip()
        if line.startswith("##XUNITS="):
            x_units = line.split("=", 1)[1].strip().upper()
            if "MICROMETER" in x_units:
                x_unit_scale = 1000.0  # microns -> nm
            elif "1/CM" in x_units or "CM-1" in x_units:
                x_is_wavenumber = True
        elif line.startswith("##YUNITS="):
            y_units = line.split("=", 1)[1].strip().upper()
            if "LOG EPSILON" in y_units or "LOG(EPSILON)" in y_units or "LOGARITHM EPSILON" in y_units:
                y_is_log_epsilon = True
        elif line.startswith("##XFACTOR="):
            try:
                x_factor = float(line.split("=", 1)[1].strip())
            except ValueError:
                pass
        elif line.startswith("##YFACTOR="):
            try:
                y_factor = float(line.split("=", 1)[1].strip())
            except ValueError:
                pass

    def to_nm(raw_x: float) -> float | None:
        """Convert a raw X value to nanometres, or None if invalid."""
        if x_is_wavenumber:
            return 1e7 / raw_x if raw_x > 0 else None
        return raw_x * x_factor * x_unit_scale

    # Extract XYDATA or XYPOINTS
    xydata_lines = []
    in_xydata = False
    for line in lines:
        if line.startswith("##XYDATA") or line.startswith("##XYPOINTS"):
            in_xydata = True
            # Might have data on same line; strip format specifier like (XY..XY) or (X++(Y..Y))
            parts = line.split("=", 1)
            if len(parts) > 1:
                data_part = parts[1].strip()
                # Remove leading format specifier
                if data_part.startswith("("):
                    # Find closing paren and take everything after
                    close_idx = data_part.find(")")
                    if close_idx != -1:
                        data_part = data_part[close_idx + 1:].strip()
                if data_part:
                    xydata_lines.append(data_part)
            continue
        if in_xydata:
            if line.startswith("##"):
                break
            # Strip $$ comments
            if "$$" in line:
                line = line.split("$$")[0]
            xydata_lines.append(line.strip())

    if not xydata_lines:
        return None

    # Parse data: split on commas, semicolons, and whitespace
    data_str = " ".join(xydata_lines)
    # NIST uses comma-separated pairs: "162.4180,3.703814"
    tokens = re.split(r"[\s,;]+", data_str)
    tokens = [t for t in tokens if t]  # remove empty strings

    if len(tokens) < 4:
        return None

    xs = []
    ys = []

    # Parse as simple pairs (x, y, x, y, ...)
    try:
        vals = [float(t) for t in tokens]
        
        if len(vals) % 2 == 0:
            for i in range(0, len(vals), 2):
                raw_x = vals[i]
                y = vals[i + 1] * y_factor
                if math.isfinite(raw_x) and math.isfinite(y):
                    x_nm = to_nm(raw_x)
                    if x_nm is not None:
                        xs.append(x_nm)
                        ys.append(y)
        else:
            # Try X++(Y..Y) format: x_start, dx, then y values
            x_start = to_nm(vals[0])
            dx = vals[1] * x_factor * x_unit_scale
            if x_is_wavenumber:
                # Can't easily handle X++(Y..Y) for wavenumber
                return None
            y_vals = [v * y_factor for v in vals[2:]]
            if x_start is not None:
                for i, y in enumerate(y_vals):
                    x = x_start + i * dx
                    if math.isfinite(x) and math.isfinite(y):
                        xs.append(x)
                        ys.append(y)

    except ValueError:
        return None

    if len(xs) < 10:
        return None

    xs_arr = np.array(xs, dtype=np.float64)
    ys_arr = np.array(ys, dtype=np.float64)

    # Sort by wavelength
    order = np.argsort(xs_arr)
    xs_arr = xs_arr[order]
    ys_arr = ys_arr[order]

    # Convert LOG EPSILON to linear absorbance and normalize to 0..1
    if y_is_log_epsilon:
        # LOG EPSILON is log10(epsilon), convert to linear then normalize
        # epsilon = 10^log_epsilon
        # For relative absorbance, we can use epsilon directly and normalize
        ys_linear = 10.0 ** ys_arr
    else:
        ys_linear = ys_arr

    return xs_arr, ys_linear


def download_nist_uvvis(out_dir: Path, force: bool = False) -> int:
    """Download curated UV-Vis spectra from NIST WebBook.

    Probes Index=0..3 for each CAS and selects the record with the widest
    standard laboratory UV-Vis range (190-1000 nm, requiring max >= 230 nm
    and span >= 40 nm) so vacuum-UV-only stubs like 162-207 nm are skipped.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== [2/2] Downloading NIST WebBook UV-Vis spectra -> {out_dir} ===")

    # Remove any old vacuum-UV-only (<230 nm) stubs from previous runs
    for old in sorted(out_dir.glob("nist_uv_*.csv")):
        try:
            vals = []
            for ln in old.read_text(encoding="utf-8", errors="ignore").splitlines():
                if not ln or ln.startswith("#") or ln.lower().startswith("wavelength"):
                    continue
                vals.append(float(ln.split(",")[0]))
            if not vals or max(vals) < 230.0 or (max(vals) - min(vals)) < 40.0:
                print(f"  🗑 Removing vacuum-UV/truncated stub: {old.name} ({min(vals) if vals else 0:.0f}-{max(vals) if vals else 0:.0f} nm)")
                old.unlink()
        except Exception:
            pass

    ok = 0
    skip = 0
    fail = 0

    for name, cas in NIST_UVVIS_TARGETS:
        dest = out_dir / f"nist_uv_{slugify(name)}.csv"
        if not force and dest.is_file() and dest.stat().st_size > 200:
            skip += 1
            continue

        cas_clean = cas.replace("-", "")
        best: tuple[tuple[int, float, int], int, np.ndarray, np.ndarray] | None = None

        for idx in (0, 1, 2, 3):
            url = f"https://webbook.nist.gov/cgi/cbook.cgi?JCAMP=C{cas_clean}&Index={idx}&Type=UVVis"
            try:
                data = http_get_retry(url, attempts=3, timeout=20)
                text = data.decode("latin-1", errors="ignore")
            except Exception:
                break
            if "##TITLE=" not in text or ("##XYDATA=" not in text and "##XYPOINTS=" not in text):
                break

            parsed = parse_jcamp_dx(text)
            if parsed is not None:
                xs_raw, ys_raw = parsed
                # Standard bench UV-Vis range >= 190 nm (avoids vacuum-UV O2 band dominating scale)
                mask = (xs_raw >= 190.0) & (xs_raw <= 1000.0) & np.isfinite(ys_raw)
                xs_f = xs_raw[mask]
                ys_f = ys_raw[mask]
                if len(xs_f) >= 30:
                    span = float(xs_f[-1] - xs_f[0])
                    x_hi = float(xs_f[-1])
                    if x_hi >= 230.0 and span >= 40.0:
                        score = (1 if x_hi >= 260.0 else 0, span, len(xs_f))
                        if best is None or score > best[0]:
                            best = (score, idx, xs_f, ys_f)
            time.sleep(0.8)

        if best is None:
            print(f"  - {name:<22} (no full >=230 nm UV-Vis record on NIST)")
            if dest.is_file():
                dest.unlink()
            fail += 1
            continue

        _, best_idx, xs_arr, ys_raw = best
        ys_norm = normalize_absorbance(xs_arr, ys_raw)
        xs_resamp, ys_resamp = resample_spectrum(xs_arr, ys_norm, max_points=800)
        idx_max = int(np.argmax(ys_resamp))
        lambda_max = float(xs_resamp[idx_max])
        write_uvvis_csv(
            dest,
            xs_resamp,
            ys_resamp,
            f"NIST WebBook (CAS {cas}, Index={best_idx}, Type=UVVis)",
            lambda_max=lambda_max,
        )
        print(f"  ✓ {name:<22} Index={best_idx} ({len(xs_resamp)} pts, {xs_resamp[0]:.0f}–{xs_resamp[-1]:.0f} nm, λmax={lambda_max:.1f} nm)")
        ok += 1

    print(f"NIST UV-Vis summary: {ok} downloaded, {skip} cached, {fail} unavailable.")
    return ok


def main() -> int:
    parser = argparse.ArgumentParser(description="Download UV-Vis reference library pack")
    parser.add_argument("--all", action="store_true", help="Download both PhotochemCAD and NIST")
    parser.add_argument("--photochemcad", action="store_true", help="Download OMLC PhotochemCAD only")
    parser.add_argument("--nist", action="store_true", help="Download NIST WebBook UV-Vis only")
    parser.add_argument("--force-nist", action="store_true", help="Re-check all NIST UV-Vis indices")
    parser.add_argument("--out-dir", default="library/uvvis", help="Output directory")
    args = parser.parse_args()

    if not (args.all or args.photochemcad or args.nist):
        args.all = True

    out_dir = Path(args.out_dir)
    total = 0

    if args.all or args.photochemcad:
        total += download_photchemcad(out_dir)
    if args.all or args.nist:
        total += download_nist_uvvis(out_dir, force=args.force_nist)

    print(f"\nTotal UV-Vis reference spectra downloaded: {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())