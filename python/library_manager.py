#!/usr/bin/env python3
"""Unified, streaming reference-library downloader for SPECTRA (stdlib + numpy)."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import io
import json
import math
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path, PurePosixPath
from typing import Callable

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
UA = "SPECTRA-Chemistry-Toolkit/0.8 (open-scientific-library-manager)"
PACKS = ("nist_ftir", "flopp_ftir", "uvvis_pack", "rruff_ir", "toffolo_ir", "infra_art", "mnhn_paris", "zenodo_plastics", "nmrshiftdb2")


def _load_sibling(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, Path(__file__).with_name(filename))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load sibling helper {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


polymer = _load_sibling("download_polymer_pack.py", "spectra_download_polymer")
mega = _load_sibling("download_mega_ir.py", "spectra_download_mega")
uvhelpers = _load_sibling("download_uvvis_pack.py", "spectra_download_uvvis")


def _emit(line: str) -> None:
    print(line.replace("\r", " ").replace("\n", " "), flush=True)


def _progress(idx: int, total: int, pack_id: str, pct: float, detail: str) -> None:
    _emit(f"PROGRESS|{idx}|{total}|{pack_id}|{max(0.0, min(100.0, pct)):.1f}|{detail.replace('|', '/').replace(chr(10), ' ')}")


def _slug(text: str, max_len: int = 54) -> str:
    value = re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_").lower()
    return (value or "spectrum")[:max_len].strip("_") or "spectrum"


def _request(url: str, timeout: int = 45, headers: dict[str, str] | None = None):
    req_headers = {"User-Agent": UA}
    if headers:
        req_headers.update(headers)
    return urllib.request.urlopen(urllib.request.Request(url, headers=req_headers), timeout=timeout)


def stream_download(
    url: str,
    dest_path: Path,
    progress_cb: Callable[[float, str], None],
    timeout: int = 240,
) -> Path:
    """Stream to a sibling .part and atomically replace the destination."""
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    part_path = dest_path.with_name(dest_path.name + ".part")
    try:
        with _request(url, timeout=timeout) as response, part_path.open("wb") as out:
            total = int(response.headers.get("Content-Length") or 0)
            done = 0
            last = 0.0
            while True:
                chunk = response.read(65536)
                if not chunk:
                    break
                out.write(chunk)
                done += len(chunk)
                now = time.monotonic()
                if now - last >= 0.15:
                    pct = 0.0 if total <= 0 else min(99.9, done * 100.0 / total)
                    progress_cb(pct, f"Downloading {done / 1e6:.1f}/{total / 1e6:.1f} MB" if total else f"Downloading {done / 1e6:.1f} MB")
                    last = now
            out.flush()
        if total and done < total:
            raise IOError(f"truncated download ({done}/{total} bytes)")
        part_path.replace(dest_path)
        progress_cb(100.0, f"Downloaded {done / 1e6:.1f} MB")
        return dest_path
    except Exception:
        part_path.unlink(missing_ok=True)
        raise


def _http_bytes(url: str, timeout: int = 60) -> bytes:
    with _request(url, timeout=timeout) as response:
        return response.read()


def _cache_download(url: str, path: Path, cb: Callable[[float, str], None], min_bytes: int = 100_000, timeout: int = 300) -> Path:
    if path.is_file() and path.stat().st_size > min_bytes:
        cb(100.0, f"Using cached {path.name} ({path.stat().st_size / 1e6:.1f} MB)")
        return path
    return stream_download(url, path, cb, timeout=timeout)


def _numeric_pairs(text: str, minimum: int = 2) -> tuple[np.ndarray, np.ndarray] | None:
    xs: list[float] = []
    ys: list[float] = []
    for line in text.splitlines():
        line = line.strip().lstrip("\ufeff")
        if not line or line.startswith(("#", "//", "$$", "##")):
            continue
        fields = re.split(r"[,;\t\s]+", line)
        if len(fields) < 2:
            continue
        try:
            x, y = float(fields[0].strip('"')), float(fields[1].strip('"'))
        except (ValueError, OverflowError):
            continue
        if math.isfinite(x) and math.isfinite(y):
            xs.append(x)
            ys.append(y)
    if len(xs) < minimum:
        return None
    x_arr = np.asarray(xs, dtype=np.float64)
    y_arr = np.asarray(ys, dtype=np.float64)
    order = np.argsort(x_arr)
    x_arr, y_arr = x_arr[order], y_arr[order]
    keep = (x_arr >= 380.0) & (x_arr <= 4050.0)
    x_arr, y_arr = x_arr[keep], y_arr[keep]
    if len(x_arr) < minimum:
        return None
    if len(x_arr) > 900:
        grid = np.linspace(float(x_arr[0]), float(x_arr[-1]), 900)
        y_arr = np.interp(grid, x_arr, y_arr)
        x_arr = grid
    return x_arr, y_arr


def _write_ftir_csv(dest: Path, x: np.ndarray, y: np.ndarray, source: str) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    ymax, ymin = float(np.max(y)), float(np.min(y))
    y_label = "Transmittance (%T)" if 10.0 < ymax <= 115.0 and ymin >= -5.0 else "Absorbance"
    with dest.open("w", encoding="utf-8", newline="") as fh:
        fh.write(f"# Source: {source}\nWavenumber (cm-1),{y_label}\n")
        for xx, yy in zip(x[::-1], y[::-1]):
            fh.write(f"{xx:.2f},{yy:.6f}\n")


def _zip_entries(blob: bytes, depth: int = 0, prefix: str = ""):
    if depth > 5:
        return
    with zipfile.ZipFile(io.BytesIO(blob)) as archive:
        for info in archive.infolist():
            if info.is_dir():
                continue
            name = f"{prefix}/{info.filename}" if prefix else info.filename
            data = archive.read(info)
            if PurePosixPath(info.filename).suffix.lower() == ".zip":
                try:
                    yield from _zip_entries(data, depth + 1, name)
                except (zipfile.BadZipFile, OSError):
                    continue
            else:
                yield name, data


def _figshare_metadata() -> dict:
    return json.loads(_http_bytes("https://api.figshare.com/v2/articles/17070059", timeout=60).decode("utf-8"))


def _xlsx_spectra(blob: bytes):
    """Read basic numeric spectra from Figshare XLSX files without openpyxl."""
    from xml.etree import ElementTree as ET

    ns = {
        "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
        "p": "http://schemas.openxmlformats.org/package/2006/relationships",
    }
    with zipfile.ZipFile(io.BytesIO(blob)) as archive:
        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared = ["".join(t.text or "" for t in si.findall(".//m:t", ns)) for si in root.findall("m:si", ns)]
        book = ET.fromstring(archive.read("xl/workbook.xml"))
        rel_root = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        rels = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rel_root.findall("p:Relationship", ns)}
        sheet_paths = []
        for sheet in book.findall("m:sheets/m:sheet", ns):
            rel_id = sheet.attrib.get("{" + ns["r"] + "}id", "")
            target = rels.get(rel_id, "")
            path = target.lstrip("/") if target.startswith("/") else "xl/" + target
            if path.startswith("xl/xl/"):
                path = path[3:]
            sheet_paths.append((sheet.attrib.get("name", "spectrum"), path))

        def col_idx(ref: str) -> int:
            letters = re.match(r"[A-Z]+", ref.upper())
            value = 0
            for char in letters.group(0) if letters else "A":
                value = value * 26 + ord(char) - 64
            return value - 1

        for sheet_name, path in sheet_paths:
            if path not in archive.namelist():
                continue
            root = ET.fromstring(archive.read(path))
            grid: list[list[str]] = []
            for row in root.findall(".//m:sheetData/m:row", ns):
                cells: list[str] = []
                for cell in row.findall("m:c", ns):
                    ix = col_idx(cell.attrib.get("r", "A1"))
                    while len(cells) <= ix:
                        cells.append("")
                    kind = cell.attrib.get("t", "")
                    if kind == "inlineStr":
                        value = "".join(t.text or "" for t in cell.findall(".//m:t", ns))
                    else:
                        v = cell.find("m:v", ns)
                        value = v.text or "" if v is not None else ""
                        if kind == "s" and value:
                            try:
                                value = shared[int(value)]
                            except (ValueError, IndexError):
                                pass
                    cells[ix] = value
                grid.append(cells)
            if not grid:
                continue
            header_idx = next((ri for ri, row in enumerate(grid[:8]) if any("data(x)" in c.lower() for c in row)), -1)
            if header_idx >= 0:
                header = [c.strip().lower() for c in grid[header_idx]]
                pairs = [(i, i + 1) for i in range(len(header) - 1) if header[i] == "data(x)" and header[i + 1] == "data(y)"]
                for pair_idx, (cx, cy) in enumerate(pairs, 1):
                    points = []
                    for row in grid[header_idx + 1:]:
                        try:
                            points.append((float(row[cx]), float(row[cy])))
                        except (ValueError, IndexError):
                            continue
                    if len(points) >= 2:
                        parsed = _numeric_pairs("\n".join(f"{x},{y}" for x, y in points))
                        if parsed is not None:
                            suffix = f"_{pair_idx}" if len(pairs) > 1 else ""
                            yield f"{sheet_name}{suffix}", parsed
                if pairs:
                    continue
            # Common two-column workbook: one X/Y spectrum per sheet.
            points = []
            for row in grid:
                try:
                    points.append((float(row[0]), float(row[1])))
                except (ValueError, IndexError):
                    continue
            if len(points) >= 2:
                parsed = _numeric_pairs("\n".join(f"{x},{y}" for x, y in points))
                if parsed is not None:
                    yield sheet_name, parsed


def _download_flopp(root: Path, idx: int, total: int) -> int:
    out = root / "ftir"
    out.mkdir(parents=True, exist_ok=True)
    cache = root / ".cache" / "flopp_figshare_17070059.zip"
    meta = _figshare_metadata()
    files = [f for f in meta.get("files", []) if str(f.get("download_url", "")).startswith("http")]
    candidates = [f for f in files if str(f.get("name", "")).lower().endswith((".zip", ".xlsx"))]
    if not candidates:
        raise RuntimeError("Figshare article 17070059 has no downloadable ZIP or workbook attachment")

    def priority(item: dict) -> tuple[int, int, int]:
        name = str(item.get("name", "")).lower()
        ext = PurePosixPath(name).suffix
        kind = 0 if ext == ".zip" and "csv" in name else 1 if ext == ".zip" and ("jdx" in name or "jcamp" in name) else 2 if ext == ".zip" and "spa" in name else 3 if ext == ".zip" else 4 if ext == ".xlsx" else 5
        return (kind, -int(item.get("size", 0) or 0), len(name))

    chosen = sorted(candidates, key=priority)[0]
    attachment = str(chosen.get("name", ""))
    cache_path = cache if attachment.lower().endswith(".zip") else cache.with_suffix(".xlsx")
    _cache_download(chosen["download_url"], cache_path, lambda p, d: _progress(idx, total, "flopp_ftir", p, d), min_bytes=1_000, timeout=300)
    records: dict[tuple[str, str], tuple[int, str, bytes, tuple[np.ndarray, np.ndarray] | None]] = {}
    allowed = {".csv", ".txt", ".dpt", ".asc", ".spa", ".jdx", ".dx"}
    if attachment.lower().endswith(".xlsx"):
        for name, parsed in _xlsx_spectra(cache_path.read_bytes()):
            low = name.lower().replace("_", "-")
            sublib = "flopp_e" if any(t in low for t in ("flopp-e", "environ")) or "flopp_e" in low else "flopp"
            stem = _slug(name, 48)
            key = (sublib, stem)
            records[key] = (0, ".csv", b"", parsed)
    else:
        with zipfile.ZipFile(cache_path) as archive:
            members = sorted((m for m in archive.infolist() if not m.is_dir()), key=lambda m: m.filename.lower())
            for mi, info in enumerate(members):
                path = info.filename
                ext = PurePosixPath(path).suffix.lower()
                if ext == ".zip":
                    try:
                        nested = _zip_entries(archive.read(info), prefix=path)
                    except (zipfile.BadZipFile, OSError):
                        nested = ()
                    entries = nested
                elif ext in allowed:
                    entries = ((path, archive.read(info)),)
                else:
                    continue
                for member, raw in entries:
                    suffix = PurePosixPath(member).suffix.lower()
                    if suffix not in allowed:
                        continue
                    low_path = member.lower().replace("_", "-")
                    sublib = "flopp_e" if any(t in low_path for t in ("flopp-e", "environ")) or "flopp_e" in member.lower() else "flopp"
                    stem = _slug(PurePosixPath(member).stem, 48)
                    key = (sublib, stem)
                    rank = 0 if suffix in {".csv", ".txt", ".dpt", ".asc"} else 1 if suffix in {".jdx", ".dx"} else 2
                    parsed = None
                    if suffix in {".csv", ".txt", ".dpt", ".asc"}:
                        parsed = _numeric_pairs(raw.decode("latin-1", errors="ignore"), minimum=2)
                        if parsed is None:
                            continue
                    candidate = (rank, suffix, raw, parsed)
                    if key not in records or rank < records[key][0]:
                        records[key] = candidate
                if mi % 50 == 0:
                    _progress(idx, total, "flopp_ftir", min(99.0, 50 + 49 * mi / max(1, len(members))), f"Extracting archive {mi}/{len(members)}")
    added = 0
    source = "FLOPP/FLOPP-e (University of Toronto, De Frond et al. Anal. Chem. 2021) | DOI: 10.1021/acs.analchem.1c02549"
    for (sublib, stem), (rank, ext, raw, parsed) in sorted(records.items()):
        dest = out / f"poly_flopp_{sublib}_{stem}.csv"
        if parsed is not None:
            if not dest.exists():
                _write_ftir_csv(dest, *parsed, source)
                added += 1
        elif ext in {".jdx", ".dx", ".spa"}:
            dest = out / f"poly_flopp_{sublib}_{stem}{ext}"
            if not dest.exists():
                dest.write_bytes(raw)
                added += 1
    _progress(idx, total, "flopp_ftir", 100, f"Extracted {added} spectra")
    return added


def _download_nist_ftir(root: Path, idx: int, total: int) -> int:
    out = root / "ftir"
    out.mkdir(parents=True, exist_ok=True)
    compounds = [(slug, cas, "monomer") for slug, cas in polymer.NIST_TARGETS]
    compounds += [(slug, cas, "organic") for slug, cas in mega.NIST_ORGANIC_MEGA]
    added = 0
    for i, (slug, cas, group) in enumerate(compounds, 1):
        filename = f"{slug}_nist.jdx" if group == "monomer" else f"org_nist_{slug}.jdx"
        dest = out / filename
        if dest.is_file() and dest.stat().st_size > 400:
            _progress(idx, total, "nist_ftir", i * 100 / len(compounds), f"{i}/{len(compounds)} {slug} (installed)")
            continue
        url = f"https://webbook.nist.gov/cgi/cbook.cgi?JCAMP=C{cas.replace('-', '')}&Index=0&Type=IR"
        try:
            raw = _http_bytes(url, timeout=35)
            text = raw.decode("latin-1", errors="ignore")
            if "##TITLE=" in text and ("##XYDATA=" in text or "##XYPOINTS=" in text):
                dest.write_bytes(raw)
                added += 1
            detail = f"{i}/{len(compounds)} {slug} ({'saved' if dest.exists() else 'no IR record'})"
        except Exception as exc:
            detail = f"{i}/{len(compounds)} {slug} ({str(exc)[:100]})"
        _progress(idx, total, "nist_ftir", i * 100 / len(compounds), detail)
        time.sleep(0.15)
    return added


def _download_uvvis(root: Path, idx: int, total: int) -> int:
    out = root / "uvvis"
    out.mkdir(parents=True, exist_ok=True)
    added = 0
    try:
        index = _http_bytes(uvhelpers.PHOTOCHEMCAD_INDEX).decode("utf-8", errors="ignore")
        names = {num: name.strip() for num, name in re.findall(r'href="html/(\d+)\.html"[^>]*>([^<]+)<', index)}
        data_index = _http_bytes(uvhelpers.PHOTOCHEMCAD_DATA_INDEX).decode("utf-8", errors="ignore")
        pcad_files = re.findall(r'href="(\d+-abs\.txt)"', data_index)
    except Exception as exc:
        raise RuntimeError(f"PhotochemCAD index unavailable: {exc}") from exc
    tasks: list[tuple[str, str, str, str]] = [("pcad", f, f"OMLC PhotochemCAD {f}", f"pcad_{f.split('-', 1)[0]}_{uvhelpers.slugify(names.get(f.split('-', 1)[0], 'compound_' + f.split('-', 1)[0]))}.csv") for f in pcad_files]
    tasks += [("nist", name, cas, f"nist_uv_{uvhelpers.slugify(name)}.csv") for name, cas in uvhelpers.NIST_UVVIS_TARGETS]
    for i, (kind, name, value, filename) in enumerate(tasks, 1):
        dest = out / filename
        if dest.is_file() and dest.stat().st_size > 200:
            _progress(idx, total, "uvvis_pack", i * 100 / len(tasks), f"{i}/{len(tasks)} {filename} (installed)")
            continue
        try:
            if kind == "pcad":
                url = uvhelpers.PHOTOCHEMCAD_DATA_INDEX + name
                text = _http_bytes(url, timeout=30).decode("utf-8", errors="ignore")
                pairs = []
                for line in text.splitlines():
                    fields = re.split(r"\s+", line.strip())
                    try:
                        x, y = float(fields[0]), float(fields[1])
                        if 180 <= x <= 1100 and math.isfinite(y):
                            pairs.append((x, y))
                    except (ValueError, IndexError):
                        continue
                if len(pairs) >= 30:
                    pairs.sort()
                    x, y = np.asarray([p[0] for p in pairs]), np.asarray([p[1] for p in pairs])
                    yn = uvhelpers.normalize_absorbance(x, y)
                    x, yn = uvhelpers.resample_spectrum(x, yn, 800)
                    uvhelpers.write_uvvis_csv(dest, x, yn, f"OMLC PhotochemCAD ({names.get(name.split('-', 1)[0], name)})")
                    added += 1
            else:
                cas = value
                best = None
                for ix in range(4):
                    url = f"https://webbook.nist.gov/cgi/cbook.cgi?JCAMP=C{cas.replace('-', '')}&Index={ix}&Type=UVVis"
                    text = _http_bytes(url, timeout=30).decode("latin-1", errors="ignore")
                    if "##TITLE=" not in text or ("##XYDATA=" not in text and "##XYPOINTS=" not in text):
                        break
                    parsed = uvhelpers.parse_jcamp_dx(text)
                    if parsed is not None:
                        x, y = parsed
                        keep = (x >= 190) & (x <= 1000) & np.isfinite(y)
                        x, y = x[keep], y[keep]
                        if len(x) >= 30 and x[-1] >= 230 and x[-1] - x[0] >= 40:
                            score = (x[-1] >= 260, float(x[-1] - x[0]), len(x))
                            if best is None or score > best[0]:
                                best = (score, ix, x, y)
                    time.sleep(0.15)
                if best:
                    _, ix, x, y = best
                    y = uvhelpers.normalize_absorbance(x, y)
                    x, y = uvhelpers.resample_spectrum(x, y, 800)
                    uvhelpers.write_uvvis_csv(dest, x, y, f"NIST WebBook (CAS {cas}, Index={ix}, Type=UVVis)")
                    added += 1
            detail = f"{i}/{len(tasks)} {filename}"
        except Exception as exc:
            detail = f"{i}/{len(tasks)} {filename}: {str(exc)[:70]}"
        _progress(idx, total, "uvvis_pack", i * 100 / len(tasks), detail)
    return added


def _download_zip_record(record_id: str, root: Path, cache_name: str, idx: int, total: int, pack_id: str, min_bytes: int = 100_000) -> tuple[Path, dict]:
    meta = json.loads(_http_bytes(f"https://zenodo.org/api/records/{record_id}", timeout=60).decode("utf-8"))
    candidates = [f for f in meta.get("files", []) if str(f.get("key", "")).lower().endswith(".zip") and f.get("links", {}).get("self")]
    if not candidates:
        raise RuntimeError(f"Zenodo record {record_id} contains no ZIP file")
    chosen = max(candidates, key=lambda f: int(f.get("size", 0) or 0))
    path = root / ".cache" / cache_name
    _cache_download(chosen["links"]["self"], path, lambda p, d: _progress(idx, total, pack_id, p, d), min_bytes=min_bytes, timeout=360)
    return path, meta


def _download_rruff(root: Path, idx: int, total: int) -> int:
    out = root / "ftir"
    cache = root / ".cache" / "rruff_infrared_RAW.zip"
    _cache_download("https://www.rruff.net/zipped_data_files/infrared/RAW.zip", cache, lambda p, d: _progress(idx, total, "rruff_ir", p, d), min_bytes=100_000, timeout=240)
    written, counts = 0, {}
    with zipfile.ZipFile(cache) as archive:
        members = sorted(m for m in archive.namelist() if m.lower().endswith(".txt"))
        for i, member in enumerate(members, 1):
            base = PurePosixPath(member).name
            mineral = base.split("__")[0] if "__" in base else PurePosixPath(member).stem
            slug = _slug(mineral, 36)
            if counts.get(slug, 0) >= 2:
                continue
            parsed = _numeric_pairs(archive.read(member).decode("latin-1", errors="ignore"), minimum=120)
            if parsed is None:
                continue
            counts[slug] = counts.get(slug, 0) + 1
            suffix = "" if counts[slug] == 1 else f"_{counts[slug]}"
            dest = out / f"inorg_rruff_{slug}{suffix}.csv"
            if not dest.exists():
                _write_ftir_csv(dest, *parsed, f"RRUFF Infrared | Mineral: {mineral}")
                written += 1
            if i % 100 == 0:
                _progress(idx, total, "rruff_ir", 5 + 95 * i / max(1, len(members)), f"Extracting spectra {i}/{len(members)}")
    return written


def _download_toffolo(root: Path, idx: int, total: int) -> int:
    cache, _ = _download_zip_record("14170891", root, "zenodo_14170891.zip", idx, total, "toffolo_ir")
    records: dict[tuple[str, str], tuple[int, str, bytes, tuple[np.ndarray, np.ndarray] | None]] = {}
    with zipfile.ZipFile(cache) as archive:
        members = [m for m in archive.namelist() if not m.endswith("/")]
        for i, member in enumerate(members, 1):
            ext = PurePosixPath(member).suffix.lower()
            if ext not in {".csv", ".spa"}:
                continue
            parts = [p.lower() for p in PurePosixPath(member).parts]
            mode = "atr" if any("atr" in p for p in parts) else "trans" if any("trans" in p or "transmission" in p for p in parts) else "trans"
            stem = _slug(PurePosixPath(member).stem, 44)
            key = (mode, stem)
            raw = archive.read(member)
            parsed = _numeric_pairs(raw.decode("latin-1", errors="ignore")) if ext == ".csv" else None
            rank = 0 if ext == ".csv" and parsed is not None else 1 if ext == ".spa" else 2
            if key not in records or rank < records[key][0]:
                records[key] = (rank, ext, raw, parsed)
            if i % 100 == 0:
                _progress(idx, total, "toffolo_ir", 5 + 94 * i / max(1, len(members)), f"Inspecting archive {i}/{len(members)}")
    added = 0
    source = "Toffolo FTIR Library (CENIEH Spain, CC-BY-4.0) | DOI: 10.5281/zenodo.14170891"
    for (mode, stem), (_, ext, raw, parsed) in sorted(records.items()):
        if parsed is not None:
            dest = root / "ftir" / f"mat_toffolo_{mode}_{stem}.csv"
            if not dest.exists():
                _write_ftir_csv(dest, *parsed, source)
                added += 1
        elif ext == ".spa":
            dest = root / "ftir" / f"mat_toffolo_{mode}_{stem}.spa"
            if not dest.exists():
                dest.write_bytes(raw)
                added += 1
    return added


def _header_key(value: str) -> str:
    return re.sub(r"[^a-z0-9]", "", value.lower())


def _infra_metadata(archive: zipfile.ZipFile) -> dict[str, str]:
    names = [n for n in archive.namelist() if PurePosixPath(n).name.lower() == "metadata.csv"]
    if not names:
        return {}
    raw = archive.read(sorted(names, key=lambda n: (n.count("/"), n.lower()))[0]).decode("utf-8-sig", errors="replace")
    sample = raw[:8192]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
    except csv.Error:
        dialect = csv.excel
    reader = csv.DictReader(io.StringIO(raw), dialect=dialect)
    headers = reader.fieldnames or []
    normalized = {_header_key(h): h for h in headers}
    id_col = next((normalized[k] for k in ("sampleid", "id", "code", "samplecode", "specimenid", "filename", "file") if k in normalized), None)
    name_col = next((normalized[k] for k in ("materialname", "name", "material", "title", "compound", "sample", "object") if k in normalized), None)
    if id_col is None:
        id_col = next((h for h in headers if "id" in _header_key(h) or "code" in _header_key(h)), None)
    if name_col is None:
        name_col = next((h for h in headers if any(word in _header_key(h) for word in ("material", "name", "title", "compound"))), None)
    if not id_col or not name_col:
        return {}
    result = {}
    for row in reader:
        sample_id = (row.get(id_col) or "").strip()
        material = (row.get(name_col) or "").strip()
        if sample_id and material:
            result[sample_id.casefold()] = material
            result[PurePosixPath(sample_id).stem.casefold()] = material
    return result


def _infra_name(member: str, metadata: dict[str, str]) -> tuple[str, str]:
    stem = PurePosixPath(member).stem
    candidates = (stem, stem.rsplit("_", 1)[0], stem.split("_", 1)[0])
    sample_id = next((c for c in candidates if c.casefold() in metadata), "")
    if not sample_id:
        match = re.search(r"(?:ma|sample)?\d+[A-Za-z0-9_-]*", stem, re.I)
        if match and match.group(0).casefold() in metadata:
            sample_id = match.group(0)
    material = metadata.get(sample_id.casefold(), "") if sample_id else ""
    return (material or stem.replace("_", " "), sample_id or stem)


def _download_infra_art(root: Path, idx: int, total: int) -> int:
    cache, _ = _download_zip_record("21993883", root, "zenodo_21993883.zip", idx, total, "infra_art")
    out = root / "ftir"
    added = 0
    seen: set[str] = set()
    with zipfile.ZipFile(cache) as archive:
        metadata = _infra_metadata(archive)
        members = [m for m in archive.namelist() if m.lower().endswith((".asc", ".txt", ".csv")) and not any(s in m.lower() for s in ("metadata", "readme", "license"))]
        members += [m for m in archive.namelist() if m.lower().endswith((".dx", ".jdx"))]
        for i, member in enumerate(sorted(members), 1):
            ext = PurePosixPath(member).suffix.lower()
            material, sample_id = _infra_name(member, metadata)
            token = _slug(sample_id, 24)
            if token in seen:
                continue
            raw = archive.read(member)
            if ext in {".asc", ".txt", ".csv"}:
                parsed = _numeric_pairs(raw.decode("latin-1", errors="ignore"), minimum=120)
                if parsed is None:
                    continue
                material_slug = _slug(material, 30)
                tail = f"{material_slug}_{token}"[:44].rstrip("_")
                dest = out / f"mat_infra_{tail}.csv"
                if not dest.exists():
                    _write_ftir_csv(dest, *parsed, f"INFRA-ART (INOE 2000, Zenodo 21993883) | Material: {material} | ID: {sample_id}")
                    added += 1
            else:
                raw_head = raw[:4000].decode("latin-1", errors="ignore")
                if "##XYDATA" not in raw_head and "##XYPOINTS" not in raw_head:
                    continue
                material_slug = _slug(material, 30)
                tail = f"{material_slug}_{token}"[:44].rstrip("_")
                dest = out / f"mat_infra_{tail}{ext}"
                if not dest.exists():
                    dest.write_bytes(raw)
                    added += 1
            seen.add(token)
            if i % 100 == 0:
                _progress(idx, total, "infra_art", 5 + 94 * i / max(1, len(members)), f"Extracting spectra {i}/{len(members)}")
    return added


def _download_mnhn(root: Path, idx: int, total: int) -> int:
    cache, _ = _download_zip_record("10658337", root, "zenodo_10658337.zip", idx, total, "mnhn_paris")
    out = root / "ftir"
    added = 0
    with zipfile.ZipFile(cache) as archive:
        members = [m for m in archive.namelist() if m.lower().endswith((".jdx", ".dx"))]
        for i, member in enumerate(sorted(members), 1):
            upper = member.upper()
            prefix = "inorg_mnhn" if "/MIN/" in upper or "MIN_" in upper else "org_mnhn" if "/ORG/" in upper else "mat_mnhn"
            dest = out / f"{prefix}_{_slug(PurePosixPath(member).stem, 40)}.jdx"
            if dest.exists():
                continue
            raw = archive.read(member)
            head = raw[:4000].decode("latin-1", errors="ignore")
            if "##XYDATA" in head or "##XYPOINTS" in head or "##TITLE" in head:
                dest.write_bytes(raw)
                added += 1
            if i % 100 == 0:
                _progress(idx, total, "mnhn_paris", 5 + 94 * i / max(1, len(members)), f"Extracting spectra {i}/{len(members)}")
    return added


def _zenodo_plastics_candidate() -> tuple[str, dict]:
    candidates: list[tuple[int, str, dict]] = []
    queries = ["https://zenodo.org/api/records/10736649/versions"]
    urls = ["10736650", "10736649", "11176740"]
    for url in queries:
        try:
            data = json.loads(_http_bytes(url, timeout=60).decode("utf-8"))
            records = data if isinstance(data, list) else data.get("hits", {}).get("hits", data.get("versions", []))
            for record in records:
                if "metadata" not in record and "files" not in record:
                    continue
                record_id = str(record.get("id") or record.get("record_id") or "")
                if record_id:
                    candidates.append((int(record_id), record_id, record))
        except Exception:
            # The concept DOI resolves to 10736649 but Zenodo currently returns
            # 404 for this non-existent /versions route; the record IDs below are
            # the verified API fallback, including the actual archive record.
            pass
    for rid in urls:
        try:
            record = json.loads(_http_bytes(f"https://zenodo.org/api/records/{rid}", timeout=60).decode("utf-8"))
            candidates.append((int(rid), rid, record))
        except Exception:
            continue
    seen = set()
    for _, rid, meta in sorted(candidates, key=lambda row: row[0], reverse=True):
        if rid in seen:
            continue
        seen.add(rid)
        files = [f for f in meta.get("files", []) if int(f.get("size", 0) or 0) > 5 * 1024 * 1024 and (str(f.get("key", "")).lower().endswith((".zip", ".csv"))) and f.get("links", {}).get("self")]
        if files:
            return rid, max(files, key=lambda f: int(f.get("size", 0) or 0))
    raise RuntimeError("No FTIR-Plastics Zenodo version has a ZIP/CSV larger than 5 MB")


def _download_plastics(root: Path, idx: int, total: int) -> int:
    record_id, file_info = _zenodo_plastics_candidate()
    key = str(file_info["key"])
    cache = root / ".cache" / f"zenodo_plastics_{_slug(key, 70)}"
    if key.lower().endswith(".zip") and cache.suffix.lower() != ".zip":
        cache = cache.with_suffix(".zip")
    if key.lower().endswith(".csv") and cache.suffix.lower() != ".csv":
        cache = cache.with_suffix(".csv")
    _cache_download(file_info["links"]["self"], cache, lambda p, d: _progress(idx, total, "zenodo_plastics", p, d), min_bytes=5 * 1024 * 1024, timeout=360)
    before = {p.name for p in (root / "ftir").glob("*_real_*.csv")}
    if key.lower().endswith(".zip"):
        with zipfile.ZipFile(cache) as archive:
            members = [m for m in archive.namelist() if m.lower().endswith(".csv")]
            for i, member in enumerate(members, 1):
                polymer._extract_zenodo_csv_to_master(archive.read(member), root / "ftir", PurePosixPath(member).stem.lower())
                _progress(idx, total, "zenodo_plastics", 5 + 94 * i / max(1, len(members)), f"Extracting CSV {i}/{len(members)}")
    else:
        polymer._extract_zenodo_csv_to_master(cache.read_bytes(), root / "ftir", Path(key).stem.lower())
    after = {p.name for p in (root / "ftir").glob("*_real_*.csv")}
    return len(after - before)


def _download_nmr(root: Path, idx: int, total: int) -> int:
    cache = root / ".cache" / polymer.NMR_SDF_NAME
    url = polymer.NMR_SDF_URL
    if not (cache.is_file() and cache.stat().st_size > polymer.NMR_SDF_MIN_BYTES):
        try:
            direct = polymer._sourceforge_mirror_url(url, timeout=90)
        except Exception:
            direct = url
        _cache_download(direct, cache, lambda p, d: _progress(idx, total, "nmrshiftdb2", p, d), min_bytes=polymer.NMR_SDF_MIN_BYTES, timeout=600)
    else:
        _progress(idx, total, "nmrshiftdb2", 100, f"Using cached SDF ({cache.stat().st_size / 1e6:.1f} MB)")
    out = root / "nmr"
    out.mkdir(parents=True, exist_ok=True)
    by_key, by_cas, by_name = polymer._nmr_target_indexes()
    best: dict[tuple[str, str], dict] = {}
    records = 0
    for lines in polymer._iter_sdf_records(cache):
        records += 1
        fields = polymer._sdf_fields(lines)
        title = lines[0].strip() if lines else ""
        slug = polymer._match_target(fields, title, by_key, by_cas, by_name)
        if slug is not None:
            for nucleus, signal_idx, field, signals in polymer._record_spectra(fields):
                solvent = polymer._indexed_field(fields.get("Solvent", ""), signal_idx)
                polymer._best_spectrum(best, slug, nucleus, fields.get("nmrshiftdb2 ID", "").strip(), solvent, field, signals)
        if records % 250 == 0:
            _progress(idx, total, "nmrshiftdb2", min(99.0, records * 100 / 58000), f"Scanning SDF {records}/58000 records")
    added = 0
    for (slug, nucleus), entry in sorted(best.items()):
        dest = out / f"{slug}_{nucleus}.csv"
        if not dest.exists():
            polymer._write_nmr_csv(dest, entry["signals"], nucleus, entry["id"], entry["solvent"])
            added += 1
    _progress(idx, total, "nmrshiftdb2", 100, f"Scanned {records} records; wrote {added} spectra")
    return added


RUNNERS = {
    "nist_ftir": _download_nist_ftir,
    "flopp_ftir": _download_flopp,
    "uvvis_pack": _download_uvvis,
    "rruff_ir": _download_rruff,
    "toffolo_ir": _download_toffolo,
    "infra_art": _download_infra_art,
    "mnhn_paris": _download_mnhn,
    "zenodo_plastics": _download_plastics,
    "nmrshiftdb2": _download_nmr,
}


def check_urls() -> int:
    checks = [
        ("nist_ftir", "https://webbook.nist.gov/cgi/cbook.cgi?JCAMP=C100425&Index=0&Type=IR", "nist"),
        ("flopp_ftir", "https://api.figshare.com/v2/articles/17070059", "json"),
        ("uvvis_pack", "https://omlc.org/spectra/PhotochemCAD/data/", "html"),
        ("rruff_ir", "https://www.rruff.net/zipped_data_files/infrared/RAW.zip", "zip"),
        ("toffolo_ir", "https://zenodo.org/api/records/14170891", "json"),
        ("infra_art", "https://zenodo.org/api/records/21993883", "json"),
        ("mnhn_paris", "https://zenodo.org/api/records/10658337", "json"),
        ("zenodo_plastics", "https://zenodo.org/api/records/10736650", "json"),
        ("nmrshiftdb2", polymer.NMR_SDF_URL, "download"),
    ]
    print("Pack              Status  HTTP  Size / detected resource")
    print("----------------  ------  ----  ----------------------------------------------")
    failed = 0
    for pack_id, url, kind in checks:
        try:
            headers = {"Range": "bytes=0-0"} if kind in {"zip", "download"} else {}
            with _request(url, timeout=12, headers=headers) as response:
                status = response.status
                final_url = response.geturl()
                size = response.headers.get("Content-Length", "unknown")
                detail = response.headers.get("Content-Disposition", "")
                if kind == "json":
                    payload = json.loads(response.read().decode("utf-8"))
                    if pack_id == "flopp_ftir":
                        detected = ", ".join(f"{f.get('name')} ({int(f.get('size', 0))/1e6:.1f} MB)" for f in payload.get("files", [])[:3]) or "article metadata"
                    elif pack_id == "zenodo_plastics":
                        found = [f for f in payload.get("files", []) if int(f.get("size", 0) or 0) > 5 * 1024 * 1024 and str(f.get("key", "")).lower().endswith((".zip", ".csv"))]
                        detected = ", ".join(f"{f.get('key')} ({int(f.get('size', 0))/1e6:.1f} MB)" for f in found) or "record metadata; no data file"
                    else:
                        files = payload.get("files", [])
                        detected = ", ".join(f"{f.get('key')} ({int(f.get('size', 0))/1e6:.1f} MB)" for f in files[:2]) or "record metadata"
                    detail = detected
                elif kind == "nist":
                    detail = f"NIST IR JCAMP response ({size} bytes)"
                elif kind == "html":
                    detail = f"PhotochemCAD data index ({size} bytes)"
                elif kind in {"zip", "download"}:
                    detail = detail or f"{final_url.rsplit('/', 1)[-1]} ({size} bytes; range/redirect checked)"
                print(f"{pack_id:<17} OK      {status:<4}  {detail[:72]}")
        except urllib.error.HTTPError as exc:
            failed += 1
            print(f"{pack_id:<17} FAIL    {exc.code:<4}  {str(exc.reason)[:72]}")
        except Exception as exc:
            failed += 1
            print(f"{pack_id:<17} FAIL    ----  {str(exc)[:72]}")
    return 1 if failed else 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Download and manage SPECTRA reference libraries")
    parser.add_argument("--packs", help="Comma-separated pack IDs")
    parser.add_argument("--lib-root", type=Path, default=ROOT / "library")
    parser.add_argument("--check-urls", action="store_true", help="Check upstream endpoint health without bulk downloads")
    args = parser.parse_args()
    if args.check_urls:
        return check_urls()
    if not args.packs:
        parser.error("--packs is required unless --check-urls is used")
    requested = [item.strip() for item in args.packs.split(",") if item.strip()]
    invalid = [item for item in requested if item not in RUNNERS]
    if invalid:
        for pack_id in invalid:
            _emit(f"ERROR|{pack_id}|Unknown library pack")
        return 2
    if len(set(requested)) != len(requested):
        requested = list(dict.fromkeys(requested))
    root = args.lib_root.expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)
    total_added = 0
    failures = []
    for idx, pack_id in enumerate(requested, 1):
        start = time.monotonic()
        _progress(idx, len(requested), pack_id, 0, "Starting…")
        try:
            count = RUNNERS[pack_id](root, idx, len(requested))
            total_added += count
            elapsed = time.monotonic() - start
            _emit(f"PACK_DONE|{pack_id}|{count}|Added {count} spectra in {elapsed:.1f}s")
        except Exception as exc:
            message = str(exc).replace("|", "/").replace("\n", " ")
            failures.append((pack_id, message))
            _emit(f"ERROR|{pack_id}|{message}")
    if failures:
        _emit(f"DONE|{total_added}|Completed with {len(failures)} error(s); added {total_added} spectra")
        return 1
    _emit(f"DONE|{total_added}|Added {total_added} spectra across {len(requested)} selected pack(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
