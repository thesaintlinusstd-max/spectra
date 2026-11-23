#!/usr/bin/env python3
"""
SPECTRA — Real Experimental Spectroscopy Library Downloader (<200 MB total)
Downloads peer-reviewed, free experimental spectra directly into library/:
  1. NIST Chemistry WebBook (.jdx): 35+ monomers, crosslinkers & polymer solvents
  2. Zenodo FTIR-Plastics (record 11176740): real Jasco FT/IR-6700 spectra of
     PET, HDPE, LDPE, PP, PS, PVC (extracts clean median + representative CSVs)
  3. nmrshiftdb2 SDF dump (nmrshiftdb2withsignals.sd, ~158 MB, cached): real
     experimental 1H & 13C shift lists for monomers, crosslinkers & solvents,
     rebuilt as clean Lorentzian CSVs (also: --search-nmr for any compound)
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import io
import json
import math
import re
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import numpy as np

UA = "SPECTRA-Chemistry-Toolkit/0.3 (academic-noncommercial)"

# Essential polymer monomers, initiators, crosslinkers & solvents (Name, CAS)
NIST_TARGETS: list[tuple[str, str]] = [
    ("styrene", "100-42-5"),
    ("methyl_methacrylate", "80-62-6"),
    ("ethyl_acrylate", "140-88-5"),
    ("butyl_acrylate", "141-32-2"),
    ("acrylonitrile", "107-13-1"),
    ("vinyl_acetate", "108-05-4"),
    ("acrylic_acid", "79-10-7"),
    ("methacrylic_acid", "79-41-4"),
    ("ethylene_glycol", "107-21-1"),
    ("terephthalic_acid", "100-21-0"),
    ("epsilon_caprolactam", "105-60-2"),
    ("adipic_acid", "124-04-9"),
    ("hexamethylenediamine", "124-09-4"),
    ("bisphenol_a", "80-05-7"),
    ("epichlorohydrin", "106-89-8"),
    ("maleic_anhydride", "108-31-6"),
    ("phthalic_anhydride", "85-44-9"),
    ("phenol", "108-95-2"),
    ("toluene", "108-88-3"),
    ("benzene", "71-43-2"),
    ("acetone", "67-64-1"),
    ("chloroform", "67-66-3"),
    ("dichloromethane", "75-09-2"),
    ("tetrahydrofuran", "109-99-9"),
    ("dmf", "68-12-2"),
    ("dmac", "127-19-5"),
    ("nmp", "872-50-4"),
    ("dmso", "67-68-5"),
    ("methanol", "67-56-1"),
    ("ethanol", "64-17-5"),
    ("isopropanol", "67-63-0"),
    ("ethyl_acetate", "141-78-6"),
    ("1_4_dioxane", "123-91-1"),
    ("cyclohexane", "110-82-7"),
    ("anisole", "100-66-3"),
    ("pyridine", "110-86-1"),
    ("triethylamine", "121-44-8"),
    ("acetonitrile", "75-05-8"),
    ("acetic_acid", "64-19-7"),
]

# Compounds to extract from the nmrshiftdb2 SDF dump.
# (slug, display name, CAS, aliases/trade names as they appear in SDfiles)
NMR_SDF_TARGETS: list[tuple[str, str, str, tuple[str, ...]]] = [
    ("styrene", "Styrene", "100-42-5", ("styrene", "ethenylbenzene", "vinylbenzene", "phenylethene", "styrol")),
    ("methyl_methacrylate", "Methyl methacrylate", "80-62-6", ("methyl methacrylate", "methyl 2-methylprop-2-enoate", "methyl 2-methylacrylate", "2-methylpropenoic acid methyl ester")),
    ("vinyl_acetate", "Vinyl acetate", "108-05-4", ("vinyl acetate", "ethenyl acetate", "1-acetoxyethylene", "acetic acid vinyl ester")),
    ("acrylonitrile", "Acrylonitrile", "107-13-1", ("acrylonitrile", "prop-2-enenitrile", "2-propenenitrile", "vinyl cyanide")),
    ("acrylic_acid", "Acrylic acid", "79-10-7", ("acrylic acid", "prop-2-enoic acid", "2-propenoic acid", "propenoic acid", "acroleic acid")),
    ("methacrylic_acid", "Methacrylic acid", "79-41-4", ("methacrylic acid", "2-methylprop-2-enoic acid", "2-methylpropenoic acid", "alpha-methylacrylic acid")),
    ("ethylene_glycol", "Ethylene glycol", "107-21-1", ("ethylene glycol", "ethane-1,2-diol", "1,2-ethanediol", "monoethylene glycol")),
    ("epsilon_caprolactam", "epsilon-Caprolactam", "105-60-2", ("epsilon-caprolactam", "caprolactam", "azepan-2-one", "hexano-6-lactam", "6-hexanelactam", "2-azepanone")),
    ("bisphenol_a", "Bisphenol A", "80-05-7", ("bisphenol a", "4,4'-isopropylidenediphenol", "2,2-bis(4-hydroxyphenyl)propane", "bpa", "bisphenol-a")),
    ("epichlorohydrin", "Epichlorohydrin", "106-89-8", ("epichlorohydrin", "2-(chloromethyl)oxirane", "1-chloro-2,3-epoxypropane", "chloromethyloxirane", "epichlorhydrin")),
    ("maleic_anhydride", "Maleic anhydride", "108-31-6", ("maleic anhydride", "furan-2,5-dione", "2,5-furandione", "maleic acid anhydride")),
    ("phenol", "Phenol", "108-95-2", ("phenol", "hydroxybenzene", "benzenol", "carbolic acid", "monohydroxybenzene")),
    ("toluene", "Toluene", "108-88-3", ("toluene", "methylbenzene", "toluol", "phenylmethane")),
    ("benzene", "Benzene", "71-43-2", ("benzene", "benzol", "cyclohexatriene")),
    ("acetone", "Acetone", "67-64-1", ("acetone", "propan-2-one", "2-propanone", "propanone", "dimethyl ketone")),
    ("chloroform", "Chloroform", "67-66-3", ("chloroform", "trichloromethane", "trichloroform")),
    ("dichloromethane", "Dichloromethane", "75-09-2", ("dichloromethane", "methylene chloride", "methylene dichloride", "dcm")),
    ("tetrahydrofuran", "Tetrahydrofuran", "109-99-9", ("tetrahydrofuran", "oxolane", "thf", "1,4-epoxybutane", "diethylene oxide")),
    ("dmf", "N,N-Dimethylformamide", "68-12-2", ("n,n-dimethylformamide", "dimethylformamide", "dmf", "formic acid dimethylamide")),
    ("dmac", "N,N-Dimethylacetamide", "127-19-5", ("n,n-dimethylacetamide", "dimethylacetamide", "dmac", "acetic acid dimethylamide")),
    ("nmp", "N-Methyl-2-pyrrolidone", "872-50-4", ("n-methyl-2-pyrrolidone", "1-methyl-2-pyrrolidone", "n-methylpyrrolidone", "nmp", "1-methylpyrrolidin-2-one")),
    ("dmso", "Dimethyl sulfoxide", "67-68-5", ("dimethyl sulfoxide", "methyl sulfoxide", "dmso", "sulfinylbismethane")),
    ("methanol", "Methanol", "67-56-1", ("methanol", "methyl alcohol", "wood alcohol", "carbinol")),
    ("ethanol", "Ethanol", "64-17-5", ("ethanol", "ethyl alcohol", "grain alcohol", "ethyl hydroxide")),
    ("isopropanol", "Isopropanol", "67-63-0", ("isopropanol", "propan-2-ol", "2-propanol", "isopropyl alcohol", "ipa")),
    ("ethyl_acetate", "Ethyl acetate", "141-78-6", ("ethyl acetate", "ethyl ethanoate", "acetic acid ethyl ester", "etoac")),
    ("1_4_dioxane", "1,4-Dioxane", "123-91-1", ("1,4-dioxane", "1,4-dioxan", "p-dioxane", "1,4-diethylene dioxide")),
    ("cyclohexane", "Cyclohexane", "110-82-7", ("cyclohexane", "hexahydrobenzene", "cyclohexan")),
    ("anisole", "Anisole", "100-66-3", ("anisole", "methoxybenzene", "methyl phenyl ether", "phenol methyl ether")),
    ("pyridine", "Pyridine", "110-86-1", ("pyridine", "azine", "pyridin")),
    ("triethylamine", "Triethylamine", "121-44-8", ("triethylamine", "n,n-diethylethanamine", "triethyl amine")),
    ("acetonitrile", "Acetonitrile", "75-05-8", ("acetonitrile", "ethanenitrile", "methyl cyanide", "ethanonitrile", "acn")),
    ("acetic_acid", "Acetic acid", "64-19-7", ("acetic acid", "ethanoic acid", "glacial acetic acid", "ethylic acid", "vinegar acid")),
    ("lactic_acid", "Lactic acid", "50-21-5", ("lactic acid", "2-hydroxypropanoic acid", "2-hydroxypropionic acid")),
]

# nmrshiftdb2 SDfiles carry no name/CAS fields (the molfile title line is blank),
# so compounds are identified by the first 14 characters of their InChIKey
# (the connectivity block — stable across stereoisomers/protonation states).
# Values taken from PubChem for the names in NMR_SDF_TARGETS.
NMR_INCHI_KEYS: dict[str, str] = {
    "styrene": "PPBRXRYQALVLMV",
    "methyl_methacrylate": "VVQNEPGJFQJSBK",
    "vinyl_acetate": "XTXRWKRVRITETP",
    "acrylonitrile": "NLHHRLWOUZZQLW",
    "acrylic_acid": "NIXOWILDQLNWCW",
    "methacrylic_acid": "CERQOIWHTDAKMF",
    "ethylene_glycol": "LYCAIKOWRPUZTN",
    "epsilon_caprolactam": "JBKVHLHDHHXQEQ",
    "bisphenol_a": "IISBACLAFKSPIT",
    "epichlorohydrin": "BRLQWZUYTZBJKN",
    "maleic_anhydride": "FPYJFEHAWHCUMM",
    "phenol": "ISWSIDIOOBJBQZ",
    "toluene": "YXFVVABEGXRONW",
    "benzene": "UHOVQNZJYSORNB",
    "acetone": "CSCPPACGZOOCGX",
    "chloroform": "HEDRZPFGACZZDS",
    "dichloromethane": "YMWUJEATGCHHMB",
    "tetrahydrofuran": "WYURNTSHIVDZCO",
    "dmf": "ZMXDDKWLCZADIW",
    "dmac": "FXHOOIRPVKKKFG",
    "nmp": "SECXISVLQFMRJM",
    "dmso": "IAZDPXIOMUYVGZ",
    "methanol": "OKKJLVBELUTLKV",
    "ethanol": "LFQSCWFLJHTTHZ",
    "isopropanol": "KFZMGEQAYNKOFK",
    "ethyl_acetate": "XEKOWRVHYACXOJ",
    "1_4_dioxane": "RYHBNJHYFVUHQT",
    "cyclohexane": "XDTMQSROBMDMFD",
    "anisole": "RDOXTESZEPMUJZ",
    "pyridine": "JUJWROOIHBZHMG",
    "triethylamine": "ZMANZCXQSJIPKH",
    "acetonitrile": "WEVYAHXRMPXWCK",
    "acetic_acid": "QTBSBXVTEAMEQO",
    "lactic_acid": "JVTAAEKCZFNVCJ",
}


def http_get(url: str, timeout: int = 25) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def http_get_retry(url: str, attempts: int = 4, timeout: int = 25) -> bytes:
    """Fetch with exponential backoff — NIST returns HTTP 429 when hammered."""
    delay = 6.0
    last: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            return http_get(url, timeout=timeout)
        except Exception as exc:  # noqa: BLE001 — urllib raises many HTTP variants
            last = exc
            if attempt == attempts:
                break
            time.sleep(delay)
            delay *= 2.0
    assert last is not None
    raise last


def download_nist_pack(out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== [1/3] Downloading NIST WebBook Experimental FTIR (.jdx) -> {out_dir} ===")
    ok, skip, fail = 0, 0, 0
    for slug, cas in NIST_TARGETS:
        dest = out_dir / f"{slug}_nist.jdx"
        if dest.is_file() and dest.stat().st_size > 400:
            skip += 1
            continue
        cas_clean = cas.replace("-", "")
        # Direct NIST JCAMP IR endpoint (Index=0)
        jdx_url = f"https://webbook.nist.gov/cgi/cbook.cgi?JCAMP=C{cas_clean}&Index=0&Type=IR"
        try:
            data = http_get_retry(jdx_url, attempts=4, timeout=20)
            text = data.decode("latin-1", errors="ignore")
            if "##TITLE=" in text and ("##XYDATA=" in text or "##XYPOINTS=" in text):
                dest.write_bytes(data)
                print(f"  ✓ {slug:<22} ({len(data)//1024} KB)")
                ok += 1
            else:
                print(f"  - {slug:<22} (no IR JCAMP on NIST)")
                fail += 1
        except Exception as exc:
            print(f"  ✗ {slug:<22} ({exc})")
            fail += 1
        time.sleep(1.2)  # Be polite to NIST servers
    print(f"NIST summary: {ok} downloaded, {skip} already cached, {fail} unavailable.")


def _write_reference_csv(dest: Path, xs: np.ndarray, ys: np.ndarray, source: str) -> None:
    order = np.argsort(xs)[::-1]  # Descending, as FTIR instruments export
    with open(dest, "w", encoding="utf-8") as fh:
        fh.write(f"# Source: {source}\n")
        fh.write("Wavenumber (cm-1),Transmittance (%T)\n")
        for x, y in zip(xs[order], ys[order]):
            fh.write(f"{x:.2f},{y:.4f}\n")


def _extract_zenodo_csv_to_master(csv_bytes: bytes, out_dir: Path, tag: str) -> int:
    """
    Parses an FTIR-Plastics CSV and writes a median master spectrum plus the two
    most representative raw spectra per polymer.

    Two layouts occur in the dataset:
      * repeated `Data(x),Data(y)` column pairs — one spectrum per row (actual
        FTIR-PLASTIC_c8.csv export);
      * one numeric wavenumber column per point (older/derived exports).
    """
    reader = csv.reader(io.StringIO(csv_bytes.decode("utf-8", errors="ignore")))
    header = next(reader, None)
    if not header:
        return 0

    clean = [h.strip() for h in header]

    wn_cols: list[tuple[int, float]] = []
    for i, h in enumerate(clean):
        try:
            val = float(h.replace("cm-1", "").strip())
        except ValueError:
            continue
        if 350.0 <= val <= 4500.0:
            wn_cols.append((i, val))

    xy_pairs = [
        (i, i + 1)
        for i in range(len(clean) - 1)
        if clean[i].lower() == "data(x)" and clean[i + 1].lower() == "data(y)"
    ]

    if len(wn_cols) >= 50:
        layout = "matrix"
    elif len(xy_pairs) >= 50:
        layout = "eav"
    else:
        return 0

    poly_col = next(
        (i for i, h in enumerate(clean) if h.lower() in ("polymer", "material", "class", "label")),
        min(1, len(clean) - 1),
    )

    by_poly: dict[str, list[tuple[np.ndarray, np.ndarray]]] = {}
    for row in reader:
        if len(row) <= poly_col:
            continue
        poly = re.sub(r"[^A-Za-z0-9_-]+", "_", row[poly_col].strip()).upper()
        if not poly or poly in ("NAN", "NONE"):
            continue

        try:
            if layout == "matrix":
                if len(row) <= max(i for i, _ in wn_cols):
                    continue
                xs = np.array([w for _, w in wn_cols], dtype=np.float64)
                ys = np.array([float(row[i]) for i, _ in wn_cols], dtype=np.float64)
            else:
                if len(row) <= max(y for _, y in xy_pairs):
                    continue
                xs = np.array(
                    [float(row[i]) if row[i].strip() else np.nan for i, _ in xy_pairs],
                    dtype=np.float64,
                )
                ys = np.array(
                    [float(row[j]) if row[j].strip() else np.nan for _, j in xy_pairs],
                    dtype=np.float64,
                )
        except ValueError:
            continue

        mask = np.isfinite(xs) & np.isfinite(ys)
        if mask.sum() < 100:
            continue
        xs, ys = xs[mask], ys[mask]
        order = np.argsort(xs)
        by_poly.setdefault(poly, []).append((xs[order], ys[order]))

    written = 0
    for poly, spectra in by_poly.items():
        ref_x = max(spectra, key=lambda s: s[0].size)[0]
        stack = np.vstack([np.interp(ref_x, x, y) for x, y in spectra])
        median_y = np.median(stack, axis=0)

        distances = np.mean(np.abs(stack - median_y), axis=1)
        picks = list(np.argsort(distances)[:2])

        source = f"Zenodo FTIR-Plastics (10.5281/zenodo.10736649) | N={len(spectra)}"
        _write_reference_csv(out_dir / f"{poly}_real_{tag}_median.csv", ref_x, median_y, source)
        written += 1
        print(f"  ✓ {poly}_real_{tag}_median.csv (median of {len(spectra)} real Jasco spectra)")

        for rank, idx in enumerate(picks, start=1):
            dest = out_dir / f"{poly}_real_{tag}_sample{rank}.csv"
            _write_reference_csv(dest, ref_x, stack[idx], f"{source} | sample {idx + 1}")
            written += 1

    return written


def download_zenodo_ftir_plastics(out_dir: Path, max_mb: float = 300.0) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== [2/3] Querying Zenodo FTIR-Plastics -> {out_dir} ===")

    # 11176740 is the metadata-only version; 10736650 carries FTIR-PLASTIC.zip.
    files: list[dict] = []
    for record in ("10736650", "11176740"):
        try:
            meta_raw = http_get(f"https://zenodo.org/api/records/{record}", timeout=25)
            meta = json.loads(meta_raw.decode("utf-8"))
        except Exception as exc:
            print(f"  ✗ Could not reach Zenodo record {record}: {exc}")
            continue
        files = [f for f in meta.get("files", []) if float(f.get("size", 0)) > 500_000]
        if files:
            print(f"  record {record}: {len(files)} data file(s)")
            break

    if not files:
        print("  - No downloadable data file found on Zenodo.")
        return

    candidates = []
    for f in files:
        key = str(f.get("key", ""))
        size_mb = float(f.get("size", 0)) / (1024 * 1024)
        url = f.get("links", {}).get("self")
        if url and (key.lower().endswith(".csv") or key.lower().endswith(".zip")) and size_mb <= max_mb:
            candidates.append((size_mb, key, url))
    candidates.sort()

    if not candidates:
        print(f"  - No CSV/ZIP file under {max_mb:.0f} MB found on Zenodo record.")
        return

    # Download the cleanest smallest CSV/ZIP dataset first
    size_mb, key, url = candidates[0]
    cache_dir = out_dir.parent / ".cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / key

    if cache_path.is_file() and cache_path.stat().st_size > 1_000_000:
        print(f"  Reusing cached {cache_path.name} ({cache_path.stat().st_size / 1048576:.1f} MB)")
        blob = cache_path.read_bytes()
    else:
        print(f"  Downloading {key} ({size_mb:.1f} MB)...")
        try:
            blob = http_get_retry(url, attempts=3, timeout=max(180, int(size_mb * 2)))
        except Exception as exc:
            print(f"  ✗ Download failed: {exc}")
            return
        cache_path.write_bytes(blob)
        print(f"  Cached -> {cache_path} ({len(blob) / 1048576:.1f} MB)")

    total_written = 0
    if key.lower().endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(blob)) as zf:
            members = [m for m in zf.namelist() if m.lower().endswith(".csv")]
            print(f"  {len(members)} CSV member(s) in archive")
            for member in members:
                data = zf.read(member)
                tag = Path(member).stem.lower()
                total_written += _extract_zenodo_csv_to_master(data, out_dir, tag)
    else:
        tag = Path(key).stem.lower()
        total_written += _extract_zenodo_csv_to_master(blob, out_dir, tag)

    print(f"Zenodo summary: {total_written} real polymer master spectra extracted.")


# ---------------------------------------------------------------------------
# nmrshiftdb2: real experimental shifts from the official SDF dump (offline)
# ---------------------------------------------------------------------------

NMR_SDF_URL = (
    "https://sourceforge.net/projects/nmrshiftdb2/files/data/"
    "nmrshiftdb2withsignals.sd/download"
)
NMR_SDF_NAME = "nmrshiftdb2withsignals.sd"
NMR_SDF_MIN_BYTES = 10 * 1024 * 1024
NMR_INDEX_NAME = "nmrshiftdb2withsignals.idx.json.gz"
NMR_QUERY_CACHE_NAME = "nmr_query_cache.json"

# nucleus -> (lo ppm, hi ppm, points, Lorentzian half-width in ppm)
NMR_WINDOWS: dict[str, tuple[float, float, int, float]] = {
    "1h": (0.0, 12.0, 3000, 0.015),
    "13c": (0.0, 220.0, 4000, 0.20),
}

_SOLVENT_PREFERENCE = (
    "chloroform",
    "benzene-d6",
    "dimethylsulphoxide",
    "dimethylsulfoxide",
    "methanol-d4",
    "water",
    "acetone-d6",
    "acetonitrile",
)


def _slug(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_").lower()
    return slug or "nmr"


def _sourceforge_mirror_url(url: str, timeout: int = 90) -> str:
    """SourceForge answers /download with an HTML interstitial + meta refresh.

    Extract the real mirror URL from it; direct file URLs pass straight through.
    """
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        ctype = (resp.headers.get("Content-Type") or "").lower()
        if "html" not in ctype:
            return url
        page = resp.read().decode("utf-8", "ignore")
    m = re.search(
        r"""<meta[^>]+http-equiv=["']refresh["'][^>]+content=["'][^"']*url=([^"']+)["']""",
        page,
        re.I,
    )
    if not m:
        m = re.search(r'data-release-url="([^"]+)"', page)
    if not m:
        raise RuntimeError("SourceForge interstitial contained no mirror link")
    return html.unescape(m.group(1))


def _cache_nmr_sdf(cache: Path) -> Path:
    """Download the 159 MB nmrshiftdb2 SDfile once and reuse it forever."""
    if cache.is_file() and cache.stat().st_size > NMR_SDF_MIN_BYTES:
        print(f"  Cached SDF -> {cache} ({cache.stat().st_size / 1e6:.1f} MB)")
        return cache

    cache.parent.mkdir(parents=True, exist_ok=True)
    part = cache.with_name(cache.name + ".part")
    delay = 6.0
    last: Exception | None = None
    for attempt in range(1, 5):
        try:
            direct = _sourceforge_mirror_url(NMR_SDF_URL)
            print(f"  Downloading {NMR_SDF_NAME} (~159 MB), attempt {attempt} ...")
            req = urllib.request.Request(direct, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=240) as resp, part.open("wb") as fh:
                total = int(resp.headers.get("Content-Length") or 0)
                done = 0
                while True:
                    chunk = resp.read(1 << 16)
                    if not chunk:
                        break
                    fh.write(chunk)
                    done += len(chunk)
            if total and done < total:
                raise IOError(f"truncated download ({done}/{total} bytes)")
            if done < NMR_SDF_MIN_BYTES:
                raise IOError(f"download too small ({done} bytes)")
            part.replace(cache)
            print(f"  Cached SDF -> {cache} ({done / 1e6:.1f} MB)")
            return cache
        except Exception as exc:  # noqa: BLE001 — urllib raises many variants
            last = exc
            part.unlink(missing_ok=True)
            if attempt == 4:
                break
            print(f"  ! {exc} — retrying in {delay:.0f}s")
            time.sleep(delay)
            delay *= 2.0
    assert last is not None
    raise last


def _iter_sdf_records(path: Path):
    """Yield records one at a time; never holds the whole 159 MB SDF in RAM."""
    buf: list[str] = []
    with path.open("r", encoding="latin-1", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("$$$$"):
                if buf:
                    yield buf
                    buf = []
                continue
            buf.append(line)
    if buf:
        yield buf


_SDF_FIELD_RE = re.compile(r"^>\s*<([^>]+)>")


def _sdf_fields(lines: list[str]) -> dict[str, str]:
    """Parse the `> <Name>` data items of a single SDF record."""
    fields: dict[str, str] = {}
    i, n = 0, len(lines)
    while i < n:
        m = _SDF_FIELD_RE.match(lines[i].strip())
        if not m:
            i += 1
            continue
        key = m.group(1).strip()
        i += 1
        vals: list[str] = []
        while i < n and lines[i].strip() and not lines[i].lstrip().startswith(">"):
            vals.append(lines[i].strip())
            i += 1
        fields.setdefault(key, " ".join(vals))
    return fields


_INDEXED_VALUE_RE = re.compile(r"(?:^|\s)(\d+):(.+?)(?=\s+\d+:|$)")


def _indexed_field(text: str, index: str) -> str:
    """Pick one entry from a multi-value field like `0:CDCl3 3:CCl4 ...`.

    Every spectrum in a record has its own solvent/temperature/field strength,
    stored as `index:value` pairs in a single data item.
    """
    if not text:
        return ""
    entries = {m.group(1): m.group(2).strip() for m in _INDEXED_VALUE_RE.finditer(text)}
    if index in entries:
        return entries[index]
    return next(iter(entries.values()), text.strip())


def _parse_shift_signals(value: str, lo: float, hi: float) -> list[tuple[float, int]]:
    """Parse an nmrshiftdb2 signal list `shift;intensity;atom|...`.

    Returns [(ppm, occurrences)]; occurrences approximate the integration
    because every equivalent atom contributes its own entry.
    """
    counts: dict[float, int] = {}
    for token in value.split("|"):
        head = token.split(";", 1)[0].strip()
        if not head:
            continue
        try:
            shift = float(head)
        except ValueError:
            continue  # NaN / textual placeholder
        if not math.isfinite(shift) or not (lo <= shift <= hi):
            continue
        key = round(shift, 4)
        counts[key] = counts.get(key, 0) + 1
    return sorted(counts.items())


def _normalise_name(text: str) -> str:
    text = text.strip().lower()
    text = re.sub(r"\((?:[rs])\)-", "", text)
    text = re.sub(r"\([^)]*\)", " ", text)
    text = text.replace("α", "alpha").replace("β", "beta").replace("γ", "gamma")
    text = text.replace("’", "'").replace("′", "'")
    text = re.sub(r"[\s\-_,;:]+", " ", text)
    return text.strip()


def _nmr_target_indexes() -> tuple[dict[str, str], dict[str, str], dict[str, str]]:
    by_key: dict[str, str] = {}
    by_cas: dict[str, str] = {}
    by_name: dict[str, str] = {}
    for slug, display, cas, aliases in NMR_SDF_TARGETS:
        key = NMR_INCHI_KEYS.get(slug)
        if key:
            by_key[key.upper()] = slug
        by_cas[cas] = slug
        for name in (display, *aliases):
            by_name[_normalise_name(name)] = slug
    return by_key, by_cas, by_name


def _match_target(
    fields: dict[str, str],
    title: str,
    by_key: dict[str, str],
    by_cas: dict[str, str],
    by_name: dict[str, str],
) -> str | None:
    """Match a record by InChIKey, then CAS, then exact case-insensitive name.

    The official nmrshiftdb2 SDfile has a blank molfile title and no Name/CAS
    data items, so the InChIKey (always present) is the primary key; the name
    and CAS branches keep the extractor working with richer SDfiles.
    """
    key = fields.get("INChI key", "").strip().upper()
    if len(key) >= 14 and key[:14] in by_key:
        return by_key[key[:14]]
    for field in ("CAS", "CAS Registry Number", "CAS Number"):
        cas = fields.get(field, "").strip()
        if cas and cas in by_cas:
            return by_cas[cas]
    names = [title]
    names += [fields.get(f, "") for f in ("Name", "CanonicalName", "Systematic Name", "CAS Name")]
    for name in names:
        if name:
            slug = by_name.get(_normalise_name(name))
            if slug:
                return slug
    return None


def _record_spectra(
    fields: dict[str, str],
) -> list[tuple[str, str, str, list[tuple[float, int]]]]:
    """Return [(nucleus, index, raw_field, signals)] for one SDF record."""
    out: list[tuple[str, str, str, list[tuple[float, int]]]] = []
    for key, value in fields.items():
        m = re.match(r"^Spectrum\s+(1H|13C)\s+(\d+)$", key)
        if not m:
            continue
        nucleus = "1h" if m.group(1) == "1H" else "13c"
        lo, hi = (NMR_WINDOWS[nucleus][0] - 1.0, NMR_WINDOWS[nucleus][1] + 4.0)
        signals = _parse_shift_signals(value, lo, hi)
        if signals:
            out.append((nucleus, m.group(2), key, signals))
    return out


def _solvent_rank(solvent: str) -> int:
    low = solvent.lower()
    for rank, fragment in enumerate(_SOLVENT_PREFERENCE):
        if fragment in low:
            return rank
    return len(_SOLVENT_PREFERENCE)


def _best_spectrum(
    best: dict[tuple[str, str], dict],
    slug: str,
    nucleus: str,
    rec_id: str,
    solvent: str,
    field: str,
    signals: list[tuple[float, int]],
) -> None:
    """Keep the spectrum with the most resolved shifts (CDCl3 wins ties)."""
    score = (
        len(signals),
        -_solvent_rank(solvent),
        -int(rec_id) if rec_id.isdigit() else 0,
    )
    current = best.get((slug, nucleus))
    if current is None or score > current["score"]:
        best[(slug, nucleus)] = {
            "score": score,
            "id": rec_id,
            "solvent": solvent,
            "field": field,
            "signals": signals,
        }


def _lorentzian_trace(
    signals: list[tuple[float, int]], nucleus: str
) -> tuple[np.ndarray, np.ndarray]:
    lo, hi, npts, gamma = NMR_WINDOWS[nucleus]
    xs = np.linspace(lo, hi, npts)
    ys = np.zeros_like(xs)
    for shift, weight in signals:
        ys += weight * (gamma * gamma) / ((xs - shift) ** 2 + gamma * gamma)
    peak = float(ys.max()) if ys.size else 0.0
    if peak > 0.0:
        ys = ys / peak * 100.0
    return xs, ys


def _write_nmr_csv(
    dest: Path,
    signals: list[tuple[float, int]],
    nucleus: str,
    rec_id: str,
    solvent: str,
) -> None:
    xs, ys = _lorentzian_trace(signals, nucleus)
    lo, hi, npts, gamma = NMR_WINDOWS[nucleus]
    label = "1H" if nucleus == "1h" else "13C"
    shift_text = " ".join(f"{s:g}" for s, _ in signals)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as fh:
        fh.write(
            "# Source: nmrshiftdb2 (nmrshiftdb2withsignals.sd)"
            f" | ID={rec_id} | shifts={shift_text}\n"
        )
        fh.write(
            f"# {label} NMR (ppm) | solvent={solvent or 'n/a'} | "
            f"Lorentzian gamma={gamma:g} ppm | {npts} points\n"
        )
        fh.write("ppm,Intensity\n")
        # NMR spectra are conventionally plotted right-to-left (descending ppm).
        for x, y in zip(xs[::-1].tolist(), ys[::-1].tolist()):
            fh.write(f"{x:.4f},{y:.4f}\n")


def download_nmr_pack(out_dir: Path, root: Path | None = None) -> None:
    root = root or out_dir.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\n=== [3/3] Extracting real nmrshiftdb2 1H/13C shifts -> {out_dir} ===")
    sdf = _cache_nmr_sdf(root / ".cache" / NMR_SDF_NAME)
    by_key, by_cas, by_name = _nmr_target_indexes()

    best: dict[tuple[str, str], dict] = {}
    records = matched = 0
    started = time.time()
    for lines in _iter_sdf_records(sdf):
        records += 1
        fields = _sdf_fields(lines)
        title = lines[0].strip() if lines else ""
        slug = _match_target(fields, title, by_key, by_cas, by_name)
        if slug is None:
            continue
        matched += 1
        rec_id = fields.get("nmrshiftdb2 ID", "").strip()
        for nucleus, idx, field, signals in _record_spectra(fields):
            solvent = _indexed_field(fields.get("Solvent", ""), idx)
            _best_spectrum(best, slug, nucleus, rec_id, solvent, field, signals)

    written = 0
    for (slug, nucleus), entry in sorted(best.items()):
        dest = out_dir / f"{slug}_{nucleus}.csv"
        _write_nmr_csv(dest, entry["signals"], nucleus, entry["id"], entry["solvent"])
        written += 1
        label = "1H" if nucleus == "1h" else "13C"
        print(
            f"  ✓ {slug:<22} {label:<3} {len(entry['signals']):>3} signals  "
            f"ID={entry['id'] or '?':<9} {entry['solvent'] or 'solvent n/a'}"
        )

    missing = sorted({t[0] for t in NMR_SDF_TARGETS} - {slug for slug, _ in best})
    if missing:
        print(f"  ! no usable shifts for: {', '.join(missing)}")
    compounds = len({slug for slug, _ in best})
    print(
        f"nmrshiftdb2 summary: {written} spectra for {compounds} compounds "
        f"({matched} matching records of {records}, {time.time() - started:.1f}s)."
    )


def _nmr_index_path(cache_dir: Path) -> Path:
    return cache_dir / NMR_INDEX_NAME


def _build_nmr_index(sdf: Path, index_path: Path) -> dict[str, list[list[int]]]:
    """Map InChIKey -> [(byte offset, length)] so lookups never rescan 159 MB."""
    print(f"  Indexing {sdf.name} once (~10 s) for instant offline lookups ...")
    index: dict[str, list[list[int]]] = {}
    records = 0
    with sdf.open("rb") as fh:
        start = fh.tell()
        key = ""
        while True:
            line = fh.readline()
            if not line:
                break
            if line.startswith(b"$$$$"):
                end = fh.tell()
                if key:
                    index.setdefault(key, []).append([start, end - start])
                    records += 1
                key = ""
                start = fh.tell()
            elif line.startswith(b"> <INChI key>"):
                value = fh.readline().decode("latin-1", "replace").strip()
                key = value[:14].upper()
    index_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(index_path, "wt", encoding="utf-8") as fh:
        json.dump({"size": sdf.stat().st_size, "records": records, "index": index}, fh)
    print(f"  Indexed {records} records -> {index_path.name}")
    return index


def _load_nmr_index(sdf: Path, cache_dir: Path) -> dict[str, list[list[int]]]:
    path = _nmr_index_path(cache_dir)
    try:
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            blob = json.load(fh)
        if blob.get("size") == sdf.stat().st_size and blob.get("records"):
            return blob["index"]
    except Exception:  # noqa: BLE001 — rebuild on any cache corruption
        pass
    return _build_nmr_index(sdf, path)


def _load_query_cache(path: Path) -> dict[str, str]:
    try:
        blob = json.loads(path.read_text(encoding="utf-8"))
        return {str(k): str(v) for k, v in blob.items()}
    except Exception:  # noqa: BLE001
        return {}


def _pubchem_inchikey(query: str, timeout: int = 30) -> str | None:
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
        f"{urllib.parse.quote(query)}/property/InChIKey/JSON"
    )
    try:
        payload = json.loads(http_get(url, timeout=timeout))
        key = payload["PropertyTable"]["Properties"][0]["InChIKey"]
        return str(key)[:14].upper()
    except Exception:  # noqa: BLE001 — offline or unknown compound
        return None


def _resolve_nmr_query(
    query: str,
    cache_dir: Path,
    by_key: dict[str, str],
    by_cas: dict[str, str],
    by_name: dict[str, str],
) -> tuple[str, str]:
    """Resolve a name/CAS/InChIKey to a 14-char InChIKey block."""
    text = query.strip()
    compact = re.sub(r"[^A-Za-z0-9]", "", text).upper()
    if re.fullmatch(r"[A-Z]{14}[A-Z0-9]{0,11}", compact):
        return compact[:14], text

    slug = by_cas.get(text) or by_name.get(_normalise_name(text))
    if slug and slug in NMR_INCHI_KEYS:
        return NMR_INCHI_KEYS[slug], slug

    cache_path = cache_dir / NMR_QUERY_CACHE_NAME
    cache = _load_query_cache(cache_path)
    key = cache.get(_normalise_name(text)) or cache.get(text.lower())
    if not key:
        key = _pubchem_inchikey(text)
        if key:
            cache[_normalise_name(text)] = key
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(json.dumps(cache, indent=0, sort_keys=True), encoding="utf-8")
            print(f"  resolved {text!r} -> InChIKey {key} (cached for offline reuse)")
    return (key, text) if key else ("", text)


def search_nmr(out_dir: Path, query: str, root: Path | None = None) -> int:
    """Extract any compound's real 1H/13C NMR from the cached nmrshiftdb2 SDF."""
    root = root or out_dir.parent
    cache_dir = root / ".cache"
    print(f"\n=== nmrshiftdb2 offline search: {query!r} ===")
    sdf = _cache_nmr_sdf(cache_dir / NMR_SDF_NAME)
    by_key, by_cas, by_name = _nmr_target_indexes()
    key, label = _resolve_nmr_query(query, cache_dir, by_key, by_cas, by_name)
    if not key:
        print(
            "ERROR: could not resolve that query. Pass an InChIKey "
            "(e.g. PPBRXRYQALVLMV-UHFFFAOYSA-N), a CAS number or a name; "
            "names need one online PubChem lookup, then they work offline.",
            file=sys.stderr,
        )
        return 1

    index = _load_nmr_index(sdf, cache_dir)
    offsets = index.get(key)
    if not offsets:
        print(f"ERROR: InChIKey {key} is not present in {NMR_SDF_NAME}.", file=sys.stderr)
        return 1

    slug = by_key.get(key) or _slug(label)
    best: dict[tuple[str, str], dict] = {}
    with sdf.open("rb") as fh:
        for offset, length in offsets:
            fh.seek(offset)
            raw = fh.read(length).decode("latin-1", "replace")
            fields = _sdf_fields([ln.rstrip("\r\n") for ln in raw.splitlines()])
            rec_id = fields.get("nmrshiftdb2 ID", "").strip()
            for nucleus, idx, field, signals in _record_spectra(fields):
                solvent = _indexed_field(fields.get("Solvent", ""), idx)
                _best_spectrum(best, slug, nucleus, rec_id, solvent, field, signals)

    written = 0
    for nucleus in ("1h", "13c"):
        entry = best.get((slug, nucleus))
        if entry is None:
            continue
        dest = out_dir / f"{slug}_{nucleus}.csv"
        _write_nmr_csv(dest, entry["signals"], nucleus, entry["id"], entry["solvent"])
        written += 1
        label_n = "1H" if nucleus == "1h" else "13C"
        print(
            f"  ✓ {dest}  ({label_n}, {len(entry['signals'])} signals, "
            f"ID={entry['id'] or '?'}, {entry['solvent'] or 'solvent n/a'})"
        )
    if not written:
        print("ERROR: matched records carry no usable 1H/13C shift lists.", file=sys.stderr)
        return 1
    print(f"saved {written} spectrum file(s) in {out_dir}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Download real FTIR & NMR libraries for SPECTRA")
    parser.add_argument("--all", action="store_true", help="Download NIST + Zenodo + nmrshiftdb2")
    parser.add_argument("--nist", action="store_true", help="Download NIST WebBook FTIR pack (~2 MB)")
    parser.add_argument("--zenodo", action="store_true", help="Download Zenodo FTIR-Plastics pack")
    parser.add_argument(
        "--nmr",
        action="store_true",
        help="Extract the nmrshiftdb2 NMR pack from the cached SDF dump (159 MB, once)",
    )
    parser.add_argument(
        "--search-nmr",
        metavar="NAME_OR_CAS",
        default=None,
        help="Extract one compound's real 1H/13C NMR offline from the cached nmrshiftdb2 SDF",
    )
    parser.add_argument(
        "--max-mb",
        type=float,
        default=300.0,
        help="Maximum Zenodo archive size in MB (default 300)",
    )
    parser.add_argument(
        "--lib-root",
        default=str(Path(__file__).resolve().parent.parent / "library"),
        help="Target library directory",
    )
    args = parser.parse_args()

    root = Path(args.lib_root)

    if args.search_nmr:
        return search_nmr(root / "nmr", args.search_nmr, root=root)

    if not (args.all or args.nist or args.zenodo or args.nmr):
        args.all = True

    if args.all or args.nist:
        download_nist_pack(root / "ftir")
    if args.all or args.zenodo:
        download_zenodo_ftir_plastics(root / "ftir", max_mb=args.max_mb)
    if args.all or args.nmr:
        download_nmr_pack(root / "nmr", root=root)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
