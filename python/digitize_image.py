#!/usr/bin/env python3
"""
SPECTRA — Offline Spectrum Image Digitizer (PNG / JPG / WebP -> CSV)
Uses only numpy, scipy, and PIL/matplotlib (100% offline, zero OpenCV required).
Auto-detects the plot frame, strips gridlines & text annotations, traces the
spectrum curve with path continuity, and calibrates to physical cm⁻¹ / ppm.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

try:
    from PIL import Image
except Exception:
    Image = None

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg


def load_rgb(path: Path) -> np.ndarray:
    if Image is not None:
        im = Image.open(path).convert("RGBA")
        arr = np.asarray(im, dtype=np.float32) / 255.0
        alpha = arr[:, :, 3:4]
        return arr[:, :, :3] * alpha + (1.0 - alpha)
    raw = mpimg.imread(str(path))
    arr = raw.astype(np.float32)
    if arr.max() > 1.5:
        arr /= 255.0
    if arr.ndim == 2:
        return np.stack([arr, arr, arr], axis=-1)
    if arr.shape[2] >= 4:
        alpha = arr[:, :, 3:4]
        return arr[:, :, :3] * alpha + (1.0 - alpha)
    return arr[:, :, :3]


def digitize_image(
    img_path: Path,
    x_left: float,
    x_right: float,
    y_bottom: float,
    y_top: float,
    bg_mode: str = "auto",
) -> tuple[np.ndarray, np.ndarray, tuple[int, int, int, int]]:
    rgb = load_rgb(img_path)
    h, w, _ = rgb.shape
    if h < 40 or w < 60:
        raise ValueError("Image resolution too small to digitize.")

    gray = 0.299 * rgb[:, :, 0] + 0.587 * rgb[:, :, 1] + 0.114 * rgb[:, :, 2]
    median_lum = float(np.median(gray))
    dark_bg = (bg_mode == "dark") or (bg_mode == "auto" and median_lum < 0.45)

    fg = np.clip(gray - median_lum, 0.0, 1.0) if dark_bg else np.clip(median_lum - gray, 0.0, 1.0)
    sat = np.max(rgb, axis=2) - np.min(rgb, axis=2)

    # 1. Auto-detect plot axes frame via long horizontal/vertical ink projections
    ink = (fg > 0.16).astype(np.float32)
    row_proj = ink.mean(axis=1)
    col_proj = ink.mean(axis=0)

    top_c = np.where(row_proj[: int(h * 0.45)] > 0.32)[0]
    bot_c = np.where(row_proj[int(h * 0.55) :] > 0.32)[0] + int(h * 0.55)
    left_c = np.where(col_proj[: int(w * 0.40)] > 0.32)[0]
    right_c = np.where(col_proj[int(w * 0.60) :] > 0.32)[0] + int(w * 0.60)

    y0 = int(top_c[-1]) + 3 if len(top_c) > 0 else int(h * 0.08)
    y1 = int(bot_c[0]) - 3 if len(bot_c) > 0 else int(h * 0.86)
    x0 = int(left_c[-1]) + 3 if len(left_c) > 0 else int(w * 0.10)
    x1 = int(right_c[0]) - 3 if len(right_c) > 0 else int(w * 0.95)

    if (y1 - y0) < h * 0.25 or (x1 - x0) < w * 0.30:
        y0, y1 = int(h * 0.08), int(h * 0.86)
        x0, x1 = int(w * 0.10), int(w * 0.95)

    crop_fg = fg[y0:y1, x0:x1]
    crop_sat = sat[y0:y1, x0:x1]
    ch, cw = crop_fg.shape

    # 2. Separate spectrum trace from light/gray gridlines
    has_color = float(np.percentile(crop_sat, 99.5)) > 0.20
    if has_color:
        trace_map = crop_sat * 1.45 + crop_fg * 0.45
    else:
        row_med = np.median(crop_fg, axis=1, keepdims=True)
        col_med = np.median(crop_fg, axis=0, keepdims=True)
        trace_map = np.clip(crop_fg - 0.85 * np.minimum(row_med, col_med), 0.0, 1.0)

    peak_ref = float(np.percentile(trace_map, 98.8))
    thresh = max(0.14, peak_ref * 0.42)

    # 3. Continuity-aware column-by-column curve extraction
    ys_raw = np.full(cw, np.nan, dtype=np.float64)
    prev_y: float | None = None
    for col in range(cw):
        col_vals = trace_map[:, col]
        mask = col_vals >= thresh
        if not np.any(mask):
            continue
        labeled, n_comp = ndimage.label(mask)
        best_y = None
        best_cost = 1e18
        for cid in range(1, n_comp + 1):
            idxs = np.where(labeled == cid)[0]
            if len(idxs) > ch * 0.40:
                continue
            weights = col_vals[idxs]
            cy = float(np.sum(idxs * weights) / max(float(np.sum(weights)), 1e-9))
            strength = float(np.max(weights))
            if prev_y is None:
                cost = -strength * 100.0
            else:
                cost = abs(cy - prev_y) - strength * 14.0
            if cost < best_cost:
                best_cost = cost
                best_y = cy
        if best_y is not None:
            ys_raw[col] = best_y
            prev_y = best_y

    valid = np.isfinite(ys_raw)
    if int(np.sum(valid)) < 12:
        raise ValueError("Could not detect a clear spectrum curve in this image.")

    cols = np.arange(cw, dtype=np.float64)
    ys_interp = np.interp(cols, cols[valid], ys_raw[valid])
    ys_clean = ndimage.median_filter(ys_interp, size=3)

    # Resample to clean uniform grid (up to 1600 points)
    n_out = min(max(cw, 400), 1600)
    u = np.linspace(0.0, 1.0, n_out)
    y_resamp = np.interp(u * (cw - 1), cols, ys_clean)

    phys_x = x_left + u * (x_right - x_left)
    phys_y = y_top - (y_resamp / max(ch - 1, 1)) * (y_top - y_bottom)

    order = np.argsort(phys_x)
    return phys_x[order], phys_y[order], (x0, y0, x1, y1)


def main() -> int:
    parser = argparse.ArgumentParser(description="SPECTRA offline spectrum image digitizer")
    parser.add_argument("--input", "-i", required=True)
    parser.add_argument("--output", "-o", required=True)
    parser.add_argument("--kind", default="ftir", choices=["ftir", "hnmr", "cnmr", "uvvis"])
    parser.add_argument("--y-unit", default="transmittance", choices=["transmittance", "absorbance", "intensity"])
    parser.add_argument("--x-left", type=float, default=4000.0)
    parser.add_argument("--x-right", type=float, default=400.0)
    parser.add_argument("--bg", default="auto", choices=["auto", "light", "dark"])
    args = parser.parse_args()

    y_bot, y_top = (0.0, 100.0) if args.y_unit == "transmittance" else (0.0, 1.0)
    try:
        xs, ys, box = digitize_image(
            Path(args.input),
            x_left=args.x_left,
            x_right=args.x_right,
            y_bottom=y_bot,
            y_top=y_top,
            bg_mode=args.bg,
        )
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        if args.kind == "uvvis":
            header = "Wavelength (nm),Absorbance" if args.y_unit == "absorbance" else "Wavelength (nm),Transmittance (%T)"
        else:
            header = (
                "Wavenumber (cm-1),Transmittance (%T)"
                if args.kind == "ftir" and args.y_unit == "transmittance"
                else "Wavenumber (cm-1),Absorbance"
                if args.kind == "ftir"
                else "Chemical Shift (ppm),Intensity"
            )
        with open(out_path, "w", encoding="utf-8") as fh:
            fh.write(header + "\n")
            # UV-Vis: wavelengths ascending (already sorted ascending by digitize_image)
            if args.kind == "uvvis":
                for x, y in zip(xs, ys):
                    fh.write(f"{x:.4f},{y:.4f}\n")
            else:
                for x, y in zip(xs[::-1], ys[::-1]):
                    fh.write(f"{x:.4f},{y:.4f}\n")
        print(f"OK:{len(xs)}:box={box}")
        return 0
    except Exception as exc:
        print(f"ERROR:{exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
