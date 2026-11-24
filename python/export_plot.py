#!/usr/bin/env python3
"""SPECTRA — publication-quality spectroscopy figure exporter.

Consumes the JSON payload produced by the Rust TUI / CLI and renders a
journal-ready figure with Matplotlib:

  * IUPAC-reversed X axis (FTIR 4000 → 400 cm⁻¹, NMR high → low ppm)
  * ACS / RSC / Everforest journal styling
  * Overlay or vertically offset (stacked) multi-trace comparisons
  * Peak annotations (pinned ★ and auto-detected)
  * PNG (600 dpi), SVG (vector) or PDF (LaTeX) output

Usage:
    python3 export_plot.py --input payload.json --output figure.png
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import textwrap

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import AutoMinorLocator, MaxNLocator  # noqa: E402

STYLES = {
    "acs": {
        "figure_facecolor": "#ffffff",
        "axes_facecolor": "#ffffff",
        "text": "#111111",
        "spine": "#111111",
        "grid": "#c8c8c8",
        "palette": ["#1f4e79", "#c00000", "#2e7d32", "#6a1b9a", "#e65100", "#00838f"],
        "figsize": (7.0, 4.3),
        "linewidth": 1.35,
        "font_size": 9.0,
        "spines": ("left", "bottom"),
        "grid_alpha": 0.35,
    },
    "rsc": {
        "figure_facecolor": "#ffffff",
        "axes_facecolor": "#ffffff",
        "text": "#101820",
        "spine": "#101820",
        "grid": "#d5d9de",
        "palette": ["#003f5c", "#bc5090", "#ffa600", "#58508d", "#ff6361", "#2f4b7c"],
        "figsize": (7.2, 4.4),
        "linewidth": 1.4,
        "font_size": 9.0,
        "spines": ("left", "bottom"),
        "grid_alpha": 0.30,
    },
    "everforest": {
        "figure_facecolor": "#2d353b",
        "axes_facecolor": "#2d353b",
        "text": "#d3c6aa",
        "spine": "#859289",
        "grid": "#475258",
        "palette": ["#a7c080", "#dbbc7f", "#7fbbb3", "#d699b6", "#e67e80", "#83c092"],
        "figsize": (7.2, 4.4),
        "linewidth": 1.45,
        "font_size": 9.0,
        "spines": ("left", "bottom", "top", "right"),
        "grid_alpha": 0.35,
    },
}

X_LABELS = {
    "ftir": r"Wavenumber (cm$^{-1}$)",
    "hnmr": r"Chemical shift $\delta$ (ppm)",
    "cnmr": r"Chemical shift $\delta$ (ppm)",
    "uvvis": r"Wavelength, $\lambda$ (nm)",
}

Y_LABELS = {
    "transmittance": "Transmittance (%)",
    "absorbance": "Absorbance (a.u.)",
    "intensity": "Intensity (a.u.)",
}


def wavelength_to_rgb(nm: float) -> tuple[float, float, float]:
    """Dan Bruton's piecewise-linear CIE wavelength -> linear RGB (380-750 nm).

    Mirrors ``src/kitty_plot.rs::wavelength_to_rgb`` so the TUI HD preview and the
    publication export agree on the visible-light rainbow strip.
    """
    if nm < 380.0 or nm > 750.0:
        return (0.0, 0.0, 0.0)

    w = float(nm)
    if w < 440.0:
        r, g, b = -(w - 440.0) / (440.0 - 380.0), 0.0, 1.0
    elif w < 490.0:
        r, g, b = 0.0, (w - 440.0) / (490.0 - 440.0), 1.0
    elif w < 510.0:
        r, g, b = 0.0, 1.0, -(w - 510.0) / (510.0 - 490.0)
    elif w < 580.0:
        r, g, b = (w - 510.0) / (580.0 - 510.0), 1.0, 0.0
    elif w < 645.0:
        r, g, b = 1.0, -(w - 645.0) / (645.0 - 580.0), 0.0
    else:
        r, g, b = 1.0, 0.0, 0.0

    if w < 420.0:
        factor = 0.3 + 0.7 * (w - 380.0) / (420.0 - 380.0)
    elif w > 700.0:
        factor = 0.3 + 0.7 * (750.0 - w) / (750.0 - 700.0)
    else:
        factor = 1.0

    gamma = 0.8
    return tuple(
        float(np.clip(c, 0.0, 1.0) * factor) ** gamma for c in (r, g, b)
    )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Render a SPECTRA JSON payload as a publication figure."
    )
    parser.add_argument("--input", required=True, help="Path to the JSON payload")
    parser.add_argument("--output", required=True, help="Output figure path (.png/.svg/.pdf)")
    return parser.parse_args(argv)


def load_payload(path):
    with open(path, "r", encoding="utf-8") as handle:
        payload = json.load(handle)

    traces = payload.get("traces") or []
    if not traces:
        raise ValueError("payload contains no traces")

    cleaned = []
    for trace in traces:
        x = np.asarray(trace.get("x") or [], dtype=float)
        y = np.asarray(trace.get("y") or [], dtype=float)
        if x.size < 2 or y.size < 2:
            continue
        n = min(x.size, y.size)
        x, y = x[:n], y[:n]
        mask = np.isfinite(x) & np.isfinite(y)
        x, y = x[mask], y[mask]
        if x.size < 2:
            continue
        order = np.argsort(x)
        cleaned.append(
            {
                "name": str(trace.get("name") or "Spectrum"),
                "x": x[order],
                "y": y[order],
            }
        )

    if not cleaned:
        raise ValueError("payload contains no usable numeric traces")

    payload["traces"] = cleaned
    return payload


def apply_style(style_name):
    cfg = STYLES.get(str(style_name).lower(), STYLES["acs"]).copy()

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": cfg["font_size"],
            "axes.linewidth": 0.9,
            "axes.edgecolor": cfg["spine"],
            "axes.labelcolor": cfg["text"],
            "axes.titlecolor": cfg["text"],
            "xtick.color": cfg["text"],
            "ytick.color": cfg["text"],
            "text.color": cfg["text"],
            "figure.facecolor": cfg["figure_facecolor"],
            "savefig.facecolor": cfg["figure_facecolor"],
            "legend.frameon": False,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    return cfg


def _wrap_peak_label(label, width=19, max_lines=2):
    lines = textwrap.wrap(label, width=width, break_long_words=False)
    if not lines:
        return label[:width]
    if len(lines) > max_lines:
        kept = lines[:max_lines]
        last = kept[-1]
        if len(last) >= width:
            last = last[: max(1, width - 1)].rstrip()
        kept[-1] = f"{last}…"
        return "\n".join(kept)
    return "\n".join(lines)


def _space_peaks_u(targets, lanes, half_widths, min_step=0.018, lane_pad=0.016, edge_pad=0.014):
    """Space normalized X positions in [0, 1] monotonically without lane overlap."""
    n = len(targets)
    if n == 0:
        return []

    lo = [hw + edge_pad for hw in half_widths]
    hi = [1.0 - hw - edge_pad for hw in half_widths]
    if n == 1:
        return [min(max(targets[0], lo[0]), hi[0])]

    prev_same = [-1] * n
    last_in_lane = {}
    for i, lane in enumerate(lanes):
        if lane in last_in_lane:
            prev_same[i] = last_in_lane[lane]
        last_in_lane[lane] = i

    next_same = [-1] * n
    for i, p in enumerate(prev_same):
        if p >= 0:
            next_same[p] = i

    fwd = [min(max(targets[i], lo[i]), hi[i]) for i in range(n)]
    for i in range(1, n):
        fwd[i] = max(fwd[i], fwd[i - 1] + min_step)
        p = prev_same[i]
        if p >= 0:
            fwd[i] = max(fwd[i], fwd[p] + half_widths[p] + half_widths[i] + lane_pad)
    if fwd[-1] > hi[-1]:
        fwd[-1] = hi[-1]
        for i in range(n - 2, -1, -1):
            fwd[i] = min(fwd[i], fwd[i + 1] - min_step)
            nx = next_same[i]
            if nx >= 0:
                fwd[i] = min(fwd[i], fwd[nx] - (half_widths[i] + half_widths[nx] + lane_pad))

    bwd = [min(max(targets[i], lo[i]), hi[i]) for i in range(n)]
    for i in range(n - 2, -1, -1):
        bwd[i] = min(bwd[i], bwd[i + 1] - min_step)
        nx = next_same[i]
        if nx >= 0:
            bwd[i] = min(bwd[i], bwd[nx] - (half_widths[i] + half_widths[nx] + lane_pad))
    if bwd[0] < lo[0]:
        bwd[0] = lo[0]
        for i in range(1, n):
            bwd[i] = max(bwd[i], bwd[i - 1] + min_step)
            p = prev_same[i]
            if p >= 0:
                bwd[i] = max(bwd[i], bwd[p] + half_widths[p] + half_widths[i] + lane_pad)

    return [min(max(0.5 * (f + b), l), h) for f, b, l, h in zip(fwd, bwd, lo, hi)]


def choose_x_limits(payload, traces):
    bounds = payload.get("x_bounds") or []
    if isinstance(bounds, (list, tuple)) and len(bounds) >= 2:
        try:
            lo, hi = float(bounds[0]), float(bounds[1])
            if np.isfinite(lo) and np.isfinite(hi) and lo != hi:
                return min(lo, hi), max(lo, hi)
        except (TypeError, ValueError):
            pass

    lo = min(float(t["x"][0]) for t in traces)
    hi = max(float(t["x"][-1]) for t in traces)
    if lo == hi:
        hi = lo + 1.0
    return min(lo, hi), max(lo, hi)


def render(payload):
    cfg = apply_style(payload.get("style", "acs"))
    traces = payload["traces"]
    mode = str(payload.get("mode", "ftir")).lower()
    y_unit = str(payload.get("y_unit", "intensity")).lower()
    stacked = bool(payload.get("stacked", False))
    stacked_multi = stacked and len(traces) > 1
    title = str(payload.get("title") or "").strip()

    try:
        active_trace = int(payload.get("active_trace", len(traces) - 1))
    except (TypeError, ValueError):
        active_trace = len(traces) - 1
    active_trace = max(0, min(active_trace, len(traces) - 1))

    # Global Y span is computed before offsets so stacked traces keep a
    # consistent separation regardless of individual trace amplitudes.
    all_y = np.concatenate([t["y"] for t in traces])
    y_min = float(np.nanmin(all_y))
    y_max = float(np.nanmax(all_y))
    y_span = max(y_max - y_min, 1e-9)

    offsets = [0.0] * len(traces)
    if stacked and len(traces) > 1:
        step = y_span * 1.12
        for i in range(len(traces)):
            offsets[i] = i * step

    shifted_max = y_max + (offsets[-1] if offsets else 0.0)
    data_span = max(shifted_max - y_min, y_span)
    pad = max(data_span * 0.06, y_span * 0.04)

    fig, ax = plt.subplots(figsize=cfg["figsize"], facecolor=cfg["figure_facecolor"])
    ax.set_facecolor(cfg["axes_facecolor"])

    x_lo, x_hi = choose_x_limits(payload, traces)

    for idx, trace in enumerate(traces):
        color = cfg["palette"][idx % len(cfg["palette"])]
        is_active = idx == active_trace
        lw = cfg["linewidth"] * (1.18 if is_active else 0.92)
        alpha = 1.0 if is_active or len(traces) == 1 else 0.88
        ax.plot(
            trace["x"],
            trace["y"] + offsets[idx],
            color=color,
            linewidth=lw,
            alpha=alpha,
            solid_capstyle="round",
            solid_joinstyle="round",
            label=trace["name"],
            zorder=3 + idx,
        )

    # UV-Vis: ascending wavelength (nm) left-to-right. Others: IUPAC reversed.
    if mode == "uvvis":
        ax.set_xlim(x_lo, x_hi)
    else:
        # IUPAC convention: X axis decreases from left to right.
        ax.set_xlim(x_hi, x_lo)

    ax.set_xlabel(X_LABELS.get(mode, X_LABELS["ftir"]), labelpad=6)
    y_label = Y_LABELS.get(y_unit, "Intensity (a.u.)")
    if stacked_multi:
        y_label = f"{y_label} (stacked)"
    ax.set_ylabel(y_label, labelpad=6)

    legend_extra = data_span * (0.08 + 0.035 * min(len(traces), 4)) if len(traces) > 1 else 0.0
    ax.set_ylim(y_min - pad, shifted_max + pad + legend_extra)

    # Baseline for NMR / absorbance traces
    if mode != "ftir" or y_unit == "absorbance":
        ax.axhline(0.0, color=cfg["grid"], linewidth=0.8, alpha=0.65, zorder=1)

    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, prune=None))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=7, prune=None))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))

    if stacked_multi:
        # Stacked traces are vertically offset for readability, so absolute
        # Y tick values (e.g. >100% transmittance) would be misleading.
        ax.set_yticks([])

    ax.grid(axis="y", color=cfg["grid"], linewidth=0.6, alpha=cfg["grid_alpha"], zorder=0)
    ax.tick_params(axis="both", which="major", direction="out", length=4.0, width=0.9, labelsize=8.2)
    ax.tick_params(axis="both", which="minor", direction="out", length=2.0, width=0.6)

    for side, spine in ax.spines.items():
        spine.set_visible(side in cfg["spines"])
        spine.set_color(cfg["spine"])
        spine.set_linewidth(0.9)

    if title:
        ax.set_title(title, fontsize=10.5, fontweight="bold", pad=10)

    # ------------------------------------------------------------------
    # Peak annotations: leaders point to peaks while labels occupy reserved
    # top/bottom lanes outside the trace area. This avoids text running over
    # the spectrum, even when several peak labels are close together.
    # ------------------------------------------------------------------
    peak_candidates = []
    for peak in payload.get("peaks") or []:
        try:
            px = float(peak.get("x"))
            py = float(peak.get("y")) + offsets[active_trace]
        except (TypeError, ValueError):
            continue
        if not (np.isfinite(px) and np.isfinite(py)):
            continue
        if px < x_lo - 1e-9 or px > x_hi + 1e-9:
            continue
        label = str(peak.get("label") or f"{px:.2f}").strip()
        if label:
            peak_candidates.append(
                (px, py, label, bool(peak.get("pinned", False)))
            )

    # Pinned annotations are intentional user choices, so keep them first.
    peak_candidates.sort(key=lambda item: not item[3])
    max_labels = 10
    visible_peaks = peak_candidates[:max_labels]
    if visible_peaks:
        # UV-Vis: ascending wavelength left-to-right. Others: reversed IUPAC.
        if mode == "uvvis":
            visible_peaks.sort(key=lambda item: item[0])  # ascending
        else:
            # Sort left-to-right on the reversed IUPAC X-axis before assigning
            # round-robin lanes so adjacent peaks land on different tiers and
            # leader lines stay nearly vertical instead of crossing.
            visible_peaks.sort(key=lambda item: item[0], reverse=True)
        wrapped_items = []
        for px, py, label, pinned in visible_peaks:
            wrapped = _wrap_peak_label(label, width=19, max_lines=2)
            wrapped_items.append((px, py, wrapped, pinned))

        n_vis = len(wrapped_items)
        lane_count = 1 if n_vis <= 2 else (2 if n_vis <= 4 else 3)
        lanes = [[] for _ in range(lane_count)]
        for idx, item in enumerate(wrapped_items):
            lanes[idx % lane_count].append(item)

        line_counts = [
            max((1 + item[2].count("\n") for item in lane), default=1)
            for lane in lanes
        ]
        lane_heights_frac = [0.068 if lc == 1 else 0.104 for lc in line_counts]
        trace_gap_frac = 0.055
        outer_frac = 0.032
        reserved_frac = trace_gap_frac + sum(lane_heights_frac) + outer_frac

        # FTIR transmittance troughs point downward, so put annotations below
        # the trace unless the user selected the top trace of a stacked plot.
        labels_below = (
            mode == "ftir"
            and y_unit == "transmittance"
            and not (stacked_multi and active_trace == len(traces) - 1)
        )

        if labels_below:
            upper_y = shifted_max + pad * 0.45 + legend_extra
            data_and_top = max(upper_y - y_min, y_span)
            total_span = data_and_top / max(1.0 - reserved_frac, 0.35)
            lower_y = upper_y - total_span
            ax.set_ylim(lower_y, upper_y)
            y_cur = y_min - trace_gap_frac * total_span
            lane_sign = -1.0
            va = "top"
        else:
            lower_y = y_min - pad * 0.45
            data_and_bottom = max(shifted_max - lower_y, y_span)
            legend_frac = (0.09 + 0.025 * min(len(traces), 4)) if len(traces) > 1 else 0.0
            total_span = data_and_bottom / max(1.0 - reserved_frac - legend_frac, 0.35)
            upper_y = lower_y + total_span
            ax.set_ylim(lower_y, upper_y)
            y_cur = shifted_max + trace_gap_frac * total_span
            lane_sign = 1.0
            va = "bottom"

        lane_ys = []
        for lane in range(lane_count):
            lane_ys.append(y_cur)
            y_cur += lane_sign * lane_heights_frac[lane] * total_span

        x_span = max(x_hi - x_lo, 1e-9)
        lane_indices = [idx % lane_count for idx in range(n_vis)]
        if mode == "uvvis":
            targets_u = [(item[0] - x_lo) / x_span for item in wrapped_items]
        else:
            targets_u = [(x_hi - item[0]) / x_span for item in wrapped_items]
        half_widths = [
            max(0.042, max(len(ln) for ln in item[2].splitlines()) * 0.0043 + 0.011)
            for item in wrapped_items
        ]
        placed_u = _space_peaks_u(targets_u, lane_indices, half_widths)

        for (px, py, wrapped, pinned), lane, u_text in zip(wrapped_items, lane_indices, placed_u):
            if mode == "uvvis":
                x_text = x_lo + u_text * x_span
            else:
                x_text = x_hi - u_text * x_span
            y_text = lane_ys[lane]
            color = cfg["palette"][active_trace % len(cfg["palette"])]
            if pinned:
                color = "#b8860b" if cfg["figure_facecolor"] == "#ffffff" else "#dbbc7f"

            # Pass 1: all leader lines first, at a low zorder, so they can never
            # slice through the opaque text pills drawn in pass 2.
            ax.plot(
                [px, x_text],
                [py, y_text],
                color=color,
                linewidth=0.60,
                alpha=0.88,
                zorder=14,
                clip_on=True,
            )

            # Pass 2: every text label on top with an opaque background pill.
            ax.text(
                x_text,
                y_text,
                wrapped,
                ha="center",
                va=va,
                fontsize=6.2,
                color=color,
                fontweight="bold" if pinned else "normal",
                bbox={
                    "boxstyle": "round,pad=0.18",
                    "facecolor": cfg["axes_facecolor"],
                    "edgecolor": "none",
                    "alpha": 0.94,
                },
                clip_on=True,
                zorder=22,
            )

    # UV-Vis: subtle 380-750 nm visible-light rainbow strip flush along the
    # bottom axis. Rendered AFTER the final set_ylim so it never drifts away
    # from the axis when the annotation lanes expand the plot.
    if mode == "uvvis":
        cur_lo, cur_hi = ax.get_ylim()
        y_bottom = cur_lo
        y_top = cur_lo + (cur_hi - cur_lo) * 0.032
        clip_lo = max(x_lo, 380.0)
        clip_hi = min(x_hi, 750.0)
        if clip_lo < clip_hi:
            wavelengths = np.linspace(clip_lo, clip_hi, 400)
            rgb_colors = np.array([wavelength_to_rgb(float(w)) for w in wavelengths])
            ax.imshow(
                rgb_colors[np.newaxis, :, :],
                aspect="auto",
                extent=[clip_lo, clip_hi, y_bottom, y_top],
                zorder=1,
                alpha=0.32,
                origin="lower",
            )

    if mode == "ftir" and y_unit == "transmittance" and not stacked_multi:
        if y_min >= -2.0 and shifted_max <= 105.0:
            ticks = [t for t in ax.get_yticks() if -1e-6 <= t <= 100.0 + 1e-6]
            if ticks:
                ax.set_yticks(ticks)

    if len(traces) > 1:
        legend = ax.legend(
            loc="upper right",
            fontsize=7.2,
            handlelength=1.8,
            borderaxespad=0.4,
            labelcolor=cfg["text"],
        )
        legend.set_zorder(30)

    fig.tight_layout(pad=0.6)
    return fig


def main(argv=None):
    args = parse_args(argv)
    try:
        payload = load_payload(args.input)
        fig = render(payload)

        out_dir = os.path.dirname(os.path.abspath(args.output))
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)

        dpi = int(payload.get("dpi", 600) or 600)
        fig.savefig(
            args.output,
            dpi=dpi,
            bbox_inches="tight",
            facecolor=fig.get_facecolor(),
            edgecolor="none",
        )
        plt.close(fig)
    except Exception as exc:  # pragma: no cover - surfaced to the Rust CLI
        print(f"SPECTRA export error: {exc}", file=sys.stderr)
        return 1

    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
