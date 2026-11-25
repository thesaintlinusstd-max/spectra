//! Pure-Rust publication-quality SVG exporter — zero extra dependencies.
//!
//! Matches the Python `export_plot.py` behaviour:
//!   * IUPAC-reversed X axis for FTIR/NMR
//!   * ACS / RSC / Everforest journal styles
//!   * Stacked or overlaid multi-trace plots
//!   * Two-pass non-crossing peak annotations with opaque pills
//!   * UV-Vis visible-light rainbow gradient (Dan Bruton CIE mapping)

use std::collections::HashMap;
use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::export::{ExportPayload, ExportStyle, ExportTrace};

const CANVAS_W: f64 = 960.0;
const PLOT_L: f64 = 78.0;
const PLOT_R: f64 = 932.0;
const PLOT_T: f64 = 54.0;
const PLOT_B: f64 = 506.0;
const TITLE_Y: f64 = 28.0;
const LEGEND_X_ANCHOR: f64 = PLOT_R - 12.0;

struct StyleCfg {
    bg: &'static str,
    axes_bg: &'static str,
    text: &'static str,
    spine: &'static str,
    grid: &'static str,
    palette: [&'static str; 6],
    all_spines: bool,
    grid_alpha: f32,
    line_width: f32,
}

fn style_cfg(style: ExportStyle) -> StyleCfg {
    match style {
        ExportStyle::Acs => StyleCfg {
            bg: "#ffffff",
            axes_bg: "#ffffff",
            text: "#111111",
            spine: "#111111",
            grid: "#c8c8c8",
            palette: ["#1f4e79", "#c00000", "#2e7d32", "#6a1b9a", "#e65100", "#00838f"],
            all_spines: false,
            grid_alpha: 0.35,
            line_width: 1.6,
        },
        ExportStyle::Rsc => StyleCfg {
            bg: "#ffffff",
            axes_bg: "#ffffff",
            text: "#101820",
            spine: "#101820",
            grid: "#d5d9de",
            palette: ["#003f5c", "#bc5090", "#ffa600", "#58508d", "#ff6361", "#2f4b7c"],
            all_spines: false,
            grid_alpha: 0.30,
            line_width: 1.65,
        },
        ExportStyle::Everforest => StyleCfg {
            bg: "#2d353b",
            axes_bg: "#2d353b",
            text: "#d3c6aa",
            spine: "#859289",
            grid: "#475258",
            palette: ["#a7c080", "#dbbc7f", "#7fbbb3", "#d699b6", "#e67e80", "#83c092"],
            all_spines: true,
            grid_alpha: 0.35,
            line_width: 1.7,
        },
    }
}

fn xml_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 16);
    for ch in s.chars() {
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&apos;"),
            c => out.push(c),
        }
    }
    out
}

fn nice_step(span: f64, target_ticks: f64) -> f64 {
    let rough = (span / target_ticks).max(1e-12);
    let mag = 10.0_f64.powf(rough.log10().floor());
    for c in [1.0, 2.0, 5.0, 10.0] {
        let step = c * mag;
        if step >= rough {
            return step;
        }
    }
    10.0 * mag
}

fn build_ticks(lo: f64, hi: f64, step: f64) -> Vec<f64> {
    let mut ticks = Vec::new();
    let start = (lo / step).ceil() * step;
    let mut v = start;
    let mut guard = 0;
    while v <= hi + 1e-9 && guard < 256 {
        ticks.push(v);
        v += step;
        guard += 1;
    }
    ticks
}

fn format_tick(v: f64) -> String {
    if v.abs() < 1e-6 {
        "0".to_string()
    } else if v.fract().abs() < 1e-6 {
        format!("{:.0}", v)
    } else {
        format!("{:.1}", v)
    }
}

fn x_label_markup(mode: &str) -> &'static str {
    match mode {
        "hnmr" | "cnmr" => "Chemical shift δ (ppm)",
        "uvvis" => "Wavelength, λ (nm)",
        _ => r#"Wavenumber (cm<tspan baseline-shift="super" font-size="75%">−1</tspan>)"#,
    }
}

fn y_label_text(y_unit: &str, stacked: bool, traces: usize) -> String {
    let base = match y_unit {
        "transmittance" => "Transmittance (%)",
        "absorbance" => "Absorbance (a.u.)",
        _ => "Intensity (a.u.)",
    };
    if stacked && traces > 1 {
        format!("{} (stacked)", base)
    } else {
        base.to_string()
    }
}

fn decimate_trace(pts: &[(f64, f64)], fx: &dyn Fn(f64) -> f64) -> Vec<(f64, f64)> {
    if pts.len() <= 3000 {
        return pts.to_vec();
    }
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < pts.len() {
        let col = fx(pts[i].0).round() as i64;
        let start = i;
        while i < pts.len() && fx(pts[i].0).round() as i64 == col {
            i += 1;
        }
        let slice = &pts[start..i];
        if slice.len() <= 2 {
            out.extend_from_slice(slice);
            continue;
        }
        let mut indices = vec![0usize, slice.len() - 1];
        let mut min_idx = 0;
        let mut max_idx = 0;
        for (j, p) in slice.iter().enumerate() {
            if p.1 < slice[min_idx].1 {
                min_idx = j;
            }
            if p.1 > slice[max_idx].1 {
                max_idx = j;
            }
        }
        indices.push(min_idx);
        indices.push(max_idx);
        indices.sort_unstable();
        indices.dedup();
        for &idx in &indices {
            out.push(slice[idx]);
        }
    }
    out
}

pub fn write_publication_svg(payload: &ExportPayload, output_path: &Path) -> Result<PathBuf> {
    if let Some(parent) = output_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }

    let style = match payload.style.as_str() {
        "rsc" => ExportStyle::Rsc,
        "everforest" => ExportStyle::Everforest,
        _ => ExportStyle::Acs,
    };
    let cfg = style_cfg(style);

    // ------------------------------------------------------ data preparation
    let mut traces: Vec<(&ExportTrace, Vec<(f64, f64)>)> = Vec::new();
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut y_lo = f64::INFINITY;
    let mut y_hi = f64::NEG_INFINITY;

    for t in &payload.traces {
        let mut pts: Vec<(f64, f64)> = t
            .x
            .iter()
            .zip(&t.y)
            .filter_map(|(&x, &y)| {
                if x.is_finite() && y.is_finite() {
                    Some((x, y))
                } else {
                    None
                }
            })
            .collect();
        if pts.len() < 2 {
            continue;
        }
        pts.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (a, b) = (pts[0].0, pts.last().unwrap().0);
        x_lo = x_lo.min(a);
        x_hi = x_hi.max(b);
        for &(_, y) in &pts {
            y_lo = y_lo.min(y);
            y_hi = y_hi.max(y);
        }
        traces.push((t, pts));
    }
    if traces.is_empty() {
        bail!("No valid traces to export");
    }

    let bounds = payload.x_bounds;
    if bounds[0].is_finite() && bounds[1].is_finite() && (bounds[0] - bounds[1]).abs() > 1e-9 {
        x_lo = bounds[0].min(bounds[1]);
        x_hi = bounds[0].max(bounds[1]);
    }

    let is_uvvis = payload.mode.eq_ignore_ascii_case("uvvis");
    let x_reversed = !is_uvvis;

    // Stacking offsets
    let y_span = (y_hi - y_lo).max(1e-9);
    let offsets: Vec<f64> = if payload.stacked && traces.len() > 1 {
        let step = y_span * 1.12;
        (0..traces.len()).map(|i| i as f64 * step).collect()
    } else {
        vec![0.0; traces.len()]
    };

    let shifted_max = y_hi + *offsets.last().unwrap_or(&0.0);
    let data_span = (shifted_max - y_lo).max(y_span);
    let pad = (data_span * 0.06).max(y_span * 0.04);

    let active_trace = payload.active_trace.min(traces.len().saturating_sub(1));
    let stacked_multi = payload.stacked && traces.len() > 1;

    // Peak candidates (pinned first, max 10)
    struct PeakCand {
        px: f64,
        py: f64,
        label: String,
        pinned: bool,
    }
    let mut peak_cands: Vec<PeakCand> = Vec::new();
    for pk in &payload.peaks {
        let px = pk.x;
        let py = pk.y + offsets[active_trace];
        if !px.is_finite() || !py.is_finite() {
            continue;
        }
        if px < x_lo - 1e-9 || px > x_hi + 1e-9 {
            continue;
        }
        let label = pk.label.trim();
        if label.is_empty() {
            continue;
        }
        peak_cands.push(PeakCand {
            px,
            py,
            label: label.to_string(),
            pinned: pk.pinned,
        });
    }
    // Pinned first
    peak_cands.sort_by_key(|a| std::cmp::Reverse(a.pinned));
    let mut visible_peaks: Vec<PeakCand> = peak_cands.into_iter().take(10).collect();
    // Sort left-to-right in screen coordinates
    visible_peaks.sort_by(|a, b| {
        if is_uvvis {
            a.px.total_cmp(&b.px)
        } else {
            b.px.total_cmp(&a.px)
        }
    });

    let n_vis = visible_peaks.len();
    let lane_count = if n_vis <= 2 {
        1
    } else if n_vis <= 4 {
        2
    } else {
        3
    };

    struct WrappedPeak {
        px: f64,
        py: f64,
        lines: Vec<String>,
        pinned: bool,
    }
    let mut wrapped: Vec<WrappedPeak> = Vec::with_capacity(n_vis);
    for p in &visible_peaks {
        let mut lines: Vec<String> = Vec::new();
        let mut cur = String::new();
        for w in p.label.split_whitespace() {
            if cur.is_empty() {
                cur.push_str(w);
            } else if cur.len() + w.len() < 19 {
                cur.push(' ');
                cur.push_str(w);
            } else {
                lines.push(cur);
                cur = w.to_string();
            }
        }
        if !cur.is_empty() {
            lines.push(cur);
        }
        if lines.len() > 2 {
            let mut kept = lines[..2].to_vec();
            let mut last = kept.pop().unwrap();
            if last.len() >= 19 {
                last = last.chars().take(18).collect();
            }
            last.push('…');
            kept.push(last);
            lines = kept;
        }
        if lines.is_empty() {
            lines.push(p.label.chars().take(19).collect());
        }
        wrapped.push(WrappedPeak {
            px: p.px,
            py: p.py,
            lines,
            pinned: p.pinned,
        });
    }

    let mut lane_heights_frac = Vec::with_capacity(lane_count);
    for lane in 0..lane_count {
        let mut max_lines = 1;
        for (i, wp) in wrapped.iter().enumerate() {
            if i % lane_count == lane {
                max_lines = max_lines.max(wp.lines.len());
            }
        }
        lane_heights_frac.push(if max_lines == 1 { 0.068 } else { 0.104 });
    }
    let trace_gap_frac = 0.055;
    let outer_frac = 0.032;
    let legend_frac = if traces.len() > 1 {
        0.09 + 0.025 * (traces.len().min(4) as f64)
    } else {
        0.0
    };
    let reserved_frac = if n_vis > 0 {
        trace_gap_frac + lane_heights_frac.iter().sum::<f64>() + outer_frac
    } else {
        0.0
    };

    let labels_below = payload.mode.eq_ignore_ascii_case("ftir")
        && payload.y_unit.eq_ignore_ascii_case("transmittance")
        && !(stacked_multi && active_trace == traces.len() - 1);

    let (y_bottom_data, y_top_data) = if n_vis == 0 {
        let legend_extra = data_span * legend_frac;
        (y_lo - pad, shifted_max + pad + legend_extra)
    } else if labels_below {
        let upper_y = shifted_max + pad * 0.45 + data_span * legend_frac;
        let data_and_top = (upper_y - y_lo).max(y_span);
        let total_span = data_and_top / (1.0 - reserved_frac).max(0.35);
        let lower_y = upper_y - total_span;
        (lower_y, upper_y)
    } else {
        let lower_y = y_lo - pad * 0.45;
        let data_and_bottom = (shifted_max - lower_y).max(y_span);
        let total_span = data_and_bottom / (1.0 - reserved_frac - legend_frac).max(0.35);
        let upper_y = lower_y + total_span;
        (lower_y, upper_y)
    };

    let data_h = (y_top_data - y_bottom_data).max(1e-9);

    let mut lane_ys = Vec::with_capacity(lane_count);
    if labels_below {
        let mut y_cur = y_lo - trace_gap_frac * data_h;
        for h_frac in &lane_heights_frac {
            lane_ys.push(y_cur);
            y_cur -= h_frac * data_h;
        }
    } else {
        let mut y_cur = shifted_max + trace_gap_frac * data_h;
        for h_frac in &lane_heights_frac {
            lane_ys.push(y_cur);
            y_cur += h_frac * data_h;
        }
    }

    let plot_w = PLOT_R - PLOT_L;
    let plot_h = PLOT_B - PLOT_T;
    let x_data_span = (x_hi - x_lo).max(1e-9);

    let fx = |x: f64| -> f64 {
        let u = (x - x_lo) / x_data_span;
        if x_reversed {
            PLOT_R - u * plot_w
        } else {
            PLOT_L + u * plot_w
        }
    };
    let fy = |y: f64| -> f64 {
        let u = (y - y_bottom_data) / data_h;
        PLOT_B - u * plot_h
    };

    let x_step = nice_step(x_data_span, 7.0);
    let x_ticks = build_ticks(x_lo, x_hi, x_step);
    let x_minor_ticks = build_ticks(x_lo, x_hi, x_step / 5.0);

    let show_y_ticks = !stacked_multi;
    let y_ticks = if show_y_ticks {
        let y_step = nice_step(data_h, 6.0);
        let mut ticks = build_ticks(y_bottom_data, y_top_data, y_step);
        if payload.y_unit.eq_ignore_ascii_case("transmittance") {
            ticks.retain(|&v| (-1e-6..=100.0 + 1e-6).contains(&v));
        }
        ticks
    } else {
        Vec::new()
    };

    let rainbow = if is_uvvis {
        let clip_lo = 380.0_f64.max(x_lo);
        let clip_hi = 750.0_f64.min(x_hi);
        if clip_lo < clip_hi {
            let n_stops = 24;
            let mut grad = String::new();
            for i in 0..n_stops {
                let u = i as f64 / (n_stops - 1) as f64;
                let nm = clip_lo + u * (clip_hi - clip_lo);
                let (r, g, b) = crate::kitty_plot::wavelength_to_rgb(nm);
                let _ = write!(
                    grad,
                    r#"<stop offset="{:.1}%" stop-color="rgb({},{},{})"/>"#,
                    u * 100.0,
                    r,
                    g,
                    b
                );
            }
            Some((clip_lo, clip_hi, grad))
        } else {
            None
        }
    } else {
        None
    };

    let lane_indices: Vec<usize> = (0..n_vis).map(|i| i % lane_count).collect();
    let mut targets_u: Vec<f64> = Vec::with_capacity(n_vis);
    let mut half_widths: Vec<f64> = Vec::with_capacity(n_vis);
    for wp in &wrapped {
        let u = (fx(wp.px) - PLOT_L) / plot_w;
        targets_u.push(u.clamp(0.0, 1.0));
        let max_len = wp.lines.iter().map(|s| s.len()).max().unwrap_or(1) as f64;
        half_widths.push((max_len * 0.0043 + 0.011).max(0.042));
    }
    let placed_u = space_peaks_u(&targets_u, &lane_indices, &half_widths);

    // ------------------------------------------------------------ SVG output
    let mut svg = String::with_capacity(16 * 1024);
    let _ = writeln!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="960" height="580" viewBox="0 0 960 580">"#
    );
    let _ = writeln!(svg, r#"<rect width="960" height="580" fill="{}"/>"#, cfg.bg);
    let _ = writeln!(svg, r#"<defs>"#);
    let _ = writeln!(
        svg,
        r#"<clipPath id="plot-clip"><rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}"/></clipPath>"#,
        PLOT_L, PLOT_T, plot_w, plot_h
    );
    if let Some((_, _, stops)) = &rainbow {
        let _ = writeln!(
            svg,
            r#"<linearGradient id="vis-rainbow" x1="0%" y1="0%" x2="100%" y2="0%">{}</linearGradient>"#,
            stops
        );
    }
    let _ = writeln!(svg, r#"</defs>"#);

    if !payload.title.trim().is_empty() {
        let _ = writeln!(
            svg,
            r#"<text x="{}" y="{}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="15" font-weight="bold" fill="{}">{}</text>"#,
            CANVAS_W / 2.0,
            TITLE_Y,
            cfg.text,
            xml_escape(payload.title.trim())
        );
    }

    // Grid
    let _ = writeln!(
        svg,
        r#"<g id="grid" stroke="{}" stroke-width="0.6" opacity="{:.2}">"#,
        cfg.grid, cfg.grid_alpha
    );
    if show_y_ticks {
        for &yv in &y_ticks {
            let y_svg = fy(yv);
            let _ = writeln!(
                svg,
                r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}"/>"#,
                PLOT_L, y_svg, PLOT_R, y_svg
            );
        }
    }
    let _ = writeln!(svg, r#"</g>"#);

    // Axes spines and ticks
    let _ = writeln!(
        svg,
        r#"<g id="axes" stroke="{}" stroke-width="1.0" fill="none">"#,
        cfg.spine
    );
    let _ = writeln!(
        svg,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}"/>"#,
        PLOT_L, PLOT_B, PLOT_R, PLOT_B
    );
    let _ = writeln!(
        svg,
        r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}"/>"#,
        PLOT_L, PLOT_T, PLOT_L, PLOT_B
    );
    if cfg.all_spines {
        let _ = writeln!(
            svg,
            r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}"/>"#,
            PLOT_L, PLOT_T, PLOT_R, PLOT_T
        );
        let _ = writeln!(
            svg,
            r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}"/>"#,
            PLOT_R, PLOT_T, PLOT_R, PLOT_B
        );
    }
    for &xv in &x_ticks {
        let x_svg = fx(xv);
        let _ = writeln!(
            svg,
            r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}"/>"#,
            x_svg,
            PLOT_B,
            x_svg,
            PLOT_B + 5.0
        );
    }
    for &xv in &x_minor_ticks {
        if x_ticks.iter().any(|&t| (t - xv).abs() < 1e-9) {
            continue;
        }
        let x_svg = fx(xv);
        let _ = writeln!(
            svg,
            r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke-width="0.6"/>"#,
            x_svg,
            PLOT_B,
            x_svg,
            PLOT_B + 2.5
        );
    }
    if show_y_ticks {
        for &yv in &y_ticks {
            let y_svg = fy(yv);
            let _ = writeln!(
                svg,
                r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}"/>"#,
                PLOT_L - 5.0,
                y_svg,
                PLOT_L,
                y_svg
            );
        }
    }
    let _ = writeln!(svg, r#"</g>"#);

    // X tick labels
    let _ = writeln!(
        svg,
        r#"<g id="x-tick-labels" font-family="Arial,Helvetica,sans-serif" font-size="10.5" fill="{}">"#,
        cfg.text
    );
    for &xv in &x_ticks {
        let x_svg = fx(xv);
        let label = format_tick(xv);
        let _ = writeln!(
            svg,
            r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" dominant-baseline="hanging">{}</text>"#,
            x_svg,
            PLOT_B + 9.0,
            xml_escape(&label)
        );
    }
    let _ = writeln!(svg, r#"</g>"#);

    // Y tick labels
    if show_y_ticks {
        let _ = writeln!(
            svg,
            r#"<g id="y-tick-labels" font-family="Arial,Helvetica,sans-serif" font-size="10.5" fill="{}">"#,
            cfg.text
        );
        for &yv in &y_ticks {
            let y_svg = fy(yv);
            let label = format_tick(yv);
            let _ = writeln!(
                svg,
                r#"<text x="{:.1}" y="{:.1}" text-anchor="end" dominant-baseline="middle">{}</text>"#,
                PLOT_L - 9.0,
                y_svg,
                xml_escape(&label)
            );
        }
        let _ = writeln!(svg, r#"</g>"#);
    }

    // Axis titles
    let _ = writeln!(
        svg,
        r#"<text x="{}" y="{}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="12" fill="{}">{}</text>"#,
        (PLOT_L + PLOT_R) / 2.0,
        PLOT_B + 46.0,
        cfg.text,
        x_label_markup(&payload.mode)
    );

    let y_label = y_label_text(&payload.y_unit, payload.stacked, traces.len());
    let _ = writeln!(
        svg,
        r#"<text x="22" y="{}" text-anchor="middle" font-family="Arial,Helvetica,sans-serif" font-size="12" fill="{}" transform="rotate(-90, 22, {})">{}</text>"#,
        (PLOT_T + PLOT_B) / 2.0,
        cfg.text,
        (PLOT_T + PLOT_B) / 2.0,
        xml_escape(&y_label)
    );

    // UV-Vis rainbow strip
    if let Some((clip_lo, clip_hi, _)) = rainbow {
        let x1 = fx(clip_lo);
        let x2 = fx(clip_hi);
        let _ = writeln!(
            svg,
            r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="10" fill="url(#vis-rainbow)" fill-opacity="0.32"/>"#,
            x1.min(x2),
            PLOT_B - 10.0,
            (x2 - x1).abs()
        );
    }

    // Zero baseline for NMR or absorbance
    if !payload.mode.eq_ignore_ascii_case("ftir") || payload.y_unit.eq_ignore_ascii_case("absorbance") {
        let y0 = fy(0.0);
        if y0 > PLOT_T && y0 < PLOT_B {
            let _ = writeln!(
                svg,
                r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="0.8" opacity="0.65"/>"#,
                PLOT_L, y0, PLOT_R, y0, cfg.grid
            );
        }
    }

    // Traces
    let _ = writeln!(svg, r#"<g id="traces" clip-path="url(#plot-clip)">"#);
    for (idx, (_, pts)) in traces.iter().enumerate() {
        let decimated = decimate_trace(pts, &fx);
        let color = cfg.palette[idx % cfg.palette.len()];
        let is_active = idx == active_trace;
        let lw = cfg.line_width * if is_active || traces.len() == 1 { 1.18 } else { 0.92 };
        let opacity = if is_active || traces.len() == 1 { 1.0 } else { 0.88 };
        let mut points_str = String::with_capacity(decimated.len() * 14);
        for (x, y) in &decimated {
            let x_svg = fx(*x);
            let y_svg = fy(*y + offsets[idx]);
            let _ = write!(points_str, "{:.2},{:.2} ", x_svg, y_svg);
        }
        let _ = writeln!(
            svg,
            r#"<polyline fill="none" stroke="{}" stroke-width="{:.2}" stroke-opacity="{:.2}" stroke-linecap="round" stroke-linejoin="round" points="{}"/>"#,
            color,
            lw,
            opacity,
            points_str.trim_end()
        );
    }
    let _ = writeln!(svg, r#"</g>"#);

    // Deconvolved sub-peaks: dashed component curves + crisp cumulative fit.
    if let Some(dc) = &payload.deconv {
        let y_off = offsets.get(active_trace).copied().unwrap_or(0.0);
        let _ = writeln!(
            svg,
            r#"<g id="deconv-components" clip-path="url(#plot-clip)">"#
        );
        for (ci, comp) in dc.components.iter().enumerate() {
            let color = cfg.palette[(ci + 1) % cfg.palette.len()];
            let mut points_str = String::with_capacity(comp.len() * 14);
            for (i, &y) in comp.iter().enumerate() {
                if i >= dc.x.len() {
                    break;
                }
                let x = dc.x[i];
                if !x.is_finite() || !y.is_finite() {
                    continue;
                }
                let _ = write!(points_str, "{:.2},{:.2} ", fx(x), fy(y + y_off));
            }
            let _ = writeln!(
                svg,
                r#"<polyline fill="none" stroke="{}" stroke-width="1.15" stroke-opacity="0.78" stroke-dasharray="5,3" stroke-linecap="round" points="{}"/>"#,
                color,
                points_str.trim_end()
            );
        }
        let _ = writeln!(svg, r#"</g>"#);

        let _ = writeln!(svg, r#"<g id="deconv-fit" clip-path="url(#plot-clip)">"#);
        let mut fit_str = String::with_capacity(dc.y_fit.len() * 14);
        for (i, &y) in dc.y_fit.iter().enumerate() {
            if i >= dc.x.len() {
                break;
            }
            let x = dc.x[i];
            if !x.is_finite() || !y.is_finite() {
                continue;
            }
            let _ = write!(fit_str, "{:.2},{:.2} ", fx(x), fy(y + y_off));
        }
        let _ = writeln!(
            svg,
            r##"<polyline fill="none" stroke="#d49b1a" stroke-width="1.35" stroke-opacity="0.92" stroke-linecap="round" stroke-linejoin="round" points="{}"/>"##,
            fit_str.trim_end()
        );
        let _ = writeln!(svg, r#"</g>"#);
    }

    // Legend
    if traces.len() > 1 {
        let _ = writeln!(
            svg,
            r#"<g id="legend" font-family="Arial,Helvetica,sans-serif" font-size="10.5" fill="{}">"#,
            cfg.text
        );
        let mut y_pos = PLOT_T + 10.0;
        for (idx, (t, _)) in traces.iter().enumerate() {
            let color = cfg.palette[idx % cfg.palette.len()];
            let label = xml_escape(&t.name);
            let text_w = t.name.chars().count() as f64 * 6.0;
            let swatch_x2 = LEGEND_X_ANCHOR - text_w - 8.0;
            let swatch_x1 = swatch_x2 - 22.0;
            let swatch_y = y_pos + 2.0;
            let _ = writeln!(
                svg,
                r#"<line x1="{:.1}" y1="{:.1}" x2="{:.1}" y2="{:.1}" stroke="{}" stroke-width="{:.2}" stroke-linecap="round"/>"#,
                swatch_x1, swatch_y, swatch_x2, swatch_y, color, cfg.line_width
            );
            let _ = writeln!(
                svg,
                r#"<text x="{:.1}" y="{:.1}" text-anchor="end" dominant-baseline="middle">{}</text>"#,
                LEGEND_X_ANCHOR,
                y_pos + 2.5,
                label
            );
            y_pos += 16.0;
        }
        let _ = writeln!(svg, r#"</g>"#);
    }

    // ------------------------------------------------------ peak annotations
    // Pass 1: leader lines
    let trace_color = cfg.palette[active_trace % cfg.palette.len()];
    let pinned_color = if cfg.bg == "#ffffff" { "#b8860b" } else { "#dbbc7f" };

    let _ = writeln!(svg, r#"<g id="peak-leaders" fill="none">"#);
    for ((wp, &lane), &u_text) in wrapped.iter().zip(&lane_indices).zip(&placed_u) {
        let x_text = PLOT_L + u_text * plot_w;
        let y_anchor = fy(lane_ys[lane]);
        let color = if wp.pinned { pinned_color } else { trace_color };
        let x_peak = fx(wp.px);
        let y_peak = fy(wp.py);
        let _ = writeln!(
            svg,
            r#"<line x1="{:.2}" y1="{:.2}" x2="{:.2}" y2="{:.2}" stroke="{}" stroke-width="0.9"/>"#,
            x_peak, y_peak, x_text, y_anchor, color
        );
    }
    let _ = writeln!(svg, r#"</g>"#);

    // Pass 2: opaque rounded pills + centered text
    let _ = writeln!(
        svg,
        r#"<g id="peak-labels" font-family="Arial,Helvetica,sans-serif" font-size="9.5">"#
    );
    const LINE_H: f64 = 11.0;
    const PAD_X: f64 = 4.5;
    const PAD_Y: f64 = 2.5;
    for ((wp, &lane), &u_text) in wrapped.iter().zip(&lane_indices).zip(&placed_u) {
        let x_text = PLOT_L + u_text * plot_w;
        let y_anchor = fy(lane_ys[lane]);
        let font_weight = if wp.pinned { "bold" } else { "normal" };
        let label_color = if wp.pinned { pinned_color } else { trace_color };

        let max_len = wp.lines.iter().map(|s| s.chars().count()).max().unwrap_or(1) as f64;
        let pill_w = max_len * 5.2 + 2.0 * PAD_X;
        let pill_h = wp.lines.len() as f64 * LINE_H + 2.0 * PAD_Y;
        let pill_x = x_text - pill_w / 2.0;
        let pill_y = if labels_below {
            y_anchor
        } else {
            y_anchor - pill_h
        };

        let _ = writeln!(
            svg,
            r#"<rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}" rx="3" ry="3" fill="{}" fill-opacity="0.94"/>"#,
            pill_x, pill_y, pill_w, pill_h, cfg.axes_bg
        );

        for (line_idx, line) in wp.lines.iter().enumerate() {
            let baseline = pill_y + PAD_Y + 8.5 + line_idx as f64 * LINE_H;
            let _ = writeln!(
                svg,
                r#"<text x="{:.1}" y="{:.1}" text-anchor="middle" fill="{}" font-weight="{}">{}</text>"#,
                x_text,
                baseline,
                label_color,
                font_weight,
                xml_escape(line)
            );
        }
    }
    let _ = writeln!(svg, r#"</g>"#);

    svg.push_str("</svg>");

    fs::write(output_path, svg)
        .with_context(|| format!("Failed to write SVG to {}", output_path.display()))?;
    Ok(output_path.to_path_buf())
}

fn space_peaks_u(targets: &[f64], lanes: &[usize], half_widths: &[f64]) -> Vec<f64> {
    let n = targets.len();
    if n == 0 {
        return Vec::new();
    }
    let mut lo = Vec::with_capacity(n);
    let mut hi = Vec::with_capacity(n);
    for hw in half_widths {
        lo.push(*hw + 0.014);
        hi.push(1.0 - *hw - 0.014);
    }
    if n == 1 {
        return vec![targets[0].clamp(lo[0], hi[0])];
    }

    let mut prev_same = vec![-1_i32; n];
    let mut last_in_lane: HashMap<usize, usize> = HashMap::new();
    for (i, &lane) in lanes.iter().enumerate() {
        if let Some(&last) = last_in_lane.get(&lane) {
            prev_same[i] = last as i32;
        }
        last_in_lane.insert(lane, i);
    }

    let mut next_same = vec![-1_i32; n];
    for (i, &p) in prev_same.iter().enumerate() {
        if p >= 0 {
            next_same[p as usize] = i as i32;
        }
    }

    const MIN_STEP: f64 = 0.018;
    const LANE_PAD: f64 = 0.016;

    let mut fwd = Vec::with_capacity(n);
    for i in 0..n {
        let mut v = targets[i].clamp(lo[i], hi[i]);
        if i > 0 {
            v = v.max(fwd[i - 1] + MIN_STEP);
        }
        let p = prev_same[i];
        if p >= 0 {
            v = v.max(fwd[p as usize] + half_widths[p as usize] + half_widths[i] + LANE_PAD);
        }
        fwd.push(v);
    }
    if fwd[n - 1] > hi[n - 1] {
        fwd[n - 1] = hi[n - 1];
        for i in (0..n - 1).rev() {
            fwd[i] = fwd[i].min(fwd[i + 1] - MIN_STEP);
            let nx = next_same[i];
            if nx >= 0 {
                fwd[i] = fwd[i].min(fwd[nx as usize] - (half_widths[i] + half_widths[nx as usize] + LANE_PAD));
            }
        }
    }

    let mut bwd: Vec<f64> = (0..n).map(|i| targets[i].clamp(lo[i], hi[i])).collect();
    for i in (0..n - 1).rev() {
        bwd[i] = bwd[i].min(bwd[i + 1] - MIN_STEP);
        let nx = next_same[i];
        if nx >= 0 {
            bwd[i] = bwd[i].min(bwd[nx as usize] - (half_widths[i] + half_widths[nx as usize] + LANE_PAD));
        }
    }
    if bwd[0] < lo[0] {
        bwd[0] = lo[0];
        for i in 1..n {
            bwd[i] = bwd[i].max(bwd[i - 1] + MIN_STEP);
            let p = prev_same[i];
            if p >= 0 {
                bwd[i] = bwd[i].max(bwd[p as usize] + half_widths[p as usize] + half_widths[i] + LANE_PAD);
            }
        }
    }

    (0..n)
        .map(|i| (0.5 * (fwd[i] + bwd[i])).clamp(lo[i], hi[i]))
        .collect()
}
