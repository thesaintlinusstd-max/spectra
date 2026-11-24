//! Pure-Rust offline spectrum image digitizer.
//!
//! Decodes a raster plot (PNG/JPEG/WebP/BMP), finds the plot frame, removes
//! gridlines with a `max(row_median, col_median)` estimate, then traces the
//! curve with a Viterbi dynamic program. No Python, no scipy, no network.

use std::cmp::Ordering;
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::data::{Spectrum, SpectrumKind, YUnit};

/// Maximum candidate runs considered per image column.
const MAX_CANDS: usize = 6;

#[derive(Clone, Copy)]
struct Candidate {
    cy: f64,
    strength: f64,
}

fn percentile(values: &mut [f32], q: f64) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let q = q.clamp(0.0, 1.0);
    let idx = ((values.len() - 1) as f64 * q).round() as usize;
    let (_, &mut v, _) = values.select_nth_unstable_by(idx, |a, b| {
        a.partial_cmp(b).unwrap_or(Ordering::Equal)
    });
    v
}

fn median(values: &mut [f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let mid = values.len() / 2;
    let (_, &mut v, _) = values.select_nth_unstable_by(mid, |a, b| {
        a.partial_cmp(b).unwrap_or(Ordering::Equal)
    });
    v
}

fn median3(a: f64, b: f64, c: f64) -> f64 {
    a.max(b).min(a.min(b).max(c))
}

/// Digitize a spectrum plot image into a calibrated `Spectrum`.
///
/// `x_left` is the physical X value at the left edge of the plot frame and
/// `x_right` the value at the right edge (so FTIR calls pass `4000, 400`).
pub fn digitize_spectrum_image(
    img_path: &Path,
    kind: SpectrumKind,
    y_unit: YUnit,
    x_left: f64,
    x_right: f64,
    bg_mode: &str,
) -> Result<Spectrum> {
    if !(x_right - x_left).is_finite() || (x_right - x_left).abs() < 1e-12 {
        bail!("X Left and X Right must be finite and distinct");
    }

    let img = image::open(img_path)
        .with_context(|| format!("Failed to decode image {}", img_path.display()))?
        .to_rgba8();
    let (w, h) = (img.width() as usize, img.height() as usize);
    if h < 40 || w < 60 {
        bail!("Image is too small to digitize ({}x{})", w, h);
    }

    // ---------------------------------------------------------------- load
    let bg_fill = if bg_mode == "dark" { 0.0_f32 } else { 1.0_f32 };
    let mut gray = vec![0.0_f32; w * h];
    let mut sat = vec![0.0_f32; w * h];
    for (i, px) in img.pixels().enumerate() {
        let p = px.0;
        let a = p[3] as f32 / 255.0;
        let rf = (p[0] as f32 / 255.0) * a + bg_fill * (1.0 - a);
        let gf = (p[1] as f32 / 255.0) * a + bg_fill * (1.0 - a);
        let bf = (p[2] as f32 / 255.0) * a + bg_fill * (1.0 - a);
        gray[i] = 0.299 * rf + 0.587 * gf + 0.114 * bf;
        sat[i] = rf.max(gf).max(bf) - rf.min(gf).min(bf);
    }

    // -------------------------------------------------- background polarity
    let total = w * h;
    let stride = (total / 20_000).max(1);
    let mut sample: Vec<f32> = gray.iter().step_by(stride).copied().collect();
    let median_lum = percentile(&mut sample, 0.5);
    let dark_bg = bg_mode == "dark" || (bg_mode == "auto" && median_lum < 0.45);

    let mut fg = vec![0.0_f32; total];
    for i in 0..total {
        fg[i] = if dark_bg {
            (gray[i] - median_lum).clamp(0.0, 1.0)
        } else {
            (median_lum - gray[i]).clamp(0.0, 1.0)
        };
    }

    // ---------------------------------------------------- plot frame search
    // Prefer dark/solid axis spines (fg > 0.30) so light-gray gridlines
    // (e.g. RGB 210,210,210 -> fg ~0.18) are not mistaken for the frame.
    // Fall back to fg > 0.16 if the plot uses very faint axis lines.
    let detect_frame = |ink_thresh: f32| -> (usize, usize, usize, usize, bool) {
        let mut row_proj = vec![0.0_f32; h];
        let mut col_proj = vec![0.0_f32; w];
        for r in 0..h {
            let mut n = 0usize;
            for c in 0..w {
                if fg[r * w + c] > ink_thresh {
                    n += 1;
                }
            }
            row_proj[r] = n as f32 / w as f32;
        }
        for c in 0..w {
            let mut n = 0usize;
            for r in 0..h {
                if fg[r * w + c] > ink_thresh {
                    n += 1;
                }
            }
            col_proj[c] = n as f32 / h as f32;
        }

        let mut y0 = h * 8 / 100;
        let mut found_t = false;
        for (r, &p) in row_proj.iter().enumerate().take(h * 45 / 100) {
            if p > 0.32 {
                y0 = r + 3;
                found_t = true;
            }
        }

        let mut y1 = h * 86 / 100;
        let mut found_b = false;
        for (r, &p) in row_proj.iter().enumerate().skip(h * 55 / 100) {
            if p > 0.32 {
                y1 = r.saturating_sub(3);
                found_b = true;
                break;
            }
        }

        let mut x0 = w * 10 / 100;
        let mut found_l = false;
        for (c, &p) in col_proj.iter().enumerate().take(w * 40 / 100) {
            if p > 0.32 {
                x0 = c + 3;
                found_l = true;
            }
        }

        let mut x1 = w * 95 / 100;
        let mut found_r = false;
        for (c, &p) in col_proj.iter().enumerate().skip(w * 60 / 100) {
            if p > 0.32 {
                x1 = c.saturating_sub(3);
                found_r = true;
                break;
            }
        }

        (y0, y1, x0, x1, found_t && found_b && found_l && found_r)
    };

    let (mut y0, mut y1, mut x0, mut x1, found_all) = detect_frame(0.30);
    if !found_all {
        let (fy0, fy1, fx0, fx1, f_all) = detect_frame(0.16);
        if f_all {
            y0 = fy0;
            y1 = fy1;
            x0 = fx0;
            x1 = fx1;
        }
    }

    let bad_frame = y1 <= y0
        || (y1 - y0) < h * 25 / 100
        || x1 <= x0
        || (x1 - x0) < w * 30 / 100;
    if bad_frame {
        y0 = h * 8 / 100;
        y1 = h * 86 / 100;
        x0 = w * 10 / 100;
        x1 = w * 95 / 100;
    }
    let y1 = y1.min(h);
    let x1 = x1.min(w);
    if y1 <= y0 + 2 || x1 <= x0 + 2 {
        bail!("Could not locate a usable plot frame in image");
    }

    let ch = y1 - y0;
    let cw = x1 - x0;

    // ------------------------------------------------- curve vs gridlines
    let mut crop_fg = Vec::with_capacity(ch * cw);
    let mut crop_sat = Vec::with_capacity(ch * cw);
    for r in y0..y1 {
        let base = r * w + x0;
        crop_fg.extend_from_slice(&fg[base..base + cw]);
        crop_sat.extend_from_slice(&sat[base..base + cw]);
    }

    let mut sat_probe = crop_sat.clone();
    let sat_p995 = percentile(&mut sat_probe, 0.995);
    let has_color = sat_p995 > 0.20;

    let mut trace_map = vec![0.0_f32; ch * cw];
    if has_color {
        for i in 0..ch * cw {
            trace_map[i] = crop_sat[i] * 1.45 + crop_fg[i] * 0.45;
        }
    } else {
        // Exact per-row and per-column medians, then subtract the stronger of
        // the two. `max` (not `min`) is what kills both horizontal AND vertical
        // gridlines without erasing the curve.
        let mut row_med = vec![0.0_f32; ch];
        let mut col_med = vec![0.0_f32; cw];
        let mut scratch = Vec::with_capacity(cw.max(ch));
        for r in 0..ch {
            scratch.clear();
            scratch.extend_from_slice(&crop_fg[r * cw..(r + 1) * cw]);
            row_med[r] = median(&mut scratch);
        }
        for c in 0..cw {
            scratch.clear();
            for r in 0..ch {
                scratch.push(crop_fg[r * cw + c]);
            }
            col_med[c] = median(&mut scratch);
        }
        for (r, &rm) in row_med.iter().enumerate().take(ch) {
            for (c, &cm) in col_med.iter().enumerate().take(cw) {
                let i = r * cw + c;
                trace_map[i] = (crop_fg[i] - 0.85 * rm.max(cm)).clamp(0.0, 1.0);
            }
        }
    }

    let mut probe = trace_map.clone();
    let peak_ref = percentile(&mut probe, 0.988);
    let thresh = (peak_ref * 0.42).max(0.14);

    // ------------------------------------------------------ Viterbi tracer
    let mut valid_cols: Vec<usize> = Vec::new();
    let mut col_cands: Vec<Vec<Candidate>> = Vec::new();
    let max_run = ch as f32 * 0.40;

    #[allow(clippy::needless_range_loop)]
    for c in 0..cw {
        let mut cands: Vec<Candidate> = Vec::new();
        let mut r = 0usize;
        while r < ch {
            if trace_map[r * cw + c] >= thresh {
                let start = r;
                let mut sum = 0.0_f64;
                let mut wsum = 0.0_f64;
                let mut strength = 0.0_f32;
                while r < ch && trace_map[r * cw + c] >= thresh {
                    let v = trace_map[r * cw + c];
                    sum += r as f64 * v as f64;
                    wsum += v as f64;
                    strength = strength.max(v);
                    r += 1;
                }
                let height = (r - start) as f32;
                if height <= max_run && wsum > 0.0 {
                    cands.push(Candidate {
                        cy: sum / wsum,
                        strength: strength as f64,
                    });
                }
            } else {
                r += 1;
            }
        }
        if !cands.is_empty() {
            cands.sort_by(|a, b| b.strength.partial_cmp(&a.strength).unwrap_or(Ordering::Equal));
            cands.truncate(MAX_CANDS);
            valid_cols.push(c);
            col_cands.push(cands);
        }
    }

    if valid_cols.len() < 12 {
        bail!("Curve trace failed: only {} usable columns found", valid_cols.len());
    }

    let n = valid_cols.len();
    let jump_cap = 0.45 * ch as f64;
    let mut dp: Vec<Vec<f64>> = vec![Vec::new(); n];
    let mut back: Vec<Vec<usize>> = vec![Vec::new(); n];

    for k in 0..n {
        let len = col_cands[k].len();
        dp[k] = vec![0.0_f64; len];
        back[k] = vec![0usize; len];
        if k == 0 {
            for j in 0..len {
                dp[0][j] = -col_cands[0][j].strength * 24.0;
            }
            continue;
        }
        let gap = (valid_cols[k] - valid_cols[k - 1]).max(1) as f64;
        let prev_len = col_cands[k - 1].len();
        for j in 0..len {
            let unary = -col_cands[k][j].strength * 24.0;
            let mut best_prev = 0usize;
            let mut best_cost = f64::INFINITY;
            for pj in 0..prev_len {
                let jump = (col_cands[k][j].cy - col_cands[k - 1][pj].cy).abs();
                let trans = jump.min(jump_cap) * (1.15 / gap.sqrt());
                let cost = dp[k - 1][pj] + trans;
                if cost < best_cost {
                    best_cost = cost;
                    best_prev = pj;
                }
            }
            dp[k][j] = unary + best_cost;
            back[k][j] = best_prev;
        }
    }

    let last = n - 1;
    let mut best_j = 0usize;
    for j in 1..dp[last].len() {
        if dp[last][j] < dp[last][best_j] {
            best_j = j;
        }
    }
    let mut chosen = vec![0.0_f64; n];
    let mut k = last;
    loop {
        chosen[k] = col_cands[k][best_j].cy;
        if k == 0 {
            break;
        }
        best_j = back[k][best_j];
        k -= 1;
    }

    let mut ys_raw = vec![f64::NAN; cw];
    for (idx, &col) in valid_cols.iter().enumerate() {
        ys_raw[col] = chosen[idx];
    }

    // ------------------------------------ interpolation + median filtering
    let first = valid_cols[0];
    let last_col = *valid_cols.last().unwrap();
    let first_v = ys_raw[first];
    let last_v = ys_raw[last_col];
    for v in ys_raw.iter_mut().take(first) {
        *v = first_v;
    }
    for v in ys_raw.iter_mut().skip(last_col + 1) {
        *v = last_v;
    }

    let mut col = 0usize;
    while col < cw {
        if ys_raw[col].is_finite() {
            col += 1;
            continue;
        }
        let gap_start = col;
        let mut gap_end = col;
        while gap_end < cw && !ys_raw[gap_end].is_finite() {
            gap_end += 1;
        }
        if gap_end >= cw {
            break;
        }
        let a = ys_raw[gap_start - 1];
        let b = ys_raw[gap_end];
        let span = (gap_end - (gap_start - 1)) as f64;
        for (j, slot) in ys_raw.iter_mut().enumerate().take(gap_end).skip(gap_start) {
            let t = (j - (gap_start - 1)) as f64 / span;
            *slot = a + (b - a) * t;
        }
        col = gap_end;
    }

    let ys_clean: Vec<f64> = (0..cw)
        .map(|i| {
            median3(
                ys_raw[i.saturating_sub(1)],
                ys_raw[i],
                ys_raw[(i + 1).min(cw - 1)],
            )
        })
        .collect();

    // ------------------------------------------------ resample + calibrate
    let n_out = cw.clamp(400, 1600);
    let (y_bottom, y_top) = if y_unit == YUnit::Transmittance {
        (0.0_f64, 100.0_f64)
    } else {
        (0.0_f64, 1.0_f64)
    };
    let denom = ch.saturating_sub(1).max(1) as f64;

    let mut points: Vec<(f64, f64)> = Vec::with_capacity(n_out);
    for k in 0..n_out {
        let u = if n_out > 1 {
            k as f64 / (n_out - 1) as f64
        } else {
            0.0
        };
        let pos = u * (cw - 1) as f64;
        let i0 = pos.floor() as usize;
        let i1 = (i0 + 1).min(cw - 1);
        let t = pos - i0 as f64;
        let y_resamp = ys_clean[i0] + (ys_clean[i1] - ys_clean[i0]) * t;
        let phys_x = x_left + u * (x_right - x_left);
        let phys_y = y_top - (y_resamp / denom) * (y_top - y_bottom);
        if phys_x.is_finite() && phys_y.is_finite() {
            points.push((phys_x, phys_y));
        }
    }

    points.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
    let mut x_out: Vec<f64> = Vec::with_capacity(points.len());
    let mut y_out: Vec<f64> = Vec::with_capacity(points.len());
    for (x, y) in points {
        if x_out.last().map(|&prev| x > prev + 1e-9).unwrap_or(true) {
            x_out.push(x);
            y_out.push(y);
        }
    }

    if x_out.len() < 2 {
        bail!("Digitized curve has too few points");
    }

    let stem = img_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("Spectrum")
        .to_string();

    Ok(Spectrum {
        name: format!("{} [Image]", stem),
        source_path: Some(img_path.display().to_string()),
        kind,
        raw_y_unit: y_unit,
        x: x_out,
        y: y_out,
        visible: true,
        color_idx: 0,
        metadata: vec![
            ("Digitizer".into(), "Native Rust (Viterbi DP)".into()),
            ("Crop Box".into(), format!("({x0}, {y0}, {x1}, {y1})")),
        ],
    })
}
