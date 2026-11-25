use std::fmt::Write as _;

use anyhow::{bail, Result};

use crate::data::{best_peak_label_in_context, SolventId, SpectrumKind, YUnit};
use crate::processing::{
    savitzky_golay, savitzky_golay_deriv2, transmittance_to_absorbance, SmoothLevel,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeakShape {
    PseudoVoigt,
    Gaussian,
    Lorentzian,
}

impl PeakShape {
    pub fn next(self) -> Self {
        match self {
            Self::PseudoVoigt => Self::Gaussian,
            Self::Gaussian => Self::Lorentzian,
            Self::Lorentzian => Self::PseudoVoigt,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::PseudoVoigt => "Pseudo-Voigt",
            Self::Gaussian => "Gaussian",
            Self::Lorentzian => "Lorentzian",
        }
    }

    pub fn prev(self) -> Self {
        match self {
            Self::PseudoVoigt => Self::Lorentzian,
            Self::Gaussian => Self::PseudoVoigt,
            Self::Lorentzian => Self::Gaussian,
        }
    }
}

#[derive(Debug, Clone)]
pub struct DeconvConfig {
    pub shape: PeakShape,
    pub max_peaks: usize,
    pub fit_baseline: bool,
    pub manual_centers: Vec<f64>,
    pub max_iter: usize,
}

impl Default for DeconvConfig {
    fn default() -> Self {
        Self {
            shape: PeakShape::PseudoVoigt,
            max_peaks: 3,
            fit_baseline: true,
            manual_centers: Vec::new(),
            max_iter: 150,
        }
    }
}

#[derive(Debug, Clone)]
pub struct FittedPeak {
    pub center: f64,
    pub err_center: f64,
    pub height: f64,
    pub err_height: f64,
    pub fwhm: f64,
    pub err_fwhm: f64,
    pub eta: f64,
    pub err_eta: f64,
    pub area: f64,
    pub area_pct: f64,
    pub label: String,
}

#[derive(Debug, Clone)]
pub struct DeconvResult {
    pub kind: SpectrumKind,
    pub shape: PeakShape,
    pub x_lo: f64,
    pub x_hi: f64,
    pub x: Vec<f64>,
    pub y_obs: Vec<f64>,
    pub y_fit: Vec<f64>,
    pub y_baseline: Vec<f64>,
    pub y_residual: Vec<f64>,
    pub components: Vec<Vec<f64>>,
    pub peaks: Vec<FittedPeak>,
    pub baseline_b0: f64,
    pub baseline_b1: f64,
    pub r_squared: f64,
    pub adj_r_squared: f64,
    pub reduced_chi2: f64,
    pub rmse: f64,
    pub iterations: usize,
    pub converged: bool,
}

/// Dense linear system solver `A x = b` via Gaussian elimination with partial pivoting.
fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    if n == 0 || a.len() != n || a.iter().any(|row| row.len() != n) {
        return None;
    }

    let mut mat = a.to_vec();
    let mut rhs = b.to_vec();

    for col in 0..n {
        let mut pivot_row = col;
        let mut pivot_val = mat[col][col].abs();
        for (r, row) in mat.iter().enumerate().skip(col + 1) {
            let v = row[col].abs();
            if v > pivot_val {
                pivot_val = v;
                pivot_row = r;
            }
        }
        if !pivot_val.is_finite() || pivot_val < 1e-18 {
            return None;
        }
        if pivot_row != col {
            mat.swap(col, pivot_row);
            rhs.swap(col, pivot_row);
        }

        let diag = mat[col][col];
        for r in (col + 1)..n {
            let factor = mat[r][col] / diag;
            if factor == 0.0 {
                continue;
            }
            mat[r][col] = 0.0;
            for c in (col + 1)..n {
                mat[r][c] -= factor * mat[col][c];
            }
            rhs[r] -= factor * rhs[col];
        }
    }

    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut sum = rhs[i];
        for c in (i + 1)..n {
            sum -= mat[i][c] * x[c];
        }
        let d = mat[i][i];
        if !d.is_finite() || d.abs() < 1e-18 {
            return None;
        }
        x[i] = sum / d;
    }

    if x.iter().all(|v| v.is_finite()) {
        Some(x)
    } else {
        None
    }
}

/// 2nd-derivative hidden-peak auto-seeder for overlapping bands and shoulders.
pub fn auto_seed_peaks_deriv2(xs: &[f64], ys: &[f64], max_peaks: usize) -> Vec<f64> {
    let n = xs.len().min(ys.len());
    if n == 0 || max_peaks == 0 {
        return Vec::new();
    }
    if n < 5 {
        let (idx, _) = ys[..n]
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .unwrap_or((0, &0.0));
        return vec![xs[idx]];
    }

    let sm = savitzky_golay(&ys[..n], SmoothLevel::Win5);
    let d2 = savitzky_golay_deriv2(&sm, n > 200);
    let neg: Vec<f64> = d2.into_iter().map(|v| -v).collect();

    let neg_max = neg.iter().copied().fold(0.0_f64, f64::max);
    let y_min = ys[..n].iter().copied().fold(f64::INFINITY, f64::min);
    let y_max = ys[..n].iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let y_span = (y_max - y_min).max(1e-12);

    let mut candidates: Vec<(f64, f64)> = Vec::new();
    for i in 1..(n - 1) {
        let is_d2_max = neg[i] >= neg[i - 1]
            && neg[i] >= neg[i + 1]
            && (neg[i] > neg[i - 1] || neg[i] > neg[i + 1]);
        let is_y_max = ys[i] >= ys[i - 1]
            && ys[i] >= ys[i + 1]
            && (ys[i] > ys[i - 1] || ys[i] > ys[i + 1]);

        let above_floor = ys[i] > y_min + 0.04 * y_span;
        if (is_d2_max && neg_max > 1e-15 && neg[i] > 0.04 * neg_max && above_floor)
            || (is_y_max && above_floor)
        {
            let neg_norm = if neg_max > 1e-15 {
                (neg[i] / neg_max).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let height_norm = ((ys[i] - y_min) / y_span).clamp(0.0, 1.0);
            let score = 0.65 * neg_norm + 0.35 * height_norm;
            candidates.push((xs[i], score));
        }
    }

    candidates.sort_by(|a, b| b.1.total_cmp(&a.1));

    let span = (xs[n - 1] - xs[0]).abs().max(1e-9);
    let dx = span / (n.saturating_sub(1).max(1)) as f64;
    let min_sep = (span / (max_peaks as f64 * 3.5)).max(dx * 2.5);

    let mut selected: Vec<f64> = Vec::new();
    for (cx, _) in candidates {
        if selected.iter().all(|&s| (s - cx).abs() >= min_sep) {
            selected.push(cx);
            if selected.len() >= max_peaks {
                break;
            }
        }
    }

    if selected.is_empty() {
        let (idx, _) = ys[..n]
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .unwrap_or((0, &0.0));
        selected.push(xs[idx]);
    }

    selected.sort_by(|a, b| a.total_cmp(b));
    selected
}

#[allow(clippy::too_many_arguments)]
pub fn deconvolve_window(
    kind: SpectrumKind,
    xs: &[f64],
    ys: &[f64],
    y_unit: YUnit,
    x_lo: f64,
    x_hi: f64,
    cfg: &DeconvConfig,
    solvent: SolventId,
) -> Result<DeconvResult> {
    let n_total = xs.len().min(ys.len());
    if n_total < 8 {
        bail!("Need at least 8 data points for peak deconvolution.");
    }

    let lo = x_lo.min(x_hi);
    let hi = x_lo.max(x_hi);

    let mut wx = Vec::new();
    let mut wy = Vec::new();
    for i in 0..n_total {
        let x = xs[i];
        let y = ys[i];
        if x >= lo && x <= hi && x.is_finite() && y.is_finite() {
            let y_lin = if y_unit == YUnit::Transmittance {
                transmittance_to_absorbance(y)
            } else {
                y
            };
            wx.push(x);
            wy.push(y_lin);
        }
    }

    let n = wx.len();
    if n < 8 {
        bail!("Zoom window [{:.2}, {:.2}] contains fewer than 8 points.", lo, hi);
    }

    let actual_lo = wx[0];
    let actual_hi = wx[n - 1];
    let dx = (actual_hi - actual_lo).max(1e-9);
    let y_min = wy.iter().copied().fold(f64::INFINITY, f64::min);
    let y_max = wy.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let dy = (y_max - y_min).max(1e-9);

    let u: Vec<f64> = wx.iter().map(|&x| (x - actual_lo) / dx).collect();
    let yn: Vec<f64> = wy.iter().map(|&y| (y - y_min) / dy).collect();

    let max_peaks = cfg.max_peaks.clamp(1, 8);
    let mut centers: Vec<f64> = if cfg.manual_centers.is_empty() {
        auto_seed_peaks_deriv2(&wx, &wy, max_peaks)
    } else {
        cfg.manual_centers
            .iter()
            .copied()
            .filter(|&c| c >= actual_lo && c <= actual_hi)
            .take(max_peaks)
            .collect()
    };
    if centers.is_empty() {
        centers = auto_seed_peaks_deriv2(&wx, &wy, max_peaks);
    }

    let k = centers.len();
    let per = match cfg.shape {
        PeakShape::PseudoVoigt => 4,
        PeakShape::Gaussian | PeakShape::Lorentzian => 3,
    };
    let n_base = if cfg.fit_baseline { 2 } else { 0 };
    let n_params = k * per + n_base;

    let du_min = (1.5 / n as f64).max(1e-4);
    let default_w = (0.8 / (k as f64 + 1.0)).clamp(du_min * 2.0, 0.35);

    let mut p = vec![0.0_f64; n_params];
    let b0_init = (yn[0].min(yn[n - 1])).clamp(-0.2, 0.5);
    for (pi, &cx) in centers.iter().enumerate() {
        let u0 = ((cx - actual_lo) / dx).clamp(0.01, 0.99);
        let idx = wx.partition_point(|&v| v < cx).min(n - 1);
        let h_init = (yn[idx] - b0_init).clamp(0.05, 1.5);
        let base = pi * per;
        p[base] = h_init;
        p[base + 1] = u0;
        p[base + 2] = default_w;
        if per == 4 {
            p[base + 3] = 0.5;
        }
    }
    if cfg.fit_baseline {
        let b1_init = (yn[n - 1] - yn[0]).clamp(-0.5, 0.5);
        let b0_mid = (0.5 * (yn[0] + yn[n - 1])).clamp(-0.2, 0.6);
        p[k * per] = b0_mid;
        p[k * per + 1] = b1_init;
    }

    let clamp_params = |params: &mut [f64]| {
        for pi in 0..k {
            let base = pi * per;
            params[base] = params[base].clamp(0.0, 2.5);
            params[base + 1] = params[base + 1].clamp(0.0, 1.0);
            params[base + 2] = params[base + 2].clamp(du_min, 1.5);
            if per == 4 {
                params[base + 3] = params[base + 3].clamp(0.0, 1.0);
            }
        }
        if cfg.fit_baseline {
            let b = k * per;
            params[b] = params[b].clamp(-0.5, 0.8);
            params[b + 1] = params[b + 1].clamp(-1.0, 1.0);
        }
    };

    let eval_model_and_jac = |params: &[f64], need_jac: bool| -> (Vec<f64>, Vec<Vec<f64>>, f64) {
        let mut y_hat = vec![0.0_f64; n];
        let mut jac = if need_jac {
            vec![vec![0.0_f64; n_params]; n]
        } else {
            Vec::new()
        };

        let (b0, b1) = if cfg.fit_baseline {
            (params[k * per], params[k * per + 1])
        } else {
            (0.0, 0.0)
        };

        let ln2 = 2.0_f64.ln();
        let mut rss = 0.0_f64;

        for i in 0..n {
            let ui = u[i];
            let mut val = b0 + b1 * (ui - 0.5);
            if need_jac && cfg.fit_baseline {
                jac[i][k * per] = 1.0;
                jac[i][k * per + 1] = ui - 0.5;
            }
            for pi in 0..k {
                let base = pi * per;
                let h = params[base];
                let u0 = params[base + 1];
                let w = params[base + 2].max(du_min);
                let eta = match cfg.shape {
                    PeakShape::PseudoVoigt => params[base + 3],
                    PeakShape::Gaussian => 0.0,
                    PeakShape::Lorentzian => 1.0,
                };
                let z = 2.0 * (ui - u0) / w;
                let z2 = z * z;
                let l = 1.0 / (1.0 + z2);
                let g = (-ln2 * z2).exp();
                let s = eta * l + (1.0 - eta) * g;
                val += h * s;

                if need_jac {
                    let ds_dz = -2.0 * eta * z * l * l - 2.0 * ln2 * (1.0 - eta) * z * g;
                    jac[i][base] = s;
                    jac[i][base + 1] = h * ds_dz * (-2.0 / w);
                    jac[i][base + 2] = h * ds_dz * (-z / w);
                    if per == 4 {
                        jac[i][base + 3] = h * (l - g);
                    }
                }
            }
            y_hat[i] = val;
            let r = yn[i] - val;
            rss += r * r;
        }
        (y_hat, jac, rss)
    };

    let mut lam = 1e-3_f64;
    let mut iterations = 0usize;
    let mut converged = false;
    for iter in 0..cfg.max_iter {
        iterations = iter + 1;
        let (y_hat, jac, rss) = eval_model_and_jac(&p, true);

        let mut jtj = vec![vec![0.0_f64; n_params]; n_params];
        let mut jtr = vec![0.0_f64; n_params];
        for i in 0..n {
            let ri = yn[i] - y_hat[i];
            for a in 0..n_params {
                let ja = jac[i][a];
                if ja == 0.0 {
                    continue;
                }
                jtr[a] += ja * ri;
                for b in a..n_params {
                    jtj[a][b] += ja * jac[i][b];
                }
            }
        }
        for a in 0..n_params {
            for b in 0..a {
                jtj[a][b] = jtj[b][a];
            }
        }

        let mut a_mat = jtj.clone();
        for a in 0..n_params {
            let diag = jtj[a][a].max(1e-8);
            a_mat[a][a] += lam * diag;
        }

        if let Some(dp) = solve_linear(&a_mat, &jtr) {
            let mut p_new: Vec<f64> = p.iter().zip(&dp).map(|(a, b)| a + b).collect();
            clamp_params(&mut p_new);
            let (_, _, rss_new) = eval_model_and_jac(&p_new, false);
            if rss_new < rss {
                let rel_impr = (rss - rss_new) / rss.max(1e-18);
                p = p_new;
                lam = (lam * 0.35).max(1e-9);
                if rel_impr < 1e-10 {
                    converged = true;
                    break;
                }
            } else {
                lam = (lam * 5.0).min(1e9);
                if lam >= 1e9 {
                    break;
                }
            }
        } else {
            lam = (lam * 10.0).min(1e9);
            if lam >= 1e9 {
                break;
            }
        }
    }

    // Final evaluation & covariance matrix for parameter standard errors
    let (y_hat_norm, jac, final_rss_norm) = eval_model_and_jac(&p, true);
    let dof = (n as isize - n_params as isize).max(1) as f64;
    let chi2_norm = final_rss_norm / dof;

    let mut jtj = vec![vec![0.0_f64; n_params]; n_params];
    for i in 0..n {
        for a in 0..n_params {
            let ja = jac[i][a];
            if ja == 0.0 {
                continue;
            }
            for b in a..n_params {
                jtj[a][b] += ja * jac[i][b];
            }
        }
    }
    for a in 0..n_params {
        for b in 0..a {
            jtj[a][b] = jtj[b][a];
        }
        jtj[a][a] += 1e-10;
    }

    let mut param_se = vec![0.0_f64; n_params];
    for col in 0..n_params {
        let mut e = vec![0.0_f64; n_params];
        e[col] = 1.0;
        if let Some(inv_col) = solve_linear(&jtj, &e) {
            param_se[col] = (inv_col[col].max(0.0) * chi2_norm).sqrt();
        }
    }

    let (b0_n, b1_n) = if cfg.fit_baseline {
        (p[k * per], p[k * per + 1])
    } else {
        (0.0, 0.0)
    };

    let y_baseline: Vec<f64> = u
        .iter()
        .map(|&ui| y_min + dy * (b0_n + b1_n * (ui - 0.5)))
        .collect();
    let y_fit: Vec<f64> = y_hat_norm.iter().map(|&yh| y_min + dy * yh).collect();
    let y_residual: Vec<f64> = wy.iter().zip(&y_fit).map(|(&o, &f)| o - f).collect();

    let baseline_b0 = y_min + dy * b0_n;
    let baseline_b1 = dy * b1_n / dx;

    let ln2 = 2.0_f64.ln();
    let lorentz_fac = std::f64::consts::FRAC_PI_2;
    let gauss_fac = 0.5 * (std::f64::consts::PI / ln2).sqrt();

    let mut peaks = Vec::with_capacity(k);
    let mut components = Vec::with_capacity(k);
    let mut total_area = 0.0_f64;

    for pi in 0..k {
        let base = pi * per;
        let h_n = p[base];
        let u0 = p[base + 1];
        let w_n = p[base + 2].max(du_min);
        let eta = match cfg.shape {
            PeakShape::PseudoVoigt => p[base + 3],
            PeakShape::Gaussian => 0.0,
            PeakShape::Lorentzian => 1.0,
        };

        let center = actual_lo + dx * u0;
        let height = dy * h_n;
        let fwhm = dx * w_n;

        let err_height = dy * param_se[base];
        let err_center = dx * param_se[base + 1];
        let err_fwhm = dx * param_se[base + 2];
        let err_eta = if per == 4 { param_se[base + 3] } else { 0.0 };

        let area = height * fwhm * (eta * lorentz_fac + (1.0 - eta) * gauss_fac);
        total_area += area;

        let mut comp = Vec::with_capacity(n);
        for i in 0..n {
            let z = 2.0 * (u[i] - u0) / w_n;
            let z2 = z * z;
            let l = 1.0 / (1.0 + z2);
            let g = (-ln2 * z2).exp();
            let s = eta * l + (1.0 - eta) * g;
            comp.push(y_baseline[i] + height * s);
        }

        peaks.push(FittedPeak {
            center,
            err_center,
            height,
            err_height,
            fwhm,
            err_fwhm,
            eta,
            err_eta,
            area,
            area_pct: 0.0,
            label: String::new(),
        });
        components.push(comp);
    }

    for pk in &mut peaks {
        pk.area_pct = if total_area > 1e-18 {
            100.0 * pk.area / total_area
        } else {
            0.0
        };
    }

    // Sort peaks & components in lockstep in display order
    let mut order: Vec<usize> = (0..k).collect();
    if kind.is_x_reversed() {
        order.sort_by(|&a, &b| peaks[b].center.total_cmp(&peaks[a].center));
    } else {
        order.sort_by(|&a, &b| peaks[a].center.total_cmp(&peaks[b].center));
    }
    let peaks_sorted: Vec<FittedPeak> = order.iter().map(|&i| peaks[i].clone()).collect();
    let comps_sorted: Vec<Vec<f64>> = order.iter().map(|&i| components[i].clone()).collect();
    let mut peaks = peaks_sorted;
    let components = comps_sorted;

    let centers_ctx: Vec<f64> = peaks.iter().map(|pk| pk.center).collect();
    for pk in &mut peaks {
        pk.label = best_peak_label_in_context(kind, pk.center, &centers_ctx, solvent);
    }

    let mean_y = wy.iter().sum::<f64>() / n as f64;
    let tss: f64 = wy.iter().map(|&y| (y - mean_y) * (y - mean_y)).sum();
    let rss_phys: f64 = y_residual.iter().map(|&r| r * r).sum();
    let r_squared = if tss > 1e-18 {
        (1.0 - rss_phys / tss).clamp(0.0, 1.0)
    } else {
        1.0
    };
    let adj_r_squared = if n > n_params + 1 && tss > 1e-18 {
        (1.0 - (1.0 - r_squared) * ((n - 1) as f64 / dof)).clamp(0.0, 1.0)
    } else {
        r_squared
    };
    let reduced_chi2 = rss_phys / dof;
    let rmse = (rss_phys / n as f64).sqrt();

    Ok(DeconvResult {
        kind,
        shape: cfg.shape,
        x_lo: actual_lo,
        x_hi: actual_hi,
        x: wx,
        y_obs: wy,
        y_fit,
        y_baseline,
        y_residual,
        components,
        peaks,
        baseline_b0,
        baseline_b1,
        r_squared,
        adj_r_squared,
        reduced_chi2,
        rmse,
        iterations,
        converged,
    })
}

pub fn format_deconv_markdown(res: &DeconvResult) -> String {
    let mut out = String::new();
    let _ = writeln!(
        out,
        "### Peak Deconvolution ({} · {})",
        res.kind.display_name(),
        res.shape.label()
    );
    let x_mid = 0.5 * (res.x_lo + res.x_hi);
    let _ = writeln!(
        out,
        "Window: {:.2}–{:.2} {} | R² = {:.6} | Adj-R² = {:.6} | RMSE = {:.5} | χ²_red = {:.6} | Baseline = {:.5} + {:.5}*(x - {:.2}) | Iters = {} ({})\n",
        res.x_lo,
        res.x_hi,
        res.kind.x_unit(),
        res.r_squared,
        res.adj_r_squared,
        res.rmse,
        res.reduced_chi2,
        res.baseline_b0,
        res.baseline_b1,
        x_mid,
        res.iterations,
        if res.converged { "converged" } else { "max-iter" }
    );
    let _ = writeln!(
        out,
        "| # | Center ± SE ({}) | FWHM ± SE | Height | Area | Area % | η | Assignment |",
        res.kind.x_unit()
    );
    let _ = writeln!(out, "|---|---|---|---|---|---|---|---|");
    for (i, pk) in res.peaks.iter().enumerate() {
        let _ = writeln!(
            out,
            "| {} | {:.2} ± {:.2} | {:.2} ± {:.2} | {:.4} | {:.4} | {:.1}% | {:.2} ± {:.2} | {} |",
            i + 1,
            pk.center,
            pk.err_center,
            pk.fwhm,
            pk.err_fwhm,
            pk.height,
            pk.area,
            pk.area_pct,
            pk.eta,
            pk.err_eta,
            pk.label
        );
    }
    out
}

pub fn format_deconv_csv(res: &DeconvResult) -> String {
    let mut out = String::new();
    let _ = writeln!(
        out,
        "# SPECTRA Peak Deconvolution Report ({}, {})",
        res.kind.display_name(),
        res.shape.label()
    );
    let _ = writeln!(
        out,
        "# Window: {:.4} to {:.4} {}, R2={:.6}, AdjR2={:.6}, RMSE={:.6}, RedChi2={:.6}",
        res.x_lo,
        res.x_hi,
        res.kind.x_unit(),
        res.r_squared,
        res.adj_r_squared,
        res.rmse,
        res.reduced_chi2
    );
    let _ = writeln!(
        out,
        "# Peak,Center,ErrCenter,FWHM,ErrFWHM,Height,ErrHeight,Area,AreaPct,Eta,Assignment"
    );
    for (i, pk) in res.peaks.iter().enumerate() {
        let _ = writeln!(
            out,
            "# {},{:.4},{:.4},{:.4},{:.4},{:.6},{:.6},{:.6},{:.2},{:.4},\"{}\"",
            i + 1,
            pk.center,
            pk.err_center,
            pk.fwhm,
            pk.err_fwhm,
            pk.height,
            pk.err_height,
            pk.area,
            pk.area_pct,
            pk.eta,
            pk.label.replace('"', "\"\"")
        );
    }

    let _ = write!(out, "X,Y_obs,Y_fit,Y_baseline,Y_residual");
    for i in 0..res.components.len() {
        let _ = write!(out, ",Comp_{}", i + 1);
    }
    let _ = writeln!(out);

    for i in 0..res.x.len() {
        let _ = write!(
            out,
            "{:.6},{:.6},{:.6},{:.6},{:.6}",
            res.x[i], res.y_obs[i], res.y_fit[i], res.y_baseline[i], res.y_residual[i]
        );
        for comp in &res.components {
            let _ = write!(out, ",{:.6}", comp[i]);
        }
        let _ = writeln!(out);
    }
    out
}
