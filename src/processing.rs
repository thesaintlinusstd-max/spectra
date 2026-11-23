use crate::data::{best_short_peak_label, Spectrum, SpectrumKind, YUnit};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SmoothLevel {
    Off,
    Win5,
    Win9,
    Win15,
}

impl SmoothLevel {
    pub fn next(self) -> Self {
        match self {
            Self::Off => Self::Win5,
            Self::Win5 => Self::Win9,
            Self::Win9 => Self::Win15,
            Self::Win15 => Self::Off,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Off => "OFF",
            Self::Win5 => "SG-5",
            Self::Win9 => "SG-9",
            Self::Win15 => "SG-15",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct ProcessConfig {
    pub ftir_absorbance: bool,
    pub baseline: bool,
    pub smooth: SmoothLevel,
    pub normalize: bool,
}

impl Default for ProcessConfig {
    fn default() -> Self {
        Self {
            ftir_absorbance: false,
            baseline: false,
            smooth: SmoothLevel::Off,
            normalize: false,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Peak {
    pub x: f64,
    pub y: f64,
    pub prominence: f64,
    pub label: String,
    pub pinned: bool,
}

pub fn transmittance_to_absorbance(pct_t: f64) -> f64 {
    2.0 - pct_t.clamp(0.01, 100.0).log10()
}

pub fn absorbance_to_transmittance(abs: f64) -> f64 {
    (100.0 * 10.0_f64.powf(-abs.max(0.0))).clamp(0.0, 100.0)
}

pub fn savitzky_golay(ys: &[f64], level: SmoothLevel) -> Vec<f64> {
    let (coeffs, norm): (&[f64], f64) = match level {
        SmoothLevel::Off => return ys.to_vec(),
        SmoothLevel::Win5 => (&[-3.0, 12.0, 17.0, 12.0, -3.0], 35.0),
        SmoothLevel::Win9 => (&[-21.0, 14.0, 39.0, 54.0, 59.0, 54.0, 39.0, 14.0, -21.0], 231.0),
        SmoothLevel::Win15 => (
            &[
                -78.0, -13.0, 42.0, 87.0, 122.0, 147.0, 162.0, 167.0, 162.0, 147.0, 122.0, 87.0,
                42.0, -13.0, -78.0,
            ],
            1105.0,
        ),
    };

    let n = ys.len();
    let half = coeffs.len() / 2;
    if n < coeffs.len() {
        return ys.to_vec();
    }

    let mut out = vec![0.0; n];
    for i in 0..n {
        let mut acc = 0.0;
        for (k, &w) in coeffs.iter().enumerate() {
            let offset = k as isize - half as isize;
            let mut idx = i as isize + offset;
            if idx < 0 {
                idx = -idx;
            }
            if idx >= n as isize {
                idx = 2 * (n as isize - 1) - idx;
            }
            let clamped = idx.clamp(0, (n - 1) as isize) as usize;
            acc += w * ys[clamped];
        }
        out[i] = acc / norm;
    }
    out
}

/// Bruker OPUS-style Rubberband (Convex Hull) Baseline Correction.
/// `upper = true` for %Transmittance (peaks point down), `false` for Absorbance/NMR (peaks point up).
pub fn rubberband_baseline(xs: &[f64], ys: &[f64], upper: bool) -> Vec<f64> {
    let n = xs.len().min(ys.len());
    if n < 3 {
        return ys.to_vec();
    }

    // Build monotone chain hull indices
    let mut hull: Vec<usize> = Vec::with_capacity(64);
    for i in 0..n {
        while hull.len() >= 2 {
            let a = hull[hull.len() - 2];
            let b = hull[hull.len() - 1];
            let cross = (xs[b] - xs[a]) * (ys[i] - ys[a]) - (ys[b] - ys[a]) * (xs[i] - xs[a]);
            let pop = if upper { cross >= 0.0 } else { cross <= 0.0 };
            if pop {
                hull.pop();
            } else {
                break;
            }
        }
        hull.push(i);
    }

    if hull.len() < 2 {
        return ys.to_vec();
    }

    // Linearly interpolate baseline between hull vertices
    let mut baseline = vec![0.0; n];
    for seg in hull.windows(2) {
        let a = seg[0];
        let b = seg[1];
        let dx = (xs[b] - xs[a]).max(1e-12);
        let dy = ys[b] - ys[a];
        for i in a..=b {
            let t = (xs[i] - xs[a]) / dx;
            baseline[i] = ys[a] + t * dy;
        }
    }

    if upper {
        // Shift %T so rubberband baseline sits cleanly at 100%T
        ys.iter()
            .zip(&baseline)
            .map(|(&y, &b)| (100.0 - (b - y)).clamp(0.0, 100.0))
            .collect()
    } else {
        // Shift Absorbance/NMR so rubberband baseline sits cleanly at 0
        ys.iter()
            .zip(&baseline)
            .map(|(&y, &b)| (y - b).max(0.0))
            .collect()
    }
}

pub fn process_trace(spec: &Spectrum, cfg: ProcessConfig) -> (Vec<f64>, Vec<f64>, YUnit) {
    let xs = spec.x.clone();
    let (mut ys, y_unit) = match spec.kind {
        SpectrumKind::Ftir => {
            if cfg.ftir_absorbance {
                let converted = match spec.raw_y_unit {
                    YUnit::Transmittance => {
                        spec.y.iter().copied().map(transmittance_to_absorbance).collect()
                    }
                    _ => spec.y.clone(),
                };
                (converted, YUnit::Absorbance)
            } else {
                let converted = match spec.raw_y_unit {
                    YUnit::Absorbance => {
                        spec.y.iter().copied().map(absorbance_to_transmittance).collect()
                    }
                    _ => spec.y.clone(),
                };
                (converted, YUnit::Transmittance)
            }
        }
        SpectrumKind::Hnmr | SpectrumKind::Cnmr => (spec.y.clone(), YUnit::Intensity),
        SpectrumKind::UvVis => {
            // Default = Absorbance; `ftir_absorbance` toggles %T <-> Abs.
            let want_transmittance = cfg.ftir_absorbance;
            if want_transmittance {
                let converted = match spec.raw_y_unit {
                    YUnit::Absorbance => {
                        spec.y.iter().copied().map(absorbance_to_transmittance).collect()
                    }
                    _ => spec.y.clone(),
                };
                (converted, YUnit::Transmittance)
            } else {
                let converted = match spec.raw_y_unit {
                    YUnit::Transmittance => {
                        spec.y.iter().copied().map(transmittance_to_absorbance).collect()
                    }
                    _ => spec.y.clone(),
                };
                (converted, YUnit::Absorbance)
            }
        }
    };

    if cfg.smooth != SmoothLevel::Off {
        ys = savitzky_golay(&ys, cfg.smooth);
    }

    if cfg.baseline {
        let upper = y_unit == YUnit::Transmittance;
        ys = rubberband_baseline(&xs, &ys, upper);
    }

    if cfg.normalize && ys.len() >= 2 {
        let y_min = ys.iter().copied().fold(f64::INFINITY, f64::min);
        let y_max = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let span = (y_max - y_min).max(1e-12);
        if y_unit == YUnit::Transmittance {
            for v in &mut ys {
                *v = 10.0 + 90.0 * (*v - y_min) / span;
            }
        } else {
            for v in &mut ys {
                *v = (*v - y_min) / span;
            }
        }
    }

    (xs, ys, y_unit)
}

pub fn interpolate_y(xs: &[f64], ys: &[f64], target_x: f64) -> Option<f64> {
    if xs.is_empty() || xs.len() != ys.len() {
        return None;
    }
    if target_x <= xs[0] {
        return Some(ys[0]);
    }
    if target_x >= xs[xs.len() - 1] {
        return Some(ys[ys.len() - 1]);
    }
    let idx = xs.partition_point(|&v| v < target_x);
    if idx == 0 {
        return Some(ys[0]);
    }
    let x0 = xs[idx - 1];
    let x1 = xs[idx];
    let t = ((target_x - x0) / (x1 - x0).max(1e-12)).clamp(0.0, 1.0);
    Some(ys[idx - 1] + t * (ys[idx] - ys[idx - 1]))
}

pub fn detect_peaks(
    xs: &[f64],
    ys: &[f64],
    kind: SpectrumKind,
    y_unit: YUnit,
    rel_prominence: f64,
    max_peaks: usize,
) -> Vec<Peak> {
    let n = xs.len().min(ys.len());
    if n < 5 {
        return Vec::new();
    }

    // Work on upward-pointing signal `s`
    let invert = y_unit == YUnit::Transmittance;
    let s: Vec<f64> = if invert {
        ys.iter().map(|&v| -v).collect()
    } else {
        ys.to_vec()
    };

    let s_min = s.iter().copied().fold(f64::INFINITY, f64::min);
    let s_max = s.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let span = (s_max - s_min).max(1e-12);
    let min_prom = rel_prominence.clamp(0.01, 0.9) * span;

    let mut candidates: Vec<(usize, f64)> = Vec::new();
    for i in 1..(n - 1) {
        if s[i] > s[i - 1] && s[i] >= s[i + 1] {
            // Left scan until higher point or edge
            let mut left_min = s[i];
            let mut l = i;
            while l > 0 {
                l -= 1;
                if s[l] < left_min {
                    left_min = s[l];
                }
                if s[l] > s[i] {
                    break;
                }
            }
            // Right scan until higher point or edge
            let mut right_min = s[i];
            let mut r = i;
            while r + 1 < n {
                r += 1;
                if s[r] < right_min {
                    right_min = s[r];
                }
                if s[r] > s[i] {
                    break;
                }
            }
            let prom = s[i] - left_min.max(right_min);
            if prom >= min_prom {
                candidates.push((i, prom));
            }
        }
    }

    // Sort strongest prominence first, enforce minimum X distance
    candidates.sort_by(|a, b| b.1.total_cmp(&a.1));
    let min_dx = match kind {
        SpectrumKind::Ftir => 14.0,
        SpectrumKind::Hnmr => 0.06,
        SpectrumKind::Cnmr => 0.7,
        SpectrumKind::UvVis => 10.0,
    };

    let mut chosen: Vec<Peak> = Vec::new();
    for (idx, prom) in candidates {
        // 3-point parabolic refinement
        let y_l = s[idx - 1];
        let y_c = s[idx];
        let y_r = s[idx + 1];
        let denom = y_l - 2.0 * y_c + y_r;
        let delta = if denom.abs() > 1e-12 {
            (0.5 * (y_l - y_r) / denom).clamp(-0.5, 0.5)
        } else {
            0.0
        };
        let dx_step = 0.5 * (xs[idx + 1] - xs[idx - 1]);
        let refined_x = xs[idx] + delta * dx_step;

        if chosen.iter().any(|p| (p.x - refined_x).abs() < min_dx) {
            continue;
        }
        chosen.push(Peak {
            x: refined_x,
            y: ys[idx],
            prominence: prom / span,
            label: best_short_peak_label(kind, refined_x),
            pinned: false,
        });
        if chosen.len() >= max_peaks {
            break;
        }
    }

    // FTIR / NMR plot descending in X (IUPAC convention); UV-Vis ascending (nm).
    if kind.is_x_reversed() {
        chosen.sort_by(|a, b| b.x.total_cmp(&a.x));
    } else {
        chosen.sort_by(|a, b| a.x.total_cmp(&b.x));
    }
    chosen
}

/// Returns (baseline-corrected integrated area, baseline-corrected peak height) on Absorbance/Intensity
pub fn integrate_region(xs: &[f64], ys: &[f64], x_lo: f64, x_hi: f64) -> Option<(f64, f64)> {
    if xs.len() < 2 || x_lo >= x_hi {
        return None;
    }
    let i0 = xs.partition_point(|&v| v < x_lo).min(xs.len() - 1);
    let i1 = xs.partition_point(|&v| v <= x_hi).saturating_sub(1).max(i0);
    if i1 <= i0 {
        return None;
    }

    let x_a = xs[i0];
    let y_a = ys[i0];
    let x_b = xs[i1];
    let y_b = ys[i1];
    let span_x = (x_b - x_a).max(1e-12);

    let mut area = 0.0;
    let mut max_h = 0.0_f64;
    for k in i0..i1 {
        let base0 = y_a + (y_b - y_a) * ((xs[k] - x_a) / span_x);
        let base1 = y_a + (y_b - y_a) * ((xs[k + 1] - x_a) / span_x);
        let h0 = (ys[k] - base0).max(0.0);
        let h1 = (ys[k + 1] - base1).max(0.0);
        max_h = max_h.max(h0).max(h1);
        area += 0.5 * (h0 + h1) * (xs[k + 1] - xs[k]);
    }
    Some((area, max_h))
}

// ============================================================================
// POLYMER CHEMISTRY CALCULATORS
// ============================================================================

#[derive(Debug, Clone, Copy)]
pub struct ConversionResult {
    pub area_conversion_pct: f64,
    pub height_conversion_pct: f64,
    pub ratio_0_area: f64,
    pub ratio_t_area: f64,
    pub ratio_0_height: f64,
    pub ratio_t_height: f64,
}

pub fn ftir_conversion_percent(
    a_rxn_0: f64,
    a_ref_0: f64,
    a_rxn_t: f64,
    a_ref_t: f64,
) -> Option<f64> {
    if a_rxn_0 <= 0.0 || a_ref_0 <= 0.0 || a_ref_t <= 0.0 || a_rxn_t < 0.0 {
        return None;
    }
    let r0 = a_rxn_0 / a_ref_0;
    let rt = a_rxn_t / a_ref_t;
    Some(((1.0 - rt / r0) * 100.0).clamp(-100.0, 100.0))
}

pub fn auto_ftir_conversion(
    monomer: &Spectrum,
    polymer: &Spectrum,
    rxn_center: f64,
    rxn_half: f64,
    ref_center: f64,
    ref_half: f64,
) -> Option<ConversionResult> {
    let abs_cfg = ProcessConfig {
        ftir_absorbance: true,
        baseline: false,
        smooth: SmoothLevel::Off,
        normalize: false,
    };
    let (mx, my, _) = process_trace(monomer, abs_cfg);
    let (px, py, _) = process_trace(polymer, abs_cfg);

    let (m_rxn_a, m_rxn_h) = integrate_region(&mx, &my, rxn_center - rxn_half, rxn_center + rxn_half)?;
    let (m_ref_a, m_ref_h) = integrate_region(&mx, &my, ref_center - ref_half, ref_center + ref_half)?;
    let (p_rxn_a, p_rxn_h) = integrate_region(&px, &py, rxn_center - rxn_half, rxn_center + rxn_half)?;
    let (p_ref_a, p_ref_h) = integrate_region(&px, &py, ref_center - ref_half, ref_center + ref_half)?;

    let area_conv = ftir_conversion_percent(m_rxn_a, m_ref_a, p_rxn_a, p_ref_a)?;
    let height_conv = ftir_conversion_percent(m_rxn_h, m_ref_h, p_rxn_h, p_ref_h)?;

    Some(ConversionResult {
        area_conversion_pct: area_conv,
        height_conversion_pct: height_conv,
        ratio_0_area: m_rxn_a / m_ref_a,
        ratio_t_area: p_rxn_a / p_ref_a,
        ratio_0_height: m_rxn_h / m_ref_h,
        ratio_t_height: p_rxn_h / p_ref_h,
    })
}

/// Returns (DP_n, M_n in g/mol)
pub fn nmr_mn_dpn(
    m0: f64,
    m_end: f64,
    n_repeat: f64,
    n_end: f64,
    i_repeat: f64,
    i_end: f64,
) -> Option<(f64, f64)> {
    if m0 <= 0.0 || n_repeat <= 0.0 || n_end <= 0.0 || i_repeat <= 0.0 || i_end <= 0.0 {
        return None;
    }
    let dpn = (i_repeat / n_repeat) / (i_end / n_end);
    let mn = dpn * m0 + m_end.max(0.0);
    Some((dpn, mn))
}

#[derive(Debug, Clone, Copy)]
pub struct TacticityResult {
    pub pct_mm: f64,
    pub pct_mr: f64,
    pub pct_rr: f64,
    pub dyad_m: f64,
    pub dyad_r: f64,
    pub rho: f64,
    pub regime: &'static str,
}

pub fn tacticity_stats(i_mm: f64, i_mr: f64, i_rr: f64) -> Option<TacticityResult> {
    let total = i_mm + i_mr + i_rr;
    if total <= 0.0 || i_mm < 0.0 || i_mr < 0.0 || i_rr < 0.0 {
        return None;
    }
    let f_mm = i_mm / total;
    let f_mr = i_mr / total;
    let f_rr = i_rr / total;
    let m = f_mm + 0.5 * f_mr;
    let r = f_rr + 0.5 * f_mr;
    let rho = if f_mr > 1e-9 { (2.0 * m * r) / f_mr } else { f64::INFINITY };
    let regime = if (rho - 1.0).abs() <= 0.12 {
        "Bernoullian (typical free-radical / chain-end control)"
    } else if rho > 1.12 {
        "Stereoblock / Enantiomorphic-site tendency (ρ > 1)"
    } else {
        "Alternating / Heterotactic tendency (ρ < 1)"
    };

    Some(TacticityResult {
        pct_mm: f_mm * 100.0,
        pct_mr: f_mr * 100.0,
        pct_rr: f_rr * 100.0,
        dyad_m: m,
        dyad_r: r,
        rho,
        regime,
    })
}

// ============================================================================
// PVDF β-PHASE FRACTION & POLYURETHANE H-BONDING INDEX
// ============================================================================

/// Gregorio (2010) β-phase fraction from FTIR absorbance at 840 and 763 cm⁻¹:
/// `F(β) = A_840 / (1.26 · A_763 + A_840) × 100%`
pub fn pvdf_beta_fraction(a_840: f64, a_763: f64) -> Option<f64> {
    if !a_840.is_finite() || !a_763.is_finite() || a_840 < 0.0 || a_763 < 0.0 {
        return None;
    }
    let denom = 1.26 * a_763 + a_840;
    if denom <= 1e-12 {
        return None;
    }
    Some((a_840 / denom * 100.0).clamp(0.0, 100.0))
}

#[derive(Debug, Clone, Copy)]
pub struct HydrogenBondIndex {
    pub a_bonded: f64,
    pub a_free: f64,
    pub hbi: f64,
    pub bonded_pct: f64,
}

/// Polyurethane carbonyl hydrogen-bonding index:
/// `HBI = A_1703 / A_1732` and bonded fraction `= HBI / (1 + HBI) × 100%`
pub fn pu_hydrogen_bond_index(a_1703: f64, a_1732: f64) -> Option<HydrogenBondIndex> {
    if !a_1703.is_finite() || !a_1732.is_finite() || a_1703 < 0.0 || a_1732 <= 1e-12 {
        return None;
    }
    let hbi = a_1703 / a_1732;
    Some(HydrogenBondIndex {
        a_bonded: a_1703,
        a_free: a_1732,
        hbi,
        bonded_pct: (hbi / (1.0 + hbi) * 100.0).clamp(0.0, 100.0),
    })
}

// ============================================================================
// UV-VIS SEMICONDUCTOR BAND GAP (TAUC) & BEER-LAMBERT QUANTITATION
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TaucMode {
    Direct,
    Indirect,
}

impl TaucMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::Direct => "Direct (n = 1/2)",
            Self::Indirect => "Indirect (n = 2)",
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct TaucResult {
    pub band_gap_ev: f64,
    pub edge_nm: f64,
    pub r_squared: f64,
    pub mode: TaucMode,
}

/// Tauc plot band gap from a UV-Vis absorption edge.
///
/// `(αhν)^(1/n)` versus `hν` is fitted in a sliding linear window on the
/// absorption edge; the x-intercept of the best-fit line is Eg.
pub fn compute_tauc_band_gap(spec: &Spectrum, mode: TaucMode) -> Option<TaucResult> {
    if spec.kind != SpectrumKind::UvVis || spec.x.len() < 10 || spec.x.len() != spec.y.len() {
        return None;
    }

    // Ensure Absorbance
    let ys: Vec<f64> = if spec.raw_y_unit == YUnit::Transmittance {
        spec.y.iter().copied().map(transmittance_to_absorbance).collect()
    } else {
        spec.y.clone()
    };

    // Baseline from the 5th percentile of absorbance
    let mut sorted = ys.clone();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let idx5 = ((sorted.len() as f64 - 1.0) * 0.05).round() as usize;
    let baseline = sorted.get(idx5).copied().unwrap_or(0.0).max(0.0);

    let mut ev: Vec<f64> = Vec::with_capacity(spec.x.len());
    let mut tauc: Vec<f64> = Vec::with_capacity(spec.x.len());
    for (&nm, &a) in spec.x.iter().zip(&ys) {
        if !(nm.is_finite() && a.is_finite()) || nm <= 0.0 {
            continue;
        }
        let e = 1239.84193 / nm;
        let alpha = (a - baseline).max(0.0);
        // Tauc ordinate: direct allowed → (αhν)²; indirect allowed → (αhν)^½
        let ord = match mode {
            TaucMode::Direct => (alpha * e) * (alpha * e),
            TaucMode::Indirect => (alpha * e).sqrt(),
        };
        ev.push(e);
        tauc.push(ord);
    }
    if ev.len() < 10 {
        return None;
    }

    let max_tauc = tauc.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max_tauc <= 1e-12 {
        return None;
    }

    let mut best: Option<(f64, f64, f64, f64)> = None; // (score, slope, intercept, r2)
    let n = ev.len();
    let win = (n / 12).clamp(5, 16);
    let step = ((n - win) / 120).max(1);

    for start in (0..=(n - win)).step_by(step) {
        let end = start + win;
        let slice_e = &ev[start..end];
        let slice_t = &tauc[start..end];
        let m = win as f64;
        let mean_x = slice_e.iter().sum::<f64>() / m;
        let mean_y = slice_t.iter().sum::<f64>() / m;
        // Restrict to the rising mid-section of the edge: this rejects both the
        // noise-dominated Urbach tail and any saturated high-energy plateau.
        if mean_y < 0.20 * max_tauc || mean_y > 0.80 * max_tauc {
            continue;
        }
        let mut sxx = 0.0;
        let mut sxy = 0.0;
        let mut syy = 0.0;
        for (&x, &y) in slice_e.iter().zip(slice_t) {
            let dx = x - mean_x;
            let dy = y - mean_y;
            sxx += dx * dx;
            sxy += dx * dy;
            syy += dy * dy;
        }
        if sxx <= 1e-12 || syy <= 1e-12 {
            continue;
        }
        let slope = sxy / sxx;
        let intercept = mean_y - slope * mean_x;
        let r2 = (sxy * sxy / (sxx * syy)).clamp(0.0, 1.0);
        if slope <= 0.0 || r2 < 0.85 {
            continue;
        }
        let eg = -intercept / slope;
        if !(0.30..=6.50).contains(&eg) {
            continue;
        }
        // Steepest, best-correlated edge segment; ties resolve to the lowest-E
        // (first) window, i.e. the true onset.
        let score = slope * r2 * r2;
        if best.map_or(true, |(bs, _, _, _)| score > bs) {
            best = Some((score, slope, intercept, r2));
        }
    }

    let (_, slope, intercept, r2) = best?;
    let eg = -intercept / slope;
    Some(TaucResult {
        band_gap_ev: eg,
        edge_nm: 1239.84193 / eg.max(1e-9),
        r_squared: r2,
        mode,
    })
}

/// Beer–Lambert law: `A = ε · b · c` → `c = A / (ε · b)`.
/// Returns concentration in mol/L when ε is in M⁻¹·cm⁻¹ and `path_cm` is cm.
pub fn beer_lambert_conc(absorbance: f64, epsilon: f64, path_cm: f64) -> Option<f64> {
    if !absorbance.is_finite() || !epsilon.is_finite() || !path_cm.is_finite() {
        return None;
    }
    if absorbance < 0.0 || epsilon <= 0.0 || path_cm <= 0.0 {
        return None;
    }
    Some(absorbance / (epsilon * path_cm))
}
