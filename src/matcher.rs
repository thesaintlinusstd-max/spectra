//! Spectral library matching engine — 100% pure Rust, zero extra crates.
//!
//! Compares an unknown FTIR spectrum against the built-in polymer fingerprint
//! library and any user-supplied reference spectra in `library/ftir/`, using a
//! Hit Quality Index (HQI) that combines baseline-sensitive and derivative
//! (baseline-invariant) Pearson correlations.

use std::path::PathBuf;

use crate::data::{
    generate_demo_preset, synthesize_polymer_reference, DemoPreset, Spectrum, SpectrumKind,
    POLYMERS_DB,
};
use crate::io::load_spectrum_file;
use crate::processing::{process_trace, rubberband_baseline, ProcessConfig, SmoothLevel};

const GRID_LO: f64 = 520.0;
const GRID_HI: f64 = 3760.0;
const GRID_STEP: f64 = 4.0;
const CO2_LO: f64 = 2280.0;
const CO2_HI: f64 = 2410.0;
const MIN_POINTS: usize = 120;
const MAX_MATCHES: usize = 25;

const UVVIS_GRID_LO: f64 = 200.0;
const UVVIS_GRID_HI: f64 = 800.0;
const UVVIS_GRID_STEP: f64 = 2.0;
const UVVIS_MIN_POINTS: usize = 30;

#[derive(Debug, Clone)]
pub struct LibraryMatch {
    pub name: String,
    pub short: String,
    pub source: String,
    pub hqi_pct: f64,
    pub r_abs: f64,
    pub r_deriv: f64,
    pub is_builtin: bool,
    pub file_path: Option<PathBuf>,
}

fn grid() -> Vec<f64> {
    let n = ((GRID_HI - GRID_LO) / GRID_STEP).round() as usize + 1;
    (0..n).map(|i| GRID_LO + i as f64 * GRID_STEP).collect()
}

fn interp(xs: &[f64], ys: &[f64], target: f64) -> Option<f64> {
    if xs.len() < 2 || xs.len() != ys.len() {
        return None;
    }
    if target < xs[0] || target > xs[xs.len() - 1] {
        return None;
    }
    let idx = xs.partition_point(|&v| v < target);
    if idx == 0 {
        return Some(ys[0]);
    }
    if idx >= xs.len() {
        return Some(ys[ys.len() - 1]);
    }
    let x0 = xs[idx - 1];
    let x1 = xs[idx];
    let t = ((target - x0) / (x1 - x0).max(1e-12)).clamp(0.0, 1.0);
    Some(ys[idx - 1] + t * (ys[idx] - ys[idx - 1]))
}

/// Absorbance → rubberband baseline → uniform grid → [0, 1] min-max normalisation.
fn prepare(spec: &Spectrum, cfg: ProcessConfig) -> Option<Vec<f64>> {
    if spec.kind != SpectrumKind::Ftir || spec.x.len() < 2 {
        return None;
    }
    let abs_cfg = ProcessConfig {
        ftir_absorbance: true,
        baseline: false,
        smooth: cfg.smooth,
        normalize: false,
    };
    let (xs, raw, _) = process_trace(spec, abs_cfg);
    if xs.len() < 2 {
        return None;
    }
    let corrected = rubberband_baseline(&xs, &raw, false);

    let g = grid();
    let mut out = vec![f64::NAN; g.len()];
    let mut n_valid = 0usize;
    for (i, &gx) in g.iter().enumerate() {
        if let Some(v) = interp(&xs, &corrected, gx) {
            if v.is_finite() {
                out[i] = v;
                n_valid += 1;
            }
        }
    }
    if n_valid < MIN_POINTS {
        return None;
    }

    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in &out {
        if v.is_finite() {
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    let span = (hi - lo).max(1e-12);
    for v in &mut out {
        if v.is_finite() {
            *v = (*v - lo) / span;
        }
    }

    // Mild 44 cm⁻¹ boxcar smoothing suppresses interpolation / peak-width
    // mismatch between experimental and synthesised references.
    let r = 4usize;
    let smoothed: Vec<f64> = (0..out.len())
        .map(|i| {
            if out[i].is_nan() {
                return f64::NAN;
            }
            let lo_i = i.saturating_sub(r);
            let hi_i = (i + r).min(out.len() - 1);
            let mut sum = 0.0;
            let mut n = 0usize;
            for &v in &out[lo_i..=hi_i] {
                if v.is_finite() {
                    sum += v;
                    n += 1;
                }
            }
            if n > 0 {
                sum / n as f64
            } else {
                f64::NAN
            }
        })
        .collect();
    Some(smoothed)
}

fn pearson(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() != b.len() || a.len() < 8 {
        return None;
    }
    let n = a.len() as f64;
    let ma = a.iter().sum::<f64>() / n;
    let mb = b.iter().sum::<f64>() / n;
    let mut num = 0.0;
    let mut da = 0.0;
    let mut db = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        let dx = x - ma;
        let dy = y - mb;
        num += dx * dy;
        da += dx * dx;
        db += dy * dy;
    }
    let den = (da * db).sqrt();
    if den <= 1e-15 {
        return None;
    }
    Some(num / den)
}

fn derivative(v: &[f64]) -> Vec<f64> {
    let mut d = vec![f64::NAN; v.len()];
    for i in 1..v.len().saturating_sub(1) {
        let (a, b) = (v[i - 1], v[i + 1]);
        if a.is_finite() && b.is_finite() {
            d[i] = b - a;
        }
    }
    d
}

fn score_pair(query: &[f64], candidate: &[f64], g: &[f64]) -> Option<(f64, f64)> {
    let mut qa = Vec::with_capacity(g.len());
    let mut ca = Vec::with_capacity(g.len());
    for (i, &x) in g.iter().enumerate() {
        if x >= CO2_LO && x <= CO2_HI {
            continue;
        }
        if query[i].is_finite() && candidate[i].is_finite() {
            qa.push(query[i]);
            ca.push(candidate[i]);
        }
    }
    if qa.len() < MIN_POINTS {
        return None;
    }
    let r_abs = pearson(&qa, &ca)?.clamp(0.0, 1.0);

    let dq = derivative(query);
    let dc = derivative(candidate);
    let mut qd = Vec::with_capacity(g.len());
    let mut cd = Vec::with_capacity(g.len());
    for (i, &x) in g.iter().enumerate() {
        if x >= CO2_LO && x <= CO2_HI {
            continue;
        }
        if dq[i].is_finite() && dc[i].is_finite() {
            qd.push(dq[i]);
            cd.push(dc[i]);
        }
    }
    let r_deriv = if qd.len() >= 8 {
        pearson(&qd, &cd).unwrap_or(0.0).clamp(0.0, 1.0)
    } else {
        0.0
    };
    Some((r_abs, r_deriv))
}

/// FTIR reference files, reusing the process-wide library index so matching
/// and name-search share a single directory scan.
fn library_files() -> Vec<PathBuf> {
    crate::io::library_index()
        .iter()
        .filter(|e| e.kind == SpectrumKind::Ftir)
        .map(|e| e.path.clone())
        .collect()
}

/// Rank the built-in polymer library plus `library/ftir/` references against `query`.
pub fn match_ftir_library(query: &Spectrum, cfg: ProcessConfig) -> Vec<LibraryMatch> {
    let Some(q) = prepare(query, cfg) else {
        return Vec::new();
    };
    let g = grid();
    let score_cfg = ProcessConfig {
        smooth: SmoothLevel::Off,
        ..cfg
    };

    let mut results: Vec<LibraryMatch> = Vec::new();

    for poly in POLYMERS_DB {
        let cand = synthesize_polymer_reference(poly, SpectrumKind::Ftir);
        let Some(c) = prepare(&cand, score_cfg) else {
            continue;
        };
        let Some((r_abs, r_deriv)) = score_pair(&q, &c, &g) else {
            continue;
        };
        results.push(LibraryMatch {
            name: poly.name.to_string(),
            short: poly.short.to_string(),
            source: "Built-in".to_string(),
            hqi_pct: (0.55 * r_abs + 0.45 * r_deriv) * 100.0,
            r_abs,
            r_deriv,
            is_builtin: true,
            file_path: None,
        });
    }

    for path in library_files() {
        let Ok(spec) = load_spectrum_file(&path, Some(SpectrumKind::Ftir)) else {
            continue;
        };
        let Some(c) = prepare(&spec, score_cfg) else {
            continue;
        };
        let Some((r_abs, r_deriv)) = score_pair(&q, &c, &g) else {
            continue;
        };
        let short = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("ref")
            .to_string();
        // Show the real origin (`library/ftir`, `library/nmr`, `samples`, …).
        let source = path
            .parent()
            .and_then(|p| p.file_name())
            .and_then(|s| s.to_str())
            .map(|s| format!("{s}/"))
            .unwrap_or_else(|| "library/".to_string());
        results.push(LibraryMatch {
            name: spec.name.clone(),
            short,
            source,
            hqi_pct: (0.55 * r_abs + 0.45 * r_deriv) * 100.0,
            r_abs,
            r_deriv,
            is_builtin: false,
            file_path: Some(path),
        });
    }

    results.sort_by(|a, b| b.hqi_pct.total_cmp(&a.hqi_pct));
    results.truncate(MAX_MATCHES);
    results
}

// ============================================================================
// UV-VIS LIBRARY MATCHING
// ============================================================================

fn grid_uvvis() -> Vec<f64> {
    let n = ((UVVIS_GRID_HI - UVVIS_GRID_LO) / UVVIS_GRID_STEP).round() as usize + 1;
    (0..n)
        .map(|i| UVVIS_GRID_LO + i as f64 * UVVIS_GRID_STEP)
        .collect()
}

/// Absorbance -> baseline -> 200..800 nm grid (2 nm) -> [0, 1] normalisation.
fn prepare_uvvis(spec: &Spectrum, cfg: ProcessConfig) -> Option<Vec<f64>> {
    if spec.kind != SpectrumKind::UvVis || spec.x.len() < 2 {
        return None;
    }
    let abs_cfg = ProcessConfig {
        ftir_absorbance: false,
        baseline: false,
        smooth: cfg.smooth,
        normalize: false,
    };
    let (xs, raw, _) = process_trace(spec, abs_cfg);
    if xs.len() < 2 {
        return None;
    }
    let corrected = rubberband_baseline(&xs, &raw, false);

    let g = grid_uvvis();
    let mut out = vec![f64::NAN; g.len()];
    let mut n_valid = 0usize;
    for (i, &gx) in g.iter().enumerate() {
        if let Some(v) = interp(&xs, &corrected, gx) {
            if v.is_finite() {
                out[i] = v;
                n_valid += 1;
            }
        }
    }
    if n_valid < UVVIS_MIN_POINTS {
        return None;
    }

    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    for &v in &out {
        if v.is_finite() {
            lo = lo.min(v);
            hi = hi.max(v);
        }
    }
    let span = (hi - lo).max(1e-12);
    for v in &mut out {
        if v.is_finite() {
            *v = (*v - lo) / span;
        }
    }

    // ~10 nm boxcar smoothing suppresses peak-width mismatch between
    // experimental spectra and synthesised / library references.
    let r = 5usize;
    let smoothed: Vec<f64> = (0..out.len())
        .map(|i| {
            if out[i].is_nan() {
                return f64::NAN;
            }
            let lo_i = i.saturating_sub(r);
            let hi_i = (i + r).min(out.len() - 1);
            let mut sum = 0.0;
            let mut n = 0usize;
            for &v in &out[lo_i..=hi_i] {
                if v.is_finite() {
                    sum += v;
                    n += 1;
                }
            }
            if n > 0 {
                sum / n as f64
            } else {
                f64::NAN
            }
        })
        .collect();
    Some(smoothed)
}

fn score_pair_uvvis(query: &[f64], candidate: &[f64]) -> Option<(f64, f64)> {
    let mut qa = Vec::with_capacity(query.len());
    let mut ca = Vec::with_capacity(candidate.len());
    for i in 0..query.len() {
        if query[i].is_finite() && candidate[i].is_finite() {
            qa.push(query[i]);
            ca.push(candidate[i]);
        }
    }
    if qa.len() < UVVIS_MIN_POINTS {
        return None;
    }
    let r_abs = pearson(&qa, &ca)?.clamp(0.0, 1.0);

    let dq = derivative(query);
    let dc = derivative(candidate);
    let mut qd = Vec::with_capacity(query.len());
    let mut cd = Vec::with_capacity(candidate.len());
    for i in 0..query.len() {
        if dq[i].is_finite() && dc[i].is_finite() {
            qd.push(dq[i]);
            cd.push(dc[i]);
        }
    }
    let r_deriv = if qd.len() >= 8 {
        pearson(&qd, &cd).unwrap_or(0.0).clamp(0.0, 1.0)
    } else {
        0.0
    };
    Some((r_abs, r_deriv))
}

/// UV-Vis reference files discovered in `library/uvvis/` (and samples).
fn library_files_uvvis() -> Vec<PathBuf> {
    crate::io::library_index()
        .iter()
        .filter(|e| e.kind == SpectrumKind::UvVis)
        .map(|e| e.path.clone())
        .collect()
}

/// Rank the built-in UV-Vis demos plus `library/uvvis/` references against `query`.
///
/// Score = `0.60 · r_abs + 0.40 · r_deriv` on a 200–800 nm grid (2 nm step),
/// requiring at least 30 overlapping finite points.
pub fn match_uvvis_library(query: &Spectrum, cfg: ProcessConfig) -> Vec<LibraryMatch> {
    let Some(q) = prepare_uvvis(query, cfg) else {
        return Vec::new();
    };
    let score_cfg = ProcessConfig {
        smooth: SmoothLevel::Off,
        ..cfg
    };

    let mut results: Vec<LibraryMatch> = Vec::new();

    for preset in DemoPreset::ALL {
        if preset.kind() != SpectrumKind::UvVis {
            continue;
        }
        for demo in generate_demo_preset(preset) {
            if demo.kind != SpectrumKind::UvVis {
                continue;
            }
            let Some(c) = prepare_uvvis(&demo, score_cfg) else {
                continue;
            };
            let Some((r_abs, r_deriv)) = score_pair_uvvis(&q, &c) else {
                continue;
            };
            results.push(LibraryMatch {
                name: demo.name.clone(),
                short: demo
                    .name
                    .split(['(', '—'])
                    .next()
                    .unwrap_or(&demo.name)
                    .trim()
                    .to_string(),
                source: "Built-in".to_string(),
                hqi_pct: (0.60 * r_abs + 0.40 * r_deriv) * 100.0,
                r_abs,
                r_deriv,
                is_builtin: true,
                file_path: None,
            });
        }
    }

    for path in library_files_uvvis() {
        let Ok(spec) = load_spectrum_file(&path, Some(SpectrumKind::UvVis)) else {
            continue;
        };
        let Some(c) = prepare_uvvis(&spec, score_cfg) else {
            continue;
        };
        let Some((r_abs, r_deriv)) = score_pair_uvvis(&q, &c) else {
            continue;
        };
        let short = path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("ref")
            .to_string();
        let source = path
            .parent()
            .and_then(|p| p.file_name())
            .and_then(|s| s.to_str())
            .map(|s| format!("{s}/"))
            .unwrap_or_else(|| "library/".to_string());
        results.push(LibraryMatch {
            name: spec.name.clone(),
            short,
            source,
            hqi_pct: (0.60 * r_abs + 0.40 * r_deriv) * 100.0,
            r_abs,
            r_deriv,
            is_builtin: false,
            file_path: Some(path),
        });
    }

    results.sort_by(|a, b| b.hqi_pct.total_cmp(&a.hqi_pct));
    results.truncate(MAX_MATCHES);
    results
}
