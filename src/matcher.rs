//! Spectral library matching engine — 100% pure Rust, zero extra crates.
//!
//! Compares an unknown FTIR, UV-Vis, ¹H NMR, or ¹³C NMR spectrum against the
//! built-in reference database and any user-supplied reference spectra in
//! `library/`, using an academic Hit Quality Index (HQI) with AsLS baseline
//! correction, Savitzky-Golay smoothing & 1st-derivative correlation, coverage
//! weighting, and an in-memory + binary disk fingerprint cache.

use std::collections::hash_map::DefaultHasher;
use std::fs;
use std::hash::{Hash, Hasher};
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{OnceLock, RwLock};
use std::time::UNIX_EPOCH;

use crate::data::{
    generate_demo_preset, synthesize_polymer_reference, DemoPreset, SolventId, Spectrum,
    SpectrumKind, POLYMERS_DB,
};
use crate::io::load_spectrum_file;
use crate::processing::{
    asls_baseline, process_trace, savitzky_golay, savitzky_golay_deriv1, BaselineMode,
    ProcessConfig, SmoothLevel, YAxisPreference,
};

// FTIR grid: 520..3760 cm⁻¹ at 4 cm⁻¹ step (811 pts)
const GRID_LO: f64 = 520.0;
const GRID_HI: f64 = 3760.0;
const GRID_STEP: f64 = 4.0;
const CO2_LO: f64 = 2280.0;
const CO2_HI: f64 = 2410.0;
const FP_HI: f64 = 1850.0;
const MIN_POINTS: usize = 120;

// UV-Vis grid: 200..800 nm at 2 nm step (301 pts)
const UVVIS_GRID_LO: f64 = 200.0;
const UVVIS_GRID_HI: f64 = 800.0;
const UVVIS_GRID_STEP: f64 = 2.0;
const UVVIS_MIN_POINTS: usize = 30;

// ¹H NMR grid: 0..12 ppm at 0.02 ppm step (601 pts)
const HNMR_GRID_LO: f64 = 0.0;
const HNMR_GRID_HI: f64 = 12.0;
const HNMR_GRID_STEP: f64 = 0.02;

// ¹³C NMR grid: 0..220 ppm at 0.5 ppm step (441 pts)
const CNMR_GRID_LO: f64 = 0.0;
const CNMR_GRID_HI: f64 = 220.0;
const CNMR_GRID_STEP: f64 = 0.5;

const NMR_MIN_POINTS: usize = 60;
const MAX_MATCHES: usize = 25;

const CACHE_MAGIC: &[u8; 8] = b"SPCQ0600";

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

// ============================================================================
// GRID & SIGNAL HELPERS
// ============================================================================

fn build_grid(lo: f64, hi: f64, step: f64) -> Vec<f64> {
    let n = ((hi - lo) / step).round() as usize + 1;
    (0..n).map(|i| lo + i as f64 * step).collect()
}

fn grid_ftir() -> Vec<f64> {
    build_grid(GRID_LO, GRID_HI, GRID_STEP)
}

fn grid_uvvis() -> Vec<f64> {
    build_grid(UVVIS_GRID_LO, UVVIS_GRID_HI, UVVIS_GRID_STEP)
}

fn grid_hnmr() -> Vec<f64> {
    build_grid(HNMR_GRID_LO, HNMR_GRID_HI, HNMR_GRID_STEP)
}

fn grid_cnmr() -> Vec<f64> {
    build_grid(CNMR_GRID_LO, CNMR_GRID_HI, CNMR_GRID_STEP)
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
        return Some(ys[xs.len() - 1]);
    }
    let x0 = xs[idx - 1];
    let x1 = xs[idx];
    let t = ((target - x0) / (x1 - x0).max(1e-12)).clamp(0.0, 1.0);
    Some(ys[idx - 1] + t * (ys[idx] - ys[idx - 1]))
}

/// Smooth & differentiate ONLY the contiguous finite slice `[first..=last]`,
/// leaving `f64::NAN` untouched outside the spectrum's measured bounds.
fn smooth_and_deriv_finite(vals: &[f64], smooth: SmoothLevel) -> (Vec<f64>, Vec<f64>) {
    let n = vals.len();
    let mut sm = vec![f64::NAN; n];
    let mut dv = vec![f64::NAN; n];

    let Some(first) = vals.iter().position(|v| v.is_finite()) else {
        return (sm, dv);
    };
    let Some(last) = vals.iter().rposition(|v| v.is_finite()) else {
        return (sm, dv);
    };
    if last <= first + 6 {
        return (sm, dv);
    }

    // Fill any isolated interior NaN linearly before SG filtering
    let mut seg: Vec<f64> = vals[first..=last].to_vec();
    let mut last_good = 0usize;
    for i in 0..seg.len() {
        if seg[i].is_finite() {
            if i > last_good + 1 {
                let y0 = seg[last_good];
                let y1 = seg[i];
                let span = (i - last_good) as f64;
                for (k, slot) in seg
                    .iter_mut()
                    .enumerate()
                    .take(i)
                    .skip(last_good + 1)
                {
                    *slot = y0 + (y1 - y0) * ((k - last_good) as f64 / span);
                }
            }
            last_good = i;
        }
    }

    let smoothed = savitzky_golay(&seg, smooth);
    let deriv = savitzky_golay_deriv1(&smoothed);
    for (offset, (&s, &d)) in smoothed.iter().zip(&deriv).enumerate() {
        sm[first + offset] = s;
        dv[first + offset] = d;
    }
    (sm, dv)
}

fn interpolate_and_normalize(
    xs: &[f64],
    corrected: &[f64],
    g: &[f64],
    min_points: usize,
) -> Option<(Vec<f64>, Vec<f64>)> {
    let mut out = vec![f64::NAN; g.len()];
    let mut n_valid = 0usize;
    for (i, &gx) in g.iter().enumerate() {
        if let Some(v) = interp(xs, corrected, gx) {
            if v.is_finite() {
                out[i] = v;
                n_valid += 1;
            }
        }
    }
    if n_valid < min_points {
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

    Some(smooth_and_deriv_finite(&out, SmoothLevel::Win5))
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

fn to_f32_vec(v: &[f64]) -> Vec<f32> {
    v.iter().map(|&x| x as f32).collect()
}

// ============================================================================
// SPECTRUM PREPARATION
// ============================================================================

/// Absorbance -> AsLS baseline -> 520..3760 cm⁻¹ grid (4 cm⁻¹) -> [0, 1] -> SG-5 + SG 1st deriv.
fn prepare_ftir(spec: &Spectrum, cfg: ProcessConfig) -> Option<(Vec<f64>, Vec<f64>)> {
    if spec.kind != SpectrumKind::Ftir || spec.x.len() < 2 {
        return None;
    }
    let abs_cfg = ProcessConfig {
        y_pref: YAxisPreference::Swapped,
        baseline: BaselineMode::Off,
        smooth: cfg.smooth,
        normalize: false,
    };
    let (xs, raw, _) = process_trace(spec, abs_cfg);
    if xs.len() < 2 {
        return None;
    }
    let corrected = asls_baseline(&raw, 1e6, 0.01, 8);
    let g = grid_ftir();
    interpolate_and_normalize(&xs, &corrected, &g, MIN_POINTS)
}

/// Absorbance -> AsLS baseline -> 200..800 nm grid (2 nm) -> [0, 1] -> SG-5 + SG 1st deriv.
fn prepare_uvvis(spec: &Spectrum, cfg: ProcessConfig) -> Option<(Vec<f64>, Vec<f64>)> {
    if spec.kind != SpectrumKind::UvVis || spec.x.len() < 2 {
        return None;
    }
    let abs_cfg = ProcessConfig {
        y_pref: YAxisPreference::Native,
        baseline: BaselineMode::Off,
        smooth: cfg.smooth,
        normalize: false,
    };
    let (xs, raw, _) = process_trace(spec, abs_cfg);
    if xs.len() < 2 {
        return None;
    }
    let corrected = asls_baseline(&raw, 5e5, 0.01, 8);
    let g = grid_uvvis();
    interpolate_and_normalize(&xs, &corrected, &g, UVVIS_MIN_POINTS)
}

/// Intensity -> AsLS baseline -> NMR ppm grid -> [0, 1] -> SG-5 + SG 1st deriv.
fn prepare_nmr(spec: &Spectrum, cfg: ProcessConfig, g: &[f64]) -> Option<(Vec<f64>, Vec<f64>)> {
    if !matches!(spec.kind, SpectrumKind::Hnmr | SpectrumKind::Cnmr) || spec.x.len() < 2 {
        return None;
    }
    let nmr_cfg = ProcessConfig {
        y_pref: YAxisPreference::Native,
        baseline: BaselineMode::Off,
        smooth: cfg.smooth,
        normalize: false,
    };
    let (xs, raw, _) = process_trace(spec, nmr_cfg);
    if xs.len() < 2 {
        return None;
    }
    let corrected = asls_baseline(&raw, 5e5, 0.005, 8);
    interpolate_and_normalize(&xs, &corrected, g, NMR_MIN_POINTS)
}

fn prepare_for_kind(spec: &Spectrum, cfg: ProcessConfig) -> Option<(Vec<f64>, Vec<f64>)> {
    match spec.kind {
        SpectrumKind::Ftir => prepare_ftir(spec, cfg),
        SpectrumKind::UvVis => prepare_uvvis(spec, cfg),
        SpectrumKind::Hnmr => prepare_nmr(spec, cfg, &grid_hnmr()),
        SpectrumKind::Cnmr => prepare_nmr(spec, cfg, &grid_cnmr()),
    }
}

// ============================================================================
// IN-MEMORY + BINARY DISK FINGERPRINT CACHE
// ============================================================================

#[derive(Clone)]
struct CachedFingerprint {
    name: String,
    short: String,
    source: String,
    path: PathBuf,
    vals: Vec<f32>,
    deriv: Vec<f32>,
}

struct KindCache {
    signature: u64,
    entries: Vec<CachedFingerprint>,
}

fn kind_slot(kind: SpectrumKind) -> usize {
    match kind {
        SpectrumKind::Ftir => 0,
        SpectrumKind::UvVis => 1,
        SpectrumKind::Hnmr => 2,
        SpectrumKind::Cnmr => 3,
    }
}

fn kind_tag(kind: SpectrumKind) -> &'static str {
    match kind {
        SpectrumKind::Ftir => "ftir",
        SpectrumKind::UvVis => "uvvis",
        SpectrumKind::Hnmr => "hnmr",
        SpectrumKind::Cnmr => "cnmr",
    }
}

fn matcher_cache_lock() -> &'static RwLock<[Option<KindCache>; 4]> {
    static CACHE: OnceLock<RwLock<[Option<KindCache>; 4]>> = OnceLock::new();
    CACHE.get_or_init(|| RwLock::new([None, None, None, None]))
}

/// Clear the in-memory fingerprint cache (called when the library index is refreshed).
pub fn invalidate_matcher_cache() {
    let mut guard = matcher_cache_lock()
        .write()
        .unwrap_or_else(|e| e.into_inner());
    *guard = [None, None, None, None];
}

fn directory_signature(paths: &[PathBuf]) -> u64 {
    let mut hasher = DefaultHasher::new();
    paths.len().hash(&mut hasher);
    for p in paths {
        p.to_string_lossy().hash(&mut hasher);
        if let Ok(meta) = fs::metadata(p) {
            meta.len().hash(&mut hasher);
            let secs = meta
                .modified()
                .ok()
                .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
                .map(|d| d.as_secs())
                .unwrap_or(0);
            secs.hash(&mut hasher);
        } else {
            0u64.hash(&mut hasher);
            0u64.hash(&mut hasher);
        }
    }
    hasher.finish()
}

fn cache_file_path(kind: SpectrumKind) -> PathBuf {
    let filename = format!("spectra_hqi_v6_{}.bin", kind_tag(kind));
    for root in crate::io::library_roots() {
        // library_roots() yields .../library/ftir, .../library/nmr, .../samples, etc.
        if let Some(parent) = root.parent() {
            if parent.file_name().and_then(|s| s.to_str()) == Some("library") {
                return parent.join(".cache").join(filename);
            }
        }
    }
    std::env::temp_dir().join(filename)
}

fn write_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn write_u64(buf: &mut Vec<u8>, v: u64) {
    buf.extend_from_slice(&v.to_le_bytes());
}

fn write_str(buf: &mut Vec<u8>, s: &str) {
    let bytes = s.as_bytes();
    write_u32(buf, bytes.len() as u32);
    buf.extend_from_slice(bytes);
}

fn write_f32_slice(buf: &mut Vec<u8>, vals: &[f32]) {
    for &v in vals {
        buf.extend_from_slice(&v.to_le_bytes());
    }
}

fn read_u32(data: &[u8], pos: &mut usize) -> Option<u32> {
    let end = pos.checked_add(4)?;
    let slice = data.get(*pos..end)?;
    *pos = end;
    Some(u32::from_le_bytes([slice[0], slice[1], slice[2], slice[3]]))
}

fn read_u64(data: &[u8], pos: &mut usize) -> Option<u64> {
    let end = pos.checked_add(8)?;
    let slice = data.get(*pos..end)?;
    *pos = end;
    Some(u64::from_le_bytes([
        slice[0], slice[1], slice[2], slice[3], slice[4], slice[5], slice[6], slice[7],
    ]))
}

fn read_str(data: &[u8], pos: &mut usize) -> Option<String> {
    let len = read_u32(data, pos)? as usize;
    let end = pos.checked_add(len)?;
    let slice = data.get(*pos..end)?;
    *pos = end;
    String::from_utf8(slice.to_vec()).ok()
}

fn read_f32_vec(data: &[u8], pos: &mut usize, len: usize) -> Option<Vec<f32>> {
    let byte_len = len.checked_mul(4)?;
    let end = pos.checked_add(byte_len)?;
    let slice = data.get(*pos..end)?;
    *pos = end;
    let mut out = Vec::with_capacity(len);
    for chunk in slice.as_chunks::<4>().0 {
        out.push(f32::from_le_bytes(*chunk));
    }
    Some(out)
}

fn load_binary_cache(path: &Path, expected_sig: u64) -> Option<Vec<CachedFingerprint>> {
    let mut file = fs::File::open(path).ok()?;
    let mut data = Vec::new();
    file.read_to_end(&mut data).ok()?;

    if data.len() < 20 || &data[0..8] != CACHE_MAGIC {
        return None;
    }
    let mut pos = 8usize;
    let sig = read_u64(&data, &mut pos)?;
    if sig != expected_sig {
        return None;
    }
    let count = read_u32(&data, &mut pos)? as usize;
    let mut entries = Vec::with_capacity(count.min(16_384));

    for _ in 0..count {
        let path_str = read_str(&data, &mut pos)?;
        let name = read_str(&data, &mut pos)?;
        let short = read_str(&data, &mut pos)?;
        let source = read_str(&data, &mut pos)?;
        let vec_len = read_u32(&data, &mut pos)? as usize;
        let vals = read_f32_vec(&data, &mut pos, vec_len)?;
        let deriv = read_f32_vec(&data, &mut pos, vec_len)?;
        entries.push(CachedFingerprint {
            name,
            short,
            source,
            path: PathBuf::from(path_str),
            vals,
            deriv,
        });
    }

    Some(entries)
}

fn save_binary_cache(path: &Path, sig: u64, entries: &[CachedFingerprint]) {
    if let Some(parent) = path.parent() {
        if fs::create_dir_all(parent).is_err() {
            return;
        }
    }

    let mut buf = Vec::with_capacity(20 + entries.len() * 4096);
    buf.extend_from_slice(CACHE_MAGIC);
    write_u64(&mut buf, sig);
    write_u32(&mut buf, entries.len() as u32);

    for e in entries {
        write_str(&mut buf, &e.path.to_string_lossy());
        write_str(&mut buf, &e.name);
        write_str(&mut buf, &e.short);
        write_str(&mut buf, &e.source);
        let vec_len = e.vals.len().min(e.deriv.len());
        write_u32(&mut buf, vec_len as u32);
        write_f32_slice(&mut buf, &e.vals[..vec_len]);
        write_f32_slice(&mut buf, &e.deriv[..vec_len]);
    }

    let tmp_name = format!(
        "{}.part.{}",
        path.file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("cache.bin"),
        std::process::id()
    );
    let tmp_path = path.with_file_name(tmp_name);
    if let Ok(mut f) = fs::File::create(&tmp_path) {
        if f.write_all(&buf).is_ok() && f.flush().is_ok() {
            let _ = fs::rename(&tmp_path, path);
            return;
        }
    }
    let _ = fs::remove_file(&tmp_path);
}

fn cached_library_entries(kind: SpectrumKind) -> Vec<CachedFingerprint> {
    let mut paths: Vec<PathBuf> = crate::io::library_index()
        .iter()
        .filter(|e| e.kind == kind)
        .map(|e| e.path.clone())
        .collect();
    paths.sort();

    let sig = directory_signature(&paths);
    let slot = kind_slot(kind);

    // 1. Check in-memory RAM cache
    {
        let guard = matcher_cache_lock()
            .read()
            .unwrap_or_else(|e| e.into_inner());
        if let Some(kc) = &guard[slot] {
            if kc.signature == sig {
                return kc.entries.clone();
            }
        }
    }

    // 2. Check binary disk cache
    let disk_path = cache_file_path(kind);
    if let Some(entries) = load_binary_cache(&disk_path, sig) {
        let mut guard = matcher_cache_lock()
            .write()
            .unwrap_or_else(|e| e.into_inner());
        guard[slot] = Some(KindCache {
            signature: sig,
            entries: entries.clone(),
        });
        return entries;
    }

    // 3. Cold build: parse & fingerprint all library files
    let prep_cfg = ProcessConfig {
        y_pref: YAxisPreference::Native,
        baseline: BaselineMode::Off,
        smooth: SmoothLevel::Off,
        normalize: false,
    };

    let mut entries = Vec::with_capacity(paths.len());
    for path in paths {
        let Ok(spec) = load_spectrum_file(&path, Some(kind)) else {
            continue;
        };
        let Some((vals, deriv)) = prepare_for_kind(&spec, prep_cfg) else {
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

        entries.push(CachedFingerprint {
            name: spec.name,
            short,
            source,
            path,
            vals: to_f32_vec(&vals),
            deriv: to_f32_vec(&deriv),
        });
    }

    // Best-effort atomic disk cache write
    save_binary_cache(&disk_path, sig, &entries);

    let mut guard = matcher_cache_lock()
        .write()
        .unwrap_or_else(|e| e.into_inner());
    guard[slot] = Some(KindCache {
        signature: sig,
        entries: entries.clone(),
    });
    entries
}

// ============================================================================
// SCORING FUNCTIONS
// ============================================================================

fn score_ftir_pair(
    q_val: &[f64],
    q_der: &[f64],
    c_val: &[f32],
    c_der: &[f32],
    g: &[f64],
) -> Option<(f64, f64, f64)> {
    let n = g
        .len()
        .min(q_val.len())
        .min(q_der.len())
        .min(c_val.len())
        .min(c_der.len());

    let mut q_valid = 0usize;
    let mut qa = Vec::with_capacity(n);
    let mut ca = Vec::with_capacity(n);
    let mut qd = Vec::with_capacity(n);
    let mut cd = Vec::with_capacity(n);
    let mut qfp = Vec::with_capacity(n);
    let mut cfp = Vec::with_capacity(n);

    for i in 0..n {
        let x = g[i];
        if (CO2_LO..=CO2_HI).contains(&x) {
            continue;
        }
        let qv = q_val[i];
        if qv.is_finite() {
            q_valid += 1;
            let cv = c_val[i] as f64;
            if cv.is_finite() {
                qa.push(qv);
                ca.push(cv);
                if x <= FP_HI {
                    qfp.push(qv);
                    cfp.push(cv);
                }
            }
        }
        let qdv = q_der[i];
        let cdv = c_der[i] as f64;
        if qdv.is_finite() && cdv.is_finite() {
            qd.push(qdv);
            cd.push(cdv);
        }
    }

    if qa.len() < MIN_POINTS {
        return None;
    }

    let r_abs = pearson(&qa, &ca)?.clamp(0.0, 1.0);
    let r_deriv = if qd.len() >= 16 {
        pearson(&qd, &cd).unwrap_or(0.0).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let r_fp = if qfp.len() >= 60 {
        pearson(&qfp, &cfp).unwrap_or(r_abs).clamp(0.0, 1.0)
    } else {
        r_abs
    };

    let overlap_ratio = (qa.len() as f64 / q_valid.max(1) as f64).clamp(0.0, 1.0);
    let cov_factor = 0.55 + 0.45 * overlap_ratio;
    let hqi_pct = (0.35 * r_abs + 0.40 * r_deriv + 0.25 * r_fp) * cov_factor * 100.0;

    Some((hqi_pct, r_abs, r_deriv))
}

fn score_uvvis_pair(
    q_val: &[f64],
    q_der: &[f64],
    c_val: &[f32],
    c_der: &[f32],
) -> Option<(f64, f64, f64)> {
    let n = q_val
        .len()
        .min(q_der.len())
        .min(c_val.len())
        .min(c_der.len());

    let mut q_valid = 0usize;
    let mut qa = Vec::with_capacity(n);
    let mut ca = Vec::with_capacity(n);
    let mut qd = Vec::with_capacity(n);
    let mut cd = Vec::with_capacity(n);

    for i in 0..n {
        let qv = q_val[i];
        if qv.is_finite() {
            q_valid += 1;
            let cv = c_val[i] as f64;
            if cv.is_finite() {
                qa.push(qv);
                ca.push(cv);
            }
        }
        let qdv = q_der[i];
        let cdv = c_der[i] as f64;
        if qdv.is_finite() && cdv.is_finite() {
            qd.push(qdv);
            cd.push(cdv);
        }
    }

    if qa.len() < UVVIS_MIN_POINTS {
        return None;
    }

    let r_abs = pearson(&qa, &ca)?.clamp(0.0, 1.0);
    let r_deriv = if qd.len() >= 8 {
        pearson(&qd, &cd).unwrap_or(0.0).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let overlap_ratio = (qa.len() as f64 / q_valid.max(1) as f64).clamp(0.0, 1.0);
    let cov_factor = 0.60 + 0.40 * overlap_ratio;
    let hqi_pct = (0.60 * r_abs + 0.40 * r_deriv) * cov_factor * 100.0;

    Some((hqi_pct, r_abs, r_deriv))
}

fn score_nmr_pair(
    q_val: &[f64],
    q_der: &[f64],
    c_val: &[f32],
    c_der: &[f32],
    masked: &[bool],
) -> Option<(f64, f64, f64)> {
    let n = masked
        .len()
        .min(q_val.len())
        .min(q_der.len())
        .min(c_val.len())
        .min(c_der.len());

    let mut q_valid = 0usize;
    let mut qa = Vec::with_capacity(n);
    let mut ca = Vec::with_capacity(n);
    let mut qd = Vec::with_capacity(n);
    let mut cd = Vec::with_capacity(n);

    for i in 0..n {
        if masked[i] {
            continue;
        }
        let qv = q_val[i];
        if qv.is_finite() {
            q_valid += 1;
            let cv = c_val[i] as f64;
            if cv.is_finite() {
                qa.push(qv);
                ca.push(cv);
            }
        }
        let qdv = q_der[i];
        let cdv = c_der[i] as f64;
        if qdv.is_finite() && cdv.is_finite() {
            qd.push(qdv);
            cd.push(cdv);
        }
    }

    if qa.len() < NMR_MIN_POINTS {
        return None;
    }

    let r_abs = pearson(&qa, &ca)?.clamp(0.0, 1.0);
    let r_deriv = if qd.len() >= 16 {
        pearson(&qd, &cd).unwrap_or(0.0).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let overlap_ratio = (qa.len() as f64 / q_valid.max(1) as f64).clamp(0.0, 1.0);
    let cov_factor = 0.60 + 0.40 * overlap_ratio;
    let hqi_pct = (0.65 * r_abs + 0.35 * r_deriv) * cov_factor * 100.0;

    Some((hqi_pct, r_abs, r_deriv))
}

// ============================================================================
// PUBLIC MATCHING ENTRY POINTS
// ============================================================================

/// Rank the built-in polymer library plus `library/ftir/` references against `query`.
pub fn match_ftir_library(query: &Spectrum, cfg: ProcessConfig) -> Vec<LibraryMatch> {
    let Some((q_val, q_der)) = prepare_ftir(query, cfg) else {
        return Vec::new();
    };
    let g = grid_ftir();
    let score_cfg = ProcessConfig {
        smooth: SmoothLevel::Off,
        ..cfg
    };

    let mut results: Vec<LibraryMatch> = Vec::new();

    for poly in POLYMERS_DB {
        let cand = synthesize_polymer_reference(poly, SpectrumKind::Ftir);
        let Some((c_val, c_der)) = prepare_ftir(&cand, score_cfg) else {
            continue;
        };
        let cv32 = to_f32_vec(&c_val);
        let cd32 = to_f32_vec(&c_der);
        let Some((hqi_pct, r_abs, r_deriv)) = score_ftir_pair(&q_val, &q_der, &cv32, &cd32, &g)
        else {
            continue;
        };
        results.push(LibraryMatch {
            name: poly.name.to_string(),
            short: poly.short.to_string(),
            source: "Built-in".to_string(),
            hqi_pct,
            r_abs,
            r_deriv,
            is_builtin: true,
            file_path: None,
        });
    }

    for entry in cached_library_entries(SpectrumKind::Ftir) {
        let Some((hqi_pct, r_abs, r_deriv)) =
            score_ftir_pair(&q_val, &q_der, &entry.vals, &entry.deriv, &g)
        else {
            continue;
        };
        results.push(LibraryMatch {
            name: entry.name,
            short: entry.short,
            source: entry.source,
            hqi_pct,
            r_abs,
            r_deriv,
            is_builtin: false,
            file_path: Some(entry.path),
        });
    }

    results.sort_by(|a, b| b.hqi_pct.total_cmp(&a.hqi_pct));
    results.truncate(MAX_MATCHES);
    results
}

/// Rank the built-in UV-Vis demos plus `library/uvvis/` references against `query`.
pub fn match_uvvis_library(query: &Spectrum, cfg: ProcessConfig) -> Vec<LibraryMatch> {
    let Some((q_val, q_der)) = prepare_uvvis(query, cfg) else {
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
            let Some((c_val, c_der)) = prepare_uvvis(&demo, score_cfg) else {
                continue;
            };
            let cv32 = to_f32_vec(&c_val);
            let cd32 = to_f32_vec(&c_der);
            let Some((hqi_pct, r_abs, r_deriv)) =
                score_uvvis_pair(&q_val, &q_der, &cv32, &cd32)
            else {
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
                hqi_pct,
                r_abs,
                r_deriv,
                is_builtin: true,
                file_path: None,
            });
        }
    }

    for entry in cached_library_entries(SpectrumKind::UvVis) {
        let Some((hqi_pct, r_abs, r_deriv)) =
            score_uvvis_pair(&q_val, &q_der, &entry.vals, &entry.deriv)
        else {
            continue;
        };
        results.push(LibraryMatch {
            name: entry.name,
            short: entry.short,
            source: entry.source,
            hqi_pct,
            r_abs,
            r_deriv,
            is_builtin: false,
            file_path: Some(entry.path),
        });
    }

    results.sort_by(|a, b| b.hqi_pct.total_cmp(&a.hqi_pct));
    results.truncate(MAX_MATCHES);
    results
}

/// Rank the built-in polymer NMR references plus `library/nmr/` files against `query`,
/// actively masking residual solvent, water/HDO, and TMS signals.
pub fn match_nmr_library(
    query: &Spectrum,
    cfg: ProcessConfig,
    solvent: SolventId,
) -> Vec<LibraryMatch> {
    if !matches!(query.kind, SpectrumKind::Hnmr | SpectrumKind::Cnmr) {
        return Vec::new();
    }
    let g = match query.kind {
        SpectrumKind::Hnmr => grid_hnmr(),
        SpectrumKind::Cnmr => grid_cnmr(),
        _ => return Vec::new(),
    };

    let Some((q_val, q_der)) = prepare_nmr(query, cfg, &g) else {
        return Vec::new();
    };

    let masked: Vec<bool> = match query.kind {
        SpectrumKind::Hnmr => {
            let res_h = solvent.residual_h();
            let water = solvent.water_h();
            g.iter()
                .map(|&x| {
                    (x - 0.0).abs() <= 0.08
                        || (x - water).abs() <= 0.12
                        || res_h.iter().any(|&s| (x - s).abs() <= 0.12)
                })
                .collect()
        }
        SpectrumKind::Cnmr => {
            let res_c = solvent.residual_c();
            g.iter()
                .map(|&x| res_c.iter().any(|&s| (x - s).abs() <= 1.5))
                .collect()
        }
        _ => vec![false; g.len()],
    };

    let score_cfg = ProcessConfig {
        smooth: SmoothLevel::Off,
        ..cfg
    };

    let mut results: Vec<LibraryMatch> = Vec::new();

    for poly in POLYMERS_DB {
        let has_peaks = match query.kind {
            SpectrumKind::Hnmr => !poly.hnmr_peaks.is_empty(),
            SpectrumKind::Cnmr => !poly.cnmr_peaks.is_empty(),
            _ => false,
        };
        if !has_peaks {
            continue;
        }
        let cand = synthesize_polymer_reference(poly, query.kind);
        let Some((c_val, c_der)) = prepare_nmr(&cand, score_cfg, &g) else {
            continue;
        };
        let cv32 = to_f32_vec(&c_val);
        let cd32 = to_f32_vec(&c_der);
        let Some((hqi_pct, r_abs, r_deriv)) =
            score_nmr_pair(&q_val, &q_der, &cv32, &cd32, &masked)
        else {
            continue;
        };
        results.push(LibraryMatch {
            name: poly.name.to_string(),
            short: poly.short.to_string(),
            source: "Built-in".to_string(),
            hqi_pct,
            r_abs,
            r_deriv,
            is_builtin: true,
            file_path: None,
        });
    }

    for entry in cached_library_entries(query.kind) {
        let Some((hqi_pct, r_abs, r_deriv)) =
            score_nmr_pair(&q_val, &q_der, &entry.vals, &entry.deriv, &masked)
        else {
            continue;
        };
        results.push(LibraryMatch {
            name: entry.name,
            short: entry.short,
            source: entry.source,
            hqi_pct,
            r_abs,
            r_deriv,
            is_builtin: false,
            file_path: Some(entry.path),
        });
    }

    results.sort_by(|a, b| b.hqi_pct.total_cmp(&a.hqi_pct));
    results.truncate(MAX_MATCHES);
    results
}
