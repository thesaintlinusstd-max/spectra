//! Multi-component spectral mixture deconvolution — pure-Rust Lawson–Hanson
//! Non-Negative Least Squares (NNLS) with residual-screened shortlisting.
//!
//! This mirrors the workflow of Thermo OMNIC Specta / Wiley KnowItAll mixture
//! analysis: a measured blend is explained as a non-negative sum of 1–3
//! reference spectra, with a composite HQI that rewards both the band envelope
//! and its first derivative.

use std::fmt::Write as _;
use std::path::PathBuf;

use crate::data::{SolventId, Spectrum, SpectrumKind};
use crate::processing::ProcessConfig;

#[derive(Debug, Clone)]
pub struct MixtureComponent {
    pub name: String,
    pub short: String,
    pub source: String,
    pub weight_pct: f64,
    #[allow(dead_code)] // absorbance-space coefficient (diagnostics / export)
    pub coeff: f64,
    pub is_builtin: bool,
    pub file_path: Option<PathBuf>,
}

#[derive(Debug, Clone)]
pub struct MixtureMatch {
    pub components: Vec<MixtureComponent>, // 1, 2, or 3 components sorted by weight_pct desc
    pub composite_hqi_pct: f64,
    pub r_abs: f64,
    pub r_deriv: f64,
    pub improvement_pct: f64, // composite_hqi_pct - best_single_hqi_pct
}

fn solve_small_linear(a: &[Vec<f64>], b: &[f64]) -> Option<Vec<f64>> {
    let n = b.len();
    if n == 0 {
        return Some(Vec::new());
    }
    let mut mat = a.to_vec();
    let mut rhs = b.to_vec();
    for col in 0..n {
        let mut piv = col;
        let mut best = mat[col][col].abs();
        for r in (col + 1)..n {
            let v = mat[r][col].abs();
            if v > best {
                best = v;
                piv = r;
            }
        }
        if !best.is_finite() || best < 1e-15 {
            return None;
        }
        if piv != col {
            mat.swap(col, piv);
            rhs.swap(col, piv);
        }
        let diag = mat[col][col];
        for r in (col + 1)..n {
            let f = mat[r][col] / diag;
            mat[r][col] = 0.0;
            for c in (col + 1)..n {
                mat[r][c] -= f * mat[col][c];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = rhs[i];
        for c in (i + 1)..n {
            s -= mat[i][c] * x[c];
        }
        x[i] = s / mat[i][i];
    }
    x.iter().all(|v| v.is_finite()).then_some(x)
}

/// Exact Lawson-Hanson NNLS minimizing 0.5 c^T G c - rhs^T c subject to c >= 0,
/// where G = A^T A (symmetric PSD) and rhs = A^T y.
pub fn nnls_gram(g: &[Vec<f64>], rhs: &[f64]) -> Vec<f64> {
    let k = rhs.len();
    let mut c = vec![0.0_f64; k];
    let mut passive = vec![false; k];

    for _ in 0..(k * 10 + 10) {
        // Negative gradient w = rhs - G c
        let mut best_j: Option<usize> = None;
        let mut best_grad = 1e-11_f64;
        for j in 0..k {
            if passive[j] {
                continue;
            }
            let mut gc = 0.0;
            for m in 0..k {
                gc += g[j][m] * c[m];
            }
            let grad = rhs[j] - gc;
            if grad > best_grad {
                best_grad = grad;
                best_j = Some(j);
            }
        }
        let Some(enter) = best_j else {
            break;
        };
        passive[enter] = true;

        // Inner feasibility loop
        for _ in 0..(k + 3) {
            let active: Vec<usize> = (0..k).filter(|&i| passive[i]).collect();
            if active.is_empty() {
                break;
            }
            let m = active.len();
            let mut gp = vec![vec![0.0; m]; m];
            let mut rp = vec![0.0; m];
            for (ri, &i) in active.iter().enumerate() {
                rp[ri] = rhs[i];
                for (ci, &j) in active.iter().enumerate() {
                    gp[ri][ci] = g[i][j] + if ri == ci { 1e-11 } else { 0.0 };
                }
            }
            let Some(sp) = solve_small_linear(&gp, &rp) else {
                passive[enter] = false;
                break;
            };
            if sp.iter().all(|&v| v > 1e-12) {
                for v in &mut c {
                    *v = 0.0;
                }
                for (&idx, &val) in active.iter().zip(&sp) {
                    c[idx] = val;
                }
                break;
            }
            let mut alpha = 1.0_f64;
            for (&idx, &sv) in active.iter().zip(&sp) {
                if sv <= 1e-12 {
                    let denom = c[idx] - sv;
                    if denom > 1e-15 {
                        alpha = alpha.min(c[idx] / denom);
                    }
                }
            }
            let mut s_full = vec![0.0; k];
            for (&idx, &val) in active.iter().zip(&sp) {
                s_full[idx] = val;
            }
            for j in 0..k {
                c[j] += alpha * (s_full[j] - c[j]);
                if passive[j] && c[j] <= 1e-12 {
                    c[j] = 0.0;
                    passive[j] = false;
                }
            }
        }
    }
    c
}

#[inline]
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

#[inline]
fn l2(a: &[f64]) -> f64 {
    a.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn normalize_unit(v: &mut [f64]) -> bool {
    let n = l2(v);
    if n <= 1e-12 {
        return false;
    }
    let inv = 1.0 / n;
    for x in v.iter_mut() {
        *x *= inv;
    }
    true
}

fn pearson_local(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len().min(b.len());
    if n < 8 {
        return 0.0;
    }
    let nf = n as f64;
    let ma = a[..n].iter().sum::<f64>() / nf;
    let mb = b[..n].iter().sum::<f64>() / nf;
    let mut num = 0.0;
    let mut da = 0.0;
    let mut db = 0.0;
    for (&x, &y) in a[..n].iter().zip(&b[..n]) {
        let dx = x - ma;
        let dy = y - mb;
        num += dx * dy;
        da += dx * dx;
        db += dy * dy;
    }
    let den = (da * db).sqrt();
    if den <= 1e-15 {
        0.0
    } else {
        num / den
    }
}

fn component_key(m: &MixtureMatch) -> String {
    let mut names: Vec<String> = m
        .components
        .iter()
        .map(|c| c.short.trim().to_ascii_lowercase())
        .collect();
    names.sort();
    names.join("|")
}

struct ProjectedRef {
    ref_idx: usize,
    norm: f64,
    ay: Vec<f64>,
    ad: Vec<f64>,
}

/// Analyze `query` as a non-negative mixture of 1–`max_components` references.
pub fn analyze_mixture(
    query: &Spectrum,
    cfg: ProcessConfig,
    solvent: SolventId,
    max_components: usize,
) -> Vec<MixtureMatch> {
    let Some(q) = crate::matcher::prepare_query_bundle(query, cfg, solvent) else {
        return Vec::new();
    };
    let refs = crate::matcher::collect_prepared_refs(query.kind, cfg);
    if refs.is_empty() {
        return Vec::new();
    }

    let active_idx: Vec<usize> = (0..q.vals.len())
        .filter(|&i| {
            !q.masked.get(i).copied().unwrap_or(true) && q.vals[i].is_finite()
        })
        .collect();
    if active_idx.len() < q.min_points {
        return Vec::new();
    }
    let n = active_idx.len();

    let mut qy = vec![0.0; n];
    for (t, &i) in active_idx.iter().enumerate() {
        qy[t] = q.vals[i];
    }
    if !normalize_unit(&mut qy) {
        return Vec::new();
    }
    let mut qd = vec![0.0; n];
    for (t, &i) in active_idx.iter().enumerate() {
        let d = q.deriv.get(i).copied().unwrap_or(0.0);
        qd[t] = if d.is_finite() { d } else { 0.0 };
    }
    normalize_unit(&mut qd);

    // Project each reference onto the active grid and normalize it to unit L2.
    // `ad` is scaled by exactly the same factor so a linear combination of `ay`
    // has a derivative equal to the same linear combination of `ad`.
    let mut projected: Vec<ProjectedRef> = Vec::new();
    for (ri, r) in refs.iter().enumerate() {
        let finite = active_idx
            .iter()
            .filter(|&&i| r.vals.get(i).map(|v| v.is_finite()).unwrap_or(false))
            .count();
        if finite * 4 < n * 3 {
            continue;
        }
        let mut ay = vec![0.0; n];
        for (t, &i) in active_idx.iter().enumerate() {
            if let Some(&v) = r.vals.get(i) {
                if v.is_finite() {
                    ay[t] = v as f64;
                }
            }
        }
        let norm = l2(&ay);
        if norm <= 1e-12 {
            continue;
        }
        let inv = 1.0 / norm;
        for v in &mut ay {
            *v *= inv;
        }
        let mut ad = vec![0.0; n];
        for (t, &i) in active_idx.iter().enumerate() {
            if let Some(&v) = r.deriv.get(i) {
                if v.is_finite() {
                    ad[t] = v as f64 * inv;
                }
            }
        }
        projected.push(ProjectedRef {
            ref_idx: ri,
            norm,
            ay,
            ad,
        });
    }
    if projected.is_empty() {
        return Vec::new();
    }

    // ------------------------------------------------------------------
    // Stage 1: residual-screened shortlist pool (<= 32 candidates)
    // ------------------------------------------------------------------
    let mut scores: Vec<(usize, f64)> = projected
        .iter()
        .enumerate()
        .map(|(p, r)| (p, dot(&qy, &r.ay)))
        .collect();
    scores.sort_by(|a, b| b.1.total_cmp(&a.1));

    let mut pool: Vec<usize> = Vec::new();
    let mut seen = vec![false; projected.len()];
    for &(p, _) in scores.iter().take(18) {
        if !seen[p] {
            seen[p] = true;
            pool.push(p);
        }
    }
    for &(top, _) in scores.iter().take(3) {
        let ay = &projected[top].ay;
        let alpha = dot(&qy, ay).max(0.0);
        let res: Vec<f64> = qy
            .iter()
            .zip(ay)
            .map(|(a, b)| (a - alpha * b).max(0.0))
            .collect();
        let mut rscores: Vec<(usize, f64)> = projected
            .iter()
            .enumerate()
            .filter(|(p, _)| !seen[*p])
            .map(|(p, r)| (p, dot(&res, &r.ay)))
            .collect();
        rscores.sort_by(|a, b| b.1.total_cmp(&a.1));
        for &(p, _) in rscores.iter().take(8) {
            if !seen[p] {
                seen[p] = true;
                pool.push(p);
            }
        }
        if pool.len() >= 32 {
            break;
        }
    }
    pool.truncate(32);
    let pool_n = pool.len();
    if pool_n == 0 {
        return Vec::new();
    }

    // ------------------------------------------------------------------
    // Stage 2: Gram matrix precomputation + exhaustive subset evaluation
    // ------------------------------------------------------------------
    let mut g = vec![vec![0.0; pool_n]; pool_n];
    let mut rhs = vec![0.0; pool_n];
    for a in 0..pool_n {
        rhs[a] = dot(&qy, &projected[pool[a]].ay);
        for b in 0..pool_n {
            g[a][b] = dot(&projected[pool[a]].ay, &projected[pool[b]].ay);
        }
    }

    let eval_subset = |subset: &[usize]| -> Option<MixtureMatch> {
        let k = subset.len();
        if k == 0 {
            return None;
        }
        let mut gg = vec![vec![0.0; k]; k];
        let mut rr = vec![0.0; k];
        for (a, &pa) in subset.iter().enumerate() {
            rr[a] = rhs[pa];
            for (b, &pb) in subset.iter().enumerate() {
                gg[a][b] = g[pa][pb];
            }
        }
        let c = nnls_gram(&gg, &rr);
        if c.iter().any(|&v| v <= 1e-9) {
            return None;
        }
        let eff: Vec<f64> = subset
            .iter()
            .enumerate()
            .map(|(a, &pa)| c[a] / projected[pool[pa]].norm.max(1e-12))
            .collect();
        let sum: f64 = eff.iter().sum();
        if sum <= 1e-12 {
            return None;
        }
        let weights: Vec<f64> = eff.iter().map(|v| v / sum * 100.0).collect();
        if weights.iter().any(|&w| w < 5.0) {
            return None;
        }

        let mut y_hat = vec![0.0; n];
        let mut d_hat = vec![0.0; n];
        for (a, &pa) in subset.iter().enumerate() {
            let pr = &projected[pool[pa]];
            for t in 0..n {
                y_hat[t] += c[a] * pr.ay[t];
                d_hat[t] += c[a] * pr.ad[t];
            }
        }
        let r_abs = pearson_local(&qy, &y_hat).clamp(0.0, 1.0);
        let r_deriv = pearson_local(&qd, &d_hat).clamp(0.0, 1.0);
        let composite_hqi_pct = (0.55 * r_abs + 0.45 * r_deriv) * 100.0;

        let mut components: Vec<MixtureComponent> = subset
            .iter()
            .enumerate()
            .map(|(a, &pa)| {
                let r = &refs[projected[pool[pa]].ref_idx];
                MixtureComponent {
                    name: r.name.clone(),
                    short: r.short.clone(),
                    source: r.source.clone(),
                    weight_pct: weights[a],
                    coeff: eff[a],
                    is_builtin: r.is_builtin,
                    file_path: r.file_path.clone(),
                }
            })
            .collect();
        components.sort_by(|a, b| b.weight_pct.total_cmp(&a.weight_pct));

        Some(MixtureMatch {
            components,
            composite_hqi_pct,
            r_abs,
            r_deriv,
            improvement_pct: 0.0,
        })
    };

    // 1-component solutions
    let mut singles: Vec<MixtureMatch> = Vec::new();
    let mut best_single_hqi = f64::NEG_INFINITY;
    for p in 0..pool_n {
        if let Some(m) = eval_subset(&[p]) {
            best_single_hqi = best_single_hqi.max(m.composite_hqi_pct);
            singles.push(m);
        }
    }
    if !best_single_hqi.is_finite() {
        best_single_hqi = 0.0;
    }
    singles.sort_by(|a, b| b.composite_hqi_pct.total_cmp(&a.composite_hqi_pct));

    // 2-component solutions (exhaustive within the shortlist)
    let mut pair_evaluated: Vec<((usize, usize), MixtureMatch)> = Vec::new();
    let mut best_pair_hqi = f64::NEG_INFINITY;
    if max_components >= 2 {
        for i in 0..pool_n {
            let si = refs[projected[pool[i]].ref_idx]
                .short
                .trim()
                .to_ascii_lowercase();
            for j in (i + 1)..pool_n {
                let sj = refs[projected[pool[j]].ref_idx]
                    .short
                    .trim()
                    .to_ascii_lowercase();
                if si == sj {
                    continue;
                }
                if let Some(m) = eval_subset(&[i, j]) {
                    best_pair_hqi = best_pair_hqi.max(m.composite_hqi_pct);
                    pair_evaluated.push(((i, j), m));
                }
            }
        }
    }
    pair_evaluated.sort_by(|a, b| b.1.composite_hqi_pct.total_cmp(&a.1.composite_hqi_pct));

    let mut top_pairs: Vec<((usize, usize), MixtureMatch)> = Vec::new();
    let mut pair_keys: Vec<String> = Vec::new();
    for (idx, m) in pair_evaluated {
        let key = component_key(&m);
        if pair_keys.iter().any(|k| k == &key) {
            continue;
        }
        pair_keys.push(key);
        top_pairs.push((idx, m));
        if top_pairs.len() >= 15 {
            break;
        }
    }
    if !best_pair_hqi.is_finite() {
        best_pair_hqi = 0.0;
    }

    // 3-component solutions: add every possible third reference to the best pairs
    let mut triples: Vec<MixtureMatch> = Vec::new();
    if max_components >= 3 {
        for ((i, j), _) in top_pairs.iter().take(12) {
            let si = refs[projected[pool[*i]].ref_idx]
                .short
                .trim()
                .to_ascii_lowercase();
            let sj = refs[projected[pool[*j]].ref_idx]
                .short
                .trim()
                .to_ascii_lowercase();
            for k in 0..pool_n {
                if k == *i || k == *j {
                    continue;
                }
                let sk = refs[projected[pool[k]].ref_idx]
                    .short
                    .trim()
                    .to_ascii_lowercase();
                if sk == si || sk == sj {
                    continue;
                }
                if let Some(m) = eval_subset(&[*i, *j, k]) {
                    triples.push(m);
                }
            }
        }
    }
    triples.sort_by(|a, b| b.composite_hqi_pct.total_cmp(&a.composite_hqi_pct));

    // ------------------------------------------------------------------
    // Anti-overfitting gate: a pure spectrum must keep its single match #1
    // ------------------------------------------------------------------
    let mut accepted: Vec<MixtureMatch> = Vec::new();
    for m in singles {
        accepted.push(m);
    }
    for (_, mut m) in top_pairs {
        if m.composite_hqi_pct >= best_single_hqi + 1.5 {
            m.improvement_pct = m.composite_hqi_pct - best_single_hqi;
            accepted.push(m);
        }
    }
    for mut m in triples {
        if m.composite_hqi_pct >= best_pair_hqi + 1.2
            && m.composite_hqi_pct >= best_single_hqi + 2.5
        {
            m.improvement_pct = m.composite_hqi_pct - best_single_hqi;
            accepted.push(m);
        }
    }

    accepted.sort_by(|a, b| b.composite_hqi_pct.total_cmp(&a.composite_hqi_pct));
    let mut out: Vec<MixtureMatch> = Vec::new();
    let mut keys: Vec<String> = Vec::new();
    for m in accepted {
        let key = component_key(&m);
        if keys.iter().any(|k| k == &key) {
            continue;
        }
        keys.push(key);
        out.push(m);
        if out.len() >= 15 {
            break;
        }
    }
    out
}

/// Markdown report suitable for clipboard export with `[c]`.
pub fn format_mixture_markdown(results: &[MixtureMatch], kind: SpectrumKind) -> String {
    let mut out = String::new();
    let _ = writeln!(
        out,
        "# SPECTRA Multi-Component Mixture Analysis — {}\n",
        kind.display_name()
    );
    let _ = writeln!(
        out,
        "| Rank | Composite HQI | ΔHQI | r_abs | r_deriv | Components |"
    );
    let _ = writeln!(
        out,
        "|-----:|--------------:|-----:|------:|--------:|:-----------|"
    );
    for (i, m) in results.iter().enumerate() {
        let comps = m
            .components
            .iter()
            .map(|c| format!("{:.1}% {}", c.weight_pct, c.short))
            .collect::<Vec<_>>()
            .join(" + ");
        let _ = writeln!(
            out,
            "| {} | {:.2}% | {:+.2}% | {:.4} | {:.4} | {} |",
            i + 1,
            m.composite_hqi_pct,
            m.improvement_pct,
            m.r_abs,
            m.r_deriv,
            comps
        );
    }
    let _ = writeln!(out);
    let _ = writeln!(
        out,
        "Composite HQI = 0.55 × absorbance correlation + 0.45 × derivative correlation."
    );
    out
}
