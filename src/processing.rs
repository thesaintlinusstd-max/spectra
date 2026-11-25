use crate::data::{best_peak_label_in_context, SolventId, Spectrum, SpectrumKind, YUnit};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BaselineMode {
    Off,
    Rubberband,
    Asls,
    Arpls,
}

impl BaselineMode {
    pub fn next(self) -> Self {
        match self {
            Self::Off => Self::Rubberband,
            Self::Rubberband => Self::Asls,
            Self::Asls => Self::Arpls,
            Self::Arpls => Self::Off,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Off => "OFF",
            Self::Rubberband => "Rubberband",
            Self::Asls => "AsLS",
            Self::Arpls => "arPLS",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NormMode {
    Off,
    MinMax,
    VectorSnv,
    Area,
}

impl NormMode {
    pub fn next(self) -> Self {
        match self {
            Self::Off => Self::MinMax,
            Self::MinMax => Self::VectorSnv,
            Self::VectorSnv => Self::Area,
            Self::Area => Self::Off,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Off => "OFF",
            Self::MinMax => "MinMax",
            Self::VectorSnv => "SNV",
            Self::Area => "Area",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum OpticsMode {
    Off,
    AtrDiamond,
    AtrZnSe,
    AtrGe,
    AtrSi,
    KubelkaMunk,
}

impl OpticsMode {
    pub fn next(self) -> Self {
        match self {
            Self::Off => Self::AtrDiamond,
            Self::AtrDiamond => Self::AtrZnSe,
            Self::AtrZnSe => Self::AtrGe,
            Self::AtrGe => Self::AtrSi,
            Self::AtrSi => Self::KubelkaMunk,
            Self::KubelkaMunk => Self::Off,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Off => "OFF",
            Self::AtrDiamond => "ATR-Dia",
            Self::AtrZnSe => "ATR-ZnSe",
            Self::AtrGe => "ATR-Ge",
            Self::AtrSi => "ATR-Si",
            Self::KubelkaMunk => "KM F(R)",
        }
    }

    pub fn atr_n1(self) -> Option<f64> {
        match self {
            Self::AtrDiamond => Some(2.400),
            Self::AtrZnSe => Some(2.403),
            Self::AtrGe => Some(4.003),
            Self::AtrSi => Some(3.422),
            Self::Off | Self::KubelkaMunk => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum DerivOrder {
    Raw,
    First,
    Second,
}

impl DerivOrder {
    pub fn next(self) -> Self {
        match self {
            Self::Raw => Self::First,
            Self::First => Self::Second,
            Self::Second => Self::Raw,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Raw => "Raw",
            Self::First => "dY/dX",
            Self::Second => "d²Y/dX²",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum YAxisPreference {
    /// Default convention per technique: FTIR -> %T, UV-Vis -> Absorbance, NMR -> Intensity
    Native,
    /// Flipped mode: FTIR -> Absorbance, UV-Vis -> %T
    Swapped,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct ProcessConfig {
    pub y_pref: YAxisPreference,
    pub baseline: BaselineMode,
    pub smooth: SmoothLevel,
    pub norm: NormMode,
    pub optics: OpticsMode,
    pub deriv: DerivOrder,
    pub atm_co2: bool,
}

impl Default for ProcessConfig {
    fn default() -> Self {
        Self {
            y_pref: YAxisPreference::Native,
            baseline: BaselineMode::Off,
            smooth: SmoothLevel::Off,
            norm: NormMode::Off,
            optics: OpticsMode::Off,
            deriv: DerivOrder::Raw,
            atm_co2: false,
        }
    }
}

impl ProcessConfig {
    pub fn ftir_absorbance() -> Self {
        Self {
            y_pref: YAxisPreference::Swapped,
            ..Self::default()
        }
    }
}

#[derive(Debug, Clone)]
pub struct Peak {
    pub x: f64,
    pub y: f64,
    pub prominence: f64,
    pub fwhm: Option<f64>,
    pub snr: f64,
    pub label: String,
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
    for (i, slot) in out.iter_mut().enumerate() {
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
        *slot = acc / norm;
    }
    out
}

/// 7-point quadratic Savitzky-Golay 1st derivative (Steinier 1972) with reflected edges.
pub fn savitzky_golay_deriv1(ys: &[f64]) -> Vec<f64> {
    const COEFFS: [f64; 7] = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
    const NORM: f64 = 28.0;
    let n = ys.len();
    if n < 7 {
        return vec![0.0; n];
    }
    let half = 3isize;
    let mut out = vec![0.0; n];
    for (i, slot) in out.iter_mut().enumerate() {
        let mut acc = 0.0;
        for (k, &w) in COEFFS.iter().enumerate() {
            let mut idx = i as isize + (k as isize - half);
            if idx < 0 {
                idx = -idx;
            }
            if idx >= n as isize {
                idx = 2 * (n as isize - 1) - idx;
            }
            let c = idx.clamp(0, (n - 1) as isize) as usize;
            acc += w * ys[c];
        }
        *slot = acc / NORM;
    }
    out
}

/// 9-point or 15-point quadratic/cubic Savitzky-Golay 2nd derivative with reflected edges.
pub fn savitzky_golay_deriv2(ys: &[f64], win15: bool) -> Vec<f64> {
    const COEFFS_9: [f64; 9] = [28.0, 7.0, -8.0, -17.0, -20.0, -17.0, -8.0, 7.0, 28.0];
    const NORM_9: f64 = 462.0;
    const COEFFS_15: [f64; 15] = [
        91.0, 52.0, 19.0, -8.0, -29.0, -44.0, -53.0, -56.0, -53.0, -44.0, -29.0, -8.0, 19.0, 52.0,
        91.0,
    ];
    const NORM_15: f64 = 6188.0;

    let n = ys.len();
    let (coeffs, norm): (&[f64], f64) = if win15 && n >= 15 {
        (&COEFFS_15, NORM_15)
    } else if n >= 9 {
        (&COEFFS_9, NORM_9)
    } else {
        return vec![0.0; n];
    };

    let half = (coeffs.len() / 2) as isize;
    let mut out = vec![0.0; n];
    for (i, slot) in out.iter_mut().enumerate() {
        let mut acc = 0.0;
        for (k, &w) in coeffs.iter().enumerate() {
            let mut idx = i as isize + (k as isize - half);
            if idx < 0 {
                idx = -idx;
            }
            if idx >= n as isize {
                idx = 2 * (n as isize - 1) - idx;
            }
            let c = idx.clamp(0, (n - 1) as isize) as usize;
            acc += w * ys[c];
        }
        *slot = acc / norm;
    }
    out
}

/// Physical Savitzky-Golay spectral derivative scaled by median |Δx|.
pub fn spectral_derivative(xs: &[f64], ys: &[f64], order: DerivOrder) -> Vec<f64> {
    if order == DerivOrder::Raw {
        return ys.to_vec();
    }
    let n = xs.len().min(ys.len());
    if n < 2 {
        return ys.to_vec();
    }
    let mut dxs: Vec<f64> = xs[..n]
        .windows(2)
        .map(|w| (w[1] - w[0]).abs())
        .filter(|&d| d > 1e-12)
        .collect();
    let dx = if dxs.is_empty() {
        1.0
    } else {
        dxs.sort_by(|a, b| a.total_cmp(b));
        dxs[dxs.len() / 2]
    };
    match order {
        DerivOrder::Raw => ys.to_vec(),
        DerivOrder::First => {
            let d1 = savitzky_golay_deriv1(&ys[..n]);
            if dx > 1e-9 {
                d1.into_iter().map(|v| v / dx).collect()
            } else {
                d1
            }
        }
        DerivOrder::Second => {
            let d2 = savitzky_golay_deriv2(&ys[..n], n > 600);
            if dx > 1e-9 {
                let dx2 = dx * dx;
                d2.into_iter().map(|v| v / dx2).collect()
            } else {
                d2
            }
        }
    }
}

/// O(N) symmetric pentadiagonal LDLᵀ solver for `(W + λ D₂ᵀD₂) z = W y`.
fn pentadiag_ldlt_solve(w: &[f64], ys: &[f64], lam: f64) -> Vec<f64> {
    let n = ys.len();
    if n < 5 {
        return ys.to_vec();
    }
    let mut p0 = vec![0.0_f64; n];
    let mut p1 = vec![0.0_f64; n.saturating_sub(1)];
    let mut p2 = vec![0.0_f64; n.saturating_sub(2)];
    for i in 0..(n - 2) {
        p0[i] += lam;
        p0[i + 1] += 4.0 * lam;
        p0[i + 2] += lam;
        p1[i] += -2.0 * lam;
        p1[i + 1] += -2.0 * lam;
        p2[i] += lam;
    }

    let mut z = vec![0.0_f64; n];
    let mut d = vec![0.0_f64; n];
    let mut l1 = vec![0.0_f64; n];
    let mut l2 = vec![0.0_f64; n];
    let mut v = vec![0.0_f64; n];

    for i in 0..n {
        let mut a0 = w[i] + p0[i];
        if i >= 1 {
            a0 -= l1[i - 1] * l1[i - 1] * d[i - 1];
        }
        if i >= 2 {
            a0 -= l2[i - 2] * l2[i - 2] * d[i - 2];
        }
        d[i] = a0.max(1e-18);

        if i + 1 < n {
            let mut a1 = p1[i];
            if i >= 1 {
                a1 -= l2[i - 1] * l1[i - 1] * d[i - 1];
            }
            l1[i] = a1 / d[i];
        }
        if i + 2 < n {
            l2[i] = p2[i] / d[i];
        }
    }
    for i in 0..n {
        let mut s = w[i] * ys[i];
        if i >= 1 {
            s -= l1[i - 1] * v[i - 1];
        }
        if i >= 2 {
            s -= l2[i - 2] * v[i - 2];
        }
        v[i] = s;
    }
    for i in (0..n).rev() {
        let mut s = v[i] / d[i];
        if i + 1 < n {
            s -= l1[i] * z[i + 1];
        }
        if i + 2 < n {
            s -= l2[i] * z[i + 2];
        }
        z[i] = s;
    }
    z
}

/// Eilers & Boelens (2005) Asymmetric Least Squares (AsLS) baseline correction.
/// Solves `(W + λ D₂ᵀD₂) z = W y` in O(N) via symmetric pentadiagonal LDLᵀ factorization.
/// `base_lambda` is auto-scaled by `(n / 1800)⁴` for sampling-density invariance.
pub fn asls_baseline(ys: &[f64], base_lambda: f64, p: f64, n_iter: usize) -> Vec<f64> {
    let n = ys.len();
    if n < 5 {
        return ys.to_vec();
    }
    let scale = (n as f64 / 1800.0).powi(4);
    let lam = (base_lambda * scale).clamp(1e2, 1e10);
    let mut w = vec![1.0_f64; n];
    let mut z = vec![0.0_f64; n];

    for _ in 0..n_iter {
        z = pentadiag_ldlt_solve(&w, ys, lam);
        for i in 0..n {
            w[i] = if ys[i] > z[i] { p } else { 1.0 - p };
        }
    }

    ys.iter().zip(&z).map(|(&y, &b)| (y - b).max(0.0)).collect()
}

/// Baek et al. (Analyst 2015, 140, 250) asymmetrically reweighted Penalized Least Squares (arPLS).
pub fn arpls_baseline(ys: &[f64], base_lambda: f64, ratio: f64, max_iter: usize) -> Vec<f64> {
    let n = ys.len();
    if n < 5 {
        return ys.to_vec();
    }
    let scale = (n as f64 / 1800.0).powi(4);
    let lam = (base_lambda * scale).clamp(1e2, 1e10);
    let mut w = vec![1.0_f64; n];
    let mut z = vec![0.0_f64; n];

    for _ in 0..max_iter {
        z = pentadiag_ldlt_solve(&w, ys, lam);
        let mut dn_sum = 0.0;
        let mut dn_cnt = 0usize;
        for i in 0..n {
            let d = ys[i] - z[i];
            if d < 0.0 {
                dn_sum += d;
                dn_cnt += 1;
            }
        }
        if dn_cnt < 2 {
            break;
        }
        let mean = dn_sum / dn_cnt as f64;
        let mut var = 0.0;
        for i in 0..n {
            let d = ys[i] - z[i];
            if d < 0.0 {
                let diff = d - mean;
                var += diff * diff;
            }
        }
        let std = (var / dn_cnt as f64).sqrt();
        if std < 1e-12 {
            break;
        }

        let mut diff_sq = 0.0;
        let mut w_sq = 0.0;
        for i in 0..n {
            let d = ys[i] - z[i];
            let w_new = if d <= 0.0 {
                1.0
            } else {
                let arg = (2.0 * (d - (2.0 * std - mean)) / std).clamp(-60.0, 60.0);
                1.0 / (1.0 + arg.exp())
            };
            let dw = w_new - w[i];
            diff_sq += dw * dw;
            w_sq += w[i] * w[i];
            w[i] = w_new;
        }
        if (diff_sq.sqrt() / w_sq.sqrt().max(1e-12)) < ratio {
            break;
        }
    }
    ys.iter().zip(&z).map(|(&y, &b)| (y - b).max(0.0)).collect()
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

#[derive(Debug, Clone, Copy)]
struct Complex {
    re: f64,
    im: f64,
}

impl Complex {
    fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }
    fn add(self, o: Self) -> Self {
        Self {
            re: self.re + o.re,
            im: self.im + o.im,
        }
    }
    fn sub(self, o: Self) -> Self {
        Self {
            re: self.re - o.re,
            im: self.im - o.im,
        }
    }
    fn mul(self, o: Self) -> Self {
        Self {
            re: self.re * o.re - self.im * o.im,
            im: self.re * o.im + self.im * o.re,
        }
    }
    fn div(self, o: Self) -> Self {
        let den = (o.re * o.re + o.im * o.im).max(1e-30);
        Self {
            re: (self.re * o.re + self.im * o.im) / den,
            im: (self.im * o.re - self.re * o.im) / den,
        }
    }
    fn norm_sq(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    fn sqrt(self) -> Self {
        let r = self.re.hypot(self.im);
        let re = ((r + self.re) * 0.5).max(0.0).sqrt();
        let im_mag = ((r - self.re) * 0.5).max(0.0).sqrt();
        let im = if self.im < 0.0 { -im_mag } else { im_mag };
        Self { re, im }
    }
}

fn fresnel_atr_absorbance_scalar(n1: f64, n2_real: f64, k: f64, sin_th: f64, cos_th: f64) -> f64 {
    let n2 = Complex::new(n2_real, k.max(0.0));
    let ratio = Complex::new(n1 * sin_th, 0.0).div(n2);
    let one_minus = Complex::new(1.0, 0.0).sub(ratio.mul(ratio));
    let mut cos_t = one_minus.sqrt();
    // Enforce evanescent decay into the sample: Im(n2 * cos_t) >= 0
    if n2.mul(cos_t).im < 0.0 {
        cos_t = Complex::new(-cos_t.re, -cos_t.im);
    }
    let n1_c = Complex::new(n1 * cos_th, 0.0);
    let n2_ct = n2.mul(cos_t);
    let rs = n1_c.sub(n2_ct).div(n1_c.add(n2_ct));

    let n2_c = n2.mul(Complex::new(cos_th, 0.0));
    let n1_ct = cos_t.mul(Complex::new(n1, 0.0));
    let rp = n2_c.sub(n1_ct).div(n2_c.add(n1_ct));

    let r_unpol = 0.5 * (rs.norm_sq() + rp.norm_sq());
    -(r_unpol.clamp(1e-6, 1.0)).log10()
}

/// Ohta & Ishida (1988) Maclaurin Cauchy principal value Kramers-Kronig transform
/// on a downsampled uniform grid (<= 450 pts for <1 ms speed).
fn kk_n_from_k(wn: &[f64], k: &[f64], n_inf: f64) -> Vec<f64> {
    let n = wn.len();
    if n < 4 || k.len() != n {
        return vec![n_inf; n];
    }
    const MAX_GRID: usize = 450;
    let (g_wn, g_k): (Vec<f64>, Vec<f64>) = if n <= MAX_GRID {
        (wn.to_vec(), k.to_vec())
    } else {
        let lo = wn[0];
        let hi = wn[n - 1];
        let mut gx = Vec::with_capacity(MAX_GRID);
        let mut gk = Vec::with_capacity(MAX_GRID);
        for i in 0..MAX_GRID {
            let x = lo + (hi - lo) * (i as f64 / (MAX_GRID - 1) as f64);
            gx.push(x);
            gk.push(interpolate_y(wn, k, x).unwrap_or(0.0));
        }
        (gx, gk)
    };

    let m = g_wn.len();
    let mut vk = vec![0.0_f64; m];
    let mut v2 = vec![0.0_f64; m];
    for j in 0..m {
        let dw = if j == 0 {
            g_wn[1] - g_wn[0]
        } else if j == m - 1 {
            g_wn[m - 1] - g_wn[m - 2]
        } else {
            0.5 * (g_wn[j + 1] - g_wn[j - 1])
        };
        vk[j] = g_wn[j] * g_k[j] * (2.0 * dw);
        v2[j] = g_wn[j] * g_wn[j];
    }

    let mut g_n = vec![n_inf; m];
    let inv_pi2 = 2.0 / std::f64::consts::PI;
    for i in 0..m {
        let vi2 = v2[i];
        let mut acc = 0.0;
        let mut j = if i % 2 == 0 { 1 } else { 0 };
        while j < m {
            let denom = v2[j] - vi2;
            if denom.abs() > 1e-9 {
                acc += vk[j] / denom;
            }
            j += 2;
        }
        g_n[i] = n_inf + inv_pi2 * acc;
    }

    if n <= MAX_GRID {
        g_n
    } else {
        wn.iter()
            .map(|&x| interpolate_y(&g_wn, &g_n, x).unwrap_or(n_inf))
            .collect()
    }
}

/// Advanced ATR-to-Transmission correction (Fresnel + Kramers-Kronig anomalous dispersion +
/// Harrick penetration depth scaling). Input & output are in linear Absorbance units.
pub fn advanced_atr_correction(
    wn: &[f64],
    abs_ys: &[f64],
    n1: f64,
    n_sample: f64,
    theta_deg: f64,
) -> Vec<f64> {
    let n = wn.len().min(abs_ys.len());
    if n < 4 || n1 <= n_sample + 0.05 {
        return abs_ys.to_vec();
    }
    let th = theta_deg.to_radians();
    let (s, c) = (th.sin(), th.cos());
    let n21 = n_sample / n1;
    let rad = (s * s - n21 * n21).max(1e-6).sqrt();
    let de_perp = n21 * c / (std::f64::consts::PI * n1 * (1.0 - n21 * n21) * rad);
    let de_par = de_perp * (2.0 * s * s - n21 * n21) / ((1.0 + n21 * n21) * s * s - n21 * n21);
    let scale_k = (4.0 * std::f64::consts::PI / 10.0_f64.ln()) * n_sample * 0.5 * (de_perp + de_par);

    let k_limit = if n1 < 3.0 { 0.80 } else { 1.40 };
    let a_single_max = fresnel_atr_absorbance_scalar(n1, n_sample, k_limit, s, c);
    let a_obs_max = abs_ys.iter().copied().fold(0.0_f64, |m, v| m.max(v.max(0.0)));
    if a_obs_max <= 1e-9 || a_single_max <= 1e-9 {
        return abs_ys.to_vec();
    }

    // Auto-estimate effective bounces if input is multi-bounce or high-absorbance
    let bounces = (a_obs_max / (0.90 * a_single_max)).ceil().max(1.0);
    let a_single: Vec<f64> = abs_ys[..n].iter().map(|&v| v.max(0.0) / bounces).collect();

    let mut k_est: Vec<f64> = a_single
        .iter()
        .map(|&a| (a / scale_k.max(1e-6)).clamp(0.0, k_limit))
        .collect();
    let n_cap = n1 * s * 0.992;

    for _ in 0..8 {
        let n_est = kk_n_from_k(&wn[..n], &k_est, n_sample);
        for i in 0..n {
            let target = a_single[i];
            let k_solved = if target < 1e-7 {
                0.0
            } else {
                let ni = n_est[i].clamp(1.05, n_cap);
                let mut lo = 0.0_f64;
                let mut hi = k_limit;
                if target >= fresnel_atr_absorbance_scalar(n1, ni, hi, s, c) {
                    hi
                } else {
                    for _ in 0..16 {
                        let mid = 0.5 * (lo + hi);
                        if fresnel_atr_absorbance_scalar(n1, ni, mid, s, c) < target {
                            lo = mid;
                        } else {
                            hi = mid;
                        }
                    }
                    0.5 * (lo + hi)
                }
            };
            k_est[i] = 0.5 * k_est[i] + 0.5 * k_solved;
        }
    }

    (0..n)
        .map(|i| bounces * (wn[i].max(1.0) / 1000.0) * k_est[i] * scale_k)
        .collect()
}

/// Kubelka-Munk remission transform F(R) = (1 - R)^2 / (2R) for DRIFTS and UV-Vis DRS.
/// Automatically detects %R (0..100) vs fractional R (0..1).
pub fn kubelka_munk_transform(ys: &[f64]) -> Vec<f64> {
    let y_max = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let is_pct = y_max > 1.5;
    ys.iter()
        .map(|&v| {
            let r = (if is_pct { v / 100.0 } else { v }).clamp(0.001, 0.9999);
            let d = 1.0 - r;
            (d * d) / (2.0 * r)
        })
        .collect()
}

/// Linearly bridge atmospheric CO₂ doublet (2280..2410 cm⁻¹) on FTIR traces.
pub fn suppress_co2_band(xs: &[f64], ys: &[f64]) -> Vec<f64> {
    let n = xs.len().min(ys.len());
    if n < 4 {
        return ys.to_vec();
    }
    let i0 = xs.partition_point(|&v| v < 2280.0);
    let i1 = xs.partition_point(|&v| v <= 2410.0);
    if i0 == 0 || i1 >= n || i1 <= i0 {
        return ys.to_vec();
    }
    let mut out = ys.to_vec();
    let (x0, y0) = (xs[i0 - 1], ys[i0 - 1]);
    let (x1, y1) = (xs[i1], ys[i1]);
    let span = (x1 - x0).max(1e-12);
    for i in i0..i1 {
        let t = ((xs[i] - x0) / span).clamp(0.0, 1.0);
        out[i] = y0 + t * (y1 - y0);
    }
    out
}

pub fn process_trace(spec: &Spectrum, cfg: ProcessConfig) -> (Vec<f64>, Vec<f64>, YUnit) {
    let xs = spec.x.clone();
    // Step A: Optical / Y-unit conversion
    let (mut ys, mut y_unit) = if cfg.optics == OpticsMode::KubelkaMunk
        && matches!(spec.kind, SpectrumKind::Ftir | SpectrumKind::UvVis)
    {
        let refl: Vec<f64> = if spec.raw_y_unit == YUnit::Absorbance {
            spec.y.iter().copied().map(absorbance_to_transmittance).collect()
        } else {
            spec.y.clone()
        };
        (kubelka_munk_transform(&refl), YUnit::Absorbance)
    } else if let Some(n1) = cfg.optics.atr_n1().filter(|_| spec.kind == SpectrumKind::Ftir) {
        let abs: Vec<f64> = match spec.raw_y_unit {
            YUnit::Transmittance => {
                spec.y.iter().copied().map(transmittance_to_absorbance).collect()
            }
            _ => spec.y.clone(),
        };
        let corr = advanced_atr_correction(&xs, &abs, n1, 1.50, 45.0);
        if cfg.y_pref == YAxisPreference::Native
            && cfg.deriv == DerivOrder::Raw
            && !matches!(cfg.norm, NormMode::VectorSnv | NormMode::Area)
        {
            (
                corr.into_iter().map(absorbance_to_transmittance).collect(),
                YUnit::Transmittance,
            )
        } else {
            (corr, YUnit::Absorbance)
        }
    } else {
        match spec.kind {
            SpectrumKind::Ftir => {
                if cfg.y_pref == YAxisPreference::Swapped {
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
                // Default = Absorbance; `YAxisPreference::Swapped` toggles %T <-> Abs.
                let want_transmittance = cfg.y_pref == YAxisPreference::Swapped;
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
        }
    };

    // Step B: Atmospheric CO₂ suppression
    if cfg.atm_co2 && spec.kind == SpectrumKind::Ftir {
        ys = suppress_co2_band(&xs, &ys);
    }

    // Step C: Smoothing
    if cfg.smooth != SmoothLevel::Off {
        ys = savitzky_golay(&ys, cfg.smooth);
    }

    // Step D: Baseline correction
    if cfg.baseline != BaselineMode::Off {
        let lam = if spec.kind == SpectrumKind::UvVis { 5e5 } else { 1e6 };
        if y_unit == YUnit::Transmittance {
            let abs: Vec<f64> = ys.iter().copied().map(transmittance_to_absorbance).collect();
            let corr = match cfg.baseline {
                BaselineMode::Rubberband => rubberband_baseline(&xs, &abs, false),
                BaselineMode::Asls => asls_baseline(&abs, lam, 0.01, 10),
                BaselineMode::Arpls => arpls_baseline(&abs, lam, 1e-3, 40),
                BaselineMode::Off => abs,
            };
            ys = corr.into_iter().map(absorbance_to_transmittance).collect();
        } else {
            ys = match cfg.baseline {
                BaselineMode::Rubberband => rubberband_baseline(&xs, &ys, false),
                BaselineMode::Asls => asls_baseline(&ys, lam, 0.01, 10),
                BaselineMode::Arpls => arpls_baseline(&ys, lam, 1e-3, 40),
                BaselineMode::Off => ys,
            };
        }
    }

    // Step E: Normalization
    if ys.len() >= 2 {
        match cfg.norm {
            NormMode::Off => {}
            NormMode::MinMax => {
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
            NormMode::VectorSnv => {
                if y_unit == YUnit::Transmittance {
                    ys = ys.into_iter().map(transmittance_to_absorbance).collect();
                    y_unit = YUnit::Absorbance;
                }
                let n_f = ys.len() as f64;
                let mean = ys.iter().sum::<f64>() / n_f;
                let var = ys.iter().map(|&v| (v - mean) * (v - mean)).sum::<f64>() / n_f;
                let std = var.sqrt().max(1e-12);
                for v in &mut ys {
                    *v = (*v - mean) / std;
                }
            }
            NormMode::Area => {
                if y_unit == YUnit::Transmittance {
                    ys = ys.into_iter().map(transmittance_to_absorbance).collect();
                    y_unit = YUnit::Absorbance;
                }
                let mut area = 0.0;
                for k in 0..(xs.len().min(ys.len()) - 1) {
                    let dx = (xs[k + 1] - xs[k]).abs();
                    area += 0.5 * (ys[k].abs() + ys[k + 1].abs()) * dx;
                }
                if area > 1e-12 {
                    for v in &mut ys {
                        *v /= area;
                    }
                }
            }
        }
    }

    // Step F: Derivative Spectroscopy
    if cfg.deriv != DerivOrder::Raw {
        if y_unit == YUnit::Transmittance {
            ys = ys.into_iter().map(transmittance_to_absorbance).collect();
        }
        ys = spectral_derivative(&xs, &ys, cfg.deriv);
        y_unit = YUnit::Intensity;
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

fn estimate_noise_mad(s: &[f64]) -> f64 {
    if s.len() < 4 {
        return 0.0;
    }
    let mut diffs: Vec<f64> = s.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    diffs.sort_by(|a, b| a.total_cmp(b));
    let med = diffs[diffs.len() / 2];
    med / (0.6745 * std::f64::consts::SQRT_2)
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
    let noise = estimate_noise_mad(&s);
    let gate = min_prom.max(3.0 * noise);

    let mut candidates: Vec<(usize, usize, usize, f64)> = Vec::new();
    let mut i = 1;
    while i < n - 1 {
        if s[i] > s[i - 1] {
            let mut j = i;
            while j + 1 < n - 1 && (s[j + 1] - s[i]).abs() <= 1e-12 {
                j += 1;
            }
            if s[j + 1] < s[i] {
                let mid = (i + j) / 2;
                // Left scan until higher point or edge
                let mut left_min = s[mid];
                let mut l = mid;
                while l > 0 {
                    l -= 1;
                    if s[l] < left_min {
                        left_min = s[l];
                    }
                    if s[l] > s[mid] {
                        break;
                    }
                }
                // Right scan until higher point or edge
                let mut right_min = s[mid];
                let mut r = mid;
                while r + 1 < n {
                    r += 1;
                    if s[r] < right_min {
                        right_min = s[r];
                    }
                    if s[r] > s[mid] {
                        break;
                    }
                }
                let prom = s[mid] - left_min.max(right_min);
                if prom >= gate {
                    candidates.push((mid, i, j, prom));
                }
                i = j + 1;
                continue;
            }
        }
        i += 1;
    }

    // Sort strongest prominence first, enforce minimum X distance
    candidates.sort_by(|a, b| b.3.total_cmp(&a.3));
    let min_dx = match kind {
        SpectrumKind::Ftir => 14.0,
        SpectrumKind::Hnmr => 0.06,
        SpectrumKind::Cnmr => 0.7,
        SpectrumKind::UvVis => 10.0,
    };

    let mut chosen: Vec<Peak> = Vec::new();
    for (idx, p_left, p_right, prom) in candidates {
        let refined_x = if p_left < p_right {
            0.5 * (xs[p_left] + xs[p_right])
        } else {
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
            xs[idx] + delta * dx_step
        };

        if chosen.iter().any(|p| (p.x - refined_x).abs() < min_dx) {
            continue;
        }

        // Compute FWHM at half-prominence height
        let half_h = s[idx] - 0.5 * prom;
        let mut x_left = None;
        let mut l = idx;
        while l > 0 {
            if s[l - 1] <= half_h {
                let denom = (s[l] - s[l - 1]).max(1e-12);
                let t = ((half_h - s[l - 1]) / denom).clamp(0.0, 1.0);
                x_left = Some(xs[l - 1] + t * (xs[l] - xs[l - 1]));
                break;
            }
            l -= 1;
        }
        let mut x_right = None;
        let mut r = idx;
        while r + 1 < n {
            if s[r + 1] <= half_h {
                let denom = (s[r] - s[r + 1]).max(1e-12);
                let t = ((s[r] - half_h) / denom).clamp(0.0, 1.0);
                x_right = Some(xs[r] + t * (xs[r + 1] - xs[r]));
                break;
            }
            r += 1;
        }
        let fwhm = match (x_left, x_right) {
            (Some(xl), Some(xr)) => Some((xr - xl).abs()),
            _ => None,
        };
        let snr = if noise > 1e-12 { prom / noise } else { 999.0 };

        chosen.push(Peak {
            x: refined_x,
            y: ys[idx],
            prominence: prom / span,
            fwhm,
            snr,
            label: String::new(),
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

    let peak_xs: Vec<f64> = chosen.iter().map(|p| p.x).collect();
    for p in &mut chosen {
        p.label = best_peak_label_in_context(kind, p.x, &peak_xs, SolventId::Cdcl3);
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
    let abs_cfg = ProcessConfig::ftir_absorbance();
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
    pub fit_e_min: f64,
    pub fit_e_max: f64,
    pub mode: TaucMode,
}

/// Tauc plot band gap from a UV-Vis absorption edge (Makuła 2018 standard).
///
/// `(αhν)^(1/n)` versus `hν` is fitted in a sliding linear window on the
/// first fundamental absorption onset; the x-intercept of the best-fit line is Eg.
#[allow(dead_code)]
pub fn compute_tauc_band_gap(spec: &Spectrum, mode: TaucMode) -> Option<TaucResult> {
    compute_tauc_band_gap_with_optics(spec, mode, OpticsMode::Off)
}

pub fn compute_tauc_band_gap_with_optics(
    spec: &Spectrum,
    mode: TaucMode,
    optics: OpticsMode,
) -> Option<TaucResult> {
    if spec.kind != SpectrumKind::UvVis || spec.x.len() < 12 || spec.x.len() != spec.y.len() {
        return None;
    }

    let ys_abs: Vec<f64> = if optics == OpticsMode::KubelkaMunk {
        let refl: Vec<f64> = if spec.raw_y_unit == YUnit::Absorbance {
            spec.y.iter().copied().map(absorbance_to_transmittance).collect()
        } else {
            spec.y.clone()
        };
        kubelka_munk_transform(&refl)
    } else if spec.raw_y_unit == YUnit::Transmittance {
        spec.y.iter().copied().map(transmittance_to_absorbance).collect()
    } else {
        spec.y.clone()
    };

    // 1. Convert λ (nm) -> photon energy E (eV) and sort ASCENDING in E.
    let mut pairs: Vec<(f64, f64)> = spec
        .x
        .iter()
        .zip(&ys_abs)
        .filter_map(|(&nm, &a)| {
            if nm.is_finite() && a.is_finite() && nm > 100.0 && nm < 2500.0 {
                Some((1239.84193 / nm, a))
            } else {
                None
            }
        })
        .collect();
    if pairs.len() < 12 {
        return None;
    }
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));

    let ev_raw: Vec<f64> = pairs.iter().map(|p| p.0).collect();
    let y_raw: Vec<f64> = pairs.iter().map(|p| p.1).collect();

    // 2. Smart sub-gap baseline: only subtract tail percentile if sub-gap tail is flat.
    let y_min = y_raw.iter().copied().fold(f64::INFINITY, f64::min);
    let y_max = y_raw.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let y_span = y_max - y_min;
    if y_span <= 1e-12 {
        return None;
    }
    let tail_n = (y_raw.len() / 6).max(5).min(y_raw.len());
    let mut tail = y_raw[..tail_n].to_vec();
    tail.sort_by(|a, b| a.total_cmp(b));
    let tail_spread = tail.last().copied().unwrap_or(0.0) - tail.first().copied().unwrap_or(0.0);
    let baseline = if tail_spread < 0.08 * y_span {
        let p15 = ((tail.len() - 1) as f64 * 0.15).round() as usize;
        tail[p15].max(0.0)
    } else {
        y_min.max(0.0)
    };
    let alpha_raw: Vec<f64> = y_raw.iter().map(|&y| (y - baseline).max(0.0)).collect();

    // 3. Resample onto uniform energy grid (ΔE ≈ 0.01 eV) so every fit window has identical eV span.
    let e_lo = ev_raw[0];
    let e_hi = *ev_raw.last()?;
    if e_hi - e_lo < 0.30 {
        return None;
    }
    let n_grid = (((e_hi - e_lo) / 0.01).round() as usize + 1).clamp(40, 600);
    let de_step = (e_hi - e_lo) / (n_grid - 1) as f64;
    let mut ev = Vec::with_capacity(n_grid);
    let mut a_grid = Vec::with_capacity(n_grid);
    let mut tauc = Vec::with_capacity(n_grid);
    for i in 0..n_grid {
        let e = e_lo + i as f64 * de_step;
        let a = interpolate_y(&ev_raw, &alpha_raw, e).unwrap_or(0.0).max(0.0);
        let ahv = a * e;
        let ord = match mode {
            TaucMode::Direct => ahv * ahv,
            TaucMode::Indirect => ahv.sqrt(),
        };
        ev.push(e);
        a_grid.push(a);
        tauc.push(ord);
    }

    let a_max = a_grid.iter().copied().fold(0.0_f64, f64::max);
    if a_max <= 1e-12 {
        return None;
    }

    // 4. Locate the FIRST fundamental absorption onset (scanning low -> high E).
    let i_onset = a_grid.iter().position(|&a| a >= 0.10 * a_max)?;
    let max_span_pts = ((0.85 / de_step).round() as usize).max(15);
    let mut i_ceil = (i_onset + max_span_pts).min(n_grid - 1);

    // Stop at first absorption band peak if a real dip follows (avoids jumping to a 220 nm deep-UV wall).
    if i_ceil > i_onset + 11 {
        for j in (i_onset + 8)..(i_ceil - 3) {
            let local_max = a_grid[j - 3..=j + 3].iter().copied().fold(0.0_f64, f64::max);
            if a_grid[j] >= 0.25 * a_max && a_grid[j] >= local_max {
                let post_end = (j + 20).min(n_grid);
                let post_min = a_grid[j..post_end].iter().copied().fold(f64::INFINITY, f64::min);
                if post_min < 0.94 * a_grid[j] {
                    i_ceil = j;
                    break;
                }
            }
        }
    }

    let win = ((0.22 / de_step).round() as usize).clamp(6, 30);
    if i_ceil < i_onset + win {
        i_ceil = (i_onset + win + 2).min(n_grid - 1);
    }
    let edge_max_tauc = tauc[i_onset..=i_ceil].iter().copied().fold(0.0_f64, f64::max);
    if edge_max_tauc <= 1e-12 {
        return None;
    }

    // 5. Slide linear fit window strictly across the first rising absorption edge.
    let mut best: Option<(f64, f64, f64, f64, f64, f64)> = None;
    let start_lo = i_onset.saturating_sub(4);
    let start_hi = i_ceil.saturating_sub(win);
    for start in start_lo..=start_hi {
        let end = start + win;
        if end > n_grid {
            break;
        }
        let se = &ev[start..end];
        let st = &tauc[start..end];
        let m = win as f64;
        let mean_y = st.iter().sum::<f64>() / m;
        if mean_y < 0.12 * edge_max_tauc || mean_y > 0.88 * edge_max_tauc {
            continue;
        }
        let mean_x = se.iter().sum::<f64>() / m;
        let mut sxx = 0.0;
        let mut sxy = 0.0;
        let mut syy = 0.0;
        for (&x, &y) in se.iter().zip(st) {
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
        if slope <= 0.0 || r2 < 0.90 {
            continue;
        }
        let eg = -intercept / slope;
        if !(0.30..=6.50).contains(&eg) || eg >= se[0] + 0.05 {
            continue;
        }
        let score = r2.powi(6) * (slope / (edge_max_tauc + 1e-12));
        if best.is_none_or(|(bs, _, _, _, _, _)| score > bs) {
            best = Some((score, eg, r2, se[0], se[win - 1], slope));
        }
    }

    let (_, eg, r2, fit_e_min, fit_e_max, _) = best?;
    Some(TaucResult {
        band_gap_ev: eg,
        edge_nm: 1239.84193 / eg.max(1e-9),
        r_squared: r2,
        fit_e_min,
        fit_e_max,
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

// ============================================================================
// INTERACTIVE SPECTRAL MATH ENGINE
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpectralOp {
    Subtract, // A - k * B
    Add,      // A + k * B
    Ratio,    // A / (k * B)
    Average,  // Mean of compatible visible traces
}

impl SpectralOp {
    pub fn next(self) -> Self {
        match self {
            Self::Subtract => Self::Add,
            Self::Add => Self::Ratio,
            Self::Ratio => Self::Average,
            Self::Average => Self::Subtract,
        }
    }

    pub fn prev(self) -> Self {
        match self {
            Self::Subtract => Self::Average,
            Self::Add => Self::Subtract,
            Self::Ratio => Self::Add,
            Self::Average => Self::Ratio,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Subtract => "A - k·B (Subtract)",
            Self::Add => "A + k·B (Add)",
            Self::Ratio => "A / (k·B) (Ratio)",
            Self::Average => "Average Visible",
        }
    }
}

fn to_linear_space(spec: &Spectrum) -> Vec<f64> {
    if spec.raw_y_unit == YUnit::Transmittance {
        spec.y
            .iter()
            .copied()
            .map(transmittance_to_absorbance)
            .collect()
    } else {
        spec.y.clone()
    }
}

/// Computes baseline-invariant optimal subtraction factor k* = sum(dA * dB) / sum(dB^2)
/// over overlapping points inside [x_lo, x_hi] in linear space (clamped to [0.0, 20.0]).
/// Falls back to full overlapping range if [x_lo, x_hi] has < 4 overlapping points.
pub fn optimal_subtraction_factor(
    a: &Spectrum,
    b: &Spectrum,
    x_lo: f64,
    x_hi: f64,
) -> Option<f64> {
    if a.kind != b.kind || a.x.len() < 4 || b.x.len() < 4 {
        return None;
    }
    let a_lin = to_linear_space(a);
    let b_lin = to_linear_space(b);
    let b_lo = *b.x.first()?;
    let b_hi = *b.x.last()?;

    let collect_overlap = |lo: f64, hi: f64| -> (Vec<f64>, Vec<f64>) {
        let mut va = Vec::new();
        let mut vb = Vec::new();
        for (&x, &ya) in a.x.iter().zip(&a_lin) {
            if x >= lo && x <= hi && x >= b_lo && x <= b_hi {
                if let Some(yb) = interpolate_y(&b.x, &b_lin, x) {
                    va.push(ya);
                    vb.push(yb);
                }
            }
        }
        (va, vb)
    };

    let (mut va, mut vb) = collect_overlap(x_lo.min(x_hi), x_lo.max(x_hi));
    if va.len() < 4 {
        (va, vb) = collect_overlap(f64::NEG_INFINITY, f64::INFINITY);
    }
    if va.len() < 4 {
        return None;
    }

    let da: Vec<f64> = va.windows(2).map(|w| w[1] - w[0]).collect();
    let db: Vec<f64> = vb.windows(2).map(|w| w[1] - w[0]).collect();
    if da.len() < 3 {
        return None;
    }
    // Removing the mean derivative rejects constant linear-baseline tilt as well as offsets.
    let mean_da = da.iter().sum::<f64>() / da.len() as f64;
    let mean_db = db.iter().sum::<f64>() / db.len() as f64;
    let mut num = 0.0;
    let mut den = 0.0;
    for (&dai, &dbi) in da.iter().zip(&db) {
        let ac = dai - mean_da;
        let bc = dbi - mean_db;
        num += ac * bc;
        den += bc * bc;
    }
    if den <= 1e-18 {
        return None;
    }
    Some((num / den).clamp(0.0, 20.0))
}

pub fn compute_spectral_math(
    a: &Spectrum,
    b: &Spectrum,
    op: SpectralOp,
    k: f64,
) -> Option<Spectrum> {
    if op == SpectralOp::Average {
        return average_spectra(&[a, b]);
    }
    if a.kind != b.kind || a.x.is_empty() || b.x.is_empty() || !k.is_finite() {
        return None;
    }
    let a_lin = to_linear_space(a);
    let b_lin = to_linear_space(b);
    let b_lo = *b.x.first()?;
    let b_hi = *b.x.last()?;

    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for (&x, &ya) in a.x.iter().zip(&a_lin) {
        if x >= b_lo && x <= b_hi {
            if let Some(yb) = interpolate_y(&b.x, &b_lin, x) {
                let val = match op {
                    SpectralOp::Subtract => ya - k * yb,
                    SpectralOp::Add => ya + k * yb,
                    SpectralOp::Ratio => ya / (k * yb).max(1e-6),
                    SpectralOp::Average => unreachable!(),
                };
                xs.push(x);
                ys.push(val);
            }
        }
    }
    if xs.len() < 2 {
        return None;
    }

    let raw_y_unit = if matches!(a.kind, SpectrumKind::Ftir | SpectrumKind::UvVis) {
        YUnit::Absorbance
    } else {
        YUnit::Intensity
    };
    let name = match op {
        SpectralOp::Subtract => format!("{} - {:.2}*{}", a.name, k, b.name),
        SpectralOp::Add => format!("{} + {:.2}*{}", a.name, k, b.name),
        SpectralOp::Ratio => format!("{} / ({:.2}*{})", a.name, k, b.name),
        SpectralOp::Average => unreachable!(),
    };
    let op_desc = match op {
        SpectralOp::Subtract => format!("{} - {:.4} * {}", a.name, k, b.name),
        SpectralOp::Add => format!("{} + {:.4} * {}", a.name, k, b.name),
        SpectralOp::Ratio => format!("{} / ({:.4} * {})", a.name, k, b.name),
        SpectralOp::Average => unreachable!(),
    };

    Some(Spectrum {
        name,
        source_path: None,
        kind: a.kind,
        raw_y_unit,
        x: xs,
        y: ys,
        visible: true,
        color_idx: a.color_idx + 1,
        metadata: vec![("Spectral Math".to_string(), op_desc)],
    })
}

pub fn average_spectra(traces: &[&Spectrum]) -> Option<Spectrum> {
    let first = *traces.first()?;
    let compat: Vec<&Spectrum> = traces
        .iter()
        .copied()
        .filter(|s| s.kind == first.kind && s.x.len() >= 2 && s.x.len() == s.y.len())
        .collect();
    if compat.is_empty() {
        return None;
    }

    let lin_vecs: Vec<Vec<f64>> = compat.iter().map(|s| to_linear_space(s)).collect();
    let mut lo = f64::NEG_INFINITY;
    let mut hi = f64::INFINITY;
    for s in &compat {
        if let (Some(&s_lo), Some(&s_hi)) = (s.x.first(), s.x.last()) {
            lo = lo.max(s_lo);
            hi = hi.min(s_hi);
        }
    }
    if lo >= hi {
        return None;
    }

    let mut xs = Vec::new();
    let mut ys = Vec::new();
    let count = compat.len() as f64;
    for &x in &first.x {
        if x >= lo && x <= hi {
            let mut sum = 0.0;
            let mut ok = true;
            for (s, lin) in compat.iter().zip(&lin_vecs) {
                if let Some(yv) = interpolate_y(&s.x, lin, x) {
                    sum += yv;
                } else {
                    ok = false;
                    break;
                }
            }
            if ok {
                let avg = sum / count;
                let out_y = if first.raw_y_unit == YUnit::Transmittance {
                    absorbance_to_transmittance(avg)
                } else {
                    avg
                };
                xs.push(x);
                ys.push(out_y);
            }
        }
    }
    if xs.len() < 2 {
        return None;
    }

    Some(Spectrum {
        name: format!("Avg ({} traces)", compat.len()),
        source_path: None,
        kind: first.kind,
        raw_y_unit: first.raw_y_unit,
        x: xs,
        y: ys,
        visible: true,
        color_idx: first.color_idx + 1,
        metadata: vec![(
            "Spectral Math".to_string(),
            format!("Average of {} {} traces", compat.len(), first.kind.display_name()),
        )],
    })
}
