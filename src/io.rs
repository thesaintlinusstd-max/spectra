use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::{RwLock, OnceLock};

use anyhow::{anyhow, Context, Result};

use crate::data::{
    mma_monomer_spectrum, pmma_hnmr_spectrum, pmma_polymer_spectrum, polystyrene_ftir_spectrum,
    tio2_zno_uvvis_spectrum, Spectrum, SpectrumKind, YUnit,
};

/// Instrument exports are frequently Latin-1 / Windows-1252, so a strict
/// UTF-8 read makes perfectly valid spectra silently unreadable.
fn decode_spectrum_bytes(bytes: Vec<u8>) -> String {
    match String::from_utf8(bytes) {
        Ok(text) => text,
        Err(err) => err
            .into_bytes()
            .into_iter()
            .map(|b| b as char)
            .collect::<String>(),
    }
}

pub fn load_spectrum_file(path: &Path, kind_override: Option<SpectrumKind>) -> Result<Spectrum> {
    let bytes = fs::read(path)
        .with_context(|| format!("Failed to read file: {}", path.display()))?;
    let name = path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("spectrum")
        .to_string();
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();

    // Native binary instrument formats must never be decoded as UTF-8 text.
    if let Some(parsed) = crate::binary_io::try_parse_binary(&bytes, &name, &ext, kind_override) {
        let mut spec = parsed?;
        spec.source_path = Some(path.display().to_string());
        return Ok(spec);
    }

    let text = decode_spectrum_bytes(bytes);
    let text = text.trim_start_matches('\u{feff}');

    let mut spec = if matches!(ext.as_str(), "jdx" | "dx" | "jcamp")
        || text.trim_start().starts_with("##TITLE")
    {
        parse_jcamp_dx(text, name, kind_override)?
    } else {
        parse_delimited(text, name, kind_override)?
    };
    spec.source_path = Some(path.display().to_string());
    Ok(spec)
}

fn infer_kind(x_min: f64, x_max: f64, hint: Option<&str>) -> SpectrumKind {
    if let Some(h) = hint {
        let low = h.to_ascii_lowercase();
        if low.contains("uvvis")
            || low.contains("uv-vis")
            || low.contains("uv/vis")
            || low.contains("wavelength")
            || low.contains("nanometer")
            || low.contains("(nm)")
            || low.contains("extinction")
            || low.contains("photochemcad")
            || low.contains("ultraviolet")
        {
            return SpectrumKind::UvVis;
        }
        if low.contains("infrared") || low.contains("ftir") || low.contains("cm-1") || low.contains("1/cm") {
            return SpectrumKind::Ftir;
        }
        if low.contains("13c") || low.contains("c13") || low.contains("carbon") {
            return SpectrumKind::Cnmr;
        }
        if low.contains("nmr") || low.contains("ppm") {
            return if x_max > 25.0 { SpectrumKind::Cnmr } else { SpectrumKind::Hnmr };
        }
    }
    if x_max > 1300.0 {
        SpectrumKind::Ftir
    } else if x_min >= 100.0 && (260.0..=1300.0).contains(&x_max) {
        SpectrumKind::UvVis
    } else if x_max > 24.0 {
        SpectrumKind::Cnmr
    } else {
        SpectrumKind::Hnmr
    }
}

fn infer_y_unit(kind: SpectrumKind, ys: &[f64], hint: Option<&str>) -> YUnit {
    if kind == SpectrumKind::UvVis {
        if let Some(h) = hint {
            let low = h.to_ascii_lowercase();
            if low.contains("%t") || low.contains("trans") {
                return YUnit::Transmittance;
            }
            if low.contains("absorb") || low.contains("abs") || low.contains("epsilon") || low.contains("log") {
                return YUnit::Absorbance;
            }
        }
        // UV-Vis defaults to Absorbance unless %T / transmittance is evident
        let y_max = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let y_min = ys.iter().copied().fold(f64::INFINITY, f64::min);
        if y_max > 15.0 && y_max <= 105.0 && y_min >= 0.0 && ys.iter().filter(|&&v| v > 20.0).count() > ys.len() / 2 {
            return YUnit::Transmittance;
        }
        return YUnit::Absorbance;
    }
    if kind != SpectrumKind::Ftir {
        return YUnit::Intensity;
    }
    if let Some(h) = hint {
        let low = h.to_ascii_lowercase();
        if low.contains("absorb") || low == "a" || low.contains("abs") {
            return YUnit::Absorbance;
        }
        if low.contains("trans") || low.contains("%t") {
            return YUnit::Transmittance;
        }
    }
    let y_max = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let y_min = ys.iter().copied().fold(f64::INFINITY, f64::min);
    if y_max > 10.0 && y_max <= 115.0 && y_min >= -5.0 {
        YUnit::Transmittance
    } else if y_max <= 1.05 && y_min >= 0.0 && ys.iter().filter(|&&v| v > 0.7).count() > ys.len() / 2 {
        // Fractional transmittance 0..1 -> caller normalizes to 0..100
        YUnit::Transmittance
    } else {
        YUnit::Absorbance
    }
}

fn finalize_points(mut pts: Vec<(f64, f64)>) -> Result<(Vec<f64>, Vec<f64>)> {
    pts.retain(|(x, y)| x.is_finite() && y.is_finite());
    if pts.len() < 2 {
        return Err(anyhow!("Not enough valid numeric (X, Y) data points found"));
    }
    pts.sort_by(|a, b| a.0.total_cmp(&b.0));

    // Deduplicate identical X points by averaging Y
    let mut xs = Vec::with_capacity(pts.len());
    let mut ys = Vec::with_capacity(pts.len());
    let mut i = 0;
    while i < pts.len() {
        let x0 = pts[i].0;
        let mut sum = pts[i].1;
        let mut cnt = 1.0;
        while i + 1 < pts.len() && (pts[i + 1].0 - x0).abs() < 1e-12 {
            i += 1;
            sum += pts[i].1;
            cnt += 1.0;
        }
        xs.push(x0);
        ys.push(sum / cnt);
        i += 1;
    }
    Ok((xs, ys))
}

fn parse_float_field(s: &str, eu_decimal: bool) -> Option<f64> {
    let trimmed = s.trim().trim_matches('"');
    if trimmed.is_empty() {
        return None;
    }
    if eu_decimal && trimmed.contains(',') && !trimmed.contains('.') {
        trimmed.replace(',', ".").parse::<f64>().ok()
    } else {
        trimmed.parse::<f64>().ok()
    }
}

pub fn parse_delimited(
    text: &str,
    name: String,
    kind_override: Option<SpectrumKind>,
) -> Result<Spectrum> {
    let mut header_hint = String::new();
    let mut pts: Vec<(f64, f64)> = Vec::new();
    let mut started_data = false;

    for raw_line in text.lines() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') || line.starts_with("//") {
            if !started_data {
                header_hint.push(' ');
                header_hint.push_str(line);
            }
            continue;
        }

        let parsed = if line.contains(';') {
            let parts: Vec<&str> = line.split(';').collect();
            if parts.len() >= 2 {
                parse_float_field(parts[0], true)
                    .and_then(|x| parse_float_field(parts[1], true).map(|y| (x, y)))
            } else {
                None
            }
        } else if line.contains('\t') {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() >= 2 {
                parse_float_field(parts[0], true)
                    .and_then(|x| parse_float_field(parts[1], true).map(|y| (x, y)))
            } else {
                None
            }
        } else if line.contains(',') {
            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() >= 2 {
                parse_float_field(parts[0], false)
                    .and_then(|x| parse_float_field(parts[1], false).map(|y| (x, y)))
            } else {
                None
            }
        } else {
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() >= 2 {
                parse_float_field(parts[0], false)
                    .and_then(|x| parse_float_field(parts[1], false).map(|y| (x, y)))
            } else {
                None
            }
        };

        if let Some((x, y)) = parsed {
            started_data = true;
            pts.push((x, y));
        } else if !started_data {
            header_hint.push(' ');
            header_hint.push_str(line);
        }
    }

    let (xs, mut ys) = finalize_points(pts)?;
    let kind = kind_override.unwrap_or_else(|| {
        infer_kind(xs[0], *xs.last().unwrap_or(&xs[0]), Some(&header_hint))
    });
    let raw_y_unit = infer_y_unit(kind, &ys, Some(&header_hint));

    // Normalize fractional transmittance (0..1) to %T (0..100)
    if kind == SpectrumKind::Ftir && raw_y_unit == YUnit::Transmittance {
        let y_max = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if y_max <= 1.05 {
            for v in &mut ys {
                *v *= 100.0;
            }
        }
    }

    Ok(Spectrum {
        name,
        source_path: None,
        kind,
        raw_y_unit,
        x: xs,
        y: ys,
        visible: true,
        color_idx: 0,
        metadata: vec![("Format".into(), "Delimited Text/CSV".into())],
    })
}

// ============================================================================
// FULL IUPAC JCAMP-DX PARSER (AFFN, PAC, SQZ, DIF, DUP, DIFDUP)
// ============================================================================

fn sqz_digit(c: char) -> Option<i64> {
    match c {
        '@' => Some(0),
        'A'..='I' => Some((c as i64) - ('A' as i64) + 1),
        'a'..='i' => Some(-((c as i64) - ('a' as i64) + 1)),
        _ => None,
    }
}

fn dif_digit(c: char) -> Option<i64> {
    match c {
        '%' => Some(0),
        'J'..='R' => Some((c as i64) - ('J' as i64) + 1),
        'j'..='r' => Some(-((c as i64) - ('j' as i64) + 1)),
        _ => None,
    }
}

fn dup_digit(c: char) -> Option<usize> {
    match c {
        'S'..='Z' => Some(((c as u8) - b'S' + 1) as usize),
        's' => Some(9),
        _ => None,
    }
}

#[derive(Debug, Clone, Copy)]
enum AsdfToken {
    Abs(f64),
    Dif(f64),
    Dup(usize),
}

fn tokenize_asdf_line(line: &str) -> Vec<AsdfToken> {
    let chars: Vec<char> = line.chars().collect();
    let n = chars.len();
    let mut i = 0;
    let mut tokens = Vec::new();

    while i < n {
        let c = chars[i];
        if c.is_whitespace() || c == ',' || c == ';' {
            i += 1;
            continue;
        }
        if c == '?' {
            tokens.push(AsdfToken::Abs(f64::NAN));
            i += 1;
            continue;
        }

        // SQZ pseudo-digit starts an absolute number
        if let Some(d0) = sqz_digit(c) {
            i += 1;
            let neg = d0 < 0;
            let mut s = String::new();
            if neg {
                s.push('-');
            }
            s.push_str(&d0.abs().to_string());
            while i < n && (chars[i].is_ascii_digit() || chars[i] == '.') {
                s.push(chars[i]);
                i += 1;
            }
            if let Ok(v) = s.parse::<f64>() {
                tokens.push(AsdfToken::Abs(v));
            }
            continue;
        }

        // DIF pseudo-digit starts a difference value
        if let Some(d0) = dif_digit(c) {
            i += 1;
            let neg = d0 < 0;
            let mut s = String::new();
            if neg {
                s.push('-');
            }
            s.push_str(&d0.abs().to_string());
            while i < n && (chars[i].is_ascii_digit() || chars[i] == '.') {
                s.push(chars[i]);
                i += 1;
            }
            if let Ok(v) = s.parse::<f64>() {
                tokens.push(AsdfToken::Dif(v));
            }
            continue;
        }

        // DUP pseudo-digit repeats previous value or previous diff (count - 1 more times)
        if let Some(d0) = dup_digit(c) {
            i += 1;
            let mut count = d0;
            while i < n && chars[i].is_ascii_digit() {
                count = count * 10 + ((chars[i] as u8 - b'0') as usize);
                i += 1;
            }
            tokens.push(AsdfToken::Dup(count));
            continue;
        }

        // Standard AFFN / PAC signed/unsigned float
        if c.is_ascii_digit() || c == '+' || c == '-' || c == '.' {
            let mut s = String::new();
            s.push(c);
            i += 1;
            while i < n {
                let ch = chars[i];
                if ch.is_ascii_digit() || ch == '.' {
                    s.push(ch);
                    i += 1;
                } else if (ch == 'E' || ch == 'e')
                    && i + 1 < n
                    && (chars[i + 1].is_ascii_digit()
                        || chars[i + 1] == '+'
                        || chars[i + 1] == '-')
                {
                    s.push(ch);
                    i += 1;
                    s.push(chars[i]);
                    i += 1;
                } else {
                    break;
                }
            }
            if let Ok(v) = s.parse::<f64>() {
                tokens.push(AsdfToken::Abs(v));
            }
            continue;
        }

        // Skip unrecognized token chars safely
        i += 1;
    }

    tokens
}

fn parse_xypoints_line(
    line: &str,
    xfactor: f64,
    yfactor: f64,
    out: &mut Vec<(f64, f64)>,
) {
    for chunk in line.split(';') {
        let nums: Vec<f64> = chunk
            .split(|c: char| c == ',' || c.is_whitespace())
            .filter(|s| !s.is_empty())
            .filter_map(|s| s.parse::<f64>().ok())
            .collect();
        let mut idx = 0;
        while idx + 1 < nums.len() {
            out.push((nums[idx] * xfactor, nums[idx + 1] * yfactor));
            idx += 2;
        }
    }
}

pub fn parse_jcamp_dx(
    text: &str,
    default_name: String,
    kind_override: Option<SpectrumKind>,
) -> Result<Spectrum> {
    let mut title = default_name;
    let mut data_type = String::new();
    let mut xunits = String::new();
    let mut yunits = String::new();
    let mut firstx: Option<f64> = None;
    let mut lastx: Option<f64> = None;
    let mut _npoints: Option<usize> = None;
    let mut deltax: Option<f64> = None;
    let mut xfactor: f64 = 1.0;
    let mut yfactor: f64 = 1.0;
    let mut observe_freq_mhz: Option<f64> = None;
    let mut metadata: Vec<(String, String)> = vec![("Format".into(), "JCAMP-DX".into())];

    #[derive(PartialEq, Eq)]
    enum DataMode {
        None,
        XyData,
        XyPoints,
    }
    let mut mode = DataMode::None;
    let mut xydata_lines: Vec<String> = Vec::new();
    let mut xypoints_pairs: Vec<(f64, f64)> = Vec::new();

    for raw in text.lines() {
        let stripped = match raw.find("$$") {
            Some(idx) => &raw[..idx],
            None => raw,
        };
        let line = stripped.trim();
        if line.is_empty() {
            continue;
        }

        if let Some(rest) = line.strip_prefix("##") {
            if let Some((k_raw, v_raw)) = rest.split_once('=') {
                let key = k_raw.trim().to_ascii_uppercase().replace([' ', '-', '_'], "");
                let val = v_raw.trim();
                mode = DataMode::None;
                match key.as_str() {
                    "TITLE" if !val.is_empty() => title = val.to_string(),
                    "DATATYPE" => {
                        data_type = val.to_string();
                        metadata.push(("Data Type".into(), val.to_string()));
                    }
                    "XUNITS" => {
                        xunits = val.to_string();
                        metadata.push(("X Units".into(), val.to_string()));
                    }
                    "YUNITS" => {
                        yunits = val.to_string();
                        metadata.push(("Y Units".into(), val.to_string()));
                    }
                    "FIRSTX" => firstx = val.parse().ok(),
                    "LASTX" => lastx = val.parse().ok(),
                    "NPOINTS" => _npoints = val.parse().ok(),
                    "DELTAX" => deltax = val.parse().ok(),
                    "XFACTOR" => {
                        if let Ok(v) = val.parse() {
                            xfactor = v;
                        }
                    }
                    "YFACTOR" => {
                        if let Ok(v) = val.parse() {
                            yfactor = v;
                        }
                    }
                    ".OBSERVEFREQUENCY" | "OBSERVEFREQUENCY" => {
                        observe_freq_mhz = val
                            .split_whitespace()
                            .next()
                            .and_then(|s| s.parse::<f64>().ok());
                        metadata.push(("Freq (MHz)".into(), val.to_string()));
                    }
                    ".SOLVENTNAME" | "SOLVENTNAME" => {
                        metadata.push(("Solvent".into(), val.to_string()));
                    }
                    "XYDATA" | "DATATABLE" => {
                        let up = val.to_ascii_uppercase();
                        if up.contains("XY..XY") {
                            mode = DataMode::XyPoints;
                            if let Some(pos) = val.find(')') {
                                let rest = val[pos + 1..].trim();
                                if !rest.is_empty() {
                                    parse_xypoints_line(rest, xfactor, yfactor, &mut xypoints_pairs);
                                }
                            }
                        } else {
                            mode = DataMode::XyData;
                        }
                    }
                    "XYPOINTS" | "PEAKTABLE" => {
                        mode = DataMode::XyPoints;
                        if let Some(pos) = val.find(')') {
                            let rest = val[pos + 1..].trim();
                            if !rest.is_empty() {
                                parse_xypoints_line(rest, xfactor, yfactor, &mut xypoints_pairs);
                            }
                        }
                    }
                    "END" if !xydata_lines.is_empty() || !xypoints_pairs.is_empty() => break,
                    _ => {}
                }
            }
            continue;
        }

        match mode {
            DataMode::None => {}
            DataMode::XyPoints => {
                parse_xypoints_line(line, xfactor, yfactor, &mut xypoints_pairs);
            }
            DataMode::XyData => {
                xydata_lines.push(line.to_string());
            }
        }
    }

    let pts = if !xypoints_pairs.is_empty() {
        xypoints_pairs
    } else if !xydata_lines.is_empty() {
        // Decompress XYDATA (X++(Y..Y)) with exact DIF checkpoint handling
        let mut raw_ys: Vec<f64> = Vec::new();
        let mut first_line_x: Option<f64> = None;
        let mut prev_line_was_dif = false;

        for line in &xydata_lines {
            let tokens = tokenize_asdf_line(line);
            if tokens.len() < 2 {
                continue;
            }
            // First token on every XYDATA line is the starting X coordinate
            if first_line_x.is_none() {
                if let AsdfToken::Abs(x0) = tokens[0] {
                    first_line_x = Some(x0 * xfactor);
                }
            }

            let mut line_y_vals: Vec<f64> = Vec::new();
            let mut last_y: f64 = 0.0;
            let mut last_dif: Option<f64> = None;
            let mut line_used_dif = false;

            for &tok in &tokens[1..] {
                match tok {
                    AsdfToken::Abs(v) => {
                        last_y = v;
                        last_dif = None;
                        line_y_vals.push(v);
                    }
                    AsdfToken::Dif(d) => {
                        line_used_dif = true;
                        last_dif = Some(d);
                        last_y += d;
                        line_y_vals.push(last_y);
                    }
                    AsdfToken::Dup(cnt) => {
                        // DUP count `cnt` means total occurrences = cnt, so append cnt - 1 more
                        for _ in 1..cnt {
                            if let Some(d) = last_dif {
                                last_y += d;
                                line_y_vals.push(last_y);
                            } else {
                                line_y_vals.push(last_y);
                            }
                        }
                    }
                }
            }

            // IUPAC DIF Checkpoint Rule: if previous line ended in DIF mode,
            // the first Y on this line duplicates the previous line's last Y -> skip it!
            let start_idx = if prev_line_was_dif && !line_y_vals.is_empty() {
                1
            } else {
                0
            };
            for &yv in &line_y_vals[start_idx..] {
                raw_ys.push(yv * yfactor);
            }
            prev_line_was_dif = line_used_dif;
        }

        if raw_ys.len() < 2 {
            return Err(anyhow!("JCAMP-DX file contained no decompressible XYDATA"));
        }

        let x0 = firstx.or(first_line_x).unwrap_or(0.0);
        let dx = if let Some(lx) = lastx {
            let denom = (raw_ys.len().max(2) - 1) as f64;
            (lx - x0) / denom
        } else {
            deltax.unwrap_or(1.0)
        };

        raw_ys
            .into_iter()
            .enumerate()
            .map(|(i, y)| (x0 + (i as f64) * dx, y))
            .collect()
    } else {
        return Err(anyhow!("No XYDATA or XYPOINTS block found in JCAMP-DX file"));
    };

    let (mut xs, mut ys) = finalize_points(pts)?;

    let xu = xunits.to_ascii_uppercase();
    let dt_up = data_type.to_ascii_uppercase();
    let is_uv_hint = dt_up.contains("UV")
        || dt_up.contains("ULTRAVIOLET")
        || kind_override == Some(SpectrumKind::UvVis);

    let is_micrometer = xu.contains("MICROMETER") || xu.contains("MICRON") || xu == "UM";
    let is_wavenumber = xu.contains("1/CM") || xu.contains("CM-1") || xu.contains("WAVENUMBER");
    let is_hz = xu == "HZ" || xu == "HERTZ" || xu.starts_with("HZ ");

    let mut need_resort = false;
    if is_micrometer {
        let x_max = xs.last().copied().unwrap_or(0.0);
        if is_uv_hint || x_max <= 2.0 {
            // UV-Vis / NIR wavelength in µm -> nm
            for v in &mut xs {
                *v *= 1000.0;
            }
        } else {
            // Infrared spectrum in µm (e.g. 2.5..15 µm) -> cm⁻¹ = 10_000 / µm
            let mut converted = Vec::with_capacity(xs.len());
            for (x, y) in xs.into_iter().zip(ys) {
                if x > 1e-6 {
                    converted.push((10_000.0 / x, y));
                }
            }
            let (nx, ny) = finalize_points(converted)?;
            xs = nx;
            ys = ny;
            need_resort = false;
        }
    } else if is_wavenumber && is_uv_hint {
        // UV-Vis stored in 1/cm -> nm = 1e7 / wn
        let mut converted = Vec::with_capacity(xs.len());
        for (x, y) in xs.into_iter().zip(ys) {
            if x > 1e-6 {
                converted.push((1.0e7 / x, y));
            }
        }
        let (nx, ny) = finalize_points(converted)?;
        xs = nx;
        ys = ny;
    } else if is_hz {
        if let Some(freq) = observe_freq_mhz.filter(|&f| f > 1.0) {
            for v in &mut xs {
                *v /= freq;
            }
            need_resort = true;
        }
    }

    if need_resort {
        let pairs: Vec<(f64, f64)> = xs.into_iter().zip(ys).collect();
        let (nx, ny) = finalize_points(pairs)?;
        xs = nx;
        ys = ny;
    }

    let mut hint = format!("{} {} {}", data_type, xunits, yunits);
    if is_hz && observe_freq_mhz.is_some() && dt_up.is_empty() {
        hint.push_str(" NMR ppm");
    }
    let kind = kind_override.unwrap_or_else(|| {
        infer_kind(xs[0], *xs.last().unwrap_or(&xs[0]), Some(&hint))
    });

    let raw_y_unit = infer_y_unit(kind, &ys, Some(&yunits));
    let yu = yunits.to_ascii_uppercase();
    let is_log_scale = yu.contains("LOG") || yu.contains("LOG10");

    if kind == SpectrumKind::UvVis && is_log_scale {
        // LOG EPSILON / LOG ABSORBANCE -> linearise then normalize to 0..1 Abs
        for v in &mut ys {
            *v = 10.0_f64.powf(*v);
        }
        let y_max = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if y_max > 1e-12 {
            for v in &mut ys {
                *v /= y_max;
            }
        }
    } else if (kind == SpectrumKind::UvVis || kind == SpectrumKind::Ftir)
        && raw_y_unit == YUnit::Transmittance
    {
        let y_max = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if y_max <= 1.05 {
            for v in &mut ys {
                *v *= 100.0;
            }
        }
    }

    Ok(Spectrum {
        name: title,
        source_path: None,
        kind,
        raw_y_unit,
        x: xs,
        y: ys,
        visible: true,
        color_idx: 0,
        metadata,
    })
}

// ============================================================================
// SAMPLE FILE GENERATOR (CSV + REAL DIFDUP JCAMP-DX)
// ============================================================================

fn encode_sqz(val: i64) -> String {
    let neg = val < 0;
    let digits = val.abs().to_string();
    let mut chs = digits.chars();
    let d0 = chs.next().unwrap_or('0') as u8 - b'0';
    let first = if d0 == 0 {
        '@'
    } else if neg {
        (b'a' + d0 - 1) as char
    } else {
        (b'A' + d0 - 1) as char
    };
    let mut out = String::new();
    out.push(first);
    out.extend(chs);
    out
}

fn encode_dif(val: i64) -> String {
    let neg = val < 0;
    let digits = val.abs().to_string();
    let mut chs = digits.chars();
    let d0 = chs.next().unwrap_or('0') as u8 - b'0';
    let first = if d0 == 0 {
        '%'
    } else if neg {
        (b'j' + d0 - 1) as char
    } else {
        (b'J' + d0 - 1) as char
    };
    let mut out = String::new();
    out.push(first);
    out.extend(chs);
    out
}

pub fn encode_jcamp_difdup(spec: &Spectrum) -> String {
    let yfactor = 0.01_f64;
    let n = spec.x.len();
    // Standard FTIR JCAMP stores 4000 -> 400 descending
    let xs: Vec<f64> = spec.x.iter().rev().copied().collect();
    let ys_int: Vec<i64> = spec
        .y
        .iter()
        .rev()
        .map(|&v| (v / yfactor).round() as i64)
        .collect();

    let mut out = String::new();
    out.push_str(&format!("##TITLE={}\n", spec.name));
    out.push_str("##JCAMP-DX=5.00 $$ Generated by SPECTRA\n");
    out.push_str("##DATA TYPE=INFRARED SPECTRUM\n");
    out.push_str("##XUNITS=1/CM\n");
    out.push_str("##YUNITS=TRANSMITTANCE\n");
    out.push_str(&format!("##FIRSTX={:.4}\n", xs[0]));
    out.push_str(&format!("##LASTX={:.4}\n", xs[n - 1]));
    out.push_str(&format!("##NPOINTS={}\n", n));
    out.push_str("##XFACTOR=1.0\n");
    out.push_str(&format!("##YFACTOR={}\n", yfactor));
    out.push_str("##XYDATA=(X++(Y..Y))\n");

    let per_line = 12;
    let mut idx = 0;
    while idx < n {
        let end = (idx + per_line).min(n);
        let mut line = format!("{:.2} ", xs[idx]);
        line.push_str(&encode_sqz(ys_int[idx]));
        let mut j = idx + 1;
        while j < end {
            let diff = ys_int[j] - ys_int[j - 1];
            let mut rep = 1usize;
            while j + rep < end && (ys_int[j + rep] - ys_int[j + rep - 1]) == diff && rep < 8 {
                rep += 1;
            }
            line.push_str(&encode_dif(diff));
            if rep > 1 {
                let dup_ch = (b'S' + (rep as u8) - 1) as char;
                line.push(dup_ch);
            }
            j += rep;
        }
        out.push_str(&line);
        out.push('\n');
        // Overlap by 1 point to emit standard IUPAC DIF checkpoint on next line!
        if end >= n {
            break;
        }
        idx = end - 1;
    }
    out.push_str("##END=\n");
    out
}

/// Write a processed trace as a clean CSV with a commented provenance header.
/// FTIR / NMR are written in the conventional descending instrument order;
/// UV-Vis stays ascending.
pub fn save_processed_csv(
    path: &Path,
    spec: &Spectrum,
    xs: &[f64],
    ys: &[f64],
    y_unit: YUnit,
    note: &str,
) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    let x_label = match spec.kind {
        SpectrumKind::Ftir => "Wavenumber (cm-1)",
        SpectrumKind::Hnmr | SpectrumKind::Cnmr => "Chemical Shift (ppm)",
        SpectrumKind::UvVis => "Wavelength (nm)",
    };
    let y_label = match y_unit {
        YUnit::Transmittance => "Transmittance (%T)",
        YUnit::Absorbance => "Absorbance (a.u.)",
        YUnit::Intensity => "Intensity (a.u.)",
    };

    let mut buf = String::new();
    buf.push_str("# SPECTRA processed spectrum\n");
    buf.push_str(&format!("# Trace: {}\n", spec.name));
    buf.push_str(&format!("# Technique: {}\n", spec.kind.display_name()));
    buf.push_str(&format!("# Y-axis: {}\n", y_unit.label()));
    for line in note.lines() {
        buf.push_str(&format!("# Processing: {}\n", line.trim()));
    }
    buf.push_str(&format!("{}\n", format!("{},{}", x_label, y_label)));

    let n = xs.len().min(ys.len());
    if spec.kind.is_x_reversed() {
        for i in (0..n).rev() {
            buf.push_str(&format!("{:.6},{:.6}\n", xs[i], ys[i]));
        }
    } else {
        for i in 0..n {
            buf.push_str(&format!("{:.6},{:.6}\n", xs[i], ys[i]));
        }
    }

    fs::write(path, buf).with_context(|| format!("Failed writing {}", path.display()))
}

/// Write a processed trace as JCAMP-DX. FTIR %T keeps the compact DIFDUP
/// encoder; everything else uses clean AFFN `##XYPOINTS=(XY..XY)` pairs.
pub fn save_processed_jcamp(
    path: &Path,
    name: &str,
    kind: SpectrumKind,
    xs: &[f64],
    ys: &[f64],
    y_unit: YUnit,
) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    let n = xs.len().min(ys.len());
    if n < 2 {
        return Err(anyhow!("Not enough points to write JCAMP-DX"));
    }

    let text = if kind == SpectrumKind::Ftir && y_unit == YUnit::Transmittance {
        let spec = Spectrum {
            name: name.to_string(),
            source_path: None,
            kind,
            raw_y_unit: y_unit,
            x: xs[..n].to_vec(),
            y: ys[..n].to_vec(),
            visible: true,
            color_idx: 0,
            metadata: Vec::new(),
        };
        encode_jcamp_difdup(&spec)
    } else {
        let (data_type, xunits, yunits) = match kind {
            SpectrumKind::Ftir => (
                "INFRARED SPECTRUM",
                "1/CM",
                match y_unit {
                    YUnit::Transmittance => "TRANSMITTANCE",
                    YUnit::Absorbance => "ABSORBANCE",
                    YUnit::Intensity => "ARBITRARY UNITS",
                },
            ),
            SpectrumKind::Hnmr | SpectrumKind::Cnmr => {
                ("NMR SPECTRUM", "PPM", "ARBITRARY UNITS")
            }
            SpectrumKind::UvVis => (
                "UV/VIS SPECTRUM",
                "NANOMETERS",
                match y_unit {
                    YUnit::Transmittance => "TRANSMITTANCE",
                    YUnit::Absorbance => "ABSORBANCE",
                    YUnit::Intensity => "ARBITRARY UNITS",
                },
            ),
        };

        let mut out = String::new();
        out.push_str(&format!("##TITLE={}\n", name));
        out.push_str("##JCAMP-DX=5.00 $$ Generated by SPECTRA\n");
        out.push_str(&format!("##DATA TYPE={}\n", data_type));
        out.push_str(&format!("##XUNITS={}\n", xunits));
        out.push_str(&format!("##YUNITS={}\n", yunits));
        let (first, last) = if kind.is_x_reversed() {
            (xs[n - 1], xs[0])
        } else {
            (xs[0], xs[n - 1])
        };
        out.push_str(&format!("##FIRSTX={:.6}\n", first));
        out.push_str(&format!("##LASTX={:.6}\n", last));
        out.push_str(&format!("##NPOINTS={}\n", n));
        out.push_str("##XYPOINTS=(XY..XY)\n");
        let mut line = String::new();
        for i in 0..n {
            let idx = if kind.is_x_reversed() { n - 1 - i } else { i };
            line.push_str(&format!("{:.6} {:.6}; ", xs[idx], ys[idx]));
            if line.len() > 96 {
                out.push_str(line.trim_end());
                out.push('\n');
                line.clear();
            }
        }
        if !line.trim().is_empty() {
            out.push_str(line.trim_end());
            out.push('\n');
        }
        out.push_str("##END=\n");
        out
    };

    fs::write(path, text).with_context(|| format!("Failed writing {}", path.display()))
}

fn write_csv(path: &Path, spec: &Spectrum, header: &str) -> Result<()> {
    let mut buf = String::new();
    buf.push_str(header);
    buf.push('\n');
    // Write descending (standard laboratory instrument export order)
    for (&x, &y) in spec.x.iter().rev().zip(spec.y.iter().rev()) {
        buf.push_str(&format!("{:.4},{:.4}\n", x, y));
    }
    fs::write(path, buf).with_context(|| format!("Failed writing {}", path.display()))
}

// ============================================================================
// SPECTRAL LIBRARY INDEX & NAME RESOLUTION (`spectra butane`, `/` search, paste)
// ============================================================================

/// One library file that can be resolved by name, e.g. `org_nist_n_butane.jdx`.
#[derive(Debug, Clone)]
pub struct LibraryEntry {
    pub path: PathBuf,
    pub stem: String,
    pub kind: SpectrumKind,
    pub score: i32,
}

pub fn is_spectrum_extension(path: &Path) -> bool {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    matches!(
        ext.as_str(),
        "csv" | "tsv" | "txt" | "dpt" | "asc" | "xy" | "jdx" | "dx" | "jcamp" | "spc" | "spa"
            | "opus"
    ) || (!ext.is_empty() && ext.chars().all(|c| c.is_ascii_digit()))
}

fn push_unique_dir(out: &mut Vec<PathBuf>, seen: &mut HashSet<PathBuf>, dir: PathBuf) {
    if let Ok(canon) = fs::canonicalize(&dir) {
        if canon.is_dir() && seen.insert(canon.clone()) {
            out.push(canon);
        }
    }
}

/// Directories scanned for reference spectra, resolved relative to the binary,
/// the manifest, `SPECTRA_HOME`, and the current working directory.
pub fn library_roots() -> Vec<PathBuf> {
    let mut bases: Vec<PathBuf> = Vec::new();
    if let Ok(home) = std::env::var("SPECTRA_HOME") {
        if !home.trim().is_empty() {
            bases.push(PathBuf::from(home));
        }
    }
    if let Ok(exe) = std::env::current_exe() {
        if let Some(root) = exe
            .parent()
            .and_then(|p| p.parent())
            .and_then(|p| p.parent())
        {
            bases.push(root.to_path_buf());
        }
    }
    bases.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")));
    if let Ok(cwd) = std::env::current_dir() {
        bases.push(cwd);
    }

    let mut roots = Vec::new();
    let mut seen = HashSet::new();
    for base in bases {
        push_unique_dir(&mut roots, &mut seen, base.join("library/ftir"));
        push_unique_dir(&mut roots, &mut seen, base.join("library/nmr"));
        push_unique_dir(&mut roots, &mut seen, base.join("library/uvvis"));
        push_unique_dir(&mut roots, &mut seen, base.join("samples"));
    }
    roots
}

fn infer_file_kind(path: &Path, stem: &str) -> SpectrumKind {
    let dir = path
        .parent()
        .and_then(|p| p.file_name())
        .and_then(|s| s.to_str())
        .unwrap_or("")
        .to_ascii_lowercase();
    let name = stem.to_ascii_lowercase();
    let uvvis = dir.contains("uvvis")
        || dir.contains("uv-vis")
        || name.contains("uvvis")
        || name.contains("uv-vis")
        || name.starts_with("pcad_")
        || name.contains("nist_uv_");
    let carbon = name.contains("13c") || name.contains("cnmr") || name.contains("_c13");
    let proton = name.contains("1h") || name.contains("hnmr") || name.contains("_h1");
    if uvvis {
        SpectrumKind::UvVis
    } else if carbon {
        SpectrumKind::Cnmr
    } else if proton || dir.contains("nmr") {
        SpectrumKind::Hnmr
    } else {
        SpectrumKind::Ftir
    }
}

fn normalize_name(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    let mut prev_space = true;
    for ch in text.chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch.to_ascii_lowercase());
            prev_space = false;
        } else if !prev_space {
            out.push(' ');
            prev_space = true;
        }
    }
    out.trim_end().to_string()
}

/// Score how well a library file stem matches a free-text query:
/// exact stem > exact token > all tokens > token prefixes > raw substring.
fn name_match_score(query_norm: &str, stem_norm: &str) -> Option<i32> {
    if query_norm.is_empty() || stem_norm.is_empty() {
        return None;
    }
    let q_tokens: Vec<&str> = query_norm.split(' ').collect();
    let s_tokens: Vec<&str> = stem_norm.split(' ').collect();
    let score = if stem_norm == query_norm {
        1000
    } else if s_tokens.contains(&query_norm) {
        900
    } else if q_tokens.iter().all(|qt| s_tokens.contains(qt)) {
        800
    } else if q_tokens
        .iter()
        .all(|qt| s_tokens.iter().any(|st| st.starts_with(qt)))
    {
        600
    } else if stem_norm.replace(' ', "").contains(&query_norm.replace(' ', "")) {
        400
    } else {
        return None;
    };
    Some(score - (s_tokens.len().min(40) as i32))
}

/// Every reference spectrum on disk, indexed lazily and refreshable.
/// Uses OnceLock<RwLock<Vec<LibraryEntry>>> so the index can be rebuilt
/// after downloading new reference packs without restarting the app.
fn library_index_lock() -> &'static RwLock<Vec<LibraryEntry>> {
    static INDEX: OnceLock<RwLock<Vec<LibraryEntry>>> = OnceLock::new();
    INDEX.get_or_init(|| RwLock::new(Vec::new()))
}

/// Build the library index from disk (scans all library roots).
fn build_library_index() -> Vec<LibraryEntry> {
    let mut out: Vec<LibraryEntry> = Vec::new();
    let mut seen: HashSet<PathBuf> = HashSet::new();
    for root in library_roots() {
        let Ok(rd) = fs::read_dir(&root) else {
            continue;
        };
        for item in rd.flatten() {
            let path = item.path();
            if !is_spectrum_extension(&path) {
                continue;
            }
            let is_file = match item.file_type() {
                Ok(ft) if ft.is_symlink() => path.is_file(),
                Ok(ft) => ft.is_file(),
                Err(_) => false,
            };
            if !is_file || !seen.insert(path.clone()) {
                continue;
            }
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("")
                .to_string();
            if stem.is_empty() {
                continue;
            }
            out.push(LibraryEntry {
                kind: infer_file_kind(&path, &stem),
                path,
                stem,
                score: 0,
            });
        }
    }
    out
}

/// Get the current library index, initializing it on first access.
pub fn library_index() -> Vec<LibraryEntry> {
    let lock = library_index_lock();
    let guard = lock.read().unwrap_or_else(|e| e.into_inner());
    if guard.is_empty() {
        drop(guard);
        let mut write_guard = lock.write().unwrap_or_else(|e| e.into_inner());
        if write_guard.is_empty() {
            *write_guard = build_library_index();
        }
    }
    lock.read().unwrap_or_else(|e| e.into_inner()).clone()
}

/// Force a re-scan of all library directories.
/// Call this after downloading new reference packs so search/matcher see them immediately.
pub fn refresh_library_index() {
    let lock = library_index_lock();
    let mut guard = lock.write().unwrap_or_else(|e| e.into_inner());
    *guard = build_library_index();
    drop(guard);
    crate::matcher::invalidate_matcher_cache();
}

/// Rank every reference spectrum whose file name matches `query`.
pub fn search_library(query: &str) -> Vec<LibraryEntry> {
    let q = normalize_name(query);
    if q.is_empty() {
        return Vec::new();
    }
    let mut out: Vec<LibraryEntry> = library_index()
        .iter()
        .filter_map(|entry| {
            let stem_norm = normalize_name(&entry.stem);
            name_match_score(&q, &stem_norm).map(|score| LibraryEntry {
                score,
                ..entry.clone()
            })
        })
        .collect();
    let nucleus_hint = q.contains("13c")
        || q.contains("1h")
        || q.contains("cnmr")
        || q.contains("hnmr")
        || q.contains("nmr");
    let kind_rank = |k: SpectrumKind| match k {
        SpectrumKind::Ftir => 0u8,
        SpectrumKind::Hnmr => 1,
        SpectrumKind::Cnmr => 2,
        SpectrumKind::UvVis => 3,
    };
    out.sort_by(|a, b| {
        // A bare name like "butane" should surface the IR reference first;
        // an explicit nucleus ("butane 13c", "1h nmr") overrides that.
        let ra = if nucleus_hint { 0 } else { kind_rank(a.kind) };
        let rb = if nucleus_hint { 0 } else { kind_rank(b.kind) };
        ra.cmp(&rb)
            .then_with(|| b.score.cmp(&a.score))
            .then_with(|| a.stem.len().cmp(&b.stem.len()))
            .then_with(|| a.path.cmp(&b.path))
    });
    out
}

/// Rank every reference spectrum whose file name matches `query`, optionally
/// restricted to a single spectroscopy mode. Restriction is strict: asking for
/// ¹H NMR never returns an FTIR file just because its name matches.
pub fn search_library_for_kind(
    query: &str,
    kind: Option<SpectrumKind>,
) -> Vec<LibraryEntry> {
    let mut hits = search_library(query);
    if let Some(kind) = kind {
        hits.retain(|entry| entry.kind == kind);
    }
    hits
}

/// Best single library file for a name, honouring an optional nucleus hint.
pub fn resolve_library_query(query: &str, prefer: Option<SpectrumKind>) -> Option<LibraryEntry> {
    search_library_for_kind(query, prefer).into_iter().next()
}

/// True when the query names a file outright rather than merely describing it:
/// `n_butane` → `org_nist_n_butane`, but `butane` → `org_nist_n_butane` is not
/// exact (the stem carries a chemical prefix). Used so an unambiguous request
/// opens immediately while a vague one goes to the candidate picker.
pub fn exact_stem_match(stem: &str, query: &str) -> bool {
    const VENDOR: [&str; 8] = [
        "org", "nist", "mnhn", "jcamp", "atr", "ftir", "ir", "sample",
    ];
    let s = normalize_name(stem);
    let q = normalize_name(query);
    if s.is_empty() || q.is_empty() {
        return false;
    }
    if s == q {
        return true;
    }
    match s.strip_suffix(&q) {
        Some(prefix) if prefix.ends_with(' ') => prefix
            .trim()
            .split(' ')
            .filter(|t| !t.is_empty())
            .all(|t| VENDOR.contains(&t)),
        _ => false,
    }
}

/// Expand `~/...` so pasted and typed paths behave like a shell would.
pub fn expand_tilde(raw: &str) -> PathBuf {
    let trimmed = raw.trim().trim_matches('"').trim_matches('\'').trim();
    if trimmed == "~" {
        if let Ok(home) = std::env::var("HOME") {
            return PathBuf::from(home);
        }
    }
    if let Some(rest) = trimmed.strip_prefix("~/") {
        if let Ok(home) = std::env::var("HOME") {
            return PathBuf::from(home).join(rest);
        }
    }
    PathBuf::from(trimmed)
}

/// Base directory for small UI state files (last browsed folder).
/// Follows the XDG state convention, with a temp fallback for odd setups.
fn state_base() -> PathBuf {
    if let Some(dir) = std::env::var_os("XDG_STATE_HOME") {
        if !dir.is_empty() {
            return PathBuf::from(dir).join("spectra");
        }
    }
    if let Some(home) = std::env::var_os("HOME") {
        if !home.is_empty() {
            return PathBuf::from(home).join(".local").join("state").join("spectra");
        }
    }
    std::env::temp_dir().join("spectra-state")
}

/// Read a remembered folder from an explicit state base (testable without
/// touching the real user profile).
pub(crate) fn read_last_dir_from(base: &Path) -> Option<PathBuf> {
    let text = fs::read_to_string(base.join("last_dir")).ok()?;
    let trimmed = text.trim();
    if trimmed.is_empty() {
        return None;
    }
    let path = PathBuf::from(trimmed);
    if path.is_dir() {
        Some(path)
    } else {
        None
    }
}

pub(crate) fn write_last_dir_to(base: &Path, dir: &Path) -> std::io::Result<()> {
    fs::create_dir_all(base)?;
    fs::write(base.join("last_dir"), format!("{}\n", dir.display()))
}

/// Folder the interactive file picker last visited, if it still exists.
pub fn remembered_picker_dir() -> Option<PathBuf> {
    read_last_dir_from(&state_base())
}

/// Persist the picker folder. Failures (read-only $HOME, full disk) are
/// deliberately ignored: remembering a directory is a convenience, never a
/// reason to interrupt the user.
pub fn remember_picker_dir(dir: &Path) {
    let _ = write_last_dir_to(&state_base(), dir);
}

/// Turn pasted / dropped text into concrete paths. Accepts `file://` URIs,
/// shell-escaped names, and one path per line (kitty drag-and-drop format).
pub fn parse_pasted_paths(text: &str) -> Vec<PathBuf> {
    let mut out = Vec::new();
    for line in text.lines() {
        let raw = line.trim();
        if raw.is_empty() {
            continue;
        }
        if let Some(rest) = raw.strip_prefix("file://") {
            // `file:///abs/path` keeps the slash; `file://host/abs/path`
            // carries an authority component that must be dropped.
            let path = if let Some(idx) = rest.find('/') {
                let host = &rest[..idx];
                if host.is_empty() {
                    rest
                } else {
                    &rest[idx..]
                }
            } else {
                rest
            };
            out.push(expand_tilde(path));
            continue;
        }
        // Drag-and-drop from some file managers encodes spaces as %20.
        let unescaped = if raw.contains("%20") {
            raw.replace("%20", " ")
        } else {
            raw.to_string()
        };
        out.push(expand_tilde(&unescaped));
    }
    out
}

pub fn write_sample_files(dir: &Path) -> Result<Vec<PathBuf>> {
    fs::create_dir_all(dir)?;
    let mut written = Vec::new();

    let p1 = dir.join("mma_monomer_ftir.csv");
    write_csv(&p1, &mma_monomer_spectrum(), "Wavenumber (cm-1),Transmittance (%T)")?;
    written.push(p1);

    let p2 = dir.join("pmma_polymer_ftir.csv");
    write_csv(&p2, &pmma_polymer_spectrum(), "Wavenumber (cm-1),Transmittance (%T)")?;
    written.push(p2);

    let p3 = dir.join("polystyrene_ftir.jdx");
    let jdx = encode_jcamp_difdup(&polystyrene_ftir_spectrum());
    fs::write(&p3, jdx)?;
    written.push(p3);

    let p4 = dir.join("pmma_1h_nmr_cdcl3.csv");
    write_csv(&p4, &pmma_hnmr_spectrum(), "Chemical Shift (ppm),Intensity")?;
    written.push(p4);

    let p5 = dir.join("tio2_zno_uvvis.csv");
    write_csv(&p5, &tio2_zno_uvvis_spectrum(), "Wavelength (nm),Absorbance")?;
    written.push(p5);

    Ok(written)
}
