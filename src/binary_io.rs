//! Zero-dependency native parsers for the three dominant laboratory binary
//! spectroscopy formats:
//!
//!   * Galactic / Thermo `.spc` (SPCHDR / SUBHDR)
//!   * Bruker OPUS `.0`, `.1`, `.0000`, `.opus` (block directory + parameter tags)
//!   * Thermo Nicolet OMNIC `.spa` (block directory with key 2 / key 3)
//!
//! Every parser returns a strictly ascending `Spectrum` so the rest of the
//! application never has to care where the numbers came from.

use anyhow::{bail, Result};

use crate::data::{Spectrum, SpectrumKind, YUnit};

/// Magic-byte dispatcher. Returns `None` when the buffer is not a known binary
/// instrument format, so `src/io.rs` can fall back to JCAMP-DX / delimited text.
pub fn try_parse_binary(
    bytes: &[u8],
    default_name: &str,
    ext_lower: &str,
    kind_override: Option<SpectrumKind>,
) -> Option<Result<Spectrum>> {
    if bytes.len() >= 24 && &bytes[0..4] == b"\x0A\x0A\xFE\xFE" {
        Some(parse_opus_bytes(bytes, default_name, kind_override))
    } else if bytes.len() >= 320 && bytes.starts_with(b"Spectral Data File") {
        Some(parse_spa_bytes(bytes, default_name, kind_override))
    } else if bytes.len() >= 544
        && (ext_lower == "spc"
            || (bytes[1] == 0x4B
                && u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]) >= 2))
    {
        Some(parse_spc_bytes(bytes, default_name, kind_override))
    } else {
        None
    }
}

fn c_string(raw: &[u8]) -> String {
    let end = raw.iter().position(|&b| b == 0).unwrap_or(raw.len());
    let text: String = raw[..end]
        .iter()
        .filter(|&&b| (0x20..=0x7E).contains(&b))
        .map(|&b| b as char)
        .collect();
    text.trim().to_string()
}

fn linear_x(first: f64, last: f64, n: usize) -> Vec<f64> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![first];
    }
    let step = (last - first) / (n as f64 - 1.0);
    (0..n).map(|i| first + step * i as f64).collect()
}

/// Sort ascending and average duplicate X values. Binary formats occasionally
/// carry repeated abscissae after unit conversion.
fn finish_pairs(mut pairs: Vec<(f64, f64)>) -> Result<(Vec<f64>, Vec<f64>)> {
    pairs.retain(|(x, y)| x.is_finite() && y.is_finite());
    if pairs.len() < 2 {
        bail!("Binary spectrum contained fewer than two finite data points");
    }
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));

    let mut xs = Vec::with_capacity(pairs.len());
    let mut ys = Vec::with_capacity(pairs.len());
    let mut i = 0;
    while i < pairs.len() {
        let x0 = pairs[i].0;
        let mut sum = pairs[i].1;
        let mut cnt = 1.0;
        while i + 1 < pairs.len() && (pairs[i + 1].0 - x0).abs() < 1e-12 {
            i += 1;
            sum += pairs[i].1;
            cnt += 1.0;
        }
        xs.push(x0);
        ys.push(sum / cnt);
        i += 1;
    }
    Ok((xs, ys))
}

fn infer_kind_range(x_min: f64, x_max: f64) -> SpectrumKind {
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

fn infer_y_unit(kind: SpectrumKind, ys: &[f64]) -> YUnit {
    match kind {
        SpectrumKind::UvVis => YUnit::Absorbance,
        SpectrumKind::Hnmr | SpectrumKind::Cnmr => YUnit::Intensity,
        SpectrumKind::Ftir => {
            let y_max = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let y_min = ys.iter().copied().fold(f64::INFINITY, f64::min);
            if y_max > 10.0 && y_max <= 115.0 && y_min >= -5.0 {
                YUnit::Transmittance
            } else {
                YUnit::Absorbance
            }
        }
    }
}

// ============================================================================
// Galactic / Thermo SPC
// ============================================================================

pub fn parse_spc_bytes(
    bytes: &[u8],
    default_name: &str,
    kind_override: Option<SpectrumKind>,
) -> Result<Spectrum> {
    if bytes.len() < 544 {
        bail!("SPC file too short ({} bytes, need at least 544)", bytes.len());
    }
    if bytes[1] != 0x4B {
        bail!("Unsupported SPC byte order (only new little-endian SPC files are supported)");
    }

    let ftflgs = bytes[0];
    let tsprec = (ftflgs & 0x01) != 0;
    let txvals = (ftflgs & 0x80) != 0;
    let fexper = bytes[2];
    let fexp = bytes[3] as i8;
    let fnpts = u32::from_le_bytes([bytes[4], bytes[5], bytes[6], bytes[7]]) as usize;
    let ffirst = f64::from_le_bytes(bytes[8..16].try_into().unwrap());
    let flast = f64::from_le_bytes(bytes[16..24].try_into().unwrap());
    let fxtype = bytes[28];
    let fytype = bytes[29];

    let memo = c_string(&bytes[48..178]);
    let title = if memo.is_empty() {
        default_name.to_string()
    } else {
        memo
    };

    let explicit_x: Option<Vec<f64>> = if txvals {
        let end = 512usize.saturating_add(fnpts.saturating_mul(4));
        if fnpts >= 2 && end <= bytes.len() {
            Some(
                (0..fnpts)
                    .map(|i| {
                        let o = 512 + i * 4;
                        f32::from_le_bytes(bytes[o..o + 4].try_into().unwrap()) as f64
                    })
                    .collect(),
            )
        } else {
            None
        }
    } else {
        None
    };

    let sub_off = if txvals {
        512usize.saturating_add(fnpts.saturating_mul(4))
    } else {
        512
    };

    let (subexp, subnpts) = if sub_off + 32 <= bytes.len() {
        (
            bytes[sub_off + 1] as i8,
            u32::from_le_bytes([
                bytes[sub_off + 16],
                bytes[sub_off + 17],
                bytes[sub_off + 18],
                bytes[sub_off + 19],
            ]) as usize,
        )
    } else {
        (0, 0)
    };

    let npts = if (2..=5_000_000).contains(&subnpts) {
        subnpts
    } else {
        fnpts
    };
    let exp = if subexp != 0 { subexp } else { fexp };

    let y_off = sub_off + 32;
    if y_off >= bytes.len() {
        bail!("SPC data section starts beyond end of file");
    }
    let bytes_per = if exp == -128 {
        4
    } else if tsprec {
        2
    } else {
        4
    };
    let available = (bytes.len() - y_off) / bytes_per;
    let n = npts.min(available);
    if n < 2 {
        bail!("SPC file contains too few Y values");
    }

    let scale = 2.0_f64.powi(exp as i32);
    let mut ys: Vec<f64> = Vec::with_capacity(n);
    for i in 0..n {
        let o = y_off + i * bytes_per;
        let v = if exp == -128 {
            f32::from_le_bytes(bytes[o..o + 4].try_into().unwrap()) as f64
        } else if tsprec {
            let raw = i16::from_le_bytes([bytes[o], bytes[o + 1]]) as f64;
            raw * scale / 65536.0
        } else {
            let raw = i32::from_le_bytes([
                bytes[o],
                bytes[o + 1],
                bytes[o + 2],
                bytes[o + 3],
            ]) as f64;
            raw * scale / 4_294_967_296.0
        };
        ys.push(v);
    }

    let mut xs = match explicit_x {
        Some(mut xv) if xv.len() >= n => {
            xv.truncate(n);
            xv
        }
        _ => linear_x(ffirst, flast, n),
    };

    // Micrometres: IR -> cm⁻¹, UV-Vis -> nm.
    if fxtype == 2 {
        let uv = fexper == 10 || kind_override == Some(SpectrumKind::UvVis);
        if uv {
            for x in &mut xs {
                *x *= 1000.0;
            }
        } else {
            for x in &mut xs {
                if *x > 0.0 {
                    *x = 10_000.0 / *x;
                } else {
                    *x = f64::NAN;
                }
            }
        }
    }

    let x_min = xs.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let kind = kind_override.unwrap_or_else(|| {
        if fxtype == 3 || fexper == 10 {
            SpectrumKind::UvVis
        } else if fxtype == 1 || fxtype == 13 || fexper == 4 {
            SpectrumKind::Ftir
        } else if fxtype == 7 || fexper == 12 {
            if x_max > 25.0 {
                SpectrumKind::Cnmr
            } else {
                SpectrumKind::Hnmr
            }
        } else {
            infer_kind_range(x_min, x_max)
        }
    });

    let y_unit = if fytype == 2 || fytype == 3 {
        YUnit::Absorbance
    } else if fytype == 7 {
        for v in &mut ys {
            *v *= 100.0;
        }
        YUnit::Transmittance
    } else if fytype == 10 {
        YUnit::Transmittance
    } else {
        infer_y_unit(kind, &ys)
    };

    let pairs: Vec<(f64, f64)> = xs.into_iter().zip(ys).collect();
    let (xs, ys) = finish_pairs(pairs)?;

    Ok(Spectrum {
        name: title,
        source_path: None,
        kind,
        raw_y_unit: y_unit,
        x: xs,
        y: ys,
        visible: true,
        color_idx: 0,
        metadata: vec![("Format".into(), "Galactic SPC (binary)".into())],
    })
}

// ============================================================================
// Bruker OPUS
// ============================================================================

#[derive(Debug, Clone, Default)]
struct OpusParams {
    npf: Option<usize>,
    fxv: Option<f64>,
    lxv: Option<f64>,
    csf: f64,
    dxu: Option<String>,
}

fn parse_opus_param_block(slice: &[u8]) -> Option<OpusParams> {
    let mut params = OpusParams {
        csf: 1.0,
        ..OpusParams::default()
    };
    let mut pos = 0usize;
    while pos + 8 <= slice.len() {
        let tag = &slice[pos..pos + 4];
        if tag.starts_with(b"END") || tag[0] == 0 {
            break;
        }
        let _ptype = u16::from_le_bytes([slice[pos + 4], slice[pos + 5]]);
        let val_bytes = u16::from_le_bytes([slice[pos + 6], slice[pos + 7]]) as usize * 2;
        let val_start = pos + 8;
        let val_end = val_start.saturating_add(val_bytes);
        if val_end > slice.len() {
            break;
        }
        let val = &slice[val_start..val_end];

        match &tag[0..3] {
            b"NPF" if val_bytes >= 4 => {
                let raw = i32::from_le_bytes(val[0..4].try_into().unwrap());
                if raw >= 0 {
                    params.npf = Some(raw as usize);
                }
            }
            b"FXV" if val_bytes >= 8 => {
                params.fxv = Some(f64::from_le_bytes(val[0..8].try_into().unwrap()));
            }
            b"LXV" if val_bytes >= 8 => {
                params.lxv = Some(f64::from_le_bytes(val[0..8].try_into().unwrap()));
            }
            b"CSF" if val_bytes >= 8 => {
                params.csf = f64::from_le_bytes(val[0..8].try_into().unwrap());
            }
            b"DXU" => {
                params.dxu = Some(c_string(val));
            }
            _ => {}
        }
        pos = val_end;
    }

    match (params.npf, params.fxv, params.lxv) {
        (Some(npf), Some(_), Some(_)) if npf >= 2 => Some(params),
        _ => None,
    }
}

fn opus_block_name(b0: u8, b1: u8) -> &'static str {
    match (b0, b1) {
        (0x0F, 0x10) => "AB (Absorbance)",
        (0x0F, 0x14) => "TR (Transmittance)",
        (0x0F, 0x1C) => "KM (Kubelka-Munk)",
        (0x0F, 0x0C) => "RFC (Reflectance)",
        (0x07, 0x04) => "ScSm (Single-channel sample)",
        _ => "Spectrum",
    }
}

fn opus_rank(b0: u8, b1: u8) -> u8 {
    match (b0, b1) {
        (0x0F, 0x10) => 7,
        (0x0F, 0x14) => 6,
        (0x0F, 0x1C) => 5,
        (0x0F, 0x0C) => 4,
        (0x0F, _) => 3,
        (0x07, 0x04) => 2,
        (0x07, _) => 1,
        _ => 0,
    }
}

pub fn parse_opus_bytes(
    bytes: &[u8],
    default_name: &str,
    kind_override: Option<SpectrumKind>,
) -> Result<Spectrum> {
    if bytes.len() < 24 || &bytes[0..4] != b"\x0A\x0A\xFE\xFE" {
        bail!("Not a Bruker OPUS file (bad magic)");
    }
    let dir_offset = u32::from_le_bytes([bytes[12], bytes[13], bytes[14], bytes[15]]) as usize;
    let ndir = u32::from_le_bytes([bytes[20], bytes[21], bytes[22], bytes[23]]) as usize;
    if dir_offset >= bytes.len() {
        bail!("OPUS directory offset is beyond end of file");
    }
    let max_entries = (bytes.len() - dir_offset) / 12;
    let n_entries = ndir.min(max_entries);

    let mut params: Vec<((u8, u8), OpusParams)> = Vec::new();
    let mut data_blocks: Vec<(u8, u8, usize, usize)> = Vec::new();

    for i in 0..n_entries {
        let o = dir_offset + i * 12;
        let block_type = u32::from_le_bytes([
            bytes[o],
            bytes[o + 1],
            bytes[o + 2],
            bytes[o + 3],
        ]);
        let len_words = u32::from_le_bytes([
            bytes[o + 4],
            bytes[o + 5],
            bytes[o + 6],
            bytes[o + 7],
        ]) as usize;
        let offset = u32::from_le_bytes([
            bytes[o + 8],
            bytes[o + 9],
            bytes[o + 10],
            bytes[o + 11],
        ]) as usize;
        let b0 = (block_type & 0xFF) as u8;
        let b1 = ((block_type >> 8) & 0xFF) as u8;
        let byte_len = len_words.saturating_mul(4);
        if offset >= bytes.len() || byte_len == 0 {
            continue;
        }
        let end = offset.saturating_add(byte_len).min(bytes.len());
        let slice = &bytes[offset..end];

        if b0 == 0x0F || b0 == 0x07 {
            if byte_len >= 8 {
                data_blocks.push((b0, b1, byte_len, offset));
            }
        } else if let Some(p) = parse_opus_param_block(slice) {
            params.push(((b0, b1), p));
        }
    }

    if data_blocks.is_empty() {
        bail!("OPUS file contains no spectrum data block");
    }

    // Prefer real spectra over interferograms when both are present.
    let filtered: Vec<(u8, u8, usize, usize)> = data_blocks
        .iter()
        .copied()
        .filter(|&(b0, b1, _, _)| !(b0 == 0x07 && (b1 == 0x08 || b1 == 0x0C)))
        .collect();
    let candidates = if filtered.is_empty() {
        data_blocks
    } else {
        filtered
    };

    let mut best = candidates[0];
    for &cand in &candidates[1..] {
        if opus_rank(cand.0, cand.1) > opus_rank(best.0, best.1) {
            best = cand;
        }
    }
    let (b0, b1, byte_len, offset) = best;

    let params_ref = params
        .iter()
        .find(|((p0, p1), _)| *p0 == b0.wrapping_add(16) && *p1 == b1)
        .map(|(_, p)| p)
        .or_else(|| params.iter().find(|((_, p1), _)| *p1 == b1).map(|(_, p)| p))
        .or_else(|| {
            params
                .iter()
                .find(|(_, p)| p.npf.map(|n| n * 4 <= byte_len).unwrap_or(false))
                .map(|(_, p)| p)
        })
        .ok_or_else(|| anyhow::anyhow!("OPUS spectrum block has no matching parameter block"))?;

    let npf = params_ref.npf.unwrap_or(0);
    let fxv = params_ref.fxv.unwrap_or(0.0);
    let lxv = params_ref.lxv.unwrap_or(0.0);
    let csf = if params_ref.csf.is_finite() && params_ref.csf != 0.0 {
        params_ref.csf
    } else {
        1.0
    };

    let n = npf.min(byte_len / 4);
    if n < 2 {
        bail!("OPUS parameter block reports fewer than two points");
    }
    let mut ys: Vec<f64> = Vec::with_capacity(n);
    for i in 0..n {
        let o = offset + i * 4;
        let v = f32::from_le_bytes(bytes[o..o + 4].try_into().unwrap()) as f64;
        ys.push(v * csf);
    }
    let xs = linear_x(fxv, lxv, n);

    let x_min = xs.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let kind = kind_override.unwrap_or_else(|| infer_kind_range(x_min, x_max));

    let y_unit = match (b0, b1) {
        (0x0F, 0x14) | (0x0F, 0x0C) => YUnit::Transmittance,
        (0x0F, 0x10) | (0x0F, 0x1C) => YUnit::Absorbance,
        _ => infer_y_unit(kind, &ys),
    };

    if y_unit == YUnit::Transmittance {
        let y_max = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if y_max <= 1.05 {
            for v in &mut ys {
                *v *= 100.0;
            }
        }
    }

    let pairs: Vec<(f64, f64)> = xs.into_iter().zip(ys).collect();
    let (xs, ys) = finish_pairs(pairs)?;

    let mut metadata = vec![("Format".into(), "Bruker OPUS (binary)".into())];
    metadata.push(("OPUS Block".into(), opus_block_name(b0, b1).into()));
    if let Some(dxu) = &params_ref.dxu {
        if !dxu.is_empty() {
            metadata.push(("OPUS X Unit".into(), dxu.clone()));
        }
    }

    Ok(Spectrum {
        name: default_name.to_string(),
        source_path: None,
        kind,
        raw_y_unit: y_unit,
        x: xs,
        y: ys,
        visible: true,
        color_idx: 0,
        metadata,
    })
}

// ============================================================================
// Thermo Nicolet OMNIC SPA
// ============================================================================

pub fn parse_spa_bytes(
    bytes: &[u8],
    default_name: &str,
    kind_override: Option<SpectrumKind>,
) -> Result<Spectrum> {
    if bytes.len() < 320 || !bytes.starts_with(b"Spectral Data File") {
        bail!("Not a Thermo Nicolet SPA file (bad signature)");
    }

    let title_raw = c_string(&bytes[30..286]);
    let title = if title_raw.is_empty() {
        default_name.to_string()
    } else {
        title_raw
    };

    let n_blocks = u16::from_le_bytes([bytes[294], bytes[295]]) as usize;
    let max_blocks = (bytes.len() - 304) / 16;
    let n_blocks = n_blocks.min(max_blocks);

    let mut header: Option<(usize, u32, u32, f64, f64)> = None;
    let mut y_data: Option<Vec<f64>> = None;

    for i in 0..n_blocks {
        let o = 304 + i * 16;
        if o + 16 > bytes.len() {
            break;
        }
        let key = u16::from_le_bytes([bytes[o], bytes[o + 1]]);
        let off = u32::from_le_bytes([bytes[o + 2], bytes[o + 3], bytes[o + 4], bytes[o + 5]]) as usize;
        let size = u32::from_le_bytes([bytes[o + 6], bytes[o + 7], bytes[o + 8], bytes[o + 9]]) as usize;
        if off >= bytes.len() {
            continue;
        }
        let end = off.saturating_add(size).min(bytes.len());

        match key {
            2 if size >= 24 && off + 24 <= bytes.len() => {
                let npts = u32::from_le_bytes([
                    bytes[off + 4],
                    bytes[off + 5],
                    bytes[off + 6],
                    bytes[off + 7],
                ]) as usize;
                let x_unit = u32::from_le_bytes([
                    bytes[off + 8],
                    bytes[off + 9],
                    bytes[off + 10],
                    bytes[off + 11],
                ]);
                let y_unit = u32::from_le_bytes([
                    bytes[off + 12],
                    bytes[off + 13],
                    bytes[off + 14],
                    bytes[off + 15],
                ]);
                let x_first = f32::from_le_bytes([
                    bytes[off + 16],
                    bytes[off + 17],
                    bytes[off + 18],
                    bytes[off + 19],
                ]) as f64;
                let x_last = f32::from_le_bytes([
                    bytes[off + 20],
                    bytes[off + 21],
                    bytes[off + 22],
                    bytes[off + 23],
                ]) as f64;
                header = Some((npts, x_unit, y_unit, x_first, x_last));
            }
            3 => {
                let count = size / 4;
                let mut vals = Vec::with_capacity(count);
                for j in 0..count {
                    let p = off + j * 4;
                    if p + 4 > end {
                        break;
                    }
                    vals.push(f32::from_le_bytes(bytes[p..p + 4].try_into().unwrap()) as f64);
                }
                if !vals.is_empty() {
                    y_data = Some(vals);
                }
            }
            _ => {}
        }
    }

    let (npts, x_unit, y_unit, x_first, x_last) =
        header.ok_or_else(|| anyhow::anyhow!("SPA file is missing its key-2 header block"))?;
    let mut ys = y_data.ok_or_else(|| anyhow::anyhow!("SPA file is missing its key-3 data block"))?;

    let n = if npts >= 2 {
        npts.min(ys.len())
    } else {
        ys.len()
    };
    if n < 2 {
        bail!("SPA file contains too few Y values");
    }
    ys.truncate(n);
    let xs = linear_x(x_first, x_last, n);

    let x_min = xs.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let kind = kind_override.unwrap_or_else(|| match x_unit {
        3 => SpectrumKind::UvVis,
        1 => SpectrumKind::Ftir,
        _ => infer_kind_range(x_min, x_max),
    });

    let y_unit = match y_unit {
        2 | 17 => YUnit::Transmittance,
        1 | 15 | 16 => YUnit::Absorbance,
        _ => infer_y_unit(kind, &ys),
    };

    let pairs: Vec<(f64, f64)> = xs.into_iter().zip(ys).collect();
    let (xs, ys) = finish_pairs(pairs)?;

    Ok(Spectrum {
        name: title,
        source_path: None,
        kind,
        raw_y_unit: y_unit,
        x: xs,
        y: ys,
        visible: true,
        color_idx: 0,
        metadata: vec![("Format".into(), "Thermo Nicolet SPA (binary)".into())],
    })
}
