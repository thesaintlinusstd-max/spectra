use std::collections::hash_map::DefaultHasher;
use std::env;
use std::fs;
use std::hash::{Hash, Hasher};
use std::io::Write;
use std::path::PathBuf;

use anyhow::Result;

use crate::app::{App, Modal};
use crate::data::{SpectrumKind, YUnit};
use crate::processing::interpolate_y;

const IMAGE_ID: u32 = 42;

// Scientific pitch-black oscilloscope palette: white-mint primary curve!
const BG_RGBA: [u8; 4] = [8, 12, 10, 255];
const GRID_RGB: (u8, u8, u8) = (26, 46, 35);
const CROSSHAIR_RGB: (u8, u8, u8) = (232, 206, 118);
const SOLVENT_RGB: (u8, u8, u8) = (214, 102, 102);
const DRAG_RGB: (u8, u8, u8) = (78, 168, 144);
const PEAK_PIN_RGB: (u8, u8, u8) = (242, 214, 114);
const PEAK_AUTO_RGB: (u8, u8, u8) = (130, 178, 148);
const PEAK_LABEL_AUTO_RGB: (u8, u8, u8) = (186, 224, 198);

// Crisp bitmap label metrics: 5×7 glyph scaled 2× → 10×14 px, 12 px advance.
const LABEL_SCALE: i32 = 2;
const GLYPH_W: i32 = 5 * LABEL_SCALE;
const GLYPH_H: i32 = 7 * LABEL_SCALE;
const LABEL_ADVANCE: i32 = 12;
const LABEL_PAD: i32 = 4;


const HD_TRACE_PALETTE: [(u8, u8, u8); 6] = [
    (226, 252, 234), // Trace 0: Luminous white-tinted mint green!
    (238, 184, 98),  // Trace 1: Warm amber
    (110, 214, 182), // Trace 2: Crisp aqua
    (194, 146, 204), // Trace 3: Soft purple
    (226, 118, 118), // Trace 4: Coral
    (236, 212, 124), // Trace 5: Yellow
];

pub fn is_kitty_terminal() -> bool {
    env::var("KITTY_WINDOW_ID").is_ok()
        || env::var("TERM")
            .map(|t| t.to_ascii_lowercase().contains("kitty"))
            .unwrap_or(false)
}

/// Approximate the physical colour of a visible wavelength (Dan Bruton's
/// piecewise-linear CIE fit, gamma 0.8). UV (<380 nm) is rendered as a dim
/// violet and NIR (>750 nm) as a dim deep red so the axis stays readable.
pub fn wavelength_to_rgb(nm: f64) -> (u8, u8, u8) {
    if !nm.is_finite() {
        return (0, 0, 0);
    }
    // Outside the visible window: dim, monotonically darker tails.
    if nm < 380.0 {
        let t = ((nm - 300.0) / 80.0).clamp(0.0, 1.0);
        let k = 0.18 + 0.42 * t;
        return (
            (120.0 * k).round() as u8,
            (36.0 * k).round() as u8,
            (190.0 * k).round() as u8,
        );
    }
    if nm > 750.0 {
        let t = ((nm - 750.0) / 150.0).clamp(0.0, 1.0);
        let k = 0.85 - 0.55 * t;
        return (
            (190.0 * k).round() as u8,
            (18.0 * k).round() as u8,
            (14.0 * k).round() as u8,
        );
    }

    let w = nm;
    let (r, g, b) = if w < 440.0 {
        (-(w - 440.0) / (440.0 - 380.0), 0.0, 1.0)
    } else if w < 490.0 {
        (0.0, (w - 440.0) / (490.0 - 440.0), 1.0)
    } else if w < 510.0 {
        (0.0, 1.0, -(w - 510.0) / (510.0 - 490.0))
    } else if w < 580.0 {
        ((w - 510.0) / (580.0 - 510.0), 1.0, 0.0)
    } else if w < 645.0 {
        (1.0, -(w - 645.0) / (645.0 - 580.0), 0.0)
    } else {
        (1.0, 0.0, 0.0)
    };

    // Intensity falloff at the edges of human vision.
    let factor = if w < 420.0 {
        0.3 + 0.7 * (w - 380.0) / (420.0 - 380.0)
    } else if w > 700.0 {
        0.3 + 0.7 * (750.0 - w) / (750.0 - 700.0)
    } else {
        1.0
    };

    let gamma = 0.8;
    let conv = |c: f64| -> u8 {
        (255.0 * (c.clamp(0.0, 1.0) * factor).powf(gamma)).round() as u8
    };
    (conv(r), conv(g), conv(b))
}

#[derive(Default)]
pub struct KittyPlotState {
    pub active_shown: bool,
    pub last_hash: u64,
    pub seq: u64,
    pub prev_tmp: Option<PathBuf>,
}

impl Drop for KittyPlotState {
    fn drop(&mut self) {
        if let Some(path) = self.prev_tmp.take() {
            if path.exists() {
                let _ = fs::remove_file(path);
            }
        }
    }
}

pub fn clear_kitty_image(stdout: &mut impl Write) {
    let _ = write!(stdout, "\x1b_Ga=d,d=I,i={},q=2\x1b\\", IMAGE_ID);
    let _ = stdout.flush();
}

fn base64_encode(input: &[u8]) -> String {
    const TABLE: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(input.len().div_ceil(3) * 4);
    for chunk in input.chunks(3) {
        let b0 = chunk[0] as u32;
        let b1 = if chunk.len() > 1 { chunk[1] as u32 } else { 0 };
        let b2 = if chunk.len() > 2 { chunk[2] as u32 } else { 0 };
        let triple = (b0 << 16) | (b1 << 8) | b2;
        out.push(TABLE[((triple >> 18) & 0x3F) as usize] as char);
        out.push(TABLE[((triple >> 12) & 0x3F) as usize] as char);
        if chunk.len() > 1 {
            out.push(TABLE[((triple >> 6) & 0x3F) as usize] as char);
        } else {
            out.push('=');
        }
        if chunk.len() > 2 {
            out.push(TABLE[(triple & 0x3F) as usize] as char);
        } else {
            out.push('=');
        }
    }
    out
}

/// Built-in 5×7 bitmap glyphs (one 5-bit row per element, MSB = leftmost pixel).
/// Only the characters needed for peak annotations are included.
fn glyph_5x7(ch: char) -> Option<[u8; 7]> {
    let rows = match ch {
        '0' => [0b01110, 0b10001, 0b10011, 0b10101, 0b11001, 0b10001, 0b01110],
        '1' => [0b00100, 0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110],
        '2' => [0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b01000, 0b11111],
        '3' => [0b11111, 0b00010, 0b00100, 0b00010, 0b00001, 0b10001, 0b01110],
        '4' => [0b00010, 0b00110, 0b01010, 0b10010, 0b11111, 0b00010, 0b00010],
        '5' => [0b11111, 0b10000, 0b11110, 0b00001, 0b00001, 0b10001, 0b01110],
        '6' => [0b00110, 0b01000, 0b10000, 0b11110, 0b10001, 0b10001, 0b01110],
        '7' => [0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b01000, 0b01000],
        '8' => [0b01110, 0b10001, 0b10001, 0b01110, 0b10001, 0b10001, 0b01110],
        '9' => [0b01110, 0b10001, 0b10001, 0b01111, 0b00001, 0b00010, 0b01100],
        '.' => [0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b01100, 0b01100],
        '-' => [0b00000, 0b00000, 0b00000, 0b11111, 0b00000, 0b00000, 0b00000],
        '*' => [0b00100, 0b10101, 0b01110, 0b11111, 0b01110, 0b10101, 0b00100],
        _ => return None,
    };
    Some(rows)
}

/// Pixel bounds of the dark pill drawn behind `text` when centred at `(cx, cy)`.
fn label_rect(cx: i32, cy: i32, text: &str, w: i32, h: i32) -> (i32, i32, i32, i32) {
    let n = text.chars().count().max(1) as i32;
    let text_w = LABEL_ADVANCE * n - (LABEL_ADVANCE - GLYPH_W);
    let mut x0 = cx - text_w / 2 - LABEL_PAD;
    let mut y0 = cy - GLYPH_H / 2 - 3;
    let mut x1 = x0 + text_w + LABEL_PAD * 2 - 1;
    let mut y1 = y0 + GLYPH_H + 5;

    // Keep the whole pill inside [2..w-2, 2..h-2].
    if x1 > w - 2 {
        let d = x1 - (w - 2);
        x0 -= d;
        x1 -= d;
    }
    if x0 < 2 {
        let d = 2 - x0;
        x0 += d;
        x1 += d;
    }
    if y1 > h - 2 {
        let d = y1 - (h - 2);
        y0 -= d;
        y1 -= d;
    }
    if y0 < 2 {
        let d = 2 - y0;
        y0 += d;
        y1 += d;
    }

    (
        x0.max(2),
        y0.max(2),
        x1.min(w - 2),
        y1.min(h - 2),
    )
}

fn rects_overlap(a: (i32, i32, i32, i32), b: (i32, i32, i32, i32)) -> bool {
    a.0 <= b.2 && b.0 <= a.2 && a.1 <= b.3 && b.1 <= a.3
}

struct PixelCanvas {
    w: usize,
    h: usize,
    buf: Vec<u8>,
}

impl PixelCanvas {
    fn new(w: usize, h: usize) -> Self {
        let mut buf = vec![0u8; w * h * 4];
        for px in buf.as_chunks_mut::<4>().0 {
            px.copy_from_slice(&BG_RGBA);
        }
        Self { w, h, buf }
    }

    #[inline]
    fn blend_px(&mut self, x: i32, y: i32, col: (u8, u8, u8), alpha: f32) {
        if x < 0 || y < 0 || x >= self.w as i32 || y >= self.h as i32 || alpha <= 0.0 {
            return;
        }
        let a = alpha.clamp(0.0, 1.0);
        let idx = ((y as usize) * self.w + (x as usize)) * 4;
        let inv = 1.0 - a;
        self.buf[idx] = (self.buf[idx] as f32 * inv + col.0 as f32 * a) as u8;
        self.buf[idx + 1] = (self.buf[idx + 1] as f32 * inv + col.1 as f32 * a) as u8;
        self.buf[idx + 2] = (self.buf[idx + 2] as f32 * inv + col.2 as f32 * a) as u8;
        self.buf[idx + 3] = 255;
    }

    fn fill_rect_alpha(&mut self, x0: i32, y0: i32, x1: i32, y1: i32, col: (u8, u8, u8), alpha: f32) {
        let xa = x0.min(x1).max(0);
        let xb = x0.max(x1).min(self.w as i32 - 1);
        let ya = y0.min(y1).max(0);
        let yb = y0.max(y1).min(self.h as i32 - 1);
        for y in ya..=yb {
            for x in xa..=xb {
                self.blend_px(x, y, col, alpha);
            }
        }
    }

    fn draw_hline_dashed(&mut self, y: i32, col: (u8, u8, u8), alpha: f32, dash: i32) {
        for x in 0..(self.w as i32) {
            if (x / dash) % 2 == 0 {
                self.blend_px(x, y, col, alpha);
            }
        }
    }

    fn draw_vline_dashed(&mut self, x: i32, y0: i32, y1: i32, col: (u8, u8, u8), alpha: f32, dash: i32) {
        let ya = y0.min(y1).max(0);
        let yb = y0.max(y1).min(self.h as i32 - 1);
        for y in ya..=yb {
            if (y / dash) % 2 == 0 {
                self.blend_px(x, y, col, alpha);
            }
        }
    }

    /// Smooth anti-aliased continuous line with optional subtle glow
    #[allow(clippy::too_many_arguments)]
    fn draw_line_aa(&mut self, x0: f64, y0: f64, x1: f64, y1: f64, col: (u8, u8, u8), alpha: f32, thick: bool) {
        let dx = x1 - x0;
        let dy = y1 - y0;
        let steps = dx.abs().max(dy.abs()).ceil().max(1.0) as usize;
        let inv = 1.0 / (steps as f64);
        for s in 0..=steps {
            let t = (s as f64) * inv;
            let fx = x0 + dx * t;
            let fy = y0 + dy * t;
            let ix = fx.floor() as i32;
            let iy = fy.floor() as i32;
            let rx = (fx - ix as f64) as f32;
            let ry = (fy - iy as f64) as f32;

            // Bilinear sub-pixel coverage for razor-sharp anti-aliasing
            self.blend_px(ix, iy, col, alpha * (1.0 - rx) * (1.0 - ry));
            self.blend_px(ix + 1, iy, col, alpha * rx * (1.0 - ry));
            self.blend_px(ix, iy + 1, col, alpha * (1.0 - rx) * ry);
            self.blend_px(ix + 1, iy + 1, col, alpha * rx * ry);

            if thick {
                let halo = alpha * 0.28;
                self.blend_px(ix - 1, iy, col, halo);
                self.blend_px(ix + 2, iy, col, halo);
                self.blend_px(ix, iy - 1, col, halo);
                self.blend_px(ix, iy + 2, col, halo);
            }
        }
    }

    fn draw_diamond(&mut self, cx: i32, cy: i32, r: i32, col: (u8, u8, u8)) {
        for dy in -r..=r {
            let rem = r - dy.abs();
            for dx in -rem..=rem {
                self.blend_px(cx + dx, cy + dy, col, 0.95);
            }
        }
    }

    /// Crisp 2× bitmap text on a translucent dark pill, centred at `(cx, cy)`.
    fn draw_label(&mut self, cx: i32, cy: i32, text: &str, col: (u8, u8, u8)) {
        if text.is_empty() || self.w < 8 || self.h < 8 {
            return;
        }
        let (x0, y0, x1, y1) = label_rect(cx, cy, text, self.w as i32, self.h as i32);
        self.fill_rect_alpha(x0, y0, x1, y1, (8, 12, 10), 0.78);

        let n = text.chars().count().max(1) as i32;
        let text_w = LABEL_ADVANCE * n - (LABEL_ADVANCE - GLYPH_W);
        let mut pen_x = (x0 + x1 + 1) / 2 - text_w / 2;
        let pen_y = (y0 + y1 + 1) / 2 - GLYPH_H / 2;

        for ch in text.chars() {
            if let Some(glyph) = glyph_5x7(ch) {
                for (row, bits) in glyph.iter().enumerate() {
                    for bit in 0..5 {
                        if bits & (1 << (4 - bit)) == 0 {
                            continue;
                        }
                        let px = pen_x + bit * LABEL_SCALE;
                        let py = pen_y + row as i32 * LABEL_SCALE;
                        for dy in 0..LABEL_SCALE {
                            for dx in 0..LABEL_SCALE {
                                self.blend_px(px + dx, py + dy, col, 1.0);
                            }
                        }
                    }
                }
            }
            pen_x += LABEL_ADVANCE;
        }
    }
}

fn compute_plot_hash(app: &App) -> u64 {
    let mut h = DefaultHasher::new();
    app.plot_rect.x.hash(&mut h);
    app.plot_rect.y.hash(&mut h);
    app.plot_rect.width.hash(&mut h);
    app.plot_rect.height.hash(&mut h);
    app.view_x_min.to_bits().hash(&mut h);
    app.view_x_max.to_bits().hash(&mut h);
    app.cursor_x.to_bits().hash(&mut h);
    app.drag_start_x.map(f64::to_bits).hash(&mut h);
    app.drag_curr_x.map(f64::to_bits).hash(&mut h);
    app.active_trace.hash(&mut h);
    app.stacked.hash(&mut h);
    app.show_peaks.hash(&mut h);
    app.show_solvent_markers.hash(&mut h);
    app.nmr_solvent.index().hash(&mut h);
    app.process_cfg.y_pref.hash(&mut h);
    app.process_cfg.baseline.hash(&mut h);
    (app.process_cfg.smooth as u8).hash(&mut h);
    app.process_cfg.normalize.hash(&mut h);
    app.pinned_peaks.len().hash(&mut h);
    for p in &app.pinned_peaks {
        p.x.to_bits().hash(&mut h);
    }
    app.peak_threshold.to_bits().hash(&mut h);
    app.auto_peaks.len().hash(&mut h);
    (app.tab as u8).hash(&mut h);
    app.traces.len().hash(&mut h);
    for (i, t) in app.traces.iter().enumerate() {
        t.visible.hash(&mut h);
        t.color_idx.hash(&mut h);
        t.name.hash(&mut h);
        if let Some(pt) = app.processed.get(i) {
            pt.x.len().hash(&mut h);
            pt.y.first().copied().map(f64::to_bits).hash(&mut h);
            pt.y.last().copied().map(f64::to_bits).hash(&mut h);
        }
    }
    h.finish()
}

pub fn sync_hd_plot(
    app: &App,
    state: &mut KittyPlotState,
    stdout: &mut impl Write,
) -> Result<()> {
    let rect = app.plot_rect;
    let should_show = app.hd_mode
        && matches!(app.modal, Modal::None)
        && app.tab.as_kind().is_some()
        && rect.width >= 12
        && rect.height >= 4;

    if !should_show {
        if state.active_shown {
            clear_kitty_image(stdout);
            state.active_shown = false;
            state.last_hash = 0;
        }
        return Ok(());
    }

    let hash = compute_plot_hash(app);
    if state.active_shown && state.last_hash == hash {
        return Ok(());
    }

    let w = (rect.width as usize * 8).clamp(320, 1280);
    let h = (rect.height as usize * 16).clamp(180, 760);
    let mut img = PixelCanvas::new(w, h);

    // 1. Subtle laboratory oscilloscope graticule grid (5×6 divisions)
    for row_div in 1..5 {
        let gy = (h * row_div / 5) as i32;
        img.draw_hline_dashed(gy, GRID_RGB, 0.75, 4);
    }
    for col_div in 1..6 {
        let gx = (w * col_div / 6) as i32;
        img.draw_vline_dashed(gx, 0, h as i32 - 1, GRID_RGB, 0.75, 4);
    }

    let x_lo = app.view_x_min;
    let x_hi = app.view_x_max;
    let x_span = (x_hi - x_lo).max(1e-9);
    let kind = app.current_kind();

    let visible_indices: Vec<usize> = app
        .traces
        .iter()
        .enumerate()
        .filter(|(idx, t)| {
            t.visible
                && app
                    .processed
                    .get(*idx)
                    .is_some_and(|p| !p.x.is_empty())
        })
        .map(|(idx, _)| idx)
        .collect();

    let y_unit = app
        .processed
        .get(app.active_trace)
        .map(|p| p.y_unit)
        .unwrap_or(YUnit::Transmittance);

    let mut global_y_min = f64::INFINITY;
    let mut global_y_max = f64::NEG_INFINITY;
    for &idx in &visible_indices {
        let pt = &app.processed[idx];
        let i0 = pt.x.partition_point(|&v| v < x_lo);
        let i1 = pt.x.partition_point(|&v| v <= x_hi);
        for &yv in &pt.y[i0..i1.min(pt.y.len())] {
            global_y_min = global_y_min.min(yv);
            global_y_max = global_y_max.max(yv);
        }
    }
    if !global_y_min.is_finite()
        || !global_y_max.is_finite()
        || (global_y_max - global_y_min).abs() < 1e-9
    {
        global_y_min = 0.0;
        global_y_max = if y_unit == YUnit::Transmittance {
            100.0
        } else {
            1.0
        };
    }

    let y_span_raw = (global_y_max - global_y_min).max(1e-9);
    let n_vis = visible_indices.len().max(1);
    let stack_step = if app.stacked && n_vis > 1 {
        y_span_raw * 1.15
    } else {
        0.0
    };
    let plot_y_min = global_y_min - y_span_raw * 0.08;
    let plot_y_max =
        global_y_max + stack_step * (n_vis.saturating_sub(1) as f64) + y_span_raw * 0.10;
    let plot_y_span = (plot_y_max - plot_y_min).max(1e-9);

    let map_x = |x: f64| -> f64 {
        let u = if kind.is_x_reversed() {
            ((x_hi - x) / x_span).clamp(0.0, 1.0)
        } else {
            ((x - x_lo) / x_span).clamp(0.0, 1.0)
        };
        u * ((w - 1) as f64)
    };
    let map_y = |y: f64| -> f64 {
        let v = ((y - plot_y_min) / plot_y_span).clamp(0.0, 1.0);
        (1.0 - v) * ((h - 1) as f64)
    };

    // 1b. UV-Vis: physical visible-light rainbow strip along the bottom edge,
    // aligned pixel-for-pixel with the [x_lo, x_hi] wavelength axis.
    if kind == SpectrumKind::UvVis && h > 12 {
        let bar_top = h as i32 - 9;
        let bar_bot = h as i32 - 1;
        for col_x in 0..w {
            let u = if w > 1 {
                col_x as f64 / ((w - 1) as f64)
            } else {
                0.0
            };
            let nm = if kind.is_x_reversed() {
                x_hi - u * x_span
            } else {
                x_lo + u * x_span
            };
            let rgb = wavelength_to_rgb(nm);
            img.fill_rect_alpha(col_x as i32, bar_top, col_x as i32, bar_bot, rgb, 0.72);
        }
    }

    // 2. Drag-to-zoom translucent highlight rectangle
    if let (Some(x1), Some(x2)) = (app.drag_start_x, app.drag_curr_x) {
        if (x2 - x1).abs() > x_span * 0.015 {
            let px1 = map_x(x1) as i32;
            let px2 = map_x(x2) as i32;
            img.fill_rect_alpha(px1, 0, px2, h as i32 - 1, DRAG_RGB, 0.22);
            img.draw_vline_dashed(px1, 0, h as i32 - 1, DRAG_RGB, 0.9, 3);
            img.draw_vline_dashed(px2, 0, h as i32 - 1, DRAG_RGB, 0.9, 3);
        }
    }

    // 3. NMR solvent / water dashed coral markers
    if app.show_solvent_markers
        && matches!(kind, SpectrumKind::Hnmr | SpectrumKind::Cnmr)
    {
        let shifts: Vec<f64> = if kind == SpectrumKind::Hnmr {
            let mut v = app.nmr_solvent.residual_h().to_vec();
            v.push(app.nmr_solvent.water_h());
            v
        } else {
            app.nmr_solvent.residual_c().to_vec()
        };
        for s in shifts {
            if s >= x_lo && s <= x_hi {
                let sx = map_x(s) as i32;
                img.draw_vline_dashed(sx, (h as f64 * 0.65) as i32, h as i32 - 1, SOLVENT_RGB, 0.85, 3);
            }
        }
    }

    // 4. Draw continuous anti-aliased spectrum curves (per-pixel min-max envelope!)
    let mut active_stack_offset = 0.0;
    for (vis_pos, &idx) in visible_indices.iter().enumerate() {
        let pt = &app.processed[idx];
        let col = HD_TRACE_PALETTE[app.traces[idx].color_idx % HD_TRACE_PALETTE.len()];
        let is_active = idx == app.active_trace;
        let offset = (vis_pos as f64) * stack_step;
        if is_active {
            active_stack_offset = offset;
        }

        let mut bucket: Vec<Option<(f64, f64, f64, f64)>> = vec![None; w];
        let i0 = pt.x.partition_point(|&v| v < x_lo).saturating_sub(1);
        let i1 = (pt.x.partition_point(|&v| v <= x_hi) + 1).min(pt.x.len());

        let mut push_pt = |k: usize| {
            let px = (map_x(pt.x[k]).round() as usize).min(w - 1);
            let py = map_y(pt.y[k] + offset);
            match &mut bucket[px] {
                None => bucket[px] = Some((py, py, py, py)),
                Some((mn, mx, _, last)) => {
                    if py < *mn { *mn = py; }
                    if py > *mx { *mx = py; }
                    *last = py;
                }
            }
        };
        if kind.is_x_reversed() {
            for k in (i0..i1).rev() {
                push_pt(k);
            }
        } else {
            for k in i0..i1 {
                push_pt(k);
            }
        }

        let alpha = if is_active || visible_indices.len() == 1 { 0.98 } else { 0.72 };
        let mut prev: Option<(f64, f64)> = None;
        for (col_x, slot) in bucket.iter().enumerate() {
            if let Some((mn, mx, first, last)) = *slot {
                let fx = col_x as f64;
                if let Some((px, py)) = prev {
                    img.draw_line_aa(px, py, fx, first, col, alpha, is_active);
                }
                if (mx - mn).abs() > 0.5 {
                    img.draw_line_aa(fx, mn, fx, mx, col, alpha, is_active);
                }
                prev = Some((fx, last));
            }
        }
    }

    // 5. Peak diamonds (auto + pinned) with crisp bitmap annotations
    let mut auto_marks: Vec<(f64, f64)> = Vec::new();
    if app.show_peaks {
        for p in app.auto_peaks.iter().take(12) {
            if p.x >= x_lo && p.x <= x_hi {
                let px = map_x(p.x) as i32;
                let py = map_y(p.y + active_stack_offset) as i32;
                img.draw_diamond(px, py, 3, PEAK_AUTO_RGB);
                auto_marks.push((p.x, p.y));
            }
        }
    }
    let mut pin_marks: Vec<(f64, f64)> = Vec::new();
    for p in &app.pinned_peaks {
        if p.x >= x_lo && p.x <= x_hi {
            let px = map_x(p.x) as i32;
            let py = map_y(p.y + active_stack_offset) as i32;
            img.draw_diamond(px, py, 5, PEAK_PIN_RGB);
            pin_marks.push((p.x, p.y));
        }
    }

    let label_for = |x: f64| -> String {
        match kind {
            SpectrumKind::Ftir => format!("{:.0}", x),
            SpectrumKind::Hnmr => format!("{:.2}", x),
            SpectrumKind::Cnmr => format!("{:.1}", x),
            SpectrumKind::UvVis => format!("{:.0}", x),
        }
    };
    let w_i = w as i32;
    let h_i = h as i32;
    let mut placed: Vec<(i32, i32, i32, i32)> = Vec::new();
    let place_label = |img: &mut PixelCanvas,
                           placed: &mut Vec<(i32, i32, i32, i32)>,
                           px: i32,
                           py: i32,
                           text: &str,
                           col: (u8, u8, u8)| {
        // Absorbance peaks point up (label above); %T troughs point down.
        let trough = y_unit == YUnit::Transmittance;
        let base = if trough { GLYPH_H / 2 + 7 } else { -(GLYPH_H / 2 + 7) };
        let away = if trough { 1 } else { -1 };
        let mut chosen: Option<(i32, i32, i32, i32)> = None;
        let mut chosen_cy = py + base;

        for attempt in 0..8 {
            let stagger = if attempt == 0 {
                0
            } else {
                let step = ((attempt + 1) / 2) * 18;
                let dir = if attempt % 2 == 1 { away } else { -away };
                step * dir
            };
            let cy = py + base + stagger;
            let rect = label_rect(px, cy, text, w_i, h_i);
            if !placed.iter().any(|b| rects_overlap(*b, rect)) {
                chosen = Some(rect);
                chosen_cy = cy;
                break;
            }
        }

        let rect = chosen.unwrap_or_else(|| label_rect(px, chosen_cy, text, w_i, h_i));
        placed.push(rect);
        img.draw_label(px, chosen_cy, text, col);
    };

    // Pinned peaks get first pick of label positions, then the top-10 auto peaks.
    for (x, y) in &pin_marks {
        let px = map_x(*x) as i32;
        let py = map_y(*y + active_stack_offset) as i32;
        let text = format!("*{}", label_for(*x));
        place_label(&mut img, &mut placed, px, py, &text, PEAK_PIN_RGB);
    }
    if app.show_peaks {
        for (x, y) in auto_marks.iter().take(10) {
            let px = map_x(*x) as i32;
            let py = map_y(*y + active_stack_offset) as i32;
            let text = label_for(*x);
            place_label(&mut img, &mut placed, px, py, &text, PEAK_LABEL_AUTO_RGB);
        }
    }

    // 6. Crosshair vertical line + glowing intersection marker on active trace
    let cx = map_x(app.cursor_x) as i32;
    img.draw_vline_dashed(cx, 0, h as i32 - 1, CROSSHAIR_RGB, 0.85, 3);
    if let Some(pt) = app.processed.get(app.active_trace) {
        if let Some(cy_val) = interpolate_y(&pt.x, &pt.y, app.cursor_x) {
            let cy = map_y(cy_val + active_stack_offset) as i32;
            img.draw_diamond(cx, cy, 4, CROSSHAIR_RGB);
        }
    }

    // 7. Transmit via Kitty Graphics Protocol temp-file transport (auto-cleaned by Kitty!)
    state.seq = state.seq.wrapping_add(1);
    if let Some(old) = state.prev_tmp.take() {
        let _ = fs::remove_file(old);
    }
    let tmp_path: PathBuf = env::temp_dir().join(format!(
        "spectra-tty-graphics-protocol-{}-{}.rgba",
        std::process::id(),
        state.seq
    ));
    fs::write(&tmp_path, &img.buf)?;
    state.prev_tmp = Some(tmp_path.clone());
    let b64_path = base64_encode(tmp_path.to_string_lossy().as_bytes());

    write!(
        stdout,
        "\x1b[{};{}H\x1b_Ga=T,f=32,s={},v={},c={},r={},i={},C=1,t=t,q=2;{}\x1b\\",
        rect.y + 1,
        rect.x + 1,
        w,
        h,
        rect.width,
        rect.height,
        IMAGE_ID,
        b64_path
    )?;
    stdout.flush()?;

    state.active_shown = true;
    state.last_hash = hash;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bitmap_glyphs_cover_peak_annotation_charset() {
        for ch in "0123456789.-*".chars() {
            let glyph = glyph_5x7(ch);
            assert!(glyph.is_some(), "missing glyph for {ch:?}");
            assert_eq!(glyph.unwrap().len(), 7);
        }
        assert!(glyph_5x7('x').is_none());
        assert_eq!(glyph_5x7('8').unwrap()[0], 0b01110);
        assert_eq!(glyph_5x7('-').unwrap()[3], 0b11111);
    }

    #[test]
    fn labels_stay_inside_canvas_and_paint_pixels() {
        let mut img = PixelCanvas::new(200, 80);
        let before = img.buf.clone();
        img.draw_label(2, 2, "1728", PEAK_PIN_RGB);
        assert!(img.buf != before, "label should paint pixels");

        let (x0, y0, x1, y1) = label_rect(2, 2, "1728", 200, 80);
        assert!(
            x0 >= 2 && y0 >= 2 && x1 <= 198 && y1 <= 78,
            "label rect escaped canvas: {:?}",
            (x0, y0, x1, y1)
        );

        // Off-canvas anchors are clamped, never panic.
        img.draw_label(-500, -500, "3.60", PEAK_PIN_RGB);
        img.draw_label(9999, 9999, "*18.7", PEAK_PIN_RGB);
    }

    #[test]
    fn staggered_labels_stop_overlapping() {
        let base = label_rect(50, 100, "1728", 400, 200);
        let same_row = label_rect(60, 100, "1728", 400, 200);
        assert!(rects_overlap(base, same_row));

        // Two 18 px steps are needed because a pill is 20 px tall.
        let staggered = label_rect(60, 100 - 36, "1728", 400, 200);
        assert!(
            !rects_overlap(base, staggered),
            "36 px stagger should clear the previous label"
        );
    }
}
