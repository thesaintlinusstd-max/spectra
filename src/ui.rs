use ratatui::layout::{Alignment, Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::symbols::Marker;
use ratatui::text::{Line as TextLine, Span};
use ratatui::widgets::canvas::{Canvas, Line as CanvasLine};
use ratatui::widgets::{Cell, Clear, Paragraph, Row, Table, Wrap};
use ratatui::Frame;

use crate::app::{ActiveTab, App, Modal, PolymerSubTab, MONOMER_PRESETS};
use crate::data::{
    match_cnmr, match_ftir, match_hnmr, match_polymers, match_solvent_impurities,
    match_uvvis, photon_energy_ev, spectral_zone_name, uvvis_color_info, DemoPreset, SolventId, SpectrumKind, YUnit, SOLVENT_IMPURITIES,
};
use crate::processing::interpolate_y;
use crate::theme::*;

pub fn render(frame: &mut Frame, app: &mut App) {
    let area = frame.area();
    if area.width < 45 || area.height < 14 {
        let msg = Paragraph::new("Terminal too small for SPECTRA (min 45×14)")
            .style(Style::default().fg(ACCENT_AMBER))
            .alignment(Alignment::Center);
        frame.render_widget(msg, area);
        return;
    }

    let root = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(1), // Top bar
            Constraint::Min(10),   // Main workspace
            Constraint::Length(2), // Bottom footer
        ])
        .split(area);

    app.topbar_rect = root[0];
    draw_topbar(frame, app, root[0]);

    match app.tab {
        ActiveTab::Ftir | ActiveTab::Hnmr | ActiveTab::Cnmr | ActiveTab::UvVis => {
            draw_spectrum_workspace(frame, app, root[1]);
        }
        ActiveTab::Solvents => draw_solvents_tab(frame, app, root[1]),
        ActiveTab::PolymerLab => draw_polymer_lab_tab(frame, app, root[1]),
    }

    draw_footer(frame, app, root[2]);

    if !matches!(app.modal, Modal::None) {
        draw_modal(frame, app, area);
    }
}

fn draw_topbar(frame: &mut Frame, app: &App, area: Rect) {
    let version = env!("CARGO_PKG_VERSION");
    let short_version = {
        let mut it = version.split('.');
        match (it.next(), it.next()) {
            (Some(a), Some(b)) => format!("v{}.{}", a, b),
            _ => format!("v{}", version),
        }
    };
    let left_spans = vec![
        Span::styled(
            format!(" ✦ SPECTRA {} ", short_version),
            Style::default()
                .fg(ACCENT_GREEN)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled("· FTIR · NMR · UV-Vis ", Style::default().fg(TEXT_MUTED)),
        tab_span("1:FTIR", app.tab == ActiveTab::Ftir),
        Span::raw(" "),
        tab_span("2:¹H NMR", app.tab == ActiveTab::Hnmr),
        Span::raw(" "),
        tab_span("3:¹³C NMR", app.tab == ActiveTab::Cnmr),
        Span::raw(" "),
        tab_span("4:UV-Vis", app.tab == ActiveTab::UvVis),
        Span::raw(" "),
        tab_span("5:Solvents", app.tab == ActiveTab::Solvents),
        Span::raw(" "),
        tab_span("6:Polymer Lab", app.tab == ActiveTab::PolymerLab),
    ];

    let right_spans = vec![
        key_badge("[f] Open", ACCENT_AQUA),
        Span::raw(" "),
        key_badge("[e] Export", ACCENT_AMBER),
        Span::raw(" "),
        key_badge("[/] Search", ACCENT_YELLOW),
        Span::raw(" "),
        key_badge("[?] Help", ACCENT_GREEN),
        Span::raw(" "),
        key_badge("[q] Quit", ACCENT_CORAL),
    ];

    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Min(40), Constraint::Length(60)])
        .split(area);

    frame.render_widget(Paragraph::new(TextLine::from(left_spans)), cols[0]);
    frame.render_widget(
        Paragraph::new(TextLine::from(right_spans)).alignment(Alignment::Right),
        cols[1],
    );
}

fn draw_spectrum_workspace(frame: &mut Frame, app: &mut App, area: Rect) {
    let panes = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(71), Constraint::Percentage(29)])
        .split(area);

    // Left pane: Spectrum Canvas + Bottom Trace/Peak Strip
    let left_split = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Percentage(78), Constraint::Percentage(22)])
        .split(panes[0]);

    draw_spectrum_canvas(frame, app, left_split[0]);
    draw_trace_and_peaks_strip(frame, app, left_split[1]);

    // Right pane: MENDEL-style Amber Inspector
    draw_inspector(frame, app, panes[1]);
}

fn draw_spectrum_canvas(frame: &mut Frame, app: &mut App, area: Rect) {
    let kind = app.current_kind();
    let y_unit = app
        .processed
        .get(app.active_trace)
        .map(|p| p.y_unit)
        .unwrap_or(YUnit::Transmittance);

    let view_tag = if app.hd_mode { "HD" } else { "Braille" };
    let title = if kind.is_x_reversed() {
        // Reversed IUPAC axis (FTIR/NMR): left = view max, right = view min.
        format!(
            "✦ {} Spectrum [{} · {} · {:.0} → {:.0} {}]",
            kind.display_name(),
            view_tag,
            y_unit.label(),
            app.view_x_max,
            app.view_x_min,
            kind.x_unit()
        )
    } else {
        // Ascending axis (UV-Vis nm): left = view min, right = view max.
        format!(
            "✦ {} Spectrum [{} · {} · {:.0} → {:.0} {}]",
            kind.display_name(),
            view_tag,
            y_unit.label(),
            app.view_x_min,
            app.view_x_max,
            kind.x_unit()
        )
    };
    let outer = panel_block(title, BORDER_LEFT);
    let inner = outer.inner(area);
    frame.render_widget(outer, area);

    if inner.width < 12 || inner.height < 6 {
        return;
    }

    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(1), // Diagnostic zones banner
            Constraint::Min(4),    // Braille Canvas
            Constraint::Length(1), // X-Axis ticks ruler
        ])
        .split(inner);

    // 1. Diagnostic Zone Banner
    let zone_bar = match kind {
        SpectrumKind::Ftir => TextLine::from(vec![
            Span::styled(" 4000 ", Style::default().fg(TEXT_MUTED)),
            Span::styled("│O–H/N–H ", Style::default().fg(ACCENT_AQUA)),
            Span::styled("│C–H ", Style::default().fg(ACCENT_GREEN)),
            Span::styled("│C≡N/NCO ", Style::default().fg(ACCENT_PURPLE)),
            Span::styled("│C=O ", Style::default().fg(ACCENT_AMBER)),
            Span::styled("│C=C/Ar ", Style::default().fg(ACCENT_CORAL)),
            Span::styled("│Fingerprint (<1450 cm⁻¹)", Style::default().fg(TEXT_MUTED)),
        ]),
        SpectrumKind::Hnmr => TextLine::from(vec![
            Span::styled(" δ ", Style::default().fg(TEXT_MUTED)),
            Span::styled("│Acid/CHO ", Style::default().fg(ACCENT_PURPLE)),
            Span::styled("│Ar–H ", Style::default().fg(ACCENT_CORAL)),
            Span::styled("│Vinyl ", Style::default().fg(ACCENT_AMBER)),
            Span::styled("│O–CH/N–CH ", Style::default().fg(ACCENT_AQUA)),
            Span::styled("│Backbone/Aliph ", Style::default().fg(ACCENT_GREEN)),
            Span::styled(
                format!("│Solvent: {} [S]", app.nmr_solvent.label()),
                Style::default().fg(ACCENT_YELLOW),
            ),
        ]),
        SpectrumKind::Cnmr => TextLine::from(vec![
            Span::styled(" ¹³C ", Style::default().fg(TEXT_MUTED)),
            Span::styled("│C=O (160–220) ", Style::default().fg(ACCENT_AMBER)),
            Span::styled("│Ar/sp² (110–160) ", Style::default().fg(ACCENT_CORAL)),
            Span::styled("│C–O/OCH₃ (50–80) ", Style::default().fg(ACCENT_AQUA)),
            Span::styled("│Aliph/Tacticity (0–50)", Style::default().fg(ACCENT_GREEN)),
        ]),
        SpectrumKind::UvVis => TextLine::from(vec![
            Span::styled(" λ ", Style::default().fg(TEXT_MUTED)),
            Span::styled("│UV (190–380 nm) ", Style::default().fg(ACCENT_PURPLE)),
            Span::styled("│Visible (380–750 nm) ", Style::default().fg(ACCENT_GREEN)),
            Span::styled("│NIR (>750 nm) ", Style::default().fg(ACCENT_AQUA)),
        ]),
    };
    frame.render_widget(Paragraph::new(zone_bar), rows[0]);

    // Save exact plot inner Rect for mouse coordinate mapping & Kitty HD graphics!
    app.plot_rect = rows[1];

    let x_lo = app.view_x_min;
    let x_hi = app.view_x_max;
    let x_span = (x_hi - x_lo).max(1e-9);
    // Direction flag: for UV-Vis (ascending nm) u=0 at left (x_lo), u=1 at right (x_hi).
    // For FTIR/NMR (reversed) u=0 at left (x_hi), u=1 at right (x_lo).
    let x_reversed = kind.is_x_reversed();
    let u_from_x = |x: f64| -> f64 {
        if x_reversed {
            ((x_hi - x) / x_span).clamp(0.0, 1.0)
        } else {
            ((x - x_lo) / x_span).clamp(0.0, 1.0)
        }
    };

    if !(app.hd_mode && matches!(app.modal, Modal::None)) {
        // 2. Precompute normalized coordinates + Min-Max Envelope Downsampling for Braille
        let visible_indices: Vec<usize> = app
        .traces
        .iter()
        .enumerate()
        .filter(|(idx, t)| t.visible && app.processed.get(*idx).map_or(false, |p| !p.x.is_empty()))
        .map(|(idx, _)| idx)
        .collect();

    let n_vis = visible_indices.len().max(1);
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
    if !global_y_min.is_finite() || !global_y_max.is_finite() || (global_y_max - global_y_min).abs() < 1e-9 {
        global_y_min = 0.0;
        global_y_max = if y_unit == YUnit::Transmittance { 100.0 } else { 1.0 };
    }
    let y_span_raw = (global_y_max - global_y_min).max(1e-9);
    let stack_step = if app.stacked && n_vis > 1 {
        y_span_raw * 1.15
    } else {
        0.0
    };

    let plot_y_min = global_y_min - y_span_raw * 0.08;
    let plot_y_max = global_y_max + stack_step * (n_vis.saturating_sub(1) as f64) + y_span_raw * 0.10;

    // Sub-pixel horizontal bins (Braille gives 2 horizontal dots per column!)
    let bins = ((rows[1].width as usize) * 2).clamp(60, 500);
    let mut segments: Vec<(f64, f64, f64, f64, Color)> = Vec::new();
    let mut active_stack_offset = 0.0;

    for (vis_pos, &idx) in visible_indices.iter().enumerate() {
        let pt = &app.processed[idx];
        let color = trace_color(app.traces[idx].color_idx);
        let offset = (vis_pos as f64) * stack_step;
        if idx == app.active_trace {
            active_stack_offset = offset;
        }

        // Bucket points in [x_lo, x_hi] into screen bins u in [0, 1]
        let mut bucket: Vec<Option<(f64, f64, f64, f64)>> = vec![None; bins]; // (min, max, first, last)
        let i0 = pt.x.partition_point(|&v| v < x_lo).saturating_sub(1);
        let i1 = (pt.x.partition_point(|&v| v <= x_hi) + 1).min(pt.x.len());

        for k in i0..i1 {
            let x = pt.x[k];
            let y = pt.y[k] + offset;
            let u = u_from_x(x);
            let b = ((u * (bins - 1) as f64).round() as usize).min(bins - 1);
            match &mut bucket[b] {
                None => bucket[b] = Some((y, y, y, y)),
                Some((mn, mx, _, last)) => {
                    if y < *mn { *mn = y; }
                    if y > *mx { *mx = y; }
                    *last = y;
                }
            }
        }

        let mut prev: Option<(f64, f64)> = None;
        for (b, slot) in bucket.iter().enumerate() {
            if let Some((mn, mx, first, last)) = *slot {
                let u = (b as f64) / ((bins - 1) as f64);
                if let Some((pu, py)) = prev {
                    segments.push((pu, py, u, first, color));
                }
                if (mx - mn).abs() > 1e-12 {
                    segments.push((u, mn, u, mx, color));
                }
                prev = Some((u, last));
            }
        }
    }

    // Cursor u-coordinate
    let u_cursor = u_from_x(app.cursor_x);

    // Drag-zoom highlight bounds
    let drag_u = match (app.drag_start_x, app.drag_curr_x) {
        (Some(x1), Some(x2)) if (x2 - x1).abs() > x_span * 0.015 => {
            let u1 = u_from_x(x1);
            let u2 = u_from_x(x2);
            Some((u1.min(u2), u1.max(u2)))
        }
        _ => None,
    };

    // Solvent markers (NMR modes)
    let mut solvent_markers: Vec<(f64, String)> = Vec::new();
    if app.show_solvent_markers && matches!(kind, SpectrumKind::Hnmr | SpectrumKind::Cnmr) {
        let shifts: Vec<(f64, &str)> = if kind == SpectrumKind::Hnmr {
            let mut v: Vec<(f64, &str)> = app
                .nmr_solvent
                .residual_h()
                .iter()
                .map(|&s| (s, app.nmr_solvent.label()))
                .collect();
            v.push((app.nmr_solvent.water_h(), "H₂O"));
            v
        } else {
            app.nmr_solvent
                .residual_c()
                .iter()
                .map(|&s| (s, app.nmr_solvent.label()))
                .collect()
        };
        for (s, lbl) in shifts {
            if s >= x_lo && s <= x_hi {
                let u = u_from_x(s);
                solvent_markers.push((u, lbl.to_string()));
            }
        }
    }

    // Peak markers (auto + pinned)
    let mut peak_markers: Vec<(f64, f64, String, bool)> = Vec::new();
    if app.show_peaks {
        for p in app.auto_peaks.iter().take(10) {
            if p.x >= x_lo && p.x <= x_hi {
                let u = u_from_x(p.x);
                let tag = match kind {
                    SpectrumKind::Ftir => format!("{:.0}", p.x),
                    SpectrumKind::Hnmr => format!("{:.2}", p.x),
                    SpectrumKind::Cnmr => format!("{:.1}", p.x),
                    SpectrumKind::UvVis => format!("{:.0}", p.x),
                };
                peak_markers.push((u, p.y + active_stack_offset, tag, false));
            }
        }
    }
    for p in &app.pinned_peaks {
        if p.x >= x_lo && p.x <= x_hi {
            let u = u_from_x(p.x);
            let tag = match kind {
                SpectrumKind::Ftir => format!("★{:.0}", p.x),
                SpectrumKind::Hnmr => format!("★{:.2}", p.x),
                SpectrumKind::Cnmr => format!("★{:.1}", p.x),
                SpectrumKind::UvVis => format!("★{:.0}", p.x),
            };
            peak_markers.push((u, p.y + active_stack_offset, tag, true));
        }
    }

    let is_trans = y_unit == YUnit::Transmittance;
    let y_total_span = (plot_y_max - plot_y_min).max(1e-9);

    let canvas = Canvas::default()
        .background_color(BG_TRANSPARENT)
        .marker(Marker::Braille)
        .x_bounds([0.0, 1.0])
        .y_bounds([plot_y_min, plot_y_max])
        .paint(move |ctx| {
            // Subtle baseline / reference grid line
            let mid_y = plot_y_min + 0.5 * y_total_span;
            ctx.draw(&CanvasLine {
                x1: 0.0,
                y1: mid_y,
                x2: 1.0,
                y2: mid_y,
                color: Color::Rgb(42, 54, 44),
            });

            // Drag-zoom selection edges
            if let Some((du1, du2)) = drag_u {
                ctx.draw(&CanvasLine {
                    x1: du1,
                    y1: plot_y_min,
                    x2: du1,
                    y2: plot_y_max,
                    color: ACCENT_AQUA,
                });
                ctx.draw(&CanvasLine {
                    x1: du2,
                    y1: plot_y_min,
                    x2: du2,
                    y2: plot_y_max,
                    color: ACCENT_AQUA,
                });
            }

            // NMR Solvent impurity reference lines
            for (su, lbl) in &solvent_markers {
                ctx.draw(&CanvasLine {
                    x1: *su,
                    y1: plot_y_min,
                    x2: *su,
                    y2: plot_y_min + y_total_span * 0.25,
                    color: ACCENT_CORAL,
                });
                ctx.print(
                    (*su - 0.02).clamp(0.0, 0.92),
                    plot_y_min + y_total_span * 0.03,
                    TextLine::from(Span::styled(lbl.clone(), Style::default().fg(ACCENT_CORAL))),
                );
            }

            // Plot spectrum traces
            for &(x1, y1, x2, y2, color) in &segments {
                ctx.draw(&CanvasLine { x1, y1, x2, y2, color });
            }

            // Active crosshair line
            ctx.draw(&CanvasLine {
                x1: u_cursor,
                y1: plot_y_min,
                x2: u_cursor,
                y2: plot_y_max,
                color: ACCENT_YELLOW,
            });

            // Peak labels
            for (pu, py, tag, pinned) in &peak_markers {
                let col = if *pinned { ACCENT_YELLOW } else { TEXT_MUTED };
                let label_y = if is_trans {
                    (*py - y_total_span * 0.06).max(plot_y_min + y_total_span * 0.02)
                } else {
                    (*py + y_total_span * 0.05).min(plot_y_max - y_total_span * 0.04)
                };
                ctx.print(
                    (*pu - 0.015).clamp(0.0, 0.93),
                    label_y,
                    TextLine::from(Span::styled(
                        tag.clone(),
                        Style::default().fg(col).add_modifier(if *pinned {
                            Modifier::BOLD
                        } else {
                            Modifier::empty()
                        }),
                    )),
                );
            }
        });

        frame.render_widget(canvas, rows[1]);
    }

    // 3. X-Axis Numeric Ticks Ruler (respects kind.is_x_reversed(): UV-Vis ascending, FTIR/NMR descending)
    let n_ticks = 6usize;
    let w = rows[2].width as usize;
    let mut ruler = vec![' '; w];
    for t in 0..=n_ticks {
        let frac = t as f64 / n_ticks as f64;
        let x_val = if x_reversed { x_hi - frac * x_span } else { x_lo + frac * x_span };
        let lbl = match kind {
            SpectrumKind::Ftir => format!("{:.0}", x_val),
            SpectrumKind::Hnmr => format!("{:.2}", x_val),
            SpectrumKind::Cnmr => format!("{:.0}", x_val),
            SpectrumKind::UvVis => format!("{:.0}", x_val),
        };
        let pos = ((frac * (w.saturating_sub(lbl.len())) as f64).round() as usize)
            .min(w.saturating_sub(lbl.len()));
        for (j, ch) in lbl.chars().enumerate() {
            if pos + j < w {
                ruler[pos + j] = ch;
            }
        }
    }
    let ruler_str: String = ruler.into_iter().collect();
    frame.render_widget(
        Paragraph::new(ruler_str).style(Style::default().fg(TEXT_MUTED)),
        rows[2],
    );
}

fn draw_trace_and_peaks_strip(frame: &mut Frame, app: &App, area: Rect) {
    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(46), Constraint::Percentage(54)])
        .split(area);

    // Left box: Loaded Traces
    let t_block = panel_block(
        "Traces [Tab:Active · Space:Hide · Del:Remove · s:Stack]",
        BORDER_LEFT,
    );
    let mut t_lines = Vec::new();
    for (idx, tr) in app.traces.iter().enumerate() {
        let col = trace_color(tr.color_idx);
        let active_mark = if idx == app.active_trace { "▶ " } else { "  " };
        let vis_mark = if tr.visible { "●" } else { "○" };
        let y_at = app
            .processed
            .get(idx)
            .and_then(|pt| interpolate_y(&pt.x, &pt.y, app.cursor_x))
            .map(|v| format!("{:.2}", v))
            .unwrap_or_else(|| "—".into());
        t_lines.push(TextLine::from(vec![
            Span::styled(active_mark, Style::default().fg(ACCENT_YELLOW)),
            Span::styled(format!("{} ", vis_mark), Style::default().fg(col)),
            Span::styled(
                tr.name.clone(),
                Style::default()
                    .fg(if idx == app.active_trace { TEXT_PRIMARY } else { TEXT_MUTED })
                    .add_modifier(if idx == app.active_trace { Modifier::BOLD } else { Modifier::empty() }),
            ),
            Span::styled(
                format!("  [{} pts | Y={}]", tr.x.len(), y_at),
                Style::default().fg(TEXT_MUTED),
            ),
        ]));
    }
    if t_lines.is_empty() {
        t_lines.push(TextLine::from(Span::styled(
            "  No traces loaded — press [f] to open or [d] for presets",
            Style::default().fg(TEXT_MUTED),
        )));
    }
    frame.render_widget(Paragraph::new(t_lines).block(t_block), cols[0]);

    // Right box: Detected & Pinned Peaks
    let p_title = format!(
        "Peaks [n/N:Jump · x:Pin · [/]:Thresh {:.0}%] {}/{}",
        app.peak_threshold * 100.0,
        if app.auto_peaks.is_empty() {
            0
        } else {
            app.selected_peak + 1
        },
        app.auto_peaks.len()
    );
    let p_block = panel_block(p_title, BORDER_LEFT);
    let inner_h = p_block.inner(cols[1]).height as usize;
    let mut p_lines = Vec::new();
    let pinned_shown = app.pinned_peaks.len().min(4);
    for pin in app.pinned_peaks.iter().take(pinned_shown) {
        p_lines.push(TextLine::from(vec![
            Span::styled(" ★ ", Style::default().fg(ACCENT_YELLOW).add_modifier(Modifier::BOLD)),
            Span::styled(pin.label.clone(), Style::default().fg(TEXT_PRIMARY)),
        ]));
    }

    // Scroll the list so the peak selected with [n]/[N] is always visible.
    let total = app.auto_peaks.len();
    let avail = inner_h.saturating_sub(pinned_shown).max(1);
    let start = if total <= avail {
        0
    } else {
        let half = avail / 2;
        let mut s = app.selected_peak.saturating_sub(half);
        if s + avail > total {
            s = total - avail;
        }
        s
    };
    for (i, pk) in app.auto_peaks.iter().enumerate().skip(start).take(avail) {
        let sel = i == app.selected_peak;
        p_lines.push(TextLine::from(vec![
            Span::styled(
                if sel { " ▸ " } else { " · " },
                Style::default().fg(if sel { ACCENT_GREEN } else { TEXT_MUTED }),
            ),
            Span::styled(
                pk.label.clone(),
                Style::default().fg(if sel { ACCENT_GREEN } else { TEXT_PRIMARY }),
            ),
            Span::styled(
                format!("  (prom {:.0}%)", pk.prominence * 100.0),
                Style::default().fg(TEXT_MUTED),
            ),
        ]));
    }
    if total == 0 {
        p_lines.push(TextLine::from(Span::styled(
            "  No peaks above threshold — adjust with [ / ]",
            Style::default().fg(TEXT_MUTED),
        )));
    }
    frame.render_widget(Paragraph::new(p_lines).block(p_block), cols[1]);
}

fn draw_inspector(frame: &mut Frame, app: &App, area: Rect) {
    let kind = app.current_kind();
    let title = format!("Inspector · {}", kind.display_name());
    let outer = panel_block(title, BORDER_RIGHT);
    let inner = outer.inner(area);
    frame.render_widget(outer, area);

    if inner.height < 8 {
        return;
    }

    let header_h = if kind == SpectrumKind::UvVis { 5 } else { 4 };
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(header_h), // Cursor Readout Card
            Constraint::Min(6),           // Functional groups + Solvent warning + Polymer matches
            Constraint::Length(4),        // MENDEL-style Copy/Export action buttons
        ])
        .split(inner);

    // 1. Cursor Readout Card
    let pos_str = match kind {
        SpectrumKind::Ftir => format!("{} = {:.1} cm⁻¹", kind.symbol(), app.cursor_x),
        SpectrumKind::Hnmr => format!("{} = {:.2} ppm", kind.symbol(), app.cursor_x),
        SpectrumKind::Cnmr => format!("{} = {:.1} ppm", kind.symbol(), app.cursor_x),
        SpectrumKind::UvVis => format!("{} = {:.1} nm", kind.symbol(), app.cursor_x),
    };
    let active_y = app
        .processed
        .get(app.active_trace)
        .and_then(|pt| {
            interpolate_y(&pt.x, &pt.y, app.cursor_x).map(|v| format!("{:.2} {}", v, pt.y_unit.label()))
        })
        .unwrap_or_else(|| "—".into());

    let zone = spectral_zone_name(kind, app.cursor_x);

    let mut card_lines = vec![
        TextLine::from(vec![
            Span::styled(
                format!(" {} ", pos_str),
                Style::default()
                    .fg(DARK_TEXT)
                    .bg(ACCENT_AMBER)
                    .add_modifier(Modifier::BOLD),
            ),
            Span::styled(format!("  Y: {}", active_y), Style::default().fg(TEXT_PRIMARY)),
        ]),
        TextLine::from(Span::styled(
            format!(" Zone: {}", zone),
            Style::default().fg(ACCENT_AQUA),
        )),
    ];

    // UV-Vis: add photon energy (eV) and absorbed/observed color info on 2 clean lines
    if kind == SpectrumKind::UvVis {
        let ev = photon_energy_ev(app.cursor_x);
        let (absorbed_region, solution_color) = uvvis_color_info(app.cursor_x);
        card_lines.push(TextLine::from(Span::styled(
            format!(" E = {:.2} eV  ·  {}", ev, absorbed_region),
            Style::default().fg(ACCENT_YELLOW),
        )));
        card_lines.push(TextLine::from(Span::styled(
            format!(" Solution: {}", solution_color),
            Style::default().fg(ACCENT_GREEN),
        )));
    }

    card_lines.push(TextLine::from(Span::styled(
        " ─────────────────────────────",
        Style::default().fg(TEXT_MUTED),
    )));

    frame.render_widget(Paragraph::new(card_lines), chunks[0]);

    // 2. Scientific Assignments + Solvent Alert + Polymer Match
    let mut info: Vec<TextLine> = Vec::new();

    // NMR Solvent/Impurity Alert Box (CRITICAL Lab Lifesaver!)
    if matches!(kind, SpectrumKind::Hnmr | SpectrumKind::Cnmr) {
        let is_c = kind == SpectrumKind::Cnmr;
        let tol = if is_c { 1.8 } else { 0.09 };
        let alerts = match_solvent_impurities(app.nmr_solvent, is_c, app.cursor_x, tol);
        if !alerts.is_empty() {
            info.push(TextLine::from(Span::styled(
                format!(" ⚠ IMPURITY / SOLVENT ALERT ({}):", app.nmr_solvent.label()),
                Style::default()
                    .fg(ACCENT_CORAL)
                    .add_modifier(Modifier::BOLD),
            )));
            for (imp, shift, diff) in alerts.iter().take(3) {
                info.push(TextLine::from(vec![
                    Span::styled(
                        format!("  • {} ", imp.name),
                        Style::default().fg(ACCENT_CORAL).add_modifier(Modifier::BOLD),
                    ),
                    Span::styled(
                        format!("{:.2} ppm (Δ{:.2})", shift, diff),
                        Style::default().fg(ACCENT_YELLOW),
                    ),
                ]));
                info.push(TextLine::from(Span::styled(
                    format!("    {} | {}", imp.mult, imp.note),
                    Style::default().fg(TEXT_MUTED),
                )));
            }
            info.push(TextLine::from(""));
        }
    }

    // Functional Group Matches
    info.push(TextLine::from(Span::styled(
        " ✦ Functional Groups:",
        Style::default().fg(ACCENT_GREEN).add_modifier(Modifier::BOLD),
    )));

    match kind {
        SpectrumKind::Ftir => {
            let matches = match_ftir(app.cursor_x);
            if matches.is_empty() {
                info.push(TextLine::from(Span::styled(
                    "  No primary band in immediate window.",
                    Style::default().fg(TEXT_MUTED),
                )));
            } else {
                for b in matches.iter().take(4) {
                    info.push(TextLine::from(vec![
                        Span::styled("  ▸ ", Style::default().fg(ACCENT_AMBER)),
                        Span::styled(b.group, Style::default().fg(TEXT_PRIMARY).add_modifier(Modifier::BOLD)),
                    ]));
                    info.push(TextLine::from(Span::styled(
                        format!("    {:.0}–{:.0} cm⁻¹ · {} · {}", b.min, b.max, b.intensity, b.mode),
                        Style::default().fg(ACCENT_YELLOW),
                    )));
                    info.push(TextLine::from(Span::styled(
                        format!("    {}", b.note),
                        Style::default().fg(TEXT_MUTED),
                    )));
                    if !b.confirm.is_empty() {
                        info.push(TextLine::from(Span::styled(
                            format!("    ↳ {}", b.confirm),
                            Style::default().fg(ACCENT_AQUA),
                        )));
                    }
                }
            }
        }
        SpectrumKind::Hnmr => {
            let matches = match_hnmr(app.cursor_x);
            if matches.is_empty() {
                info.push(TextLine::from(Span::styled(
                    "  Move cursor over 0–13 ppm to inspect.",
                    Style::default().fg(TEXT_MUTED),
                )));
            } else {
                for b in matches.iter().take(4) {
                    info.push(TextLine::from(vec![
                        Span::styled("  ▸ ", Style::default().fg(ACCENT_AMBER)),
                        Span::styled(b.group, Style::default().fg(TEXT_PRIMARY).add_modifier(Modifier::BOLD)),
                    ]));
                    info.push(TextLine::from(Span::styled(
                        format!("    {:.2}–{:.2} ppm · {}", b.min, b.max, b.example),
                        Style::default().fg(ACCENT_YELLOW),
                    )));
                    info.push(TextLine::from(Span::styled(
                        format!("    {}", b.note),
                        Style::default().fg(TEXT_MUTED),
                    )));
                    if !b.j_hz.is_empty() {
                        info.push(TextLine::from(Span::styled(
                            format!("    J: {}", b.j_hz),
                            Style::default().fg(ACCENT_AQUA),
                        )));
                    }
                }
            }
        }
        SpectrumKind::Cnmr => {
            let matches = match_cnmr(app.cursor_x);
            if matches.is_empty() {
                info.push(TextLine::from(Span::styled(
                    "  Move cursor over 0–220 ppm to inspect.",
                    Style::default().fg(TEXT_MUTED),
                )));
            } else {
                for b in matches.iter().take(4) {
                    info.push(TextLine::from(vec![
                        Span::styled("  ▸ ", Style::default().fg(ACCENT_AMBER)),
                        Span::styled(b.group, Style::default().fg(TEXT_PRIMARY).add_modifier(Modifier::BOLD)),
                    ]));
                    info.push(TextLine::from(Span::styled(
                        format!("    {:.0}–{:.0} ppm · {}", b.min, b.max, b.example),
                        Style::default().fg(ACCENT_YELLOW),
                    )));
                    info.push(TextLine::from(Span::styled(
                        format!("    {}", b.note),
                        Style::default().fg(TEXT_MUTED),
                    )));
                }
            }
        }
        SpectrumKind::UvVis => {
            let bands = match_uvvis(app.cursor_x);
            if bands.is_empty() {
                info.push(TextLine::from(Span::styled(
                    "  Move cursor over 190–900 nm to inspect.",
                    Style::default().fg(TEXT_MUTED),
                )));
            } else {
                for b in bands.iter().take(4) {
                    info.push(TextLine::from(vec![
                        Span::styled("  ▸ ", Style::default().fg(ACCENT_AMBER)),
                        Span::styled(b.group, Style::default().fg(TEXT_PRIMARY).add_modifier(Modifier::BOLD)),
                    ]));
                    info.push(TextLine::from(Span::styled(
                        format!("    {:.0}–{:.0} nm · {}", b.min, b.max, b.transition),
                        Style::default().fg(ACCENT_YELLOW),
                    )));
                    info.push(TextLine::from(Span::styled(
                        format!("    ε {} · {}", b.epsilon, b.note),
                        Style::default().fg(TEXT_MUTED),
                    )));
                }
            }
        }
    }

    // Polymer Fingerprints
    let poly_hits = match_polymers(kind, app.cursor_x);
    if !poly_hits.is_empty() {
        info.push(TextLine::from(""));
        info.push(TextLine::from(Span::styled(
            " ✦ Polymer Signatures:",
            Style::default().fg(ACCENT_PURPLE).add_modifier(Modifier::BOLD),
        )));
        for (p, _, desc) in poly_hits.iter().take(4) {
            info.push(TextLine::from(vec![
                Span::styled(
                    format!("  {} ", p.short),
                    Style::default().fg(ACCENT_AQUA).add_modifier(Modifier::BOLD),
                ),
                Span::styled(desc.to_string(), Style::default().fg(TEXT_PRIMARY)),
            ]));
        }
    }

    frame.render_widget(
        Paragraph::new(info).wrap(Wrap { trim: false }),
        chunks[1],
    );

    // 3. MENDEL-style Quick Action Buttons
    let btn_lines = vec![
        TextLine::from(Span::styled(
            " Quick Actions (Wayland / Python):",
            Style::default().fg(TEXT_MUTED),
        )),
        TextLine::from(vec![
            Span::raw(" "),
            key_badge("[c] Copy Peak", ACCENT_YELLOW),
            Span::raw(" "),
            key_badge("[x] Pin Peak", ACCENT_GREEN),
        ]),
        TextLine::from(vec![
            Span::raw(" "),
            key_badge("[e] Export Py", ACCENT_AMBER),
            Span::raw(" "),
            key_badge("[d] Presets", ACCENT_AQUA),
        ]),
    ];
    frame.render_widget(Paragraph::new(btn_lines), chunks[2]);
}

// ============================================================================
// TAB 4: FULMER 2010 NMR SOLVENTS & IMPURITIES EXPLORER
// ============================================================================
fn draw_solvents_tab(frame: &mut Frame, app: &App, area: Rect) {
    let outer = panel_block(
        "NMR Trace Impurities & Residual Solvents (Fulmer 2010 / Gottlieb 1997)",
        BORDER_LEFT,
    );
    let inner = outer.inner(area);
    frame.render_widget(outer, area);

    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(2), // Solvent pills
            Constraint::Length(3), // Filter bar
            Constraint::Min(6),    // Impurity table
        ])
        .split(inner);

    // Solvent pills
    let mut pills = vec![Span::styled(" Solvent [h/l]: ", Style::default().fg(TEXT_MUTED))];
    for s in SolventId::ALL {
        pills.push(tab_span(s.label(), s == app.solvents_tab.solvent));
        pills.push(Span::raw(" "));
    }
    frame.render_widget(Paragraph::new(TextLine::from(pills)).wrap(Wrap { trim: false }), rows[0]);

    // Filter bar
    let st = &app.solvents_tab;
    let nuc = if st.is_carbon { "¹³C NMR" } else { "¹H NMR" };
    let filter_line = TextLine::from(vec![
        Span::styled(" Nucleus [n]: ", Style::default().fg(TEXT_MUTED)),
        Span::styled(
            format!(" {} ", nuc),
            Style::default().fg(DARK_TEXT).bg(ACCENT_AQUA).add_modifier(Modifier::BOLD),
        ),
        Span::styled("  Tolerance [t]: ", Style::default().fg(TEXT_MUTED)),
        Span::styled(
            format!(" {} ", st.tol_label()),
            Style::default().fg(ACCENT_YELLOW).bg(SURFACE_DARK),
        ),
        Span::styled("  Query shift or name [s/Enter]: ", Style::default().fg(TEXT_MUTED)),
        Span::styled(
            if st.query.is_empty() && !st.editing {
                " <all> ".to_string()
            } else if st.editing {
                format!(" {}█ ", st.query)
            } else {
                format!(" {} ", st.query)
            },
            Style::default()
                .fg(if st.editing { ACCENT_GREEN } else { TEXT_PRIMARY })
                .bg(SURFACE_DARK)
                .add_modifier(Modifier::BOLD),
        ),
        Span::styled("  [x] Clear", Style::default().fg(TEXT_MUTED)),
    ]);
    frame.render_widget(Paragraph::new(filter_line), rows[1]);

    // Filter rows
    let solv_idx = st.solvent.index();
    let q_trim = st.query.trim();
    let q_num = q_trim.parse::<f64>().ok();
    let q_low = q_trim.to_ascii_lowercase();
    let tol = st.tolerance();

    let mut items: Vec<(&str, &str, &str, f64, Option<f64>, &str)> = Vec::new();
    for imp in SOLVENT_IMPURITIES {
        let Some(shift) = (if st.is_carbon { imp.c_shifts[solv_idx] } else { imp.h_shifts[solv_idx] }) else {
            continue;
        };
        if let Some(target) = q_num {
            let diff = (shift - target).abs();
            if diff <= tol {
                items.push((imp.name, imp.group, imp.mult, shift, Some(diff), imp.note));
            }
        } else if q_low.is_empty()
            || imp.name.to_ascii_lowercase().contains(&q_low)
            || imp.group.to_ascii_lowercase().contains(&q_low)
            || imp.note.to_ascii_lowercase().contains(&q_low)
        {
            items.push((imp.name, imp.group, imp.mult, shift, None, imp.note));
        }
    }
    if q_num.is_some() {
        items.sort_by(|a, b| a.4.unwrap_or(999.0).total_cmp(&b.4.unwrap_or(999.0)));
    } else {
        items.sort_by(|a, b| b.3.total_cmp(&a.3));
    }

    let sel = st.selected.min(items.len().saturating_sub(1));
    let table_rows: Vec<Row> = items
        .iter()
        .enumerate()
        .map(|(i, (name, grp, mult, shift, diff, note))| {
            let is_sel = i == sel;
            let row_style = if is_sel {
                Style::default().bg(SURFACE_ALT).fg(ACCENT_GREEN).add_modifier(Modifier::BOLD)
            } else if name.contains("Residual") || name.contains("Water") || name.contains("Grease") {
                Style::default().fg(ACCENT_AMBER)
            } else {
                Style::default().fg(TEXT_PRIMARY)
            };
            let diff_str = diff.map(|d| format!("±{:.2}", d)).unwrap_or_else(|| "—".into());
            Row::new(vec![
                Cell::from(*name),
                Cell::from(*grp),
                Cell::from(*mult),
                Cell::from(format!("{:.2} ppm", shift)),
                Cell::from(diff_str),
                Cell::from(*note),
            ])
            .style(row_style)
        })
        .collect();

    let header = Row::new(vec!["Impurity / Solvent", "Group", "Mult", "δ (ppm)", "Δppm", "Polymer Lab Note"])
        .style(Style::default().fg(ACCENT_AMBER).add_modifier(Modifier::BOLD))
        .bottom_margin(1);

    let table = Table::new(
        table_rows,
        [
            Constraint::Length(26),
            Constraint::Length(14),
            Constraint::Length(10),
            Constraint::Length(12),
            Constraint::Length(8),
            Constraint::Min(25),
        ],
    )
    .header(header);

    frame.render_widget(table, rows[2]);
}

// ============================================================================
// TAB 5: POLYMER CHEMISTRY LAB CALCULATORS
// ============================================================================
fn draw_polymer_lab_tab(frame: &mut Frame, app: &App, area: Rect) {
    let outer = panel_block("Polymer Chemistry Toolkit [h/l:Switch Tool · j/k:Field · Enter:Edit · c:Copy]", BORDER_LEFT);
    let inner = outer.inner(area);
    frame.render_widget(outer, area);

    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(2), Constraint::Min(10)])
        .split(inner);

    let pl = &app.polymer_lab;
    let sub_tabs = TextLine::from(vec![
        Span::raw(" "),
        tab_span("1. FTIR Polymerization Conversion", pl.sub_tab == PolymerSubTab::Conversion),
        Span::raw("  "),
        tab_span("2. ¹H NMR End-Group M_n / DP_n", pl.sub_tab == PolymerSubTab::NmrMn),
        Span::raw("  "),
        tab_span("3. PMMA / Vinyl Tacticity & ρ", pl.sub_tab == PolymerSubTab::Tacticity),
        Span::raw("  "),
        tab_span("4. PVDF F(β) & PU H-Bond", pl.sub_tab == PolymerSubTab::PhaseIndex),
        Span::raw("  "),
        tab_span("5. UV-Vis Tauc & Beer-Lambert", pl.sub_tab == PolymerSubTab::UvVisLab),
    ]);
    frame.render_widget(Paragraph::new(sub_tabs), rows[0]);

    let split = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(48), Constraint::Percentage(52)])
        .split(rows[1]);

    let field_line = |idx: usize, label: &str, val: &str| -> TextLine {
        let active = pl.field_idx == idx;
        let editing = active && pl.editing;
        TextLine::from(vec![
            Span::styled(
                if active { " ▸ " } else { "   " },
                Style::default().fg(ACCENT_GREEN),
            ),
            Span::styled(format!("{:<26}", label), Style::default().fg(TEXT_PRIMARY)),
            Span::styled(
                if editing { format!(" {}█ ", val) } else { format!(" {} ", val) },
                Style::default()
                    .fg(if editing { DARK_TEXT } else if active { ACCENT_YELLOW } else { TEXT_PRIMARY })
                    .bg(if editing { ACCENT_GREEN } else { SURFACE_DARK })
                    .add_modifier(if active { Modifier::BOLD } else { Modifier::empty() }),
            ),
        ])
    };

    match pl.sub_tab {
        PolymerSubTab::Conversion => {
            let mut left = vec![
                TextLine::from(Span::styled(
                    format!(
                        " Mode [m]: {} ",
                        if pl.manual_mode { "Manual Absorbance Inputs" } else { "Auto-Pull from Loaded Trace 1 (t=0) & Trace 2 (t)" }
                    ),
                    Style::default().fg(ACCENT_AQUA).add_modifier(Modifier::BOLD),
                )),
                TextLine::from(""),
            ];
            if pl.manual_mode {
                left.push(field_line(0, "A_reactive (t=0, C=C):", &pl.man_rxn_0));
                left.push(field_line(1, "A_reference (t=0, C=O):", &pl.man_ref_0));
                left.push(field_line(2, "A_reactive (t, C=C):", &pl.man_rxn_t));
                left.push(field_line(3, "A_reference (t, C=O):", &pl.man_ref_t));
            } else {
                left.push(field_line(0, "Reactive Peak Center (cm⁻¹):", &pl.rxn_center));
                left.push(field_line(1, "Reactive Half-Width (±cm⁻¹):", &pl.rxn_half));
                left.push(field_line(2, "Ref Peak Center (cm⁻¹):", &pl.ref_center));
                left.push(field_line(3, "Ref Half-Width (±cm⁻¹):", &pl.ref_half));
                left.push(TextLine::from(""));
                left.push(TextLine::from(Span::styled(
                    " Tip: Load monomer as Trace 1 & polymer as Trace 2",
                    Style::default().fg(TEXT_MUTED),
                )));
                left.push(TextLine::from(Span::styled(
                    " (default demo already has MMA + PMMA loaded!).",
                    Style::default().fg(TEXT_MUTED),
                )));
            }
            frame.render_widget(Paragraph::new(left).block(panel_block("Inputs", BORDER_LEFT)), split[0]);

            let mut right = vec![
                TextLine::from(Span::styled(
                    " Formula: Conv(%) = [1 − (A_rxn/A_ref)_t / (A_rxn/A_ref)_0] × 100",
                    Style::default().fg(TEXT_MUTED),
                )),
                TextLine::from(""),
            ];
            if pl.manual_mode {
                if let Some(conv) = app.compute_manual_conversion() {
                    right.push(TextLine::from(Span::styled(
                        format!("  Conversion = {:.2} %", conv),
                        Style::default().fg(ACCENT_GREEN).add_modifier(Modifier::BOLD),
                    )));
                } else {
                    right.push(TextLine::from(Span::styled(
                        "  Enter positive numeric absorbances.",
                        Style::default().fg(ACCENT_CORAL),
                    )));
                }
            } else if let Some(res) = app.compute_live_conversion() {
                right.push(TextLine::from(Span::styled(
                    format!("  Trace 1 (t=0): {}", app.traces[0].name),
                    Style::default().fg(ACCENT_GREEN),
                )));
                right.push(TextLine::from(Span::styled(
                    format!("  Trace 2 (t)  : {}", app.traces[1].name),
                    Style::default().fg(ACCENT_AMBER),
                )));
                right.push(TextLine::from(""));
                right.push(TextLine::from(vec![
                    Span::styled("  Integrated Area Conversion : ", Style::default().fg(TEXT_PRIMARY)),
                    Span::styled(
                        format!(" {:.2} % ", res.area_conversion_pct),
                        Style::default().fg(DARK_TEXT).bg(ACCENT_GREEN).add_modifier(Modifier::BOLD),
                    ),
                ]));
                right.push(TextLine::from(vec![
                    Span::styled("  Peak Height Conversion     : ", Style::default().fg(TEXT_PRIMARY)),
                    Span::styled(
                        format!(" {:.2} % ", res.height_conversion_pct),
                        Style::default().fg(ACCENT_YELLOW).add_modifier(Modifier::BOLD),
                    ),
                ]));
                right.push(TextLine::from(""));
                right.push(TextLine::from(Span::styled(
                    format!("  Area ratios   : R₀ = {:.4}  →  Rₜ = {:.4}", res.ratio_0_area, res.ratio_t_area),
                    Style::default().fg(TEXT_MUTED),
                )));
                right.push(TextLine::from(Span::styled(
                    format!("  Height ratios : R₀ = {:.4}  →  Rₜ = {:.4}", res.ratio_0_height, res.ratio_t_height),
                    Style::default().fg(TEXT_MUTED),
                )));
                right.push(TextLine::from(Span::styled(
                    "  (Beer-Lambert %T → Absorbance applied automatically)",
                    Style::default().fg(ACCENT_AQUA),
                )));
            } else {
                right.push(TextLine::from(Span::styled(
                    "  Need ≥ 2 FTIR traces loaded (or press [m] for Manual mode).",
                    Style::default().fg(ACCENT_CORAL),
                )));
            }
            frame.render_widget(Paragraph::new(right).block(panel_block("Result [c:Copy]", BORDER_RIGHT)), split[1]);
        }
        PolymerSubTab::NmrMn => {
            let preset_name = MONOMER_PRESETS[pl.monomer_preset].0;
            let left = vec![
                TextLine::from(Span::styled(
                    format!(" Monomer Preset [p]: {}", preset_name),
                    Style::default().fg(ACCENT_AQUA).add_modifier(Modifier::BOLD),
                )),
                TextLine::from(""),
                field_line(0, "M₀ repeat unit (g/mol):", &pl.m0),
                field_line(1, "M_end group (g/mol):", &pl.m_end),
                field_line(2, "n_repeat (H per repeat):", &pl.n_repeat),
                field_line(3, "n_end (H per end-group):", &pl.n_end),
                field_line(4, "I_repeat (integral):", &pl.i_repeat),
                field_line(5, "I_end (integral):", &pl.i_end),
            ];
            frame.render_widget(Paragraph::new(left).block(panel_block("Inputs", BORDER_LEFT)), split[0]);

            let mut right = vec![
                TextLine::from(Span::styled(
                    " DP_n = (I_repeat / n_repeat) / (I_end / n_end)",
                    Style::default().fg(TEXT_MUTED),
                )),
                TextLine::from(Span::styled(
                    " M_n  = DP_n · M₀ + M_end",
                    Style::default().fg(TEXT_MUTED),
                )),
                TextLine::from(""),
            ];
            if let Some((dpn, mn)) = app.compute_nmr_mn() {
                right.push(TextLine::from(vec![
                    Span::styled("  Degree of Polymerization (DP_n): ", Style::default().fg(TEXT_PRIMARY)),
                    Span::styled(
                        format!(" {:.1} ", dpn),
                        Style::default().fg(ACCENT_YELLOW).add_modifier(Modifier::BOLD),
                    ),
                ]));
                right.push(TextLine::from(vec![
                    Span::styled("  Number-Average M_n             : ", Style::default().fg(TEXT_PRIMARY)),
                    Span::styled(
                        format!(" {:.0} g/mol ({:.2} kDa) ", mn, mn / 1000.0),
                        Style::default().fg(DARK_TEXT).bg(ACCENT_GREEN).add_modifier(Modifier::BOLD),
                    ),
                ]));
                right.push(TextLine::from(""));
                right.push(TextLine::from(Span::styled(
                    "  Note: Valid for controlled polymerization (RAFT/ATRP/ROP)",
                    Style::default().fg(TEXT_MUTED),
                )));
                right.push(TextLine::from(Span::styled(
                    "  up to ~25–30 kDa where end-group S/N remains accurate.",
                    Style::default().fg(TEXT_MUTED),
                )));
            } else {
                right.push(TextLine::from(Span::styled(
                    "  Enter positive numeric values.",
                    Style::default().fg(ACCENT_CORAL),
                )));
            }
            frame.render_widget(Paragraph::new(right).block(panel_block("Result [c:Copy]", BORDER_RIGHT)), split[1]);
        }
        PolymerSubTab::Tacticity => {
            let left = vec![
                TextLine::from(Span::styled(
                    " PMMA α-CH₃ Triads in CDCl₃: mm ~1.21 | mr ~1.02 | rr ~0.85 ppm",
                    Style::default().fg(ACCENT_AQUA),
                )),
                TextLine::from(""),
                field_line(0, "Integral I(mm) isotactic:", &pl.i_mm),
                field_line(1, "Integral I(mr) heterotactic:", &pl.i_mr),
                field_line(2, "Integral I(rr) syndiotactic:", &pl.i_rr),
            ];
            frame.render_widget(Paragraph::new(left).block(panel_block("Triad Integrals", BORDER_LEFT)), split[0]);

            let mut right = Vec::new();
            if let Some(t) = app.compute_tacticity() {
                let bar_w = 34usize;
                let w_mm = ((t.pct_mm / 100.0) * bar_w as f64).round() as usize;
                let w_mr = ((t.pct_mr / 100.0) * bar_w as f64).round() as usize;
                let w_rr = bar_w.saturating_sub(w_mm + w_mr);

                right.push(TextLine::from(vec![
                    Span::styled("  mm (isotactic)    : ", Style::default().fg(ACCENT_CORAL)),
                    Span::styled(format!("{:5.1} %", t.pct_mm), Style::default().fg(TEXT_PRIMARY).add_modifier(Modifier::BOLD)),
                ]));
                right.push(TextLine::from(vec![
                    Span::styled("  mr (heterotactic) : ", Style::default().fg(ACCENT_AMBER)),
                    Span::styled(format!("{:5.1} %", t.pct_mr), Style::default().fg(TEXT_PRIMARY).add_modifier(Modifier::BOLD)),
                ]));
                right.push(TextLine::from(vec![
                    Span::styled("  rr (syndiotactic) : ", Style::default().fg(ACCENT_GREEN)),
                    Span::styled(format!("{:5.1} %", t.pct_rr), Style::default().fg(TEXT_PRIMARY).add_modifier(Modifier::BOLD)),
                ]));
                right.push(TextLine::from(""));
                right.push(TextLine::from(vec![
                    Span::raw("  ["),
                    Span::styled("█".repeat(w_mm), Style::default().fg(ACCENT_CORAL)),
                    Span::styled("█".repeat(w_mr), Style::default().fg(ACCENT_AMBER)),
                    Span::styled("█".repeat(w_rr), Style::default().fg(ACCENT_GREEN)),
                    Span::raw("]"),
                ]));
                right.push(TextLine::from(""));
                right.push(TextLine::from(Span::styled(
                    format!("  Dyads: m = {:.3}   r = {:.3}   (P_m = {:.1}%)", t.dyad_m, t.dyad_r, t.dyad_m * 100.0),
                    Style::default().fg(ACCENT_AQUA),
                )));
                right.push(TextLine::from(vec![
                    Span::styled("  Bernoullian persistence ρ = 2mr/[mr]: ", Style::default().fg(TEXT_PRIMARY)),
                    Span::styled(
                        format!("{:.2}", t.rho),
                        Style::default().fg(ACCENT_YELLOW).add_modifier(Modifier::BOLD),
                    ),
                ]));
                right.push(TextLine::from(Span::styled(
                    format!("  → {}", t.regime),
                    Style::default().fg(ACCENT_GREEN),
                )));
            } else {
                right.push(TextLine::from(Span::styled(
                    "  Enter non-negative triad integrals.",
                    Style::default().fg(ACCENT_CORAL),
                )));
            }
            frame.render_widget(Paragraph::new(right).block(panel_block("Microstructure [c:Copy]", BORDER_RIGHT)), split[1]);
        }
        PolymerSubTab::PhaseIndex => {
            let left = vec![
                TextLine::from(Span::styled(
                    format!(
                        " Mode [m]: {} ",
                        if pl.phase_manual {
                            "Manual Absorbance Inputs"
                        } else {
                            "Auto-Pull from Active FTIR Trace (baseline-corrected Abs)"
                        }
                    ),
                    Style::default().fg(ACCENT_AQUA).add_modifier(Modifier::BOLD),
                )),
                TextLine::from(""),
                field_line(0, "A_840 (β-phase, ±15 cm⁻¹):", &pl.a840),
                field_line(1, "A_763 (α-phase, ±15 cm⁻¹):", &pl.a763),
                field_line(2, "A_1703 (bonded C=O, ±15):", &pl.a1703),
                field_line(3, "A_1732 (free C=O, ±15):", &pl.a1732),
                TextLine::from(""),
                TextLine::from(Span::styled(
                    " Gregorio: F(β) = A_840 / (1.26·A_763 + A_840)",
                    Style::default().fg(TEXT_MUTED),
                )),
                TextLine::from(Span::styled(
                    " PU: HBI = A_1703 / A_1732; Bonded% = HBI/(1+HBI)",
                    Style::default().fg(TEXT_MUTED),
                )),
            ];
            frame.render_widget(
                Paragraph::new(left).block(panel_block("Inputs", BORDER_LEFT)),
                split[0],
            );

            let mut right = Vec::new();
            match app.compute_pvdf_beta() {
                Some((a840, a763, f_beta)) => {
                    let bar_w = 34usize;
                    let w_beta = ((f_beta / 100.0) * bar_w as f64).round() as usize;
                    let w_alpha = bar_w.saturating_sub(w_beta);
                    right.push(TextLine::from(vec![
                        Span::styled("  PVDF β-phase F(β) : ", Style::default().fg(TEXT_PRIMARY)),
                        Span::styled(
                            format!(" {:5.1} % ", f_beta),
                            Style::default()
                                .fg(DARK_TEXT)
                                .bg(ACCENT_GREEN)
                                .add_modifier(Modifier::BOLD),
                        ),
                    ]));
                    right.push(TextLine::from(vec![
                        Span::raw("  ["),
                        Span::styled("█".repeat(w_beta), Style::default().fg(ACCENT_GREEN)),
                        Span::styled("█".repeat(w_alpha), Style::default().fg(ACCENT_CORAL)),
                        Span::raw("]"),
                    ]));
                    right.push(TextLine::from(Span::styled(
                        format!("  A_840 = {:.4}   A_763 = {:.4}   (piezoelectric β)", a840, a763),
                        Style::default().fg(ACCENT_AQUA),
                    )));
                }
                None => right.push(TextLine::from(Span::styled(
                    "  Load an FTIR trace or switch to Manual mode.",
                    Style::default().fg(ACCENT_CORAL),
                ))),
            }

            right.push(TextLine::from(""));
            match app.compute_pu_hbi() {
                Some(h) => {
                    let bar_w = 34usize;
                    let w_bond = ((h.bonded_pct / 100.0) * bar_w as f64).round() as usize;
                    let w_free = bar_w.saturating_sub(w_bond);
                    right.push(TextLine::from(vec![
                        Span::styled("  PU H-bonding index : ", Style::default().fg(TEXT_PRIMARY)),
                        Span::styled(
                            format!(" {:.3} ", h.hbi),
                            Style::default()
                                .fg(DARK_TEXT)
                                .bg(ACCENT_AMBER)
                                .add_modifier(Modifier::BOLD),
                        ),
                    ]));
                    right.push(TextLine::from(vec![
                        Span::styled("  Bonded C=O        : ", Style::default().fg(TEXT_PRIMARY)),
                        Span::styled(
                            format!("{:5.1} %", h.bonded_pct),
                            Style::default().fg(ACCENT_YELLOW).add_modifier(Modifier::BOLD),
                        ),
                    ]));
                    right.push(TextLine::from(vec![
                        Span::raw("  ["),
                        Span::styled("█".repeat(w_bond), Style::default().fg(ACCENT_AMBER)),
                        Span::styled("█".repeat(w_free), Style::default().fg(ACCENT_AQUA)),
                        Span::raw("]"),
                    ]));
                    right.push(TextLine::from(Span::styled(
                        format!(
                            "  A_1703 = {:.4}   A_1732 = {:.4}",
                            h.a_bonded, h.a_free
                        ),
                        Style::default().fg(ACCENT_AQUA),
                    )));
                }
                None => right.push(TextLine::from(Span::styled(
                    "  Enter A_1732 > 0 to compute the H-bonding index.",
                    Style::default().fg(ACCENT_CORAL),
                ))),
            }
            frame.render_widget(
                Paragraph::new(right).block(panel_block("Phase Indices [c:Copy]", BORDER_RIGHT)),
                split[1],
            );
        }
        PolymerSubTab::UvVisLab => {
            let tauc_lbl = if pl.tauc_direct { "Direct" } else { "Indirect" };
            let left = vec![
                TextLine::from(Span::styled(
                    " UV-Vis Laboratory [t:Tauc · m:Manual]",
                    Style::default().fg(ACCENT_AQUA).add_modifier(Modifier::BOLD),
                )),
                TextLine::from(""),
                field_line(0, "Absorbance (A):", &pl.uv_abs),
                field_line(1, "Molar ε (L·mol⁻¹·cm⁻¹):", &pl.uv_eps),
                field_line(2, "Path (cm):", &pl.uv_path_cm),
                TextLine::from(""),
                TextLine::from(Span::styled(
                    format!("  Tauc mode: {}  (toggle with [t])", tauc_lbl),
                    Style::default().fg(TEXT_MUTED),
                )),
                TextLine::from(Span::styled(
                    format!(
                        "  Manual absorbance input: {}  (toggle with [m])",
                        if pl.uv_manual { "ON" } else { "OFF" }
                    ),
                    Style::default().fg(TEXT_MUTED),
                )),
            ];
            frame.render_widget(Paragraph::new(left).block(panel_block("Inputs", BORDER_LEFT)), split[0]);

            let mut right = vec![
                TextLine::from(Span::styled(
                    " Beer–Lambert: c = A / (ε · l)",
                    Style::default().fg(TEXT_MUTED),
                )),
                TextLine::from(""),
            ];
            if let Some(conc) = app.compute_beer_lambert() {
                let (unit_label, unit_val) = if conc < 1e-6 {
                    ("nM", conc * 1e9)
                } else if conc < 1e-3 {
                    ("µM", conc * 1e6)
                } else if conc < 1.0 {
                    ("mM", conc * 1e3)
                } else {
                    ("M", conc)
                };
                right.push(TextLine::from(Span::styled(
                    format!("  Concentration = {:.3e} M  ({:.2} {})", conc, unit_val, unit_label),
                    Style::default().fg(ACCENT_GREEN).add_modifier(Modifier::BOLD),
                )));
            } else {
                right.push(TextLine::from(Span::styled(
                    "  Enter A, ε and path to compute concentration.",
                    Style::default().fg(ACCENT_CORAL),
                )));
            }
            if let Some(tauc) = app.compute_tauc() {
                right.push(TextLine::from(""));
                right.push(TextLine::from(Span::styled(
                    format!(
                        "  Tauc band gap: {:.3} eV  (edge ≈ {:.0} nm, R² = {:.3})",
                        tauc.band_gap_ev, tauc.edge_nm, tauc.r_squared
                    ),
                    Style::default().fg(ACCENT_AMBER).add_modifier(Modifier::BOLD),
                )));
            }
            frame.render_widget(
                Paragraph::new(right).block(panel_block("Results [c:Copy]", BORDER_RIGHT)),
                split[1],
            );
        }
    }
}

// ============================================================================
// FOOTER BAR
// ============================================================================
fn draw_footer(frame: &mut Frame, app: &App, area: Rect) {
    let y_status = app
        .processed
        .get(app.active_trace)
        .map(|pt| pt.y_unit.label().to_string())
        .unwrap_or_else(|| "Int".to_string());
    let flags = TextLine::from(vec![
        Span::styled(" ● ", Style::default().fg(ACCENT_GREEN)),
        Span::styled(app.status_msg.clone(), Style::default().fg(TEXT_PRIMARY)),
        Span::styled("   │ ", Style::default().fg(TEXT_MUTED)),
        Span::styled(
            format!(
                "Y:{}  View:{}  Base:{}  Smooth:{}  Norm:{}  Stack:{}",
                y_status,
                if app.hd_mode { "HD" } else { "Braille" },
                if app.process_cfg.baseline { "ON" } else { "OFF" },
                app.process_cfg.smooth.label(),
                if app.process_cfg.normalize { "ON" } else { "OFF" },
                if app.stacked { "ON" } else { "OFF" },
            ),
            Style::default().fg(ACCENT_AMBER),
        ),
    ]);

    let hints = TextLine::from(Span::styled(
        " Click/Drag:Cursor/Zoom | Wheel/+/-:Zoom | h/l:Move | g:GoTo | v:HD/Dots | a:%T↔Abs | b:Base | m:Smooth | w:Norm | s:Stack | p:Peaks | M:Library | n/N:Jump | x:Pin | c:Copy | e:Export",
        Style::default().fg(TEXT_MUTED),
    ));

    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(1), Constraint::Length(1)])
        .split(area);
    frame.render_widget(Paragraph::new(flags), rows[0]);
    frame.render_widget(Paragraph::new(hints), rows[1]);
}

// ============================================================================
// MODALS (GoTo, Search, FilePicker, DemoPicker, Publication Export, Help)
// ============================================================================
fn centered_rect(pct_x: u16, pct_y: u16, area: Rect) -> Rect {
    let v = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage((100 - pct_y) / 2),
            Constraint::Percentage(pct_y),
            Constraint::Percentage((100 - pct_y) / 2),
        ])
        .split(area);
    Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage((100 - pct_x) / 2),
            Constraint::Percentage(pct_x),
            Constraint::Percentage((100 - pct_x) / 2),
        ])
        .split(v[1])[1]
}

fn draw_modal(frame: &mut Frame, app: &App, area: Rect) {
    match &app.modal {
        Modal::None => {}
        Modal::GoTo(buf) => {
            let rect = centered_rect(42, 22, area);
            frame.render_widget(Clear, rect);
            let unit = app.current_kind().x_unit();
            let lines = vec![
                TextLine::from(""),
                TextLine::from(Span::styled(
                    format!("  Enter exact {} value ({}):", app.current_kind().symbol(), unit),
                    Style::default().fg(TEXT_MUTED),
                )),
                TextLine::from(""),
                TextLine::from(Span::styled(
                    format!("  > {}█", buf),
                    Style::default().fg(ACCENT_GREEN).add_modifier(Modifier::BOLD),
                )),
                TextLine::from(""),
                TextLine::from(Span::styled("  [Enter] Jump   [Esc] Cancel", Style::default().fg(TEXT_MUTED))),
            ];
            frame.render_widget(
                Paragraph::new(lines).block(modal_block("Go To Position", ACCENT_YELLOW)),
                rect,
            );
        }
        Modal::DemoPicker(sel) => {
            let rect = centered_rect(58, 42, area);
            frame.render_widget(Clear, rect);
            let mut lines = vec![
                TextLine::from(Span::styled(
                    "  Choose a built-in Polymer / Organic Reference Spectrum:",
                    Style::default().fg(TEXT_MUTED),
                )),
                TextLine::from(""),
            ];
            for (i, preset) in DemoPreset::ALL.iter().enumerate() {
                let active = i == *sel;
                lines.push(TextLine::from(Span::styled(
                    format!("  {} {}", if active { "▶" } else { " " }, preset.label()),
                    if active {
                        Style::default().fg(DARK_TEXT).bg(ACCENT_GREEN).add_modifier(Modifier::BOLD)
                    } else {
                        Style::default().fg(TEXT_PRIMARY)
                    },
                )));
            }
            lines.push(TextLine::from(""));
            lines.push(TextLine::from(Span::styled(
                "  [j/k] Navigate   [Enter] Load   [Esc] Close",
                Style::default().fg(TEXT_MUTED),
            )));
            frame.render_widget(
                Paragraph::new(lines).block(modal_block("Built-in Demo Presets", ACCENT_AQUA)),
                rect,
            );
        }
        Modal::Search(st) => {
            let rect = centered_rect(72, 68, area);
            frame.render_widget(Clear, rect);
            let mut lines = vec![
                TextLine::from(vec![
                    Span::styled("  Search: ", Style::default().fg(ACCENT_AMBER).add_modifier(Modifier::BOLD)),
                    Span::styled(format!("{}█", st.query), Style::default().fg(ACCENT_GREEN).add_modifier(Modifier::BOLD)),
                    Span::styled(format!("  ({} matches)", st.results.len()), Style::default().fg(TEXT_MUTED)),
                ]),
                TextLine::from(Span::styled(
                    "  ────────────────────────────────────────────────────────────",
                    Style::default().fg(TEXT_MUTED),
                )),
            ];
            let max_vis = rect.height.saturating_sub(6) as usize;
            let start = st.selected.saturating_sub(max_vis.saturating_sub(1));
            for (idx, item) in st.results.iter().enumerate().skip(start).take(max_vis) {
                let active = idx == st.selected;
                lines.push(TextLine::from(vec![
                    Span::styled(
                        format!(" [{}] ", item.category),
                        Style::default().fg(if active { DARK_TEXT } else { ACCENT_AQUA }).bg(if active { ACCENT_GREEN } else { SURFACE_DARK }),
                    ),
                    Span::styled(
                        format!(" {}", item.title),
                        Style::default()
                            .fg(if active { ACCENT_GREEN } else { TEXT_PRIMARY })
                            .add_modifier(if active { Modifier::BOLD } else { Modifier::empty() }),
                    ),
                ]));
                if active && !item.detail.is_empty() {
                    lines.push(TextLine::from(Span::styled(
                        format!("      {}", item.detail),
                        Style::default().fg(TEXT_MUTED),
                    )));
                }
            }
            frame.render_widget(
                Paragraph::new(lines).block(modal_block("Global Search (FTIR / NMR / Polymers / Solvents)", ACCENT_YELLOW)),
                rect,
            );
        }
        Modal::FilePicker(st) => {
            let rect = centered_rect(74, 72, area);
            frame.render_widget(Clear, rect);
            let total = st.all_entries.len();
            let mut lines = vec![
                TextLine::from(vec![
                    Span::styled("  Dir: ", Style::default().fg(TEXT_MUTED)),
                    Span::styled(st.dir.display().to_string(), Style::default().fg(ACCENT_AQUA)),
                ]),
                TextLine::from(vec![
                    Span::styled("  Load Mode [Ctrl+A]: ", Style::default().fg(TEXT_MUTED)),
                    Span::styled(
                        if st.append_mode { "APPEND OVERLAY" } else { "REPLACE WORKSPACE" },
                        Style::default().fg(ACCENT_AMBER).add_modifier(Modifier::BOLD),
                    ),
                    Span::styled(
                        format!("   {} entries on disk", total),
                        Style::default().fg(TEXT_MUTED),
                    ),
                ]),
            ];
            if let Some(buf) = &st.path_input {
                lines.push(TextLine::from(vec![
                    Span::styled(
                        "  Path: ",
                        Style::default()
                            .fg(ACCENT_AMBER)
                            .add_modifier(Modifier::BOLD),
                    ),
                    Span::styled(
                        format!("{}█", buf),
                        Style::default()
                            .fg(ACCENT_GREEN)
                            .add_modifier(Modifier::BOLD),
                    ),
                    Span::styled(
                        "   [Enter] Go / open   [Esc] Back to search",
                        Style::default().fg(TEXT_MUTED),
                    ),
                ]));
            } else {
                let counter = if st.filter.trim().is_empty() {
                    String::new()
                } else {
                    format!("   {}/{} shown", st.entries.len(), total)
                };
                lines.push(TextLine::from(vec![
                    Span::styled(
                        "  Filter: ",
                        Style::default()
                            .fg(ACCENT_AMBER)
                            .add_modifier(Modifier::BOLD),
                    ),
                    Span::styled(
                        format!("{}█", st.filter),
                        Style::default()
                            .fg(ACCENT_GREEN)
                            .add_modifier(Modifier::BOLD),
                    ),
                    Span::styled(counter, Style::default().fg(TEXT_MUTED)),
                ]));
            }
            lines.push(TextLine::from(""));

            let max_vis = rect.height.saturating_sub(7) as usize;
            let start = if st.selected >= max_vis {
                st.selected + 1 - max_vis
            } else {
                0
            };
            for (idx, entry) in st.entries.iter().enumerate().skip(start).take(max_vis) {
                let active = idx == st.selected;
                lines.push(TextLine::from(Span::styled(
                    format!("  {} {}", if active { "▶" } else { " " }, entry.display),
                    if active {
                        Style::default().fg(DARK_TEXT).bg(ACCENT_GREEN).add_modifier(Modifier::BOLD)
                    } else if entry.is_dir {
                        Style::default().fg(ACCENT_AQUA)
                    } else {
                        Style::default().fg(TEXT_PRIMARY)
                    },
                )));
            }
            if st.entries.is_empty() {
                lines.push(TextLine::from(Span::styled(
                    "  No files match — keep typing or [Esc] clears the search",
                    Style::default().fg(ACCENT_CORAL),
                )));
            }
            lines.push(TextLine::from(Span::styled(
                "  type=search · ↑/↓ move · Enter open · ⌫ up · / path · Ctrl+A append · Esc clear",
                Style::default().fg(TEXT_MUTED),
            )));
            frame.render_widget(
                Paragraph::new(lines).block(modal_block(
                    "Open Spectrum Data or Image — live search (.csv/.jdx/.png/.jpg)",
                    ACCENT_AQUA,
                )),
                rect,
            );
        }
        Modal::ImportConfig(st) => {
            let rect = centered_rect(68, if st.is_image { 56 } else { 42 }, area);
            frame.render_widget(Clear, rect);
            let file_name = st
                .path
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("file");

            let kind_labels: &[(&str, usize)] = &[
                ("[1] FTIR %T", 0),
                ("[2] FTIR Abs", 1),
                ("[3] ¹H NMR", 2),
                ("[4] ¹³C NMR", 3),
                ("[5] Auto", 4),
                ("[6] UV-Vis Abs", 5),
                ("[7] UV-Vis %T", 6),
            ];

            let mut kind_spans = vec![
                Span::styled(
                    if st.field_idx == 0 { "  ▶ " } else { "    " },
                    Style::default().fg(ACCENT_GREEN),
                ),
                Span::styled(format!("{:<16}", "Spectrum Type:"), Style::default().fg(TEXT_MUTED)),
            ];
            for &(lbl, idx) in kind_labels {
                let sel = st.kind_choice == idx;
                kind_spans.push(Span::styled(
                    format!(" {} ", lbl),
                    if sel {
                        Style::default()
                            .fg(DARK_TEXT)
                            .bg(if st.field_idx == 0 { ACCENT_GREEN } else { ACCENT_AMBER })
                            .add_modifier(Modifier::BOLD)
                    } else {
                        Style::default().fg(TEXT_PRIMARY).bg(SURFACE_ALT)
                    },
                ));
                kind_spans.push(Span::raw(" "));
            }

            let row_line = |idx: usize, label: &str, val: String, editable: bool| -> TextLine {
                let active = st.field_idx == idx;
                TextLine::from(vec![
                    Span::styled(
                        if active { "  ▶ " } else { "    " },
                        Style::default().fg(ACCENT_GREEN),
                    ),
                    Span::styled(format!("{:<16}", label), Style::default().fg(TEXT_MUTED)),
                    Span::styled(
                        if active && editable {
                            format!(" {}█ ", val)
                        } else {
                            format!(" {} ", val)
                        },
                        Style::default()
                            .fg(if active { DARK_TEXT } else { TEXT_PRIMARY })
                            .bg(if active { ACCENT_AMBER } else { SURFACE_ALT })
                            .add_modifier(if active { Modifier::BOLD } else { Modifier::empty() }),
                    ),
                ])
            };

            let mut lines = vec![
                TextLine::from(vec![
                    Span::styled("  File: ", Style::default().fg(TEXT_MUTED)),
                    Span::styled(file_name, Style::default().fg(TEXT_PRIMARY).add_modifier(Modifier::BOLD)),
                    Span::raw("  "),
                    Span::styled(
                        if st.is_image { " [Image Digitizer] " } else { " [Data File] " },
                        Style::default()
                            .fg(DARK_TEXT)
                            .bg(if st.is_image { ACCENT_PURPLE } else { ACCENT_AQUA })
                            .add_modifier(Modifier::BOLD),
                    ),
                ]),
                TextLine::from(""),
                TextLine::from(kind_spans),
                row_line(
                    1,
                    "Load Mode:",
                    if st.append_mode {
                        "Append Overlay".into()
                    } else {
                        "Replace Workspace".into()
                    },
                    false,
                ),
            ];

            if st.is_image {
                let bg_lbl = match st.bg_idx {
                    1 => "Light (Paper/White)",
                    2 => "Dark (Screen/Black)",
                    _ => "Auto-Detect",
                };
                lines.push(row_line(2, "X Left (axis):", st.x_left.clone(), true));
                lines.push(row_line(3, "X Right (axis):", st.x_right.clone(), true));
                lines.push(row_line(4, "Background:", bg_lbl.into(), false));
            }

            lines.push(TextLine::from(""));
            lines.push(TextLine::from(Span::styled(
                "  [↑/↓/Tab] Field   [←/→/Space] Toggle   [1..7] Quick Type   [Enter] Import   [Esc] Cancel",
                Style::default().fg(TEXT_MUTED),
            )));

            if let Some(ref err) = st.error {
                lines.push(TextLine::from(""));
                lines.push(TextLine::from(Span::styled(
                    format!("  ✖ {}", err),
                    Style::default().fg(ACCENT_CORAL).add_modifier(Modifier::BOLD),
                )));
            }

            frame.render_widget(
                Paragraph::new(lines).block(modal_block("Import Spectrum Configuration", ACCENT_GREEN)),
                rect,
            );
        }
        Modal::CandidatePicker(st) => {
            let rect = centered_rect(84, 66, area);
            frame.render_widget(Clear, rect);
            let mut lines = vec![
                TextLine::from(Span::styled(
                    "  Several library spectra match your query",
                    Style::default().fg(ACCENT_AQUA),
                )),
                TextLine::from(Span::styled(
                    "  Pick one with ↑/↓ or j/k — modes are never mixed",
                    Style::default().fg(TEXT_MUTED),
                )),
                TextLine::from(""),
            ];
            if st.entries.is_empty() {
                lines.push(TextLine::from(Span::styled(
                    "  No candidates.",
                    Style::default().fg(ACCENT_CORAL),
                )));
            }
            let max_vis = (rect.height.saturating_sub(8) as usize).max(1);
            let start = if st.selected >= max_vis {
                st.selected + 1 - max_vis
            } else {
                0
            };
            for (i, entry) in st
                .entries
                .iter()
                .enumerate()
                .skip(start)
                .take(max_vis)
            {
                let active = i == st.selected;
                let kind_color = match entry.kind {
                    SpectrumKind::Ftir => ACCENT_AQUA,
                    SpectrumKind::Hnmr => ACCENT_GREEN,
                    SpectrumKind::Cnmr => ACCENT_PURPLE,
                    SpectrumKind::UvVis => ACCENT_CORAL,
                };
                let source = entry
                    .path
                    .parent()
                    .and_then(|p| p.file_name())
                    .and_then(|s| s.to_str())
                    .unwrap_or("library");
                let name: String = entry.stem.replace('_', " ");
                lines.push(TextLine::from(vec![
                    Span::styled(
                        if active { "  ▶ " } else { "    " },
                        Style::default().fg(ACCENT_YELLOW),
                    ),
                    Span::styled(
                        format!("{:<8}", format!("[{}]", entry.kind.display_name())),
                        Style::default().fg(kind_color).add_modifier(Modifier::BOLD),
                    ),
                    Span::styled(
                        format!(" {:<32}", name.chars().take(32).collect::<String>()),
                        Style::default()
                            .fg(if active { DARK_TEXT } else { TEXT_PRIMARY })
                            .bg(if active { ACCENT_AMBER } else { SURFACE_ALT }),
                    ),
                    Span::styled(
                        format!("  {}/", source),
                        Style::default().fg(TEXT_MUTED),
                    ),
                ]));
            }
            lines.push(TextLine::from(""));
            if st.entries.len() > max_vis {
                lines.push(TextLine::from(Span::styled(
                    format!(
                        "  {}/{}  (scroll with ↑/↓)",
                        st.selected + 1,
                        st.entries.len()
                    ),
                    Style::default().fg(TEXT_MUTED),
                )));
            }
            lines.push(TextLine::from(Span::styled(
                "  [↑/↓ or j/k] Move   [Enter] Open   [Esc] Dismiss",
                Style::default().fg(TEXT_MUTED),
            )));
            frame.render_widget(
                Paragraph::new(lines)
                    .block(modal_block("Select Reference Spectrum", ACCENT_YELLOW)),
                rect,
            );
        }
        Modal::LibraryMatch(st) => {
            let rect = centered_rect(86, 74, area);
            frame.render_widget(Clear, rect);
            let is_uv = app.current_kind() == SpectrumKind::UvVis;
            let (sub1, sub2) = if is_uv {
                (
                    "  UV-Vis Hit Quality Index (HQI) — absorbance + derivative correlation",
                    "  200–800 nm grid (2 nm step) · Built-in + library/uvvis/ references",
                )
            } else {
                (
                    "  FTIR Hit Quality Index (HQI) — absorbance + derivative correlation",
                    "  CO₂ purge region 2280–2410 cm⁻¹ excluded · 811-point grid 520–3760 cm⁻¹",
                )
            };
            let mut lines = vec![
                TextLine::from(Span::styled(sub1, Style::default().fg(ACCENT_AQUA))),
                TextLine::from(Span::styled(sub2, Style::default().fg(TEXT_MUTED))),
                TextLine::from(Span::styled(
                    "  Press [/] to search all library files by name",
                    Style::default().fg(ACCENT_YELLOW),
                )),
                TextLine::from(""),
            ];

            if st.results.is_empty() {
                lines.push(TextLine::from(Span::styled(
                    "  No matches found. Load an FTIR or UV-Vis trace and try again.",
                    Style::default().fg(ACCENT_CORAL),
                )));
            }

            // Two rendered lines per hit (score row + name row): scroll the
            // window so the j/k selection is always visible.
            let max_vis = (rect.height.saturating_sub(10) as usize / 2).max(1);
            let start = st.selected.saturating_sub(max_vis.saturating_sub(1));
            for (i, m) in st
                .results
                .iter()
                .enumerate()
                .skip(start)
                .take(max_vis)
            {
                let active = i == st.selected;
                let bar_w = 20usize;
                let filled = ((m.hqi_pct / 100.0) * bar_w as f64).round() as usize;
                let color = if m.hqi_pct >= 88.0 {
                    ACCENT_GREEN
                } else if m.hqi_pct >= 70.0 {
                    ACCENT_YELLOW
                } else {
                    ACCENT_CORAL
                };

                let mut spans = vec![
                    Span::styled(
                        if active { " ▶ " } else { "   " },
                        Style::default().fg(ACCENT_GREEN),
                    ),
                    Span::styled(
                        format!("{:<12}", m.short),
                        Style::default()
                            .fg(if active { DARK_TEXT } else { TEXT_PRIMARY })
                            .bg(if active { ACCENT_GREEN } else { SURFACE_DARK })
                            .add_modifier(Modifier::BOLD),
                    ),
                    Span::raw(" "),
                    Span::styled("█".repeat(filled), Style::default().fg(color)),
                    Span::styled("░".repeat(bar_w - filled), Style::default().fg(TEXT_MUTED)),
                    Span::styled(
                        format!(" {:5.1}%", m.hqi_pct),
                        Style::default().fg(color).add_modifier(Modifier::BOLD),
                    ),
                    Span::raw("  "),
                    Span::styled(
                        format!("[{}]", m.source),
                        Style::default()
                            .fg(if m.is_builtin { ACCENT_AQUA } else { ACCENT_AMBER }),
                    ),
                ];
                if active {
                    spans.push(Span::styled(
                        format!("  r_abs {:.3} · r_deriv {:.3}", m.r_abs, m.r_deriv),
                        Style::default().fg(ACCENT_YELLOW),
                    ));
                }
                lines.push(TextLine::from(spans));
                lines.push(TextLine::from(Span::styled(
                    format!("     {}", m.name),
                    Style::default().fg(TEXT_MUTED),
                )));
            }

            lines.push(TextLine::from(""));
            if st.results.len() > max_vis {
                lines.push(TextLine::from(Span::styled(
                    format!(
                        "  {}/{} matches shown — keep pressing [j] to scroll",
                        st.selected + 1,
                        st.results.len()
                    ),
                    Style::default().fg(ACCENT_YELLOW),
                )));
            }
            lines.push(TextLine::from(Span::styled(
                "  [j/k] Navigate   [Enter] Overlay Reference   [/] Search by name   [Esc] Close",
                Style::default().fg(TEXT_MUTED),
            )));

            frame.render_widget(
                Paragraph::new(lines).block(modal_block("Spectral Library Matcher [M]", ACCENT_GREEN)),
                rect,
            );
        }
        Modal::Export(st) => {
            let rect = centered_rect(68, 56, area);
            frame.render_widget(Clear, rect);
            let row_line = |idx: usize, label: &str, val: String, editable: bool| -> TextLine {
                let active = st.field_idx == idx;
                TextLine::from(vec![
                    Span::styled(if active { "  ▶ " } else { "    " }, Style::default().fg(ACCENT_GREEN)),
                    Span::styled(format!("{:<18}", label), Style::default().fg(TEXT_MUTED)),
                    Span::styled(
                        if active && editable { format!(" {}█ ", val) } else { format!(" {} ", val) },
                        Style::default()
                            .fg(if active { DARK_TEXT } else { TEXT_PRIMARY })
                            .bg(if active { ACCENT_AMBER } else { SURFACE_ALT })
                            .add_modifier(if active { Modifier::BOLD } else { Modifier::empty() }),
                    ),
                ])
            };
            let mut lines = vec![
                TextLine::from(Span::styled(
                    "  Powered by Python Matplotlib — Publication-Grade Reversed-Axis Figures",
                    Style::default().fg(ACCENT_AQUA),
                )),
                TextLine::from(""),
                row_line(0, "Format [Space]:", st.format.label().into(), false),
                row_line(1, "Journal Style:", st.style.label().into(), false),
                row_line(2, "Layout:", if st.stacked { "Stacked (Offset)" } else { "Overlay" }.into(), false),
                row_line(3, "X-Axis Range:", if st.use_zoom_range { "Current Zoom Window" } else { "Full Spectrum" }.into(), false),
                row_line(4, "Peak Labels:", st.peak_mode.label().into(), false),
                row_line(5, "Figure Title:", st.title.clone(), true),
                row_line(6, "Output Path:", st.output_path.clone(), true),
                TextLine::from(""),
                TextLine::from(Span::styled(
                    "  [↑/↓] Field   [Space/←/→] Toggle   [Enter] Export Figure   [Esc] Close",
                    Style::default().fg(TEXT_MUTED),
                )),
            ];
            if let Some(ref msg) = st.feedback {
                lines.push(TextLine::from(""));
                lines.push(TextLine::from(Span::styled(
                    format!("  ✦ {}", msg),
                    Style::default().fg(ACCENT_GREEN).add_modifier(Modifier::BOLD),
                )));
            }
            frame.render_widget(
                Paragraph::new(lines).block(modal_block("Publication Export (Python / Matplotlib)", ACCENT_AMBER)),
                rect,
            );
        }
        Modal::Busy(busy) => {
            let rect = centered_rect(56, 18, area);
            frame.render_widget(Clear, rect);
            let elapsed = busy.started.elapsed().as_secs_f64();
            let frame_chars = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏'];
            let spin = frame_chars[(elapsed * 10.0) as usize % frame_chars.len()];
            let lines = vec![
                TextLine::from(""),
                TextLine::from(Span::styled(
                    format!("   {}  {}", spin, busy.label),
                    Style::default().fg(ACCENT_GREEN).add_modifier(Modifier::BOLD),
                )),
                TextLine::from(""),
                TextLine::from(Span::styled(
                    format!("   elapsed {:.1}s  ·  Esc hides this panel", elapsed),
                    Style::default().fg(TEXT_MUTED),
                )),
            ];
            frame.render_widget(
                Paragraph::new(lines).block(modal_block("Working…", ACCENT_AMBER)),
                rect,
            );
        }
        Modal::Help => {
            let rect = centered_rect(88, 84, area);
            let block = modal_block(
                format!(
                    "SPECTRA v{} · Beginner's Guide — What It Is & How To Use It",
                    env!("CARGO_PKG_VERSION")
                ),
                ACCENT_GREEN,
            );
            let inner = block.inner(rect);
            frame.render_widget(Clear, rect);
            frame.render_widget(block, rect);
            if inner.height < 3 {
                return;
            }

            let lines = wrap_help_lines(help_guide_lines(), inner.width as usize);
            let total = lines.len();
            // Reserve one blank row and one footer row inside the block.
            let list_h = inner.height.saturating_sub(2) as usize;
            let start = app.help_scroll.min(total.saturating_sub(list_h));
            let visible: Vec<TextLine> = lines.into_iter().skip(start).take(list_h).collect();
            frame.render_widget(
                Paragraph::new(visible),
                Rect {
                    height: list_h as u16,
                    ..inner
                },
            );

            let footer = "  ↑/↓ or j/k scroll · PgUp/PgDn · Home/End · Esc/?/Enter close";
            frame.render_widget(
                Paragraph::new(TextLine::from(Span::styled(
                    footer,
                    Style::default().fg(TEXT_MUTED),
                ))),
                Rect {
                    y: inner.y + inner.height - 1,
                    height: 1,
                    ..inner
                },
            );
        }
    }
}

/// Wrap guide rows to the available terminal width, preserving each row's
/// style while preventing text from being clipped in narrower terminals.
fn wrap_help_lines(lines: Vec<TextLine<'static>>, width: usize) -> Vec<TextLine<'static>> {
    let width = width.max(1);
    let mut wrapped = Vec::new();
    for line in lines {
        let style = line
            .spans
            .first()
            .map(|span| span.style)
            .unwrap_or_default();
        let text = line
            .spans
            .into_iter()
            .map(|span| span.content.into_owned())
            .collect::<String>();
        if text.trim().is_empty() {
            wrapped.push(TextLine::from(""));
            continue;
        }

        let indent = text.chars().take_while(|c| *c == ' ').count().min(width - 1);
        let mut current = " ".repeat(indent);
        for word in text.split_whitespace() {
            let word_len = word.chars().count();
            if current.trim().is_empty() {
                current.push_str(word);
            } else if current.chars().count() + 1 + word_len <= width {
                current.push(' ');
                current.push_str(word);
            } else {
                wrapped.push(TextLine::from(Span::styled(current, style)));
                current = " ".repeat(indent);
                current.push_str(word);
            }
        }
        wrapped.push(TextLine::from(Span::styled(current, style)));
    }
    wrapped
}

/// Content of the `?` guide.
fn help_guide_lines() -> Vec<TextLine<'static>> {
    let section = |title: &'static str| {
        TextLine::from(Span::styled(
            title,
            Style::default()
                .fg(ACCENT_AMBER)
                .add_modifier(Modifier::BOLD),
        ))
    };
    let credit = TextLine::from(Span::styled(
        "   Created by Muhammad Musavi",
        Style::default()
            .fg(ACCENT_GREEN)
            .add_modifier(Modifier::BOLD),
    ));
    vec![
        section("  WHAT IS SPECTRA?"),
        TextLine::from("   An offline terminal workbench for molecular spectroscopy: FTIR,"),
        TextLine::from("   ¹H NMR, ¹³C NMR and UV-Vis. It reads instrument exports, draws interactive"),
        TextLine::from("   spectra, labels peaks and functional groups, and identifies unknown"),
        TextLine::from("   samples against a built-in polymer database plus 1,800+ reference"),
        TextLine::from("   spectra in library/. It also calculates polymerization conversion,"),
        TextLine::from("   NMR molar mass, tacticity and phase indices."),
        TextLine::from("   Core analysis runs locally; only optional reference downloads need network."),
        TextLine::from(""),
        section("  THE SCREEN AT A GLANCE"),
        TextLine::from("   Top bar     Tabs [1..6] and quick actions; colour marks the active tab."),
        TextLine::from("   Plot        Reversed IUPAC axis (FTIR 4000→400 cm⁻¹, NMR ppm→0); UV-Vis uses ascending nm (190→900 left→right)."),
        TextLine::from("               Mouse-drag a box to zoom, scroll to zoom around the cursor,"),
        TextLine::from("               right click or [0]/[r] to reset."),
        TextLine::from("   Inspector   Cursor position, spectral zone and functional-group notes."),
        TextLine::from("               UV-Vis: shows photon energy (eV) + absorbed/observed color."),
        TextLine::from("   Peaks       Auto-detected peaks: n/N jump, x pins, [ ] threshold."),
        TextLine::from("   Traces      Loaded spectra: Tab = active, Space = hide, Del = remove,"),
        TextLine::from("               s = stack, w = normalize heights."),
        TextLine::from("   Status bar  One-line feedback for every action you take."),
        TextLine::from(""),
        section("  THE SIX TABS"),
        TextLine::from("   1 FTIR         Wavenumber cm⁻¹, %T or absorbance, band assignments,"),
        TextLine::from("                  MMA→PMMA conversion and the HQI library matcher [M]."),
        TextLine::from("   2 ¹H NMR       Chemical shifts (ppm), solvent/impurity markers,"),
        TextLine::from("                  tacticity (mm/mr/rr) and polymer proton assignments."),
        TextLine::from("   3 ¹³C NMR      Carbon shifts, CDCl₃ triplet and polymer assignments."),
        TextLine::from("   4 UV-Vis       Wavelength nm, absorbance/%T, Tauc band gap, Beer–Lambert"),
        TextLine::from("                  (auto-scaling M/mM/µM/nM), photon energy (eV), color info."),
        TextLine::from("                  Matched against UV-Vis library [M] (PhotochemCAD + NIST)."),
        TextLine::from("   5 Solvents     Searchable solvent/impurity database with editable"),
        TextLine::from("                  tolerance — separates real peaks from contamination."),
        TextLine::from("   6 Polymer Lab  Conversion, NMR Mn/DPn, tacticity, PVDF F(β), PU H-bond."),
        TextLine::from("                  These tools support analysis; verify results against your method."),
        TextLine::from(""),
        section("  MOVING AROUND"),
        TextLine::from("   ←/→ or h/l       Move the cursor (Shift = coarse jump)"),
        TextLine::from("   Mouse            Click to place cursor; drag a rectangle to zoom"),
        TextLine::from("                    Scroll zooms around the pointer; right-click resets view."),
        TextLine::from("   g                Jump to an exact cm⁻¹ / ppm value"),
        TextLine::from("   + / = / z        Zoom in; - / Z zoom out around the cursor"),
        TextLine::from("   0 or r           Reset zoom to the full spectrum"),
        TextLine::from("   Tab              Cycle loaded traces in the spectrum tabs"),
        TextLine::from("   ?                This guide — ↑/↓ or PgUp/PgDn scroll, Esc closes"),
        TextLine::from(""),
        section("  SIGNAL PROCESSING (SPECTRUM TABS)"),
        TextLine::from("   a                FTIR %Transmittance ↔ Absorbance (Beer-Lambert)"),
        TextLine::from("   b                Rubberband baseline correction"),
        TextLine::from("   m                Savitzky-Golay smoothing: OFF → 5 → 9 → 15 points"),
        TextLine::from("   w                Min-max normalization (helps when overlaying traces)"),
        TextLine::from("   s                Overlay ↔ stacked-offset layout"),
        TextLine::from("   v                HD Kitty graphics ↔ portable Braille rendering"),
        TextLine::from("   p                Show/hide auto peak markers"),
        TextLine::from("   [ / ]            Lower / raise peak detection threshold (2–60%)"),
        TextLine::from(""),
        section("  UV-VIS SPECIFIC (Tab 4)"),
        TextLine::from("   Ascending X-axis (190→900 nm left→right) unlike FTIR/NMR reversed axes."),
        TextLine::from("   Inspector shows: photon energy (eV), absorbed region, solution color."),
        TextLine::from("   Beer-Lambert: concentration auto-scales M / mM / µM / nM."),
        TextLine::from("   Tauc plot: direct/indirect band gap (eV) with R² fit quality."),
        TextLine::from("   Library matcher [M]: HQI on 200–800 nm grid (2 nm step)."),
        TextLine::from("   Export: visible-light rainbow strip (380–750 nm) at bottom."),
        TextLine::from("   Digitize image: --kind uvvis (ascending wavelengths out)."),
        TextLine::from("   Ref library: spectra fetch pack --uvvis (PhotochemCAD + NIST)."),
        TextLine::from(""),
        section("  PEAKS & NMR HELPERS"),
        TextLine::from("   n / N            Jump to the next / previous detected peak"),
        TextLine::from("   x                Pin or unpin a peak at the cursor"),
        TextLine::from("   c                Copy the cursor assignment to the clipboard"),
        TextLine::from("   S                Cycle the NMR reference solvent (CDCl₃, DMSO, D₂O…)"),
        TextLine::from("   i                Toggle solvent / impurity markers (¹H/¹³C NMR only)"),
        TextLine::from("   Solvents tab     n or Tab switches ¹H/¹³C lists; ←/→ solvent,"),
        TextLine::from("                    j/k select a row; s/Enter searches; t tolerance; x clears."),
        TextLine::from("   Polymer Lab     h/l switches calculator; j/k selects an input; Enter edits,"),
        TextLine::from("                    type numbers then Enter; m changes mode; c copies result."),
        TextLine::from("                    NMR Mn: p cycles monomer preset; Phase Index: m manual/auto."),
        TextLine::from("                    UV-Vis Lab: t toggles Tauc direct/indirect; m manual A."),
        TextLine::from("   Global          Ctrl+C quits too; Tab cycles loaded traces in spectrum tabs."),
        TextLine::from(""),
        section("  FILES, LIBRARY & SEARCH"),
        TextLine::from("   f / F            File browser (replace workspace / append overlay)"),
        TextLine::from("   In the browser   type to live-filter the current folder; ↑/↓ move,"),
        TextLine::from("                    Enter opens, ⌫ goes up, / types a path, Ctrl+A append,"),
        TextLine::from("                    Ctrl+U clears, Esc clears the search then closes."),
        TextLine::from("                    The last folder you visited is remembered."),
        TextLine::from("   d                Built-in demo presets (FTIR, ¹H, ¹³C, UV-Vis examples)"),
        TextLine::from("   /                Global search: files, library spectra, polymers,"),
        TextLine::from("                    solvents and tools in one place."),
        TextLine::from("   M                Library matcher: ranks references by Hit Quality"),
        TextLine::from("                    Index (HQI). FTIR: 520–3760 cm⁻¹; UV-Vis: 200–800 nm."),
        TextLine::from("                    Top 25 hits; Enter overlays the chosen one."),
        TextLine::from("   e                Publication export: ACS/RSC/Everforest, stacked or"),
        TextLine::from("                    overlay, PNG/SVG/PDF up to 600 DPI (via Python)."),
        TextLine::from("                    UV-Vis export includes visible-light rainbow strip."),
        TextLine::from("   Import image     Pick a PNG/JPG scan to digitize a printed spectrum:"),
        TextLine::from("                    set the X-axis bounds and background, then Enter."),
        TextLine::from("                    --kind uvvis for UV-Vis scans (ascending λ output)."),
        TextLine::from(""),
        section("  COMMAND LINE (NO TUI)"),
        TextLine::from("   spectra FILE...                 Open data files directly"),
        TextLine::from("   spectra butane                  Look up a name in the library"),
        TextLine::from("   spectra ftir butane             Strict mode search (also nmr / cnmr)"),
        TextLine::from("   spectra nmr                     Open the ¹H tab with a demo preset"),
        TextLine::from("   spectra export FILE -o out.png  Headless publication figure"),
        TextLine::from("   spectra fetch nmr-search NAME   Download extra NMR references (network)"),
        TextLine::from("   spectra fetch uvvis CAS         Download NIST UV-Vis by CAS (e.g. 71-43-2)"),
        TextLine::from("   spectra fetch pack --uvvis      Download full UV-Vis library (PhotochemCAD + NIST)"),
        TextLine::from("   spectra init-samples            Write example CSV spectra"),
        TextLine::from("   Optional fetch/export tools need Python; fetched references need network."),
        TextLine::from("   spectra --help                  Full usage summary"),
        TextLine::from(""),
        section("  BEGINNER TIPS"),
        TextLine::from("   • Press [d] first: the labelled demo spectra need no files."),
        TextLine::from("   • Move the cursor — the Inspector explains each spectral region."),
        TextLine::from("   • In FTIR, press [a] for absorbance; conversion math needs it."),
        TextLine::from("   • Run [M] on a real FTIR trace to see how library matching works."),
        TextLine::from("   • [f] then a few letters is the fastest way to find a data file."),
        TextLine::from("   • Everything is offline; extra references can be fetched later."),
        TextLine::from("   • UV-Vis: press [t] in Polymer Lab to toggle Tauc direct/indirect."),
        TextLine::from("   • UV-Vis: [M] matches against 189 reference spectra (PhotochemCAD + NIST)."),
        TextLine::from("   • UV-Vis: export adds a visible-light color bar (380–750 nm)."),
        TextLine::from(""),
        section("  CREDITS"),
        TextLine::from("   SPECTRA is written mostly in Rust (ratatui + crossterm for the TUI);"),
        TextLine::from("   a small amount of Python is used only for optional Matplotlib"),
        TextLine::from("   figures and image digitizing. Free, offline and terminal-first."),
        credit,
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{DemoPreset, SpectrumKind};
    use ratatui::backend::TestBackend;
    use ratatui::Terminal;

    #[test]
    fn beginner_guide_wraps_and_includes_credit_and_language_note() {
        let guide = wrap_help_lines(help_guide_lines(), 44);
        assert!(guide.iter().all(|line| {
            line.spans
                .iter()
                .map(|span| span.content.chars().count())
                .sum::<usize>()
                <= 44
        }));
        let text = guide
            .iter()
            .flat_map(|line| line.spans.iter())
            .map(|span| span.content.as_ref())
            .collect::<String>();
        assert!(text.contains("Created by Muhammad Musavi"));
        assert!(text.contains("Rust"));
        assert!(text.contains("WHAT IS SPECTRA?"));
    }

    #[test]
    fn selected_peak_scrolls_into_view() {
        let mut app = App::new(Vec::new(), Some(SpectrumKind::Hnmr));
        app.load_demo(DemoPreset::PmmaHnmr);

        let target = app
            .auto_peaks
            .iter()
            .position(|p| (p.x - 0.85).abs() < 0.03)
            .expect("0.85 ppm rr triad must be detected");
        app.selected_peak = target;

        let backend = TestBackend::new(180, 40);
        let mut terminal = Terminal::new(backend).unwrap();
        terminal.draw(|frame| render(frame, &mut app)).unwrap();

        let buffer = terminal.backend().buffer();
        let width = buffer.area.width as usize;
        let cells = buffer.content();
        let mut selected_line = None;
        let mut full_text = String::new();
        for row in 0..buffer.area.height as usize {
            let line: String = cells[row * width..(row + 1) * width]
                .iter()
                .map(|c| c.symbol())
                .collect();
            if line.contains("prom") && line.contains("0.85") {
                selected_line = Some(line.clone());
            }
            full_text.push_str(&line);
            full_text.push('\n');
        }

        let line = selected_line.expect("0.85 ppm peak must be rendered in the Peaks strip");
        assert!(
            line.contains('▸'),
            "selected 0.85 ppm peak must show the ▸ marker: {line:?}"
        );

        assert!(
            full_text.contains("Y:Int"),
            "NMR footer must report intensity, not transmittance"
        );
        assert!(
            full_text.contains("[q] Quit"),
            "top-bar quit hint must not be clipped"
        );
    }
}
