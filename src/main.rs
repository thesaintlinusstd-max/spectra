mod app;
mod clipboard;
mod data;
mod deconv;
mod digitize;
mod export;
mod io;
mod kitty_plot;
mod matcher;
mod mixture;
mod processing;
mod svg_export;
mod theme;
mod ui;

use std::env;
use std::io as stdio;
use std::path::PathBuf;
use std::time::Duration;

use anyhow::{bail, Result};
use crossterm::event::{
    self, DisableBracketedPaste, DisableMouseCapture, EnableBracketedPaste, EnableMouseCapture,
    Event, KeyEventKind,
};
use crossterm::execute;
use crossterm::terminal::{
    disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen,
};
use ratatui::backend::CrosstermBackend;
use ratatui::Terminal;

use crate::app::App;
use crate::data::{SpectrumKind, YUnit};
use crate::export::{find_python_script, run_export, ExportPayload, ExportPeak, ExportStyle, ExportTrace};
use crate::io::{
    load_spectrum_file, resolve_library_query, search_library_for_kind, write_sample_files,
    LibraryEntry,
};
use crate::processing::{
    detect_peaks, process_trace, BaselineMode, DerivOrder, NormMode, OpticsMode, ProcessConfig,
    SmoothLevel, YAxisPreference,
};

/// Positional CLI arguments are usually file paths, but `spectra butane` is a
/// natural thing to type. Bare words that name a technique (`nmr`, `ftir`,
/// `13c`, …) select a mode instead of becoming a spectrum query.
fn mode_keyword(s: &str) -> Option<SpectrumKind> {
    match s.trim().to_ascii_lowercase().as_str() {
        "ftir" | "ir" => Some(SpectrumKind::Ftir),
        "nmr" | "hnmr" | "1h" | "1hnmr" | "proton" => Some(SpectrumKind::Hnmr),
        "cnmr" | "13c" | "13cnmr" | "carbon" => Some(SpectrumKind::Cnmr),
        "uv" | "uvvis" | "uv-vis" => Some(SpectrumKind::UvVis),
        _ => None,
    }
}

/// What the positional CLI arguments resolved to.
#[derive(Debug, Default)]
struct CliTargets {
    /// Files that exist on disk (or resolved library paths) to open.
    files: Vec<PathBuf>,
    /// Ambiguous library hits: show the picker instead of guessing.
    candidates: Vec<LibraryEntry>,
    /// Human-readable explanation appended to the status line.
    notes: Vec<String>,
    /// Original query, for the picker title.
    query: Option<String>,
}

/// Message shown when a bare name matches nothing in the requested mode.
/// Never silently falls back to the demo preset: the user asked for a specific
/// molecule and must learn why it is not there.
fn no_match_message(query: &str, mode: Option<SpectrumKind>) -> String {
    let wanted = match mode {
        Some(kind) => format!("{} spectrum", kind.display_name()),
        None => "spectrum".to_string(),
    };
    let mut msg = format!(
        "No {} matching '{}' found in library.",
        wanted, query
    );

    let others = search_library_for_kind(query, None);
    let mut kinds: Vec<SpectrumKind> = Vec::new();
    for hit in &others {
        if Some(hit.kind) == mode || kinds.contains(&hit.kind) {
            continue;
        }
        kinds.push(hit.kind);
    }
    if !kinds.is_empty() {
        msg.push_str("\n\nMatches exist in other modes:");
        for kind in kinds.iter().take(3) {
            let example = others
                .iter()
                .find(|h| h.kind == *kind)
                .map(|h| h.stem.as_str())
                .unwrap_or(query);
            msg.push_str(&format!(
                "\n  Try: spectra {} {}    ({})",
                kind.as_str(),
                query,
                example
            ));
        }
    }
    msg.push_str(&format!(
        "\n\nFetch 58k NMR dump: spectra fetch nmr-search {}",
        query
    ));
    msg
}

/// Resolve positional arguments: existing paths open directly, bare words are
/// looked up strictly inside the selected mode.
fn resolve_cli_targets(raw: Vec<PathBuf>, mode: Option<SpectrumKind>) -> Result<CliTargets> {
    let mut files = Vec::new();
    let mut words: Vec<String> = Vec::new();
    for path in raw {
        if path.exists() {
            files.push(path);
        } else {
            words.push(path.to_string_lossy().to_string());
        }
    }

    let mut out = CliTargets {
        files,
        ..CliTargets::default()
    };

    if words.is_empty() {
        return Ok(out);
    }

    if out.files.is_empty() {
        // Every positional word is a name query: `spectra nmr butane`.
        let query = words.join(" ");
        if query.contains('/')
            || query.starts_with('~')
            || crate::io::is_spectrum_extension(std::path::Path::new(&query))
        {
            bail!("File not found: {}", query);
        }
        let hits = search_library_for_kind(&query, mode);
        if hits.is_empty() {
            bail!("{}", no_match_message(&query, mode));
        }
        let exact: Vec<&LibraryEntry> = hits
            .iter()
            .filter(|h| crate::io::exact_stem_match(&h.stem, &query))
            .collect();
        if hits.len() == 1 {
            let hit = &hits[0];
            out.notes.push(format!(
                "{} → {} ({})",
                query,
                hit.stem,
                hit.kind.display_name()
            ));
            out.files.push(hit.path.clone());
        } else if exact.len() == 1 {
            let hit = exact[0];
            out.notes.push(format!(
                "{} → {} ({})",
                query,
                hit.stem,
                hit.kind.display_name()
            ));
            out.files.push(hit.path.clone());
        } else {
            out.notes
                .push(format!("{} ambiguous: {} candidates", query, hits.len()));
            out.query = Some(query);
            out.candidates = hits;
        }
    } else {
        // Mixed invocation: real files plus names to append from the library.
        for word in words {
            match resolve_library_query(&word, mode) {
                Some(hit) => {
                    out.notes.push(format!(
                        "{} → {} ({})",
                        word,
                        hit.stem,
                        hit.kind.display_name()
                    ));
                    out.files.push(hit.path);
                }
                None => out.notes.push(no_match_message(&word, mode)),
            }
        }
    }

    Ok(out)
}

struct TerminalGuard;

impl TerminalGuard {
    fn enter() -> Result<Self> {
        enable_raw_mode()?;
        let mut stdout = stdio::stdout();
        execute!(
            stdout,
            EnterAlternateScreen,
            EnableMouseCapture,
            EnableBracketedPaste
        )?;
        Ok(Self)
    }
}

impl Drop for TerminalGuard {
    fn drop(&mut self) {
        let _ = disable_raw_mode();
        let mut stdout = stdio::stdout();
        kitty_plot::clear_kitty_image(&mut stdout);
        let _ = execute!(
            stdout,
            DisableBracketedPaste,
            DisableMouseCapture,
            LeaveAlternateScreen
        );
    }
}

fn install_panic_hook() {
    let orig = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        let _ = disable_raw_mode();
        let mut stdout = stdio::stdout();
        kitty_plot::clear_kitty_image(&mut stdout);
        let _ = execute!(
            stdout,
            DisableBracketedPaste,
            DisableMouseCapture,
            LeaveAlternateScreen
        );
        orig(info);
    }));
}

fn parse_mode(s: &str) -> Result<SpectrumKind> {
    match s.trim().to_ascii_lowercase().as_str() {
        "ftir" | "ir" => Ok(SpectrumKind::Ftir),
        "nmr" | "hnmr" | "1h" | "1hnmr" | "proton" => Ok(SpectrumKind::Hnmr),
        "cnmr" | "13c" | "13cnmr" | "carbon" => Ok(SpectrumKind::Cnmr),
        "uv" | "uvvis" | "uv-vis" => Ok(SpectrumKind::UvVis),
        other => bail!("Unknown mode '{}' (use ftir, hnmr, cnmr, or uvvis)", other),
    }
}

fn parse_style(s: &str) -> Result<ExportStyle> {
    match s.trim().to_ascii_lowercase().as_str() {
        "acs" => Ok(ExportStyle::Acs),
        "rsc" => Ok(ExportStyle::Rsc),
        "everforest" | "dark" => Ok(ExportStyle::Everforest),
        other => bail!("Unknown style '{}' (use acs, rsc, or everforest)", other),
    }
}

fn print_help() {
    println!(
        "SPECTRA {} — FTIR & NMR Spectroscopy Explorer & Polymer Chemistry Toolkit\n\n\
         USAGE:\n  \
           spectra [FILES...] [--mode ftir|hnmr|cnmr|uvvis]\n  \
           spectra [ftir|nmr|cnmr|uvvis] [QUERY...]     (e.g. spectra ftir butane)\n  \
           spectra init-samples [--dir samples]\n  \
           spectra fetch <nist|nmr|nmr-search|import|list|pack> [ARGS...]\n  \
           spectra export <FILES...> -o <OUT.svg|png|pdf> [--style acs|rsc|everforest]\n                 \
                   [--stacked] [--absorbance] [--baseline (AsLS) | --rubberband | --arpls]\n                 \
                   [--atr diamond|znse|ge|si | --km] [--deriv 1|2]\n                 \
                   [--norm minmax|snv|area] [--no-co2]\n                 \
                   [--smooth 0|5|9|15] [--title \"...\"] [--dpi 600]\n\n\
         SVG figures are rendered natively in Rust; PNG/PDF use the Python helper.\n",
        env!("CARGO_PKG_VERSION")
    );
}

#[allow(clippy::too_many_arguments)]
fn run_headless_export(
    files: Vec<PathBuf>,
    output: PathBuf,
    style: ExportStyle,
    stacked: bool,
    absorbance: bool,
    baseline: bool,
    rubberband: bool,
    arpls: bool,
    optics: OpticsMode,
    deriv: DerivOrder,
    norm: NormMode,
    atm_co2: bool,
    smooth: usize,
    title: String,
    dpi: u32,
) -> Result<()> {
    if files.is_empty() {
        bail!("spectra export requires at least one input spectrum file");
    }
    let smooth_lvl = match smooth {
        5 => SmoothLevel::Win5,
        9 => SmoothLevel::Win9,
        15 => SmoothLevel::Win15,
        _ => SmoothLevel::Off,
    };
    let cfg = ProcessConfig {
        y_pref: if absorbance {
            YAxisPreference::Swapped
        } else {
            YAxisPreference::Native
        },
        baseline: if rubberband {
            BaselineMode::Rubberband
        } else if arpls {
            BaselineMode::Arpls
        } else if baseline {
            BaselineMode::Asls
        } else {
            BaselineMode::Off
        },
        smooth: smooth_lvl,
        norm,
        optics,
        deriv,
        atm_co2,

    };

    let mut traces = Vec::new();
    let mut kind = SpectrumKind::Ftir;
    let mut y_unit = YUnit::Transmittance;
    let mut x_lo = f64::INFINITY;
    let mut x_hi = f64::NEG_INFINITY;
    let mut peaks_out = Vec::new();

    for (idx, path) in files.iter().enumerate() {
        let spec = load_spectrum_file(path, None)?;
        let (x, y, yu) = process_trace(&spec, cfg);
        if idx == 0 {
            kind = spec.kind;
            y_unit = yu;
            for p in detect_peaks(&x, &y, spec.kind, yu, 0.10, 8) {
                peaks_out.push(ExportPeak {
                    x: p.x,
                    y: p.y,
                    label: p.label,
                    pinned: false,
                });
            }
        }
        if let (Some(&a), Some(&b)) = (x.first(), x.last()) {
            x_lo = x_lo.min(a);
            x_hi = x_hi.max(b);
        }
        traces.push(ExportTrace {
            name: spec.name,
            x,
            y,
        });
    }

    let payload = ExportPayload {
        mode: kind.as_str().to_string(),
        y_unit: y_unit.as_str().to_string(),
        style: style.as_str().to_string(),
        stacked,
        dpi,
        title,
        x_bounds: [x_lo, x_hi],
        active_trace: 0,
        traces,
        peaks: peaks_out,
    };

    let saved = run_export(&payload, &output)?;
    println!("Exported publication figure: {}", saved.display());
    Ok(())
}

fn run_tui(targets: CliTargets, mode: Option<SpectrumKind>) -> Result<()> {
    install_panic_hook();
    let _guard = TerminalGuard::enter()?;
    let backend = CrosstermBackend::new(stdio::stdout());
    let mut terminal = Terminal::new(backend)?;
    // No `terminal.clear()` here: ratatui's clear() round-trips a cursor-position
    // query (`ESC[6n`) and aborts after ~2 s on terminals that do not answer it.
    // The alternate screen is already blank and the first draw repaints everything.

    let mut app = App::new(targets.files, mode);
    if !targets.candidates.is_empty() {
        app.open_candidates(
            targets.query.clone().unwrap_or_default(),
            targets.candidates,
        );
    }
    if !targets.notes.is_empty() {
        app.status_msg = format!("{} · {}", app.status_msg, targets.notes.join(" · "));
    }
    let mut hd_state = kitty_plot::KittyPlotState::default();

    while !app.should_quit {
        app.poll_background();
        terminal.draw(|frame| ui::render(frame, &mut app))?;
        let _ = kitty_plot::sync_hd_plot(&app, &mut hd_state, &mut stdio::stdout());

        if event::poll(Duration::from_millis(40))? {
            match event::read()? {
                Event::Key(key)
                    if key.kind == KeyEventKind::Press || key.kind == KeyEventKind::Repeat =>
                {
                    app.handle_key(key);
                }
                Event::Mouse(mouse) => app.handle_mouse(mouse),
                // Kitty / terminal drag-and-drop and Ctrl+Shift+V arrive as a
                // single Paste event: treat them as a path (or a name query)
                // instead of replaying every character as a keystroke.
                Event::Paste(text) => app.handle_paste(&text),
                _ => {}
            }
        }
    }

    Ok(())
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().skip(1).collect();

    if args.iter().any(|a| a == "-h" || a == "--help") {
        print_help();
        return Ok(());
    }
    if args.len() == 1 && (args[0] == "-V" || args[0] == "--version" || args[0] == "version") {
        println!("spectra {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }

    if args.first().map(|s| s.as_str()) == Some("init-samples") {
        let mut dir = PathBuf::from("samples");
        let mut i = 1;
        while i < args.len() {
            match args[i].as_str() {
                "-d" | "--dir" => {
                    i += 1;
                    if i < args.len() {
                        dir = PathBuf::from(&args[i]);
                    }
                }
                other if !other.starts_with('-') => {
                    dir = PathBuf::from(other);
                }
                _ => {}
            }
            i += 1;
        }
        let paths = write_sample_files(&dir)?;
        println!("Generated {} sample files:", paths.len());
        for p in paths {
            println!("  - {}", p.display());
        }
        return Ok(());
    }

    if args.first().map(|s| s.as_str()) == Some("fetch") {
        let (python, script) = find_python_script("fetch_library.py")?;
        let status = std::process::Command::new(&python)
            .arg(&script)
            .args(&args[1..])
            .status()
            .map_err(|e| anyhow::anyhow!("Failed to launch fetch_library.py: {}", e))?;
        if !status.success() {
            bail!("fetch_library.py exited with status {}", status);
        }
        return Ok(());
    }

    if args.first().map(|s| s.as_str()) == Some("export") {
        let mut files = Vec::new();
        let mut output = PathBuf::from("exports/figure.png");
        let mut style = ExportStyle::Acs;
        let mut stacked = false;
        let mut absorbance = false;
        let mut baseline = false;
        let mut rubberband = false;
        let mut arpls = false;
        let mut optics = OpticsMode::Off;
        let mut deriv = DerivOrder::Raw;
        let mut norm = NormMode::Off;
        let mut atm_co2 = false;
        let mut smooth: usize = 0;
        let mut title = String::new();
        let mut dpi: u32 = 600;

        let mut i = 1;
        while i < args.len() {
            match args[i].as_str() {
                "-o" | "--output" => {
                    i += 1;
                    if i >= args.len() {
                        bail!("Missing value after --output");
                    }
                    output = PathBuf::from(&args[i]);
                }
                "--style" => {
                    i += 1;
                    if i >= args.len() {
                        bail!("Missing value after --style");
                    }
                    style = parse_style(&args[i])?;
                }
                "--stacked" => stacked = true,
                "--absorbance" => absorbance = true,
                "--baseline" => baseline = true,
                "--rubberband" => rubberband = true,
                "--arpls" => arpls = true,
                "--km" => optics = OpticsMode::KubelkaMunk,
                "--atr" => {
                    i += 1;
                    if i >= args.len() {
                        bail!("Missing value after --atr");
                    }
                    optics = match args[i].to_ascii_lowercase().as_str() {
                        "diamond" | "dia" => OpticsMode::AtrDiamond,
                        "znse" => OpticsMode::AtrZnSe,
                        "ge" | "germanium" => OpticsMode::AtrGe,
                        "si" | "silicon" => OpticsMode::AtrSi,
                        other => bail!("Unknown ATR crystal '{}' (diamond|znse|ge|si)", other),
                    };
                }
                "--deriv" => {
                    i += 1;
                    if i >= args.len() {
                        bail!("Missing value after --deriv");
                    }
                    deriv = match args[i].as_str() {
                        "1" => DerivOrder::First,
                        "2" => DerivOrder::Second,
                        other => bail!("Unknown derivative order '{}' (1|2)", other),
                    };
                }
                "--norm" => {
                    i += 1;
                    if i >= args.len() {
                        bail!("Missing value after --norm");
                    }
                    norm = match args[i].to_ascii_lowercase().as_str() {
                        "minmax" => NormMode::MinMax,
                        "snv" | "vector" => NormMode::VectorSnv,
                        "area" => NormMode::Area,
                        "off" | "none" => NormMode::Off,
                        other => bail!("Unknown normalization '{}' (minmax|snv|area)", other),
                    };
                }
                "--no-co2" => atm_co2 = true,
                "--smooth" => {
                    i += 1;
                    if i >= args.len() {
                        bail!("Missing value after --smooth");
                    }
                    smooth = args[i].parse()?;
                }
                "--title" => {
                    i += 1;
                    if i >= args.len() {
                        bail!("Missing value after --title");
                    }
                    title = args[i].clone();
                }
                "--dpi" => {
                    i += 1;
                    if i >= args.len() {
                        bail!("Missing value after --dpi");
                    }
                    dpi = args[i].parse()?;
                }
                other if !other.starts_with('-') => files.push(PathBuf::from(other)),
                flag => bail!("Unknown export flag: {}", flag),
            }
            i += 1;
        }

        return run_headless_export(
            files, output, style, stacked, absorbance, baseline, rubberband, arpls, optics,
            deriv, norm, atm_co2, smooth, title, dpi,
        );
    }

    // Default: launch interactive TUI
    let mut files = Vec::new();
    let mut mode: Option<SpectrumKind> = None;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-m" | "--mode" | "--kind" => {
                i += 1;
                if i >= args.len() {
                    bail!("Missing value after --mode");
                }
                mode = Some(parse_mode(&args[i])?);
            }
            other if !other.starts_with('-') => {
                let path = PathBuf::from(other);
                // `spectra nmr butane` / `spectra 13c butane`: a bare mode word
                // selects the technique unless a file of that name exists.
                match (path.exists(), mode_keyword(other)) {
                    (false, Some(kind)) if mode.is_none() => mode = Some(kind),
                    _ => files.push(path),
                }
            }
            flag => bail!("Unknown option: {}", flag),
        }
        i += 1;
    }

    let targets = resolve_cli_targets(files, mode)?;
    run_tui(targets, mode)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{
        match_ftir, match_solvent_impurities, match_uvvis, mma_monomer_spectrum,
        pmma_polymer_spectrum, polystyrene_ftir_spectrum, DemoPreset, SolventId,
        SpectrumKind, Spectrum, YUnit,
    };
    use crate::io::{encode_jcamp_difdup, parse_delimited, parse_jcamp_dx};
    use crate::processing::{
        absorbance_to_transmittance, auto_ftir_conversion, beer_lambert_conc,
        compute_tauc_band_gap, nmr_mn_dpn, pu_hydrogen_bond_index, pvdf_beta_fraction,
        rubberband_baseline, savitzky_golay, tacticity_stats, transmittance_to_absorbance,
        ProcessConfig, TaucMode,
    };
    use crate::processing::{
        advanced_atr_correction, arpls_baseline, compute_spectral_math, kubelka_munk_transform,
        optimal_subtraction_factor, NormMode, SpectralOp,
    };
    use crate::deconv::{deconvolve_window, DeconvConfig, PeakShape};

    #[test]
    fn test_csv_eu_decimal_and_headers() {
        let sample = "# Instrument: Bruker\nWavenumber;Transmittance\n4000,0;98,5\n1725,5;15,2\n400,0;95,0\nEND\n";
        let spec = parse_delimited(sample, "test".into(), None).unwrap();
        assert_eq!(spec.kind, SpectrumKind::Ftir);
        assert_eq!(spec.raw_y_unit, YUnit::Transmittance);
        assert_eq!(spec.x.len(), 3);
        assert!((spec.x[0] - 400.0).abs() < 1e-6);
        assert!((spec.x[2] - 4000.0).abs() < 1e-6);
        assert!((spec.y[1] - 15.2).abs() < 1e-6);
    }

    #[test]
    fn test_jcamp_difdup_roundtrip_and_checkpoint_dedup() {
        let orig = polystyrene_ftir_spectrum();
        let jdx_text = encode_jcamp_difdup(&orig);
        let decoded = parse_jcamp_dx(&jdx_text, "ps".into(), None).unwrap();
        assert_eq!(decoded.kind, SpectrumKind::Ftir);
        assert_eq!(decoded.x.len(), orig.x.len());
        let max_err = orig
            .y
            .iter()
            .zip(&decoded.y)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!(max_err < 0.05, "Max DIFDUP round-trip error too large: {}", max_err);
    }

    #[test]
    fn test_beer_lambert_roundtrip() {
        for pct_t in [10.0, 25.0, 50.0, 88.5, 100.0] {
            let a = transmittance_to_absorbance(pct_t);
            let t2 = absorbance_to_transmittance(a);
            assert!((pct_t - t2).abs() < 1e-9);
        }
    }

    #[test]
    fn test_dsp_smoothing_baseline_and_peak_detection() {
        let pmma = pmma_polymer_spectrum();
        let sm = savitzky_golay(&pmma.y, SmoothLevel::Win9);
        assert_eq!(sm.len(), pmma.y.len());
        assert!(sm.iter().all(|v| v.is_finite()));

        let base = rubberband_baseline(&pmma.x, &sm, true);
        assert!(base.iter().all(|&v| v.is_finite() && (0.0..=100.0).contains(&v)));

        let peaks = detect_peaks(&pmma.x, &base, SpectrumKind::Ftir, YUnit::Transmittance, 0.12, 12);
        assert!(
            peaks.iter().any(|p| (p.x - 1728.0).abs() < 4.0),
            "Expected PMMA ester C=O peak near 1728 cm⁻¹"
        );
    }

    #[test]
    fn test_polymer_calculators() {
        let mono = mma_monomer_spectrum();
        let poly = pmma_polymer_spectrum();
        let conv = auto_ftir_conversion(&mono, &poly, 1638.0, 25.0, 1725.0, 35.0).unwrap();
        assert!(
            conv.area_conversion_pct > 85.0 && conv.area_conversion_pct < 98.0,
            "Unexpected MMA->PMMA conversion: {:.2}%",
            conv.area_conversion_pct
        );

        let (dpn, mn) = nmr_mn_dpn(100.12, 221.3, 3.0, 2.0, 240.0, 2.0).unwrap();
        assert!((dpn - 80.0).abs() < 1e-6);
        assert!((mn - (80.0 * 100.12 + 221.3)).abs() < 1e-3);

        let tact = tacticity_stats(6.25, 37.5, 56.25).unwrap();
        assert!((tact.rho - 1.0).abs() < 1e-4);
    }

    #[test]
    fn test_fulmer_solvent_database() {
        let cdcl3_hits = match_solvent_impurities(SolventId::Cdcl3, false, 7.26, 0.02);
        assert!(!cdcl3_hits.is_empty());
        let dmso_water = match_solvent_impurities(SolventId::Dmso, false, 3.33, 0.02);
        assert!(dmso_water.iter().any(|(e, _, _)| e.name.contains("Water")));
        let grease = match_solvent_impurities(SolventId::Cdcl3, false, 0.07, 0.02);
        assert!(grease.iter().any(|(e, _, _)| e.name.contains("Silicone")));
    }

    #[test]
    fn test_offline_json_serializer() {
        let payload = ExportPayload {
            mode: "ftir".into(),
            y_unit: "transmittance".into(),
            style: "acs".into(),
            stacked: true,
            dpi: 600,
            title: "MMA \"vs\" PMMA".into(),
            x_bounds: [400.0, 4000.0],
            active_trace: 0,
            traces: vec![ExportTrace {
                name: "MMA".into(),
                x: vec![400.0, 1725.0],
                y: vec![95.0, 12.0],
            }],
            peaks: vec![ExportPeak {
                x: 1725.0,
                y: 12.0,
                label: "1725 C=O".into(),
                pinned: true,
            }],
        };
        let js = payload.to_json_string();
        assert!(js.contains("\"title\":\"MMA \\\"vs\\\" PMMA\""));
        assert!(js.contains("\"stacked\":true"));
        assert!(js.contains("\"peaks\":[{\"x\":1725.000000"));
    }

    #[test]
    fn test_ftir_specificity_ranking() {
        let ch2 = match_ftir(2916.0);
        assert!(!ch2.is_empty(), "expected FTIR assignments at 2916 cm⁻¹");
        assert!(
            ch2[0].group.contains("CH₂"),
            "narrow CH₂ stretch must outrank broad O–H, got: {}",
            ch2[0].group
        );
        assert!(
            !ch2[0].group.contains("O–H"),
            "broad carboxylic O–H must not rank first at 2916 cm⁻¹"
        );

        let pvdf = crate::data::match_ftir_with_context(840.0, &[1402.0, 1275.0, 1175.0, 840.0, 510.0]);
        assert!(
            pvdf.iter().any(|b| b.band.group.contains("PVDF β-phase")),
            "expected PVDF β-phase assignment at 840 cm⁻¹, got {:?}",
            pvdf.iter().map(|b| b.band.group).collect::<Vec<_>>()
        );
        assert!(
            pvdf[0].band.group.contains("PVDF β-phase"),
            "PVDF β-phase should rank first at 840 cm⁻¹ with PVDF context, got: {}",
            pvdf[0].band.group
        );
        assert!(pvdf[0].confirmed, "PVDF β-phase should be confirmed in PVDF context");
    }

    #[test]
    fn test_library_matcher_identifies_pmma_and_ps() {
        let pmma = pmma_polymer_spectrum();
        let pm = crate::matcher::match_ftir_library(&pmma, ProcessConfig::default());
        assert!(!pm.is_empty(), "matcher returned no PMMA candidates");
        assert!(
            pm[0].short.to_uppercase().contains("PMMA")
                || pm[0].name.to_uppercase().contains("PMMA"),
            "PMMA should rank #1, got {} ({})",
            pm[0].short,
            pm[0].name
        );
        assert!(
            pm[0].hqi_pct > 88.0,
            "PMMA HQI too low: {:.1}% (r_abs={:.3}, r_deriv={:.3})",
            pm[0].hqi_pct, pm[0].r_abs, pm[0].r_deriv
        );

        let ps = polystyrene_ftir_spectrum();
        let psm = crate::matcher::match_ftir_library(&ps, ProcessConfig::default());
        assert!(!psm.is_empty(), "matcher returned no PS candidates");
        assert!(
            psm[0].short.to_uppercase().contains("PS")
                || psm[0].name.to_uppercase().contains("POLYSTYRENE"),
            "PS should rank #1, got {} ({})",
            psm[0].short,
            psm[0].name
        );
        assert!(
            psm[0].hqi_pct > 88.0,
            "PS HQI too low: {:.1}% (r_abs={:.3}, r_deriv={:.3})",
            psm[0].hqi_pct, psm[0].r_abs, psm[0].r_deriv
        );
    }

    #[test]
    fn test_pvdf_beta_and_pu_hbi_calculators() {
        let f_beta = pvdf_beta_fraction(1.26, 1.0).expect("valid Gregorio inputs");
        assert!(
            (f_beta - 50.0).abs() < 1e-6,
            "Gregorio F(β) expected 50.0%, got {:.4}%",
            f_beta
        );
        assert!(pvdf_beta_fraction(0.0, 0.0).is_none());

        let hbi = pu_hydrogen_bond_index(0.75, 0.25).expect("valid HBI inputs");
        assert!((hbi.hbi - 3.0).abs() < 1e-9, "HBI expected 3.0, got {:.4}", hbi.hbi);
        assert!(
            (hbi.bonded_pct - 75.0).abs() < 1e-6,
            "bonded fraction expected 75.0%, got {:.3}%",
            hbi.bonded_pct
        );
        assert!(pu_hydrogen_bond_index(0.5, 0.0).is_none());
    }

    #[test]
    fn test_extended_solvents_dmac_nmp_aibn() {
        let dmso = |ppm: f64| {
            match_solvent_impurities(SolventId::Dmso, false, ppm, 0.06)
                .into_iter()
                .map(|(e, _, _)| e.name)
                .collect::<Vec<_>>()
        };
        assert!(
            dmso(2.73).iter().any(|n| n.contains("DMF")),
            "DMF second N–CH₃ (2.73 in DMSO-d₆) not resolved"
        );
        assert!(
            dmso(2.94).iter().any(|n| n.contains("DMAc")),
            "DMAc N–CH₃ (2.94 in DMSO-d₆) not resolved"
        );
        assert!(
            dmso(3.31).iter().any(|n| n.contains("NMP")),
            "NMP NCH₂ (3.31 in DMSO-d₆) not resolved"
        );
        assert!(
            dmso(1.68).iter().any(|n| n.contains("AIBN")),
            "AIBN C(CH₃)₂ (1.68 in DMSO-d₆) not resolved"
        );

        let cdcl3 = |ppm: f64| {
            match_solvent_impurities(SolventId::Cdcl3, false, ppm, 0.06)
                .into_iter()
                .map(|(e, _, _)| e.name)
                .collect::<Vec<_>>()
        };
        assert!(
            cdcl3(2.85).iter().any(|n| n.contains("NMP")),
            "NMP NCH₃ (2.85 in CDCl₃) not resolved"
        );
        assert!(
            cdcl3(1.70).iter().any(|n| n.contains("AIBN")),
            "AIBN C(CH₃)₂ (1.70 in CDCl₃) not resolved"
        );
        assert!(
            cdcl3(3.02).iter().any(|n| n.contains("DMAc")),
            "DMAc N–CH₃ (3.02 in CDCl₃) not resolved"
        );
    }

    /// Smoke test for the real `library/ftir/` pack: every downloadable reference
    /// must parse, and querying a reference against the library must identify
    /// itself at ~100% HQI. Skips silently when the pack has not been fetched.
    #[test]
    fn test_real_library_pack_self_identification() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("library/ftir");
        let Ok(rd) = std::fs::read_dir(&dir) else {
            return;
        };
        let mut files: Vec<std::path::PathBuf> = rd
            .flatten()
            .map(|e| e.path())
            .filter(|p| p.is_file())
            .collect();
        files.sort();
        if files.is_empty() {
            return;
        }

        let mut probe: Option<(&std::path::PathBuf, crate::data::Spectrum)> = None;
        for path in &files {
            match load_spectrum_file(path, Some(SpectrumKind::Ftir)) {
                Ok(spec) => {
                    assert!(
                        spec.x.len() > 50,
                        "{} parsed to only {} points",
                        path.display(),
                        spec.x.len()
                    );
                    assert!(
                        spec.x.windows(2).all(|w| w[0] <= w[1]),
                        "{} is not sorted ascending",
                        path.display()
                    );
                    if probe.is_none() {
                        probe = Some((path, spec));
                    }
                }
                Err(err) => eprintln!("skipping unreadable library file {}: {err}", path.display()),
            }
        }

        let Some((path, spec)) = probe else {
            return;
        };
        let matches = crate::matcher::match_ftir_library(&spec, ProcessConfig::default());
        assert!(!matches.is_empty(), "library matcher returned no candidates");
        let top = &matches[0];
        assert_eq!(
            top.file_path.as_deref(),
            Some(path.as_path()),
            "self-match should rank {} first, got {} ({:.1}%)",
            path.display(),
            top.short,
            top.hqi_pct
        );
        assert!(
            top.hqi_pct > 99.0,
            "self-match HQI should be ~100%, got {:.1}%",
            top.hqi_pct
        );
    }

    /// Smoke test for the real `library/nmr/` pack produced from the cached
    /// nmrshiftdb2 SDfile: every CSV must parse, infer the right nucleus from
    /// its ppm range, carry a Lorentzian peak, and keep ascending abscissae.
    #[test]
    fn test_real_library_nmr_pack_parses() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("library/nmr");
        let Ok(rd) = std::fs::read_dir(&dir) else {
            return;
        };
        let mut files: Vec<std::path::PathBuf> = rd
            .flatten()
            .map(|e| e.path())
            .filter(|p| p.extension().and_then(|e| e.to_str()) == Some("csv"))
            .collect();
        files.sort();
        if files.is_empty() {
            return;
        }

        for path in &files {
            let stem = path.file_stem().and_then(|s| s.to_str()).unwrap_or("");
            let spec = load_spectrum_file(path, None)
                .unwrap_or_else(|e| panic!("{} failed to parse: {e}", path.display()));
            let expected = if stem.ends_with("_1h") {
                SpectrumKind::Hnmr
            } else {
                SpectrumKind::Cnmr
            };
            assert_eq!(
                spec.kind,
                expected,
                "{} inferred {:?}, expected {:?}",
                path.display(),
                spec.kind,
                expected
            );
            assert!(
                spec.x.len() >= 1000,
                "{} parsed to only {} points",
                path.display(),
                spec.x.len()
            );
            assert!(
                spec.x.windows(2).all(|w| w[0] <= w[1]),
                "{} is not sorted ascending",
                path.display()
            );
            let y_max = spec.y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            assert!(
                y_max > 90.0,
                "{} has no normalized Lorentzian peak (y_max={y_max:.1})",
                path.display()
            );
        }
    }

    #[test]
    fn test_bare_name_resolves_to_library_spectrum() {
        let hit = crate::io::resolve_library_query("butane", Some(SpectrumKind::Ftir))
            .expect("`butane` should resolve to an FTIR library file");
        assert!(
            hit.stem.to_ascii_lowercase().contains("butane"),
            "unexpected library hit: {}",
            hit.stem
        );
        assert_eq!(hit.kind, SpectrumKind::Ftir);
        let spec = load_spectrum_file(&hit.path, None).expect("resolved file should load");
        assert!(spec.x.len() > 50);
    }

    #[test]
    fn test_cli_bare_word_becomes_library_candidates() {
        let targets = resolve_cli_targets(
            vec![std::path::PathBuf::from("butane")],
            Some(SpectrumKind::Ftir),
        )
        .expect("butane exists in the FTIR library");
        // Several butane references exist (n-butane, isobutane, bromobutane…):
        // the CLI must offer the picker rather than guessing one.
        assert!(
            !targets.candidates.is_empty(),
            "notes: {:?}",
            targets.notes
        );
        assert!(
            targets.candidates.iter().all(|c| c.kind == SpectrumKind::Ftir),
            "strict mode filter leaked a non-FTIR reference"
        );
        assert!(
            targets
                .candidates
                .iter()
                .any(|c| c.stem.to_ascii_lowercase().contains("butane"))
        );
    }

    #[test]
    fn test_nmr_mode_keyword_and_strict_kind_filter() {
        assert_eq!(mode_keyword("nmr"), Some(SpectrumKind::Hnmr));
        assert_eq!(mode_keyword("HNMR"), Some(SpectrumKind::Hnmr));
        assert_eq!(mode_keyword("13c"), Some(SpectrumKind::Cnmr));
        assert_eq!(mode_keyword("ftir"), Some(SpectrumKind::Ftir));
        // Real molecule names are never swallowed as mode words.
        assert_eq!(mode_keyword("butane"), None);

        // `spectra nmr butane`: strict filtering must never return the IR file.
        let hnmr = crate::io::search_library_for_kind("butane", Some(SpectrumKind::Hnmr));
        assert!(hnmr.iter().all(|h| h.kind == SpectrumKind::Hnmr));
        assert!(
            !hnmr.iter().any(|h| h.stem.contains("org_nist_n_butane")),
            "FTIR butane leaked into the ¹H NMR results"
        );

        let err = resolve_cli_targets(
            vec![std::path::PathBuf::from("butane")],
            Some(SpectrumKind::Hnmr),
        )
        .expect_err("there is no ¹H NMR butane reference");
        let msg = format!("{err}");
        assert!(msg.contains("¹H NMR"), "{msg}");
        assert!(msg.contains("spectra ftir butane"), "{msg}");
        assert!(msg.contains("nmr-search butane"), "{msg}");
    }

    #[test]
    fn test_alkane_ch_stretch_ranks_first() {
        let hits = match_ftir(2906.7);
        assert!(!hits.is_empty(), "2906.7 cm⁻¹ matched nothing");
        let top = hits[0];
        assert!(
            top.group.contains("C–H") && top.group.contains("sp³"),
            "expected C–H stretch first, got {} ({})",
            top.group,
            top.mode
        );
    }

    #[test]
    fn test_remove_active_trace_is_safe() {
        let mut app = App::new(vec![], None);
        let n = app.traces.len();
        assert!(n > 0, "demo workspace should start with traces");
        for _ in 0..n {
            app.remove_active_trace();
        }
        assert!(app.traces.is_empty());
        assert!(app.pinned_peaks.is_empty());
        assert!(app.processed.is_empty());
        // One extra delete on an empty workspace must not panic.
        app.remove_active_trace();
        assert!(app.traces.is_empty());
    }

    #[test]
    fn test_file_picker_live_filter() {
        use crate::app::{filter_file_entries, FileEntry};
        let mk = |name: &str, is_dir: bool| FileEntry {
            path: std::path::PathBuf::from("/tmp/demo").join(name),
            display: format!("  {name}"),
            is_dir,
        };
        let entries = vec![
            mk("../ (parent directory)", true),
            mk("samples", true),
            mk("isobutane.jdx", false),
            mk("butane_13c.csv", false),
            mk("butane_1h.csv", false),
            mk("pmma_ftir.csv", false),
        ];

        let hits = filter_file_entries(&entries, "butane");
        assert_eq!(hits.len(), 3, "hits: {:?}", hits.iter().map(|e| &e.display).collect::<Vec<_>>());
        assert!(hits.iter().all(|e| !e.is_dir));

        let dirs = filter_file_entries(&entries, "sam");
        assert_eq!(dirs.len(), 1);
        assert!(dirs[0].is_dir);

        // Tokens are ANDed: "butane csv" drops the JCAMP file.
        let multi = filter_file_entries(&entries, "butane csv");
        assert_eq!(multi.len(), 2);
        assert!(multi.iter().all(|e| e.display.ends_with(".csv")));

        assert!(filter_file_entries(&entries, "zzz").is_empty());
        // Clearing the query restores the full listing, "../" included.
        assert_eq!(filter_file_entries(&entries, "").len(), entries.len());
    }

    #[test]
    fn test_last_picker_dir_roundtrip() {
        let base =
            std::env::temp_dir().join(format!("spectra-state-test-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&base);
        let target = std::env::temp_dir();
        crate::io::write_last_dir_to(&base, &target).expect("state write");
        let read = crate::io::read_last_dir_from(&base).expect("state read");
        assert_eq!(
            read.canonicalize().unwrap(),
            target.canonicalize().unwrap()
        );

        // A remembered folder that no longer exists must not be returned.
        let gone = base.join("deleted-folder");
        crate::io::write_last_dir_to(&base, &gone).expect("state write");
        assert!(crate::io::read_last_dir_from(&base).is_none());
        let _ = std::fs::remove_dir_all(&base);
    }

    #[test]
    fn test_version_is_0_7() {
        assert_eq!(env!("CARGO_PKG_VERSION"), "0.7.0");
    }

    #[test]
    fn test_uvvis_key_a_and_baseline_cycle_and_pinned_y_refresh() {
        use crossterm::event::{KeyCode, KeyEvent, KeyModifiers};

        let mut app = App::new(vec![], Some(SpectrumKind::UvVis));
        assert_eq!(app.processed[0].y_unit, YUnit::Absorbance);

        let key_a = KeyEvent::new(KeyCode::Char('a'), KeyModifiers::NONE);
        app.handle_key(key_a);
        assert_eq!(app.processed[0].y_unit, YUnit::Transmittance);

        let key_b = KeyEvent::new(KeyCode::Char('b'), KeyModifiers::NONE);
        app.handle_key(key_b);
        assert_eq!(app.process_cfg.baseline, crate::processing::BaselineMode::Rubberband);
        app.handle_key(key_b);
        assert_eq!(app.process_cfg.baseline, crate::processing::BaselineMode::Asls);
        app.handle_key(key_b);
        assert_eq!(app.process_cfg.baseline, crate::processing::BaselineMode::Arpls);
        app.handle_key(key_b);
        assert_eq!(app.process_cfg.baseline, crate::processing::BaselineMode::Off);

        // FTIR: pinned %T y-value must be refreshed to absorbance after [a].
        let mut ftir = App::new(vec![], Some(SpectrumKind::Ftir));
        ftir.load_demo(DemoPreset::MmaVsPmmaFtir);
        ftir.cursor_x = 1722.0;
        ftir.toggle_pin_at_cursor();
        assert!(
            ftir.pinned_peaks[0].y > 1.0,
            "expected %T-scale pinned y > 1.0, got {}",
            ftir.pinned_peaks[0].y
        );
        ftir.handle_key(key_a);
        assert!(
            ftir.pinned_peaks[0].y < 2.0,
            "expected absorbance-scale refreshed pinned y < 2.0, got {}",
            ftir.pinned_peaks[0].y
        );
    }

    #[test]
    fn test_dropped_image_opens_import_config() {
        let tmp = std::env::temp_dir().join(format!(
            "spectra-dropped-{}-{}.png",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .map(|d| d.as_millis())
                .unwrap_or(0)
        ));
        let img = image::RgbaImage::from_pixel(80, 60, image::Rgba([255, 255, 255, 255]));
        img.save(&tmp).expect("write tiny PNG");

        let mut app = App::new(vec![], Some(SpectrumKind::UvVis));
        app.handle_paste(&tmp.display().to_string());
        assert!(
            matches!(app.modal, crate::app::Modal::ImportConfig(ref st) if st.is_image),
            "dropped PNG should open the ImportConfig image dialog"
        );

        let _ = std::fs::remove_file(&tmp);
    }

    #[test]
    fn test_latin1_jcamp_dx_loads() {
        // Many MNHN exports are ISO-8859-1; a strict UTF-8 read used to make
        // these files look corrupt.
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("library/ftir/org_mnhn_wood_mnhn_170031_atr_0.jdx");
        if !path.exists() {
            return;
        }
        let spec = load_spectrum_file(&path, Some(SpectrumKind::Ftir))
            .expect("Latin-1 JCAMP-DX should be decoded");
        assert!(spec.x.len() > 50);
    }

    #[test]
    fn test_uvvis_kind_detection_no_ftir_conflict() {
        // UV-Vis data (nm) must be detected as UvVis, never confused with FTIR.
        let uvvis = "Wavelength (nm),Absorbance\n200,0.10\n300,0.20\n400,0.30\n500,0.20\n600,0.10\n700,0.05\n800,0.01\n";
        let spec = parse_delimited(uvvis, "uvvis".into(), None).unwrap();
        assert_eq!(spec.kind, SpectrumKind::UvVis);
        assert!(spec.x.iter().any(|&v| v > 190.0 && v < 900.0));

        // FTIR data (cm⁻¹) stays FTIR; no UV-Vis / FTIR conflict.
        let ftir = "Wavenumber,Transmittance\n4000,98\n1725,15\n400,95\n";
        let ftir_spec = parse_delimited(ftir, "ftir".into(), None).unwrap();
        assert_eq!(ftir_spec.kind, SpectrumKind::Ftir);

        // match_uvvis returns bands around a known wavelength.
        let bands = match_uvvis(260.0);
        assert!(!bands.is_empty());
    }

    #[test]
    fn test_uvvis_tauc_band_gap_and_beer_lambert() {
        // Beer–Lambert: c = A / (ε · l).
        assert!((beer_lambert_conc(0.5, 50_000.0, 1.0).unwrap() - 1e-5).abs() < 1e-12);
        assert!(beer_lambert_conc(-1.0, 50_000.0, 1.0).is_none());
        assert!(beer_lambert_conc(0.5, 0.0, 1.0).is_none());
        assert!(beer_lambert_conc(0.5, 50_000.0, 0.0).is_none());

        // Synthetic direct-allowed absorption edge at Eg ≈ 3.0 eV (≈ 413 nm).
        let mut x = Vec::new();
        let mut y = Vec::new();
        for nm in (300..=413).map(|i| i as f64) {
            let e = 1239.84193 / nm;
            let alpha = if e > 3.0 { (100.0 * (e - 3.0)).sqrt() / e } else { 0.0 };
            x.push(nm);
            y.push(alpha);
        }
        let spec = Spectrum {
            name: "uvvis-edge".into(),
            source_path: None,
            kind: SpectrumKind::UvVis,
            raw_y_unit: YUnit::Absorbance,
            x,
            y,
            visible: true,
            color_idx: 0,
            metadata: vec![],
        };
        let result = compute_tauc_band_gap(&spec, TaucMode::Direct).unwrap();
        assert!(result.band_gap_ev > 0.0);
        assert!((2.8..=3.4).contains(&result.band_gap_ev));
        assert!(result.r_squared > 0.85);

        // Non-UV-Vis kind yields no Tauc result.
        let ftir_spec = Spectrum {
            name: "ftir".into(),
            source_path: None,
            kind: SpectrumKind::Ftir,
            raw_y_unit: YUnit::Transmittance,
            x: spec.x.clone(),
            y: spec.y.clone(),
            visible: true,
            color_idx: 0,
            metadata: vec![],
        };
        assert!(compute_tauc_band_gap(&ftir_spec, TaucMode::Direct).is_none());
    }

    #[test]
    fn test_uvvis_cli_mode() {
        assert_eq!(mode_keyword("uv"), Some(SpectrumKind::UvVis));
        assert_eq!(mode_keyword("uvvis"), Some(SpectrumKind::UvVis));
        assert_eq!(mode_keyword("uv-vis"), Some(SpectrumKind::UvVis));
        assert_eq!(parse_mode("uvvis").unwrap(), SpectrumKind::UvVis);
        assert!(parse_mode("notamode").is_err());
    }

    #[test]
    fn test_uvvis_library_matching_and_inspector() {
        // Verify match_uvvis_library returns scored hits for a synthetic spectrum
        // with a peak around 260 nm (benzene-like).
        let mut x = Vec::new();
        let mut y = Vec::new();
        for nm in (190..=400).map(|i| i as f64) {
            let peak = (-(nm - 260.0).powi(2) / (2.0 * 15.0_f64.powi(2))).exp();
            x.push(nm);
            y.push(peak);
        }
        let spec = Spectrum {
            name: "test_uvvis".into(),
            source_path: None,
            kind: SpectrumKind::UvVis,
            raw_y_unit: YUnit::Absorbance,
            x,
            y,
            visible: true,
            color_idx: 0,
            metadata: vec![],
        };
        let matches = crate::matcher::match_uvvis_library(&spec, crate::processing::ProcessConfig::default());
        assert!(!matches.is_empty(), "match_uvvis_library should return hits");
        // Should contain known compounds
        let names: Vec<_> = matches.iter().map(|m| m.name.as_str()).collect();
        assert!(names.iter().any(|n| n.contains("benzene") || n.contains("Benzene")),
            "Should match benzene-like spectrum: {:?}", names);

        // Verify photon_energy_ev and uvvis_color_info work for UV-Vis range
        let ev_260 = crate::data::photon_energy_ev(260.0);
        assert!((ev_260 - 4.768).abs() < 0.01, "photon_energy_ev(260) ≈ 4.77 eV");
        let (abs_region, sol_color) = crate::data::uvvis_color_info(260.0);
        assert!(abs_region.contains("UV"));
        assert!(!sol_color.is_empty());

        let ev_500 = crate::data::photon_energy_ev(500.0);
        assert!((ev_500 - 2.479).abs() < 0.01, "photon_energy_ev(500) ≈ 2.48 eV");
        let (abs_region, sol_color) = crate::data::uvvis_color_info(500.0);
        assert!(abs_region.contains("Green"));
        assert!(sol_color.contains("Purple") || sol_color.contains("Magenta"));
    }

    #[test]
    fn test_mma_monomer_context_aware_labels_no_false_positives() {
        let mma = crate::data::mma_monomer_spectrum();
        let (xs, ys, yu) = process_trace(&mma, ProcessConfig::default());
        let peaks = detect_peaks(&xs, &ys, mma.kind, yu, 0.08, 16);
        assert!(!peaks.is_empty(), "MMA monomer should yield detected peaks");

        for p in &peaks {
            for bad in ["Polyimide", "Calcium", "CaCO", "Sulfone", "POM", "Si–CH₃"] {
                assert!(
                    !p.label.contains(bad),
                    "MMA peak {} should not be mislabeled with '{}'",
                    p.label,
                    bad
                );
            }
        }

        let p1722 = peaks
            .iter()
            .min_by(|a, b| (a.x - 1722.0).abs().total_cmp(&(b.x - 1722.0).abs()))
            .expect("expected ~1722 peak");
        assert!(
            p1722.label.contains("C=O"),
            "~1722 cm⁻¹ peak should be labeled with C=O, got: {}",
            p1722.label
        );

        let p1443 = peaks
            .iter()
            .min_by(|a, b| (a.x - 1443.0).abs().total_cmp(&(b.x - 1443.0).abs()))
            .expect("expected ~1443 peak");
        assert!(
            p1443.label.contains("CH"),
            "~1443 cm⁻¹ peak should be labeled with CH, got: {}",
            p1443.label
        );
    }

    #[test]
    fn test_confirmed_specific_materials_still_identify() {
        let pi = crate::data::match_ftir_with_context(1720.0, &[1778.0, 1720.0, 1375.0, 725.0]);
        assert!(
            pi[0].confirmed && pi[0].band.group.contains("Polyimide"),
            "Polyimide should be confirmed at 1720 cm⁻¹ with PI context, got: {} (confirmed={})",
            pi[0].band.group,
            pi[0].confirmed
        );

        let caco3 = crate::data::match_ftir_with_context(1425.0, &[1425.0, 875.0, 712.0]);
        assert!(
            caco3[0].confirmed && caco3[0].band.group.contains("Calcium Carbonate"),
            "CaCO₃ should be confirmed at 1425 cm⁻¹ with calcite context, got: {} (confirmed={})",
            caco3[0].band.group,
            caco3[0].confirmed
        );

        let pdms_795 = crate::data::match_ftir_with_context(795.0, &[2962.0, 1260.0, 1080.0, 1015.0, 795.0]);
        assert!(
            pdms_795[0].confirmed && pdms_795[0].band.group.contains("Si–CH₃"),
            "PDMS Si–CH₃ should be confirmed at 795 cm⁻¹, got: {} (confirmed={})",
            pdms_795[0].band.group,
            pdms_795[0].confirmed
        );

        let pdms_1260 = crate::data::match_ftir_with_context(1260.0, &[2962.0, 1260.0, 1080.0, 1015.0, 795.0]);
        assert!(
            pdms_1260[0].confirmed && pdms_1260[0].band.group.contains("Si–CH₃"),
            "PDMS Si–CH₃ should be confirmed at 1260 cm⁻¹, got: {} (confirmed={})",
            pdms_1260[0].band.group,
            pdms_1260[0].confirmed
        );
    }

    #[test]
    fn test_asls_recovers_curved_baseline() {
        let n = 1000;
        let mut baseline = Vec::with_capacity(n);
        let mut signal = Vec::with_capacity(n);
        let mut raw = Vec::with_capacity(n);
        for i in 0..n {
            let t = i as f64 / (n - 1) as f64;
            let b = 0.25 + 0.45 * (t - 0.5).powi(2) + 0.10 * t;
            let p1 = 0.90 / (1.0 + ((i as f64 - 300.0) / 8.0).powi(2));
            let p2 = 0.70 / (1.0 + ((i as f64 - 700.0) / 10.0).powi(2));
            let s = p1 + p2;
            baseline.push(b);
            signal.push(s);
            raw.push(b + s);
        }

        let corr = crate::processing::asls_baseline(&raw, 1e6, 0.01, 10);
        let mut max_off_err = 0.0_f64;
        for (i, &c) in corr.iter().enumerate().take(n) {
            if (i as f64 - 300.0).abs() > 45.0 && (i as f64 - 700.0).abs() > 55.0 {
                max_off_err = max_off_err.max(c.abs());
            }
        }
        assert!(
            max_off_err < 0.025,
            "AsLS off-peak residual should be < 0.025, got {:.4}",
            max_off_err
        );
        assert!(
            corr[300] > 0.90 * signal[300],
            "AsLS should preserve >90% of peak 1 height: {:.4} vs {:.4}",
            corr[300],
            signal[300]
        );
        assert!(
            corr[700] > 0.90 * signal[700],
            "AsLS should preserve >90% of peak 2 height: {:.4} vs {:.4}",
            corr[700],
            signal[700]
        );
    }

    #[test]
    fn test_tauc_ignores_deep_uv_second_band() {
        let mut x = Vec::new();
        let mut y = Vec::new();
        let eg_true = 3.20_f64;
        for nm in (200..=600).step_by(2) {
            let nm_f = nm as f64;
            let hv = 1239.84193 / nm_f;
            // Direct-allowed TiO2 onset around 3.20 eV + strong 220 nm deep-UV Gaussian band
            let edge = if hv > eg_true {
                0.75 * ((hv - eg_true) / 0.35).sqrt().min(1.3)
            } else {
                0.01 * ((hv - eg_true) / 0.12).exp()
            };
            let deep_uv = 2.4 * (-0.5 * ((nm_f - 220.0) / 14.0).powi(2)).exp();
            x.push(nm_f);
            y.push(0.02 + edge + deep_uv);
        }

        let spec = Spectrum {
            name: "tio2_deep_uv".into(),
            source_path: None,
            kind: SpectrumKind::UvVis,
            raw_y_unit: YUnit::Absorbance,
            x,
            y,
            visible: true,
            color_idx: 0,
            metadata: vec![],
        };

        let res = crate::processing::compute_tauc_band_gap(
            &spec,
            crate::processing::TaucMode::Direct,
        )
        .expect("should compute Tauc band gap");
        assert!(
            (3.10..=3.30).contains(&res.band_gap_ev),
            "Tauc Eg should be in [3.10, 3.30] eV, got {:.3} eV",
            res.band_gap_ev
        );
        assert!(
            res.r_squared > 0.95,
            "Tauc R² should be > 0.95, got {:.4}",
            res.r_squared
        );
    }

    #[test]
    fn test_detect_peaks_plateau_and_fwhm() {
        let xs: Vec<f64> = (0..9).map(|i| 100.0 + i as f64 * 10.0).collect();
        let ys = vec![0.0, 0.2, 0.5, 1.0, 1.0, 1.0, 0.5, 0.2, 0.0];
        let peaks = detect_peaks(&xs, &ys, SpectrumKind::UvVis, YUnit::Absorbance, 0.1, 5);
        assert_eq!(peaks.len(), 1, "expected 1 plateau peak, got {:?}", peaks);
        assert!(
            (peaks[0].x - 140.0).abs() < 1e-6,
            "plateau peak should center on midpoint 140.0, got {}",
            peaks[0].x
        );
        let fwhm = peaks[0].fwhm.expect("expected finite FWHM");
        assert!(fwhm.is_finite() && fwhm > 0.0, "expected positive finite FWHM, got {}", fwhm);
    }

    #[test]
    fn test_jcamp_micrometers_ir_and_hz_nmr() {
        let ir_um = "\
##TITLE=IR Micrometers Test
##JCAMP-DX=4.24
##DATA TYPE=INFRARED SPECTRUM
##XUNITS=MICROMETERS
##YUNITS=TRANSMITTANCE
##XYPOINTS=(XY..XY)
2.5 95.0; 5.0 80.0; 10.0 30.0; 15.0 90.0
##END=
";
        let ir_spec = parse_jcamp_dx(ir_um, "ir_um".into(), None).unwrap();
        assert_eq!(ir_spec.kind, SpectrumKind::Ftir);
        assert!(ir_spec.x.windows(2).all(|w| w[0] <= w[1]));
        assert!(
            (ir_spec.x[0] - 666.666667).abs() < 0.05,
            "expected ~666.67 cm⁻¹, got {}",
            ir_spec.x[0]
        );
        assert!(
            (*ir_spec.x.last().unwrap() - 4000.0).abs() < 0.05,
            "expected ~4000.0 cm⁻¹, got {}",
            ir_spec.x.last().unwrap()
        );

        let nmr_hz = "\
##TITLE=1H NMR Hz Test
##JCAMP-DX=5.00
##DATA TYPE=NMR SPECTRUM
##XUNITS=HZ
##YUNITS=ARBITRARY UNITS
##.OBSERVE FREQUENCY=400.0
##XYPOINTS=(XY..XY)
0.0 10.0; 1440.0 85.0; 4000.0 5.0
##END=
";
        let nmr_spec = parse_jcamp_dx(nmr_hz, "nmr_hz".into(), None).unwrap();
        assert_eq!(nmr_spec.kind, SpectrumKind::Hnmr);
        assert!(nmr_spec.x.windows(2).all(|w| w[0] <= w[1]));
        assert!(
            (nmr_spec.x[0] - 0.0).abs() < 1e-6,
            "expected 0.0 ppm, got {}",
            nmr_spec.x[0]
        );
        assert!(
            (*nmr_spec.x.last().unwrap() - 10.0).abs() < 1e-6,
            "expected 10.0 ppm, got {}",
            nmr_spec.x.last().unwrap()
        );

        let uv_wn = "\
##TITLE=UV-Vis Wavenumber Test
##JCAMP-DX=4.24
##DATA TYPE=UV/VIS SPECTRUM
##XUNITS=1/CM
##YUNITS=ABSORBANCE
##XYPOINTS=(XY..XY)
20000.0 0.10; 35000.0 0.80; 50000.0 0.25
##END=
";
        let uv_spec = parse_jcamp_dx(uv_wn, "uv_wn".into(), None).unwrap();
        assert_eq!(uv_spec.kind, SpectrumKind::UvVis);
        assert!(uv_spec.x.windows(2).all(|w| w[0] <= w[1]));
        assert!(
            (uv_spec.x[0] - 200.0).abs() < 1e-6,
            "expected 200.0 nm, got {}",
            uv_spec.x[0]
        );
        assert!(
            (*uv_spec.x.last().unwrap() - 500.0).abs() < 1e-6,
            "expected 500.0 nm, got {}",
            uv_spec.x.last().unwrap()
        );
    }

    #[test]
    fn test_jcamp_multiblock_link_parses() {
        let jdx = "\
##TITLE=Compound Link Container
##JCAMP-DX=5.01
##DATA TYPE=LINK
##BLOCKS=2
##END=
##TITLE=Real Spectrum
##JCAMP-DX=5.01
##DATA TYPE=INFRARED SPECTRUM
##XUNITS=1/CM
##YUNITS=ABSORBANCE
##XYPOINTS=(XY..XY) 500.0 0.05; 1725.0 0.90; 4000.0 0.02
##END=
";
        let spec = parse_jcamp_dx(jdx, "fallback".into(), None)
            .expect("LINK container should parse the inner spectrum block");
        assert_eq!(spec.name, "Real Spectrum");
        assert_eq!(spec.kind, SpectrumKind::Ftir);
        assert_eq!(spec.x.len(), 3);
        assert!((spec.x[1] - 1725.0).abs() < 1e-6);
    }

    #[test]
    fn test_nmr_library_matcher_identifies_pmma_hnmr() {
        let pmma = crate::data::pmma_hnmr_spectrum();
        let matches = crate::matcher::match_nmr_library(
            &pmma,
            ProcessConfig::default(),
            SolventId::Cdcl3,
        );
        assert!(!matches.is_empty(), "NMR matcher returned no candidates");
        assert!(
            matches
                .iter()
                .take(3)
                .any(|m| (m.short.eq_ignore_ascii_case("PMMA") || m.name.contains("PMMA"))
                    && m.hqi_pct > 70.0),
            "Expected PMMA in top 3 NMR matches with HQI > 70%, got: {:?}",
            matches
                .iter()
                .take(3)
                .map(|m| (&m.short, m.hqi_pct))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_native_rust_digitizer_recovers_curve() {
        let w = 400u32;
        let h = 240u32;
        let mut img = image::RgbImage::from_pixel(w, h, image::Rgb([255, 255, 255]));

        // Light-gray horizontal & vertical gridlines inside box (x=40..380, y=20..200)
        for gy in [50u32, 90, 130, 170] {
            for x in 40..=380 {
                img.put_pixel(x, gy, image::Rgb([210, 210, 210]));
            }
        }
        for gx in [100u32, 170, 240, 310] {
            for y in 20..=200 {
                img.put_pixel(gx, y, image::Rgb([210, 210, 210]));
            }
        }

        // Black axis box (x=40..380, y=20..200)
        for x in 40..=380 {
            img.put_pixel(x, 20, image::Rgb([0, 0, 0]));
            img.put_pixel(x, 200, image::Rgb([0, 0, 0]));
        }
        for y in 20..=200 {
            img.put_pixel(40, y, image::Rgb([0, 0, 0]));
            img.put_pixel(380, y, image::Rgb([0, 0, 0]));
        }

        // Small black text/legend rectangle (x=310..360, y=30..45)
        for y in 30..=45 {
            for x in 310..=360 {
                img.put_pixel(x, y, image::Rgb([0, 0, 0]));
            }
        }

        // Solid blue spectrum curve (25, 80, 180), 3px thick:
        // flat 95%T baseline with Lorentzian trough down to 15%T centered at 1725 cm^-1
        // inside the crop box (x=43..377, y=23..197).
        let x0 = 43.0_f64;
        let x1 = 377.0_f64;
        let y0 = 23.0_f64;
        let y1 = 197.0_f64;
        for px in 43u32..=377 {
            let u = (px as f64 - x0) / (x1 - x0);
            let wn = 4000.0 + u * (400.0 - 4000.0);
            let dx = (wn - 1725.0) / 45.0;
            let t_pct = 95.0 - 80.0 / (1.0 + dx * dx);
            let py_f = y0 + (100.0 - t_pct) / 100.0 * (y1 - y0);
            let py = py_f.round() as i32;
            for dy in -1..=1 {
                let yy = py + dy;
                if (21..=199).contains(&yy) {
                    img.put_pixel(px, yy as u32, image::Rgb([25, 80, 180]));
                }
            }
        }

        let tmp = std::env::temp_dir().join(format!(
            "spectra_test_digitize_{}.png",
            std::process::id()
        ));
        img.save(&tmp).expect("save synthetic PNG");

        let spec = crate::digitize::digitize_spectrum_image(
            &tmp,
            SpectrumKind::Ftir,
            YUnit::Transmittance,
            4000.0,
            400.0,
            "auto",
        )
        .expect("digitize synthetic spectrum");
        let _ = std::fs::remove_file(&tmp);

        assert_eq!(spec.kind, SpectrumKind::Ftir);
        assert!(spec.x.len() >= 300, "expected >= 300 points, got {}", spec.x.len());
        assert!(spec.x.windows(2).all(|w| w[0] < w[1]), "x must be strictly ascending");

        let (min_idx, _) = spec
            .y
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| a.total_cmp(b))
            .unwrap();
        let min_x = spec.x[min_idx];
        assert!(
            (min_x - 1725.0).abs() <= 40.0,
            "expected trough near 1725 cm^-1, got {}",
            min_x
        );
    }

    #[test]
    fn test_native_svg_export_ftir_and_uvvis() {
        let tmp_ftir = std::env::temp_dir().join(format!(
            "spectra_test_ftir_{}.svg",
            std::process::id()
        ));
        let tmp_uv = std::env::temp_dir().join(format!(
            "spectra_test_uvvis_{}.svg",
            std::process::id()
        ));

        let x_ir: Vec<f64> = (0..100).map(|i| 400.0 + i as f64 * 36.0).collect();
        let y_ir1: Vec<f64> = x_ir
            .iter()
            .map(|&x| 95.0 - 60.0 / (1.0 + ((x - 1725.0) / 40.0).powi(2)))
            .collect();
        let y_ir2: Vec<f64> = x_ir
            .iter()
            .map(|&x| 92.0 - 50.0 / (1.0 + ((x - 2950.0) / 60.0).powi(2)))
            .collect();

        let ftir_payload = ExportPayload {
            mode: "ftir".into(),
            y_unit: "transmittance".into(),
            style: "acs".into(),
            stacked: true,
            dpi: 600,
            title: "MMA <vs> PMMA & Co.".into(),
            x_bounds: [400.0, 4000.0],
            active_trace: 0,
            traces: vec![
                ExportTrace {
                    name: "PMMA".into(),
                    x: x_ir.clone(),
                    y: y_ir1,
                },
                ExportTrace {
                    name: "MMA".into(),
                    x: x_ir,
                    y: y_ir2,
                },
            ],
            peaks: vec![
                ExportPeak { x: 1725.0, y: 35.0, label: "C=O 1725".into(), pinned: true },
                ExportPeak { x: 2950.0, y: 78.0, label: "C-H 2950".into(), pinned: false },
                ExportPeak { x: 1450.0, y: 82.0, label: "CH2 1450".into(), pinned: false },
                ExportPeak { x: 1150.0, y: 80.0, label: "C-O 1150".into(), pinned: false },
            ],
        };

        let saved_ftir = crate::export::run_export(&ftir_payload, &tmp_ftir)
            .expect("export FTIR SVG");
        assert!(saved_ftir.is_file());
        let ftir_svg = std::fs::read_to_string(&saved_ftir).expect("read FTIR SVG");
        assert!(ftir_svg.starts_with("<svg"));
        assert!(ftir_svg.ends_with("</svg>"));
        assert!(ftir_svg.contains("&amp;"));
        assert!(ftir_svg.contains("&lt;"));
        assert!(ftir_svg.contains("peak-leaders"));
        assert!(ftir_svg.contains("peak-labels"));

        let x_uv: Vec<f64> = (200..=800).step_by(5).map(|v| v as f64).collect();
        let y_uv: Vec<f64> = x_uv
            .iter()
            .map(|&x| 0.05 + 0.9 * (-((x - 450.0) / 40.0).powi(2)).exp())
            .collect();
        let uv_payload = ExportPayload {
            mode: "uvvis".into(),
            y_unit: "absorbance".into(),
            style: "everforest".into(),
            stacked: false,
            dpi: 600,
            title: "UV-Vis Dye".into(),
            x_bounds: [200.0, 800.0],
            active_trace: 0,
            traces: vec![ExportTrace {
                name: "Dye".into(),
                x: x_uv,
                y: y_uv,
            }],
            peaks: vec![ExportPeak {
                x: 450.0,
                y: 0.95,
                label: "λmax 450 nm".into(),
                pinned: true,
            }],
        };

        let saved_uv = crate::export::run_export(&uv_payload, &tmp_uv)
            .expect("export UV-Vis SVG");
        let uv_svg = std::fs::read_to_string(&saved_uv).expect("read UV-Vis SVG");
        assert!(uv_svg.contains("vis-rainbow"));
        assert!(uv_svg.contains("Wavelength"));

        let _ = std::fs::remove_file(&tmp_ftir);
        let _ = std::fs::remove_file(&tmp_uv);
    }

    fn test_spectrum(
        name: &str,
        kind: SpectrumKind,
        raw_y_unit: YUnit,
        x: Vec<f64>,
        y: Vec<f64>,
    ) -> Spectrum {
        Spectrum {
            name: name.to_string(),
            source_path: None,
            kind,
            raw_y_unit,
            x,
            y,
            visible: true,
            color_idx: 0,
            metadata: Vec::new(),
        }
    }

    fn lorentz(x: f64, center: f64, height: f64, fwhm: f64) -> f64 {
        let z = 2.0 * (x - center) / fwhm.max(1e-9);
        height / (1.0 + z * z)
    }

    #[test]
    fn test_deconv_recovers_overlapping_voigt_peaks() {
        let n = 420usize;
        let x_lo = 1660.0;
        let x_hi = 1775.0;
        let xs: Vec<f64> = (0..n)
            .map(|i| x_lo + (x_hi - x_lo) * i as f64 / (n - 1) as f64)
            .collect();
        let specs = [(1703.0, 0.62, 9.0), (1717.0, 0.30, 11.0), (1732.0, 0.47, 8.5)];
        let ys: Vec<f64> = xs
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                let mut v = 0.035 + 0.00006 * (x - 1700.0);
                for (c, h, w) in specs {
                    v += lorentz(x, c, h, w);
                }
                // deterministic pseudo-noise, no RNG dependency
                let noise = 0.0006 * (((i * 1103515245 + 12345) % 1000) as f64 / 500.0 - 1.0);
                v + noise
            })
            .collect();

        let cfg = DeconvConfig {
            shape: PeakShape::PseudoVoigt,
            max_peaks: 3,
            fit_baseline: true,
            manual_centers: Vec::new(),
            max_iter: 150,
        };
        let res = deconvolve_window(
            SpectrumKind::Ftir,
            &xs,
            &ys,
            YUnit::Absorbance,
            x_lo,
            x_hi,
            &cfg,
            SolventId::Cdcl3,
        )
        .expect("deconvolution failed");

        assert!(res.r_squared > 0.995, "R² too low: {:.6}", res.r_squared);
        assert_eq!(res.peaks.len(), 3, "expected 3 peaks, got {}", res.peaks.len());
        for target in [1703.0, 1717.0, 1732.0] {
            assert!(
                res.peaks.iter().any(|p| (p.center - target).abs() <= 1.5),
                "peak at {} not recovered; fitted centers: {:?}",
                target,
                res.peaks.iter().map(|p| p.center).collect::<Vec<_>>()
            );
        }
        let pct: f64 = res.peaks.iter().map(|p| p.area_pct).sum();
        assert!((pct - 100.0).abs() < 0.1, "area percentages sum to {:.3}", pct);
    }

    #[test]
    fn test_deconv_nmr_tacticity_lorentzian() {
        let n = 640usize;
        let x_lo = 0.65;
        let x_hi = 1.40;
        let xs: Vec<f64> = (0..n)
            .map(|i| x_lo + (x_hi - x_lo) * i as f64 / (n - 1) as f64)
            .collect();
        let specs = [(1.21, 0.55, 0.030), (1.02, 0.55, 0.028), (0.85, 0.55, 0.032)];
        let ys: Vec<f64> = xs
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                let mut v = 0.002;
                for (c, h, w) in specs {
                    v += lorentz(x, c, h, w);
                }
                v + 0.0004 * (((i * 2654435761 + 7) % 997) as f64 / 498.0 - 1.0)
            })
            .collect();

        let cfg = DeconvConfig {
            shape: PeakShape::Lorentzian,
            max_peaks: 3,
            fit_baseline: true,
            manual_centers: Vec::new(),
            max_iter: 150,
        };
        let res = deconvolve_window(
            SpectrumKind::Hnmr,
            &xs,
            &ys,
            YUnit::Intensity,
            x_lo,
            x_hi,
            &cfg,
            SolventId::Cdcl3,
        )
        .expect("NMR deconvolution failed");

        assert!(res.r_squared > 0.995, "NMR R² too low: {:.6}", res.r_squared);
        assert_eq!(res.peaks.len(), 3, "expected 3 triads, got {}", res.peaks.len());
        for target in [1.21, 1.02, 0.85] {
            assert!(
                res.peaks.iter().any(|p| (p.center - target).abs() <= 0.01),
                "triad at {} ppm not recovered; got {:?}",
                target,
                res.peaks.iter().map(|p| p.center).collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn test_arpls_noisy_baseline_recovery() {
        let n = 900usize;
        let xs: Vec<f64> = (0..n)
            .map(|i| 400.0 + (4000.0 - 400.0) * i as f64 / (n - 1) as f64)
            .collect();
        let truth_baseline = |x: f64| 0.06 + 0.20 * ((x - 2200.0) / 1800.0).powi(2);
        let ys: Vec<f64> = xs
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                let peak1 = lorentz(x, 1728.0, 0.65, 26.0);
                let peak2 = lorentz(x, 2940.0, 0.42, 45.0);
                let noise = 0.004 * (((i * 1103515245 + 12345) % 10000) as f64 / 5000.0 - 1.0);
                truth_baseline(x) + peak1 + peak2 + noise
            })
            .collect();

        let corrected = arpls_baseline(&ys, 1e6, 1e-3, 40);
        // Off-peak residual must be essentially flat and near zero.
        for (i, &x) in xs.iter().enumerate() {
            let off_peak = (x < 1600.0) || (x > 3200.0);
            if off_peak {
                assert!(
                    corrected[i].abs() < 0.04,
                    "off-peak residual {:.4} at {:.0} cm⁻¹ exceeds 0.04",
                    corrected[i],
                    x
                );
            }
        }
        let max_near = |center: f64| {
            let i0 = xs.partition_point(|&v| v < center - 20.0);
            let i1 = xs.partition_point(|&v| v <= center + 20.0);
            corrected[i0..i1.max(i0 + 1)]
                .iter()
                .copied()
                .fold(0.0_f64, f64::max)
        };
        assert!(
            max_near(1728.0) > 0.85 * 0.65,
            "carbonyl peak lost: {:.4}",
            max_near(1728.0)
        );
        assert!(
            max_near(2940.0) > 0.85 * 0.42,
            "C-H peak lost: {:.4}",
            max_near(2940.0)
        );
    }

    #[test]
    fn test_advanced_atr_kk_recovers_dispersion_shift() {
        #[derive(Clone, Copy)]
        struct C {
            re: f64,
            im: f64,
        }
        impl C {
            fn new(re: f64, im: f64) -> Self {
                Self { re, im }
            }
            fn add(self, o: Self) -> Self {
                Self::new(self.re + o.re, self.im + o.im)
            }
            fn sub(self, o: Self) -> Self {
                Self::new(self.re - o.re, self.im - o.im)
            }
            fn mul(self, o: Self) -> Self {
                Self::new(
                    self.re * o.re - self.im * o.im,
                    self.re * o.im + self.im * o.re,
                )
            }
            fn div(self, o: Self) -> Self {
                let den = (o.re * o.re + o.im * o.im).max(1e-30);
                Self::new(
                    (self.re * o.re + self.im * o.im) / den,
                    (self.im * o.re - self.re * o.im) / den,
                )
            }
            fn norm_sq(self) -> f64 {
                self.re * self.re + self.im * self.im
            }
            fn sqrt(self) -> Self {
                let r = self.re.hypot(self.im);
                let re = ((r + self.re) * 0.5).max(0.0).sqrt();
                let im = ((r - self.re) * 0.5).max(0.0).sqrt() * if self.im < 0.0 { -1.0 } else { 1.0 };
                Self::new(re, im)
            }
        }
        fn fresnel_abs(n1: f64, n2r: f64, k: f64, s: f64, c: f64) -> f64 {
            let n2 = C::new(n2r, k.max(0.0));
            let ratio = C::new(n1 * s, 0.0).div(n2);
            let one_minus = C::new(1.0, 0.0).sub(ratio.mul(ratio));
            let mut ct = one_minus.sqrt();
            if n2.mul(ct).im < 0.0 {
                ct = C::new(-ct.re, -ct.im);
            }
            let n1c = C::new(n1 * c, 0.0);
            let n2ct = n2.mul(ct);
            let rs = n1c.sub(n2ct).div(n1c.add(n2ct));
            let n2c = n2.mul(C::new(c, 0.0));
            let n1ct = ct.mul(C::new(n1, 0.0));
            let rp = n2c.sub(n1ct).div(n2c.add(n1ct));
            let r = (0.5 * (rs.norm_sq() + rp.norm_sq())).clamp(1e-6, 1.0);
            -r.log10()
        }

        let n = 1800usize;
        let wn: Vec<f64> = (0..n)
            .map(|i| 400.0 + (4000.0 - 400.0) * i as f64 / (n - 1) as f64)
            .collect();
        let th = 45.0_f64.to_radians();
        let (s, c) = (th.sin(), th.cos());
        let center = 1728.0_f64;
        let k0 = 0.28_f64;
        let fwhm = 24.0_f64;

        // Exact Lorentz oscillator anomalous dispersion: n(ν) = n_inf - k0 * z / (1 + z^2)
        let raw: Vec<f64> = wn
            .iter()
            .map(|&x| {
                let z = 2.0 * (x - center) / fwhm;
                let denom = 1.0 + z * z;
                let k = k0 / denom;
                let n2_real = 1.50 - k0 * z / denom;
                fresnel_abs(2.4, n2_real, k, s, c)
            })
            .collect();

        let (raw_idx, _) = raw
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .expect("raw ATR empty");
        let raw_peak = wn[raw_idx];
        assert!(
            raw_peak < center - 2.0,
            "Expected synthetic Fresnel ATR peak to red-shift below {:.1} cm⁻¹, got {:.1}",
            center - 2.0,
            raw_peak
        );

        let corrected = advanced_atr_correction(&wn, &raw, 2.4, 1.50, 45.0);
        let (idx, _) = corrected
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.total_cmp(b))
            .expect("ATR correction returned empty");
        let recovered = wn[idx];
        assert!(
            (recovered - center).abs() <= 4.0,
            "ATR-corrected peak at {:.1} cm⁻¹ (raw was {:.1}), expected {:.1} ± 4",
            recovered,
            raw_peak,
            center
        );
    }

    #[test]
    fn test_kubelka_munk_and_snv_normalization() {
        let km = kubelka_munk_transform(&[50.0]);
        assert!(
            (km[0] - 0.25).abs() < 1e-12,
            "KM F(50%R) should be 0.25, got {:.12}",
            km[0]
        );

        let n = 300usize;
        let xs: Vec<f64> = (0..n)
            .map(|i| 190.0 + (900.0 - 190.0) * i as f64 / (n - 1) as f64)
            .collect();
        let ys: Vec<f64> = xs
            .iter()
            .map(|&x| 0.04 + 0.8 * (-((x - 420.0) / 55.0).powi(2)).exp())
            .collect();
        let spec = test_spectrum(
            "SNV test",
            SpectrumKind::UvVis,
            YUnit::Absorbance,
            xs,
            ys,
        );
        let (_, yn, yu) = process_trace(
            &spec,
            ProcessConfig {
                norm: NormMode::VectorSnv,
                ..ProcessConfig::default()
            },
        );
        assert_eq!(yu, YUnit::Absorbance);
        let mean = yn.iter().sum::<f64>() / yn.len() as f64;
        let var = yn.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / yn.len() as f64;
        assert!(mean.abs() < 1e-9, "SNV mean = {:.3e}", mean);
        assert!((var.sqrt() - 1.0).abs() < 1e-9, "SNV std = {:.6}", var.sqrt());
    }

    #[test]
    fn test_spectral_subtraction_optimal_k_derivative() {
        let n = 600usize;
        let xs: Vec<f64> = (0..n)
            .map(|i| 170.0 + (700.0 - 170.0) * i as f64 / (n - 1) as f64)
            .collect();
        let b_y: Vec<f64> = xs
            .iter()
            .map(|&x| 0.05 + 0.7 * (-((x - 300.0) / 30.0).powi(2)).exp())
            .collect();
        let additive: Vec<f64> = xs
            .iter()
            .map(|&x| 0.22 * (-((x - 470.0) / 20.0).powi(2)).exp())
            .collect();
        let a_y: Vec<f64> = xs
            .iter()
            .enumerate()
            .map(|(i, &x)| {
                0.72 * b_y[i] + additive[i] + 0.00025 * (x - 435.0)
            })
            .collect();
        let a = test_spectrum(
            "Sample A",
            SpectrumKind::UvVis,
            YUnit::Absorbance,
            xs.clone(),
            a_y,
        );
        let b = test_spectrum(
            "Component B",
            SpectrumKind::UvVis,
            YUnit::Absorbance,
            xs.clone(),
            b_y,
        );

        let k = optimal_subtraction_factor(&a, &b, 170.0, 700.0).expect("no optimal k");
        assert!(
            (k - 0.72).abs() <= 0.02,
            "optimal k* should be ~0.72, got {:.4}",
            k
        );

        let diff = compute_spectral_math(&a, &b, SpectralOp::Subtract, k).expect("spectral math");
        assert_eq!(diff.kind, SpectrumKind::UvVis);
        assert_eq!(diff.raw_y_unit, YUnit::Absorbance);
        let i_peak = diff
            .x
            .iter()
            .position(|&x| x >= 470.0)
            .expect("isolation peak missing");
        assert!(
            diff.y[i_peak] > 0.15,
            "subtraction did not isolate additive peak (y = {:.4})",
            diff.y[i_peak]
        );
        let i_flat = diff
            .x
            .iter()
            .position(|&x| x >= 250.0)
            .expect("flat region missing");
        assert!(
            diff.y[i_flat].abs() < 0.05,
            "baseline not cancelled off-peak (y = {:.4})",
            diff.y[i_flat]
        );
    }

    fn find_polymer_by_short(short: &str) -> &'static crate::data::PolymerRef {
        crate::data::POLYMERS_DB
            .iter()
            .find(|p| p.short.eq_ignore_ascii_case(short))
            .unwrap_or_else(|| panic!("Polymer '{}' not found in POLYMERS_DB", short))
    }

    fn blend_ftir_abs(parts: &[(&crate::data::PolymerRef, f64)]) -> Spectrum {
        let mut specs: Vec<(Vec<f64>, Vec<f64>, f64)> = Vec::new();
        for &(poly, w) in parts {
            let s = crate::data::synthesize_polymer_reference(poly, SpectrumKind::Ftir);
            let abs_raw: Vec<f64> = if s.raw_y_unit == YUnit::Transmittance {
                s.y.iter().copied().map(transmittance_to_absorbance).collect()
            } else {
                s.y.clone()
            };
            let min_v = abs_raw.iter().copied().fold(f64::INFINITY, f64::min);
            let max_v = abs_raw.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let span = (max_v - min_v).max(1e-12);
            let norm: Vec<f64> = abs_raw.iter().map(|&v| (v - min_v) / span).collect();
            specs.push((s.x, norm, w));
        }
        let xs = specs[0].0.clone();
        let mut ys = vec![0.0; xs.len()];
        for (sx, sy, w) in &specs {
            for (i, &x) in xs.iter().enumerate() {
                if let Some(v) = crate::processing::interpolate_y(sx, sy, x) {
                    ys[i] += w * v;
                }
            }
        }
        Spectrum {
            name: "Synthetic FTIR Blend".to_string(),
            source_path: None,
            kind: SpectrumKind::Ftir,
            raw_y_unit: YUnit::Absorbance,
            x: xs,
            y: ys,
            visible: true,
            color_idx: 0,
            metadata: Vec::new(),
        }
    }

    #[test]
    fn test_nnls_mixture_identifies_binary_and_ternary_blends() {
        let pmma = find_polymer_by_short("PMMA");
        let ps = find_polymer_by_short("PS");
        let pet = find_polymer_by_short("PET");
        let pa66 = crate::data::POLYMERS_DB
            .iter()
            .find(|p| {
                let s = p.short.to_ascii_uppercase();
                s == "PA66" || s == "PA6,6" || s.contains("NYLON")
            })
            .expect("PA66 / Nylon polymer not found in POLYMERS_DB");

        // Binary blend: 60% PMMA + 40% PS
        let bin_blend = blend_ftir_abs(&[(pmma, 0.60), (ps, 0.40)]);
        let bin_res = crate::mixture::analyze_mixture(
            &bin_blend,
            ProcessConfig::default(),
            SolventId::Cdcl3,
            3,
        );
        assert!(!bin_res.is_empty(), "binary mixture returned no matches");
        let top_bin = &bin_res[0];
        assert!(
            top_bin.composite_hqi_pct > 92.0,
            "binary composite HQI too low: {:.2}%",
            top_bin.composite_hqi_pct
        );
        let find_w = |m: &crate::mixture::MixtureMatch, poly: &crate::data::PolymerRef| -> Option<f64> {
            let target = poly.name.to_ascii_lowercase();
            m.components
                .iter()
                .find(|c| {
                    c.short.eq_ignore_ascii_case(poly.short)
                        || c.name.to_ascii_lowercase().contains(&target)
                        || target.contains(&c.name.to_ascii_lowercase())
                })
                .map(|c| c.weight_pct)
        };
        let w_pmma = find_w(top_bin, pmma).unwrap_or_else(|| {
            panic!("PMMA missing from top binary match: {:?}", top_bin.components)
        });
        let w_ps = find_w(top_bin, ps).unwrap_or_else(|| {
            panic!("PS missing from top binary match: {:?}", top_bin.components)
        });
        assert!(
            (w_pmma - 60.0).abs() <= 12.0 && (w_ps - 40.0).abs() <= 12.0,
            "binary weights out of ±12% tolerance: PMMA={:.1}%, PS={:.1}%",
            w_pmma,
            w_ps
        );

        // Ternary blend: 50% PET + 30% PA66 + 20% PS
        let tern_blend = blend_ftir_abs(&[(pet, 0.50), (pa66, 0.30), (ps, 0.20)]);
        let tern_res = crate::mixture::analyze_mixture(
            &tern_blend,
            ProcessConfig::default(),
            SolventId::Cdcl3,
            3,
        );
        assert!(!tern_res.is_empty(), "ternary mixture returned no matches");
        let top_tern = &tern_res[0];
        assert!(
            top_tern.composite_hqi_pct > 92.0,
            "ternary composite HQI too low: {:.2}%",
            top_tern.composite_hqi_pct
        );
        let w_pet = find_w(top_tern, pet).unwrap_or_else(|| {
            panic!("PET missing from top ternary match: {:?}", top_tern.components)
        });
        let w_pa = find_w(top_tern, pa66).unwrap_or_else(|| {
            panic!("PA66 missing from top ternary match: {:?}", top_tern.components)
        });
        let w_ps3 = find_w(top_tern, ps).unwrap_or_else(|| {
            panic!("PS missing from top ternary match: {:?}", top_tern.components)
        });
        assert!(
            (w_pet - 50.0).abs() <= 12.0
                && (w_pa - 30.0).abs() <= 12.0
                && (w_ps3 - 20.0).abs() <= 12.0,
            "ternary weights out of ±12% tolerance: PET={:.1}%, PA66={:.1}%, PS={:.1}%",
            w_pet,
            w_pa,
            w_ps3
        );
    }

    #[test]
    fn test_nnls_does_not_overfit_pure_spectrum() {
        let pmma = pmma_polymer_spectrum();
        let res = crate::mixture::analyze_mixture(
            &pmma,
            ProcessConfig::default(),
            SolventId::Cdcl3,
            3,
        );
        assert!(!res.is_empty(), "pure PMMA returned no matches");
        let top = &res[0];
        assert_eq!(
            top.components.len(),
            1,
            "pure PMMA overfitted into {} components: {:?}",
            top.components.len(),
            top.components
        );
        assert!(
            top.components[0].short.eq_ignore_ascii_case("PMMA")
                || top.components[0].name.to_ascii_uppercase().contains("PMMA")
                || top.components[0].name.to_ascii_lowercase().contains("methyl methacrylate"),
            "unexpected top single match: {}",
            top.components[0].name
        );
        assert!(
            top.composite_hqi_pct > 88.0,
            "pure PMMA HQI too low: {:.2}%",
            top.composite_hqi_pct
        );
    }

    #[test]
    fn test_nnls_gram_enforces_non_negativity() {
        // SPD Gram matrix whose unconstrained solution is [1.0, -0.6, 0.4]
        let g = vec![
            vec![2.0, 1.0, 0.5],
            vec![1.0, 2.0, 0.5],
            vec![0.5, 0.5, 2.0],
        ];
        let rhs = vec![1.6, 0.0, 1.0];
        let c = crate::mixture::nnls_gram(&g, &rhs);
        assert_eq!(c.len(), 3);
        for (j, &cj) in c.iter().enumerate() {
            assert!(cj >= 0.0, "c[{}] = {} is negative", j, cj);
        }
        assert!(c[1].abs() < 1e-9, "expected c[1] == 0, got {}", c[1]);
        assert!((c[0] - 0.72).abs() < 1e-5, "c[0] = {}", c[0]);
        assert!((c[2] - 0.32).abs() < 1e-5, "c[2] = {}", c[2]);

        // KKT optimality: grad = rhs - G c
        for j in 0..3 {
            let gc: f64 = (0..3).map(|m| g[j][m] * c[m]).sum();
            let grad = rhs[j] - gc;
            if c[j] <= 1e-12 {
                assert!(grad <= 1e-6, "KKT inactive violation at {}: grad={}", j, grad);
            } else {
                assert!(grad.abs() <= 1e-6, "KKT active violation at {}: grad={}", j, grad);
            }
        }
    }

    #[test]
    fn test_integration_regions_relative_ratio() {
        let spec = crate::data::pmma_hnmr_spectrum();
        let (a_och3, h_och3) = crate::processing::integrate_spectrum_window(
            &spec.x,
            &spec.y,
            spec.raw_y_unit,
            3.35,
            3.85,
        )
        .expect("OCH3 integration failed");
        let (a_ch3, h_ch3) = crate::processing::integrate_spectrum_window(
            &spec.x,
            &spec.y,
            spec.raw_y_unit,
            0.65,
            1.35,
        )
        .expect("alpha-CH3 integration failed");

        assert!(a_och3 > 0.0 && a_och3.is_finite() && h_och3 > 0.0);
        assert!(a_ch3 > 0.0 && a_ch3.is_finite() && h_ch3 > 0.0);
        let ratio = a_och3 / a_ch3;
        assert!(ratio > 0.2 && ratio < 5.0 && ratio.is_finite(), "unexpected ratio {:.4}", ratio);
    }

    #[test]
    fn test_x_axis_calibration_shifts_peaks() {
        let orig = polystyrene_ftir_spectrum();
        let peak_in = |s: &Spectrum, lo: f64, hi: f64| -> f64 {
            let mut best_x = lo;
            let mut best_a = f64::NEG_INFINITY;
            for (&x, &y) in s.x.iter().zip(&s.y) {
                if x >= lo && x <= hi {
                    let a = if s.raw_y_unit == YUnit::Transmittance {
                        transmittance_to_absorbance(y)
                    } else {
                        y
                    };
                    if a > best_a {
                        best_a = a;
                        best_x = x;
                    }
                }
            }
            best_x
        };

        let orig_peak = peak_in(&orig, 1580.0, 1620.0);
        let mut shifted = orig.clone();
        for x in &mut shifted.x {
            *x -= 4.5;
        }
        crate::processing::apply_x_calibration(&mut shifted, 4.5);

        let restored_peak = peak_in(&shifted, 1580.0, 1620.0);
        assert!(
            (restored_peak - orig_peak).abs() <= 0.1,
            "calibrated peak {:.3} != original {:.3}",
            restored_peak,
            orig_peak
        );
        assert!(
            shifted.metadata.iter().any(|(k, _)| k == "X Calibration"),
            "X Calibration metadata not recorded"
        );
    }

    #[test]
    fn test_hnmr_j_coupling_doublet_triplet_detection() {
        let mk_peak = |x: f64| crate::processing::Peak {
            x,
            y: 1.0,
            prominence: 1.0,
            fwhm: Some(0.01),
            snr: 50.0,
            label: String::new(),
        };

        // Doublet: J = 7.2 Hz at 400 MHz -> spacing = 0.018 ppm
        let doublet = vec![mk_peak(4.100), mk_peak(4.118)];
        let d_info = crate::processing::analyze_hnmr_multiplet(&doublet, 4.109, 400.0)
            .expect("doublet not detected");
        assert_eq!(d_info.multiplicity, "d");
        assert_eq!(d_info.line_count, 2);
        assert!(
            (d_info.j_hz[0] - 7.2).abs() <= 0.1,
            "expected doublet J=7.2 Hz, got {:.3}",
            d_info.j_hz[0]
        );

        // Triplet: J = 7.0 Hz at 400 MHz -> spacing = 0.0175 ppm
        let triplet = vec![mk_peak(1.2000), mk_peak(1.2175), mk_peak(1.2350)];
        let t_info = crate::processing::analyze_hnmr_multiplet(&triplet, 1.218, 400.0)
            .expect("triplet not detected");
        assert_eq!(t_info.multiplicity, "t");
        assert_eq!(t_info.line_count, 3);
        assert!(
            (t_info.j_hz[0] - 7.0).abs() <= 0.1,
            "expected triplet J=7.0 Hz, got {:.3}",
            t_info.j_hz[0]
        );
    }
}

