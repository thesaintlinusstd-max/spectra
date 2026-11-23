mod app;
mod clipboard;
mod data;
mod export;
mod io;
mod kitty_plot;
mod matcher;
mod processing;
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
use crate::export::{find_python_script, run_python_export, ExportPayload, ExportPeak, ExportStyle, ExportTrace};
use crate::io::{
    load_spectrum_file, resolve_library_query, search_library_for_kind, write_sample_files,
    LibraryEntry,
};
use crate::processing::{detect_peaks, process_trace, ProcessConfig, SmoothLevel};

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
        let _ = execute!(stdout, DisableMouseCapture, LeaveAlternateScreen);
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
           spectra [FILES...] [--mode ftir|hnmr|cnmr]\n  \
           spectra [ftir|nmr|cnmr] [QUERY...]     (e.g. spectra ftir butane)\n  \
           spectra init-samples [--dir samples]\n  \
           spectra fetch <nist|nmr|nmr-search|import|list|pack> [ARGS...]\n  \
           spectra export <FILES...> -o <OUT.png|svg|pdf> [--style acs|rsc|everforest]\n                 \
                   [--stacked] [--absorbance] [--baseline] [--smooth 0|5|9|15]\n                 \
                   [--title \"...\"] [--dpi 600]\n",
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
        ftir_absorbance: absorbance,
        baseline,
        smooth: smooth_lvl,
        normalize: false,
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

    let saved = run_python_export(&payload, &output)?;
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
            files, output, style, stacked, absorbance, baseline, smooth, title, dpi,
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
        pmma_polymer_spectrum, polystyrene_ftir_spectrum, SolventId, SpectrumKind, Spectrum,
        YUnit,
    };
    use crate::io::{encode_jcamp_difdup, parse_delimited, parse_jcamp_dx};
    use crate::processing::{
        absorbance_to_transmittance, auto_ftir_conversion, beer_lambert_conc,
        compute_tauc_band_gap, nmr_mn_dpn, pu_hydrogen_bond_index, pvdf_beta_fraction,
        rubberband_baseline, savitzky_golay, tacticity_stats, transmittance_to_absorbance,
        TaucMode,
    };

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

        let pvdf = match_ftir(840.0);
        assert!(
            pvdf.iter().any(|b| b.group.contains("PVDF β-phase")),
            "expected PVDF β-phase assignment at 840 cm⁻¹, got {:?}",
            pvdf.iter().map(|b| b.group).collect::<Vec<_>>()
        );
        assert!(
            pvdf[0].group.contains("PVDF β-phase"),
            "PVDF β-phase should rank first at 840 cm⁻¹, got: {}",
            pvdf[0].group
        );
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
            is_image: false,
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
    fn test_version_is_0_5() {
        assert_eq!(env!("CARGO_PKG_VERSION"), "0.5.1");
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
}

