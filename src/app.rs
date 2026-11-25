use std::fs;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver};
use std::thread;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

use crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseButton, MouseEvent, MouseEventKind};
use ratatui::layout::Rect;

use crate::clipboard::{copy_png_file, copy_text};
use crate::data::{
    best_peak_label_in_context, find_polymer_by_short,
    generate_demo_preset, match_cnmr, match_ftir, match_hnmr, match_solvent_impurities,
    match_uvvis, synthesize_polymer_reference, uvvis_color_info, DemoPreset, SolventId, Spectrum,
    SpectrumKind, YUnit, CNMR_DB, FTIR_DB, HNMR_DB, POLYMERS_DB, SOLVENT_IMPURITIES, UVVIS_DB,
};
use crate::deconv::{
    deconvolve_window, format_deconv_csv, format_deconv_markdown, DeconvConfig, DeconvResult,
    PeakShape,
};
use crate::export::{
    open_exported_file, run_export, ExportFormat, ExportPayload, ExportPeak, ExportPeakMode,
    ExportStyle, ExportTrace,
};
use crate::io::{expand_tilde, load_spectrum_file, parse_pasted_paths, search_library};
use crate::processing::{
    apply_x_calibration, auto_ftir_conversion, average_spectra, beer_lambert_conc,
    compute_spectral_math, compute_tauc_band_gap_with_optics, detect_peaks,
    ftir_conversion_percent, integrate_spectrum_window, interpolate_y, nmr_mn_dpn,
    optimal_subtraction_factor, process_trace, pu_hydrogen_bond_index, pvdf_beta_fraction,
    tacticity_stats, BaselineMode, ConversionResult, HydrogenBondIndex, IntegrationRegion,
    OpticsMode, Peak, ProcessConfig, SmoothLevel, SpectralOp, TacticityResult, TaucMode,
    TaucResult, YAxisPreference,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActiveTab {
    Ftir,
    Hnmr,
    Cnmr,
    UvVis,
    Solvents,
    PolymerLab,
}

impl ActiveTab {
    pub fn as_kind(self) -> Option<SpectrumKind> {
        match self {
            Self::Ftir => Some(SpectrumKind::Ftir),
            Self::Hnmr => Some(SpectrumKind::Hnmr),
            Self::Cnmr => Some(SpectrumKind::Cnmr),
            Self::UvVis => Some(SpectrumKind::UvVis),
            _ => None,
        }
    }

    pub fn from_kind(k: SpectrumKind) -> Self {
        match k {
            SpectrumKind::Ftir => Self::Ftir,
            SpectrumKind::Hnmr => Self::Hnmr,
            SpectrumKind::Cnmr => Self::Cnmr,
            SpectrumKind::UvVis => Self::UvVis,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PolymerSubTab {
    Conversion,
    NmrMn,
    Tacticity,
    PhaseIndex,
    UvVisLab,
}

#[derive(Debug, Clone)]
pub struct ProcessedTrace {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub y_unit: YUnit,
}

#[derive(Debug, Clone)]
pub struct SearchItem {
    pub category: &'static str,
    pub target_tab: ActiveTab,
    pub jump_x: Option<f64>,
    pub title: String,
    pub detail: String,
    /// Set only for Polymer entries, so Enter can load a real reference spectrum.
    pub polymer_short: Option<&'static str>,
    /// Set for reference files in `library/` (or samples), so Enter loads them.
    pub file_path: Option<PathBuf>,
}

#[derive(Debug, Clone)]
pub struct SearchState {
    pub query: String,
    pub results: Vec<SearchItem>,
    pub selected: usize,
}

#[derive(Debug, Clone)]
pub struct FileEntry {
    pub path: PathBuf,
    pub display: String,
    pub is_dir: bool,
}

#[derive(Debug, Clone)]
pub struct FilePickerState {
    pub dir: PathBuf,
    /// Filtered view (`filter` applied); `selected` indexes this vector.
    pub entries: Vec<FileEntry>,
    /// Full directory listing, kept so typing a query is instant.
    pub all_entries: Vec<FileEntry>,
    pub selected: usize,
    pub append_mode: bool,
    /// Inline path entry (`/` or `Ctrl+L`), so files can be typed or pasted.
    pub path_input: Option<String>,
    /// Live search typed directly in the picker (substring / token match).
    pub filter: String,
}

/// Score one file name against a live-search query.
/// Returns `None` when any query token is missing; lower ranks are better:
/// 0 = filename prefix, 1 = exact word, 2 = word prefix, 3 = plain substring.
fn filter_rank(name: &str, tokens: &[&str]) -> Option<u8> {
    if !tokens.iter().all(|t| name.contains(t)) {
        return None;
    }
    let words: Vec<&str> = name
        .split(|c: char| !c.is_alphanumeric())
        .filter(|w| !w.is_empty())
        .collect();
    let mut best = 4u8;
    for token in tokens {
        let rank = if name.starts_with(token) {
            0
        } else if words.contains(token) {
            1
        } else if words.iter().any(|w| w.starts_with(token)) {
            2
        } else {
            3
        };
        best = best.min(rank);
    }
    Some(best)
}

/// Apply the live search to a directory listing: directories stay grouped
/// first, then prefix matches, word-start matches, and plain substrings.
pub(crate) fn filter_file_entries(entries: &[FileEntry], filter: &str) -> Vec<FileEntry> {
    let query = filter.trim().to_lowercase();
    if query.is_empty() {
        return entries.to_vec();
    }
    let tokens: Vec<&str> = query.split_whitespace().collect();
    let mut scored: Vec<(u8, usize, FileEntry)> = Vec::new();
    for (order, entry) in entries.iter().enumerate() {
        // The synthetic "../" row only makes sense in the unfiltered listing.
        if entry.display.starts_with("..") {
            continue;
        }
        let name = entry
            .path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_lowercase();
        if name.is_empty() {
            continue;
        }
        if let Some(rank) = filter_rank(&name, &tokens) {
            scored.push((rank, order, entry.clone()));
        }
    }
    scored.sort_by(|a, b| {
        let group = |e: &FileEntry| if e.is_dir { 0u8 } else { 1 };
        group(&a.2)
            .cmp(&group(&b.2))
            .then_with(|| a.0.cmp(&b.0))
            .then_with(|| a.1.cmp(&b.1))
    });
    scored.into_iter().map(|(_, _, entry)| entry).collect()
}

impl FilePickerState {
    pub fn new(dir: PathBuf, append_mode: bool) -> Self {
        let all_entries = scan_directory(&dir);
        let mut state = Self {
            dir,
            entries: Vec::new(),
            all_entries,
            selected: 0,
            append_mode,
            path_input: None,
            filter: String::new(),
        };
        state.refresh();
        crate::io::remember_picker_dir(&state.dir);
        state
    }

    /// Recompute the filtered view and keep the selection in range.
    fn refresh(&mut self) {
        self.entries = filter_file_entries(&self.all_entries, &self.filter);
        if self.selected >= self.entries.len() {
            self.selected = self.entries.len().saturating_sub(1);
        }
    }

    fn change_filter(&mut self) {
        self.selected = 0;
        self.refresh();
    }

    /// Enter a directory, reset the search and remember it for next launch.
    fn set_dir(&mut self, dir: PathBuf) {
        self.dir = dir;
        self.all_entries = scan_directory(&self.dir);
        self.filter.clear();
        self.path_input = None;
        self.selected = 0;
        self.refresh();
        crate::io::remember_picker_dir(&self.dir);
    }
}

#[derive(Debug, Clone)]
pub struct ImportModalState {
    pub path: PathBuf,
    pub is_image: bool,
    pub append_mode: bool,
    pub field_idx: usize,
    pub kind_choice: usize, // 0=FTIR %T, 1=FTIR Abs, 2=¹H NMR, 3=¹³C NMR, 4=Auto (data only), 5=UV-Vis Abs, 6=UV-Vis %T
    pub x_left: String,
    pub x_right: String,
    pub bg_idx: usize, // 0=Auto, 1=Light, 2=Dark
    pub error: Option<String>,
}

#[derive(Debug, Clone)]
pub struct ExportModalState {
    pub field_idx: usize,
    pub format: ExportFormat,
    pub style: ExportStyle,
    pub stacked: bool,
    pub use_zoom_range: bool,
    pub peak_mode: ExportPeakMode,
    pub title: String,
    pub output_path: String,
    pub last_exported: Option<PathBuf>,
    pub feedback: Option<String>,
}

#[derive(Debug, Clone)]
pub struct LibraryMatchState {
    pub results: Vec<crate::matcher::LibraryMatch>,
    pub selected: usize,
    pub mixture_mode: bool,
    pub mixture_results: Vec<crate::mixture::MixtureMatch>,
    pub mixture_selected: usize,
}

#[derive(Debug, Clone)]
pub struct IntegrationModalState {
    pub selected: usize,
    pub feedback: Option<String>,
}

#[derive(Debug, Clone)]
pub struct CalibrationModalState {
    pub observed_x: f64,
    pub target_buf: String,
    pub feedback: Option<String>,
}

/// Shown when a CLI query matches several reference spectra (e.g. `butane`
/// matches n-butane, isobutane, 1-chlorobutane, …). The user picks one.
#[derive(Debug, Clone)]
pub struct CandidatePickerState {
    pub entries: Vec<crate::io::LibraryEntry>,
    pub selected: usize,
}

/// State for the Reference Library Downloader modal (Shift+D / D).
#[derive(Debug, Clone)]
pub struct LibrarySetupState {
    pub selected: usize,          // 0..4 index of current row
    pub checks: [bool; 5],        // checkbox state per pack
    pub feedback: Option<String>, // status message after download
}

impl Default for LibrarySetupState {
    fn default() -> Self {
        // Default: NIST FTIR (~2 MB) + UV-Vis (~2 MB) checked for quick start.
        // Large packs (Zenodo 253 MB, nmrshiftdb2 159 MB, Mega FTIR 158 MB) unchecked.
        Self {
            selected: 0,
            checks: [true, true, false, false, false],
            feedback: None,
        }
    }
}

/// Long-running work that must not block the UI thread.
#[derive(Debug, Clone)]
pub struct BusyState {
    pub label: String,
    pub started: Instant,
}

#[derive(Debug, Clone)]
pub struct DeconvOverlay {
    pub trace_idx: usize,
    pub result: DeconvResult,
    pub show_components: bool,
    pub show_residual: bool,
}

#[derive(Debug, Clone)]
pub struct DeconvModalState {
    pub field_idx: usize,
    pub shape: PeakShape,
    pub max_peaks: usize,
    pub fit_baseline: bool,
    pub x_lo: String,
    pub x_hi: String,
    pub result: Option<DeconvResult>,
    pub feedback: Option<String>,
}

#[derive(Debug, Clone)]
pub struct SpectralMathState {
    pub field_idx: usize,
    pub op: SpectralOp,
    pub trace_a: usize,
    pub trace_b: usize,
    pub factor_k: f64,
    pub feedback: Option<String>,
}

pub enum TaskOutcome {
    Loaded {
        specs: Vec<Spectrum>,
        errors: Vec<String>,
        append: bool,
    },
    LibraryMatch(Vec<crate::matcher::LibraryMatch>),
    MixtureMatch(Vec<crate::mixture::MixtureMatch>),
    Export(std::result::Result<PathBuf, String>),
    LibraryDownload(std::result::Result<String, String>),
    Deconv(std::result::Result<DeconvResult, String>),
}

#[derive(Debug, Clone)]
pub enum Modal {
    None,
    Help,
    GoTo(String),
    Search(SearchState),
    FilePicker(FilePickerState),
    ImportConfig(ImportModalState),
    DemoPicker(usize),
    Export(ExportModalState),
    LibraryMatch(LibraryMatchState),
    CandidatePicker(CandidatePickerState),
    Busy(BusyState),
    LibrarySetup(LibrarySetupState),
    Deconv(DeconvModalState),
    SpectralMath(SpectralMathState),
    Integration(IntegrationModalState),
    Calibration(CalibrationModalState),
}

#[derive(Debug, Clone)]
pub struct SolventsTabState {
    pub solvent: SolventId,
    pub is_carbon: bool,
    pub query: String,
    pub editing: bool,
    pub tol_idx: usize, // 0: ±0.05, 1: ±0.12, 2: ±0.30, 3: All
    pub selected: usize,
}

impl SolventsTabState {
    pub fn tolerance(&self) -> f64 {
        let base = match self.tol_idx {
            0 => 0.05,
            1 => 0.12,
            2 => 0.30,
            _ => 999.0,
        };
        if self.is_carbon && self.tol_idx < 3 {
            base * 20.0
        } else {
            base
        }
    }

    pub fn tol_label(&self) -> &'static str {
        if self.is_carbon {
            match self.tol_idx {
                0 => "±1.0 ppm",
                1 => "±2.4 ppm",
                2 => "±6.0 ppm",
                _ => "All",
            }
        } else {
            match self.tol_idx {
                0 => "±0.05 ppm",
                1 => "±0.12 ppm",
                2 => "±0.30 ppm",
                _ => "All",
            }
        }
    }
}

#[derive(Debug, Clone)]
pub struct PolymerLabState {
    pub sub_tab: PolymerSubTab,
    pub field_idx: usize,
    pub editing: bool,
    // Sub 1: FTIR Conversion
    pub rxn_center: String,
    pub rxn_half: String,
    pub ref_center: String,
    pub ref_half: String,
    pub manual_mode: bool,
    pub man_rxn_0: String,
    pub man_ref_0: String,
    pub man_rxn_t: String,
    pub man_ref_t: String,
    // Sub 2: NMR Mn / DPn
    pub monomer_preset: usize,
    pub m0: String,
    pub m_end: String,
    pub n_repeat: String,
    pub n_end: String,
    pub i_repeat: String,
    pub i_end: String,
    // Sub 3: Tacticity
    pub i_mm: String,
    pub i_mr: String,
    pub i_rr: String,
    // Sub 4: PVDF β-phase & PU H-bonding
    pub phase_manual: bool,
    pub a840: String,
    pub a763: String,
    pub a1703: String,
    pub a1732: String,
    // Sub 5: UV-Vis Tauc band gap & Beer-Lambert
    /// true = direct allowed transition (αhν)², false = indirect (αhν)^½
    pub tauc_direct: bool,
    /// true = manual Beer-Lambert inputs, false = pull A at the cursor
    pub uv_manual: bool,
    pub uv_abs: String,
    pub uv_eps: String,
    pub uv_path_cm: String,
}

impl Default for PolymerLabState {
    fn default() -> Self {
        Self {
            sub_tab: PolymerSubTab::Conversion,
            field_idx: 0,
            editing: false,
            rxn_center: "1638".into(),
            rxn_half: "25".into(),
            ref_center: "1725".into(),
            ref_half: "35".into(),
            manual_mode: false,
            man_rxn_0: "0.58".into(),
            man_ref_0: "0.92".into(),
            man_rxn_t: "0.046".into(),
            man_ref_t: "0.95".into(),
            monomer_preset: 0,
            m0: "100.12".into(),
            m_end: "221.3".into(),
            n_repeat: "3".into(),
            n_end: "2".into(),
            i_repeat: "240.0".into(),
            i_end: "2.0".into(),
            i_mm: "14.0".into(),
            i_mr: "38.0".into(),
            i_rr: "62.0".into(),
            phase_manual: false,
            a840: "1.26".into(),
            a763: "1.00".into(),
            a1703: "0.72".into(),
            a1732: "0.45".into(),
            tauc_direct: true,
            uv_manual: false,
            uv_abs: "0.50".into(),
            uv_eps: "50000".into(),
            uv_path_cm: "1.0".into(),
        }
    }
}

pub const MONOMER_PRESETS: &[(&str, f64, &str)] = &[
    ("MMA (PMMA)", 100.12, "3"),
    ("Styrene (PS)", 104.15, "5"),
    ("EO (PEG)", 44.05, "4"),
    ("LA (PLA)", 72.06, "1"),
    ("ε-CL (PCL)", 114.14, "2"),
    ("VAc (PVAc)", 86.09, "1"),
    ("NIPAM (PNIPAM)", 113.16, "1"),
    ("Acrylamide", 71.08, "1"),
];

pub struct App {
    pub should_quit: bool,
    pub tab: ActiveTab,
    pub modal: Modal,

    pub traces: Vec<Spectrum>,
    pub processed: Vec<ProcessedTrace>,
    pub active_trace: usize,

    pub view_x_min: f64,
    pub view_x_max: f64,
    pub cursor_x: f64,

    pub process_cfg: ProcessConfig,
    pub stacked: bool,
    pub show_peaks: bool,
    pub peak_threshold: f64,
    pub auto_peaks: Vec<Peak>,
    pub pinned_peaks: Vec<Peak>,
    pub selected_peak: usize,

    pub nmr_solvent: SolventId,
    pub show_solvent_markers: bool,

    pub solvents_tab: SolventsTabState,
    pub polymer_lab: PolymerLabState,

    pub status_msg: String,
    pub hd_mode: bool,
    /// Scroll offset (in lines) for the beginner guide shown with `?`.
    pub help_scroll: usize,

    /// Receiver for the currently running background task (loading, matching, export).
    pub pending: Option<Receiver<TaskOutcome>>,
    /// Export settings preserved while the Python figure is generated.
    pub pending_export: Option<ExportModalState>,
    /// Active multi-peak deconvolution overlay on the spectrum plot.
    pub deconv_overlay: Option<DeconvOverlay>,
    /// Deconvolution modal state preserved while LM fitting runs in background.
    pub pending_deconv: Option<DeconvModalState>,
    /// Library-match modal state preserved while the NNLS mixture task runs.
    pub pending_match_state: Option<LibraryMatchState>,

    /// Quantitative integration regions (Mnova-style relative integrals).
    pub integration_regions: Vec<IntegrationRegion>,
    pub integration_ref_idx: usize,
    pub integration_ref_val: f64,
    pub next_region_id: usize,

    // Cached layout geometry for accurate mouse interactions
    pub topbar_rect: Rect,
    pub tab_hitboxes: Vec<(u16, u16, ActiveTab)>,
    pub plot_rect: Rect,
    pub drag_start_x: Option<f64>,
    pub drag_curr_x: Option<f64>,
}

impl App {
    pub fn new(files: Vec<PathBuf>, kind_override: Option<SpectrumKind>) -> Self {
        let mut traces: Vec<Spectrum> = Vec::new();
        let mut warnings: Vec<String> = Vec::new();
        let status_msg;

        for path in files.into_iter().take(8) {
            // A bare word such as `spectra butane` is a library name, not a path.
            let loaded = load_spectrum_file(&path, kind_override).or_else(|first_err| {
                let stem = path
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or_default();
                crate::io::resolve_library_query(stem, kind_override)
                    .and_then(|hit| load_spectrum_file(&hit.path, kind_override).ok())
                    .ok_or(first_err)
            });
            match loaded {
                Ok(mut s) => {
                    // FTIR uses cm⁻¹ (4000→400) while NMR uses ppm (12→0):
                    // overlaying them produces a meaningless 3966→0 axis, so an
                    // incompatible file replaces the workspace instead.
                    if let Some(first) = traces.first() {
                        if first.kind != s.kind {
                            warnings.push(format!(
                                "switched to {} — {} and {} axes cannot be overlaid",
                                s.kind.display_name(),
                                first.kind.display_name(),
                                s.kind.display_name()
                            ));
                            traces.clear();
                        }
                    }
                    s.color_idx = traces.len();
                    traces.push(s);
                }
                Err(e) => {
                    warnings.push(format!("{}: {}", path.display(), e));
                }
            }
        }

        let initial_tab = if let Some(first) = traces.first() {
            let base = format!("Loaded {} trace(s) · {}", traces.len(), first.name);
            status_msg = if warnings.is_empty() {
                base
            } else {
                format!("{} — ⚠ {}", base, warnings.join("; "))
            };
            ActiveTab::from_kind(first.kind)
        } else {
            // `spectra nmr` / `spectra cnmr` / `spectra ftir` with no query:
            // land on the matching tab with a representative demo preset.
            let preset = match kind_override {
                Some(SpectrumKind::Hnmr) => DemoPreset::PmmaHnmr,
                Some(SpectrumKind::Cnmr) => DemoPreset::PmmaCnmr,
                Some(SpectrumKind::UvVis) => DemoPreset::Tio2ZnoNanocompositeUvVis,
                _ => DemoPreset::MmaVsPmmaFtir,
            };
            traces = generate_demo_preset(preset);
            if warnings.is_empty() {
                status_msg = format!(
                    "Demo preset: {} · Press [?] for shortcuts, [d] Demo presets, [f] Open files",
                    preset.label()
                );
            } else {
                status_msg = format!(
                    "Could not load {} — showing demo preset instead",
                    warnings.join("; ")
                );
            }
            ActiveTab::from_kind(preset.kind())
        };

        let mut app = Self {
            should_quit: false,
            tab: initial_tab,
            modal: Modal::None,
            traces,
            processed: Vec::new(),
            active_trace: 0,
            view_x_min: 400.0,
            view_x_max: 4000.0,
            cursor_x: 1725.0,
            process_cfg: ProcessConfig::default(),
            stacked: false,
            show_peaks: true,
            peak_threshold: 0.08,
            auto_peaks: Vec::new(),
            pinned_peaks: Vec::new(),
            selected_peak: 0,
            nmr_solvent: SolventId::Cdcl3,
            show_solvent_markers: true,
            solvents_tab: SolventsTabState {
                solvent: SolventId::Cdcl3,
                is_carbon: false,
                query: String::new(),
                editing: false,
                tol_idx: 1,
                selected: 0,
            },
            polymer_lab: PolymerLabState::default(),
            status_msg,
            hd_mode: crate::kitty_plot::is_kitty_terminal(),
            help_scroll: 0,
            pending: None,
            pending_export: None,
            deconv_overlay: None,
            pending_deconv: None,
            pending_match_state: None,
            integration_regions: Vec::new(),
            integration_ref_idx: 0,
            integration_ref_val: 1.0,
            next_region_id: 0,
            topbar_rect: Rect::default(),
            tab_hitboxes: Vec::new(),
            plot_rect: Rect::default(),
            drag_start_x: None,
            drag_curr_x: None,
        };

        app.recompute_processing();
        app.reset_zoom();

        // If library is empty (fresh clone with no downloaded packs), show a helpful hint.
        // Only do this when there was no CLI error (warnings empty) to avoid spamming tests.
        if warnings.is_empty() {
            let lib_count = crate::io::library_index()
                .iter()
                .filter(|e| e.path.to_string_lossy().contains("library"))
                .count();
            if lib_count == 0 {
                app.status_msg.push_str(" · Library empty — press [D] to download reference spectra");
            }
        }

        app
    }

    pub fn current_kind(&self) -> SpectrumKind {
        self.tab
            .as_kind()
            .or_else(|| self.traces.get(self.active_trace).map(|t| t.kind))
            .unwrap_or(SpectrumKind::Ftir)
    }

    pub fn recompute_processing(&mut self) {
        self.processed = self
            .traces
            .iter()
            .map(|spec| {
                let (x, y, y_unit) = process_trace(spec, self.process_cfg);
                ProcessedTrace { x, y, y_unit }
            })
            .collect();
        self.recompute_peaks();
        self.recompute_integrations();
    }

    /// Recompute every integration region against the active processed trace so
    /// areas stay synchronized after baseline / ATR / smoothing / derivative changes.
    pub fn recompute_integrations(&mut self) {        let Some(pt) = self.processed.get(self.active_trace).cloned() else {
            return;
        };
        for r in &mut self.integration_regions {
            match integrate_spectrum_window(&pt.x, &pt.y, pt.y_unit, r.x_lo, r.x_hi) {
                Some((area, h)) => {
                    r.area = area;
                    r.peak_height = h;
                }
                None => {
                    r.area = 0.0;
                    r.peak_height = 0.0;
                }
            }
        }
    }

    fn clear_integrations(&mut self) {
        self.integration_regions.clear();
        self.integration_ref_idx = 0;
        self.integration_ref_val = 1.0;
        self.next_region_id = 0;
    }

    pub fn recompute_peaks(&mut self) {
        self.auto_peaks.clear();
        if let (Some(spec), Some(proc_t)) = (
            self.traces.get(self.active_trace),
            self.processed.get(self.active_trace),
        ) {
            let kind = spec.kind;
            self.auto_peaks = detect_peaks(
                &proc_t.x,
                &proc_t.y,
                kind,
                proc_t.y_unit,
                self.peak_threshold,
                16,
            );
            let peak_xs: Vec<f64> = self.auto_peaks.iter().map(|p| p.x).collect();
            for peak in &mut self.auto_peaks {
                peak.label =
                    best_peak_label_in_context(kind, peak.x, &peak_xs, self.nmr_solvent);

                // NMR: label detected solvent/impurity peaks explicitly, but never
                // hijack a genuine polymer signal that sits even closer.
                if matches!(kind, SpectrumKind::Hnmr | SpectrumKind::Cnmr) {
                    let is_carbon = kind == SpectrumKind::Cnmr;
                    let tol = if is_carbon { 0.8 } else { 0.06 };
                    if let Some((imp, _shift, imp_dist)) =
                        match_solvent_impurities(self.nmr_solvent, is_carbon, peak.x, tol)
                            .first()
                    {
                        let poly_dist = POLYMERS_DB
                            .iter()
                            .flat_map(|p| {
                                let peaks = match kind {
                                    SpectrumKind::Ftir => p.ftir_peaks,
                                    SpectrumKind::Hnmr => p.hnmr_peaks,
                                    SpectrumKind::Cnmr => p.cnmr_peaks,
                                    SpectrumKind::UvVis => &[],
                                };
                                peaks.iter().map(|(x, _)| (peak.x - *x).abs())
                            })
                            .fold(f64::INFINITY, f64::min);
                        if *imp_dist <= poly_dist {
                            let tag = if imp.name.contains("Residual") {
                                self.nmr_solvent.label()
                            } else {
                                imp.name.split(" (").next().unwrap_or(imp.name)
                            };
                            peak.label = if is_carbon {
                                format!("{:.1} {} [solvent/imp]", peak.x, tag)
                            } else {
                                format!("{:.2} {} [solvent/imp]", peak.x, tag)
                            };
                        }
                    }
                }
            }
            for p in &mut self.pinned_peaks {
                if let Some(new_y) = interpolate_y(&proc_t.x, &proc_t.y, p.x) {
                    p.y = new_y;
                }
            }
        }
        if self.selected_peak >= self.auto_peaks.len() {
            self.selected_peak = 0;
        }
    }

    pub fn full_bounds(&self) -> (f64, f64) {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for t in &self.traces {
            if t.visible && !t.x.is_empty() {
                let (a, b) = t.x_range();
                lo = lo.min(a);
                hi = hi.max(b);
            }
        }
        if lo.is_finite() && hi.is_finite() && lo < hi {
            (lo, hi)
        } else {
            self.current_kind().default_bounds()
        }
    }

    pub fn reset_zoom(&mut self) {
        let (lo, hi) = self.full_bounds();
        self.view_x_min = lo;
        self.view_x_max = hi;
        if self.cursor_x < lo || self.cursor_x > hi {
            self.cursor_x = match self.current_kind() {
                SpectrumKind::Ftir => 1725.0_f64.clamp(lo, hi),
                SpectrumKind::Hnmr => 3.60_f64.clamp(lo, hi),
                SpectrumKind::Cnmr => 77.16_f64.clamp(lo, hi),
                SpectrumKind::UvVis => 300.0_f64.clamp(lo, hi),
            };
        }
    }

    pub fn switch_tab(&mut self, next: ActiveTab) {
        self.tab = next;
        // When switching to an analytical tab whose kind differs from currently loaded demo traces,
        // auto-load a matching built-in reference if no external file was loaded.
        if let Some(target_kind) = next.as_kind() {
            let has_matching = self.traces.iter().any(|t| t.kind == target_kind);
            let all_builtin = self.traces.iter().all(|t| t.source_path.is_none());
            if !has_matching && all_builtin {
                let preset = match target_kind {
                    SpectrumKind::Ftir => DemoPreset::MmaVsPmmaFtir,
                    SpectrumKind::Hnmr => DemoPreset::PmmaHnmr,
                    SpectrumKind::Cnmr => DemoPreset::PmmaCnmr,
                    SpectrumKind::UvVis => DemoPreset::Tio2ZnoNanocompositeUvVis,
                };
                self.load_demo(preset);
            } else {
                self.reset_zoom();
            }
        }
    }

    pub fn load_demo(&mut self, preset: DemoPreset) {
        self.traces = generate_demo_preset(preset);
        self.active_trace = 0;
        self.pinned_peaks.clear();
        self.deconv_overlay = None;
        self.clear_integrations();
        self.tab = ActiveTab::from_kind(preset.kind());
        self.recompute_processing();
        self.reset_zoom();
        self.status_msg = format!("Loaded preset: {}", preset.label());
    }

    // ========================================================================
    // BACKGROUND TASKS (loading / matching / exporting never block the UI)
    // ========================================================================
    fn start_task(
        &mut self,
        label: impl Into<String>,
        work: impl FnOnce() -> TaskOutcome + Send + 'static,
    ) {
        if self.pending.is_some() {
            self.status_msg = "A background task is still running — please wait…".into();
            return;
        }
        let (tx, rx) = mpsc::channel();
        thread::spawn(move || {
            let _ = tx.send(work());
        });
        self.pending = Some(rx);
        self.modal = Modal::Busy(BusyState {
            label: label.into(),
            started: Instant::now(),
        });
    }

    /// Called from the main loop; applies a finished background task, if any.
    pub fn poll_background(&mut self) {
        let Some(rx) = &self.pending else {
            return;
        };
        match rx.try_recv() {
            Ok(outcome) => {
                self.pending = None;
                self.apply_task_outcome(outcome);
            }
            Err(mpsc::TryRecvError::Empty) => {}
            Err(mpsc::TryRecvError::Disconnected) => {
                self.pending = None;
                if matches!(self.modal, Modal::Busy(_)) {
                    self.modal = Modal::None;
                }
                self.status_msg = "Background task failed unexpectedly".into();
            }
        }
    }

    /// Merge freshly loaded spectra into the workspace while enforcing a single
    /// spectroscopy mode per plot. FTIR (cm⁻¹) and NMR (ppm) axes cannot be
    /// overlaid, so an incompatible append replaces the workspace instead.
    /// Returns (loaded names, notes for the status line).
    fn absorb_specs(
        &mut self,
        specs: Vec<Spectrum>,
        append: bool,
    ) -> (Vec<String>, Vec<String>) {
        let mut notes = Vec::new();
        let mut names = Vec::new();
        if specs.is_empty() {
            return (names, notes);
        }

        let appending = append && !self.traces.is_empty();
        if appending {
            let base_kind = self.traces[0].kind;
            let mut same: Vec<Spectrum> = Vec::new();
            let mut other: Vec<Spectrum> = Vec::new();
            for spec in specs {
                if spec.kind == base_kind {
                    same.push(spec);
                } else {
                    other.push(spec);
                }
            }

            if same.is_empty() {
                // Nothing compatible with the open axes: switch mode cleanly.
                let new_kind = other[0].kind;
                self.traces.clear();
                self.pinned_peaks.clear();
                self.deconv_overlay = None;
                self.clear_integrations();
                for (k, mut spec) in other
                    .into_iter()
                    .filter(|s| s.kind == new_kind)
                    .enumerate()
                {
                    names.push(spec.name.clone());
                    spec.color_idx = k;
                    self.traces.push(spec);
                }
                self.active_trace = 0;
                notes.push(format!(
                    "switched to {} — {} and {} axes cannot be overlaid",
                    new_kind.display_name(),
                    base_kind.display_name(),
                    new_kind.display_name()
                ));
            } else {
                for mut spec in same {
                    if let Some(path) = spec.source_path.as_deref() {
                        if let Some(idx) = self
                            .traces
                            .iter()
                            .position(|t| t.source_path.as_deref() == Some(path))
                        {
                            self.active_trace = idx;
                            names.push(spec.name.clone());
                            notes.push(format!("already open: {}", spec.name));
                            continue;
                        }
                    }
                    names.push(spec.name.clone());
                    spec.color_idx = self.traces.len();
                    self.traces.push(spec);
                    self.active_trace = self.traces.len() - 1;
                }
                if !other.is_empty() {
                    notes.push(format!(
                        "skipped {} {} trace(s): cannot overlay {} and {} axes",
                        other.len(),
                        other[0].kind.display_name(),
                        base_kind.display_name(),
                        other[0].kind.display_name()
                    ));
                }
            }
        } else {
            let kind = specs[0].kind;
            let total = specs.len();
            let mut replacement: Vec<Spectrum> = Vec::new();
            for mut spec in specs {
                if spec.kind != kind {
                    continue;
                }
                names.push(spec.name.clone());
                spec.color_idx = replacement.len();
                replacement.push(spec);
            }
            let skipped = total.saturating_sub(replacement.len());
            self.traces = replacement;
            self.active_trace = 0;
            self.pinned_peaks.clear();
            self.deconv_overlay = None;
            self.clear_integrations();
            if skipped > 0 {
                notes.push(format!(
                    "skipped {} trace(s) of a different spectroscopy mode",
                    skipped
                ));
            }
        }

        if self.traces.is_empty() {
            self.auto_peaks.clear();
            self.clear_integrations();
        } else {
            let active = self.active_trace.min(self.traces.len() - 1);
            self.active_trace = active;
            self.tab = ActiveTab::from_kind(self.traces[active].kind);
            self.recompute_processing();
            self.reset_zoom();
        }
        (names, notes)
    }

    fn apply_task_outcome(&mut self, outcome: TaskOutcome) {
        match outcome {
            TaskOutcome::Loaded {
                specs,
                errors,
                append,
            } => {
                self.modal = Modal::None;
                if specs.is_empty() {
                    self.status_msg = if errors.is_empty() {
                        "Nothing to load".into()
                    } else {
                        format!("Load failed — {}", errors.join("; "))
                    };
                    return;
                }

                let (names, mut notes) = self.absorb_specs(specs, append);
                notes.extend(errors);
                let base = format!(
                    "Loaded {} · {}",
                    if names.len() == 1 { "spectrum" } else { "spectra" },
                    names.join(", ")
                );
                self.status_msg = if notes.is_empty() {
                    base
                } else {
                    format!("{} — ⚠ {}", base, notes.join("; "))
                };
            }
            TaskOutcome::LibraryMatch(results) => {
                if results.is_empty() {
                    self.modal = Modal::None;
                    self.status_msg =
                        "No library matches — check that the active FTIR trace has data".into();
                    return;
                }
                self.status_msg = format!("Library match: {} candidates ranked by HQI", results.len());
                self.modal = Modal::LibraryMatch(LibraryMatchState {
                    results,
                    selected: 0,
                    mixture_mode: false,
                    mixture_results: Vec::new(),
                    mixture_selected: 0,
                });
            }
            TaskOutcome::MixtureMatch(mix) => {
                let mut st = self.pending_match_state.take().unwrap_or(LibraryMatchState {
                    results: Vec::new(),
                    selected: 0,
                    mixture_mode: true,
                    mixture_results: Vec::new(),
                    mixture_selected: 0,
                });
                if mix.is_empty() {
                    self.modal = Modal::LibraryMatch(st);
                    self.status_msg =
                        "No multi-component mixture fit — try a wider spectral window".into();
                    return;
                }
                self.status_msg = format!(
                    "Mixture NNLS: {} candidate blend(s) — top {:.1}% composite HQI",
                    mix.len(),
                    mix[0].composite_hqi_pct
                );
                st.mixture_mode = true;
                st.mixture_results = mix;
                st.mixture_selected = 0;
                self.modal = Modal::LibraryMatch(st);
            }
            TaskOutcome::Export(result) => {
                let mut state = self.pending_export.take();
                match (result, state.as_mut()) {
                    (Ok(saved), Some(st)) => {
                        st.output_path = saved.display().to_string();
                        st.last_exported = Some(saved.clone());
                        st.feedback = Some(format!(
                            "Saved: {}  ([y] Copy PNG | [o] Open)",
                            saved.display()
                        ));
                        self.status_msg = format!("Exported figure -> {}", saved.display());
                        self.modal = Modal::Export(st.clone());
                    }
                    (Err(e), Some(st)) => {
                        st.feedback = Some(format!("Error: {}", e));
                        self.modal = Modal::Export(st.clone());
                    }
                    (Ok(saved), None) => {
                        self.modal = Modal::None;
                        self.status_msg = format!("Exported figure -> {}", saved.display());
                    }
                    (Err(e), None) => {
                        self.modal = Modal::None;
                        self.status_msg = format!("Export failed: {}", e);
                    }
                }
            }
            TaskOutcome::LibraryDownload(result) => {
                match result {
                    Ok(summary) => {
                        // Refresh the library index so new files are immediately visible
                        crate::io::refresh_library_index();
                        self.modal = Modal::LibrarySetup(LibrarySetupState {
                            selected: 0,
                            checks: [false, false, false, false, false], // uncheck all after success
                            feedback: Some(summary),
                        });
                        self.status_msg = "Library download complete — press [M] to match or [/] to search".into();
                    }
                    Err(e) => {
                        self.modal = Modal::LibrarySetup(LibrarySetupState {
                            selected: 0,
                            checks: [true, true, false, false, false],
                            feedback: Some(format!("Error: {}", e)),
                        });
                        self.status_msg = format!("Library download failed: {}", e);
                    }
                }
            }
            TaskOutcome::Deconv(result) => {
                let mut st = self.pending_deconv.take().unwrap_or(DeconvModalState {
                    field_idx: 0,
                    shape: PeakShape::PseudoVoigt,
                    max_peaks: 3,
                    fit_baseline: true,
                    x_lo: format!("{:.2}", self.view_x_min),
                    x_hi: format!("{:.2}", self.view_x_max),
                    result: None,
                    feedback: None,
                });
                match result {
                    Ok(res) => {
                        let summary = format!(
                            "Fitted {} peak(s) | R² = {:.5} | RMSE = {:.4}",
                            res.peaks.len(),
                            res.r_squared,
                            res.rmse
                        );
                        st.result = Some(res.clone());
                        st.feedback = Some(summary.clone());
                        self.deconv_overlay = Some(DeconvOverlay {
                            trace_idx: self.active_trace,
                            result: res,
                            show_components: true,
                            show_residual: false,
                        });
                        self.status_msg = format!("Deconvolution: {}", summary);
                        self.modal = Modal::Deconv(st);
                    }
                    Err(e) => {
                        st.feedback = Some(format!("Fit error: {}", e));
                        self.status_msg = format!("Deconvolution failed: {}", e);
                        self.modal = Modal::Deconv(st);
                    }
                }
            }
        }
    }

    /// Load one or more spectrum files on a background thread.
    pub fn load_paths(
        &mut self,
        mut paths: Vec<PathBuf>,
        kind: Option<SpectrumKind>,
        append: bool,
        label: impl Into<String>,
    ) {
        if paths.is_empty() {
            return;
        }
        if append {
            // Appending the same file twice would create an invisible duplicate;
            // focus the trace the user already has instead.
            let mut focused: Option<usize> = None;
            paths.retain(|p| {
                if let Some(idx) = self.traces.iter().position(|t| {
                    t.source_path
                        .as_deref()
                        .map(|s| std::path::Path::new(s) == p.as_path())
                        .unwrap_or(false)
                }) {
                    focused = Some(idx);
                    false
                } else {
                    true
                }
            });
            if paths.is_empty() {
                if let Some(idx) = focused {
                    self.active_trace = idx;
                    self.recompute_peaks();
                    self.status_msg = format!("Trace already open: {}", self.traces[idx].name);
                }
                return;
            }
        }
        self.start_task(label, move || {
            let mut specs = Vec::new();
            let mut errors = Vec::new();
            for p in paths {
                match load_spectrum_file(&p, kind) {
                    Ok(spec) => specs.push(spec),
                    Err(e) => errors.push(format!("{}: {}", p.display(), e)),
                }
            }
            TaskOutcome::Loaded {
                specs,
                errors,
                append,
            }
        });
    }

    fn import_state_for(path: PathBuf, append_mode: bool) -> Modal {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        let is_image = matches!(
            ext.as_str(),
            "png" | "jpg" | "jpeg" | "webp" | "bmp" | "tif" | "tiff"
        );
        Modal::ImportConfig(ImportModalState {
            path,
            is_image,
            append_mode,
            field_idx: 0,
            kind_choice: if is_image { 0 } else { 4 },
            x_left: "4000".into(),
            x_right: "400".into(),
            bg_idx: 0,
            error: None,
        })
    }

    // ========================================================================
    // PASTE / DRAG-AND-DROP
    // ========================================================================
    /// Kitty drag-and-drop and terminal paste arrive as one Paste event.
    /// Without this, a dropped path is replayed as keystrokes: the leading `/`
    /// opens Global Search and the rest of the path becomes a useless query.
    pub fn handle_paste(&mut self, text: &str) {
        let pasted = text.trim();
        if pasted.is_empty() {
            return;
        }

        match &mut self.modal {
            Modal::Search(st) => {
                st.query.push_str(pasted);
                st.results = build_search_results(&st.query);
                st.selected = 0;
                return;
            }
            Modal::GoTo(buf) => {
                for c in pasted.chars() {
                    if c.is_ascii_digit() || c == '.' || c == '-' {
                        buf.push(c);
                    }
                }
                return;
            }
            Modal::Export(st) => {
                if st.field_idx == 5 {
                    st.title.push_str(pasted);
                } else if st.field_idx == 6 {
                    st.output_path = pasted.to_string();
                }
                return;
            }
            Modal::Deconv(_) | Modal::SpectralMath(_) | Modal::Integration(_) | Modal::Calibration(_) => {
                return;
            }
            Modal::FilePicker(_) => {
                self.paste_into_file_picker(pasted);
                return;
            }
            Modal::None => {}
            Modal::Help
            | Modal::DemoPicker(_)
            | Modal::ImportConfig(_)
            | Modal::LibraryMatch(_)
            | Modal::CandidatePicker(_)
            | Modal::LibrarySetup(_)
            | Modal::Busy(_) => return,
        }

        self.open_pasted_spectra(pasted);
    }

    fn open_pasted_spectra(&mut self, text: &str) {
        let candidates = parse_pasted_paths(text);
        let mut files: Vec<PathBuf> = candidates.iter().filter(|p| p.is_file()).cloned().collect();

        if files.is_empty() {
            // Not an existing path: allow pasting a library name such as "butane".
            let query = text
                .trim()
                .trim_matches('"')
                .trim_matches('\'')
                .to_string();
            if let Some(hit) = crate::io::resolve_library_query(&query, None) {
                files.push(hit.path);
            } else {
                let shown = candidates
                    .iter()
                    .map(|p| p.display().to_string())
                    .collect::<Vec<_>>()
                    .join(", ");
                self.status_msg = format!("Dropped / pasted item not found: {}", shown);
                return;
            }
        }

        if files.len() == 1 {
            if let Some(ext) = files[0].extension().and_then(|s| s.to_str()) {
                if matches!(
                    ext.to_ascii_lowercase().as_str(),
                    "png" | "jpg" | "jpeg" | "webp" | "bmp" | "tif" | "tiff"
                ) {
                    self.modal = Self::import_state_for(files[0].clone(), false);
                    return;
                }
            }
        }

        // Dropping a file that is already open should focus it, not reload it.
        if files.len() == 1 {
            let target = files[0].as_path();
            if let Some(idx) = self.traces.iter().position(|t| {
                t.source_path
                    .as_deref()
                    .map(|s| std::path::Path::new(s) == target)
                    .unwrap_or(false)
            }) {
                self.active_trace = idx;
                self.recompute_peaks();
                self.status_msg = format!("Trace already open: {}", self.traces[idx].name);
                return;
            }
        }

        let label = if files.len() == 1 {
            "Opening dropped spectrum…".to_string()
        } else {
            format!("Opening {} dropped spectra…", files.len())
        };
        self.load_paths(files, None, false, label);
    }

    fn paste_into_file_picker(&mut self, text: &str) {
        let append_mode = match &self.modal {
            Modal::FilePicker(st) => st.append_mode,
            _ => return,
        };
        let Some(target) = parse_pasted_paths(text).into_iter().next() else {
            return;
        };
        if target.is_file() {
            if let Some(parent) = target.parent() {
                crate::io::remember_picker_dir(parent);
            }
            self.modal = Self::import_state_for(target, append_mode);
        } else if target.is_dir() {
            if let Modal::FilePicker(st) = &mut self.modal {
                st.set_dir(target);
            }
        } else if let Some(hit) = crate::io::resolve_library_query(text.trim(), None) {
            self.modal = Self::import_state_for(hit.path, append_mode);
        } else if let Modal::FilePicker(st) = &mut self.modal {
            if text.contains('/') || text.trim_start().starts_with('~') {
                st.path_input = Some(target.display().to_string());
            } else {
                // Plain text becomes the live search, same as typing it.
                st.filter = text.trim().to_string();
                st.path_input = None;
                st.change_filter();
            }
        }
    }

    pub fn zoom_around_cursor(&mut self, factor: f64) {
        let (full_lo, full_hi) = self.full_bounds();
        let full_span = (full_hi - full_lo).max(1e-6);
        let span = (self.view_x_max - self.view_x_min).max(1e-6);
        let new_span = (span * factor).clamp(full_span * 0.01, full_span);
        let frac = ((self.cursor_x - self.view_x_min) / span).clamp(0.05, 0.95);
        let mut lo = self.cursor_x - frac * new_span;
        let mut hi = lo + new_span;
        if lo < full_lo {
            lo = full_lo;
            hi = (lo + new_span).min(full_hi);
        }
        if hi > full_hi {
            hi = full_hi;
            lo = (hi - new_span).max(full_lo);
        }
        self.view_x_min = lo;
        self.view_x_max = hi;
    }

    pub fn move_cursor(&mut self, screen_dir: f64, coarse: bool) {
        // FTIR/NMR plot with the X axis reversed (IUPAC), so moving right (+1)
        // DECREASES wavenumber / ppm. UV-Vis runs ascending in nm.
        let span = (self.view_x_max - self.view_x_min).max(1e-6);
        let step_frac = if coarse { 0.04 } else { 0.008 };
        let delta = if self.current_kind().is_x_reversed() {
            -screen_dir * span * step_frac
        } else {
            screen_dir * span * step_frac
        };
        let (full_lo, full_hi) = self.full_bounds();
        self.cursor_x = (self.cursor_x + delta).clamp(full_lo, full_hi);

        // Pan window smoothly if cursor exits current view
        if self.cursor_x < self.view_x_min {
            let shift = self.view_x_min - self.cursor_x;
            self.view_x_min = (self.view_x_min - shift).max(full_lo);
            self.view_x_max = (self.view_x_min + span).min(full_hi);
        } else if self.cursor_x > self.view_x_max {
            let shift = self.cursor_x - self.view_x_max;
            self.view_x_max = (self.view_x_max + shift).min(full_hi);
            self.view_x_min = (self.view_x_max - span).max(full_lo);
        }
    }

    pub fn jump_to_x(&mut self, target: f64) {
        let (lo, hi) = self.full_bounds();
        let clamped = target.clamp(lo, hi);
        self.cursor_x = clamped;
        if clamped < self.view_x_min || clamped > self.view_x_max {
            let span = (self.view_x_max - self.view_x_min).max((hi - lo) * 0.25);
            let mut nlo = (clamped - 0.5 * span).max(lo);
            let nhi = (nlo + span).min(hi);
            if nhi - nlo < span {
                nlo = (nhi - span).max(lo);
            }
            self.view_x_min = nlo;
            self.view_x_max = nhi;
        }
    }

    pub fn toggle_pin_at_cursor(&mut self) {
        let kind = self.current_kind();
        let tol = match kind {
            SpectrumKind::Ftir => 12.0,
            SpectrumKind::Hnmr => 0.08,
            SpectrumKind::Cnmr => 1.0,
            SpectrumKind::UvVis => 8.0,
        };
        if let Some(pos) = self.pinned_peaks.iter().position(|p| (p.x - self.cursor_x).abs() <= tol) {
            let removed = self.pinned_peaks.remove(pos);
            self.status_msg = format!("Unpinned peak at {:.2} {}", removed.x, kind.x_unit());
            return;
        }
        let y = self
            .processed
            .get(self.active_trace)
            .and_then(|pt| interpolate_y(&pt.x, &pt.y, self.cursor_x))
            .unwrap_or(0.0);
        let peak_xs: Vec<f64> = self.auto_peaks.iter().map(|p| p.x).collect();
        let label = best_peak_label_in_context(kind, self.cursor_x, &peak_xs, self.nmr_solvent);
        self.pinned_peaks.push(Peak {
            x: self.cursor_x,
            y,
            prominence: 1.0,
            fwhm: None,
            snr: 0.0,
            label: label.clone(),
        });
        if kind.is_x_reversed() {
            self.pinned_peaks.sort_by(|a, b| b.x.total_cmp(&a.x));
        } else {
            self.pinned_peaks.sort_by(|a, b| a.x.total_cmp(&b.x));
        }
        self.status_msg = format!("Pinned ★ {}", label);
    }

    pub fn copy_cursor_summary(&mut self) {
        let kind = self.current_kind();
        let text = match kind {
            SpectrumKind::Ftir => {
                let top = match_ftir(self.cursor_x)
                    .first()
                    .map(|b| format!("{} ({}, {})", b.group, b.intensity, b.mode))
                    .unwrap_or_else(|| "Unassigned band".into());
                format!("FTIR: {:.1} cm⁻¹ — {}", self.cursor_x, top)
            }
            SpectrumKind::Hnmr => {
                let top = match_hnmr(self.cursor_x)
                    .first()
                    .map(|b| format!("{} [{}]", b.group, b.example))
                    .unwrap_or_else(|| "Unassigned".into());
                let solv = match_solvent_impurities(self.nmr_solvent, false, self.cursor_x, 0.08)
                    .first()
                    .map(|(e, s, _)| format!(" | Possible impurity: {} ({:.2} ppm)", e.name, s))
                    .unwrap_or_default();
                format!("¹H NMR ({}): δ {:.2} ppm — {}{}", self.nmr_solvent.label(), self.cursor_x, top, solv)
            }
            SpectrumKind::Cnmr => {
                let top = match_cnmr(self.cursor_x)
                    .first()
                    .map(|b| format!("{} [{}]", b.group, b.example))
                    .unwrap_or_else(|| "Unassigned".into());
                format!("¹³C NMR ({}): δ {:.1} ppm — {}", self.nmr_solvent.label(), self.cursor_x, top)
            }
            SpectrumKind::UvVis => {
                let ev = crate::data::photon_energy_ev(self.cursor_x);
                let (region, color) = uvvis_color_info(self.cursor_x);
                let top = match_uvvis(self.cursor_x)
                    .first()
                    .map(|b| format!("{} ({}, ε {})", b.group, b.transition, b.epsilon))
                    .unwrap_or_else(|| "Unassigned".into());
                format!(
                    "UV-Vis: λ {:.1} nm (E = {:.2} eV) — {} · observed {}. {}",
                    self.cursor_x, ev, region, color, top
                )
            }
        };
        match copy_text(&text) {
            Ok(tool) => self.status_msg = format!("Copied via {}: {}", tool, text),
            Err(e) => self.status_msg = format!("Copy failed: {}", e),
        }
    }

    pub fn open_file_picker(&mut self, append_mode: bool) {
        // Resume where the user last browsed; on a first run fall back to the
        // folder of the active trace, then to the process working directory.
        let cwd = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
        let active_dir = self
            .traces
            .get(self.active_trace)
            .and_then(|t| t.source_path.as_deref())
            .map(PathBuf::from)
            .and_then(|p| p.parent().map(|p| p.to_path_buf()));
        let dir = crate::io::remembered_picker_dir()
            .filter(|p| p.is_dir())
            .or_else(|| active_dir.filter(|p| p.is_dir()))
            .unwrap_or(cwd);
        self.modal = Modal::FilePicker(FilePickerState::new(dir, append_mode));
    }

    pub fn open_export_modal(&mut self) {
        let stamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_secs() % 100_000)
            .unwrap_or(0);
        let stem = self
            .traces
            .get(self.active_trace)
            .map(|t| {
                t.name
                    .chars()
                    .map(|c| if c.is_ascii_alphanumeric() { c.to_ascii_lowercase() } else { '_' })
                    .collect::<String>()
            })
            .unwrap_or_else(|| "spectrum".into());
        let title = if self.traces.len() > 1 {
            format!("{} Comparison", self.current_kind().display_name())
        } else {
            self.traces
                .first()
                .map(|t| t.name.clone())
                .unwrap_or_else(|| "Spectrum".into())
        };
        self.modal = Modal::Export(ExportModalState {
            field_idx: 0,
            format: ExportFormat::Png,
            style: ExportStyle::Acs,
            stacked: self.stacked || self.traces.len() > 1,
            use_zoom_range: true,
            peak_mode: ExportPeakMode::PinnedAndAuto,
            title,
            output_path: format!("exports/{}_{}.png", stem.trim_matches('_'), stamp),
            last_exported: None,
            feedback: None,
        });
    }

    pub fn open_search_modal(&mut self) {
        let results = build_search_results("");
        self.modal = Modal::Search(SearchState {
            query: String::new(),
            results,
            selected: 0,
        });
    }

    /// Ambiguous CLI query (e.g. `spectra butane`): let the user pick the exact
    /// reference instead of silently guessing one.
    pub fn open_candidates(&mut self, query: String, entries: Vec<crate::io::LibraryEntry>) {
        if entries.is_empty() {
            return;
        }
        self.status_msg = format!(
            "{} reference spectra match '{}' — ↑/↓ select · Enter open · Esc dismiss",
            entries.len(),
            query
        );
        self.modal = Modal::CandidatePicker(CandidatePickerState {
            entries,
            selected: 0,
        });
    }

    pub fn open_library_match(&mut self) {
        let Some(spec) = self.traces.get(self.active_trace) else {
            self.status_msg = "Load a spectrum before running library matching".into();
            return;
        };
        let kind = spec.kind;
        // Matching parses the whole reference pack on cold start, which takes seconds
        // on a slow/external disk — never do it on the UI thread.
        let spec = spec.clone();
        let cfg = self.process_cfg;
        let solvent = self.nmr_solvent;
        self.start_task("Matching against the reference library…", move || {
            TaskOutcome::LibraryMatch(match kind {
                SpectrumKind::Ftir => crate::matcher::match_ftir_library(&spec, cfg),
                SpectrumKind::UvVis => crate::matcher::match_uvvis_library(&spec, cfg),
                SpectrumKind::Hnmr | SpectrumKind::Cnmr => {
                    crate::matcher::match_nmr_library(&spec, cfg, solvent)
                }
            })
        });
    }

    pub fn apply_library_match(&mut self, m: &crate::matcher::LibraryMatch) {
        let kind = self.current_kind();
        let loaded = if m.is_builtin {
            match kind {
                SpectrumKind::UvVis => {
                    // Built-in UV-Vis references are demo presets, not polymers.
                    DemoPreset::ALL
                        .iter()
                        .filter(|p| p.kind() == SpectrumKind::UvVis)
                        .flat_map(|p| generate_demo_preset(*p))
                        .find(|s| s.name == m.name)
                }
                SpectrumKind::Ftir | SpectrumKind::Hnmr | SpectrumKind::Cnmr => {
                    find_polymer_by_short(&m.short)
                        .map(|p| synthesize_polymer_reference(p, kind))
                }
            }
        } else {
            m.file_path
                .as_deref()
                .and_then(|p| load_spectrum_file(p, Some(kind)).ok())
        };

        match loaded {
            Some(mut spec) => {
                spec.visible = true;
                let label = spec.name.clone();
                let unit = spec.raw_y_unit;
                if kind == SpectrumKind::UvVis {
                    // UV-Vis: Swapped renders %T, Native renders Abs.
                    self.process_cfg.y_pref = if unit == YUnit::Transmittance {
                        YAxisPreference::Swapped
                    } else {
                        YAxisPreference::Native
                    };
                } else if kind == SpectrumKind::Ftir && unit == YUnit::Transmittance {
                    self.process_cfg.y_pref = YAxisPreference::Native;
                }
                let (_, notes) = self.absorb_specs(vec![spec], true);
                self.status_msg = format!(
                    "Overlaid library match: {} ({:.1}% HQI)",
                    label, m.hqi_pct
                );
                if !notes.is_empty() {
                    self.status_msg =
                        format!("{} — ⚠ {}", self.status_msg, notes.join("; "));
                }
            }
            None => {
                self.status_msg = format!("Could not load overlay for {}", m.name);
            }
        }
    }

    /// Overlay every component of a multi-component NNLS mixture onto the workspace.
    fn overlay_mixture_match(&mut self, m: &crate::mixture::MixtureMatch) {
        let mut parts: Vec<String> = Vec::new();
        for c in &m.components {
            parts.push(format!("{:.0}% {}", c.weight_pct, c.short));
            let lm = crate::matcher::LibraryMatch {
                name: c.name.clone(),
                short: c.short.clone(),
                source: c.source.clone(),
                hqi_pct: m.composite_hqi_pct,
                r_abs: m.r_abs,
                r_deriv: m.r_deriv,
                is_builtin: c.is_builtin,
                file_path: c.file_path.clone(),
            };
            self.apply_library_match(&lm);
        }
        self.status_msg = format!(
            "Overlaid mixture: {} ({:.1}% HQI)",
            parts.join(" + "),
            m.composite_hqi_pct
        );
    }

    fn integration_ref_area(&self) -> Option<f64> {
        self.integration_regions
            .get(self.integration_ref_idx)
            .map(|r| r.area)
            .filter(|a| a.is_finite() && *a > 1e-12)
    }

    pub fn integration_relative(&self, idx: usize) -> Option<f64> {
        let area = self.integration_regions.get(idx)?.area;
        let ref_area = self.integration_ref_area()?;
        Some(area / ref_area * self.integration_ref_val)
    }

    fn add_integration_region(&mut self, x_lo: f64, x_hi: f64) -> Option<usize> {
        let pt = self.processed.get(self.active_trace).cloned()?;
        let lo = x_lo.min(x_hi);
        let hi = x_lo.max(x_hi);
        if !(hi - lo).is_finite() || hi - lo < 1e-9 {
            return None;
        }
        let (area, peak_height) =
            integrate_spectrum_window(&pt.x, &pt.y, pt.y_unit, lo, hi)?;
        let kind = self.current_kind();
        let center = 0.5 * (lo + hi);
        let peak_xs: Vec<f64> = self.auto_peaks.iter().map(|p| p.x).collect();
        let label = best_peak_label_in_context(kind, center, &peak_xs, self.nmr_solvent);
        let id = self.next_region_id;
        self.next_region_id += 1;
        let was_empty = self.integration_regions.is_empty();
        self.integration_regions.push(IntegrationRegion {
            id,
            x_lo: lo,
            x_hi: hi,
            area,
            peak_height,
            label,
        });
        let idx = self.integration_regions.len() - 1;
        if was_empty {
            self.integration_ref_idx = idx;
            if kind == SpectrumKind::Hnmr {
                self.integration_ref_val = 3.0;
            }
        }
        Some(idx)
    }

    /// Push the current relative integrals into Polymer Lab. Returns a status
    /// message when the active technique has a matching calculator.
    fn push_integrations_to_polymer_lab(&mut self) -> Option<String> {
        let kind = self.current_kind();
        let len = self.integration_regions.len();
        if len < 2 {
            return None;
        }
        let mut order: Vec<usize> = (0..len).collect();
        order.sort_by(|&a, &b| {
            let ca = 0.5 * (self.integration_regions[a].x_lo + self.integration_regions[a].x_hi);
            let cb = 0.5 * (self.integration_regions[b].x_lo + self.integration_regions[b].x_hi);
            cb.total_cmp(&ca)
        });

        match kind {
            SpectrumKind::Hnmr if len >= 3 => {
                let vals: Vec<f64> = order
                    .iter()
                    .take(3)
                    .map(|&i| self.integration_relative(i).unwrap_or(0.0))
                    .collect();
                self.polymer_lab.i_mm = format!("{:.4}", vals[0]);
                self.polymer_lab.i_mr = format!("{:.4}", vals[1]);
                self.polymer_lab.i_rr = format!("{:.4}", vals[2]);
                self.polymer_lab.sub_tab = PolymerSubTab::Tacticity;
                self.tab = ActiveTab::PolymerLab;
                Some(format!(
                    "Pushed 3 integrals to Tacticity: mm={:.3}, mr={:.3}, rr={:.3}",
                    vals[0], vals[1], vals[2]
                ))
            }
            SpectrumKind::Hnmr => {
                let a = self.integration_relative(order[0]).unwrap_or(0.0);
                let b = self.integration_relative(order[1]).unwrap_or(0.0);
                let (repeat, end) = if a >= b { (a, b) } else { (b, a) };
                self.polymer_lab.i_repeat = format!("{:.4}", repeat);
                self.polymer_lab.i_end = format!("{:.4}", end);
                self.polymer_lab.sub_tab = PolymerSubTab::NmrMn;
                self.tab = ActiveTab::PolymerLab;
                Some(format!(
                    "Pushed 2 integrals to NMR Mn: repeat={:.3}, end={:.3}",
                    repeat, end
                ))
            }
            SpectrumKind::Ftir => {
                let a = self.integration_regions[order[0]].area;
                let b = self.integration_regions[order[1]].area;
                if b <= 1e-12 {
                    return Some("FTIR ratio unavailable: second region area is zero".into());
                }
                self.polymer_lab.man_rxn_t = format!("{:.6}", a);
                self.polymer_lab.man_ref_t = format!("{:.6}", b);
                Some(format!("Area₁ / Area₂ = {:.6}", a / b))
            }
            _ => None,
        }
    }

    fn integration_table_text(&self) -> String {
        let kind = self.current_kind();
        let mut out = format!(
            "# SPECTRA Integration Regions — {}\n# Ref = region {} (×{:.2})\n\n",
            kind.display_name(),
            self.integration_ref_idx + 1,
            self.integration_ref_val
        );
        out.push_str("| # | Range | Center | Abs Area | Rel Integral | Ref | Assignment |\n");
        out.push_str("|--:|:------|-------:|---------:|-------------:|:---:|:-----------|\n");
        for (i, r) in self.integration_regions.iter().enumerate() {
            let rel = self
                .integration_relative(i)
                .map(|v| format!("{:.4}", v))
                .unwrap_or_else(|| "—".into());
            let mark = if i == self.integration_ref_idx { "★" } else { "" };
            out.push_str(&format!(
                "| {} | {:.3} – {:.3} {} | {:.3} | {:.5} | {} | {} | {} |\n",
                i + 1,
                r.x_lo.min(r.x_hi),
                r.x_lo.max(r.x_hi),
                kind.x_unit(),
                0.5 * (r.x_lo + r.x_hi),
                r.area,
                rel,
                mark,
                r.label
            ));
        }
        out
    }

    fn apply_calibration_target(&mut self, observed: f64, target: f64) -> String {
        let delta = target - observed;
        let kind = self.current_kind();
        if !delta.is_finite() || delta.abs() < 1e-12 {
            return "Target equals observed X — no calibration needed".into();
        }
        if let Some(spec) = self.traces.get_mut(self.active_trace) {
            apply_x_calibration(spec, delta);
        }
        self.cursor_x += delta;
        self.view_x_min += delta;
        self.view_x_max += delta;
        for p in &mut self.pinned_peaks {
            p.x += delta;
        }
        for r in &mut self.integration_regions {
            r.x_lo += delta;
            r.x_hi += delta;
        }
        self.deconv_overlay = None;
        self.recompute_processing();
        format!("Calibrated X-axis by {:+.3} {}", delta, kind.x_unit())
    }

    pub fn compute_live_conversion(&self) -> Option<ConversionResult> {
        if self.traces.len() < 2
            || self.traces[0].kind != SpectrumKind::Ftir
            || self.traces[1].kind != SpectrumKind::Ftir
        {
            return None;
        }
        let rxn_c: f64 = self.polymer_lab.rxn_center.parse().ok()?;
        let rxn_h: f64 = self.polymer_lab.rxn_half.parse().ok()?;
        let ref_c: f64 = self.polymer_lab.ref_center.parse().ok()?;
        let ref_h: f64 = self.polymer_lab.ref_half.parse().ok()?;
        auto_ftir_conversion(&self.traces[0], &self.traces[1], rxn_c, rxn_h, ref_c, ref_h)
    }

    pub fn compute_manual_conversion(&self) -> Option<f64> {
        ftir_conversion_percent(
            self.polymer_lab.man_rxn_0.parse().ok()?,
            self.polymer_lab.man_ref_0.parse().ok()?,
            self.polymer_lab.man_rxn_t.parse().ok()?,
            self.polymer_lab.man_ref_t.parse().ok()?,
        )
    }

    pub fn compute_nmr_mn(&self) -> Option<(f64, f64)> {
        nmr_mn_dpn(
            self.polymer_lab.m0.parse().ok()?,
            self.polymer_lab.m_end.parse().ok()?,
            self.polymer_lab.n_repeat.parse().ok()?,
            self.polymer_lab.n_end.parse().ok()?,
            self.polymer_lab.i_repeat.parse().ok()?,
            self.polymer_lab.i_end.parse().ok()?,
        )
    }

    pub fn compute_tacticity(&self) -> Option<TacticityResult> {
        tacticity_stats(
            self.polymer_lab.i_mm.parse().ok()?,
            self.polymer_lab.i_mr.parse().ok()?,
            self.polymer_lab.i_rr.parse().ok()?,
        )
    }

    /// Baseline-corrected maximum absorbance within `center ± half` on the active trace.
    fn pull_peak_absorbance(&self, center: f64, half: f64) -> Option<f64> {
        let spec = self.traces.get(self.active_trace)?;
        if spec.kind != SpectrumKind::Ftir {
            return None;
        }
        let cfg = ProcessConfig {
            y_pref: YAxisPreference::Swapped,
            baseline: BaselineMode::Asls,
            smooth: SmoothLevel::Off,
            ..ProcessConfig::default()
        };
        let (x, y, _) = process_trace(spec, cfg);
        let i0 = x.partition_point(|&v| v < center - half);
        let i1 = x.partition_point(|&v| v <= center + half);
        if i1 <= i0 {
            return None;
        }
        let peak = y[i0..i1].iter().copied().fold(f64::NEG_INFINITY, f64::max);
        peak.is_finite().then_some(peak.max(0.0))
    }

    /// Returns (A_840, A_763, F(β) %) using manual inputs or the active trace.
    pub fn compute_pvdf_beta(&self) -> Option<(f64, f64, f64)> {
        let (a840, a763) = if self.polymer_lab.phase_manual {
            (
                self.polymer_lab.a840.parse().ok()?,
                self.polymer_lab.a763.parse().ok()?,
            )
        } else {
            (
                self.pull_peak_absorbance(840.0, 15.0)?,
                self.pull_peak_absorbance(763.0, 15.0)?,
            )
        };
        let f_beta = pvdf_beta_fraction(a840, a763)?;
        Some((a840, a763, f_beta))
    }

    /// Returns the PU carbonyl hydrogen-bonding index using manual inputs or the active trace.
    pub fn compute_pu_hbi(&self) -> Option<HydrogenBondIndex> {
        let (a1703, a1732) = if self.polymer_lab.phase_manual {
            (
                self.polymer_lab.a1703.parse().ok()?,
                self.polymer_lab.a1732.parse().ok()?,
            )
        } else {
            (
                self.pull_peak_absorbance(1703.0, 15.0)?,
                self.pull_peak_absorbance(1732.0, 15.0)?,
            )
        };
        pu_hydrogen_bond_index(a1703, a1732)
    }

    /// Tauc band gap of the active UV-Vis trace in the selected transition mode.
    pub fn compute_tauc(&self) -> Option<TaucResult> {
        let spec = self.traces.get(self.active_trace)?;
        if spec.kind != SpectrumKind::UvVis {
            return None;
        }
        let mode = if self.polymer_lab.tauc_direct {
            TaucMode::Direct
        } else {
            TaucMode::Indirect
        };
        compute_tauc_band_gap_with_optics(spec, mode, self.process_cfg.optics)
    }

    /// Absorbance for Beer-Lambert: manual input, else the processed value at
    /// the cursor (auto-converted from %T when the trace is displayed as %T).
    pub fn compute_uv_absorbance(&self) -> Option<f64> {
        if self.polymer_lab.uv_manual {
            return self.polymer_lab.uv_abs.parse().ok();
        }
        let pt = self.processed.get(self.active_trace)?;
        let y = interpolate_y(&pt.x, &pt.y, self.cursor_x)?;
        if pt.y_unit == YUnit::Transmittance {
            Some(crate::processing::transmittance_to_absorbance(y))
        } else {
            Some(y)
        }
    }

    /// Beer-Lambert concentration in mol/L from the current A, ε and path length.
    pub fn compute_beer_lambert(&self) -> Option<f64> {
        let a = self.compute_uv_absorbance()?;
        let eps: f64 = self.polymer_lab.uv_eps.parse().ok()?;
        let path: f64 = self.polymer_lab.uv_path_cm.parse().ok()?;
        beer_lambert_conc(a, eps, path)
    }

    fn execute_export(&mut self) {
        let state = match &self.modal {
            Modal::Export(s) => s.clone(),
            _ => return,
        };

        let mut traces_out = Vec::new();
        let mut y_unit = YUnit::Transmittance;
        let mut exported_active = 0usize;
        for (idx, spec) in self.traces.iter().enumerate() {
            if !spec.visible {
                continue;
            }
            if let Some(pt) = self.processed.get(idx) {
                if traces_out.is_empty() {
                    y_unit = pt.y_unit;
                }
                if idx == self.active_trace {
                    exported_active = traces_out.len();
                }
                traces_out.push(ExportTrace {
                    name: spec.name.clone(),
                    x: pt.x.clone(),
                    y: pt.y.clone(),
                });
            }
        }
        if traces_out.is_empty() {
            if let Modal::Export(s) = &mut self.modal {
                s.feedback = Some("No visible traces to export!".into());
            }
            return;
        }

        let bounds = if state.use_zoom_range {
            [self.view_x_min, self.view_x_max]
        } else {
            let (lo, hi) = self.full_bounds();
            [lo, hi]
        };

        let mut peaks_out: Vec<ExportPeak> = Vec::new();
        let peak_min_dx = if self.current_kind() == SpectrumKind::UvVis {
            10.0
        } else {
            15.0
        };
        match state.peak_mode {
            ExportPeakMode::None => {}
            ExportPeakMode::PinnedOnly => {
                for p in &self.pinned_peaks {
                    peaks_out.push(ExportPeak {
                        x: p.x,
                        y: p.y,
                        label: p.label.clone(),
                        pinned: true,
                    });
                }
            }
            ExportPeakMode::PinnedAndAuto => {
                for p in &self.pinned_peaks {
                    peaks_out.push(ExportPeak {
                        x: p.x,
                        y: p.y,
                        label: p.label.clone(),
                        pinned: true,
                    });
                }
                for p in self.auto_peaks.iter().take(8) {
                    if peaks_out.iter().all(|q| (q.x - p.x).abs() > peak_min_dx) {
                        peaks_out.push(ExportPeak {
                            x: p.x,
                            y: p.y,
                            label: p.label.clone(),
                            pinned: false,
                        });
                    }
                }
            }
        }

        let mut out_path = PathBuf::from(&state.output_path);
        out_path.set_extension(state.format.ext());

        let payload = ExportPayload {
            mode: self.current_kind().as_str().to_string(),
            y_unit: y_unit.as_str().to_string(),
            style: state.style.as_str().to_string(),
            stacked: state.stacked,
            dpi: 600,
            title: state.title.clone(),
            x_bounds: bounds,
            active_trace: exported_active.min(traces_out.len().saturating_sub(1)),
            traces: traces_out,
            peaks: peaks_out,
        };

        self.pending_export = Some(state.clone());
        self.start_task(
            format!("Rendering {} …", state.format.label()),
            move || TaskOutcome::Export(run_export(&payload, &out_path).map_err(|e| e.to_string())),
        );
    }

    // ========================================================================
    // KEYBOARD EVENT DISPATCHER
    // ========================================================================
    pub fn handle_key(&mut self, key: KeyEvent) {
        if key.modifiers.contains(KeyModifiers::CONTROL) && key.code == KeyCode::Char('c') {
            self.should_quit = true;
            return;
        }

        if !matches!(self.modal, Modal::None) {
            self.handle_modal_key(key);
            return;
        }

        // Tab-specific input capture when editing fields
        if self.tab == ActiveTab::Solvents && self.solvents_tab.editing {
            self.handle_solvents_key(key);
            return;
        }
        if self.tab == ActiveTab::PolymerLab && self.polymer_lab.editing {
            self.handle_polymer_lab_key(key);
            return;
        }

        // Global shortcuts
        match key.code {
            KeyCode::Char('q') => {
                self.should_quit = true;
                return;
            }
            KeyCode::Char('?') => {
                self.help_scroll = 0;
                self.modal = Modal::Help;
                return;
            }
            KeyCode::Char('/') => {
                self.open_search_modal();
                return;
            }
            KeyCode::Char('1') => {
                self.switch_tab(ActiveTab::Ftir);
                return;
            }
            KeyCode::Char('2') => {
                self.switch_tab(ActiveTab::Hnmr);
                return;
            }
            KeyCode::Char('3') => {
                self.switch_tab(ActiveTab::Cnmr);
                return;
            }
            KeyCode::Char('4') => {
                self.switch_tab(ActiveTab::UvVis);
                return;
            }
            KeyCode::Char('5') => {
                self.switch_tab(ActiveTab::Solvents);
                return;
            }
            KeyCode::Char('6') => {
                self.switch_tab(ActiveTab::PolymerLab);
                return;
            }
            KeyCode::Char('f') => {
                self.open_file_picker(false);
                return;
            }
            KeyCode::Char('F') => {
                self.open_file_picker(true);
                return;
            }
            KeyCode::Char('D') => {
                self.modal = Modal::LibrarySetup(LibrarySetupState::default());
                return;
            }
            // Also allow Shift+d (explicit) as some terminals report it differently
            KeyCode::Char('d') if key.modifiers.contains(KeyModifiers::SHIFT) => {
                self.modal = Modal::LibrarySetup(LibrarySetupState::default());
                return;
            }
            KeyCode::Char('d') => {
                self.modal = Modal::DemoPicker(0);
                return;
            }
            KeyCode::Char('e') => {
                self.open_export_modal();
                return;
            }
            KeyCode::Char('M') => {
                self.open_library_match();
                return;
            }
            // CapsLock + Shift is often reported as lowercase 'm' with SHIFT set.
            KeyCode::Char('m') if key.modifiers.contains(KeyModifiers::SHIFT) => {
                self.open_library_match();
                return;
            }
            _ => {}
        }

        match self.tab {
            ActiveTab::Solvents => self.handle_solvents_key(key),
            ActiveTab::PolymerLab => self.handle_polymer_lab_key(key),
            ActiveTab::Ftir | ActiveTab::Hnmr | ActiveTab::Cnmr | ActiveTab::UvVis => {
                self.handle_spectrum_key(key)
            }
        }
    }

    /// Remove the active trace from the workspace; safe all the way down to 0.
    pub fn remove_active_trace(&mut self) {
        if self.traces.is_empty() {
            self.status_msg = "No traces to remove".into();
            return;
        }
        let idx = self.active_trace.min(self.traces.len() - 1);
        let removed = self.traces.remove(idx);
        self.active_trace = idx.min(self.traces.len().saturating_sub(1));
        self.deconv_overlay = None;

        if self.traces.is_empty() {
            self.pinned_peaks.clear();
            self.auto_peaks.clear();
            self.processed.clear();
            self.selected_peak = 0;
            self.clear_integrations();
            self.reset_zoom();
        } else {
            self.recompute_processing();
            self.reset_zoom();
        }
        self.status_msg = format!("Removed trace: {}", removed.name);
    }

    fn handle_spectrum_key(&mut self, key: KeyEvent) {
        let shift = key.modifiers.contains(KeyModifiers::SHIFT);
        match key.code {            KeyCode::Left | KeyCode::Char('h') => self.move_cursor(-1.0, shift),
            KeyCode::Right | KeyCode::Char('l') => self.move_cursor(1.0, shift),
            KeyCode::Char('H') => self.move_cursor(-1.0, true),
            KeyCode::Char('L') => self.move_cursor(1.0, true),
            KeyCode::Char('g') => self.modal = Modal::GoTo(String::new()),
            KeyCode::Char('+') | KeyCode::Char('=') | KeyCode::Char('z') => {
                self.zoom_around_cursor(0.65);
            }
            KeyCode::Char('-') | KeyCode::Char('Z') => {
                self.zoom_around_cursor(1.55);
            }
            KeyCode::Char('0') | KeyCode::Char('r') => {
                self.reset_zoom();
                self.status_msg = "Reset zoom to full spectrum".into();
            }
            KeyCode::Char('A') | KeyCode::Char('a') if shift || key.code == KeyCode::Char('A') => {
                match self.current_kind() {
                    SpectrumKind::Ftir => {
                        self.process_cfg.optics = self.process_cfg.optics.next();
                        self.recompute_processing();
                        self.status_msg = format!(
                            "Optical Physics Correction: {}",
                            self.process_cfg.optics.label()
                        );
                    }
                    SpectrumKind::UvVis => {
                        self.process_cfg.optics = if self.process_cfg.optics == OpticsMode::KubelkaMunk {
                            OpticsMode::Off
                        } else {
                            OpticsMode::KubelkaMunk
                        };
                        self.recompute_processing();
                        self.status_msg = format!(
                            "Optical Physics Correction: {}",
                            self.process_cfg.optics.label()
                        );
                    }
                    SpectrumKind::Hnmr | SpectrumKind::Cnmr => {
                        self.status_msg = "ATR / Kubelka-Munk optics apply to FTIR and UV-Vis".into();
                    }
                }
            }
            KeyCode::Char('a') => {
                if matches!(
                    self.current_kind(),
                    SpectrumKind::Ftir | SpectrumKind::UvVis
                ) {
                    self.process_cfg.y_pref = match self.process_cfg.y_pref {
                        YAxisPreference::Native => YAxisPreference::Swapped,
                        YAxisPreference::Swapped => YAxisPreference::Native,
                    };
                    self.recompute_processing();
                    self.status_msg = match (self.current_kind(), self.process_cfg.y_pref) {
                        (SpectrumKind::Ftir, YAxisPreference::Swapped) => {
                            "FTIR Y-Axis: Absorbance (A)".into()
                        }
                        (SpectrumKind::Ftir, YAxisPreference::Native) => {
                            "FTIR Y-Axis: Transmittance (%T)".into()
                        }
                        (_, YAxisPreference::Swapped) => {
                            "UV-Vis Y-Axis: Transmittance (%T)".into()
                        }
                        (_, YAxisPreference::Native) => {
                            "UV-Vis Y-Axis: Absorbance (A)".into()
                        }
                    };
                }
            }
            KeyCode::Char('b') => {
                self.process_cfg.baseline = self.process_cfg.baseline.next();
                self.recompute_processing();
                self.status_msg = format!(
                    "Baseline Correction: {}",
                    self.process_cfg.baseline.label()
                );
            }
            KeyCode::Char('m') => {
                self.process_cfg.smooth = self.process_cfg.smooth.next();
                self.recompute_processing();
                self.status_msg = format!("Savitzky-Golay Smoothing: {}", self.process_cfg.smooth.label());
            }
            KeyCode::Char('w') => {
                self.process_cfg.norm = self.process_cfg.norm.next();
                self.recompute_processing();
                self.status_msg = format!("Normalization: {}", self.process_cfg.norm.label());
            }
            KeyCode::Char('t') => {
                self.process_cfg.deriv = self.process_cfg.deriv.next();
                self.recompute_processing();
                self.status_msg = format!("Derivative View: {}", self.process_cfg.deriv.label());
            }
            KeyCode::Char('s') => {
                self.stacked = !self.stacked;
                self.status_msg = format!("Layout: {}", if self.stacked { "Stacked Offset" } else { "Overlay" });
            }
            KeyCode::Char('v') => {
                self.hd_mode = !self.hd_mode;
                self.status_msg = format!(
                    "Display Mode: {}",
                    if self.hd_mode {
                        "HD Oscilloscope (Kitty)"
                    } else {
                        "Braille (Classic)"
                    }
                );
            }
            KeyCode::Char('P') | KeyCode::Char('p') if shift || key.code == KeyCode::Char('P') => {
                if self.traces.get(self.active_trace).is_some() {
                    let prev_res = self
                        .deconv_overlay
                        .as_ref()
                        .filter(|ov| ov.trace_idx == self.active_trace)
                        .map(|ov| ov.result.clone());
                    let shape = prev_res
                        .as_ref()
                        .map(|r| r.shape)
                        .unwrap_or(PeakShape::PseudoVoigt);
                    let max_peaks = prev_res
                        .as_ref()
                        .map(|r| r.peaks.len().clamp(1, 8))
                        .unwrap_or(3);
                    self.modal = Modal::Deconv(DeconvModalState {
                        field_idx: 0,
                        shape,
                        max_peaks,
                        fit_baseline: true,
                        x_lo: format!("{:.2}", self.view_x_min),
                        x_hi: format!("{:.2}", self.view_x_max),
                        result: prev_res,
                        feedback: None,
                    });
                } else {
                    self.status_msg = "Load a spectrum before opening Peak Deconvolution Studio".into();
                }
            }
            KeyCode::Char('p') => {
                self.show_peaks = !self.show_peaks;
                self.status_msg = format!("Auto Peak Markers: {}", if self.show_peaks { "ON" } else { "OFF" });
            }
            KeyCode::Char('[') => {
                self.peak_threshold = (self.peak_threshold - 0.02).clamp(0.02, 0.60);
                self.recompute_peaks();
                self.status_msg = format!("Peak threshold: {:.0}% ({} peaks)", self.peak_threshold * 100.0, self.auto_peaks.len());
            }
            KeyCode::Char(']') => {
                self.peak_threshold = (self.peak_threshold + 0.02).clamp(0.02, 0.60);
                self.recompute_peaks();
                self.status_msg = format!("Peak threshold: {:.0}% ({} peaks)", self.peak_threshold * 100.0, self.auto_peaks.len());
            }
            KeyCode::Char('n') => {
                if !self.auto_peaks.is_empty() {
                    self.selected_peak = (self.selected_peak + 1) % self.auto_peaks.len();
                    let px = self.auto_peaks[self.selected_peak].x;
                    self.jump_to_x(px);
                }
            }
            KeyCode::Char('N') => {
                if !self.auto_peaks.is_empty() {
                    self.selected_peak = (self.selected_peak + self.auto_peaks.len() - 1) % self.auto_peaks.len();
                    let px = self.auto_peaks[self.selected_peak].x;
                    self.jump_to_x(px);
                }
            }
            KeyCode::Char('X') | KeyCode::Char('x') if shift || key.code == KeyCode::Char('X') => {
                if self.traces.get(self.active_trace).is_some() {
                    let observed = self.cursor_x;
                    self.modal = Modal::Calibration(CalibrationModalState {
                        observed_x: observed,
                        target_buf: format!("{:.4}", observed),
                        feedback: None,
                    });
                    self.status_msg = format!(
                        "X-Axis Calibration — observed {:.3} {} · type target or use [t]/[s]/[w]",
                        observed,
                        self.current_kind().x_unit()
                    );
                } else {
                    self.status_msg = "Load a spectrum before calibrating the X-axis".into();
                }
            }
            KeyCode::Char('x') => self.toggle_pin_at_cursor(),
            KeyCode::Char('C') | KeyCode::Char('c') if shift || key.code == KeyCode::Char('C') => {
                if self.current_kind() == SpectrumKind::Ftir {
                    self.process_cfg.atm_co2 = !self.process_cfg.atm_co2;
                    self.recompute_processing();
                    self.status_msg = format!(
                        "Atmospheric CO₂ Suppression (2280–2410 cm⁻¹): {}",
                        if self.process_cfg.atm_co2 { "ON" } else { "OFF" }
                    );
                } else {
                    self.status_msg = "CO₂ suppression applies to FTIR spectra".into();
                }
            }
            KeyCode::Char('c') => self.copy_cursor_summary(),
            KeyCode::Char('o') => {
                if self.traces.len() >= 2 {
                    let ta = self.active_trace.min(self.traces.len() - 1);
                    let tb = (ta + 1) % self.traces.len();
                    let k_opt = optimal_subtraction_factor(
                        &self.traces[ta],
                        &self.traces[tb],
                        self.view_x_min,
                        self.view_x_max,
                    )
                    .unwrap_or(1.0);
                    self.modal = Modal::SpectralMath(SpectralMathState {
                        field_idx: 0,
                        op: SpectralOp::Subtract,
                        trace_a: ta,
                        trace_b: tb,
                        factor_k: k_opt,
                        feedback: None,
                    });
                } else {
                    self.status_msg = "Spectral Math requires at least 2 traces in workspace".into();
                }
            }
            KeyCode::Char('S') => {
                self.nmr_solvent = self.nmr_solvent.next();
                self.solvents_tab.solvent = self.nmr_solvent;
                self.recompute_peaks();
                self.status_msg = format!("NMR Reference Solvent: {}", self.nmr_solvent.label());
            }
            KeyCode::Char('I') | KeyCode::Char('i') if shift || key.code == KeyCode::Char('I') => {
                if self.traces.get(self.active_trace).is_some() {
                    let selected = self.integration_regions.len().saturating_sub(1);
                    self.modal = Modal::Integration(IntegrationModalState {
                        selected,
                        feedback: None,
                    });
                    self.status_msg =
                        "Quantitative Integration Studio — [z] zoom window · [a] cursor · [j/k] navigate"
                            .into();
                } else {
                    self.status_msg = "Load a spectrum before opening the Integration Studio".into();
                }
            }
            KeyCode::Char('i') => {
                self.show_solvent_markers = !self.show_solvent_markers;
                self.status_msg = format!(
                    "NMR Solvent Markers: {}",
                    if self.show_solvent_markers { "ON" } else { "OFF" }
                );
            }
            KeyCode::Tab => {
                if !self.traces.is_empty() {
                    self.active_trace = (self.active_trace + 1) % self.traces.len();
                    self.recompute_peaks();
                    self.status_msg = format!("Active trace: {}", self.traces[self.active_trace].name);
                }
            }
            KeyCode::Char(' ') => {
                if let Some(t) = self.traces.get_mut(self.active_trace) {
                    t.visible = !t.visible;
                    self.status_msg = format!(
                        "{}: {}",
                        t.name,
                        if t.visible { "Visible" } else { "Hidden" }
                    );
                }
            }
            KeyCode::Delete | KeyCode::Backspace => self.remove_active_trace(),
            _ => {}
        }
    }

    fn handle_solvents_key(&mut self, key: KeyEvent) {
        if self.solvents_tab.editing {
            match key.code {
                KeyCode::Esc | KeyCode::Enter => self.solvents_tab.editing = false,
                KeyCode::Backspace => {
                    self.solvents_tab.query.pop();
                    self.solvents_tab.selected = 0;
                }
                KeyCode::Char(c) => {
                    self.solvents_tab.query.push(c);
                    self.solvents_tab.selected = 0;
                }
                _ => {}
            }
            return;
        }

        match key.code {
            KeyCode::Left | KeyCode::Char('h') => {
                self.solvents_tab.solvent = self.solvents_tab.solvent.prev();
                self.nmr_solvent = self.solvents_tab.solvent;
                self.recompute_peaks();
            }
            KeyCode::Right | KeyCode::Char('l') => {
                self.solvents_tab.solvent = self.solvents_tab.solvent.next();
                self.nmr_solvent = self.solvents_tab.solvent;
                self.recompute_peaks();
            }
            KeyCode::Up | KeyCode::Char('k') => {
                self.solvents_tab.selected = self.solvents_tab.selected.saturating_sub(1);
            }
            KeyCode::Down | KeyCode::Char('j') => {
                self.solvents_tab.selected = self.solvents_tab.selected.saturating_add(1);
            }
            KeyCode::Char('n') | KeyCode::Tab => {
                self.solvents_tab.is_carbon = !self.solvents_tab.is_carbon;
                self.solvents_tab.selected = 0;
            }
            KeyCode::Char('t') => {
                self.solvents_tab.tol_idx = (self.solvents_tab.tol_idx + 1) % 4;
            }
            KeyCode::Char('s') | KeyCode::Enter => {
                self.solvents_tab.editing = true;
            }
            KeyCode::Char('x') => {
                self.solvents_tab.query.clear();
                self.solvents_tab.selected = 0;
            }
            _ => {}
        }
    }

    fn active_polymer_field_mut(&mut self) -> Option<&mut String> {
        let p = &mut self.polymer_lab;
        match p.sub_tab {
            PolymerSubTab::Conversion => {
                if p.manual_mode {
                    match p.field_idx {
                        0 => Some(&mut p.man_rxn_0),
                        1 => Some(&mut p.man_ref_0),
                        2 => Some(&mut p.man_rxn_t),
                        3 => Some(&mut p.man_ref_t),
                        _ => None,
                    }
                } else {
                    match p.field_idx {
                        0 => Some(&mut p.rxn_center),
                        1 => Some(&mut p.rxn_half),
                        2 => Some(&mut p.ref_center),
                        3 => Some(&mut p.ref_half),
                        _ => None,
                    }
                }
            }
            PolymerSubTab::NmrMn => match p.field_idx {
                0 => Some(&mut p.m0),
                1 => Some(&mut p.m_end),
                2 => Some(&mut p.n_repeat),
                3 => Some(&mut p.n_end),
                4 => Some(&mut p.i_repeat),
                5 => Some(&mut p.i_end),
                _ => None,
            },
            PolymerSubTab::Tacticity => match p.field_idx {
                0 => Some(&mut p.i_mm),
                1 => Some(&mut p.i_mr),
                2 => Some(&mut p.i_rr),
                _ => None,
            },
            PolymerSubTab::PhaseIndex => match p.field_idx {
                0 => Some(&mut p.a840),
                1 => Some(&mut p.a763),
                2 => Some(&mut p.a1703),
                3 => Some(&mut p.a1732),
                _ => None,
            },
            PolymerSubTab::UvVisLab => match p.field_idx {
                0 => Some(&mut p.uv_abs),
                1 => Some(&mut p.uv_eps),
                2 => Some(&mut p.uv_path_cm),
                _ => None,
            },
        }
    }

    fn polymer_field_count(&self) -> usize {
        match self.polymer_lab.sub_tab {
            PolymerSubTab::Conversion => 4,
            PolymerSubTab::NmrMn => 6,
            PolymerSubTab::Tacticity => 3,
            PolymerSubTab::PhaseIndex => 4,
            PolymerSubTab::UvVisLab => 3,
        }
    }

    fn handle_polymer_lab_key(&mut self, key: KeyEvent) {
        if self.polymer_lab.editing {
            match key.code {
                KeyCode::Esc | KeyCode::Enter => self.polymer_lab.editing = false,
                KeyCode::Backspace => {
                    if let Some(f) = self.active_polymer_field_mut() {
                        f.pop();
                    }
                }
                KeyCode::Char(c) if c.is_ascii_digit() || c == '.' || c == '-' => {
                    if let Some(f) = self.active_polymer_field_mut() {
                        f.push(c);
                    }
                }
                _ => {}
            }
            return;
        }

        let cnt = self.polymer_field_count();
        match key.code {
            KeyCode::Tab | KeyCode::Right | KeyCode::Char('l') => {
                self.polymer_lab.sub_tab = match self.polymer_lab.sub_tab {
                    PolymerSubTab::Conversion => PolymerSubTab::NmrMn,
                    PolymerSubTab::NmrMn => PolymerSubTab::Tacticity,
                    PolymerSubTab::Tacticity => PolymerSubTab::PhaseIndex,
                    PolymerSubTab::PhaseIndex => PolymerSubTab::UvVisLab,
                    PolymerSubTab::UvVisLab => PolymerSubTab::Conversion,
                };
                self.polymer_lab.field_idx = 0;
            }
            KeyCode::BackTab | KeyCode::Left | KeyCode::Char('h') => {
                self.polymer_lab.sub_tab = match self.polymer_lab.sub_tab {
                    PolymerSubTab::Conversion => PolymerSubTab::UvVisLab,
                    PolymerSubTab::NmrMn => PolymerSubTab::Conversion,
                    PolymerSubTab::Tacticity => PolymerSubTab::NmrMn,
                    PolymerSubTab::PhaseIndex => PolymerSubTab::Tacticity,
                    PolymerSubTab::UvVisLab => PolymerSubTab::PhaseIndex,
                };
                self.polymer_lab.field_idx = 0;
            }
            KeyCode::Down | KeyCode::Char('j') => {
                self.polymer_lab.field_idx = (self.polymer_lab.field_idx + 1) % cnt;
            }
            KeyCode::Up | KeyCode::Char('k') => {
                self.polymer_lab.field_idx = (self.polymer_lab.field_idx + cnt - 1) % cnt;
            }
            KeyCode::Enter => {
                self.polymer_lab.editing = true;
            }
            KeyCode::Char('m') if self.polymer_lab.sub_tab == PolymerSubTab::Conversion => {
                self.polymer_lab.manual_mode = !self.polymer_lab.manual_mode;
                self.polymer_lab.field_idx = 0;
            }
            KeyCode::Char('m') if self.polymer_lab.sub_tab == PolymerSubTab::PhaseIndex => {
                self.polymer_lab.phase_manual = !self.polymer_lab.phase_manual;
                self.polymer_lab.field_idx = 0;
            }
            KeyCode::Char('m') if self.polymer_lab.sub_tab == PolymerSubTab::UvVisLab => {
                self.polymer_lab.uv_manual = !self.polymer_lab.uv_manual;
                self.polymer_lab.field_idx = 0;
                self.status_msg = format!(
                    "UV-Vis absorbance: {}",
                    if self.polymer_lab.uv_manual {
                        "Manual entry"
                    } else {
                        "Auto (A at cursor)"
                    }
                );
            }
            KeyCode::Char('t') if self.polymer_lab.sub_tab == PolymerSubTab::UvVisLab => {
                self.polymer_lab.tauc_direct = !self.polymer_lab.tauc_direct;
                self.status_msg = format!(
                    "Tauc transition: {}",
                    if self.polymer_lab.tauc_direct {
                        "Direct allowed — (αhν)²"
                    } else {
                        "Indirect allowed — (αhν)^½"
                    }
                );
            }
            KeyCode::Char('p') if self.polymer_lab.sub_tab == PolymerSubTab::NmrMn => {
                let next = (self.polymer_lab.monomer_preset + 1) % MONOMER_PRESETS.len();
                self.polymer_lab.monomer_preset = next;
                let (name, m0, nr) = MONOMER_PRESETS[next];
                self.polymer_lab.m0 = format!("{:.2}", m0);
                self.polymer_lab.n_repeat = nr.to_string();
                self.status_msg = format!("Monomer preset: {} (M₀ = {:.2} g/mol)", name, m0);
            }
            KeyCode::Char('c') => {
                let summary = match self.polymer_lab.sub_tab {
                    PolymerSubTab::Conversion => {
                        if !self.polymer_lab.manual_mode {
                            if let Some(res) = self.compute_live_conversion() {
                                format!(
                                    "FTIR Conversion: {:.1}% (area) / {:.1}% (height)",
                                    res.area_conversion_pct, res.height_conversion_pct
                                )
                            } else {
                                "Load 2 FTIR traces (monomer & polymer) first".into()
                            }
                        } else if let Some(c) = self.compute_manual_conversion() {
                            format!("FTIR Conversion: {:.1}%", c)
                        } else {
                            "Invalid conversion inputs".into()
                        }
                    }
                    PolymerSubTab::NmrMn => {
                        if let Some((dpn, mn)) = self.compute_nmr_mn() {
                            format!("¹H NMR End-Group Analysis: DP_n = {:.1}, M_n = {:.0} g/mol ({:.2} kDa)", dpn, mn, mn / 1000.0)
                        } else {
                            "Invalid Mn inputs".into()
                        }
                    }
                    PolymerSubTab::Tacticity => {
                        if let Some(t) = self.compute_tacticity() {
                            format!(
                                "Tacticity: mm={:.1}%, mr={:.1}%, rr={:.1}% | m={:.3}, r={:.3}, ρ={:.2}",
                                t.pct_mm, t.pct_mr, t.pct_rr, t.dyad_m, t.dyad_r, t.rho
                            )
                        } else {
                            "Invalid tacticity inputs".into()
                        }
                    }
                    PolymerSubTab::PhaseIndex => match (self.compute_pvdf_beta(), self.compute_pu_hbi()) {
                        (Some((a840, a763, f)), Some(h)) => format!(
                            "PVDF F(β) = {:.1}% (A840={:.3}, A763={:.3}) | PU HBI = {:.3}, bonded = {:.1}%",
                            f, a840, a763, h.hbi, h.bonded_pct
                        ),
                        (Some((a840, a763, f)), None) => {
                            format!("PVDF F(β) = {:.1}% (A840={:.3}, A763={:.3})", f, a840, a763)
                        }
                        (None, Some(h)) => {
                            format!("PU HBI = {:.3}, bonded C=O = {:.1}%", h.hbi, h.bonded_pct)
                        }
                        (None, None) => "Invalid phase-index inputs".into(),
                    },
                    PolymerSubTab::UvVisLab => {
                        let mode = if self.polymer_lab.tauc_direct {
                            "direct"
                        } else {
                            "indirect"
                        };
                        let tauc = match self.compute_tauc() {
                            Some(t) => format!(
                                "Tauc E_g = {:.3} eV ({:.1} nm, {} , R² = {:.4})",
                                t.band_gap_ev, t.edge_nm, mode, t.r_squared
                            ),
                            None => "Tauc: load a UV-Vis spectrum on the active trace".into(),
                        };
                        let conc = match self.compute_beer_lambert() {
                            Some(c) => {
                                let a = self.compute_uv_absorbance().unwrap_or(0.0);
                                format!(
                                    "Beer-Lambert: c = {:.4e} mol/L (A = {:.3}, ε = {}, b = {} cm)",
                                    c, a, self.polymer_lab.uv_eps, self.polymer_lab.uv_path_cm
                                )
                            }
                            None => "Beer-Lambert: invalid A / ε / path".into(),
                        };
                        format!("{} | {}", tauc, conc)
                    }
                };
                match copy_text(&summary) {
                    Ok(tool) => self.status_msg = format!("Copied via {}: {}", tool, summary),
                    Err(e) => self.status_msg = format!("Copy failed: {}", e),
                }
            }
            _ => {}
        }
    }

    fn handle_modal_key(&mut self, key: KeyEvent) {
        if key.code == KeyCode::Esc {
            // File picker: Esc backs out one layer at a time instead of
            // discarding the whole browser (path entry → live search → close).
            if let Modal::FilePicker(st) = &mut self.modal {
                if st.path_input.is_some() {
                    st.path_input = None;
                    return;
                }
                if !st.filter.is_empty() {
                    st.filter.clear();
                    st.change_filter();
                    return;
                }
            }
            // Dismissing the ambiguous-query picker should explain what to do
            // next rather than leave a stale "N candidates" status behind.
            if matches!(self.modal, Modal::CandidatePicker(_)) {
                self.status_msg =
                    "Candidate picker dismissed — [f] opens files, [?] shows shortcuts".into();
            } else if matches!(self.modal, Modal::Busy(_)) {
                self.status_msg = "Working in the background — results appear when ready".into();
            }
            self.modal = Modal::None;
            return;
        }

        let current = std::mem::replace(&mut self.modal, Modal::None);
        match current {
            Modal::Help => {
                match key.code {
                    KeyCode::Down | KeyCode::Char('j') => {
                        self.help_scroll = self.help_scroll.saturating_add(1);
                    }
                    KeyCode::Up | KeyCode::Char('k') => {
                        self.help_scroll = self.help_scroll.saturating_sub(1);
                    }
                    KeyCode::PageDown | KeyCode::Char('d') | KeyCode::Char(' ') => {
                        self.help_scroll = self.help_scroll.saturating_add(12);
                    }
                    KeyCode::PageUp | KeyCode::Char('u') => {
                        self.help_scroll = self.help_scroll.saturating_sub(12);
                    }
                    KeyCode::Left | KeyCode::Right => {}
                    KeyCode::Home => self.help_scroll = 0,
                    KeyCode::End => self.help_scroll = usize::MAX / 2,
                    KeyCode::Char('q') | KeyCode::Char('?') | KeyCode::Enter => {
                        self.modal = Modal::None;
                        return;
                    }
                    _ => {}
                }
                self.modal = Modal::Help;
            }
            Modal::GoTo(mut buf) => match key.code {
                KeyCode::Enter => {
                    if let Ok(val) = buf.trim().parse::<f64>() {
                        self.jump_to_x(val);
                        self.status_msg = format!(
                            "Jumped to {:.2} {}",
                            self.cursor_x,
                            self.current_kind().x_unit()
                        );
                    }
                    self.modal = Modal::None;
                }
                KeyCode::Backspace => {
                    buf.pop();
                    self.modal = Modal::GoTo(buf);
                }
                KeyCode::Char(c) if c.is_ascii_digit() || c == '.' || c == '-' => {
                    buf.push(c);
                    self.modal = Modal::GoTo(buf);
                }
                _ => self.modal = Modal::GoTo(buf),
            },
            Modal::DemoPicker(mut idx) => match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    idx = (idx + DemoPreset::ALL.len() - 1) % DemoPreset::ALL.len();
                    self.modal = Modal::DemoPicker(idx);
                }
                KeyCode::Down | KeyCode::Char('j') => {
                    idx = (idx + 1) % DemoPreset::ALL.len();
                    self.modal = Modal::DemoPicker(idx);
                }
                KeyCode::Enter => {
                    self.load_demo(DemoPreset::ALL[idx]);
                    self.modal = Modal::None;
                }
                _ => self.modal = Modal::DemoPicker(idx),
            },
            Modal::CandidatePicker(mut st) => match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    if !st.entries.is_empty() {
                        st.selected = (st.selected + st.entries.len() - 1) % st.entries.len();
                    }
                    self.modal = Modal::CandidatePicker(st);
                }
                KeyCode::Down | KeyCode::Char('j') => {
                    if !st.entries.is_empty() {
                        st.selected = (st.selected + 1) % st.entries.len();
                    }
                    self.modal = Modal::CandidatePicker(st);
                }
                KeyCode::Enter => {
                    if let Some(entry) = st.entries.get(st.selected).cloned() {
                        let label = format!(
                            "Opening {} ({}) …",
                            entry.stem.replace('_', " "),
                            entry.kind.display_name()
                        );
                        self.load_paths(vec![entry.path], Some(entry.kind), false, label);
                    }
                }
                KeyCode::Esc | KeyCode::Char('q') => {
                    self.modal = Modal::None;
                }
                _ => self.modal = Modal::CandidatePicker(st),
            },
            Modal::LibraryMatch(mut st) => match key.code {
                KeyCode::Char('m') => {
                    st.mixture_mode = !st.mixture_mode;
                    if st.mixture_mode && st.mixture_results.is_empty() {
                        if let Some(spec) = self.traces.get(self.active_trace).cloned() {
                            let cfg = self.process_cfg;
                            let solvent = self.nmr_solvent;
                            self.pending_match_state = Some(st);
                            self.start_task(
                                "Analyzing 2–3 component mixtures (NNLS)…",
                                move || {
                                    TaskOutcome::MixtureMatch(crate::mixture::analyze_mixture(
                                        &spec, cfg, solvent, 3,
                                    ))
                                },
                            );
                        } else {
                            self.status_msg =
                                "Load a spectrum before mixture analysis".into();
                            self.modal = Modal::LibraryMatch(st);
                        }
                    } else {
                        self.modal = Modal::LibraryMatch(st);
                    }
                }
                KeyCode::Up | KeyCode::Char('k') => {
                    if st.mixture_mode {
                        if !st.mixture_results.is_empty() {
                            st.mixture_selected = (st.mixture_selected
                                + st.mixture_results.len()
                                - 1)
                                % st.mixture_results.len();
                        }
                    } else if !st.results.is_empty() {
                        st.selected = (st.selected + st.results.len() - 1) % st.results.len();
                    }
                    self.modal = Modal::LibraryMatch(st);
                }
                KeyCode::Down | KeyCode::Char('j') => {
                    if st.mixture_mode {
                        if !st.mixture_results.is_empty() {
                            st.mixture_selected =
                                (st.mixture_selected + 1) % st.mixture_results.len();
                        }
                    } else if !st.results.is_empty() {
                        st.selected = (st.selected + 1) % st.results.len();
                    }
                    self.modal = Modal::LibraryMatch(st);
                }
                KeyCode::Char('/') => self.open_search_modal(),
                KeyCode::Char('c') if st.mixture_mode => {
                    if st.mixture_results.is_empty() {
                        self.status_msg = "No mixture report to copy yet".into();
                    } else {
                        let report = crate::mixture::format_mixture_markdown(
                            &st.mixture_results,
                            self.current_kind(),
                        );
                        match copy_text(&report) {
                            Ok(tool) => {
                                self.status_msg =
                                    format!("Copied mixture report via {}", tool);
                            }
                            Err(e) => self.status_msg = format!("Copy failed: {}", e),
                        }
                    }
                    self.modal = Modal::LibraryMatch(st);
                }
                KeyCode::Enter => {
                    if st.mixture_mode {
                        if let Some(m) = st.mixture_results.get(st.mixture_selected).cloned() {
                            self.overlay_mixture_match(&m);
                        }
                    } else if let Some(m) = st.results.get(st.selected).cloned() {
                        self.apply_library_match(&m);
                    }
                    self.modal = Modal::None;
                }
                _ => self.modal = Modal::LibraryMatch(st),
            },
            Modal::Search(mut st) => match key.code {
                KeyCode::Down => {
                    if !st.results.is_empty() {
                        st.selected = (st.selected + 1) % st.results.len();
                    }
                    self.modal = Modal::Search(st);
                }
                KeyCode::Up => {
                    if !st.results.is_empty() {
                        st.selected = (st.selected + st.results.len() - 1) % st.results.len();
                    }
                    self.modal = Modal::Search(st);
                }
                KeyCode::Enter => {
                    if let Some(item) = st.results.get(st.selected).cloned() {
                        if let Some(path) = item.file_path.clone() {
                            // Reference spectrum from `library/`: open it directly.
                            let label = format!("Opening {} …", item.title);
                            self.load_paths(vec![path], None, false, label);
                            return;
                        }
                        if let Some(short) = item.polymer_short {
                            if let Some(poly) = find_polymer_by_short(short) {
                                let mut spec =
                                    synthesize_polymer_reference(poly, self.current_kind());
                                let kind = spec.kind;
                                let has_external =
                                    self.traces.iter().any(|t| t.source_path.is_some());

                                if has_external
                                    && self.traces.iter().all(|t| t.kind == kind)
                                {
                                    // Keep the user's experimental file and overlay the
                                    // built-in reference spectrum on top of it.
                                    self.traces.retain(|t| !t.name.ends_with("[Ref]"));
                                    spec.color_idx = self.traces.len();
                                    self.traces.push(spec);
                                    self.active_trace = self.traces.len() - 1;
                                    self.recompute_processing();
                                    self.reset_zoom();
                                    self.status_msg = format!(
                                        "Overlaid {} reference ({}, {})",
                                        poly.short, poly.repeat_unit, poly.tg_tm
                                    );
                                } else {
                                    // No matching lab file: draw the polymer's own
                                    // reference spectrum instead of leaving a stale one.
                                    spec.color_idx = 0;
                                    self.tab = ActiveTab::from_kind(kind);
                                    self.traces = vec![spec];
                                    self.active_trace = 0;
                                    self.pinned_peaks.clear();
                                    self.recompute_processing();
                                    self.reset_zoom();
                                    let jx = match kind {
                                        SpectrumKind::Ftir => item
                                            .jump_x
                                            .or_else(|| poly.ftir_peaks.last().map(|(x, _)| *x)),
                                        SpectrumKind::Hnmr => {
                                            poly.hnmr_peaks.first().map(|(x, _)| *x)
                                        }
                                        SpectrumKind::Cnmr => {
                                            poly.cnmr_peaks.first().map(|(x, _)| *x)
                                        }
                                        SpectrumKind::UvVis => item.jump_x,
                                    };
                                    if let Some(jx) = jx {
                                        self.jump_to_x(jx);
                                    }
                                    self.status_msg = format!(
                                        "Loaded {} reference: {} | {} | {}",
                                        poly.short, poly.name, poly.repeat_unit, poly.tg_tm
                                    );
                                }
                                self.modal = Modal::None;
                                return;
                            }
                        }

                        self.switch_tab(item.target_tab);
                        if let Some(jx) = item.jump_x {
                            self.jump_to_x(jx);
                        }
                        self.status_msg = format!("Selected: {}", item.title);
                    }
                    self.modal = Modal::None;
                }
                KeyCode::Backspace => {
                    st.query.pop();
                    st.results = build_search_results(&st.query);
                    st.selected = 0;
                    self.modal = Modal::Search(st);
                }
                KeyCode::Char(c) => {
                    st.query.push(c);
                    st.results = build_search_results(&st.query);
                    st.selected = 0;
                    self.modal = Modal::Search(st);
                }
                _ => self.modal = Modal::Search(st),
            },
            Modal::FilePicker(mut st) => {
                // While a path is being typed/pasted, the picker behaves like a
                // text field so files outside the current directory are reachable.
                if st.path_input.is_some() {
                    match key.code {
                        KeyCode::Char(c)
                            if !key.modifiers.contains(KeyModifiers::CONTROL) =>
                        {
                            if let Some(buf) = st.path_input.as_mut() {
                                buf.push(c);
                            }
                        }
                        KeyCode::Backspace => {
                            if let Some(buf) = st.path_input.as_mut() {
                                buf.pop();
                            }
                        }
                        KeyCode::Char('u') if key.modifiers.contains(KeyModifiers::CONTROL) => {
                            if let Some(buf) = st.path_input.as_mut() {
                                buf.clear();
                            }
                        }
                        KeyCode::Enter => {
                            let input = st.path_input.clone().unwrap_or_default();
                            let target = expand_tilde(&input);
                            if target.is_file() {
                                if let Some(parent) = target.parent() {
                                    crate::io::remember_picker_dir(parent);
                                }
                                self.modal = Self::import_state_for(target, st.append_mode);
                                return;
                            } else if target.is_dir() {
                                st.set_dir(target);
                            } else if let Some(hit) =
                                crate::io::resolve_library_query(&input, None)
                            {
                                self.modal = Self::import_state_for(hit.path, st.append_mode);
                                return;
                            } else {
                                // Not a path: reuse the text as a live search.
                                st.filter = input;
                                st.path_input = None;
                                st.change_filter();
                            }
                        }
                        KeyCode::Esc => st.path_input = None,
                        _ => {}
                    }
                    self.modal = Modal::FilePicker(st);
                    return;
                }

                match key.code {
                    KeyCode::Up => {
                        if !st.entries.is_empty() {
                            st.selected = (st.selected + st.entries.len() - 1) % st.entries.len();
                        }
                    }
                    KeyCode::Down => {
                        if !st.entries.is_empty() {
                            st.selected = (st.selected + 1) % st.entries.len();
                        }
                    }
                    KeyCode::PageUp => {
                        st.selected = st.selected.saturating_sub(10);
                    }
                    KeyCode::PageDown => {
                        st.selected = (st.selected + 10).min(st.entries.len().saturating_sub(1));
                    }
                    KeyCode::Home => st.selected = 0,
                    KeyCode::End => {
                        st.selected = st.entries.len().saturating_sub(1);
                    }
                    KeyCode::Char('a') if key.modifiers.contains(KeyModifiers::CONTROL) => {
                        st.append_mode = !st.append_mode;
                    }
                    KeyCode::Char('l') if key.modifiers.contains(KeyModifiers::CONTROL) => {
                        st.path_input = Some(String::new());
                    }
                    KeyCode::Char('u') if key.modifiers.contains(KeyModifiers::CONTROL) => {
                        st.filter.clear();
                        st.change_filter();
                    }
                    KeyCode::Char('/') | KeyCode::Char('~') => {
                        let seed = if key.code == KeyCode::Char('/') {
                            "/".to_string()
                        } else {
                            "~".to_string()
                        };
                        st.path_input = Some(seed);
                    }
                    KeyCode::Delete => {
                        st.filter.clear();
                        st.change_filter();
                    }
                    KeyCode::Backspace => {
                        if !st.filter.is_empty() {
                            st.filter.pop();
                            st.change_filter();
                        } else if let Some(parent) = st.dir.parent() {
                            st.set_dir(parent.to_path_buf());
                        }
                    }
                    KeyCode::Left => {
                        if let Some(parent) = st.dir.parent() {
                            st.set_dir(parent.to_path_buf());
                        }
                    }
                    KeyCode::Enter | KeyCode::Right => {
                        if let Some(entry) = st.entries.get(st.selected).cloned() {
                            if entry.is_dir {
                                st.set_dir(entry.path);
                            } else {
                                if let Some(parent) = entry.path.parent() {
                                    crate::io::remember_picker_dir(parent);
                                }
                                self.modal = Self::import_state_for(entry.path, st.append_mode);
                                return;
                            }
                        }
                    }
                    KeyCode::Char(c)
                        if !key.modifiers.contains(KeyModifiers::CONTROL)
                            && !c.is_control() =>
                    {
                        // Every printable key edits the live search immediately.
                        st.filter.push(c);
                        st.change_filter();
                    }
                    _ => {}
                }
                self.modal = Modal::FilePicker(st);
            }
            Modal::ImportConfig(mut st) => {
                let n_fields = if st.is_image { 5 } else { 2 };
                let text_field = st.is_image && (st.field_idx == 2 || st.field_idx == 3);
                let sync_x_bounds = |st: &mut ImportModalState| match st.kind_choice {
                    0 | 1 => {
                        st.x_left = "4000".into();
                        st.x_right = "400".into();
                    }
                    2 => {
                        st.x_left = "12.0".into();
                        st.x_right = "0.0".into();
                    }
                    3 => {
                        st.x_left = "220".into();
                        st.x_right = "0".into();
                    }
                    5 | 6 => {
                        st.x_left = "200".into();
                        st.x_right = "800".into();
                    }
                    _ => {}
                };

                match key.code {
                    KeyCode::Down | KeyCode::Tab => {
                        st.field_idx = (st.field_idx + 1) % n_fields;
                        self.modal = Modal::ImportConfig(st);
                    }
                    KeyCode::Up | KeyCode::BackTab => {
                        st.field_idx = (st.field_idx + n_fields - 1) % n_fields;
                        self.modal = Modal::ImportConfig(st);
                    }
                    KeyCode::Left if !text_field => {
                        match st.field_idx {
                            0 => {
                                let n_kinds = 7;
                                st.kind_choice = (st.kind_choice + n_kinds - 1) % n_kinds;
                                sync_x_bounds(&mut st);
                            }
                            1 => st.append_mode = !st.append_mode,
                            4 if st.is_image => st.bg_idx = (st.bg_idx + 2) % 3,
                            _ => {}
                        }
                        st.error = None;
                        self.modal = Modal::ImportConfig(st);
                    }
                    KeyCode::Right | KeyCode::Char(' ') if !text_field => {
                        match st.field_idx {
                            0 => {
                                let n_kinds = 7;
                                st.kind_choice = (st.kind_choice + 1) % n_kinds;
                                sync_x_bounds(&mut st);
                            }
                            1 => st.append_mode = !st.append_mode,
                            4 if st.is_image => st.bg_idx = (st.bg_idx + 1) % 3,
                            _ => {}
                        }
                        st.error = None;
                        self.modal = Modal::ImportConfig(st);
                    }
                    KeyCode::Char(c @ '1'..='7') if !text_field => {
                        st.kind_choice = (c as usize) - ('1' as usize);
                        sync_x_bounds(&mut st);
                        st.error = None;
                        self.modal = Modal::ImportConfig(st);
                    }
                    KeyCode::Backspace if text_field => {
                        if st.field_idx == 2 {
                            st.x_left.pop();
                        } else {
                            st.x_right.pop();
                        }
                        st.error = None;
                        self.modal = Modal::ImportConfig(st);
                    }
                    KeyCode::Char(c) if text_field && (c.is_ascii_digit() || c == '.' || c == '-') => {
                        if st.field_idx == 2 {
                            st.x_left.push(c);
                        } else {
                            st.x_right.push(c);
                        }
                        st.error = None;
                        self.modal = Modal::ImportConfig(st);
                    }
                    KeyCode::Enter => {
                        let kind_override = match st.kind_choice {
                            0 | 1 => Some(SpectrumKind::Ftir),
                            2 => Some(SpectrumKind::Hnmr),
                            3 => Some(SpectrumKind::Cnmr),
                            5 | 6 => Some(SpectrumKind::UvVis),
                            _ => None,
                        };
                        if !st.is_image {
                            match st.kind_choice {
                                0 => self.process_cfg.y_pref = YAxisPreference::Native,
                                1 => self.process_cfg.y_pref = YAxisPreference::Swapped,
                                // UV-Vis: Abs (5) keeps the default; %T (6) flips the toggle.
                                5 => self.process_cfg.y_pref = YAxisPreference::Native,
                                6 => self.process_cfg.y_pref = YAxisPreference::Swapped,
                                _ => {}
                            }
                            let path = st.path.clone();
                            let append = st.append_mode;
                            self.load_paths(
                                vec![path],
                                kind_override,
                                append,
                                "Opening spectrum file…",
                            );
                        } else {
                            let xl = st.x_left.trim().parse::<f64>();
                            let xr = st.x_right.trim().parse::<f64>();
                            match (xl, xr) {
                                (Ok(l), Ok(r)) if (l - r).abs() > 1e-6 => {
                                    let kind = kind_override.unwrap_or(SpectrumKind::Ftir);
                                    let y_unit = match st.kind_choice {
                                        1 | 5 => YUnit::Absorbance,
                                        2 | 3 => YUnit::Intensity,
                                        6 => YUnit::Transmittance,
                                        _ => YUnit::Transmittance,
                                    };
                                    let bg_mode = match st.bg_idx {
                                        1 => "light",
                                        2 => "dark",
                                        _ => "auto",
                                    }
                                    .to_string();
                                    match st.kind_choice {
                                        0 => self.process_cfg.y_pref = YAxisPreference::Native,
                                        1 => self.process_cfg.y_pref = YAxisPreference::Swapped,
                                        5 => self.process_cfg.y_pref = YAxisPreference::Native,
                                        6 => self.process_cfg.y_pref = YAxisPreference::Swapped,
                                        _ => {}
                                    }
                                    // Pure-Rust offline digitizer (Viterbi DP):
                                    // fast, but still runs off the UI thread so
                                    // large scans never block input.
                                    let path = st.path.clone();
                                    let append = st.append_mode;
                                    self.start_task("Digitizing spectrum image…", move || {
                                        match crate::digitize::digitize_spectrum_image(
                                            &path, kind, y_unit, l, r, &bg_mode,
                                        ) {
                                            Ok(spec) => TaskOutcome::Loaded {
                                                specs: vec![spec],
                                                errors: Vec::new(),
                                                append,
                                            },
                                            Err(e) => TaskOutcome::Loaded {
                                                specs: Vec::new(),
                                                errors: vec![format!(
                                                    "{}: {}",
                                                    path.display(),
                                                    e
                                                )],
                                                append,
                                            },
                                        }
                                    });
                                }
                                _ => {
                                    st.error = Some(
                                        "Enter valid distinct X Left and X Right numbers"
                                            .to_string(),
                                    );
                                    self.modal = Modal::ImportConfig(st);
                                }
                            }
                        }
                    }
                    _ => self.modal = Modal::ImportConfig(st),
                }
            }
            Modal::Export(mut st) => {
                // Fields: 0=Format, 1=Style, 2=Layout, 3=Range, 4=Peaks, 5=Title, 6=OutputPath
                let text_field = st.field_idx == 5 || st.field_idx == 6;
                match key.code {
                    KeyCode::Down | KeyCode::Tab => {
                        st.field_idx = (st.field_idx + 1) % 7;
                        self.modal = Modal::Export(st);
                    }
                    KeyCode::Up | KeyCode::BackTab => {
                        st.field_idx = (st.field_idx + 6) % 7;
                        self.modal = Modal::Export(st);
                    }
                    KeyCode::Left | KeyCode::Right | KeyCode::Char(' ') if !text_field => {
                        match st.field_idx {
                            0 => {
                                st.format = st.format.next();
                                let mut p = PathBuf::from(&st.output_path);
                                p.set_extension(st.format.ext());
                                st.output_path = p.display().to_string();
                            }
                            1 => st.style = st.style.next(),
                            2 => st.stacked = !st.stacked,
                            3 => st.use_zoom_range = !st.use_zoom_range,
                            4 => st.peak_mode = st.peak_mode.next(),
                            _ => {}
                        }
                        self.modal = Modal::Export(st);
                    }
                    KeyCode::Char('y') if !text_field => {
                        if let Some(ref p) = st.last_exported {
                            if p.extension().and_then(|e| e.to_str()) == Some("png") {
                                match copy_png_file(p) {
                                    Ok(()) => st.feedback = Some("Copied PNG image to Wayland clipboard!".into()),
                                    Err(e) => st.feedback = Some(format!("Clipboard error: {}", e)),
                                }
                            } else {
                                st.feedback = Some("Clipboard image copy requires PNG format".into());
                            }
                        }
                        self.modal = Modal::Export(st);
                    }
                    KeyCode::Char('o') if !text_field => {
                        if let Some(ref p) = st.last_exported {
                            match open_exported_file(p) {
                                Ok(()) => st.feedback = Some(format!("Opened {}", p.display())),
                                Err(e) => st.feedback = Some(format!("Open error: {}", e)),
                            }
                        }
                        self.modal = Modal::Export(st);
                    }
                    KeyCode::Backspace if text_field => {
                        if st.field_idx == 5 {
                            st.title.pop();
                        } else {
                            st.output_path.pop();
                        }
                        self.modal = Modal::Export(st);
                    }
                    KeyCode::Char(c) if text_field => {
                        if st.field_idx == 5 {
                            st.title.push(c);
                        } else {
                            st.output_path.push(c);
                        }
                        self.modal = Modal::Export(st);
                    }
                    KeyCode::Enter => {
                        self.modal = Modal::Export(st);
                        self.execute_export();
                    }
                    _ => self.modal = Modal::Export(st),
                }
            }
            Modal::LibrarySetup(mut st) => {
                match key.code {
                    KeyCode::Up | KeyCode::Char('k') => {
                        st.selected = (st.selected + st.checks.len() - 1) % st.checks.len();
                        st.feedback = None;
                        self.modal = Modal::LibrarySetup(st);
                    }
                    KeyCode::Down | KeyCode::Char('j') => {
                        st.selected = (st.selected + 1) % st.checks.len();
                        st.feedback = None;
                        self.modal = Modal::LibrarySetup(st);
                    }
                    KeyCode::Char(' ') => {
                        st.checks[st.selected] = !st.checks[st.selected];
                        st.feedback = None;
                        self.modal = Modal::LibrarySetup(st);
                    }
                    KeyCode::Char('a') => {
                        let all_checked = st.checks.iter().all(|&c| c);
                        st.checks.fill(!all_checked);
                        st.feedback = None;
                        self.modal = Modal::LibrarySetup(st);
                    }
                    KeyCode::Enter => {
                        // Build the list of (script, args) for checked packs
                        let packs: [(&str, &str); 5] = [
                            ("download_polymer_pack.py", "--nist"),      // 0: NIST FTIR (~2 MB)
                            ("download_uvvis_pack.py", "--all"),         // 1: UV-Vis PhotochemCAD + NIST (~2 MB)
                            ("download_polymer_pack.py", "--zenodo"),    // 2: Zenodo FTIR-Plastics (~253 MB)
                            ("download_polymer_pack.py", "--nmr"),       // 3: nmrshiftdb2 NMR (~159 MB)
                            ("download_mega_ir.py", "--all"),            // 4: Mega FTIR RRUFF+INFRA-ART+MNHN+NIST (~158 MB)
                        ];
                        let selected: Vec<(&str, &str)> = st
                            .checks
                            .iter()
                            .enumerate()
                            .filter(|(_, &c)| c)
                            .map(|(i, _)| packs[i])
                            .collect();
                        if selected.is_empty() {
                            st.feedback = Some("No packs selected — press [Space] to check at least one".into());
                            self.modal = Modal::LibrarySetup(st);
                            return;
                        }
                        // `start_task` already opens the Busy modal.
                        self.start_task("Downloading open reference libraries…", move || {
                            let result = crate::export::run_library_downloads(&selected);
                            TaskOutcome::LibraryDownload(result.map_err(|e| e.to_string()))
                        });
                    }
                    _ => {
                        st.feedback = None;
                        self.modal = Modal::LibrarySetup(st);
                    }
                }
            }
            Modal::Deconv(mut st) => {
                let editing = (3..=4).contains(&st.field_idx);
                match key.code {
                    KeyCode::Up | KeyCode::BackTab | KeyCode::Char('k') if !key.modifiers.contains(KeyModifiers::CONTROL) => {
                        st.field_idx = (st.field_idx + 4) % 5;
                        self.modal = Modal::Deconv(st);
                    }
                    KeyCode::Down | KeyCode::Tab | KeyCode::Char('j') => {
                        st.field_idx = (st.field_idx + 1) % 5;
                        self.modal = Modal::Deconv(st);
                    }
                    KeyCode::Left | KeyCode::Right | KeyCode::Char(' ') if st.field_idx < 3 => {
                        let forward = matches!(key.code, KeyCode::Right | KeyCode::Char(' '));
                        match st.field_idx {
                            0 => {
                                st.shape = if forward { st.shape.next() } else { st.shape.prev() };
                            }
                            1 => {
                                st.max_peaks = if forward {
                                    (st.max_peaks + 1).min(8)
                                } else {
                                    st.max_peaks.saturating_sub(1).max(1)
                                };
                            }
                            2 => st.fit_baseline = !st.fit_baseline,
                            _ => {}
                        }
                        self.modal = Modal::Deconv(st);
                    }
                    KeyCode::Char(c) if editing && (c.is_ascii_digit() || c == '.' || c == '-') => {
                        let buf = if st.field_idx == 3 { &mut st.x_lo } else { &mut st.x_hi };
                        if c != '-' || buf.is_empty() {
                            buf.push(c);
                        }
                        self.modal = Modal::Deconv(st);
                    }
                    KeyCode::Backspace if editing => {
                        if st.field_idx == 3 {
                            st.x_lo.pop();
                        } else {
                            st.x_hi.pop();
                        }
                        self.modal = Modal::Deconv(st);
                    }
                    KeyCode::Enter => {
                        let Some(spec) = self.traces.get(self.active_trace).cloned() else {
                            st.feedback = Some("No active spectrum to fit".into());
                            self.modal = Modal::Deconv(st);
                            return;
                        };
                        let Some(proc_t) = self.processed.get(self.active_trace).cloned() else {
                            st.feedback = Some("No processed data available".into());
                            self.modal = Modal::Deconv(st);
                            return;
                        };
                        let (Ok(x_lo), Ok(x_hi)) = (st.x_lo.trim().parse::<f64>(), st.x_hi.trim().parse::<f64>()) else {
                            st.feedback = Some("Enter valid numeric X window bounds".into());
                            self.modal = Modal::Deconv(st);
                            return;
                        };
                        let cfg = DeconvConfig {
                            shape: st.shape,
                            max_peaks: st.max_peaks.clamp(1, 8),
                            fit_baseline: st.fit_baseline,
                            manual_centers: Vec::new(),
                            max_iter: 150,
                        };
                        let solvent = self.nmr_solvent;
                        self.pending_deconv = Some(st);
                        self.start_task("Running Levenberg-Marquardt peak deconvolution…", move || {
                            let result = deconvolve_window(
                                spec.kind,
                                &proc_t.x,
                                &proc_t.y,
                                proc_t.y_unit,
                                x_lo,
                                x_hi,
                                &cfg,
                                solvent,
                            )
                            .map_err(|e| e.to_string());
                            TaskOutcome::Deconv(result)
                        });
                    }
                    KeyCode::Char('c') if st.field_idx < 3 => {
                        if let Some(res) = &st.result {
                            match copy_text(&format_deconv_markdown(res)) {
                                Ok(_) => st.feedback = Some("Copied fit report to clipboard".into()),
                                Err(e) => st.feedback = Some(format!("Clipboard error: {}", e)),
                            }
                        }
                        self.modal = Modal::Deconv(st);
                    }
                    KeyCode::Char('s') if st.field_idx < 3 => {
                        if let Some(res) = &st.result {
                            let saved = (|| -> std::io::Result<PathBuf> {
                                let dir = PathBuf::from("exports");
                                fs::create_dir_all(&dir)?;
                                let stamp = SystemTime::now()
                                    .duration_since(UNIX_EPOCH)
                                    .unwrap_or_default()
                                    .as_secs();
                                let path = dir.join(format!("deconv_fit_{}.csv", stamp));
                                fs::write(&path, format_deconv_csv(res))?;
                                Ok(path)
                            })();
                            st.feedback = Some(match saved {
                                Ok(path) => format!("Saved CSV: {}", path.display()),
                                Err(e) => format!("CSV save failed: {}", e),
                            });
                        }
                        self.modal = Modal::Deconv(st);
                    }
                    KeyCode::Char('r') if st.field_idx < 3 => {
                        match &mut self.deconv_overlay {
                            Some(ov) => {
                                ov.show_residual = !ov.show_residual;
                                st.feedback = Some(format!(
                                    "Residual inset {}",
                                    if ov.show_residual { "shown" } else { "hidden" }
                                ));
                            }
                            None => st.feedback = Some("Fit the window first to show a residual".into()),
                        }
                        self.modal = Modal::Deconv(st);
                    }
                    KeyCode::Char('x') if st.field_idx < 3 => {
                        self.deconv_overlay = None;
                        st.result = None;
                        st.feedback = Some("Deconvolution overlay cleared".into());
                        self.modal = Modal::Deconv(st);
                    }
                    _ => self.modal = Modal::Deconv(st),
                }
            }
            Modal::SpectralMath(mut st) => {
                match key.code {
                    KeyCode::Up | KeyCode::BackTab => {
                        st.field_idx = (st.field_idx + 3) % 4;
                        self.modal = Modal::SpectralMath(st);
                    }
                    KeyCode::Down | KeyCode::Tab => {
                        st.field_idx = (st.field_idx + 1) % 4;
                        self.modal = Modal::SpectralMath(st);
                    }
                    KeyCode::Left | KeyCode::Right => {
                        let forward = key.code == KeyCode::Right;
                        match st.field_idx {
                            0 => st.op = if forward { st.op.next() } else { st.op.prev() },
                            1 if !self.traces.is_empty() => {
                                let len = self.traces.len();
                                st.trace_a = if forward {
                                    (st.trace_a + 1) % len
                                } else {
                                    (st.trace_a + len - 1) % len
                                };
                            }
                            2 if !self.traces.is_empty() => {
                                let len = self.traces.len();
                                st.trace_b = if forward {
                                    (st.trace_b + 1) % len
                                } else {
                                    (st.trace_b + len - 1) % len
                                };
                            }
                            3 => {
                                let step = if key.modifiers.contains(KeyModifiers::SHIFT) { 0.01 } else { 0.05 };
                                st.factor_k = (st.factor_k + if forward { step } else { -step }).clamp(0.0, 20.0);
                            }
                            _ => {}
                        }
                        self.modal = Modal::SpectralMath(st);
                    }
                    KeyCode::Char('k') => {
                        if let (Some(a), Some(b)) = (self.traces.get(st.trace_a), self.traces.get(st.trace_b)) {
                            if let Some(k) = optimal_subtraction_factor(a, b, self.view_x_min, self.view_x_max) {
                                st.factor_k = k;
                                st.feedback = Some(format!("Auto-optimal k* = {:.4}", k));
                            } else {
                                st.feedback = Some("Could not estimate k* in the selected overlap".into());
                            }
                        }
                        self.modal = Modal::SpectralMath(st);
                    }
                    KeyCode::Enter => {
                        let result = if st.op == SpectralOp::Average {
                            let refs: Vec<&Spectrum> = self.traces.iter().filter(|t| t.visible).collect();
                            average_spectra(&refs)
                        } else {
                            self.traces.get(st.trace_a).zip(self.traces.get(st.trace_b))
                                .and_then(|(a, b)| compute_spectral_math(a, b, st.op, st.factor_k))
                        };
                        if let Some(mut spec) = result {
                            spec.color_idx = self.traces.len();
                            let name = spec.name.clone();
                            self.traces.push(spec);
                            self.active_trace = self.traces.len() - 1;
                            self.deconv_overlay = None;
                            self.recompute_processing();
                            self.reset_zoom();
                            self.modal = Modal::None;
                            self.status_msg = format!("Spectral Math result added: {}", name);
                        } else {
                            st.feedback = Some("No compatible overlapping data for this operation".into());
                            self.modal = Modal::SpectralMath(st);
                        }
                    }
                    _ => self.modal = Modal::SpectralMath(st),
                }
            }
            Modal::Integration(mut st) => {
                match key.code {
                    KeyCode::Char('z') => {
                        let idx = self.add_integration_region(self.view_x_min, self.view_x_max);
                        match idx {
                            Some(i) => {
                                st.selected = i;
                                st.feedback = Some(format!("Added zoom-window region #{}", i + 1));
                            }
                            None => st.feedback = Some("No spectral data in the zoom window".into()),
                        }
                        self.modal = Modal::Integration(st);
                    }
                    KeyCode::Char('a') => {
                        let half = match self.current_kind() {
                            SpectrumKind::Ftir => 25.0,
                            SpectrumKind::Hnmr => 0.08,
                            SpectrumKind::Cnmr => 1.5,
                            SpectrumKind::UvVis => 20.0,
                        };
                        let idx = self
                            .add_integration_region(self.cursor_x - half, self.cursor_x + half);
                        match idx {
                            Some(i) => {
                                st.selected = i;
                                st.feedback =
                                    Some(format!("Added cursor region #{} (±{:.2})", i + 1, half));
                            }
                            None => st.feedback = Some("No spectral data at the cursor".into()),
                        }
                        self.modal = Modal::Integration(st);
                    }
                    KeyCode::Up | KeyCode::Char('k') => {
                        st.selected = st.selected.saturating_sub(1);
                        self.modal = Modal::Integration(st);
                    }
                    KeyCode::Down | KeyCode::Char('j') => {
                        if !self.integration_regions.is_empty() {
                            st.selected = (st.selected + 1).min(self.integration_regions.len() - 1);
                        }
                        self.modal = Modal::Integration(st);
                    }
                    KeyCode::Char('r') => {
                        if !self.integration_regions.is_empty() {
                            st.selected = st.selected.min(self.integration_regions.len() - 1);
                            self.integration_ref_idx = st.selected;
                            st.feedback =
                                Some(format!("Region #{} is now the reference ★", st.selected + 1));
                        }
                        self.modal = Modal::Integration(st);
                    }
                    KeyCode::Left | KeyCode::Right => {
                        let step = if key.code == KeyCode::Right { 0.5 } else { -0.5 };
                        self.integration_ref_val =
                            (self.integration_ref_val + step).clamp(0.1, 10_000.0);
                        st.feedback = Some(format!(
                            "Reference integral = {:.2}",
                            self.integration_ref_val
                        ));
                        self.modal = Modal::Integration(st);
                    }
                    KeyCode::Char('d') => {
                        if !self.integration_regions.is_empty() {
                            let idx = st.selected.min(self.integration_regions.len() - 1);
                            self.integration_regions.remove(idx);
                            if self.integration_regions.is_empty() {
                                self.integration_ref_idx = 0;
                            } else {
                                self.integration_ref_idx =
                                    self.integration_ref_idx.min(self.integration_regions.len() - 1);
                                st.selected = idx.min(self.integration_regions.len() - 1);
                            }
                            st.feedback = Some(format!("Deleted region #{}", idx + 1));
                        }
                        self.modal = Modal::Integration(st);
                    }
                    KeyCode::Char('x') => {
                        self.clear_integrations();
                        st.selected = 0;
                        st.feedback = Some("Cleared all integration regions".into());
                        self.modal = Modal::Integration(st);
                    }
                    KeyCode::Char('c') => {
                        if self.integration_regions.is_empty() {
                            st.feedback = Some("No regions to copy".into());
                        } else {
                            let text = self.integration_table_text();
                            match copy_text(&text) {
                                Ok(tool) => st.feedback = Some(format!("Copied table via {}", tool)),
                                Err(e) => st.feedback = Some(format!("Copy failed: {}", e)),
                            }
                        }
                        self.modal = Modal::Integration(st);
                    }
                    KeyCode::Char('p') => {
                        match self.push_integrations_to_polymer_lab() {
                            Some(msg) => {
                                self.status_msg = msg.clone();
                                st.feedback = Some(msg);
                                if self.tab == ActiveTab::PolymerLab {
                                    self.modal = Modal::None;
                                    return;
                                }
                            }
                            None => {
                                st.feedback = Some(
                                    "Need at least 2 regions (3 for ¹H tacticity push)".into(),
                                );
                            }
                        }
                        self.modal = Modal::Integration(st);
                    }
                    _ => self.modal = Modal::Integration(st),
                }
            }
            Modal::Calibration(mut st) => {
                match key.code {
                    KeyCode::Char(c) if c.is_ascii_digit() || c == '.' || c == '-' => {
                        st.target_buf.push(c);
                        st.feedback = None;
                        self.modal = Modal::Calibration(st);
                    }
                    KeyCode::Backspace => {
                        st.target_buf.pop();
                        st.feedback = None;
                        self.modal = Modal::Calibration(st);
                    }
                    KeyCode::Char('t') => {
                        let target = match self.current_kind() {
                            SpectrumKind::Ftir => 1601.4,
                            SpectrumKind::UvVis => 361.5,
                            SpectrumKind::Hnmr | SpectrumKind::Cnmr => 0.0,
                        };
                        st.target_buf = format!("{:.4}", target);
                        st.feedback = Some(match self.current_kind() {
                            SpectrumKind::Ftir => "Polystyrene 1601.4 cm⁻¹".into(),
                            SpectrumKind::UvVis => "Holmium oxide 361.5 nm".into(),
                            _ => "TMS 0.00 ppm".into(),
                        });
                        self.modal = Modal::Calibration(st);
                    }
                    KeyCode::Char('s') => {
                        let (target, label) = match self.current_kind() {
                            SpectrumKind::Ftir => (3026.0, "Polystyrene 3026.0 cm⁻¹"),
                            SpectrumKind::UvVis => (536.3, "Holmium oxide 536.3 nm"),
                            SpectrumKind::Hnmr => (
                                self.nmr_solvent.residual_h().first().copied().unwrap_or(7.26),
                                "Active solvent residual ¹H",
                            ),
                            SpectrumKind::Cnmr => (
                                self.nmr_solvent.residual_c().first().copied().unwrap_or(77.16),
                                "Active solvent residual ¹³C",
                            ),
                        };
                        st.target_buf = format!("{:.4}", target);
                        st.feedback = Some(label.into());
                        self.modal = Modal::Calibration(st);
                    }
                    KeyCode::Char('w') => {
                        let (target, label) = match self.current_kind() {
                            SpectrumKind::Ftir => (2349.0, "Atmospheric CO₂ 2349.0 cm⁻¹"),
                            SpectrumKind::Hnmr => {
                                (self.nmr_solvent.water_h(), "Residual H₂O in active solvent")
                            }
                            SpectrumKind::Cnmr => (
                                self.nmr_solvent.residual_c().first().copied().unwrap_or(77.16),
                                "Active solvent residual ¹³C",
                            ),
                            SpectrumKind::UvVis => (400.0, "UV-Vis 400.0 nm"),
                        };
                        st.target_buf = format!("{:.4}", target);
                        st.feedback = Some(label.into());
                        self.modal = Modal::Calibration(st);
                    }
                    KeyCode::Enter => match st.target_buf.trim().parse::<f64>() {
                        Ok(target) => {
                            let msg = self.apply_calibration_target(st.observed_x, target);
                            self.status_msg = msg;
                            self.modal = Modal::None;
                            return;
                        }
                        Err(_) => {
                            st.feedback = Some("Enter a valid numeric target X".into());
                            self.modal = Modal::Calibration(st);
                        }
                    },
                    _ => self.modal = Modal::Calibration(st),
                }
            }
            Modal::Busy(busy) => {
                // Swallow input while a background task runs; the outcome
                // arrives through `poll_background`.
                self.modal = Modal::Busy(busy);
            }
            Modal::None => {}
        }
    }

    // ========================================================================
    // MOUSE EVENT HANDLER (Click Cursor, Drag-to-Zoom, Wheel Zoom, Top Tabs)
    // ========================================================================
    pub fn handle_mouse(&mut self, mouse: MouseEvent) {
        if !matches!(self.modal, Modal::None) {
            return;
        }

        let col = mouse.column;
        let row = mouse.row;

        // Topbar click -> switch tabs or trigger actions
        if matches!(mouse.kind, MouseEventKind::Down(MouseButton::Left))
            && row == self.topbar_rect.y
        {
            // Tab strip shown by `draw_topbar`; hitboxes are computed there.
            if let Some(&(_, _, tab)) = self
                .tab_hitboxes
                .iter()
                .find(|&&(x0, x1, _)| col >= x0 && col < x1)
            {
                self.switch_tab(tab);
                return;
            }
        }

        if self.tab.as_kind().is_none() || self.plot_rect.width < 4 {
            return;
        }

        let inside_plot = col >= self.plot_rect.x
            && col < self.plot_rect.x + self.plot_rect.width
            && row >= self.plot_rect.y
            && row < self.plot_rect.y + self.plot_rect.height;

        let map_col_to_x = |c: u16| -> f64 {
            let rel = (c.saturating_sub(self.plot_rect.x) as f64)
                / (self.plot_rect.width.saturating_sub(1).max(1) as f64);
            let u = rel.clamp(0.0, 1.0);
            if self.current_kind().is_x_reversed() {
                // Reversed IUPAC X axis: left=view_x_max, right=view_x_min
                self.view_x_max - u * (self.view_x_max - self.view_x_min)
            } else {
                // Ascending axis (UV-Vis nm): left=view_x_min, right=view_x_max
                self.view_x_min + u * (self.view_x_max - self.view_x_min)
            }
        };

        match mouse.kind {
            MouseEventKind::Down(MouseButton::Left) if inside_plot => {
                let x = map_col_to_x(col);
                self.cursor_x = x;
                self.drag_start_x = Some(x);
                self.drag_curr_x = Some(x);
            }
            MouseEventKind::Drag(MouseButton::Left) => {
                if self.drag_start_x.is_some() {
                    let x = map_col_to_x(col);
                    self.cursor_x = x;
                    self.drag_curr_x = Some(x);
                }
            }
            MouseEventKind::Up(MouseButton::Left) => {
                if let (Some(x1), Some(x2)) = (self.drag_start_x.take(), self.drag_curr_x.take()) {
                    let span = (self.view_x_max - self.view_x_min).max(1e-6);
                    if (x2 - x1).abs() > span * 0.02 {
                        self.view_x_min = x1.min(x2);
                        self.view_x_max = x1.max(x2);
                        self.cursor_x = 0.5 * (x1 + x2);
                        let (left_lbl, right_lbl) = if self.current_kind().is_x_reversed() {
                            (self.view_x_max, self.view_x_min)
                        } else {
                            (self.view_x_min, self.view_x_max)
                        };
                        self.status_msg = format!(
                            "Zoomed: {:.1} – {:.1} {}",
                            left_lbl,
                            right_lbl,
                            self.current_kind().x_unit()
                        );
                    }
                }
            }
            MouseEventKind::Down(MouseButton::Right) if inside_plot => {
                self.reset_zoom();
                self.status_msg = "Right-click: reset zoom".into();
            }
            MouseEventKind::ScrollUp if inside_plot => {
                self.cursor_x = map_col_to_x(col);
                self.zoom_around_cursor(0.75);
            }
            MouseEventKind::ScrollDown if inside_plot => {
                self.cursor_x = map_col_to_x(col);
                self.zoom_around_cursor(1.35);
            }
            _ => {}
        }
    }
}

fn scan_directory(dir: &Path) -> Vec<FileEntry> {
    let mut entries = Vec::new();
    if let Some(parent) = dir.parent() {
        entries.push(FileEntry {
            path: parent.to_path_buf(),
            display: "../ (parent directory)".into(),
            is_dir: true,
        });
    }
    if let Ok(rd) = fs::read_dir(dir) {
        let mut dirs = Vec::new();
        let mut files = Vec::new();
        for item in rd.flatten() {
            let path = item.path();
            let name = item.file_name().to_string_lossy().to_string();
            if name.starts_with('.') {
                continue;
            }
            // `DirEntry::file_type` reads `d_type` straight from the directory
            // listing, unlike `path.is_dir()`, which issues one blocking stat
            // per entry — painfully slow on stale/network mounts.
            let is_dir = match item.file_type() {
                Ok(ft) if ft.is_symlink() => path.is_dir(),
                Ok(ft) => ft.is_dir(),
                Err(_) => continue,
            };
            if is_dir {
                dirs.push(FileEntry {
                    path,
                    display: format!("📁 {}/", name),
                    is_dir: true,
                });
            } else {
                let ext = path
                    .extension()
                    .and_then(|e| e.to_str())
                    .unwrap_or("")
                    .to_ascii_lowercase();
                if matches!(ext.as_str(), "csv" | "tsv" | "txt" | "dpt" | "asc" | "xy" | "jdx" | "dx" | "jcamp") {
                    files.push(FileEntry {
                        path,
                        display: format!("📈 {}", name),
                        is_dir: false,
                    });
                } else if matches!(ext.as_str(), "png" | "jpg" | "jpeg" | "webp" | "bmp" | "tif" | "tiff") {
                    files.push(FileEntry {
                        path,
                        display: format!("🖼 {}", name),
                        is_dir: false,
                    });
                }
            }
        }
        dirs.sort_by(|a, b| a.display.cmp(&b.display));
        files.sort_by(|a, b| a.display.cmp(&b.display));
        entries.extend(dirs);
        entries.extend(files);
    }
    entries
}

fn build_search_results(query: &str) -> Vec<SearchItem> {
    let q = query.trim().to_ascii_lowercase();
    let mut out = Vec::new();

    if !q.is_empty() {
        // A dropped / pasted path should show up as an openable result.
        let candidate = expand_tilde(query);
        if candidate.is_file() {
            out.push(SearchItem {
                category: "File",
                target_tab: ActiveTab::Ftir,
                jump_x: None,
                title: candidate
                    .file_name()
                    .and_then(|s| s.to_str())
                    .unwrap_or("spectrum")
                    .to_string(),
                detail: candidate.display().to_string(),
                polymer_short: None,
                file_path: Some(candidate),
            });
        }
        // Reference spectra shipped in `library/ftir` and `library/nmr`.
        for hit in search_library(query).into_iter().take(15) {
            out.push(SearchItem {
                category: match hit.kind {
                    SpectrumKind::Ftir => "Library · FTIR",
                    SpectrumKind::Hnmr => "Library · ¹H",
                    SpectrumKind::Cnmr => "Library · ¹³C",
                    SpectrumKind::UvVis => "Library · UV-Vis",
                },
                target_tab: ActiveTab::from_kind(hit.kind),
                jump_x: None,
                title: hit.stem.replace('_', " "),
                detail: hit.path.display().to_string(),
                polymer_short: None,
                file_path: Some(hit.path),
            });
        }
    }

    for b in FTIR_DB {
        let hay = format!("{} {} {} {:.0}-{:.0}", b.group, b.mode, b.note, b.min, b.max).to_ascii_lowercase();
        if q.is_empty() || hay.contains(&q) {
            out.push(SearchItem {
                category: "FTIR",
                target_tab: ActiveTab::Ftir,
                jump_x: Some(0.5 * (b.min + b.max)),
                title: format!("{} ({:.0}–{:.0} cm⁻¹)", b.group, b.min, b.max),
                detail: b.note.to_string(),
                polymer_short: None,
                file_path: None,
            });
        }
    }

    for b in HNMR_DB {
        let hay = format!("{} {} {}", b.group, b.example, b.note).to_ascii_lowercase();
        if q.is_empty() || hay.contains(&q) {
            out.push(SearchItem {
                category: "¹H NMR",
                target_tab: ActiveTab::Hnmr,
                jump_x: Some(0.5 * (b.min + b.max)),
                title: format!("{} ({:.2}–{:.2} ppm)", b.group, b.min, b.max),
                detail: format!("{} · {}", b.example, b.note),
                polymer_short: None,
                file_path: None,
            });
        }
    }

    for b in CNMR_DB {
        let hay = format!("{} {} {}", b.group, b.example, b.note).to_ascii_lowercase();
        if q.is_empty() || hay.contains(&q) {
            out.push(SearchItem {
                category: "¹³C NMR",
                target_tab: ActiveTab::Cnmr,
                jump_x: Some(0.5 * (b.min + b.max)),
                title: format!("{} ({:.0}–{:.0} ppm)", b.group, b.min, b.max),
                detail: format!("{} · {}", b.example, b.note),
                polymer_short: None,
                file_path: None,
            });
        }
    }

    for b in UVVIS_DB {
        let hay = format!(
            "{} {} {} {} {:.0}-{:.0}",
            b.group, b.transition, b.note, b.epsilon, b.min, b.max
        )
        .to_ascii_lowercase();
        if q.is_empty() || hay.contains(&q) {
            out.push(SearchItem {
                category: "UV-Vis",
                target_tab: ActiveTab::UvVis,
                jump_x: Some(0.5 * (b.min + b.max)),
                title: format!("{} ({:.0}–{:.0} nm)", b.group, b.min, b.max),
                detail: format!("{} · ε {} · {}", b.transition, b.epsilon, b.note),
                polymer_short: None,
                file_path: None,
            });
        }
    }

    for p in POLYMERS_DB {
        let hay = format!("{} {} {}", p.short, p.name, p.repeat_unit).to_ascii_lowercase();
        if q.is_empty() || hay.contains(&q) {
            let main_peak = p.ftir_peaks.last().map(|(x, _)| *x);
            out.push(SearchItem {
                category: "Polymer",
                target_tab: ActiveTab::Ftir,
                jump_x: main_peak,
                title: format!("{} — {} (M₀={:.2})", p.short, p.name, p.m0),
                detail: format!("{} | {}", p.repeat_unit, p.tg_tm),
                polymer_short: Some(p.short),
                file_path: None,
            });
        }
    }

    for imp in SOLVENT_IMPURITIES {
        let hay = format!("{} {} {}", imp.name, imp.group, imp.note).to_ascii_lowercase();
        if q.is_empty() || hay.contains(&q) {
            out.push(SearchItem {
                category: "Solvent",
                target_tab: ActiveTab::Solvents,
                jump_x: imp.h_shifts[0],
                title: format!("{} ({}) — CDCl₃: {}",
                    imp.name,
                    imp.group,
                    imp.h_shifts[0].map(|v| format!("{:.2} ppm", v)).unwrap_or_else(|| "—".into())
                ),
                detail: imp.note.to_string(),
                polymer_short: None,
                file_path: None,
            });
        }
    }

    out.truncate(80);
    out
}
