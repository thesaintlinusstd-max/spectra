use std::env;
use std::fmt::Write as _;
use std::fs;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::time::{SystemTime, UNIX_EPOCH};

use anyhow::{anyhow, Context, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportStyle {
    Acs,
    Rsc,
    Everforest,
}

impl ExportStyle {
    pub fn next(self) -> Self {
        match self {
            Self::Acs => Self::Rsc,
            Self::Rsc => Self::Everforest,
            Self::Everforest => Self::Acs,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Acs => "ACS / Macromolecules",
            Self::Rsc => "RSC / Polymer Chemistry",
            Self::Everforest => "Dark Everforest (Slides)",
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Acs => "acs",
            Self::Rsc => "rsc",
            Self::Everforest => "everforest",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportFormat {
    Png,
    Svg,
    Pdf,
}

impl ExportFormat {
    pub fn next(self) -> Self {
        match self {
            Self::Png => Self::Svg,
            Self::Svg => Self::Pdf,
            Self::Pdf => Self::Png,
        }
    }

    pub fn ext(self) -> &'static str {
        match self {
            Self::Png => "png",
            Self::Svg => "svg",
            Self::Pdf => "pdf",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Png => "PNG (600 DPI)",
            Self::Svg => "SVG (Vector)",
            Self::Pdf => "PDF (Journal/LaTeX)",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportPeakMode {
    PinnedAndAuto,
    PinnedOnly,
    None,
}

impl ExportPeakMode {
    pub fn next(self) -> Self {
        match self {
            Self::PinnedAndAuto => Self::PinnedOnly,
            Self::PinnedOnly => Self::None,
            Self::None => Self::PinnedAndAuto,
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::PinnedAndAuto => "Pinned + Top Auto Peaks",
            Self::PinnedOnly => "Pinned Peaks Only",
            Self::None => "No Annotations",
        }
    }
}

#[derive(Debug, Clone)]
pub struct ExportTrace {
    pub name: String,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct ExportPeak {
    pub x: f64,
    pub y: f64,
    pub label: String,
    pub pinned: bool,
}

/// Deconvolved sub-peaks rendered alongside the measured trace.
#[derive(Debug, Clone)]
pub struct ExportDeconv {
    pub x: Vec<f64>,
    pub y_fit: Vec<f64>,
    pub components: Vec<Vec<f64>>,
}

#[derive(Debug, Clone)]
pub struct ExportPayload {
    pub mode: String,
    pub y_unit: String,
    pub style: String,
    pub stacked: bool,
    pub dpi: u32,
    pub title: String,
    pub x_bounds: [f64; 2],
    pub active_trace: usize,
    pub traces: Vec<ExportTrace>,
    pub peaks: Vec<ExportPeak>,
    pub deconv: Option<ExportDeconv>,
}

fn push_json_str(out: &mut String, s: &str) {
    out.push('"');
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

fn push_json_f64(out: &mut String, v: f64) {
    if v.is_finite() {
        let _ = write!(out, "{:.6}", v);
    } else {
        out.push_str("0.0");
    }
}

fn push_json_f64_slice(out: &mut String, vals: &[f64]) {
    out.push('[');
    for (i, &v) in vals.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        push_json_f64(out, v);
    }
    out.push(']');
}

impl ExportPayload {
    pub fn to_json_string(&self) -> String {
        let mut out = String::with_capacity(16384);
        out.push('{');
        out.push_str("\"mode\":");
        push_json_str(&mut out, &self.mode);
        out.push_str(",\"y_unit\":");
        push_json_str(&mut out, &self.y_unit);
        out.push_str(",\"style\":");
        push_json_str(&mut out, &self.style);
        let _ = write!(
            out,
            ",\"stacked\":{},\"dpi\":{},\"title\":",
            self.stacked, self.dpi
        );
        push_json_str(&mut out, &self.title);
        out.push_str(",\"x_bounds\":[");
        push_json_f64(&mut out, self.x_bounds[0]);
        out.push(',');
        push_json_f64(&mut out, self.x_bounds[1]);
        let _ = write!(out, "],\"active_trace\":{},\"traces\":[", self.active_trace);
        for (i, tr) in self.traces.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str("{\"name\":");
            push_json_str(&mut out, &tr.name);
            out.push_str(",\"x\":");
            push_json_f64_slice(&mut out, &tr.x);
            out.push_str(",\"y\":");
            push_json_f64_slice(&mut out, &tr.y);
            out.push('}');
        }
        out.push_str("],\"peaks\":[");
        for (i, pk) in self.peaks.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            out.push_str("{\"x\":");
            push_json_f64(&mut out, pk.x);
            out.push_str(",\"y\":");
            push_json_f64(&mut out, pk.y);
            out.push_str(",\"label\":");
            push_json_str(&mut out, &pk.label);
            let _ = write!(out, ",\"pinned\":{}}}", pk.pinned);
        }
        out.push(']');
        if let Some(dc) = &self.deconv {
            out.push_str(",\"deconv\":{\"x\":");
            push_json_f64_slice(&mut out, &dc.x);
            out.push_str(",\"y_fit\":");
            push_json_f64_slice(&mut out, &dc.y_fit);
            out.push_str(",\"components\":[");
            for (i, comp) in dc.components.iter().enumerate() {
                if i > 0 {
                    out.push(',');
                }
                push_json_f64_slice(&mut out, comp);
            }
            out.push_str("]}");
        }
        out.push('}');
        out
    }
}

pub fn find_python_script(script_name: &str) -> Result<(PathBuf, PathBuf)> {
    let mut candidates: Vec<PathBuf> = Vec::new();
    let rel = PathBuf::from("python").join(script_name);

    if let Ok(home) = env::var("SPECTRA_HOME") {
        candidates.push(PathBuf::from(home).join(&rel));
    }
    candidates.push(rel.clone());
    if let Ok(exe) = env::current_exe() {
        if let Some(root) = exe.parent().and_then(|p| p.parent()).and_then(|p| p.parent()) {
            candidates.push(root.join(&rel));
        }
        if let Some(dir) = exe.parent() {
            candidates.push(dir.join(&rel));
        }
    }
    candidates.push(PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(&rel));

    let script = candidates
        .into_iter()
        .find(|p| p.is_file())
        .ok_or_else(|| anyhow!("Could not locate python/{}", script_name))?;

    let venv_py = script
        .parent()
        .and_then(|p| p.parent())
        .map(|root| root.join(".venv/bin/python"));

    let python = match venv_py {
        Some(p) if p.is_file() => p,
        _ => PathBuf::from("python3"),
    };

    Ok((python, script))
}

pub fn find_python_and_script() -> Result<(PathBuf, PathBuf)> {
    find_python_script("export_plot.py")
}

pub fn run_export(payload: &ExportPayload, output_path: &Path) -> Result<PathBuf> {
    let is_svg = output_path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.eq_ignore_ascii_case("svg"))
        .unwrap_or(false);
    if is_svg {
        return crate::svg_export::write_publication_svg(payload, output_path);
    }
    run_python_export(payload, output_path)
}

pub fn run_python_export(payload: &ExportPayload, output_path: &Path) -> Result<PathBuf> {
    if let Some(parent) = output_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }

    let (python, script) = find_python_and_script()?;
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    let tmp_json = env::temp_dir().join(format!(
        "spectra_export_{}_{}.json",
        std::process::id(),
        stamp
    ));
    fs::write(&tmp_json, payload.to_json_string())?;

    let out = Command::new(&python)
        .arg(&script)
        .arg("--input")
        .arg(&tmp_json)
        .arg("--output")
        .arg(output_path)
        .output()
        .with_context(|| format!("Failed executing {} {}", python.display(), script.display()))?;

    let _ = fs::remove_file(&tmp_json);

    if !out.status.success() {
        let stderr = String::from_utf8_lossy(&out.stderr);
        return Err(anyhow!("Python export failed: {}", stderr.trim()));
    }

    if !output_path.is_file() {
        return Err(anyhow!("Exporter finished but output file was not created"));
    }

    Ok(output_path.to_path_buf())
}

pub fn open_exported_file(path: &Path) -> Result<()> {
    Command::new("xdg-open")
        .arg(path)
        .stdin(Stdio::null())
        .stdout(Stdio::null())
        .stderr(Stdio::null())
        .spawn()
        .context("Failed to launch xdg-open")?;
    Ok(())
}

// ============================================================================
// Library Downloader Integration
// ============================================================================

#[allow(dead_code)]
/// Run selected library download packs and return a summary string.
/// Items are tuples of (script_name, args_as_single_string).
/// Returns a clean summary with total files downloaded per pack.
pub fn run_library_downloads(items: &[(&str, &str)]) -> Result<String> {
    use std::process::Command;

    let mut summary_lines: Vec<String> = Vec::new();
    let before = crate::io::library_index().len();

    for (script_name, args_str) in items {
        let (python, script) = find_python_script(script_name)?;
        let args: Vec<&str> = if args_str.is_empty() {
            vec![]
        } else {
            args_str.split_whitespace().collect()
        };
        let mut cmd = Command::new(&python);
        cmd.arg(&script).args(&args);

        let output = cmd
            .output()
            .with_context(|| format!("Failed executing {} {}", python.display(), script.display()))?;

        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);

        if !output.status.success() {
            let err_msg = stderr.trim();
            if err_msg.is_empty() {
                summary_lines.push(format!("✗ {}: failed (exit code {})", script_name, output.status));
            } else {
                summary_lines.push(format!("✗ {}: {}", script_name, err_msg));
            }
            continue;
        }

        // Scripts may print an itemised summary, but the authoritative count is
        // the on-disk library index refreshed after all packs finish.
        let tail = stdout
            .lines()
            .rev()
            .find(|l| !l.trim().is_empty())
            .map(|l| l.trim().to_string());
        if let Some(line) = tail {
            summary_lines.push(format!("  {}", line));
        } else {
            summary_lines.push(format!("✓ {}: completed", script_name));
        }
    }

    if summary_lines.is_empty() {
        return Ok("No packs selected or all failed.".into());
    }

    crate::io::refresh_library_index();
    let after = crate::io::library_index().len();
    let added = after.saturating_sub(before);

    let mut out = String::new();
    out.push_str(&format!(
        "Library updated: +{} new spectra ({} total indexed).\n",
        added, after
    ));
    out.push_str(&summary_lines.join("\n"));
    Ok(out)
}

/// One live update from `python/library_manager.py`.
#[derive(Debug, Clone)]
pub enum DownloadEvent {
    Progress {
        pack_idx: usize,
        total_packs: usize,
        pack_id: String,
        pct: f64,
        detail: String,
    },
    Finished(std::result::Result<String, String>),
}

/// Parse `PROGRESS|idx|total|pack_id|pct|detail` into its components.
pub fn parse_download_progress_line(line: &str) -> Option<(usize, usize, String, f64, String)> {
    let rest = line.trim().strip_prefix("PROGRESS|")?;
    let mut parts = rest.splitn(5, '|');
    let pack_idx = parts.next()?.trim().parse::<usize>().ok()?;
    let total_packs = parts.next()?.trim().parse::<usize>().ok()?;
    let pack_id = parts.next()?.trim().to_string();
    let pct = parts.next()?.trim().parse::<f64>().ok()?;
    let detail = parts.next().unwrap_or("").trim().to_string();
    Some((pack_idx, total_packs, pack_id, pct, detail))
}

/// Spawn `python/library_manager.py` and stream its machine-readable progress
/// lines to the TUI without blocking the main thread.
pub fn spawn_library_download_stream(
    pack_ids: Vec<String>,
) -> std::sync::mpsc::Receiver<DownloadEvent> {
    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let result = (|| -> Result<String> {
            let (python, script) = find_python_script("library_manager.py")?;
            let mut child = Command::new(&python)
                .arg("-u")
                .arg(&script)
                .arg("--packs")
                .arg(pack_ids.join(","))
                .arg("--lib-root")
                .arg(crate::io::library_base_dir())
                .stdout(Stdio::piped())
                .stderr(Stdio::piped())
                .spawn()
                .with_context(|| {
                    format!("Failed executing {} {}", python.display(), script.display())
                })?;

            let stdout = child
                .stdout
                .take()
                .ok_or_else(|| anyhow!("library_manager.py produced no stdout"))?;
            let stderr = child.stderr.take();
            let stderr_handle = std::thread::spawn(move || {
                let mut buf = String::new();
                if let Some(mut err) = stderr {
                    let _ = err.read_to_string(&mut buf);
                }
                buf
            });

            let mut summary: Option<String> = None;
            let mut errors: Vec<String> = Vec::new();
            for line in BufReader::new(stdout).lines() {
                let Ok(line) = line else { break };
                let trimmed = line.trim();
                if trimmed.is_empty() {
                    continue;
                }
                if let Some((pack_idx, total_packs, pack_id, pct, detail)) =
                    parse_download_progress_line(trimmed)
                {
                    let _ = tx.send(DownloadEvent::Progress {
                        pack_idx,
                        total_packs,
                        pack_id,
                        pct,
                        detail,
                    });
                } else if let Some(rest) = trimmed.strip_prefix("DONE|") {
                    let mut parts = rest.splitn(2, '|');
                    let _added = parts.next().unwrap_or("0");
                    summary = Some(
                        parts
                            .next()
                            .unwrap_or("Library download complete")
                            .to_string(),
                    );
                } else if let Some(rest) = trimmed.strip_prefix("ERROR|") {
                    let mut parts = rest.splitn(2, '|');
                    let _pack = parts.next().unwrap_or("");
                    errors.push(parts.next().unwrap_or("unknown error").to_string());
                }
            }

            let status = child
                .wait()
                .context("Failed waiting for library_manager.py")?;
            let stderr_text = stderr_handle.join().unwrap_or_default();
            crate::io::refresh_library_index();

            if !status.success() {
                let msg = if !errors.is_empty() {
                    errors.join("; ")
                } else if !stderr_text.trim().is_empty() {
                    stderr_text.trim().to_string()
                } else {
                    format!("library_manager.py exited with {}", status)
                };
                return Err(anyhow!(msg));
            }
            if !errors.is_empty() {
                return Err(anyhow!(errors.join("; ")));
            }
            Ok(summary.unwrap_or_else(|| "Library download complete".into()))
        })();

        let _ = tx.send(DownloadEvent::Finished(result.map_err(|e| e.to_string())));
    });
    rx
}

#[allow(dead_code)]
pub fn run_python_digitize(
    img_path: &Path,
    kind: crate::data::SpectrumKind,
    y_unit: crate::data::YUnit,
    x_left: f64,
    x_right: f64,
    bg_mode: &str,
) -> Result<crate::data::Spectrum> {
    let (python, script) = find_python_script("digitize_image.py")?;
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    let tmp_csv = env::temp_dir().join(format!(
        "spectra_digitize_{}_{}.csv",
        std::process::id(),
        stamp
    ));

    let out = Command::new(&python)
        .arg(&script)
        .arg("--input")
        .arg(img_path)
        .arg("--output")
        .arg(&tmp_csv)
        .arg("--kind")
        .arg(kind.as_str())
        .arg("--y-unit")
        .arg(y_unit.as_str())
        .arg("--x-left")
        .arg(format!("{}", x_left))
        .arg("--x-right")
        .arg(format!("{}", x_right))
        .arg("--bg")
        .arg(bg_mode)
        .output()
        .with_context(|| format!("Failed executing {} {}", python.display(), script.display()))?;

    if !out.status.success() {
        let _ = fs::remove_file(&tmp_csv);
        let stderr = String::from_utf8_lossy(&out.stderr);
        return Err(anyhow!("Image digitizer failed: {}", stderr.trim()));
    }

    let mut spec = crate::io::load_spectrum_file(&tmp_csv, Some(kind))?;
    let _ = fs::remove_file(&tmp_csv);

    let stem = img_path
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("spectrum");
    spec.raw_y_unit = y_unit;
    spec.name = format!("{} [Image]", stem);
    spec.source_path = Some(img_path.display().to_string());
    Ok(spec)
}
