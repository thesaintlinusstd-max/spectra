use std::env;
use std::fmt::Write as _;
use std::fs;
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
        out.push_str("]}");
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

/// Run selected library download packs and return a summary string.
/// Items are tuples of (script_name, args_as_single_string).
/// Returns a clean summary with total files downloaded per pack.
pub fn run_library_downloads(items: &[(&str, &str)]) -> Result<String> {
    use std::process::Command;

    let mut summary_lines: Vec<String> = Vec::new();
    let mut total_files = 0;

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

        // Extract "Total ... downloaded: N" or similar from stdout
        let mut pack_total = 0;
        for line in stdout.lines() {
            let line = line.trim();
            if line.starts_with("Total ") && line.contains("downloaded:") {
                if let Some(num_str) = line.split("downloaded:").nth(1) {
                    if let Ok(n) = num_str.trim().parse::<usize>() {
                        pack_total = n;
                    }
                }
            } else if line.contains("summary:") && line.contains("downloaded") {
                // PhotochemCAD/NIST UV-Vis summary lines
                if let Some(idx) = line.find("downloaded") {
                    let before = &line[..idx];
                    if let Some(num_str) = before.split_whitespace().last() {
                        if let Ok(n) = num_str.parse::<usize>() {
                            pack_total = n;
                        }
                    }
                }
            }
        }
        if pack_total > 0 {
            summary_lines.push(format!("✓ {}: {} spectra", script_name, pack_total));
            total_files += pack_total;
        } else if !stdout.trim().is_empty() {
            // Fallback: show last few lines
            let last_lines: Vec<&str> = stdout.lines().rev().take(3).collect();
            for l in last_lines.iter().rev() {
                if !l.trim().is_empty() {
                    summary_lines.push(format!("  {}", l.trim()));
                }
            }
        } else {
            summary_lines.push(format!("✓ {}: completed", script_name));
        }
    }

    if summary_lines.is_empty() {
        return Ok("No packs selected or all failed.".into());
    }

    let mut out = String::new();
    out.push_str(&format!("Downloaded {} reference spectra total.\n", total_files));
    out.push_str(&summary_lines.join("\n"));
    Ok(out)
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
