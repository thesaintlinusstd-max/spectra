# SPECTRA v0.7.2

**Offline Terminal Workbench for FTIR, ¹H NMR, ¹³C NMR & UV-Vis Spectroscopy and Polymer Chemistry**

Built in **Rust** (`ratatui` + `crossterm`) with optional **Python** (`matplotlib` / `scipy` / `numpy`) tools for publication figures, printed-spectrum digitization, and open-library fetching.

---

## What's New in v0.7.0

- **Native binary instrument parsers** — pure Rust, zero new dependencies: Galactic/Thermo `.spc` (SPCHDR, explicit or linear X, float32/int32/scaled Y), Bruker **OPUS** `.0`/`.1`/`.0000`/`.opus` (block directory + NPF/FXV/LXV/CSF tags), and Thermo Nicolet OMNIC **`.spa`** (key-2 header / key-3 data). All are auto-detected by magic bytes.
- **Peak Deconvolution Studio (`P`)** — Pseudo-Voigt / Gaussian / Lorentzian multi-peak Levenberg–Marquardt fit inside the zoom window, with optional baseline and overlay of components + residual.
- **Deconvolution in publication figures** — dashed sub-peak component curves and a cumulative fit are written into the native SVG (`deconv-components` / `deconv-fit` groups) and Matplotlib PNG/PDF exports.
- **Save Processed Data & Report (`W`)** — export the processed trace as CSV (with a processing provenance header), as JCAMP-DX, or write a full Markdown **analytical report** (metadata, pipeline, peaks, deconvolution fit, integration regions).
- **Kramers–Kronig ATR correction & Kubelka–Munk (`A`)** — diamond / ZnSe / Ge / Si penetration-depth correction plus diffuse-reflectance transforms.
- **arPLS baselines, SNV normalization, derivatives** — `b` cycles OFF → Rubberband → AsLS → arPLS; `w` cycles OFF → MinMax → SNV → Area; `t` cycles Raw → 1st → 2nd derivative.
- **Interactive Spectral Math (`o`)** — weighted subtract, add, ratio and average with automatic X-grid alignment.
- **Multi-component NNLS mixture analysis (`M` then `m`)** — Lawson–Hanson non-negative least squares over 1–3 references with a composite HQI.
- **Quantitative Integration Studio (`I`)** — Mnova-style regions with relative integrals and one-key push into Polymer Lab.
- **X-Axis Calibration (`X`)** — snap an observed peak to TMS, solvent residual, polystyrene 1601.4 cm⁻¹, or holmium oxide.
- **Headless CLI pipelines** — `spectra deconv`, `spectra mix`, `spectra convert`.

---

## Screenshots

| Feature | Screenshot |
|---------|------------|
| FTIR HD Oscilloscope (Kitty) | ![FTIR HD](assets/ftir_hd.png) |
| UV-Vis HD (with Visible Rainbow Strip) | ![UV-Vis HD](assets/uvvis_hd.png) |
| NMR HD (`CDCl₃` & Impurity Markers) | ![NMR HD](assets/nmr_hd.png) |
| Publication Export (600 DPI ACS Style) | ![Export ACS](assets/export_acs.png) |

---

## Beginner Quick-Start (5 Minutes)

1. **Look at something real** — run `spectra` and press **`d`**. Nine labelled demo spectra (FTIR, ¹H NMR, ¹³C NMR, UV-Vis) open instantly; no files required. Move the cursor with `←`/`→` and read the Inspector to learn each spectral zone.
2. **Open your own data** — press **`f`** and pick a file. `.spc`, OPUS `.0`, `.spa`, JCAMP-DX (`.jdx`/`.dx`), CSV, TSV and TXT are all recognised automatically.
3. **Clean it up** — `b` baseline (OFF → Rubberband → AsLS → arPLS), `m` Savitzky–Golay smoothing, `A` ATR/Kubelka–Munk optics, `t` derivatives, `w` normalization.
4. **Assign and identify** — `p` toggles peak labels, `x` pins a peak, `c` copies the assignment, `M` runs the HQI library matcher (press `m` inside it for NNLS mixture analysis).
5. **Dig into overlaps** — drag a box around a band to zoom, press **`P`**, then **`Enter`** to fit 1–8 sub-peaks (and **`r`** to toggle the residual). The result stays overlaid on the plot.
6. **Quantify** — `I` opens the Integration Studio; `X` calibrates the X axis; Polymer Lab (tab `6`) turns ratios into Mn/DPn, tacticity and phase indices.
7. **Save or publish** — **`W`** writes CSV/JCAMP-DX or a Markdown analytical report into `exports/`; **`e`** renders a 600 DPI PNG, vector SVG or PDF.

> Processing is **non-destructive**: nothing is written to disk until you use `W` (save files/report) or `X` (calibration).

---

## Core Features

### 6 Interactive Tabs
| Key | Tab | Description |
|-----|-----|-------------|
| `1` | **FTIR** | `4000 → 400 cm⁻¹` (IUPAC reversed), `%T` / `Abs`, 122 functional-group & filler bands with confirmation peaks |
| `2` | **¹H NMR** | `14.0 → -0.5 ppm` (IUPAC reversed), J-coupling notes, residual solvent & water/grease/stabilizer alerts |
| `3` | **¹³C NMR** | `220 → 0 ppm` (IUPAC reversed), residual ¹³C solvent markers & polymer tacticity/backbone shifts |
| `4` | **UV-Vis** | **Ascending `190 → 900 nm`**, live photon energy (`eV`), complementary solution color, 67 chromophore/NP bands |
| `5` | **Solvents** | Interactive ¹H / ¹³C NMR trace-impurity & residual-solvent explorer across 10 deuterated solvents |
| `6` | **Polymer Lab** | 5 specialized calculators for polymer & materials characterization |

### Dual Terminal Rendering (`v` toggles)
- **HD Kitty Oscilloscope** — Anti-aliased continuous curves on a pitch-black oscilloscope canvas via the Kitty Graphics Protocol, with non-overlapping bitmap peak labels and a **380–750 nm visible-light rainbow strip** in UV-Vis mode.
- **Portable Braille Canvas** — Automatic fallback in any 24-bit terminal (`alacritty`, `foot`, `tmux`, `konsole`, VS Code, SSH).

### Instrument Parsers & Image Digitizer (`f` / `F`)
- **Native binary formats (pure Rust, zero dependencies)**:
  - **Galactic / Thermo `.spc`** — new little-endian SPCHDR, explicit (`TXVALS`) or evenly spaced X, IEEE float32 / int32 / scaled 16-bit Y, µm→cm⁻¹ and µm→nm conversion.
  - **Bruker OPUS** `.0`, `.1`, `.0000`, `.opus` — block directory with parameter tags (`NPF`, `FXV`, `LXV`, `CSF`, `DXU`), AB/TR/KM/RFC spectra preferred over interferograms.
  - **Thermo Nicolet OMNIC `.spa`** — block directory, key-2 header (points, X/Y units, first/last X) and key-3 float data.
  - Detection is by magic bytes, so the file extension never has to be right.
- **Delimited exports**: `.csv`, `.tsv`, `.txt`, `.dpt`, `.asc`, `.xy` — handles multi-line vendor headers, European decimal commas (`,` → `.`), and Latin-1 / Windows-1252 encodings.
- **IUPAC JCAMP-DX** (`.jdx`, `.dx`, `.jcamp`) — Full decompression of `AFFN`, `PAC`, `SQZ`, `DIF`, `DUP`, and `DIFDUP` (with DIF checkpoint handling).
- **Offline Image Digitizer** — Load a `.png`/`.jpg` spectrum scan directly from the file browser (`f`), calibrate X-axis bounds and background, and trace the curve into an interactive spectrum.

### Signal Processing, Optics & Analysis Studios
| Key | Studio / Operation |
|-----|--------------------|
| `b` | Baseline: **OFF → Rubberband → AsLS → arPLS** |
| `m` | **Savitzky–Golay smoothing** (`OFF → SG-5 → SG-9 → SG-15`) |
| `w` | Normalization: **OFF → Min–Max → SNV → Area** |
| `t` | Derivative spectroscopy: **Raw → 1st → 2nd** |
| `A` | **ATR Kramers–Kronig** (diamond / ZnSe / Ge / Si) and **Kubelka–Munk** |
| `C` | Atmospheric **CO₂ suppression** (`2280–2410 cm⁻¹`) |
| `P` | **Peak Deconvolution Studio** — Pseudo-Voigt/Gaussian/Lorentzian LM fit |
| `o` | **Interactive Spectral Math** — subtract / add / ratio / average |
| `I` | **Quantitative Integration Studio** — relative integrals → Polymer Lab |
| `X` | **X-Axis Calibration** — TMS, solvent residual, polystyrene 1601.4 cm⁻¹, holmium oxide |
| `M` → `m` | **Library matcher** ↔ **NNLS mixture analysis** (1–3 components) |

### Spectral Library Matcher (`Shift+M`)
- **Hit Quality Index (HQI)** combining absorbance Pearson correlation and 1st-derivative correlation (robust against baseline offset and film thickness).
- **FTIR**: Ranks 39 built-in polymer fingerprints + all spectra in `library/ftir/` (`520–3760 cm⁻¹` grid, atmospheric `CO₂` `2280–2410 cm⁻¹` masked).
- **UV-Vis**: Ranks built-in demos + all spectra in `library/uvvis/` (`200–800 nm` grid).
- Displays top 25 scrollable hits; press `Enter` to overlay a matched reference onto your sample.
- **Mixture mode (`m` inside the matcher)** — Lawson–Hanson **NNLS** explains the sample as a non-negative sum of 1–3 references, ranked by a composite HQI; `Enter` overlays every component and `c` copies the Markdown report.

### Save Processed Data & Analytical Reports (`Shift+W`)
- **Processed trace → CSV** with a commented provenance header (trace, technique, Y-axis, full processing pipeline) and instrument-conventional X ordering.
- **Processed trace → JCAMP-DX** (`.jdx`) — compact `DIFDUP` for FTIR `%T`, standard `AFFN` `XYPOINTS` for everything else.
- **Full Markdown analytical report** — trace metadata, processing pipeline, detected/pinned peaks with prominence and FWHM, the deconvolution fit table, and quantitative integration regions.
- Files are written to `exports/`; `c` copies the Markdown report straight to the clipboard.

### Polymer & Materials Chemistry Lab (Tab 6)
1. **FTIR Polymerization Conversion (%)** — Automated or manual `C=C` consumption vs internal carbonyl reference (`%T → Abs` applied automatically).
2. **¹H NMR End-Group Analysis** — Computes *DPₙ* and *Mₙ* (`g/mol` & `kDa`) with built-in monomer presets (`MMA`, `Styrene`, `EO`, `LA`, `ε-CL`, `VAc`, `NIPAM`, `Acrylamide`).
3. **Vinyl / PMMA Tacticity** — *mm / mr / rr* triads, dyad fractions (*m*, *r*), and Bernoullian persistence ratio (*ρ*).
4. **PVDF Phase Fraction & PU H-Bonding** — Gregorio *F(β)* (`840 / 763 cm⁻¹`) and Polyurethane carbonyl H-bonding index (`1703 / 1732 cm⁻¹`).
5. **UV-Vis Tauc Band Gap & Beer–Lambert** — Direct `(αhν)²` / Indirect `(αhν)¹/²` Tauc onset *E_g* (`eV`) with regression *R²* (`t` toggles), plus auto-scaled concentration (`M` / `mM` / `µM` / `nM`).

### Publication-Grade Export (`e`)
- **600 DPI PNG / Vector SVG / PDF** — SVG is rendered natively in Rust; PNG/PDF use the Python Matplotlib helper.
- **Styles**: `ACS`, `RSC`, and `Dark Everforest`.
- **Layout**: `Overlay` or `Stacked` (with non-overlapping leader-line peak labels and visible-light rainbow strip for UV-Vis).
- **Deconvolution overlays**: dashed sub-peak component curves plus a cumulative fit are included automatically whenever a deconvolution result is displayed.
- Copy exported PNG directly to Wayland clipboard (`y` via `wl-copy`) or open in viewer (`o`).

---

## Installation & Quick Start

### Prerequisites
- **Rust** (stable, edition 2021) — [rustup.rs](https://rustup.rs)
- **Python 3.10+** with `numpy`, `scipy`, `matplotlib` *(only needed for PNG/PDF publication figures and `fetch`; SVG export and printed-image digitizing are pure Rust)*
  - Arch / CachyOS: `sudo pacman -S --needed base-devel wl-clipboard python python-numpy python-scipy python-matplotlib`
  - Debian / Ubuntu: `pip install numpy scipy matplotlib pillow`

### Build & Install
```bash
git clone https://github.com/thesaintlinusstd-max/spectra.git
cd spectra
cargo build --release
mkdir -p ~/.local/bin
ln -sf "$PWD/target/release/spectra" ~/.local/bin/spectra
```
*(After the first build, `cargo build --release --offline` and `cargo test --release --offline` work completely offline.)*

### Quick Usage
```bash
# Launch interactive TUI (opens built-in reference workspace immediately)
spectra

# Open your own instrument or image files
spectra monomer.csv polymer.csv
spectra sample.jdx
spectra sample.spc                # native binary formats auto-detected
spectra background.0              # Bruker OPUS
spectra film.spa                  # Thermo Nicolet OMNIC

# Library lookup by compound name & mode
spectra butane               # opens interactive picker across matches
spectra ftir butane          # strict FTIR search
spectra cnmr butane          # strict ¹³C NMR search
spectra uvvis coumarin       # strict UV-Vis search

# Generate sample datasets
spectra init-samples

# Headless publication-grade export
spectra export samples/mma_monomer_ftir.csv samples/pmma_polymer_ftir.csv \
  -o exports/figure1.png --style acs --stacked --title "MMA vs PMMA"
```

### Headless Analysis Pipelines

| Command | Purpose |
|---------|---------|
| `spectra deconv FILE [--range LO:HI] [--shape voigt\|gauss\|lorentz] [--peaks N] [--no-baseline] [-o fit.csv]` | Deconvolute overlapping bands in a window and print a Markdown fit table (optionally write the sub-peak CSV). |
| `spectra mix FILE [--components 1\|2\|3] [--solvent cdcl3\|dmso\|d2o\|…]` | NNLS multi-component mixture analysis against the reference library. |
| `spectra convert INPUT -o OUT.csv\|OUT.jdx [--atr diamond\|znse\|ge\|si] [--km] [--baseline\|--rubberband\|--arpls] [--smooth 0\|5\|9\|15] [--absorbance]` | Apply the full processing pipeline and write a processed CSV or JCAMP-DX. |
| `spectra export FILES… -o OUT.svg\|png\|pdf [--style acs\|rsc\|everforest] [--stacked] [--absorbance] [--baseline\|--rubberband\|--arpls] [--atr …\|--km] [--deriv 1\|2] [--norm minmax\|snv\|area] [--smooth …] [--dpi 600]` | Headless publication figure (SVG is native Rust). |
| `spectra fetch …` | Download optional open reference libraries. |
| `spectra init-samples` | Write example CSV spectra into `samples/`. |

```bash
# Deconvolute the PMMA carbonyl region and save the components
spectra deconv samples/pmma_polymer_ftir.csv --range 1680:1780 --peaks 3 -o exports/carbonyl.csv

# Identify a blend, then export the processed trace as JCAMP-DX
spectra mix polymer_blend.csv --components 2
spectra convert polymer_blend.csv -o exports/blend_processed.jdx --arpls --smooth 9
```

---

## Fetching Open Reference Libraries (`Shift+D` or CLI)

The core binary embeds 280+ correlation bands, 63 solvent impurity entries, 39 polymer fingerprints, and 9 demo presets with zero downloads required.

To download optional experimental reference spectra into `library/`:
- **Inside the TUI**: Press **`Shift+D`** to open the interactive **Reference Library Downloader** and toggle the packs you want.
- **From CLI**:

| Command | Source | Content | Approx. Download |
|---------|--------|---------|------------------|
| `spectra fetch pack --nist` | NIST Chemistry WebBook | ~38 polymer monomers & solvents (FTIR `.jdx`) | ~0.4 MB |
| `spectra fetch pack --uvvis` | OMLC PhotochemCAD + NIST WebBook | ~187 experimental UV-Vis spectra | ~2 MB |
| `spectra fetch pack --zenodo` | Zenodo FTIR-Plastics | 36 Jasco ATR-FTIR spectra (`PET, HDPE, LDPE, PP, PS, PVC`) | ~253 MB zip |
| `spectra fetch pack --nmr` | nmrshiftdb2 SDF dump | 60 monomer/solvent NMR spectra + 58k offline search cache | ~159 MB SDF |
| `python3 python/download_mega_ir.py --all` | RRUFF + INFRA-ART + MNHN + NIST | ~1,800+ mineral, pigment, resin & organic FTIR spectra | ~70 MB |
| `spectra fetch nmr-search <name>` | Local cached nmrshiftdb2 SDF | Extract ¹H/¹³C shifts for any of 58,183 compounds | offline after cache |
| `spectra fetch uvvis <CAS>` | NIST WebBook (`Type=UVVis`) | Fetch a single UV-Vis JCAMP-DX by CAS | <50 KB |
| `spectra fetch import <dir>` | Local CSV/JCAMP folder | Index external packs (e.g., FLOPP / Open Specy) | local |

**Open Data Acknowledgments:**
- **NIST Chemistry WebBook** (U.S. Secretary of Commerce) — experimental IR & UV-Vis JCAMP-DX data
- **OMLC PhotochemCAD** (Oregon Medical Laser Center) — measured UV-Vis absorption spectra
- **Zenodo FTIR-Plastics** (`10.5281/zenodo.10736649`), **INFRA-ART** (`10.5281/zenodo.21993883`), and **MNHN Paris ATR-FTIR** (`10.5281/zenodo.10658337`)
- **RRUFF Project** (University of Arizona) — mineral & inorganic infrared spectra
- **nmrshiftdb2** (University of Cologne) — open peer-reviewed NMR chemical shift database

---

## Keyboard & Mouse Shortcuts

### Global & Navigation
| Key / Mouse | Action |
|-------------|--------|
| `1`–`6` | Switch tabs (`FTIR` \| `¹H NMR` \| `¹³C NMR` \| `UV-Vis` \| `Solvents` \| `Polymer Lab`) |
| `←`/`→` or `h`/`l` | Move cursor (`Shift` = coarse jump) |
| Drag box / Scroll / Right-click | Zoom window / zoom at cursor / reset zoom (`0` or `r`) |
| `g` | Jump to exact `cm⁻¹`, `ppm`, or `nm` |
| `Tab` / `Space` / `Del` | Switch active trace / hide trace / remove active trace |
| `s` / `v` | Toggle Stacked vs Overlay / toggle HD Kitty vs Braille |
| `f` / `F` | Open file browser (replace / append overlay) |
| `/` / `d` | Global search / built-in demo presets |
| `a` / `b` / `m` / `w` / `t` | `%T↔Abs` / baseline / smoothing / normalization / derivative |
| `A` / `C` | ATR Kramers–Kronig or Kubelka–Munk / atmospheric CO₂ suppression |
| `P` / `o` / `I` / `X` | Deconvolution / Spectral Math / Integration / X-axis calibration studios |
| `Shift+M` / `Shift+D` | Spectral library matcher (HQI + NNLS mixture) / interactive library downloader |
| `Shift+W` / `e` / `?` / `q` | Save processed data & report / export figure / user guide / quit |

### NMR & Polymer Lab
| Key | Action |
|-----|--------|
| `Shift+S` / `i` | Cycle deuterated solvent / toggle NMR solvent & water markers |
| `h` / `l` or `Tab` (Tab 6) | Switch between the 5 Polymer Lab calculators |
| `j` / `k` / `Enter` (Tab 6) | Select and edit input field |
| `m` / `t` / `p` / `c` (Tab 6) | Toggle Auto/Manual (`m`), Direct/Indirect Tauc (`t`), Monomer preset (`p`), Copy (`c`) |

---

## Scientific References & Disclaimer

### References
- **Fulmer, G. R. et al.** *Organometallics* **2010**, *29*, 2176–2179 (¹H/¹³C NMR trace impurities).
- **Gottlieb, H. E.; Kotlyar, V.; Nudelman, A.** *J. Org. Chem.* **1997**, *62*, 7512–7515.
- **Babij, N. R. et al.** *Org. Process Res. Dev.* **2016**, *20*, 661–667.
- **Gregorio, R., Jr.; Cestari, M.** *J. Polym. Sci., Part B: Polym. Phys.* **1994**, *32*, 859–870 (PVDF β-phase quantification).
- **Tauc, J.; Grigorovici, R.; Vancu, A.** *Phys. Status Solidi (b)* **1966**, *15*, 627–637 (optical band gap).
- **Socrates, G.** *Infrared and Raman Characteristic Group Frequencies: Tables and Charts*, 3rd ed.; Wiley: Chichester, **2001**.

### Scope & Scientific Limitations
- **Built-in `[Ref]` polymer spectra are synthetic educational references** generated from literature peak positions rather than raw measured scans; use `library/` experimental spectra (or your own instrument exports) for definitive comparison.
- **Image digitization** from PNG/JPG scans is intended for qualitative peak inspection and quick study; **always use raw numeric instrument exports (`.spc`, `.spa`, OPUS, `.csv`/`.jdx`/`.dpt`) for quantitative conversion, *Mₙ*, or band-gap analysis**.
- **HQI library matching** provides ranked similarity candidates, not automated proof of structure — verify assignments using orthogonal characterization.

---

## License

**MIT License** — Copyright (c) 2026 Muhammad Musavi. See [LICENSE](LICENSE).
