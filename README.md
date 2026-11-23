# SPECTRA v0.5.0

**Offline Terminal Workbench for FTIR, ¹H NMR, ¹³C NMR & UV-Vis Spectroscopy and Polymer Chemistry**

Built in **Rust** (`ratatui` + `crossterm`) with optional **Python** (`matplotlib` / `scipy` / `numpy`) tools for publication figures, printed-spectrum digitization, and open-library fetching.

---

## Screenshots

| Feature | Screenshot |
|---------|------------|
| FTIR HD Oscilloscope (Kitty) | ![FTIR HD](assets/ftir_hd.png) |
| UV-Vis HD (with Visible Rainbow Strip) | ![UV-Vis HD](assets/uvvis_hd.png) |
| NMR HD (`CDCl₃` & Impurity Markers) | ![NMR HD](assets/nmr_hd.png) |
| Publication Export (600 DPI ACS Style) | ![Export ACS](assets/export_acs.png) |

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
- **Delimited exports**: `.csv`, `.tsv`, `.txt`, `.dpt`, `.asc`, `.xy` — handles multi-line vendor headers, European decimal commas (`,` → `.`), and Latin-1 / Windows-1252 encodings.
- **IUPAC JCAMP-DX** (`.jdx`, `.dx`, `.jcamp`) — Full decompression of `AFFN`, `PAC`, `SQZ`, `DIF`, `DUP`, and `DIFDUP` (with DIF checkpoint handling).
- **Offline Image Digitizer** — Load a `.png`/`.jpg` spectrum scan directly from the file browser (`f`), calibrate X-axis bounds and background, and trace the curve into an interactive spectrum.

### Pure-Rust Signal Processing (No SciPy Required)
| Key | Operation |
|-----|-----------|
| `a` | Beer–Lambert **`%Transmittance ↔️ Absorbance`** toggle (FTIR & UV-Vis) |
| `b` | **Rubberband baseline** correction (Bruker OPUS-style convex hull) |
| `m` | **Savitzky–Golay smoothing** (`OFF → SG-5 → SG-9 → SG-15`) |
| `w` | **Min–Max normalization** for trace comparison |
| `p` | **Prominence peak detection** with 3-point parabolic refinement |
| `n` / `N` | Jump to next / previous detected peak |
| `x` | Pin / unpin (`★`) a peak at the cursor |
| `[` / `]` | Lower / raise peak prominence threshold (`2–60%`) |

### Spectral Library Matcher (`Shift+M`)
- **Hit Quality Index (HQI)** combining absorbance Pearson correlation and 1st-derivative correlation (robust against baseline offset and film thickness).
- **FTIR**: Ranks 39 built-in polymer fingerprints + all spectra in `library/ftir/` (`520–3760 cm⁻¹` grid, atmospheric `CO₂` `2280–2410 cm⁻¹` masked).
- **UV-Vis**: Ranks built-in demos + all spectra in `library/uvvis/` (`200–800 nm` grid).
- Displays top 25 scrollable hits; press `Enter` to overlay a matched reference onto your sample.

### Polymer & Materials Chemistry Lab (Tab 6)
1. **FTIR Polymerization Conversion (%)** — Automated or manual `C=C` consumption vs internal carbonyl reference (`%T → Abs` applied automatically).
2. **¹H NMR End-Group Analysis** — Computes *DPₙ* and *Mₙ* (`g/mol` & `kDa`) with built-in monomer presets (`MMA`, `Styrene`, `EO`, `LA`, `ε-CL`, `VAc`, `NIPAM`, `Acrylamide`).
3. **Vinyl / PMMA Tacticity** — *mm / mr / rr* triads, dyad fractions (*m*, *r*), and Bernoullian persistence ratio (*ρ*).
4. **PVDF Phase Fraction & PU H-Bonding** — Gregorio *F(β)* (`840 / 763 cm⁻¹`) and Polyurethane carbonyl H-bonding index (`1703 / 1732 cm⁻¹`).
5. **UV-Vis Tauc Band Gap & Beer–Lambert** — Direct `(αhν)²` / Indirect `(αhν)¹/²` Tauc onset *E_g* (`eV`) with regression *R²* (`t` toggles), plus auto-scaled concentration (`M` / `mM` / `µM` / `nM`).

### Publication-Grade Export (`e`)
- **600 DPI PNG / Vector SVG / PDF** rendered via Python Matplotlib.
- **Styles**: `ACS`, `RSC`, and `Dark Everforest`.
- **Layout**: `Overlay` or `Stacked` (with non-overlapping leader-line peak labels and visible-light rainbow strip for UV-Vis).
- Copy exported PNG directly to Wayland clipboard (`y` via `wl-copy`) or open in viewer (`o`).

---

## Installation & Quick Start

### Prerequisites
- **Rust** (stable, edition 2021) — [rustup.rs](https://rustup.rs)
- **Python 3.10+** with `numpy`, `scipy`, `matplotlib` *(only needed for `[e]` export, image digitization, and `fetch`)*
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
| `Shift+M` / `Shift+D` | Spectral library matcher (HQI) / interactive library downloader |
| `e` / `?` / `q` | Export Matplotlib figure / user guide / quit |

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
- **Image digitization** from PNG/JPG scans is intended for qualitative peak inspection and quick study; **always use raw numeric instrument exports (`.csv`/`.jdx`/`.dpt`) for quantitative conversion, *Mₙ*, or band-gap analysis**.
- **HQI library matching** provides ranked similarity candidates, not automated proof of structure — verify assignments using orthogonal characterization.

---

## License

**MIT License** — Copyright (c) 2026 Muhammad Musavi. See [LICENSE](LICENSE).
