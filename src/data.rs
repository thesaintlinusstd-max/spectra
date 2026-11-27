#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpectrumKind {
    Ftir,
    Hnmr,
    Cnmr,
    UvVis,
}

impl SpectrumKind {
    pub fn x_unit(self) -> &'static str {
        match self {
            Self::Ftir => "cm⁻¹",
            Self::Hnmr | Self::Cnmr => "ppm",
            Self::UvVis => "nm",
        }
    }

    pub fn symbol(self) -> &'static str {
        match self {
            Self::Ftir => "ν",
            Self::Hnmr | Self::Cnmr => "δ",
            Self::UvVis => "λ",
        }
    }

    pub fn default_bounds(self) -> (f64, f64) {
        match self {
            Self::Ftir => (400.0, 4000.0),
            Self::Hnmr => (-0.5, 13.0),
            Self::Cnmr => (0.0, 220.0),
            Self::UvVis => (190.0, 900.0),
        }
    }

    pub fn display_name(self) -> &'static str {
        match self {
            Self::Ftir => "FTIR",
            Self::Hnmr => "¹H NMR",
            Self::Cnmr => "¹³C NMR",
            Self::UvVis => "UV-Vis",
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Ftir => "ftir",
            Self::Hnmr => "hnmr",
            Self::Cnmr => "cnmr",
            Self::UvVis => "uvvis",
        }
    }

    pub fn is_x_reversed(self) -> bool {
        !matches!(self, Self::UvVis)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum YUnit {
    Transmittance,
    Absorbance,
    Intensity,
}

impl YUnit {
    pub fn label(self) -> &'static str {
        match self {
            Self::Transmittance => "%T",
            Self::Absorbance => "Abs",
            Self::Intensity => "Int",
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Transmittance => "transmittance",
            Self::Absorbance => "absorbance",
            Self::Intensity => "intensity",
        }
    }
}

#[derive(Debug, Clone)]
pub struct Spectrum {
    pub name: String,
    pub source_path: Option<String>,
    pub kind: SpectrumKind,
    pub raw_y_unit: YUnit,
    /// Strictly ascending X values
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub visible: bool,
    pub color_idx: usize,
    pub metadata: Vec<(String, String)>,
}

impl Spectrum {
    pub fn x_range(&self) -> (f64, f64) {
        match (self.x.first(), self.x.last()) {
            (Some(&lo), Some(&hi)) if lo < hi => (lo, hi),
            _ => self.kind.default_bounds(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct FtirBand {
    pub min: f64,
    pub max: f64,
    pub group: &'static str,
    pub mode: &'static str,
    pub intensity: &'static str,
    pub note: &'static str,
    pub confirm: &'static str,
}

#[derive(Debug, Clone, Copy)]
pub struct NmrBand {
    pub min: f64,
    pub max: f64,
    pub group: &'static str,
    pub example: &'static str,
    pub note: &'static str,
    pub j_hz: &'static str,
}

#[derive(Debug, Clone, Copy)]
pub struct UvVisBand {
    pub min: f64,
    pub max: f64,
    pub group: &'static str,
    pub transition: &'static str,
    pub epsilon: &'static str,
    pub note: &'static str,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolventId {
    Cdcl3 = 0,
    Dmso = 1,
    Cd3od = 2,
    D2o = 3,
    Acetone = 4,
    Thf = 5,
    C6d6 = 6,
    Toluene = 7,
    Cd2cl2 = 8,
    Cd3cn = 9,
}

impl SolventId {
    pub const ALL: [SolventId; 10] = [
        SolventId::Cdcl3,
        SolventId::Dmso,
        SolventId::Cd3od,
        SolventId::D2o,
        SolventId::Acetone,
        SolventId::Thf,
        SolventId::C6d6,
        SolventId::Toluene,
        SolventId::Cd2cl2,
        SolventId::Cd3cn,
    ];

    pub fn index(self) -> usize {
        self as usize
    }

    pub fn next(self) -> Self {
        Self::ALL[(self.index() + 1) % Self::ALL.len()]
    }

    pub fn prev(self) -> Self {
        Self::ALL[(self.index() + Self::ALL.len() - 1) % Self::ALL.len()]
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Cdcl3 => "CDCl₃",
            Self::Dmso => "DMSO-d₆",
            Self::Cd3od => "CD₃OD",
            Self::D2o => "D₂O",
            Self::Acetone => "Acetone-d₆",
            Self::Thf => "THF-d₈",
            Self::C6d6 => "C₆D₆",
            Self::Toluene => "Toluene-d₈",
            Self::Cd2cl2 => "CD₂Cl₂",
            Self::Cd3cn => "CD₃CN",
        }
    }

    /// Primary residual solvent ¹H signal(s) (Fulmer 2010)
    pub fn residual_h(self) -> &'static [f64] {
        match self {
            Self::Cdcl3 => &[7.26],
            Self::Dmso => &[2.50],
            Self::Cd3od => &[3.31],
            Self::D2o => &[4.79],
            Self::Acetone => &[2.05],
            Self::Thf => &[1.72, 3.58],
            Self::C6d6 => &[7.16],
            Self::Toluene => &[2.08, 6.97, 7.01, 7.09],
            Self::Cd2cl2 => &[5.32],
            Self::Cd3cn => &[1.94],
        }
    }

    /// Water/HDO ¹H shift in this solvent (Fulmer 2010)
    pub fn water_h(self) -> f64 {
        match self {
            Self::Cdcl3 => 1.56,
            Self::Dmso => 3.33,
            Self::Cd3od => 4.87,
            Self::D2o => 4.79,
            Self::Acetone => 2.84,
            Self::Thf => 2.46,
            Self::C6d6 => 0.40,
            Self::Toluene => 0.43,
            Self::Cd2cl2 => 1.52,
            Self::Cd3cn => 2.13,
        }
    }

    /// Primary residual solvent ¹³C signal(s) (Fulmer 2010)
    pub fn residual_c(self) -> &'static [f64] {
        match self {
            Self::Cdcl3 => &[77.16],
            Self::Dmso => &[39.52],
            Self::Cd3od => &[49.00],
            Self::D2o => &[],
            Self::Acetone => &[29.84, 206.26],
            Self::Thf => &[25.31, 67.21],
            Self::C6d6 => &[128.06],
            Self::Toluene => &[20.43, 125.13, 127.96, 128.87, 137.48],
            Self::Cd2cl2 => &[53.84],
            Self::Cd3cn => &[1.32, 118.26],
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SolventImpurity {
    pub name: &'static str,
    pub group: &'static str,
    pub mult: &'static str,
    pub h_shifts: [Option<f64>; 10],
    pub c_shifts: [Option<f64>; 10],
    pub note: &'static str,
}

#[derive(Debug, Clone, Copy)]
pub struct PolymerRef {
    pub short: &'static str,
    pub name: &'static str,
    pub repeat_unit: &'static str,
    pub m0: f64,
    pub tg_tm: &'static str,
    pub ftir_peaks: &'static [(f64, &'static str)],
    pub hnmr_peaks: &'static [(f64, &'static str)],
    pub cnmr_peaks: &'static [(f64, &'static str)],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DemoPreset {
    MmaVsPmmaFtir,
    PolystyreneFtir,
    PetVsNylonFtir,
    PmmaHnmr,
    PolystyreneHnmr,
    PmmaCnmr,
    PaniUvVis,
    Tio2ZnoNanocompositeUvVis,
    AgAuNanoparticlesUvVis,
}

impl DemoPreset {
    pub const ALL: [DemoPreset; 9] = [
        Self::MmaVsPmmaFtir,
        Self::PolystyreneFtir,
        Self::PetVsNylonFtir,
        Self::PmmaHnmr,
        Self::PolystyreneHnmr,
        Self::PmmaCnmr,
        Self::PaniUvVis,
        Self::Tio2ZnoNanocompositeUvVis,
        Self::AgAuNanoparticlesUvVis,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Self::MmaVsPmmaFtir => "FTIR: MMA Monomer vs PMMA Polymer (Conversion Demo)",
            Self::PolystyreneFtir => "FTIR: Polystyrene (Aromatic & Fingerprint)",
            Self::PetVsNylonFtir => "FTIR: PET Ester vs Nylon-6,6 Amide Overlay",
            Self::PmmaHnmr => "¹H NMR: PMMA in CDCl₃ (Tacticity mm/mr/rr + Impurities)",
            Self::PolystyreneHnmr => "¹H NMR: Polystyrene in CDCl₃",
            Self::PmmaCnmr => "¹³C NMR: PMMA in CDCl₃ (Backbone & Tacticity)",
            Self::PaniUvVis => "UV-Vis: Polyaniline (Emeraldine Base vs H⁺-Doped Salt)",
            Self::Tio2ZnoNanocompositeUvVis => "UV-Vis: TiO₂ / ZnO Nanocomposite (Tauc Band Gap ~3.2–3.3 eV)",
            Self::AgAuNanoparticlesUvVis => "UV-Vis: Ag (412 nm) & Au (524 nm) Nanoparticles LSPR",
        }
    }

    pub fn kind(self) -> SpectrumKind {
        match self {
            Self::MmaVsPmmaFtir | Self::PolystyreneFtir | Self::PetVsNylonFtir => SpectrumKind::Ftir,
            Self::PmmaHnmr | Self::PolystyreneHnmr => SpectrumKind::Hnmr,
            Self::PmmaCnmr => SpectrumKind::Cnmr,
            Self::PaniUvVis
            | Self::Tio2ZnoNanocompositeUvVis
            | Self::AgAuNanoparticlesUvVis => SpectrumKind::UvVis,
        }
    }
}

// ============================================================================
// COMPREHENSIVE FTIR FUNCTIONAL GROUP DATABASE (Organic + Polymer Chemistry)
// ============================================================================
pub const FTIR_DB: &[FtirBand] = &[
    FtirBand { min: 3580.0, max: 3650.0, group: "O–H (free alcohol/phenol)", mode: "ν(O–H) stretch", intensity: "Sharp, Med", note: "Unbonded hydroxyl; epoxy/PHEMA end-groups", confirm: "Confirm sharp 3640; compare H-bonded O–H 3200–3550 cm⁻¹" },
    FtirBand { min: 3200.0, max: 3550.0, group: "O–H (H-bonded alcohol/moisture)", mode: "ν(O–H) stretch", intensity: "Very Broad, Strong", note: "PVA, PEG end-groups, cured epoxy, absorbed H₂O", confirm: "Confirm C–O stretch 1050–1150 cm⁻¹" },
    FtirBand { min: 3300.0, max: 3500.0, group: "N–H (amide / amine / urethane)", mode: "ν(N–H) stretch", intensity: "Med–Strong", note: "Nylon-6,6 (~3300), Polyurethane (~3330), primary amine doublet", confirm: "Confirm Amide I 1630–1685 & Amide II 1515–1570 cm⁻¹" },
    FtirBand {
        min: 1625.0,
        max: 1652.0,
        group: "H₂O (adsorbed moisture bend)",
        mode: "δ(H–O–H) bend",
        intensity: "Med, Broad",
        note: "Adsorbed water in KBr pellets, wood, cellulose, clays & hygroscopic polymers",
        confirm: "Confirm broad O–H stretch 3200–3550 cm⁻¹",
    },
    FtirBand { min: 3260.0, max: 3330.0, group: "≡C–H (terminal alkyne)", mode: "ν(≡C–H) stretch", intensity: "Strong, Sharp", note: "Click-chemistry propargyl polymer end-groups", confirm: "Confirm C≡C 2100–2160 cm⁻¹" },
    FtirBand { min: 3000.0, max: 3110.0, group: "=C–H (aromatic / vinyl sp²)", mode: "ν(=C–H) stretch", intensity: "Weak–Med", note: "Polystyrene (3026, 3060, 3082), PET, vinyl monomers", confirm: "Confirm out-of-plane Ar–H 690–760 cm⁻¹" },
    FtirBand { min: 2945.0, max: 2995.0, group: "C–H (CH₃ asym / sp³)", mode: "ν_as(CH₃)", intensity: "Med–Strong", note: "PMMA, PP, PDMS (2962 cm⁻¹)", confirm: "Confirm δ_s(CH₃) umbrella 1365–1395 cm⁻¹" },
    FtirBand { min: 2900.0, max: 2948.0, group: "C–H (CH₂ asym / sp³)", mode: "ν_as(CH₂)", intensity: "Strong", note: "PE (2916), PP, PS, Nylon, PCL aliphatic backbone", confirm: "Confirm δ(CH₂) scissor 1465–1475 cm⁻¹" },
    FtirBand { min: 2865.0, max: 2902.0, group: "C–H (CH₃ sym / sp³)", mode: "ν_s(CH₃)", intensity: "Medium", note: "Polypropylene, PMMA, branched alkanes", confirm: "Confirm CH₃ asym 2945–2995 cm⁻¹" },
    FtirBand { min: 2840.0, max: 2865.0, group: "C–H (CH₂ sym / sp³)", mode: "ν_s(CH₂)", intensity: "Strong", note: "PE (2848), aliphatic polymer chains", confirm: "Confirm CH₂ asym 2910–2945 cm⁻¹" },
    FtirBand { min: 2700.0, max: 2830.0, group: "O=C–H (aldehyde Fermi doublet)", mode: "ν(C–H)", intensity: "Weak–Med (2 bands)", note: "Characteristic doublet near 2820 & 2720 cm⁻¹", confirm: "Confirm aldehyde C=O 1710–1740 cm⁻¹" },
    FtirBand { min: 2500.0, max: 3300.0, group: "O–H (carboxylic acid dimer)", mode: "ν(O–H) stretch", intensity: "Very Broad", note: "Poly(acrylic acid), poly(methacrylic acid)", confirm: "Confirm acid C=O 1705–1725 cm⁻¹" },
    FtirBand { min: 2550.0, max: 2600.0, group: "S–H (thiol)", mode: "ν(S–H) stretch", intensity: "Weak", note: "Thiol-ene click chemistry monomers & chain transfer agents", confirm: "Confirm C–S stretch 600–700 cm⁻¹" },
    FtirBand { min: 2250.0, max: 2285.0, group: "N=C=O (isocyanate)", mode: "ν_as(N=C=O)", intensity: "Very Strong", note: "Polyurethane synthesis! Disappears as NCO reacts with OH", confirm: "Confirm disappearance of O–H 3300–3500 cm⁻¹ on cure" },
    FtirBand { min: 2210.0, max: 2260.0, group: "C≡N (nitrile)", mode: "ν(C≡N) stretch", intensity: "Med–Strong, Sharp", note: "PAN (2243 cm⁻¹), SAN, ABS, NBR rubber", confirm: "Confirm CH₂ scissor 1454 cm⁻¹" },
    FtirBand { min: 2100.0, max: 2160.0, group: "C≡C (alkyne) / N₃ (azide ~2100)", mode: "ν(C≡C) / ν_as(N₃)", intensity: "Weak / Very Strong", note: "Azide-terminated polymers (ATRP → CuAAC click at ~2100)", confirm: "Very strong sharp ~2100; confirm 1250–1300 cm⁻¹" },
    FtirBand { min: 1800.0, max: 1850.0, group: "C=O (anhydride sym / acyl chloride)", mode: "ν(C=O) stretch", intensity: "Strong", note: "Maleic anhydride grafted polymers (1855/1780 doublet)", confirm: "Confirm second anhydride C=O ~1780 cm⁻¹" },
    FtirBand { min: 1765.0, max: 1795.0, group: "C=O (carbonate / vinyl ester / anhydride)", mode: "ν(C=O) stretch", intensity: "Very Strong", note: "Polycarbonate (PC ~1770), maleic anhydride (1780)", confirm: "Confirm C–O–C 1190–1260 cm⁻¹" },
    FtirBand { min: 1735.0, max: 1760.0, group: "C=O (saturated aliphatic ester)", mode: "ν(C=O) stretch", intensity: "Very Strong", note: "PLA (~1750), PVAc (~1738), saturated acrylates", confirm: "Confirm C–O–C 1130–1210 cm⁻¹" },
    FtirBand { min: 1715.0, max: 1738.0, group: "C=O (ester / methacrylate / PET / PCL)", mode: "ν(C=O) stretch", intensity: "Very Strong", note: "PMMA (~1728), PCL (~1725), PET (~1715), MMA monomer (~1722)", confirm: "Confirm C–O–C 1148 & 1192 cm⁻¹" },
    FtirBand { min: 1700.0, max: 1730.0, group: "C=O (urethane / carbamate)", mode: "ν(C=O) stretch", intensity: "Strong", note: "Polyurethane free (~1730) & H-bonded (~1705)", confirm: "Compare free vs H-bonded C=O at 1732 / 1703 cm⁻¹" },
    FtirBand { min: 1705.0, max: 1725.0, group: "C=O (aliphatic ketone / acid)", mode: "ν(C=O) stretch", intensity: "Very Strong", note: "PAA/PMAA carboxylic acid (~1705), aliphatic ketones (1715)", confirm: "Confirm broad acid O–H 2500–3300 cm⁻¹" },
    FtirBand { min: 1665.0, max: 1705.0, group: "C=O (α,β-unsaturated / aryl ketone)", mode: "ν(C=O) stretch", intensity: "Strong", note: "Conjugation lowers C=O by 20–40 cm⁻¹; PEEK (~1655)", confirm: "Confirm C=C 1620–1660 cm⁻¹" },
    FtirBand { min: 1630.0, max: 1685.0, group: "C=O (Amide I)", mode: "ν(C=O) amide", intensity: "Very Strong", note: "Nylon-6/6,6 (~1640), polyacrylamide (~1660), PNIPAM (~1645)", confirm: "Confirm Amide II 1515–1570 cm⁻¹" },
    FtirBand { min: 1620.0, max: 1660.0, group: "C=C (vinyl / acrylate monomer!)", mode: "ν(C=C) stretch", intensity: "Med–Strong", note: "KEY polymerization peak! MMA/acrylate 1638, styrene 1630 — disappears on cure!", confirm: "Key polymerization peak; confirm =C–H 905–990 cm⁻¹" },
    FtirBand { min: 1580.0, max: 1620.0, group: "C=C (aromatic ring skeletal)", mode: "ν(C=C) ring", intensity: "Med, Sharp", note: "PS (1601), PET (1578/1612), Epoxy/PC (1608)", confirm: "Confirm ring 1485–1520 cm⁻¹" },
    FtirBand { min: 1515.0, max: 1570.0, group: "N–H bend + C–N (Amide II)", mode: "δ(N–H) + ν(C–N)", intensity: "Strong", note: "Nylon-6,6 (~1540), Polyurethane (~1530), PNIPAM (~1545)", confirm: "Confirm Amide I 1630–1685 cm⁻¹" },
    FtirBand { min: 1485.0, max: 1520.0, group: "C=C (aromatic ring skeletal)", mode: "ν(C=C) ring", intensity: "Strong, Sharp", note: "Polystyrene (1493), Bisphenol-A Epoxy/PC (1508)", confirm: "Confirm 1580–1620 cm⁻¹" },
    FtirBand { min: 1430.0, max: 1480.0, group: "CH₂ scissor / CH₃ asym bend", mode: "δ(CH₂) / δ_as(CH₃)", intensity: "Med–Strong", note: "PE (1462/1472 crystalline doublet), PS (1452), PP (1456), PMMA (1448)", confirm: "Confirm CH₂ rock 720–730 cm⁻¹" },
    FtirBand { min: 1400.0, max: 1435.0, group: "=CH₂ in-plane scissor / C–H α-C=O", mode: "δ(=CH₂)", intensity: "Medium", note: "Acrylate/vinyl monomers (~1408), PVC (~1427)", confirm: "Confirm =C–H out-of-plane 905–990 cm⁻¹" },
    FtirBand { min: 1365.0, max: 1395.0, group: "CH₃ symmetric bend (umbrella)", mode: "δ_s(CH₃)", intensity: "Medium", note: "PP (1376 sharp), PMMA (1388), PVAc (1372), t-Bu doublet (1395/1365)", confirm: "Confirm CH₃ asym 2945–2995 cm⁻¹" },
    FtirBand { min: 1300.0, max: 1360.0, group: "C–N (aromatic amine) / SO₂ asym / CH wag", mode: "ν(C–N) / ν_as(SO₂)", intensity: "Med–Strong", note: "PET (1340 trans), PSU/PES sulfone (~1320), PVC (1330)", confirm: "Confirm ν_s(SO₂) 1145–1165 cm⁻¹" },
    FtirBand { min: 1250.0, max: 1275.0, group: "Si–CH₃ symmetric deformation", mode: "δ_s(Si–CH₃)", intensity: "Very Strong, Sharp", note: "Diagnostic PDMS/silicone peak at 1260 cm⁻¹!", confirm: "Confirm Si–CH₃ rock 795 cm⁻¹" },
    FtirBand { min: 1215.0, max: 1290.0, group: "C–O (ester / aryl ether / epoxy)", mode: "ν_as(C–O–C)", intensity: "Very Strong", note: "PET (1245), PMMA (1242/1270), Epoxy (1245), PVAc (1235), PC (1225)", confirm: "Confirm 1130–1210 cm⁻¹" },
    FtirBand { min: 1285.0, max: 1335.0, group: "C–O (α,β-unsaturated ester / acrylate)", mode: "ν(C–O) conj.", intensity: "Strong", note: "MMA/acrylate conjugated C–O stretch (1300/1324 cm⁻¹)", confirm: "Confirm ester C=O 1715–1738 & C=C 1620–1660 cm⁻¹" },
    FtirBand { min: 1140.0, max: 1220.0, group: "C–F / CF₂ stretch (fluoropolymer)", mode: "ν(CF₂)", intensity: "Very Strong", note: "PTFE doublet at 1200 & 1146 cm⁻¹, PVDF (1175)", confirm: "Confirm CF₂ wag 500–645 cm⁻¹" },
    FtirBand { min: 1130.0, max: 1210.0, group: "C–O–C (ester / aliphatic ether)", mode: "ν(C–O)", intensity: "Very Strong", note: "PMMA (1148/1192), PCL (1175), PLA (1182), Polyurethane", confirm: "Confirm 1215–1290 cm⁻¹" },
    FtirBand { min: 1050.0, max: 1150.0, group: "C–O–C (aliphatic ether / PEG / PVA)", mode: "ν_as(C–O–C)", intensity: "Very Strong", note: "PEG/PEO (1105), PET (1096), PVA (1092), cellulose (1055)", confirm: "Confirm strong PEG 1105 cm⁻¹" },
    FtirBand { min: 1000.0, max: 1100.0, group: "Si–O–Si (siloxane backbone)", mode: "ν_as(Si–O–Si)", intensity: "Very Strong, Broad", note: "PDMS characteristic doublet at ~1080 & ~1015 cm⁻¹", confirm: "Confirm Si–CH₃ 1260 & 795 cm⁻¹" },
    FtirBand { min: 960.0, max: 1005.0, group: "trans-C=C–H / vinyl =C–H OOP", mode: "γ(=C–H)", intensity: "Strong", note: "trans-Polybutadiene (966), PP isotacticity (998/973), vinyl monomers (990)", confirm: "Confirm C=C 1620–1660 cm⁻¹" },
    FtirBand { min: 905.0, max: 950.0, group: "Epoxide oxirane ring / =CH₂ wag", mode: "ν(C–O) ring / γ(=CH₂)", intensity: "Medium", note: "DGEBA Epoxy uncured ring (915 cm⁻¹ — tracks epoxy cure!); MMA =CH₂ (940)", confirm: "Confirm oxirane 830–850 cm⁻¹" },
    FtirBand { min: 820.0, max: 880.0, group: "para-Aromatic C–H OOP / Si–C / =CH₂", mode: "γ(Ar–H)", intensity: "Strong", note: "PET (872), Epoxy/PC (830 para-disubstituted), PMMA (842), MMA (814)", confirm: "Confirm ring 1500–1610 cm⁻¹" },
    FtirBand { min: 780.0, max: 820.0, group: "Si–CH₃ rock / trisubstituted C=CH", mode: "ρ(Si–CH₃)", intensity: "Very Strong", note: "PDMS (795 cm⁻¹), MMA monomer (814), PMMA (750)", confirm: "Confirm Si–O–Si 1000–1100 cm⁻¹" },
    FtirBand { min: 800.0, max: 825.0, group: "=CH₂ twist (vinyl / acrylate monomer)", mode: "τ(=CH₂)", intensity: "Med–Strong", note: "Vinylidene =CH₂ out-of-plane twist (MMA 814 cm⁻¹)", confirm: "Confirm C=C 1620–1660 & =CH₂ wag 935–950 cm⁻¹" },
    FtirBand { min: 735.0, max: 775.0, group: "Mono-/ortho-Aromatic C–H OOP", mode: "γ(Ar–H)", intensity: "Strong", note: "Polystyrene (756 cm⁻¹), PMMA (750), PET (725)", confirm: "Confirm ring bend 685–710 cm⁻¹" },
    FtirBand { min: 710.0, max: 735.0, group: "–(CH₂)_n– long-chain rocking (n≥4)", mode: "ρ(CH₂)", intensity: "Med–Strong", note: "PE crystalline doublet (720/730 cm⁻¹), Nylon-6,6 (725), PCL (732)", confirm: "Confirm CH₂ scissor 1465–1475 cm⁻¹" },
    FtirBand { min: 685.0, max: 710.0, group: "Mono-substituted phenyl ring bend", mode: "φ(ring)", intensity: "Very Strong", note: "Polystyrene hallmark peak at 698 cm⁻¹!", confirm: "Confirm Ar out-of-plane 735–775 cm⁻¹" },
    FtirBand { min: 600.0, max: 700.0, group: "C–Cl stretch (chlorinated polymer)", mode: "ν(C–Cl)", intensity: "Strong", note: "PVC syndiotactic/isotactic bands (615, 638, 690 cm⁻¹)", confirm: "Confirm CH–Cl 1250 cm⁻¹" },
    FtirBand { min: 500.0, max: 645.0, group: "CF₂ wag/bend / C–Br / aromatic ring", mode: "δ(CF₂) / ν(C–Br)", intensity: "Med–Strong", note: "PTFE (638, 554, 503), PS (540)", confirm: "Confirm CF₂ 1140–1220 cm⁻¹" },

    // ── Atmospheric / instrument artefacts ──────────────────────────────────
    FtirBand { min: 2330.0, max: 2375.0, group: "CO₂ (atmospheric, asym stretch)", mode: "ν_as(O=C=O)", intensity: "Sharp doublet", note: "Unpurged instrument artefact — NOT a sample band", confirm: "Sharp doublet ~2349; confirm CO₂ bend at 667 cm⁻¹" },
    FtirBand { min: 660.0, max: 675.0, group: "CO₂ (atmospheric, bend)", mode: "δ(O=C=O)", intensity: "Sharp", note: "Confirms an unpurged CO₂ artefact", confirm: "Confirm asymmetric stretch doublet 2330–2375 cm⁻¹" },
    // ── PVDF crystalline polymorphs ─────────────────────────────────────────
    FtirBand { min: 755.0, max: 772.0, group: "PVDF α-phase (TTTG)", mode: "γ(CH₂) + CF₂ rock", intensity: "Strong", note: "Non-polar α crystal; reference band at 763 cm⁻¹", confirm: "Confirm 795 & 976 cm⁻¹" },
    FtirBand { min: 833.0, max: 848.0, group: "PVDF β-phase (piezoelectric)", mode: "ν_s(CF₂) + CH₂ rock", intensity: "Very Strong", note: "All-trans polar β crystal — ferroelectric phase", confirm: "Confirm 1275 & 510 cm⁻¹; Gregorio F(β)" },
    FtirBand { min: 1268.0, max: 1285.0, group: "PVDF β-phase CF₂", mode: "ν_as(CF₂)", intensity: "Strong", note: "Second diagnostic β-phase band at 1275 cm⁻¹", confirm: "Confirm 840 & 510 cm⁻¹" },
    // ── Polyimide imidization tracking ──────────────────────────────────────
    FtirBand { min: 1770.0, max: 1788.0, group: "Polyimide Imide I (sym C=O)", mode: "ν_s(C=O) imide", intensity: "Strong, Sharp", note: "Kapton/PI imide ring; grows on thermal imidization", confirm: "Confirm 1720, 1375 & 725 cm⁻¹" },
    FtirBand { min: 1712.0, max: 1730.0, group: "Polyimide Imide II (asym C=O)", mode: "ν_as(C=O) imide", intensity: "Very Strong", note: "Imide carbonyl; polyamic acid C=O at ~1660 disappears", confirm: "Confirm Imide I 1770–1788 cm⁻¹" },
    FtirBand { min: 1365.0, max: 1388.0, group: "Polyimide Imide III (C–N stretch)", mode: "ν(C–N) imide ring", intensity: "Strong", note: "Imide C–N stretch at ~1375 cm⁻¹", confirm: "Confirm Imide IV 718–734 cm⁻¹" },
    FtirBand { min: 718.0, max: 734.0, group: "Polyimide Imide IV (ring bend)", mode: "δ(C=O) imide", intensity: "Medium", note: "Imide ring out-of-plane bend at ~725 cm⁻¹", confirm: "Confirm Imide I 1770–1788 cm⁻¹" },
    // ── Polyurethane microphase separation ──────────────────────────────────
    FtirBand { min: 1725.0, max: 1738.0, group: "Free Urethane C=O", mode: "ν(C=O) non-H-bonded", intensity: "Strong", note: "Soft-segment dissolved carbonyl", confirm: "Compare with H-bonded C=O at ~1703 cm⁻¹" },
    FtirBand { min: 1695.0, max: 1712.0, group: "H-bonded Urethane C=O", mode: "ν(C=O) H-bonded", intensity: "Strong", note: "Hard-segment ordered carbonyl — microphase separation", confirm: "Compare free C=O 1725–1738 cm⁻¹" },
    FtirBand { min: 3430.0, max: 3460.0, group: "Free N–H (urethane/amide)", mode: "ν(N–H) free", intensity: "Weak, Sharp", note: "Non-hydrogen-bonded N–H", confirm: "Compare H-bonded N–H 3300–3345 cm⁻¹" },
    FtirBand { min: 3300.0, max: 3345.0, group: "H-bonded N–H (urethane/amide)", mode: "ν(N–H) H-bonded", intensity: "Broad, Strong", note: "Ordered hard-segment N–H", confirm: "Compare free N–H 3430–3460 cm⁻¹" },
    // ── Polyethylene crystallinity & branching ──────────────────────────────
    FtirBand { min: 1372.0, max: 1382.0, group: "LDPE CH₃ branching", mode: "δ_s(CH₃)", intensity: "Weak", note: "Branch-point methyl; absent in linear HDPE", confirm: "Compare 1378 vs 1465 cm⁻¹" },
    FtirBand { min: 1468.0, max: 1476.0, group: "PE crystalline CH₂ scissor", mode: "δ(CH₂) orthorhombic", intensity: "Strong", note: "Crystalline splitting marker", confirm: "Confirm CH₂ rock 726–734 cm⁻¹" },
    FtirBand { min: 726.0, max: 734.0, group: "PE crystalline CH₂ rock", mode: "ρ(CH₂) crystalline", intensity: "Strong", note: "Orthorhombic crystalline PE; 720 amorphous", confirm: "Confirm CH₂ scissor 1468–1476 cm⁻¹" },
    // ── RAFT / ATRP / click end-groups ──────────────────────────────────────
    FtirBand { min: 2090.0, max: 2125.0, group: "Azide –N₃ (click end-group)", mode: "ν_as(N₃)", intensity: "Very Strong, Sharp", note: "Telechelic azide for CuAAC; vanishes on click", confirm: "Very strong sharp ~2100 cm⁻¹" },
    FtirBand { min: 1705.0, max: 1718.0, group: "Maleimide C=O (ring)", mode: "ν(C=O) imide", intensity: "Strong", note: "Maleimide/BMAm click handle", confirm: "Confirm maleimide ring 690–702 cm⁻¹" },
    FtirBand { min: 690.0, max: 702.0, group: "Maleimide ring (C–H OOP)", mode: "γ(C=C–H)", intensity: "Medium", note: "Maleimide enone ring", confirm: "Confirm C=O 1705–1718 cm⁻¹" },
    FtirBand { min: 1045.0, max: 1075.0, group: "Thiocarbonylthio C=S (RAFT)", mode: "ν(C=S)", intensity: "Weak–Med", note: "Trithiocarbonate / xanthate RAFT agent", confirm: "Confirm trithiocarbonate –SCH₂– at 3.35 ppm (¹H)" },
    FtirBand { min: 540.0, max: 650.0, group: "Aliphatic C–Br (ATRP initiator)", mode: "ν(C–Br)", intensity: "Medium", note: "α-Bromoester / –CH(Br)– ATRP end-group", confirm: "Confirm –CH(Br)– ¹H NMR at ~4.3 ppm" },
    // ── Nanocomposite fillers ───────────────────────────────────────────────
    FtirBand { min: 1060.0, max: 1130.0, group: "Silica (SiO₂) ν_as(Si–O–Si)", mode: "ν_as(Si–O–Si)", intensity: "Very Strong, Broad", note: "Fumed/colloidal silica filler", confirm: "Very strong broad ~1100; confirm 800 & 470 cm⁻¹" },
    FtirBand { min: 790.0, max: 815.0, group: "Silica (SiO₂) ν_s(Si–O)", mode: "ν_s(Si–O–Si)", intensity: "Strong", note: "Silica symmetric Si–O–Si", confirm: "Confirm 1100 & 470 cm⁻¹" },
    FtirBand { min: 455.0, max: 485.0, group: "Silica (SiO₂) δ(Si–O–Si)", mode: "δ(Si–O–Si)", intensity: "Strong", note: "Silica bending mode", confirm: "Confirm 1100 & 800 cm⁻¹" },
    FtirBand { min: 3615.0, max: 3640.0, group: "Montmorillonite (MMT) Al–OH", mode: "ν(Al–OH) structural", intensity: "Sharp, Weak", note: "Structural clay hydroxyl", confirm: "Sharp structural OH ~3625; confirm Si–O 1040 cm⁻¹" },
    FtirBand { min: 1025.0, max: 1055.0, group: "Montmorillonite (MMT) Si–O", mode: "ν(Si–O) clay", intensity: "Very Strong, Broad", note: "Silicate layer Si–O stretch", confirm: "Confirm Al–OH 3615–3640 cm⁻¹" },
    FtirBand { min: 510.0, max: 535.0, group: "Montmorillonite (MMT) Al–O–Si", mode: "δ(Al–O–Si)", intensity: "Medium", note: "Clay octahedral–tetrahedral bending", confirm: "Confirm Si–O 1025–1055 cm⁻¹" },
    FtirBand { min: 1405.0, max: 1455.0, group: "Calcium Carbonate (CaCO₃) CO₃²⁻", mode: "ν₃(CO₃²⁻) asym", intensity: "Very Strong, Broad", note: "Common mineral filler", confirm: "Very strong broad ~1420; confirm sharp 875 & 712 cm⁻¹" },
    FtirBand { min: 870.0, max: 882.0, group: "CaCO₃ out-of-plane CO₃", mode: "ν₂(CO₃²⁻)", intensity: "Sharp, Strong", note: "Calcite out-of-plane bend", confirm: "Confirm 1420 & 712 cm⁻¹" },
    FtirBand { min: 708.0, max: 718.0, group: "CaCO₃ in-plane CO₃", mode: "ν₄(CO₃²⁻)", intensity: "Sharp", note: "Calcite in-plane bend", confirm: "Confirm 1420 & 875 cm⁻¹" },
    FtirBand { min: 480.0, max: 700.0, group: "TiO₂ nanoparticles (Ti–O–Ti)", mode: "ν(Ti–O–Ti)", intensity: "Very Broad, Strong", note: "Anatase/rutile metal-oxide lattice", confirm: "Very broad strong metal-oxide band" },
    FtirBand { min: 420.0, max: 550.0, group: "ZnO nanoparticles", mode: "ν(Zn–O)", intensity: "Broad, Strong", note: "Wurtzite ZnO lattice", confirm: "Confirm broad 480–550 cm⁻¹" },
    FtirBand { min: 1720.0, max: 1740.0, group: "Graphene Oxide (GO) C=O", mode: "ν(C=O) carboxyl", intensity: "Medium", note: "GO edge carboxyl", confirm: "Confirm C=C / H₂O 1610–1635 cm⁻¹" },
    FtirBand { min: 1610.0, max: 1635.0, group: "Graphene Oxide (GO) C=C / H₂O", mode: "ν(C=C) + δ(H₂O)", intensity: "Strong", note: "Graphitic domain plus intercalated water", confirm: "Confirm C=O 1720–1740 cm⁻¹" },
    FtirBand { min: 1210.0, max: 1235.0, group: "Graphene Oxide (GO) epoxy C–O", mode: "ν(C–O–C) epoxy", intensity: "Strong", note: "Basal-plane epoxide", confirm: "Confirm alkoxy 1040–1060 cm⁻¹" },
    FtirBand { min: 1040.0, max: 1060.0, group: "Graphene Oxide (GO) alkoxy C–O", mode: "ν(C–O) alkoxy", intensity: "Strong", note: "GO basal-plane alkoxy", confirm: "Confirm epoxy 1210–1235 cm⁻¹" },
    // ── Polysulfone / PES / PEEK / cellulose / chitosan ─────────────────────
    FtirBand { min: 1290.0, max: 1335.0, group: "Sulfone ν_as(SO₂) (PSU/PES)", mode: "ν_as(SO₂)", intensity: "Strong", note: "Diaryl sulfone asymmetric stretch", confirm: "Confirm ν_s(SO₂) at 1145–1165 cm⁻¹" },
    FtirBand { min: 1140.0, max: 1170.0, group: "Sulfone ν_s(SO₂) (PSU/PES)", mode: "ν_s(SO₂)", intensity: "Strong", note: "Diaryl sulfone symmetric stretch", confirm: "Confirm ν_as(SO₂) 1290–1335 cm⁻¹" },
    FtirBand { min: 1230.0, max: 1255.0, group: "Diaryl ether Ar–O–Ar (PEEK/PSU)", mode: "ν_as(Ar–O–Ar)", intensity: "Very Strong", note: "High-performance engineering thermoplastic", confirm: "Confirm sulfone 1290–1335 cm⁻¹" },
    FtirBand { min: 1150.0, max: 1165.0, group: "Chitosan/Cellulose β-glycosidic C–O–C", mode: "ν_as(C–O–C)", intensity: "Strong", note: "Polysaccharide glycosidic linkage", confirm: "Confirm 1025–1075 cm⁻¹" },
    FtirBand { min: 1025.0, max: 1075.0, group: "Cellulose/Chitosan C–O ring", mode: "ν(C–O) pyranose", intensity: "Very Strong, Broad", note: "Dominant polysaccharide band", confirm: "Confirm 895 glycosidic cm⁻¹" },
    FtirBand { min: 885.0, max: 905.0, group: "β-glycosidic linkage (cellulose)", mode: "ν(β-(1→4) C–O–C)", intensity: "Medium", note: "Anomeric β-linkage", confirm: "Confirm 1025–1075 cm⁻¹" },
    // ── Conducting polymers ─────────────────────────────────────────────────
    FtirBand { min: 1575.0, max: 1600.0, group: "Polyaniline (PANI) quinoid ring", mode: "ν(C=N/C=C) quinoid", intensity: "Strong", note: "Emeraldine-base quinoid segment", confirm: "Confirm Benzenoid ring at 1490–1505 cm⁻¹" },
    FtirBand { min: 1485.0, max: 1510.0, group: "Polyaniline (PANI) benzenoid ring", mode: "ν(C=C) benzenoid", intensity: "Strong", note: "Emeraldine-base benzenoid segment", confirm: "Confirm quinoid 1575–1600 cm⁻¹" },
    FtirBand { min: 1120.0, max: 1140.0, group: "PANI doped (emeraldine salt)", mode: "ν(C–N⁺)", intensity: "Broad, Strong", note: "Conducting protonated form", confirm: "Confirm quinoid 1575–1600 cm⁻¹" },
    // ── Additional engineering & commodity polymers ─────────────────────────
    FtirBand { min: 1645.0, max: 1665.0, group: "PEEK aryl ketone C=O", mode: "ν(C=O) diaryl ketone", intensity: "Strong", note: "Conjugated diaryl ketone carbonyl", confirm: "Confirm diaryl ether 1230–1255 cm⁻¹" },
    FtirBand { min: 1185.0, max: 1215.0, group: "PEEK Ar–O–Ar", mode: "ν(Ar–O–Ar)", intensity: "Very Strong", note: "PEEK ether linkage", confirm: "Confirm C=O 1645–1665 cm⁻¹" },
    FtirBand { min: 1080.0, max: 1120.0, group: "POM C–O–C (polyoxymethylene)", mode: "ν_as(C–O–C)", intensity: "Very Strong", note: "Acetal backbone", confirm: "Confirm CH₂ rock 930–950 cm⁻¹" },
    FtirBand { min: 930.0, max: 950.0, group: "POM CH₂ rock", mode: "ρ(CH₂)", intensity: "Strong", note: "Polyoxymethylene helix", confirm: "Confirm 1080–1120 cm⁻¹" },
    FtirBand { min: 1630.0, max: 1650.0, group: "PA6 Amide I (Nylon-6)", mode: "ν(C=O) amide", intensity: "Very Strong", note: "Nylon-6 α-form", confirm: "Confirm Amide II 1530–1550 cm⁻¹" },
    FtirBand { min: 1260.0, max: 1300.0, group: "Nylon Amide III", mode: "δ(N–H) + ν(C–N)", intensity: "Medium", note: "Amide III combination", confirm: "Confirm Amide I / II 1630–1550 cm⁻¹" },
    FtirBand { min: 1705.0, max: 1725.0, group: "PBT ester C=O", mode: "ν(C=O) aryl ester", intensity: "Very Strong", note: "Poly(butylene terephthalate)", confirm: "Confirm p-Ar out-of-plane 815–835 cm⁻¹" },
    FtirBand { min: 815.0, max: 835.0, group: "PBT/PET p-Ar out-of-plane", mode: "γ(Ar–H)", intensity: "Strong", note: "para-Disubstituted terephthalate", confirm: "Confirm ester C=O 1705–1725 cm⁻¹" },
    FtirBand { min: 2235.0, max: 2250.0, group: "ABS/SAN nitrile", mode: "ν(C≡N)", intensity: "Medium, Sharp", note: "Acrylonitrile component", confirm: "Confirm butadiene C=C 960–980 cm⁻¹" },
    FtirBand { min: 960.0, max: 980.0, group: "ABS butadiene C=C", mode: "γ(=C–H) trans", intensity: "Medium", note: "Polybutadiene rubber phase", confirm: "Confirm nitrile 2235–2250 cm⁻¹" },
    FtirBand { min: 1640.0, max: 1660.0, group: "PNIPAM Amide I", mode: "ν(C=O) amide", intensity: "Very Strong", note: "Thermoresponsive LCST polymer", confirm: "Confirm Amide II 1540–1560 cm⁻¹" },
    FtirBand { min: 1540.0, max: 1560.0, group: "PNIPAM Amide II", mode: "δ(N–H) + ν(C–N)", intensity: "Strong", note: "Isopropylamide N–H bend", confirm: "Confirm Amide I 1640–1660 cm⁻¹" },
    FtirBand { min: 3300.0, max: 3450.0, group: "PHEMA O–H", mode: "ν(O–H) H-bonded", intensity: "Broad, Strong", note: "Hydroxyethyl methacrylate hydrogel", confirm: "Confirm ester C=O 1720–1735 cm⁻¹" },
    FtirBand { min: 1715.0, max: 1735.0, group: "PHEMA ester C=O", mode: "ν(C=O) methacrylate", intensity: "Very Strong", note: "Hydrogel ester carbonyl", confirm: "Confirm O–H 3300–3450 cm⁻¹" },
    FtirBand { min: 1700.0, max: 1720.0, group: "PAA carboxylic acid C=O", mode: "ν(C=O) acid", intensity: "Very Strong", note: "Poly(acrylic acid)", confirm: "Confirm broad acid O–H 2500–3300 cm⁻¹" },
    FtirBand { min: 1655.0, max: 1685.0, group: "PVP lactam C=O", mode: "ν(C=O) N-vinylpyrrolidone", intensity: "Very Strong", note: "Poly(vinylpyrrolidone)", confirm: "Confirm C–N 1270–1300 cm⁻¹" },
    FtirBand { min: 1270.0, max: 1300.0, group: "PVP lactam C–N", mode: "ν(C–N) lactam", intensity: "Strong", note: "PVP ring C–N", confirm: "Confirm lactam C=O 1655–1685 cm⁻¹" },
    FtirBand { min: 3200.0, max: 3500.0, group: "PVA O–H (polyvinyl alcohol)", mode: "ν(O–H) H-bonded", intensity: "Very Broad, Strong", note: "Partially hydrolysed PVAc", confirm: "Confirm C–O 1080–1120 cm⁻¹" },
    FtirBand { min: 1080.0, max: 1120.0, group: "PVA C–O", mode: "ν(C–O) secondary alcohol", intensity: "Strong", note: "PVA backbone hydroxyl", confirm: "Confirm O–H 3200–3500 cm⁻¹" },
    FtirBand { min: 1725.0, max: 1745.0, group: "PBA ester C=O (acrylate)", mode: "ν(C=O) acrylate", intensity: "Very Strong", note: "Poly(butyl acrylate) soft segment", confirm: "Confirm C–O–C 1150–1200 cm⁻¹" },
    FtirBand { min: 1745.0, max: 1765.0, group: "PLGA ester C=O", mode: "ν(C=O) aliphatic ester", intensity: "Very Strong", note: "Lactide/glycolide copolymer", confirm: "Confirm C–O–C 1080–1180 cm⁻¹" },
    FtirBand { min: 1715.0, max: 1735.0, group: "PHB ester C=O (polyhydroxybutyrate)", mode: "ν(C=O) bacterial polyester", intensity: "Very Strong", note: "Biopolyester carbonyl", confirm: "Confirm crystalline 1215–1235 cm⁻¹" },
    FtirBand { min: 1215.0, max: 1235.0, group: "PHB crystalline C–O–C", mode: "ν(C–O–C) crystalline", intensity: "Strong", note: "PHB crystallinity marker", confirm: "Confirm C=O 1715–1735 cm⁻¹" },
    FtirBand { min: 1635.0, max: 1660.0, group: "PBD 1,4-cis C=C", mode: "ν(C=C) cis", intensity: "Medium", note: "cis-1,4-polybutadiene microstructure", confirm: "Confirm =C–H out-of-plane 690–740 cm⁻¹" },
    FtirBand { min: 960.0, max: 975.0, group: "PBD 1,4-trans C=C", mode: "γ(=C–H) trans", intensity: "Strong", note: "trans-1,4-polybutadiene", confirm: "Confirm C=C 1635–1660 cm⁻¹" },
    FtirBand { min: 905.0, max: 925.0, group: "PBD 1,2-vinyl", mode: "γ(=CH₂) vinyl", intensity: "Medium", note: "1,2-polybutadiene pendant vinyl", confirm: "Confirm 990–1000 cm⁻¹" },
    FtirBand { min: 880.0, max: 895.0, group: "Polyisoprene 3,4-unit", mode: "γ(=CH₂) 3,4", intensity: "Medium", note: "3,4-polyisoprene microstructure", confirm: "Confirm 1640 cm⁻¹ C=C" },
    FtirBand { min: 835.0, max: 850.0, group: "Polyisoprene 1,4-unit", mode: "γ(=C–H) 1,4", intensity: "Medium", note: "1,4-polyisoprene (natural rubber)", confirm: "Confirm 1665 cm⁻¹ C=C" },
    FtirBand { min: 3350.0, max: 3500.0, group: "Cured epoxy O–H", mode: "ν(O–H) network", intensity: "Broad", note: "Grows as oxirane ring opens", confirm: "Confirm oxirane 915 cm⁻¹ disappearance" },
    FtirBand { min: 910.0, max: 920.0, group: "Epoxy oxirane ring", mode: "ν(C–O) epoxide", intensity: "Medium", note: "Consumed during cure", confirm: "Confirm p-Ar out-of-plane 830 cm⁻¹" },
];

// ============================================================================
// ¹H NMR FUNCTIONAL GROUP & POLYMER DATABASE
// ============================================================================
pub const HNMR_DB: &[NmrBand] = &[
    NmrBand { min: 10.5, max: 13.0, group: "R–COOH (carboxylic acid)", example: "PAA, PMAA, PLA acid end-group", note: "Broad, exchangeable with D₂O", j_hz: "br (exchangeable)" },
    NmrBand { min: 9.4, max: 10.2, group: "R–CHO (aldehyde)", example: "Formyl end-group, PMMA–CHO", note: "Sharp triplet from coupling to α-CH₂", j_hz: "J = 7–8 Hz (t)" },
    NmrBand { min: 8.8, max: 9.2, group: "H–CO–O (formate ester)", example: "RAFT formate end-group", note: "Singlet, no α-protons", j_hz: "s (singlet)" },
    NmrBand { min: 8.2, max: 8.8, group: "Electron-poor heteroaryl Ar–H", example: "Pyridine α-H (8.62), PET Ar–H (8.10)", note: "Deshielded by electronegative ring", j_hz: "d, J = 7–9 Hz (ortho)" },
    NmrBand { min: 7.8, max: 8.2, group: "Amide / urethane N–H", example: "Nylon NH (7.80), PU NH (7.35–7.9)", note: "Broad; shifts strongly with H-bonding", j_hz: "br s (exchangeable)" },
    NmrBand { min: 7.5, max: 7.8, group: "ortho-Substituted / diaryl Ar–H", example: "PET Ar–H (8.10), PEEK Ar–H", note: "Two ortho substituents", j_hz: "d, J = 7–9 Hz (ortho)" },
    NmrBand { min: 7.2, max: 7.5, group: "Mono-substituted Ar–H (m/p)", example: "PS m/p (7.08), Toluene m/p (7.25)", note: "meta/para ring protons", j_hz: "m (multiplet)" },
    NmrBand { min: 6.9, max: 7.2, group: "Mono-substituted Ar–H (ortho)", example: "PS ortho (6.55), BHT ArH (6.98)", note: "Shielded by neighbouring rings", j_hz: "d, J = 7–9 Hz (ortho)" },
    NmrBand { min: 6.4, max: 6.9, group: "PS ortho / vinyl aromatic =CH", example: "PS o-ArH (6.55), Styrene =CH (6.72)", note: "Ring-current shielded ortho protons", j_hz: "dd, J = 17.6, 10.8 Hz (styrene)" },
    NmrBand { min: 5.9, max: 6.4, group: "Vinyl =CH₂ (cis to ester)", example: "MMA (6.10), acrylates (6.2–6.4)", note: "KEY monomer conversion marker!", j_hz: "dd, J = 17.4, 1.4 Hz" },
    NmrBand { min: 5.4, max: 5.9, group: "Vinyl =CH₂ (trans to ester)", example: "MMA (5.56), styrene (5.75)", note: "Second vinyl proton of monomer", j_hz: "dd, J = 10.8, 1.4 Hz" },
    NmrBand { min: 5.0, max: 5.4, group: "PLA –CH– / PBD 1,4-cis", example: "PLA CH (5.16), PBD cis (5.40)", note: "PLA methine quartet", j_hz: "q, J = 7.0 Hz (PLA)" },
    NmrBand { min: 4.8, max: 5.2, group: "PBD 1,2-vinyl –CH= / PVAc CH", example: "PBD 1,2 CH (5.05), PVAc CH (4.87)", note: "Pendant vinyl microstructure", j_hz: "br m" },
    NmrBand { min: 4.4, max: 4.8, group: "PET –OCH₂– / PVC –CHCl–", example: "PET (4.70), PVC CHCl (4.45)", note: "Also cellulose anomeric H (4.5)", j_hz: "s (singlet)" },
    NmrBand { min: 4.1, max: 4.4, group: "ATRP –CH(Br)– / PCL –OCH₂–", example: "–CH(Br)– (4.3), PCL (4.06)", note: "α-Bromoester ATRP end-group", j_hz: "t, J = 6.7 Hz (PCL)" },
    NmrBand { min: 3.9, max: 4.1, group: "PVA –CH(OH)– / PLA end-group", example: "PVA (3.95), lactic acid end (4.1)", note: "Hydroxyl-bearing methine", j_hz: "br m" },
    NmrBand { min: 3.7, max: 3.9, group: "PEG/PEO –OCH₂–", example: "PEG (3.64–3.75)", note: "Very strong singlet of PEG backbone", j_hz: "s (singlet)" },
    NmrBand { min: 3.5, max: 3.7, group: "PMMA –O–CH₃ (methoxy)", example: "PMMA (3.60), MMA monomer (3.74)", note: "Strongest singlet in PMMA", j_hz: "s (singlet)" },
    NmrBand { min: 3.3, max: 3.5, group: "RAFT trithiocarbonate –SCH₂–", example: "CTA –CH₂–S–C(=S)– (3.30–3.40)", note: "Trithiocarbonate RAFT end-group", j_hz: "t, J = 7.2 Hz" },
    NmrBand { min: 3.1, max: 3.3, group: "N–CH₂ (nylon) / PAN –CH(CN)–", example: "Nylon (3.20), PAN (3.15, DMSO-d₆)", note: "Electronegative N/CN deshielding", j_hz: "m (multiplet)" },
    NmrBand { min: 2.9, max: 3.1, group: "DMF / DMAc N–CH₃", example: "DMF (2.96), DMAc (3.02)", note: "Polymerization-solvent impurity", j_hz: "s (singlet)" },
    NmrBand { min: 2.6, max: 2.9, group: "DMSO-d₆ residual / NMP N–CH₃", example: "DMSO (2.54), NMP (2.85 in CDCl₃)", note: "DMSO residual is a quintet", j_hz: "quintet, J = 1.9 Hz (DMSO-d₆)" },
    NmrBand { min: 2.4, max: 2.6, group: "Et₃N N–CH₂", example: "Triethylamine (2.53 CDCl₃)", note: "Paired with CH₃ triplet at 1.03", j_hz: "q, J = 7.1 Hz" },
    NmrBand { min: 2.2, max: 2.4, group: "PCL α-CH₂ / Ar–CH₃ (toluene)", example: "PCL (2.31), Toluene (2.36)", note: "α-Carbonyl methylene", j_hz: "t, J = 7.4 Hz" },
    NmrBand { min: 2.0, max: 2.2, group: "PVAc –OCOCH₃ / allylic CH₂", example: "PVAc (2.04), PBD allylic (2.05)", note: "Acetate singlet", j_hz: "s (singlet)" },
    NmrBand { min: 1.8, max: 2.0, group: "PMMA –CH₂– / PS –CH– backbone", example: "PMMA (1.81–1.95), PS (1.85)", note: "Broadened by tacticity", j_hz: "br m" },
    NmrBand { min: 1.6, max: 1.8, group: "PLA –CH₃ / PE –CH₂–", example: "PLA (1.58), PE (1.29)", note: "PLA methyl doublet", j_hz: "d, J = 7.0 Hz (PLA)" },
    NmrBand { min: 1.4, max: 1.6, group: "PBD 1,4-trans CH₂ / aliphatic", example: "PBD trans (1.45), PS CH₂ (1.43)", note: "Methylene backbone", j_hz: "br m" },
    NmrBand { min: 1.2, max: 1.4, group: "PCL/PE –CH₂– backbone", example: "PCL (1.38), PE (1.29)", note: "Long aliphatic methylene envelope", j_hz: "m (multiplet)" },
    NmrBand { min: 1.1, max: 1.25, group: "PMMA mm α-CH₃ (isotactic)", example: "PMMA mm (1.21)", note: "Isotactic triad marker", j_hz: "br s" },
    NmrBand { min: 0.98, max: 1.08, group: "PMMA mr α-CH₃ (heterotactic)", example: "PMMA mr (1.02), EtO–CH₃ (1.25)", note: "Heterotactic triad marker", j_hz: "br s" },
    NmrBand { min: 0.78, max: 0.90, group: "PMMA rr α-CH₃ (syndiotactic)", example: "PMMA rr (0.85), PP CH₃ (0.88)", note: "Free-radical PMMA is predominantly rr!", j_hz: "br s" },
    NmrBand { min: 0.6, max: 0.78, group: "Terminal aliphatic –CH₃", example: "n-hexane CH₃ (0.88), grease CH₃", note: "Chain-end methyl triplet", j_hz: "t, J = 6.8 Hz" },
    NmrBand { min: -0.2, max: 0.3, group: "Si–CH₃ (PDMS / TMS / silicone)", example: "PDMS (0.07), TMS (0.00)", note: "Classic silicone-grease singlet", j_hz: "s (singlet)" },
    NmrBand { min: 4.60, max: 4.75, group: "Cellulose anomeric H-1", example: "Cellulose (4.5–4.6, DMSO-d₆)", note: "β-(1→4) glycosidic anomeric proton", j_hz: "d, J = 7–8 Hz" },
    NmrBand { min: 3.4, max: 4.0, group: "Chitosan/cellulose ring H-2…H-6", example: "Chitosan (3.4–3.9)", note: "Overlapping pyranose ring protons", j_hz: "br m" },
    NmrBand { min: 6.0, max: 7.5, group: "Acrylamide / PNIPAM N–H", example: "PNIPAM NH (7.0–7.5 in D₂O)", note: "Exchangeable amide proton", j_hz: "br (exchangeable)" },
    NmrBand { min: 1.0, max: 1.15, group: "PNIPAM –CH(CH₃)₂", example: "PNIPAM isopropyl (1.10)", note: "Six-proton doublet", j_hz: "d, J = 6.8 Hz" },
    NmrBand { min: 7.2, max: 7.4, group: "Benzyl RAFT/ATRP end-group Ar–H", example: "PhCH₂– (7.25–7.35)", note: "Aromatic end-group", j_hz: "m (multiplet)" },
    NmrBand { min: 5.0, max: 5.2, group: "Benzyl end-group –PhCH₂–", example: "Benzyl ester/ether (5.05–5.15)", note: "Benzylic methylene end-group", j_hz: "s (singlet)" },
    NmrBand { min: 1.6, max: 1.75, group: "AIBN / cumyl –C(CH₃)₂CN", example: "AIBN (1.70), cumyl end (1.68)", note: "Radical-initiator fragment", j_hz: "s (singlet)" },
    NmrBand { min: 1.5, max: 2.0, group: "Thiol –S–H (thiol-ene)", example: "Thiol chain-transfer agent", note: "Weak, exchangeable", j_hz: "t, J = 7–8 Hz (exchangeable)" },
    NmrBand { min: 5.3, max: 5.5, group: "PBD 1,4-trans –CH=CH–", example: "trans-PBD (5.38)", note: "Internal trans alkene", j_hz: "br m" },
    NmrBand { min: 5.5, max: 5.7, group: "PBD 1,2-vinyl =CH₂", example: "1,2-PBD (5.55)", note: "Pendant vinyl methylene", j_hz: "dd, J = 10.9, 17.2 Hz" },
    NmrBand { min: 4.6, max: 4.8, group: "Polyisoprene 3,4-unit =CH₂", example: "3,4-PI (4.68)", note: "Pendant isopropenyl", j_hz: "br s" },
    NmrBand { min: 5.0, max: 5.2, group: "Polyisoprene 1,4-unit –CH=", example: "1,4-PI (5.10)", note: "Internal trisubstituted alkene", j_hz: "br m" },
    NmrBand { min: 4.8, max: 4.95, group: "Vinyl acetate monomer =CH₂", example: "VAc (4.87)", note: "Monomer conversion marker", j_hz: "dd, J = 6.7, 14.0 Hz" },
    NmrBand { min: 7.25, max: 7.35, group: "Vinyl acetate monomer =CH–", example: "VAc (7.28)", note: "Monomer methine", j_hz: "dd, J = 6.7, 14.0 Hz" },
    NmrBand { min: 2.1, max: 2.3, group: "PVP backbone –CH– / –CH₂–", example: "PVP (2.0–2.4 in D₂O)", note: "Broad lactam backbone", j_hz: "br m" },
    NmrBand { min: 3.6, max: 3.8, group: "PVP lactam N–CH₂", example: "PVP (3.6–3.8 in D₂O)", note: "Ring N-methylene", j_hz: "br m" },
    NmrBand { min: 7.0, max: 7.4, group: "PANI aromatic N–H / ArH", example: "Polyaniline (7.0–7.3, DMSO-d₆)", note: "Exchangeable amine plus aromatic H", j_hz: "br m" },
];

// ============================================================================
// ¹³C NMR FUNCTIONAL GROUP & POLYMER DATABASE
// ============================================================================
pub const CNMR_DB: &[NmrBand] = &[
    NmrBand { min: 190.0, max: 220.0, group: "C=O (ketone / aldehyde)", example: "Acetone (206.3), PEEK (195)", note: "Most downfield carbons", j_hz: "s (¹³C)" },
    NmrBand { min: 170.0, max: 185.0, group: "C=O (ester / acid / amide)", example: "PMMA (177.8), PCL (173.6), PLA (169.6)", note: "PMMA C=O splits into pentads", j_hz: "s (¹³C)" },
    NmrBand { min: 160.0, max: 170.0, group: "C=O (aryl ester / urethane)", example: "PET (166.0), PU (155)", note: "Heteroatom-conjugated carbonyl", j_hz: "s (¹³C)" },
    NmrBand { min: 152.0, max: 160.0, group: "C=O (carbonate) / Ar C–O", example: "PC (152.1), PU urethane (154)", note: "Carbonate and carbamate carbonyls", j_hz: "s (¹³C)" },
    NmrBand { min: 140.0, max: 152.0, group: "Aromatic ipso-C", example: "PS ipso (145.3), PC C–O (148.9)", note: "PS ipso-C is tacticity-sensitive", j_hz: "s (¹³C)" },
    NmrBand { min: 126.0, max: 140.0, group: "Aromatic CH (o/m/p)", example: "PS (128.0/125.6), PET (129.6)", note: "Protonated aromatic carbons", j_hz: "d, ¹J_CH ≈ 160 Hz" },
    NmrBand { min: 118.0, max: 126.0, group: "Aromatic CH / C≡N", example: "PAN (120.1), PS p (125.6)", note: "Nitrile and para-aromatic overlap", j_hz: "s (¹³C)" },
    NmrBand { min: 112.0, max: 118.0, group: "C≡N (nitrile) / CF₂", example: "PAN (120.1), PTFE CF₂ (112)", note: "Quaternary sp carbon", j_hz: "s (¹³C)" },
    NmrBand { min: 74.0, max: 80.0, group: "CDCl₃ triplet / C–O", example: "CDCl₃ (77.16), PEG (70.6)", note: "Solvent triplet sits here", j_hz: "t, ¹J_CD ≈ 32 Hz" },
    NmrBand { min: 66.0, max: 74.0, group: "C–O (ether / ester / alcohol)", example: "PEG (70.6), PLA CH (69.0), PCL OCH₂ (64.1)", note: "Oxygenated sp³ carbons", j_hz: "s (¹³C)" },
    NmrBand { min: 58.0, max: 66.0, group: "C–O (PET OCH₂ / cellulose C-6)", example: "PET (63.2), cellulose (60–65)", note: "Methylene adjacent to oxygen", j_hz: "s (¹³C)" },
    NmrBand { min: 50.0, max: 58.0, group: "PMMA –OCH₃ / PVC –CHCl–", example: "PMMA OCH₃ (51.8), PVC CHCl (56.5)", note: "Methoxy and chloromethine", j_hz: "s (¹³C)" },
    NmrBand { min: 42.0, max: 50.0, group: "Backbone –CH– / PMMA Cq", example: "PMMA Cq (44.6), PS CH (40.5)", note: "Quaternary and methine backbone", j_hz: "s (¹³C)" },
    NmrBand { min: 36.0, max: 42.0, group: "PS –CH₂– / PCL α-CH₂", example: "PS CH₂ (44.0), PCL (34.1)", note: "Backbone methylene", j_hz: "s (¹³C)" },
    NmrBand { min: 28.0, max: 36.0, group: "Aliphatic –CH₂–", example: "PE (30.0), PP CH (28.8), PCL (28.4)", note: "Crystalline PE at 32.8, amorphous 30.0", j_hz: "s (¹³C)" },
    NmrBand { min: 20.0, max: 28.0, group: "Aliphatic –CH₃ / PP CH₃", example: "PP (21.6), PCL (25.5)", note: "Methyl carbons", j_hz: "s (¹³C)" },
    NmrBand { min: 14.0, max: 20.0, group: "PMMA α-CH₃ triads / PLA CH₃", example: "PMMA mm (21.1)/mr (18.7)/rr (16.5), PLA (16.7)", note: "Classic ¹³C tacticity region", j_hz: "s (¹³C)" },
    NmrBand { min: -2.0, max: 5.0, group: "Si–CH₃ (PDMS / TMS)", example: "PDMS (1.0), TMS (0.0)", note: "Silicone grease at ~1.0 ppm", j_hz: "s (¹³C)" },
    NmrBand { min: 175.0, max: 182.0, group: "PMMA C=O pentads", example: "PMMA mm (178.6)/mr (177.8)/rr (176.9)", note: "Pentad-resolution carbonyl", j_hz: "s (¹³C)" },
    NmrBand { min: 44.0, max: 46.0, group: "PMMA quaternary C (tacticity)", example: "PMMA mm (45.0)/rr (44.3)", note: "Cq triad splitting", j_hz: "s (¹³C)" },
    NmrBand { min: 16.0, max: 22.0, group: "PMMA α-CH₃ triads (mm/mr/rr)", example: "PMMA mm (21.1)/mr (18.7)/rr (16.5)", note: "Syndiotactic rr dominates radical PMMA", j_hz: "s (¹³C)" },
    NmrBand { min: 165.0, max: 174.0, group: "PLA C=O / PCL C=O", example: "PLA (169.6), PCL (173.6)", note: "Aliphatic polyester carbonyl", j_hz: "s (¹³C)" },
    NmrBand { min: 66.0, max: 70.0, group: "PLA –CH– (tetrads)", example: "PLA CH (69.0, tetrad splitting)", note: "Methine tetrad microstructure", j_hz: "s (¹³C)" },
    NmrBand { min: 26.0, max: 36.0, group: "PBD 1,4-cis / 1,4-trans CH₂", example: "cis (27.6), trans (32.9)", note: "Diene microstructure", j_hz: "s (¹³C)" },
    NmrBand { min: 125.0, max: 136.0, group: "PBD 1,4-cis/trans –CH=", example: "cis (129.6), trans (130.4)", note: "Internal alkene carbons", j_hz: "s (¹³C)" },
    NmrBand { min: 142.0, max: 150.0, group: "PBD 1,2-vinyl –CH= (pendant)", example: "1,2-PBD (142.8)", note: "Pendant vinyl methine", j_hz: "s (¹³C)" },
    NmrBand { min: 112.0, max: 116.0, group: "PBD 1,2-vinyl =CH₂", example: "1,2-PBD (114.0)", note: "Pendant vinyl methylene", j_hz: "s (¹³C)" },
    NmrBand { min: 134.0, max: 136.0, group: "Polyisoprene 1,4 –C(CH₃)=", example: "1,4-PI (135.0)", note: "Trisubstituted alkene", j_hz: "s (¹³C)" },
    NmrBand { min: 124.0, max: 126.0, group: "Polyisoprene 1,4 –CH=", example: "1,4-PI (125.0)", note: "Internal alkene CH", j_hz: "s (¹³C)" },
    NmrBand { min: 106.0, max: 112.0, group: "Polyisoprene 3,4 =CH₂", example: "3,4-PI (109.0)", note: "Pendant isopropenyl", j_hz: "s (¹³C)" },
    NmrBand { min: 60.0, max: 68.0, group: "PVA –CH(OH)–", example: "PVA (64–68, DMSO-d₆)", note: "Tacticity-sensitive methine", j_hz: "s (¹³C)" },
    NmrBand { min: 20.0, max: 24.0, group: "PVAc –OCOCH₃", example: "PVAc acetate CH₃ (21.1)", note: "Acetate methyl", j_hz: "s (¹³C)" },
    NmrBand { min: 168.0, max: 172.0, group: "PVAc C=O", example: "PVAc (170.5)", note: "Acetate carbonyl", j_hz: "s (¹³C)" },
    NmrBand { min: 30.0, max: 34.0, group: "PAN –CH₂– / –CH(CN)–", example: "PAN CH₂ (33.0), CH (27.4)", note: "Nitrile polymer backbone", j_hz: "s (¹³C)" },
    NmrBand { min: 172.0, max: 178.0, group: "Nylon amide C=O", example: "Nylon-6,6 (173.5)", note: "Amide carbonyl", j_hz: "s (¹³C)" },
    NmrBand { min: 102.0, max: 108.0, group: "Cellulose C-1 (anomeric)", example: "Cellulose (105, DMSO-d₆)", note: "Anomeric acetal carbon", j_hz: "s (¹³C)" },
    NmrBand { min: 80.0, max: 90.0, group: "Cellulose C-4 (crystalline)", example: "Cellulose (88 cryst / 84 amorph)", note: "Crystallinity marker", j_hz: "s (¹³C)" },
    NmrBand { min: 71.0, max: 80.0, group: "Cellulose C-2/C-3/C-5", example: "Cellulose (72–75)", note: "Ring carbons", j_hz: "s (¹³C)" },
    NmrBand { min: 55.0, max: 62.0, group: "Chitosan C-2 (N-acetyl glucosamine)", example: "Chitosan (56–60)", note: "Glucosamine ring carbon", j_hz: "s (¹³C)" },
    NmrBand { min: 140.0, max: 160.0, group: "PANI quinoid C=N/C=C", example: "PANI (156, 148)", note: "Emeraldine-base ring carbons", j_hz: "s (¹³C)" },
    NmrBand { min: 118.0, max: 128.0, group: "PANI benzenoid aromatic CH", example: "PANI (122–126)", note: "Benzenoid ring carbons", j_hz: "s (¹³C)" },
    NmrBand { min: 154.0, max: 158.0, group: "PU urethane C=O", example: "Polyurethane (155)", note: "Carbamate carbonyl", j_hz: "s (¹³C)" },
    NmrBand { min: 64.0, max: 72.0, group: "PU/polyether –OCH₂–", example: "PU (64.0), PEG soft segment (70.5)", note: "Polyether soft segment", j_hz: "s (¹³C)" },
];

// ============================================================================
// COMPREHENSIVE UV-VIS CHROMOPHORE, POLYMER & NANOPARTICLE DATABASE
// ============================================================================
pub const UVVIS_DB: &[UvVisBand] = &[
    // --- Simple & Conjugated Organic Chromophores ---
    UvVisBand { min: 175.0, max: 198.0, group: "C=C isolated alkene", transition: "π → π*", epsilon: "10,000–15,000", note: "Unconjugated vinyl/alkene; residual monomer end-group" },
    UvVisBand { min: 215.0, max: 258.0, group: "Conjugated diene (C=C–C=C)", transition: "π → π*", epsilon: "15,000–25,000", note: "1,3-butadiene (~217 nm), cyclic dienes (235–258 nm); Woodward-Fieser" },
    UvVisBand { min: 258.0, max: 315.0, group: "Conjugated triene / short polyene", transition: "π → π*", epsilon: "30,000–55,000", note: "1,3,5-hexatriene (~268 nm); thermal PVC dehydrochlorination marker" },
    UvVisBand { min: 210.0, max: 248.0, group: "α,β-Unsaturated C=O (enone/acrylate)", transition: "π → π*", epsilon: "8,000–18,000", note: "MMA/acrylate vinyl-ester conjugation (~212–235 nm), enones (~220–248 nm)" },
    UvVisBand { min: 310.0, max: 355.0, group: "α,β-Unsaturated C=O (enone)", transition: "n → π*", epsilon: "20–100", note: "Forbidden carbonyl n→π* shifted red by C=C conjugation" },
    UvVisBand { min: 200.0, max: 216.0, group: "Ester / Carboxylic acid C=O", transition: "n → π*", epsilon: "40–100", note: "PMMA, aliphatic polyesters, acetic acid" },
    UvVisBand { min: 190.0, max: 210.0, group: "Amide / Peptide bond (–CONH–)", transition: "π → π*", epsilon: "7,000–12,000", note: "Nylon, polyacrylamide, protein backbone (190–205 nm)" },
    UvVisBand { min: 212.0, max: 230.0, group: "Amide / Peptide C=O", transition: "n → π*", epsilon: "50–200", note: "Weak n→π* shoulder on amide backbone" },
    UvVisBand { min: 270.0, max: 298.0, group: "Ketone / Aldehyde C=O", transition: "n → π*", epsilon: "15–35", note: "Acetone (279 nm), aliphatic ketones & photo-oxidation carbonyls" },
    UvVisBand { min: 196.0, max: 210.0, group: "Benzene (E₂ band)", transition: "π → π* (¹A₁g → ¹B₁u)", epsilon: "7,000–9,000", note: "Primary UV aromatic ring absorption" },
    UvVisBand { min: 248.0, max: 264.0, group: "Benzene (B band, vibronic)", transition: "π → π* (¹A₁g → ¹B₂u)", epsilon: "200–300", note: "Characteristic 5-finger vibronic fine structure (~254 nm)" },
    UvVisBand { min: 256.0, max: 270.0, group: "Toluene / Alkylbenzene", transition: "π → π* (B band)", epsilon: "250–400", note: "Bathochromic shift +4–8 nm vs benzene from hyperconjugation" },
    UvVisBand { min: 210.0, max: 224.0, group: "Polystyrene / Phenyl ring (E₂)", transition: "π → π*", epsilon: "8,000–12,000", note: "Polystyrene (PS) strong UV band (~215 nm)" },
    UvVisBand { min: 244.0, max: 254.0, group: "Styrene monomer (vinyl-phenyl)", transition: "π → π* (K band)", epsilon: "14,000–16,000", note: "Conjugated C=C–Ph (~248 nm); disappears upon polymerization to PS" },
    UvVisBand { min: 254.0, max: 268.0, group: "Polystyrene (PS) aromatic B-band", transition: "π → π* (¹A₁g → ¹B₂u)", epsilon: "220–320", note: "PS phenyl vibronic bands at 254, 261, 268 nm" },
    UvVisBand { min: 280.0, max: 294.0, group: "Styrene monomer (B band)", transition: "π → π*", epsilon: "500–800", note: "Residual styrene monomer marker (~282 & 291 nm)" },
    UvVisBand { min: 208.0, max: 220.0, group: "Phenol (neutral, E₂)", transition: "π → π*", epsilon: "6,200–7,500", note: "Phenolic OH auxochrome shifts benzene E₂ red" },
    UvVisBand { min: 266.0, max: 280.0, group: "Phenol (neutral, B band)", transition: "π → π*", epsilon: "1,400–2,000", note: "Phenol ~270 nm; shifts to 287 nm upon deprotonation (phenolate)" },
    UvVisBand { min: 282.0, max: 296.0, group: "Phenolate anion (ArO⁻)", transition: "π → π* (CT)", epsilon: "2,400–3,200", note: "Alkaline phenol / bisphenol-A pH-sensitive red shift" },
    UvVisBand { min: 226.0, max: 238.0, group: "Aniline (E₂ band)", transition: "π → π*", epsilon: "8,000–9,500", note: "Primary aromatic amine (~230 nm); shifts back to 203 nm as PhNH₃⁺" },
    UvVisBand { min: 276.0, max: 290.0, group: "Aniline (B band)", transition: "π → π* / n→π*", epsilon: "1,400–1,800", note: "Aniline monomer (~280 nm); PANI precursor" },
    UvVisBand { min: 250.0, max: 263.0, group: "Pyridine / Aza-aromatic", transition: "π → π* / n→π*", epsilon: "2,000–3,000", note: "Pyridine (~257 nm); slight blue shift in polar protic solvents" },
    UvVisBand { min: 255.0, max: 275.0, group: "Nitroaromatic (Ar–NO₂)", transition: "π → π* (CT)", epsilon: "7,000–10,000", note: "Nitrobenzene (~260–268 nm); strong intramolecular CT" },
    // --- UV Stabilizers & PAHs ---
    UvVisBand { min: 248.0, max: 266.0, group: "Benzophenone UV absorber (Band I)", transition: "π → π*", epsilon: "16,000–22,000", note: "2-Hydroxybenzophenone polymer UV stabilizer (~254–260 nm)" },
    UvVisBand { min: 322.0, max: 352.0, group: "Benzophenone UV absorber (Band II)", transition: "π → π* / ESIPT", epsilon: "10,000–15,000", note: "Photoprotective UV-A absorption (~328–345 nm)" },
    UvVisBand { min: 295.0, max: 318.0, group: "Benzotriazole UV stabilizer (UV-B)", transition: "π → π*", epsilon: "13,000–17,000", note: "Tinuvin P / 328 polymer light stabilizer (~300–310 nm)" },
    UvVisBand { min: 335.0, max: 360.0, group: "Benzotriazole UV stabilizer (UV-A)", transition: "π → π* (H-bond)", epsilon: "14,000–18,000", note: "Intramolecular H-bond ESIPT UV-A band (~340–355 nm)" },
    UvVisBand { min: 218.0, max: 228.0, group: "Naphthalene (β-band)", transition: "π → π* (¹B_b)", epsilon: "95,000–120,000", note: "Intense naphthalene UV-C band (~221 nm)" },
    UvVisBand { min: 268.0, max: 290.0, group: "Naphthalene (p-band)", transition: "π → π* (¹L_a)", epsilon: "5,000–7,000", note: "Vibronic progression 275, 286 nm; PET/PEN naphthalene analogues" },
    UvVisBand { min: 246.0, max: 258.0, group: "Anthracene (β-band)", transition: "π → π*", epsilon: "180,000–220,000", note: "Very intense PAH UV band (~252 nm)" },
    UvVisBand { min: 338.0, max: 382.0, group: "Anthracene (p-band vibronic)", transition: "π → π* (¹L_a)", epsilon: "7,000–9,500", note: "Characteristic 4-finger vibronic peaks (340, 356, 375 nm)" },
    UvVisBand { min: 235.0, max: 246.0, group: "Pyrene (S₃/S₄ band)", transition: "π → π*", epsilon: "60,000–85,000", note: "Pyrene fluorescent probe (~241 nm)" },
    UvVisBand { min: 318.0, max: 342.0, group: "Pyrene (S₂ band vibronic)", transition: "π → π*", epsilon: "45,000–55,000", note: "Strong vibronic bands at 320 & 335 nm; excimer probe" },
    UvVisBand { min: 312.0, max: 328.0, group: "trans-Azobenzene", transition: "π → π*", epsilon: "20,000–24,000", note: "Photochromic trans isomer (~318–322 nm); drops under UV irradiation" },
    UvVisBand { min: 432.0, max: 458.0, group: "Azobenzene (cis/trans)", transition: "n → π*", epsilon: "400–1,500", note: "Visible band (~440–445 nm); intensifies in cis-azobenzene" },
    // --- Conducting & Commodity Polymers ---
    UvVisBand { min: 315.0, max: 340.0, group: "PANI benzenoid ring", transition: "π → π*", epsilon: "High (film/soln)", note: "Polyaniline EB & ES benzenoid absorption (~325 nm)" },
    UvVisBand { min: 418.0, max: 448.0, group: "PANI polaron → π*", transition: "Polaron band", epsilon: "Medium–High", note: "H⁺-doped Emeraldine Salt (ES) conductive polaron (~430 nm)" },
    UvVisBand { min: 605.0, max: 655.0, group: "PANI quinoid exciton", transition: "Exciton (benzenoid → quinoid)", epsilon: "High", note: "Undoped Emeraldine Base (EB) blue chromophore (~630 nm)" },
    UvVisBand { min: 775.0, max: 870.0, group: "PANI polaron / bipolaron (NIR)", transition: "π → polaron (free carrier)", epsilon: "Broad, High", note: "Doped Emeraldine Salt delocalized polaron tail (~815 nm)" },
    UvVisBand { min: 218.0, max: 236.0, group: "PEDOT:PSS aromatic / thiophene", transition: "π → π*", epsilon: "High", note: "PSS phenyl & PEDOT backbone UV absorption (~225 nm)" },
    UvVisBand { min: 750.0, max: 950.0, group: "PEDOT:PSS polaron / bipolaron NIR", transition: "Polaron / Bipolaron", epsilon: "Broad NIR", note: "Conductive oxidized PEDOT free-carrier NIR absorption" },
    UvVisBand { min: 440.0, max: 565.0, group: "P3HT conjugated backbone", transition: "π → π*", epsilon: "High", note: "Regioregular P3HT (~450 nm solution; 520–560 nm + 602 nm π-stack film)" },
    UvVisBand { min: 360.0, max: 415.0, group: "Polypyrrole (PPy) backbone", transition: "π → π*", epsilon: "High", note: "Conducting polypyrrole neutral/oxidized transition (~385 nm)" },
    UvVisBand { min: 760.0, max: 890.0, group: "Polypyrrole (PPy) bipolaron", transition: "Bipolaron NIR", epsilon: "Broad", note: "Doped polypyrrole free-carrier NIR band (~820 nm)" },
    UvVisBand { min: 410.0, max: 465.0, group: "PPV (Poly(p-phenylene vinylene))", transition: "π → π*", epsilon: "High", note: "Conjugated electroluminescent polymer (~430 nm)" },
    UvVisBand { min: 238.0, max: 256.0, group: "PET terephthalate ester-aryl", transition: "π → π* (CT)", epsilon: "14,000–18,000", note: "PET/PBT aromatic ester band (~244 nm)" },
    UvVisBand { min: 285.0, max: 304.0, group: "PET / PEN aromatic UV-B", transition: "π → π*", epsilon: "1,800–2,500", note: "PET absorption edge (~292–300 nm)" },
    // --- Semiconductors, Carbon Nanomaterials & Plasmonic NPs ---
    UvVisBand { min: 365.0, max: 398.0, group: "Anatase TiO₂ band edge", transition: "O 2p → Ti 3d (indirect Eg)", epsilon: "Semiconductor edge", note: "Absorption edge ~385 nm (Eg ≈ 3.20–3.22 eV); photocatalytic UV-A" },
    UvVisBand { min: 396.0, max: 418.0, group: "Rutile TiO₂ band edge", transition: "O 2p → Ti 3d (Eg)", epsilon: "Semiconductor edge", note: "Absorption edge ~410 nm (Eg ≈ 3.02 eV)" },
    UvVisBand { min: 360.0, max: 386.0, group: "ZnO wurtzite band edge", transition: "Direct band-edge / exciton", epsilon: "Semiconductor edge", note: "Sharp excitonic onset ~375 nm (Eg ≈ 3.30–3.37 eV)" },
    UvVisBand { min: 440.0, max: 515.0, group: "CdS quantum dots / semiconductor", transition: "1S_e–1S_h / band edge", epsilon: "Quantum confined", note: "Size-tunable Eg ~2.4–2.8 eV (450–510 nm)" },
    UvVisBand { min: 224.0, max: 238.0, group: "Graphene Oxide (GO) C=C", transition: "π → π*", epsilon: "High", note: "Aromatic sp² C=C domains in GO (~230 nm)" },
    UvVisBand { min: 295.0, max: 315.0, group: "Graphene Oxide (GO) C=O shoulder", transition: "n → π*", epsilon: "Shoulder", note: "Epoxide/carbonyl n→π* shoulder (~305 nm); vanishes on reduction" },
    UvVisBand { min: 262.0, max: 278.0, group: "Reduced Graphene Oxide (rGO)", transition: "π → π* (restored sp²)", epsilon: "High", note: "Red-shift from 230→270 nm confirms electronic conjugation recovery" },
    UvVisBand { min: 395.0, max: 432.0, group: "Ag NPs (Silver nanoparticles LSPR)", transition: "Surface Plasmon (LSPR)", epsilon: "10⁸–10⁹ M⁻¹cm⁻¹", note: "Spherical Ag NPs (10–40 nm) dipole LSPR centered ~412 nm" },
    UvVisBand { min: 514.0, max: 542.0, group: "Au NPs (Gold nanoparticles LSPR)", transition: "Surface Plasmon (LSPR)", epsilon: "10⁸–10⁹ M⁻¹cm⁻¹", note: "Spherical Au NPs (12–30 nm) ruby-red LSPR centered ~524 nm" },
    UvVisBand { min: 650.0, max: 850.0, group: "Au Nanorods (longitudinal LSPR)", transition: "Longitudinal LSPR", epsilon: "10⁹–10¹⁰ M⁻¹cm⁻¹", note: "Aspect-ratio tunable NIR plasmon band" },
    // --- Common Analytical Dyes, Pigments & Fullerenes ---
    UvVisBand { min: 484.0, max: 502.0, group: "Fluorescein (dianion, pH > 7)", transition: "π → π* (S₀ → S₁)", epsilon: "75,000–88,000", note: "Sharp visible band at ~490–494 nm; bright green fluorescence" },
    UvVisBand { min: 520.0, max: 538.0, group: "Rhodamine 6G", transition: "π → π* (S₀ → S₁)", epsilon: "100,000–116,000", note: "Xanthene laser dye (~526–530 nm in EtOH/H₂O)" },
    UvVisBand { min: 546.0, max: 562.0, group: "Rhodamine B", transition: "π → π* (S₀ → S₁)", epsilon: "95,000–110,000", note: "Photocatalysis probe dye (~554 nm in water)" },
    UvVisBand { min: 654.0, max: 674.0, group: "Methylene Blue (monomer)", transition: "π → π* (S₀ → S₁)", epsilon: "74,000–95,000", note: "Photocatalytic degradation standard (~664 nm; dimer shoulder ~612 nm)" },
    UvVisBand { min: 444.0, max: 486.0, group: "β-Carotene (conjugated polyene)", transition: "π → π* (S₀ → S₂)", epsilon: "125,000–145,000", note: "11 conjugated C=C; vibronic peaks ~452 & 480 nm" },
    UvVisBand { min: 395.0, max: 426.0, group: "Porphyrin Soret band (B band)", transition: "π → π* (S₀ → S₂)", epsilon: "200,000–450,000", note: "Ultra-intense macrocycle Soret peak (~400–420 nm)" },
    UvVisBand { min: 505.0, max: 650.0, group: "Porphyrin Q-bands (I–IV)", transition: "π → π* (S₀ → S₁)", epsilon: "5,000–25,000", note: "4 Q-bands (free base) collapsing to 2 bands upon metallation" },
    UvVisBand { min: 424.0, max: 440.0, group: "Chlorophyll a (Soret band)", transition: "π → π* (B band)", epsilon: "110,000–130,000", note: "Blue Soret peak (~430 nm in acetone/ether)" },
    UvVisBand { min: 658.0, max: 672.0, group: "Chlorophyll a (Q_y band)", transition: "π → π* (Q_y)", epsilon: "75,000–90,000", note: "Red chlorophyll a diagnostic peak (~662–665 nm)" },
    UvVisBand { min: 252.0, max: 266.0, group: "Fullerene C₆₀ (UV band)", transition: "π → π*", epsilon: "140,000–175,000", note: "Buckminsterfullerene intense UV band (~257 nm) + 330 nm" },
];

// ============================================================================
// FULMER 2010 / GOTTLIEB 1997 NMR SOLVENT & IMPURITY DATABASE
// Columns: [CDCl3, DMSO, CD3OD, D2O, Acetone, THF, C6D6, Toluene, CD2Cl2, CD3CN]
// ============================================================================
pub const SOLVENT_IMPURITIES: &[SolventImpurity] = &[
    SolventImpurity {
        name: "Residual Solvent", group: "Solvent lock", mult: "m",
        h_shifts: [Some(7.26), Some(2.50), Some(3.31), Some(4.79), Some(2.05), Some(3.58), Some(7.16), Some(2.08), Some(5.32), Some(1.94)],
        c_shifts: [Some(77.16), Some(39.52), Some(49.00), None, Some(29.84), Some(67.21), Some(128.06), Some(20.43), Some(53.84), Some(1.32)],
        note: "Primary residual proton/carbon peak of deuterated solvent",
    },
    SolventImpurity {
        name: "Residual Solvent (2)", group: "Solvent lock #2", mult: "m",
        h_shifts: [None, None, None, None, None, Some(1.72), None, Some(7.01), None, None],
        c_shifts: [None, None, None, None, Some(206.26), Some(25.31), None, Some(128.87), None, Some(118.26)],
        note: "Secondary residual peak (THF-d₈, Acetone-d₆, Toluene-d₈, CD₃CN)",
    },
    SolventImpurity {
        name: "Water (H₂O / HDO)", group: "OH", mult: "s / br",
        h_shifts: [Some(1.56), Some(3.33), Some(4.87), Some(4.79), Some(2.84), Some(2.46), Some(0.40), Some(0.43), Some(1.52), Some(2.13)],
        c_shifts: [None; 10],
        note: "Very common moisture peak! Shift varies slightly with concentration",
    },
    SolventImpurity {
        name: "Silicone Grease (PDMS)", group: "Si–CH₃", mult: "s",
        h_shifts: [Some(0.07), Some(0.04), Some(0.10), None, Some(0.13), Some(0.07), Some(0.29), Some(0.26), Some(0.08), Some(0.08)],
        c_shifts: [Some(1.04), None, Some(2.10), None, Some(1.40), Some(1.19), Some(1.38), Some(1.37), Some(1.22), None],
        note: "Joint grease / syringe lubricant — sharp singlet near 0.07 ppm!",
    },
    SolventImpurity {
        name: "H-Grease (CH₂)", group: "CH₂", mult: "br s",
        h_shifts: [Some(1.26), Some(1.24), Some(1.29), None, Some(1.29), Some(1.29), Some(1.36), Some(1.36), Some(1.27), Some(1.27)],
        c_shifts: [Some(29.76), Some(29.20), None, None, Some(30.73), Some(30.45), Some(30.21), Some(30.15), Some(30.14), None],
        note: "Hydrocarbon grease / vacuum oil / PE-like aliphatic contaminant",
    },
    SolventImpurity {
        name: "H-Grease (CH₃)", group: "CH₃", mult: "m",
        h_shifts: [Some(0.86), Some(0.85), Some(0.88), None, Some(0.87), Some(0.88), Some(0.92), Some(0.92), Some(0.87), Some(0.87)],
        c_shifts: [Some(14.14), None, None, None, None, None, None, None, None, None],
        note: "Terminal CH₃ of hydrocarbon grease",
    },
    SolventImpurity {
        name: "Acetone", group: "CH₃", mult: "s",
        h_shifts: [Some(2.17), Some(2.09), Some(2.15), Some(2.22), Some(2.09), Some(2.05), Some(1.55), Some(1.55), Some(2.12), Some(2.08)],
        c_shifts: [Some(30.92), Some(30.56), Some(30.67), Some(30.89), Some(30.60), Some(30.17), Some(30.14), Some(30.04), Some(31.00), Some(30.91)],
        note: "Glassware rinse solvent! Also C=O at ~206–207 ppm in ¹³C",
    },
    SolventImpurity {
        name: "Dichloromethane (DCM)", group: "CH₂", mult: "s",
        h_shifts: [Some(5.30), Some(5.76), Some(5.49), None, Some(5.63), Some(5.51), Some(4.27), Some(4.32), Some(5.33), Some(5.44)],
        c_shifts: [Some(53.52), Some(54.84), Some(54.78), None, Some(54.95), Some(54.67), Some(53.46), Some(53.47), Some(54.24), Some(55.32)],
        note: "Sharp singlet at 5.30 ppm in CDCl₃; hard to pump off polymers",
    },
    SolventImpurity {
        name: "Chloroform (CHCl₃)", group: "CH", mult: "s",
        h_shifts: [Some(7.26), Some(8.32), Some(7.90), None, Some(8.02), Some(7.89), Some(6.15), Some(6.10), Some(7.32), Some(7.58)],
        c_shifts: [Some(77.36), Some(79.16), Some(79.44), None, Some(79.19), Some(79.24), Some(77.79), Some(77.89), Some(77.99), Some(79.17)],
        note: "Residual CHCl₃ trapped in precipitated polymer film",
    },
    SolventImpurity {
        name: "THF (CH₂ 2,5)", group: "O–CH₂", mult: "m",
        h_shifts: [Some(3.76), Some(3.60), Some(3.71), Some(3.74), Some(3.63), Some(3.62), Some(3.57), Some(3.54), Some(3.69), Some(3.64)],
        c_shifts: [Some(67.97), Some(67.03), Some(68.83), Some(68.68), Some(68.07), Some(67.96), Some(67.80), Some(67.75), Some(68.16), Some(68.33)],
        note: "Common polymer/GPC solvent — watch for BHT stabilizer too!",
    },
    SolventImpurity {
        name: "THF (CH₂ 3,4)", group: "CH₂", mult: "m",
        h_shifts: [Some(1.85), Some(1.76), Some(1.87), Some(1.88), Some(1.79), Some(1.78), Some(1.40), Some(1.43), Some(1.82), Some(1.80)],
        c_shifts: [Some(25.62), Some(25.14), Some(26.48), Some(25.67), Some(26.15), Some(26.19), Some(25.72), Some(25.79), Some(25.98), Some(26.27)],
        note: "Upfield multiplet of THF",
    },
    SolventImpurity {
        name: "BHT (ArH) [Stabilizer]", group: "Ar–H", mult: "s",
        h_shifts: [Some(6.98), Some(6.87), Some(6.92), None, Some(6.96), Some(6.90), Some(7.05), Some(7.02), Some(6.98), Some(6.97)],
        c_shifts: [Some(125.55), Some(124.76), Some(126.11), None, Some(126.03), Some(125.76), Some(125.83), Some(125.75), Some(125.88), Some(126.38)],
        note: "THF/monomer antioxidant! Also shows t-Bu (1.43), ArCH₃ (2.27), OH (5.01)",
    },
    SolventImpurity {
        name: "BHT (t-Bu) [Stabilizer]", group: "C(CH₃)₃", mult: "s",
        h_shifts: [Some(1.43), Some(1.36), Some(1.40), None, Some(1.41), Some(1.39), Some(1.38), Some(1.40), Some(1.42), Some(1.39)],
        c_shifts: [Some(30.33), Some(30.33), Some(30.77), None, Some(30.64), Some(30.72), Some(30.34), Some(30.40), Some(30.54), Some(30.77)],
        note: "Huge 18H singlet if THF wasn't distilled!",
    },
    SolventImpurity {
        name: "Methanol (MeOH)", group: "CH₃", mult: "s",
        h_shifts: [Some(3.49), Some(3.16), Some(3.34), Some(3.34), Some(3.31), Some(3.26), Some(3.07), Some(3.07), Some(3.42), Some(3.28)],
        c_shifts: [Some(50.41), Some(48.59), Some(49.86), Some(49.50), Some(49.77), Some(49.64), Some(49.97), Some(49.90), Some(50.45), Some(49.90)],
        note: "Standard precipitation non-solvent for PMMA/PS!",
    },
    SolventImpurity {
        name: "Ethanol (CH₂)", group: "O–CH₂", mult: "q, J=7",
        h_shifts: [Some(3.72), Some(3.44), Some(3.60), Some(3.65), Some(3.57), Some(3.54), Some(3.34), Some(3.32), Some(3.66), Some(3.54)],
        c_shifts: [Some(58.28), Some(56.07), Some(58.26), Some(58.05), Some(57.72), Some(57.60), Some(57.86), Some(57.81), Some(58.57), Some(57.96)],
        note: "Paired with CH₃ triplet at 1.25 ppm (CDCl₃)",
    },
    SolventImpurity {
        name: "Ethanol (CH₃)", group: "CH₃", mult: "t, J=7",
        h_shifts: [Some(1.25), Some(1.06), Some(1.19), Some(1.17), Some(1.12), Some(1.11), Some(0.96), Some(0.96), Some(1.21), Some(1.12)],
        c_shifts: [Some(18.41), Some(18.51), Some(18.40), Some(17.47), Some(18.89), Some(18.90), Some(18.72), Some(18.71), Some(18.69), Some(18.80)],
        note: "Paired with O–CH₂ quartet at 3.72 ppm (CDCl₃)",
    },
    SolventImpurity {
        name: "Ethyl Acetate (CH₂)", group: "O–CH₂", mult: "q, J=7",
        h_shifts: [Some(4.12), Some(4.03), Some(4.09), Some(4.14), Some(4.05), Some(4.05), Some(3.89), Some(3.89), Some(4.08), Some(4.06)],
        c_shifts: [Some(60.49), Some(59.74), Some(61.50), Some(62.32), Some(60.56), Some(60.30), Some(60.21), Some(60.08), Some(60.63), Some(60.98)],
        note: "Also CH₃CO (s) 2.05 and CH₃ (t) 1.26 ppm in CDCl₃",
    },
    SolventImpurity {
        name: "Ethyl Acetate (CH₃CO)", group: "CH₃C=O", mult: "s",
        h_shifts: [Some(2.05), Some(1.99), Some(2.01), Some(2.07), Some(1.97), Some(1.96), Some(1.65), Some(1.67), Some(2.01), Some(1.97)],
        c_shifts: [Some(21.04), Some(20.68), Some(20.88), Some(21.15), Some(20.83), Some(20.45), Some(20.56), Some(20.46), Some(21.15), Some(21.16)],
        note: "Paired with 4.12 (q) and 1.26 (t)",
    },
    SolventImpurity {
        name: "Diethyl Ether (CH₂)", group: "O–CH₂", mult: "q, J=7",
        h_shifts: [Some(3.48), Some(3.38), Some(3.49), Some(3.56), Some(3.41), Some(3.38), Some(3.26), Some(3.26), Some(3.45), Some(3.42)],
        c_shifts: [Some(65.91), Some(62.05), Some(66.88), Some(66.42), Some(66.12), Some(66.23), Some(65.94), Some(65.96), Some(66.11), Some(66.32)],
        note: "Common PEG/polymer precipitation solvent; CH₃ (t) at 1.21 (CDCl₃)",
    },
    SolventImpurity {
        name: "Diethyl Ether (CH₃)", group: "CH₃", mult: "t, J=7",
        h_shifts: [Some(1.21), Some(1.09), Some(1.18), Some(1.17), Some(1.11), Some(1.11), Some(1.11), Some(1.10), Some(1.17), Some(1.12)],
        c_shifts: [Some(15.20), Some(15.12), Some(15.46), Some(14.77), Some(15.78), Some(15.49), Some(15.46), Some(15.47), Some(15.44), Some(15.63)],
        note: "Paired with O–CH₂ quartet at 3.48 ppm (CDCl₃)",
    },
    SolventImpurity {
        name: "DMF (CHO)", group: "H–C=O", mult: "s",
        h_shifts: [Some(8.02), Some(7.95), Some(7.97), Some(7.92), Some(7.96), Some(7.89), Some(7.63), Some(7.57), Some(7.96), Some(7.92)],
        c_shifts: [Some(162.62), Some(162.29), Some(164.73), Some(165.53), Some(162.79), Some(162.00), Some(162.13), Some(161.93), Some(162.57), Some(163.31)],
        note: "RAFT/ATRP polymerization solvent! Paired with N–CH₃ at 2.96 & 2.88",
    },
    SolventImpurity {
        name: "DMF (N–CH₃)", group: "N(CH₃)₂", mult: "2×s",
        h_shifts: [Some(2.96), Some(2.89), Some(2.99), Some(3.01), Some(2.94), Some(2.86), Some(2.36), Some(2.36), Some(2.91), Some(2.89)],
        c_shifts: [Some(36.50), Some(35.73), Some(36.89), Some(37.54), Some(36.15), Some(35.65), Some(35.25), Some(35.22), Some(36.56), Some(36.57)],
        note: "Second N–CH₃ singlet at 2.88 (CDCl₃) / 2.73 (DMSO)",
    },
    SolventImpurity {
        name: "Toluene (CH₃)", group: "Ar–CH₃", mult: "s",
        h_shifts: [Some(2.36), Some(2.30), Some(2.33), None, Some(2.32), Some(2.31), Some(2.11), Some(2.11), Some(2.34), Some(2.33)],
        c_shifts: [Some(21.46), Some(20.99), Some(21.50), None, Some(21.46), Some(21.37), Some(21.10), Some(21.37), Some(21.53), Some(21.50)],
        note: "Classic radical polymerization solvent; Ar–H at 7.17–7.25",
    },
    SolventImpurity {
        name: "1,4-Dioxane", group: "O–CH₂", mult: "s",
        h_shifts: [Some(3.71), Some(3.57), Some(3.66), Some(3.75), Some(3.59), Some(3.57), Some(3.35), Some(3.36), Some(3.65), Some(3.60)],
        c_shifts: [Some(67.14), Some(66.36), Some(68.11), Some(67.19), Some(67.60), Some(67.65), Some(67.16), Some(67.17), Some(67.47), Some(67.72)],
        note: "Very common RAFT polymerization solvent — sharp 8H singlet",
    },
    SolventImpurity {
        name: "n-Hexane (CH₂/CH₃)", group: "CH₂ / CH₃", mult: "m / t",
        h_shifts: [Some(1.26), Some(1.25), Some(1.29), None, Some(1.28), Some(1.29), Some(1.24), Some(1.24), Some(1.27), Some(1.27)],
        c_shifts: [Some(31.64), Some(30.95), Some(32.76), None, Some(32.30), Some(32.45), Some(31.96), Some(32.02), Some(31.99), Some(32.36)],
        note: "Polymer precipitation solvent; CH₃ triplet at 0.88 (CDCl₃)",
    },
    SolventImpurity {
        name: "DMSO (Residual)", group: "S(O)CH₃", mult: "s",
        h_shifts: [Some(2.62), Some(2.54), Some(2.65), Some(2.71), Some(2.52), Some(2.49), Some(1.68), Some(1.68), Some(2.55), Some(2.50)],
        c_shifts: [Some(40.76), Some(40.45), Some(40.45), Some(39.39), Some(41.23), Some(41.05), Some(40.03), Some(40.24), Some(41.31), Some(41.31)],
        note: "High-boiling solvent — hard to remove by rotary evaporation",
    },
    SolventImpurity {
        name: "Acetic Acid", group: "CH₃", mult: "s",
        h_shifts: [Some(2.10), Some(1.91), Some(1.99), Some(2.08), Some(1.96), Some(1.94), Some(1.55), Some(1.58), Some(2.06), Some(1.96)],
        c_shifts: [Some(20.81), Some(20.95), Some(20.56), Some(21.03), Some(20.51), Some(20.13), Some(20.37), Some(20.27), Some(20.91), Some(20.73)],
        note: "C=O at ~176–178 ppm in ¹³C",
    },
    SolventImpurity {
        name: "Triethylamine (Et₃N)", group: "N–CH₂", mult: "q, J=7",
        h_shifts: [Some(2.53), Some(2.43), Some(2.58), Some(2.57), Some(2.45), Some(2.44), Some(2.40), Some(2.39), Some(2.48), Some(2.45)],
        c_shifts: [Some(46.25), Some(45.74), Some(46.96), Some(47.19), Some(47.07), Some(47.14), Some(46.77), Some(46.82), Some(46.75), Some(47.10)],
        note: "Paired with CH₃ triplet at 1.03 (CDCl₃)",
    },

    SolventImpurity {
        name: "DMF (N–CH₃ #2)", group: "N(CH₃)₂", mult: "s",
        h_shifts: [Some(2.88), Some(2.73), Some(2.92), Some(2.96), Some(2.86), Some(2.80), Some(2.30), Some(2.30), Some(2.84), Some(2.82)],
        c_shifts: [Some(31.30), Some(30.82), Some(31.62), Some(32.10), Some(31.08), Some(30.70), Some(30.20), Some(30.16), Some(31.34), Some(31.34)],
        note: "Second N–CH₃ singlet; rotates slowly, two signals (2.96 & 2.88 in CDCl₃)",
    },
    SolventImpurity {
        name: "Ethyl Acetate (CH₃ t)", group: "CH₃", mult: "t, J=7",
        h_shifts: [Some(1.26), Some(1.17), Some(1.24), Some(1.25), Some(1.20), Some(1.19), Some(0.99), Some(1.00), Some(1.24), Some(1.24)],
        c_shifts: [Some(14.20), Some(14.08), Some(14.16), Some(14.20), Some(14.12), Some(14.05), Some(14.05), Some(14.00), Some(14.22), Some(14.20)],
        note: "Upfield triplet; paired with O–CH₂ (4.12 q) and CH₃CO (2.05 s)",
    },
    SolventImpurity {
        name: "Toluene (Ar–H)", group: "Ar–H", mult: "m/d",
        h_shifts: [Some(7.25), Some(7.25), Some(7.18), None, Some(7.20), Some(7.18), Some(7.09), Some(7.01), Some(7.22), Some(7.22)],
        c_shifts: [Some(129.20), Some(129.10), Some(129.30), None, Some(129.20), Some(129.10), Some(129.30), Some(129.10), Some(129.24), Some(129.24)],
        note: "Aromatic envelope; paired with Ar–CH₃ at 2.36 (CDCl₃)",
    },
    SolventImpurity {
        name: "BHT (Ar–CH₃)", group: "Ar–CH₃", mult: "s",
        h_shifts: [Some(2.27), Some(2.18), Some(2.22), None, Some(2.24), Some(2.21), Some(2.24), Some(2.22), Some(2.25), Some(2.22)],
        c_shifts: [Some(21.00), Some(20.60), Some(20.90), None, Some(20.80), Some(20.70), Some(21.10), Some(21.00), Some(21.00), Some(20.90)],
        note: "Ring methyl of BHT antioxidant; also t-Bu (1.43) & ArH (6.98)",
    },
    SolventImpurity {
        name: "BHT (OH)", group: "Ar–OH", mult: "s",
        h_shifts: [Some(5.01), Some(6.65), None, None, None, None, None, None, None, None],
        c_shifts: [None; 10],
        note: "Hindered phenolic OH; shifts strongly with solvent & concentration",
    },
    SolventImpurity {
        name: "Triethylamine (CH₃ t)", group: "CH₃", mult: "t, J=7",
        h_shifts: [Some(1.03), Some(0.93), Some(0.99), Some(0.96), Some(0.98), Some(0.97), Some(0.92), Some(0.93), Some(1.01), Some(0.99)],
        c_shifts: [Some(12.30), Some(12.10), Some(12.20), Some(12.10), Some(12.30), Some(12.30), Some(12.50), Some(12.40), Some(12.30), Some(12.30)],
        note: "Base/acid scavenger; paired with N–CH₂ quartet at 2.53 (CDCl₃)",
    },
    SolventImpurity {
        name: "DMAc — CH₃CO", group: "CH₃C=O", mult: "s",
        h_shifts: [Some(2.09), Some(1.96), Some(2.02), Some(2.06), Some(1.99), Some(1.97), Some(1.72), Some(1.74), Some(2.05), Some(2.01)],
        c_shifts: [Some(21.50), Some(21.00), Some(21.30), Some(21.60), Some(21.20), Some(21.10), Some(21.00), Some(21.00), Some(21.50), Some(21.50)],
        note: "N,N-Dimethylacetamide acetyl methyl; polymer solvent",
    },
    SolventImpurity {
        name: "DMAc — N(CH₃)₂ (#1)", group: "N(CH₃)₂", mult: "s",
        h_shifts: [Some(3.02), Some(2.94), Some(2.98), Some(3.00), Some(2.96), Some(2.92), Some(2.40), Some(2.40), Some(2.98), Some(2.94)],
        c_shifts: [Some(38.30), Some(37.80), Some(38.10), Some(38.60), Some(38.00), Some(37.70), Some(37.50), Some(37.40), Some(38.20), Some(38.10)],
        note: "One of the two magnetically inequivalent N-methyl singlets",
    },
    SolventImpurity {
        name: "DMAc — N(CH₃)₂ (#2)", group: "N(CH₃)₂", mult: "s",
        h_shifts: [Some(2.94), Some(2.78), Some(2.90), Some(2.94), Some(2.88), Some(2.84), Some(2.32), Some(2.32), Some(2.90), Some(2.86)],
        c_shifts: [Some(35.60), Some(35.10), Some(35.40), Some(35.90), Some(35.30), Some(35.00), Some(34.80), Some(34.70), Some(35.50), Some(35.40)],
        note: "Second N-methyl singlet (restricted rotation)",
    },
    SolventImpurity {
        name: "NMP — NCH₂", group: "N–CH₂", mult: "t",
        h_shifts: [Some(3.39), Some(3.31), Some(3.35), Some(3.38), Some(3.34), Some(3.30), Some(2.95), Some(2.95), Some(3.36), Some(3.32)],
        c_shifts: [Some(48.20), Some(47.80), Some(48.00), Some(48.40), Some(48.00), Some(47.70), Some(47.80), Some(47.70), Some(48.10), Some(48.00)],
        note: "N-Methyl-2-pyrrolidone ring N-methylene",
    },
    SolventImpurity {
        name: "NMP — NCH₃", group: "N–CH₃", mult: "s",
        h_shifts: [Some(2.85), Some(2.70), Some(2.78), Some(2.82), Some(2.76), Some(2.72), Some(2.38), Some(2.40), Some(2.82), Some(2.78)],
        c_shifts: [Some(29.80), Some(29.40), Some(29.60), Some(30.00), Some(29.60), Some(29.30), Some(29.50), Some(29.40), Some(29.80), Some(29.70)],
        note: "NMP N-methyl singlet; strong polymer solvent",
    },
    SolventImpurity {
        name: "NMP — CH₂CO", group: "CH₂C=O", mult: "t",
        h_shifts: [Some(2.38), Some(2.18), Some(2.28), Some(2.32), Some(2.26), Some(2.22), Some(1.95), Some(1.96), Some(2.34), Some(2.26)],
        c_shifts: [Some(30.90), Some(30.40), Some(30.60), Some(31.00), Some(30.60), Some(30.30), Some(30.40), Some(30.30), Some(30.90), Some(30.80)],
        note: "NMP α-carbonyl methylene",
    },
    SolventImpurity {
        name: "Anisole — OCH₃", group: "O–CH₃", mult: "s",
        h_shifts: [Some(3.79), Some(3.74), Some(3.76), None, Some(3.78), Some(3.76), Some(3.32), Some(3.32), Some(3.78), Some(3.74)],
        c_shifts: [Some(55.10), Some(55.00), Some(55.00), None, Some(55.20), Some(55.10), Some(54.60), Some(54.60), Some(55.20), Some(55.10)],
        note: "Anisole methoxy; aromatic 6.91–6.93",
    },
    SolventImpurity {
        name: "Anisole — ArH", group: "Ar–H", mult: "m",
        h_shifts: [Some(6.91), Some(6.93), Some(6.88), None, Some(6.90), Some(6.88), Some(6.78), Some(6.80), Some(6.90), Some(6.88)],
        c_shifts: [Some(129.50), Some(129.40), Some(129.20), None, Some(129.30), Some(129.20), Some(129.00), Some(128.90), Some(129.40), Some(129.30)],
        note: "Anisole aromatic protons",
    },
    SolventImpurity {
        name: "HFIP — CH (sept)", group: "CH(CF₃)₂", mult: "sept, J=6",
        h_shifts: [Some(4.41), Some(5.04), Some(4.38), Some(4.30), Some(4.40), Some(4.38), None, None, Some(4.42), Some(4.36)],
        c_shifts: [Some(68.00), Some(66.50), Some(67.50), Some(67.00), Some(67.80), Some(67.60), None, None, Some(67.80), Some(67.60)],
        note: "Hexafluoroisopropanol — strongly H-bond donating NMR solvent",
    },
    SolventImpurity {
        name: "AIBN (Initiator) — C(CH₃)₂", group: "C(CH₃)₂CN", mult: "s",
        h_shifts: [Some(1.70), Some(1.68), Some(1.72), None, Some(1.70), Some(1.68), Some(1.68), Some(1.66), Some(1.70), Some(1.68)],
        c_shifts: [Some(25.00), Some(24.80), Some(25.10), None, Some(25.00), Some(24.90), Some(24.90), Some(24.80), Some(25.00), Some(24.90)],
        note: "Azobisisobutyronitrile radical initiator fragment (12H singlet)",
    },
    SolventImpurity {
        name: "MEK (2-Butanone) — CH₂", group: "CH₂", mult: "q, J=7",
        h_shifts: [Some(2.45), Some(2.43), Some(2.44), Some(2.46), Some(2.44), Some(2.42), Some(2.28), Some(2.28), Some(2.44), Some(2.42)],
        c_shifts: [Some(36.80), Some(36.20), Some(36.50), Some(36.90), Some(36.50), Some(36.20), Some(36.40), Some(36.30), Some(36.80), Some(36.70)],
        note: "Methyl ethyl ketone α-methylene",
    },
    SolventImpurity {
        name: "MEK (2-Butanone) — CH₃", group: "CH₃", mult: "s / t",
        h_shifts: [Some(2.14), Some(2.07), Some(2.10), Some(2.12), Some(2.09), Some(2.07), Some(1.96), Some(1.96), Some(2.12), Some(2.09)],
        c_shifts: [Some(29.10), Some(28.60), Some(28.90), Some(29.20), Some(28.90), Some(28.60), Some(28.80), Some(28.70), Some(29.10), Some(29.00)],
        note: "Acetyl methyl singlet plus terminal CH₃ triplet at 1.05",
    },
    SolventImpurity {
        name: "Chlorobenzene — ArH", group: "Ar–H", mult: "m",
        h_shifts: [Some(7.31), Some(7.41), Some(7.28), None, Some(7.36), Some(7.30), Some(6.93), Some(6.94), Some(7.32), Some(7.34)],
        c_shifts: [Some(129.80), Some(130.20), Some(129.50), None, Some(129.80), Some(129.50), Some(129.20), Some(129.10), Some(129.90), Some(129.90)],
        note: "High-boiling halogenated solvent",
    },
    SolventImpurity {
        name: "MEHQ (Monomer Inhibitor) — OCH₃", group: "O–CH₃", mult: "s",
        h_shifts: [Some(3.75), Some(3.68), Some(3.72), None, Some(3.72), Some(3.70), Some(3.30), Some(3.32), Some(3.74), Some(3.70)],
        c_shifts: [Some(55.90), Some(55.60), Some(55.70), None, Some(55.80), Some(55.70), Some(55.30), Some(55.20), Some(55.90), Some(55.80)],
        note: "Monomethyl ether of hydroquinone; polymerization inhibitor in MMA",
    },
    SolventImpurity {
        name: "MEHQ (Monomer Inhibitor) — ArH", group: "Ar–H", mult: "m",
        h_shifts: [Some(6.76), Some(6.69), Some(6.72), None, Some(6.74), Some(6.70), Some(6.80), Some(6.82), Some(6.74), Some(6.70)],
        c_shifts: [Some(114.50), Some(114.20), Some(114.30), None, Some(114.40), Some(114.20), Some(114.60), Some(114.50), Some(114.50), Some(114.40)],
        note: "MEHQ aromatic protons; also OH (4.5–5.5, br)",
    },
    SolventImpurity {
        name: "Styrene Monomer — =CH", group: "Ar–CH=", mult: "dd",
        h_shifts: [Some(6.72), Some(6.75), Some(6.70), None, Some(6.72), Some(6.70), Some(6.60), Some(6.62), Some(6.72), Some(6.72)],
        c_shifts: [Some(137.40), Some(137.30), Some(137.20), None, Some(137.30), Some(137.20), Some(137.60), Some(137.50), Some(137.40), Some(137.40)],
        note: "Unreacted styrene monomer; disappears on polymerization",
    },
    SolventImpurity {
        name: "Styrene Monomer — =CH₂ (trans)", group: "=CH₂", mult: "d",
        h_shifts: [Some(5.75), Some(5.78), Some(5.72), None, Some(5.74), Some(5.72), Some(5.62), Some(5.64), Some(5.74), Some(5.74)],
        c_shifts: [Some(112.90), Some(112.80), Some(112.70), None, Some(112.80), Some(112.70), Some(112.40), Some(112.50), Some(112.90), Some(112.90)],
        note: "Styrene terminal vinyl proton, trans to phenyl",
    },
    SolventImpurity {
        name: "Styrene Monomer — =CH₂ (cis)", group: "=CH₂", mult: "d",
        h_shifts: [Some(5.24), Some(5.27), Some(5.22), None, Some(5.24), Some(5.22), Some(5.12), Some(5.14), Some(5.24), Some(5.24)],
        c_shifts: [Some(112.90), Some(112.80), Some(112.70), None, Some(112.80), Some(112.70), Some(112.40), Some(112.50), Some(112.90), Some(112.90)],
        note: "Styrene terminal vinyl proton, cis to phenyl",
    },
    SolventImpurity {
        name: "MMA Monomer — =CH₂ (cis)", group: "=CH₂", mult: "dd",
        h_shifts: [Some(6.10), Some(6.12), Some(6.08), None, Some(6.08), Some(6.06), Some(5.98), Some(6.00), Some(6.08), Some(6.06)],
        c_shifts: [Some(125.20), Some(125.10), Some(125.00), None, Some(125.10), Some(125.00), Some(125.40), Some(125.30), Some(125.20), Some(125.20)],
        note: "Unreacted MMA monomer; KEY conversion marker",
    },
    SolventImpurity {
        name: "MMA Monomer — =CH₂ (trans)", group: "=CH₂", mult: "dd",
        h_shifts: [Some(5.56), Some(5.58), Some(5.54), None, Some(5.54), Some(5.52), Some(5.44), Some(5.46), Some(5.56), Some(5.54)],
        c_shifts: [Some(125.20), Some(125.10), Some(125.00), None, Some(125.10), Some(125.00), Some(125.40), Some(125.30), Some(125.20), Some(125.20)],
        note: "Second MMA vinyl proton",
    },
    SolventImpurity {
        name: "Pyridine — ArH (α)", group: "Ar–H", mult: "d",
        h_shifts: [Some(8.62), Some(8.58), Some(8.60), None, Some(8.60), Some(8.56), Some(8.20), Some(8.22), Some(8.60), Some(8.60)],
        c_shifts: [Some(149.90), Some(149.80), Some(149.60), None, Some(149.80), Some(149.60), Some(149.20), Some(149.30), Some(149.90), Some(149.90)],
        note: "Pyridine α-proton, most deshielded",
    },
    SolventImpurity {
        name: "Pyridine — ArH (γ)", group: "Ar–H", mult: "t",
        h_shifts: [Some(7.68), Some(7.79), Some(7.70), None, Some(7.70), Some(7.66), Some(7.10), Some(7.12), Some(7.68), Some(7.68)],
        c_shifts: [Some(136.00), Some(136.10), Some(135.80), None, Some(135.90), Some(135.80), Some(135.20), Some(135.30), Some(136.00), Some(136.00)],
        note: "Pyridine γ-proton",
    },
    SolventImpurity {
        name: "Pyridine — ArH (β)", group: "Ar–H", mult: "m",
        h_shifts: [Some(7.29), Some(7.39), Some(7.30), None, Some(7.32), Some(7.28), Some(6.98), Some(7.00), Some(7.30), Some(7.30)],
        c_shifts: [Some(123.60), Some(123.70), Some(123.50), None, Some(123.60), Some(123.50), Some(123.20), Some(123.30), Some(123.60), Some(123.60)],
        note: "Pyridine β-protons",
    },
    SolventImpurity {
        name: "2-Propanol (IPA) — CH", group: "O–CH", mult: "sept, J=6",
        h_shifts: [Some(4.04), Some(3.78), Some(3.92), None, Some(3.94), Some(3.90), Some(3.64), Some(3.66), Some(3.98), Some(3.90)],
        c_shifts: [Some(64.50), Some(63.50), Some(64.00), None, Some(64.20), Some(64.00), Some(64.50), Some(64.40), Some(64.60), Some(64.40)],
        note: "Isopropanol methine septet; paired with CH₃ doublet",
    },
    SolventImpurity {
        name: "2-Propanol (IPA) — CH₃", group: "CH₃", mult: "d, J=6",
        h_shifts: [Some(1.22), Some(1.04), Some(1.16), None, Some(1.14), Some(1.12), Some(0.96), Some(0.98), Some(1.18), Some(1.14)],
        c_shifts: [Some(25.20), Some(25.00), Some(25.10), None, Some(25.10), Some(25.00), Some(25.30), Some(25.20), Some(25.20), Some(25.10)],
        note: "Isopropanol methyl doublet",
    },
    SolventImpurity {
        name: "1,2-Dichloroethane (DCE)", group: "ClCH₂CH₂Cl", mult: "s",
        h_shifts: [Some(3.73), Some(3.90), Some(3.80), None, Some(3.82), Some(3.78), Some(3.22), Some(3.24), Some(3.80), Some(3.78)],
        c_shifts: [Some(43.50), Some(43.80), Some(43.60), None, Some(43.60), Some(43.40), Some(43.90), Some(43.80), Some(43.60), Some(43.50)],
        note: "DCE polymer solvent; 4H singlet",
    },
    SolventImpurity {
        name: "Cyclohexane", group: "CH₂", mult: "s",
        h_shifts: [Some(1.43), Some(1.40), Some(1.42), None, Some(1.42), Some(1.40), Some(1.38), Some(1.38), Some(1.42), Some(1.42)],
        c_shifts: [Some(27.10), Some(26.90), Some(27.00), None, Some(27.00), Some(26.90), Some(27.10), Some(27.00), Some(27.10), Some(27.10)],
        note: "Reference aliphatic solvent",
    },
    SolventImpurity {
        name: "Benzene", group: "Ar–H", mult: "s",
        h_shifts: [Some(7.36), Some(7.37), Some(7.34), None, Some(7.36), Some(7.34), Some(7.16), Some(7.16), Some(7.36), Some(7.36)],
        c_shifts: [Some(128.37), Some(128.30), Some(128.50), None, Some(128.40), Some(128.30), Some(128.00), Some(128.00), Some(128.40), Some(128.40)],
        note: "Non-deuterated benzene contaminant",
    },
    SolventImpurity {
        name: "Acetonitrile (CH₃)", group: "CH₃CN", mult: "s",
        h_shifts: [Some(2.10), Some(2.07), Some(2.08), Some(2.10), Some(2.08), Some(2.06), Some(1.86), Some(1.88), Some(2.08), Some(2.06)],
        c_shifts: [Some(1.32), Some(1.10), Some(1.20), Some(1.40), Some(1.20), Some(1.10), Some(1.50), Some(1.40), Some(1.30), Some(1.20)],
        note: "Acetonitrile impurity / solvent",
    },
];

// ============================================================================
// 16 MAJOR POLYMERS FINGERPRINT LIBRARY
// ============================================================================
pub const POLYMERS_DB: &[PolymerRef] = &[
    PolymerRef {
        short: "PMMA", name: "Poly(methyl methacrylate)", repeat_unit: "-[CH₂-C(CH₃)(COOCH₃)]-", m0: 100.12, tg_tm: "Tg ≈ 105 °C (amorphous)",
        ftir_peaks: &[(2995.0, "medium CH₃ asym"), (2950.0, "medium CH₂ asym"), (2844.0, "weak CH₂ sym"), (1728.0, "very strong ester C=O"), (1638.0, "weak C=C vinyl residual"), (1484.0, "weak CH₂ scissor"), (1448.0, "medium CH₃ bend"), (1388.0, "weak CH₃ umbrella"), (1270.0, "medium C-O-C"), (1242.0, "strong C-O-C"), (1192.0, "strong C-O-C"), (1148.0, "very strong C-O-C"), (988.0, "medium CH₂ rock"), (965.0, "weak CH₂ rock"), (842.0, "weak CH₂ rock"), (750.0, "medium C=O OOP")],
        hnmr_peaks: &[(3.60, "OCH₃"), (1.85, "CH₂ backbone"), (1.21, "α-CH₃ mm"), (1.02, "α-CH₃ mr"), (0.85, "α-CH₃ rr")],
        cnmr_peaks: &[(177.8, "C=O"), (54.4, "CH₂"), (51.8, "OCH₃"), (44.6, "Cq"), (18.7, "α-CH₃")],
    },
    PolymerRef {
        short: "PS", name: "Polystyrene", repeat_unit: "-[CH₂-CH(Ph)]-", m0: 104.15, tg_tm: "Tg ≈ 100 °C (amorphous)",
        ftir_peaks: &[(3082.0, "weak Ar C-H"), (3060.0, "weak Ar C-H"), (3026.0, "medium Ar C-H"), (2924.0, "strong CH₂ asym"), (2850.0, "medium CH₂ sym"), (1942.0, "weak Ar overtone"), (1870.0, "weak Ar overtone"), (1802.0, "weak Ar overtone"), (1601.0, "medium Ar ring"), (1583.0, "weak Ar ring"), (1493.0, "strong Ar ring"), (1452.0, "strong CH₂/Ar"), (1372.0, "weak CH bend"), (1182.0, "weak C-H"), (1154.0, "weak C-H"), (1069.0, "weak C-H"), (1028.0, "medium Ar C-H"), (906.0, "weak Ar C-H"), (842.0, "weak Ar C-H"), (756.0, "strong Ar OOP"), (698.0, "very strong phenyl bend"), (540.0, "medium ring deformation")],
        hnmr_peaks: &[(7.08, "m/p Ar-H"), (6.55, "o Ar-H"), (1.85, "CH"), (1.43, "CH₂")],
        cnmr_peaks: &[(145.3, "Ar ipso"), (128.0, "Ar o/m"), (125.6, "Ar p"), (44.0, "CH₂"), (40.5, "CH")],
    },
    PolymerRef {
        short: "PE", name: "Polyethylene (HDPE/LDPE)", repeat_unit: "-[CH₂-CH₂]-", m0: 28.05, tg_tm: "Tg ≈ -120 °C, Tm ≈ 115–135 °C",
        ftir_peaks: &[(2916.0, "CH₂ asym"), (2848.0, "CH₂ sym"), (1472.0, "CH₂ cryst"), (1462.0, "CH₂ amorph"), (730.0, "CH₂ rock cryst"), (720.0, "CH₂ rock")],
        hnmr_peaks: &[(1.29, "CH₂ backbone")],
        cnmr_peaks: &[(30.0, "CH₂ backbone")],
    },
    PolymerRef {
        short: "PP", name: "Polypropylene", repeat_unit: "-[CH₂-CH(CH₃)]-", m0: 42.08, tg_tm: "Tg ≈ -10 °C, Tm ≈ 165 °C",
        ftir_peaks: &[(2950.0, "CH₃ asym"), (2917.0, "CH₂ asym"), (2868.0, "CH₃ sym"), (1456.0, "CH₂/CH₃"), (1376.0, "CH₃ umbrella"), (998.0, "isotactic helix"), (973.0, "CH₃ rock"), (841.0, "isotactic")],
        hnmr_peaks: &[(1.62, "CH"), (1.28, "CH₂"), (0.88, "CH₃")],
        cnmr_peaks: &[(46.5, "CH₂"), (28.8, "CH"), (21.6, "CH₃ (isotactic)")],
    },
    PolymerRef {
        short: "PET", name: "Poly(ethylene terephthalate)", repeat_unit: "-[O-CH₂CH₂-O-CO-C₆H₄-CO]-", m0: 192.17, tg_tm: "Tg ≈ 75 °C, Tm ≈ 255 °C",
        ftir_peaks: &[(1715.0, "Aryl ester C=O"), (1408.0, "Ar ring"), (1340.0, "CH₂ wag trans"), (1245.0, "C-O-C"), (1096.0, "C-O-C"), (872.0, "p-Ar OOP"), (725.0, "C=O/Ar OOP")],
        hnmr_peaks: &[(8.10, "Ar-H"), (4.70, "O-CH₂CH₂-O")],
        cnmr_peaks: &[(166.0, "C=O"), (134.0, "Ar Cq"), (129.6, "Ar CH"), (63.2, "O-CH₂")],
    },
    PolymerRef {
        short: "PA66", name: "Nylon-6,6 (Polyamide)", repeat_unit: "-[NH-(CH₂)₆-NH-CO-(CH₂)₄-CO]-", m0: 226.32, tg_tm: "Tg ≈ 50 °C, Tm ≈ 265 °C",
        ftir_peaks: &[(3300.0, "N-H H-bonded"), (2932.0, "CH₂ asym"), (2858.0, "CH₂ sym"), (1640.0, "Amide I C=O"), (1540.0, "Amide II N-H"), (1275.0, "Amide III"), (690.0, "N-H OOP")],
        hnmr_peaks: &[(7.80, "Amide NH"), (3.20, "N-CH₂"), (2.20, "CH₂-CO"), (1.55, "CH₂")],
        cnmr_peaks: &[(173.5, "C=O"), (39.5, "N-CH₂"), (36.2, "CH₂-CO"), (26.5, "CH₂")],
    },
    PolymerRef {
        short: "PVC", name: "Poly(vinyl chloride)", repeat_unit: "-[CH₂-CH(Cl)]-", m0: 62.50, tg_tm: "Tg ≈ 82 °C",
        ftir_peaks: &[(2970.0, "CH-Cl stretch"), (2912.0, "CH₂"), (1427.0, "CH₂ bend"), (1254.0, "CH-Cl rock"), (960.0, "CH₂ rock"), (690.0, "C-Cl"), (615.0, "C-Cl syndio")],
        hnmr_peaks: &[(4.45, "CH-Cl"), (2.15, "CH₂")],
        cnmr_peaks: &[(56.5, "CH-Cl"), (46.0, "CH₂")],
    },
    PolymerRef {
        short: "PDMS", name: "Poly(dimethylsiloxane)", repeat_unit: "-[Si(CH₃)₂-O]-", m0: 74.15, tg_tm: "Tg ≈ -125 °C (elastomer)",
        ftir_peaks: &[(2962.0, "CH₃ asym"), (1260.0, "Si-CH₃ sharp"), (1080.0, "Si-O-Si"), (1015.0, "Si-O-Si"), (795.0, "Si-CH₃ rock")],
        hnmr_peaks: &[(0.07, "Si-CH₃")],
        cnmr_peaks: &[(1.0, "Si-CH₃")],
    },
    PolymerRef {
        short: "PEG", name: "Poly(ethylene glycol) / PEO", repeat_unit: "-[CH₂-CH₂-O]-", m0: 44.05, tg_tm: "Tg ≈ -60 °C, Tm ≈ 65 °C",
        ftir_peaks: &[(3450.0, "OH end-group"), (2885.0, "CH₂ stretch"), (1467.0, "CH₂ scissor"), (1342.0, "CH₂ wag"), (1105.0, "C-O-C very strong"), (962.0, "CH₂ rock"), (842.0, "C-O-C")],
        hnmr_peaks: &[(3.64, "O-CH₂CH₂-O")],
        cnmr_peaks: &[(70.6, "O-CH₂")],
    },
    PolymerRef {
        short: "PLA", name: "Poly(lactic acid)", repeat_unit: "-[O-CH(CH₃)-CO]-", m0: 72.06, tg_tm: "Tg ≈ 60 °C, Tm ≈ 175 °C",
        ftir_peaks: &[(2995.0, "CH₃"), (1750.0, "Ester C=O"), (1452.0, "CH₃ bend"), (1182.0, "C-O-C"), (1085.0, "C-O-C"), (870.0, "amorphous"), (755.0, "crystalline")],
        hnmr_peaks: &[(5.16, "CH quartet"), (1.58, "CH₃ doublet")],
        cnmr_peaks: &[(169.6, "C=O"), (69.0, "CH"), (16.7, "CH₃")],
    },
    PolymerRef {
        short: "PCL", name: "Poly(ε-caprolactone)", repeat_unit: "-[O-(CH₂)₅-CO]-", m0: 114.14, tg_tm: "Tg ≈ -60 °C, Tm ≈ 60 °C",
        ftir_peaks: &[(2945.0, "CH₂ asym"), (2865.0, "CH₂ sym"), (1725.0, "Ester C=O"), (1294.0, "C-O/C-C cryst"), (1240.0, "C-O-C asym"), (1175.0, "C-O-C sym"), (732.0, "CH₂ rock")],
        hnmr_peaks: &[(4.06, "O-CH₂"), (2.31, "CH₂-CO"), (1.65, "CH₂"), (1.38, "CH₂")],
        cnmr_peaks: &[(173.6, "C=O"), (64.1, "O-CH₂"), (34.1, "CH₂-CO"), (28.4, "CH₂"), (25.5, "CH₂")],
    },
    PolymerRef {
        short: "PU", name: "Polyurethane", repeat_unit: "-[NH-CO-O-R-O-CO-NH-R']-", m0: 250.0, tg_tm: "Segmented elastomer/foam",
        ftir_peaks: &[(3330.0, "N-H urethane"), (2270.0, "Unreacted NCO!"), (1730.0, "C=O free"), (1705.0, "C=O H-bonded"), (1530.0, "Amide II"), (1220.0, "C-O-C"), (1100.0, "Polyether C-O-C")],
        hnmr_peaks: &[(7.35, "Ar/NH"), (4.10, "O-CH₂-OCONH"), (3.60, "Polyether O-CH₂")],
        cnmr_peaks: &[(154.0, "Urethane C=O"), (70.5, "Polyether C-O"), (64.0, "O-CH₂")],
    },
    PolymerRef {
        short: "PC", name: "Polycarbonate (BPA-PC)", repeat_unit: "-[O-C₆H₄-C(CH₃)₂-C₆H₄-O-CO]-", m0: 254.28, tg_tm: "Tg ≈ 147 °C (amorphous)",
        ftir_peaks: &[(1770.0, "Carbonate C=O"), (1505.0, "p-Ar ring"), (1225.0, "C-O-C"), (1190.0, "C-O"), (1160.0, "C-O"), (1014.0, "Ar"), (830.0, "p-Ar OOP")],
        hnmr_peaks: &[(7.25, "Ar-H"), (7.16, "Ar-H"), (1.68, "C(CH₃)₂")],
        cnmr_peaks: &[(152.1, "Carbonate C=O"), (148.9, "Ar C-O"), (128.0, "Ar CH"), (120.3, "Ar CH"), (42.5, "Cq"), (30.9, "CH₃")],
    },
    PolymerRef {
        short: "PTFE", name: "Polytetrafluoroethylene (Teflon)", repeat_unit: "-[CF₂-CF₂]-", m0: 100.02, tg_tm: "Tm ≈ 327 °C",
        ftir_peaks: &[(1200.0, "CF₂ asym"), (1146.0, "CF₂ sym"), (638.0, "CF₂ wag"), (554.0, "CF₂ bend"), (503.0, "CF₂ rock")],
        hnmr_peaks: &[],
        cnmr_peaks: &[(112.0, "CF₂ (solid-state)")],
    },
    PolymerRef {
        short: "Epoxy", name: "DGEBA Epoxy Resin", repeat_unit: "Bisphenol-A diglycidyl ether network", m0: 340.41, tg_tm: "Thermoset network",
        ftir_peaks: &[(3400.0, "OH (increases on cure)"), (1608.0, "Ar ring"), (1508.0, "p-Ar ring"), (1245.0, "Ar-O-C"), (1035.0, "C-O-C"), (915.0, "Oxirane ring (tracks cure!)"), (830.0, "p-Ar OOP")],
        hnmr_peaks: &[(7.12, "Ar-H"), (6.82, "Ar-H"), (3.34, "Oxirane CH"), (2.88, "Oxirane CH₂"), (1.63, "C(CH₃)₂")],
        cnmr_peaks: &[(156.4, "Ar C-O"), (127.8, "Ar CH"), (114.0, "Ar CH"), (68.8, "O-CH₂"), (50.2, "Oxirane CH"), (44.7, "Oxirane CH₂")],
    },
    PolymerRef {
        short: "PAN", name: "Polyacrylonitrile", repeat_unit: "-[CH₂-CH(CN)]-", m0: 53.06, tg_tm: "Tg ≈ 95 °C",
        ftir_peaks: &[(2940.0, "CH₂ stretch"), (2243.0, "C≡N sharp nitrile"), (1454.0, "CH₂ scissor"), (1360.0, "CH bend"), (1070.0, "C-CN")],
        hnmr_peaks: &[(3.15, "CH-CN (DMSO)"), (2.06, "CH₂ (DMSO)")],
        cnmr_peaks: &[(120.1, "C≡N"), (33.0, "CH₂"), (27.4, "CH-CN")],
    },
    PolymerRef {
        short: "PVAc", name: "Poly(vinyl acetate)", repeat_unit: "-[CH₂-CH(OCOCH₃)]-", m0: 86.09, tg_tm: "Tg ≈ 32 °C",
        ftir_peaks: &[(1738.0, "Acetate C=O"), (1435.0, "CH₃/CH₂"), (1372.0, "CH₃ umbrella"), (1235.0, "C-O-C very strong"), (1022.0, "C-O")],
        hnmr_peaks: &[(4.87, "CH-OAc"), (2.04, "OCOCH₃"), (1.78, "CH₂")],
        cnmr_peaks: &[(170.5, "C=O"), (67.0, "CH-O"), (39.5, "CH₂"), (21.1, "CH₃")],
    },

    PolymerRef {
        short: "PVDF", name: "Poly(vinylidene fluoride)", repeat_unit: "-[CH₂-CF₂]-", m0: 64.03, tg_tm: "Tg ≈ -40 °C, Tm ≈ 170 °C (α/β polymorphs)",
        ftir_peaks: &[(3025.0, "CH₂ asym"), (2985.0, "CH₂ sym"), (1402.0, "CH₂ scissor"), (1275.0, "β CF₂ asym"), (1210.0, "CF₂ asym"), (1170.0, "CF₂ sym"), (976.0, "α CH₂ rock"), (840.0, "β-phase CF₂ rock"), (795.0, "α CF₂ rock"), (763.0, "α-phase 763"), (510.0, "β CF₂ bend")],
        hnmr_peaks: &[(2.90, "CH₂ (head-tail)"), (2.30, "CH₂ (defect)")],
        cnmr_peaks: &[(120.0, "CF₂"), (33.0, "CH₂")],
    },
    PolymerRef {
        short: "PEEK", name: "Poly(ether ether ketone)", repeat_unit: "-[O-C₆H₄-O-C₆H₄-CO-C₆H₄]-", m0: 288.30, tg_tm: "Tg ≈ 143 °C, Tm ≈ 343 °C",
        ftir_peaks: &[(1655.0, "Diaryl ketone C=O"), (1595.0, "Ar ring"), (1490.0, "Ar ring"), (1230.0, "Ar-O-Ar asym"), (1185.0, "Ar-O-Ar sym"), (1010.0, "Ar"), (835.0, "p-Ar OOP")],
        hnmr_peaks: &[(7.75, "Ar-H (ketone)"), (7.05, "Ar-H (ether)")],
        cnmr_peaks: &[(195.0, "C=O"), (161.0, "Ar C-O"), (130.0, "Ar CH"), (118.0, "Ar CH")],
    },
    PolymerRef {
        short: "PI", name: "Polyimide (Kapton PMDA-ODA)", repeat_unit: "-[N(CO)₂-C₆H₂-O-C₆H₄-O-]", m0: 382.30, tg_tm: "Tg ≈ 360–410 °C (thermostable film)",
        ftir_peaks: &[(1778.0, "Imide I sym C=O"), (1720.0, "Imide II asym C=O"), (1500.0, "Ar ring"), (1375.0, "Imide III C-N"), (1230.0, "Ar-O-Ar"), (725.0, "Imide IV ring bend")],
        hnmr_peaks: &[(7.90, "Ar-H (PMDA)"), (7.20, "Ar-H (ODA)")],
        cnmr_peaks: &[(166.0, "Imide C=O"), (156.0, "Ar C-O"), (128.0, "Ar CH")],
    },
    PolymerRef {
        short: "PSU", name: "Polysulfone (Bisphenol-A PSU)", repeat_unit: "-[O-C₆H₄-C(CH₃)₂-C₆H₄-O-C₆H₄-SO₂-C₆H₄]-", m0: 442.50, tg_tm: "Tg ≈ 185 °C (amorphous)",
        ftir_peaks: &[(1585.0, "Ar ring"), (1500.0, "Ar ring"), (1320.0, "SO₂ asym"), (1290.0, "Ar-O-Ar"), (1240.0, "Ar-O-C"), (1150.0, "SO₂ sym"), (1015.0, "Ar"), (835.0, "p-Ar OOP")],
        hnmr_peaks: &[(7.85, "Ar-H (sulfone)"), (7.25, "Ar-H"), (6.95, "Ar-H (ether)"), (1.70, "C(CH₃)₂")],
        cnmr_peaks: &[(162.0, "Ar C-O"), (152.0, "Ar C-SO₂"), (128.0, "Ar CH"), (42.5, "Cq"), (31.0, "CH₃")],
    },
    PolymerRef {
        short: "PES", name: "Poly(ether sulfone)", repeat_unit: "-[O-C₆H₄-SO₂-C₆H₄]-", m0: 232.26, tg_tm: "Tg ≈ 220 °C",
        ftir_peaks: &[(1575.0, "Ar ring"), (1485.0, "Ar ring"), (1320.0, "SO₂ asym"), (1240.0, "Ar-O-Ar"), (1150.0, "SO₂ sym"), (1010.0, "Ar"), (840.0, "p-Ar OOP")],
        hnmr_peaks: &[(7.90, "Ar-H (sulfone)"), (7.05, "Ar-H (ether)")],
        cnmr_peaks: &[(161.0, "Ar C-O"), (136.0, "Ar C-SO₂"), (129.0, "Ar CH")],
    },
    PolymerRef {
        short: "POM", name: "Polyoxymethylene (Acetal)", repeat_unit: "-[CH₂-O]-", m0: 30.03, tg_tm: "Tg ≈ -60 °C, Tm ≈ 175 °C",
        ftir_peaks: &[(2960.0, "CH₂ asym"), (2920.0, "CH₂ sym"), (1470.0, "CH₂ scissor"), (1240.0, "CH₂ wag"), (1090.0, "C-O-C asym"), (940.0, "C-O-C sym"), (620.0, "C-O bend")],
        hnmr_peaks: &[(4.80, "O-CH₂-O")],
        cnmr_peaks: &[(95.0, "O-CH₂-O")],
    },
    PolymerRef {
        short: "PA6", name: "Nylon-6 (Polycaprolactam)", repeat_unit: "-[NH-(CH₂)₅-CO]-", m0: 113.16, tg_tm: "Tg ≈ 50 °C, Tm ≈ 220 °C",
        ftir_peaks: &[(3300.0, "N-H H-bonded"), (2930.0, "CH₂ asym"), (2855.0, "CH₂ sym"), (1640.0, "Amide I C=O"), (1540.0, "Amide II N-H"), (1265.0, "Amide III"), (690.0, "N-H OOP")],
        hnmr_peaks: &[(7.80, "Amide NH"), (3.20, "N-CH₂"), (2.20, "CH₂-CO"), (1.55, "CH₂")],
        cnmr_peaks: &[(174.0, "C=O"), (39.5, "N-CH₂"), (36.2, "CH₂-CO"), (26.5, "CH₂")],
    },
    PolymerRef {
        short: "PBT", name: "Poly(butylene terephthalate)", repeat_unit: "-[O-(CH₂)₄-O-CO-C₆H₄-CO]-", m0: 220.22, tg_tm: "Tg ≈ 45 °C, Tm ≈ 225 °C",
        ftir_peaks: &[(1715.0, "Aryl ester C=O"), (1470.0, "CH₂ scissor"), (1395.0, "Ar ring"), (1270.0, "C-O-C"), (1100.0, "C-O-C"), (1015.0, "Ar"), (815.0, "p-Ar OOP")],
        hnmr_peaks: &[(8.10, "Ar-H"), (4.40, "O-CH₂"), (1.95, "CH₂")],
        cnmr_peaks: &[(166.0, "C=O"), (134.0, "Ar Cq"), (129.6, "Ar CH"), (64.5, "O-CH₂"), (25.6, "CH₂")],
    },
    PolymerRef {
        short: "ABS", name: "Acrylonitrile-Butadiene-Styrene", repeat_unit: "SAN grafted onto polybutadiene", m0: 211.00, tg_tm: "Tg ≈ 105 °C (SAN) + -80 °C (PBD)",
        ftir_peaks: &[(3060.0, "Ar C-H"), (3026.0, "Ar C-H"), (2920.0, "CH₂ asym"), (2238.0, "C≡N nitrile"), (1601.0, "Ar ring"), (1493.0, "Ar ring"), (966.0, "trans C=C"), (756.0, "Ar OOP"), (698.0, "Phenyl bend")],
        hnmr_peaks: &[(7.30, "Ar-H (SAN)"), (7.10, "Ar-H"), (5.40, "PBD CH=CH"), (3.10, "CH-CN"), (1.90, "CH₂")],
        cnmr_peaks: &[(145.3, "PS ipso"), (128.0, "Ar CH"), (120.1, "C≡N"), (44.0, "CH₂"), (40.5, "CH")],
    },
    PolymerRef {
        short: "PNIPAM", name: "Poly(N-isopropylacrylamide)", repeat_unit: "-[CH₂-CH(CONH-iPr)]-", m0: 113.16, tg_tm: "Tg ≈ 130 °C, LCST ≈ 32 °C",
        ftir_peaks: &[(3300.0, "N-H stretch"), (2970.0, "CH₃ asym"), (2935.0, "CH₂ asym"), (1650.0, "Amide I C=O"), (1550.0, "Amide II N-H"), (1465.0, "CH₃ bend"), (1170.0, "C-N")],
        hnmr_peaks: &[(7.20, "Amide NH"), (3.90, "CH(CH₃)₂"), (2.10, "CH backbone"), (1.55, "CH₂ backbone"), (1.10, "CH(CH₃)₂")],
        cnmr_peaks: &[(175.0, "C=O"), (42.0, "CH"), (35.0, "CH₂"), (22.5, "CH₃")],
    },
    PolymerRef {
        short: "PHEMA", name: "Poly(2-hydroxyethyl methacrylate)", repeat_unit: "-[CH₂-C(CH₃)(COOCH₂CH₂OH)]-", m0: 130.14, tg_tm: "Tg ≈ 55 °C (hydrogel)",
        ftir_peaks: &[(3400.0, "O-H H-bonded"), (2950.0, "CH₃/CH₂"), (1725.0, "Ester C=O"), (1450.0, "CH₃ bend"), (1160.0, "C-O-C"), (1075.0, "C-OH"), (750.0, "C=O OOP")],
        hnmr_peaks: &[(4.85, "OH"), (4.10, "O-CH₂"), (3.80, "CH₂-OH"), (1.85, "CH₂ backbone"), (1.10, "α-CH₃")],
        cnmr_peaks: &[(177.0, "C=O"), (67.0, "O-CH₂"), (60.0, "CH₂-OH"), (54.0, "Cq"), (45.0, "CH₂"), (18.0, "CH₃")],
    },
    PolymerRef {
        short: "PAA", name: "Poly(acrylic acid)", repeat_unit: "-[CH₂-CH(COOH)]-", m0: 72.06, tg_tm: "Tg ≈ 106 °C",
        ftir_peaks: &[(3300.0, "O-H acid broad"), (2940.0, "CH₂ asym"), (2860.0, "CH₂ sym"), (1705.0, "Acid C=O"), (1450.0, "CH₂ scissor"), (1240.0, "C-O"), (800.0, "CH₂ rock")],
        hnmr_peaks: &[(12.0, "COOH"), (2.30, "CH backbone"), (1.65, "CH₂ backbone")],
        cnmr_peaks: &[(178.0, "COOH"), (42.0, "CH"), (35.0, "CH₂")],
    },
    PolymerRef {
        short: "PVP", name: "Poly(vinylpyrrolidone)", repeat_unit: "-[CH₂-CH(NC₄H₆O)]-", m0: 111.14, tg_tm: "Tg ≈ 175 °C (amorphous)",
        ftir_peaks: &[(2950.0, "CH₂ asym"), (1670.0, "Lactam C=O"), (1460.0, "CH₂ scissor"), (1420.0, "CH₂ bend"), (1290.0, "C-N lactam"), (1015.0, "C-N"), (650.0, "N-C=O bend")],
        hnmr_peaks: &[(3.70, "N-CH₂"), (3.20, "N-CH"), (2.20, "CH/CH₂ backbone"), (1.60, "CH₂")],
        cnmr_peaks: &[(175.0, "C=O"), (45.0, "CH"), (35.0, "CH₂"), (31.0, "N-CH₂")],
    },
    PolymerRef {
        short: "PVA", name: "Poly(vinyl alcohol)", repeat_unit: "-[CH₂-CH(OH)]-", m0: 44.05, tg_tm: "Tg ≈ 85 °C, Tm ≈ 230 °C",
        ftir_peaks: &[(3400.0, "O-H H-bonded"), (2940.0, "CH₂ asym"), (2890.0, "CH₂ sym"), (1430.0, "CH₂ scissor"), (1330.0, "CH-OH bend"), (1090.0, "C-O stretch"), (850.0, "CH₂ rock")],
        hnmr_peaks: &[(4.80, "OH"), (3.95, "CH-OH"), (1.60, "CH₂")],
        cnmr_peaks: &[(66.0, "CH-OH"), (45.0, "CH₂")],
    },
    PolymerRef {
        short: "PBA", name: "Poly(butyl acrylate)", repeat_unit: "-[CH₂-CH(COOC₄H₉)]-", m0: 128.17, tg_tm: "Tg ≈ -50 °C (soft segment)",
        ftir_peaks: &[(2960.0, "CH₃ asym"), (2875.0, "CH₂ sym"), (1735.0, "Ester C=O"), (1460.0, "CH₂ scissor"), (1250.0, "C-O-C asym"), (1160.0, "C-O-C sym"), (1060.0, "C-O")],
        hnmr_peaks: &[(4.05, "O-CH₂"), (2.30, "CH backbone"), (1.60, "CH₂ backbone"), (1.40, "CH₂"), (0.95, "CH₃")],
        cnmr_peaks: &[(175.0, "C=O"), (64.5, "O-CH₂"), (41.0, "CH"), (35.0, "CH₂"), (30.5, "CH₂"), (19.2, "CH₂"), (13.7, "CH₃")],
    },
    PolymerRef {
        short: "PLGA", name: "Poly(lactic-co-glycolic acid)", repeat_unit: "-[O-CH(CH₃)-CO]ₘ-[O-CH₂-CO]ₙ-", m0: 130.14, tg_tm: "Tg ≈ 45–55 °C (composition dependent)",
        ftir_peaks: &[(2995.0, "CH₃ asym"), (2945.0, "CH₂ asym"), (1755.0, "Ester C=O"), (1450.0, "CH₃ bend"), (1180.0, "C-O-C asym"), (1090.0, "C-O-C sym"), (820.0, "C-C")],
        hnmr_peaks: &[(5.20, "PLA CH"), (4.80, "PGA CH₂"), (1.58, "PLA CH₃")],
        cnmr_peaks: &[(169.6, "PLA C=O"), (166.5, "PGA C=O"), (69.0, "PLA CH"), (61.0, "PGA CH₂"), (16.7, "PLA CH₃")],
    },
    PolymerRef {
        short: "PHB", name: "Poly(3-hydroxybutyrate)", repeat_unit: "-[O-CH(CH₃)-CH₂-CO]-", m0: 86.09, tg_tm: "Tg ≈ 5 °C, Tm ≈ 175 °C (bacterial polyester)",
        ftir_peaks: &[(2970.0, "CH₃ asym"), (2935.0, "CH₂ asym"), (1725.0, "Ester C=O"), (1455.0, "CH₃ bend"), (1380.0, "CH₃ umbrella"), (1280.0, "C-O-C"), (1225.0, "C-O-C cryst"), (1180.0, "C-O")],
        hnmr_peaks: &[(5.25, "CH"), (2.55, "CH₂"), (1.28, "CH₃")],
        cnmr_peaks: &[(169.5, "C=O"), (67.5, "CH"), (40.8, "CH₂"), (19.8, "CH₃")],
    },
    PolymerRef {
        short: "Chitosan", name: "Chitosan (deacetylated chitin)", repeat_unit: "-[GlcN]-[GlcNAc]-", m0: 161.16, tg_tm: "Tg ≈ 200 °C (biopolymer)",
        ftir_peaks: &[(3450.0, "O-H/N-H broad"), (2870.0, "CH stretch"), (1650.0, "Amide I C=O"), (1590.0, "NH₂ bend"), (1380.0, "CH₃ bend"), (1155.0, "C-O-C asym"), (1030.0, "C-O pyranose"), (895.0, "β-glycosidic")],
        hnmr_peaks: &[(4.60, "H-1 anomeric"), (3.60, "H-2…H-6 ring"), (2.00, "NHAc CH₃")],
        cnmr_peaks: &[(174.0, "NHAc C=O"), (104.0, "C-1"), (82.0, "C-4"), (75.0, "C-3/C-5"), (60.0, "C-6"), (56.0, "C-2"), (23.0, "CH₃")],
    },
    PolymerRef {
        short: "Cellulose", name: "Cellulose", repeat_unit: "-[β-D-glucopyranose-(1→4)]-", m0: 162.14, tg_tm: "Td ≈ 300 °C (semicrystalline biopolymer)",
        ftir_peaks: &[(3400.0, "O-H broad"), (2900.0, "CH₂ asym"), (1640.0, "Adsorbed H₂O"), (1425.0, "CH₂ scissor"), (1160.0, "C-O-C asym"), (1060.0, "C-O ring"), (895.0, "β-glycosidic")],
        hnmr_peaks: &[(4.50, "H-1 anomeric"), (3.70, "H-2…H-6"), (3.20, "H-2")],
        cnmr_peaks: &[(105.0, "C-1"), (88.0, "C-4 cryst"), (75.0, "C-3/C-5"), (72.0, "C-2"), (61.0, "C-6")],
    },
    PolymerRef {
        short: "PBD", name: "Polybutadiene", repeat_unit: "-[CH₂-CH=CH-CH₂]-", m0: 54.09, tg_tm: "Tg ≈ -105 °C (elastomer)",
        ftir_peaks: &[(2915.0, "CH₂ asym"), (2845.0, "CH₂ sym"), (1638.0, "C=C cis"), (1450.0, "CH₂ scissor"), (966.0, "trans C=C"), (911.0, "1,2-vinyl"), (730.0, "cis CH₂ rock")],
        hnmr_peaks: &[(5.40, "1,4-cis CH=CH"), (5.05, "1,2-vinyl CH"), (2.05, "allylic CH₂")],
        cnmr_peaks: &[(130.0, "1,4-cis/trans CH="), (114.0, "1,2-vinyl =CH₂"), (32.9, "trans CH₂"), (27.6, "cis CH₂")],
    },
    PolymerRef {
        short: "PIs", name: "Polyisoprene (Natural Rubber)", repeat_unit: "-[CH₂-C(CH₃)=CH-CH₂]-", m0: 68.12, tg_tm: "Tg ≈ -70 °C (elastomer)",
        ftir_peaks: &[(2960.0, "CH₃ asym"), (2920.0, "CH₂ asym"), (1665.0, "C=C 1,4"), (1450.0, "CH₂ scissor"), (1375.0, "CH₃ umbrella"), (890.0, "3,4-unit"), (840.0, "1,4-unit")],
        hnmr_peaks: &[(5.10, "1,4 CH="), (4.68, "3,4 =CH₂"), (2.05, "allylic CH₂"), (1.68, "CH₃")],
        cnmr_peaks: &[(135.0, "1,4 C(CH₃)="), (125.0, "1,4 CH="), (109.0, "3,4 =CH₂"), (32.0, "CH₂"), (23.4, "CH₃")],
    },
    PolymerRef {
        short: "PANI", name: "Polyaniline (Emeraldine)", repeat_unit: "-[C₆H₄-NH-]-[C₆H₄-N=]-", m0: 93.13, tg_tm: "Tg ≈ 250 °C (conducting polymer)",
        ftir_peaks: &[(1588.0, "Quinoid C=N"), (1495.0, "Benzenoid C=C"), (1300.0, "C-N stretch"), (1140.0, "Doped C-N⁺"), (820.0, "p-Ar OOP"), (505.0, "Ring deformation")],
        hnmr_peaks: &[(7.20, "Ar-H/NH"), (6.90, "Ar-H")],
        cnmr_peaks: &[(156.0, "Quinoid C"), (148.0, "Benzenoid C-N"), (128.0, "Ar CH")],
    },
];

// ============================================================================
// DATABASE LOOKUP HELPERS
// ============================================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BandTier {
    /// Universal functional-group vibrational mode (C=O ester, δ(CH₂/CH₃), C–O–C, etc.)
    Generic,
    /// Material/polymer/mineral-specific marker (Polyimide, CaCO₃, PVDF, Sulfone, POM, Si–CH₃, GO, etc.)
    Specific,
}

#[derive(Debug, Clone, Copy)]
pub struct FtirBandRule {
    pub tier: BandTier,
    pub short_label: &'static str,
    /// Companion wavenumber windows that confirm this band in a full spectrum.
    pub requires: &'static [(f64, f64)],
    /// Contradicting wavenumber windows that suppress this assignment if present.
    pub forbids: &'static [(f64, f64)],
}

#[derive(Debug, Clone, Copy)]
pub struct ScoredFtirBand {
    pub band: &'static FtirBand,
    pub short_label: &'static str,
    pub tier: BandTier,
    pub confirmed: bool,
    pub score: f64,
}

pub fn ftir_band_rule(band: &FtirBand) -> FtirBandRule {
    let g = band.group;
    // --- Material / polymer / mineral specific markers with co-occurrence rules ---
    if g.starts_with("Polyimide Imide I ") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PI Imide I", requires: &[(1710.0, 1732.0), (1362.0, 1390.0)], forbids: &[] };
    }
    if g.starts_with("Polyimide Imide II") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PI Imide II", requires: &[(1768.0, 1790.0), (1362.0, 1390.0)], forbids: &[] };
    }
    if g.starts_with("Polyimide Imide III") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PI Imide III", requires: &[(1768.0, 1790.0), (1710.0, 1732.0)], forbids: &[] };
    }
    if g.starts_with("Polyimide Imide IV") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PI Imide IV", requires: &[(1768.0, 1790.0), (1362.0, 1390.0)], forbids: &[] };
    }
    if g.starts_with("Calcium Carbonate") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "CaCO₃ (CO₃²⁻)", requires: &[(868.0, 885.0), (705.0, 722.0)], forbids: &[] };
    }
    if g.starts_with("CaCO₃ out-of-plane") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "CaCO₃ ν₂", requires: &[(1400.0, 1465.0), (705.0, 722.0)], forbids: &[] };
    }
    if g.starts_with("CaCO₃ in-plane") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "CaCO₃ ν₄", requires: &[(1400.0, 1465.0), (868.0, 885.0)], forbids: &[] };
    }
    if g.starts_with("Sulfone ν_as") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "SO₂ asym", requires: &[(1140.0, 1172.0), (1480.0, 1615.0)], forbids: &[] };
    }
    if g.starts_with("Sulfone ν_s") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "SO₂ sym", requires: &[(1288.0, 1338.0), (1480.0, 1615.0)], forbids: &[] };
    }
    if g.starts_with("Diaryl ether") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "Ar–O–Ar", requires: &[(1480.0, 1615.0)], forbids: &[] };
    }
    if g.starts_with("POM C–O–C") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "POM C–O–C", requires: &[(925.0, 955.0)], forbids: &[(1695.0, 1775.0)] };
    }
    if g.starts_with("POM CH₂ rock") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "POM CH₂ rock", requires: &[(1080.0, 1120.0)], forbids: &[(1695.0, 1775.0)] };
    }
    if g.starts_with("Si–CH₃ symmetric") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "Si–CH₃", requires: &[(1000.0, 1105.0), (780.0, 820.0)], forbids: &[] };
    }
    if g.starts_with("Si–CH₃ rock") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "Si–CH₃ rock", requires: &[(1248.0, 1278.0), (1000.0, 1105.0)], forbids: &[] };
    }
    if g.starts_with("Si–O–Si") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "Si–O–Si", requires: &[(1248.0, 1278.0)], forbids: &[] };
    }
    if g.starts_with("PVDF β-phase (piezo") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PVDF β-phase", requires: &[(1140.0, 1225.0), (1265.0, 1290.0)], forbids: &[] };
    }
    if g.starts_with("PVDF β-phase CF₂") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PVDF β CF₂", requires: &[(830.0, 850.0), (1140.0, 1225.0)], forbids: &[] };
    }
    if g.starts_with("PVDF α-phase") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PVDF α-phase", requires: &[(1140.0, 1225.0), (965.0, 990.0)], forbids: &[] };
    }
    if g.starts_with("C–F / CF₂") {
        return FtirBandRule {
            tier: BandTier::Specific,
            short_label: "C–F / CF₂",
            requires: &[(500.0, 645.0)],
            forbids: &[(1700.0, 1765.0)],
        };
    }
    if g.starts_with("Graphene Oxide (GO) C=O") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "GO C=O", requires: &[(3200.0, 3560.0), (1605.0, 1638.0), (1035.0, 1068.0)], forbids: &[] };
    }
    if g.starts_with("Graphene Oxide") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "GO", requires: &[(3200.0, 3560.0), (1715.0, 1745.0)], forbids: &[] };
    }
    if g.starts_with("PHEMA") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PHEMA", requires: &[(3250.0, 3550.0), (1710.0, 1740.0), (1060.0, 1105.0)], forbids: &[] };
    }
    if g.starts_with("PHB ") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PHB ester", requires: &[(1715.0, 1735.0), (1215.0, 1240.0), (1270.0, 1295.0)], forbids: &[] };
    }
    if g.starts_with("Silica (SiO₂)") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "SiO₂", requires: &[(1055.0, 1135.0), (785.0, 820.0)], forbids: &[] };
    }
    if g.starts_with("Montmorillonite") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "MMT clay", requires: &[(3610.0, 3645.0), (1020.0, 1060.0)], forbids: &[] };
    }
    if g.starts_with("Polyaniline") || g.starts_with("PANI ") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PANI ring", requires: &[(1570.0, 1605.0), (1485.0, 1515.0)], forbids: &[] };
    }
    if g.starts_with("PBT/PET p-Ar") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "p-Ar OOP", requires: &[(1700.0, 1730.0), (1235.0, 1280.0), (718.0, 738.0)], forbids: &[] };
    }
    if g.starts_with("PA6 ") || g.starts_with("PNIPAM ") || g.starts_with("Nylon ") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "Amide (PA/PNIPAM)", requires: &[(1628.0, 1665.0), (1525.0, 1565.0), (3250.0, 3420.0)], forbids: &[] };
    }
    if g.starts_with("PEEK ") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PEEK", requires: &[(1645.0, 1665.0), (1220.0, 1255.0), (1480.0, 1605.0)], forbids: &[] };
    }
    if g.starts_with("PVP ") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PVP lactam", requires: &[(1655.0, 1685.0), (1270.0, 1300.0)], forbids: &[] };
    }
    if g.starts_with("PVA ") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PVA", requires: &[(3200.0, 3520.0), (1080.0, 1120.0)], forbids: &[(1700.0, 1765.0)] };
    }
    if g.starts_with("Chitosan") || g.starts_with("Cellulose") || g.starts_with("β-glycosidic") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "Polysaccharide", requires: &[(1020.0, 1080.0), (3200.0, 3550.0), (885.0, 908.0)], forbids: &[] };
    }
    if g.starts_with("Maleimide") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "Maleimide", requires: &[(1702.0, 1722.0), (688.0, 704.0)], forbids: &[] };
    }
    if g.starts_with("Free Urethane") || g.starts_with("H-bonded Urethane") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "Urethane C=O", requires: &[(3280.0, 3465.0), (1515.0, 1565.0)], forbids: &[] };
    }
    if g.starts_with("Free N–H") || g.starts_with("H-bonded N–H") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "ν(N–H)", requires: &[(1515.0, 1570.0)], forbids: &[] };
    }
    if g.starts_with("O–H (H-bonded") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "O–H (H-bonded)", requires: &[(1000.0, 1160.0), (1620.0, 1655.0)], forbids: &[] };
    }
    if g.starts_with("N–H (amide") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "N–H", requires: &[(1515.0, 1575.0)], forbids: &[] };
    }
    if g.starts_with("H₂O (adsorbed") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "δ(H₂O)", requires: &[(3200.0, 3550.0)], forbids: &[(800.0, 995.0)] };
    }
    if g.starts_with("≡C–H") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "≡C–H", requires: &[(2100.0, 2160.0)], forbids: &[] };
    }
    if g.starts_with("PE crystalline") || g.starts_with("LDPE ") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "PE CH₂", requires: &[(2905.0, 2935.0), (2840.0, 2865.0), (715.0, 735.0)], forbids: &[(1680.0, 1775.0)] };
    }
    if g.starts_with("ABS") {
        return FtirBandRule { tier: BandTier::Specific, short_label: "ABS/SAN", requires: &[(2230.0, 2255.0), (690.0, 710.0)], forbids: &[] };
    }
    if g.starts_with("PBT ") || g.starts_with("PAA ") || g.starts_with("PBA ") || g.starts_with("PLGA")
        || g.starts_with("PBD ") || g.starts_with("Polyisoprene") || g.starts_with("Cured epoxy")
        || g.starts_with("Epoxy oxirane") || g.starts_with("TiO₂") || g.starts_with("ZnO")
        || g.starts_with("Thiocarbonylthio") || g.starts_with("Aliphatic C–Br") || g.starts_with("Azide ")
    {
        let short = g.split('(').next().unwrap_or(g).trim();
        return FtirBandRule { tier: BandTier::Specific, short_label: short, requires: &[], forbids: &[] };
    }

    // --- Universal functional-group bands (Generic tier) ---
    if g.starts_with("–(CH₂)_n–") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "ρ(CH₂) rock", requires: &[(2840.0, 2945.0)], forbids: &[] };
    }
    if g.starts_with("CH₂ scissor") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "δ(CH₂/CH₃)", requires: &[(2835.0, 3005.0)], forbids: &[] };
    }
    if g.starts_with("CH₃ symmetric") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "δ_s(CH₃)", requires: &[(2940.0, 3005.0)], forbids: &[] };
    }
    if g.starts_with("C=O (ester / methacrylate") || g.starts_with("C=O (saturated aliphatic ester)") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "C=O (ester)", requires: &[(1120.0, 1320.0)], forbids: &[] };
    }
    if g.starts_with("C=O (Amide I)") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "Amide I C=O", requires: &[(1515.0, 1575.0)], forbids: &[] };
    }
    if g.starts_with("N–H bend + C–N (Amide II)") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "Amide II", requires: &[(1625.0, 1685.0)], forbids: &[] };
    }
    if g.starts_with("C=C (vinyl / acrylate") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "C=C (vinyl)", requires: &[(800.0, 995.0)], forbids: &[] };
    }
    if g.starts_with("C–O (α,β-unsaturated ester") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "C–O (conj. ester)", requires: &[(1700.0, 1745.0), (1615.0, 1660.0)], forbids: &[] };
    }
    if g.starts_with("C–O (ester / aryl ether") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "C–O (ester)", requires: &[(1700.0, 1785.0)], forbids: &[] };
    }
    if g.starts_with("C–O–C (ester / aliphatic ether)") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "C–O–C (ester)", requires: &[(1700.0, 1785.0)], forbids: &[] };
    }
    if g.starts_with("C–O–C (aliphatic ether") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "C–O–C (ether)", requires: &[], forbids: &[] };
    }
    if g.starts_with("Epoxide oxirane ring / =CH₂ wag") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "=CH₂ wag / epoxide", requires: &[(1620.0, 1660.0)], forbids: &[] };
    }
    if g.starts_with("=CH₂ twist") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "=CH₂ twist", requires: &[(1620.0, 1660.0)], forbids: &[] };
    }
    if g.starts_with("C–H (CH₃ asym") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "ν_as(CH₃)", requires: &[(1360.0, 1480.0)], forbids: &[] };
    }
    if g.starts_with("C–H (CH₂ asym") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "ν_as(CH₂)", requires: &[(1430.0, 1480.0)], forbids: &[] };
    }
    if g.starts_with("C–H (CH₃ sym") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "ν_s(CH₃)", requires: &[(2940.0, 3000.0)], forbids: &[] };
    }
    if g.starts_with("C–H (CH₂ sym") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "ν_s(CH₂)", requires: &[(2900.0, 2948.0)], forbids: &[] };
    }
    if g.starts_with("Mono-substituted phenyl") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "Ar ring bend", requires: &[(735.0, 775.0), (1480.0, 1615.0)], forbids: &[] };
    }
    if g.starts_with("Mono-/ortho-Aromatic") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "γ(Ar–H) OOP", requires: &[(1480.0, 1615.0)], forbids: &[] };
    }
    if g.starts_with("CO₂ (atmospheric") {
        return FtirBandRule { tier: BandTier::Generic, short_label: "CO₂ (atm)", requires: &[], forbids: &[] };
    }

    let short = g.split('(').next().unwrap_or(g).trim();
    let short_label = if short.is_empty() || short == "–" { g } else { short };
    FtirBandRule { tier: BandTier::Generic, short_label, requires: &[], forbids: &[] }
}

pub fn match_ftir_with_context(wn: f64, spectrum_peaks: &[f64]) -> Vec<ScoredFtirBand> {
    let has_peak = |lo: f64, hi: f64| spectrum_peaks.iter().any(|&p| p >= lo && p <= hi);
    let mut exact: Vec<ScoredFtirBand> = FTIR_DB
        .iter()
        .filter(|b| wn >= b.min && wn <= b.max)
        .map(|b| {
            let rule = ftir_band_rule(b);
            let width = (b.max - b.min).max(10.0);
            let center_dist = (wn - 0.5 * (b.min + b.max)).abs() / (0.5 * width);
            let base = center_dist + 0.0012 * width;

            let (confirmed, mut score) = if !spectrum_peaks.is_empty() {
                match rule.tier {
                    BandTier::Specific => {
                        let forb_hit = rule.forbids.iter().any(|&(lo, hi)| has_peak(lo, hi));
                        let all_req = !rule.requires.is_empty()
                            && rule.requires.iter().all(|&(lo, hi)| has_peak(lo, hi));
                        if !forb_hit && all_req {
                            (true, base - 0.45)
                        } else {
                            (false, base + 2.50)
                        }
                    }
                    BandTier::Generic => {
                        let n_supp = rule
                            .requires
                            .iter()
                            .filter(|&&(lo, hi)| has_peak(lo, hi))
                            .count();
                        (true, base - 0.30 * (n_supp as f64))
                    }
                }
            } else {
                // Context-free lookup: slightly prefer universal functional groups
                match rule.tier {
                    BandTier::Generic => (true, base),
                    BandTier::Specific => (false, base + 0.35),
                }
            };

            // Without Amide II evidence, broad H-bonded O–H is the natural
            // assignment; with Amide II present the N–H stretch is supported.
            if rule.short_label == "O–H (H-bonded)" && !spectrum_peaks.is_empty() {
                if has_peak(1515.0, 1575.0) {
                    score += 0.50;
                } else {
                    score -= 0.15;
                }
            }

            ScoredFtirBand {
                band: b,
                short_label: rule.short_label,
                tier: rule.tier,
                confirmed,
                score,
            }
        })
        .collect();

    exact.sort_by(|a, b| a.score.total_cmp(&b.score));
    if !exact.is_empty() {
        return exact;
    }

    // Nearest within ±45 cm⁻¹
    let mut near: Vec<(ScoredFtirBand, f64)> = FTIR_DB
        .iter()
        .filter_map(|b| {
            let d = if wn < b.min { b.min - wn } else { wn - b.max };
            if d <= 45.0 {
                let rule = ftir_band_rule(b);
                let penalty = if rule.tier == BandTier::Specific { 15.0 } else { 0.0 };
                Some((
                    ScoredFtirBand {
                        band: b,
                        short_label: rule.short_label,
                        tier: rule.tier,
                        confirmed: rule.tier == BandTier::Generic,
                        score: d + penalty,
                    },
                    d + penalty,
                ))
            } else {
                None
            }
        })
        .collect();
    near.sort_by(|a, b| a.1.total_cmp(&b.1));
    near.into_iter().take(3).map(|(s, _)| s).collect()
}

pub fn match_ftir(wn: f64) -> Vec<&'static FtirBand> {
    match_ftir_with_context(wn, &[])
        .into_iter()
        .map(|s| s.band)
        .collect()
}

pub fn match_hnmr(ppm: f64) -> Vec<&'static NmrBand> {
    let mut exact: Vec<&'static NmrBand> = HNMR_DB
        .iter()
        .filter(|b| ppm >= b.min && ppm <= b.max)
        .collect();
    exact.sort_by(|a, b| {
        let da = (ppm - 0.5 * (a.min + a.max)).abs();
        let db = (ppm - 0.5 * (b.min + b.max)).abs();
        da.total_cmp(&db)
    });
    exact
}

pub fn match_cnmr(ppm: f64) -> Vec<&'static NmrBand> {
    let mut exact: Vec<&'static NmrBand> = CNMR_DB
        .iter()
        .filter(|b| ppm >= b.min && ppm <= b.max)
        .collect();
    exact.sort_by(|a, b| {
        let da = (ppm - 0.5 * (a.min + a.max)).abs();
        let db = (ppm - 0.5 * (b.min + b.max)).abs();
        da.total_cmp(&db)
    });
    exact
}

pub fn match_uvvis(nm: f64) -> Vec<&'static UvVisBand> {
    let mut exact: Vec<&'static UvVisBand> = UVVIS_DB
        .iter()
        .filter(|b| nm >= b.min && nm <= b.max)
        .collect();
    exact.sort_by(|a, b| {
        let da = (nm - 0.5 * (a.min + a.max)).abs();
        let db = (nm - 0.5 * (b.min + b.max)).abs();
        da.total_cmp(&db)
    });
    exact
}

pub fn photon_energy_ev(nm: f64) -> f64 {
    if nm > 0.0 && nm.is_finite() {
        1239.84193 / nm
    } else {
        0.0
    }
}

/// Returns `(absorbed_region, complementary_observed_solution_color)`.
pub fn uvvis_color_info(nm: f64) -> (&'static str, &'static str) {
    if nm < 280.0 {
        ("UV-C (<280 nm)", "Colorless (transparent)")
    } else if nm < 315.0 {
        ("UV-B (280–315 nm)", "Colorless (transparent)")
    } else if nm < 380.0 {
        ("UV-A (315–380 nm)", "Colorless / faint pale yellow")
    } else if nm < 430.0 {
        ("Violet (380–430 nm)", "Yellow-Green")
    } else if nm < 480.0 {
        ("Blue (430–480 nm)", "Yellow-Orange")
    } else if nm < 500.0 {
        ("Blue-Green (480–500 nm)", "Orange-Red")
    } else if nm < 560.0 {
        ("Green (500–560 nm)", "Purple / Magenta")
    } else if nm < 580.0 {
        ("Yellow-Green (560–580 nm)", "Violet")
    } else if nm < 600.0 {
        ("Yellow (580–600 nm)", "Blue-Violet")
    } else if nm < 650.0 {
        ("Orange (600–650 nm)", "Cyan / Blue")
    } else if nm < 750.0 {
        ("Red (650–750 nm)", "Blue-Green / Cyan")
    } else {
        ("Near-Infrared (>750 nm)", "Colorless (NIR-absorbing)")
    }
}

pub fn match_solvent_impurities(
    solvent: SolventId,
    is_carbon: bool,
    ppm: f64,
    tol: f64,
) -> Vec<(&'static SolventImpurity, f64, f64)> {
    let idx = solvent.index();
    let mut out: Vec<(&'static SolventImpurity, f64, f64)> = SOLVENT_IMPURITIES
        .iter()
        .filter_map(|entry| {
            let shift = if is_carbon { entry.c_shifts[idx] } else { entry.h_shifts[idx] }?;
            let diff = (ppm - shift).abs();
            if diff <= tol {
                Some((entry, shift, diff))
            } else {
                None
            }
        })
        .collect();
    out.sort_by(|a, b| a.2.total_cmp(&b.2));
    out
}

pub fn match_polymers(kind: SpectrumKind, x: f64) -> Vec<(&'static PolymerRef, f64, &'static str)> {
    if kind == SpectrumKind::UvVis {
        return Vec::new();
    }
    let tol = match kind {
        SpectrumKind::Ftir => 18.0,
        SpectrumKind::Hnmr => 0.14,
        SpectrumKind::Cnmr => 2.5,
        SpectrumKind::UvVis => return Vec::new(),
    };
    let mut out = Vec::new();
    for p in POLYMERS_DB {
        let peaks = match kind {
            SpectrumKind::Ftir => p.ftir_peaks,
            SpectrumKind::Hnmr => p.hnmr_peaks,
            SpectrumKind::Cnmr => p.cnmr_peaks,
            SpectrumKind::UvVis => &[],
        };
        let mut best: Option<(f64, &'static str)> = None;
        for &(px, desc) in peaks {
            let d = (x - px).abs();
            if d <= tol && best.is_none_or(|(bd, _)| d < bd) {
                best = Some((d, desc));
            }
        }
        if let Some((d, desc)) = best {
            out.push((p, d, desc));
        }
    }
    out.sort_by(|a, b| a.1.total_cmp(&b.1));
    out
}

pub fn spectral_zone_name(kind: SpectrumKind, x: f64) -> &'static str {
    match kind {
        SpectrumKind::Ftir => {
            if x >= 3200.0 { "X–H Stretch (O–H / N–H)" }
            else if x >= 2750.0 { "C–H Stretch Region" }
            else if x >= 2000.0 { "Triple / Cumulated (C≡N / NCO)" }
            else if x >= 1630.0 { "Carbonyl (C=O) Region" }
            else if x >= 1450.0 { "Double Bond / Aromatic / Bend" }
            else { "Fingerprint Region (<1450)" }
        }
        SpectrumKind::Hnmr => {
            if x >= 9.0 { "Aldehyde / Carboxylic Acid" }
            else if x >= 6.5 { "Aromatic Region" }
            else if x >= 4.5 { "Vinyl / Alkene Region" }
            else if x >= 3.0 { "Heteroatom (O–CH / N–CH)" }
            else if x >= 1.8 { "α-Carbonyl / Allylic / Backbone" }
            else if x >= 0.4 { "Aliphatic / Tacticity CH₃" }
            else { "Silane / TMS Region" }
        }
        SpectrumKind::Cnmr => {
            if x >= 160.0 { "Carbonyl (C=O) Region" }
            else if x >= 110.0 { "Aromatic / sp² / Nitrile" }
            else if x >= 50.0 { "C–O / C–Cl / OCH₃ Region" }
            else { "Aliphatic / Backbone / Tacticity" }
        }
        SpectrumKind::UvVis => {
            if x < 200.0 { "Far-UV / σ→σ* & Isolated C=C (<200 nm)" }
            else if x < 280.0 { "UV-C / π→π* Chromophores (200–280 nm)" }
            else if x < 315.0 { "UV-B / Conjugated Aromatics (280–315 nm)" }
            else if x < 380.0 { "UV-A / n→π* & Band Edge (315–380 nm)" }
            else if x < 500.0 { "Visible Short (Violet–Blue, 380–500 nm)" }
            else if x < 620.0 { "Visible Mid (Green–Orange, 500–620 nm)" }
            else if x < 750.0 { "Visible Long (Red, 620–750 nm)" }
            else { "Near-Infrared / Polaron NIR (>750 nm)" }
        }
    }
}

fn silicone_grease_shift(solvent: SolventId, is_carbon: bool) -> Option<f64> {
    SOLVENT_IMPURITIES
        .iter()
        .find(|e| e.name.contains("Silicone"))
        .and_then(|e| {
            if is_carbon {
                e.c_shifts[solvent.index()]
            } else {
                e.h_shifts[solvent.index()]
            }
        })
}

pub fn best_peak_label_in_context(
    kind: SpectrumKind,
    x: f64,
    spectrum_peaks: &[f64],
    solvent: SolventId,
) -> String {
    match kind {
        SpectrumKind::Ftir => {
            if let Some(hit) = match_ftir_with_context(x, spectrum_peaks).first() {
                let label = if hit.short_label == "=CH₂ wag / epoxide" {
                    if spectrum_peaks.iter().any(|&p| (1618.0..=1662.0).contains(&p)) {
                        "=CH₂ wag"
                    } else {
                        "Epoxide / =CH₂"
                    }
                } else {
                    hit.short_label
                };
                format!("{:.0} {}", x, label)
            } else {
                format!("{:.0}", x)
            }
        }
        SpectrumKind::Hnmr | SpectrumKind::Cnmr | SpectrumKind::UvVis => {
            best_short_peak_label_with_solvent(kind, x, solvent)
        }
    }
}

/// Short peak label that recognizes the active NMR solvent and the most
/// common trace impurities before falling back to functional-group assignment.
/// Impurity labels are only used when the impurity is closer than any polymer
/// or functional-group match (avoids false positives on PMMA O–CH₃ at 3.60 ppm
/// and PMMA C=O at ~177 ppm).
pub fn best_short_peak_label_with_solvent(
    kind: SpectrumKind,
    x: f64,
    solvent: SolventId,
) -> String {
    match kind {
        SpectrumKind::Ftir => best_peak_label_in_context(SpectrumKind::Ftir, x, &[], solvent),
        SpectrumKind::Hnmr | SpectrumKind::Cnmr => {
            let is_carbon = kind == SpectrumKind::Cnmr;
            let (decimals, tol) = if is_carbon { (1usize, 0.55) } else { (2usize, 0.035) };

            // Compute closest polymer peak distance (from POLYMERS_DB)
            let poly_dist = POLYMERS_DB
                .iter()
                .flat_map(|p| {
                    let peaks = if is_carbon { p.cnmr_peaks } else { p.hnmr_peaks };
                    peaks.iter().map(|(px, _)| (x - *px).abs())
                })
                .fold(f64::INFINITY, f64::min);

            // Compute closest functional-group match distance (from NMR_DB)
            let fg_dist = if is_carbon {
                match_cnmr(x)
            } else {
                match_hnmr(x)
            }
            .first()
            .map(|b| (x - 0.5 * (b.min + b.max)).abs())
            .unwrap_or(f64::INFINITY);

            let min_signal_dist = poly_dist.min(fg_dist);

            // 1) Residual solvent signal (e.g. CDCl₃ at 7.26 ppm)
            let residual = if is_carbon {
                solvent.residual_c()
            } else {
                solvent.residual_h()
            };
            if let Some(&s) = residual.iter().min_by(|a, b| {
                (x - **a).abs().total_cmp(&(x - **b).abs())
            }) {
                let imp_dist = (x - s).abs();
                if imp_dist <= tol && imp_dist <= min_signal_dist {
                    return format!("{:.*} {} (solvent)", decimals, x, solvent.label());
                }
            }

            // 2) Water / HDO impurity (¹H only)
            if !is_carbon {
                let water_shift = solvent.water_h();
                let imp_dist = (x - water_shift).abs();
                if imp_dist <= tol && imp_dist <= min_signal_dist {
                    return format!("{:.*} H₂O (impurity)", decimals, x);
                }
            }

            // 3) Silicone grease (PDMS) impurity
            if let Some(shift) = silicone_grease_shift(solvent, is_carbon) {
                let imp_dist = (x - shift).abs();
                if imp_dist <= tol && imp_dist <= min_signal_dist {
                    return format!("{:.*} Silicone grease (impurity)", decimals, x);
                }
            }

            // 4) Normal functional-group assignment
            let fallback = if is_carbon {
                match_cnmr(x)
            } else {
                match_hnmr(x)
            };
            if let Some(b) = fallback.first() {
                let short = b.group.split('(').next().unwrap_or(b.group).trim();
                format!("{:.*} {}", decimals, x, short)
            } else {
                format!("{:.*}", decimals, x)
            }
        }
        SpectrumKind::UvVis => {
            if let Some(b) = match_uvvis(x).first() {
                let short = b.group.split('(').next().unwrap_or(b.group).trim();
                format!("{:.0} nm {}", x, short)
            } else {
                format!("{:.0} nm", x)
            }
        }
    }
}

// ============================================================================
// SYNTHETIC SCIENTIFIC DEMO SPECTRA GENERATOR
// ============================================================================

fn lorentz(x: f64, x0: f64, h: f64, gamma: f64) -> f64 {
    let d = (x - x0) / gamma;
    h / (1.0 + d * d)
}

fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

fn pseudo_noise(i: usize, amp: f64) -> f64 {
    let h1 = splitmix64((i as u64).wrapping_mul(0xD1B5_4A32_D192_ED03).wrapping_add(1));
    let h2 = splitmix64(h1.wrapping_add(0x94D0_49BB_1331_11EB));
    let u1 = ((h1 >> 11) as f64) * (1.0 / ((1u64 << 53) as f64)) - 0.5;
    let u2 = ((h2 >> 11) as f64) * (1.0 / ((1u64 << 53) as f64)) - 0.5;
    (u1 + u2) * amp
}

fn synth_ftir(
    name: &str,
    color_idx: usize,
    peaks: &[(f64, f64, f64)],
    slope: f64,
    offset: f64,
    noise: f64,
) -> Spectrum {
    let n = 1801;
    let mut x = Vec::with_capacity(n);
    let mut y = Vec::with_capacity(n);
    for i in 0..n {
        let wn = 400.0 + (i as f64) * 2.0; // 400..4000 step 2 cm⁻¹
        let mut abs = offset + slope * (wn - 400.0) / 3600.0 + pseudo_noise(i + color_idx * 997, noise);
        for &(x0, h, g) in peaks {
            abs += lorentz(wn, x0, h, g);
        }
        let abs_clamped = abs.max(0.005);
        let pct_t = (100.0 * 10.0_f64.powf(-abs_clamped)).clamp(0.5, 100.0);
        x.push(wn);
        y.push(pct_t);
    }
    Spectrum {
        name: name.to_string(),
        source_path: None,
        kind: SpectrumKind::Ftir,
        raw_y_unit: YUnit::Transmittance,
        x,
        y,
        visible: true,
        color_idx,
        metadata: vec![("Origin".into(), "SPECTRA Built-in Reference".into())],
    }
}

#[allow(clippy::too_many_arguments)]
fn synth_nmr(
    name: &str,
    kind: SpectrumKind,
    color_idx: usize,
    x_lo: f64,
    x_hi: f64,
    n: usize,
    peaks: &[(f64, f64, f64)],
    noise: f64,
) -> Spectrum {
    let mut x = Vec::with_capacity(n);
    let mut y = Vec::with_capacity(n);
    let step = (x_hi - x_lo) / (n - 1) as f64;
    for i in 0..n {
        let ppm = x_lo + i as f64 * step;
        let mut val = pseudo_noise(i + color_idx * 521, noise).abs() * 0.5;
        for &(x0, h, g) in peaks {
            val += lorentz(ppm, x0, h, g);
        }
        x.push(ppm);
        y.push(val);
    }
    Spectrum {
        name: name.to_string(),
        source_path: None,
        kind,
        raw_y_unit: YUnit::Intensity,
        x,
        y,
        visible: true,
        color_idx,
        metadata: vec![
            ("Origin".into(), "SPECTRA Built-in Reference".into()),
            ("Solvent".into(), "CDCl₃ (400 MHz)".into()),
        ],
    }
}

fn gaussian(x: f64, x0: f64, h: f64, sigma: f64) -> f64 {
    let d = (x - x0) / sigma;
    h * (-0.5 * d * d).exp()
}

fn synth_uvvis(
    name: &str,
    color_idx: usize,
    peaks: &[(f64, f64, f64)],
    baseline: f64,
    noise: f64,
) -> Spectrum {
    let n = 711;
    let mut x = Vec::with_capacity(n);
    let mut y = Vec::with_capacity(n);
    for i in 0..n {
        let nm = 190.0 + i as f64;
        let mut val = baseline + pseudo_noise(i + color_idx * 419, noise).abs() * 0.4;
        for &(x0, h, s) in peaks {
            val += gaussian(nm, x0, h, s);
        }
        x.push(nm);
        y.push(val.max(0.0));
    }
    Spectrum {
        name: name.to_string(),
        source_path: None,
        kind: SpectrumKind::UvVis,
        raw_y_unit: YUnit::Absorbance,
        x,
        y,
        visible: true,
        color_idx,
        metadata: vec![("Origin".into(), "SPECTRA Built-in Reference".into())],
    }
}

fn synth_semiconductor_uvvis(
    name: &str,
    color_idx: usize,
    eg_ev: f64,
    direct: bool,
    amp: f64,
    exciton: Option<(f64, f64, f64)>,
    noise: f64,
) -> Spectrum {
    let n = 711;
    let mut x = Vec::with_capacity(n);
    let mut y = Vec::with_capacity(n);
    for i in 0..n {
        let nm = 190.0 + i as f64;
        let e = 1239.84193 / nm;
        let alpha = if e > eg_ev {
            let de = e - eg_ev;
            // Tauc-consistent: direct α·hν ∝ √(hν−Eg); indirect α·hν ∝ (hν−Eg)²
            if direct {
                amp * de.sqrt() / e
            } else {
                amp * de * de / e
            }
        } else {
            // Urbach (excitonic) tail below the band edge
            0.012 * ((e - eg_ev) / 0.06).exp()
        };
        let exc = exciton
            .map(|(x0, h, s)| gaussian(nm, x0, h, s))
            .unwrap_or(0.0);
        let val = 0.02 + alpha + exc + pseudo_noise(i + color_idx * 313, noise).abs() * 0.3;
        x.push(nm);
        y.push(val.max(0.0));
    }
    Spectrum {
        name: name.to_string(),
        source_path: None,
        kind: SpectrumKind::UvVis,
        raw_y_unit: YUnit::Absorbance,
        x,
        y,
        visible: true,
        color_idx,
        metadata: vec![("Origin".into(), "SPECTRA Built-in Reference".into())],
    }
}

pub fn tio2_zno_uvvis_spectrum() -> Spectrum {
    // Anatase TiO₂ direct edge ~385 nm (3.22 eV) + ZnO exciton shoulder ~366 nm
    let mut spec = synth_semiconductor_uvvis(
        "TiO₂/ZnO Nanocomposite (UV-Vis)",
        0,
        3.22,
        true,
        5.4,
        Some((366.0, 0.12, 10.0)),
        0.004,
    );
    spec.metadata.push((
        "Assignment".into(),
        "Anatase TiO₂ (~385 nm, 3.22 eV) + ZnO (~375 nm, 3.31 eV)".into(),
    ));
    spec
}

pub fn mma_monomer_spectrum() -> Spectrum {
    synth_ftir(
        "MMA Monomer (t=0)",
        0,
        &[
            (3105.0, 0.08, 12.0),
            (2994.0, 0.28, 14.0),
            (2954.0, 0.34, 16.0),
            (1722.0, 0.92, 14.0), // α,β-unsat C=O
            (1638.0, 0.58, 10.0), // Strong C=C vinyl!
            (1450.0, 0.36, 15.0),
            (1438.0, 0.28, 12.0),
            (1324.0, 0.42, 14.0),
            (1300.0, 0.48, 14.0),
            (1200.0, 0.55, 16.0),
            (1162.0, 0.85, 22.0),
            (1015.0, 0.24, 14.0),
            (940.0, 0.44, 12.0),  // =CH₂ wag
            (814.0, 0.36, 10.0),  // =CH₂ twist
            (652.0, 0.18, 12.0),
        ],
        0.03,
        0.02,
        0.004,
    )
}

pub fn pmma_polymer_spectrum() -> Spectrum {
    synth_ftir(
        "PMMA Polymer (92% conv)",
        1,
        &[
            (2995.0, 0.32, 18.0),
            (2950.0, 0.40, 20.0),
            (2844.0, 0.12, 16.0),
            (1728.0, 0.95, 16.0), // Saturated ester C=O
            (1638.0, 0.046, 10.0), // ~92% depleted C=C!
            (1484.0, 0.22, 15.0),
            (1448.0, 0.44, 18.0),
            (1388.0, 0.25, 14.0),
            (1270.0, 0.38, 16.0),
            (1242.0, 0.54, 18.0),
            (1192.0, 0.58, 16.0),
            (1148.0, 0.88, 26.0), // Strong C-O-C
            (988.0, 0.34, 16.0),
            (965.0, 0.22, 14.0),
            (842.0, 0.26, 14.0),
            (750.0, 0.32, 14.0),
        ],
        0.06, // Slight baseline tilt to demonstrate Rubberband baseline!
        0.035,
        0.006,
    )
}

pub fn polystyrene_ftir_spectrum() -> Spectrum {
    synth_ftir(
        "Polystyrene (PS)",
        0,
        &[
            (3082.0, 0.20, 12.0),
            (3060.0, 0.28, 12.0),
            (3026.0, 0.44, 14.0),
            (2924.0, 0.62, 18.0),
            (2850.0, 0.36, 16.0),
            (1942.0, 0.07, 16.0),
            (1870.0, 0.06, 16.0),
            (1802.0, 0.06, 16.0),
            (1601.0, 0.45, 10.0),
            (1583.0, 0.16, 10.0),
            (1493.0, 0.72, 12.0),
            (1452.0, 0.78, 14.0),
            (1372.0, 0.16, 14.0),
            (1182.0, 0.14, 12.0),
            (1154.0, 0.14, 12.0),
            (1069.0, 0.22, 12.0),
            (1028.0, 0.34, 12.0),
            (906.0, 0.20, 12.0),
            (842.0, 0.12, 12.0),
            (756.0, 0.84, 14.0),
            (698.0, 1.05, 14.0),
            (540.0, 0.40, 16.0),
        ],
        0.02,
        0.02,
        0.003,
    )
}

pub fn pmma_hnmr_spectrum() -> Spectrum {
    synth_nmr(
        "PMMA (CDCl₃, 400 MHz)",
        SpectrumKind::Hnmr,
        0,
        0.0,
        10.0,
        4000,
        &[
            (7.26, 0.36, 0.012), // Residual CDCl₃
            (6.10, 0.08, 0.014), // Trace MMA monomer vinyl H
            (5.56, 0.08, 0.014), // Trace MMA monomer vinyl H
            (3.60, 1.00, 0.026), // PMMA -OCH₃ (3H)
            (1.91, 0.30, 0.045), // PMMA backbone -CH₂-
            (1.81, 0.36, 0.048), // PMMA backbone -CH₂-
            (1.56, 0.22, 0.018), // H₂O impurity in CDCl₃!
            (1.21, 0.14, 0.022), // PMMA α-CH₃ mm (isotactic)
            (1.02, 0.38, 0.025), // PMMA α-CH₃ mr (heterotactic)
            (0.85, 0.62, 0.025), // PMMA α-CH₃ rr (syndiotactic)
            (0.07, 0.10, 0.010), // Silicone grease trace
            (0.00, 0.15, 0.009), // TMS reference
        ],
        0.006,
    )
}

pub fn generate_demo_preset(preset: DemoPreset) -> Vec<Spectrum> {
    match preset {
        DemoPreset::MmaVsPmmaFtir => vec![mma_monomer_spectrum(), pmma_polymer_spectrum()],
        DemoPreset::PolystyreneFtir => vec![polystyrene_ftir_spectrum()],
        DemoPreset::PetVsNylonFtir => vec![
            synth_ftir(
                "PET (Polyester)",
                0,
                &[
                    (2968.0, 0.22, 16.0),
                    (1715.0, 0.96, 16.0),
                    (1578.0, 0.14, 12.0),
                    (1505.0, 0.22, 12.0),
                    (1408.0, 0.34, 12.0),
                    (1340.0, 0.40, 14.0),
                    (1245.0, 0.92, 24.0),
                    (1096.0, 0.86, 22.0),
                    (1018.0, 0.45, 14.0),
                    (872.0, 0.42, 12.0),
                    (725.0, 0.82, 14.0),
                ],
                0.02,
                0.02,
                0.003,
            ),
            synth_ftir(
                "Nylon-6,6 (Polyamide)",
                1,
                &[
                    (3300.0, 0.68, 42.0),
                    (3075.0, 0.18, 22.0),
                    (2932.0, 0.64, 20.0),
                    (2858.0, 0.48, 18.0),
                    (1640.0, 0.95, 18.0),
                    (1540.0, 0.84, 20.0),
                    (1464.0, 0.36, 16.0),
                    (1372.0, 0.24, 14.0),
                    (1275.0, 0.38, 16.0),
                    (1198.0, 0.26, 14.0),
                    (936.0, 0.22, 14.0),
                    (690.0, 0.45, 28.0),
                ],
                0.025,
                0.025,
                0.004,
            ),
        ],
        DemoPreset::PmmaHnmr => vec![pmma_hnmr_spectrum()],
        DemoPreset::PolystyreneHnmr => vec![synth_nmr(
            "Polystyrene (CDCl₃, 400 MHz)",
            SpectrumKind::Hnmr,
            2,
            0.0,
            10.0,
            4000,
            &[
                (7.26, 0.34, 0.012), // CDCl₃
                (7.08, 0.85, 0.095), // m/p Ar-H (3H)
                (6.55, 0.56, 0.095), // ortho Ar-H (2H)
                (1.85, 0.38, 0.075), // backbone -CH- (1H)
                (1.56, 0.18, 0.018), // H₂O
                (1.43, 0.72, 0.085), // backbone -CH₂- (2H)
                (0.00, 0.14, 0.010), // TMS
            ],
            0.006,
        )],
        DemoPreset::PmmaCnmr => vec![synth_nmr(
            "PMMA (¹³C NMR, CDCl₃)",
            SpectrumKind::Cnmr,
            2,
            0.0,
            210.0,
            4200,
            &[
                (177.8, 0.45, 0.40), // C=O rr
                (177.0, 0.28, 0.35), // C=O mr
                (77.48, 0.62, 0.16), // CDCl₃ triplet
                (77.16, 0.68, 0.16),
                (76.84, 0.62, 0.16),
                (54.4, 0.64, 0.45),  // CH₂ backbone
                (51.8, 0.95, 0.25),  // OCH₃
                (44.9, 0.44, 0.30),  // Cq mr
                (44.5, 0.62, 0.30),  // Cq rr
                (21.1, 0.18, 0.30),  // α-CH₃ mm
                (18.7, 0.46, 0.32),  // α-CH₃ mr
                (16.5, 0.72, 0.32),  // α-CH₃ rr
            ],
            0.006,
        )],
        DemoPreset::PaniUvVis => vec![
            synth_uvvis(
                "Polyaniline Emeraldine Base (EB)",
                0,
                &[(325.0, 0.90, 28.0), (630.0, 0.74, 48.0)],
                0.02,
                0.004,
            ),
            synth_uvvis(
                "Polyaniline Emeraldine Salt (ES, H⁺-doped)",
                1,
                &[
                    (325.0, 0.68, 26.0),
                    (430.0, 0.64, 32.0),
                    (815.0, 0.92, 85.0),
                ],
                0.03,
                0.004,
            ),
        ],
        DemoPreset::Tio2ZnoNanocompositeUvVis => vec![
            tio2_zno_uvvis_spectrum(),
            synth_semiconductor_uvvis(
                "ZnO Nanoparticles (Direct Edge ~375 nm)",
                1,
                3.31,
                true,
                5.0,
                Some((366.0, 0.22, 9.0)),
                0.004,
            ),
        ],
        DemoPreset::AgAuNanoparticlesUvVis => vec![
            synth_uvvis(
                "Ag Nanoparticles (LSPR 412 nm)",
                0,
                &[(240.0, 0.32, 35.0), (412.0, 0.96, 22.0)],
                0.02,
                0.004,
            ),
            synth_uvvis(
                "Au Nanoparticles (LSPR 524 nm)",
                1,
                &[(280.0, 0.38, 45.0), (524.0, 0.88, 24.0)],
                0.02,
                0.004,
            ),
        ],
    }
}

// ============================================================================
// POLYMER REFERENCE SYNTHESIS (used by the "/" search palette)
// ============================================================================

/// Look up a polymer by its short code (e.g. "PVC", "PTFE", "Epoxy").
pub fn find_polymer_by_short(short: &str) -> Option<&'static PolymerRef> {
    POLYMERS_DB
        .iter()
        .find(|p| p.short.eq_ignore_ascii_case(short))
}

/// Build a synthetic reference spectrum from the polymer fingerprint database.
/// Falls back to FTIR when the requested NMR kind has no tabulated peaks.
pub fn synthesize_polymer_reference(
    poly: &PolymerRef,
    requested: SpectrumKind,
) -> Spectrum {
    let kind = match requested {
        SpectrumKind::Hnmr if !poly.hnmr_peaks.is_empty() => SpectrumKind::Hnmr,
        SpectrumKind::Cnmr if !poly.cnmr_peaks.is_empty() => SpectrumKind::Cnmr,
        _ => SpectrumKind::Ftir,
    };

    let title = format!("{} — {} [Ref]", poly.short, poly.name);
    let mut spec = match kind {
        SpectrumKind::Ftir => {
            let peaks: Vec<(f64, f64, f64)> = poly
                .ftir_peaks
                .iter()
                .map(|&(wn, desc)| {
                    let d = desc.to_ascii_lowercase();
                    let h = if d.contains("very strong") {
                        0.95
                    } else if d.contains("weak") {
                        0.22
                    } else if d.contains("medium") {
                        0.50
                    } else if d.contains("strong")
                        || d.contains("c=o")
                        || d.contains("c-o-c")
                        || d.contains("cf")
                        || d.contains("si-o")
                        || d.contains("phenyl")
                        || d.contains("c-cl")
                        || d.contains("amide i")
                        || d.contains("c≡n")
                    {
                        0.75
                    } else {
                        0.48
                    };
                    let gamma = if d.contains("broad")
                        || d.contains("o-h")
                        || d.contains("n-h")
                    {
                        30.0
                    } else if d.contains("sharp")
                        || d.contains("c≡n")
                        || d.contains("nco")
                        || d.contains("si-ch")
                    {
                        11.0
                    } else {
                        16.0
                    };
                    (wn, h, gamma)
                })
                .collect();
            synth_ftir(&title, 2, &peaks, 0.015, 0.02, 0.002)
        }
        SpectrumKind::Hnmr => {
            let mut peaks: Vec<(f64, f64, f64)> = poly
                .hnmr_peaks
                .iter()
                .enumerate()
                .map(|(i, &(ppm, _))| (ppm, (0.95 - (i as f64) * 0.12).max(0.35), 0.035))
                .collect();
            peaks.push((7.26, 0.25, 0.012)); // CDCl3 ref
            peaks.push((0.00, 0.12, 0.009)); // TMS
            synth_nmr(&title, SpectrumKind::Hnmr, 2, 0.0, 10.0, 3500, &peaks, 0.004)
        }
        SpectrumKind::Cnmr => {
            let mut peaks: Vec<(f64, f64, f64)> = poly
                .cnmr_peaks
                .iter()
                .enumerate()
                .map(|(i, &(ppm, _))| (ppm, (0.92 - (i as f64) * 0.08).max(0.40), 0.35))
                .collect();
            peaks.push((77.48, 0.45, 0.16));
            peaks.push((77.16, 0.50, 0.16));
            peaks.push((76.84, 0.45, 0.16));
            synth_nmr(&title, SpectrumKind::Cnmr, 2, 0.0, 210.0, 3800, &peaks, 0.004)
        }
        SpectrumKind::UvVis => unreachable!("polymer references carry no UV-Vis fingerprint"),
    };

    spec.metadata = vec![
        ("Polymer".into(), poly.name.into()),
        ("Repeat Unit".into(), poly.repeat_unit.into()),
        ("M₀".into(), format!("{:.2} g/mol", poly.m0)),
        ("Thermal".into(), poly.tg_tm.into()),
    ];
    spec
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solvent_and_impurity_labels_take_priority() {
        let label =
            best_short_peak_label_with_solvent(SpectrumKind::Hnmr, 7.26, SolventId::Cdcl3);
        assert!(label.contains("CDCl₃ (solvent)"), "got {label:?}");

        let label =
            best_short_peak_label_with_solvent(SpectrumKind::Hnmr, 1.56, SolventId::Cdcl3);
        assert!(label.contains("H₂O (impurity)"), "got {label:?}");

        let label =
            best_short_peak_label_with_solvent(SpectrumKind::Hnmr, 0.07, SolventId::Cdcl3);
        assert!(label.contains("Silicone grease (impurity)"), "got {label:?}");

        // A genuine polymer signal must not be hijacked by the impurity rules.
        let label =
            best_short_peak_label_with_solvent(SpectrumKind::Hnmr, 3.60, SolventId::Cdcl3);
        assert!(label.contains("O–CH₃"), "got {label:?}");
    }

    #[test]
    fn polymer_reference_uses_its_own_fingerprint() {
        let pvc = find_polymer_by_short("pvc").expect("PVC must exist");
        let spec = synthesize_polymer_reference(pvc, SpectrumKind::Ftir);
        assert!(spec.name.contains("PVC"), "got {:?}", spec.name);
        assert_eq!(spec.kind, SpectrumKind::Ftir);
        // C–Cl fingerprint bands at 615 / 690 cm⁻¹ must be present in the data.
        assert!(spec.x.len() > 100);
        let lowest = spec.x.first().copied().unwrap_or(0.0);
        let highest = spec.x.last().copied().unwrap_or(0.0);
        assert!(lowest <= 615.0 && highest >= 690.0);
    }

    #[test]
    fn polymer_reference_falls_back_to_ftir_without_nmr_peaks() {
        let ptfe = find_polymer_by_short("PTFE").expect("PTFE must exist");
        let spec = synthesize_polymer_reference(ptfe, SpectrumKind::Hnmr);
        assert_eq!(spec.kind, SpectrumKind::Ftir);
    }
}
