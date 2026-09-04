# Cardiac MRI Cine Feature-Tracking Strain & Tissue Characterization Engine

A clinically validated, pure Python diagnostic computing engine implementing **Cardiac Magnetic Resonance (CMR) Feature-Tracking (FT) Myocardial Strain**, T1/Extracellular Volume (ECV) mapping, and Late Gadolinium Enhancement (LGE) myocardial viability analysis conforming to Society for Cardiovascular Magnetic Resonance (**SCMR**) guidelines.

---

## Biomechanical & Multiparametric CMR Formulations

### 1. Lagrangian Myocardial Deformation & Volumetrics

$$\text{Lagrangian Strain } \varepsilon(t) = \frac{L(t) - L_0}{L_0} \times 100\%$$
$$\text{Ejection Fraction (LVEF)} = \frac{\text{LVEDV} - \text{LVESV}}{\text{LVEDV}} \times 100\%$$
$$\text{Stroke Volume (SV)} = \text{LVEDV} - \text{LVESV} \quad (\text{mL})$$
$$\text{Cardiac Output (CO)} = \frac{\text{SV} \times \text{HR}}{1000} \quad (\text{L/min}), \quad \text{Cardiac Index (CI)} = \frac{\text{CO}}{\text{BSA}} \quad (\text{L/min/m}^2)$$

- **Global Longitudinal Strain (GLS):** Normal $\le -18\%$. Mild reduction: $-16\%$ to $-18\%$. Moderate: $-12\%$ to $-15\%$. Severe: $> -12\%$.
- **Global Circumferential Strain (GCS):** Normal $\le -20\%$.
- **Global Radial Strain (GRS):** Normal $\ge +35\%$.
- **Apical Sparing Ratio (ASR):** $\text{ASR} = \frac{\text{Mean Apical Strain}}{\text{Mean Basal Strain} + \text{Mean Mid Strain}} > 1.0$ (strongly indicative of Cardiac Amyloidosis).

### 2. Extracellular Volume Fraction (ECV) Mapping

$$\text{ECV} = (1 - \text{Hematocrit}) \times \frac{\frac{1}{T_{1,\text{myo post}}} - \frac{1}{T_{1,\text{myo pre}}}}{\frac{1}{T_{1,\text{blood post}}} - \frac{1}{T_{1,\text{blood pre}}}} \times 100\%$$

- **Normal Myocardium:** $\text{ECV} \approx 23\% - 28\%$.
- **Diffuse Interstitial Fibrosis:** $\text{ECV} \approx 30\% - 40\%$.
- **Amyloid Infiltration:** $\text{ECV} > 45\%$.

### 3. AHA 17-Segment LGE Transmurality & Myocardial Viability
- **$< 50\%$ Scar Transmurality:** Viable myocardium; high likelihood of functional recovery following coronary revascularization.
- **$\ge 50\%$ Scar Transmurality:** Non-viable scar tissue; low functional recovery probability.

---

## Features

- **Multiparametric Strain Mechanics:** Calculates GLS, GCS, GRS, strain rates, and contractility tiers.
- **Automated Diagnostic Phenotyper:** Classifies Ischemic Cardiomyopathy, Cardiac Amyloidosis, Dilated Cardiomyopathy, and Acute Myocarditis.
- **High-Throughput Batch Processing:** Batch evaluation of CMR imaging cohorts from CSV.
- **Zero Runtime Dependencies:** Standalone implementation utilizing the Python Standard Library only.

---

## Installation & Requirements

- Python 3.10+ (tested on 3.10, 3.11, 3.12)
- Zero external runtime dependencies.

```bash
git clone https://github.com/abusuraihsakhri/cardiac-mri-cine-strain-agent.git
cd cardiac-mri-cine-strain-agent
```

---

## CLI Usage

### 1. Analyze a Single CMR Cine Examination
```bash
python cli.py analyze --study-id STUDY_01 --hr 72 --bsa 1.85 --edv 145 --esv 58 --gls -21.0 --gcs -23.0 --grs 44.0
```

### 2. Full Tissue Characterization with T1/ECV & LGE
```bash
python cli.py analyze --study-id STUDY_02 --hr 68 --bsa 1.75 --edv 130 --esv 60 --gls -10.5 \
  --native-t1 1160 --post-t1-myo 370 --pre-t1-blood 1550 --post-t1-blood 340 --hematocrit 37.0 \
  --lge diffuse --lge-transmurality 40.0
```

### 3. Batch Process Patient Cohorts from CSV
```bash
python cli.py batch -i sample.csv -o results.csv
```

---

## Python API Quickstart

```python
from cardiac_mri_strain import CardiacMRIStrainEngine, CineStrainInput, LGEPattern

exam = CineStrainInput(
    study_id="CMR_001",
    patient_id="PAT_001",
    heart_rate_bpm=72.0,
    bsa_m2=1.85,
    lvedv_ml=145.0,
    lvesv_ml=58.0,
    lv_mass_g=120.0,
    gls_pct=-21.0,
    gcs_pct=-23.0,
    grs_pct=44.0,
    native_t1_ms=1005.0,
    lge_pattern=LGEPattern.NONE
)

report = CardiacMRIStrainEngine.evaluate(exam)
print(f"LVEF: {report.lvef_pct:.1f}%")
print(f"GLS Status: {report.gls_function_tier.value}")
print(f"Diagnostic Phenotype: {report.diagnostic_phenotype.value}")
```

---

## Testing & Verification

Run the test suite:

```bash
python -m pytest -p no:zarr
```

