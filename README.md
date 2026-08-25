# Cardiac MRI Cine Feature-Tracking Strain & Tissue Characterization Agent

A zero-dependency Python clinical biophysics platform for **Cardiovascular Magnetic Resonance (CMR) Cine Feature-Tracking (FT) Myocardial Strain Analysis**, **Native $T_1$ & Extracellular Volume (ECV) Mapping**, **Late Gadolinium Enhancement (LGE) Scar Transmurality**, and **Cardiomyopathy Phenotyping**.

---

## Myocardial Biomechanics & Mathematical Modeling

```
               [ Cine Long-Axis & Short-Axis MRI ]
                                |
             ----------------------------------------
            |                   |                    |
     [ Longitudinal ]    [ Circumferential ]     [ Radial ]
        GLS (%)             GCS (%)              GRS (%)
   (Normal: -18 to -25) (Normal: -19 to -26) (Normal: +35 to +55)
            |                   |                    |
             ----------------------------------------
                                |
      + [ Parametric Mapping: Native T1, ECV Fraction, LGE ]
                                |
            [ Multiparametric Diagnostic Phenotyping ]
```

### 1. Lagrangian Myocardial Deformation & Strain
- **Peak Systolic Strain ($\epsilon$):**
  $$\epsilon(t) = \frac{L(t) - L_0}{L_0} \times 100\%$$
- **Global Longitudinal Strain (GLS):** Shortening along the base-to-apex long axis (normally $-18\%$ to $-25\%$).
- **Global Circumferential Strain (GCS):** Myocardial circumferential shortening along the short-axis curvature (normally $-19\%$ to $-26\%$).
- **Global Radial Strain (GRS):** Transmural myocardial thickening (normally $+35\%$ to $+55\%$).

### 2. Extracellular Volume Fraction ($ECV$) Mapping
Quantifies diffuse interstitial myocardial fibrosis and amyloid infiltration:
$$ECV = (1 - \text{Hematocrit}) \times \frac{\left(\frac{1}{T_{1,\text{myo,post}}} - \frac{1}{T_{1,\text{myo,pre}}}\right)}{\left(\frac{1}{T_{1,\text{blood,post}}} - \frac{1}{T_{1,\text{blood,pre}}}\right)} \times 100\%$$
- **Normal:** $23\% - 28\%$
- **Diffuse Interstitial Fibrosis:** $30\% - 38\%$
- **Cardiac Amyloidosis:** $\ge 40\%$

### 3. Left Ventricular Hemodynamics & Volumetrics
- **Ejection Fraction (LVEF):** $\frac{\text{LVEDV} - \text{LVESV}}{\text{LVEDV}} \times 100\%$
- **Stroke Volume (SV):** $\text{LVEDV} - \text{LVESV} \quad [\text{mL}]$
- **Cardiac Output (CO):** $\frac{\text{SV} \times \text{HR}}{1000} \quad [\text{L/min}]$
- **Cardiac Index (CI):** $\frac{\text{CO}}{\text{BSA}} \quad [\text{L/min/m}^2]$

### 4. Apical Sparing Ratio (ASR)
Identifies cardiac amyloidosis ("cherry-on-top" bullseye plot):
$$\text{ASR} = \frac{|\text{Average Apical Longitudinal Strain}|}{\frac{|\text{Average Basal LS}| + |\text{Average Mid LS}|}{2}}$$
- $\text{ASR} > 1.0$ is highly sensitive and specific for Cardiac Amyloidosis.

---

## Installation & Quick Start

```bash
# Clone the repository
git clone https://github.com/example/cardiac-mri-cine-strain-agent.git
cd cardiac-mri-cine-strain-agent
```

### Python API Example
```python
from cardiac_mri_strain import CineStrainInput, evaluate_cine_strain, LGEPattern

study = CineStrainInput(
    study_id="CMR-AMYLOID-01",
    patient_id="PAT-102",
    lvedv_ml=130.0,
    lvesv_ml=60.0,
    lv_mass_g=180.0,
    gls_pct=-10.5,
    gcs_pct=-12.0,
    grs_pct=18.0,
    apical_ls_pct=-22.0,
    mid_ls_pct=-8.0,
    basal_ls_pct=-6.0,
    native_t1_ms=1160.0,
    post_contrast_t1_myo_ms=370.0,
    pre_contrast_t1_blood_ms=1550.0,
    post_contrast_t1_blood_ms=340.0,
    hematocrit_pct=37.0,
    lge_pattern=LGEPattern.DIFFUSE_SUBENDOCARDIAL,
)

report = evaluate_cine_strain(study)
print(f"Phenotype          : {report.diagnostic_phenotype}")
print(f"LVEF               : {report.lvef_pct:.1f}%")
print(f"GLS Function       : {report.gls_function_tier}")
print(f"Calculated ECV     : {report.calculated_ecv_pct:.1f}%")
print(f"Apical Sparing ASR : {report.apical_sparing_ratio:.2f}")
```

---

## CLI Usage

### 1. Evaluate Single Study
```bash
python cli.py evaluate \
    --study-id "CMR-EXP-01" \
    --edv 150 \
    --esv 60 \
    --mass 120 \
    --gls -20.5 \
    --gcs -22.0 \
    --grs 42.0 \
    --native-t1 1005 \
    --lge-pattern none
```

### 2. JSON Output Mode
```bash
python cli.py evaluate --edv 150 --esv 60 --gls -20.5 --json
```

### 3. Interactive Clinical Wizard
```bash
python cli.py interactive
```

### 4. Batch CSV Processing
```bash
python cli.py batch -i sample.csv -o cmr_results.csv
```

### 5. Normal Reference Ranges
```bash
python cli.py norms
```

---

## Test Suite Execution

Run the complete 28-case unit test suite:
```bash
python -m unittest test_cardiac_mri_strain.py
```

All 28 tests confirm 100% pass rates across myocardial strain classifications, ECV formulas, amyloidosis apical sparing detection, ischemic scar patterns, and CLI workflows.

---

## References
1. Scatteia A, et al. Comprehensive review of CMR feature tracking for myocardial strain analysis: clinical applications and technical considerations. *European Heart Journal - Cardiovascular Imaging*. 2017;18(12):1303–1314.
2. Messroghli DR, et al. Clinical recommendations for cardiovascular magnetic resonance mapping of T1, T2, T2* and extracellular volume: A consensus statement by the SCMR. *J Cardiovasc Magn Reson*. 2017;19(1):75.
3. Phelan D, et al. Relative apical sparing of longitudinal strain using two-dimensional speckle-tracking echocardiography is both sensitive and specific for the diagnosis of cardiac amyloidosis. *Heart*. 2012;98(19):1442–1448.
