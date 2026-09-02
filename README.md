# Cardiac MRI Cine Strain Agent

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Cardiac MRI Cine Feature-Tracking (FT) Strain Analysis Engine
============================================================
Comprehensive quantitative myocardial biomechanics and multiparametric CMR
characterization engine for Global Longitudinal Strain (GLS), Global Circumferential
Strain (GCS), Global Radial Strain (GRS), strain rates (SRE/SRA), T1/ECV mapping,
LGE scar transmurality, and diagnostic phenotyping.

Mathematical Formulations:
  1. Lagrangian Myocardial Strain:
     epsilon(t) = (L(t) - L_0) / L_0 * 100%
  2. Extracellular Volume Fraction (ECV):
     ECV = (1 - Hematocrit) * [ (1/T1_myo_post - 1/T1_myo_pre) / (1/T1_blood_post - 1/T1_blood_pre) ] * 100%
  3. Left Ventricular Volumetrics:
     LVEF = (LVEDV - LVESV) / LVEDV * 100%
     Cardiac Output = (LVEDV - LVESV) * HeartRate / 1000  [L/min]
     Cardiac Index = Cardiac Output / BSA  [L/min/m^2]
  4. Apical Sparing Ratio (Relative Apical Strain):
     ASR = Average Apical Longitudinal Strain / (Average Basal Strain + Average Mid Strain)

LGETransmuralityAgent: late gadolinium enhancement transmurality and viability.

Per AHA 17-segment model, each segment's scar transmurality classifies
viability for revascularization decisions:

    0%            -> normal myocardium
    1-49%         -> hibernating / viable (revascularization benefit likely)
    >50% (>=50%)  -> non-viable scar

Segments map to coronary territories:
    LAD : 1,2,7,8,13,14,17 (+ apical 15 in wrap variants)
    RCA : 3,4,9,10,15
    LCx : 5,6,11,12,16

Outputs a per-territory viability report and an ASCII bull's-eye map.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`MyocardialFunctionTier`** — dedicated module for myocardial function tier evaluation and state verification.
- **`LGEPattern`** — dedicated module for l g e pattern evaluation and state verification.
- **`DiagnosticPhenotype`** — dedicated module for diagnostic phenotype evaluation and state verification.
- **`CineStrainInput`**: Input features for Cine CMR strain and multiparametric evaluation.
- **`CineStrainReport`**: Output dossier for cardiac MRI cine strain and tissue characterization.
- **`SegmentScar`** — dedicated module for segment scar evaluation and state verification.

---

## 📐 Mathematical Formulation & Logic

```text
  Mathematical Formulations:
  Cardiac Index = Cardiac Output / BSA  [L/min/m^2]
  calculated_ecv_pct: Optional[float]
  calculated_ecv = None
  ecv_val = calculate_ecv(
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --- <value> --study-id <value> --patient-id <value> --hr <value>
```

### Parameter Reference
- `---`: Specifies input measurement or parameter value.
- `--study-id`: Specifies input measurement or parameter value.
- `--patient-id`: Specifies input measurement or parameter value.
- `--hr`: Specifies input measurement or parameter value.
- `--bsa`: Specifies input measurement or parameter value.
- `--edv`: Specifies input measurement or parameter value.
- `--esv`: Specifies input measurement or parameter value.
- `--mass`: Specifies input measurement or parameter value.
- `--gls`: Specifies input measurement or parameter value.
- `--gcs`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `study_id` | Parameter / observation metric | Required |
| `patient_id` | Parameter / observation metric | Required |
| `heart_rate_bpm` | Parameter / observation metric | Required |
| `bsa_m2` | Parameter / observation metric | Required |
| `lvedv_ml` | Parameter / observation metric | Required |
| `lvesv_ml` | Parameter / observation metric | Required |
| `lv_mass_g` | Parameter / observation metric | Required |
| `gls_pct` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t cardiac-mri-cine-strain-agent .
docker run -p 8000:8000 cardiac-mri-cine-strain-agent
```
