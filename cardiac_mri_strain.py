#!/usr/bin/env python3
"""
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
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union


class MyocardialFunctionTier(str, Enum):
    NORMAL = "Normal"
    MILDLY_REDUCED = "Mildly Reduced"
    MODERATELY_REDUCED = "Moderately Reduced"
    SEVERELY_REDUCED = "Severely Reduced"


class LGEPattern(str, Enum):
    NONE = "none"
    SUBENDOCARDIAL = "subendocardial"      # Ischemic / CAD
    TRANSMURAL = "transmural"              # Ischemic transmural infarct
    MID_WALL = "mid_wall"                  # Non-ischemic / DCM / Fibrosis
    SUBEPICARDIAL = "subepicardial"        # Myocarditis / Sarcoidosis
    DIFFUSE_SUBENDOCARDIAL = "diffuse"     # Cardiac Amyloidosis


class DiagnosticPhenotype(str, Enum):
    NORMAL = "Normal Myocardial Mechanics"
    ISCHEMIC_INFARCT = "Ischemic Cardiomyopathy / Myocardial Infarction"
    CARDIAC_AMYLOIDOSIS = "Cardiac Amyloidosis (Apical Sparing Strain Pattern)"
    HYPERTROPHIC_CM = "Hypertrophic Cardiomyopathy (HCM)"
    DILATED_CM = "Dilated Cardiomyopathy (DCM)"
    ACUTE_MYOCARDITIS = "Acute Myocarditis"
    DIASTOLIC_DYSFUNCTION = "Isolated Diastolic Impairment"


@dataclass
class CineStrainInput:
    """Input features for Cine CMR strain and multiparametric evaluation."""
    study_id: str = "CMR-EXAM-001"
    patient_id: Optional[str] = None
    heart_rate_bpm: float = 72.0
    bsa_m2: float = 1.85
    
    # Volumetrics (mL)
    lvedv_ml: float = 145.0
    lvesv_ml: float = 58.0
    lv_mass_g: float = 120.0
    
    # Global Peak Strains (%)
    gls_pct: float = -20.5                # Global Longitudinal Strain (-25% to -18% normal)
    gcs_pct: float = -22.0                # Global Circumferential Strain (-26% to -19% normal)
    grs_pct: float = 42.0                 # Global Radial Strain (+35% to +55% normal)
    
    # Segmental Longitudinal Strains for Apical Sparing Ratio (optional, %)
    apical_ls_pct: Optional[float] = None
    mid_ls_pct: Optional[float] = None
    basal_ls_pct: Optional[float] = None
    
    # Strain Rates (s^-1)
    peak_systolic_sr_s1: float = -1.25     # Normal <= -1.0 s^-1
    early_diastolic_sr_e_s1: float = 1.45  # Normal >= 1.2 s^-1
    late_diastolic_sr_a_s1: float = 0.95   # Normal ~0.8-1.1 s^-1
    
    # Mechanical Dyssynchrony (Standard deviation of time-to-peak strain in ms)
    sd_ttp_ms: float = 28.0                # Normal < 35 ms, dyssynchrony > 45 ms
    
    # Tissue Characterization (T1 & ECV)
    native_t1_ms: float = 1005.0          # Normal 950-1060 ms at 1.5T
    post_contrast_t1_myo_ms: Optional[float] = 480.0
    pre_contrast_t1_blood_ms: Optional[float] = 1580.0
    post_contrast_t1_blood_ms: Optional[float] = 320.0
    hematocrit_pct: float = 42.0          # Hematocrit (35-50%)
    
    # Late Gadolinium Enhancement (LGE)
    lge_pattern: Union[LGEPattern, str] = LGEPattern.NONE
    lge_transmurality_pct: float = 0.0     # 0 - 100%
    lge_scar_mass_pct: float = 0.0         # % of LV mass


@dataclass
class CineStrainReport:
    """Output dossier for cardiac MRI cine strain and tissue characterization."""
    study_id: str
    patient_id: Optional[str]
    lvef_pct: float
    stroke_volume_ml: float
    cardiac_output_l_min: float
    cardiac_index_l_min_m2: float
    lv_mass_index_g_m2: float
    
    # Strain Analysis
    gls_pct: float
    gcs_pct: float
    grs_pct: float
    gls_function_tier: str
    gcs_function_tier: str
    grs_function_tier: str
    overall_contractility_tier: str
    
    # Diastolic & Dyssynchrony
    sre_sra_ratio: float
    diastolic_function_grade: str
    mechanical_dyssynchrony: str
    sd_ttp_ms: float
    
    # Tissue Characterization
    native_t1_ms: float
    t1_status: str
    calculated_ecv_pct: Optional[float]
    ecv_status: str
    lge_summary: str
    
    # Apical Sparing & Diagnostic Phenotyping
    apical_sparing_ratio: Optional[float]
    diagnostic_phenotype: str
    clinical_findings: List[str] = field(default_factory=list)
    remediation_recommendations: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def calculate_ecv(
    t1_myo_pre: float,
    t1_myo_post: float,
    t1_blood_pre: float,
    t1_blood_post: float,
    hematocrit_pct: float,
) -> float:
    """
    Calculate Extracellular Volume Fraction (ECV) from pre- and post-T1 values and hematocrit.
    """
    if t1_myo_pre <= 0 or t1_myo_post <= 0 or t1_blood_pre <= 0 or t1_blood_post <= 0:
        raise ValueError("T1 relaxation times must be positive non-zero values.")
    if hematocrit_pct <= 10.0 or hematocrit_pct >= 70.0:
        raise ValueError(f"Hematocrit ({hematocrit_pct}%) is outside plausible range [10, 70].")

    delta_r1_myo = (1.0 / (t1_myo_post / 1000.0)) - (1.0 / (t1_myo_pre / 1000.0))
    delta_r1_blood = (1.0 / (t1_blood_post / 1000.0)) - (1.0 / (t1_blood_pre / 1000.0))

    if delta_r1_blood <= 0:
        raise ValueError("Post-contrast blood T1 shortening is insufficient to compute ECV.")

    hct_fraction = hematocrit_pct / 100.0
    ecv = (1.0 - hct_fraction) * (delta_r1_myo / delta_r1_blood) * 100.0
    return ecv


def classify_longitudinal_strain(gls: float) -> MyocardialFunctionTier:
    """Classify peak GLS (magnitudes more negative indicate superior deformation)."""
    # GLS is negative, e.g. -20%
    if gls <= -18.0:
        return MyocardialFunctionTier.NORMAL
    elif gls <= -15.0:
        return MyocardialFunctionTier.MILDLY_REDUCED
    elif gls <= -11.0:
        return MyocardialFunctionTier.MODERATELY_REDUCED
    else:
        return MyocardialFunctionTier.SEVERELY_REDUCED


def classify_circumferential_strain(gcs: float) -> MyocardialFunctionTier:
    """Classify peak GCS."""
    if gcs <= -19.0:
        return MyocardialFunctionTier.NORMAL
    elif gcs <= -15.0:
        return MyocardialFunctionTier.MILDLY_REDUCED
    elif gcs <= -10.0:
        return MyocardialFunctionTier.MODERATELY_REDUCED
    else:
        return MyocardialFunctionTier.SEVERELY_REDUCED


def classify_radial_strain(grs: float) -> MyocardialFunctionTier:
    """Classify peak GRS (positive percentage)."""
    if grs >= 35.0:
        return MyocardialFunctionTier.NORMAL
    elif grs >= 25.0:
        return MyocardialFunctionTier.MILDLY_REDUCED
    elif grs >= 15.0:
        return MyocardialFunctionTier.MODERATELY_REDUCED
    else:
        return MyocardialFunctionTier.SEVERELY_REDUCED


def evaluate_cine_strain(inp: CineStrainInput) -> CineStrainReport:
    """
    Comprehensive pipeline evaluating CMR Cine Feature-Tracking strain and tissue maps.
    """
    # Input validation
    if inp.lvedv_ml <= 0 or inp.lvesv_ml <= 0:
        raise ValueError("Left ventricular volumes must be positive non-zero values.")
    if inp.lvedv_ml <= inp.lvesv_ml:
        raise ValueError(f"LVEDV ({inp.lvedv_ml} mL) must be strictly greater than LVESV ({inp.lvesv_ml} mL).")
    if inp.bsa_m2 <= 0.5 or inp.bsa_m2 >= 4.0:
        raise ValueError(f"BSA ({inp.bsa_m2} m^2) is outside physiological range [0.5, 4.0].")

    # 1. Volumetrics & Hemodynamics
    stroke_vol = inp.lvedv_ml - inp.lvesv_ml
    lvef = (stroke_vol / inp.lvedv_ml) * 100.0
    cardiac_out = (stroke_vol * inp.heart_rate_bpm) / 1000.0
    cardiac_idx = cardiac_out / inp.bsa_m2
    lv_mass_idx = inp.lv_mass_g / inp.bsa_m2

    # 2. Strain Classifications
    gls_tier = classify_longitudinal_strain(inp.gls_pct)
    gcs_tier = classify_circumferential_strain(inp.gcs_pct)
    grs_tier = classify_radial_strain(inp.grs_pct)

    # Composite contractility tier
    tier_ranks = {
        MyocardialFunctionTier.NORMAL: 0,
        MyocardialFunctionTier.MILDLY_REDUCED: 1,
        MyocardialFunctionTier.MODERATELY_REDUCED: 2,
        MyocardialFunctionTier.SEVERELY_REDUCED: 3,
    }
    max_rank = max(tier_ranks[gls_tier], tier_ranks[gcs_tier], tier_ranks[grs_tier])
    overall_contractility = [k for k, v in tier_ranks.items() if v == max_rank][0].value

    # 3. Diastolic Strain Rates & Dyssynchrony
    sre_sra = round(inp.early_diastolic_sr_e_s1 / max(0.01, inp.late_diastolic_sr_a_s1), 2)
    if inp.early_diastolic_sr_e_s1 >= 1.20:
        diastolic_grade = "Normal Active Diastolic Relaxation"
    elif inp.early_diastolic_sr_e_s1 >= 0.85:
        diastolic_grade = "Grade 1 (Impaired Relaxation)"
    else:
        diastolic_grade = "Grade 2/3 (Advanced Diastolic Dysfunction)"

    if inp.sd_ttp_ms > 45.0:
        dyssynchrony_status = "Significant Left Ventricular Mechanical Dyssynchrony"
    elif inp.sd_ttp_ms > 35.0:
        dyssynchrony_status = "Borderline Mechanical Dyssynchrony"
    else:
        dyssynchrony_status = "Synchronous Myocardial Activation"

    # 4. T1 and ECV Mapping
    if inp.native_t1_ms > 1070.0:
        t1_status = f"Significantly Elevated Native T1 ({inp.native_t1_ms:.0f} ms) - Extracellular expansion / edema"
    elif inp.native_t1_ms < 930.0:
        t1_status = f"Low Native T1 ({inp.native_t1_ms:.0f} ms) - Suggests Myocardial Iron Overload or Fabry Disease"
    else:
        t1_status = f"Normal Native T1 ({inp.native_t1_ms:.0f} ms)"

    calculated_ecv = None
    ecv_status = "ECV not assessed"
    if (
        inp.post_contrast_t1_myo_ms is not None
        and inp.pre_contrast_t1_blood_ms is not None
        and inp.post_contrast_t1_blood_ms is not None
    ):
        try:
            ecv_val = calculate_ecv(
                t1_myo_pre=inp.native_t1_ms,
                t1_myo_post=inp.post_contrast_t1_myo_ms,
                t1_blood_pre=inp.pre_contrast_t1_blood_ms,
                t1_blood_post=inp.post_contrast_t1_blood_ms,
                hematocrit_pct=inp.hematocrit_pct,
            )
            calculated_ecv = round(ecv_val, 1)
            if ecv_val >= 40.0:
                ecv_status = f"Markedly Elevated ECV ({ecv_val:.1f}%) - Highly characteristic of Cardiac Amyloidosis"
            elif ecv_val >= 30.0:
                ecv_status = f"Elevated ECV ({ecv_val:.1f}%) - Diffuse Interstitial Myocardial Fibrosis"
            elif ecv_val < 22.0:
                ecv_status = f"Low ECV ({ecv_val:.1f}%)"
            else:
                ecv_status = f"Normal ECV ({ecv_val:.1f}%)"
        except ValueError as e:
            ecv_status = f"ECV Calculation Error: {e}"

    # 5. LGE Pattern Description
    lge_pat = inp.lge_pattern.value if isinstance(inp.lge_pattern, LGEPattern) else str(inp.lge_pattern).lower()
    if lge_pat in ("none", "no", "false", ""):
        lge_summary = "No late gadolinium enhancement detected (no focal scar)."
    else:
        lge_summary = f"LGE Present: Pattern [{lge_pat}], Transmurality [{inp.lge_transmurality_pct:.0f}%], Scar Mass [{inp.lge_scar_mass_pct:.1f}% LV mass]."

    # 6. Apical Sparing Ratio (ASR)
    # ASR = Apical LS / (Basal LS + Mid LS)
    asr = None
    if inp.apical_ls_pct is not None and inp.basal_ls_pct is not None and inp.mid_ls_pct is not None:
        denom = abs(inp.basal_ls_pct) + abs(inp.mid_ls_pct)
        if denom > 0:
            asr = round(abs(inp.apical_ls_pct) / (denom / 2.0), 2)

    # 7. Diagnostic Phenotyping & Findings
    findings: List[str] = []
    recs: List[str] = []

    # Phenotyping heuristic
    if lge_pat in ("subendocardial", "transmural") or (inp.lge_transmurality_pct > 25.0 and lge_pat != "diffuse"):
        phenotype = DiagnosticPhenotype.ISCHEMIC_INFARCT.value
        findings.append(f"Ischemic territorial scar pattern with {inp.lge_transmurality_pct:.0f}% transmurality.")
        recs.append("Assess coronary revascularization viability if transmurality < 50%.")
    elif (asr is not None and asr > 1.0) or (calculated_ecv is not None and calculated_ecv >= 40.0) or lge_pat == "diffuse":
        phenotype = DiagnosticPhenotype.CARDIAC_AMYLOIDOSIS.value
        findings.append("Prominent apical sparing strain pattern and severe extracellular matrix expansion.")
        recs.append("Recommend 99mTc-PYP scintigraphy and serum free light chain electrophoresis for Amyloidosis typing.")
    elif lv_mass_idx > 115.0 and gls_tier in (MyocardialFunctionTier.MODERATELY_REDUCED, MyocardialFunctionTier.SEVERELY_REDUCED):
        phenotype = DiagnosticPhenotype.HYPERTROPHIC_CM.value
        findings.append("Marked LV hypertrophy with reduced regional longitudinal strain.")
        recs.append("Genetic testing and sudden cardiac death risk stratification.")
    elif lvef < 50.0 and gls_tier != MyocardialFunctionTier.NORMAL:
        phenotype = DiagnosticPhenotype.DILATED_CM.value
        findings.append("Global systolic dysfunction with diminished multi-directional strain indices.")
        recs.append("Guideline-directed medical therapy (GDMT) optimization for HFrEF.")
    elif lge_pat == "subepicardial" or (inp.native_t1_ms > 1080.0 and lvef >= 50.0):
        phenotype = DiagnosticPhenotype.ACUTE_MYOCARDITIS.value
        findings.append("Subepicardial LGE with elevated native T1 compatible with acute inflammatory myocarditis.")
        recs.append("Activity restriction for 3-6 months and serial CMR follow-up.")
    elif diastolic_grade != "Normal Active Diastolic Relaxation" and lvef >= 50.0:
        phenotype = DiagnosticPhenotype.DIASTOLIC_DYSFUNCTION.value
        findings.append("Preserved LVEF with impaired early diastolic strain rate relaxation.")
        recs.append("Evaluate for HFpEF and manage afterload / volume status.")
    else:
        phenotype = DiagnosticPhenotype.NORMAL.value
        findings.append("Normal left ventricular volumes, systolic ejection fraction, and multi-directional strain.")
        recs.append("Routine clinical follow-up.")

    # Additional alerts
    if lvef < 35.0:
        findings.append(f"WARNING: Severe LVEF depression ({lvef:.1f}%) - ICD eligibility evaluation indicated.")
    if inp.sd_ttp_ms > 45.0 and lvef < 35.0:
        recs.append("Consider cardiac resynchronization therapy (CRT) evaluation.")

    return CineStrainReport(
        study_id=inp.study_id,
        patient_id=inp.patient_id,
        lvef_pct=round(lvef, 1),
        stroke_volume_ml=round(stroke_vol, 1),
        cardiac_output_l_min=round(cardiac_out, 2),
        cardiac_index_l_min_m2=round(cardiac_idx, 2),
        lv_mass_index_g_m2=round(lv_mass_idx, 1),
        gls_pct=round(inp.gls_pct, 1),
        gcs_pct=round(inp.gcs_pct, 1),
        grs_pct=round(inp.grs_pct, 1),
        gls_function_tier=gls_tier.value,
        gcs_function_tier=gcs_tier.value,
        grs_function_tier=grs_tier.value,
        overall_contractility_tier=overall_contractility,
        sre_sra_ratio=sre_sra,
        diastolic_function_grade=diastolic_grade,
        mechanical_dyssynchrony=dyssynchrony_status,
        sd_ttp_ms=round(inp.sd_ttp_ms, 1),
        native_t1_ms=round(inp.native_t1_ms, 0),
        t1_status=t1_status,
        calculated_ecv_pct=calculated_ecv,
        ecv_status=ecv_status,
        lge_summary=lge_summary,
        apical_sparing_ratio=asr,
        diagnostic_phenotype=phenotype,
        clinical_findings=findings,
        remediation_recommendations=recs,
    )


def calculate_metrics(**kwargs) -> Dict[str, Any]:
    """
    Standard top-level interface for single evaluation and testing.
    """
    def _float(key: str, default: float) -> float:
        val = kwargs.get(key)
        if val is None:
            return default
        try:
            return float(val)
        except (ValueError, TypeError):
            return default

    def _str(key: str, default: str) -> str:
        val = kwargs.get(key)
        return str(val) if val is not None else default

    study_id = _str("study_id", _str("id", "CMR-001"))
    patient_id = kwargs.get("patient_id") or kwargs.get("Patient")
    
    hr = _float("heart_rate_bpm", _float("heart_rate", 72.0))
    bsa = _float("bsa_m2", _float("bsa", 1.85))
    edv = _float("lvedv_ml", _float("edv", _float("primary_metric", 145.0)))
    esv = _float("lvesv_ml", _float("esv", _float("secondary_metric", 58.0)))
    mass = _float("lv_mass_g", _float("lv_mass", 120.0))
    
    gls = _float("gls_pct", _float("gls", -20.5))
    gcs = _float("gcs_pct", _float("gcs", -22.0))
    grs = _float("grs_pct", _float("grs", 42.0))
    
    apical_ls = kwargs.get("apical_ls_pct")
    mid_ls = kwargs.get("mid_ls_pct")
    basal_ls = kwargs.get("basal_ls_pct")

    t1 = _float("native_t1_ms", _float("native_t1", 1005.0))
    t1_myo_post = kwargs.get("post_contrast_t1_myo_ms") or kwargs.get("t1_myo_post")
    t1_blood_pre = kwargs.get("pre_contrast_t1_blood_ms") or kwargs.get("t1_blood_pre")
    t1_blood_post = kwargs.get("post_contrast_t1_blood_ms") or kwargs.get("t1_blood_post")
    hct = _float("hematocrit_pct", _float("hct", 42.0))
    
    lge_pattern = _str("lge_pattern", "none")
    lge_trans = _float("lge_transmurality_pct", 0.0)
    lge_scar = _float("lge_scar_mass_pct", 0.0)

    inp = CineStrainInput(
        study_id=study_id,
        patient_id=str(patient_id) if patient_id is not None else None,
        heart_rate_bpm=hr,
        bsa_m2=bsa,
        lvedv_ml=edv,
        lvesv_ml=esv,
        lv_mass_g=mass,
        gls_pct=gls,
        gcs_pct=gcs,
        grs_pct=grs,
        apical_ls_pct=float(apical_ls) if apical_ls is not None else None,
        mid_ls_pct=float(mid_ls) if mid_ls is not None else None,
        basal_ls_pct=float(basal_ls) if basal_ls is not None else None,
        native_t1_ms=t1,
        post_contrast_t1_myo_ms=float(t1_myo_post) if t1_myo_post is not None else 480.0,
        pre_contrast_t1_blood_ms=float(t1_blood_pre) if t1_blood_pre is not None else 1580.0,
        post_contrast_t1_blood_ms=float(t1_blood_post) if t1_blood_post is not None else 320.0,
        hematocrit_pct=hct,
        lge_pattern=lge_pattern,
        lge_transmurality_pct=lge_trans,
        lge_scar_mass_pct=lge_scar,
    )

    report = evaluate_cine_strain(inp)
    res = report.to_dict()
    res["tool"] = "cardiac-mri-cine-strain-agent"
    res["score"] = report.lvef_pct
    res["classification"] = report.diagnostic_phenotype
    res["clinical_recommendation"] = "; ".join(report.remediation_recommendations)
    return res


def process_batch(input_csv: str, output_csv: str) -> int:
    """
    Batch process patient CMR examinations from CSV file.
    """
    with open(input_csv, mode="r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    out_fields = fieldnames + [
        "lvef_pct",
        "stroke_volume_ml",
        "cardiac_index_l_min_m2",
        "gls_pct",
        "gls_function_tier",
        "calculated_ecv_pct",
        "diagnostic_phenotype",
        "overall_contractility_tier",
    ]
    dedup_fields = []
    for fn in out_fields:
        if fn not in dedup_fields:
            dedup_fields.append(fn)

    out_rows = []
    for r in rows:
        calc_res = calculate_metrics(**r)
        row_dict = dict(r)
        row_dict["lvef_pct"] = calc_res["lvef_pct"]
        row_dict["stroke_volume_ml"] = calc_res["stroke_volume_ml"]
        row_dict["cardiac_index_l_min_m2"] = calc_res["cardiac_index_l_min_m2"]
        row_dict["gls_pct"] = calc_res["gls_pct"]
        row_dict["gls_function_tier"] = calc_res["gls_function_tier"]
        row_dict["calculated_ecv_pct"] = calc_res["calculated_ecv_pct"]
        row_dict["diagnostic_phenotype"] = calc_res["diagnostic_phenotype"]
        row_dict["overall_contractility_tier"] = calc_res["overall_contractility_tier"]
        out_rows.append(row_dict)

    with open(output_csv, mode="w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=dedup_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    return len(out_rows)
