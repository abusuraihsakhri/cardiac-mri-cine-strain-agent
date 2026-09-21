#!/usr/bin/env python3
"""Quantitative cardiac MRI cine/strain calculations.

This module accepts already-measured CMR values and derives volumetric,
hemodynamic, strain-summary, and extracellular-volume metrics. It does not
segment images, perform feature tracking, make diagnoses, or recommend therapy.
"""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Union


REFERENCE_NOTE = (
    "Strain, native T1, ECV, and related reference intervals depend on scanner, "
    "field strength, acquisition sequence, post-processing software, and local "
    "validation. Interpret quantitative values against an appropriate local or "
    "method-specific reference interval."
)


class MyocardialFunctionTier(str, Enum):
    """Legacy illustrative strain bands retained for API compatibility."""

    NORMAL = "Within illustrative range"
    MILDLY_REDUCED = "Mildly outside illustrative range"
    MODERATELY_REDUCED = "Moderately outside illustrative range"
    SEVERELY_REDUCED = "Markedly outside illustrative range"


class LGEPattern(str, Enum):
    NONE = "none"
    SUBENDOCARDIAL = "subendocardial"
    TRANSMURAL = "transmural"
    MID_WALL = "mid_wall"
    SUBEPICARDIAL = "subepicardial"
    DIFFUSE_SUBENDOCARDIAL = "diffuse"


class DiagnosticPhenotype(str, Enum):
    """Compatibility enum; this package no longer infers disease phenotypes."""

    NOT_INFERRED = "Not inferred from quantitative inputs"
    NORMAL = "Normal Myocardial Mechanics"
    ISCHEMIC_INFARCT = "Ischemic Cardiomyopathy / Myocardial Infarction"
    CARDIAC_AMYLOIDOSIS = "Cardiac Amyloidosis (Apical Sparing Strain Pattern)"
    HYPERTROPHIC_CM = "Hypertrophic Cardiomyopathy (HCM)"
    DILATED_CM = "Dilated Cardiomyopathy (DCM)"
    ACUTE_MYOCARDITIS = "Acute Myocarditis"
    DIASTOLIC_DYSFUNCTION = "Isolated Diastolic Impairment"


@dataclass
class CineStrainInput:
    """User-supplied quantitative CMR measurements."""

    study_id: str = "CMR-EXAM-001"
    patient_id: Optional[str] = None
    heart_rate_bpm: float = 72.0
    bsa_m2: float = 1.85
    lvedv_ml: float = 145.0
    lvesv_ml: float = 58.0
    lv_mass_g: float = 120.0
    gls_pct: float = -20.5
    gcs_pct: float = -22.0
    grs_pct: float = 42.0
    apical_ls_pct: Optional[float] = None
    mid_ls_pct: Optional[float] = None
    basal_ls_pct: Optional[float] = None
    peak_systolic_sr_s1: float = -1.25
    early_diastolic_sr_e_s1: float = 1.45
    late_diastolic_sr_a_s1: float = 0.95
    sd_ttp_ms: float = 28.0
    native_t1_ms: Optional[float] = 1005.0
    post_contrast_t1_myo_ms: Optional[float] = 480.0
    pre_contrast_t1_blood_ms: Optional[float] = 1580.0
    post_contrast_t1_blood_ms: Optional[float] = 320.0
    hematocrit_pct: Optional[float] = 42.0
    lge_pattern: Union[LGEPattern, str] = LGEPattern.NONE
    lge_transmurality_pct: float = 0.0
    lge_scar_mass_pct: float = 0.0


@dataclass
class CineStrainReport:
    study_id: str
    patient_id: Optional[str]
    lvef_pct: float
    stroke_volume_ml: float
    lvedv_index_ml_m2: float
    lvesv_index_ml_m2: float
    stroke_volume_index_ml_m2: float
    cardiac_output_l_min: float
    cardiac_index_l_min_m2: float
    lv_mass_index_g_m2: float
    gls_pct: float
    gcs_pct: float
    grs_pct: float
    gls_function_tier: str
    gcs_function_tier: str
    grs_function_tier: str
    overall_contractility_tier: str
    sre_sra_ratio: Optional[float]
    diastolic_function_grade: str
    mechanical_dyssynchrony: str
    sd_ttp_ms: float
    native_t1_ms: Optional[float]
    t1_status: str
    calculated_ecv_pct: Optional[float]
    ecv_status: str
    lge_summary: str
    apical_sparing_ratio: Optional[float]
    diagnostic_phenotype: str
    clinical_findings: List[str] = field(default_factory=list)
    remediation_recommendations: List[str] = field(default_factory=list)
    reference_note: str = REFERENCE_NOTE

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def _require_positive(name: str, value: float) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be greater than zero.")


def calculate_ecv(
    t1_myo_pre: float,
    t1_myo_post: float,
    t1_blood_pre: float,
    t1_blood_post: float,
    hematocrit_pct: float,
) -> float:
    """Calculate ECV (%) from paired myocardial/blood T1 values and hematocrit."""

    for name, value in (
        ("Pre-contrast myocardial T1", t1_myo_pre),
        ("Post-contrast myocardial T1", t1_myo_post),
        ("Pre-contrast blood T1", t1_blood_pre),
        ("Post-contrast blood T1", t1_blood_post),
    ):
        _require_positive(name, value)
    if not 0.0 < hematocrit_pct < 100.0:
        raise ValueError("Hematocrit must be between 0 and 100 percent.")

    delta_r1_myo = (1000.0 / t1_myo_post) - (1000.0 / t1_myo_pre)
    delta_r1_blood = (1000.0 / t1_blood_post) - (1000.0 / t1_blood_pre)
    if delta_r1_blood <= 0:
        raise ValueError("Post-contrast blood R1 must exceed pre-contrast blood R1.")

    ecv = (1.0 - hematocrit_pct / 100.0) * (delta_r1_myo / delta_r1_blood) * 100.0
    if ecv < 0:
        raise ValueError("Calculated ECV is negative; check T1 inputs and acquisition pairing.")
    return ecv


def classify_longitudinal_strain(gls: float) -> MyocardialFunctionTier:
    """Legacy illustrative banding; not a universal CMR reference standard."""
    if gls <= -18.0:
        return MyocardialFunctionTier.NORMAL
    if gls <= -15.0:
        return MyocardialFunctionTier.MILDLY_REDUCED
    if gls <= -11.0:
        return MyocardialFunctionTier.MODERATELY_REDUCED
    return MyocardialFunctionTier.SEVERELY_REDUCED


def classify_circumferential_strain(gcs: float) -> MyocardialFunctionTier:
    """Legacy illustrative banding; not a universal CMR reference standard."""
    if gcs <= -19.0:
        return MyocardialFunctionTier.NORMAL
    if gcs <= -15.0:
        return MyocardialFunctionTier.MILDLY_REDUCED
    if gcs <= -10.0:
        return MyocardialFunctionTier.MODERATELY_REDUCED
    return MyocardialFunctionTier.SEVERELY_REDUCED


def classify_radial_strain(grs: float) -> MyocardialFunctionTier:
    """Legacy illustrative banding; not a universal CMR reference standard."""
    if grs >= 35.0:
        return MyocardialFunctionTier.NORMAL
    if grs >= 25.0:
        return MyocardialFunctionTier.MILDLY_REDUCED
    if grs >= 15.0:
        return MyocardialFunctionTier.MODERATELY_REDUCED
    return MyocardialFunctionTier.SEVERELY_REDUCED


def _overall_tier(*tiers: MyocardialFunctionTier) -> str:
    order = {
        MyocardialFunctionTier.NORMAL: 0,
        MyocardialFunctionTier.MILDLY_REDUCED: 1,
        MyocardialFunctionTier.MODERATELY_REDUCED: 2,
        MyocardialFunctionTier.SEVERELY_REDUCED: 3,
    }
    return max(tiers, key=order.get).value


def evaluate_cine_strain(inp: CineStrainInput) -> CineStrainReport:
    """Derive quantitative metrics from already-measured CMR inputs."""

    _require_positive("LVEDV", inp.lvedv_ml)
    _require_positive("LVESV", inp.lvesv_ml)
    if inp.lvedv_ml <= inp.lvesv_ml:
        raise ValueError("LVEDV must be greater than LVESV.")
    _require_positive("Body surface area", inp.bsa_m2)
    _require_positive("Heart rate", inp.heart_rate_bpm)
    if inp.lv_mass_g < 0:
        raise ValueError("LV mass cannot be negative.")
    if inp.sd_ttp_ms < 0:
        raise ValueError("SD time-to-peak cannot be negative.")
    if not 0.0 <= inp.lge_transmurality_pct <= 100.0:
        raise ValueError("LGE transmurality must be between 0 and 100 percent.")
    if not 0.0 <= inp.lge_scar_mass_pct <= 100.0:
        raise ValueError("LGE scar mass must be between 0 and 100 percent.")

    stroke_vol = inp.lvedv_ml - inp.lvesv_ml
    lvef = stroke_vol / inp.lvedv_ml * 100.0
    cardiac_out = stroke_vol * inp.heart_rate_bpm / 1000.0

    gls_tier = classify_longitudinal_strain(inp.gls_pct)
    gcs_tier = classify_circumferential_strain(inp.gcs_pct)
    grs_tier = classify_radial_strain(inp.grs_pct)

    sre_sra: Optional[float] = None
    if inp.late_diastolic_sr_a_s1 != 0:
        sre_sra = round(inp.early_diastolic_sr_e_s1 / inp.late_diastolic_sr_a_s1, 2)

    ecv_value: Optional[float] = None
    ecv_fields = (
        inp.native_t1_ms,
        inp.post_contrast_t1_myo_ms,
        inp.pre_contrast_t1_blood_ms,
        inp.post_contrast_t1_blood_ms,
        inp.hematocrit_pct,
    )
    if all(v is not None for v in ecv_fields):
        ecv_value = calculate_ecv(
            float(inp.native_t1_ms),
            float(inp.post_contrast_t1_myo_ms),
            float(inp.pre_contrast_t1_blood_ms),
            float(inp.post_contrast_t1_blood_ms),
            float(inp.hematocrit_pct),
        )
    elif any(v is not None for v in ecv_fields):
        raise ValueError("ECV calculation requires all four T1 values and hematocrit.")

    apical_ratio: Optional[float] = None
    if any(v is not None for v in (inp.apical_ls_pct, inp.mid_ls_pct, inp.basal_ls_pct)):
        if not all(v is not None for v in (inp.apical_ls_pct, inp.mid_ls_pct, inp.basal_ls_pct)):
            raise ValueError("Relative apical ratio requires apical, mid, and basal strain values.")
        denominator = (abs(float(inp.basal_ls_pct)) + abs(float(inp.mid_ls_pct))) / 2.0
        if denominator == 0:
            raise ValueError("Mid/basal strain magnitudes cannot both be zero.")
        apical_ratio = abs(float(inp.apical_ls_pct)) / denominator

    lge_pattern = inp.lge_pattern.value if isinstance(inp.lge_pattern, LGEPattern) else str(inp.lge_pattern).strip().lower()
    if lge_pattern in ("", "no", "false"):
        lge_pattern = "none"
    allowed_lge = {item.value for item in LGEPattern}
    if lge_pattern not in allowed_lge:
        raise ValueError(f"Unsupported LGE pattern: {lge_pattern!r}.")

    if lge_pattern == "none":
        lge_summary = "No LGE pattern supplied."
    else:
        lge_summary = (
            f"Reported LGE pattern: {lge_pattern}; transmurality "
            f"{inp.lge_transmurality_pct:.1f}%; scar mass {inp.lge_scar_mass_pct:.1f}% of LV mass."
        )

    findings = [
        "Quantitative summary only; no disease phenotype is inferred.",
        "Illustrative strain bands are retained for compatibility and require method-specific validation.",
    ]
    if apical_ratio is not None:
        findings.append(f"Relative apical strain ratio: {apical_ratio:.2f}; interpret in clinical context.")

    return CineStrainReport(
        study_id=inp.study_id,
        patient_id=inp.patient_id,
        lvef_pct=round(lvef, 1),
        stroke_volume_ml=round(stroke_vol, 1),
        lvedv_index_ml_m2=round(inp.lvedv_ml / inp.bsa_m2, 1),
        lvesv_index_ml_m2=round(inp.lvesv_ml / inp.bsa_m2, 1),
        stroke_volume_index_ml_m2=round(stroke_vol / inp.bsa_m2, 1),
        cardiac_output_l_min=round(cardiac_out, 2),
        cardiac_index_l_min_m2=round(cardiac_out / inp.bsa_m2, 2),
        lv_mass_index_g_m2=round(inp.lv_mass_g / inp.bsa_m2, 1),
        gls_pct=round(inp.gls_pct, 1),
        gcs_pct=round(inp.gcs_pct, 1),
        grs_pct=round(inp.grs_pct, 1),
        gls_function_tier=gls_tier.value,
        gcs_function_tier=gcs_tier.value,
        grs_function_tier=grs_tier.value,
        overall_contractility_tier=_overall_tier(gls_tier, gcs_tier, grs_tier),
        sre_sra_ratio=sre_sra,
        diastolic_function_grade="Not classified from strain-rate values alone",
        mechanical_dyssynchrony="Not classified from SD time-to-peak alone",
        sd_ttp_ms=round(inp.sd_ttp_ms, 1),
        native_t1_ms=round(inp.native_t1_ms, 1) if inp.native_t1_ms is not None else None,
        t1_status="Interpret with a site-, scanner-, sequence-, and method-specific reference interval",
        calculated_ecv_pct=round(ecv_value, 1) if ecv_value is not None else None,
        ecv_status=(
            "Calculated; interpret with an appropriate local/method-specific reference interval"
            if ecv_value is not None
            else "Not calculated"
        ),
        lge_summary=lge_summary,
        apical_sparing_ratio=round(apical_ratio, 2) if apical_ratio is not None else None,
        diagnostic_phenotype=DiagnosticPhenotype.NOT_INFERRED.value,
        clinical_findings=findings,
        remediation_recommendations=[],
    )


def _optional_float(value: Any) -> Optional[float]:
    if value is None:
        return None
    if isinstance(value, str) and not value.strip():
        return None
    return float(value)


def _float_alias(kwargs: Dict[str, Any], names: List[str], default: float) -> float:
    for name in names:
        value = kwargs.get(name)
        if value is not None and not (isinstance(value, str) and not value.strip()):
            return float(value)
    return default


def calculate_metrics(**kwargs: Any) -> Dict[str, Any]:
    """Compatibility wrapper accepting common CSV/CLI aliases."""

    study_id = str(kwargs.get("study_id") or kwargs.get("id") or "CMR-001")
    patient = kwargs.get("patient_id") or kwargs.get("Patient")
    inp = CineStrainInput(
        study_id=study_id,
        patient_id=str(patient) if patient not in (None, "") else None,
        heart_rate_bpm=_float_alias(kwargs, ["heart_rate_bpm", "heart_rate", "hr"], 72.0),
        bsa_m2=_float_alias(kwargs, ["bsa_m2", "bsa"], 1.85),
        lvedv_ml=_float_alias(kwargs, ["lvedv_ml", "edv", "primary_metric"], 145.0),
        lvesv_ml=_float_alias(kwargs, ["lvesv_ml", "esv", "secondary_metric"], 58.0),
        lv_mass_g=_float_alias(kwargs, ["lv_mass_g", "lv_mass", "mass"], 120.0),
        gls_pct=_float_alias(kwargs, ["gls_pct", "gls"], -20.5),
        gcs_pct=_float_alias(kwargs, ["gcs_pct", "gcs"], -22.0),
        grs_pct=_float_alias(kwargs, ["grs_pct", "grs"], 42.0),
        apical_ls_pct=_optional_float(kwargs.get("apical_ls_pct")),
        mid_ls_pct=_optional_float(kwargs.get("mid_ls_pct")),
        basal_ls_pct=_optional_float(kwargs.get("basal_ls_pct")),
        native_t1_ms=_optional_float(kwargs.get("native_t1_ms", kwargs.get("native_t1", 1005.0))),
        post_contrast_t1_myo_ms=_optional_float(kwargs.get("post_contrast_t1_myo_ms", kwargs.get("t1_myo_post", 480.0))),
        pre_contrast_t1_blood_ms=_optional_float(kwargs.get("pre_contrast_t1_blood_ms", kwargs.get("t1_blood_pre", 1580.0))),
        post_contrast_t1_blood_ms=_optional_float(kwargs.get("post_contrast_t1_blood_ms", kwargs.get("t1_blood_post", 320.0))),
        hematocrit_pct=_optional_float(kwargs.get("hematocrit_pct", kwargs.get("hct", 42.0))),
        lge_pattern=str(kwargs.get("lge_pattern") or "none"),
        lge_transmurality_pct=_float_alias(kwargs, ["lge_transmurality_pct", "lge_transmurality"], 0.0),
        lge_scar_mass_pct=_float_alias(kwargs, ["lge_scar_mass_pct", "lge_scar_mass"], 0.0),
    )
    report = evaluate_cine_strain(inp)
    result = report.to_dict()
    result["tool"] = "cardiac-mri-cine-strain-agent"
    result["score"] = report.lvef_pct
    result["classification"] = report.diagnostic_phenotype
    result["clinical_recommendation"] = ""
    return result


def process_batch(input_csv: str, output_csv: str) -> int:
    """Batch-process a CSV while preserving its original columns."""

    with open(input_csv, mode="r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError("Input CSV must contain a header row.")
        rows = list(reader)
        original_fields = list(reader.fieldnames)

    derived_fields = [
        "lvef_pct",
        "stroke_volume_ml",
        "lvedv_index_ml_m2",
        "lvesv_index_ml_m2",
        "stroke_volume_index_ml_m2",
        "cardiac_index_l_min_m2",
        "lv_mass_index_g_m2",
        "calculated_ecv_pct",
        "diagnostic_phenotype",
    ]
    output_fields = list(dict.fromkeys(original_fields + derived_fields))
    output_rows: List[Dict[str, Any]] = []

    for index, row in enumerate(rows, start=2):
        try:
            result = calculate_metrics(**row)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid data on CSV row {index}: {exc}") from exc
        merged = dict(row)
        for field_name in derived_fields:
            merged[field_name] = result[field_name]
        output_rows.append(merged)

    with open(output_csv, mode="w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields)
        writer.writeheader()
        writer.writerows(output_rows)
    return len(output_rows)
