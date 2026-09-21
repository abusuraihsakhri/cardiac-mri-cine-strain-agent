import csv
import json
import subprocess
import sys
from pathlib import Path

import pytest

from cardiac_mri_strain import (
    CineStrainInput,
    DiagnosticPhenotype,
    calculate_ecv,
    calculate_metrics,
    evaluate_cine_strain,
    process_batch,
)


def test_volumetric_and_index_calculations():
    report = evaluate_cine_strain(CineStrainInput(lvedv_ml=160, lvesv_ml=60, heart_rate_bpm=70, bsa_m2=2.0, lv_mass_g=120))
    assert report.stroke_volume_ml == 100.0
    assert report.lvef_pct == 62.5
    assert report.lvedv_index_ml_m2 == 80.0
    assert report.lvesv_index_ml_m2 == 30.0
    assert report.stroke_volume_index_ml_m2 == 50.0
    assert report.cardiac_output_l_min == 7.0
    assert report.cardiac_index_l_min_m2 == 3.5
    assert report.lv_mass_index_g_m2 == 60.0


def test_ecv_formula():
    ecv = calculate_ecv(1000, 500, 1600, 350, 40)
    assert ecv == pytest.approx(26.88, abs=0.05)


def test_no_disease_phenotype_or_treatment_is_inferred():
    report = evaluate_cine_strain(CineStrainInput(lge_pattern="transmural", lge_transmurality_pct=80, gls_pct=-8))
    assert report.diagnostic_phenotype == DiagnosticPhenotype.NOT_INFERRED.value
    assert report.remediation_recommendations == []
    assert "no disease phenotype" in report.clinical_findings[0].lower()


def test_partial_ecv_inputs_rejected():
    with pytest.raises(ValueError, match="requires all four T1 values"):
        evaluate_cine_strain(CineStrainInput(native_t1_ms=1000, post_contrast_t1_myo_ms=None, pre_contrast_t1_blood_ms=1500, post_contrast_t1_blood_ms=350, hematocrit_pct=40))


def test_invalid_volumes_rejected():
    with pytest.raises(ValueError, match="LVEDV must be greater"):
        evaluate_cine_strain(CineStrainInput(lvedv_ml=80, lvesv_ml=90))


def test_compatibility_wrapper():
    result = calculate_metrics(primary_metric=160, secondary_metric=60, heart_rate=70, bsa=2)
    assert result["tool"] == "cardiac-mri-cine-strain-agent"
    assert result["score"] == 62.5
    assert result["classification"] == DiagnosticPhenotype.NOT_INFERRED.value


def test_batch_round_trip(tmp_path: Path):
    source = tmp_path / "input.csv"
    target = tmp_path / "output.csv"
    source.write_text("study_id,lvedv_ml,lvesv_ml,bsa_m2,heart_rate_bpm\nS1,150,60,2,70\n", encoding="utf-8")
    assert process_batch(str(source), str(target)) == 1
    rows = list(csv.DictReader(target.open(encoding="utf-8")))
    assert rows[0]["study_id"] == "S1"
    assert rows[0]["lvef_pct"] == "60.0"
    assert rows[0]["diagnostic_phenotype"] == DiagnosticPhenotype.NOT_INFERRED.value


def test_cli_json_smoke():
    proc = subprocess.run([sys.executable, "cli.py", "evaluate", "--edv", "150", "--esv", "60", "--json"], check=True, capture_output=True, text=True)
    payload = json.loads(proc.stdout)
    assert payload["lvef_pct"] == 60.0
