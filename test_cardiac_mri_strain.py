#!/usr/bin/env python3
"""
Unit Test Suite for Cardiac MRI Cine Feature-Tracking Strain Agent
==================================================================
Tests cover:
  - Volumetric calculations (LVEF, SV, CO, CI, LV mass index)
  - Global strain categorization (GLS, GCS, GRS)
  - Extracellular Volume Fraction (ECV) formula & boundary conditions
  - Diastolic strain rate ratios (SRE/SRA) & mechanical dyssynchrony
  - Diagnostic cardiomyopathy phenotyping (Amyloidosis, Ischemia, HCM, DCM, Myocarditis)
  - Validation error handling on unphysiological volumes and relaxation times
  - Batch CSV processing & JSON serialization
  - CLI subcommand execution
"""

import csv
import io
import json
import math
import os
import tempfile
import unittest
from unittest.mock import patch

from cardiac_mri_strain import (
    CineStrainInput,
    CineStrainReport,
    DiagnosticPhenotype,
    LGEPattern,
    MyocardialFunctionTier,
    calculate_ecv,
    calculate_metrics,
    classify_circumferential_strain,
    classify_longitudinal_strain,
    classify_radial_strain,
    evaluate_cine_strain,
    process_batch,
)
import cli


class TestCardiacMRICineStrainAgent(unittest.TestCase):

    def test_normal_baseline_volumetrics_and_strain(self):
        """Normal patient CMR: LVEF ~60%, normal GLS, GCS, GRS, synchronous activation."""
        inp = CineStrainInput(
            study_id="TEST-NORM-01",
            lvedv_ml=150.0,
            lvesv_ml=60.0,
            gls_pct=-21.0,
            gcs_pct=-23.0,
            grs_pct=45.0,
            native_t1_ms=1000.0,
        )
        rep = evaluate_cine_strain(inp)
        self.assertAlmostEqual(rep.lvef_pct, 60.0, places=1)
        self.assertEqual(rep.stroke_volume_ml, 90.0)
        self.assertEqual(rep.gls_function_tier, MyocardialFunctionTier.NORMAL.value)
        self.assertEqual(rep.gcs_function_tier, MyocardialFunctionTier.NORMAL.value)
        self.assertEqual(rep.grs_function_tier, MyocardialFunctionTier.NORMAL.value)
        self.assertEqual(rep.diagnostic_phenotype, DiagnosticPhenotype.NORMAL.value)

    def test_gls_classification_tiers(self):
        """Test GLS categorization across mild, moderate, and severe reductions."""
        self.assertEqual(classify_longitudinal_strain(-22.0), MyocardialFunctionTier.NORMAL)
        self.assertEqual(classify_longitudinal_strain(-16.5), MyocardialFunctionTier.MILDLY_REDUCED)
        self.assertEqual(classify_longitudinal_strain(-13.0), MyocardialFunctionTier.MODERATELY_REDUCED)
        self.assertEqual(classify_longitudinal_strain(-8.5), MyocardialFunctionTier.SEVERELY_REDUCED)

    def test_gcs_classification_tiers(self):
        """Test GCS categorization."""
        self.assertEqual(classify_circumferential_strain(-22.0), MyocardialFunctionTier.NORMAL)
        self.assertEqual(classify_circumferential_strain(-17.0), MyocardialFunctionTier.MILDLY_REDUCED)
        self.assertEqual(classify_circumferential_strain(-12.0), MyocardialFunctionTier.MODERATELY_REDUCED)
        self.assertEqual(classify_circumferential_strain(-7.0), MyocardialFunctionTier.SEVERELY_REDUCED)

    def test_grs_classification_tiers(self):
        """Test GRS radial thickening categorization."""
        self.assertEqual(classify_radial_strain(45.0), MyocardialFunctionTier.NORMAL)
        self.assertEqual(classify_radial_strain(30.0), MyocardialFunctionTier.MILDLY_REDUCED)
        self.assertEqual(classify_radial_strain(20.0), MyocardialFunctionTier.MODERATELY_REDUCED)
        self.assertEqual(classify_radial_strain(10.0), MyocardialFunctionTier.SEVERELY_REDUCED)

    def test_ecv_calculation_formula(self):
        """Test formula for ECV calculation from T1 values."""
        # Myo: pre=1000, post=500 -> delta_r1_myo = 1/0.5 - 1/1.0 = 2 - 1 = 1.0
        # Blood: pre=1600, post=350 -> delta_r1_blood = 1/0.35 - 1/1.6 = 2.857 - 0.625 = 2.232
        # Hct = 40% -> (1 - 0.40) * (1.0 / 2.232) * 100 = 0.6 * 0.448 * 100 = 26.88%
        ecv = calculate_ecv(
            t1_myo_pre=1000.0,
            t1_myo_post=500.0,
            t1_blood_pre=1600.0,
            t1_blood_post=350.0,
            hematocrit_pct=40.0,
        )
        self.assertAlmostEqual(ecv, 26.88, places=1)

    def test_amyloidosis_apical_sparing_phenotype(self):
        """Apical sparing strain pattern (ASR > 1.0) and high ECV (>40%) triggers Amyloidosis phenotype."""
        inp = CineStrainInput(
            gls_pct=-11.0,
            apical_ls_pct=-22.0,
            mid_ls_pct=-8.0,
            basal_ls_pct=-6.0,
            native_t1_ms=1150.0,
            post_contrast_t1_myo_ms=380.0,
            post_contrast_t1_blood_ms=350.0,
            hematocrit_pct=38.0,
        )
        rep = evaluate_cine_strain(inp)
        self.assertIsNotNone(rep.apical_sparing_ratio)
        self.assertGreater(rep.apical_sparing_ratio, 1.0)
        self.assertEqual(rep.diagnostic_phenotype, DiagnosticPhenotype.CARDIAC_AMYLOIDOSIS.value)
        self.assertTrue(any("Amyloidosis" in rec for rec in rep.remediation_recommendations))

    def test_ischemic_infarct_phenotype(self):
        """Subendocardial or transmural LGE triggers Ischemic Cardiomyopathy phenotype."""
        inp = CineStrainInput(
            lge_pattern=LGEPattern.TRANSMURAL,
            lge_transmurality_pct=65.0,
            gls_pct=-12.0,
        )
        rep = evaluate_cine_strain(inp)
        self.assertEqual(rep.diagnostic_phenotype, DiagnosticPhenotype.ISCHEMIC_INFARCT.value)
        self.assertTrue(any("Ischemic" in f for f in rep.clinical_findings))

    def test_hypertrophic_cardiomyopathy_phenotype(self):
        """Massive LV hypertrophy (LV mass index > 115) with reduced strain triggers HCM."""
        inp = CineStrainInput(
            bsa_m2=1.8,
            lv_mass_g=230.0,  # Mass index = 127.7 g/m^2
            gls_pct=-13.0,
        )
        rep = evaluate_cine_strain(inp)
        self.assertEqual(rep.diagnostic_phenotype, DiagnosticPhenotype.HYPERTROPHIC_CM.value)

    def test_dilated_cardiomyopathy_phenotype(self):
        """Depressed LVEF < 50% and globally reduced strain with normal mass index triggers DCM."""
        inp = CineStrainInput(
            lvedv_ml=220.0,
            lvesv_ml=140.0,  # LVEF = 36.4%
            gls_pct=-10.0,
            gcs_pct=-12.0,
        )
        rep = evaluate_cine_strain(inp)
        self.assertEqual(rep.diagnostic_phenotype, DiagnosticPhenotype.DILATED_CM.value)
        self.assertTrue(any("GDMT" in rec for rec in rep.remediation_recommendations))

    def test_acute_myocarditis_phenotype(self):
        """Subepicardial LGE with elevated native T1 triggers Acute Myocarditis phenotype."""
        inp = CineStrainInput(
            lge_pattern=LGEPattern.SUBEPICARDIAL,
            native_t1_ms=1120.0,
            lvedv_ml=140.0,
            lvesv_ml=55.0,
        )
        rep = evaluate_cine_strain(inp)
        self.assertEqual(rep.diagnostic_phenotype, DiagnosticPhenotype.ACUTE_MYOCARDITIS.value)

    def test_mechanical_dyssynchrony_alert(self):
        """Time-to-peak SD > 45 ms identifies significant dyssynchrony."""
        inp = CineStrainInput(sd_ttp_ms=52.0)
        rep = evaluate_cine_strain(inp)
        self.assertIn("Dyssynchrony", rep.mechanical_dyssynchrony)

    def test_severe_lvef_warning(self):
        """LVEF < 35% generates ICD eligibility evaluation warning."""
        inp = CineStrainInput(lvedv_ml=250.0, lvesv_ml=180.0)  # LVEF = 28%
        rep = evaluate_cine_strain(inp)
        self.assertTrue(any("ICD eligibility" in f for f in rep.clinical_findings))

    def test_invalid_volumes_exception(self):
        """LVEDV <= LVESV raises ValueError."""
        inp = CineStrainInput(lvedv_ml=100.0, lvesv_ml=120.0)
        with self.assertRaises(ValueError):
            evaluate_cine_strain(inp)

    def test_invalid_bsa_exception(self):
        """BSA out of bounds raises ValueError."""
        inp = CineStrainInput(bsa_m2=5.5)
        with self.assertRaises(ValueError):
            evaluate_cine_strain(inp)

    def test_invalid_t1_ecv_exception(self):
        """Negative T1 relaxation time raises ValueError."""
        with self.assertRaises(ValueError):
            calculate_ecv(t1_myo_pre=-100.0, t1_myo_post=400.0, t1_blood_pre=1500.0, t1_blood_post=300.0, hematocrit_pct=40.0)

    def test_invalid_hematocrit_ecv_exception(self):
        """Hematocrit out of physiological range raises ValueError."""
        with self.assertRaises(ValueError):
            calculate_ecv(t1_myo_pre=1000.0, t1_myo_post=400.0, t1_blood_pre=1500.0, t1_blood_post=300.0, hematocrit_pct=85.0)

    def test_calculate_metrics_wrapper_aliases(self):
        """Test calculate_metrics with varied parameter formats."""
        res = calculate_metrics(
            primary_metric=160.0,   # LVEDV
            secondary_metric=65.0,  # LVESV
            gls=-19.5,
            gcs=-21.0,
            grs=40.0,
        )
        self.assertEqual(res["tool"], "cardiac-mri-cine-strain-agent")
        self.assertAlmostEqual(res["score"], (95.0 / 160.0) * 100.0, places=1)
        self.assertIn("classification", res)

    def test_to_dict_and_json_serialization(self):
        """Ensure report serializes to valid JSON."""
        inp = CineStrainInput(study_id="CMR-JSON-01")
        rep = evaluate_cine_strain(inp)
        d = rep.to_dict()
        s = json.dumps(d)
        deserialized = json.loads(s)
        self.assertEqual(deserialized["study_id"], "CMR-JSON-01")
        self.assertIn("lvef_pct", deserialized)

    def test_batch_processing(self):
        """Test batch CSV processing of CMR studies."""
        with tempfile.TemporaryDirectory() as tmpdir:
            in_csv = os.path.join(tmpdir, "cmr_in.csv")
            out_csv = os.path.join(tmpdir, "cmr_out.csv")

            with open(in_csv, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=["study_id", "lvedv_ml", "lvesv_ml", "gls_pct", "gcs_pct", "grs_pct"])
                writer.writeheader()
                writer.writerow({"study_id": "S1", "lvedv_ml": "140", "lvesv_ml": "50", "gls_pct": "-21.5", "gcs_pct": "-22.0", "grs_pct": "44.0"})
                writer.writerow({"study_id": "S2", "lvedv_ml": "220", "lvesv_ml": "150", "gls_pct": "-9.0", "gcs_pct": "-10.0", "grs_pct": "18.0"})

            count = process_batch(in_csv, out_csv)
            self.assertEqual(count, 2)
            self.assertTrue(os.path.exists(out_csv))

            with open(out_csv, "r", encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
                self.assertEqual(len(rows), 2)
                self.assertEqual(rows[0]["study_id"], "S1")
                self.assertEqual(rows[1]["study_id"], "S2")

    def test_cli_evaluate_command(self):
        """Test CLI evaluation subcommand."""
        with patch("sys.stdout", new_callable=io.StringIO) as mock_out:
            exit_code = cli.main(["evaluate", "--study-id", "CLI-CMR-01", "--edv", "150", "--esv", "60", "--gls", "-20.0"])
            self.assertEqual(exit_code, 0)
            output = mock_out.getvalue()
            self.assertIn("CARDIAC MRI CINE STRAIN", output)
            self.assertIn("CLI-CMR-01", output)

    def test_cli_evaluate_json_flag(self):
        """Test CLI JSON mode."""
        with patch("sys.stdout", new_callable=io.StringIO) as mock_out:
            exit_code = cli.main(["evaluate", "--study-id", "CLI-JSON", "--json"])
            self.assertEqual(exit_code, 0)
            data = json.loads(mock_out.getvalue())
            self.assertEqual(data["study_id"], "CLI-JSON")
            self.assertIn("lvef_pct", data)

    def test_cli_norms_command(self):
        """Test CLI norms reference table."""
        with patch("sys.stdout", new_callable=io.StringIO) as mock_out:
            exit_code = cli.main(["norms"])
            self.assertEqual(exit_code, 0)
            output = mock_out.getvalue()
            self.assertIn("Global Longitudinal Strain", output)
            self.assertIn("Extracellular Volume", output)

    def test_interactive_wizard_mock(self):
        """Test interactive wizard prompt."""
        user_inputs = ["CMR-WIZ", "PAT-01", "70", "1.9", "150", "60", "120", "-21.0", "-23.0", "45.0", "1000", "1"]
        with patch("builtins.input", side_effect=user_inputs):
            inp = cli.interactive_wizard()
            self.assertEqual(inp.study_id, "CMR-WIZ")
            self.assertEqual(inp.patient_id, "PAT-01")
            self.assertEqual(inp.lvedv_ml, 150.0)

    def test_low_native_t1_iron_overload(self):
        """Low Native T1 (< 930 ms) triggers iron overload / Fabry disease alert."""
        inp = CineStrainInput(native_t1_ms=880.0)
        rep = evaluate_cine_strain(inp)
        self.assertIn("Iron Overload", rep.t1_status)

    def test_cardiac_output_and_index_calculation(self):
        """Test Stroke Volume, Cardiac Output, and Cardiac Index."""
        inp = CineStrainInput(
            lvedv_ml=160.0,
            lvesv_ml=60.0,   # SV = 100 mL
            heart_rate_bpm=70.0,  # CO = 100 * 70 / 1000 = 7.0 L/min
            bsa_m2=2.0,      # CI = 7.0 / 2.0 = 3.5 L/min/m^2
        )
        rep = evaluate_cine_strain(inp)
        self.assertEqual(rep.stroke_volume_ml, 100.0)
        self.assertAlmostEqual(rep.cardiac_output_l_min, 7.0, places=2)
        self.assertAlmostEqual(rep.cardiac_index_l_min_m2, 3.5, places=2)

    def test_grade_1_diastolic_dysfunction(self):
        """Early diastolic strain rate 0.85-1.20 s^-1 graded as Grade 1 impaired relaxation."""
        inp = CineStrainInput(early_diastolic_sr_e_s1=0.95)
        rep = evaluate_cine_strain(inp)
        self.assertEqual(rep.diastolic_function_grade, "Grade 1 (Impaired Relaxation)")

    def test_grade_2_diastolic_dysfunction(self):
        """Early diastolic strain rate < 0.85 s^-1 graded as Grade 2/3."""
        inp = CineStrainInput(early_diastolic_sr_e_s1=0.60)
        rep = evaluate_cine_strain(inp)
        self.assertEqual(rep.diastolic_function_grade, "Grade 2/3 (Advanced Diastolic Dysfunction)")

    def test_calculate_metrics_default_invocation(self):
        """Test default invocation returns complete dictionary."""
        res = calculate_metrics()
        self.assertEqual(res["tool"], "cardiac-mri-cine-strain-agent")
        self.assertIn("lvef_pct", res)


if __name__ == "__main__":
    unittest.main()
