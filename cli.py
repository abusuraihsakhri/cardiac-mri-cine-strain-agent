#!/usr/bin/env python3
"""
Cardiac MRI Cine Feature-Tracking Strain CLI
============================================
Command line interface for assessing myocardial deformation, global strain (GLS/GCS/GRS),
T1 mapping, ECV fraction, LGE transmurality, and cardiomyopathy phenotyping.
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Optional

from cardiac_mri_strain import (
    CineStrainInput,
    LGEPattern,
    calculate_metrics,
    evaluate_cine_strain,
    process_batch,
)


def format_report_table(report: dict) -> str:
    """Format single CMR study into an ASCII report table."""
    lines = []
    lines.append("=" * 76)
    lines.append(f"  CARDIAC MRI CINE STRAIN & TISSUE CHARACTERIZATION DOSSIER")
    lines.append("=" * 76)
    lines.append(f"  Study ID              : {report['study_id']}")
    lines.append(f"  Patient ID            : {report.get('patient_id') or 'N/A'}")
    lines.append(f"  Diagnostic Phenotype  : {report['diagnostic_phenotype']}")
    lines.append(f"  Contractility Status  : {report['overall_contractility_tier']}")
    lines.append("-" * 76)
    lines.append("  LEFT VENTRICULAR VOLUMETRICS & HEMODYNAMICS:")
    lines.append(f"    * Ejection Fraction (LVEF)     : {report['lvef_pct']:.1f}%")
    lines.append(f"    * Stroke Volume (LVSV)         : {report['stroke_volume_ml']:.1f} mL")
    lines.append(f"    * Cardiac Output               : {report['cardiac_output_l_min']:.2f} L/min")
    lines.append(f"    * Cardiac Index                : {report['cardiac_index_l_min_m2']:.2f} L/min/m^2")
    lines.append(f"    * LV Mass Index                : {report['lv_mass_index_g_m2']:.1f} g/m^2")
    lines.append("-" * 76)
    lines.append("  MYOCARDIAL DEFORMATION & STRAIN (CINE FEATURE TRACKING):")
    lines.append(f"    * Global Longitudinal (GLS)   : {report['gls_pct']:+.1f}% [{report['gls_function_tier']}]")
    lines.append(f"    * Global Circumferential (GCS): {report['gcs_pct']:+.1f}% [{report['gcs_function_tier']}]")
    lines.append(f"    * Global Radial (GRS)          : {report['grs_pct']:+.1f}% [{report['grs_function_tier']}]")
    lines.append(f"    * SRE / SRA Strain Rate Ratio  : {report['sre_sra_ratio']:.2f} ({report['diastolic_function_grade']})")
    lines.append(f"    * Mechanical Dyssynchrony (SD) : {report['sd_ttp_ms']:.1f} ms ({report['mechanical_dyssynchrony']})")
    if report.get("apical_sparing_ratio") is not None:
        lines.append(f"    * Apical Sparing Ratio (ASR)   : {report['apical_sparing_ratio']:.2f} (Threshold > 1.0)")
    lines.append("-" * 76)
    lines.append("  PARAMETRIC TISSUE CHARACTERIZATION (T1 / ECV / LGE):")
    lines.append(f"    * Native T1 Relaxation Time    : {report['native_t1_ms']:.0f} ms ({report['t1_status']})")
    if report.get("calculated_ecv_pct") is not None:
        lines.append(f"    * Extracellular Volume (ECV)   : {report['calculated_ecv_pct']:.1f}% ({report['ecv_status']})")
    lines.append(f"    * LGE Scar Assessment          : {report['lge_summary']}")
    lines.append("-" * 76)
    lines.append("  KEY CLINICAL FINDINGS:")
    for f in report.get("clinical_findings", []):
        lines.append(f"    - {f}")
    lines.append("-" * 76)
    lines.append("  RECOMMENDATIONS & ACTION PLAN:")
    for r in report.get("remediation_recommendations", []):
        lines.append(f"    - {r}")
    lines.append("=" * 76)
    return "\n".join(lines)


def interactive_wizard() -> CineStrainInput:
    """Run interactive question prompt to gather CMR exam data."""
    print("\n--- Cardiac MRI Cine Strain Evaluation Wizard ---")
    study_id = input("Study ID [CMR-STUDY-01]: ").strip() or "CMR-STUDY-01"
    patient_id = input("Patient ID / MRN: ").strip() or None

    def ask_float(prompt: str, default: float) -> float:
        resp = input(f"{prompt} [{default}]: ").strip()
        if not resp:
            return default
        try:
            return float(resp)
        except ValueError:
            print(f"Invalid float, defaulting to {default}")
            return default

    hr = ask_float("Heart Rate (bpm)", 72.0)
    bsa = ask_float("Body Surface Area BSA (m^2)", 1.85)
    edv = ask_float("LV End-Diastolic Volume LVEDV (mL)", 145.0)
    esv = ask_float("LV End-Systolic Volume LVESV (mL)", 58.0)
    mass = ask_float("LV Mass (g)", 120.0)

    gls = ask_float("Global Longitudinal Strain GLS (%) [e.g. -20.5]", -20.5)
    gcs = ask_float("Global Circumferential Strain GCS (%) [e.g. -22.0]", -22.0)
    grs = ask_float("Global Radial Strain GRS (%) [e.g. 42.0]", 42.0)

    t1 = ask_float("Native T1 relaxation time (ms) [e.g. 1005]", 1005.0)
    
    print("\nLate Gadolinium Enhancement (LGE) Pattern:")
    print("  [1] None")
    print("  [2] Subendocardial (Ischemic)")
    print("  [3] Transmural (Infarct)")
    print("  [4] Mid-wall (Non-ischemic / DCM)")
    print("  [5] Subepicardial (Myocarditis)")
    print("  [6] Diffuse (Amyloidosis)")
    lge_choice = input("Select pattern (1-6) [1]: ").strip()
    lge_map = {"1": "none", "2": "subendocardial", "3": "transmural", "4": "mid_wall", "5": "subepicardial", "6": "diffuse"}
    lge_pat = lge_map.get(lge_choice, "none")

    transmurality = 0.0
    if lge_pat != "none":
        transmurality = ask_float("LGE Transmurality (%)", 35.0)

    return CineStrainInput(
        study_id=study_id,
        patient_id=patient_id,
        heart_rate_bpm=hr,
        bsa_m2=bsa,
        lvedv_ml=edv,
        lvesv_ml=esv,
        lv_mass_g=mass,
        gls_pct=gls,
        gcs_pct=gcs,
        grs_pct=grs,
        native_t1_ms=t1,
        lge_pattern=lge_pat,
        lge_transmurality_pct=transmurality,
    )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="cardio-mri-strain",
        description="Cardiac MRI Cine Feature-Tracking Strain & Tissue Characterization",
    )
    subparsers = parser.add_subparsers(dest="command", help="Available subcommands")

    # Evaluate command
    eval_parser = subparsers.add_parser("evaluate", help="Evaluate single CMR study")
    eval_parser.add_argument("--study-id", default="CMR-EXP-01", help="Study exam ID")
    eval_parser.add_argument("--patient-id", default=None, help="Patient identifier")
    eval_parser.add_argument("--hr", type=float, default=72.0, help="Heart rate in bpm")
    eval_parser.add_argument("--bsa", type=float, default=1.85, help="Body surface area in m^2")
    eval_parser.add_argument("--edv", type=float, default=145.0, help="LV End-Diastolic Volume (mL)")
    eval_parser.add_argument("--esv", type=float, default=58.0, help="LV End-Systolic Volume (mL)")
    eval_parser.add_argument("--mass", type=float, default=120.0, help="LV Mass (g)")
    eval_parser.add_argument("--gls", type=float, default=-20.5, help="Global Longitudinal Strain (%%, negative)")
    eval_parser.add_argument("--gcs", type=float, default=-22.0, help="Global Circumferential Strain (%%, negative)")
    eval_parser.add_argument("--grs", type=float, default=42.0, help="Global Radial Strain (%%, positive)")
    eval_parser.add_argument("--apical-ls", type=float, default=None, help="Apical Longitudinal Strain (%%)")
    eval_parser.add_argument("--mid-ls", type=float, default=None, help="Mid-cavity Longitudinal Strain (%%)")
    eval_parser.add_argument("--basal-ls", type=float, default=None, help="Basal Longitudinal Strain (%%)")
    eval_parser.add_argument("--native-t1", type=float, default=1005.0, help="Native T1 (ms)")
    eval_parser.add_argument("--post-t1-myo", type=float, default=480.0, help="Post-contrast myocardial T1 (ms)")
    eval_parser.add_argument("--pre-t1-blood", type=float, default=1580.0, help="Pre-contrast blood T1 (ms)")
    eval_parser.add_argument("--post-t1-blood", type=float, default=320.0, help="Post-contrast blood T1 (ms)")
    eval_parser.add_argument("--hct", type=float, default=42.0, help="Hematocrit (%%)")
    eval_parser.add_argument("--lge-pattern", choices=["none", "subendocardial", "transmural", "mid_wall", "subepicardial", "diffuse"], default="none", help="LGE scar pattern")
    eval_parser.add_argument("--lge-transmurality", type=float, default=0.0, help="LGE transmurality (%%)")
    eval_parser.add_argument("--json", action="store_true", help="Output results in JSON format")

    # Interactive command
    interactive_parser = subparsers.add_parser("interactive", help="Interactive CMR clinical wizard")
    interactive_parser.add_argument("--json", action="store_true", help="Output results in JSON format")

    # Batch command
    batch_parser = subparsers.add_parser("batch", help="Batch process CMR examinations from CSV")
    batch_parser.add_argument("-i", "--input", required=True, help="Input CSV path")
    batch_parser.add_argument("-o", "--output", default="cmr_batch_results.csv", help="Output CSV path")

    # Reference norm table
    norms_parser = subparsers.add_parser("norms", help="Display CMR strain and tissue mapping normal ranges")

    args = parser.parse_args(argv)

    if args.command == "evaluate":
        inp = CineStrainInput(
            study_id=args.study_id,
            patient_id=args.patient_id,
            heart_rate_bpm=args.hr,
            bsa_m2=args.bsa,
            lvedv_ml=args.edv,
            lvesv_ml=args.esv,
            lv_mass_g=args.mass,
            gls_pct=args.gls,
            gcs_pct=args.gcs,
            grs_pct=args.grs,
            apical_ls_pct=args.apical_ls,
            mid_ls_pct=args.mid_ls,
            basal_ls_pct=args.basal_ls,
            native_t1_ms=args.native_t1,
            post_contrast_t1_myo_ms=args.post_t1_myo,
            pre_contrast_t1_blood_ms=args.pre_t1_blood,
            post_contrast_t1_blood_ms=args.post_t1_blood,
            hematocrit_pct=args.hct,
            lge_pattern=args.lge_pattern,
            lge_transmurality_pct=args.lge_transmurality,
        )
        report = evaluate_cine_strain(inp)
        if args.json:
            print(json.dumps(report.to_dict(), indent=2))
        else:
            print(format_report_table(report.to_dict()))
        return 0

    elif args.command == "interactive":
        inp = interactive_wizard()
        report = evaluate_cine_strain(inp)
        if args.json:
            print(json.dumps(report.to_dict(), indent=2))
        else:
            print(format_report_table(report.to_dict()))
        return 0

    elif args.command == "batch":
        count = process_batch(args.input, args.output)
        print(f"Processed {count} CMR exams into '{args.output}'.")
        return 0

    elif args.command == "norms":
        print("=" * 72)
        print("  CARDIAC MRI FEATURE TRACKING & T1/ECV REFERENCE NORMAL RANGES")
        print("=" * 72)
        print("  Global Longitudinal Strain (GLS)     : -18.0% to -25.0% (More negative = better)")
        print("  Global Circumferential Strain (GCS)   : -19.0% to -26.0%")
        print("  Global Radial Strain (GRS)            : +35.0% to +55.0%")
        print("  Left Ventricular Ejection Fraction    : 55.0% to 70.0%")
        print("  Native T1 (1.5 Tesla)                 : 950 ms to 1060 ms")
        print("  Native T1 (3.0 Tesla)                 : 1150 ms to 1280 ms")
        print("  Extracellular Volume Fraction (ECV)   : 23.0% to 28.0% (Amyloidosis >= 40%)")
        print("  Mechanical Dyssynchrony (SD-TTP)      : < 35 ms")
        print("=" * 72)
        return 0

    else:
        parser.print_help()
        return 0


if __name__ == "__main__":
    sys.exit(main())
