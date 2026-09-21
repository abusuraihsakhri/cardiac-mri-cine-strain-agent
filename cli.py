#!/usr/bin/env python3
"""Command-line interface for quantitative cardiac MRI calculations."""
from __future__ import annotations

import argparse
import json
import sys

from cardiac_mri_strain import REFERENCE_NOTE, CineStrainInput, evaluate_cine_strain, process_batch


def format_report_table(report: dict) -> str:
    lines = [
        "=" * 72,
        "CARDIAC MRI QUANTITATIVE SUMMARY",
        "=" * 72,
        f"Study: {report['study_id']}",
        f"LVEF: {report['lvef_pct']:.1f}% | SV: {report['stroke_volume_ml']:.1f} mL | CI: {report['cardiac_index_l_min_m2']:.2f} L/min/m²",
        f"EDVI: {report['lvedv_index_ml_m2']:.1f} | ESVI: {report['lvesv_index_ml_m2']:.1f} | LVMI: {report['lv_mass_index_g_m2']:.1f}",
        f"GLS: {report['gls_pct']:+.1f}% | GCS: {report['gcs_pct']:+.1f}% | GRS: {report['grs_pct']:+.1f}%",
    ]
    if report.get("calculated_ecv_pct") is not None:
        lines.append(f"ECV: {report['calculated_ecv_pct']:.1f}%")
    lines.extend([report["lge_summary"], "-" * 72, REFERENCE_NOTE, "=" * 72])
    return "\n".join(lines)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cardio-mri-strain",
        description="Derive quantitative CMR metrics from user-supplied measurements.",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    ev = sub.add_parser("evaluate", help="Evaluate one set of measurements")
    ev.add_argument("--study-id", default="CMR-001")
    ev.add_argument("--hr", type=float, default=72.0)
    ev.add_argument("--bsa", type=float, default=1.85)
    ev.add_argument("--edv", type=float, default=145.0)
    ev.add_argument("--esv", type=float, default=58.0)
    ev.add_argument("--mass", type=float, default=120.0)
    ev.add_argument("--gls", type=float, default=-20.5)
    ev.add_argument("--gcs", type=float, default=-22.0)
    ev.add_argument("--grs", type=float, default=42.0)
    ev.add_argument("--native-t1", type=float, default=1005.0)
    ev.add_argument("--post-t1-myo", type=float, default=480.0)
    ev.add_argument("--pre-t1-blood", type=float, default=1580.0)
    ev.add_argument("--post-t1-blood", type=float, default=320.0)
    ev.add_argument("--hct", type=float, default=42.0)
    ev.add_argument("--lge-pattern", default="none", choices=["none", "subendocardial", "transmural", "mid_wall", "subepicardial", "diffuse"])
    ev.add_argument("--lge-transmurality", type=float, default=0.0)
    ev.add_argument("--json", action="store_true")
    batch = sub.add_parser("batch", help="Process a CSV file")
    batch.add_argument("-i", "--input", required=True)
    batch.add_argument("-o", "--output", default="results.csv")
    sub.add_parser("norms", help="Explain reference-range handling")
    return parser


def main(argv=None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "evaluate":
            report = evaluate_cine_strain(CineStrainInput(
                study_id=args.study_id,
                heart_rate_bpm=args.hr,
                bsa_m2=args.bsa,
                lvedv_ml=args.edv,
                lvesv_ml=args.esv,
                lv_mass_g=args.mass,
                gls_pct=args.gls,
                gcs_pct=args.gcs,
                grs_pct=args.grs,
                native_t1_ms=args.native_t1,
                post_contrast_t1_myo_ms=args.post_t1_myo,
                pre_contrast_t1_blood_ms=args.pre_t1_blood,
                post_contrast_t1_blood_ms=args.post_t1_blood,
                hematocrit_pct=args.hct,
                lge_pattern=args.lge_pattern,
                lge_transmurality_pct=args.lge_transmurality,
            ))
            print(json.dumps(report.to_dict(), indent=2) if args.json else format_report_table(report.to_dict()))
            return 0
        if args.command == "batch":
            count = process_batch(args.input, args.output)
            print(f"Processed {count} row(s) -> {args.output}")
            return 0
        if args.command == "norms":
            print(REFERENCE_NOTE)
            return 0
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
