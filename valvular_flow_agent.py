#!/usr/bin/env python3
"""
ValvularFlowAgent: phase-contrast MRI regurgitation quantification.

Through-plane velocity-encoded cine data are integrated over the cardiac
cycle to yield forward and reverse stroke volumes:

    V_phase [mL] = velocity [m/s] * area [cm2] * dt [ms] / 10

    RegurgitantVolume = integral of reverse flow
    RegurgitantFraction = RegurgitantVolume / ForwardVolume * 100

CMR severity bands for chronic regurgitation (regurgitant fraction):
    trivial <15% | mild 16-25% | moderate 26-48% | severe >48%
Severe aortic regurgitation additionally implies regurgitant volume >60 mL.

Also cross-checks forward-flow stroke volume against the cine-derived LV
stroke volume and flags velocity-encoding (VENC) aliasing when peak velocity
approaches the encoding limit.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional


@dataclass
class FlowSeries:
    phase_intervals_ms: List[float]     # dt between phases (same length as velocities)
    net_velocities_m_s: List[float]     # signed through-plane mean velocities
    roi_area_cm2: float
    venc_cm_s: float = 150.0
    peak_velocity_m_s: float = 1.2
    cine_lv_stroke_volume_ml: Optional[float] = None


SEVERITY_BANDS = [
    ("trivial", 0.0, 15.0),
    ("mild", 15.0, 25.0),
    ("moderate", 25.0, 48.0),
    ("severe", 48.0, 1e9),
]


def _band_for_rf(rf_pct: float) -> str:
    for name, lo, hi in SEVERITY_BANDS:
        if lo <= rf_pct < hi:
            return name
    return "severe"


def quantify_regurgitation(f: FlowSeries, valve: str = "aortic") -> Dict[str, Any]:
    if len(f.phase_intervals_ms) != len(f.net_velocities_m_s):
        raise ValueError("phase intervals and velocities length mismatch")

    volumes_ml = [v * f.roi_area_cm2 * dt / 10.0
                  for v, dt in zip(f.net_velocities_m_s, f.phase_intervals_ms)]
    forward = sum(v for v in volumes_ml if v > 0)
    reverse = -sum(v for v in volumes_ml if v < 0)

    if forward <= 0:
        raise ValueError("no net forward flow detected")
    rv = reverse
    rf_pct = 100.0 * reverse / forward
    band = _band_for_rf(rf_pct)

    checks: List[str] = []
    if valve.lower().startswith("aort") and band == "severe" and rv <= 60.0:
        checks.append("RF severe but RVol <= 60 mL: reconcile before calling severe")
    if f.peak_velocity_m_s * 100.0 > 0.8 * f.venc_cm_s:
        checks.append("peak velocity > 80% of VENC: possible aliasing; repeat "
                      "with higher VENC or correct for wrapped phase")
    if f.cine_lv_stroke_volume_ml is not None:
        mismatch_pct = 100.0 * abs(forward - f.cine_lv_stroke_volume_ml) / \
                       f.cine_lv_stroke_volume_ml
        if mismatch_pct > 10.0:
            checks.append(f"forward-flow SV differs from cine SV by "
                          f"{mismatch_pct:.1f}%: verify ROI placement")
        else:
            checks.append(f"forward-flow SV agrees with cine SV "
                          f"(difference {mismatch_pct:.1f}%)")

    return {
        "valve": valve,
        "forward_volume_ml": round(forward, 1),
        "regurgitant_volume_ml": round(rv, 1),
        "regurgitant_fraction_pct": round(rf_pct, 1),
        "severity": band,
        "quality_checks": checks,
    }


if __name__ == "__main__":
    # Synthetic moderate AR waveform through the aortic root ROI
    rng_dt = [40.0] * 25
    velocities = [
        1.05, 1.18, 1.12, 0.98, 0.78, 0.52, 0.25, 0.06,   # systolic ejection
        0.0, 0.0,
        -0.22, -0.24, -0.22, -0.19, -0.16,                # diastolic regurgitation
        -0.13, -0.11, -0.09, -0.07, -0.05, -0.04, -0.03, -0.02, 0.0, 0.0,
    ]
    series = FlowSeries(rng_dt, velocities, roi_area_cm2=3.5, venc_cm_s=150,
                        peak_velocity_m_s=1.18, cine_lv_stroke_volume_ml=83.0)
    rep = quantify_regurgitation(series)
    print("Phase-contrast aortic valve analysis")
    print("-" * 50)
    print(f"forward SV      : {rep['forward_volume_ml']} mL")
    print(f"regurgitant vol : {rep['regurgitant_volume_ml']} mL")
    print(f"RF              : {rep['regurgitant_fraction_pct']}% -> {rep['severity'].upper()}")
    for c in rep["quality_checks"]:
        print(f"  - {c}")
