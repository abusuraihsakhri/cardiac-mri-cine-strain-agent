#!/usr/bin/env python3
"""
QuantitativeMappingAgent: MOLLI/ShMOLLI T1 mapping and ECV computation.

MOLLI inversion-recovery signal model fitted per pixel/ROI:

    s(TI) = | A - B * exp(-TI / T1*) |

with the true longitudinal relaxation time recovered as

    T1 = T1* * (B/A - 1)

Fitting is done by a coarse grid search over T1* followed by closed-form
linear least squares for (A, B) at each candidate and golden-section
refinement. Extracellular volume fraction follows the SCMR consensus:

    ECV = dR1_myocardium / dR1_blood * (1 - Hct),   R1 = 1/T1

Diffuse fibrosis threshold: ECV > 27% (upper normal limit ~25-28% at 1.5T).
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Sequence


@dataclass
class MOLLIAcquisition:
    ti_times_ms: List[float]
    signals: List[float]          # magnitude images, arbitrary units


def _fit_linear(x: Sequence[float], y: Sequence[float]):
    n = len(x)
    sx, sy = sum(x), sum(y)
    sxx = sum(v * v for v in x)
    sxy = sum(a * b for a, b in zip(x, y))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-12:
        return 0.0, 0.0
    b = (n * sxy - sx * sy) / denom
    a = (sy - b * sx) / n
    return a, b


def _sse_for_params(t1_star: float, ratio_b_over_a: float, acq: MOLLIAcquisition):
    """Magnitude-model fit: s(TI) = |A - B*exp(-TI/T1*)|.

    For fixed T1* and fixed r = B/A the prediction is linear in the single
    amplitude A, so its optimal value has a closed form; only T1* and r
    need searching. Returns (sse, A, B).
    """
    m = [abs(1.0 - ratio_b_over_a * math.exp(-ti / t1_star))
         for ti in acq.ti_times_ms]
    denom = sum(v * v for v in m)
    if denom < 1e-12:
        return float("inf"), 0.0, 0.0
    A = sum(s * mv for s, mv in zip(acq.signals, m)) / denom
    if A <= 0:
        return float("inf"), A, A * ratio_b_over_a
    B = A * ratio_b_over_a
    sse = sum((s - abs(A - B * math.exp(-ti / t1_star))) ** 2
              for s, ti in zip(acq.signals, acq.ti_times_ms))
    return sse, A, B


def fit_molli_t1(acq: MOLLIAcquisition) -> Dict[str, float]:
    """Coarse 2-D grid (T1* x B/A ratio) + local refinement; magnitude-safe."""
    best = (float("inf"), None, None)
    for t1_star in range(50, 3001, 25):
        r = 1.2
        while r <= 12.0:
            sse, _, _ = _sse_for_params(float(t1_star), r, acq)
            if sse < best[0]:
                best = (sse, float(t1_star), r)
            r += 0.1
    # local refinement around the grid winner
    _, bt, br = best
    step_t, step_r = 10.0, 0.02
    improved = True
    while improved:
        improved = False
        for dt, dr in ((step_t, 0), (-step_t, 0), (0, step_r), (0, -step_r)):
            cand = _sse_for_params(max(20.0, bt + dt), max(1.01, br + dr), acq)[0]
            if cand < best[0]:
                best = (cand, max(20.0, bt + dt), max(1.01, br + dr))
                bt, br = best[1], best[2]
                improved = True

    t1_star, ratio = best[1], best[2]
    xs = [math.exp(-ti / t1_star) for ti in acq.ti_times_ms]
    m = [abs(1.0 - ratio * x) for x in xs]
    denom = sum(v * v for v in m)
    A = sum(s * mv for s, mv in zip(acq.signals, m)) / denom
    B = A * ratio
    true_t1 = t1_star * (ratio - 1.0)
    residuals = [s - abs(A - B * x) for s, x in zip(acq.signals, xs)]
    rmse = math.sqrt(sum(r * r for r in residuals) / len(residuals))

    native_ok = 200.0 < true_t1 < 3500.0
    return {
        "t1_star_ms": round(t1_star, 1),
        "native_T1_ms": round(true_t1, 1),
        "amplitude_A": round(A, 3),
        "amplitude_B": round(B, 3),
        "fit_rmse": round(rmse, 4),
        "physiologic_fit": native_ok,
    }


def compute_ecv(native_t1_myo: float, post_contrast_t1_myo: float,
                native_t1_blood: float, post_contrast_t1_blood: float,
                hematocrit_fraction: float) -> Dict[str, float]:
    """SCMR-consensus extracellular volume fraction."""
    r1_native_myo = 1000.0 / native_t1_myo
    r1_post_myo = 1000.0 / post_contrast_t1_myo
    r1_native_blood = 1000.0 / native_t1_blood
    r1_post_blood = 1000.0 / post_contrast_t1_blood
    delta_r1_myo = r1_post_myo - r1_native_myo
    delta_r1_blood = r1_post_blood - r1_native_blood
    if delta_r1_blood <= 0:
        raise ValueError("post-contrast blood R1 must exceed native R1")
    ecv = (delta_r1_myo / delta_r1_blood) * (1.0 - hematocrit_fraction)
    return {
        "ecv_percent": round(ecv * 100.0, 1),
        "delta_r1_myo": round(delta_r1_myo, 4),
        "delta_r1_blood": round(delta_r1_blood, 4),
        "diffuse_fibrosis_suspected": ecv > 0.27,
        "threshold_reference": "ECV >27% suggests diffuse fibrosis (SCMR)",
    }


if __name__ == "__main__":
    # Synthetic MOLLI curve: myocardium with T1* = 640 ms and B/A chosen
    # so that T1 = T1*(B/A - 1) equals ~980 ms (normal myocardium at 1.5T)
    true_t1star = 640.0
    A_true, B_true = 1200.0, 3040.0
    tis = [100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0]
    sig = [abs(A_true - B_true * math.exp(-t / true_t1star)) for t in tis]
    fit = fit_molli_t1(MOLLIAcquisition(tis, sig))
    print("MOLLI fit (synthetic myocardium)")
    print("-" * 46)
    print(f"true T1*={true_t1star} ms  -> fitted {fit['t1_star_ms']} ms")
    print(f"derived native T1 = {fit['native_T1_ms']} ms "
          f"(rmse {fit['fit_rmse']}, physiologic={fit['physiologic_fit']})")

    ecv = compute_ecv(native_t1_myo=fit["native_T1_ms"],
                      post_contrast_t1_myo=520.0,
                      native_t1_blood=1650.0,
                      post_contrast_t1_blood=420.0,
                      hematocrit_fraction=0.40)
    print("\nECV:", ecv["ecv_percent"], "% | diffuse fibrosis:",
          ecv["diffuse_fibrosis_suspected"])
