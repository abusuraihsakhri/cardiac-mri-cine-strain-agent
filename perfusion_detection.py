#!/usr/bin/env python3
"""
Perfusion Defect Detection for Cardiac MRI Cine Strain Agent.
Detects myocardial perfusion defects from first-pass perfusion imaging.
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass


@dataclass
class PerfusionSegment:
    """Myocardial perfusion segment."""
    segment_name: str  # e.g., "LAD_1", "LCx_4", "RCA_7"
    territory: str
    signal_intensity: float
    time_to_peak: float
    upslope: float
    defect: bool = False


def detect_perfusion_defects(segments: List[PerfusionSegment]) -> Dict[str, Any]:
    """Detect myocardial perfusion defects from segmental data."""
    if not segments:
        return {"error": "No perfusion segments provided"}

    mean_si = sum(s.signal_intensity for s in segments) / len(segments)
    mean_ttp = sum(s.time_to_peak for s in segments) / len(segments)
    mean_upslope = sum(s.upslope for s in segments) / len(segments)

    defects = []
    for seg in segments:
        si_ratio = seg.signal_intensity / max(mean_si, 0.01)
        ttp_delay = seg.time_to_peak - mean_ttp
        upslope_ratio = seg.upslope / max(mean_upslope, 0.01)

        defect_score = 0
        if si_ratio < 0.7:
            defect_score += 2
        elif si_ratio < 0.85:
            defect_score += 1

        if ttp_delay > 2.0:
            defect_score += 2
        elif ttp_delay > 1.0:
            defect_score += 1

        if upslope_ratio < 0.5:
            defect_score += 2
        elif upslope_ratio < 0.7:
            defect_score += 1

        seg.defect = defect_score >= 2
        if seg.defect:
            defects.append({
                "segment": seg.segment_name,
                "territory": seg.territory,
                "si_ratio": round(si_ratio, 2),
                "ttp_delay_sec": round(ttp_delay, 1),
                "upslope_ratio": round(upslope_ratio, 2),
                "severity": "severe" if defect_score >= 4 else "moderate",
            })

    total_segments = len(segments)
    defect_count = len(defects)
    defect_pct = (defect_count / total_segments * 100) if total_segments > 0 else 0

    affected_territories = list(set(d["territory"] for d in defects))

    if defect_pct >= 30:
        perfusion_status = "severely_reduced"
        recommendation = "High likelihood of flow-limiting stenosis. Coronary angiography recommended."
    elif defect_pct >= 15:
        perfusion_status = "moderately_reduced"
        recommendation = "Moderate perfusion defects. Consider stress testing and angiography."
    elif defect_pct > 0:
        perfusion_status = "mildly_reduced"
        recommendation = "Focal defects detected. Clinical correlation recommended."
    else:
        perfusion_status = "normal"
        recommendation = "No significant perfusion defects detected."

    return {
        "total_segments": total_segments,
        "defect_count": defect_count,
        "defect_percentage": round(defect_pct, 1),
        "defects": defects,
        "affected_territories": affected_territories,
        "perfusion_status": perfusion_status,
        "recommendation": recommendation,
        "mean_si": round(mean_si, 2),
    }


class PerfusionAgent:
    """Sub-agent for perfusion defect detection."""

    def __init__(self):
        self.agent_name = "PerfusionAgent"

    def evaluate(self, segments: List[PerfusionSegment]) -> Dict[str, Any]:
        """Evaluate perfusion segments."""
        result = detect_perfusion_defects(segments)
        alerts = []

        if "error" in result:
            return {"perfusion_result": result, "alerts": [{"type": "NO_DATA", "severity": "ERROR",
                    "message": result["error"], "recommendation": "Provide perfusion data."}]}

        if result["perfusion_status"] in ("severely_reduced", "moderately_reduced"):
            alerts.append({
                "type": "PERFUSION_DEFECT", "severity": "WARNING",
                "message": f"{result['defect_count']} perfusion defect(s) in "
                           f"{len(result['affected_territories'])} territory(ies).",
                "recommendation": result["recommendation"]
            })

        return {"perfusion_result": result, "alerts": alerts}
