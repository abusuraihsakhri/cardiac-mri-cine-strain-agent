#!/usr/bin/env python3
"""
Tissue Tracking Agent for Cardiac MRI Cine Strain Agent.
Tracks myocardial tissue deformation and strain patterns across cardiac phases.
"""

from typing import Dict, Any, List, Optional
from dataclasses import dataclass


@dataclass
class TissuePoint:
    """Single tissue tracking point."""
    phase: int
    x: float
    y: float
    z: float
    displacement: float = 0.0
    strain: float = 0.0


def calculate_tissue_strain(points: List[TissuePoint], reference_phase: int = 0) -> Dict[str, Any]:
    """Calculate tissue strain from tracking points."""
    reference = [p for p in points if p.phase == reference_phase]
    if not reference:
        return {"error": "Reference phase not found"}

    ref_center = (
        sum(p.x for p in reference) / len(reference),
        sum(p.y for p in reference) / len(reference),
    )

    phase_strains = {}
    for phase in set(p.phase for p in points if p.phase != reference_phase):
        phase_points = [p for p in points if p.phase == phase]
        if phase_points:
            curr_center = (
                sum(p.x for p in phase_points) / len(phase_points),
                sum(p.y for p in phase_points) / len(phase_points),
            )
            dx = curr_center[0] - ref_center[0]
            dy = curr_center[1] - ref_center[1]
            displacement = (dx**2 + dy**2) ** 0.5
            ref_dist = (ref_center[0]**2 + ref_center[1]**2) ** 0.5
            strain = displacement / max(ref_dist, 0.01)
            phase_strains[phase] = {
                "displacement": round(displacement, 3),
                "strain": round(strain, 4),
            }

    strains = [v["strain"] for v in phase_strains.values()]
    peak_strain = max(strains) if strains else 0.0
    peak_phase = max(phase_strains, key=lambda k: phase_strains[k]["strain"]) if phase_strains else 0

    if peak_strain >= 0.20:
        function_category = "normal"
    elif peak_strain >= 0.14:
        function_category = "mildly_reduced"
    elif peak_strain >= 0.08:
        function_category = "moderately_reduced"
    else:
        function_category = "severely_reduced"

    return {
        "phase_strains": phase_strains,
        "peak_strain": round(peak_strain, 4),
        "peak_phase": peak_phase,
        "function_category": function_category,
        "tracking_points": len(points),
    }


class TissueTrackingAgent:
    """Sub-agent for tissue tracking."""

    def __init__(self):
        self.agent_name = "TissueTrackingAgent"

    def evaluate(self, points: List[TissuePoint], reference_phase: int = 0) -> Dict[str, Any]:
        """Evaluate tissue tracking."""
        result = calculate_tissue_strain(points, reference_phase)
        alerts = []

        if "error" in result:
            return {"tracking_result": result, "alerts": [{"type": "NO_DATA", "severity": "ERROR",
                    "message": result["error"], "recommendation": "Provide tracking points."}]}

        if result["function_category"] in ("moderately_reduced", "severely_reduced"):
            alerts.append({
                "type": "REDUCED_MYOCARDIAL_FUNCTION", "severity": "WARNING",
                "message": f"Peak strain {result['peak_strain']:.4f} ({result['function_category']}).",
                "recommendation": "Assess for ischemia or cardiomyopathy. Consider stress testing."
            })

        return {"tracking_result": result, "alerts": alerts}
