#!/usr/bin/env python3
"""
LGETransmuralityAgent: late gadolinium enhancement transmurality and viability.

Per AHA 17-segment model, each segment's scar transmurality classifies
viability for revascularization decisions:

    0%            -> normal myocardium
    1-49%         -> hibernating / viable (revascularization benefit likely)
    >50% (>=50%)  -> non-viable scar

Segments map to coronary territories:
    LAD : 1,2,7,8,13,14,17 (+ apical 15 in wrap variants)
    RCA : 3,4,9,10,15
    LCx : 5,6,11,12,16

Outputs a per-territory viability report and an ASCII bull's-eye map.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional


TERRITORY_MAP = {
    "LAD": {1, 2, 7, 8, 13, 14, 17},
    "RCA": {3, 4, 9, 10, 15},
    "LCx": {5, 6, 11, 12, 16},
}
SEGMENT_TO_TERRITORY = {s: t for t, segs in TERRITORY_MAP.items() for s in segs}


@dataclass
class SegmentScar:
    segment_number: int                 # AHA 1-17
    scar_fraction: float                # 0.0-1.0 transmural extent


def classify_viability(scar_fraction: float) -> str:
    if scar_fraction <= 0.001:
        return "normal"
    if scar_fraction < 0.50:
        return "viable_hibernating"
    return "non_viable"


def territory_of(segment: int) -> Optional[str]:
    return SEGMENT_TO_TERRITORY.get(segment)


def analyze_lge(segments: List[SegmentScar]) -> Dict[str, Any]:
    if not segments:
        raise ValueError("no segments provided")

    def label(seg: int) -> str:
        return f"{['','basal','mid','apical'][min((seg-1)//4+1 if seg<17 else 3,3)]}" \
               f" seg{seg}" if seg != 17 else "apex seg17"

    per_segment = []
    for s in segments:
        if s.scar_fraction < 0 or s.scar_fraction > 1:
            raise ValueError(f"segment {s.segment_number}: scar fraction out of range")
        cls = classify_viability(s.scar_fraction)
        per_segment.append({
            "segment": s.segment_number,
            "territory": territory_of(s.segment_number),
            "scar_pct": round(s.scar_fraction * 100.0),
            "viability": cls,
        })

    territories: Dict[str, Dict[str, Any]] = {}
    for t_name, seg_ids in TERRITORY_MAP.items():
        rows = [p for p in per_segment if p["segment"] in seg_ids]
        if not rows:
            continue
        nonviable = [p["segment"] for p in rows if p["viability"] == "non_viable"]
        hibernating = [p["segment"] for p in rows if p["viability"] == "viable_hibernating"]
        territories[t_name] = {
            "segments_scored": len(rows),
            "non_viable_segments": nonviable,
            "hibernating_segments": hibernating,
            "revascularization_value": (
                "high - predominantly viable myocardium" if len(nonviable) <= 1 and hibernating
                else "limited - predominantly non-viable scar" if len(nonviable) >= 2
                else "intermediate"),
        }

    total_scar_pct = round(sum(p["scar_pct"] for p in per_segment)
                           / len(per_segment), 1)
    guidance: List[str] = []
    if total_scar_pct > 35:
        guidance.append("Global scar burden high; assess with viability imaging before CABG/PCI")
    lad_rows = territories.get("LAD", {})
    if len(lad_rows.get("non_viable_segments", [])) >= 2:
        guidance.append("Extensive LAD territory scar: consider ICD evaluation")
    if any(len(t.get("hibernating_segments", [])) > 0 for t in territories.values()):
        guidance.append("Hibernating segments present: ischemia evaluation indicated")

    return {
        "per_segment": per_segment,
        "territories": territories,
        "mean_transmurality_pct": total_scar_pct,
        "clinical_guidance": guidance,
        "bullseye": render_bullseye(per_segment),
    }


def render_bullseye(per_segment: List[Dict[str, Any]]) -> List[str]:
    """Compact ASCII bull's-eye: basal/mid/apical rings plus apex."""
    by_seg = {p["segment"]: p["scar_pct"] for p in per_segment}

    def glyph(seg: int) -> str:
        pct = by_seg.get(seg)
        if pct is None:
            return "?"
        if pct == 0:
            return "."
        if pct < 25:
            return "+"
        if pct < 50:
            return "#"
        return "@"

    ring_labels = ["basal", "mid", "apical"]
    lines = ["LGE bull's-eye ( . none | + mild | # sub-50% | @ transmural )"]
    for r, lab in enumerate(ring_labels):
        base = r * 6 + 1
        cells = " ".join(glyph(base + i) for i in range(6))
        lines.append(f"{lab:7s} [{cells}]")
    apex_cells = " ".join(glyph(13 + i) for i in range(4))
    lines.append(f"{'apical4':7s} [{apex_cells}]")
    lines.append(f"{'apex':7s} [  {glyph(17)}  ]")
    return lines


if __name__ == "__main__":
    scars = [
        SegmentScar(1, 0.0), SegmentScar(2, 0.0), SegmentScar(3, 0.9),
        SegmentScar(4, 0.8), SegmentScar(5, 0.0), SegmentScar(6, 0.0),
        SegmentScar(7, 0.3), SegmentScar(8, 0.45), SegmentScar(9, 0.95),
        SegmentScar(10, 0.85), SegmentScar(11, 0.0), SegmentScar(12, 0.0),
        SegmentScar(13, 0.2), SegmentScar(14, 0.0), SegmentScar(15, 0.7),
        SegmentScar(16, 0.0), SegmentScar(17, 0.1),
    ]
    report = analyze_lge(scars)
    print("LGE transmurality & viability report")
    print("-" * 52)
    print(f"mean global transmurality: {report['mean_transmurality_pct']}%")
    for line in report["bullseye"]:
        print(line)
    print("\nTerritory summary:")
    for t, info in report["territories"].items():
        print(f"  {t}: non-viable={info['non_viable_segments']} "
              f"hibernating={info['hibernating_segments']} "
              f"-> {info['revascularization_value']}")
    print("\nGuidance:")
    for g in report["clinical_guidance"]:
        print(f"  * {g}")
