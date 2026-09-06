"""Produce the blind MSLN benchmark output and freeze it.

Reads no answers. The literature panel is revealed only after this output is
frozen and committed, and the guard refuses to freeze if the answers are
present on disk or if any held-out value has reached the algorithm.
"""

from __future__ import annotations

import json
import sys

from benchmark import blind
from car_pipeline.api import pipeline, server
from car_pipeline.stages import binder_check, binder_profile, binder_ranking

TARGET = "MSLN"
INDICATION = "Pancreatic Ductal Adenocarcinoma"


def build() -> dict:
    """The benchmark output, in the shape section 8 specifies."""
    run = pipeline.run(INDICATION, progress=lambda *a: None)
    rows = server._checked_rows(run, server._binder_rows(run))
    mine = [b for b in rows if b["target_id"] == TARGET]
    ranked = binder_ranking.rank(mine)
    by_id = {b["binder_id"]: b for b in ranked["binders"]}

    target = binder_profile.normalise_target(
        (run.get("surface_by_gene") or {}).get(TARGET))

    binders = []
    for row in mine:
        scored = by_id.get(row["binder_id"], {})
        sanitised = binder_profile.sanitise(
            row.get("heavy_sequence"), domain_declared=bool(row.get("format")),
            note=row.get("format") or "")
        binders.append({
            "binder": row["identifier"],
            "binder_id": row["binder_id"],
            "target_match": row["target_match"],
            "wrong_antigen_flag": row.get("wrong_antigen_flag"),
            "recorded_antigens": row.get("recorded_antigens"),
            "origin": row["origin"],
            "route": row["route"],
            "affinity": {"value": None, "unit": "nM", "type": "UNKNOWN",
                         "confidence": None},
            "epitope": {"location": None, "membrane_proximity": None,
                        "confidence": None},
            "sequence_original": sanitised["sequence_original"],
            "sequence_clean": sanitised["sequence_clean"],
            "modifications": sanitised["modifications"],
            "annotation": binder_profile.annotate(
                row.get("heavy_sequence"), row.get("light_sequence"),
                row.get("format")),
            "developability": binder_profile.developability(
                row.get("heavy_sequence")),
            "structural_evidence": row.get("structural_evidence"),
            "specificity_score": None,
            "binding_rank": scored.get("binding_rank"),
            "car_suitability_rank": scored.get("car_suitability_rank"),
            "binding": scored.get("binding"),
            "suitability": scored.get("suitability"),
            "on_front": scored.get("on_front"),
        })

    verdicts: dict[str, int] = {}
    for row in mine:
        verdicts[row["target_match"]] = verdicts.get(row["target_match"], 0) + 1

    return {
        "target": TARGET,
        "source_coverage": {
            "affinity_values_present": 0,
            "retrieved_candidates_pool_wide": len(rows),
            "statement": (
                "Affinity is UNKNOWN on every retrieved binder. No connected "
                "evidence release carries an affinity, KD or free-energy "
                "column. Fold error against a literature value and rank "
                "correlation against a literature ordering are therefore "
                "structurally uncomputable rather than computed and failed: "
                "there is no value on either side of the comparison. This is "
                "a source-coverage fact and should be read before any result."),
        },
        "retrieval": {
            "total_binders_found": len(mine),
            "structure_route": sum(1 for b in mine if b["route"] == "structure"),
            "sequence_route": sum(1 for b in mine if b["route"] == "sequence"),
            "sequences_available": sum(1 for b in mine if b["sequence_available"]),
            "target_match_counts": verdicts,
            "source_coverage": [
                "deposited structural complexes, searched by target accession",
                "named therapeutic antibodies, matched by recorded target",
            ],
        },
        "normalized_target": target,
        "binders": binders,
        "ranking": {
            "status": ranked["status"],
            "weight_versions": ranked["weight_versions"],
            "front": ranked["front"],
            "front_axes": ranked["front_axes"],
            "front_discriminates": ranked["front_discriminates"],
            "front_note": ranked["front_note"],
            "excluded_binders": ranked["excluded_binders"],
        },
        "validation": {
            "benchmark_revealed_after_prediction": True,
            "binder_retrieval_recall": None,
            "target_match_accuracy": None,
            "affinity_rank_correlation": None,
            "reason": (
                "Every validation figure is null because each compares against "
                "the literature panel, which has not been read. They are filled "
                "only after this output is frozen and committed."),
        },
    }


def main() -> int:
    """Build the blind output, guard it, freeze it."""
    print("guarding the blind rule", flush=True)
    state = blind.guard()
    print(f"  layer 1 literals            {state['layer_1_literals'] or 'clean'}")
    print(f"  layer 2 answers reachable   "
          f"{state['layer_2_answers_reachable'] or 'clean'}")
    print(f"  layer 3 ordering            {state['layer_3_ordering']['reason']}")

    print("\nbuilding the blind output", flush=True)
    payload = build()
    record = blind.freeze(payload)
    print(f"  frozen to {blind.FROZEN.relative_to(blind.ROOT).as_posix()}")
    print(f"  fingerprint {record['fingerprint']}")

    again = blind.verify_frozen(build())
    print(f"  re-running reproduces it: {again['reproduces']}")
    if not again["reproduces"]:
        print(f"  *** {again} ***")
        return 2

    print()
    print(payload["source_coverage"]["statement"])
    print()
    print(f"  retrieval: {payload['retrieval']['total_binders_found']} binder(s), "
          f"{payload['retrieval']['target_match_counts']}")
    print(f"  ranking:   {payload['ranking']['status']}, "
          f"front {len(payload['ranking']['front'])}, "
          f"discriminates {payload['ranking']['front_discriminates']}")
    print()
    print("  Commit this artifact before reading the literature panel.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
