"""Reproducible whole-unit portioning / inspection / sampling demo.

Run ``python3 -m prepquant.units_demo`` (or ``python3 -m prepquant.demo`` which
also calls :func:`units_demo`).  Everything is exact ``Fraction`` arithmetic and
is asserted against hand-computed values, so the demo is a self-check as well as
an illustration.

Scenarios
---------
1. **Portion and inspect.**  A parent holds 4 beads, 2 colored.  Split off A (2
   units); the remainder P keeps 2.  The joint law of (A, P) is printed to show the
   anti-correlation.  Then inspect A to be at least 1 colored and watch the joint
   renormalise.
2. **The two-stage contract from the statement.**  4 beads / 2 colored split into
   A and B (2 each); draw 1 from A, continue into B, accept iff the running total
   is <= 1.  Release probability 5/6, expected draws 3/2.  Tightening ``b`` to 0
   drops release to 1/2.  A one-colored parent releases with certainty (draws 5/4).
3. **Re-merging split-off boxes.**  Two portions recombined into M keep the joint
   anti-correlation with the remainder.
4. **Continuing in the same box.**  Draw from A, then keep drawing A with its
   *remaining* units -- not a fresh full box.  The report distinguishes the two.
5. **A contradictory inspection** is refused with an explicit message.
"""

from __future__ import annotations

import json
from pathlib import Path

from .contracts import ContractError, evaluate_contract, format_report
from .units import UnitArchive, UnitContradiction


class Failure(Exception):
    pass


def check(condition: bool, message: str) -> None:
    if not condition:
        raise Failure(message)


def _fresh() -> UnitArchive:
    return UnitArchive()  # in-memory


def scenario_portion_and_inspect(report: dict) -> None:
    arc = _fresh()
    arc.apply_batch(
        [
            {"op": "beads", "box": "P", "colored": 2, "colorless": 2},
            {"op": "transfer", "source": "P", "units": 2, "destination": "A"},
        ]
    )
    joint = arc.joint_distribution("A", "P")
    check(
        {k: str(v) for k, v in joint.items()}
        == {(0, 2): "1/6", (1, 1): "2/3", (2, 0): "1/6"},
        f"split joint wrong: {joint}",
    )
    record = {"step": "portion", "boxes": arc.boxes(), "joint": {f"{a},{b}": str(p) for (a, b), p in joint.items()}}
    # A non-destructive scan: A shows at least one colored bead.
    arc.apply_batch([{"op": "inspect", "box": "A", "min": 1, "max": 2}])
    conditioned = arc.joint_distribution("A", "P")
    check(
        {k: str(v) for k, v in conditioned.items()} == {(1, 1): "4/5", (2, 0): "1/5"},
        f"inspection renormalisation wrong: {conditioned}",
    )
    record["after_inspect"] = {f"{a},{b}": str(p) for (a, b), p in conditioned.items()}
    report["steps"].append(record)


def scenario_statement_contract(report: dict) -> None:
    arc = _fresh()
    arc.apply_batch(
        [
            {"op": "beads", "box": "P", "colored": 2, "colorless": 2},
            {"op": "transfer", "source": "P", "units": 2, "destination": "A"},
        ]
    )
    product = arc.joint_product("A", "P")
    contract = {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "P"}
    result = evaluate_contract(product, contract)
    check(result["accept_prob"] == "5/6", f"release expected 5/6, got {result['accept_prob']}")
    check(result["expected_draws"] == "3/2", f"draws expected 3/2, got {result['expected_draws']}")
    report["steps"].append(
        {
            "step": "contract_b1",
            "report": result,
            "text": format_report(result),
        }
    )

    tight = evaluate_contract(product, {**contract, "b": 0})
    check(tight["accept_prob"] == "1/2", f"b=0 release expected 1/2, got {tight['accept_prob']}")
    report["steps"].append({"step": "contract_b0", "accept_prob": tight["accept_prob"],
                            "expected_draws": tight["expected_draws"]})

    # One colored of four: release is certain, expected draws 5/4.
    slight = _fresh()
    slight.apply_batch(
        [
            {"op": "beads", "box": "P", "colored": 1, "colorless": 3},
            {"op": "transfer", "source": "P", "units": 2, "destination": "A"},
        ]
    )
    r = evaluate_contract(slight.joint_product("A", "P"), contract)
    check(r["accept_prob"] == "1", f"one-colored release expected 1, got {r['accept_prob']}")
    check(r["expected_draws"] == "5/4", f"one-colored draws expected 5/4, got {r['expected_draws']}")
    report["steps"].append({"step": "contract_one_colored", "accept_prob": r["accept_prob"],
                            "expected_draws": r["expected_draws"]})


def scenario_merge(report: dict) -> None:
    arc = _fresh()
    arc.apply_batch(
        [
            {"op": "beads", "box": "P", "colored": 3, "colorless": 3},
            {"op": "transfer", "source": "P", "units": 2, "destination": "A"},
            {"op": "transfer", "source": "P", "units": 2, "destination": "B"},
            {"op": "mix", "sources": [{"box": "A"}, {"box": "B"}], "destination": "M"},
        ]
    )
    check(arc.boxes() == [{"box": "M", "units": 4}, {"box": "P", "units": 2}], f"merge boxes wrong: {arc.boxes()}")
    joint = arc.joint_distribution("M", "P")
    for (cm, cp) in joint:
        check(cm + cp == 3, f"merged+remainder colored should sum to 3, got {cm}+{cp}")
    check(sum(joint.values()) == 1, "merged joint must sum to one")
    report["steps"].append(
        {"step": "merge", "boxes": arc.boxes(), "joint": {f"{a},{b}": str(p) for (a, b), p in joint.items()}}
    )


def scenario_same_box(report: dict) -> None:
    # 4 beads / 2 colored, A holds 2. Draw 1 from A; on a colored draw continue in
    # A with its single remaining unit.
    arc = _fresh()
    arc.apply_batch(
        [
            {"op": "beads", "box": "P", "colored": 2, "colorless": 2},
            {"op": "transfer", "source": "P", "units": 2, "destination": "A"},
        ]
    )
    product = arc.joint_product("A", "P")
    got = evaluate_contract(product, {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "continue"})
    check(got["accept_prob"] == "5/6", f"same-box release expected 5/6, got {got['accept_prob']}")
    check(got["contract"]["same_box"] is True, "same-box continuation should be flagged")
    # A naive fresh-full-box redraw would report 2/3; the remainder logic must not.
    check(got["accept_prob"] != "2/3", "same-box continuation must use the remainder, not a fresh box")
    # Asking for more than the remainder holds on a live branch is refused.
    refused = None
    try:
        evaluate_contract(product, {"n": 1, "m": 2, "a": 0, "r": 2, "b": 2, "second_box": "continue"})
    except ContractError as exc:
        refused = str(exc)
    check(refused is not None, "over-drawing the remainder on a live branch must be refused")
    report["steps"].append(
        {"step": "same_box", "accept_prob": got["accept_prob"], "expected_draws": got["expected_draws"],
         "text": format_report(got), "refusal": refused}
    )


def scenario_contradiction(report: dict) -> None:
    message = None
    try:
        arc = _fresh()
        arc.apply_batch(
            [
                {"op": "beads", "box": "P", "colored": 1, "colorless": 3},
                {"op": "inspect", "box": "P", "min": 4, "max": 4},
            ]
        )
    except UnitContradiction as exc:
        message = str(exc)
    check(message is not None, "an impossible inspection must be refused")
    report["steps"].append({"step": "contradiction", "refused": message})


def _json_safe(obj):
    """Drop the ``exact`` Fraction cache (kept for programmatic callers) so the
    report serialises; everything user-facing is already a qstr."""
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items() if k != "exact"}
    if isinstance(obj, list):
        return [_json_safe(v) for v in obj]
    return obj


def units_demo(out_dir: Path | None = None) -> dict:
    report = {"steps": [], "checks": []}
    scenario_portion_and_inspect(report)
    report["checks"].append("portion joint (0,2)=1/6,(1,1)=2/3,(2,0)=1/6; inspection renormalises to 4/5,1/5")
    scenario_statement_contract(report)
    report["checks"].append("statement contract 5/6 & draws 3/2; b=0 -> 1/2; one-colored -> 1 & 5/4")
    scenario_merge(report)
    report["checks"].append("re-merge keeps joint anti-correlation with remainder")
    scenario_same_box(report)
    report["checks"].append("same-box continuation uses the remainder (5/6, not 2/3) and refuses over-draw")
    scenario_contradiction(report)
    report["checks"].append("impossible inspection refused explicitly")
    if out_dir is not None:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "units_demo_report.json").write_text(
            json.dumps(_json_safe(report), indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
    return report


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(prog="prepquant.units_demo")
    parser.add_argument("--out", default=str(Path(__file__).resolve().parent.parent / "out"))
    args = parser.parse_args(argv)
    try:
        report = units_demo(Path(args.out))
    except Failure as exc:
        print(f"UNITS DEMO FAILED: {exc}")
        return 1
    print("units demo OK")
    for line in report["checks"]:
        print(f"  - {line}")
    print(f"report written to {Path(args.out) / 'units_demo_report.json'}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
