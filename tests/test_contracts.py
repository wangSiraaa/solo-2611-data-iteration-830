import json
from fractions import Fraction

import pytest

from prepquant.contracts import (
    ContractError,
    evaluate_contract,
    format_report,
    parse_contract,
)
from prepquant.units import UnitArchive


def _archive_split(total, colored, first_units, first="A", parent="P"):
    """Declare ``parent`` with ``colored`` of ``total`` and move ``first_units`` to ``first``."""
    arc = UnitArchive()
    arc.apply_batch(
        [
            {"op": "beads", "box": parent, "colored": colored, "colorless": total - colored},
            {"op": "transfer", "source": parent, "units": first_units, "destination": first},
        ]
    )
    return arc


# --------------------------------------------------------------- the examples
def test_statement_example_two_of_four():
    # 4 beads, 2 colored; A=2, B=2; n=1 from A, a=0 r=2; m=1 from B, b=1.
    arc = _archive_split(4, 2, 2)
    product = arc.joint_product("A", "P")
    report = evaluate_contract(product, {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "P"})
    assert report["accept_prob"] == "5/6"
    assert report["expected_draws"] == "3/2"
    # per branch: x=0 accepts immediately (1/2); x=1 continues (1/2), accepts 2/3 of it.
    branches = {b["x"]: b for b in report["branches"]}
    assert branches[0]["status"] == "accept" and branches[0]["prob"] == "1/2"
    assert branches[1]["status"] == "continue" and branches[1]["prob"] == "1/2"
    assert branches[1]["accept_prob"] == "1/3"


def test_statement_example_b_zero():
    arc = _archive_split(4, 2, 2)
    product = arc.joint_product("A", "P")
    report = evaluate_contract(product, {"n": 1, "m": 1, "a": 0, "r": 2, "b": 0, "second_box": "P"})
    assert report["accept_prob"] == "1/2"
    assert report["expected_draws"] == "3/2"


def test_statement_example_one_colored_always_accepts():
    # 4 beads, 1 colored; A=2, B=2; same contract. Release is certain.
    arc = _archive_split(4, 1, 2)
    product = arc.joint_product("A", "P")
    report = evaluate_contract(product, {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "P"})
    assert report["accept_prob"] == "1"
    assert report["expected_draws"] == "5/4"


# ------------------------------------------------------- stage-1 branch detail
def test_early_reject_branch_is_exact():
    # Draw both units of a box that can hold 2 colored; r=2 => reject on x=2.
    arc = _archive_split(4, 2, 2)
    product = arc.joint_product("A", "P")
    report = evaluate_contract(product, {"n": 2, "m": 1, "a": 0, "r": 2, "b": 2, "second_box": "P"})
    branches = {b["x"]: b for b in report["branches"]}
    assert branches[2]["status"] == "reject"
    assert branches[2]["accept_prob"] == "0"
    # branch split sums to one
    s = report["stage1"]
    assert Fraction(s["accept"]) + Fraction(s["reject"]) + Fraction(s["continue"]) == 1


def test_expected_draws_uses_continue_mass():
    arc = _archive_split(4, 2, 2)
    product = arc.joint_product("A", "P")
    report = evaluate_contract(product, {"n": 1, "m": 2, "a": 0, "r": 2, "b": 1, "second_box": "P"})
    cont = Fraction(report["stage1"]["continue"])
    assert Fraction(report["expected_draws"]) == 1 + 2 * cont


# -------------------------------------------------- continue in the same box
def test_same_box_continuation_uses_remainder():
    # 4 beads/2 colored, A holds 2. Draw 1 from A; if it is colored continue in A,
    # which now has exactly 1 bead left.  Accept iff total colored <= 1.
    arc = _archive_split(4, 2, 2)
    product = arc.joint_product("A", "P")
    report = evaluate_contract(product, {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "continue"})
    assert report["contract"]["same_box"] is True
    assert report["accept_prob"] == "5/6"
    # A naive "redraw from a fresh full A" would give 2/3; the remainder gives 5/6.
    assert report["accept_prob"] != "2/3"


def test_same_box_continuation_refuses_overdraw_on_live_branch():
    arc = _archive_split(4, 2, 2)
    product = arc.joint_product("A", "P")
    # stage 1 can show x=1 with positive probability -> continue -> needs 2 more
    # units but A only has 1 left.
    with pytest.raises(ContractError):
        evaluate_contract(product, {"n": 1, "m": 2, "a": 0, "r": 2, "b": 2, "second_box": "continue"})


def test_no_continue_branch_allows_arbitrary_stage_two_size():
    arc = _archive_split(4, 2, 2)
    product = arc.joint_product("A", "P")
    # a >= n means every stage-1 outcome is an early accept; stage 2 never runs.
    report = evaluate_contract(product, {"n": 1, "m": 99, "a": 1, "r": 2, "b": 1, "second_box": "continue"})
    assert report["accept_prob"] == "1"
    assert report["expected_draws"] == "1"


def test_named_second_box_overdraw_refused():
    arc = _archive_split(4, 2, 2)
    product = arc.joint_product("A", "P")
    with pytest.raises(ContractError):
        evaluate_contract(product, {"n": 1, "m": 3, "a": 0, "r": 2, "b": 1, "second_box": "P"})


# ------------------------------------------------------------- input hygiene
def test_parse_contract_validates():
    with pytest.raises(ContractError):
        parse_contract({"n": 1, "m": 1, "a": 2, "r": 1, "b": 1})  # a < r violated
    with pytest.raises(ContractError):
        parse_contract({"n": 1, "m": 1, "a": 0, "r": 2})  # missing b
    with pytest.raises(ContractError):
        parse_contract({"n": "1/2", "m": 1, "a": 0, "r": 2, "b": 1})
    with pytest.raises(ContractError):
        parse_contract({"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "surprise": 1})


def test_named_second_box_must_be_in_product():
    arc = _archive_split(4, 2, 2)
    product = arc.joint_product("A", "P")
    with pytest.raises(ContractError):
        evaluate_contract(product, {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "Q"})


def test_stage_one_overdraw_refused():
    arc = _archive_split(4, 2, 2)
    product = arc.joint_product("A", "P")
    with pytest.raises(ContractError):
        evaluate_contract(product, {"n": 3, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "P"})


def test_malformed_product_rejected():
    with pytest.raises(ContractError):
        evaluate_contract({"boxes": [{"box": "A", "units": 2}], "joint": [{"counts": [1], "prob": "1/2"}]},
                          {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "continue"})
    # probabilities that do not sum to one
    with pytest.raises(ContractError):
        evaluate_contract(
            {"boxes": [{"box": "A", "units": 2}], "joint": [{"counts": [1], "prob": "1/2"}]},
            {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "continue"},
        )


def test_contract_reads_only_product_and_contract():
    # Evaluate a hand-written product with no archive at all.
    product = {
        "boxes": [{"box": "A", "units": 2}, {"box": "B", "units": 2}],
        "joint": [
            {"counts": [0, 2], "prob": "1/6"},
            {"counts": [1, 1], "prob": "2/3"},
            {"counts": [2, 0], "prob": "1/6"},
        ],
    }
    report = evaluate_contract(product, {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "B"})
    assert report["accept_prob"] == "5/6"


# -------------------------------------------------------- inspection of product
def test_contract_reflects_upstream_inspection():
    # Baseline (no inspection): A's early-accept mass for x=0 is 1/2.
    base = _archive_split(4, 2, 2)
    base_report = evaluate_contract(
        base.joint_product("A", "P"), {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "P"}
    )
    assert base_report["stage1"]["accept"] == "1/2"

    # After an inspection that already says A holds at least 1 colored, drawing a
    # single colorless bead is less likely: x=0 mass drops from 1/2 to 2/5.
    arc = _archive_split(4, 2, 2)
    arc.apply_batch([{"op": "inspect", "box": "A", "min": 1, "max": 2}])
    report = evaluate_contract(
        arc.joint_product("A", "P"), {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "P"}
    )
    assert report["stage1"]["accept"] == "2/5"
    assert report["accept_prob"] != base_report["accept_prob"]


def test_report_round_trips_through_json():
    arc = _archive_split(4, 2, 2)
    product = json.loads(json.dumps(arc.joint_product("A", "P")))
    report = evaluate_contract(product, {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "P"})
    text = format_report(report)
    assert "release probability: 5/6" in text
    assert "expected units drawn: 3/2" in text
