"""Exact two-stage acceptance sampling of a whole-unit joint product.

An inspector reads only a **joint product** payload (the exact joint law of two
boxes' colored counts, produced by :meth:`prepquant.units.UnitArchive.joint_product`)
and a sampling contract.  Nothing else is consulted.

Contract shape (all counts are whole numbers; ``a < r``):

* Stage 1: draw ``n`` units **without replacement** from the first box.  Let ``X``
  be the colored count.
  - ``X <= a``  -> early **accept** (release).
  - ``X >= r``  -> early **reject** (return).
  - ``a < X < r`` -> **continue** to stage 2.
* Stage 2 (only on the continue branch): draw ``m`` units without replacement
  from a second box.  By default that box is the *first* box continuing with its
  remaining units; the contract may instead name the other box of the product.
  Naming the first box again (``second_box`` equal to the first box's id) is the
  same as "continue".  Let ``Y`` be the new colored count; the final decision is
  **accept** iff ``X + Y <= b``.

Everything is computed as an exact :class:`fractions.Fraction`.  For each possible
stage-1 colored count ``X = x`` the report gives the probability mass, the
early-accept / early-reject / continue probabilities, and the eventual total
accept probability through that branch.  The overall release probability and the
expected number of units drawn are exact.

A contract that would have a stage-1 continue branch ask for more units than the
chosen stage-2 box actually holds *on a positive-probability branch* is refused:
that branch can really occur, so promising to draw more than exists is invalid.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Mapping, Sequence

from .rational import ZERO, qstr
from .units import hypergeom_cdf, hypergeom_pmf


class ContractError(Exception):
    """The contract cannot be honoured by the product as described."""


# ------------------------------------------------------------------ parsing
def _whole(value, field: str, minimum: int = 0) -> int:
    if isinstance(value, bool):
        raise ContractError(f"{field} must be a whole number, not a bool")
    if isinstance(value, int):
        result = value
    elif isinstance(value, str):
        try:
            frac = Fraction(value.strip())
        except (ValueError, ZeroDivisionError) as exc:
            raise ContractError(f"{field}={value!r} is not an exact number") from exc
        if frac.denominator != 1:
            raise ContractError(f"{field}={value!r} is not a whole number")
        result = frac.numerator
    else:
        raise ContractError(f"{field} must be an int or integer string, got {type(value).__name__}")
    if result < minimum:
        raise ContractError(f"{field} must be >= {minimum}, got {result}")
    return result


def parse_contract(contract: Mapping) -> dict:
    """Validate a contract mapping and return a normalised dict of ints."""
    if not isinstance(contract, Mapping):
        raise ContractError("a contract must be an object")
    allowed = {"n", "m", "a", "r", "b", "second_box", "continue_box"}
    extra = set(contract) - allowed
    if extra:
        raise ContractError(f"contract has unknown field(s): {', '.join(sorted(extra))}")
    for field in ("n", "m", "a", "r", "b"):
        if field not in contract:
            raise ContractError(f"contract is missing '{field}'")
    n = _whole(contract["n"], "n", minimum=0)
    m = _whole(contract["m"], "m", minimum=0)
    a = _whole(contract["a"], "a", minimum=0)
    r = _whole(contract["r"], "r", minimum=1)
    b = _whole(contract["b"], "b", minimum=0)
    if not a < r:
        raise ContractError(f"need a < r; got a={a}, r={r}")
    # Which box stage 2 draws from.  ``second_box`` names a box id; the sentinel
    # "continue" (default) means keep drawing the first box with its remainder.
    second_box = contract.get("second_box", contract.get("continue_box", "continue"))
    if not isinstance(second_box, str) or not second_box:
        raise ContractError("second_box must be a box id or 'continue'")
    if second_box == "continue":
        second_box = None
    return {"n": n, "m": m, "a": a, "r": r, "b": b, "second_box": second_box}


# --------------------------------------------------------------- evaluation
def _joint_rows(payload: Mapping) -> tuple[list[dict], list[dict]]:
    try:
        boxes = payload["boxes"]
        rows = payload["joint"]
    except (KeyError, TypeError) as exc:
        raise ContractError(f"malformed joint product: missing {exc}") from exc
    if len(boxes) not in (1, 2):
        raise ContractError("a joint product must describe one or two boxes")
    normalised = []
    for box in boxes:
        cid = box.get("box")
        if not isinstance(cid, str) or not cid:
            raise ContractError("each joint product box needs a non-empty 'box' id")
        units = _whole(box.get("units"), "box units", minimum=1)
        normalised.append({"box": cid, "units": units})
    if len({b["box"] for b in normalised}) != len(normalised):
        raise ContractError("joint product repeats a box id")
    parsed_rows = []
    total = ZERO
    for row in rows:
        counts = row.get("counts")
        if not isinstance(counts, Sequence) or len(counts) != len(normalised):
            raise ContractError("each joint row must give one count per box")
        parsed_counts = tuple(_whole(c, "joint count") for c in counts)
        for box, value in zip(normalised, parsed_counts):
            if value > box["units"]:
                raise ContractError(
                    f"joint row says box {box['box']!r} has {value} colored but it holds {box['units']}"
                )
        try:
            prob = Fraction(str(row.get("prob")))
        except (ValueError, ZeroDivisionError) as exc:
            raise ContractError(f"joint row has a non-numeric probability: {row.get('prob')!r}") from exc
        if prob < 0:
            raise ContractError("joint probabilities must be non-negative")
        parsed_rows.append({"counts": parsed_counts, "prob": prob})
        total += prob
    if total != 1:
        raise ContractError(f"joint product probabilities sum to {qstr(total)}, expected 1")
    return normalised, parsed_rows


def evaluate_contract(joint_product: Mapping, contract: Mapping) -> dict:
    """Assess a two-stage sampling contract exactly against a joint product."""
    boxes, rows = _joint_rows(joint_product)
    spec = parse_contract(contract)
    n, m, a, r, b = spec["n"], spec["m"], spec["a"], spec["r"], spec["b"]

    first_box = boxes[0]["box"]
    first_size = boxes[0]["units"]
    if n > first_size:
        raise ContractError(
            f"stage 1 asks for {n} units from {first_box!r}, which holds only {first_size}"
        )

    # Resolve the stage-2 box: either the same first box (with its post-stage-1
    # remainder) or a named other box.  Naming the first box explicitly is the
    # same thing as "continue".
    if spec["second_box"] is None or spec["second_box"] == first_box:
        second_box = first_box
        same_box = True
        second_size = first_size
    else:
        same_box = False
        second_box = spec["second_box"]
        matches = [bx for bx in boxes if bx["box"] == second_box]
        if not matches:
            raise ContractError(
                f"stage 2 names box {second_box!r}, which is not part of this joint product"
            )
        second_size = matches[0]["units"]

    # Per stage-1 outcome record.  A branch through x contributes when a < x < r.
    per_x: list[dict] = []
    prob_by_x: dict[int, Fraction] = {}
    for row in rows:
        c1 = row["counts"][0]
        for x in range(0, min(n, c1) + 1):
            weight = hypergeom_pmf(first_size, c1, n, x)
            if weight == 0:
                continue
            prob_by_x[x] = prob_by_x.get(x, ZERO) + row["prob"] * weight
    # Every row/pair is accounted for; now build branch statistics.
    for x in sorted(prob_by_x):
        px = prob_by_x[x]
        if x <= a:
            status, through = "accept", ZERO
        elif x >= r:
            status, through = "reject", ZERO
        else:
            status, through = "continue", ZERO
            # Sum over every joint row the eventual conditional acceptance
            # probability, continuing without replacement.
            for row in rows:
                c1 = row["counts"][0]
                weight = hypergeom_pmf(first_size, c1, n, x)
                if weight == 0:
                    continue
                branch = row["prob"] * weight
                if same_box:
                    remaining = first_size - n
                    residual_colored = c1 - x
                    if remaining < m:
                        raise ContractError(
                            f"continuing in {first_box!r} at x={x} would draw {m} units "
                            f"but only {remaining} remain"
                        )
                    p = hypergeom_cdf(remaining, residual_colored, m, b - x)
                else:
                    c2 = row["counts"][1]
                    if second_size < m:
                        raise ContractError(
                            f"stage 2 asks for {m} units from {second_box!r}, "
                            f"which holds only {second_size}"
                        )
                    p = hypergeom_cdf(second_size, c2, m, b - x)
                through += branch * p
        per_x.append(
            {
                "x": x,
                "prob": px,
                "status": status,
                "accept_prob": (px if status == "accept" else through),
            }
        )

    # Overall exact figures.
    total_accept = sum((entry["accept_prob"] for entry in per_x), ZERO)
    cont_prob = sum((entry["prob"] for entry in per_x if entry["status"] == "continue"), ZERO)
    expected_draws = Fraction(n) + Fraction(m) * cont_prob

    # Verify the split of stage-1 mass into accept / reject / continue is exact.
    stage1_accept = sum((e["prob"] for e in per_x if e["status"] == "accept"), ZERO)
    stage1_reject = sum((e["prob"] for e in per_x if e["status"] == "reject"), ZERO)
    if stage1_accept + stage1_reject + cont_prob != 1:
        raise ContractError("internal error: stage-1 branch probabilities do not sum to 1")

    return {
        "contract": {
            "n": n, "m": m, "a": a, "r": r, "b": b,
            "first_box": first_box,
            "second_box": second_box,
            "same_box": same_box,
        },
        "branches": [
            {
                "x": entry["x"],
                "prob": qstr(entry["prob"]),
                "status": entry["status"],
                "accept_prob": qstr(entry["accept_prob"]),
            }
            for entry in per_x
        ],
        "stage1": {
            "accept": qstr(stage1_accept),
            "reject": qstr(stage1_reject),
            "continue": qstr(cont_prob),
        },
        "accept_prob": qstr(total_accept),
        "expected_draws": qstr(expected_draws),
        "exact": {
            "accept_prob": total_accept,
            "expected_draws": expected_draws,
        },
    }


def format_report(report: Mapping) -> str:
    """Human-readable rendering of :func:`evaluate_contract` output."""
    c = report["contract"]
    lines = [
        f"stage 1: draw {c['n']} from {c['first_box']!r}; x<= {c['a']} accept, x>= {c['r']} reject, else continue",
        f"stage 2: draw {c['m']} from "
        + (f"{c['first_box']!r} (continue same box)" if c["same_box"] else f"{c['second_box']!r}")
        + f"; final accept iff total <= {c['b']}",
        "  x   P(X=x)   branch     eventual accept",
    ]
    for br in report["branches"]:
        lines.append(
            f"  {br['x']:>1}   {br['prob']:>7}   {br['status']:<9}  {br['accept_prob']}"
        )
    lines.append(
        f"stage-1 split: accept {report['stage1']['accept']}, "
        f"reject {report['stage1']['reject']}, continue {report['stage1']['continue']}"
    )
    lines.append(f"release probability: {report['accept_prob']}")
    lines.append(f"expected units drawn: {report['expected_draws']}")
    return "\n".join(lines)
