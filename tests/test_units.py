import os
import tempfile
from fractions import Fraction

import pytest

from prepquant.units import (
    UnitArchive,
    UnitContradiction,
    UnitError,
    beads_from_ledger,
    hypergeom_cdf,
    hypergeom_pmf,
)


@pytest.fixture()
def archive():
    d = tempfile.mkdtemp()
    A = UnitArchive(os.path.join(d, "units.db"))
    yield A
    A.close()


def op(**kw):
    return kw


# ---------------------------------------------------------------- hypergeometric
def test_hypergeom_basic():
    # draw 1 from 4 beads, 2 colored -> 1/2 each
    assert hypergeom_pmf(4, 2, 1, 0) == Fraction(1, 2)
    assert hypergeom_pmf(4, 2, 1, 1) == Fraction(1, 2)
    # draw 2 from 4, 2 colored -> 0,1,2 = 1/6, 4/6, 1/6
    assert [hypergeom_pmf(4, 2, 2, d) for d in range(3)] == [Fraction(1, 6), Fraction(4, 6), Fraction(1, 6)]
    assert sum(hypergeom_pmf(6, 3, 2, d) for d in range(3)) == 1
    # out of support is zero, not an error
    assert hypergeom_pmf(4, 2, 2, 3) == 0
    assert hypergeom_cdf(4, 2, 2, 1) == Fraction(5, 6)


def test_hypergeom_rejects_impossible_draws():
    with pytest.raises(ValueError):
        hypergeom_pmf(2, 1, 3, 1)  # draw more than the box holds
    with pytest.raises(ValueError):
        hypergeom_pmf(4, 5, 1, 0)  # more colored than total


# ------------------------------------------------------------- split / joint
def test_split_joint_is_anticorrelated(archive):
    archive.apply_batch(
        [
            op(op="beads", box="P", colored=2, colorless=2),
            op(op="transfer", source="P", units=2, destination="A"),
        ]
    )
    assert archive.boxes() == [{"box": "A", "units": 2}, {"box": "P", "units": 2}]
    joint = archive.joint_distribution("A", "P")
    # both boxes always hold 2 beads; the colored counts sum to 2
    assert joint == {
        (0, 2): Fraction(1, 6),
        (1, 1): Fraction(4, 6),
        (2, 0): Fraction(1, 6),
    }
    assert archive.marginal("A") == {0: Fraction(1, 6), 1: Fraction(4, 6), 2: Fraction(1, 6)}


def test_repeated_split_conditions_on_current_box(archive):
    # 8 beads, 4 colored; A gets 3, B gets 2 from the remainder, 3 stay.
    archive.apply_batch(
        [
            op(op="beads", box="P", colored=4, colorless=4),
            op(op="transfer", source="P", units=3, destination="A"),
            op(op="transfer", source="P", units=2, destination="B"),
        ]
    )
    assert archive.boxes() == [
        {"box": "A", "units": 3},
        {"box": "B", "units": 2},
        {"box": "P", "units": 3},
    ]
    # Every joint row sums colored over all three boxes to 4.
    for key, p in archive.joint.items():
        assert p > 0
        assert sum(count for _, count in key) == 4
    # marginal of B is the second draw from the depleted remainder, not fresh.
    marg = archive.marginal("B")
    assert sum(marg.values()) == 1


def test_sequential_splits_match_multivariate_hypergeometric():
    # Two conditional draws from one parent must equal the closed-form
    # multivariate hypergeometric partition law over the three final boxes.
    from math import comb

    archive = UnitArchive()
    archive.apply_batch(
        [
            op(op="beads", box="P", colored=3, colorless=3),
            op(op="transfer", source="P", units=2, destination="A"),
            op(op="transfer", source="P", units=2, destination="B"),
        ]
    )
    sizes = {"A": 2, "B": 2, "P": 2}
    total_colored = 3
    total = 6

    def weight(counts):
        w = 1
        for cid, c in counts.items():
            w *= comb(sizes[cid], c)
        return Fraction(w, comb(total, total_colored))

    brute = {}
    for a in range(3):
        for b in range(3):
            for p in range(3):
                if a + b + p == total_colored:
                    brute[(a, b, p)] = weight({"A": a, "B": b, "P": p})
    got = {}
    for key, prob in archive.joint.items():
        d = dict(key)
        got[(d["A"], d["B"], d["P"])] = prob
    assert got == brute


def test_mix_is_deterministic_but_keeps_correlation(archive):
    archive.apply_batch(
        [
            op(op="beads", box="P", colored=2, colorless=2),
            op(op="transfer", source="P", units=2, destination="A"),
            op(op="mix", sources=[{"box": "A"}], destination="C"),
        ]
    )
    # C is just A moved wholesale; the joint with the remainder P is unchanged.
    assert archive.boxes() == [{"box": "C", "units": 2}, {"box": "P", "units": 2}]
    assert archive.joint_distribution("C", "P") == archive.joint_distribution("C", "P")
    assert sum(p for p in archive.joint_distribution("C", "P").values()) == 1


def test_merge_two_boxes_recombine(archive):
    archive.apply_batch(
        [
            op(op="beads", box="P", colored=3, colorless=3),
            op(op="transfer", source="P", units=2, destination="A"),
            op(op="transfer", source="P", units=2, destination="B"),
            op(op="mix", sources=[{"box": "A"}, {"box": "B"}], destination="M"),
        ]
    )
    # A and B merged back into M (4 units); P has the remaining 2.
    assert archive.boxes() == [{"box": "M", "units": 4}, {"box": "P", "units": 2}]
    joint = archive.joint_distribution("M", "P")
    for (cm, cp), p in joint.items():
        assert cm + cp == 3  # 3 colored split between M and P


# ---------------------------------------------------------------- inspection
def test_inspection_conditions_and_renormalises(archive):
    archive.apply_batch(
        [
            op(op="beads", box="P", colored=2, colorless=2),
            op(op="transfer", source="P", units=2, destination="A"),
            op(op="inspect", box="A", min=2, max=2),
        ]
    )
    # Knowing A has both colored pins P to zero colored.
    assert archive.marginal("A") == {2: Fraction(1)}
    assert archive.joint_distribution("A", "P") == {(2, 0): Fraction(1)}
    assert sum(archive.joint.values()) == 1


def test_inspection_interval_keeps_mass(archive):
    archive.apply_batch(
        [
            op(op="beads", box="P", colored=2, colorless=2),
            op(op="transfer", source="P", units=2, destination="A"),
            op(op="inspect", box="A", min=1, max=2),
        ]
    )
    # Removes the (0,2) row; renormalise the rest.
    joint = archive.joint_distribution("A", "P")
    assert (0, 2) not in joint
    assert (1, 1) in joint and (2, 0) in joint
    assert sum(joint.values()) == 1
    # original masses 4/6 and 1/6 -> normalised 4/5 and 1/5
    assert joint[(1, 1)] == Fraction(4, 5)
    assert joint[(2, 0)] == Fraction(1, 5)


def test_contradictory_inspection_is_reported(archive):
    archive.apply_batch(
        [op(op="beads", box="P", colored=1, colorless=3)]
    )
    with pytest.raises(UnitContradiction):
        archive.apply_batch([op(op="inspect", box="P", min=4, max=4)])
    with pytest.raises(UnitContradiction):
        archive.apply_batch([op(op="inspect", box="P", min=5, max=9)])  # more than exist
    # impossible interval at declaration time is a plain error
    with pytest.raises(UnitError):
        archive.apply_batch([op(op="inspect", box="P", min=3, max=2)])


def test_inspection_zero_probability_after_split(archive):
    with pytest.raises(UnitContradiction):
        archive.apply_batch(
            [
                op(op="beads", box="P", colored=1, colorless=3),
                op(op="transfer", source="P", units=2, destination="A"),
                op(op="inspect", box="A", min=2, max=2),  # impossible: only 1 colored exists
            ]
        )


# ----------------------------------------------------------- rules / limits
def test_batch_rolls_back_on_error(archive):
    archive.apply_batch([op(op="beads", box="P", colored=2, colorless=2)])
    before_ops = [dict(o) for o in archive.ops]
    before_sizes = dict(archive.sizes)
    before_joint = dict(archive.joint)
    with pytest.raises(UnitContradiction):
        archive.apply_batch(
            [
                op(op="transfer", source="P", units=2, destination="X"),
                op(op="inspect", box="X", min=2, max=2),
                op(op="inspect", box="P", min=2, max=2),  # now impossible
            ]
        )
    assert archive.ops == before_ops
    assert archive.sizes == before_sizes
    assert archive.joint == before_joint


def test_overdraw_and_identity_errors(archive):
    archive.apply_batch([op(op="beads", box="P", colored=1, colorless=1)])
    with pytest.raises(UnitError):
        archive.apply_batch([op(op="transfer", source="P", units=3, destination="A")])
    with pytest.raises(UnitError):
        archive.apply_batch([op(op="beads", box="P", colored=1, colorless=1)])
    # destination must be empty
    archive.apply_batch([op(op="transfer", source="P", units=1, destination="A")])
    with pytest.raises(UnitError):
        archive.apply_batch([op(op="transfer", source="P", units=1, destination="A")])


def test_limits_beads_and_boxes(archive):
    with pytest.raises(UnitError):
        archive.apply_batch([op(op="beads", box="P", colored=10, colorless=10)])  # > 16 beads
    # five boxes exceeds the four-box limit
    with pytest.raises(UnitError):
        archive.apply_batch(
            [
                op(op="beads", box="P", colored=5, colorless=0),
                op(op="transfer", source="P", units=1, destination="A"),
                op(op="transfer", source="P", units=1, destination="B"),
                op(op="transfer", source="P", units=1, destination="C"),
                op(op="transfer", source="P", units=1, destination="D"),
            ]
        )


def test_fractional_units_rejected(archive):
    with pytest.raises(UnitError):
        archive.apply_batch([op(op="transfer", source="P", units="1/2", destination="A")])
    with pytest.raises(UnitError):
        archive.apply_batch([op(op="beads", box="P", colored="1.5", colorless=1)])


def test_string_integer_units_accepted(archive):
    archive.apply_batch(
        [
            op(op="beads", box="P", colored="2", colorless="2"),
            op(op="transfer", source="P", units="2", destination="A"),
        ]
    )
    assert archive.sizes == {"P": 2, "A": 2}


# -------------------------------------------------------- ledger integration
def test_beads_from_real_ledger_records():
    from prepquant.ledger import Ledger

    d = tempfile.mkdtemp()
    L = Ledger(os.path.join(d, "ledger.db"))
    L.apply_batch(
        [
            {"op": "sample", "container": "P", "volume": "4", "sample": "s1"},
            {"op": "transfer", "source": "P", "volume": "2", "destination": "A"},
        ]
    )
    records = beads_from_ledger(L, {"P": (2, 2)})
    A = UnitArchive()
    A.apply_batch(records)
    assert A.sizes == {"P": 2, "A": 2}
    assert A.joint_distribution("A", "P")[(1, 1)] == Fraction(2, 3)
    L.close()


def test_beads_from_ledger_rejects_non_unit_ops():
    from prepquant.ledger import Ledger

    d = tempfile.mkdtemp()
    L = Ledger(os.path.join(d, "ledger.db"))
    L.apply_batch(
        [
            {"op": "sample", "container": "P", "volume": "4", "sample": "s1"},
            {"op": "transfer", "source": "P", "volume": "2", "destination": "A"},
            {"op": "solvent", "container": "A", "volume": "1"},
        ]
    )
    with pytest.raises(UnitError):
        beads_from_ledger(L, {"P": (2, 2)})
    L.close()


# -------------------------------------------------------------- persistence
def test_persistence_roundtrip():
    d = tempfile.mkdtemp()
    path = os.path.join(d, "units.db")
    A = UnitArchive(path)
    A.apply_batch(
        [
            op(op="beads", box="P", colored=2, colorless=2),
            op(op="transfer", source="P", units=2, destination="A"),
            op(op="inspect", box="A", min=1, max=2),
        ]
    )
    joint_before = dict(A.joint)
    A.close()
    B = UnitArchive(path)
    assert B.joint == joint_before
    assert B.sizes == {"P": 2, "A": 2}
    B.close()
