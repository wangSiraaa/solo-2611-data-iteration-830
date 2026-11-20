import tempfile
from pathlib import Path

from prepquant.contracts import evaluate_contract
from prepquant.ledger import Ledger
from prepquant.units import UnitArchive, beads_from_ledger
from prepquant.units_demo import units_demo


def test_units_demo_runs_and_checks():
    d = tempfile.mkdtemp()
    report = units_demo(Path(d))
    assert len(report["checks"]) == 5
    assert (Path(d) / "units_demo_report.json").exists()


def test_ledger_to_contract_end_to_end():
    """A real dilution ledger's transfer record drives the bead archive and contract."""
    d = tempfile.mkdtemp()
    L = Ledger(Path(d) / "ledger.db")
    L.apply_batch(
        [
            {"op": "sample", "container": "P", "volume": "4", "sample": "s1"},
            {"op": "transfer", "source": "P", "volume": "2", "destination": "A"},
        ]
    )
    archive = UnitArchive()
    archive.apply_batch(beads_from_ledger(L, {"P": (2, 2)}))
    product = archive.joint_product("A", "P")
    report = evaluate_contract(product, {"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "P"})
    assert report["accept_prob"] == "5/6"
    assert report["expected_draws"] == "3/2"
    L.close()
