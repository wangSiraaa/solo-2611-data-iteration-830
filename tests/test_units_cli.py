import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def run_units(tmp_path, *args, input_text=None):
    cmd = [
        sys.executable,
        "-m",
        "prepquant.cli",
        "--units",
        str(tmp_path / "units.db"),
        *args,
    ]
    return subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, input=input_text)


def test_units_cli_full_flow(tmp_path):
    ops = [
        {"op": "beads", "box": "P", "colored": 2, "colorless": 2},
        {"op": "transfer", "source": "P", "units": 2, "destination": "A"},
    ]
    ops_file = tmp_path / "uops.json"
    ops_file.write_text(json.dumps(ops), encoding="utf-8")
    r = run_units(tmp_path, "units-prep", str(ops_file))
    assert r.returncode == 0, r.stderr

    r = run_units(tmp_path, "units-boxes", "--json")
    assert r.returncode == 0, r.stderr
    boxes = json.loads(r.stdout)
    assert boxes == [{"box": "A", "units": 2}, {"box": "P", "units": 2}]

    prod = tmp_path / "prod.json"
    r = run_units(tmp_path, "units-joint", "A", "P", "--out", str(prod))
    assert r.returncode == 0, r.stderr
    payload = json.loads(prod.read_text(encoding="utf-8"))
    assert len(payload["joint"]) == 3

    contract_file = tmp_path / "contract.json"
    contract_file.write_text(
        json.dumps({"n": 1, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "P"}), encoding="utf-8"
    )
    r = run_units(tmp_path, "contract", "--product", str(prod), "--contract", str(contract_file))
    assert r.returncode == 0, r.stderr
    assert "release probability: 5/6" in r.stdout

    r = run_units(tmp_path, "contract", "--product", str(prod), "--contract", str(contract_file), "--json")
    report = json.loads(r.stdout)
    assert report["accept_prob"] == "5/6" and report["expected_draws"] == "3/2"


def test_units_cli_contradiction_is_an_error(tmp_path):
    ops = [
        {"op": "beads", "box": "P", "colored": 1, "colorless": 3},
        {"op": "inspect", "box": "P", "min": 4, "max": 4},
    ]
    f = tmp_path / "bad.json"
    f.write_text(json.dumps(ops), encoding="utf-8")
    r = run_units(tmp_path, "units-prep", str(f))
    assert r.returncode != 0
    assert "zero conditional probability" in r.stderr or "impossible" in r.stderr
    # nothing committed
    r = run_units(tmp_path, "units-boxes", "--json")
    assert json.loads(r.stdout) == []


def test_units_cli_bad_contract(tmp_path):
    ops = [{"op": "beads", "box": "P", "colored": 2, "colorless": 2}]
    f = tmp_path / "uops.json"
    f.write_text(json.dumps(ops), encoding="utf-8")
    assert run_units(tmp_path, "units-prep", str(f)).returncode == 0
    prod = tmp_path / "prod.json"
    assert run_units(tmp_path, "units-joint", "P", "--out", str(prod)).returncode == 0
    contract = tmp_path / "c.json"
    contract.write_text(json.dumps({"n": 5, "m": 1, "a": 0, "r": 2, "b": 1, "second_box": "continue"}),
                        encoding="utf-8")
    r = run_units(tmp_path, "contract", "--product", str(prod), "--contract", str(contract))
    assert r.returncode != 0
    assert "asks for 5" in r.stderr
