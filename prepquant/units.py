"""Exact whole-unit (colored / colorless bead) portioning archive.

The dilution ledger tracks continuous volume and solute mass.  Next to it lives a
deliberately *discrete* archive for boxes that hold whole units: every volume unit
is exactly one bead, and every bead is either **colored** or **colorless**.  A
parent box declares how many beads of each kind it holds; a portioning step draws
``k`` beads uniformly **without replacement** and leaves the rest behind.

Design choices and their consequences
-------------------------------------
* **Joint state, not per-box marginals.**  Two boxes split from the same parent
  are negatively correlated: if the first took many colored beads the second must
  have taken few.  Collapsing to marginals would destroy that.  The exact state is
  therefore the joint distribution over the *vector of colored counts of every
  non-empty box*, keyed by ``((cid, count), ...)`` sorted by container id.  With
  at most 16 beads and at most four non-empty boxes this stays tiny, and every
  probability is a :class:`fractions.Fraction`.
* **Exchangeability makes the count vector sufficient.**  Which particular beads
  sit in a box never matters, only how many are colored.  Drawing ``k`` from a box
  of ``(t, c)`` yields ``d`` colored with the hypergeometric law
  ``C(c,d) C(t-c,k-d) / C(t,k)`` and leaves ``c-d`` behind.  Conditioning a later
  draw on "the box as it actually is now" is exactly conditioning on the current
  ``c`` in the joint state, so repeated splits are handled by repeated
  hypergeometric convolution -- no hand-rolled independence assumption.
* **Merges are deterministic.**  Re-combining split-off boxes moves every bead, so
  the colored counts simply add.  No randomness enters, but the joint state still
  remembers how the merged box co-varies with the rest.
* **Inspections are conditions; contradictions are first class.**  A
  non-destructive inspection only narrows a box's colored count to an interval.
  The matching mass of the joint is renormalised.  An inspection that leaves zero
  conditional probability (or that asks for more beads than exist, or with
  ``min > max``) raises :class:`UnitContradiction` so the conflict is never
  silently dropped.
* **Replay, atomic batches.**  The durable truth is the ordered log of
  ``beads / transfer / mix / inspect`` records -- the same transfer and mix record
  shapes as the dilution ledger -- and any batch that would break a rule leaves the
  stored state untouched.

The module also produces a self-contained **joint product** payload
(:meth:`UnitArchive.joint_product`) so a separate inspector can evaluate a
sampling contract from nothing but that payload and the contract itself;
see :mod:`prepquant.contracts`.
"""

from __future__ import annotations

import json
import math
import sqlite3
from fractions import Fraction
from pathlib import Path
from typing import Mapping, Sequence

from .rational import ZERO, qstr

MAX_BEADS = 16
MAX_BOXES = 4

UNIT_OP_KINDS = ("beads", "transfer", "mix", "inspect")


class UnitError(Exception):
    """A whole-unit record violated a rule; nothing was committed."""


class UnitContradiction(UnitError):
    """An inspection has zero conditional probability (or is itself impossible)."""


# --------------------------------------------------------------------- helpers
def _comb(a: int, b: int) -> int:
    return math.comb(a, b)


def hypergeom_pmf(total: int, colored: int, draw: int, hits: int) -> Fraction:
    """``P(exactly hits colored when draw beads come from (total, colored))``.

    Drawing ``draw`` beads uniformly without replacement from a box of ``total``
    beads of which ``colored`` are colored.  Out-of-support arguments return
    ``Fraction(0)``.
    """
    if draw < 0 or total < 0 or colored < 0 or colored > total:
        raise ValueError("invalid hypergeometric parameters")
    if draw > total:
        raise ValueError("cannot draw more units than the box holds")
    if hits < 0 or hits > draw or hits > colored or draw - hits > total - colored:
        return ZERO
    return Fraction(_comb(colored, hits) * _comb(total - colored, draw - hits), _comb(total, draw))


def hypergeom_cdf(total: int, colored: int, draw: int, limit: int) -> Fraction:
    """``P(at most ``limit`` colored in ``draw`` draws from (total, colored))``."""
    if limit < 0:
        return ZERO
    limit = min(limit, draw)
    return sum((hypergeom_pmf(total, colored, draw, y) for y in range(limit + 1)), ZERO)


def _int(value, field: str, minimum: int = 0) -> int:
    """Parse a whole number from an int or an exact integer string (no floats)."""
    if isinstance(value, bool):
        raise UnitError(f"{field} must be a whole number, not a bool")
    if isinstance(value, int):
        result = value
    elif isinstance(value, str):
        text = value.strip()
        try:
            frac = Fraction(text)
        except (ValueError, ZeroDivisionError) as exc:
            raise UnitError(f"{field}={value!r} is not an exact number") from exc
        if frac.denominator != 1:
            raise UnitError(f"{field}={value!r} is not a whole number of units")
        result = frac.numerator
    else:
        raise UnitError(f"{field} must be an int or integer string, got {type(value).__name__}")
    if result < minimum:
        raise UnitError(f"{field} must be >= {minimum}, got {result}")
    return result


def _key(pairs: Mapping[str, int]) -> tuple[tuple[str, int], ...]:
    """Canonical joint-state key: ``((box, colored_count), ...)`` sorted by box id."""
    return tuple(sorted(((cid, int(count)) for cid, count in pairs.items()), key=lambda item: item[0]))


def _sum_prob(joint: Mapping[tuple, Fraction]) -> Fraction:
    return sum(joint.values(), ZERO)


# ------------------------------------------------------------- op normalisation
def _normalise(op: Mapping, seq: int | None) -> dict:
    """Canonicalise one record, accepting the dilution ledger's field aliases."""

    def bad(msg: str) -> "UnitError":
        where = f"seq {seq}: " if seq is not None else ""
        return UnitError(where + msg)

    if not isinstance(op, Mapping):
        raise bad("a whole-unit record must be an object")
    kind = op.get("op")
    if kind not in UNIT_OP_KINDS:
        raise bad(f"unknown whole-unit operation {kind!r}; expected one of {', '.join(UNIT_OP_KINDS)}")

    label = op.get("label")

    if kind == "beads":
        box = op.get("box", op.get("container"))
        if not isinstance(box, str) or not box:
            raise bad("beads record needs a non-empty 'box'")
        allowed = {"op", "box", "container", "colored", "colorless", "label"}
        extra = set(op) - allowed
        if extra:
            raise bad(f"beads record has unknown field(s): {', '.join(sorted(extra))}")
        if "colored" not in op or "colorless" not in op:
            raise bad("beads record needs 'colored' and 'colorless' counts")
        colored = _int(op["colored"], "colored")
        colorless = _int(op["colorless"], "colorless")
        if colored + colorless == 0:
            raise bad("a declared box must hold at least one bead")
        return {"op": "beads", "box": box, "colored": colored, "colorless": colorless, "label": label}

    if kind == "transfer":
        source = op.get("source")
        destination = op.get("destination")
        if not isinstance(source, str) or not source:
            raise bad("transfer needs a non-empty 'source'")
        if not isinstance(destination, str) or not destination:
            raise bad("transfer needs a non-empty 'destination'")
        if source == destination:
            raise bad("cannot transfer into the source box itself")
        if ("units" in op) == ("volume" in op):
            raise bad("transfer needs exactly one of 'units' or 'volume'")
        units = _int(op.get("units", op.get("volume")), "units", minimum=1)
        allowed = {"op", "source", "destination", "units", "volume", "label"}
        extra = set(op) - allowed
        if extra:
            raise bad(f"transfer record has unknown field(s): {', '.join(sorted(extra))}")
        return {"op": "transfer", "source": source, "destination": destination, "units": units, "label": label}

    if kind == "mix":
        sources = op.get("sources")
        destination = op.get("destination")
        if not isinstance(sources, list) or not sources:
            raise bad("mix needs a non-empty 'sources' list")
        if not isinstance(destination, str) or not destination:
            raise bad("mix needs a non-empty 'destination'")
        allowed = {"op", "sources", "destination", "label"}
        extra = set(op) - allowed
        if extra:
            raise bad(f"mix record has unknown field(s): {', '.join(sorted(extra))}")
        parsed: list[dict] = []
        seen: set[str] = set()
        for item in sources:
            if isinstance(item, str):
                cid, units = item, None
            elif isinstance(item, Mapping):
                cid = item.get("box", item.get("container"))
                if not isinstance(cid, str) or not cid:
                    raise bad("each mix source object needs a non-empty 'box'")
                if ("units" in item) == ("volume" in item) and ("units" in item or "volume" in item):
                    raise bad("each mix source needs at most one of 'units' or 'volume'")
                raw = item.get("units", item.get("volume"))
                units = None if raw is None else _int(raw, "mix source units", minimum=1)
            else:
                raise bad("each mix source must be a box id or an object")
            if cid == destination:
                raise bad("destination cannot also be a mix source")
            if cid in seen:
                raise bad(f"box {cid!r} listed twice in one mix")
            seen.add(cid)
            parsed.append({"box": cid, "units": units})
        return {"op": "mix", "sources": parsed, "destination": destination, "label": label}

    # inspect
    box = op.get("box", op.get("container"))
    if not isinstance(box, str) or not box:
        raise bad("inspect needs a non-empty 'box'")
    if ("min" in op) or ("max" in op):
        if ("lo" in op) or ("hi" in op):
            raise bad("inspect needs either min/max or lo/hi, not both")
        lo, hi = op.get("min"), op.get("max")
    else:
        lo, hi = op.get("lo"), op.get("hi")
    if lo is None or hi is None:
        raise bad("inspect needs a lower and an upper bound")
    allowed = {"op", "box", "container", "min", "max", "lo", "hi", "label"}
    extra = set(op) - allowed
    if extra:
        raise bad(f"inspect record has unknown field(s): {', '.join(sorted(extra))}")
    lo_i = _int(lo, "inspect min")
    hi_i = _int(hi, "inspect max")
    if lo_i > hi_i:
        raise bad(f"inspect interval is empty: min {lo_i} > max {hi_i}")
    return {"op": "inspect", "box": box, "min": lo_i, "max": hi_i, "label": label}


# -------------------------------------------------------------------- archive
class UnitArchive:
    """Replayable log of whole-unit portioning records with an exact joint state."""

    def __init__(self, db_path: str | Path | None = None) -> None:
        self.db_path = None if db_path is None else str(db_path)
        if self.db_path is not None:
            self.conn: sqlite3.Connection | None = sqlite3.connect(self.db_path)
            self.conn.row_factory = sqlite3.Row
            self._create_schema()
        else:
            self.conn = None
        self.ops: list[dict] = []
        self.sizes: dict[str, int] = {}
        self.joint: dict[tuple[tuple[str, int], ...], Fraction] = {}
        self._reload()

    # -- persistence -------------------------------------------------------
    def _create_schema(self) -> None:
        assert self.conn is not None
        with self.conn:
            self.conn.executescript(
                """
                CREATE TABLE IF NOT EXISTS unit_ops(
                    seq INTEGER PRIMARY KEY, kind TEXT NOT NULL, payload TEXT NOT NULL
                );
                """
            )

    def _reload(self) -> None:
        if self.conn is None:
            self._replay([], persist=False)
            return
        rows = self.conn.execute("SELECT seq, kind, payload FROM unit_ops ORDER BY seq").fetchall()
        self.ops = [{"seq": r["seq"], "kind": r["kind"], "op": json.loads(r["payload"])} for r in rows]
        self._replay(self.ops, persist=False)

    def close(self) -> None:
        if self.conn is not None:
            self.conn.close()

    # -- replay ------------------------------------------------------------
    def _replay(self, ops: Sequence[Mapping], persist: bool) -> None:
        sizes: dict[str, int] = {}
        joint: dict[tuple[tuple[str, int], ...], Fraction] = {}
        for record in ops:
            seq = record["seq"]
            op = record["op"]
            kind = op.get("op")
            if kind != record.get("kind"):
                raise UnitError(f"seq {seq}: stored kind does not match record payload")
            try:
                self._apply_one(seq, op, sizes, joint)
            except UnitContradiction as exc:
                message = str(exc)
                raise UnitContradiction(
                    message if message.startswith(f"seq {seq}:") else f"seq {seq}: {message}"
                ) from None
            except UnitError as exc:
                # Re-raise with the seq prefix if the low-level message lacks it.
                message = str(exc)
                raise UnitError(message if message.startswith(f"seq {seq}:") else f"seq {seq}: {message}") from None
            if len(sizes) > MAX_BOXES:
                raise UnitError(
                    f"seq {seq}: at most {MAX_BOXES} non-empty boxes are supported, got {len(sizes)} "
                    f"({', '.join(sorted(sizes))})"
                )
            if sum(sizes.values()) > MAX_BEADS:
                raise UnitError(f"seq {seq}: at most {MAX_BEADS} beads are supported, got {sum(sizes.values())}")
        total = _sum_prob(joint)
        if joint and total != 1:
            raise UnitError(f"internal error: joint mass is {qstr(total)}, expected 1")
        self.sizes = sizes
        self.joint = joint
        if persist:
            self._persist(ops)

    def _apply_one(self, seq: int, op: Mapping, sizes: dict[str, int], joint: dict[tuple, Fraction]) -> None:
        kind = op.get("op")
        if kind == "beads":
            self._apply_beads(op, sizes, joint)
        elif kind == "transfer":
            self._apply_transfer(op, sizes, joint)
        elif kind == "mix":
            self._apply_mix(op, sizes, joint)
        elif kind == "inspect":
            self._apply_inspect(op, sizes, joint)
        else:  # pragma: no cover - guarded by normalise
            raise UnitError(f"unknown whole-unit operation {kind!r}")

    def _apply_beads(self, op: Mapping, sizes: dict[str, int], joint: dict[tuple, Fraction]) -> None:
        box = op["box"]
        if box in sizes:
            raise UnitError(f"box {box!r} already holds beads; box identity is stable (discard is not modelled)")
        colored = op["colored"]
        total = colored + op["colorless"]
        if sum(sizes.values()) + total > MAX_BEADS:
            raise UnitError(
                f"declaring {total} beads in {box!r} would exceed {MAX_BEADS} beads "
                f"(already {sum(sizes.values())})"
            )
        # The declared colored count is known, so the new box is deterministic.
        pairs = sorted(sizes.items())
        if not joint and not sizes:
            joint[((box, colored),)] = Fraction(1)
        else:
            new_joint: dict[tuple, Fraction] = {}
            for key, p in joint.items():
                merged = dict(key)
                merged[box] = colored
                nk = _key(merged)
                new_joint[nk] = new_joint.get(nk, ZERO) + p
            joint.clear()
            joint.update(new_joint)
        sizes[box] = total

    def _apply_transfer(self, op: Mapping, sizes: dict[str, int], joint: dict[tuple, Fraction]) -> None:
        source = op["source"]
        dest = op["destination"]
        units = op["units"]
        if source not in sizes:
            raise UnitError(f"source box {source!r} does not exist or is empty")
        if dest in sizes:
            raise UnitError(f"destination {dest!r} is not empty")
        ts = sizes[source]
        if units > ts:
            raise UnitError(f"cannot move {units} units out of {source!r} which holds only {ts}")
        new_ts = ts - units
        new_joint: dict[tuple, Fraction] = {}
        for key, p in joint.items():
            current = dict(key)
            cs = current[source]
            for moved in range(0, min(units, cs) + 1):
                weight = hypergeom_pmf(ts, cs, units, moved)
                if weight == 0:
                    continue
                assigned = {cid: cnt for cid, cnt in current.items() if cid != source}
                if new_ts > 0:
                    assigned[source] = cs - moved
                assigned[dest] = moved
                nk = _key(assigned)
                new_joint[nk] = new_joint.get(nk, ZERO) + p * weight
        joint.clear()
        joint.update(new_joint)
        if new_ts > 0:
            sizes[source] = new_ts
        else:
            del sizes[source]
        sizes[dest] = units

    def _apply_mix(self, op: Mapping, sizes: dict[str, int], joint: dict[tuple, Fraction]) -> None:
        sources = [item["box"] for item in op["sources"]]
        dest = op["destination"]
        if dest in sizes:
            raise UnitError(f"destination {dest!r} is not empty")
        for item in op["sources"]:
            box = item["box"]
            if box not in sizes:
                raise UnitError(f"mix source {box!r} does not exist or is empty")
            full = sizes[box]
            if item["units"] is not None and item["units"] != full:
                raise UnitError(
                    f"mix moves whole boxes: {box!r} holds {full} units but {item['units']} were requested"
                )
        source_set = set(sources)
        new_joint: dict[tuple, Fraction] = {}
        for key, p in joint.items():
            current = dict(key)
            moved = sum(current[box] for box in sources)
            assigned = {cid: cnt for cid, cnt in current.items() if cid not in source_set}
            assigned[dest] = moved
            nk = _key(assigned)
            new_joint[nk] = new_joint.get(nk, ZERO) + p
        joint.clear()
        joint.update(new_joint)
        total = 0
        for box in sources:
            total += sizes.pop(box)
        sizes[dest] = total

    def _apply_inspect(self, op: Mapping, sizes: dict[str, int], joint: dict[tuple, Fraction]) -> None:
        box = op["box"]
        lo, hi = op["min"], op["max"]
        if box not in sizes:
            raise UnitContradiction(f"inspect refers to box {box!r} which does not exist or is empty")
        size = sizes[box]
        if lo > size or hi < 0:
            raise UnitContradiction(
                f"inspect {box!r} ∈ [{lo}, {hi}] is impossible: box holds only {size} beads"
            )
        kept: dict[tuple, Fraction] = {}
        for key, p in joint.items():
            c = dict(key)[box]
            if lo <= c <= hi:
                kept[key] = kept.get(key, ZERO) + p
        total = _sum_prob(kept)
        if total == 0:
            raise UnitContradiction(
                f"inspect {box!r} ∈ [{lo}, {hi}] has zero conditional probability "
                f"given the archive and earlier inspections"
            )
        joint.clear()
        joint.update({key: p / total for key, p in kept.items()})

    # -- batch API ---------------------------------------------------------
    def _persist(self, ops: Sequence[Mapping]) -> None:
        if self.conn is None:
            return
        with self.conn:
            self.conn.execute("DELETE FROM unit_ops")
            for record in ops:
                self.conn.execute(
                    "INSERT INTO unit_ops(seq, kind, payload) VALUES (?,?,?)",
                    (record["seq"], record["kind"], json.dumps(record["op"], sort_keys=True)),
                )

    def apply_batch(self, ops: Sequence[Mapping]) -> list[int]:
        """Validate and append a batch of records atomically; contradictions abort it."""
        if not ops:
            return []
        candidate = [dict(o) for o in self.ops]
        next_seq = len(candidate) + 1
        for raw in ops:
            normalised = _normalise(raw, next_seq)
            candidate.append({"seq": next_seq, "kind": normalised["op"], "op": normalised})
            next_seq += 1
        # Full replay validates every rule; any UnitError leaves stored state
        # intact because the replay builds its state locally and only assigns at
        # the very end.  Persistence happens after that assignment.
        self._replay(candidate, persist=True)
        self.ops = candidate
        return [record["seq"] for record in candidate[len(candidate) - len(ops):]]

    # -- queries -----------------------------------------------------------
    def latest_seq(self) -> int:
        return len(self.ops)

    def boxes(self) -> list[dict]:
        return [{"box": cid, "units": self.sizes[cid]} for cid in sorted(self.sizes)]

    def _require(self, box: str) -> int:
        if box not in self.sizes:
            raise UnitError(f"box {box!r} does not exist or is empty")
        return self.sizes[box]

    def marginal(self, box: str) -> dict[int, Fraction]:
        self._require(box)
        out: dict[int, Fraction] = {}
        for key, p in self.joint.items():
            out[dict(key)[box]] = out.get(dict(key)[box], ZERO) + p
        return dict(sorted(out.items()))

    def joint_distribution(self, first: str, second: str) -> dict[tuple[int, int], Fraction]:
        """Exact joint law of two *different* current boxes' colored counts.

        Marginals alone would lose the anti-correlation between boxes split from a
        shared parent; this keeps it.
        """
        if first == second:
            raise UnitError("joint distribution needs two different boxes; use marginal() for one")
        self._require(first)
        self._require(second)
        out: dict[tuple[int, int], Fraction] = {}
        for key, p in self.joint.items():
            current = dict(key)
            pair = (current[first], current[second])
            out[pair] = out.get(pair, ZERO) + p
        return dict(sorted(out.items()))

    def joint_product(self, first: str, second: str | None = None) -> dict:
        """Self-contained payload a contract inspector may read.

        Contains only the involved boxes' sizes and the exact joint (or marginal)
        law, so the inspector needs no access to the archive's private log.
        """
        self._require(first)
        if second is None or second == first:
            counts = sorted(self.marginal(first).items())
            boxes = [{"box": first, "units": self.sizes[first]}]
            joint_rows = [{"counts": [c], "prob": qstr(p)} for c, p in counts]
        else:
            self._require(second)
            joint = self.joint_distribution(first, second)
            boxes = [
                {"box": first, "units": self.sizes[first]},
                {"box": second, "units": self.sizes[second]},
            ]
            joint_rows = [{"counts": [a, b], "prob": qstr(p)} for (a, b), p in joint.items()]
        return {"boxes": boxes, "joint": joint_rows}

    def log(self) -> list[dict]:
        return [{"seq": r["seq"], "kind": r["kind"], "op": r["op"]} for r in self.ops]


def beads_from_ledger(
    ledger,
    initial: Mapping[str, tuple[int, int]],
    upto_seq: int | None = None,
) -> list[dict]:
    """Translate an existing dilution ledger's real records into bead records.

    ``initial`` maps each *parent* bead box to ``(colored, colorless)`` -- the two
    known counts the ledger itself never stores.  The function then walks the
    ledger log in order and, for every ``transfer`` / ``mix`` that moves a bead
    box, emits the matching whole-unit record (the volume is read as a whole
    number of units).  A destination of such a transfer becomes a bead box too.

    Only genuine ledger operations drive the result; nothing is invented.  A bead
    box touched by a non-whole-unit operation (``solvent``/``evaporate``/``discard``)
    or by a transfer/mix that also involves a non-bead container raises
    :class:`UnitError`, so the discrete archive can never silently diverge from the
    dilution ledger it mirrors.
    """
    records = [{"op": "beads", "box": box, "colored": int(c), "colorless": int(w)} for box, (c, w) in initial.items()]
    bead_set = set(initial)
    for record in ledger.log():
        seq = record["seq"]
        if upto_seq is not None and seq > upto_seq:
            break
        op = record["op"]
        kind = record["kind"]
        if kind == "transfer":
            src, dst = op["source"], op["destination"]
            if src not in bead_set and dst not in bead_set:
                continue
            if src not in bead_set:
                raise UnitError(
                    f"seq {seq}: {dst!r} is a bead box but the transfer source {src!r} is not"
                )
            bead_set.add(dst)
            records.append(
                {"op": "transfer", "source": src, "units": op["volume"], "destination": dst}
            )
        elif kind == "mix":
            ids = [s["container"] for s in op["sources"]]
            dst = op["destination"]
            if not any(cid in bead_set for cid in ids) and dst not in bead_set:
                continue
            outsiders = [cid for cid in ids if cid not in bead_set]
            if outsiders:
                raise UnitError(
                    f"seq {seq}: mix into a bead box pulls in non-bead containers {outsiders}"
                )
            bead_set.add(dst)
            records.append(
                {
                    "op": "mix",
                    "sources": [{"box": s["container"], "units": s["volume"]} for s in op["sources"]],
                    "destination": dst,
                }
            )
        elif kind in ("solvent", "evaporate", "discard"):
            cid = op.get("container") or op.get("source")
            if cid in bead_set:
                raise UnitError(
                    f"seq {seq}: {kind} on bead box {cid!r} has no whole-unit meaning"
                )
        # standard/sample on non-bead containers is irrelevant to the archive.
    return records
