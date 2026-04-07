"""
two_node_db.py — SQLite database + CLI for two-node quiver universality classes.

Build the database once (takes ~3 min), then query instantly.

Usage:
  python two_node_db.py build                   # scan & store all theories
  python two_node_db.py classes                 # list all universality classes
  python two_node_db.py classes --pair SU-SU    # filter by gauge pair
  python two_node_db.py show 7                  # show all theories in class 7
  python two_node_db.py search --matter0 adj    # search theories
  python two_node_db.py stats                   # database summary
"""

from __future__ import annotations

import argparse
import signal
import sqlite3
import sys
import time
from datetime import datetime
from fractions import Fraction

from quiver_generation import (
    Quiver, enumerate_quivers, enumerate_quivers_mixed_rank,
    chiral_excess_coeffs, nf_bound, nf_bound_str,
)
from a_maximization_large_N import (
    a_maximize_large_N_fast, a_maximize_large_N_fast_full,
    a_maximize_large_N, build_fields_large_N,
    _fmt_matter, _fmt_edges, _fmt_linear, _fmt_expr, _fmt_R,
)

DEFAULT_DB = "quivers.db"
TOL = 1e-5


# ── Schema ────────────────────────────────────────────────────────────────────

_DDL = """
CREATE TABLE IF NOT EXISTS universality_class (
    class_id       INTEGER PRIMARY KEY AUTOINCREMENT,
    gauge_pair     TEXT    NOT NULL,
    rank0_mult     INTEGER NOT NULL DEFAULT 1,
    rank1_mult     INTEGER NOT NULL DEFAULT 1,
    a_over_N2      REAL    NOT NULL,
    n_theories     INTEGER NOT NULL,
    rep_theory_id  INTEGER,
    a_exact        TEXT,
    c_exact        TEXT,
    R_exact        TEXT
);

CREATE TABLE IF NOT EXISTS theory (
    theory_id      INTEGER PRIMARY KEY AUTOINCREMENT,
    class_id       INTEGER NOT NULL,
    gauge_pair     TEXT    NOT NULL,
    gauge0         TEXT    NOT NULL,
    gauge1         TEXT    NOT NULL,
    rank0_mult     INTEGER NOT NULL DEFAULT 1,
    rank1_mult     INTEGER NOT NULL DEFAULT 1,
    matter0        TEXT    NOT NULL,
    matter1        TEXT    NOT NULL,
    edges          TEXT    NOT NULL,
    delta0         TEXT,
    delta1         TEXT,
    delta0_a       INTEGER,
    delta0_b       INTEGER,
    delta1_a       INTEGER,
    delta1_b       INTEGER,
    nf_bound0      TEXT    NOT NULL,
    nf_bound1      TEXT    NOT NULL,
    a_over_N2      REAL    NOT NULL,
    c_over_N2      REAL,
    R_numerical    TEXT,
    FOREIGN KEY (class_id) REFERENCES universality_class(class_id)
);

CREATE TABLE IF NOT EXISTS build_info (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"""

_INDEXES = """
CREATE INDEX IF NOT EXISTS idx_theory_class  ON theory(class_id);
CREATE INDEX IF NOT EXISTS idx_theory_gpair  ON theory(gauge_pair);
CREATE INDEX IF NOT EXISTS idx_theory_a      ON theory(a_over_N2);
CREATE INDEX IF NOT EXISTS idx_class_gpair_a ON universality_class(gauge_pair, a_over_N2);
"""

GAUGE_PAIR_LABEL = {
    "SU-SU": "SU(N)×SU(N)",
    "SU-SO": "SU(N)×SO(N)",
    "SO-SU": "SO(N)×SU(N)",
    "SU-Sp": "SU(N)×Sp(N)",
    "Sp-SU": "Sp(N)×SU(N)",
    "SO-SO": "SO(N)×SO(N)",
    "SO-Sp": "SO(N)×Sp(N)",
    "Sp-SO": "Sp(N)×SO(N)",
    "Sp-Sp": "Sp(N)×Sp(N)",
}

PAIR_ORDER = ["SU-SU", "SU-SO", "SO-SU", "SU-Sp", "Sp-SU",
              "SO-SO", "SO-Sp", "Sp-SO", "Sp-Sp"]


def _gauge_pair_label(g0: str, g1: str, m0: int = 1, m1: int = 1) -> str:
    """Human-readable gauge pair label, e.g. 'SU(2N)×SU(N)'."""
    def _grp(g: str, m: int) -> str:
        rank = f"{m}N" if m > 1 else "N"
        return f"{g}({rank})"
    return f"{_grp(g0, m0)}×{_grp(g1, m1)}"


# ── Helpers ───────────────────────────────────────────────────────────────────

def _gauge_pair_key(g0: str, g1: str) -> str:
    return f"{g0}-{g1}"


def _delta_for_node(q: Quiver, node: int) -> tuple[str | None, int | None, int | None]:
    """Return (display_str, a_coeff, b_coeff) for chiral excess at a node."""
    if q.gauge_types[node] != "SU":
        return None, None, None
    a, b = chiral_excess_coeffs(q, node)
    return _fmt_linear(a, b), a, b


def _nf_bound_str_for(q: Quiver, node: int) -> str:
    """Format the N_f AF bound for a node."""
    res = nf_bound(q, node)
    if res is None:
        return "—"
    alpha, gamma = res
    return nf_bound_str(alpha, gamma)


def _field_signature(fields) -> tuple:
    """Hashable signature of a quiver's leading-order field structure."""
    return tuple(sorted(
        (tuple(sorted(f.T_lead.items())), f.dim_lead)
        for f in fields
    ))


class _Timeout(Exception):
    pass


def _exact_with_timeout(quiver: Quiver, timeout: int = 30):
    """Run exact a-maximization with a timeout. Returns LargeNResult or None."""
    def _handler(signum, frame):
        raise _Timeout()
    old = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(timeout)
    try:
        result = a_maximize_large_N(quiver)
        signal.alarm(0)
        return result
    except (_Timeout, Exception):
        signal.alarm(0)
        return None
    finally:
        signal.signal(signal.SIGALRM, old)


def _format_R_charges(result, quiver) -> str:
    """Format exact R-charges as a compact display string."""
    parts = []
    for label, R in result.R_charges.items():
        short = _fmt_R(label)
        # Add node index subscript for 2-node theories
        if label.startswith("node"):
            node_idx = label.split("_")[0].replace("node", "")
            rep_name = "_".join(label.split("_")[1:])
            short = f"{_fmt_R('node0_' + rep_name)}_{node_idx}"
        elif label.startswith("edge"):
            short = "bif"
        parts.append(f"R_{short}={_fmt_expr(R)}")
    return ",  ".join(parts) if parts else "—"


def _cluster(a_vals: list[float]) -> list[tuple[float, list[int]]]:
    """
    Cluster sorted float values into groups within TOL.
    Returns list of (centroid, [indices_into_a_vals]) sorted by centroid.
    """
    if not a_vals:
        return []
    import numpy as np
    order = sorted(range(len(a_vals)), key=lambda i: a_vals[i])
    clusters: list[tuple[float, list[int]]] = []
    cur_center = a_vals[order[0]]
    cur_members = [order[0]]
    for idx in order[1:]:
        if abs(a_vals[idx] - cur_center) < TOL:
            cur_members.append(idx)
            cur_center = sum(a_vals[i] for i in cur_members) / len(cur_members)
        else:
            clusters.append((cur_center, cur_members))
            cur_center = a_vals[idx]
            cur_members = [idx]
    clusters.append((cur_center, cur_members))
    return sorted(clusters, key=lambda x: x[0])


def _rep_idx(rows: list[dict], members: list[int]) -> int:
    """Pick the representative: fewest edges, then alphabetical on matter."""
    return min(
        members,
        key=lambda i: (len(rows[i]["edges"]), rows[i]["matter0"], rows[i]["matter1"]),
    )


# ── Build ─────────────────────────────────────────────────────────────────────

def _format_R_numerical(fast_result) -> str | None:
    """Format fast numerical R-charges as a compact string for DB storage."""
    if fast_result is None or not fast_result.R_charges:
        return None
    parts = []
    for label, r in fast_result.R_charges.items():
        short = _fmt_R(label)
        parts.append(f"R_{short}={r:.5f}")
    return ",  ".join(parts)


def cmd_build(args: argparse.Namespace) -> None:
    import os
    db_path = args.db
    max_a = args.max_a
    exact_timeout = args.exact_timeout

    if os.path.exists(db_path):
        if not args.force:
            ans = input(f"'{db_path}' already exists. Overwrite? [y/N] ").strip().lower()
            if ans != "y":
                print("Aborted.")
                return
        os.remove(db_path)

    # ── Phase 1: enumerate and fast numerical scan ────────────────────────────
    print("Phase 1: enumerating 2-node quivers (equal + mixed rank)...", flush=True)
    t0 = time.time()

    # Equal-rank quivers (all nodes rank N)
    equal_rank = enumerate_quivers(
        n_nodes=2, max_multiedge=4, min_multiedge=1, require_connected=True,
    )
    # Mixed-rank quivers: [2,1] and [1,2]
    mixed_21 = enumerate_quivers_mixed_rank(
        n_nodes=2, rank_multipliers=[2, 1],
        max_multiedge=4, min_multiedge=1, require_connected=True,
    )
    mixed_12 = enumerate_quivers_mixed_rank(
        n_nodes=2, rank_multipliers=[1, 2],
        max_multiedge=4, min_multiedge=1, require_connected=True,
    )
    all_quivers = (
        [(q, bounds, (1, 1)) for q, bounds in equal_rank] +
        [(q, bounds, (2, 1)) for q, bounds in mixed_21] +
        [(q, bounds, (1, 2)) for q, bounds in mixed_12]
    )
    total = len(all_quivers)
    print(f"  {len(equal_rank)} equal-rank + {len(mixed_21)} [2,1] + {len(mixed_12)} [1,2] "
          f"= {total} quivers found in {time.time()-t0:.1f}s", flush=True)

    print("Phase 1: fast a-maximization scan...", flush=True)
    t1 = time.time()

    from collections import defaultdict
    # Bucket key: (gauge_pair, m0, m1) to keep equal-rank and mixed-rank separate
    buckets: dict[tuple, list[dict]] = defaultdict(list)
    n_none = 0
    n_filtered = 0

    for i, (q, bounds, (m0, m1)) in enumerate(all_quivers):
        if i % 200 == 0:
            print(f"  {i}/{total}\r", end="", flush=True)
        fast_res = a_maximize_large_N_fast_full(q)
        if fast_res is None:
            n_none += 1
            continue
        a_val = fast_res.a_over_N2
        if abs(a_val) > max_a:
            n_filtered += 1
            continue

        g0, g1 = q.gauge_types[0], q.gauge_types[1]
        pair = _gauge_pair_key(g0, g1)
        bucket_key = (pair, m0, m1)

        d0_str, d0_a, d0_b = _delta_for_node(q, 0)
        d1_str, d1_a, d1_b = _delta_for_node(q, 1)
        nf0 = _nf_bound_str_for(q, 0)
        nf1 = _nf_bound_str_for(q, 1)

        fields = build_fields_large_N(q)

        buckets[bucket_key].append({
            "gauge_pair": pair, "gauge0": g0, "gauge1": g1,
            "rank0_mult": m0, "rank1_mult": m1,
            "matter0": _fmt_matter(q.node_matter[0]),
            "matter1": _fmt_matter(q.node_matter[1]),
            "edges": _fmt_edges(q),
            "delta0": d0_str, "delta1": d1_str,
            "delta0_a": d0_a, "delta0_b": d0_b,
            "delta1_a": d1_a, "delta1_b": d1_b,
            "nf_bound0": nf0, "nf_bound1": nf1,
            "a_over_N2": a_val,
            "c_over_N2": fast_res.c_over_N2,
            "R_numerical": _format_R_numerical(fast_res),
            "quiver": q,
            "n_fields": len(fields),
        })

    n_stored = sum(len(v) for v in buckets.values())
    print(f"  Done in {time.time()-t1:.1f}s. "
          f"{n_none} diverged, {n_filtered} filtered (|a|>{max_a}), "
          f"{n_stored} stored.", flush=True)

    # ── Phase 2: cluster into classes ─────────────────────────────────────────
    class_list: list[dict] = []

    for (pair, m0, m1), rows in buckets.items():
        if not rows:
            continue
        a_vals = [r["a_over_N2"] for r in rows]
        clusters = _cluster(a_vals)
        for centroid, members in clusters:
            class_list.append({
                "pair": pair, "m0": m0, "m1": m1,
                "centroid": centroid,
                "members": members,
                "rows": rows,
            })

    total_classes = len(class_list)
    print(f"  {total_classes} universality classes identified.", flush=True)

    # ── Phase 3: exact symbolic a-maximization per class ──────────────────────
    print(f"Phase 3: exact a-maximization ({total_classes} classes, "
          f"timeout={exact_timeout}s each)...", flush=True)
    t2 = time.time()
    n_exact_ok = 0
    n_exact_fail = 0

    for ci, cls in enumerate(class_list):
        rows = cls["rows"]
        members = cls["members"]

        best_idx = min(members, key=lambda i: rows[i]["n_fields"])
        rep_q = rows[best_idx]["quiver"]
        n_f = rows[best_idx]["n_fields"]

        print(f"  [{ci+1}/{total_classes}] {cls['pair']}({cls['m0']},{cls['m1']}) "
              f"a≈{cls['centroid']:.6f} ({n_f} fields)\r", end="", flush=True)

        result = _exact_with_timeout(rep_q, timeout=exact_timeout)
        if result is not None:
            try:
                exact_float = float(result.a_over_N2)
            except (TypeError, ValueError):
                exact_float = None
            if exact_float is not None and abs(exact_float - cls["centroid"]) < 0.001:
                cls["a_exact"] = _fmt_expr(result.a_over_N2)
                cls["c_exact"] = _fmt_expr(result.c_over_N2)
                cls["R_exact"] = _format_R_charges(result, rep_q)
                n_exact_ok += 1
            else:
                cls["a_exact"] = None
                cls["c_exact"] = None
                cls["R_exact"] = None
                n_exact_fail += 1
        else:
            cls["a_exact"] = None
            cls["c_exact"] = None
            cls["R_exact"] = None
            n_exact_fail += 1

    print(f"\n  Done in {time.time()-t2:.1f}s. "
          f"{n_exact_ok} exact, {n_exact_fail} failed.", flush=True)

    # ── Phase 4: write to DB ──────────────────────────────────────────────────
    con = sqlite3.connect(db_path)
    con.executescript(_DDL)
    con.executescript(_INDEXES)

    total_theories = 0
    for ci, cls in enumerate(class_list):
        cid = ci + 1
        rows = cls["rows"]
        members = cls["members"]
        m0, m1 = cls["m0"], cls["m1"]

        with con:
            con.execute(
                "INSERT INTO universality_class "
                "(class_id, gauge_pair, rank0_mult, rank1_mult, "
                " a_over_N2, n_theories, a_exact, c_exact, R_exact) "
                "VALUES (?,?,?,?,?,?,?,?,?)",
                (cid, cls["pair"], m0, m1, cls["centroid"], len(members),
                 cls["a_exact"], cls["c_exact"], cls["R_exact"]),
            )

            theory_ids = []
            for idx in members:
                r = rows[idx]
                cur = con.execute(
                    "INSERT INTO theory "
                    "(class_id, gauge_pair, gauge0, gauge1, rank0_mult, rank1_mult, "
                    " matter0, matter1, edges, delta0, delta1, "
                    " delta0_a, delta0_b, delta1_a, delta1_b, "
                    " nf_bound0, nf_bound1, a_over_N2, c_over_N2, R_numerical) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (cid, r["gauge_pair"], r["gauge0"], r["gauge1"],
                     r["rank0_mult"], r["rank1_mult"],
                     r["matter0"], r["matter1"], r["edges"],
                     r["delta0"], r["delta1"],
                     r["delta0_a"], r["delta0_b"],
                     r["delta1_a"], r["delta1_b"],
                     r["nf_bound0"], r["nf_bound1"],
                     r["a_over_N2"], r["c_over_N2"], r["R_numerical"]),
                )
                theory_ids.append(cur.lastrowid)

            rep_local = _rep_idx(rows, members)
            rep_tid = theory_ids[members.index(rep_local)]
            con.execute(
                "UPDATE universality_class SET rep_theory_id=? WHERE class_id=?",
                (rep_tid, cid),
            )
            total_theories += len(members)

    with con:
        con.executemany(
            "INSERT OR REPLACE INTO build_info (key, value) VALUES (?,?)",
            [
                ("built_at",      datetime.now().isoformat(timespec="seconds")),
                ("n_theories",    str(total_theories)),
                ("n_classes",     str(total_classes)),
                ("n_exact",       str(n_exact_ok)),
                ("n_exact_fail",  str(n_exact_fail)),
                ("max_a_filter",  str(max_a)),
                ("tolerance",     str(TOL)),
                ("max_multiedge", "4"),
                ("exact_timeout", str(exact_timeout)),
            ],
        )

    con.close()
    print(f"\nDatabase written to '{db_path}'.")
    print(f"  {total_theories} theories, {total_classes} classes "
          f"({n_exact_ok} exact, {n_exact_fail} numerical).")


# ── Classes ───────────────────────────────────────────────────────────────────

def cmd_classes(args: argparse.Namespace) -> None:
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row

    pair_filter = args.pair.upper().replace("X", "-") if args.pair else None
    rank_filter = args.rank if hasattr(args, "rank") and args.rank else None

    query = """
        SELECT uc.class_id, uc.gauge_pair, uc.rank0_mult, uc.rank1_mult,
               uc.a_over_N2, uc.n_theories,
               uc.a_exact, uc.R_exact,
               t.matter0, t.matter1, t.edges
        FROM universality_class uc
        LEFT JOIN theory t ON t.theory_id = uc.rep_theory_id
        WHERE (:pair IS NULL OR uc.gauge_pair = :pair)
          AND (:min_a IS NULL OR uc.a_over_N2 >= :min_a)
          AND (:max_a IS NULL OR uc.a_over_N2 <= :max_a)
        ORDER BY uc.gauge_pair, uc.rank0_mult, uc.rank1_mult, uc.a_over_N2
    """
    rows = con.execute(query, {
        "pair":  pair_filter,
        "min_a": args.min_a,
        "max_a": args.max_a,
    }).fetchall()
    con.close()

    if rank_filter:
        m0f, m1f = rank_filter
        rows = [r for r in rows if r["rank0_mult"] == m0f and r["rank1_mult"] == m1f]

    if not rows:
        print("No classes found.")
        return

    # Compute column widths for exact a/N²
    a_exact_w = max(
        (len(r["a_exact"] or f"{r['a_over_N2']:.6f}") for r in rows),
        default=10,
    )
    a_exact_w = max(a_exact_w, 10)

    current_group = None
    total = 0
    for r in rows:
        group_key = (r["gauge_pair"], r["rank0_mult"], r["rank1_mult"])
        if group_key != current_group:
            current_group = group_key
            label = _gauge_pair_label(r["gauge_pair"].split("-")[0],
                                      r["gauge_pair"].split("-")[1],
                                      r["rank0_mult"], r["rank1_mult"])
            print(f"\n  {label}")
            print(f"  {'─'*90}")
            print(f"  {'ID':>5}  {'a/N² (exact)':<{a_exact_w}}  {'≈':>10}  {'#th':>4}  "
                  f"Representative theory")
            print(f"  {'─'*90}")
        m0 = r["matter0"] or "—"
        m1 = r["matter1"] or "—"
        e  = r["edges"] or "—"
        rep = f"({m0})  |  ({m1})  |  {e}"
        a_ex = r["a_exact"] or f"≈ {r['a_over_N2']:.6f}"
        print(f"  {r['class_id']:>5}  {a_ex:<{a_exact_w}}  "
              f"{r['a_over_N2']:>10.6f}  {r['n_theories']:>4}  {rep}")
        total += 1

    print(f"\n  Total: {total} classes")


# ── Show ──────────────────────────────────────────────────────────────────────

def cmd_show(args: argparse.Namespace) -> None:
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row

    cls = con.execute(
        "SELECT * FROM universality_class WHERE class_id=?", (args.class_id,)
    ).fetchone()
    if cls is None:
        print(f"Class {args.class_id} not found.")
        con.close()
        return

    theories = con.execute(
        "SELECT * FROM theory WHERE class_id=? ORDER BY length(edges), matter0, matter1",
        (args.class_id,),
    ).fetchall()
    con.close()

    label = _gauge_pair_label(cls["gauge_pair"].split("-")[0],
                              cls["gauge_pair"].split("-")[1],
                              cls["rank0_mult"], cls["rank1_mult"])
    print(f"\nUniversality class #{cls['class_id']}  "
          f"{label}  "
          f"({cls['n_theories']} theories)")

    a_ex = cls["a_exact"] or f"≈ {cls['a_over_N2']:.8f}  (numerical)"
    c_ex = cls["c_exact"] or "(numerical only)"
    R_ex = cls["R_exact"] or "(numerical only)"
    print(f"  a/N² = {a_ex}  ≈ {cls['a_over_N2']:.8f}")
    print(f"  c/N² = {c_ex}")
    print(f"  R-charges (large N): {R_ex}")
    print("─" * 100)

    # Column widths
    m0w = max(len(t["matter0"] or "—") for t in theories)
    m1w = max(len(t["matter1"] or "—") for t in theories)
    ew  = max(len(t["edges"]   or "—") for t in theories)
    d0w = max(len(t["delta0"]  or "—") for t in theories)
    d1w = max(len(t["delta1"]  or "—") for t in theories)
    m0w = max(m0w, 8);  m1w = max(m1w, 8)
    ew  = max(ew,  5);  d0w = max(d0w, 7);  d1w = max(d1w, 7)

    hdr = (f"  {'ID':>5}  "
           f"{'Matter(0)':<{m0w}}  "
           f"{'Matter(1)':<{m1w}}  "
           f"{'Edges':<{ew}}  "
           f"{'delta(0)':>{d0w}}  "
           f"{'delta(1)':>{d1w}}  "
           f"{'N_f bound(0)':<14}  "
           f"{'N_f bound(1)':<14}  "
           f"{'a/N²':>10}")
    print(hdr)
    print("  " + "─" * (len(hdr) - 2))

    for t in theories:
        m0 = t["matter0"] or "—"
        m1 = t["matter1"] or "—"
        e  = t["edges"]   or "—"
        d0 = t["delta0"]  or "—"
        d1 = t["delta1"]  or "—"
        print(f"  {t['theory_id']:>5}  "
              f"{m0:<{m0w}}  "
              f"{m1:<{m1w}}  "
              f"{e:<{ew}}  "
              f"{d0:>{d0w}}  "
              f"{d1:>{d1w}}  "
              f"{t['nf_bound0']:<14}  "
              f"{t['nf_bound1']:<14}  "
              f"{t['a_over_N2']:>10.6f}")


# ── Search ────────────────────────────────────────────────────────────────────

def cmd_search(args: argparse.Namespace) -> None:
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row

    conditions = []
    params: dict = {}

    if args.pair:
        conditions.append("gauge_pair = :pair")
        params["pair"] = args.pair.upper().replace("X", "-")
    if args.matter0:
        conditions.append("matter0 LIKE :m0")
        params["m0"] = f"%{args.matter0}%"
    if args.matter1:
        conditions.append("matter1 LIKE :m1")
        params["m1"] = f"%{args.matter1}%"
    if args.delta0 is not None:
        conditions.append("delta0_a = :d0a")
        params["d0a"] = args.delta0
    if args.delta1 is not None:
        conditions.append("delta1_a = :d1a")
        params["d1a"] = args.delta1
    if args.min_a is not None:
        conditions.append("a_over_N2 >= :min_a")
        params["min_a"] = args.min_a
    if args.max_a is not None:
        conditions.append("a_over_N2 <= :max_a")
        params["max_a"] = args.max_a

    where = ("WHERE " + " AND ".join(conditions)) if conditions else ""
    limit_clause = f"LIMIT {args.limit}" if args.limit else ""

    rows = con.execute(
        f"SELECT * FROM theory {where} ORDER BY gauge_pair, a_over_N2, matter0 {limit_clause}",
        params,
    ).fetchall()
    con.close()

    if not rows:
        print("No theories found.")
        return

    m0w = max((len(r["matter0"] or "—") for r in rows), default=8)
    m1w = max((len(r["matter1"] or "—") for r in rows), default=8)
    ew  = max((len(r["edges"]   or "—") for r in rows), default=5)
    d0w = max((len(r["delta0"]  or "—") for r in rows), default=7)
    d1w = max((len(r["delta1"]  or "—") for r in rows), default=7)
    m0w = max(m0w, 9);  m1w = max(m1w, 9)

    hdr = (f"  {'ID':>5}  {'Cl':>4}  {'Pair':<7}  "
           f"{'Matter(0)':<{m0w}}  "
           f"{'Matter(1)':<{m1w}}  "
           f"{'Edges':<{ew}}  "
           f"{'delta(0)':>{d0w}}  "
           f"{'delta(1)':>{d1w}}  "
           f"{'a/N²':>10}")
    print(hdr)
    print("  " + "─" * (len(hdr) - 2))

    for r in rows:
        m0 = r["matter0"] or "—"
        m1 = r["matter1"] or "—"
        e  = r["edges"]   or "—"
        d0 = r["delta0"]  or "—"
        d1 = r["delta1"]  or "—"
        print(f"  {r['theory_id']:>5}  {r['class_id']:>4}  {r['gauge_pair']:<7}  "
              f"{m0:<{m0w}}  "
              f"{m1:<{m1w}}  "
              f"{e:<{ew}}  "
              f"{d0:>{d0w}}  "
              f"{d1:>{d1w}}  "
              f"{r['a_over_N2']:>10.6f}")

    print(f"\n  {len(rows)} theories found.")


# ── Stats ─────────────────────────────────────────────────────────────────────

def cmd_stats(args: argparse.Namespace) -> None:
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row

    info = {r["key"]: r["value"] for r in con.execute("SELECT * FROM build_info")}
    print("\nBuild info:")
    for k, v in info.items():
        print(f"  {k:<16}: {v}")

    print("\nTheories and classes by gauge pair and rank:")
    print(f"  {'Gauge pair':<22}  {'#theories':>9}  {'#classes':>8}")
    print(f"  {'─'*44}")
    rows_stat = con.execute(
        "SELECT gauge_pair, rank0_mult, rank1_mult, "
        "COUNT(*) n_th FROM theory GROUP BY gauge_pair, rank0_mult, rank1_mult"
    ).fetchall()
    cls_stat = con.execute(
        "SELECT gauge_pair, rank0_mult, rank1_mult, "
        "COUNT(*) n_cl FROM universality_class GROUP BY gauge_pair, rank0_mult, rank1_mult"
    ).fetchall()
    cls_map = {(r["gauge_pair"], r["rank0_mult"], r["rank1_mult"]): r["n_cl"]
               for r in cls_stat}
    for r in sorted(rows_stat, key=lambda x: (PAIR_ORDER.index(x["gauge_pair"])
                                               if x["gauge_pair"] in PAIR_ORDER else 99,
                                               x["rank0_mult"], x["rank1_mult"])):
        label = _gauge_pair_label(r["gauge_pair"].split("-")[0],
                                  r["gauge_pair"].split("-")[1],
                                  r["rank0_mult"], r["rank1_mult"])
        n_cl = cls_map.get((r["gauge_pair"], r["rank0_mult"], r["rank1_mult"]), 0)
        print(f"  {label:<22}  {r['n_th']:>9}  {n_cl:>8}")

    a_stats = con.execute(
        "SELECT MIN(a_over_N2) mn, MAX(a_over_N2) mx, AVG(a_over_N2) av "
        "FROM universality_class"
    ).fetchone()
    print(f"\na/N² distribution (all classes):")
    print(f"  min = {a_stats['mn']:.6f}   max = {a_stats['mx']:.6f}   mean = {a_stats['av']:.6f}")

    con.close()


# ── CLI entry point ───────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Browse two-node quiver universality classes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--db", default=DEFAULT_DB, help="SQLite database path")
    sub = parser.add_subparsers(dest="cmd", required=True)

    # build
    p = sub.add_parser("build", help="Build/rebuild the database from scratch")
    p.add_argument("--max-a", type=float, default=2.0,
                   help="Exclude theories with |a/N²| > this (default 2.0)")
    p.add_argument("--force", action="store_true", help="Overwrite existing DB without asking")
    p.add_argument("--exact-timeout", type=int, default=30, dest="exact_timeout",
                   help="Timeout in seconds per exact solve (default 30)")

    def _parse_rank(s):
        parts = s.split(",")
        if len(parts) != 2:
            raise argparse.ArgumentTypeError("--rank expects 'm0,m1' e.g. '2,1'")
        return (int(parts[0]), int(parts[1]))

    # classes
    p = sub.add_parser("classes", help="List universality classes")
    p.add_argument("--pair", help="Filter by gauge pair, e.g. SU-SU or SU×SU")
    p.add_argument("--min-a", type=float, dest="min_a")
    p.add_argument("--max-a", type=float, dest="max_a")
    p.add_argument("--rank", type=_parse_rank,
                   help="Filter by rank multipliers, e.g. '2,1' for SU(2N)×SU(N)")

    # show
    p = sub.add_parser("show", help="Show all theories in a universality class")
    p.add_argument("class_id", type=int, help="Class ID (from 'classes' output)")

    # search
    p = sub.add_parser("search", help="Filter theories by properties")
    p.add_argument("--pair",    help="Gauge pair, e.g. SU-SU")
    p.add_argument("--matter0", help="Substring match on matter at node 0")
    p.add_argument("--matter1", help="Substring match on matter at node 1")
    p.add_argument("--delta0",  type=int,
                   help="Exact match on the N-coefficient of delta at node 0")
    p.add_argument("--delta1",  type=int,
                   help="Exact match on the N-coefficient of delta at node 1")
    p.add_argument("--min-a", type=float, dest="min_a")
    p.add_argument("--max-a", type=float, dest="max_a")
    p.add_argument("--limit", type=int, help="Maximum number of results")

    # stats
    sub.add_parser("stats", help="Database summary and build metadata")

    args = parser.parse_args()

    if args.cmd != "build":
        try:
            open(args.db).close()
        except FileNotFoundError:
            print(f"Database '{args.db}' not found. Run 'build' first.")
            sys.exit(1)

    dispatch = {
        "build":   cmd_build,
        "classes": cmd_classes,
        "show":    cmd_show,
        "search":  cmd_search,
        "stats":   cmd_stats,
    }
    dispatch[args.cmd](args)


if __name__ == "__main__":
    main()
