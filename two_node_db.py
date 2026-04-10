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
import re
import signal
import sqlite3
import sys
import time
from datetime import datetime
from fractions import Fraction

from quiver_generation import (
    Quiver, enumerate_quivers, enumerate_quivers_mixed_rank,
    chiral_excess_coeffs, nf_bound, nf_bound_str,
    _dedup_symmetries,
)
from a_maximization_large_N import (
    a_maximize_large_N_fast, a_maximize_large_N_fast_full,
    a_maximize_large_N, a_maximize_batch_mathematica,
    build_fields_large_N,
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
    R_exact        TEXT,
    a_over_c       REAL,
    veneziano_any  INTEGER NOT NULL DEFAULT 0,
    veneziano_all  INTEGER NOT NULL DEFAULT 0
);

CREATE TABLE IF NOT EXISTS morphology_class (
    morph_id   INTEGER PRIMARY KEY AUTOINCREMENT,
    N_rank2_0  REAL    NOT NULL,
    N_rank2_1  REAL    NOT NULL,
    N_bif      INTEGER NOT NULL,
    N_fund_0   INTEGER NOT NULL,
    N_fund_1   INTEGER NOT NULL,
    n_theories INTEGER NOT NULL DEFAULT 0,
    n_classes  INTEGER NOT NULL DEFAULT 0,
    UNIQUE(N_rank2_0, N_rank2_1, N_bif, N_fund_0, N_fund_1)
);

CREATE TABLE IF NOT EXISTS theory (
    theory_id      INTEGER PRIMARY KEY AUTOINCREMENT,
    class_id       INTEGER,
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
    a_over_N2      REAL,
    c_over_N2      REAL,
    R_numerical    TEXT,
    a_over_c       REAL,
    veneziano      INTEGER NOT NULL DEFAULT 0,
    N_rank2_0      REAL,
    N_rank2_1      REAL,
    N_bif          INTEGER,
    N_fund_0       INTEGER,
    N_fund_1       INTEGER,
    morph_id       INTEGER REFERENCES morphology_class(morph_id),
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
CREATE INDEX IF NOT EXISTS idx_theory_morph  ON theory(morph_id);
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


# ── Morphology helpers ────────────────────────────────────────────────────────

_RANK2_WEIGHT: dict[str, float] = {
    "adj": 1.0, "S": 0.5, "Sbar": 0.5, "A": 0.5, "Abar": 0.5,
    # display aliases used in matter text column
    "S̄": 0.5, "Ā": 0.5,
}


def _parse_N_rank2(matter_str: str) -> float:
    """Parse matter text (e.g. '2adj + S̄') → effective rank-2 tensor count."""
    if matter_str == "—":
        return 0.0
    total = 0.0
    for part in matter_str.split(" + "):
        part = part.strip()
        m = re.match(r"^(\d+)?(.*)", part)
        count = int(m.group(1)) if m and m.group(1) else 1
        rep = m.group(2).strip() if m else part
        total += count * _RANK2_WEIGHT.get(rep, 0.0)
    return total


def _parse_N_bif(edges_str: str) -> int:
    """Parse edges text (e.g. '2×(0→1,+-)') → total bifundamental count."""
    if edges_str == "—":
        return 0
    total = 0
    for m in re.finditer(r"(\d+)×\(|\(", edges_str):
        grp = m.group()
        total += int(grp.rstrip("×(")) if "×" in grp else 1
    return total


def _morph_vec_from_quiver(q: "Quiver", d0_a, d1_a) -> tuple:
    """Compute (N_rank2_0, N_rank2_1, N_bif, N_fund_0, N_fund_1) from a Quiver."""
    nm0, nm1 = q.node_matter[0], q.node_matter[1]
    N_rank2_0 = sum(
        cnt * (1.0 if rep == "adj" else 0.5)
        for rep, cnt in nm0.items()
        if rep in _RANK2_WEIGHT
    )
    N_rank2_1 = sum(
        cnt * (1.0 if rep == "adj" else 0.5)
        for rep, cnt in nm1.items()
        if rep in _RANK2_WEIGHT
    )
    return (N_rank2_0, N_rank2_1, len(q.edges), abs(d0_a or 0), abs(d1_a or 0))


def _morph_vec_from_text(matter0: str, matter1: str, edges: str,
                         delta0_a, delta1_a) -> tuple:
    """Compute morphology vector from stored text columns (for migration)."""
    return (
        _parse_N_rank2(matter0),
        _parse_N_rank2(matter1),
        _parse_N_bif(edges),
        abs(delta0_a or 0),
        abs(delta1_a or 0),
    )


def _get_or_create_morph_id(con: sqlite3.Connection, vec: tuple) -> int:
    """INSERT OR IGNORE the morphology vector, then return its morph_id."""
    N_r0, N_r1, N_bif, N_f0, N_f1 = vec
    con.execute(
        "INSERT OR IGNORE INTO morphology_class "
        "(N_rank2_0, N_rank2_1, N_bif, N_fund_0, N_fund_1, n_theories, n_classes) "
        "VALUES (?,?,?,?,?,0,0)",
        (N_r0, N_r1, N_bif, N_f0, N_f1),
    )
    row = con.execute(
        "SELECT morph_id FROM morphology_class "
        "WHERE N_rank2_0=? AND N_rank2_1=? AND N_bif=? AND N_fund_0=? AND N_fund_1=?",
        (N_r0, N_r1, N_bif, N_f0, N_f1),
    ).fetchone()
    return row[0]


def cmd_build(args: argparse.Namespace) -> None:
    import os
    db_path = args.db
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

    # Equal-rank quivers (all nodes rank N) — already deduped by symmetries internally
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
    # Cross-set dedup: [2,1] and [1,2] quivers may be node-swaps of each other
    all_qb = equal_rank + mixed_21 + mixed_12
    all_qb = _dedup_symmetries(all_qb)
    total = len(all_qb)
    print(f"  {len(equal_rank)} equal-rank + {len(mixed_21)} [2,1] + {len(mixed_12)} [1,2] "
          f"→ {total} after cross-set dedup  ({time.time()-t0:.1f}s)", flush=True)

    print("Phase 1: Mathematica NSolve batch scan...", flush=True)
    t1 = time.time()

    mma_results = a_maximize_batch_mathematica(
        [q for q, _ in all_qb], timeout=7200,
    )
    print(f"  Mathematica scan done in {time.time()-t1:.1f}s.", flush=True)

    from collections import defaultdict
    buckets: dict[tuple, list[dict]] = defaultdict(list)
    diverged_rows: list[dict] = []
    n_diverged = 0

    for (q, bounds), fast_res in zip(all_qb, mma_results):
        m0, m1 = q.rank_multipliers[0], q.rank_multipliers[1]
        g0, g1 = q.gauge_types[0], q.gauge_types[1]
        pair = _gauge_pair_key(g0, g1)

        d0_str, d0_a, d0_b = _delta_for_node(q, 0)
        d1_str, d1_a, d1_b = _delta_for_node(q, 1)
        nf0 = _nf_bound_str_for(q, 0)
        nf1 = _nf_bound_str_for(q, 1)
        veneziano = int(
            (d0_a is not None and d0_a != 0) or
            (d1_a is not None and d1_a != 0)
        )

        morph_vec = _morph_vec_from_quiver(q, d0_a, d1_a)

        if fast_res is None:
            n_diverged += 1
            diverged_rows.append({
                "gauge_pair": pair, "gauge0": g0, "gauge1": g1,
                "rank0_mult": m0, "rank1_mult": m1,
                "matter0": _fmt_matter(q.node_matter[0]),
                "matter1": _fmt_matter(q.node_matter[1]),
                "edges": _fmt_edges(q),
                "delta0": d0_str, "delta1": d1_str,
                "delta0_a": d0_a, "delta0_b": d0_b,
                "delta1_a": d1_a, "delta1_b": d1_b,
                "nf_bound0": nf0, "nf_bound1": nf1,
                "a_over_N2": None, "c_over_N2": None,
                "R_numerical": None, "a_over_c": None,
                "veneziano": veneziano,
                "morph_vec": morph_vec,
            })
            continue

        a_val = fast_res.a_over_N2
        c_val = fast_res.c_over_N2
        a_over_c = (a_val / c_val) if (c_val is not None and abs(c_val) > 1e-12) else None

        fields = build_fields_large_N(q)
        buckets[(pair, m0, m1)].append({
            "gauge_pair": pair, "gauge0": g0, "gauge1": g1,
            "rank0_mult": m0, "rank1_mult": m1,
            "matter0": _fmt_matter(q.node_matter[0]),
            "matter1": _fmt_matter(q.node_matter[1]),
            "edges": _fmt_edges(q),
            "delta0": d0_str, "delta1": d1_str,
            "delta0_a": d0_a, "delta0_b": d0_b,
            "delta1_a": d1_a, "delta1_b": d1_b,
            "nf_bound0": nf0, "nf_bound1": nf1,
            "a_over_N2": a_val, "c_over_N2": c_val,
            "R_numerical": _format_R_numerical(fast_res),
            "a_over_c": a_over_c,
            "veneziano": veneziano,
            "quiver": q,
            "n_fields": len(fields),
            "morph_vec": morph_vec,
        })

    n_stored = sum(len(v) for v in buckets.values())
    print(f"  {n_diverged} diverged (stored as class NULL), "
          f"{n_stored} clusterable.", flush=True)

    # ── Phase 2: cluster into classes ─────────────────────────────────────────
    class_list: list[dict] = []

    for (pair, m0, m1), rows in buckets.items():
        if not rows:
            continue
        a_vals = [r["a_over_N2"] for r in rows]
        clusters = _cluster(a_vals)
        for centroid, members in clusters:
            flags = [rows[i]["veneziano"] for i in members]
            rep_local = _rep_idx(rows, members)
            class_list.append({
                "pair": pair, "m0": m0, "m1": m1,
                "centroid": centroid,
                "members": members,
                "rows": rows,
                "veneziano_any": int(any(flags)),
                "veneziano_all": int(all(flags)),
                "a_over_c": rows[rep_local].get("a_over_c"),
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
                " a_over_N2, n_theories, a_exact, c_exact, R_exact, "
                " a_over_c, veneziano_any, veneziano_all) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                (cid, cls["pair"], m0, m1, cls["centroid"], len(members),
                 cls["a_exact"], cls["c_exact"], cls["R_exact"],
                 cls["a_over_c"], cls["veneziano_any"], cls["veneziano_all"]),
            )

            theory_ids = []
            for idx in members:
                r = rows[idx]
                mv = r["morph_vec"]
                mid = _get_or_create_morph_id(con, mv)
                cur = con.execute(
                    "INSERT INTO theory "
                    "(class_id, gauge_pair, gauge0, gauge1, rank0_mult, rank1_mult, "
                    " matter0, matter1, edges, delta0, delta1, "
                    " delta0_a, delta0_b, delta1_a, delta1_b, "
                    " nf_bound0, nf_bound1, a_over_N2, c_over_N2, R_numerical, "
                    " a_over_c, veneziano, "
                    " N_rank2_0, N_rank2_1, N_bif, N_fund_0, N_fund_1, morph_id) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (cid, r["gauge_pair"], r["gauge0"], r["gauge1"],
                     r["rank0_mult"], r["rank1_mult"],
                     r["matter0"], r["matter1"], r["edges"],
                     r["delta0"], r["delta1"],
                     r["delta0_a"], r["delta0_b"],
                     r["delta1_a"], r["delta1_b"],
                     r["nf_bound0"], r["nf_bound1"],
                     r["a_over_N2"], r["c_over_N2"], r["R_numerical"],
                     r["a_over_c"], r["veneziano"],
                     mv[0], mv[1], mv[2], mv[3], mv[4], mid),
                )
                theory_ids.append(cur.lastrowid)

            rep_local = _rep_idx(rows, members)
            rep_tid = theory_ids[members.index(rep_local)]
            con.execute(
                "UPDATE universality_class SET rep_theory_id=? WHERE class_id=?",
                (rep_tid, cid),
            )
            total_theories += len(members)

    # Insert diverged theories (class_id = NULL)
    with con:
        for r in diverged_rows:
            mv = r["morph_vec"]
            mid = _get_or_create_morph_id(con, mv)
            con.execute(
                "INSERT INTO theory "
                "(class_id, gauge_pair, gauge0, gauge1, rank0_mult, rank1_mult, "
                " matter0, matter1, edges, delta0, delta1, "
                " delta0_a, delta0_b, delta1_a, delta1_b, "
                " nf_bound0, nf_bound1, a_over_N2, c_over_N2, R_numerical, "
                " a_over_c, veneziano, "
                " N_rank2_0, N_rank2_1, N_bif, N_fund_0, N_fund_1, morph_id) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (None, r["gauge_pair"], r["gauge0"], r["gauge1"],
                 r["rank0_mult"], r["rank1_mult"],
                 r["matter0"], r["matter1"], r["edges"],
                 r["delta0"], r["delta1"],
                 r["delta0_a"], r["delta0_b"],
                 r["delta1_a"], r["delta1_b"],
                 r["nf_bound0"], r["nf_bound1"],
                 r["a_over_N2"], r["c_over_N2"], r["R_numerical"],
                 r["a_over_c"], r["veneziano"],
                 mv[0], mv[1], mv[2], mv[3], mv[4], mid),
            )

    with con:
        con.executemany(
            "INSERT OR REPLACE INTO build_info (key, value) VALUES (?,?)",
            [
                ("built_at",      datetime.now().isoformat(timespec="seconds")),
                ("n_theories",    str(total_theories)),
                ("n_diverged",    str(n_diverged)),
                ("n_classes",     str(total_classes)),
                ("n_exact",       str(n_exact_ok)),
                ("n_exact_fail",  str(n_exact_fail)),
                ("tolerance",     str(TOL)),
                ("max_multiedge", "4"),
                ("exact_timeout", str(exact_timeout)),
            ],
        )

    # Update morphology_class aggregate counts
    with con:
        con.execute("""
            UPDATE morphology_class SET n_theories = (
                SELECT COUNT(*) FROM theory WHERE theory.morph_id = morphology_class.morph_id
            )
        """)
        con.execute("""
            UPDATE morphology_class SET n_classes = (
                SELECT COUNT(DISTINCT class_id) FROM theory
                WHERE theory.morph_id = morphology_class.morph_id AND class_id IS NOT NULL
            )
        """)

    con.close()
    print(f"\nDatabase written to '{db_path}'.")
    print(f"  {total_theories} theories in {total_classes} classes "
          f"({n_exact_ok} exact, {n_exact_fail} numerical); "
          f"{n_diverged} diverged (class_id NULL).")


# ── Classes ───────────────────────────────────────────────────────────────────

def cmd_classes(args: argparse.Namespace) -> None:
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row

    pair_filter = args.pair.upper().replace("X", "-") if args.pair else None
    rank_filter = args.rank if hasattr(args, "rank") and args.rank else None
    veneziano_filter = getattr(args, "veneziano", None)

    veneziano_cond = ""
    if veneziano_filter is True:
        veneziano_cond = "AND uc.veneziano_all = 1"
    elif veneziano_filter is False:
        veneziano_cond = "AND uc.veneziano_all = 0"

    query = f"""
        SELECT uc.class_id, uc.gauge_pair, uc.rank0_mult, uc.rank1_mult,
               uc.a_over_N2, uc.n_theories,
               uc.a_exact, uc.R_exact, uc.a_over_c,
               uc.veneziano_any, uc.veneziano_all,
               t.matter0, t.matter1, t.edges
        FROM universality_class uc
        LEFT JOIN theory t ON t.theory_id = uc.rep_theory_id
        WHERE (:pair IS NULL OR uc.gauge_pair = :pair)
          AND (:min_a IS NULL OR uc.a_over_N2 >= :min_a)
          AND (:max_a IS NULL OR uc.a_over_N2 <= :max_a)
          {veneziano_cond}
        ORDER BY uc.gauge_pair, uc.rank0_mult, uc.rank1_mult, uc.a_over_N2
    """
    rows = con.execute(query, {
        "pair":  pair_filter,
        "min_a": args.min_a,
        "max_a": args.max_a,
    }).fetchall()

    n_diverged = con.execute(
        "SELECT COUNT(*) FROM theory WHERE class_id IS NULL"
    ).fetchone()[0]
    con.close()

    if rank_filter:
        m0f, m1f = rank_filter
        rows = [r for r in rows if r["rank0_mult"] == m0f and r["rank1_mult"] == m1f]

    if not rows:
        print("No classes found.")
        if n_diverged:
            print(f"  Diverged / no convergence: {n_diverged} theories  (use 'show null')")
        return

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
            print(f"  {'─'*100}")
            print(f"  {'ID':>5}  {'a/N² (exact)':<{a_exact_w}}  {'≈':>10}  "
                  f"{'#th':>4}  {'Ven':>3}  Representative theory")
            print(f"  {'─'*100}")
        m0 = r["matter0"] or "—"
        m1 = r["matter1"] or "—"
        e  = r["edges"] or "—"
        rep = f"({m0})  |  ({m1})  |  {e}"
        a_ex = r["a_exact"] or f"≈ {r['a_over_N2']:.6f}"
        ven = "Y" if r["veneziano_any"] else "N"
        print(f"  {r['class_id']:>5}  {a_ex:<{a_exact_w}}  "
              f"{r['a_over_N2']:>10.6f}  {r['n_theories']:>4}  {ven:>3}  {rep}")
        total += 1

    print(f"\n  Total: {total} classes")
    if n_diverged:
        print(f"  Diverged / no convergence: {n_diverged} theories  (use 'show null')")


# ── Show ──────────────────────────────────────────────────────────────────────

def _print_theory_table(theories) -> None:
    """Print a table of theories (shared by cmd_show and _show_diverged)."""
    if not theories:
        print("  (none)")
        return

    m0w = max(len(t["matter0"] or "—") for t in theories)
    m1w = max(len(t["matter1"] or "—") for t in theories)
    ew  = max(len(t["edges"]   or "—") for t in theories)
    d0w = max(len(t["delta0"]  or "—") for t in theories)
    d1w = max(len(t["delta1"]  or "—") for t in theories)
    rw  = max(len(t["R_numerical"] or "—") for t in theories)
    m0w = max(m0w, 8);  m1w = max(m1w, 8)
    ew  = max(ew,  5);  d0w = max(d0w, 7);  d1w = max(d1w, 7)
    rw  = min(rw, 60)   # cap R_numerical display width

    hdr = (f"  {'ID':>5}  "
           f"{'Matter(0)':<{m0w}}  "
           f"{'Matter(1)':<{m1w}}  "
           f"{'Edges':<{ew}}  "
           f"{'delta(0)':>{d0w}}  "
           f"{'delta(1)':>{d1w}}  "
           f"{'N_f bound(0)':<14}  "
           f"{'N_f bound(1)':<14}  "
           f"{'a/N²':>10}  "
           f"{'c/N²':>10}  "
           f"{'a/c':>8}  "
           f"{'Ven':>3}  "
           f"{'R-charges':<{rw}}")
    print(hdr)
    print("  " + "─" * (len(hdr) - 2))

    for t in theories:
        m0  = t["matter0"]     or "—"
        m1  = t["matter1"]     or "—"
        e   = t["edges"]       or "—"
        d0  = t["delta0"]      or "—"
        d1  = t["delta1"]      or "—"
        a_s = f"{t['a_over_N2']:.6f}" if t["a_over_N2"] is not None else "—"
        c_s = f"{t['c_over_N2']:.6f}" if t["c_over_N2"] is not None else "—"
        ac_s = f"{t['a_over_c']:.5f}" if t["a_over_c"] is not None else "—"
        ven = "Y" if t["veneziano"] else "N"
        r_s = (t["R_numerical"] or "—")[:rw]
        print(f"  {t['theory_id']:>5}  "
              f"{m0:<{m0w}}  "
              f"{m1:<{m1w}}  "
              f"{e:<{ew}}  "
              f"{d0:>{d0w}}  "
              f"{d1:>{d1w}}  "
              f"{t['nf_bound0']:<14}  "
              f"{t['nf_bound1']:<14}  "
              f"{a_s:>10}  "
              f"{c_s:>10}  "
              f"{ac_s:>8}  "
              f"{ven:>3}  "
              f"{r_s:<{rw}}")


def _show_diverged(args: argparse.Namespace) -> None:
    """Show all theories with class_id IS NULL (optimizer diverged)."""
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    theories = con.execute(
        "SELECT * FROM theory WHERE class_id IS NULL "
        "ORDER BY gauge_pair, matter0, matter1"
    ).fetchall()
    con.close()
    print(f"\nDiverged theories (class_id IS NULL): {len(theories)}")
    print("─" * 100)
    _print_theory_table(theories)


def cmd_show(args: argparse.Namespace) -> None:
    if args.class_id is None:
        _show_diverged(args)
        return

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

    a_ex  = cls["a_exact"] or f"≈ {cls['a_over_N2']:.8f}  (numerical)"
    c_ex  = cls["c_exact"] or "(numerical only)"
    R_ex  = cls["R_exact"] or "(numerical only)"
    ac_s  = f"{cls['a_over_c']:.6f}" if cls["a_over_c"] is not None else "(n/a)"
    ven_s = f"any={cls['veneziano_any']} all={cls['veneziano_all']}"
    print(f"  a/N² = {a_ex}  ≈ {cls['a_over_N2']:.8f}")
    print(f"  c/N² = {c_ex}")
    print(f"  a/c  = {ac_s}")
    print(f"  R-charges (large N): {R_ex}")
    print(f"  Veneziano: {ven_s}")
    print("─" * 100)

    _print_theory_table(theories)


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
    veneziano_filter = getattr(args, "veneziano", None)
    if veneziano_filter is True:
        conditions.append("veneziano = 1")
    elif veneziano_filter is False:
        conditions.append("veneziano = 0")

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


# ── Morphology build (migration) ──────────────────────────────────────────────

def cmd_morph_build(args: argparse.Namespace) -> None:
    """Populate morphology columns on an existing DB without recomputing R-charges."""
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row

    # Add new columns to theory (idempotent)
    new_cols = [
        ("N_rank2_0", "REAL"),
        ("N_rank2_1", "REAL"),
        ("N_bif",     "INTEGER"),
        ("N_fund_0",  "INTEGER"),
        ("N_fund_1",  "INTEGER"),
        ("morph_id",  "INTEGER"),
    ]
    for col, typ in new_cols:
        try:
            con.execute(f"ALTER TABLE theory ADD COLUMN {col} {typ}")
        except sqlite3.OperationalError:
            pass  # column already exists

    # Create morphology_class table and index
    con.executescript("""
        CREATE TABLE IF NOT EXISTS morphology_class (
            morph_id   INTEGER PRIMARY KEY AUTOINCREMENT,
            N_rank2_0  REAL    NOT NULL,
            N_rank2_1  REAL    NOT NULL,
            N_bif      INTEGER NOT NULL,
            N_fund_0   INTEGER NOT NULL,
            N_fund_1   INTEGER NOT NULL,
            n_theories INTEGER NOT NULL DEFAULT 0,
            n_classes  INTEGER NOT NULL DEFAULT 0,
            UNIQUE(N_rank2_0, N_rank2_1, N_bif, N_fund_0, N_fund_1)
        );
        CREATE INDEX IF NOT EXISTS idx_theory_morph ON theory(morph_id);
    """)

    rows = con.execute(
        "SELECT theory_id, matter0, matter1, edges, delta0_a, delta1_a FROM theory"
    ).fetchall()
    print(f"Processing {len(rows)} theories...", flush=True)

    with con:
        for row in rows:
            vec = _morph_vec_from_text(
                row["matter0"], row["matter1"], row["edges"],
                row["delta0_a"], row["delta1_a"],
            )
            mid = _get_or_create_morph_id(con, vec)
            con.execute(
                "UPDATE theory SET "
                "N_rank2_0=?, N_rank2_1=?, N_bif=?, N_fund_0=?, N_fund_1=?, morph_id=? "
                "WHERE theory_id=?",
                (vec[0], vec[1], vec[2], vec[3], vec[4], mid, row["theory_id"]),
            )

        # Update aggregate counts
        con.execute("""
            UPDATE morphology_class SET n_theories = (
                SELECT COUNT(*) FROM theory WHERE theory.morph_id = morphology_class.morph_id
            )
        """)
        con.execute("""
            UPDATE morphology_class SET n_classes = (
                SELECT COUNT(DISTINCT class_id) FROM theory
                WHERE theory.morph_id = morphology_class.morph_id AND class_id IS NOT NULL
            )
        """)

    n_morph = con.execute("SELECT COUNT(*) FROM morphology_class").fetchone()[0]
    con.close()
    print(f"Done. {len(rows)} theories assigned to {n_morph} morphology classes.")


# ── Morphologies display ───────────────────────────────────────────────────────

def cmd_morphologies(args: argparse.Namespace) -> None:
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row

    if args.show is not None:
        # Show universality classes within a morphology
        mid = args.show
        morph = con.execute(
            "SELECT * FROM morphology_class WHERE morph_id=?", (mid,)
        ).fetchone()
        if morph is None:
            print(f"Morphology {mid} not found.")
            con.close()
            return
        print(f"\nMorphology {mid}: "
              f"N_r2₀={morph['N_rank2_0']}, N_r2₁={morph['N_rank2_1']}, "
              f"N_bif={morph['N_bif']}, N_f₀={morph['N_fund_0']}, "
              f"N_f₁={morph['N_fund_1']}  "
              f"({morph['n_theories']} theories, {morph['n_classes']} classes)")
        cls_rows = con.execute(
            "SELECT DISTINCT uc.class_id, uc.gauge_pair, uc.a_over_N2, uc.n_theories "
            "FROM universality_class uc "
            "JOIN theory t ON t.class_id = uc.class_id "
            "WHERE t.morph_id=? ORDER BY uc.a_over_N2",
            (mid,),
        ).fetchall()
        print(f"\n  {'CID':>5}  {'Gauge pair':<12}  {'a/N²':>12}  {'#th':>4}")
        print(f"  {'─'*40}")
        for r in cls_rows:
            print(f"  {r['class_id']:>5}  {r['gauge_pair']:<12}  {r['a_over_N2']:>12.6f}  {r['n_theories']:>4}")
    else:
        morph_rows = con.execute(
            "SELECT * FROM morphology_class "
            "ORDER BY (N_rank2_0+N_rank2_1), N_bif, N_fund_0, N_fund_1, morph_id"
        ).fetchall()
        print(f"\n  {'MorphID':>7}  {'N_r2₀':>6}  {'N_r2₁':>6}  {'N_bif':>5}  "
              f"{'N_f₀':>4}  {'N_f₁':>4}  {'#th':>5}  {'#cls':>5}")
        print(f"  {'─'*60}")
        for r in morph_rows:
            print(f"  {r['morph_id']:>7}  {r['N_rank2_0']:>6.1f}  {r['N_rank2_1']:>6.1f}  "
                  f"{r['N_bif']:>5}  {r['N_fund_0']:>4}  {r['N_fund_1']:>4}  "
                  f"{r['n_theories']:>5}  {r['n_classes']:>5}")
        print(f"\n  {len(morph_rows)} morphology classes total.")

    con.close()


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
    p.add_argument("--force", action="store_true", help="Overwrite existing DB without asking")
    p.add_argument("--exact-timeout", type=int, default=30, dest="exact_timeout",
                   help="Timeout in seconds per exact solve (default 30)")

    def _parse_rank(s):
        parts = s.split(",")
        if len(parts) != 2:
            raise argparse.ArgumentTypeError("--rank expects 'm0,m1' e.g. '2,1'")
        return (int(parts[0]), int(parts[1]))

    def _parse_class_id(s):
        if s.lower() == "null":
            return None
        return int(s)

    # classes
    p = sub.add_parser("classes", help="List universality classes")
    p.add_argument("--pair", help="Filter by gauge pair, e.g. SU-SU or SU×SU")
    p.add_argument("--min-a", type=float, dest="min_a")
    p.add_argument("--max-a", type=float, dest="max_a")
    p.add_argument("--rank", type=_parse_rank,
                   help="Filter by rank multipliers, e.g. '2,1' for SU(2N)×SU(N)")
    vg = p.add_mutually_exclusive_group()
    vg.add_argument("--veneziano",    dest="veneziano", action="store_true",  default=None,
                    help="Show only classes where all theories require Veneziano limit")
    vg.add_argument("--no-veneziano", dest="veneziano", action="store_false",
                    help="Show only classes where no theory requires Veneziano limit")

    # show
    p = sub.add_parser("show", help="Show all theories in a universality class")
    p.add_argument("class_id", type=_parse_class_id,
                   help="Class ID (from 'classes') or 'null' for diverged theories")

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
    vg = p.add_mutually_exclusive_group()
    vg.add_argument("--veneziano",    dest="veneziano", action="store_true",  default=None,
                    help="Show only theories requiring Veneziano limit")
    vg.add_argument("--no-veneziano", dest="veneziano", action="store_false",
                    help="Show only theories not requiring Veneziano limit")

    # stats
    sub.add_parser("stats", help="Database summary and build metadata")

    # morph-build
    sub.add_parser("morph-build",
                   help="Populate morphology columns on existing DB (no Mathematica needed)")

    # morphologies
    p = sub.add_parser("morphologies", help="List secondary morphology classes")
    p.add_argument("--show", type=int, metavar="MORPH_ID",
                   help="Show universality classes within a morphology")

    args = parser.parse_args()

    if args.cmd != "build":
        try:
            open(args.db).close()
        except FileNotFoundError:
            print(f"Database '{args.db}' not found. Run 'build' first.")
            sys.exit(1)

    dispatch = {
        "build":        cmd_build,
        "classes":      cmd_classes,
        "show":         cmd_show,
        "search":       cmd_search,
        "stats":        cmd_stats,
        "morph-build":  cmd_morph_build,
        "morphologies": cmd_morphologies,
    }
    dispatch[args.cmd](args)


if __name__ == "__main__":
    main()
