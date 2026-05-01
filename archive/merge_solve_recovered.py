"""Re-run a-max on the 207 class_id-NULL theories using the new Solve
fallback, and assign each newly-converged theory to an existing
universality class (matched by a/N² within tolerance) or create a
fresh class if no match exists. Preserves existing class_ids."""

from __future__ import annotations

import sqlite3
import time
from math import gcd
from pathlib import Path

from quiver_generation import (
    enumerate_quivers, enumerate_quivers_mixed_rank, _dedup_symmetries,
)
from a_maximization_large_N import (
    a_maximize_batch_mathematica, _fmt_matter, _fmt_edges,
)
from two_node_db import _format_R_numerical

DB = Path(__file__).parent / "quivers.db"
TOL = 1e-5


def _sig(pair, m0, m1, matter0, matter1, edges):
    return (pair, m0, m1, matter0, matter1, edges)


def _sig_q(q):
    g0, g1 = q.gauge_types
    return _sig(f"{g0}-{g1}", q.rank_multipliers[0], q.rank_multipliers[1],
                _fmt_matter(q.node_matter[0]),
                _fmt_matter(q.node_matter[1]),
                _fmt_edges(q))


def main():
    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row

    diverged = con.execute(
        "SELECT theory_id, gauge_pair, rank0_mult, rank1_mult, "
        "matter0, matter1, edges, veneziano "
        "FROM theory WHERE class_id IS NULL"
    ).fetchall()
    print(f"Diverged theories in DB: {len(diverged)}")

    wanted = {_sig(r["gauge_pair"], r["rank0_mult"], r["rank1_mult"],
                   r["matter0"], r["matter1"], r["edges"]): r["theory_id"]
              for r in diverged}

    # Enumerate quivers to find matching Quiver objects
    t0 = time.time()
    equal = enumerate_quivers(n_nodes=2, max_multiedge=4, min_multiedge=1,
                              require_connected=True)
    mixed = []
    for m0 in range(1, 5):
        for m1 in range(1, 5):
            if (m0, m1) == (1, 1) or gcd(m0, m1) != 1:
                continue
            mixed.extend(enumerate_quivers_mixed_rank(
                n_nodes=2, rank_multipliers=[m0, m1],
                max_multiedge=4, min_multiedge=1, require_connected=True))
    all_qb = _dedup_symmetries(equal + mixed)
    print(f"Enumerated {len(all_qb)} quivers in {time.time()-t0:.1f}s")

    matched = []
    for q, _ in all_qb:
        tid = wanted.get(_sig_q(q))
        if tid is not None:
            matched.append((tid, q))
    print(f"Matched {len(matched)}/{len(wanted)} signatures")

    t1 = time.time()
    results = a_maximize_batch_mathematica([q for _, q in matched], timeout=3600)
    print(f"Mathematica re-scan done in {time.time()-t1:.1f}s")

    # ── Partition results ────────────────────────────────────────────────────
    converged = []  # list of (theory_id, quiver, fast_result)
    still_none = 0
    for (tid, q), res in zip(matched, results):
        if res is None:
            still_none += 1
        else:
            converged.append((tid, q, res))
    print(f"Converged: {len(converged)} / Still diverged: {still_none}\n")

    if not converged:
        print("Nothing to merge.")
        return

    # ── Load existing classes (with a/N²) grouped by (pair, m0, m1) ─────────
    existing_classes = con.execute(
        "SELECT class_id, gauge_pair, rank0_mult, rank1_mult, a_over_N2, "
        "veneziano_any, veneziano_all, n_theories "
        "FROM universality_class WHERE a_over_N2 IS NOT NULL"
    ).fetchall()
    by_sector: dict[tuple, list] = {}
    for c in existing_classes:
        key = (c["gauge_pair"], c["rank0_mult"], c["rank1_mult"])
        by_sector.setdefault(key, []).append(dict(c))

    max_cid = con.execute("SELECT MAX(class_id) FROM universality_class").fetchone()[0]

    # ── Assign each converged theory ─────────────────────────────────────────
    updates = []           # list of (class_id, a, c, a_over_c, R_str, theory_id)
    class_delta: dict = {} # class_id -> extra veneziano flags list
    new_classes = []       # list of (class_id, pair, m0, m1, centroid, veneziano)

    for tid, q, res in converged:
        pair = "-".join(q.gauge_types)
        m0, m1 = q.rank_multipliers
        key = (pair, m0, m1)
        a = res.a_over_N2
        c = res.c_over_N2
        a_over_c = a / c if c else None
        R_str = _format_R_numerical(res)
        ven = next((r["veneziano"] for r in diverged if r["theory_id"] == tid), 0)

        # Find existing class in this sector within TOL
        best_cid = None
        for ec in by_sector.get(key, []):
            if abs(ec["a_over_N2"] - a) < TOL:
                best_cid = ec["class_id"]
                break

        if best_cid is None:
            max_cid += 1
            best_cid = max_cid
            new_classes.append((best_cid, pair, m0, m1, a, ven))
            by_sector.setdefault(key, []).append({
                "class_id": best_cid, "gauge_pair": pair,
                "rank0_mult": m0, "rank1_mult": m1,
                "a_over_N2": a, "veneziano_any": ven, "veneziano_all": ven,
                "n_theories": 0,
            })

        updates.append((best_cid, a, c, a_over_c, R_str, ven, tid))
        class_delta.setdefault(best_cid, []).append(ven)

    print(f"Assigned to existing classes: "
          f"{sum(1 for u in updates if u[0] <= 913)}")
    print(f"New classes created: {len(new_classes)}")
    for cid, pair, m0, m1, a, ven in new_classes:
        print(f"  +class {cid}  {pair}({m0},{m1})  a/N²={a:.6f}")

    # ── Write all updates atomically ─────────────────────────────────────────
    with con:
        for cid, pair, m0, m1, a, ven in new_classes:
            con.execute(
                "INSERT INTO universality_class "
                "(class_id, gauge_pair, rank0_mult, rank1_mult, "
                " a_over_N2, n_theories, a_exact, c_exact, R_exact, "
                " a_over_c, veneziano_any, veneziano_all) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                (cid, pair, m0, m1, a, 0, None, None, None, None, ven, ven),
            )

        for cid, a, c, aoc, R_str, ven, tid in updates:
            con.execute(
                "UPDATE theory SET class_id=?, a_over_N2=?, c_over_N2=?, "
                "a_over_c=?, R_numerical=? WHERE theory_id=?",
                (cid, a, c, aoc, R_str, tid),
            )

        # Recompute n_theories and veneziano flags for every affected class
        for cid in set(u[0] for u in updates):
            row = con.execute(
                "SELECT COUNT(*) AS n, "
                "MAX(veneziano) AS ven_any, MIN(veneziano) AS ven_all "
                "FROM theory WHERE class_id=?", (cid,)
            ).fetchone()
            con.execute(
                "UPDATE universality_class "
                "SET n_theories=?, veneziano_any=?, veneziano_all=? "
                "WHERE class_id=?",
                (row["n"], row["ven_any"], row["ven_all"], cid),
            )

    # ── Final check ──────────────────────────────────────────────────────────
    n_null_after = con.execute(
        "SELECT COUNT(*) FROM theory WHERE class_id IS NULL"
    ).fetchone()[0]
    n_classes_after = con.execute(
        "SELECT COUNT(*) FROM universality_class"
    ).fetchone()[0]
    print(f"\nAfter merge: {n_null_after} theories still NULL  |  "
          f"{n_classes_after} total classes")
    con.close()


if __name__ == "__main__":
    main()
