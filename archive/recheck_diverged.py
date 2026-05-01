"""Re-run Mathematica a-maximization on all class_id-NULL (diverged) theories
and report whether any now converge to a bounded maximum."""

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

DB = Path(__file__).parent / "quivers.db"


def _signature(gauge_pair, m0, m1, matter0, matter1, edges):
    return (gauge_pair, m0, m1, matter0, matter1, edges)


def _sig_from_quiver(q):
    g0, g1 = q.gauge_types[0], q.gauge_types[1]
    pair = f"{g0}-{g1}"
    m0, m1 = q.rank_multipliers[0], q.rank_multipliers[1]
    return _signature(pair, m0, m1,
                      _fmt_matter(q.node_matter[0]),
                      _fmt_matter(q.node_matter[1]),
                      _fmt_edges(q))


def main():
    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row
    diverged = con.execute(
        "SELECT theory_id, gauge_pair, rank0_mult, rank1_mult, "
        "matter0, matter1, edges "
        "FROM theory WHERE class_id IS NULL"
    ).fetchall()
    con.close()

    wanted = {
        _signature(r["gauge_pair"], r["rank0_mult"], r["rank1_mult"],
                   r["matter0"], r["matter1"], r["edges"]): r["theory_id"]
        for r in diverged
    }
    print(f"Diverged theories in DB: {len(diverged)}  (unique sigs: {len(wanted)})")

    # Re-enumerate exactly as build does
    t0 = time.time()
    equal = enumerate_quivers(n_nodes=2, max_multiedge=4, min_multiedge=1,
                              require_connected=True)
    MAX_RANK = 4
    mixed_sets = []
    for m0 in range(1, MAX_RANK + 1):
        for m1 in range(1, MAX_RANK + 1):
            if (m0, m1) == (1, 1) or gcd(m0, m1) != 1:
                continue
            mixed_sets.extend(enumerate_quivers_mixed_rank(
                n_nodes=2, rank_multipliers=[m0, m1],
                max_multiedge=4, min_multiedge=1, require_connected=True,
            ))
    all_qb = _dedup_symmetries(equal + mixed_sets)
    print(f"Enumerated {len(all_qb)} quivers in {time.time()-t0:.1f}s")

    # Match to wanted signatures
    matched: list = []  # list of (theory_id, quiver)
    unmatched_ids = set(wanted.values())
    for q, _bounds in all_qb:
        sig = _sig_from_quiver(q)
        tid = wanted.get(sig)
        if tid is not None and tid in unmatched_ids:
            matched.append((tid, q))
            unmatched_ids.discard(tid)
    print(f"Matched {len(matched)}/{len(wanted)} "
          f"(unmatched DB rows: {len(unmatched_ids)})")

    if not matched:
        return

    # Re-run a-maximization
    t1 = time.time()
    quivers = [q for _, q in matched]
    results = a_maximize_batch_mathematica(quivers, timeout=3600)
    print(f"Mathematica re-scan done in {time.time()-t1:.1f}s")

    n_none = 0
    converged = []
    for (tid, q), res in zip(matched, results):
        if res is None:
            n_none += 1
        else:
            converged.append((tid, q, res))

    print(f"\nResult: {n_none}/{len(matched)} still diverge, "
          f"{len(converged)} now converge.")
    if converged:
        print("\nTheories that now converge:")
        for tid, q, res in converged:
            sig = _sig_from_quiver(q)
            print(f"  theory_id={tid}  {sig[0]}({sig[1]},{sig[2]})  "
                  f"matter=[{sig[3]} | {sig[4]}]  edges={sig[5]}  "
                  f"a/N²={res.a_over_N2:.6f}  c/N²={res.c_over_N2:.6f}")


if __name__ == "__main__":
    main()
