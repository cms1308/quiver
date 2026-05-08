"""Consolidate `quivers.db` by merging theories related by per-node charge
conjugation (in addition to global conjugation + node swap, which were already
handled at enumeration time).

The dedup signature group is (Z_2)^k_conj × S_k_swap; for a 2-node quiver this
is 8 elements. Existing _dedup_symmetries inside quiver_generation.py only
considered the order-4 subgroup {e, c_global, swap, swap·c_global}, missing the
single-node conjugations c_0 and c_1 (and their swap composites).

This script:
  1. Backs up the DB to quivers.db.bak-pre-conj-dedup.
  2. Computes the orbit canonical signature for every theory row.
  3. Picks one survivor per orbit (the row whose own _quiver_signature equals
     the canonical, ties broken by min theory_id).
  4. Sanity-checks that orbit-mates agree on physically-derived columns
     (class_id, a_over_N2, c_over_N2, veneziano).
  5. Promotes OR'd has_flat_direction onto the survivor.
  6. Repoints universality_class.rep_theory_id off any deleted theory.
  7. Deletes the loser rows.
  8. Refreshes n_theories on universality_class and morphology_class.
  9. VACUUMs.

Usage:
  python3 scripts/dedup_charge_conjugation.py --dry-run
  python3 scripts/dedup_charge_conjugation.py --apply
"""

from __future__ import annotations

import argparse
import os
import shutil
import sqlite3
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from marginal_operators import quiver_from_row
from quiver_generation import (
    _orbit_canonical_signature,
    _quiver_signature,
)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DB = os.path.join(REPO_ROOT, "quivers.db")
BACKUP_PATH = os.path.join(REPO_ROOT, "quivers.db.bak-pre-conj-dedup")

# Columns that must agree across orbit mates (physical observables).
CONSISTENCY_COLS = (
    "class_id", "a_over_N2", "c_over_N2", "veneziano",
    "gauge_pair", "rank0_mult", "rank1_mult",
)
# Columns OR'd onto the survivor.
OR_COLS = ("has_flat_direction",)


def compute_orbits(rows: list[sqlite3.Row]) -> dict[tuple, list[int]]:
    by_orbit: dict[tuple, list[int]] = defaultdict(list)
    for r in rows:
        q = quiver_from_row(dict(r))
        canon = _orbit_canonical_signature(q)
        by_orbit[canon].append(r["theory_id"])
    return by_orbit


def pick_survivor(tids: list[int], rows_by_id: dict[int, sqlite3.Row],
                  canon: tuple) -> int:
    """Survivor is the row whose own signature is the canonical (so the kept
    row's matter0/matter1/edges strings are already in canonical form). Ties
    broken by min theory_id."""
    matching = []
    for tid in tids:
        q = quiver_from_row(dict(rows_by_id[tid]))
        if _quiver_signature(q) == canon:
            matching.append(tid)
    if not matching:
        # Fall back: lowest theory_id (existing row will be left as-is — the
        # canonical signature differs from the existing serialization, but
        # symmetry-equivalent so still physically the same theory).
        return min(tids)
    return min(matching)


def check_consistency(tids: list[int], rows_by_id: dict[int, sqlite3.Row]) -> list[str]:
    if len(tids) <= 1:
        return []
    issues: list[str] = []
    ref = rows_by_id[tids[0]]
    for col in CONSISTENCY_COLS:
        ref_val = ref[col]
        for tid in tids[1:]:
            v = rows_by_id[tid][col]
            if v != ref_val:
                if col in ("a_over_N2", "c_over_N2") and ref_val is not None and v is not None:
                    if abs(float(v) - float(ref_val)) < 1e-6:
                        continue
                issues.append(
                    f"orbit {tids}: {col} mismatch ({tids[0]}={ref_val!r} vs {tid}={v!r})"
                )
    return issues


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=DEFAULT_DB)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--show-orbits", type=int, default=10,
                    help="number of orbit examples to print")
    args = ap.parse_args()

    if not args.dry_run and not args.apply:
        ap.error("specify --dry-run or --apply")

    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    rows = con.execute("SELECT * FROM theory ORDER BY theory_id").fetchall()
    rows_by_id = {r["theory_id"]: r for r in rows}
    total = len(rows)
    print(f"Total theories: {total}")

    by_orbit = compute_orbits(rows)
    n_orbits = len(by_orbit)
    print(f"Orbits found:   {n_orbits}")
    print(f"Duplicates:     {total - n_orbits}")

    survivors: dict[tuple, int] = {}
    losers: list[int] = []
    or_promote: dict[int, dict[str, int]] = {}
    consistency_issues: list[str] = []
    orbit_size_dist: dict[int, int] = defaultdict(int)
    pair_loser_count: dict[str, int] = defaultdict(int)

    for canon, tids in by_orbit.items():
        orbit_size_dist[len(tids)] += 1
        survivors[canon] = pick_survivor(tids, rows_by_id, canon)
        for t in tids:
            if t != survivors[canon]:
                losers.append(t)
                pair_loser_count[rows_by_id[t]["gauge_pair"]] += 1
        if len(tids) > 1:
            consistency_issues.extend(check_consistency(tids, rows_by_id))
            for col in OR_COLS:
                vals = [rows_by_id[t][col] or 0 for t in tids]
                if any(vals) and not all(vals):
                    or_promote.setdefault(survivors[canon], {})[col] = 1

    print()
    print("Orbit-size distribution:")
    for sz in sorted(orbit_size_dist):
        print(f"  size {sz}: {orbit_size_dist[sz]} orbits")
    print()
    print("Losers by gauge_pair:")
    for pair, n in sorted(pair_loser_count.items()):
        print(f"  {pair}: {n}")
    print()
    print(f"Annotation conflicts: {len(consistency_issues)}")
    for issue in consistency_issues[:10]:
        print(f"  {issue}")
    if len(consistency_issues) > 10:
        print(f"  ... and {len(consistency_issues) - 10} more")

    print()
    print(f"Orbit examples (up to {args.show_orbits}):")
    multi = [(c, t) for c, t in by_orbit.items() if len(t) > 1]
    for canon, tids in multi[: args.show_orbits]:
        print(f"  orbit (canonical gauge_pair={rows_by_id[survivors[canon]]['gauge_pair']}):")
        for t in tids:
            r = rows_by_id[t]
            tag = "[KEEP]" if t == survivors[canon] else "[DROP]"
            print(f"    {tag} tid={t} class={r['class_id']} matter=({r['matter0']!r},{r['matter1']!r}) edges={r['edges']!r}")

    print()
    print(f"Will delete {len(losers)} theories; {total - len(losers)} survive.")
    print(f"OR-promotions on survivor: {len(or_promote)}")

    if consistency_issues:
        print()
        print("ABORT: annotation conflicts detected. Resolve upstream before --apply.")
        return 2

    if args.dry_run:
        print()
        print("Dry run complete. Re-run with --apply to commit.")
        return 0

    # ── Apply ─────────────────────────────────────────────────────────────────
    print()
    print(f"Backing up DB → {BACKUP_PATH}")
    shutil.copy(args.db, BACKUP_PATH)

    with con:
        # Promote OR'd flags onto survivors
        for sid, cols in or_promote.items():
            for col, val in cols.items():
                con.execute(
                    f"UPDATE theory SET {col}=? WHERE theory_id=?", (val, sid)
                )

        # Repoint rep_theory_id off losers
        loser_set = set(losers)
        repoint_count = 0
        rep_rows = con.execute(
            "SELECT class_id, rep_theory_id FROM universality_class"
        ).fetchall()
        for rr in rep_rows:
            if rr["rep_theory_id"] in loser_set:
                # Find a surviving theory in this class
                surviving = con.execute(
                    "SELECT theory_id FROM theory WHERE class_id=? "
                    "AND theory_id NOT IN ({}) "
                    "ORDER BY length(edges), matter0, matter1, theory_id "
                    "LIMIT 1".format(",".join(str(x) for x in losers)),
                    (rr["class_id"],),
                ).fetchone()
                if surviving is not None:
                    con.execute(
                        "UPDATE universality_class SET rep_theory_id=? WHERE class_id=?",
                        (surviving["theory_id"], rr["class_id"]),
                    )
                    repoint_count += 1

        print(f"  Repointed rep_theory_id on {repoint_count} classes")

        # Delete losers (chunked to avoid SQLite parameter limit)
        CHUNK = 500
        for i in range(0, len(losers), CHUNK):
            batch = losers[i : i + CHUNK]
            con.execute(
                f"DELETE FROM theory WHERE theory_id IN ({','.join('?'*len(batch))})",
                batch,
            )
        print(f"  Deleted {len(losers)} theory rows")

        # Refresh n_theories
        con.execute(
            "UPDATE universality_class SET n_theories = ("
            "  SELECT COUNT(*) FROM theory WHERE class_id = universality_class.class_id"
            ")"
        )
        con.execute(
            "UPDATE morphology_class SET n_theories = ("
            "  SELECT COUNT(*) FROM theory WHERE morph_id = morphology_class.morph_id"
            ")"
        )
        # n_classes per morphology: unchanged in principle (each orbit lives in
        # one class), but recompute defensively from current theory rows.
        con.execute(
            "UPDATE morphology_class SET n_classes = ("
            "  SELECT COUNT(DISTINCT class_id) FROM theory "
            "   WHERE morph_id = morphology_class.morph_id"
            "     AND class_id IS NOT NULL"
            ")"
        )

    con.execute("VACUUM")
    con.close()

    # Verify
    con = sqlite3.connect(args.db)
    con.row_factory = sqlite3.Row
    surviving = con.execute("SELECT COUNT(*) FROM theory").fetchone()[0]
    n_classes = con.execute("SELECT COUNT(*) FROM universality_class").fetchone()[0]
    n_classes_used = con.execute(
        "SELECT COUNT(DISTINCT class_id) FROM theory WHERE class_id IS NOT NULL"
    ).fetchone()[0]
    n_scft = con.execute(
        "SELECT COUNT(*) FROM theory WHERE class_id IS NOT NULL"
    ).fetchone()[0]
    n_nonscft = con.execute(
        "SELECT COUNT(*) FROM theory WHERE class_id IS NULL"
    ).fetchone()[0]
    con.close()

    print()
    print(f"After dedup: {surviving} theories ({n_scft} SCFT in {n_classes_used} classes, {n_nonscft} non-SCFT)")
    print(f"Universality_class rows: {n_classes}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
