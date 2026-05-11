"""Renumber theory_id and morph_id to be 1-based contiguous (no gaps).

After charge-conjugation dedup, primary keys are no longer contiguous:
  theory_id: 12539 rows but max=21315 (8776 gaps)
  morph_id:  1353 rows but max=21315 (19962 gaps)
  class_id:  876 rows, already 1..876

Mapping is computed by sorting current IDs ASC and assigning ROW_NUMBER().
Foreign-key references (theory.class_id, theory.morph_id, universality_class.rep_theory_id)
are updated atomically using a temporary high offset to avoid PK collisions
mid-update.

Usage:
    python3 scripts/renumber_pks.py --dry-run
    python3 scripts/renumber_pks.py --apply
"""

from __future__ import annotations

import argparse
import os
import shutil
import sqlite3
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DB = os.path.join(REPO_ROOT, "quivers.db")
BACKUP_PATH = os.path.join(REPO_ROOT, "quivers.db.bak-pre-pk-renumber")
OFFSET = 100_000_000  # safely above any current PK


def renumber_theory(con: sqlite3.Connection) -> tuple[int, int]:
    """Assign theory_id 1..N grouped by class_id (NULLs at the tail), with
    existing theory_id as the intra-class tiebreaker. This keeps theories in
    one universality class contiguous."""
    rows = con.execute(
        "SELECT theory_id, "
        "ROW_NUMBER() OVER ("
        "  ORDER BY (class_id IS NULL), class_id, theory_id"
        ") AS new_id "
        "FROM theory"
    ).fetchall()
    n_real = sum(1 for old, new in rows if old != new)
    max_new = max(new for _, new in rows) if rows else 0

    # Build TEMP map table
    con.execute("DROP TABLE IF EXISTS _tid_map")
    con.execute("CREATE TEMP TABLE _tid_map (old_id INTEGER PRIMARY KEY, new_id INTEGER NOT NULL)")
    con.executemany("INSERT INTO _tid_map(old_id, new_id) VALUES (?, ?)", rows)

    # Step 1: shift all theory_id and rep_theory_id off the live range
    con.execute("UPDATE theory SET theory_id = theory_id + ?", (OFFSET,))
    con.execute("UPDATE universality_class SET rep_theory_id = rep_theory_id + ?", (OFFSET,))

    # Step 2: write new IDs from the map
    con.execute(
        "UPDATE theory SET theory_id = "
        "  (SELECT new_id FROM _tid_map WHERE old_id = theory_id - ?)",
        (OFFSET,),
    )
    con.execute(
        "UPDATE universality_class SET rep_theory_id = "
        "  (SELECT new_id FROM _tid_map WHERE old_id = rep_theory_id - ?)",
        (OFFSET,),
    )

    con.execute("DROP TABLE _tid_map")

    # Reset AUTOINCREMENT counter
    con.execute(
        "UPDATE sqlite_sequence SET seq = ? WHERE name = 'theory'", (max_new,)
    )
    return n_real, max_new


def renumber_morph(con: sqlite3.Connection) -> tuple[int, int]:
    """Assign morph_id 1..M in order of current morph_id ASC."""
    rows = con.execute(
        "SELECT morph_id, ROW_NUMBER() OVER (ORDER BY morph_id) AS new_id "
        "FROM morphology_class"
    ).fetchall()
    n_real = sum(1 for old, new in rows if old != new)
    max_new = max(new for _, new in rows) if rows else 0

    con.execute("DROP TABLE IF EXISTS _mid_map")
    con.execute("CREATE TEMP TABLE _mid_map (old_id INTEGER PRIMARY KEY, new_id INTEGER NOT NULL)")
    con.executemany("INSERT INTO _mid_map(old_id, new_id) VALUES (?, ?)", rows)

    con.execute("UPDATE morphology_class SET morph_id = morph_id + ?", (OFFSET,))
    con.execute("UPDATE theory SET morph_id = morph_id + ?", (OFFSET,))

    con.execute(
        "UPDATE morphology_class SET morph_id = "
        "  (SELECT new_id FROM _mid_map WHERE old_id = morph_id - ?)",
        (OFFSET,),
    )
    con.execute(
        "UPDATE theory SET morph_id = "
        "  (SELECT new_id FROM _mid_map WHERE old_id = morph_id - ?)",
        (OFFSET,),
    )

    con.execute("DROP TABLE _mid_map")

    con.execute(
        "UPDATE sqlite_sequence SET seq = ? WHERE name = 'morphology_class'",
        (max_new,),
    )
    return n_real, max_new


def report_gaps(con: sqlite3.Connection, label: str) -> None:
    for table, col in (
        ("theory", "theory_id"),
        ("universality_class", "class_id"),
        ("morphology_class", "morph_id"),
    ):
        n, lo, hi = con.execute(
            f"SELECT COUNT(*), MIN({col}), MAX({col}) FROM {table}"
        ).fetchone()
        gaps = (hi or 0) - (n or 0)
        print(f"  [{label}] {table}.{col}: {n} rows, min={lo}, max={hi}, gaps={gaps}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=DEFAULT_DB)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    if not args.dry_run and not args.apply:
        ap.error("specify --dry-run or --apply")

    con = sqlite3.connect(args.db)
    print("Before:")
    report_gaps(con, "before")

    if args.dry_run:
        # Use a transaction we will rollback to preview
        con.execute("BEGIN")
        try:
            n_t, max_t = renumber_theory(con)
            n_m, max_m = renumber_morph(con)
            print()
            print(f"Would remap {n_t} theory_id values (new max = {max_t})")
            print(f"Would remap {n_m} morph_id values (new max = {max_m})")
            print()
            print("Preview (after):")
            report_gaps(con, "preview")
        finally:
            con.execute("ROLLBACK")
        return 0

    # Apply
    print()
    print(f"Backing up DB → {BACKUP_PATH}")
    shutil.copy(args.db, BACKUP_PATH)

    con.execute("BEGIN")
    try:
        n_t, max_t = renumber_theory(con)
        n_m, max_m = renumber_morph(con)
        con.execute("COMMIT")
    except Exception:
        con.execute("ROLLBACK")
        raise

    con.execute("VACUUM")
    con.close()
    con = sqlite3.connect(args.db)
    print()
    print("After:")
    report_gaps(con, "after")
    print()
    print(f"Remapped {n_t} theory_id, {n_m} morph_id values.")

    # Sanity: FK integrity
    orphans = con.execute(
        "SELECT COUNT(*) FROM universality_class "
        "WHERE rep_theory_id NOT IN (SELECT theory_id FROM theory)"
    ).fetchone()[0]
    print(f"Orphan rep_theory_id refs: {orphans}")
    return 0 if orphans == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
