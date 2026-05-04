"""Sweep quivers.db, find always-marginal flavor-singlet operators per theory.

For each theory:
  1. Reconstruct Quiver from DB row.
  2. Run finite-N a-maximization at N ∈ {10, 20, 30}.
  3. Enumerate single-trace gauge-invariant operators up to max-degree.
  4. Filter to operators with |R − 2| < 1e-6 at every tested N AND flavor singlet.
  5. Dump to a TeX longtable.

A flavor-singlet operator marginal at every tested N is a candidate exactly
marginal deformation; the count exceeds the number of independent gauge
anomaly conditions iff there are non-trivial conformal manifold directions.
"""

from __future__ import annotations

import argparse
import sqlite3
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from marginal_operators import (
    CandidateOp,
    _flavor_multiplets,
    enumerate_candidates,
    is_flavor_singlet,
    is_marginal_at_all_N,
    op_R_at_N,
    quiver_from_row,
    r_values_at_N,
)

DB = Path(__file__).parent.parent / "quivers.db"
OUT = Path(__file__).parent.parent / "paper/sections/generated/marginal_operators.tex"

DEFAULT_N_LIST = (10, 20, 30)
DEFAULT_MAX_DEGREE = 6
DEFAULT_TOL = 1e-6

# ── Operator pretty-printing ──────────────────────────────────────────────────

_REP_TEX = {
    "adj": r"\mathrm{adj}",
    "S": "S", "Sbar": r"\bar{S}",
    "A": "A", "Abar": r"\bar{A}",
    "V": "V", "fund": r"\square", "antifund": r"\bar{\square}",
}


def _label_to_tex(label: str) -> str:
    if label.startswith("node"):
        node, rep = label[4:].split("_", 1)
        sym = _REP_TEX.get(rep, rep)
        return f"{sym}_{{{node}}}"
    parts = label.split("_")
    src, dst, rep = parts[1], parts[2], "_".join(parts[3:])
    rep_tex = {
        "pm": r"\square\bar{\square}",
        "pp": r"\square\square",
        "mm": r"\bar{\square}\bar{\square}",
        "p":  r"\square",
        "m":  r"\bar{\square}",
        "std": "V",
    }.get(rep, rep)
    return f"Q^{{{rep_tex}}}_{{{src}{dst}}}"


def _op_to_tex(op: CandidateOp) -> str:
    if op.kind == "bifund-loop" and op.word is not None:
        parts = [_label_to_tex(lbl) for lbl in op.word]
        return r"\mathrm{tr}\!\left(" + r"\,".join(parts) + r"\right)"
    parts = []
    for label, mult in op.factors:
        tex = _label_to_tex(label)
        parts.append(tex if mult == 1 else f"{tex}^{{{mult}}}")
    return r"\mathrm{tr}\!\left(" + r"\,".join(parts) + r"\right)"


def _factors_str(op: CandidateOp) -> str:
    return " · ".join(f"{lbl}^{m}" if m > 1 else lbl for lbl, m in op.factors)


# ── Sweep ──────────────────────────────────────────────────────────────────────

def sweep(
    con: sqlite3.Connection,
    where: str = "",
    params: tuple = (),
    N_list: tuple[int, ...] = DEFAULT_N_LIST,
    max_degree: int = DEFAULT_MAX_DEGREE,
    tol: float = DEFAULT_TOL,
    verbose: bool = False,
) -> list[dict]:
    con.row_factory = sqlite3.Row
    sql = "SELECT * FROM theory"
    if where:
        sql += " WHERE " + where
    sql += " ORDER BY theory_id"
    rows = list(con.execute(sql, params))

    results: list[dict] = []
    n_failed = 0
    n_processed = 0
    t0 = time.time()
    for row in rows:
        try:
            q = quiver_from_row(dict(row))
        except Exception as e:
            n_failed += 1
            if verbose:
                print(f"  theory_id={row['theory_id']}: parse failure: {e}", file=sys.stderr)
            continue
        try:
            R_per_N = {N: r_values_at_N(q, N) for N in N_list}
        except Exception as e:
            n_failed += 1
            if verbose:
                print(f"  theory_id={row['theory_id']}: a-max failure: {e}", file=sys.stderr)
            continue
        cands = enumerate_candidates(q, max_degree=max_degree)
        multiplets = _flavor_multiplets(q)
        marginals = []
        for op in cands:
            if not is_marginal_at_all_N(op, R_per_N, tol=tol):
                continue
            singlet = is_flavor_singlet(op, multiplets)
            marginals.append((op, singlet))
        n_processed += 1
        if marginals:
            singlet_marginals = [op for op, s in marginals if s]
            results.append({
                "theory_id": row["theory_id"],
                "class_id": row["class_id"],
                "gauge_pair": row["gauge_pair"],
                "matter0": row["matter0"],
                "matter1": row["matter1"],
                "edges": row["edges"],
                "n_marginal_total": len(marginals),
                "n_marginal_singlet": len(singlet_marginals),
                "marginals": marginals,
                "R_per_N": {N: {k: float(v) for k, v in r.items()} for N, r in R_per_N.items()},
            })
        if verbose and n_processed % 200 == 0:
            print(f"  ... processed {n_processed}/{len(rows)} ({time.time()-t0:.1f}s)", file=sys.stderr)

    print(f"Sweep complete: {n_processed} processed, {n_failed} failed, "
          f"{len(results)} with marginals  ({time.time()-t0:.1f}s)", file=sys.stderr)
    return results


# ── TeX output ────────────────────────────────────────────────────────────────

def write_tex(results: list[dict], out_path: Path, N_list: tuple[int, ...]) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as f:
        f.write("% Auto-generated from quivers.db by scripts/dump_marginal_operators.py\n")
        f.write("% Requires: \\usepackage{booktabs,longtable,amssymb,amsmath}\n")
        f.write("\\section*{Always-marginal flavor-singlet operators (finite $N$)}\n\n")

        f.write("Operators built from chiral matter fields with $R[O]=2$ at every $N\\in")
        f.write(f"\\{{{','.join(str(N) for N in N_list)}\\}}$ ")
        f.write("(numerical $a$-maximization, $|R-2|<10^{-6}$). Flavor singlets ")
        f.write("are candidates for exactly marginal deformations; non-singlets ")
        f.write("would generate flavor-current beta functions.\n\n")
        f.write(f"Total theories with at least one always-marginal flavor singlet: "
                f"\\textbf{{{sum(1 for r in results if r['n_marginal_singlet'] > 0)}}} of "
                f"\\textbf{{{len(results)}}} reporting any marginal.\n\n")

        # Group by class_id
        by_class: dict = {}
        for r in results:
            by_class.setdefault(r["class_id"], []).append(r)

        f.write("\\begin{longtable}{rrlllr}\n")
        f.write("\\toprule\n")
        f.write("class & theory & gauge & matter & operator & deg \\\\\n")
        f.write("\\midrule\n\\endhead\n")
        for cls_id in sorted(by_class.keys(), key=lambda x: (x is None, x)):
            for r in by_class[cls_id]:
                singlets = [op for op, s in r["marginals"] if s]
                if not singlets:
                    continue
                matter_str = f"{r['matter0']} \\| {r['matter1']}"
                cls_str = str(r["class_id"]) if r["class_id"] is not None else "—"
                for i, op in enumerate(singlets):
                    cls_cell = cls_str if i == 0 else ""
                    th_cell = str(r["theory_id"]) if i == 0 else ""
                    gp_cell = r["gauge_pair"] if i == 0 else ""
                    m_cell = matter_str if i == 0 else ""
                    f.write(
                        f"{cls_cell} & {th_cell} & {gp_cell} & {m_cell} & "
                        f"${_op_to_tex(op)}$ & {op.degree} \\\\\n"
                    )
        f.write("\\bottomrule\n\\end{longtable}\n")
    print(f"Wrote {out_path}", file=sys.stderr)


def write_summary_text(results: list[dict], N_list: tuple[int, ...]) -> None:
    """Print a summary to stdout."""
    n_with_singlet = sum(1 for r in results if r["n_marginal_singlet"] > 0)
    print(f"\nTheories with marginal operators: {len(results)}")
    print(f"Theories with flavor-singlet marginal operators: {n_with_singlet}")
    print(f"Tested N values: {N_list}")
    print()
    # Top 20 by class
    print("First 20 theories with flavor-singlet marginals:")
    print(f"{'theory':>6} {'class':>5} {'pair':<6} {'matter (0|1)':<30} {'#sing':>5} {'op example'}")
    shown = 0
    for r in results:
        if r["n_marginal_singlet"] == 0:
            continue
        singlets = [op for op, s in r["marginals"] if s]
        first = singlets[0]
        m_str = f"{r['matter0']} | {r['matter1']}"[:30]
        print(f"{r['theory_id']:>6} {str(r['class_id']):>5} {r['gauge_pair']:<6} "
              f"{m_str:<30} {r['n_marginal_singlet']:>5}  "
              f"deg={first.degree} {_factors_str(first)}")
        shown += 1
        if shown >= 20:
            break


# ── CLI ────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    sel = ap.add_mutually_exclusive_group()
    sel.add_argument("--theory-id", type=int, help="single theory to test")
    sel.add_argument("--class-id", type=int, help="all theories in this class")
    sel.add_argument("--all", action="store_true", help="sweep entire DB")
    ap.add_argument("--max-degree", type=int, default=DEFAULT_MAX_DEGREE)
    ap.add_argument("--N-list", type=str, default=",".join(str(N) for N in DEFAULT_N_LIST),
                    help="comma-separated N values, default 10,20,30")
    ap.add_argument("--tol", type=float, default=DEFAULT_TOL)
    ap.add_argument("--tex", action="store_true", help=f"write TeX to {OUT}")
    ap.add_argument("--out", type=str, default=str(OUT), help="output TeX path")
    ap.add_argument("--verbose", "-v", action="store_true")
    args = ap.parse_args()

    N_list = tuple(int(x) for x in args.N_list.split(","))
    con = sqlite3.connect(DB)
    if args.theory_id is not None:
        results = sweep(con, "theory_id=?", (args.theory_id,), N_list, args.max_degree, args.tol, args.verbose)
    elif args.class_id is not None:
        results = sweep(con, "class_id=?", (args.class_id,), N_list, args.max_degree, args.tol, args.verbose)
    elif args.all:
        results = sweep(con, "", (), N_list, args.max_degree, args.tol, args.verbose)
    else:
        ap.error("must specify --theory-id, --class-id, or --all")

    write_summary_text(results, N_list)
    if args.tex or args.all:
        write_tex(results, Path(args.out), N_list)


if __name__ == "__main__":
    main()
