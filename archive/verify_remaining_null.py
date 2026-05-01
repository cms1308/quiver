"""For each of the 115 class_id-NULL theories, re-run a-max symbolically and
report the reason each one fails — either the anomaly constraint has no
solution, or every critical point has a Hessian with a positive eigenvalue.

Writes a JSON report: for each theory, the failure category and (when a
critical point exists) the Hessian eigenvalues at the best candidate."""

from __future__ import annotations

import json
import os
import subprocess
import sqlite3
import tempfile
import time
from math import gcd
from pathlib import Path

from quiver_generation import (
    enumerate_quivers, enumerate_quivers_mixed_rank, _dedup_symmetries,
)
from a_maximization_large_N import (
    _fmt_matter, _fmt_edges, build_fields_large_N,
    _anomaly_matrix_exact, dim_group_lead,
)

DB = Path(__file__).parent / "quivers.db"


def _f(x): return f"{x.numerator}/{x.denominator}"
def _v(v): return "{" + ",".join(_f(x) for x in v) + "}"
def _m(M): return "{" + ",".join(_v(r) for r in M) + "}"


def main():
    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row
    nulls = con.execute(
        "SELECT theory_id, gauge_pair, rank0_mult, rank1_mult, "
        "matter0, matter1, edges "
        "FROM theory WHERE class_id IS NULL"
    ).fetchall()
    con.close()
    print(f"Null theories: {len(nulls)}")

    wanted = {(r['gauge_pair'], r['rank0_mult'], r['rank1_mult'],
               r['matter0'], r['matter1'], r['edges']): r['theory_id']
              for r in nulls}

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

    matched: list = []
    for q, _ in _dedup_symmetries(equal + mixed):
        g0, g1 = q.gauge_types
        sig = (f"{g0}-{g1}", q.rank_multipliers[0], q.rank_multipliers[1],
               _fmt_matter(q.node_matter[0]),
               _fmt_matter(q.node_matter[1]),
               _fmt_edges(q))
        tid = wanted.get(sig)
        if tid is not None:
            matched.append((tid, q))
    print(f"Matched {len(matched)}/{len(wanted)}")

    # Build batch Mathematica script
    theory_data = []
    for tid, q in matched:
        fields = build_fields_large_N(q)
        A, b = _anomaly_matrix_exact(fields, q)
        dims = [f.dim_lead for f in fields]
        mults = q.rank_multipliers
        gaugino = sum(dim_group_lead(g) * m * m
                      for g, m in zip(q.gauge_types, mults))
        theory_data.append({
            "tid": tid, "A": _m(A), "b": _v(b),
            "dims": _v(dims), "gaugino": _f(gaugino),
        })

    assoc = ",\n".join(
        f'<|"tid"->{d["tid"]},"A"->{d["A"]},"b"->{d["b"]},'
        f'"dims"->{d["dims"]},"gaugino"->{d["gaugino"]}|>'
        for d in theory_data
    )

    tmpdir = tempfile.mkdtemp()
    script_path = os.path.join(tmpdir, "verify.wl")
    result_path = os.path.join(tmpdir, "verify.json")

    script = """
LaunchKernels[];
processOne[th_] := Module[
  {tid, Amat, bvec, dims, gaugino, Fnull, nfree, Fmat, R0, svars, R, r1,
   trR, trR3, aF, grad, hess, solsSym, freeParams, report, eigsPerSol},
  tid = th["tid"]; Amat = th["A"]; bvec = th["b"];
  dims = th["dims"]; gaugino = th["gaugino"];
  R0 = Quiet[LinearSolve[Amat, bvec]];
  If[!ListQ[R0],
    Return[<|"tid" -> tid, "reason" -> "LinearSolve_nosol"|>]];
  Fnull = NullSpace[Amat]; nfree = Length[Fnull]; Fmat = Transpose[Fnull];
  If[nfree == 0,
    Return[<|"tid" -> tid, "reason" -> "fully_constrained",
             "a" -> N[(3/32)*(3*(gaugino + dims.((R0 - 1)^3)) - (gaugino + dims.(R0 - 1))), 20]|>]];
  svars = Array[s, nfree];
  R = R0 + Fmat.svars; r1 = R - 1;
  trR = gaugino + dims.r1; trR3 = gaugino + dims.(r1^3);
  aF = (3/32)*(3*trR3 - trR);
  grad = Table[D[aF, svars[[j]]], {j, nfree}];
  hess = Outer[D[aF, #1, #2]&, svars, svars];
  solsSym = TimeConstrained[Quiet[Solve[Thread[grad == 0], svars]], 120, $Failed];
  If[solsSym === $Failed,
    Return[<|"tid" -> tid, "reason" -> "Solve_timeout", "nfree" -> nfree|>]];
  If[!ListQ[solsSym] || Length[solsSym] == 0,
    Return[<|"tid" -> tid, "reason" -> "Solve_no_solutions", "nfree" -> nfree|>]];
  freeParams = Variables[svars /. solsSym];
  If[Length[freeParams] > 0,
    solsSym = (Join[#, Thread[freeParams -> 0]]&) /@ solsSym];
  (* Keep only solutions with all-real svars *)
  solsSym = Select[solsSym,
    Function[rules,
      AllTrue[N[svars /. rules, 30],
        NumericQ[#] && Abs[Im[#]] < 10^-9 &]]];
  If[Length[solsSym] == 0,
    Return[<|"tid" -> tid, "reason" -> "no_real_critical_points", "nfree" -> nfree|>]];
  (* Collect Hessian eigenvalues at each real critical point *)
  eigsPerSol = Table[
    Sort[Eigenvalues[N[hess /. solsSym[[k]], 30]]], {k, Length[solsSym]}];
  (* Does ANY critical point have a NegSemiDef Hessian? *)
  anyNSD = AnyTrue[eigsPerSol, AllTrue[#, # <= 10^-9 &]&];
  <|"tid" -> tid, "reason" -> If[anyNSD, "unexpected_should_have_passed", "indefinite_at_all_critical_points"],
    "nfree" -> nfree,
    "n_real_critpts" -> Length[solsSym],
    "eigs" -> N[eigsPerSol, 6],
    "avals" -> N[aF /. solsSym, 8]|>
];

theories = {__ASSOC__};
results = ParallelMap[processOne, theories];
Export["__PATH__", results, "JSON"];
""".replace("__ASSOC__", assoc).replace("__PATH__", result_path)

    with open(script_path, "w") as fh:
        fh.write(script)
    print("Running Mathematica verification...")
    t0 = time.time()
    proc = subprocess.run(["wolframscript", "-file", script_path],
                   timeout=3600, capture_output=True, text=True)
    print(f"Done in {time.time()-t0:.1f}s")
    print("STDOUT:", proc.stdout[-2000:])
    print("STDERR:", proc.stderr[-2000:])
    print("script:", script_path)

    # Parse
    import re as _re
    with open(result_path) as fh:
        raw_text = fh.read()
    raw_text = _re.sub(r'(\d)\.(?!\d)', r'\1.0', raw_text)
    raw = json.loads(raw_text)

    # Summarize
    by_reason: dict = {}
    for r in raw:
        by_reason.setdefault(r["reason"], []).append(r["tid"])
    print(f"\nFailure reasons across {len(raw)} null theories:")
    for reason, tids in sorted(by_reason.items(), key=lambda x: -len(x[1])):
        print(f"  {reason}: {len(tids)}")
        if len(tids) <= 5:
            print(f"    tids: {tids}")

    # Any alarming results?
    unexpected = by_reason.get("unexpected_should_have_passed", [])
    if unexpected:
        print(f"\n⚠ {len(unexpected)} theories claim a NegSemiDef Hessian but came back None!")
        print(f"  tids: {unexpected}")

    # Save for reference
    report_path = Path(__file__).parent / "remaining_null_verification.json"
    with open(report_path, "w") as fh:
        json.dump(raw, fh, indent=2)
    print(f"\nFull report → {report_path}")


if __name__ == "__main__":
    main()
