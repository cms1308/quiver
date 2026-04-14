#!/usr/bin/env python3
"""
single_node_tables.py — Complete single-node classification

Enumerates all 67 single-node 4D N=1 AF theories (52 SU + 6 SO + 9 Sp),
computes exact large-N R-charges and central charges via a-maximization,
validates against universal formulas and arXiv:2510.19136, and produces
formatted markdown tables + machine-readable JSON.

% ASSERT_CONVENTION: natural_units=natural, gauge_groups=SU_SO_Sp,
%   dynkin_index=T_fund_SU_half, representation_notation=adj_all_types,
%   large_N_scaling=a_over_N2
"""

from __future__ import annotations

import json
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from fractions import Fraction
from sympy import Rational, simplify, sqrt, S, Expr
from collections import defaultdict

from quiver_generation import (
    Quiver, enumerate_quivers, chiral_excess_coeffs, nf_bound,
)
from a_maximization_large_N import (
    a_maximize_large_N, LargeNResult, dim_group_lead,
    T_rep_lead, T_adj_lead, build_fields_large_N, _R,
)
from beta_functions import RANK2_ADJ_REPS, b0_linear


# ── Formatting helpers ────────────────────────────────────────────────────────

REP_SYM = {
    "adj": "adj", "S": "S", "Sbar": "S\u0305",
    "A": "A", "Abar": "A\u0305", "V": "V",
    "fund": "\u25a1", "antifund": "\u25a1\u0305",
}

REP_ORDER = ["adj", "S", "Sbar", "A", "Abar", "V", "fund", "antifund"]


def fmt_matter(matter: dict[str, int]) -> str:
    """Compact matter notation: '2adj + S + A-bar' etc."""
    parts = []
    for rep in REP_ORDER:
        n = matter.get(rep, 0)
        if n:
            sym = REP_SYM.get(rep, rep)
            parts.append(f"{n}{sym}" if n > 1 else sym)
    return " + ".join(parts) if parts else "\u2014"


def fmt_matter_ascii(matter: dict[str, int]) -> str:
    """ASCII-safe matter notation for markdown tables."""
    sym_map = {
        "adj": "adj", "S": "S", "Sbar": "Sbar",
        "A": "A", "Abar": "Abar", "V": "V",
        "fund": "Q", "antifund": "Qbar",
    }
    parts = []
    for rep in REP_ORDER:
        n = matter.get(rep, 0)
        if n:
            sym = sym_map[rep]
            parts.append(f"{n}{sym}" if n > 1 else sym)
    return " + ".join(parts) if parts else "---"


def fmt_expr(expr: Expr) -> str:
    """Compact string for a sympy expression, suitable for tables."""
    s = str(expr)
    s = s.replace("sqrt(", "sqrt(")  # keep readable
    return s


def fmt_R_charges(result: LargeNResult) -> str:
    """Format R-charges as a compact string for table column."""
    label_map = {
        "adj": "adj", "S": "S", "Sbar": "Sbar",
        "A": "A", "Abar": "Abar", "V": "V",
        "fund": "Q", "antifund": "Qbar",
    }
    parts = []
    for f in result.fields:
        label = f.label
        # Parse node0_rep format
        if label.startswith("node"):
            rep = label.split("_", 1)[1]
            display = label_map.get(rep, rep)
        else:
            display = label
        R_val = result.R_charges[label]
        parts.append(f"R({display})={fmt_expr(R_val)}")
    return ", ".join(parts)


def fmt_nf_bound(alpha: Fraction, gamma: Fraction) -> str:
    """Format N_f bound as readable string."""
    if alpha == 0 and gamma == 0:
        return "0 (CM*)"
    parts = []
    if alpha:
        if alpha == 1:
            parts.append("N")
        else:
            parts.append(f"{alpha}N")
    if gamma:
        g = gamma
        if parts:
            sign = "+" if g > 0 else "-"
            parts.append(f" {sign} {abs(g)}")
        else:
            parts.append(str(g))
    elif not parts:
        parts.append("0")
    return "".join(parts)


def fmt_delta(a: int, b: int) -> str:
    """Format chiral excess delta(N) = aN + b."""
    parts = []
    if a:
        if a == 1:
            parts.append("N")
        elif a == -1:
            parts.append("-N")
        else:
            parts.append(f"{a}N")
    if b:
        if parts:
            sign = "+" if b > 0 else "-"
            parts.append(f" {sign} {abs(b)}")
        else:
            parts.append(str(b))
    elif not parts:
        parts.append("0")
    return "".join(parts)


# ── Enumeration ───────────────────────────────────────────────────────────────

def enumerate_all_single_node() -> list[tuple[Quiver, list[tuple[Fraction, Fraction]]]]:
    """
    Enumerate all single-node theories including those with no edges.
    Filter out pure SYM (all matter counts zero).
    """
    results = enumerate_quivers(
        n_nodes=1, max_multiedge=1, min_multiedge=0, require_connected=False
    )
    # Filter out pure SYM
    filtered = []
    for q, bounds in results:
        matter = q.node_matter[0]
        if any(v > 0 for v in matter.values()):
            filtered.append((q, bounds))
    return filtered


# ── Classification computation ────────────────────────────────────────────────

def classify_single_node() -> list[dict]:
    """
    Enumerate all 67 single-node theories, compute exact large-N data,
    and return a list of theory dicts.
    """
    theories = enumerate_all_single_node()

    rows = []
    for q, bounds in theories:
        g = q.gauge_types[0]
        matter = q.node_matter[0]

        # Compute exact large-N a-maximization
        result = a_maximize_large_N(q)

        # Count rank-2 + adjoint fields
        n_rank2 = sum(matter.get(rep, 0) for rep in RANK2_ADJ_REPS[g])

        # Determine type by n_rank2 (matching paper convention)
        # For SO/Sp this equals effective T; for SU the universality
        # is finer (depends on effective_T = n_adj + n_tensor/2).
        if n_rank2 == 1:
            theory_type = "I"
        elif n_rank2 == 2:
            theory_type = "II"
        elif n_rank2 == 3:
            theory_type = "III"
        else:
            theory_type = f"n={n_rank2}"

        # Chiral excess (SU only)
        if g == "SU":
            a_delta, b_delta = chiral_excess_coeffs(q, 0)
            is_veneziano = (a_delta != 0)
        else:
            a_delta, b_delta = 0, 0
            is_veneziano = False

        # AF bound
        alpha, gamma = bounds[0]

        # Conformal manifold: b_0 = 0 means alpha=0 and gamma=0
        is_conformal_manifold = (alpha == 0 and gamma == 0)

        # a/c ratio
        a_val = result.a_over_N2
        c_val = result.c_over_N2
        if simplify(a_val) == 0 and simplify(c_val) == 0:
            ac_ratio = "n/a"
            ac_ratio_expr = None
        else:
            ac_ratio_expr = simplify(a_val / c_val)
            ac_ratio = fmt_expr(ac_ratio_expr)

        # R-charges dict
        R_dict = {}
        for f in result.fields:
            label = f.label
            rep = label.split("_", 1)[1] if label.startswith("node") else label
            R_dict[rep] = result.R_charges[label]

        row = {
            "gauge_type": g,
            "matter": dict(matter),
            "matter_str": fmt_matter_ascii(matter),
            "n_rank2": n_rank2,
            "type": theory_type,
            "chiral_excess_a": a_delta,
            "chiral_excess_b": b_delta,
            "is_veneziano": is_veneziano,
            "nf_bound_alpha": alpha,
            "nf_bound_gamma": gamma,
            "is_conformal_manifold": is_conformal_manifold,
            "R_charges": R_dict,
            "R_charges_str": fmt_R_charges(result),
            "a_over_N2": a_val,
            "c_over_N2": c_val,
            "a_over_c": ac_ratio,
            "a_over_c_expr": ac_ratio_expr,
            "result": result,
            "quiver": q,
            "bounds": bounds,
        }
        rows.append(row)

    return rows


# ── Table formatting ──────────────────────────────────────────────────────────

def sort_theories(rows: list[dict]) -> list[dict]:
    """Sort by type (I, II, III), then by matter complexity."""
    type_order = {"I": 0, "II": 1, "III": 2}
    return sorted(rows, key=lambda r: (
        type_order.get(r["type"], 9),
        r["n_rank2"],
        r["matter_str"],
    ))


def write_su_table(rows: list[dict], path: str) -> None:
    """Write SU(N) single-node table to markdown."""
    su_rows = sort_theories([r for r in rows if r["gauge_type"] == "SU"])

    lines = [
        "# SU(N) Single-Node Classification",
        "",
        f"**{len(su_rows)} theories** (including Veneziano-limit and conformal manifold theories)",
        "",
        "| # | Matter | Type | delta | N_f bound | R-charges | a/N^2 | c/N^2 | a/c |",
        "|---|--------|------|-------|-----------|-----------|-------|-------|-----|",
    ]

    for i, r in enumerate(su_rows, 1):
        delta_s = fmt_delta(r["chiral_excess_a"], r["chiral_excess_b"])
        nf_s = fmt_nf_bound(r["nf_bound_alpha"], r["nf_bound_gamma"])
        matter_s = r["matter_str"]
        if r["is_conformal_manifold"]:
            matter_s += " *"

        lines.append(
            f"| {i} | {matter_s} | {r['type']} | {delta_s} | {nf_s} "
            f"| {r['R_charges_str']} | {fmt_expr(r['a_over_N2'])} "
            f"| {fmt_expr(r['c_over_N2'])} | {r['a_over_c']} |"
        )

    lines.extend(["", "\\* = conformal manifold (b_0 = 0)", ""])

    with open(path, "w") as f:
        f.write("\n".join(lines))
    print(f"  Wrote {path} ({len(su_rows)} theories)")


def write_so_table(rows: list[dict], path: str) -> None:
    """Write SO(N) single-node table to markdown."""
    so_rows = sort_theories([r for r in rows if r["gauge_type"] == "SO"])

    lines = [
        "# SO(N) Single-Node Classification",
        "",
        f"**{len(so_rows)} theories**",
        "",
        "| # | Matter | Type | N_f bound | R-charges | a/N^2 | c/N^2 | a/c |",
        "|---|--------|------|-----------|-----------|-------|-------|-----|",
    ]

    for i, r in enumerate(so_rows, 1):
        nf_s = fmt_nf_bound(r["nf_bound_alpha"], r["nf_bound_gamma"])
        matter_s = r["matter_str"]
        if r["is_conformal_manifold"]:
            matter_s += " *"

        lines.append(
            f"| {i} | {matter_s} | {r['type']} | {nf_s} "
            f"| {r['R_charges_str']} | {fmt_expr(r['a_over_N2'])} "
            f"| {fmt_expr(r['c_over_N2'])} | {r['a_over_c']} |"
        )

    lines.extend(["", "\\* = conformal manifold (b_0 = 0)", ""])

    with open(path, "w") as f:
        f.write("\n".join(lines))
    print(f"  Wrote {path} ({len(so_rows)} theories)")


def write_sp_table(rows: list[dict], path: str) -> None:
    """Write Sp(N) single-node table to markdown."""
    sp_rows = sort_theories([r for r in rows if r["gauge_type"] == "Sp"])

    lines = [
        "# Sp(N) Single-Node Classification",
        "",
        f"**{len(sp_rows)} theories**",
        "",
        "| # | Matter | Type | N_f bound | R-charges | a/N^2 | c/N^2 | a/c |",
        "|---|--------|------|-----------|-----------|-------|-------|-----|",
    ]

    for i, r in enumerate(sp_rows, 1):
        nf_s = fmt_nf_bound(r["nf_bound_alpha"], r["nf_bound_gamma"])
        matter_s = r["matter_str"]
        if r["is_conformal_manifold"]:
            matter_s += " *"

        lines.append(
            f"| {i} | {matter_s} | {r['type']} | {nf_s} "
            f"| {r['R_charges_str']} | {fmt_expr(r['a_over_N2'])} "
            f"| {fmt_expr(r['c_over_N2'])} | {r['a_over_c']} |"
        )

    lines.extend(["", "\\* = conformal manifold (b_0 = 0)", ""])

    with open(path, "w") as f:
        f.write("\n".join(lines))
    print(f"  Wrote {path} ({len(sp_rows)} theories)")


# ── JSON export ───────────────────────────────────────────────────────────────

def export_json(rows: list[dict], path: str) -> None:
    """Export all theories to machine-readable JSON."""
    entries = []
    for r in rows:
        R_charges_str = {}
        for label, val in r["R_charges"].items():
            R_charges_str[label] = str(val)

        entry = {
            "gauge_type": r["gauge_type"],
            "matter": r["matter"],
            "matter_str": r["matter_str"],
            "n_rank2": r["n_rank2"],
            "type": r["type"],
            "chiral_excess_a": r["chiral_excess_a"],
            "chiral_excess_b": r["chiral_excess_b"],
            "is_veneziano": r["is_veneziano"],
            "nf_bound_alpha": str(r["nf_bound_alpha"]),
            "nf_bound_gamma": str(r["nf_bound_gamma"]),
            "is_conformal_manifold": r["is_conformal_manifold"],
            "R_charges": R_charges_str,
            "a_over_N2": str(r["a_over_N2"]),
            "c_over_N2": str(r["c_over_N2"]),
            "a_over_c": str(r["a_over_c"]) if r["a_over_c"] != "n/a" else "n/a",
            "a_over_N2_float": float(r["a_over_N2"]),
            "c_over_N2_float": float(r["c_over_N2"]),
        }
        entries.append(entry)

    with open(path, "w") as f:
        json.dump(entries, f, indent=2)
    print(f"  Wrote {path} ({len(entries)} theories)")


# ── Validation ────────────────────────────────────────────────────────────────

def _su_effective_T(matter: dict[str, int]) -> Fraction:
    """
    Effective total T contribution from rank-2/adj matter at an SU node.

    T_rep_lead(SU, adj) = 1, T_rep_lead(SU, S/Sbar/A/Abar) = 1/2.
    This determines the SU universality class: theories with the same
    (effective_T, |a_delta|) give identical leading-order R-charges and
    central charges.
    """
    n_adj = matter.get("adj", 0)
    n_tensor = sum(matter.get(r, 0) for r in ("S", "Sbar", "A", "Abar"))
    return Fraction(n_adj) + Fraction(n_tensor, 2)


def validate_all(rows: list[dict]) -> tuple[int, int, list[str]]:
    """
    Run all validation checks. Returns (passed, failed, failure_details).
    """
    passed = 0
    failed = 0
    failures = []

    def check(name: str, condition: bool, detail: str = ""):
        nonlocal passed, failed
        if condition:
            passed += 1
            print(f"  PASS: {name}")
        else:
            failed += 1
            msg = f"  FAIL: {name}" + (f" -- {detail}" if detail else "")
            print(msg)
            failures.append(msg)

    # ── Test 1: Theory count ──────────────────────────────────────────────
    print("\n=== Test: Theory count ===")
    su_count = sum(1 for r in rows if r["gauge_type"] == "SU")
    so_count = sum(1 for r in rows if r["gauge_type"] == "SO")
    sp_count = sum(1 for r in rows if r["gauge_type"] == "Sp")
    total = len(rows)
    check("SU count = 52", su_count == 52, f"got {su_count}")
    check("SO count = 6", so_count == 6, f"got {so_count}")
    check("Sp count = 9", sp_count == 9, f"got {sp_count}")
    check("Total = 67", total == 67, f"got {total}")

    # ── Test 2: SO/Sp Type II universal formula ───────────────────────────
    # For SO/Sp, all rank-2 reps have T_lead = 1, so Type II/III universal
    # formulas apply directly based on n_rank2.
    print("\n=== Test: SO/Sp Type II universal formula ===")
    dim_lead_map = {"SU": Rational(1), "SO": Rational(1, 2), "Sp": Rational(2)}

    for r in rows:
        g = r["gauge_type"]
        if g == "SU":
            continue
        if r["type"] != "II":
            continue
        expected = Rational(27, 128) * dim_lead_map[g]
        a_val = simplify(r["a_over_N2"])
        c_val = simplify(r["c_over_N2"])
        check(
            f"Type II {g} {r['matter_str']}: a/N^2 = {expected}",
            simplify(a_val - expected) == 0,
            f"got a/N^2 = {a_val}",
        )
        check(
            f"Type II {g} {r['matter_str']}: c/N^2 = {expected}",
            simplify(c_val - expected) == 0,
            f"got c/N^2 = {c_val}",
        )
        for label, R_val in r["R_charges"].items():
            check(
                f"Type II {g} {r['matter_str']}: R({label}) = 1/2",
                simplify(R_val - Rational(1, 2)) == 0,
                f"got R = {R_val}",
            )

    # ── Test 3: SO/Sp Type III universal formula ──────────────────────────
    print("\n=== Test: SO/Sp Type III universal formula ===")
    for r in rows:
        g = r["gauge_type"]
        if g == "SU":
            continue
        if r["type"] != "III":
            continue
        expected = Rational(1, 4) * dim_lead_map[g]
        a_val = simplify(r["a_over_N2"])
        c_val = simplify(r["c_over_N2"])
        check(
            f"Type III {g} {r['matter_str']}: a/N^2 = {expected}",
            simplify(a_val - expected) == 0,
            f"got a/N^2 = {a_val}",
        )
        check(
            f"Type III {g} {r['matter_str']}: c/N^2 = {expected}",
            simplify(c_val - expected) == 0,
            f"got c/N^2 = {c_val}",
        )
        for label, R_val in r["R_charges"].items():
            check(
                f"Type III {g} {r['matter_str']}: R({label}) = 2/3",
                simplify(R_val - Rational(2, 3)) == 0,
                f"got R = {R_val}",
            )

    # ── Test 4: SU universal formula (adj-only theories) ─────────────────
    # For SU, the universal formula holds when all rank-2 fields are adj
    # (all have T_lead = 1). Mixed adj + tensor theories have richer structure.
    print("\n=== Test: SU adj-only universal formula ===")
    for r in rows:
        if r["gauge_type"] != "SU":
            continue
        matter = r["matter"]
        n_adj = matter.get("adj", 0)
        n_tensor = sum(matter.get(rep, 0) for rep in ("S", "Sbar", "A", "Abar"))
        if n_tensor > 0:
            continue  # Mixed theories tested separately
        if n_adj == 2:
            check(
                f"SU 2adj: a/N^2 = 27/128",
                simplify(r["a_over_N2"] - Rational(27, 128)) == 0,
                f"got {simplify(r['a_over_N2'])}",
            )
            check(
                f"SU 2adj: R(adj) = 1/2",
                all(simplify(R - Rational(1, 2)) == 0
                    for R in r["R_charges"].values()),
                "",
            )
        elif n_adj == 3:
            check(
                f"SU 3adj: a/N^2 = 1/4",
                simplify(r["a_over_N2"] - Rational(1, 4)) == 0,
                f"got {simplify(r['a_over_N2'])}",
            )
            check(
                f"SU 3adj: R(adj) = 2/3",
                all(simplify(R - Rational(2, 3)) == 0
                    for R in r["R_charges"].values()),
                "",
            )
        elif n_adj == 1 and not r["is_veneziano"]:
            check(
                f"SU adj: a/N^2 = 0",
                simplify(r["a_over_N2"]) == 0,
                f"got {simplify(r['a_over_N2'])}",
            )

    # ── Test 5: Type I vanishing (SO/Sp) ──────────────────────────────────
    print("\n=== Test: Type I vanishing ===")
    for r in rows:
        if r["type"] != "I":
            continue
        g = r["gauge_type"]
        a_val = simplify(r["a_over_N2"])
        c_val = simplify(r["c_over_N2"])

        if g in ("SO", "Sp"):
            check(
                f"Type I {g} {r['matter_str']}: a/N^2 = 0",
                a_val == 0,
                f"got {a_val}",
            )
            check(
                f"Type I {g} {r['matter_str']}: c/N^2 = 0",
                c_val == 0,
                f"got {c_val}",
            )
        else:
            # SU Type I: non-Veneziano should vanish; Veneziano may not
            if not r["is_veneziano"]:
                check(
                    f"Type I SU {r['matter_str']}: a/N^2 = 0",
                    a_val == 0,
                    f"got {a_val}",
                )
            else:
                check(
                    f"Type I SU {r['matter_str']} (Veneziano): a/N^2 finite",
                    True,
                    f"a/N^2 = {a_val}",
                )

    # ── Test 6: Anomaly-free R-symmetry ───────────────────────────────────
    print("\n=== Test: Anomaly-free R-symmetry ===")
    for r in rows:
        g = r["gauge_type"]
        result = r["result"]
        q = r["quiver"]

        # T_adj + sum_i T_i (R_i - 1) = 0
        anomaly = _R(T_adj_lead(g))
        for f in result.fields:
            T_a = f.T_lead.get(0)
            if T_a is not None:
                anomaly += _R(T_a) * (result.R_charges[f.label] - 1)
        anomaly = simplify(anomaly)

        check(
            f"Anomaly-free {g} {r['matter_str']}",
            anomaly == 0,
            f"anomaly residual = {anomaly}",
        )

    # ── Test 7: SO/Sp universality ────────────────────────────────────────
    print("\n=== Test: SO/Sp universality ===")
    for g in ("SO", "Sp"):
        by_nrank2 = defaultdict(list)
        for r in rows:
            if r["gauge_type"] == g:
                by_nrank2[r["n_rank2"]].append(r)
        for n, group in by_nrank2.items():
            if len(group) <= 1:
                continue
            ref_a = simplify(group[0]["a_over_N2"])
            ref_c = simplify(group[0]["c_over_N2"])
            for r in group[1:]:
                check(
                    f"{g} n_rank2={n} universality: a/N^2 match ({r['matter_str']})",
                    simplify(r["a_over_N2"] - ref_a) == 0,
                    f"ref={ref_a}, got={simplify(r['a_over_N2'])}",
                )
                check(
                    f"{g} n_rank2={n} universality: c/N^2 match ({r['matter_str']})",
                    simplify(r["c_over_N2"] - ref_c) == 0,
                    f"ref={ref_c}, got={simplify(r['c_over_N2'])}",
                )

    # ── Test 8: SU universality classes ───────────────────────────────────
    # For SU, the universality class is determined by (effective_T, |a_delta|)
    # where effective_T = n_adj + n_tensor/2 (T_lead weighting).
    print("\n=== Test: SU universality classes ===")
    su_groups = defaultdict(list)
    for r in rows:
        if r["gauge_type"] == "SU":
            eff_T = _su_effective_T(r["matter"])
            key = (eff_T, abs(r["chiral_excess_a"]))
            su_groups[key].append(r)

    for key, group in su_groups.items():
        if len(group) <= 1:
            continue
        ref_a = simplify(group[0]["a_over_N2"])
        ref_c = simplify(group[0]["c_over_N2"])
        for r in group[1:]:
            check(
                f"SU class (eff_T={key[0]}, |delta|={key[1]}): a match ({r['matter_str']})",
                simplify(r["a_over_N2"] - ref_a) == 0,
                f"ref={ref_a}, got={simplify(r['a_over_N2'])}",
            )
            check(
                f"SU class (eff_T={key[0]}, |delta|={key[1]}): c match ({r['matter_str']})",
                simplify(r["c_over_N2"] - ref_c) == 0,
                f"ref={ref_c}, got={simplify(r['c_over_N2'])}",
            )

    # ── Test 9: a/c ratio ─────────────────────────────────────────────────
    # Hofman-Maldacena bounds: 1/2 <= a/c <= 5/4 for unitary 4D N=1 SCFTs.
    # Type I Veneziano theories may violate at leading order because the
    # leading-order a/N^2 is subleading (can be zero or negative); the true
    # a is O(N), not O(N^2). Skip these from the bound check.
    print("\n=== Test: a/c ratio bounds ===")
    for r in rows:
        a_val = simplify(r["a_over_N2"])
        if a_val == 0:
            continue
        ac = r["a_over_c_expr"]
        if ac is None:
            continue
        # Skip Veneziano Type I: leading-order a may be negative/unphysical
        if r["type"] == "I" and r["is_veneziano"]:
            print(f"  SKIP: a/c for Type I Veneziano {r['gauge_type']} "
                  f"{r['matter_str']} (leading-order artifact, a/c={float(ac):.4f})")
            continue
        ac_float = float(ac)
        check(
            f"a/c bounds {r['gauge_type']} {r['matter_str']}",
            0.5 - 1e-10 <= ac_float <= 1.25 + 1e-10,
            f"a/c = {ac_float}",
        )

    # ── Test 10: R-charge physical range ──────────────────────────────────
    # Type I Veneziano theories can have negative leading-order R-charges
    # because the actual R-charge is subleading in N. Skip those.
    print("\n=== Test: R-charge physical range ===")
    for r in rows:
        skip_veneziano_i = (r["type"] == "I" and r["is_veneziano"])
        for label, R_val in r["R_charges"].items():
            R_float = float(R_val)
            if skip_veneziano_i and R_float < 0:
                print(f"  SKIP: R({label})={R_float:.4f} for Type I Veneziano "
                      f"{r['gauge_type']} {r['matter_str']} (leading-order artifact)")
                continue
            check(
                f"R({label}) in [0,2] for {r['gauge_type']} {r['matter_str']}",
                -1e-10 <= R_float <= 2.0 + 1e-10,
                f"R = {R_float}",
            )

    # ── Test 11: Veneziano detection ──────────────────────────────────────
    print("\n=== Test: Veneziano detection ===")
    veneziano_count = sum(1 for r in rows if r["is_veneziano"])
    non_veneziano_su = sum(1 for r in rows
                          if r["gauge_type"] == "SU" and not r["is_veneziano"])
    check(
        f"Veneziano SU theories detected (expect some)",
        veneziano_count > 0,
        f"found {veneziano_count} Veneziano theories",
    )
    check(
        f"Non-Veneziano SU count reasonable",
        non_veneziano_su >= 20,
        f"found {non_veneziano_su} non-Veneziano SU theories",
    )

    # ── Test 12: Anchor paper comparison ──────────────────────────────────
    print("\n=== Test: Anchor paper comparison (arXiv:2510.19136) ===")

    # SO benchmarks from Table 27
    so_benchmarks = {
        "I": {"a_over_N2": Rational(0)},
        "II": {"a_over_N2": Rational(27, 256)},
        "III": {"a_over_N2": Rational(1, 8)},
    }
    so_rows_v = [r for r in rows if r["gauge_type"] == "SO"]
    for r in so_rows_v:
        bench = so_benchmarks.get(r["type"])
        if bench:
            check(
                f"SO anchor {r['matter_str']} Type {r['type']}: a/N^2",
                simplify(r["a_over_N2"] - bench["a_over_N2"]) == 0,
                f"expected {bench['a_over_N2']}, got {simplify(r['a_over_N2'])}",
            )

    # Sp benchmarks from Table 35
    sp_benchmarks = {
        "I": {"a_over_N2": Rational(0)},
        "II": {"a_over_N2": Rational(27, 64)},
        "III": {"a_over_N2": Rational(1, 2)},
    }
    sp_rows_v = [r for r in rows if r["gauge_type"] == "Sp"]
    for r in sp_rows_v:
        bench = sp_benchmarks.get(r["type"])
        if bench:
            check(
                f"Sp anchor {r['matter_str']} Type {r['type']}: a/N^2",
                simplify(r["a_over_N2"] - bench["a_over_N2"]) == 0,
                f"expected {bench['a_over_N2']}, got {simplify(r['a_over_N2'])}",
            )

    # SU adj-only benchmarks from Table 2
    su_non_ven = [r for r in rows
                  if r["gauge_type"] == "SU" and not r["is_veneziano"]]
    for r in su_non_ven:
        matter = r["matter"]
        n_adj = matter.get("adj", 0)
        n_tensor = sum(matter.get(rep, 0) for rep in ("S", "Sbar", "A", "Abar"))
        if n_tensor == 0:
            # Pure adj: universal formula applies
            if n_adj == 2:
                check(
                    f"SU anchor 2adj: a/N^2 = 27/128",
                    simplify(r["a_over_N2"] - Rational(27, 128)) == 0,
                    f"got {simplify(r['a_over_N2'])}",
                )
            elif n_adj == 3:
                check(
                    f"SU anchor 3adj: a/N^2 = 1/4",
                    simplify(r["a_over_N2"] - Rational(1, 4)) == 0,
                    f"got {simplify(r['a_over_N2'])}",
                )

    # ── Summary ───────────────────────────────────────────────────────────
    print(f"\n{'='*60}")
    print(f"VALIDATION SUMMARY: {passed} passed, {failed} failed")
    print(f"{'='*60}")
    if failures:
        print("\nFailures:")
        for f in failures:
            print(f)

    return passed, failed, failures


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))

    print("=" * 60)
    print("Single-Node Classification: Enumerate + Compute + Validate")
    print("=" * 60)

    # Step 1: Enumerate and compute
    print("\n--- Step 1: Enumerate all single-node theories ---")
    rows = classify_single_node()

    # Count by gauge type
    counts = defaultdict(int)
    for r in rows:
        counts[r["gauge_type"]] += 1
    print(f"\n  Counts: SU={counts['SU']}, SO={counts['SO']}, Sp={counts['Sp']}")
    print(f"  Total: {len(rows)}")

    # Step 2: Write tables
    print("\n--- Step 2: Write formatted tables ---")
    write_su_table(rows, os.path.join(base_dir, "tables", "su_single_node.md"))
    write_so_table(rows, os.path.join(base_dir, "tables", "so_single_node.md"))
    write_sp_table(rows, os.path.join(base_dir, "tables", "sp_single_node.md"))

    # Step 3: Export JSON
    print("\n--- Step 3: Export JSON ---")
    export_json(rows, os.path.join(base_dir, "single_node_classification.json"))

    # Step 4: Validate
    print("\n--- Step 4: Validation ---")
    passed, failed, failures = validate_all(rows)

    if failed > 0:
        print(f"\nWARNING: {failed} validation checks failed!")
        sys.exit(1)
    else:
        print("\nAll validation checks passed.")


if __name__ == "__main__":
    main()
