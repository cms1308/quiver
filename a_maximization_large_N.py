"""
a_maximization_large_N.py — Module 3 (large N limit, exact)

Exact leading-order R-charges and central charges in the large N limit
via symbolic a-maximization using sympy.

At large N (fixed N_f), dividing anomaly constraints by N and traces by N^2
gives N-independent rational coefficients. The anomaly constraints are linear
in R_i; after solving, the a-function is a cubic polynomial in the remaining
free parameters. Its critical points satisfy quadratic equations, yielding
exact closed-form R-charges (rational or algebraic numbers involving sqrt).

Reference: module3_a_maximization.md
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

from sympy import (
    Rational, Symbol, symbols, Matrix, solve, diff,
    simplify, sqrt, oo, S,
)
from sympy.core.expr import Expr

from beta_functions import GaugeType
from quiver_generation import Quiver, chiral_excess_coeffs, quiver_summary, enumerate_quivers


# ── Fraction → sympy Rational ─────────────────────────────────────────────────

def _R(f: Fraction | int) -> Rational:
    """Convert a fractions.Fraction (or int) to sympy Rational."""
    if isinstance(f, int):
        return Rational(f)
    return Rational(f.numerator, f.denominator)


# ── Leading-order coefficient functions ───────────────────────────────────────

def dim_group_lead(gauge_type: GaugeType) -> Fraction:
    """Leading N^2 coefficient of dim(G): lim_{N→∞} dim(G)/N^2."""
    if gauge_type == "SU":
        return Fraction(1)       # (N^2-1)/N^2 → 1
    if gauge_type == "SO":
        return Fraction(1, 2)    # N(N-1)/2 / N^2 → 1/2
    if gauge_type == "Sp":
        return Fraction(2)       # N(2N+1)/N^2 → 2
    raise ValueError(f"Unknown gauge type: {gauge_type!r}")


def T_adj_lead(gauge_type: GaugeType) -> Fraction:
    """Leading N coefficient of T(adj): lim_{N→∞} T(adj)/N. Always 1."""
    if gauge_type in ("SU", "SO", "Sp"):
        return Fraction(1)
    raise ValueError(f"Unknown gauge type: {gauge_type!r}")


def T_rep_lead(gauge_type: GaugeType, rep: str) -> Fraction:
    """Leading N coefficient of T(rep): lim_{N→∞} T(rep)/N."""
    if gauge_type == "SU":
        if rep in ("fund", "antifund"):
            return Fraction(0)
        if rep == "adj":
            return Fraction(1)
        if rep in ("S", "Sbar"):
            return Fraction(1, 2)
        if rep in ("A", "Abar"):
            return Fraction(1, 2)
    elif gauge_type == "SO":
        if rep == "V":
            return Fraction(0)
        if rep == "adj":
            return Fraction(1)
        if rep == "S":
            return Fraction(1)
    elif gauge_type == "Sp":
        if rep == "fund":
            return Fraction(0)
        if rep == "adj":
            return Fraction(1)
        if rep == "A":
            return Fraction(1)
    raise ValueError(f"Unknown rep {rep!r} for gauge type {gauge_type!r}")


def dim_rep_lead(gauge_type: GaugeType, rep: str) -> Fraction:
    """Leading N^2 coefficient of dim(rep): lim_{N→∞} dim(rep)/N^2."""
    if gauge_type == "SU":
        if rep in ("fund", "antifund"):
            return Fraction(0)
        if rep == "adj":
            return Fraction(1)
        if rep in ("S", "Sbar"):
            return Fraction(1, 2)
        if rep in ("A", "Abar"):
            return Fraction(1, 2)
    elif gauge_type == "SO":
        if rep == "V":
            return Fraction(0)
        if rep == "adj":
            return Fraction(1, 2)
        if rep == "S":
            return Fraction(1, 2)
    elif gauge_type == "Sp":
        if rep == "fund":
            return Fraction(0)
        if rep == "adj":
            return Fraction(2)
        if rep == "A":
            return Fraction(2)
    raise ValueError(f"Unknown rep {rep!r} for gauge type {gauge_type!r}")


def T_bifund_lead(ga: GaugeType, gb: GaugeType) -> tuple[Fraction, Fraction]:
    """Leading N coefficient of T_bifund at each endpoint: lim_{N→∞} T/N."""
    pair = (ga, gb)
    if pair == ("SU", "SU"):
        return Fraction(1, 2), Fraction(1, 2)
    if pair == ("SU", "SO"):
        return Fraction(1, 2), Fraction(1)
    if pair == ("SO", "SU"):
        return Fraction(1), Fraction(1, 2)
    if pair == ("SU", "Sp"):
        return Fraction(1), Fraction(1, 2)
    if pair == ("Sp", "SU"):
        return Fraction(1, 2), Fraction(1)
    if pair == ("SO", "SO"):
        return Fraction(1), Fraction(1)
    if pair == ("SO", "Sp"):
        return Fraction(2), Fraction(1, 2)
    if pair == ("Sp", "SO"):
        return Fraction(1, 2), Fraction(2)
    if pair == ("Sp", "Sp"):
        return Fraction(1), Fraction(1)
    raise ValueError(f"Unknown gauge pair: ({ga!r}, {gb!r})")


def dim_bifund_lead(ga: GaugeType, gb: GaugeType) -> Fraction:
    """Leading N^2 coefficient of dim(bifund) = dim_src*dim_dst / N^2."""
    pair = (ga, gb)
    if pair == ("SU", "SU"):
        return Fraction(1)
    if pair in (("SU", "SO"), ("SO", "SU")):
        return Fraction(1)
    if pair in (("SU", "Sp"), ("Sp", "SU")):
        return Fraction(2)
    if pair == ("SO", "SO"):
        return Fraction(1)
    if pair in (("SO", "Sp"), ("Sp", "SO")):
        return Fraction(2)
    if pair == ("Sp", "Sp"):
        return Fraction(4)
    raise ValueError(f"Unknown gauge pair: ({ga!r}, {gb!r})")


# ── Leading-order field descriptor ─────────────────────────────────────────────

@dataclass
class LeadField:
    """One chiral multiplet type at leading order in N."""
    label: str                           # e.g. "node0_adj", "edge_0_1_pm"
    R_index: int                         # index into R-charge vector
    dim_lead: Fraction                   # coefficient of N^2 in dim
    T_lead: dict[int, Fraction]          # node → coefficient of N in T


# ── Build leading-order fields ─────────────────────────────────────────────────

def build_fields_large_N(quiver: Quiver) -> list[LeadField]:
    """
    Build the leading-order chiral fields for large N.

    No N or N_f parameter: fund-like fields with fixed N_f are subleading.
    Fund-like matter appears only when chiral excess has nonzero N-coefficient
    (anomaly cancellation forces O(N) extra fundamentals).

    rank_multipliers[i] = m_i means node i has gauge group G(m_i * N).
    All leading-order coefficients are scaled accordingly:
      - T_adj_lead(g_a) * m_a    (T_adj ~ m_a * N)
      - T_rep_lead(g, rep) * m_a  (T_rep ~ m_a * N for rank-2/adj)
      - dim_rep_lead(g, rep) * m_a²  (dim ~ m_a² * N²)
      - T_bifund_lead at node a from node b: T_bifund_lead(g_a, g_b) * m_b
      - dim_bifund_lead(g_a, g_b) * m_a * m_b
    """
    fields: list[LeadField] = []
    idx = 0
    mults = quiver.rank_multipliers

    # 1. Rank-2/adj node matter
    for i, (g, matter) in enumerate(zip(quiver.gauge_types, quiver.node_matter)):
        m_i = mults[i]
        for rep, count in matter.items():
            if count == 0:
                continue
            t = count * T_rep_lead(g, rep) * m_i
            d = count * dim_rep_lead(g, rep) * m_i * m_i
            if t == 0 and d == 0:
                continue
            fields.append(LeadField(
                label=f"node{i}_{rep}",
                R_index=idx,
                dim_lead=d,
                T_lead={i: t},
            ))
            idx += 1

    # 2. Fund-like matter at SU nodes with nonzero chiral excess coefficient
    for i, g in enumerate(quiver.gauge_types):
        if g != "SU":
            continue
        m_i = mults[i]
        a_delta, _ = chiral_excess_coeffs(quiver, i)
        if a_delta == 0:
            continue
        # |a_delta|*N extra fundamentals of SU(m_i*N): dim = |a_delta|*m_i per N²,
        # T = |a_delta|/2 per N (T_fund = 1/2 exact)
        dim_lead = Fraction(abs(a_delta)) * m_i
        T_lead_val = Fraction(abs(a_delta), 2)
        if a_delta > 0:
            fields.append(LeadField(
                label=f"node{i}_fund",
                R_index=idx,
                dim_lead=dim_lead,
                T_lead={i: T_lead_val},
            ))
        else:
            fields.append(LeadField(
                label=f"node{i}_antifund",
                R_index=idx,
                dim_lead=dim_lead,
                T_lead={i: T_lead_val},
            ))
        idx += 1

    # 3. Bifundamental edges — merge identical (src, dst, rep) groups.
    # Identical edges must share the same R-charge at an SCFT, so treating
    # them as one field with count×weights is the physically correct choice.
    # It also keeps n_free minimal, avoiding spurious symmetry-breaking
    # critical points in the NSolve polynomial system.
    from collections import Counter
    edge_counts: dict[tuple, int] = {}
    seen_edge_keys: list[tuple] = []
    for e in quiver.edges:
        key = (e.src, e.dst, e.rep)
        if key not in edge_counts:
            edge_counts[key] = 0
            seen_edge_keys.append(key)
        edge_counts[key] += 1
    for key in seen_edge_keys:
        count = edge_counts[key]
        i, j, rep = key
        gi, gj = quiver.gauge_types[i], quiver.gauge_types[j]
        m_i, m_j = mults[i], mults[j]
        T_i_base, T_j_base = T_bifund_lead(gi, gj)
        T_i = T_i_base * m_j * count   # T at node i scales with neighbor's rank mult
        T_j = T_j_base * m_i * count
        d = dim_bifund_lead(gi, gj) * m_i * m_j * count
        rep_label = rep.replace("+", "p").replace("-", "m")
        fields.append(LeadField(
            label=f"edge_{i}_{j}_{rep_label}",
            R_index=idx,
            dim_lead=d,
            T_lead={i: T_i, j: T_j},
        ))
        idx += 1

    return fields


# ── Symbolic anomaly constraints ──────────────────────────────────────────────

def _build_anomaly_system(
    fields: list[LeadField],
    quiver: Quiver,
    R: list[Symbol],
) -> list[Expr]:
    """
    Build anomaly-free constraint equations (one per gauge node).

    Each equation: T_adj_lead * m_a + Σ T_lead(field) * (R_field - 1) = 0
    where m_a = rank_multipliers[a].
    """
    mults = quiver.rank_multipliers
    eqs: list[Expr] = []
    for a, g in enumerate(quiver.gauge_types):
        m_a = mults[a]
        eq: Expr = _R(T_adj_lead(g)) * m_a
        for f in fields:
            T_a = f.T_lead.get(a)
            if T_a is not None:
                eq = eq + _R(T_a) * (R[f.R_index] - 1)
        eqs.append(eq)
    return eqs


# ── Symbolic a and c functions ────────────────────────────────────────────────

def _symbolic_traces(
    fields: list[LeadField],
    quiver: Quiver,
    R: list[Expr],
) -> tuple[Expr, Expr]:
    """Return symbolic (Tr R / N^2, Tr R^3 / N^2) at leading order.
    Gaugino contribution scales as dim_group_lead(g) * m_a²."""
    mults = quiver.rank_multipliers
    gaugino = sum(_R(dim_group_lead(g)) * m * m
                  for g, m in zip(quiver.gauge_types, mults))
    tr_R: Expr = gaugino
    tr_R3: Expr = gaugino
    for f in fields:
        d = _R(f.dim_lead)
        r = R[f.R_index]
        tr_R = tr_R + d * (r - 1)
        tr_R3 = tr_R3 + d * (r - 1) ** 3
    return tr_R, tr_R3


def _symbolic_a(
    fields: list[LeadField],
    quiver: Quiver,
    R: list[Expr],
) -> Expr:
    """Symbolic a/N^2 = (3/32)(3 Tr R^3 - Tr R)."""
    tr_R, tr_R3 = _symbolic_traces(fields, quiver, R)
    return Rational(3, 32) * (3 * tr_R3 - tr_R)


def _symbolic_c(
    fields: list[LeadField],
    quiver: Quiver,
    R: list[Expr],
) -> Expr:
    """Symbolic c/N^2 = (1/32)(9 Tr R^3 - 5 Tr R)."""
    tr_R, tr_R3 = _symbolic_traces(fields, quiver, R)
    return Rational(1, 32) * (9 * tr_R3 - 5 * tr_R)


# ── Gauge-invariant operator check ────────────────────────────────────────────

def _gauge_invariant_ops_symbolic(
    fields: list[LeadField],
    quiver: Quiver,
    R: list[Expr],
) -> dict[str, Expr]:
    """
    Return gauge-invariant operators and their symbolic R-charges.

    Operators: meson (fund+antifund), bilinears (adj Tr Φ², S·Sbar, A·Abar,
    V·V, Sp f·Ω·f), bifundamental mesons (Q·Qtilde for +- pairs).
    """
    by_label: dict[str, LeadField] = {f.label: f for f in fields}
    ops: dict[str, Expr] = {}

    for i, g in enumerate(quiver.gauge_types):
        fl, afl = f"node{i}_fund", f"node{i}_antifund"
        if g == "SU" and fl in by_label and afl in by_label:
            ops[f"meson_{i}"] = R[by_label[fl].R_index] + R[by_label[afl].R_index]

        bilinear_pairs = [
            ("adj",  "adj",  f"adj_bilinear_{i}"),
            ("S",    "Sbar", f"SSbar_{i}"),
            ("A",    "Abar", f"AAbar_{i}"),
            ("V",    "V",    f"VV_{i}"),
        ]
        if g == "Sp":
            bilinear_pairs.append(("fund", "fund", f"Sp_ff_{i}"))

        for rep1, rep2, op_label in bilinear_pairs:
            l1, l2 = f"node{i}_{rep1}", f"node{i}_{rep2}"
            if rep1 == rep2:
                if l1 in by_label:
                    ops[op_label] = 2 * R[by_label[l1].R_index]
            else:
                if l1 in by_label and l2 in by_label:
                    ops[op_label] = R[by_label[l1].R_index] + R[by_label[l2].R_index]

    # Bifundamental pairs: +- edges in opposite directions
    edge_fields = [f for f in fields if f.label.startswith("edge_")]
    for f1 in edge_fields:
        parts1 = f1.label.split("_")
        src1, dst1, rep1 = parts1[1], parts1[2], "_".join(parts1[3:])
        if rep1 != "pm":
            continue
        for f2 in edge_fields:
            parts2 = f2.label.split("_")
            src2, dst2, rep2 = parts2[1], parts2[2], "_".join(parts2[3:])
            if rep2 != "pm" or src2 != dst1 or dst2 != src1:
                continue
            key = f"bif_meson_{min(src1,src2)}_{max(src1,src2)}"
            if key not in ops:
                ops[key] = R[f1.R_index] + R[f2.R_index]

    return ops


# ── Result dataclass ───────────────────────────────────────────────────────────

@dataclass
class LargeNResult:
    """
    Exact result of leading-order large N a-maximization.

    R-charges and central charges are sympy expressions (exact algebraic numbers).
    """
    fields: list[LeadField]
    R_charges: dict[str, Expr]           # field.label → exact R-charge
    a_over_N2: Expr                      # exact a/N^2
    c_over_N2: Expr                      # exact c/N^2
    unitarity_ok: bool                   # all gauge-invariant ops have R ≥ 2/3
    gauge_invariant_R: dict[str, Expr] = field(default_factory=dict)  # op → R[O]

    def print_summary(self) -> None:
        """Print a human-readable summary with both exact and numeric values."""
        print("  R-charges (leading order in N):")
        for label, R in self.R_charges.items():
            print(f"    {label}: {R}  ≈ {float(R):.6f}")
        print(f"  a/N^2 = {self.a_over_N2}  ≈ {float(self.a_over_N2):.6f}")
        print(f"  c/N^2 = {self.c_over_N2}  ≈ {float(self.c_over_N2):.6f}")
        if self.gauge_invariant_R:
            print("  Gauge-invariant operators:")
            for op, R in self.gauge_invariant_R.items():
                print(f"    R[{op}] = {R}  ≈ {float(R):.6f}")
        print(f"  unitarity_ok = {self.unitarity_ok}")


# ── Core exact a-maximization ─────────────────────────────────────────────────

def a_maximize_large_N(quiver: Quiver) -> LargeNResult:
    """
    Exact leading-order a-maximization at large N using sympy.

    Algorithm:
    1. Build leading-order fields (T_lead, dim_lead as exact rationals)
    2. Write anomaly constraints as symbolic linear equations in R_i
    3. Solve constraints: express dependent R_i in terms of free parameters
    4. Substitute into a_trial (cubic polynomial in free parameters)
    5. Differentiate, solve quadratic critical-point equations exactly
    6. Select the maximum among critical points
    7. Return exact R-charges, a/N^2, c/N^2 as sympy expressions
    """
    fields = build_fields_large_N(quiver)
    n_R = len(fields)

    if n_R == 0:
        # No leading-order matter: gaugino-only contribution.
        mults = quiver.rank_multipliers
        gaugino_sum = sum(_R(dim_group_lead(g)) * m * m
                          for g, m in zip(quiver.gauge_types, mults))
        a_val = Rational(3, 32) * (2 * gaugino_sum)
        c_val = Rational(1, 32) * (4 * gaugino_sum)
        return LargeNResult(
            fields=[], R_charges={}, a_over_N2=a_val, c_over_N2=c_val,
            unitarity_ok=True,
        )

    return _solve_a_max(fields, quiver)


def a_maximize_large_N_with_decoupling(quiver: Quiver) -> LargeNResult:
    """
    Exact leading-order a-maximization with iterative unitarity decoupling.

    When a gauge-invariant operator O has R[O] < 2/3, add the constraint
    Σ R[fields in O] = 2/3 and redo. Iterate until all R[O] ≥ 2/3.
    """
    fields = build_fields_large_N(quiver)
    n_R = len(fields)

    if n_R == 0:
        return a_maximize_large_N(quiver)

    extra_constraints: list[tuple[list[int], Rational]] = []  # (R_indices, rhs)

    for _ in range(n_R + 1):
        result = _solve_a_max(fields, quiver, extra_constraints)

        violating = [
            (op, r_expr)
            for op, r_expr in result.gauge_invariant_R.items()
            if simplify(r_expr - Rational(2, 3)) < 0
        ]
        if not violating:
            return result

        already = {
            tuple(sorted(idxs)) for idxs, _ in extra_constraints
        }
        added = False
        for op, _ in violating:
            # Find R-indices for this operator from the symbolic expression
            r_indices = _op_r_indices(op, fields, quiver)
            key = tuple(sorted(set(r_indices)))
            if key in already:
                continue
            extra_constraints.append((r_indices, Rational(2, 3)))
            already.add(key)
            added = True
        if not added:
            break

    return result


def _op_r_indices(op_label: str, fields: list[LeadField], quiver: Quiver) -> list[int]:
    """Extract the R_indices involved in a gauge-invariant operator."""
    by_label = {f.label: f for f in fields}

    if op_label.startswith("meson_"):
        node = int(op_label.split("_")[1])
        return [by_label[f"node{node}_fund"].R_index,
                by_label[f"node{node}_antifund"].R_index]
    if op_label.startswith("adj_bilinear_"):
        node = int(op_label.split("_")[2])
        idx = by_label[f"node{node}_adj"].R_index
        return [idx, idx]
    if op_label.startswith("SSbar_"):
        node = int(op_label.split("_")[1])
        return [by_label[f"node{node}_S"].R_index,
                by_label[f"node{node}_Sbar"].R_index]
    if op_label.startswith("AAbar_"):
        node = int(op_label.split("_")[1])
        return [by_label[f"node{node}_A"].R_index,
                by_label[f"node{node}_Abar"].R_index]
    if op_label.startswith("VV_"):
        node = int(op_label.split("_")[1])
        idx = by_label[f"node{node}_V"].R_index
        return [idx, idx]
    if op_label.startswith("Sp_ff_"):
        node = int(op_label.split("_")[2])
        idx = by_label[f"node{node}_fund"].R_index
        return [idx, idx]
    if op_label.startswith("bif_meson_"):
        parts = op_label.split("_")
        s, d = parts[2], parts[3]
        edge_fields = [f for f in fields if f.label.startswith("edge_")]
        idxs = []
        for f in edge_fields:
            fp = f.label.split("_")
            if fp[3:] == ["pm"] and ((fp[1] == s and fp[2] == d) or
                                      (fp[1] == d and fp[2] == s)):
                idxs.append(f.R_index)
        return idxs
    return []


def _solve_a_max(
    fields: list[LeadField],
    quiver: Quiver,
    extra_constraints: list[tuple[list[int], Rational]] | None = None,
) -> LargeNResult:
    """
    Core symbolic solver: solve anomaly constraints + optional extra constraints,
    parameterize free directions, maximize a symbolically.
    """
    n_R = len(fields)
    R_syms = symbols(f"R:{n_R}", real=True)

    # Build anomaly constraint equations
    eqs = _build_anomaly_system(fields, quiver, list(R_syms))

    # Add extra (decoupling) constraints: Σ R[k] = rhs
    if extra_constraints:
        for r_indices, rhs in extra_constraints:
            eq = -rhs
            for k in r_indices:
                eq = eq + R_syms[k]
            eqs.append(eq)

    # Build constraint matrix A @ R = b using sympy Matrix
    A_rows = []
    b_vals = []
    for eq in eqs:
        row = [Rational(0)] * n_R
        const = Rational(0)
        # eq is linear in R_syms: eq = Σ a_i R_i + c = 0
        for i, Ri in enumerate(R_syms):
            coeff = eq.coeff(Ri)
            row[i] = coeff
        const = eq.subs({Ri: 0 for Ri in R_syms})
        A_rows.append(row)
        b_vals.append(-const)

    A = Matrix(A_rows)
    b = Matrix(b_vals)

    # Null space of A
    ns_vecs = A.nullspace()
    n_free = len(ns_vecs)

    # Particular solution: A^T (A A^T)^{-1} b (minimum-norm)
    AAT = A * A.T
    if AAT.det() != 0:
        R0 = A.T * AAT.inv() * b
    else:
        # Fall back: row-reduce augmented matrix
        aug = A.row_join(b)
        rref_mat, pivots = aug.rref()
        R0 = Matrix([Rational(0)] * n_R)
        for row_i, p in enumerate(pivots):
            if p < n_R:
                R0[p] = rref_mat[row_i, n_R]

    # Parameterize: R = R0 + Σ s_i * ns_i
    if n_free > 0:
        s_syms = symbols(f"s:{n_free}", real=True)
        perturbation = Matrix([Rational(0)] * n_R)
        for i in range(n_free):
            perturbation = perturbation + s_syms[i] * ns_vecs[i]
        R_param = R0 + perturbation
    else:
        s_syms = ()
        R_param = R0

    R_list = [R_param[i] for i in range(n_R)]

    # Build symbolic a-function
    a_expr = _symbolic_a(fields, quiver, R_list)

    if n_free > 0:
        # Differentiate and solve
        grad = [diff(a_expr, si) for si in s_syms]
        solutions = solve(grad, s_syms, dict=True)

        if solutions:
            # Evaluate a at each critical point, pick the maximum
            # Use float comparison to identify the maximum
            # Skip complex solutions (arise from spurious branches)
            best_sol = None
            best_a = None
            for sol in solutions:
                a_val = simplify(a_expr.subs(sol))
                try:
                    a_float = float(a_val)
                except (TypeError, ValueError):
                    continue  # complex or non-numeric solution
                if best_a is None or a_float > float(best_a):
                    best_a = a_val
                    best_sol = sol

            R_opt = [simplify(r.subs(best_sol)) for r in R_list]
            a_opt = simplify(best_a)
        else:
            # No critical points — use R0 (boundary or degenerate case)
            R_opt = [simplify(R0[i]) for i in range(n_R)]
            a_opt = simplify(a_expr.subs({si: 0 for si in s_syms}))
    else:
        R_opt = [simplify(R0[i]) for i in range(n_R)]
        a_opt = simplify(a_expr)

    # Compute c at the optimum
    c_opt = simplify(_symbolic_c(fields, quiver, R_opt))

    # Build result
    R_charges = {f.label: R_opt[f.R_index] for f in fields}

    # Check gauge-invariant operators
    gi_ops = _gauge_invariant_ops_symbolic(fields, quiver, R_opt)
    gi_ops_simplified = {op: simplify(r) for op, r in gi_ops.items()}
    unitarity_ok = all(
        simplify(r - Rational(2, 3)) >= 0
        for r in gi_ops_simplified.values()
    )

    return LargeNResult(
        fields=fields,
        R_charges=R_charges,
        a_over_N2=a_opt,
        c_over_N2=c_opt,
        unitarity_ok=unitarity_ok,
        gauge_invariant_R=gi_ops_simplified,
    )


# ── Numerical result container and anomaly matrix ────────────────────────────

@dataclass
class FastNumericalResult:
    """Numerical result of large-N a-maximization (from Mathematica NSolve)."""
    a_over_N2: float
    c_over_N2: float
    R_charges: dict[str, float]   # field label → R-charge


def _anomaly_matrix_exact(fields, quiver) -> tuple[list[list], list]:
    """
    Build the anomaly matrix A and RHS vector b as exact Fraction values.

    Anomaly-free condition per gauge node a:
        Σ_i T_lead_a(field_i) × (R_i − 1) = −T_adj_lead(g_a) × m_a

    Rewritten as A·R = b with:
        A[a, i] = T_lead_a(field_i)
        b[a]    = −T_adj_lead(g_a)×m_a + Σ_i T_lead_a(field_i)
    """
    from fractions import Fraction
    mults = quiver.rank_multipliers
    n_nodes = quiver.n_nodes
    n_R = len(fields)
    A: list[list] = [[Fraction(0)] * n_R for _ in range(n_nodes)]
    b: list       = [Fraction(0)] * n_nodes
    for a, g in enumerate(quiver.gauge_types):
        m_a = mults[a]
        b[a] = -T_adj_lead(g) * m_a
        for f in fields:
            T_a = f.T_lead.get(a)
            if T_a is not None:
                A[a][f.R_index] += T_a
                b[a] += T_a
    return A, b


# ── Mathematica batch NSolve scanner ─────────────────────────────────────────

def a_maximize_batch_mathematica(
    quivers: list[Quiver],
    timeout: int = 7200,
    working_precision: int = 30,
) -> list[FastNumericalResult | None]:
    """
    Numerical large-N a-maximization for a list of quivers via a single
    Mathematica NSolve session.

    Finds ALL real critical points of the a-function (stationarity equations are
    quadratic in the free parameters s), keeps only local maxima (negative-
    semidefinite Hessian), and returns the maximum.  Returns None for theories
    with no bounded maximum (unbounded a-function or no real critical points).

    A single wolframscript process handles the entire batch, avoiding per-theory
    subprocess startup overhead.
    """
    import json
    import subprocess
    import tempfile
    import os
    from fractions import Fraction as Frac

    def _frac_to_mma(f: Frac) -> str:
        return f"{f.numerator}/{f.denominator}"

    def _fvec_to_mma(v: list) -> str:
        return "{" + ",".join(_frac_to_mma(x) for x in v) + "}"

    def _fmat_to_mma(M: list) -> str:
        return "{" + ",".join(_fvec_to_mma(row) for row in M) + "}"

    # ── Serialise each quiver's data ─────────────────────────────────────────
    theory_data = []
    for idx, q in enumerate(quivers):
        fields = build_fields_large_N(q)
        n_R = len(fields)
        mults = q.rank_multipliers
        gaugino = sum(
            dim_group_lead(g) * m * m
            for g, m in zip(q.gauge_types, mults)
        )

        if n_R == 0:
            # No leading-order matter: a and c are fixed by gaugino alone
            a_val = Frac(3, 32) * 2 * gaugino
            c_val = Frac(1, 32) * 4 * gaugino
            theory_data.append({
                "idx": idx,
                "trivial": True,
                "a": _frac_to_mma(a_val),
                "c": _frac_to_mma(c_val),
            })
            continue

        # Exact anomaly matrix and RHS (all Fraction entries)
        A, b = _anomaly_matrix_exact(fields, q)
        dims = [f.dim_lead for f in fields]  # Fraction values

        theory_data.append({
            "idx":         idx,
            "trivial":     False,
            "A_mma":       _fmat_to_mma(A),
            "b_mma":       _fvec_to_mma(b),
            "dims_mma":    _fvec_to_mma(dims),
            "gaugino_mma": _frac_to_mma(gaugino),
            "labels":      [f.label for f in fields],
        })

    # ── Write Mathematica batch script (parallel) ─────────────────────────────
    tmpdir = tempfile.mkdtemp()
    script_path = os.path.join(tmpdir, "mma_batch.wl")
    result_path = os.path.join(tmpdir, "mma_results.json")

    wp = working_precision

    # Build the Mathematica list of Association literals, one per theory
    assoc_parts = []
    for td in theory_data:
        idx = td["idx"]
        if td.get("trivial"):
            assoc_parts.append(
                f'<|"idx"->{idx},"trivial"->True,"a0"->{td["a"]},"c0"->{td["c"]}|>'
            )
        else:
            assoc_parts.append(
                f'<|"idx"->{idx},"trivial"->False,'
                f'"A"->{td["A_mma"]},'
                f'"b"->{td["b_mma"]},'
                f'"dims"->{td["dims_mma"]},'
                f'"gaugino"->{td["gaugino_mma"]}|>'
            )

    theories_literal = "{\n" + ",\n".join(assoc_parts) + "\n}"

    script = f"""\
(* Parallel a-maximization via NSolve — anomaly constraint enforced exactly *)
LaunchKernels[];

processTheory[th_] := Module[
  {{idx, Amat, bvec, Fnull, nfree, Fmat, R0, dims, gaugino,
    svars, R, r1, trR, trR3, a, c, grad, hess, sols, solsSym,
    freeParams, avals, best, Ropt}},
  idx = th["idx"];
  If[th["trivial"],
    Return[<|"idx" -> idx, "a" -> th["a0"], "c" -> th["c0"], "R" -> {{}}|>]
  ];
  (* Exact anomaly constraint: A.R = b  (all entries are exact rationals) *)
  Amat  = th["A"];
  bvec  = th["b"];
  Fnull = NullSpace[Amat];      (* exact rational null vectors, shape n_free × n_R *)
  nfree = Length[Fnull];
  Fmat  = Transpose[Fnull];    (* n_R × n_free, so R = R0 + Fmat.svars *)
  R0    = LinearSolve[Amat, bvec];  (* exact particular solution *)
  If[!ListQ[R0],
    Return[<|"idx" -> idx, "a" -> "NONE"|>]
  ];
  dims    = th["dims"];
  gaugino = th["gaugino"];
  If[nfree == 0,
    (* R fully determined by anomaly constraint — no free parameters *)
    r1   = R0 - 1;
    trR  = gaugino + dims.r1;
    trR3 = gaugino + dims.(r1^3);
    Return[<|"idx" -> idx,
             "a" -> N[(3/32)*(3*trR3 - trR), {wp}],
             "c" -> N[(1/32)*(9*trR3 - 5*trR), {wp}],
             "R" -> N[R0, {wp}]|>]
  ];
  (* nfree > 0: maximize a over the anomaly-free subspace *)
  svars = Array[s, nfree];
  R     = R0 + Fmat.svars;    (* satisfies A.R = b for all svars by construction *)
  r1    = R - 1;
  trR   = gaugino + dims.r1;
  trR3  = gaugino + dims.(r1^3);
  a     = (3/32)*(3*trR3 - trR);
  c     = (1/32)*(9*trR3 - 5*trR);
  grad  = Table[D[a, svars[[j]]], {{j, nfree}}];
  hess  = Outer[D[a, #1, #2]&, svars, svars];
  sols  = Quiet[NSolve[Thread[grad == 0], svars, Reals,
                       WorkingPrecision -> {wp}]];
  If[!ListQ[sols], sols = {{}}];
  (* Keep only local maxima: Hessian is negative definite. *)
  sols  = Select[sols,
    NegativeDefiniteMatrixQ[N[hess /. #, {wp}]] &];
  If[Length[sols] == 0,
    (* Fallback: symbolic Solve — handles flat directions NSolve misses.
       Use complex Solve (much faster than Reals domain for quadratic systems);
       any remaining free parameters (flat directions) are pinned to 0. Then
       filter out complex solutions and accept negative-semidefinite Hessians
       (zero eigenvalues come from flat directions). *)
    solsSym = TimeConstrained[
      Quiet[Solve[Thread[grad == 0], svars]],
      120, {{}}];
    If[ListQ[solsSym] && Length[solsSym] > 0,
      (* Residual variables after substitution: unsolved svars (flat directions
         within svars) AND external parameters like C[_] from Solve. *)
      freeParams = Variables[svars /. solsSym];
      If[Length[freeParams] > 0,
        solsSym = (Join[#, Thread[freeParams -> 0]]&) /@ solsSym];
      solsSym = Select[solsSym,
        Function[rules,
          AllTrue[N[svars /. rules, {wp}],
            NumericQ[#] && Abs[Im[#]] < 10^-10 &]]];
      sols = Select[solsSym,
        Function[rules,
          NegativeSemidefiniteMatrixQ[N[hess /. rules, {wp}]]]];
    ];
  ];
  If[Length[sols] == 0,
    Return[<|"idx" -> idx, "a" -> "NONE"|>]
  ];
  avals = N[a /. sols, {wp}];
  best  = First[Ordering[avals, -1]];
  Ropt  = N[R /. sols[[best]], {wp}];
  <|"idx" -> idx,
    "a" -> avals[[best]],
    "c" -> N[c /. sols[[best]], {wp}],
    "R" -> Ropt|>
];

theories = {theories_literal};
results  = ParallelMap[processTheory, theories];
Export["{result_path}", results, "JSON"];
"""

    with open(script_path, "w") as fh:
        fh.write(script)

    # ── Run wolframscript ─────────────────────────────────────────────────────
    subprocess.run(
        ["wolframscript", "-file", script_path],
        timeout=timeout,
        check=True,
    )

    # ── Parse results ─────────────────────────────────────────────────────────
    # Mathematica can emit e.g. "0.e-30" (no digit after decimal) which is
    # invalid JSON.  Fix by inserting a zero: "0.e-30" → "0.0e-30".
    import re as _re
    with open(result_path) as fh:
        raw_text = fh.read()
    # Mathematica emits bare decimals like "0.", "0.e7", "1.e-30" — all invalid
    # JSON.  Any digit followed by "." not followed by another digit gets a "0"
    # appended, e.g. "0.e7" → "0.0e7", "0." → "0.0".
    raw_text = _re.sub(r'(\d)\.(?!\d)', r'\1.0', raw_text)
    raw = json.loads(raw_text)

    by_idx: dict[int, dict] = {int(r["idx"]): r for r in raw}

    out: list[FastNumericalResult | None] = []
    for idx, q in enumerate(quivers):
        r = by_idx.get(idx)
        if r is None or r.get("a") == "NONE":
            out.append(None)
            continue
        a_val = float(r["a"])
        c_val = float(r["c"])
        R_list = [float(x) for x in r["R"]]
        fields = build_fields_large_N(q)
        R_charges = {f.label: R_list[f.R_index] for f in fields} if R_list else {}
        out.append(FastNumericalResult(
            a_over_N2=a_val, c_over_N2=c_val, R_charges=R_charges,
        ))
    return out


# ── Classification table ──────────────────────────────────────────────────────

def _fmt_matter(matter: dict[str, int]) -> str:
    """Compact matter notation: '2adj + S + Ā' etc."""
    rep_sym = {
        "adj": "adj", "S": "S", "Sbar": "S̄",
        "A": "A", "Abar": "Ā", "V": "V",
        "fund": "□", "antifund": "□̄",
    }
    order = ["adj", "S", "Sbar", "A", "Abar", "V", "fund", "antifund"]
    parts = []
    for rep in order:
        n = matter.get(rep, 0)
        if n:
            sym = rep_sym.get(rep, rep)
            parts.append(f"{n}{sym}" if n > 1 else sym)
    return " + ".join(parts) if parts else "—"


def _fmt_R(label: str) -> str:
    """Compact field-name from a ChiralField label."""
    # "node0_adj" → "adj", "node0_antifund" → "□̄", "edge_0_1_pm" → "Q₀₁"
    rep_sym = {
        "adj": "adj", "S": "S", "Sbar": "S̄",
        "A": "A", "Abar": "Ā", "V": "V",
        "fund": "□", "antifund": "□̄",
    }
    if label.startswith("node"):
        _, rep = label.split("_", 1)
        return rep_sym.get(rep, rep)
    return label


def _fmt_expr(expr) -> str:
    """Compact string for a sympy expression."""
    from sympy import Rational as Rat, nsimplify, sqrt
    s = str(expr)
    # Make it a bit shorter for common forms
    s = s.replace("sqrt(", "√(").replace("**2", "²").replace("**3", "³")
    return s


def _fmt_linear(a: int, b: int) -> str:
    """Format an integer linear expression aN+b as e.g. 'N-4', '2N', '-N+3', '0'."""
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
            parts.append(f"+{b}" if b > 0 else str(b))
        else:
            parts.append(str(b))
    return "".join(parts) if parts else "0"


def _fmt_af(bounds: list) -> str:
    """
    Format the AF upper bound on N_f.

    bounds is a list of (alpha, gamma) Fraction tuples such that
    AF holds for all nodes iff N_f < alpha*N + gamma at each node.
    For a single-node theory there is one element; take the tightest.
    """
    if not bounds:
        return "?"
    alpha, gamma = bounds[0]   # single-node: one constraint
    parts = []
    if alpha:
        parts.append(f"{alpha}N" if alpha != 1 else "N")
    if gamma or not parts:
        g = int(gamma)
        if g >= 0:
            parts.append(f"+{g}" if parts else str(g))
        else:
            parts.append(str(g))
    return "".join(parts)


def _classify_single_node() -> None:
    """
    Print a classification table of all single-node quiver theories
    by their large N leading-order R-charges and central charges.
    """
    from quiver_generation import enumerate_quivers
    from sympy import Rational as Rat

    # Pure SYM a/N² per gauge type (gaugino only, R=1 always)
    _a_sym = {
        "SU": Rational(3, 16),   # (3/32)*2*1
        "SO": Rational(3, 32),   # (3/32)*2*(1/2)
        "Sp": Rational(3, 8),    # (3/32)*2*2
    }

    quivers = enumerate_quivers(
        n_nodes=1, max_multiedge=1, min_multiedge=1, require_connected=True
    )

    # Collect results, skipping pure SYM (no matter)
    rows = []
    for q, bounds in quivers:
        g = q.gauge_types[0]
        matter = q.node_matter[0]
        if not any(matter.values()):
            continue  # skip pure SYM

        res = a_maximize_large_N(q)

        # AF bound
        af_str = _fmt_af(bounds)

        # R-charges
        r_strs = [
            f"R_{_fmt_R(lbl)}={_fmt_expr(R)}"
            for lbl, R in res.R_charges.items()
        ]
        r_col = ",  ".join(r_strs) if r_strs else "(gaugino only)"

        a_val = float(res.a_over_N2)

        # Chiral excess: full aN+b form from anomaly cancellation (SU nodes only)
        if g == "SU":
            a_delta, b_delta = chiral_excess_coeffs(q, 0)
            delta_str = _fmt_linear(a_delta, b_delta)
        else:
            delta_str = "—"

        rows.append({
            "g":        g,
            "matter":   _fmt_matter(matter),
            "af":       af_str,
            "R":        r_col,
            "a":        _fmt_expr(res.a_over_N2),
            "a_float":  a_val,
            "a_num":    f"{a_val:.6f}",
            "delta":    delta_str,
        })

    # Print one sub-table per gauge group, sorted by a/N²
    for g_type in ("SU", "SO", "Sp"):
        sub = [r for r in rows if r["g"] == g_type]
        sub.sort(key=lambda r: r["a_float"])

        col_w = {
            "matter": max(len(r["matter"]) for r in sub),
            "delta":  max(len(r["delta"])  for r in sub),
            "af":     max(len(r["af"])     for r in sub),
            "R":      max(len(r["R"])      for r in sub),
            "a":      max(len(r["a"])      for r in sub),
            "a_num":  max(len(r["a_num"])  for r in sub),
        }
        col_w = {k: max(v, len(k)) for k, v in col_w.items()}

        title = f"  {g_type}(N) — single-node theories  ({len(sub)} theories)"
        print()
        print(title)
        print("  " + "─" * (len(title) - 2))

        hdr = (
            f"  {'Matter':<{col_w['matter']}}  "
            f"{'delta':>{col_w['delta']}}  "
            f"{'AF: N_f <':<{col_w['af']}}  "
            f"{'R-charges (large N)':<{col_w['R']}}  "
            f"{'a/N²':<{col_w['a']}}  "
            f"{'a/N² (num)':<{col_w['a_num']}}"
        )
        print(hdr)
        print("  " + "─" * (len(hdr) - 2))

        for r in sub:
            line = (
                f"  {r['matter']:<{col_w['matter']}}  "
                f"{r['delta']:>{col_w['delta']}}  "
                f"{r['af']:<{col_w['af']}}  "
                f"{r['R']:<{col_w['R']}}  "
                f"{r['a']:<{col_w['a']}}  "
                f"{r['a_num']:<{col_w['a_num']}}"
            )
            print(line)


def _fmt_edges(quiver: Quiver) -> str:
    """Compact edge description: '2(+-)  (++)'  etc."""
    from collections import Counter
    c: Counter = Counter()
    for e in quiver.edges:
        c[(e.src, e.dst, e.rep)] += 1
    parts = []
    for (src, dst, rep), n in sorted(c.items()):
        mult = f"{n}×" if n > 1 else ""
        parts.append(f"{mult}({src}→{dst},{rep})")
    return "  ".join(parts) if parts else "—"


def _classify_two_node() -> None:
    """
    Print a classification table of all two-node quiver theories
    by their large N a/N² value (Mathematica NSolve scan).
    """
    import time
    import numpy as np

    print("Enumerating 2-node quivers...", flush=True)
    t0 = time.time()
    quivers_bounds = enumerate_quivers(
        n_nodes=2, max_multiedge=4, min_multiedge=1, require_connected=True
    )
    print(f"  {len(quivers_bounds)} quivers found in {time.time()-t0:.1f}s", flush=True)

    print("Running large-N a-maximization (Mathematica NSolve)...", flush=True)
    t1 = time.time()
    mma_results = a_maximize_batch_mathematica([q for q, _ in quivers_bounds])
    rows = []
    n_none = 0
    for (q, bounds), fast_res in zip(quivers_bounds, mma_results):
        if fast_res is None:
            n_none += 1
            continue
        g0, g1 = q.gauge_types[0], q.gauge_types[1]
        rows.append({
            "gauge_pair": (g0, g1),
            "matter0":    _fmt_matter(q.node_matter[0]),
            "matter1":    _fmt_matter(q.node_matter[1]),
            "edges":      _fmt_edges(q),
            "a":          fast_res.a_over_N2,
        })
    print(f"  Done in {time.time()-t1:.1f}s. "
          f"{n_none} no-SCFT (diverged), {len(rows)} physical theories.", flush=True)

    # Cluster into universality classes: same a/N² up to tolerance 1e-5
    TOL = 1e-5

    def _cluster(vals: list[float]) -> list[tuple[float, list[int]]]:
        """Return (representative_a, [indices]) clusters sorted by a."""
        if not vals:
            return []
        order = np.argsort(vals)
        clusters: list[tuple[float, list[int]]] = []
        cur_center = vals[order[0]]
        cur_members = [order[0]]
        for idx in order[1:]:
            if abs(vals[idx] - cur_center) < TOL:
                cur_members.append(idx)
                cur_center = np.mean([vals[i] for i in cur_members])
            else:
                clusters.append((cur_center, cur_members))
                cur_center = vals[idx]
                cur_members = [idx]
        clusters.append((cur_center, cur_members))
        return sorted(clusters, key=lambda x: x[0])

    # Group by gauge pair
    pair_order = [("SU", "SU"), ("SO", "SU"), ("SU", "Sp"),
                  ("SO", "SO"), ("SO", "Sp"), ("Sp", "Sp")]
    pair_labels = {
        ("SU", "SU"): "SU(N)×SU(N)",
        ("SO", "SU"): "SO(N)×SU(N)",
        ("SU", "Sp"): "SU(N)×Sp(N)",
        ("SO", "SO"): "SO(N)×SO(N)",
        ("SO", "Sp"): "SO(N)×Sp(N)",
        ("Sp", "Sp"): "Sp(N)×Sp(N)",
    }

    total_classes = 0
    for pair in pair_order:
        sub = [r for r in rows if r["gauge_pair"] == pair]
        if not sub:
            continue
        a_vals = [r["a"] for r in sub]
        clusters = _cluster(a_vals)
        total_classes += len(clusters)

        title = f"  {pair_labels[pair]}  ({len(sub)} theories, {len(clusters)} classes)"
        print()
        print(title)
        print("  " + "─" * (len(title) - 2))
        print(f"  {'a/N²':>9}  {'# theories':>11}  Representative theory")
        print("  " + "─" * 80)

        for a_center, members in clusters:
            n_theories = len(members)
            # Pick a representative: the one with fewest edges, then fewest matter
            rep_idx = min(
                members,
                key=lambda i: (len(sub[i]["edges"]), sub[i]["matter0"], sub[i]["matter1"]),
            )
            rep = sub[rep_idx]
            m0 = rep["matter0"] or "—"
            m1 = rep["matter1"] or "—"
            e  = rep["edges"]
            theory_str = f"node0: {m0}  |  node1: {m1}  |  edges: {e}"
            print(f"  {a_center:>9.6f}  {n_theories:>11}  {theory_str}")

    print()
    print(f"  Total universality classes (|a/N²| ≤ {max_a}): {total_classes}")


# ── Validation ─────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    from quiver_generation import Edge
    from a_maximization import a_maximize

    def _convergence_check(quiver: Quiver, N_values: list[int] = [50, 100, 200]) -> None:
        """Show convergence of numeric a/N^2 to the exact large N result."""
        for N in N_values:
            res = a_maximize(quiver, N=N, N_f=0)
            print(f"    N={N:>4}:  a/N^2 = {res.a/N**2:.6f},  c/N^2 = {res.c/N**2:.6f}")

    # # ── Example 1: SU(N) + 1 adj ──────────────────────────────────────────────
    # # Anomaly: 1 + 1*(R-1) = 0  ⟹  R = 0  (fully determined, no free params)
    # print("=== SU(N) + 1 adj ===")
    # q1 = Quiver(["SU"], node_matter=[{"adj": 1}])
    # res1 = a_maximize_large_N(q1)
    # res1.print_summary()
    # print("  Numeric convergence:")
    # _convergence_check(q1)

    # print()

    # # ── Example 2: SU(N)^2 circular (balanced +-) ─────────────────────────────
    # # Anomaly at each node: 1 + ½(R01-1) + ½(R10-1) = 0  ⟹  R01 + R10 = 0
    # # a-max by symmetry: R01 = R10 = 0
    # print("=== SU(N)^2 circular (balanced +-) ===")
    # q2 = Quiver(["SU", "SU"], [Edge(0, 1, "+-"), Edge(1, 0, "+-")])
    # res2 = a_maximize_large_N(q2)
    # res2.print_summary()
    # print("  Numeric convergence:")
    # _convergence_check(q2)

    # print()

    # # ── Example 3: SU(N) + 1 S + 1 Sbar ──────────────────────────────────────
    # # Chiral excess: n_S - n_Sbar = 0, so a_delta = 0. No fund-like fields.
    # # Anomaly: 1 + ½(R_S-1) + ½(R_Sbar-1) = 0  ⟹  R_S + R_Sbar = 0
    # # Free direction R_S - R_Sbar: a-max gives R_S = R_Sbar = 0
    # print("=== SU(N) + 1 S + 1 Sbar ===")
    # q3 = Quiver(["SU"], node_matter=[{"S": 1, "Sbar": 1}])
    # res3 = a_maximize_large_N(q3)
    # res3.print_summary()
    # print("  Numeric convergence:")
    # _convergence_check(q3)

    # print()

    # # ── Example 4: SU(N) + 1 S (chiral) ──────────────────────────────────────
    # # Chiral excess: a_delta = -1 (need ~N antifundamentals)
    # # Fields: S (dim_lead=1/2, T_lead=1/2) + antifund (dim_lead=1, T_lead=1/2)
    # # Anomaly: 1 + ½(R_S-1) + ½(R_af-1) = 0  ⟹  R_S + R_af = 0
    # # One free param: nontrivial a-max with exact algebraic solution
    # print("=== SU(N) + 1 S (chiral) ===")
    # q4 = Quiver(["SU"], node_matter=[{"S": 1}])
    # res4 = a_maximize_large_N(q4)
    # res4.print_summary()
    # print("  Numeric convergence:")
    # _convergence_check(q4)

    # print()

    # # ── Example 5: SU(N) + 1 adj + 1 S + 1 Sbar ─────────────────────────────
    # # Three fields, one constraint ⟹ two free parameters
    # print("=== SU(N) + 1 adj + 1 S + 1 Sbar ===")
    # q5 = Quiver(["SU"], node_matter=[{"adj": 1, "S": 1, "Sbar": 1}])
    # res5 = a_maximize_large_N(q5)
    # res5.print_summary()
    # print("  Numeric convergence:")
    # _convergence_check(q5)


    # ── Classification table: all single-node theories ───────────────────────
    # _classify_single_node()

    # ── Classification table: all two-node theories ───────────────────────────
    _classify_two_node()
