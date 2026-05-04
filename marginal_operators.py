"""
marginal_operators.py — Module 6

Enumerate gauge-invariant single-trace chiral operators in two-node quiver
gauge theories, evaluate their R-charges at finite N via numerical
a-maximization, and identify operators that are marginal (R = 2) at every N
in a chosen test set, restricted to flavor singlets.

A flavor-singlet operator with R = 2 across multiple N is a candidate for an
exactly marginal deformation in the Leigh–Strassler sense: turning it on
preserves the flavor symmetry, so the only running comes from gauge couplings.
Marginality robust to 1/N corrections (rather than coincidental at one N) is a
necessary condition for a genuine conformal-manifold direction.

Conventions follow CLAUDE.md: SU(N) reps {□, □̄, adj, S, S̄, A, Ā}; SO(N)
{V, adj, S}; Sp(N) {f, adj, A}; bifundamentals only between distinct nodes.
"""

from __future__ import annotations

import re
import unicodedata
from collections import Counter
from dataclasses import dataclass, field
from itertools import product

import numpy as np
from scipy.linalg import null_space
from scipy.optimize import minimize as _minimize

from quiver_generation import Edge, Quiver, _nf_bound_at_N
from a_maximization import (
    a_maximize, build_fields, ChiralField,
    anomaly_matrix, a_trial, dim_rep,
)
from beta_functions import T_rep


# ── Field/half-edge bookkeeping ────────────────────────────────────────────────

# Half-edge counts (upper, lower) per node-rep label
_NODE_REP_BALANCE = {
    "adj":      (1, 1),
    "S":        (2, 0),
    "Sbar":     (0, 2),
    "A":        (2, 0),
    "Abar":     (0, 2),
    "fund":     (1, 0),
    "antifund": (0, 1),
    "V":        (1, 1),    # SO vector — δ contracts upper to lower
}


def _edge_label(e: Edge) -> str:
    """Match a_maximization.build_fields label format."""
    rep = e.rep.replace("+", "p").replace("-", "m")
    return f"edge_{e.src}_{e.dst}_{rep}"


def _half_edge_balance(quiver: Quiver, node: int, label: str) -> tuple[int, int]:
    """Return (upper, lower) index-slot count contributed at `node` by one
    copy of the field referenced by `label`. SU nodes use upper/lower
    distinction (mesonic singlet ↔ upper == lower). SO/Sp nodes use total
    count (singlet ↔ even, since δ/Ω contracts any pair); we report the
    total in the upper slot so callers can sum upper+lower uniformly."""
    g = quiver.gauge_types[node]
    if label.startswith("node"):
        prefix, rep = label.split("_", 1)
        if int(prefix[4:]) != node:
            return (0, 0)
        if g == "SU":
            return _NODE_REP_BALANCE.get(rep, (0, 0))
        if rep == "V":          # SO vector: 1 index
            return (1, 0)
        if rep == "adj":        # 2 indices
            return (2, 0)
        if rep in ("S", "A"):   # 2 indices (symmetric or antisymmetric)
            return (2, 0)
        if rep == "fund":       # Sp fundamental: 1 index
            return (1, 0)
        return (0, 0)
    parts = label.split("_")
    src, dst, rep = int(parts[1]), int(parts[2]), "_".join(parts[3:])
    if node not in (src, dst):
        return (0, 0)
    side = "src" if node == src else "dst"
    if g == "SU":
        if rep == "pm":
            return (1, 0) if side == "src" else (0, 1)
        if rep == "pp":
            return (1, 0)
        if rep == "mm":
            return (0, 1)
        if rep == "p":          # SU-side fund of SU-SO/Sp edge
            return (1, 0)
        if rep == "m":          # SU-side antifund
            return (0, 1)
    # SO or Sp endpoint of any edge: 1 vector/fundamental index
    return (1, 0)


def _index_balance_node(
    quiver: Quiver, node: int, multiset: dict[str, int]
) -> bool:
    """Gauge-invariance proxy at leading-N: SU nodes need upper == lower
    (mesonic), SO/Sp nodes need (upper + lower) even (δ/Ω contractions)."""
    upper = lower = 0
    for label, mult in multiset.items():
        u, l = _half_edge_balance(quiver, node, label)
        upper += mult * u
        lower += mult * l
    g = quiver.gauge_types[node]
    if g == "SU":
        return upper == lower
    return (upper + lower) % 2 == 0


# ── Flavor multiplets ──────────────────────────────────────────────────────────

def _flavor_multiplets(quiver: Quiver) -> dict[str, int]:
    """Map multiplet_id → multiplicity n (size of U(n) flavor group acting).

    Multiplets:
      - node{i}_{rep} for each rank-2/adj rep with count ≥ 2
      - edge_{src}_{dst}_{rep} for each bifundamental group with multiplicity ≥ 2
    Multiplicity-1 fields contribute no flavor symmetry (no flavor charge to
    worry about).
    """
    mults: dict[str, int] = {}
    for i, matter in enumerate(quiver.node_matter):
        for rep, count in matter.items():
            if count >= 2:
                mults[f"node{i}_{rep}"] = count
    edge_groups: Counter = Counter()
    for e in quiver.edges:
        edge_groups[(e.src, e.dst, e.rep)] += 1
    for (src, dst, rep), n in edge_groups.items():
        if n >= 2:
            rep_lbl = rep.replace("+", "p").replace("-", "m")
            mults[f"edge_{src}_{dst}_{rep_lbl}"] = n
    return mults


def _flavor_signature(
    factors: tuple[tuple[str, int], ...], multiplets: dict[str, int]
) -> dict[str, int]:
    """Count copies of fields belonging to each flavor multiplet."""
    sig: dict[str, int] = {}
    for label, mult in factors:
        if label in multiplets:
            sig[label] = sig.get(label, 0) + mult
    return sig


# ── Candidate operator dataclass ──────────────────────────────────────────────

@dataclass(frozen=True)
class CandidateOp:
    factors: tuple[tuple[str, int], ...]      # canonical-sorted ((label, mult), ...)
    kind: str                                  # "single-node" | "bifund-loop" | "mixed"
    word: tuple[str, ...] | None               # canonical cyclic word for loops
    degree: int
    flavor_signature: dict[str, int] = field(default_factory=dict, compare=False, hash=False)


def _canonical_factors(multiset: dict[str, int]) -> tuple[tuple[str, int], ...]:
    return tuple(sorted((label, mult) for label, mult in multiset.items() if mult > 0))


# ── Enumerators ────────────────────────────────────────────────────────────────

def _multisets_up_to(labels: list[str], max_degree: int) -> list[dict[str, int]]:
    """All multisets over `labels` with 1 ≤ total degree ≤ max_degree."""
    if not labels:
        return []
    out: list[dict[str, int]] = []
    n = len(labels)

    def rec(i: int, remaining: int, cur: dict[str, int]) -> None:
        if i == n:
            if cur and sum(cur.values()) > 0:
                out.append(dict(cur))
            return
        for k in range(remaining + 1):
            if k > 0:
                cur[labels[i]] = k
            rec(i + 1, remaining - k, cur)
            if k > 0:
                del cur[labels[i]]

    rec(0, max_degree, {})
    return [m for m in out if sum(m.values()) >= 1]


_LABEL_PROBE_N = (10, 20, 30)


def _all_labels_via_build(quiver: Quiver) -> set[str]:
    """Union of build_fields labels across multiple N probes — captures
    fund/antifund identity changes when chiral excess flips sign with N."""
    out: set[str] = set()
    for N in _LABEL_PROBE_N:
        try:
            for f in build_fields(quiver, N=N, N_f=0):
                out.add(f.label)
        except Exception:
            pass
    return out


def _node_intra_labels(quiver: Quiver, node: int,
                       all_labels: set[str] | None = None) -> list[str]:
    """Labels for fields living entirely at `node`. Includes rank-2/adj from
    `node_matter` plus fund/antifund/V derived by `build_fields` when chiral
    excess or N_f forces them. Pass `all_labels` to override the discovery
    (e.g. include max-N_f fund/antifund/V)."""
    if all_labels is None:
        all_labels = _all_labels_via_build(quiver)
    out = []
    for label in all_labels:
        if not label.startswith("node"):
            continue
        if int(label[4:].split("_", 1)[0]) == node:
            out.append(label)
    return sorted(out)


def _all_field_labels(quiver: Quiver,
                      all_labels: set[str] | None = None) -> list[str]:
    """All field labels in the quiver (intra-node + edges)."""
    if all_labels is None:
        all_labels = _all_labels_via_build(quiver)
    return sorted(all_labels)


def _enumerate_single_node(quiver: Quiver, max_degree: int,
                           all_labels: set[str] | None = None) -> list[CandidateOp]:
    out: list[CandidateOp] = []
    multiplets = _flavor_multiplets(quiver)
    for i in range(quiver.n_nodes):
        labels = _node_intra_labels(quiver, i, all_labels=all_labels)
        if not labels:
            continue
        for ms in _multisets_up_to(labels, max_degree):
            if not _index_balance_node(quiver, i, ms):
                continue
            # Other nodes: no fields used → trivially balanced.
            ok = all(_index_balance_node(quiver, j, {}) for j in range(quiver.n_nodes) if j != i)
            if not ok:
                continue
            factors = _canonical_factors(ms)
            sig = _flavor_signature(factors, multiplets)
            out.append(CandidateOp(
                factors=factors, kind="single-node", word=None,
                degree=sum(ms.values()), flavor_signature=sig,
            ))
    return out


def _enumerate_bifund_loops(quiver: Quiver, max_degree: int) -> list[CandidateOp]:
    """Closed walks in the bifundamental multigraph.

    For SU-SU edges, orientation is dictated by `rep`: +- has 1 upper at src
    and 1 lower at dst; ++ has 2 upper; -- has 2 lower. SO/Sp endpoints
    contract via δ/Ω.

    We enumerate by trying each ordered cyclic sequence of edges (with possible
    repetition) of length 2..max_degree and check whether it gives a
    gauge-invariant operator via _index_balance_node at every node. Then dedup
    by canonical cyclic rotation + reversal.
    """
    edges = quiver.edges
    if not edges:
        return []
    out: dict[tuple, CandidateOp] = {}
    multiplets = _flavor_multiplets(quiver)

    edge_labels = [_edge_label(e) for e in edges]
    # The cyclic word uses edge indices as alphabet so that distinct copies of
    # identical bifundamentals yield the same canonical form (they share label).
    seen_canonical: set[tuple[str, ...]] = set()

    for length in range(2, max_degree + 1):
        for combo in product(range(len(edges)), repeat=length):
            ms: dict[str, int] = {}
            for idx in combo:
                lbl = edge_labels[idx]
                ms[lbl] = ms.get(lbl, 0) + 1
            # Check gauge invariance at every node simultaneously
            if not all(_index_balance_node(quiver, n, ms) for n in range(quiver.n_nodes)):
                continue
            word = tuple(edge_labels[i] for i in combo)
            canon = _canonical_cycle(word)
            if canon in seen_canonical:
                continue
            seen_canonical.add(canon)
            factors = _canonical_factors(ms)
            sig = _flavor_signature(factors, multiplets)
            out[(factors, canon)] = CandidateOp(
                factors=factors, kind="bifund-loop", word=canon,
                degree=length, flavor_signature=sig,
            )
    return list(out.values())


def _canonical_cycle(word: tuple[str, ...]) -> tuple[str, ...]:
    """Lex-minimum over rotations and reversal."""
    n = len(word)
    rotations = [word[i:] + word[:i] for i in range(n)]
    rev = tuple(reversed(word))
    rotations += [rev[i:] + rev[:i] for i in range(n)]
    return min(rotations)


def _enumerate_mixed(quiver: Quiver, max_degree: int,
                     all_labels: set[str] | None = None) -> list[CandidateOp]:
    """General enumerator: multisets drawn from intra-node + edge labels with
    gauge invariance at every node. Subsumes single-node and bifund-loop;
    deduplicated against those by (factors, kind=…)."""
    labels = _all_field_labels(quiver, all_labels=all_labels)
    if not labels:
        return []
    out: list[CandidateOp] = []
    multiplets = _flavor_multiplets(quiver)
    for ms in _multisets_up_to(labels, max_degree):
        if not all(_index_balance_node(quiver, n, ms) for n in range(quiver.n_nodes)):
            continue
        # Skip trivially single-node (already covered by _enumerate_single_node)
        edge_count = sum(m for lbl, m in ms.items() if lbl.startswith("edge_"))
        intra_nodes = {int(lbl[4]) for lbl in ms if lbl.startswith("node")}
        if edge_count == 0 and len(intra_nodes) <= 1:
            continue  # single-node, already enumerated
        if len(intra_nodes) == 0:
            continue  # pure bifund loop, already enumerated
        factors = _canonical_factors(ms)
        sig = _flavor_signature(factors, multiplets)
        out.append(CandidateOp(
            factors=factors, kind="mixed", word=None,
            degree=sum(ms.values()), flavor_signature=sig,
        ))
    return out


def enumerate_candidates(quiver: Quiver, max_degree: int = 6,
                         all_labels: set[str] | None = None) -> list[CandidateOp]:
    """Enumerate distinct gauge-invariant single-trace candidate operators
    of degree ≥ 2. Degree-1 trace operators are either identically zero
    (Tr T^a = 0 for all gauge generators) or non-singlet (any rank > 0
    representation needs index contraction).

    `all_labels` overrides automatic label discovery — pass the union of
    labels from a custom field-builder (e.g. `_all_labels_max_Nf`) to include
    fund/antifund/V/f added by N_f saturation."""
    cands: dict[tuple, CandidateOp] = {}
    for op in _enumerate_single_node(quiver, max_degree, all_labels=all_labels):
        if op.degree < 2:
            continue
        cands[(op.kind, op.factors, op.word)] = op
    for op in _enumerate_bifund_loops(quiver, max_degree):
        if op.degree < 2:
            continue
        cands[(op.kind, op.factors, op.word)] = op
    for op in _enumerate_mixed(quiver, max_degree, all_labels=all_labels):
        if op.degree < 2:
            continue
        cands[(op.kind, op.factors, op.word)] = op
    return list(cands.values())


def _enumerate_with_labels(quiver: Quiver, max_degree: int,
                           all_labels: set[str]) -> list[CandidateOp]:
    return enumerate_candidates(quiver, max_degree=max_degree, all_labels=all_labels)


# ── Flavor singlet check ──────────────────────────────────────────────────────

def is_flavor_singlet(op: CandidateOp, multiplets: dict[str, int]) -> bool:
    """Singlet iff every multiplet of size ≥ 2 is used with count ≠ 1.
    A lone field from an n ≥ 2 multiplet carries uncontracted flavor charge."""
    for mult_id, n in multiplets.items():
        if n < 2:
            continue
        if op.flavor_signature.get(mult_id, 0) == 1:
            return False
    return True


# ── R-charge evaluation at finite N ───────────────────────────────────────────

def r_values_at_N(quiver: Quiver, N: int, N_f: int = 0) -> dict[str, float]:
    """Run numerical a-maximization at this N and return label → R-charge."""
    res = a_maximize(quiver, N=N, N_f=N_f)
    return res.R_charges


def op_R_at_N(op: CandidateOp, R_values: dict[str, float]) -> float | None:
    """Sum constituent R-charges. Returns None if any field is missing
    (e.g. fund/antifund count is zero at this N)."""
    total = 0.0
    for label, mult in op.factors:
        if label not in R_values:
            return None
        total += mult * R_values[label]
    return total


def is_marginal_at_all_N(
    op: CandidateOp,
    R_per_N: dict[int, dict[str, float]],
    tol: float = 1e-6,
) -> bool:
    for R_values in R_per_N.values():
        r = op_R_at_N(op, R_values)
        if r is None or abs(r - 2.0) > tol:
            return False
    return True


def _r_values_sane(R: dict[str, float]) -> bool:
    """Reject obviously diverged a-max output (R outside a generous bound)."""
    return all(-0.5 <= v <= 2.5 for v in R.values())


def find_marginal_ops(
    quiver: Quiver,
    N_list: tuple[int, ...] = (10, 20, 30),
    max_degree: int = 6,
    tol: float = 1e-6,
    N_f: int = 0,
) -> list[CandidateOp]:
    """Enumerate candidates, evaluate R at each N, return ops with R=2 at
    every tested N. Returns [] if a-max fails to converge at any tested N
    (detected as R-charges outside the physical [0,2] range)."""
    R_per_N: dict[int, dict[str, float]] = {}
    for N in N_list:
        R = r_values_at_N(quiver, N, N_f=N_f)
        if not _r_values_sane(R):
            return []
        R_per_N[N] = R
    candidates = enumerate_candidates(quiver, max_degree=max_degree)
    return [op for op in candidates if is_marginal_at_all_N(op, R_per_N, tol=tol)]


def find_marginal_singlet_ops(
    quiver: Quiver,
    N_list: tuple[int, ...] = (10, 20, 30),
    max_degree: int = 6,
    tol: float = 1e-6,
    N_f: int = 0,
) -> list[CandidateOp]:
    """Same as find_marginal_ops but additionally filtered by the heuristic
    flavor-singlet rule (lone field from a multiplet of size ≥ 2 → not singlet).
    The check is incomplete (misses U(1) phases on multiplicity-1 fields and
    finer rep-theory of multi-field flavor tensors)."""
    R_per_N = {N: r_values_at_N(quiver, N, N_f=N_f) for N in N_list}
    candidates = enumerate_candidates(quiver, max_degree=max_degree)
    multiplets = _flavor_multiplets(quiver)
    return [
        op for op in candidates
        if is_marginal_at_all_N(op, R_per_N, tol=tol)
        and is_flavor_singlet(op, multiplets)
    ]


# ── Max-N_f mode: saturate b_0 = 0 at each node ───────────────────────────────

def nf_max_per_node(quiver: Quiver, N: int) -> list[int]:
    """Max N_f at each node such that b_0 ≤ 0 there. Floor of the AF bound;
    clamped at 0 (no negative N_f). For nodes already at b_0 = 0 with the
    anomaly-required matter alone, this returns 0."""
    out: list[int] = []
    for i in range(quiver.n_nodes):
        f = _nf_bound_at_N(quiver, i, N)
        out.append(max(0, int(f)))
    return out


def build_fields_max_Nf(
    quiver: Quiver, N: int, Nf_per_node: list[int]
) -> list[ChiralField]:
    """Build chiral fields with per-node N_f added on top of the anomaly-required
    fund/antifund/V/f matter. Mirrors a_maximization.build_fields but allows
    each gauge node to carry a different N_f (so each node can independently
    saturate b_0 = 0)."""
    fields = build_fields(quiver, N=N, N_f=0)
    by_label = {f.label: f for f in fields}
    next_idx = max((f.R_index for f in fields), default=-1) + 1
    mults = quiver.rank_multipliers

    def _bump(node: int, rep: str, count: int) -> None:
        nonlocal next_idx
        if count <= 0:
            return
        g = quiver.gauge_types[node]
        N_i = mults[node] * N
        d_extra = count * dim_rep(g, rep, N_i)
        T_extra = count * float(T_rep(g, rep, N_i))
        lbl = f"node{node}_{rep}"
        f = by_label.get(lbl)
        if f is None:
            f = ChiralField(label=lbl, R_index=next_idx, dim=d_extra,
                            T_contributions={node: T_extra})
            fields.append(f)
            by_label[lbl] = f
            next_idx += 1
        else:
            f.dim += d_extra
            f.T_contributions[node] = f.T_contributions.get(node, 0) + T_extra

    for i, g in enumerate(quiver.gauge_types):
        Nf = Nf_per_node[i]
        if Nf <= 0:
            continue
        if g == "SU":
            _bump(i, "fund", Nf)
            _bump(i, "antifund", Nf)
        elif g == "SO":
            _bump(i, "V", Nf)
        elif g == "Sp":
            _bump(i, "fund", 2 * Nf)
    return fields


def r_values_max_Nf(quiver: Quiver, N: int) -> tuple[dict[str, float], list[int]]:
    """Run a-max with N_f saturated to b_0 = 0 at each node independently.
    Returns (R_charges by label, Nf_per_node)."""
    Nf_per_node = nf_max_per_node(quiver, N)
    fields = build_fields_max_Nf(quiver, N, Nf_per_node)
    A, b = anomaly_matrix(fields, quiver, N)
    R0, *_ = np.linalg.lstsq(A, b, rcond=None)
    F = null_space(A)
    n_free = F.shape[1]
    if n_free == 0:
        s_opt = np.zeros(0)
    else:
        def neg_a(s: np.ndarray) -> float:
            R = R0 + F @ s
            return -a_trial(R, fields, quiver, N)
        s_opt = _minimize(neg_a, x0=np.zeros(n_free), method="BFGS").x
    R_opt = R0 + F @ s_opt
    return {f.label: float(R_opt[f.R_index]) for f in fields}, Nf_per_node


def _all_labels_max_Nf(quiver: Quiver) -> set[str]:
    """Union of build_fields_max_Nf labels across N probes; captures
    fund/antifund identity changes as chiral excess flips sign with N."""
    out: set[str] = set()
    for N in _LABEL_PROBE_N:
        try:
            Nf_per_node = nf_max_per_node(quiver, N)
            for f in build_fields_max_Nf(quiver, N, Nf_per_node):
                out.add(f.label)
        except Exception:
            pass
    return out


def enumerate_candidates_max_Nf(quiver: Quiver, max_degree: int = 6) -> list[CandidateOp]:
    """Same enumeration as `enumerate_candidates`, but the available label
    set is taken from build_fields_max_Nf so that fund/antifund/V/f added by
    saturating b_0 = 0 are eligible to appear in operators."""
    return _enumerate_with_labels(quiver, max_degree, _all_labels_max_Nf(quiver))


def find_marginal_ops_max_Nf(
    quiver: Quiver,
    N_list: tuple[int, ...] = (10, 20, 30),
    max_degree: int = 6,
    tol: float = 1e-6,
) -> tuple[list[CandidateOp], dict[int, list[int]]]:
    """Like find_marginal_ops, but with N_f saturating b_0 = 0 per node.
    Returns (marginal ops, per-N Nf_per_node mapping)."""
    R_per_N: dict[int, dict[str, float]] = {}
    Nf_per_N: dict[int, list[int]] = {}
    for N in N_list:
        R, Nf_pn = r_values_max_Nf(quiver, N)
        if not _r_values_sane(R):
            return [], {}
        R_per_N[N] = R
        Nf_per_N[N] = Nf_pn
    candidates = enumerate_candidates_max_Nf(quiver, max_degree=max_degree)
    return [op for op in candidates if is_marginal_at_all_N(op, R_per_N, tol=tol)], Nf_per_N


# ── Compact pretty-printing ───────────────────────────────────────────────────

_REP_UNI = {
    "adj": "adj", "S": "S", "Sbar": "S̄", "A": "A", "Abar": "Ā",
    "V": "V", "fund": "□", "antifund": "□̄",
}
_EDGE_REP_UNI = {
    "pm": "(+−)", "pp": "(++)", "mm": "(−−)",
    "p": "(+)", "m": "(−)", "std": "",
}
_SUPER = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
_SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def _label_short(label: str) -> str:
    if label.startswith("node"):
        node, rep = label[4:].split("_", 1)
        return _REP_UNI.get(rep, rep) + node.translate(_SUB)
    parts = label.split("_")
    src, dst, rep = parts[1], parts[2], "_".join(parts[3:])
    return f"Q{_EDGE_REP_UNI.get(rep, rep)}{src.translate(_SUB)}{dst.translate(_SUB)}"


def op_short(op: CandidateOp) -> str:
    """Compact one-line operator notation, e.g. 'tr(adj₁³)' or 'tr(Q(++)₀₁ Q(−−)₀₁ adj₁²)'."""
    if op.kind == "bifund-loop" and op.word is not None:
        parts = [_label_short(lbl) for lbl in op.word]
        return "tr(" + " ".join(parts) + ")"
    parts = []
    for label, mult in op.factors:
        s = _label_short(label)
        if mult > 1:
            s += str(mult).translate(_SUPER)
        parts.append(s)
    return "tr(" + " ".join(parts) + ")"


def ops_short_summary(ops: list[CandidateOp], max_chars: int = 60) -> str:
    """Render a list of ops as a semicolon-separated string, truncated."""
    if not ops:
        return "—"
    strs = [op_short(op) for op in ops]
    out = "; ".join(strs)
    if len(out) > max_chars:
        # truncate keeping count visible
        kept = []
        used = 0
        for s in strs:
            if used + len(s) + 2 > max_chars - 6:
                break
            kept.append(s)
            used += len(s) + 2
        out = "; ".join(kept) + f"; …+{len(strs) - len(kept)}"
    return out


# ── DB row → Quiver reconstruction ────────────────────────────────────────────

def _nfc(s: str) -> str:
    return unicodedata.normalize("NFC", s)


_REP_UNICODE_MAP = {_nfc(k): v for k, v in {
    "adj": "adj", "S": "S", "Sbar": "Sbar", "S̄": "Sbar",
    "A": "A", "Abar": "Abar", "Ā": "Abar",
    "V": "V", "fund": "fund", "antifund": "antifund",
    "□": "fund", "□̄": "antifund",
}.items()}


def parse_matter(matter_str: str) -> dict[str, int]:
    """Parse the Unicode matter string (e.g. 'S + 2Ā', 'adj + S̄') back to a
    rep→count dict. Inverse of a_maximization_large_N._fmt_matter."""
    matter_str = _nfc((matter_str or "").strip())
    if not matter_str or matter_str == "—":
        return {}
    out: dict[str, int] = {}
    for term in matter_str.split("+"):
        term = term.strip()
        if not term:
            continue
        m = re.match(r"^(\d*)(.+)$", term)
        n = int(m.group(1)) if m.group(1) else 1
        sym = _nfc(m.group(2).strip())
        rep = _REP_UNICODE_MAP.get(sym)
        if rep is None:
            raise ValueError(f"Unknown matter symbol: {sym!r} in {matter_str!r}")
        out[rep] = out.get(rep, 0) + n
    return out


def parse_edges(edges_str: str) -> list[Edge]:
    """Parse the edges string (e.g. '2×(0→1,+-)  (0→1,++)') back to Edges.
    Inverse of a_maximization_large_N._fmt_edges."""
    edges_str = (edges_str or "").strip()
    if not edges_str or edges_str == "—":
        return []
    edges: list[Edge] = []
    pattern = re.compile(r"(?:(\d+)×)?\((\d+)→(\d+),([^)]+)\)")
    for m in pattern.finditer(edges_str):
        mult = int(m.group(1)) if m.group(1) else 1
        src, dst, rep = int(m.group(2)), int(m.group(3)), m.group(4)
        for _ in range(mult):
            edges.append(Edge(src=src, dst=dst, rep=rep))
    return edges


def quiver_from_row(row: dict) -> Quiver:
    """Reconstruct a Quiver from a theory-table row (dict-like)."""
    return Quiver(
        gauge_types=[row["gauge0"], row["gauge1"]],
        edges=parse_edges(row["edges"]),
        node_matter=[parse_matter(row["matter0"]), parse_matter(row["matter1"])],
        rank_multipliers=[row["rank0_mult"], row["rank1_mult"]],
    )
