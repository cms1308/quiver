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

from quiver_generation import Edge, Quiver
from a_maximization import a_maximize, build_fields


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


def _node_intra_labels(quiver: Quiver, node: int) -> list[str]:
    """Labels for fields living entirely at `node` (rank-2/adj only)."""
    return [f"node{node}_{rep}" for rep, c in quiver.node_matter[node].items() if c > 0]


def _all_field_labels(quiver: Quiver) -> list[str]:
    """All field labels in the quiver (intra-node + edges)."""
    out: list[str] = []
    for i in range(quiver.n_nodes):
        out.extend(_node_intra_labels(quiver, i))
    for e in quiver.edges:
        out.append(_edge_label(e))
    # dedup edges of identical (src,dst,rep) — they are flavor-degenerate
    return sorted(set(out))


def _enumerate_single_node(quiver: Quiver, max_degree: int) -> list[CandidateOp]:
    out: list[CandidateOp] = []
    multiplets = _flavor_multiplets(quiver)
    for i in range(quiver.n_nodes):
        labels = _node_intra_labels(quiver, i)
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


def _enumerate_mixed(quiver: Quiver, max_degree: int) -> list[CandidateOp]:
    """General enumerator: multisets drawn from intra-node + edge labels with
    gauge invariance at every node. Subsumes single-node and bifund-loop;
    deduplicated against those by (factors, kind=…)."""
    labels = _all_field_labels(quiver)
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


def enumerate_candidates(quiver: Quiver, max_degree: int = 6) -> list[CandidateOp]:
    """Enumerate distinct gauge-invariant single-trace candidate operators
    of degree ≥ 2. Degree-1 trace operators are either identically zero
    (Tr T^a = 0 for all gauge generators) or non-singlet (any rank > 0
    representation needs index contraction)."""
    cands: dict[tuple, CandidateOp] = {}
    for op in _enumerate_single_node(quiver, max_degree):
        if op.degree < 2:
            continue
        cands[(op.kind, op.factors, op.word)] = op
    for op in _enumerate_bifund_loops(quiver, max_degree):
        if op.degree < 2:
            continue
        cands[(op.kind, op.factors, op.word)] = op
    for op in _enumerate_mixed(quiver, max_degree):
        if op.degree < 2:
            continue
        cands[(op.kind, op.factors, op.word)] = op
    return list(cands.values())


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


def find_marginal_singlet_ops(
    quiver: Quiver,
    N_list: tuple[int, ...] = (10, 20, 30),
    max_degree: int = 6,
    tol: float = 1e-6,
    N_f: int = 0,
) -> list[CandidateOp]:
    """End-to-end: enumerate, evaluate R at each N, filter to flavor-singlet
    operators that are marginal at every N. Returns filtered list."""
    R_per_N = {N: r_values_at_N(quiver, N, N_f=N_f) for N in N_list}
    candidates = enumerate_candidates(quiver, max_degree=max_degree)
    multiplets = _flavor_multiplets(quiver)
    return [
        op for op in candidates
        if is_marginal_at_all_N(op, R_per_N, tol=tol)
        and is_flavor_singlet(op, multiplets)
    ]


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
