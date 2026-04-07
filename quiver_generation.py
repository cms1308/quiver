"""
quiver_generation.py — Module 2

Quiver data structure, N_f-parametrized AF bounds, anomaly constraints,
and enumeration of valid quiver gauge theories.

N_f is defined as:
  SU(N): min(n_f, n_fbar) — number of complete (fund, antifund) pairs
  SO(N): number of vectors
  Sp(N): half the number of fundamentals

For each theory (gauge types, edges, rank-2/adj matter), the one-loop AF
condition is N_f < alpha * N + gamma, where (alpha, gamma) = nf_bound(q, node).

Reference: module2_quiver_generation.md
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from itertools import product
from typing import Iterator

from beta_functions import (
    GaugeType,
    N_MIN,
    RANK2_ADJ_REPS,
    NodeSpec,
    compute_b0,
    b0_linear,
    is_af_all_N,
)

# ── Constants ──────────────────────────────────────────────────────────────────

# Valid edge representation types per gauge group pair.
# Sign refers to the SU-side chirality: "+" = □ (fund), "-" = □̄ (antifund).
# SO/Sp reps (V, f) are real/pseudo-real, so the non-SU side is always self-conjugate.
EDGE_REPS: dict[frozenset[str], list[str]] = {
    frozenset(["SU", "SU"]): ["+-", "++", "--"],  # (□,□̄) directed; (□,□) and (□̄,□̄) undirected
    frozenset(["SU", "SO"]): ["+", "-"],           # (□,V) and (□̄,V)
    frozenset(["SU", "Sp"]): ["+", "-"],           # (□,f) and (□̄,f)
    frozenset(["SO", "SO"]): ["std"],              # (V,V) self-conjugate
    frozenset(["SO", "Sp"]): ["std"],              # (V,f) self-conjugate
    frozenset(["Sp", "Sp"]): ["std"],              # (f,f) self-conjugate
}

# Rank-2/adj reps only — fund-like reps are absorbed into N_f.
NODE_REPS: dict[str, list[str]] = RANK2_ADJ_REPS


# ── Data structures ────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class Edge:
    """A single bifundamental edge in the quiver."""
    src: int   # source node index
    dst: int   # destination node index
    rep: str   # SU-SU: "+-" (□,□̄), "++" (□,□), "--" (□̄,□̄)
               # SU-SO/Sp: "+" (□,V/f), "-" (□̄,V/f)
               # SO/Sp-SO/Sp: "std" (self-conjugate)


@dataclass
class Quiver:
    """
    A quiver gauge theory Q = (V, E, gauge_types, edges, node_matter).

    node_matter contains only rank-2/adj representations; fund-like matter
    (fundamentals, antifundamentals, vectors) is parametrized by N_f.

    rank_multipliers[i] = m_i means node i has gauge group G(m_i * N).
    Defaults to [1, 1, ...] (all nodes share the same rank N).
    """
    gauge_types: list[GaugeType]
    edges: list[Edge] = field(default_factory=list)
    node_matter: list[dict[str, int]] = field(default_factory=list)
    rank_multipliers: list[int] = field(default_factory=list)

    def __post_init__(self) -> None:
        n = len(self.gauge_types)
        if len(self.node_matter) == 0:
            self.node_matter = [{} for _ in range(n)]
        assert len(self.node_matter) == n
        if len(self.rank_multipliers) == 0:
            self.rank_multipliers = [1] * n
        assert len(self.rank_multipliers) == n

    @property
    def n_nodes(self) -> int:
        return len(self.gauge_types)

    def degree(self, node: int) -> int:
        """Total number of bifundamental edges incident to node (in + out)."""
        return sum(1 for e in self.edges if e.src == node or e.dst == node)

    def incident_edges(self, node: int) -> list[tuple[Edge, str]]:
        """Return (edge, side) pairs where side is 'src' or 'dst'."""
        return [
            (e, "src" if e.src == node else "dst")
            for e in self.edges
            if e.src == node or e.dst == node
        ]

    def bifund_neighbor_types(self, node: int) -> list[GaugeType]:
        """List of gauge types of bifundamental neighbors (with multiplicity)."""
        result = []
        for e in self.edges:
            if e.src == node:
                result.append(self.gauge_types[e.dst])
            elif e.dst == node:
                result.append(self.gauge_types[e.src])
        return result

    def bifund_neighbor_mults(self, node: int) -> list[int]:
        """List of rank multipliers of bifundamental neighbors (parallel to bifund_neighbor_types)."""
        result = []
        for e in self.edges:
            if e.src == node:
                result.append(self.rank_multipliers[e.dst])
            elif e.dst == node:
                result.append(self.rank_multipliers[e.src])
        return result

    def node_spec(self, node: int, N: int) -> NodeSpec:
        """Build a NodeSpec for node at base rank N (uses rank_multipliers[node] * N)."""
        m = self.rank_multipliers[node]
        return NodeSpec(
            gauge_type=self.gauge_types[node],
            N=m * N,
            matter=self.node_matter[node],
            bifund_neighbors=self.bifund_neighbor_types(node),
        )

    def is_connected(self) -> bool:
        """Return True if the underlying undirected graph is connected."""
        if self.n_nodes == 1:
            return True
        adj: list[set[int]] = [set() for _ in range(self.n_nodes)]
        for e in self.edges:
            adj[e.src].add(e.dst)
            adj[e.dst].add(e.src)
        visited = {0}
        queue = [0]
        while queue:
            v = queue.pop()
            for u in adj[v]:
                if u not in visited:
                    visited.add(u)
                    queue.append(u)
        return len(visited) == self.n_nodes


# ── N_f bound computation ──────────────────────────────────────────────────────

def chiral_excess_coeffs(quiver: Quiver, node: int) -> tuple[int, int]:
    """
    For an SU(m_a*N) node, return (a, b) such that the chiral excess
    delta(N) = a*N + b, where delta = n_f - n_fbar is fixed by
    gauge anomaly cancellation:

        delta + (n_S - n_Sbar)(m_a*N+4) + (n_A - n_Abar)(m_a*N-4) + bif_anomaly = 0

    Rank-2 anomaly contributions at SU(m_a*N):
      (n_S - n_Sb)*(m_a*N+4) + (n_A - n_Ab)*(m_a*N-4)
      → N-coeff: [(n_S-n_Sb) + (n_A-n_Ab)] * m_a
      → const:   4*(n_S-n_Sb) - 4*(n_A-n_Ab)

    Bifundamental anomaly from neighbor of type g_nb with mult m_b
    (A_fund × dim_fund(g_nb, m_b*N)):
      SU-SU "+-" src: +m_b*N;  dst: -m_b*N
      SU-SU "++": +m_b*N;  "--": -m_b*N
      SU-SO "+": +m_b*N;  "-": -m_b*N  (dim_V = m_b*N)
      SU-Sp "+": +2*m_b*N;  "-": -2*m_b*N  (dim_f = 2*m_b*N)
    """
    assert quiver.gauge_types[node] == "SU"
    m_a = quiver.rank_multipliers[node]
    matter = quiver.node_matter[node]
    n_S  = matter.get("S",    0)
    n_Sb = matter.get("Sbar", 0)
    n_A  = matter.get("A",    0)
    n_Ab = matter.get("Abar", 0)

    # Rank-2 anomaly: scaled by m_a
    r2_coeff_N     = ((n_S - n_Sb) + (n_A - n_Ab)) * m_a
    r2_coeff_const = 4 * (n_S - n_Sb) - 4 * (n_A - n_Ab)

    # Bifundamental anomaly: scaled by neighbor's m_b
    bif_coeff_N = 0
    for e, side in quiver.incident_edges(node):
        neighbor = e.dst if side == "src" else e.src
        g_nb = quiver.gauge_types[neighbor]
        m_b  = quiver.rank_multipliers[neighbor]
        if g_nb == "SU":
            if e.rep == "+-":
                bif_coeff_N += m_b * (+1 if side == "src" else -1)
            elif e.rep == "++":
                bif_coeff_N += m_b
            elif e.rep == "--":
                bif_coeff_N -= m_b
        elif g_nb == "SO":
            if e.rep == "+":
                bif_coeff_N += m_b
            elif e.rep == "-":
                bif_coeff_N -= m_b
        elif g_nb == "Sp":
            if e.rep == "+":
                bif_coeff_N += 2 * m_b
            elif e.rep == "-":
                bif_coeff_N -= 2 * m_b

    # delta = -(rank-2 anomaly + bif anomaly)
    return -(r2_coeff_N + bif_coeff_N), -r2_coeff_const


def _nf_bound_at_N(quiver: Quiver, node: int, N: int) -> Fraction:
    """
    Return f(N) such that the node is AF iff N_f < f(N).

    f(N) = b_0(rank-2/adj matter, N_f=0) - |delta(N)| / 2

    where delta(N) is the SU chiral excess (zero for SO/Sp nodes).
    """
    b0_no_fund = compute_b0(quiver.node_spec(node, N))

    if quiver.gauge_types[node] != "SU":
        return b0_no_fund

    a, b = chiral_excess_coeffs(quiver, node)
    delta = a * N + b
    return b0_no_fund - Fraction(abs(delta), 2)


def nf_bound(quiver: Quiver, node: int) -> tuple[Fraction, Fraction] | None:
    """
    Return (alpha, gamma) such that the one-loop AF condition is N_f < alpha*N + gamma
    for all large N.

    For SU nodes, f(N) = b0_rank2_adj(N) - |delta(N)|/2 is piecewise linear (delta can
    change sign). We compute the large-N asymptotic linear form directly:
        For large N: |delta(N)| ~ |a|*N + sgn(a)*b
        => large-N slope: alpha = alpha_b0 - |a|/2

    Validity requires: large-N alpha >= 0 with positive constant.
    Returns None if the condition cannot be satisfied.
    """
    g = quiver.gauge_types[node]
    m = quiver.rank_multipliers[node]
    neighbors = quiver.bifund_neighbor_types(node)
    n_mults   = quiver.bifund_neighbor_mults(node)

    if g != "SU":
        # No chiral excess: f = b0_rank2_adj is exactly linear
        alpha, gamma = b0_linear(g, quiver.node_matter[node], neighbors,
                                 rank_mult=m, neighbor_mults=n_mults)
        if alpha < 0 or (alpha == 0 and gamma <= 0):
            return None
        return alpha, gamma

    # SU: compute large-N asymptotic form of f(N) = b0_rank2_adj(N) - |delta(N)|/2
    a, b = chiral_excess_coeffs(quiver, node)
    alpha_b0, beta_b0 = b0_linear(g, quiver.node_matter[node], neighbors,
                                   rank_mult=m, neighbor_mults=n_mults)

    alpha = alpha_b0 - Fraction(abs(a), 2)
    if a == 0:
        gamma = beta_b0 - Fraction(abs(b), 2)
    elif a > 0:
        gamma = beta_b0 - Fraction(b, 2)
    else:  # a < 0: for large N, |delta| = -a*N - b
        gamma = beta_b0 + Fraction(b, 2)

    if alpha < 0 or (alpha == 0 and gamma <= 0):
        return None
    return alpha, gamma


# ── Anomaly checks ─────────────────────────────────────────────────────────────

def check_sp_witten(quiver: Quiver, node: int) -> bool:
    """
    Check the Sp(N) Witten anomaly in the N_f parametrization.

    N_f pairs contribute 2*N_f fundamentals (always even). Each bifundamental
    edge contributes one more fundamental. Witten condition: degree must be even.
    """
    assert quiver.gauge_types[node] == "Sp"
    return quiver.degree(node) % 2 == 0


def check_anomalies(quiver: Quiver) -> bool:
    """Check Sp(N) Witten anomaly at all Sp nodes."""
    for i, g in enumerate(quiver.gauge_types):
        if g == "Sp" and not check_sp_witten(quiver, i):
            return False
    return True


# ── Enumeration ────────────────────────────────────────────────────────────────

def _enumerate_node_matter(
    gauge_type: GaugeType,
    bifund_neighbors: list[GaugeType],
    rank_mult: int = 1,
    neighbor_mults: list[int] | None = None,
) -> Iterator[dict[str, int]]:
    """
    Yield all rank-2/adj matter dicts satisfying b_0(rank-2/adj, N_f=0) > 0
    for all valid N, using backtracking pruning.

    Pruning: b_0 decreases monotonically with each rep count, so if count=k
    violates AF, all count > k also do — break immediately.
    """
    reps = NODE_REPS[gauge_type]
    if neighbor_mults is None:
        neighbor_mults = [rank_mult] * len(bifund_neighbors)

    def _gen(idx: int, current: dict[str, int]) -> Iterator[dict[str, int]]:
        if idx == len(reps):
            yield dict(current)
            return
        rep = reps[idx]
        # count = 0: always recurse (omitting a rep cannot worsen AF)
        yield from _gen(idx + 1, current)
        # count > 0: recurse while AF holds, break on first failure
        for count in range(1, 100):
            current[rep] = count
            if not is_af_all_N(gauge_type, current, bifund_neighbors,
                                rank_mult=rank_mult, neighbor_mults=neighbor_mults):
                del current[rep]
                break
            yield from _gen(idx + 1, current)
        if rep in current:
            del current[rep]

    yield from _gen(0, {})


def _enumerate_edges(
    gauge_types: list[GaugeType],
    max_multiedge: int = 2,
    min_multiedge: int = 0,
) -> Iterator[list[Edge]]:
    """
    Generate all edge multisets for a fixed gauge type assignment.
    No self-loops (enforced by i < j strictly).

    min_multiedge: minimum total bifundamental edges per node pair.
    max_multiedge: maximum per individual slot (edge type) per node pair.
    """
    n = len(gauge_types)

    # Group slots by node pair (i, j)
    pair_slots: dict[tuple[int, int], list[str]] = {}
    for i in range(n):
        for j in range(i + 1, n):
            gi, gj = gauge_types[i], gauge_types[j]
            key = frozenset([gi, gj])
            if key not in EDGE_REPS:
                continue
            slots: list[str] = []
            for rep in EDGE_REPS[key]:
                if rep == "+-":
                    slots.append("+-_fwd")
                    slots.append("+-_bwd")
                else:
                    slots.append(rep)
            pair_slots[(i, j)] = slots

    if not pair_slots:
        yield []
        return

    # For each pair, enumerate valid count-tuples (total >= min_multiedge, each <= max_multiedge)
    pairs = list(pair_slots.keys())
    pair_choices: list[list[tuple[int, ...]]] = []
    for pair in pairs:
        n_slots = len(pair_slots[pair])
        choices = [
            combo for combo in product(range(max_multiedge + 1), repeat=n_slots)
            if sum(combo) >= min_multiedge
        ]
        pair_choices.append(choices)

    def _make_edges(pair: tuple[int, int], rep_key: str, count: int) -> list[Edge]:
        i, j = pair
        result = []
        for _ in range(count):
            if rep_key == "+-_fwd":
                result.append(Edge(i, j, "+-"))
            elif rep_key == "+-_bwd":
                result.append(Edge(j, i, "+-"))
            else:
                result.append(Edge(i, j, rep_key))
        return result

    for combo in product(*pair_choices):
        edges: list[Edge] = []
        for pair, counts in zip(pairs, combo):
            for rep_key, count in zip(pair_slots[pair], counts):
                edges.extend(_make_edges(pair, rep_key, count))
        yield edges


def enumerate_quivers(
    n_nodes: int,
    max_multiedge: int = 2,
    min_multiedge: int = 0,
    require_connected: bool = True,
) -> list[tuple[Quiver, list[tuple[Fraction, Fraction]]]]:
    """
    Enumerate all quivers with exactly n_nodes gauge nodes satisfying:
    - N_f bound f(N) > 0 and non-decreasing at every node
    - Sp(N) Witten anomaly: even degree at every Sp node

    Returns list of (quiver, bounds) where bounds[i] = (alpha, gamma) gives
    the AF condition N_f < alpha*N + gamma at node i.

    Note: does not deduplicate isomorphic quivers.
    """
    gauge_group_types: list[GaugeType] = ["SU", "SO", "Sp"]
    results: list[tuple[Quiver, list[tuple[Fraction, Fraction]]]] = []

    for gauge_types in product(gauge_group_types, repeat=n_nodes):
        gauge_types = list(gauge_types)

        for edges in _enumerate_edges(gauge_types, max_multiedge, min_multiedge):
            q_partial = Quiver(gauge_types, edges)

            if require_connected and not q_partial.is_connected():
                continue

            # Witten anomaly depends only on degree — prune early
            if not check_anomalies(q_partial):
                continue

            # Enumerate rank-2/adj matter for each node
            matter_choices: list[list[dict[str, int]]] = []
            feasible = True
            for i in range(n_nodes):
                neighbors = q_partial.bifund_neighbor_types(i)
                choices = list(_enumerate_node_matter(gauge_types[i], neighbors))
                if not choices:
                    feasible = False
                    break
                matter_choices.append(choices)

            if not feasible:
                continue

            for matter_combo in product(*matter_choices):
                q = Quiver(gauge_types, edges, list(matter_combo))
                # Compute N_f bounds; all nodes must have valid bounds
                bounds: list[tuple[Fraction, Fraction]] = []
                valid = True
                for i in range(n_nodes):
                    b = nf_bound(q, i)
                    if b is None:
                        valid = False
                        break
                    bounds.append(b)
                if valid:
                    results.append((q, bounds))

    return _dedup_conjugation(results)


def enumerate_quivers_mixed_rank(
    n_nodes: int,
    rank_multipliers: list[int],
    max_multiedge: int = 2,
    min_multiedge: int = 0,
    require_connected: bool = True,
) -> list[tuple[Quiver, list[tuple[Fraction, Fraction]]]]:
    """
    Enumerate all quivers with the given fixed rank multipliers satisfying:
    - N_f bound f(N) > 0 and non-decreasing at every node
    - Sp(N) Witten anomaly: even degree at every Sp node

    Returns list of (quiver, bounds) where bounds[i] = (alpha, gamma) gives
    the AF condition N_f < alpha*N + gamma at node i.

    Note: does not deduplicate isomorphic quivers; caller should handle pairs
    like [[2,1],[1,2]] to cover both orderings.
    """
    assert len(rank_multipliers) == n_nodes
    gauge_group_types: list[GaugeType] = ["SU", "SO", "Sp"]
    results: list[tuple[Quiver, list[tuple[Fraction, Fraction]]]] = []

    for gauge_types in product(gauge_group_types, repeat=n_nodes):
        gauge_types = list(gauge_types)

        for edges in _enumerate_edges(gauge_types, max_multiedge, min_multiedge):
            q_partial = Quiver(gauge_types, edges, rank_multipliers=rank_multipliers)

            if require_connected and not q_partial.is_connected():
                continue

            if not check_anomalies(q_partial):
                continue

            matter_choices: list[list[dict[str, int]]] = []
            feasible = True
            for i in range(n_nodes):
                neighbors    = q_partial.bifund_neighbor_types(i)
                n_mults      = q_partial.bifund_neighbor_mults(i)
                choices = list(_enumerate_node_matter(
                    gauge_types[i], neighbors,
                    rank_mult=rank_multipliers[i],
                    neighbor_mults=n_mults,
                ))
                if not choices:
                    feasible = False
                    break
                matter_choices.append(choices)

            if not feasible:
                continue

            for matter_combo in product(*matter_choices):
                q = Quiver(gauge_types, edges, list(matter_combo),
                           rank_multipliers=list(rank_multipliers))
                bounds: list[tuple[Fraction, Fraction]] = []
                valid = True
                for i in range(n_nodes):
                    b = nf_bound(q, i)
                    if b is None:
                        valid = False
                        break
                    bounds.append(b)
                if valid:
                    results.append((q, bounds))

    return results


# ── Conjugation ────────────────────────────────────────────────────────────────

_CONJ_REP: dict[str, str] = {
    "S": "Sbar", "Sbar": "S", "A": "Abar", "Abar": "A", "adj": "adj",
}


def _conjugate_matter(matter: dict[str, int]) -> dict[str, int]:
    """Swap S↔Sbar, A↔Abar (charge conjugation of SU(N) reps)."""
    return {_CONJ_REP[rep]: count for rep, count in matter.items()}


def _quiver_signature(q: Quiver) -> tuple:
    """Hashable canonical signature of a quiver."""
    edges = tuple(sorted((e.src, e.dst, e.rep) for e in q.edges))
    matter = tuple(tuple(sorted(m.items())) for m in q.node_matter)
    return (tuple(q.gauge_types), edges, matter)


def _conjugate_signature(q: Quiver) -> tuple:
    """Signature of the charge-conjugate quiver."""
    conj_edges = []
    for e in q.edges:
        if e.rep == "+-":
            conj_edges.append((e.dst, e.src, "+-"))
        elif e.rep == "++":
            conj_edges.append((e.src, e.dst, "--"))
        elif e.rep == "--":
            conj_edges.append((e.src, e.dst, "++"))
        elif e.rep == "+":
            conj_edges.append((e.src, e.dst, "-"))
        elif e.rep == "-":
            conj_edges.append((e.src, e.dst, "+"))
        else:
            conj_edges.append((e.src, e.dst, e.rep))
    edges = tuple(sorted(conj_edges))
    matter = tuple(tuple(sorted(_conjugate_matter(m).items())) for m in q.node_matter)
    return (tuple(q.gauge_types), edges, matter)


def _dedup_conjugation(
    results: list[tuple[Quiver, list[tuple[Fraction, Fraction]]]],
) -> list[tuple[Quiver, list[tuple[Fraction, Fraction]]]]:
    """Remove duplicate theories related by charge conjugation."""
    seen: set[tuple] = set()
    unique: list[tuple[Quiver, list[tuple[Fraction, Fraction]]]] = []
    for q, bounds in results:
        sig = _quiver_signature(q)
        conj_sig = _conjugate_signature(q)
        canonical = min(sig, conj_sig)
        if canonical not in seen:
            seen.add(canonical)
            unique.append((q, bounds))
    return unique


# ── Utilities ──────────────────────────────────────────────────────────────────

def nf_bound_str(alpha: Fraction, gamma: Fraction) -> str:
    """Format N_f bound as a readable string."""
    if alpha == 0:
        return f"N_f < {gamma}"
    if gamma == 0:
        return f"N_f < {alpha}*N"
    sign = "+" if gamma > 0 else "-"
    return f"N_f < {alpha}*N {sign} {abs(gamma)}"


def delta_str(a: int, b: int) -> str:
    """Format chiral excess delta(N) = a*N + b as a readable string."""
    if a == 0:
        return str(b)
    n_term = ("" if a == 1 else ("-" if a == -1 else f"{a}*")) + "N"
    if b == 0:
        return n_term
    sign = "+" if b > 0 else "-"
    return f"{n_term}{sign}{abs(b)}"


def quiver_summary(
    q: Quiver,
    bounds: list[tuple[Fraction, Fraction]] | None = None,
) -> str:
    """Human-readable summary of a quiver."""
    lines = []
    lines.append(f"Nodes: {' '.join(f'{g}({i})' for i, g in enumerate(q.gauge_types))}")
    for e in q.edges:
        lines.append(f"  Edge {e.src}→{e.dst}  rep={e.rep}")
    for i, (g, matter) in enumerate(zip(q.gauge_types, q.node_matter)):
        non_zero = {k: v for k, v in matter.items() if v > 0}
        if non_zero:
            lines.append(f"  Node {i} matter: {non_zero}")
        if g == "SU":
            a, b = chiral_excess_coeffs(q, i)
            lines.append(f"  Node {i} delta: {delta_str(a, b)}")
    if bounds:
        for i, (alpha, gamma) in enumerate(bounds):
            lines.append(f"  Node {i} AF: {nf_bound_str(alpha, gamma)}")
    return "\n".join(lines)


# ── Quick validation ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    # SU(N), no rank-2 matter: delta=0, f(N) = 3N => N_f < 3N
    print("=== SU(N), no rank-2 matter ===")
    q = Quiver(["SU"])
    b = nf_bound(q, 0)
    print(f"  {nf_bound_str(*b)}")  # expect N_f < 3*N

    print()

    # SU(N) + 1 S: delta = -(N+4), f(N) = 2N - 3
    print("=== SU(N) + 1 S ===")
    q = Quiver(["SU"], node_matter=[{"S": 1}])
    b = nf_bound(q, 0)
    print(f"  {nf_bound_str(*b)}")  # expect N_f < 2*N - 3

    print()

    # SU(N) + 1 A: delta = -(N-4), f(N) = 2N + 3
    print("=== SU(N) + 1 A ===")
    q = Quiver(["SU"], node_matter=[{"A": 1}])
    b = nf_bound(q, 0)
    print(f"  {nf_bound_str(*b)}")  # expect N_f < 2*N + 3

    print()

    # SU(N) + 3 A: constant bound N_f < 9
    print("=== SU(N) + 3 A ===")
    q = Quiver(["SU"], node_matter=[{"A": 3}])
    b = nf_bound(q, 0)
    print(f"  {nf_bound_str(*b) if b else 'not valid'}")  # expect N_f < 9

    print()

    # SU(N) + 4 A: not AF for all N (alpha < 0)
    print("=== SU(N) + 4 A ===")
    q = Quiver(["SU"], node_matter=[{"A": 4}])
    b = nf_bound(q, 0)
    print(f"  {'not valid (alpha < 0)' if b is None else nf_bound_str(*b)}")  # expect not valid

    print()

    # Two SU(N) nodes, balanced +- (one each direction)
    # b_0 per node = 3N - N/2 = 5N/2, delta=0 => N_f < 5N/2
    print("=== Two SU(N) nodes, balanced +- ===")
    q = Quiver(["SU", "SU"], [Edge(0, 1, "+-"), Edge(1, 0, "+-")])
    b0 = nf_bound(q, 0)
    b1 = nf_bound(q, 1)
    print(f"  Node 0: {nf_bound_str(*b0)}")  # expect N_f < 5/2*N
    print(f"  Node 1: {nf_bound_str(*b1)}")

    print()

    # Sp(N) node with odd degree: Witten violation
    print("=== Sp(N), one SU neighbor (degree=1, odd) ===")
    q = Quiver(["SU", "Sp"], [Edge(0, 1, "std")])
    print(f"  Witten OK={check_sp_witten(q, 1)}")  # expect False

    print()

    # Sp(N) node with even degree: Witten OK
    print("=== Sp(N), two SU neighbors (degree=2, even) ===")
    q = Quiver(["SU", "Sp"], [Edge(0, 1, "std"), Edge(0, 1, "std")])
    print(f"  Witten OK={check_sp_witten(q, 1)}")  # expect True

    print()

    # 1-node enumeration with N_f bounds
    print("=== Enumerate 2-node quivers (N_f parametrized) ===")
    results = enumerate_quivers(n_nodes=3, max_multiedge=1, min_multiedge=1, require_connected=True)
    print(f"  Found {len(results)} valid 2-node quivers")

    for q, bounds in results[:14]:
        summary = quiver_summary(q, bounds)
        print(f"  {summary}")