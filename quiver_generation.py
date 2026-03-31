"""
quiver_generation.py — Module 2

Quiver data structure, anomaly constraints, and enumeration of valid quiver
gauge theories.

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
    NodeSpec,
    T_bifund,
    compute_b0,
    is_af_all_N,
)

# ── Constants (§2, §3) ─────────────────────────────────────────────────────────

# Valid edge representation types per gauge group pair (§3 bifundamental table).
# SU-SU admits two distinct types; all other pairs have exactly one.
EDGE_REPS: dict[frozenset[str], list[str]] = {
    frozenset(["SU", "SU"]): ["+-", "++"],  # "+-": (□,□̄) directed; "++": (□,□) undirected
    frozenset(["SU", "SO"]): ["std"],        # (□, V)
    frozenset(["SU", "Sp"]): ["std"],        # (□, f)
    frozenset(["SO", "SO"]): ["std"],        # (V, V)
    frozenset(["SO", "Sp"]): ["std"],        # (V, f)
    frozenset(["Sp", "Sp"]): ["std"],        # (f, f)
}

# Valid node-level matter representations per gauge group type (§3 node matter table).
NODE_REPS: dict[str, list[str]] = {
    "SU": ["fund", "antifund", "adj", "S", "Sbar", "A", "Abar"],
    "SO": ["V", "adj", "S"],
    "Sp": ["fund", "adj", "A"],
}

# Fund-like reps per gauge group type (those enumerated up to max_fund).
FUND_REPS: dict[str, list[str]] = {
    "SU": ["fund", "antifund"],
    "SO": ["V"],
    "Sp": ["fund"],
}


# ── Data structures (§1) ───────────────────────────────────────────────────────

@dataclass(frozen=True)
class Edge:
    """A single bifundamental edge in the quiver."""
    src: int   # source node index
    dst: int   # destination node index (src != dst; no self-loops)
    rep: str   # "+-" (□,□̄) | "++" (□,□) | "std" (unique type for non-SU-SU pairs)


@dataclass
class Quiver:
    """
    A quiver gauge theory Q = (V, E, gauge_types, edges, node_matter).

    Attributes
    ----------
    gauge_types : list of GaugeType, one per node
    edges : list of Edge (directed bifundamental edges, with multiplicity; no self-loops)
    node_matter : list of dicts mapping rep name -> count, one per node
    """
    gauge_types: list[GaugeType]
    edges: list[Edge] = field(default_factory=list)
    node_matter: list[dict[str, int]] = field(default_factory=list)

    def __post_init__(self) -> None:
        n = len(self.gauge_types)
        if len(self.node_matter) == 0:
            self.node_matter = [{} for _ in range(n)]
        assert len(self.node_matter) == n

    @property
    def n_nodes(self) -> int:
        return len(self.gauge_types)

    def degree(self, node: int) -> int:
        """Total number of bifundamental edges incident to node (in + out, with multiplicity)."""
        return sum(1 for e in self.edges if e.src == node or e.dst == node)

    def incident_edges(self, node: int) -> list[tuple[Edge, str]]:
        """Return (edge, side) pairs where side is 'src' or 'dst'."""
        return [
            (e, "src" if e.src == node else "dst")
            for e in self.edges
            if e.src == node or e.dst == node
        ]

    def bifund_neighbor_types(self, node: int) -> list[GaugeType]:
        """
        List of gauge types of bifundamental neighbors of node (with multiplicity).
        Suitable for NodeSpec.bifund_neighbors.
        """
        result = []
        for e in self.edges:
            if e.src == node:
                result.append(self.gauge_types[e.dst])
            elif e.dst == node:
                result.append(self.gauge_types[e.src])
        return result

    def node_spec(self, node: int, N: int) -> NodeSpec:
        """Build a NodeSpec for node at rank N."""
        return NodeSpec(
            gauge_type=self.gauge_types[node],
            N=N,
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


# ── Anomaly checks (§4) ────────────────────────────────────────────────────────

def su_anomaly_polynomial(quiver: Quiver, node: int) -> tuple[int, int]:
    """
    Return (coeff_N, coeff_const) for the SU(N) cubic gauge anomaly at node,
    expressed as a polynomial: anomaly(N) = coeff_N * N + coeff_const.

    The anomaly-free condition for all N requires both coefficients to vanish.

    Contributions (§4a):
    - Node matter:     (n_S - n_Sbar) + (n_A - n_Abar) to coeff_N
                       (n_f - n_fbar) + 4*(n_S-n_Sbar) - 4*(n_A-n_Abar) to coeff_const
    - SU-SU "+-" (node→b): +N;  (b→node): -N
    - SU-SU "++" edge:  +N  (to both endpoints)
    - SU-SO edge:  +N   [dim(V_SO) = N]
    - SU-Sp edge:  +2N  [dim(f_Sp) = 2N]
    """
    assert quiver.gauge_types[node] == "SU"
    m = quiver.node_matter[node]

    n_f   = m.get("fund", 0)
    n_fb  = m.get("antifund", 0)
    n_S   = m.get("S", 0)
    n_Sb  = m.get("Sbar", 0)
    n_A   = m.get("A", 0)
    n_Ab  = m.get("Abar", 0)

    coeff_N     = (n_S - n_Sb) + (n_A - n_Ab)
    coeff_const = (n_f - n_fb) + 4 * (n_S - n_Sb) - 4 * (n_A - n_Ab)

    for e, side in quiver.incident_edges(node):
        neighbor   = e.dst if side == "src" else e.src
        g_neighbor = quiver.gauge_types[neighbor]

        if g_neighbor == "SU":
            if e.rep == "+-":
                coeff_N += +1 if side == "src" else -1
            elif e.rep == "++":
                coeff_N += +1
        elif g_neighbor == "SO":
            coeff_N += +1   # dim(V_SO) = N
        elif g_neighbor == "Sp":
            coeff_N += +2   # dim(f_Sp) = 2N

    return coeff_N, coeff_const


def check_su_anomaly(quiver: Quiver, node: int) -> bool:
    """Return True if the SU(N) gauge anomaly vanishes at node for all N."""
    coeff_N, coeff_const = su_anomaly_polynomial(quiver, node)
    return coeff_N == 0 and coeff_const == 0


def sp_fund_count(quiver: Quiver, node: int) -> int:
    """
    Total number of fundamental chiral multiplets at an Sp(N) node:
    node-level funds + one per bifundamental edge (§4b).
    """
    assert quiver.gauge_types[node] == "Sp"
    return quiver.node_matter[node].get("fund", 0) + quiver.degree(node)


def check_sp_witten(quiver: Quiver, node: int) -> bool:
    """Return True if the Sp(N) Witten anomaly vanishes (fund count even, §4b)."""
    return sp_fund_count(quiver, node) % 2 == 0


def check_anomalies(quiver: Quiver) -> bool:
    """Check all gauge and global anomaly conditions for every node."""
    for i, g in enumerate(quiver.gauge_types):
        if g == "SU" and not check_su_anomaly(quiver, i):
            return False
        if g == "Sp" and not check_sp_witten(quiver, i):
            return False
    return True


# ── Asymptotic freedom (§5) ────────────────────────────────────────────────────

def check_af(quiver: Quiver) -> bool:
    """Return True if b_0 > 0 for all valid N at every node."""
    for i, g in enumerate(quiver.gauge_types):
        if not is_af_all_N(g, quiver.node_matter[i], quiver.bifund_neighbor_types(i)):
            return False
    return True


# ── Enumeration helpers ────────────────────────────────────────────────────────

def _matter_bounds(gauge_type: GaugeType) -> dict[str, int]:
    """
    Upper bounds on rank-2 matter multiplicities, derived from requiring the
    leading-N coefficient of b_0 to be non-negative (necessary for AF for all N).
    Fund-like reps are bounded separately by max_fund.
    """
    if gauge_type == "SU":
        # leading coeff of b_0: 3 - n_adj - (n_S+n_Sbar)/2 - (n_A+n_Abar)/2 >= 0
        return {"fund": 0, "antifund": 0, "adj": 2, "S": 5, "Sbar": 5, "A": 5, "Abar": 5}
    if gauge_type == "SO":
        # leading coeff: 3 - n_adj - n_S >= 0
        return {"V": 0, "adj": 2, "S": 2}
    if gauge_type == "Sp":
        # leading coeff: 3 - n_adj - n_A >= 0
        return {"fund": 0, "adj": 2, "A": 2}
    raise ValueError(gauge_type)


def _enumerate_node_matter(
    gauge_type: GaugeType,
    bifund_neighbors: list[GaugeType],
    max_fund: int = 6,
) -> Iterator[dict[str, int]]:
    """
    Yield all matter dicts for a node that satisfy b_0 > 0 for all valid N.

    Reps are enumerated in the order given by NODE_REPS. For each rep:
    - count = 0: always allowed (cannot worsen AF); recurse
    - count > 0: recurse only if AF still holds; break on first failure
      (valid since b_0 is monotonically decreasing in each rep count)
    """
    reps       = NODE_REPS[gauge_type]
    bounds     = _matter_bounds(gauge_type)
    fund_reps  = FUND_REPS[gauge_type]

    def _gen(idx: int, current: dict[str, int]) -> Iterator[dict[str, int]]:
        if idx == len(reps):
            if is_af_all_N(gauge_type, current, bifund_neighbors):
                yield dict(current)
            return
        rep = reps[idx]
        hi  = max_fund if rep in fund_reps else bounds[rep]
        for count in range(hi + 1):
            if count > 0:
                current[rep] = count
                if not is_af_all_N(gauge_type, current, bifund_neighbors):
                    del current[rep]
                    break  # b_0 only decreases with more matter; no need to try higher
            yield from _gen(idx + 1, current)
        if rep in current:
            del current[rep]

    yield from _gen(0, {})


def _enumerate_edges(
    gauge_types: list[GaugeType],
    max_multiedge: int = 2,
) -> Iterator[list[Edge]]:
    """
    Yield all bifundamental edge multisets for a fixed gauge type assignment.

    Iterates over ordered pairs (i, j) with i < j (no self-loops).
    For each pair:
    - SU-SU: three independent counts — "+-" i→j, "+-" j→i, "++"
    - Other pairs: one count for "std" edges (unordered)
    Each count ranges from 0 to max_multiedge.
    """
    n = len(gauge_types)
    # Build list of (i, j, slot_key) for all valid slots
    slots: list[tuple[int, int, str]] = []
    for i in range(n):
        for j in range(i + 1, n):   # strict i < j: no self-loops
            gi, gj = gauge_types[i], gauge_types[j]
            key = frozenset([gi, gj])
            if key not in EDGE_REPS:
                continue
            for rep in EDGE_REPS[key]:
                if rep == "+-":
                    slots.append((i, j, "+-_fwd"))   # i→j
                    slots.append((i, j, "+-_bwd"))   # j→i
                else:
                    slots.append((i, j, rep))

    for counts in product(range(max_multiedge + 1), repeat=len(slots)):
        edges: list[Edge] = []
        for (i, j, slot_key), count in zip(slots, counts):
            for _ in range(count):
                if slot_key == "+-_fwd":
                    edges.append(Edge(i, j, "+-"))
                elif slot_key == "+-_bwd":
                    edges.append(Edge(j, i, "+-"))
                else:
                    edges.append(Edge(i, j, slot_key))
        yield edges


# ── Main enumeration (§6) ──────────────────────────────────────────────────────

def enumerate_quivers(
    max_nodes: int = 2,
    max_fund: int = 6,
    max_multiedge: int = 2,
    require_connected: bool = True,
) -> list[Quiver]:
    """
    Enumerate all quiver gauge theories satisfying (§6 steps 1–6):
    - UV asymptotic freedom: b_0 > 0 for all valid N at every node
    - SU(N) gauge anomaly cancellation
    - Sp(N) Witten anomaly cancellation
    - Connected (if require_connected=True)

    Does not deduplicate graph-isomorphic quivers.

    Parameters
    ----------
    max_nodes     : maximum number of gauge nodes
    max_fund      : maximum multiplicity for fund-like representations per node
    max_multiedge : maximum multiplicity for each edge slot between a pair of nodes
    """
    gauge_group_types: list[GaugeType] = ["SU", "SO", "Sp"]
    valid: list[Quiver] = []

    for n_nodes in range(1, max_nodes + 1):
        # Step 3: enumerate gauge type assignments
        for gauge_types in product(gauge_group_types, repeat=n_nodes):
            gauge_types = list(gauge_types)

            # Steps 1–2 + 4: enumerate edge configurations
            for edges in _enumerate_edges(gauge_types, max_multiedge):
                q_partial = Quiver(gauge_types, edges)

                # Connectivity filter
                if require_connected and not q_partial.is_connected():
                    continue

                # Step 5: enumerate node matter (AF-filtered per node)
                matter_per_node: list[list[dict[str, int]]] = []
                for i in range(n_nodes):
                    neighbors = q_partial.bifund_neighbor_types(i)
                    choices   = list(_enumerate_node_matter(gauge_types[i], neighbors, max_fund))
                    matter_per_node.append(choices)

                # Step 6: apply all constraints
                for matter_combo in product(*matter_per_node):
                    q = Quiver(gauge_types, edges, list(matter_combo))
                    if check_af(q) and check_anomalies(q):
                        valid.append(q)

    return valid


# ── Utilities ──────────────────────────────────────────────────────────────────

def quiver_summary(q: Quiver) -> str:
    """Human-readable summary of a quiver."""
    lines = [f"Nodes: {' '.join(f'{g}({i})' for i, g in enumerate(q.gauge_types))}"]
    for e in q.edges:
        lines.append(f"  Edge {e.src}→{e.dst}  rep={e.rep}")
    for i, matter in enumerate(q.node_matter):
        non_zero = {k: v for k, v in matter.items() if v > 0}
        if non_zero:
            lines.append(f"  Node {i} matter: {non_zero}")
    return "\n".join(lines)


# ── Validation ─────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Single SU(N) node: SQCD with N_f fund + N_f antifund
    # b_0 = 3N - N_f > 0 for all N>=2  =>  N_f <= 5
    # Anomaly: n_f - n_fbar = 0  =>  ok
    print("=== Single SU(N) node, SQCD-like ===")
    for n_f in range(7):
        q   = Quiver(["SU"], [], [{"fund": n_f, "antifund": n_f}])
        print(f"  N_f={n_f}: AF={check_af(q)}  anomaly_free={check_anomalies(q)}")

    print()

    # Two SU(N) nodes, one "+-" bifundamental (unbalanced anomaly)
    print("=== Two SU(N) nodes, one +- edge, no extra matter ===")
    q = Quiver(["SU", "SU"], [Edge(0, 1, "+-")])
    print(f"  AF={check_af(q)}  anomaly_free={check_anomalies(q)}")
    coeff = su_anomaly_polynomial(q, 0)
    print(f"  SU anomaly at node 0: {coeff[0]}*N + {coeff[1]}")

    print()

    # Two SU(N) nodes, one +- each direction: balanced
    # b_0 = 3N - N = 2N > 0;  anomaly: +N - N = 0
    print("=== Two SU(N) nodes, one +- each direction (balanced) ===")
    q = Quiver(["SU", "SU"], [Edge(0, 1, "+-"), Edge(1, 0, "+-")])
    print(f"  AF={check_af(q)}  anomaly_free={check_anomalies(q)}")

    print()

    # Sp(N) Witten anomaly
    print("=== Sp(N) node, 2 fundamentals ===")
    q = Quiver(["Sp"], [], [{"fund": 2}])
    print(f"  AF={check_af(q)}  Witten OK={check_sp_witten(q, 0)}")

    print()

    print("=== Sp(N) node, 3 fundamentals ===")
    q = Quiver(["Sp"], [], [{"fund": 3}])
    print(f"  AF={check_af(q)}  Witten OK={check_sp_witten(q, 0)}")

    print()

    # 1-node enumeration
    print("=== Enumerate 1-node quivers ===")
    results = enumerate_quivers(max_nodes=1, max_fund=6, max_multiedge=0)
    print(f"  Found {len(results)} valid 1-node quivers")

    # Verify every result passes both checks
    for q in results:
        assert check_af(q),        f"AF failed: {quiver_summary(q)}"
        assert check_anomalies(q), f"Anomaly failed: {quiver_summary(q)}"
    print("  All pass check_af and check_anomalies.")
