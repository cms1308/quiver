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

# ── Data structures ────────────────────────────────────────────────────────────

# Valid edge representation types per gauge group pair.
# SU-SU allows two types; all other pairs have exactly one.
EDGE_REPS: dict[frozenset[str], list[str]] = {
    frozenset(["SU", "SU"]): ["+-", "++"],  # (□,□̄) directed; (□,□) undirected
    frozenset(["SU", "SO"]): ["std"],
    frozenset(["SU", "Sp"]): ["std"],
    frozenset(["SO", "SO"]): ["std"],
    frozenset(["SO", "Sp"]): ["std"],
    frozenset(["Sp", "Sp"]): ["std"],
}

# Valid node-level matter representations per gauge group type
NODE_REPS: dict[str, list[str]] = {
    "SU": ["fund", "antifund", "adj", "S", "Sbar", "A", "Abar"],
    "SO": ["V", "adj", "S"],
    "Sp": ["fund", "adj", "A"],
}


@dataclass(frozen=True)
class Edge:
    """A single bifundamental edge in the quiver."""
    src: int   # source node index
    dst: int   # destination node index
    rep: str   # "+-" (□,□̄); "++" (□,□); "std" (unique type for non-SU-SU pairs)


@dataclass
class Quiver:
    """
    A quiver gauge theory Q = (V, E, gauge_types, edges, node_matter).

    Attributes
    ----------
    gauge_types : list of GaugeType, one per node
    edges : list of Edge (directed bifundamental edges, with multiplicity)
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
        """
        List of gauge types of bifundamental neighbors of node (with multiplicity),
        suitable for passing to NodeSpec.bifund_neighbors.
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


# ── Anomaly checks ─────────────────────────────────────────────────────────────

def su_anomaly_polynomial(quiver: Quiver, node: int) -> tuple[int, int]:
    """
    Return (coeff_N, coeff_const) for the SU(N) cubic gauge anomaly at node,
    as a polynomial anomaly(N) = coeff_N * N + coeff_const.

    The anomaly-free condition requires both coefficients to vanish.

    Sources:
    - Node-level matter: A(r) * 1 (no other group factor for node matter)
    - SU-SU "+-" edge (node→b): A(□_node) * dim(□̄_b) = +N
    - SU-SU "+-" edge (b→node): A(□̄_node) * dim(□_b) = -N
    - SU-SU "++" edge: A(□_node) * dim(□_b) = +N
    - SU-SO edge: A(□_node) * dim(V_SO) = +N
    - SU-Sp edge: A(□_node) * dim(f_Sp) = +2N
    """
    assert quiver.gauge_types[node] == "SU"
    matter = quiver.node_matter[node]

    n_f    = matter.get("fund", 0)
    n_fb   = matter.get("antifund", 0)
    n_S    = matter.get("S", 0)
    n_Sb   = matter.get("Sbar", 0)
    n_A    = matter.get("A", 0)
    n_Ab   = matter.get("Abar", 0)

    # Node matter contributions
    coeff_N     = (n_S - n_Sb) + (n_A - n_Ab)
    coeff_const = (n_f - n_fb) + 4 * (n_S - n_Sb) - 4 * (n_A - n_Ab)

    # Bifundamental contributions
    for e, side in quiver.incident_edges(node):
        neighbor = e.dst if side == "src" else e.src
        g_neighbor = quiver.gauge_types[neighbor]

        if g_neighbor == "SU":
            if e.rep == "+-":
                # (□_node, □̄_neighbor) if node is src; (□_neighbor, □̄_node) if dst
                coeff_N += +1 if side == "src" else -1
            elif e.rep == "++":
                coeff_N += +1  # always +N to both endpoints
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
    Count the total number of fundamental chiral multiplets at an Sp(N) node.
    Each bifundamental edge contributes 1; node-level funds contribute their count.

    The Witten anomaly condition: sp_fund_count ≡ 0 (mod 2).
    """
    assert quiver.gauge_types[node] == "Sp"
    n_f = quiver.node_matter[node].get("fund", 0)
    return n_f + quiver.degree(node)


def check_sp_witten(quiver: Quiver, node: int) -> bool:
    """Return True if the Sp(N) Witten anomaly vanishes at node."""
    return sp_fund_count(quiver, node) % 2 == 0


def check_anomalies(quiver: Quiver) -> bool:
    """Check all anomaly conditions for every node."""
    for i, g in enumerate(quiver.gauge_types):
        if g == "SU" and not check_su_anomaly(quiver, i):
            return False
        if g == "Sp" and not check_sp_witten(quiver, i):
            return False
    return True


# ── Asymptotic freedom ─────────────────────────────────────────────────────────

def check_af(quiver: Quiver) -> bool:
    """Return True if b_0 > 0 for all valid N at every node."""
    for i, g in enumerate(quiver.gauge_types):
        if not is_af_all_N(g, quiver.node_matter[i], quiver.bifund_neighbor_types(i)):
            return False
    return True


# ── Enumeration ────────────────────────────────────────────────────────────────

def _enumerate_node_matter(
    gauge_type: GaugeType,
    bifund_neighbors: list[GaugeType],
    max_rep: int = 6,
) -> Iterator[dict[str, int]]:
    """
    Yield all matter dicts for a node that satisfy b_0 > 0 for all valid N.

    Each representation count ranges from 0 to max_rep. Pruning: since b_0
    decreases monotonically in each count, if AF fails at some count > 0 we
    skip all higher counts of that rep.
    """
    reps = NODE_REPS[gauge_type]

    def _gen(idx: int, current: dict[str, int]) -> Iterator[dict[str, int]]:
        if idx == len(reps):
            if is_af_all_N(gauge_type, current, bifund_neighbors):
                yield dict(current)
            return
        rep = reps[idx]
        for count in range(max_rep + 1):
            if count > 0:
                current[rep] = count
                if not is_af_all_N(gauge_type, current, bifund_neighbors):
                    del current[rep]
                    break
            yield from _gen(idx + 1, current)
        if rep in current:
            del current[rep]

    yield from _gen(0, {})


def _enumerate_edges(
    gauge_types: list[GaugeType],
    max_multiedge: int = 2,
) -> Iterator[list[Edge]]:
    """
    Generate all edge multisets for a fixed gauge type assignment,
    up to max_multiedge edges of each type between each ordered pair of nodes.

    Yields lists of Edge objects (may be empty).
    """
    n = len(gauge_types)
    # Collect all (i, j, rep) slots where i <= j
    slots: list[tuple[int, int, str]] = []
    for i in range(n):
        for j in range(i, n):
            gi, gj = gauge_types[i], gauge_types[j]
            key = frozenset([gi, gj])
            if key not in EDGE_REPS:
                continue
            for rep in EDGE_REPS[key]:
                if rep == "+-":
                    # Directed: separate counts for i→j and j→i
                    slots.append((i, j, "+-_fwd"))  # i→j
                    if i != j:
                        slots.append((i, j, "+-_bwd"))  # j→i
                else:
                    slots.append((i, j, rep))

    # Enumerate multiplicities for each slot
    ranges = [range(max_multiedge + 1)] * len(slots)
    for counts in product(*ranges):
        edges: list[Edge] = []
        for (i, j, rep_key), count in zip(slots, counts):
            for _ in range(count):
                if rep_key == "+-_fwd":
                    edges.append(Edge(i, j, "+-"))
                elif rep_key == "+-_bwd":
                    edges.append(Edge(j, i, "+-"))
                else:
                    edges.append(Edge(i, j, rep_key))
        yield edges


def enumerate_quivers(
    max_nodes: int = 2,
    max_rep: int = 6,
    max_multiedge: int = 2,
    require_connected: bool = True,
) -> list[Quiver]:
    """
    Enumerate all quivers satisfying:
    - UV asymptotic freedom (b_0 > 0) for all valid N at every node
    - SU(N) gauge anomaly cancellation
    - Sp(N) Witten anomaly cancellation
    - Connected (if require_connected=True)

    Note: does not deduplicate isomorphic quivers. For small max_nodes the
    output can be filtered using a graph isomorphism library (e.g. networkx).

    Parameters
    ----------
    max_nodes : maximum number of gauge nodes
    max_rep : maximum multiplicity for each representation per node
    max_multiedge : maximum multiplicity of each edge type between a pair of nodes
    """
    gauge_group_types: list[GaugeType] = ["SU", "SO", "Sp"]
    valid: list[Quiver] = []

    for n_nodes in range(1, max_nodes + 1):
        # Enumerate all gauge type assignments
        for gauge_types in product(gauge_group_types, repeat=n_nodes):
            gauge_types = list(gauge_types)

            # Enumerate all edge configurations
            for edges in _enumerate_edges(gauge_types, max_multiedge):
                # Build partial quiver with no node matter to get neighbor lists
                q_partial = Quiver(gauge_types, edges)

                if require_connected and not q_partial.is_connected():
                    continue

                # Enumerate node matter for each node independently,
                # then combine (with constraint checking)
                matter_choices: list[list[dict[str, int]]] = []
                feasible = True
                for i in range(n_nodes):
                    neighbors = q_partial.bifund_neighbor_types(i)
                    choices = list(_enumerate_node_matter(gauge_types[i], neighbors, max_rep))
                    if not choices:
                        feasible = False
                        break
                    matter_choices.append(choices)

                if not feasible:
                    continue

                for matter_combo in product(*matter_choices):
                    q = Quiver(gauge_types, edges, list(matter_combo))
                    if check_anomalies(q):
                        valid.append(q)

    return valid


# ── Utilities ──────────────────────────────────────────────────────────────────

def quiver_summary(q: Quiver) -> str:
    """Human-readable summary of a quiver."""
    lines = []
    lines.append(f"Nodes: {' '.join(f'{g}({i})' for i, g in enumerate(q.gauge_types))}")
    for e in q.edges:
        lines.append(f"  Edge {e.src}→{e.dst}  rep={e.rep}")
    for i, matter in enumerate(q.node_matter):
        non_zero = {k: v for k, v in matter.items() if v > 0}
        if non_zero:
            lines.append(f"  Node {i} matter: {non_zero}")
    return "\n".join(lines)


# ── Quick validation ───────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Single SU(N) node: SQCD with N_f fund + N_f antifund
    # AF: 3N - N_f > 0 for all N>=2 => N_f <= 5
    # Anomaly: n_f - n_fbar = 0 => ok
    print("=== Single SU(N) node, SQCD-like ===")
    for n_f in range(7):
        q = Quiver(["SU"], [], [{"fund": n_f, "antifund": n_f}])
        af  = check_af(q)
        ano = check_anomalies(q)
        print(f"  N_f={n_f}: AF={af}  anomaly_free={ano}")

    print()

    # Two SU(N) nodes, one "+-" bifundamental
    # AF: b_0 = 3N - N/2 = 5N/2 > 0 ✓; anomaly: +N at node 0, -N at node 1 → not balanced
    print("=== Two SU(N) nodes, one +- edge, no extra matter ===")
    q = Quiver(["SU", "SU"], [Edge(0, 1, "+-")])
    print(f"  AF={check_af(q)}  anomaly_free={check_anomalies(q)}")
    coeff = su_anomaly_polynomial(q, 0)
    print(f"  SU anomaly at node 0: {coeff[0]}*N + {coeff[1]}")

    print()

    # Two SU(N) nodes, two +- edges in opposite directions: balanced
    # b_0 = 3N - N = 2N > 0 ✓; anomaly: +N - N = 0 ✓
    print("=== Two SU(N) nodes, one +- each direction (balanced) ===")
    q = Quiver(["SU", "SU"], [Edge(0, 1, "+-"), Edge(1, 0, "+-")])
    print(f"  AF={check_af(q)}  anomaly_free={check_anomalies(q)}")

    print()

    # Sp(N) node with 2 fund: Witten OK
    print("=== Sp(N) node, 2 fundamentals ===")
    q = Quiver(["Sp"], [], [{"fund": 2}])
    print(f"  AF={check_af(q)}  Witten OK={check_sp_witten(q, 0)}")

    print()

    # Sp(N) node with 3 fund: Witten violation
    print("=== Sp(N) node, 3 fundamentals ===")
    q = Quiver(["Sp"], [], [{"fund": 3}])
    print(f"  AF={check_af(q)}  Witten OK={check_sp_witten(q, 0)}")

    print()

    # Small enumeration: 1 node, no multiedge
    print("=== Enumerate 1-node quivers ===")
    results = enumerate_quivers(max_nodes=1, max_rep=6, max_multiedge=0)
    print(f"  Found {len(results)} valid 1-node quivers")
