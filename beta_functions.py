"""
beta_functions.py — Module 1

One-loop beta function coefficients for 4D N=1 quiver gauge theories.

Reference: module1_beta_functions.md
"""

from fractions import Fraction
from dataclasses import dataclass, field
from typing import Literal

GaugeType = Literal["SU", "SO", "Sp"]

# Minimum valid N for each gauge group type
N_MIN: dict[str, int] = {"SU": 5, "SO": 5, "Sp": 3}

# Valid node-level matter representations per gauge group type (all reps)
VALID_REPS: dict[str, set[str]] = {
    "SU": {"fund", "antifund", "adj", "S", "Sbar", "A", "Abar"},
    "SO": {"V", "adj", "S"},
    "Sp": {"fund", "adj", "A"},
}

# Rank-2 and adjoint representations only (fund-like reps absorbed into N_f)
RANK2_ADJ_REPS: dict[str, list[str]] = {
    "SU": ["adj", "S", "Sbar", "A", "Abar"],
    "SO": ["adj", "S"],
    "Sp": ["adj", "A"],
}


# ── Dynkin indices ─────────────────────────────────────────────────────────────

def T_adj(gauge_type: GaugeType, N: int) -> Fraction:
    """Dynkin index of the adjoint representation."""
    if gauge_type == "SU":
        return Fraction(N)
    if gauge_type == "SO":
        return Fraction(N - 2)
    if gauge_type == "Sp":
        return Fraction(N + 1)
    raise ValueError(f"Unknown gauge type: {gauge_type!r}")


def T_rep(gauge_type: GaugeType, rep: str, N: int) -> Fraction:
    """Dynkin index T(r) for a node-level matter representation."""
    if gauge_type == "SU":
        if rep in ("fund", "antifund"):
            return Fraction(1, 2)
        if rep == "adj":
            return Fraction(N)
        if rep in ("S", "Sbar"):
            return Fraction(N + 2, 2)
        if rep in ("A", "Abar"):
            return Fraction(N - 2, 2)
    elif gauge_type == "SO":
        if rep == "V":
            return Fraction(1)
        if rep == "adj":
            return Fraction(N - 2)
        if rep == "S":
            return Fraction(N + 2)
    elif gauge_type == "Sp":
        if rep == "fund":
            return Fraction(1, 2)
        if rep == "adj":
            return Fraction(N + 1)
        if rep == "A":
            return Fraction(N - 1)
    raise ValueError(f"Unknown rep {rep!r} for gauge type {gauge_type!r}")


def _fund_dim(gauge_type: GaugeType, N: int) -> int:
    """Dimension of the fundamental/vector/spinor representation at each endpoint."""
    if gauge_type == "SU":
        return N          # □ of SU(N)
    if gauge_type == "SO":
        return N          # V of SO(N)
    if gauge_type == "Sp":
        return 2 * N      # f of USp(2N)
    raise ValueError(f"Unknown gauge type: {gauge_type!r}")


def _fund_T(gauge_type: GaugeType) -> Fraction:
    """Dynkin index T(fund-like rep) — exact, independent of N."""
    if gauge_type == "SU":
        return Fraction(1, 2)
    if gauge_type == "SO":
        return Fraction(1)
    if gauge_type == "Sp":
        return Fraction(1, 2)
    raise ValueError(f"Unknown gauge type: {gauge_type!r}")


def T_bifund(
    gauge_a: GaugeType, gauge_b: GaugeType, N: int,
    N_a: int | None = None, N_b: int | None = None,
) -> tuple[Fraction, Fraction]:
    """
    Bifundamental Dynkin index contributions (T_a, T_b) to b_0 at each endpoint
    for one bifundamental edge between a node of type gauge_a (rank N_a) and
    one of gauge_b (rank N_b).

    When N_a and N_b are omitted both default to N (equal-rank case).
    Formula: T_a = T_fund(gauge_a) × dim_fund(gauge_b, N_b), and vice versa.
    """
    if N_a is None:
        N_a = N
    if N_b is None:
        N_b = N
    T_a = _fund_T(gauge_a) * _fund_dim(gauge_b, N_b)
    T_b = _fund_T(gauge_b) * _fund_dim(gauge_a, N_a)
    return T_a, T_b


# ── Anomaly coefficients (SU(N) only) ─────────────────────────────────────────

def A_rep(rep: str, N: int) -> int:
    """
    Cubic anomaly coefficient A(r) for SU(N) representations.
    The gauge anomaly at an SU(N) node vanishes iff sum_i A(r_i) * dim(other) = 0.
    """
    if rep == "fund":
        return +1
    if rep == "antifund":
        return -1
    if rep == "adj":
        return 0
    if rep == "S":
        return +(N + 4)
    if rep == "Sbar":
        return -(N + 4)
    if rep == "A":
        return +(N - 4)
    if rep == "Abar":
        return -(N - 4)
    raise ValueError(f"Unknown SU(N) rep for anomaly: {rep!r}")


# ── Beta function computation ──────────────────────────────────────────────────

@dataclass
class NodeSpec:
    """
    Specification of a single gauge node.

    Attributes
    ----------
    gauge_type : "SU", "SO", or "Sp"
    N : gauge group rank parameter
    matter : dict mapping representation name to multiplicity
        SU(N): "fund", "antifund", "adj", "S", "Sbar", "A", "Abar"
        SO(N): "V", "adj", "S"
        Sp(N): "fund", "adj", "A"
    bifund_neighbors : list of gauge types of neighboring nodes,
        one entry per bifundamental edge (with multiplicity)
    """
    gauge_type: GaugeType
    N: int
    matter: dict[str, int] = field(default_factory=dict)
    bifund_neighbors: list[GaugeType] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.N < N_MIN[self.gauge_type]:
            raise ValueError(
                f"{self.gauge_type}(N) requires N >= {N_MIN[self.gauge_type]}, got N={self.N}"
            )
        for rep in self.matter:
            if rep not in VALID_REPS[self.gauge_type]:
                raise ValueError(
                    f"Rep {rep!r} is not valid for {self.gauge_type}(N). "
                    f"Valid: {VALID_REPS[self.gauge_type]}"
                )


def compute_b0(node: NodeSpec) -> Fraction:
    """
    Compute the one-loop beta function coefficient:

        b_0 = 3 T(adj) - sum_i T(r_i)

    where the sum runs over all chiral multiplets charged under the node,
    including bifundamental contributions from neighboring nodes.

    Asymptotic freedom requires b_0 > 0.
    """
    N, g = node.N, node.gauge_type
    b0 = 3 * T_adj(g, N)
    for rep, count in node.matter.items():
        b0 -= count * T_rep(g, rep, N)
    for neighbor in node.bifund_neighbors:
        T_a, _ = T_bifund(g, neighbor, N)
        b0 -= T_a
    return b0


def b0_linear(
    gauge_type: GaugeType,
    matter: dict[str, int],
    bifund_neighbors: list[GaugeType],
    rank_mult: int = 1,
    neighbor_mults: list[int] | None = None,
) -> tuple[Fraction, Fraction]:
    """
    Return (alpha, beta) such that b_0 = alpha*N_base + beta, where N_base is
    the shared large-N parameter and each node has rank rank_mult * N_base.

    neighbor_mults[i] is the rank multiplier for the i-th bifundamental neighbor;
    defaults to rank_mult for each neighbor (equal-rank case).

    Since b_0 is linear in N_base, two evaluations suffice.
    """
    if neighbor_mults is None:
        neighbor_mults = [rank_mult] * len(bifund_neighbors)
    N_min_base = N_MIN[gauge_type]

    def _b0_at(N_base: int) -> Fraction:
        N_a = rank_mult * N_base
        b0 = 3 * T_adj(gauge_type, N_a)
        for rep, count in matter.items():
            b0 -= count * T_rep(gauge_type, rep, N_a)
        for g_nb, m_nb in zip(bifund_neighbors, neighbor_mults):
            N_b = m_nb * N_base
            T_a, _ = T_bifund(gauge_type, g_nb, N_a, N_a, N_b)
            b0 -= T_a
        return b0

    b0_lo = _b0_at(N_min_base)
    b0_hi = _b0_at(N_min_base + 1)
    alpha = b0_hi - b0_lo
    beta  = b0_lo - alpha * N_min_base
    return alpha, beta


def is_af(node: NodeSpec) -> bool:
    """Return True if b_0 > 0 at the given N."""
    return compute_b0(node) > 0


def is_af_all_N(
    gauge_type: GaugeType,
    matter: dict[str, int],
    bifund_neighbors: list[GaugeType],
    rank_mult: int = 1,
    neighbor_mults: list[int] | None = None,
) -> bool:
    """
    Return True if b_0 is positive for all sufficiently large N (large-N classification).

    For large N, b_0(N) = alpha*N + beta > 0 requires:
      - alpha > 0, or
      - alpha == 0 and beta > 0  (constant positive b_0)
    """
    alpha, beta = b0_linear(gauge_type, matter, bifund_neighbors, rank_mult, neighbor_mults)
    return alpha > 0 or (alpha == 0 and beta >= 0)


# ── Examples ───────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Example A: SU(N) SQCD with N_f fundamentals + N_f anti-fundamentals
    # b_0 = 3N - N_f  =>  AF iff N_f < 3N
    for N_f in [2, 5, 9]:
        node = NodeSpec("SU", 5, matter={"fund": N_f, "antifund": N_f})
        print(f"SQCD N=3, N_f={N_f}: b_0 = {compute_b0(node)}  AF={is_af(node)}")

    print()

    # Example B: N=2 SQCD — SU(N) with 1 adjoint + N_f fund + N_f antifund
    # b_0 = 2N - N_f  =>  AF iff N_f < 2N
    for N_f in [1, 4, 6]:
        node = NodeSpec("SU", 5, matter={"adj": 1, "fund": N_f, "antifund": N_f})
        print(f"N=2 SQCD N=4, N_f={N_f}: b_0 = {compute_b0(node)}  AF={is_af(node)}")

    print()

    # AF for all N check: SU with 1 adjoint (b_0 = 2N > 0 for all N >= 2)
    print("SU + adj, AF all N:", is_af_all_N("SU", {"adj": 1}, []))

    # SU with 3 antifund: b_0 = 3N - 3/2, alpha=3 > 0, b_0(2)=9/2 > 0 => AF all N
    print("SU + 3 antifund, AF all N:", is_af_all_N("SU", {"antifund": 3}, []))

    # Bifundamental: SU(N)-SU(N) circular quiver (2 SU neighbors each)
    # b_0 = 3N - 2*(N/2) = 2N > 0 for all N >= 2
    alpha, beta = b0_linear("SU", {}, ["SU", "SU"])
    print(f"SU with 2 SU bifund neighbors: b_0 = {alpha}*N + {beta}")
