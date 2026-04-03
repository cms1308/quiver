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


def T_bifund(gauge_a: GaugeType, gauge_b: GaugeType, N: int) -> tuple[Fraction, Fraction]:
    """
    Bifundamental Dynkin index contributions (T_a, T_b) to b_0 at each endpoint
    for one bifundamental edge between a node of type gauge_a and one of gauge_b.
    """
    pair = (gauge_a, gauge_b)
    if pair == ("SU", "SU"):
        return Fraction(N, 2), Fraction(N, 2)
    if pair in (("SU", "SO"), ("SO", "SU")):
        T_SU, T_SO = Fraction(N, 2), Fraction(N)
        return (T_SU, T_SO) if gauge_a == "SU" else (T_SO, T_SU)
    if pair in (("SU", "Sp"), ("Sp", "SU")):
        T_SU, T_Sp = Fraction(N), Fraction(N, 2)
        return (T_SU, T_Sp) if gauge_a == "SU" else (T_Sp, T_SU)
    if pair == ("SO", "SO"):
        return Fraction(N), Fraction(N)
    if pair in (("SO", "Sp"), ("Sp", "SO")):
        T_SO, T_Sp = Fraction(2 * N), Fraction(N, 2)
        return (T_SO, T_Sp) if gauge_a == "SO" else (T_Sp, T_SO)
    if pair == ("Sp", "Sp"):
        return Fraction(N), Fraction(N)
    raise ValueError(f"Unknown gauge pair: ({gauge_a!r}, {gauge_b!r})")


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
) -> tuple[Fraction, Fraction]:
    """
    Return (alpha, beta) such that b_0(N) = alpha*N + beta.

    Since b_0 is linear in N, two evaluations at N_min and N_min+1 suffice.
    """
    N_min = N_MIN[gauge_type]
    b0_lo = compute_b0(NodeSpec(gauge_type, N_min,     matter, bifund_neighbors))
    b0_hi = compute_b0(NodeSpec(gauge_type, N_min + 1, matter, bifund_neighbors))
    alpha = b0_hi - b0_lo
    beta  = b0_lo - alpha * N_min
    return alpha, beta


def is_af(node: NodeSpec) -> bool:
    """Return True if b_0 > 0 at the given N."""
    return compute_b0(node) > 0


def is_af_all_N(
    gauge_type: GaugeType,
    matter: dict[str, int],
    bifund_neighbors: list[GaugeType],
) -> bool:
    """
    Return True if b_0 is positive for all sufficiently large N (large-N classification).

    For large N, b_0(N) = alpha*N + beta > 0 requires:
      - alpha > 0, or
      - alpha == 0 and beta > 0  (constant positive b_0)
    """
    alpha, beta = b0_linear(gauge_type, matter, bifund_neighbors)
    return alpha > 0 or (alpha == 0 and beta > 0)


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
