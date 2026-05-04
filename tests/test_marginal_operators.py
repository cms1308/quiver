"""Acceptance tests for marginal_operators module."""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from quiver_generation import Edge, Quiver
from marginal_operators import (
    CandidateOp,
    _flavor_multiplets,
    _index_balance_node,
    enumerate_candidates,
    is_flavor_singlet,
    is_marginal_at_all_N,
    op_R_at_N,
    parse_edges,
    parse_matter,
    quiver_from_row,
    r_values_at_N,
)


# ── Index balance ──────────────────────────────────────────────────────────────

def test_index_balance_su_S_Sbar_balanced():
    q = Quiver(["SU"], node_matter=[{"S": 1, "Sbar": 1}])
    assert _index_balance_node(q, 0, {"node0_S": 1, "node0_Sbar": 1})


def test_index_balance_su_unbalanced():
    q = Quiver(["SU"], node_matter=[{"S": 1, "fund": 1}])
    # S contributes (2,0), fund contributes (1,0) → 3 upper, 0 lower
    assert not _index_balance_node(q, 0, {"node0_S": 1, "node0_fund": 2})


def test_index_balance_su_adj_balanced():
    q = Quiver(["SU"], node_matter=[{"adj": 1}])
    assert _index_balance_node(q, 0, {"node0_adj": 3})


def test_index_balance_so_V_pair():
    q = Quiver(["SO"], node_matter=[{"V": 1}])
    # 2 V's: total 2 indices → even → balanced
    assert _index_balance_node(q, 0, {"node0_V": 2})


def test_index_balance_so_V_single_fail():
    q = Quiver(["SO"], node_matter=[{"V": 1}])
    # 1 V: total 1 → odd → not balanced
    assert not _index_balance_node(q, 0, {"node0_V": 1})


def test_index_balance_sp_fund_odd_fail():
    q = Quiver(["Sp"], node_matter=[{"fund": 1}])
    assert not _index_balance_node(q, 0, {"node0_fund": 1})
    assert _index_balance_node(q, 0, {"node0_fund": 2})


# ── Flavor multiplets / singlet ────────────────────────────────────────────────

def test_flavor_multiplets_bif_count():
    q = Quiver(
        ["SU", "SU"],
        edges=[Edge(0, 1, "+-"), Edge(0, 1, "+-")],
        node_matter=[{}, {}],
    )
    mults = _flavor_multiplets(q)
    assert mults == {"edge_0_1_pm": 2}


def test_flavor_singlet_rule():
    multiplets = {"edge_0_1_pm": 2}

    op_one = CandidateOp(
        factors=(("edge_0_1_pm", 1),), kind="bifund-loop", word=None,
        degree=1, flavor_signature={"edge_0_1_pm": 1},
    )
    assert not is_flavor_singlet(op_one, multiplets)

    op_two = CandidateOp(
        factors=(("edge_0_1_pm", 2),), kind="bifund-loop", word=None,
        degree=2, flavor_signature={"edge_0_1_pm": 2},
    )
    assert is_flavor_singlet(op_two, multiplets)

    op_zero = CandidateOp(
        factors=(("node0_adj", 2),), kind="single-node", word=None,
        degree=2, flavor_signature={},
    )
    assert is_flavor_singlet(op_zero, multiplets)


# ── Parsers ────────────────────────────────────────────────────────────────────

def test_parse_matter_examples():
    assert parse_matter("") == {}
    assert parse_matter("—") == {}
    assert parse_matter("Ā") == {"Abar": 1}
    assert parse_matter("S + 2Ā") == {"S": 1, "Abar": 2}
    assert parse_matter("adj + S̄") == {"adj": 1, "Sbar": 1}
    assert parse_matter("2adj + S + Ā") == {"adj": 2, "S": 1, "Abar": 1}


def test_parse_edges_examples():
    assert parse_edges("") == []
    assert parse_edges("(0→1,++)") == [Edge(0, 1, "++")]
    assert parse_edges("2×(0→1,+-)") == [Edge(0, 1, "+-"), Edge(0, 1, "+-")]
    parsed = parse_edges("(0→1,+-)  2×(0→1,++)")
    assert parsed == [Edge(0, 1, "+-"), Edge(0, 1, "++"), Edge(0, 1, "++")]


# ── End-to-end smoke test on DB ────────────────────────────────────────────────

def test_db_roundtrip_theory_1():
    import sqlite3
    db_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "quivers.db",
    )
    if not os.path.exists(db_path):
        return  # skip when DB unavailable
    con = sqlite3.connect(db_path)
    con.row_factory = sqlite3.Row
    row = con.execute("SELECT * FROM theory WHERE theory_id=1").fetchone()
    q = quiver_from_row(dict(row))
    assert q.gauge_types == ["SU", "SU"]
    assert q.node_matter == [{"Abar": 1}, {"Abar": 1}]
    assert q.edges == [Edge(0, 1, "++")]


def test_op_R_evaluation():
    # SU-SU theory with adj at node 1 + 2 bifunds (++, --)
    q = Quiver(
        ["SU", "SU"],
        edges=[Edge(0, 1, "++"), Edge(0, 1, "--")],
        node_matter=[{"A": 1, "Abar": 1}, {"adj": 1}],
    )
    cands = enumerate_candidates(q, max_degree=4)
    assert len(cands) > 0
    R10 = r_values_at_N(q, 10)
    # Every candidate must evaluate (no missing labels)
    for op in cands:
        r = op_R_at_N(op, R10)
        assert r is not None, f"{op.factors} produced None"


def test_max_nf_saturation():
    # SU-SU with one ++ bifund and Ā at each node: not conformal at default
    # N_f=0; saturating each node to b_0=0 should yield R=2/3 for all matter
    # (free SQCD-like) and produce cubic mesonic marginals.
    q = Quiver(
        ["SU", "SU"],
        edges=[Edge(0, 1, "++")],
        node_matter=[{"Abar": 1}, {"Abar": 1}],
    )
    from marginal_operators import (
        find_marginal_ops_max_Nf, nf_max_per_node, r_values_max_Nf,
    )
    Nf = nf_max_per_node(q, 10)
    assert all(n > 0 for n in Nf), f"expected positive N_f at every node, got {Nf}"
    R, _ = r_values_max_Nf(q, 10)
    # All matter should have R ≈ 2/3 at the b_0=0 saturation
    for label, r in R.items():
        assert abs(r - 2/3) < 1e-6, f"{label}: R={r}, expected 2/3"
    ops, _ = find_marginal_ops_max_Nf(q, max_degree=4)
    assert len(ops) >= 1, "expected at least one marginal operator at max-N_f"


def test_marginal_filter_consistency():
    # The gauge-anomaly Konishi-like operator at any conformal SU node has
    # multiplicities (n_i) ∝ T_a(rep_i) and should be marginal at every N
    # by construction. Verify the filter catches at least one such candidate.
    q = Quiver(
        ["SU", "SU"],
        edges=[Edge(0, 1, "++"), Edge(0, 1, "--")],
        node_matter=[{"A": 1, "Abar": 1}, {"adj": 1}],
    )
    R_per_N = {N: r_values_at_N(q, N) for N in (10, 20, 30)}
    cands = enumerate_candidates(q, max_degree=4)
    marginal_count = sum(1 for op in cands if is_marginal_at_all_N(op, R_per_N))
    assert marginal_count >= 1


if __name__ == "__main__":
    import sys
    failed = 0
    for name in sorted(globals()):
        if not name.startswith("test_"):
            continue
        fn = globals()[name]
        try:
            fn()
            print(f"  PASS  {name}")
        except Exception as e:
            print(f"  FAIL  {name}: {e}")
            failed += 1
    sys.exit(1 if failed else 0)
