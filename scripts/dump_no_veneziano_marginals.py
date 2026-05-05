"""SU(N)-SU(N) no-Veneziano theories: matter (with delta-induced fund/antifund),
AF condition at O point, default-mode marginal operators per theory.

The matter-richer node is placed in the second column. AF bounds use ≤ when
the theory admits at least one always-marginal operator at b_0 = 0
(max-N_f mode), strictly < otherwise.
"""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from marginal_operators import (
    quiver_from_row, find_marginal_ops, find_marginal_ops_max_Nf,
    op_short, _r_values_sane, r_values_max_Nf,
)
from quiver_generation import Edge, Quiver

DB = Path(__file__).parent.parent / "quivers.db"
OUT = Path(__file__).parent.parent / "paper/sections/generated/no_veneziano_marginals.tex"

N_LIST = (10, 20, 30)
MAX_DEGREE = 6


# ── matter formatting ─────────────────────────────────────────────────────────

_REP_TEX = {
    "adj": r"\mathrm{adj}",
    "S":   "S",   "Sbar": r"\bar{S}",
    "A":   "A",   "Abar": r"\bar{A}",
}
_REP_ORDER = ("adj", "S", "Sbar", "A", "Abar")


def _matter_dict_to_tex(matter: dict[str, int]) -> list[str]:
    parts = []
    for rep in _REP_ORDER:
        n = matter.get(rep, 0)
        if n == 0:
            continue
        sym = _REP_TEX[rep]
        parts.append(sym if n == 1 else f"{n}\\,{sym}")
    return parts


def _delta_mag_tex(a: int, b: int) -> tuple[str, int]:
    """Return (magnitude TeX, sign) for |a*N + b| at large N.
    sign is +1 (fund), -1 (antifund), 0 (no extras)."""
    if a == 0 and b == 0:
        return "", 0
    if a > 0 or (a == 0 and b > 0):
        sign = +1
    else:
        sign = -1
        a, b = -a, -b
    parts = []
    if a:
        parts.append("N" if a == 1 else f"{a}N")
    if b:
        if not parts:
            parts.append(str(b))
        elif b > 0:
            parts.append(f"+{b}")
        else:
            parts.append(str(b))
    return "".join(parts), sign


def _matter_with_delta_tex(matter: dict[str, int], delta_a: int, delta_b: int) -> str:
    parts = _matter_dict_to_tex(matter)
    mag, sign = _delta_mag_tex(delta_a, delta_b)
    if sign != 0:
        rep = r"\square" if sign == +1 else r"\bar{\square}"
        # If mag has no N (purely numeric) and equals "1", drop it
        prefix = "" if mag in ("", "1") else mag + r"\,"
        # If mag contains "+" or "-" inside, wrap in parens
        if any(c in mag for c in "+-") and mag and mag != "":
            prefix = f"({mag})\\,"
        parts.append(f"{prefix}{rep}")
    return " + ".join(parts) if parts else "-"


def _matter_count(matter: dict[str, int], delta_a: int, delta_b: int) -> int:
    """Sort key: total field count at the node (rank-2 + |delta| at large N)."""
    base = sum(matter.values())
    # |delta| as N-coeff weighted (larger a means more fundamentals at large N)
    extra = abs(delta_a) * 100 + abs(delta_b)  # weight a heavily, tie-break on b
    return base * 1000 + extra


# ── edges formatting ──────────────────────────────────────────────────────────

_EDGE_REP_TEX = {
    "+-": r"\square\bar{\square}",
    "++": r"\square\square",
    "--": r"\bar{\square}\bar{\square}",
}


def _edges_tex(edges: list[Edge], swap: bool) -> str:
    """Format edges. If swap, relabel node 0↔1 (and flip src/dst)."""
    from collections import Counter
    c: Counter = Counter()
    for e in edges:
        src, dst = (e.dst, e.src) if swap else (e.src, e.dst)
        c[(src, dst, e.rep)] += 1
    parts = []
    for (src, dst, rep), n in sorted(c.items()):
        rep_t = _EDGE_REP_TEX.get(rep, rep)
        edge_t = f"({src+1}\\!\\to\\!{dst+1},\\,{rep_t})"
        parts.append(edge_t if n == 1 else f"{n}\\!\\times\\!{edge_t}")
    return ",\\;".join(parts) if parts else "-"


# ── nf bound formatting ──────────────────────────────────────────────────────

def _bound_tex_from_db(nf_bound_str: str, comparison: str) -> str:
    """Convert DB string 'N_f ≤ 2*N - 1' → '$N_f {≤,<} 2N-1$'."""
    if not nf_bound_str:
        return "-"
    import re as _re
    s = nf_bound_str.strip()
    s = s.replace("≤", comparison)
    s = s.replace("*N", "N").replace("* N", "N")
    s = _re.sub(r"\b1N\b", "N", s)
    s = s.replace("N_f", "N_f")
    return f"${s}$"


# ── operator pretty-printing ──────────────────────────────────────────────────

def _op_tex(op) -> str:
    """Compact LaTeX for a CandidateOp."""
    from marginal_operators import _label_short
    if op.kind == "bifund-loop" and op.word is not None:
        parts = [_op_label_tex(lbl) for lbl in op.word]
        return r"\mathrm{tr}(" + r"\,".join(parts) + ")"
    parts = []
    for label, mult in op.factors:
        s = _op_label_tex(label)
        if mult > 1:
            s += f"^{{{mult}}}"
        parts.append(s)
    return r"\mathrm{tr}(" + r"\,".join(parts) + ")"


def _op_label_tex(label: str) -> str:
    if label.startswith("node"):
        node_idx, rep = label[4:].split("_", 1)
        node = int(node_idx) + 1
        rep_t = {
            "adj": r"\mathrm{adj}", "S": "S", "Sbar": r"\bar{S}",
            "A": "A", "Abar": r"\bar{A}", "fund": r"\square",
            "antifund": r"\bar{\square}",
        }.get(rep, rep)
        return f"{rep_t}_{{{node}}}"
    parts = label.split("_")
    src, dst, rep = parts[1], parts[2], "_".join(parts[3:])
    rep_t = _EDGE_REP_TEX.get(rep, rep)
    return f"Q^{{{rep_t}}}_{{{int(src)+1}{int(dst)+1}}}"


def _ops_tex(ops, swap: bool, max_chars: int = 80) -> str:
    if not ops:
        return "-"
    if swap:
        ops = [_swap_op(op) for op in ops]
    strs = [_op_tex(op) for op in ops]
    # cap the number shown
    out = ";\\; ".join(f"${s}$" for s in strs)
    if len(out) > max_chars:
        kept, used = [], 0
        for s in strs:
            if used + len(s) + 4 > max_chars - 12:
                break
            kept.append(s)
            used += len(s) + 4
        out = ";\\; ".join(f"${s}$" for s in kept) + f";\\; \\ldots+{len(strs) - len(kept)}"
    return out


def _swap_op(op):
    """Return a CandidateOp with all node indices in labels swapped 0↔1."""
    from marginal_operators import CandidateOp
    def swap_label(lbl):
        if lbl.startswith("node"):
            node, rest = lbl[4:].split("_", 1)
            return f"node{1 - int(node)}_{rest}"
        # edge_{src}_{dst}_{rep}
        parts = lbl.split("_")
        src, dst = 1 - int(parts[1]), 1 - int(parts[2])
        return f"edge_{src}_{dst}_{'_'.join(parts[3:])}"
    new_factors = tuple(sorted((swap_label(lbl), m) for lbl, m in op.factors))
    new_word = tuple(swap_label(l) for l in op.word) if op.word else None
    return CandidateOp(
        factors=new_factors, kind=op.kind, word=new_word,
        degree=op.degree, flavor_signature=op.flavor_signature,
    )


# ── max-N_f marginal existence ────────────────────────────────────────────────

def _max_nf_has_marginal(q: Quiver) -> bool:
    """Return True iff at least one always-marginal operator exists at b_0=0
    saturation. Returns False if max-N_f a-max diverges at any N."""
    try:
        ops, _ = find_marginal_ops_max_Nf(q, N_list=N_LIST, max_degree=MAX_DEGREE)
        return bool(ops)
    except Exception:
        return False


# ── main loop ────────────────────────────────────────────────────────────────

def main() -> None:
    con = sqlite3.connect(DB)
    con.row_factory = sqlite3.Row
    rows = list(con.execute(
        "SELECT * FROM theory WHERE gauge_pair='SU-SU' AND veneziano=0 AND class_id IS NOT NULL "
        "ORDER BY class_id, theory_id"
    ))
    print(f"Processing {len(rows)} SU-SU no-Veneziano theories...", file=sys.stderr)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w") as f:
        f.write("% Auto-generated by scripts/dump_no_veneziano_marginals.py\n")
        f.write("% Requires: \\usepackage{booktabs,longtable,amssymb,amsmath}\n")
        f.write("\\section*{SU(N)-SU(N) non-Veneziano theories: marginal operators}\n\n")
        f.write(
            "Matter content (with anomaly-required $\\square$/$\\bar{\\square}$ from $\\delta$), "
            "edges, asymptotic-freedom condition at $\\mathcal{O}=(g_1=0,g_2=0)$, "
            "and always-marginal operators at the IR fixed point (default mode, $|R-2|<10^{-6}$ "
            "at $N\\in\\{10,20,30\\}$). The node with more matter is placed second. "
            "AF bound uses $\\le$ when a marginal operator exists at $b_0=0$ saturation, "
            "$<$ otherwise.\n\n"
        )
        f.write("\\begin{longtable}{rrlllllp{0.30\\linewidth}}\n")
        f.write("\\toprule\n")
        f.write("class & th. & matter (1) & matter (2) & edges & AF (1) & AF (2) & marginal ops \\\\\n")
        f.write("\\midrule\n\\endhead\n")

        n_processed = 0
        for row in rows:
            r = dict(row)
            try:
                q = quiver_from_row(r)
            except Exception:
                continue

            d0a = r["delta0_a"] or 0
            d0b = r["delta0_b"] or 0
            d1a = r["delta1_a"] or 0
            d1b = r["delta1_b"] or 0

            count0 = _matter_count(q.node_matter[0], d0a, d0b)
            count1 = _matter_count(q.node_matter[1], d1a, d1b)
            swap = count0 > count1

            if swap:
                m_first  = _matter_with_delta_tex(q.node_matter[1], d1a, d1b)
                m_second = _matter_with_delta_tex(q.node_matter[0], d0a, d0b)
                af_first  = r["nf_bound1"]
                af_second = r["nf_bound0"]
            else:
                m_first  = _matter_with_delta_tex(q.node_matter[0], d0a, d0b)
                m_second = _matter_with_delta_tex(q.node_matter[1], d1a, d1b)
                af_first  = r["nf_bound0"]
                af_second = r["nf_bound1"]

            edges_tex = _edges_tex(q.edges, swap)
            comparison = "\\le" if _max_nf_has_marginal(q) else "<"
            af1_tex = _bound_tex_from_db(af_first, comparison)
            af2_tex = _bound_tex_from_db(af_second, comparison)

            try:
                ops = find_marginal_ops(q, N_list=N_LIST, max_degree=MAX_DEGREE)
            except Exception:
                ops = []
            ops_tex = _ops_tex(ops, swap)

            f.write(
                f"{r['class_id']} & {r['theory_id']} & "
                f"${m_first}$ & ${m_second}$ & ${edges_tex}$ & "
                f"{af1_tex} & {af2_tex} & {ops_tex} \\\\\n"
            )
            n_processed += 1
            if n_processed % 100 == 0:
                print(f"  ... {n_processed}/{len(rows)}", file=sys.stderr)
        f.write("\\bottomrule\n\\end{longtable}\n")

    con.close()
    print(f"Wrote {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
