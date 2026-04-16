"""Dump non-Veneziano universality classes from quivers.db into TeX tables, grouped by gauge pair."""

import re
import sqlite3
from pathlib import Path

DB = Path(__file__).parent / "quivers.db"
OUT = Path(__file__).parent / "no_veneziano_tables.tex"

GAUGE_PAIRS = ["SU-SU", "SU-SO", "SU-Sp", "SO-SO", "SO-Sp", "Sp-Sp"]


def _fractionize(s: str) -> str:
    """Replace p/q patterns with \\frac{p}{q}, handling sqrt terms in numerators."""
    def _repl(m):
        sign, num, den = m.group(1) or "", m.group(2), m.group(3)
        return f"{sign}\\frac{{{num}}}{{{den}}}"
    s = re.sub(
        r"([-+]\s*)?([0-9]*\s*\\sqrt\{[^}]*\}|[0-9]+)\s*/\s*([0-9]+)",
        _repl, s
    )
    return s


def texify_expr(s: str) -> str:
    if s is None or s == "":
        return "---"
    s = s.replace("*", "")
    s = re.sub(r"√\(([^()]*)\)", r"\\sqrt{\1}", s)
    s = _fractionize(s)
    return f"${s}$"


_RANK2_REPS = {"A", "Ā", "S", "S̄", r"\bar{A}", r"\bar{S}", "adj"}


def texify_R(s: str) -> str:
    if not s:
        return "---"
    s = s.replace("*", "")
    # Relabel nodes: 0->1, 1->2
    s = re.sub(r"_0\b", "_@@TEMP@@", s)
    s = re.sub(r"_1\b", "_2", s)
    s = s.replace("_@@TEMP@@", "_1")

    # Parse into (name, node, value) triples before TeX conversion
    raw_parts = [p.strip() for p in s.split(",")]
    entries = []  # (rep, node, value)  where node is "1","2" or None (bif)
    for p in raw_parts:
        m = re.match(r"R_(\S+?)_([12])\s*=\s*(.*)", p)
        if m:
            entries.append((m.group(1), m.group(2), m.group(3)))
        else:
            m2 = re.match(r"R_(\S+?)\s*=\s*(.*)", p)
            if m2:
                entries.append((m2.group(1), None, m2.group(2)))
            else:
                entries.append((None, None, p))

    # Group rank-2 reps by node; always write as R_{r2,node}
    nodes_done = set()
    output_parts = []
    for rep, node, val in entries:
        if rep in _RANK2_REPS and node is not None:
            if node not in nodes_done:
                nodes_done.add(node)
                output_parts.append(f"R_{{r2,{node}}}={val}")
        else:
            if node is not None:
                output_parts.append(f"R_{{{rep}_{node}}}={val}")
            elif rep is not None:
                output_parts.append(f"R_{{{rep}}}={val}")
            else:
                output_parts.append(val)

    # Now apply TeX transformations
    tex_parts = []
    for p in output_parts:
        p = re.sub(r"√\(([^()]*)\)", r"\\sqrt{\1}", p)
        p = p.replace("Ā", r"\bar{A}")
        p = p.replace("S̄", r"\bar{S}")
        p = _fractionize(p)
        tex_parts.append(p)

    return "$\\begin{array}{l}" + " \\\\ ".join(tex_parts) + "\\end{array}$"


def texify_matter(matter0, matter1, gauge0, gauge1):
    """Format matter content for display."""
    parts = []
    if matter0 and matter0 != "—":
        parts.append(matter0)
    if matter1 and matter1 != "—":
        parts.append(matter1)
    if not parts:
        return "---"
    return ", ".join(parts)


def mult_str(r0, r1):
    if r0 == 1 and r1 == 1:
        return "(1,1)"
    return f"({r0},{r1})"


def main():
    con = sqlite3.connect(DB)
    cur = con.cursor()

    with OUT.open("w") as f:
        f.write("% Auto-generated from quivers.db by dump_no_veneziano_tables.py\n")
        f.write("% Requires: \\usepackage{booktabs,longtable,amssymb,amsmath}\n")
        f.write("\\section*{Non-Veneziano universality classes}\n\n")

        # Include classes that are not entirely Veneziano (veneziano_all=0).
        # For mixed classes (veneziano_any=1, veneziano_all=0), count only
        # non-Veneziano theories.
        n_nv = cur.execute(
            "SELECT COUNT(*) FROM universality_class WHERE veneziano_all=0"
        ).fetchone()[0]
        n_nv_scft = cur.execute(
            "SELECT COUNT(*) FROM universality_class WHERE veneziano_all=0 AND a_exact IS NOT NULL"
        ).fetchone()[0]
        n_nv_bw = n_nv - n_nv_scft
        n_th = cur.execute(
            "SELECT COUNT(*) FROM theory WHERE veneziano=0"
        ).fetchone()[0]
        f.write(
            f"{n_nv} non-Veneziano universality classes "
            f"({n_nv_scft} with IR fixed points, {n_nv_bw} below the conformal window), "
            f"covering {n_th} theories.\n\n"
        )

        for pair in GAUGE_PAIRS:
            rank_pairs = cur.execute(
                "SELECT DISTINCT rank0_mult, rank1_mult FROM universality_class "
                "WHERE gauge_pair=? AND veneziano_all=0 "
                "ORDER BY rank0_mult, rank1_mult",
                (pair,),
            ).fetchall()
            if not rank_pairs:
                continue

            f.write(f"\\subsection*{{{pair}}}\n\n")

            for r0, r1 in rank_pairs:
                classes = cur.execute(
                    """SELECT uc.class_id, uc.a_over_N2,
                              (SELECT COUNT(*) FROM theory t2
                               WHERE t2.class_id=uc.class_id AND t2.veneziano=0) AS n_nonven,
                              uc.a_exact, uc.c_exact, uc.R_exact, uc.a_over_c,
                              t.gauge0, t.gauge1, t.matter0, t.matter1, t.edges,
                              t.N_bif, t.N_fund_0, t.N_fund_1, t.N_rank2_0, t.N_rank2_1
                       FROM universality_class uc
                       JOIN theory t ON t.theory_id = uc.rep_theory_id
                       WHERE uc.gauge_pair=? AND uc.rank0_mult=? AND uc.rank1_mult=?
                             AND uc.veneziano_all=0
                       ORDER BY uc.a_over_N2 NULLS LAST""",
                    (pair, r0, r1),
                ).fetchall()
                if not classes:
                    continue

                scft_classes = [r for r in classes if r[3] is not None]
                bw_classes = [r for r in classes if r[3] is None]

                if scft_classes:
                    f.write(
                        f"\\subsubsection*{{{pair}, $(m_1,m_2)={mult_str(r0,r1)}$: "
                        f"IR fixed points}}\n"
                    )
                    f.write("\\begin{longtable}{rrrrrll}\n")
                    f.write("\\toprule\n")
                    f.write(
                        "ID & $\\#$ & $N^{(1)}_{r2}$ & $N^{(2)}_{r2}$ & $N_{\\text{bif}}$ "
                        "& $a/N^2 = c/N^2$ & $R$-charges \\\\\n"
                    )
                    f.write("\\midrule\n\\endhead\n")
                    for r in scft_classes:
                        (cid, a_num, nth, a_ex, c_ex, R_ex, aoc,
                         g0, g1, m0, m1, edges, nbif, nf0, nf1, nr0, nr1) = r
                        f.write(
                            f"{cid} & {nth} & {nr0:g} & {nr1:g} & {nbif} "
                            f"& {texify_expr(a_ex)} "
                            f"& {texify_R(R_ex)} \\\\\n"
                        )
                    f.write("\\bottomrule\n\\end{longtable}\n\n")

                if bw_classes:
                    f.write(
                        f"\\subsubsection*{{{pair}, $(m_1,m_2)={mult_str(r0,r1)}$: "
                        f"below conformal window}}\n"
                    )
                    f.write("\\begin{longtable}{rrrrrr}\n")
                    f.write("\\toprule\n")
                    f.write(
                        "ID & $\\#$ & $N^{(1)}_{r2}$ & $N^{(2)}_{r2}$ & $N_{\\text{bif}}$ "
                        "& $N_f$ \\\\\n"
                    )
                    f.write("\\midrule\n\\endhead\n")
                    for r in bw_classes:
                        (cid, a_num, nth, a_ex, c_ex, R_ex, aoc,
                         g0, g1, m0, m1, edges, nbif, nf0, nf1, nr0, nr1) = r
                        morphs = cur.execute(
                            "SELECT mc.N_rank2_0, mc.N_rank2_1, mc.N_bif, "
                            "mc.N_fund_0, mc.N_fund_1, COUNT(*) "
                            "FROM theory t "
                            "JOIN morphology_class mc ON t.morph_id = mc.morph_id "
                            "WHERE t.class_id = ? AND t.veneziano = 0 "
                            "GROUP BY mc.morph_id "
                            "ORDER BY mc.N_rank2_0, mc.N_rank2_1, mc.N_bif, mc.N_fund_0, mc.N_fund_1",
                            (cid,),
                        ).fetchall()
                        first = True
                        for mr0, mr1, mbif, mf0, mf1, mth in morphs:
                            id_col = str(cid) if first else ""
                            nth_col = str(nth) if first else ""
                            nf_total = mf0 + mf1
                            f.write(
                                f"{id_col} & {nth_col} & {mr0:g} & {mr1:g} & {mbif} "
                                f"& {nf_total} \\\\\n"
                            )
                            first = False
                        if len(morphs) > 1:
                            f.write("\\addlinespace\n")
                    f.write("\\bottomrule\n\\end{longtable}\n\n")

    con.close()
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
