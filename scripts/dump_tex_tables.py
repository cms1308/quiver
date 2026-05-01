"""Dump universality classes and morphology classes from quivers.db into TeX tables, grouped by gauge pair."""

import re
import sqlite3
from pathlib import Path

DB = Path(__file__).parent.parent / "quivers.db"
OUT = Path(__file__).parent.parent / "paper/sections/generated/classification_tables.tex"

GAUGE_PAIRS = ["SU-SU", "SU-SO", "SU-Sp", "SO-SO", "SO-Sp", "Sp-Sp"]


def texify_expr(s: str) -> str:
    if s is None or s == "":
        return ""
    s = s.replace("*", "")
    s = re.sub(r"√\(([^()]*)\)", r"\\sqrt{\1}", s)
    return f"${s}$"


def texify_R(s: str) -> str:
    if not s:
        return ""
    s = s.replace("*", "")
    s = re.sub(r"√\(([^()]*)\)", r"\\sqrt{\1}", s)
    s = s.replace("□̄", r"\bar{\square}")
    s = s.replace("□", r"\square")
    s = s.replace("Ā", r"\bar{A}")
    s = s.replace("S̄", r"\bar{S}")
    parts = [p.strip() for p in s.split(",")]
    parts = [re.sub(r"R_([^=]+)=", lambda m: f"R_{{{m.group(1).strip()}}}=", p) for p in parts]
    return "$" + ",\\;".join(parts) + "$"


def mult_str(r0, r1):
    if r0 == 1 and r1 == 1:
        return "(1,1)"
    return f"({r0},{r1})"


def write_class_table(f, pair, r0, r1, rows):
    f.write(f"\\subsubsection*{{Universality classes: {pair}, $(m_1,m_2)={mult_str(r0,r1)}$}}\n")
    f.write("\\begin{longtable}{rrrll}\n")
    f.write("\\toprule\n")
    f.write("ID & $\\#$ & $a/N^2$ & $a_{\\text{exact}}$ & $c_{\\text{exact}}$ \\\\\n")
    f.write("\\midrule\n\\endhead\n")
    for r in rows:
        cid, _, _, _, a_num, n, _, a_ex, c_ex, _, _, _, _ = r
        a_str = f"{a_num:.5f}" if a_num is not None else "---"
        f.write(f"{cid} & {n} & {a_str} & {texify_expr(a_ex)} & {texify_expr(c_ex)} \\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n\n")


def write_morph_table(f, pair, r0, r1, rows):
    f.write(f"\\subsubsection*{{Morphology classes: {pair}, $(m_1,m_2)={mult_str(r0,r1)}$}}\n")
    f.write("\\begin{longtable}{rrrrrrrr}\n")
    f.write("\\toprule\n")
    f.write("ID & $N^{(1)}_{r2}$ & $N^{(2)}_{r2}$ & $N_{\\text{bif}}$ & $N^{(1)}_f$ & $N^{(2)}_f$ & classes & theories \\\\\n")
    f.write("\\midrule\n\\endhead\n")
    for r in rows:
        mid, _, _, _, n20, n21, nb, nf0, nf1, nth, ncl = r
        f.write(f"{mid} & {n20:g} & {n21:g} & {nb} & {nf0} & {nf1} & {ncl} & {nth} \\\\\n")
    f.write("\\bottomrule\n\\end{longtable}\n\n")


def main():
    con = sqlite3.connect(DB)
    cur = con.cursor()

    with OUT.open("w") as f:
        f.write("% Auto-generated from quivers.db by dump_tex_tables.py\n")
        f.write("% Requires: \\usepackage{booktabs,longtable,amssymb,amsmath}\n")
        f.write("\\section*{Two-node quiver classification tables}\n\n")

        n_cls_total = cur.execute("SELECT COUNT(*) FROM universality_class").fetchone()[0]
        n_morph_total = cur.execute("SELECT COUNT(*) FROM morphology_class").fetchone()[0]
        n_th_total = cur.execute("SELECT COUNT(*) FROM theory").fetchone()[0]
        f.write(f"Database summary: {n_th_total} theories, {n_cls_total} universality classes, {n_morph_total} morphology classes.\n\n")

        for pair in GAUGE_PAIRS:
            f.write(f"\\section*{{{pair}}}\n\n")
            rank_pairs = cur.execute(
                "SELECT DISTINCT rank0_mult, rank1_mult FROM theory "
                "WHERE gauge_pair=? ORDER BY rank0_mult, rank1_mult", (pair,)
            ).fetchall()
            for r0, r1 in rank_pairs:
                morphs = cur.execute(
                    "SELECT morph_id, gauge_pair, rank0_mult, rank1_mult, N_rank2_0, N_rank2_1, "
                    "N_bif, N_fund_0, N_fund_1, n_theories, n_classes "
                    "FROM morphology_class WHERE gauge_pair=? AND rank0_mult=? AND rank1_mult=? "
                    "ORDER BY morph_id", (pair, r0, r1)
                ).fetchall()
                classes = cur.execute(
                    "SELECT class_id, gauge_pair, rank0_mult, rank1_mult, a_over_N2, n_theories, "
                    "rep_theory_id, a_exact, c_exact, R_exact, a_over_c, veneziano_any, veneziano_all "
                    "FROM universality_class WHERE gauge_pair=? AND rank0_mult=? AND rank1_mult=? "
                    "ORDER BY a_over_N2", (pair, r0, r1)
                ).fetchall()
                if not morphs and not classes:
                    continue
                write_morph_table(f, pair, r0, r1, morphs)
                write_class_table(f, pair, r0, r1, classes)

    con.close()
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
