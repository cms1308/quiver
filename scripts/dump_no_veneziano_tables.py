"""Dump non-Veneziano universality classes from quivers.db into TeX tables, grouped by gauge pair."""

import re
import sqlite3
import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from two_node_db import _node_types

DB = Path(__file__).parent.parent / "quivers.db"
OUT = Path(__file__).parent.parent / "paper/sections/generated/no_veneziano_tables.tex"

GAUGE_PAIRS = ["SU-SU", "SU-SO", "SU-Sp", "SO-SO", "SO-Sp", "Sp-Sp"]


def _type_tex(t: Fraction) -> str:
    if t.denominator == 1:
        return f"${t.numerator}$"
    sign = "-" if t.numerator < 0 else ""
    return f"${sign}\\frac{{{abs(t.numerator)}}}{{{t.denominator}}}$"


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


_RANK2_REPS = {"A", "Ā", "S", "S̄", "adj"}


def _parse_matter_reps(matter_str: str):
    if not matter_str or matter_str == "—":
        return []
    out = []
    for part in matter_str.split(" + "):
        part = part.strip()
        m = re.match(r"^(\d+)?(.+)", part)
        if m:
            cnt = int(m.group(1)) if m.group(1) else 1
            out.append((cnt, m.group(2).strip()))
    return out


def _val_to_tex(val):
    if val is None:
        return "?"
    f = Fraction(val).limit_denominator(1000)
    if abs(float(f) - val) < 1e-4:
        if f.denominator == 1:
            return str(f.numerator)
        s = f"\\frac{{{abs(f.numerator)}}}{{{f.denominator}}}"
        return ("-" + s) if f.numerator < 0 else s
    return f"{val:.4f}"


def _R_exact_kinds(R_exact):
    """Parse R_exact → {'r2_1': sym, 'r2_2': sym, 'bif': sym} (first occurrence only)."""
    kinds = {}
    if not R_exact:
        return kinds
    for e in R_exact.split(","):
        e = e.strip()
        m = re.match(r"R_(\S+?)_0\s*=\s*(.+)", e)
        if m and m.group(1) in _RANK2_REPS and "r2_1" not in kinds:
            kinds["r2_1"] = m.group(2).strip()
            continue
        m = re.match(r"R_(\S+?)_1\s*=\s*(.+)", e)
        if m and m.group(1) in _RANK2_REPS and "r2_2" not in kinds:
            kinds["r2_2"] = m.group(2).strip()
            continue
        m = re.match(r"R_bif\s*=\s*(.+)", e)
        if m and "bif" not in kinds:
            kinds["bif"] = m.group(1).strip()
    return kinds


def _R_numerical_kinds(R_num, matter0, matter1):
    """Parse R_numerical → (r2_1_val, r2_2_val, bif_val) as floats (or None)."""
    if not R_num:
        return (None, None, None)
    parsed = []
    for e in R_num.split(","):
        m = re.match(r"\s*R_(\S+?)\s*=\s*([-0-9.eE+]+)", e)
        if m:
            parsed.append((m.group(1), float(m.group(2))))
    reps0 = _parse_matter_reps(matter0)
    reps1 = _parse_matter_reps(matter1)
    n0, n1 = len(reps0), len(reps1)
    node0 = parsed[:n0]
    node1 = parsed[n0:n0 + n1]
    edges = parsed[n0 + n1:]
    r2_1 = next((v for nm, v in node0 if nm in _RANK2_REPS), None)
    r2_2 = next((v for nm, v in node1 if nm in _RANK2_REPS), None)
    bif = edges[0][1] if edges else None
    return (r2_1, r2_2, bif)


def _fmt_sym(s):
    if s is None:
        return "?"
    s = s.replace("*", "")
    s = re.sub(r"√\(([^()]*)\)", r"\\sqrt{\1}", s)
    s = s.replace("Ā", r"\bar{A}")
    s = s.replace("S̄", r"\bar{S}")
    s = _fractionize(s)
    return s


def class_R_tex(R_exact, R_num, matter0, matter1,
                has_r2_1, has_r2_2, has_bif):
    """Return TeX array with up to three rows: R_{r2,1}, R_{r2,2}, R_{bif}."""
    kinds = _R_exact_kinds(R_exact)
    need_fallback = ((has_r2_1 and "r2_1" not in kinds)
                     or (has_r2_2 and "r2_2" not in kinds)
                     or (has_bif and "bif" not in kinds))
    if need_fallback:
        num_vals = _R_numerical_kinds(R_num, matter0, matter1)
        if has_r2_1 and "r2_1" not in kinds and num_vals[0] is not None:
            kinds["r2_1"] = _val_to_tex(num_vals[0])
        if has_r2_2 and "r2_2" not in kinds and num_vals[1] is not None:
            kinds["r2_2"] = _val_to_tex(num_vals[1])
        if has_bif and "bif" not in kinds and num_vals[2] is not None:
            kinds["bif"] = _val_to_tex(num_vals[2])

    rows = []
    if has_r2_1:
        rows.append(f"R_{{r2,1}}={_fmt_sym(kinds.get('r2_1'))}")
    if has_r2_2:
        rows.append(f"R_{{r2,2}}={_fmt_sym(kinds.get('r2_2'))}")
    if has_bif:
        rows.append(f"R_{{bif}}={_fmt_sym(kinds.get('bif'))}")
    if not rows:
        return "---"
    return "$\\begin{array}{l}" + " \\\\ ".join(rows) + "\\end{array}$"


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

    # Now apply TeX transformations, dedup identical entries (e.g. R_bif=X listed twice)
    tex_parts = []
    seen = set()
    for p in output_parts:
        p = re.sub(r"√\(([^()]*)\)", r"\\sqrt{\1}", p)
        p = p.replace("Ā", r"\bar{A}")
        p = p.replace("S̄", r"\bar{S}")
        p = _fractionize(p)
        if p in seen:
            continue
        seen.add(p)
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


def gauge_label(pair, r0, r1):
    g0, g1 = pair.split("-")
    def _grp(g, m):
        rank = f"{m}N" if m > 1 else "N"
        return f"{g}({rank})"
    return f"{_grp(g0, r0)}$\\times${_grp(g1, r1)}"


def main():
    con = sqlite3.connect(DB)
    cur = con.cursor()

    with OUT.open("w") as f:
        f.write("% Auto-generated from quivers.db by dump_no_veneziano_tables.py\n")
        f.write("% Requires: \\usepackage{booktabs,longtable,amssymb,amsmath,multirow}\n")
        f.write("\\renewcommand{\\arraystretch}{1.5}\n")
        f.write("\\section*{Non-Veneziano universality classes}\n\n")

        n_nv = cur.execute(
            "SELECT COUNT(*) FROM universality_class WHERE veneziano_all=0 AND a_exact IS NOT NULL"
        ).fetchone()[0]
        n_th = cur.execute(
            """SELECT COUNT(*) FROM theory t
               JOIN universality_class uc ON t.class_id=uc.class_id
               WHERE uc.veneziano_all=0 AND uc.a_exact IS NOT NULL AND t.veneziano=0"""
        ).fetchone()[0]
        f.write(
            f"{n_nv} non-Veneziano universality classes with IR fixed points, "
            f"covering {n_th} theories.\n\n"
        )

        # Summary table: classes and theories per gauge pair × rank combination
        summary = cur.execute(
            """SELECT uc.gauge_pair, uc.rank0_mult, uc.rank1_mult,
                      COUNT(DISTINCT uc.class_id),
                      SUM(CASE WHEN t.veneziano=0 THEN 1 ELSE 0 END)
               FROM universality_class uc
               JOIN theory t ON t.class_id=uc.class_id
               WHERE uc.veneziano_all=0 AND uc.a_exact IS NOT NULL
               GROUP BY uc.gauge_pair, uc.rank0_mult, uc.rank1_mult
               ORDER BY uc.gauge_pair, uc.rank0_mult, uc.rank1_mult"""
        ).fetchall()

        # Group by gauge pair for multirow
        from itertools import groupby
        pair_groups = {}
        for row in summary:
            pair_groups.setdefault(row[0], []).append(row)

        f.write("\\begin{longtable}{l|l|r|r}\n")
        f.write("\\toprule\n")
        f.write("Gauge pair & Rank structure & Classes & Theories \\\\\n")
        f.write("\\midrule\n\\endhead\n")
        for pair in GAUGE_PAIRS:
            rows = pair_groups.get(pair, [])
            if not rows:
                continue
            n_rows = len(rows)
            for i, (gp, r0, r1, nc, nth_s) in enumerate(rows):
                pair_cell = (
                    f"\\multirow{{{n_rows}}}{{*}}{{{gp}}}" if i == 0 else ""
                )
                f.write(
                    f"{pair_cell} & {gauge_label(gp, r0, r1)} & {nc} & {nth_s} \\\\\n"
                )
            if pair != [p for p in GAUGE_PAIRS if p in pair_groups][-1]:
                f.write("\\hline\n")
        f.write("\\bottomrule\n\\end{longtable}\n\n")

        for pair in GAUGE_PAIRS:
            rank_pairs = cur.execute(
                "SELECT DISTINCT rank0_mult, rank1_mult FROM universality_class "
                "WHERE gauge_pair=? AND veneziano_all=0 AND a_exact IS NOT NULL "
                "ORDER BY rank0_mult, rank1_mult",
                (pair,),
            ).fetchall()
            if not rank_pairs:
                continue

            pair_total_classes = cur.execute(
                "SELECT COUNT(*) FROM universality_class "
                "WHERE gauge_pair=? AND veneziano_all=0 AND a_exact IS NOT NULL",
                (pair,),
            ).fetchone()[0]
            pair_total_theories = cur.execute(
                """SELECT COUNT(*) FROM theory t
                   JOIN universality_class uc ON t.class_id=uc.class_id
                   WHERE uc.gauge_pair=? AND uc.veneziano_all=0
                         AND uc.a_exact IS NOT NULL AND t.veneziano=0""",
                (pair,),
            ).fetchone()[0]
            f.write(
                f"\\subsection*{{{pair} "
                f"({pair_total_classes} classes, {pair_total_theories} theories)}}\n\n"
            )

            for r0, r1 in rank_pairs:
                classes = cur.execute(
                    """SELECT uc.class_id, uc.a_over_N2,
                              (SELECT COUNT(*) FROM theory t2
                               WHERE t2.class_id=uc.class_id AND t2.veneziano=0) AS n_nonven,
                              uc.a_exact, uc.c_exact, uc.R_exact, uc.a_over_c,
                              t.gauge0, t.gauge1, t.matter0, t.matter1, t.edges,
                              t.N_bif, t.N_fund_0, t.N_fund_1, t.N_rank2_0, t.N_rank2_1,
                              t.R_numerical
                       FROM universality_class uc
                       JOIN theory t ON t.theory_id = uc.rep_theory_id
                       WHERE uc.gauge_pair=? AND uc.rank0_mult=? AND uc.rank1_mult=?
                             AND uc.veneziano_all=0 AND uc.a_exact IS NOT NULL
                       ORDER BY uc.a_over_N2""",
                    (pair, r0, r1),
                ).fetchall()
                if not classes:
                    continue

                class_has = {}
                for r in classes:
                    cid = r[0]
                    row = cur.execute(
                        "SELECT MAX(N_rank2_0), MAX(N_rank2_1), MAX(N_bif) "
                        "FROM theory WHERE class_id=? AND veneziano=0",
                        (cid,),
                    ).fetchone()
                    class_has[cid] = (
                        (row[0] or 0) > 0,
                        (row[1] or 0) > 0,
                        (row[2] or 0) > 0,
                    )

                f.write(f"\\subsubsection*{{{gauge_label(pair, r0, r1)}}}\n")
                f.write("\\begin{longtable}{cc|ccc|l|l|r}\n")
                f.write("\\toprule\n")
                f.write(
                    "Type$^{(1)}$ & Type$^{(2)}$ "
                    "& $N^{(1)}_{r2}$ & $N^{(2)}_{r2}$ & $N_{\\text{bif}}$ "
                    "& $a/N^2 = c/N^2$ & $R$-charges & $\\#$ \\\\\n"
                )
                f.write("\\midrule\n\\endhead\n")
                for i, r in enumerate(classes):
                    (cid, a_num, nth, a_ex, c_ex, R_ex, aoc,
                     g0, g1, m0, m1, edges, nbif, nf0, nf1, nr0, nr1,
                     R_num) = r
                    morphs = cur.execute(
                        "SELECT mc.N_rank2_0, mc.N_rank2_1, mc.N_bif, COUNT(*) "
                        "FROM theory t "
                        "JOIN morphology_class mc ON t.morph_id = mc.morph_id "
                        "WHERE t.class_id = ? AND t.veneziano = 0 "
                        "GROUP BY mc.morph_id "
                        "ORDER BY mc.N_rank2_0, mc.N_rank2_1, mc.N_bif",
                        (cid,),
                    ).fetchall()
                    if r0 == r1:
                        canon = {}
                        for mr0, mr1, mnb, mc_cnt in morphs:
                            key = (min(mr0, mr1), max(mr0, mr1), mnb)
                            canon[key] = canon.get(key, 0) + mc_cnt
                        morphs = [(k[0], k[1], k[2], canon[k]) for k in sorted(canon)]
                    typed = [(_node_types(pair, r0, r1, mr0, mr1, mnb), mr0, mr1, mnb)
                             for mr0, mr1, mnb, _ in morphs]
                    n_sub = len(typed)
                    uniform = len({ts for ts, *_ in typed}) == 1
                    for j, ((t1, t2), mr0, mr1, mnb) in enumerate(typed):
                        if uniform and n_sub > 1:
                            t1c = (f"\\multirow{{{n_sub}}}{{*}}{{{_type_tex(t1)}}}"
                                   if j == 0 else "")
                            t2c = (f"\\multirow{{{n_sub}}}{{*}}{{{_type_tex(t2)}}}"
                                   if j == 0 else "")
                        else:
                            t1c = _type_tex(t1)
                            t2c = _type_tex(t2)
                        first = (j == 0)
                        if first:
                            h0, h1, hb = class_has.get(cid, (True, True, True))
                            a_raw = texify_expr(a_ex)
                            R_raw = class_R_tex(R_ex, R_num, m0, m1, h0, h1, hb)
                            nth_raw = str(nth)
                            if n_sub > 1:
                                a_col = f"\\multirow{{{n_sub}}}{{*}}{{{a_raw}}}"
                                R_col = f"\\multirow{{{n_sub}}}{{*}}{{{R_raw}}}"
                                nth_col = f"\\multirow{{{n_sub}}}{{*}}{{{nth_raw}}}"
                            else:
                                a_col, R_col, nth_col = a_raw, R_raw, nth_raw
                        else:
                            a_col = R_col = nth_col = ""
                        f.write(
                            f"{t1c} & {t2c} "
                            f"& {mr0:g} & {mr1:g} & {mnb} "
                            f"& {a_col} & {R_col} & {nth_col} \\\\\n"
                        )
                    if i < len(classes) - 1:
                        f.write("\\hline\n")
                f.write("\\bottomrule\n\\end{longtable}\n\n")

    con.close()
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
