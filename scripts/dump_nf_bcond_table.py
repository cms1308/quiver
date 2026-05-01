"""
Generate nf_bcond_table.tex: paper section on the N_f conformal window
for SU(N)×SU(N) non-Veneziano classes.

Notation:
  O = (g_1=0, g_2=0)  : origin / free UV fixed point
  A = (g_1=g_1*, g_2=0): node-1 fixed point, node-2 free
  B = (g_1=0, g_2=g_2*): node-1 free, node-2 fixed point
"""

import re
import sqlite3
from fractions import Fraction
from pathlib import Path

DB  = Path(__file__).parent.parent / "quivers.db"
OUT = Path(__file__).parent.parent / "paper/sections/generated/nf_bcond_table.tex"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fractionize(s):
    def _repl(m):
        sign = m.group(1) or ""
        return f"{sign}\\frac{{{m.group(2)}}}{{{m.group(3)}}}"
    return re.sub(r"([-+]\s*)?([0-9]+)\s*/\s*([0-9]+)", _repl, s)

def texify_a(s):
    if s is None:
        return "---"
    s = s.replace("*", "")
    s = re.sub(r"√\(([^()]*)\)", r"\\sqrt{\1}", s)
    s = _fractionize(s)
    return f"${s}$"

def _parse_alpha(bound_str):
    """Extract coefficient α of N from a bound string."""
    if not bound_str:
        return "?"
    if "conformal" in bound_str:
        return "0"
    m = re.match(r"N_f\s*[≤=]\s*(-?[0-9./]+)\s*\*\s*N", bound_str)
    if m:
        return m.group(1).strip()
    m2 = re.match(r"N_f\s*[≤=]\s*-?[0-9]+", bound_str)
    if m2:
        return "0"
    return "?"

def _alpha_tex(a):
    """Format α as exact fraction (rational) or 4-digit decimal (irrational)."""
    a = a.strip().lstrip("+")
    if not a or a == "?":
        return a
    if "." in a:
        try:
            v = float(a)
            return "0" if abs(v) < 1e-4 else f"{v:.4f}"
        except Exception:
            return a
    try:
        f = Fraction(a)
        if abs(f) < Fraction(1, 10000):
            return "0"
        if f.denominator == 1:
            return str(f.numerator)
        sign = "-" if f.numerator < 0 else ""
        return f"{sign}\\frac{{{abs(f.numerator)}}}{{{f.denominator}}}"
    except Exception:
        return a

def _morph_tex(nr0, nr1, nb):
    def _h(v):
        f = Fraction(v).limit_denominator(10)
        return str(f.numerator) if f.denominator == 1 \
               else f"\\tfrac{{{f.numerator}}}{{{f.denominator}}}"
    return (f"$({_h(nr0)},\\,{_h(nr1)},\\,{int(nb)})$")

def _morph_full(nr0, nr1, nb):
    return (f"$(N_{{r2}}^{{(1)}}, N_{{r2}}^{{(2)}}, N_{{\\rm bif}})"
            f" = {_morph_tex(nr0, nr1, nb)}$")


# ---------------------------------------------------------------------------
# Physics prose for each (class_id, nr0, nr1, nb) morphology
# ---------------------------------------------------------------------------

MORPH_PROSE = {

  # ---- Class 1 ----
  (1, 0.5, 0.5, 1): (
    r"With one rank-2 tensor on each node and a single bifundamental, "
    r"the $\mathcal{O}$-point beta function gives $\alpha_{\mathcal{O}} = 2$ on both nodes. "
    r"At $a/N^2 = 0$ all matter R-charges vanish at the fixed point, "
    r"so the bifundamental enters the $\mathcal{B}$-point (and $\mathcal{A}$-point) "
    r"beta function with R-charge zero, exactly as at $\mathcal{O}$. "
    r"The difference between the two bounds arises solely from the "
    r"rank-2 tensor contribution: at $\mathcal{B}$ the rank-2 tensor on node 1 "
    r"is free ($R = 0$) while the bifundamental contribution from node 2 "
    r"is unchanged, halving the net coefficient. "
    r"The result is a uniform compression $\alpha_{\mathcal{O}} = 2 \to \alpha_{\mathcal{B}} = 1$ "
    r"on both nodes."
  ),

  # ---- Class 2 ----
  (2, 0.5, 1.5, 1): (
    r"Node 1 carries one rank-2 tensor ($\alpha_{\mathcal{O}}^{(1)} = 2$) "
    r"and node 2 carries three ($\alpha_{\mathcal{O}}^{(2)} = 1$). "
    r"At the interacting fixed point the bifundamental acquires R-charge "
    r"$R_{\rm bif} \approx 0.146 < \tfrac{2}{3}$, giving a negative anomalous "
    r"dimension $\gamma_{\rm bif} = 3R_{\rm bif} - 2 \approx -1.562$. "
    r"This increased depletion of the one-loop coefficient tightens both "
    r"boundary conditions by the same amount "
    r"$\Delta\alpha = \tfrac{1}{2}|\gamma_{\rm bif}| \approx 0.781$: "
    r"$\alpha_{\mathcal{B}}^{(1)} \approx 1.219 < \alpha_{\mathcal{O}}^{(1)} = 2$ "
    r"and $\alpha_{\mathcal{A}}^{(2)} \approx 0.219 < \alpha_{\mathcal{O}}^{(2)} = 1$. "
    r"The binding constraint on the conformal window is $\alpha_{\mathcal{A}}^{(2)} \approx 0.219$."
  ),
  (2, 1.5, 0.5, 1): (
    r"The node-exchange mirror of the previous morphology. "
    r"The same $\Delta\alpha \approx 0.781$ tightening applies: "
    r"$\alpha_{\mathcal{B}}^{(1)} \approx 0.219 < \alpha_{\mathcal{O}}^{(1)} = 1$ "
    r"and $\alpha_{\mathcal{A}}^{(2)} \approx 1.219 < \alpha_{\mathcal{O}}^{(2)} = 2$. "
    r"The binding constraint is now $\alpha_{\mathcal{B}}^{(1)} \approx 0.219$."
  ),

  # ---- Class 4 ----
  (4, 2.5, 0.5, 1): (
    r"Node 1 carries five half-units of rank-2 matter "
    r"($N_{r2}^{(1)} = \tfrac{5}{2}$), saturating the $\mathcal{O}$-point "
    r"beta function: $\alpha_{\mathcal{O}}^{(1)} = 0$. "
    r"The bifundamental R-charge at the $\sqrt{13}$ fixed point is "
    r"$R_{\rm bif} \approx 0.162$, giving $\Delta\alpha \approx 0.757$. "
    r"Since $\alpha_{\mathcal{O}}^{(1)} = 0$, the tightening pushes node 1's "
    r"coefficient below zero: $\alpha_{\mathcal{B}}^{(1)} \approx -0.757 < 0$, "
    r"so $b_0^{(1)}|_{\mathcal{B}} < 0$ is automatically satisfied for any "
    r"$N_f \geq 0$ at large $N$. "
    r"Node 2, with one rank-2 tensor ($\alpha_{\mathcal{O}}^{(2)} = 2$), "
    r"is tightened to $\alpha_{\mathcal{A}}^{(2)} \approx 1.243$: "
    r"the conformal window is controlled by $b_0^{(2)}|_{\mathcal{A}}$."
  ),

  # ---- Class 9, 5 morphologies ----
  (9, 0.0, 0.0, 4): (
    r"Both nodes carry no rank-2 tensors and are connected by four bifundamentals. "
    r"At $\mathcal{O}$, each node has $b_0 = N$ giving $\alpha_{\mathcal{O}} = 1$. "
    r"At $\mathcal{B}$, the four bifundamentals carry their fixed-point R-charges "
    r"from node 2; since $N_{\rm bif} = 4$ saturates the formula "
    r"$\alpha_{\mathcal{B}} = 1 - N_{\rm bif}/4$, the coefficient collapses to zero. "
    r"The conformal window for additional fundamentals is thus "
    r"$N$-independent at both boundaries, anticipating the behaviour of class~44."
  ),
  (9, 0.0, 1.0, 4): (
    r"Node 1 has no rank-2 tensors ($\alpha_{\mathcal{O}}^{(1)} = 1$) while "
    r"node 2 carries one adjoint unit ($N_{r2}^{(2)} = 1$, "
    r"$\alpha_{\mathcal{O}}^{(2)} = 0$, UV threshold). "
    r"With $\Delta\alpha = N_{\rm bif}/4 = 1$: "
    r"node 1's coefficient drops to $\alpha_{\mathcal{B}}^{(1)} = 0$ (conformal threshold at $\mathcal{B}$), "
    r"and node 2's UV-saturated condition is pushed to "
    r"$\alpha_{\mathcal{A}}^{(2)} = -1 < 0$, automatically satisfied. "
    r"No fundamental deformation with $N_f \sim N$ is possible for either node."
  ),
  (9, 0.5, 0.5, 3): (
    r"One rank-2 tensor on each node and three bifundamentals. "
    r"The formula $\alpha_{\mathcal{B}} = 1 - N_{\rm bif}/4 = \tfrac{1}{4}$ applies "
    r"to both nodes by symmetry. "
    r"The conformal window is compressed to one-quarter of the $\mathcal{O}$-point estimate."
  ),
  (9, 1.0, 1.0, 2): (
    r"Two rank-2 equivalent units per node and two bifundamentals. "
    r"The formula gives $\alpha_{\mathcal{B}} = 1 - 2/4 = \tfrac{1}{2}$ on both nodes. "
    r"The conformal window is compressed to one-half of the $\mathcal{O}$-point estimate."
  ),
  (9, 1.5, 1.5, 1): (
    r"Three rank-2 equivalent units per node and one bifundamental. "
    r"The formula gives $\alpha_{\mathcal{B}} = 1 - 1/4 = \tfrac{3}{4}$ on both nodes: "
    r"the mildest compression within class~9, consistent with the single "
    r"bifundamental carrying only a small fixed-point R-charge."
  ),

  # ---- Class 16 ----
  (16, 1.5, 0.5, 3): (
    r"Node 1 carries $\tfrac{3}{2}$ rank-2 equivalent units with three bifundamentals "
    r"($\alpha_{\mathcal{O}}^{(1)} = 0$, at the $\mathcal{O}$-point threshold). "
    r"$\Delta\alpha = \tfrac{3}{2}|\gamma_{\rm bif}| \approx 0.573$: "
    r"$\alpha_{\mathcal{B}}^{(1)} \approx -0.573 < 0$ "
    r"(automatically satisfying $b_0^{(1)}|_{\mathcal{B}} < 0$ for large $N$), "
    r"and node 2's condition tightens to "
    r"$\alpha_{\mathcal{A}}^{(2)} \approx 0.427 < \alpha_{\mathcal{O}}^{(2)} = 1$. "
    r"The conformal window is entirely controlled by $b_0^{(2)}|_{\mathcal{A}}$."
  ),

  # ---- Class 21, 2 morphologies ----
  (21, 1.0, 2.0, 2): (
    r"Node 1 carries one rank-2 equivalent unit ($\alpha_{\mathcal{O}}^{(1)} = 1$) "
    r"and node 2 carries two ($\alpha_{\mathcal{O}}^{(2)} = 0$, at the threshold). "
    r"$\Delta\alpha = N_{\rm bif}|\gamma_{\rm bif}|/2 \approx 0.325$: "
    r"node 1's condition tightens to "
    r"$\alpha_{\mathcal{B}}^{(1)} \approx 0.675 < \alpha_{\mathcal{O}}^{(1)} = 1$, "
    r"and node 2's UV-saturated condition is pushed to "
    r"$\alpha_{\mathcal{A}}^{(2)} \approx -0.325 < 0$, automatically satisfied. "
    r"The conformal window is controlled by $b_0^{(1)}|_{\mathcal{B}}$."
  ),
  (21, 2.0, 1.0, 2): (
    r"The mirror of the previous morphology with nodes 1 and 2 exchanged. "
    r"Node 1 ($N_{r2}^{(1)} = 2$, $\alpha_{\mathcal{O}}^{(1)} = 0$) "
    r"is at the $\mathcal{O}$-point threshold: "
    r"$\alpha_{\mathcal{B}}^{(1)} \approx -0.325 < 0$, automatically satisfied. "
    r"Node 2's condition tightens to "
    r"$\alpha_{\mathcal{A}}^{(2)} \approx 0.675 < \alpha_{\mathcal{O}}^{(2)} = 1$: "
    r"the conformal window is controlled by $b_0^{(2)}|_{\mathcal{A}}$."
  ),

  # ---- Class 23, 2 morphologies ----
  (23, 1.5, 2.5, 1): (
    r"Node 1 carries $\tfrac{3}{2}$ rank-2 equivalent units ($\alpha_{\mathcal{O}}^{(1)} = 1$) "
    r"and node 2 carries $\tfrac{5}{2}$ ($\alpha_{\mathcal{O}}^{(2)} = 0$, at the threshold). "
    r"$\Delta\alpha = \tfrac{1}{2}|\gamma_{\rm bif}| \approx 0.147$: "
    r"node 1's condition tightens to "
    r"$\alpha_{\mathcal{B}}^{(1)} \approx 0.853 < \alpha_{\mathcal{O}}^{(1)} = 1$, "
    r"the largest tightened coefficient among all irrational classes, "
    r"and node 2's UV-saturated condition is pushed to "
    r"$\alpha_{\mathcal{A}}^{(2)} \approx -0.147 < 0$, automatically satisfied. "
    r"The conformal window is controlled by $b_0^{(1)}|_{\mathcal{B}}$."
  ),
  (23, 2.5, 1.5, 1): (
    r"The mirror morphology: node 1 carries $\tfrac{5}{2}$ rank-2 equivalent units "
    r"($\alpha_{\mathcal{O}}^{(1)} = 0$, at threshold): "
    r"$\alpha_{\mathcal{B}}^{(1)} \approx -0.147 < 0$, automatically satisfied. "
    r"Node 2 ($\tfrac{3}{2}$ rank-2 units, $\alpha_{\mathcal{O}}^{(2)} = 1$) sees "
    r"its condition tighten to "
    r"$\alpha_{\mathcal{A}}^{(2)} \approx 0.853$: "
    r"the conformal window is controlled by $b_0^{(2)}|_{\mathcal{A}}$."
  ),

  # ---- Class 44, 5 morphologies ----
  (44, 0.0, 0.0, 6): (
    r"Six bifundamentals, no rank-2 tensors on either node. "
    r"Both nodes are already at the $\mathcal{O}$-point conformal threshold "
    r"($\alpha_{\mathcal{O}} = 0$). "
    r"At both $\mathcal{A}$ and $\mathcal{B}$, $\alpha = 0$ as well: "
    r"no fundamental deformation is possible in the large-$N$ limit."
  ),
  (44, 1.0, 1.0, 4): (
    r"One rank-2 equivalent unit per node and four bifundamentals. "
    r"At $\mathcal{O}$, node 1 admits up to a constant number of fundamentals "
    r"(independent of $N$) while node 2 is at the threshold. "
    r"Both $\mathcal{B}$- and $\mathcal{A}$-point bounds have $\alpha = 0$: "
    r"the conformal window size is $\mathcal{O}(1)$ and $N$-independent."
  ),
  (44, 1.5, 1.5, 3): (
    r"Three rank-2 equivalent units per node with three bifundamentals. "
    r"Both nodes admit at most $\mathcal{O}(1)$ fundamentals at $\mathcal{O}$, "
    r"and $\alpha_{\mathcal{B}} = \alpha_{\mathcal{A}} = 0$: $N$-independence is confirmed."
  ),
  (44, 2.0, 2.0, 2): (
    r"Two rank-2 equivalent units per node and two bifundamentals. "
    r"At $\mathcal{O}$ some configurations of the rank-2 matter place node 2 "
    r"at the conformal threshold while node 1 can admit a constant number of fundamentals. "
    r"In all cases $\alpha_{\mathcal{B}} = \alpha_{\mathcal{A}} = 0$: "
    r"the conformal window is $N$-independent."
  ),
  (44, 2.5, 2.5, 1): (
    r"Five rank-2 equivalent units per node and one bifundamental. "
    r"Both nodes are near the $\mathcal{O}$-point threshold and admit "
    r"at most $\mathcal{O}(1)$ fundamentals. "
    r"The $\alpha = 0$ property at both $\mathcal{A}$ and $\mathcal{B}$ "
    r"confirms that this morphology, like all others in class~44, "
    r"lies beyond the reach of large-$N$ fundamental deformations."
  ),
}

CLASS_INTRO = {
  1: (
    r"This class has $a/N^2 = 0$, the lowest possible value, "
    r"corresponding to a fixed point at which all matter R-charges "
    r"vanish identically. "
    r"The $\mathcal{B}$-condition and the $\mathcal{O}$-condition therefore "
    r"differ only in the rank-2 tensor contribution, "
    r"yielding a simple uniform compression of the conformal window "
    r"by a factor of two."
  ),
  2: (
    r"The fixed point involves irrational R-charges proportional to $\sqrt{5}$. "
    r"The class has two morphologies related by node exchange. "
    r"In both, the same shift $\Delta\alpha \approx 0.781$ tightens both boundary "
    r"conditions relative to $\mathcal{O}$. "
    r"The binding constraint comes from the heavier node "
    r"(larger $N_{r2}$, smaller $\alpha_{\mathcal{O}}$), "
    r"whose tightened coefficient $\approx 0.219$ controls the conformal window."
  ),
  4: (
    r"This class has irrational R-charges proportional to $\sqrt{13}$ "
    r"and a single asymmetric morphology. "
    r"The heavier node saturates the $\mathcal{O}$-point beta function "
    r"($\alpha_{\mathcal{O}}^{(1)} = 0$): "
    r"the same tightening shift $\Delta\alpha \approx 0.757$ "
    r"pushes $\alpha_{\mathcal{B}}^{(1)} \approx -0.757 < 0$, "
    r"automatically satisfying $b_0^{(1)}|_{\mathcal{B}} < 0$ for large $N$. "
    r"The conformal window is controlled by the lighter node's "
    r"tightened coefficient $\alpha_{\mathcal{A}}^{(2)} \approx 1.243$."
  ),
  9: (
    r"This is the class with the simplest rational structure. "
    r"All morphologies share $a/N^2 = 27/64$ and have $\alpha_{\mathcal{O}} = 1$ "
    r"on both nodes. The $\mathcal{B}$- (and $\mathcal{A}$-) point coefficient "
    r"follows the exact rational formula\vspace{-2pt} "
    "\n"
    r"$$\alpha_{\mathcal{B}} = 1 - \frac{N_{\rm bif}}{4}$$"
    "\n"
    r"for the symmetric morphologies ($N_{r2}^{(1)} = N_{r2}^{(2)}$). "
    r"The formula counts how much the four bifundamental R-charges "
    r"--- each equal to $\tfrac{1}{2}$ at the $a/N^2 = 27/64$ fixed point "
    r"--- deplete the boundary beta function relative to the free-field value."
  ),
  16: (
    r"Irrational R-charges proportional to $\sqrt{5}$, "
    r"a single asymmetric morphology with three bifundamentals. "
    r"The structure is analogous to class~4: "
    r"the UV-saturated node's condition is automatically satisfied at $\mathcal{B}$, "
    r"and the lighter node's condition tightens to "
    r"$\alpha_{\mathcal{A}}^{(2)} \approx 0.427$ at $\mathcal{A}$."
  ),
  21: (
    r"Irrational R-charges proportional to $\sqrt{10}$, "
    r"two morphologies related by node exchange. "
    r"In both, the UV-saturated node's boundary condition is automatically "
    r"satisfied (pushed to $\alpha < 0$), "
    r"and the conformal window is controlled by the other node's tightened coefficient."
  ),
  23: (
    r"The highest-$a$ irrational class, with R-charges proportional to $\sqrt{17}$. "
    r"Two morphologies related by node exchange. "
    r"The tightened boundary coefficient $\alpha \approx 0.853$ is the "
    r"largest among all non-Veneziano irrational classes "
    r"(at $\mathcal{B}$ or $\mathcal{A}$ depending on the morphology), "
    r"approaching unity as $a/N^2 \to \tfrac{1}{2}$."
  ),
  44: (
    r"At the maximum value $a/N^2 = \tfrac{1}{2}$, "
    r"the $\mathcal{B}$- and $\mathcal{A}$-point coefficients vanish on "
    r"\emph{both} nodes for every morphology: $\alpha_{\mathcal{B}} = \alpha_{\mathcal{A}} = 0$. "
    r"The conformal window for additional fundamental flavours is "
    r"entirely $N$-independent, bounded only by the $\mathcal{O}(1)$ "
    r"constant $\beta$ that depends on the S/A character of the "
    r"rank-2 matter. "
    r"Class~44 is the endpoint of the formula $\alpha_{\mathcal{B}} = 1 - N_{\rm bif}/4$ "
    r"extrapolated to $N_{\rm bif} = 4$, generalised to all morphologies."
  ),
}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    con = sqlite3.connect(DB)
    cur = con.cursor()

    # Fetch α coefficients per (class, morphology) — representative row per morph
    morph_rows = cur.execute(
        """SELECT t.class_id, uc.a_exact, uc.a_over_N2,
                  t.N_rank2_0, t.N_rank2_1, t.N_bif,
                  t.nf_bound0, t.nf_bound1,
                  t.B_cond_A, t.B_cond_B
           FROM theory t
           JOIN universality_class uc ON t.class_id = uc.class_id
           WHERE t.gauge_pair = 'SU-SU' AND t.veneziano = 0
                 AND uc.veneziano_all = 0 AND uc.a_exact IS NOT NULL
                 AND uc.rank0_mult = 1 AND uc.rank1_mult = 1
                 AND t.delta0 = -4 AND (t.delta1 = -4 OR t.delta1 = 0)
           GROUP BY t.class_id, t.N_rank2_0, t.N_rank2_1, t.N_bif
           ORDER BY uc.a_over_N2,
                    t.N_rank2_0, t.N_rank2_1, t.N_bif"""
    ).fetchall()

    # Also need class 9 (0,0,4) and (0,1,4) which have delta=0
    # Already covered by delta0=-4 OR delta1=0 condition above,
    # but let's also explicitly get them
    extra = cur.execute(
        """SELECT t.class_id, uc.a_exact, uc.a_over_N2,
                  t.N_rank2_0, t.N_rank2_1, t.N_bif,
                  t.nf_bound0, t.nf_bound1,
                  t.B_cond_A, t.B_cond_B
           FROM theory t
           JOIN universality_class uc ON t.class_id = uc.class_id
           WHERE t.gauge_pair = 'SU-SU' AND t.veneziano = 0
                 AND uc.veneziano_all = 0 AND uc.a_exact IS NOT NULL
                 AND uc.rank0_mult = 1 AND uc.rank1_mult = 1
                 AND t.delta0 = 0 AND t.delta1 = 0
           GROUP BY t.class_id, t.N_rank2_0, t.N_rank2_1, t.N_bif
           ORDER BY uc.a_over_N2,
                    t.N_rank2_0, t.N_rank2_1, t.N_bif"""
    ).fetchall()

    # Merge, deduplicate by (class_id, nr0, nr1, nb)
    seen = set()
    all_rows = []
    for r in list(morph_rows) + list(extra):
        key = (r[0], r[3], r[4], r[5])
        if key not in seen:
            seen.add(key)
            all_rows.append(r)
    all_rows.sort(key=lambda r: (r[2], r[3], r[4], r[5]))

    # Group by class
    classes = {}
    for row in all_rows:
        classes.setdefault(row[0], []).append(row)

    with OUT.open("w") as f:

        f.write("% Auto-generated by dump_nf_bcond_table.py\n")
        f.write("% Requires: booktabs, longtable, amsmath, amssymb\n\n")
        f.write("\\renewcommand{\\arraystretch}{1.4}\n\n")

        # ---- Section header + introduction --------------------------------
        f.write("\\section*{Asymptotic freedom beyond the UV}\n\n")
        f.write(
            "For each non-Veneziano class we study the effect of adding "
            "$N_f$ fundamental chiral multiplets to one or both gauge nodes "
            "of the $\\mathrm{SU}(N)\\times\\mathrm{SU}(N)$ quiver. "
            "The asymptotic-freedom condition $b_0^{(a)} > 0$ on each node "
            "must hold not only at the free UV fixed point "
            "but also at points where one coupling has flowed to its interacting "
            "fixed point while the other remains free. "
            "At such points the bifundamental R-charges are no longer free-field "
            "values but are fixed by the interacting node's dynamics, "
            "modifying the one-loop coefficient and tightening the upper bound on $N_f$. "
            "We label three relevant points in the two-dimensional "
            "$(g_1, g_2)$ coupling space as follows:\n"
            "\\begin{itemize}\n"
            "\\item $\\mathcal{O} = (g_1=0,\\,g_2=0)$: the free UV fixed point. "
            "The one-loop AF condition $b_0^{(a)} > 0$ at $\\mathcal{O}$ "
            "gives the standard asymptotic-freedom upper bound\n"
            "  $N_f < N_f^{\\mathcal{O},(a)} \\equiv \\alpha_{\\mathcal{O}}^{(a)}\\,N + \\beta^{(a)},$\n"
            "  where $\\alpha_{\\mathcal{O}}^{(a)}$ is determined by the rank-2 and "
            "bifundamental content of node $a$.\n"
            "\\item $\\mathcal{B} = (g_1=0,\\,g_2=g_2^*)$: node 2 is at its "
            "interacting fixed point, node 1 is free. "
            "The R-charges of the bifundamentals are now fixed by node 2's "
            "fixed point. The condition $b_0^{(1)}|_{\\mathcal{B}} > 0$ gives\n"
            "  $N_f < N_f^{\\mathcal{B},(1)} \\equiv \\alpha_{\\mathcal{B}}^{(1)}\\,N + \\beta^{(1)}.$\n"
            "\\item $\\mathcal{A} = (g_1=g_1^*,\\,g_2=0)$: node 1 is at its "
            "interacting fixed point, node 2 is free. "
            "The condition $b_0^{(2)}|_{\\mathcal{A}} > 0$ gives\n"
            "  $N_f < N_f^{\\mathcal{A},(2)} \\equiv \\alpha_{\\mathcal{A}}^{(2)}\\,N + \\beta^{(2)}.$\n"
            "\\end{itemize}\n\n"
            "The theory flows to a non-trivial IR fixed point only if $N_f$ "
            "satisfies all three bounds simultaneously. "
            "The coefficient $\\alpha$ of $N$ is fixed by the morphology "
            "$(N_{r2}^{(1)}, N_{r2}^{(2)}, N_{\\rm bif})$; "
            "the constant $\\beta$ depends on the S/A character of the "
            "rank-2 tensors but does not affect the large-$N$ scaling "
            "and is not tabulated here. "
            "Each table below lists the coefficients $\\alpha$ for all "
            "morphologies within a class; negative values of "
            "$\\alpha_{\\mathcal{B}}$ or $\\alpha_{\\mathcal{A}}$ indicate "
            "that the corresponding boundary condition is automatically "
            "satisfied for large $N$ and imposes no constraint on the "
            "conformal window.\n\n"
        )

        # ---- Per-class sections -------------------------------------------
        for cid, crows in classes.items():
            a_exact = crows[0][1]
            n_total = cur.execute(
                "SELECT COUNT(*) FROM theory WHERE class_id=? AND veneziano=0",
                (cid,)
            ).fetchone()[0]
            a_tex = texify_a(a_exact)

            f.write(
                f"\\subsection*{{Class~{cid}: "
                f"$a/N^2 = {a_tex[1:-1]}$"
                f"\\quad ({n_total} theories)}}\n\n"
            )

            intro = CLASS_INTRO.get(cid, "")
            if intro:
                f.write(intro + "\n\n")

            # Per-morphology prose
            for row in crows:
                nr0, nr1, nb = row[3], row[4], row[5]
                key = (cid, nr0, nr1, nb)
                prose = MORPH_PROSE.get(key, "")
                if prose:
                    f.write(
                        f"\\textit{{Morphology {_morph_tex(nr0,nr1,nb)}.}}\\;"
                        + prose + "\n\n"
                    )

            # Table for this class: one row per morphology, α coefficients only
            label = f"tab:nf_cl{cid}"
            morph_label = (
                f"$(N_{{r2}}^{{(1)}},\\,N_{{r2}}^{{(2)}},\\,N_{{\\rm bif}})$"
            )
            caption = (
                f"Coefficients $\\alpha$ of $N$ in the upper bounds "
                f"$N_f < \\alpha\\,N + \\beta$ at each boundary point, "
                f"for all morphologies of class~{cid}. "
                r"The constant $\beta$ (not shown) depends on the S/A "
                r"character of the rank-2 matter. "
                r"Negative $\alpha_{\mathcal{B}}$ or $\alpha_{\mathcal{A}}$ "
                r"means the boundary condition is automatically satisfied "
                r"at large $N$; `conf.' means the node is already at the "
                r"conformal threshold at $\mathcal{O}$."
            )

            f.write(f"\\begin{{longtable}}{{l|cc|cc}}\n")
            f.write(f"\\caption{{{caption}}}\\label{{{label}}}\\\\\n")
            f.write("\\toprule\n")
            f.write(
                f"{morph_label}"
                r" & $\alpha_{\mathcal{O}}^{(1)}$"
                r" & $\alpha_{\mathcal{B}}^{(1)}$"
                r" & $\alpha_{\mathcal{O}}^{(2)}$"
                r" & $\alpha_{\mathcal{A}}^{(2)}$"
                " \\\\\n"
            )
            f.write("\\midrule\n\\endhead\n")

            for row in crows:
                nr0, nr1, nb = row[3], row[4], row[5]
                aO1 = _alpha_tex(_parse_alpha(row[6]))
                aO2 = _alpha_tex(_parse_alpha(row[7]))
                aB1 = _alpha_tex(_parse_alpha(row[9]))  # B_cond_B = α_B^(1)
                aA2 = _alpha_tex(_parse_alpha(row[8]))  # B_cond_A = α_A^(2)

                def _cell(v):
                    if v == "0":
                        return r"$\mathrm{conf.}$"
                    return f"${v}$"

                f.write(
                    f"{_morph_tex(nr0,nr1,nb)}"
                    f" & {_cell(aO1)}"
                    f" & {_cell(aB1)}"
                    f" & {_cell(aO2)}"
                    f" & {_cell(aA2)}"
                    " \\\\\n"
                )

            f.write("\\bottomrule\n\\end{longtable}\n\n")

        # ---- Summary section ----------------------------------------------
        f.write(
            "\\subsection*{Summary}\n\n"
            "The results above reveal three qualitatively distinct "
            "behaviours of the conformal window under fundamental deformations, "
            "summarised in Table~\\ref{tab:nf_summary}.\n\n"
            "\\begin{enumerate}\n\n"
            "\\item \\textbf{Rational compression (classes 1 and 9).} "
            "Both classes have rational R-charges at the fixed point, "
            "leading to rational $\\alpha_{\\mathcal{B}}$ and $\\alpha_{\\mathcal{A}}$. "
            "For class~1, the uniform compression factor is exactly $\\tfrac{1}{2}$: "
            "$\\alpha_{\\mathcal{B}} = 1$ against $\\alpha_{\\mathcal{O}} = 2$. "
            "For class~9, the compression follows the exact formula "
            "$\\alpha_{\\mathcal{B}} = 1 - N_{\\rm bif}/4$ across all symmetric morphologies, "
            "reflecting the fixed-point bifundamental R-charge "
            "$R_{\\rm bif} = \\tfrac{1}{2}$ of the $a/N^2 = 27/64$ solution. "
            "The constant offset $\\beta$ is identical at $\\mathcal{O}$ "
            "and $\\mathcal{B}$, so the compression acts only on the "
            "$N$-scaling of the conformal window.\n\n"
            "\\item \\textbf{$N$-independent conformal window (class 44).} "
            "At the maximum $a/N^2 = \\tfrac{1}{2}$, "
            "$\\alpha_{\\mathcal{B}} = \\alpha_{\\mathcal{A}} = 0$ "
            "for every morphology. "
            "The conformal window is bounded solely by the $\\mathcal{O}(1)$ "
            "constant $\\beta$, which can be at most 4 (for the $(2,2,2)$ "
            "morphology with both nodes carrying adjoint-type rank-2 tensors). "
            "No fundamental deformation with $N_f \\sim N$ is possible "
            "in this class, regardless of the matter configuration. "
            "Class~44 is the endpoint $N_{\\rm bif} = 4$ of the "
            "class~9 formula, extended to all morphologies.\n\n"
            "\\item \\textbf{Universal tightening (classes 2, 4, 16, 21, 23).} "
            "For all classes with irrational R-charges, the bifundamental "
            "R-charges at the fixed point satisfy $R_{\\rm bif} < \\tfrac{2}{3}$, "
            "giving a negative anomalous dimension "
            "$\\gamma_{\\rm bif} = 3R_{\\rm bif} - 2 < 0$. "
            "The resulting shift $\\Delta\\alpha = N_{\\rm bif}|\\gamma_{\\rm bif}|/2$ "
            "tightens \\emph{both} boundary conditions relative to $\\mathcal{O}$: "
            "nodes with $\\alpha_{\\mathcal{O}} > 0$ see their coefficient reduced, "
            "while nodes at the UV threshold ($\\alpha_{\\mathcal{O}} = 0$) acquire "
            "$\\alpha_{\\mathcal{B/A}} = -\\Delta\\alpha < 0$, "
            "automatically satisfying their boundary condition for large $N$. "
            "The binding (positive) boundary coefficient "
            "increases monotonically with $a/N^2$: "
            "$\\approx 0.219$ (class~2, $a/N^2 \\approx 0.152$), "
            "$\\approx 0.427$ (class~16, $a/N^2 \\approx 0.443$), "
            "$\\approx 0.675$ (class~21, $a/N^2 \\approx 0.453$), "
            "$\\approx 0.853$ (class~23, $a/N^2 \\approx 0.458$), "
            "approaching unity as $a/N^2 \\to \\tfrac{1}{2}$. "
            "Class~4 ($a/N^2 \\approx 0.168$) has the largest rank-2 content "
            "($N_{r2}^{(1)} = \\tfrac{5}{2}$) among irrational classes, "
            "giving the lighter node a wide tightened window "
            "$\\alpha_{\\mathcal{A}}^{(2)} \\approx 1.243$. "
            "In all cases the conformal window is controlled by the "
            "$\\mathcal{B}$ or $\\mathcal{A}$ condition, not by the "
            "asymptotic-freedom condition at $\\mathcal{O}$.\n\n"
            "\\end{enumerate}\n\n"
            "Across all three regimes, the constant $\\beta$ in the "
            "$\\mathcal{B}$/$\\mathcal{A}$ bounds equals the constant "
            "in the $\\mathcal{O}$ bound. "
            "The IR dynamics thus affect only the large-$N$ scaling "
            "of the conformal window, leaving the $\\mathcal{O}(1)$ "
            "structure --- which encodes the S/A character of the "
            "rank-2 matter --- unchanged.\n\n"
        )

        # Summary table
        f.write(
            "\\begin{table}[h]\n"
            "\\centering\\small\n"
            "\\renewcommand{\\arraystretch}{1.3}\n"
            "\\caption{Coefficients of $N$ in the three upper bounds "
            "on the number of additional fundamental flavours, "
            "for representative morphologies of each "
            "$\\mathrm{SU}(N)\\times\\mathrm{SU}(N)$ non-Veneziano class. "
            "For asymmetric classes with a UV-saturated node ($\\alpha_{\\mathcal{O}}=0$), "
            "that node is listed as node~1; for class~2 (no UV-saturated node), "
            "node~1 is the lighter node. "
            "Irrational values are given to four decimal places.}\n"
            "\\label{tab:nf_summary}\n"
            "\\begin{tabular}{c|c|ccc|ccc}\n"
            "\\toprule\n"
            "Class & $a/N^2$ & "
            "$\\alpha_{\\mathcal{O}}^{(1)}$ & $\\alpha_{\\mathcal{B}}^{(1)}$ & "
            "$\\alpha_{\\mathcal{A}}^{(1)}$ & "
            "$\\alpha_{\\mathcal{O}}^{(2)}$ & $\\alpha_{\\mathcal{B}}^{(2)}$ & "
            "$\\alpha_{\\mathcal{A}}^{(2)}$ \\\\\n"
            "\\midrule\n"
        )
        summary_rows = [
          # cid, a_tex,   aO1,  aB1,    aA1,   aO2,  aB2,   aA2
          (1,  r"0",           r"2", r"1",     r"1",     r"2", r"1",      r"1"),
          (2,  r"-\frac{297}{32}+\frac{135\sqrt{5}}{32}",
                               r"2", r"1.2188",r"—",  r"1", r"—",    r"0.2188"),
          (4,  r"-\frac{229}{32}+\frac{65\sqrt{13}}{32}",
                               r"0", r"-0.7569",r"—", r"2", r"—",    r"1.2431"),
          (9,  r"\frac{27}{64}\ (N_{\rm bif}=1)",
                               r"1", r"\frac{3}{4}",r"\frac{3}{4}",
                                                       r"1", r"\frac{3}{4}",r"\frac{3}{4}"),
          (9,  r"(N_{\rm bif}=2)", r"1", r"\frac{1}{2}",r"\frac{1}{2}",
                                                       r"1", r"\frac{1}{2}",r"\frac{1}{2}"),
          (9,  r"(N_{\rm bif}=3)", r"1", r"\frac{1}{4}",r"\frac{1}{4}",
                                                       r"1", r"\frac{1}{4}",r"\frac{1}{4}"),
          (9,  r"(N_{\rm bif}=4)", r"1", r"0",r"0",   r"1", r"0",    r"0"),
          (16, r"\frac{3}{32}+\frac{5\sqrt{5}}{32}",
                               r"0", r"-0.5730",r"—", r"1", r"—",    r"0.4270"),
          (21, r"-\frac{7}{2}+\frac{5\sqrt{10}}{4}",
                               r"0", r"-0.3246",r"—", r"1", r"—",    r"0.6754"),
          (23, r"\frac{2295\sqrt{17}}{86528}+\frac{30159}{86528}",
                               r"0", r"-0.1468",r"—", r"1", r"—",    r"0.8532"),
          (44, r"\frac{1}{2}", r"—", r"0",     r"0",  r"—", r"0",    r"0"),
        ]
        for row in summary_rows:
            cid_s = f"${row[0]}$"
            cols  = " & ".join(f"${v}$" for v in row[1:])
            f.write(f"{cid_s} & {cols} \\\\\n")
        f.write(
            "\\bottomrule\n"
            "\\end{tabular}\n"
            "\\end{table}\n"
        )

    con.close()
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
