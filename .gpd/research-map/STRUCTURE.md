# Project Structure

**Analysis Date:** 2026-04-14

## Directory Layout

```
quiver/
|-- beta_functions.py            # Module 1: Dynkin indices, b_0, AF checks
|-- quiver_generation.py         # Module 2: Quiver dataclass, enumeration, anomaly checks
|-- a_maximization.py            # Module 3: numerical a-max at fixed (N, N_f)
|-- a_maximization_large_N.py    # Module 3: exact symbolic large-N a-max + Mathematica batch
|-- two_node_db.py               # Module 5: SQLite DB builder + CLI
|-- quivers.db                   # Pre-built database (~2.5 MB, 7757 theories)
|-- README.md                    # Pipeline overview and key references
|-- CLAUDE.md                    # AI assistant instructions
|-- plan.md                      # Original project plan
|-- module1_beta_functions.md    # Module 1 specification
|-- module2_quiver_generation.md # Module 2 specification
|-- module3_a_maximization.md    # Module 3 specification
|-- module4_ir_fixed_points.md   # Module 4 specification (not yet implemented)
|-- module5_classification.md    # Module 5 specification
|-- __pycache__/                 # Python bytecode cache (generated)
|-- .claude/                     # Claude settings
|-- .gpd/                        # GPD analysis directory
|   +-- research-map/            # Research mapping documents
|-- .git/                        # Git repository
```

## Directory Purposes

**Root directory (`quiver/`):**
- Purpose: Flat project layout; all source code, documentation, and data live at root level
- Contains: `.py` (source), `.md` (documentation/specs), `.db` (SQLite database)
- No subdirectories for code, data, or docs -- everything is at the top level

**`__pycache__/`:**
- Purpose: Python bytecode cache
- Generated: Yes (automatic)
- Committed: No (should be in .gitignore but currently tracked)

## Key File Locations

**Theory / Specifications:**
- `module1_beta_functions.md`: one-loop b_0 computation formulas, Dynkin index tables, large-N AF constraints
- `module2_quiver_generation.md`: quiver enumeration algorithm, anomaly cancellation rules, edge rep types
- `module3_a_maximization.md`: a-maximization formalism, anomaly-free constraints, unitarity bounds
- `module4_ir_fixed_points.md`: IR fixed point boundary analysis via B_a = Tr[R^(a) G_a^2] < 0
- `module5_classification.md`: full pipeline specification, output format, classification criteria

**Computation / Implementation:**
- `beta_functions.py` (288 lines): Dynkin indices (`T_adj`, `T_rep`, `T_bifund`), anomaly coefficients (`A_rep`), beta function (`compute_b0`, `b0_linear`), AF checks (`is_af`, `is_af_all_N`)
- `quiver_generation.py` (~780 lines): `Quiver` and `Edge` dataclasses, N_f bound computation, chiral excess, Witten anomaly check, combinatorial enumeration with backtracking, deduplication by conjugation + node swap
- `a_maximization.py` (~530 lines): `ChiralField` dataclass, field builder at fixed (N, N_f), anomaly matrix construction, numerical a-maximization via SciPy BFGS, iterative unitarity decoupling, gauge-invariant operator catalog
- `a_maximization_large_N.py` (~1000 lines): `LeadField` dataclass, leading-order field builder, symbolic anomaly constraints, exact a-maximization via SymPy, Mathematica batch NSolve via wolframscript subprocess, formatting utilities
- `two_node_db.py` (~1100 lines): SQLite schema, DB build pipeline (4 phases), CLI commands (build, classes, show, search, stats, morph-build), clustering, morphology classification

**Data / Results:**
- `quivers.db` (SQLite, ~2.5 MB): pre-built database containing 7757 theories in 326 universality classes and 580 morphology classes. Built 2026-04-13.
  - Table `theory`: individual quiver theories with matter content, edges, R-charges, central charges
  - Table `universality_class`: classes grouped by a/N^2 value (182 exact, 144 numerical-only)
  - Table `morphology_class`: structural classification by (N_rank2, N_bif, N_fund) vector
  - Table `build_info`: build metadata (timestamp, parameters, counts)

## Computation Dependencies

**Import Graph:**
```
beta_functions.py            (no internal imports)
       ^
       |
quiver_generation.py         (imports from beta_functions)
       ^
       |
a_maximization.py            (imports from beta_functions, quiver_generation)
       ^
       |
a_maximization_large_N.py    (imports from beta_functions, quiver_generation)
       ^
       |
two_node_db.py               (imports from quiver_generation, a_maximization_large_N)
```

**External Library Dependencies:**
```
beta_functions.py:           fractions (stdlib)
quiver_generation.py:        fractions, itertools (stdlib)
a_maximization.py:           numpy, scipy.linalg, scipy.optimize
a_maximization_large_N.py:   sympy, subprocess, json, tempfile (+ wolframscript on PATH)
two_node_db.py:              sqlite3, argparse, signal, time, re, numpy (via imports)
```

**External Tool Dependencies:**
- `wolframscript` (Mathematica): required for `a_maximize_batch_mathematica` in `a_maximization_large_N.py`. Invoked via `subprocess.run`. Must be on system PATH.
- No `requirements.txt`, `pyproject.toml`, `setup.py`, or `Makefile` exists. Dependencies are implicit.

## Naming Conventions

**Files:**
- Python modules: `{physics_concept}.py` or `{physics_concept}_{variant}.py`
  - e.g., `beta_functions.py`, `a_maximization.py`, `a_maximization_large_N.py`
- Specification documents: `module{N}_{topic}.md`
  - e.g., `module1_beta_functions.md`, `module3_a_maximization.md`
- Database: `quivers.db` (SQLite)
- CLI entry point: `two_node_db.py` (Module 5, doubles as library and CLI)

**Variables in Code:**
- Gauge types: `"SU"`, `"SO"`, `"Sp"` (string literals, typed as `GaugeType = Literal["SU", "SO", "Sp"]`)
- Representations: `"fund"`, `"antifund"`, `"adj"`, `"S"`, `"Sbar"`, `"A"`, `"Abar"`, `"V"` (string literals)
- Edge representations: `"+-"`, `"++"`, `"--"`, `"+"`, `"-"`, `"std"` (string literals)
- Beta function coefficient: `b0`, `b_0`
- Dynkin index: `T_adj`, `T_rep`, `T_bifund`
- R-charges: `R`, `R_opt`, `R_values`, `R_charges`
- Central charges: `a_val`, `c_val`, `a_over_N2`, `c_over_N2`
- Null space matrix: `F`, `Fmat`, `Fnull`
- Free parameters: `s`, `s_opt`, `svars`
- Anomaly matrix: `A` (overloaded -- also used for antisymmetric rep label)
- Node rank parameter: `N` (shared base), `N_i = mults[i] * N` (effective rank)
- Rank multipliers: `mults`, `rank_multipliers`, `m_a`, `m_b`, `m0`, `m1`

**Database Columns:**
- Snake_case: `gauge_pair`, `a_over_N2`, `class_id`, `theory_id`, `morph_id`
- Display text stored as formatted strings: `matter0`, `matter1`, `edges`, `delta0`, `R_numerical`

## Where to Add New Content

**New Gauge Group Type:**
1. Add to `GaugeType` literal in `beta_functions.py` (line 13)
2. Add to `N_MIN`, `VALID_REPS`, `RANK2_ADJ_REPS` dicts in `beta_functions.py`
3. Implement `T_adj`, `T_rep` cases in `beta_functions.py`
4. Add edge rep types in `EDGE_REPS` dict in `quiver_generation.py` (line 40)
5. Add dimension formulas in `a_maximization.py` (`dim_group`, `dim_rep`)
6. Add leading-order coefficients in `a_maximization_large_N.py`

**New Representation Type:**
1. Add to `VALID_REPS` and `RANK2_ADJ_REPS` in `beta_functions.py`
2. Add `T_rep` case in `beta_functions.py`
3. Add `dim_rep` case in `a_maximization.py`
4. Add leading-order coefficients in `a_maximization_large_N.py`
5. Add gauge-invariant operator if applicable in `gauge_invariant_ops` (`a_maximization.py`, line 266)

**New Quiver Topology (e.g., 3-node quivers):**
1. Use existing `enumerate_quivers(n_nodes=3, ...)` -- already supports arbitrary n_nodes
2. Extend `two_node_db.py` or create `three_node_db.py` for the classification database
3. Node-swap deduplication currently hardcodes 2-node swap; generalize `_node_swap_signature` for n > 2

**Module 4 Implementation (IR Fixed Points):**
1. Create `ir_fixed_points.py` at project root
2. Implement B_a = Tr[R^(a) G_a^2] computation using R-charges from a-maximization
3. Check B_a < 0 at the boundary g_a -> 0 for each node
4. Integrate into `two_node_db.py` build pipeline as an additional filter
5. Specification: `module4_ir_fixed_points.md`

**New Limiting Case / Cross-Check:**
- Each `.py` file has `if __name__ == "__main__"` validation blocks at the bottom
- `beta_functions.py` (line 261): SQCD examples with known b_0 values
- `quiver_generation.py` (line 710): N_f bound validation against known results
- `a_maximization.py` (line 467): SQCD R-charges against analytic R = 1 - N/N_f
- Add new validation examples in the `__main__` block of the relevant module

## Build and Execution

**Running the Classification Pipeline:**
```bash
# Build the full database (~3 min, requires Mathematica)
python two_node_db.py build

# Build with custom parameters
python two_node_db.py build --force --exact-timeout 60

# Rebuild on existing DB (--force to overwrite)
python two_node_db.py build --force
```

**Querying the Database:**
```bash
# List all universality classes
python two_node_db.py classes

# Filter by gauge pair
python two_node_db.py classes --pair SU-SU

# Filter by a/N^2 range
python two_node_db.py classes --min-a 0.5 --max-a 1.0

# Show all theories in a class
python two_node_db.py show 7

# Show diverged theories (no convergence)
python two_node_db.py show

# Search theories by properties
python two_node_db.py search --matter0 adj --pair SU-SU

# Database summary statistics
python two_node_db.py stats
```

**Running Individual Module Validations:**
```bash
python beta_functions.py          # Module 1 examples
python quiver_generation.py       # Module 2 examples + 3-node enumeration
python a_maximization.py          # Module 3 SQCD validation
```

## Input/Output Formats

**Input:**
- No external input files. All theory parameters are generated programmatically.
- Command-line arguments for DB build configuration (max_multiedge, exact_timeout, etc.)

**Output:**
- `quivers.db`: SQLite database (primary output)
- Console: formatted tables for CLI commands (classes, show, search, stats)
- Temporary files: Mathematica script (`mma_batch.wl`) and results (`mma_results.json`) written to `tempfile.mkdtemp()`, not cleaned up

**Intermediate Data Formats:**
- Mathematica batch: JSON file written by wolframscript containing per-theory results `{idx, a, c, R}`. Requires post-processing for Mathematica's non-standard decimal format (e.g., `"0."` -> `"0.0"`).
  - File: `a_maximization_large_N.py` (lines 886-896, JSON fix regex)

## Special Directories

**`__pycache__/`:**
- Purpose: Python bytecode cache
- Generated: Yes (automatic)
- Committed: Should not be (currently appears in git status as untracked)

**`.gpd/research-map/`:**
- Purpose: GPD analysis documents
- Generated: Yes (by GPD tools)
- Committed: No

---

_Structure analysis: 2026-04-14_
