# ── 00 · RUN THE FULL PIPELINE ─────────────────────────────────────────────────
# Sources the split scripts in dependency order. 01 builds the cleaned data, 02
# estimates and persists every matrix, and 03–07 are independent consumers that
# read the two .rds artifacts. 10–12 are robustness re-cuts of the main analysis
# on the same derived frame. Run from the project root.

source("code/01-prepare-data.R")      # → data/derived/gss_clean.rds
source("code/02-estimate-matrices.R") # → data/derived/matrices.rds
source("code/03-diagnostics.R")       # console tables
source("code/04-homogeneity.R")       # homogeneity figures + LaTeX table
source("code/05-memory-measures.R")   # national memory/mobility/MTE figures
source("code/06-stratified-matrices.R") # binary/nativity/sex/political figures

# ── ROBUSTNESS CHECKS (consume data/derived/gss_clean.rds) ───────────────────
source("code/10-robustness-nonblack.R")   # non-Black sample + difference matrices
source("code/11-robustness-gss-decade.R") # GSS survey-period stratification
source("code/12-robustness-bp.R")         # 6-state space (Black Protestant separate)

# ── RUN MANUALLY (not sourced here) ──────────────────────────────────────────
# 08-age-standardization-explore.R  and  09-period-cohort-decomp-explore.R are
#   exploratory diagnostics, run on demand.
# robustness-add-health.R requires the restricted Add Health files in add-health/
#   (not on the gssr package), so it is deliberately kept out of the auto-run.
