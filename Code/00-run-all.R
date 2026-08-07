# ── 00 · RUN THE FULL PIPELINE ─────────────────────────────────────────────────
# Sources all numbered scripts (01–20) in dependency order. 01 builds the
# cleaned data, 02 estimates and persists every matrix, and 03–20 are consumers
# of the derived artifacts. Run from the project root.
#
# robustness-add-health.R requires the restricted Add Health files in add-health/
# (not on the gssr package), so it is deliberately kept out of the auto-run.

source("code/01-prepare-data.R")      # → data/derived/gss_clean.rds
source("code/02-estimate-matrices.R") # → data/derived/matrices.rds
source("code/03-diagnostics.R")       # console tables
source("code/04-homogeneity.R")       # homogeneity figures + LaTeX table
source("code/05-memory-measures.R")   # national memory/mobility/MTE figures
source("code/06-movement-measures.R") # OM/SM/EM/SSM + diagonal exit figures
source("code/07-lambda2-bootstrap.R") # λ₂ trend with bootstrapped CIs

# ── EXPLORATORY DIAGNOSTICS ───────────────────────────────────────────────────
source("code/08-age-standardization-explore.R")  # age-standardized λ₂/mean-diagonal sensitivity
source("code/09-period-cohort-decomp-explore.R") # age-period-cohort decomposition

# ── ROBUSTNESS CHECKS (consume data/derived/gss_clean.rds) ───────────────────
source("code/10-robustness-nonblack.R")     # non-Black sample + difference matrices
source("code/11-robustness-gss-decade.R")   # GSS survey-period stratification
source("code/12-robustness-bp.R")           # 6-state space (Black Protestant separate)
source("code/13-robustness-nativity.R")     # nativity-stratified decadal matrices
source("code/14-robustness-sex.R")          # sex-stratified decadal matrices
source("code/15-robustness-political.R")    # political stratification matrices
source("code/16-robustness-attend12.R")     # childhood attendance-stratified matrices
source("code/17-robustness-demographics.R") # demographic stratification (siblings & children)

# ── DIFFERENCE FIGURES ────────────────────────────────────────────────────────
source("code/18-diff-figures-sex.R")        # P_male − P_female + π comparison
source("code/19-diff-figures-nonblack.R")   # P_nonblack − P_full + π comparison
source("code/20-diff-figures-political.R")  # P_dem − P_rep, P_lib − P_con + π comparison
