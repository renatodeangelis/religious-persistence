# ── 00 · RUN THE FULL PIPELINE ─────────────────────────────────────────────────
# Sources the split scripts in dependency order. 01 builds the cleaned data, 02
# estimates and persists every matrix, and 03–07 are independent consumers that
# read the two .rds artifacts. Run from the project root.

source("code/01-prepare-data.R")      # → data/derived/gss_clean.rds
source("code/02-estimate-matrices.R") # → data/derived/matrices.rds
source("code/03-diagnostics.R")       # console tables
source("code/04-homogeneity.R")       # homogeneity figures + LaTeX table
source("code/05-memory-measures.R")   # national memory/mobility/MTE figures
source("code/06-stratified-matrices.R") # region/binary/nativity/sex/political figures
source("code/07-attitudes.R")         # attitude matrices + cohort trend figures
