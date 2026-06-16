# Religious Persistence

This project applies **Markov-chain memory measures** from the intergenerational mobility literature to measure **intergenerational persistence of religious affiliation** in the United States. The framework adapts Blume, Cholli, Durlauf & Lukina (forthcoming, *SMR*) and Wodtke, Wang, Butaeva & Durlauf (NBER w34800, 2026): rather than asking "what is the stationary memory of a pooled Markov chain?", the project estimates cohort-specific one-step transition matrices and computes implied memory curves — tracking how rapidly the influence of religious origin dissipates — across birth cohorts from 1900 to the present.

## Data

**Primary**: GSS cumulative cross-section (1974–2024), using `relig` (current affiliation) and `relig16` (childhood affiliation) to construct one-step transition matrices. Two loading paths:
- `summary-stats.R` loads via a Dropbox `.dta` URL (Stata extract)
- `transition-matrices.R` loads via the `gssr` R package (`data(gss_all)`)

**Calibration/validation**:
- **Add Health W1/W3**: Parent's Wave 1 religion (`PA22`) vs. adult child's Wave 3 retrospective recall (`H3RE26`) — tests recall bias in `matrix-validation.R`
- **Add Health W1/W4**: Parent RELTRAD split matrices, testing whether parent religiosity moderates transmission — `transition-matrices-add-health.R`

Key variables: `relig`, `relig16`, `denom`, `denom16`, `reliten`, `cohort`, `year`, `region`, `educ`

## Methods

Cohort-specific transition matrices P_t are estimated for 5-, 10-, and 20-year birth cohort windows using the RELTRAD classification (evangelical, mainline, Catholic, other, none). From each P_t:

- **Diagonal persistence rates** by origin state (most robust to small N)
- **λ₂** (second eigenvalue) — scalar memory summary
- **Individual memory curves**: log TV distance from π\* at steps t = 0–4 for each origin state
- **Overall, exchange, and structural mobility** decomposition
- **Binary (affiliated/unaffiliated) 2×2 matrices** as a fallback

A measurement section benchmarks the GSS RELIG16-based matrices against Add Health and tests for endogenous recall bias (Hayward et al. 2012).

## Code

| File | Purpose |
|------|---------|
| `Code/summary-stats.R` | Hout (2016) replication and extension, Figs 1–6 |
| `Code/transition-matrices.R` | Main matrix estimation, memory curves, mobility decomposition |
| `Code/transition-matrices-add-health.R` | Add Health W1/W4 religiosity-split matrices |
| `Code/matrix-validation.R` | Recall bias analysis (Add Health W1/W3) |
| `Code/utils.R` | Shared functions: matrix math, memory measures, plotting |
| `Code/presentation.qmd` | Quarto slide deck |

## Figures

### Hout replication (`Code/summary-stats.R` → `Figures/`)
- **Figure 1**: Share currently vs. raised Catholic, by survey year
- **Figure 2**: Decomposition of religious trajectories for those raised Catholic, by survey year
- **Figure 3**: Same decomposition by birth cohort
- **Figure 4**: Persistence and switching across Catholic, Mainline Protestant, and Conservative Protestant by birth cohort
- **Figure 5**: Diagonal persistence by education and birth cohort
- **Figure 6**: Currently vs. raised Catholic by education and survey year

### Transition matrix outputs (`Code/transition-matrices.R` → `output/figures/`)
- `output/figures/`: National cohort heatmaps (5-, 10-, 20-year bins)
- `output/figures/region/`: Regional cohort heatmaps and mobility trends
- `output/figures/binary/`: 2×2 affiliated/unaffiliated matrices and persistence trends
- `output/figures/attitude/`: Religiosity-split matrices (evolution, abortion, homosexuality, etc.)
- `im_memory_10yr.png`, `im_memory_20yr.png`: Memory curves by cohort
- `mobility_pooled.png`: Overall mobility trend
- `em_sm_10yr.png`: Exchange vs. structural mobility by cohort
