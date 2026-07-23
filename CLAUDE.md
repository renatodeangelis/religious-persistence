# Religious Persistence Project — Claude Context

## Project Overview

This project adapts **Markov-chain memory measures** from the class mobility literature — specifically the "memory curve" framework from Blume, Cholli, Durlauf & Lukina (forthcoming, *SMR*; NBER w33166) and its empirical application in Wodtke, Wang, Butaeva & Durlauf (NBER w34800, February 2026) — to measure **intergenerational persistence of religious affiliation** in the United States using the **General Social Survey (GSS)**.

**Core question**: How rapidly does the influence of religious origin dissipate across (synthetic) generations, and how has that rate changed across birth cohorts?

**Reframed research question**: Instead of "what is the memory of religious origin in a stationary Markov chain?" → *How has the one-step transition matrix changed across birth cohorts, and what do those changes imply for projected memory under each cohort-specific regime?* Each cohort-specific P_t yields implied memory curves "as if" that regime persisted. The comparison across cohorts is the finding.

---

## Methodological Concerns (Revised Status)

### Resolved by design
- **Non-stationarity** (was #3): Not an issue. The paper never pools across cohorts or computes a single P. Each cohort-specific P_t is a self-contained object; memory curves are conditional on that cohort's regime. The rise of the nones is what the cohort comparison *measures*, not a confound.
- **Ill-conditioned minority-group matrices** (was #6): Not an issue. Minority categories (Jewish, Muslim, Buddhist, Hindu) are residualized into "other" in the 5-state scheme. Thin-cell instability is not a problem for the main estimation sample.

### Framing issues (not methodological flaws)
- **Nominal state space** (was #2): The math requires no ordinality — TV distance, Altham index, and λ₂ are all well-defined on nominal spaces. The real concern is that λ₂ is a scalar that collapses heterogeneous transition types (Catholic→None vs. Catholic→Evangelical). *Paper response*: lead with origin-state-specific memory curves and diagonal persistence; treat λ₂ as a secondary summary. Brief acknowledgment that the scalar aggregates transitions with different sociological content.
- **Ergodic concept** (was #5): π\* is a property of the matrix, not a prediction about the social world. *Paper response*: frame π\* explicitly as "the implied long-run destination under cohort t's regime, if that regime were to persist" — not an equilibrium prediction. Present it as a conditional summary alongside the memory curves, not a headline result. If it's redundant with the curves, cut it.
- **Conceptual equivalence** (was #7): The Markovian framework embeds no ordinal assumptions and requires no "structural mobility" correction. The US religious context has no analogue to the agricultural-to-service transition or postwar Japan. *Paper response*: this is a strength, not a weakness — state briefly that the model may be more appropriate for religious affiliation than for occupational mobility precisely because no structural adjustment is needed. The source of persistence (socialization and identity formation) differs from occupational gatekeeping, but the formal structure (discrete-state transition) is valid in both cases.

### Genuine concerns requiring explicit handling
- **First-order Markov assumption** (was #1): Bengtson's LSOG data show grandparental religiosity affects grandchildren net of parents. Each cohort-specific P_t captures only the single parent-to-child transmission step; the memory curves are therefore a *lower bound* on true multi-generational religious persistence. The grandparent channel is not captured. The trend comparison across cohorts is valid conditional on the grandparent effect being approximately stable across cohorts — a testable claim using LSOG. *Paper response*: explicitly state that memory curves are lower bounds on true persistence; note that the trend is valid under stable grandparent effects; use LSOG to assess whether this assumption holds.
- **RELIG16 endogenous recall bias** (was #4): The most serious surviving concern. Bias is not random — people who have left religion are more likely to recall a secular upbringing, inflating every diagonal. Critically, the bias may *trend with the cohort variable*: younger cohorts have left religion at higher rates and thus have more leavers revising their childhood recall downward. The paper's core finding (declining persistence across cohorts) could be partly a trend in bias rather than a genuine trend in transmission. *Paper response*: (a) present Add Health W1/W3 validation result prominently — it establishes magnitude; (b) discuss explicitly whether the bias trends with cohort; (c) run sensitivity analysis under estimated misclassification rates (Kuha & Skinner 1997). Diagonal persistence is more robust to this bias than λ₂ and should be treated as the primary fallback.

---

## Data Architecture

### Primary estimation sample
**GSS RELIG16 → RELIG** (full time series, 1974–present). Maximum statistical power. Transition matrices stratified by birth cohort in 5- or 10-year windows.

### Calibration/diagnostic samples

| Dataset | Purpose | Birth cohort covered |
|---------|---------|----------------------|
| **Add Health** | Prospective transition matrix for birth cohort ~1976–1983; compare cell-by-cell to GSS RELIG16 matrix for same cohort | ~1976–1983 |
| **NSYR** | Wave 1 parent + youth separately → Waves 3–4 young adult religion | ~1984–1990 |
| **GSS panel (2006–2010)** | Test-retest stability of RELIG16; instability correlated with religious switching = evidence of identity contamination | — |
| **MIDUS twin subsample** | Twins who diverged in adult religiosity reporting on same childhood household | — |

**Note on Add Health recall data**: Add Health Wave 3 includes `H3RE26`, a retrospective item asking the adult child what religion they were raised in — this is used in `matrix-validation.R` to test recall bias against `PA22` (parent's Wave 1 religion). Wave 4 does not include an equivalent retrospective item. NSYR does not appear to have a retrospective childhood religion variable in later waves. The paper should note the Wave 4 gap and the absence of a comparable item in NSYR.

---

## Measures to Report

- Diagonal persistence rates by origin category across cohort windows (most robust)
- λ₂ across cohort windows with bootstrapped CIs (headline memory summary)
- Altham index d(P_t, I) and d(P_t, P_{t+1}) (does not require ordinality)
- Memory curves (P_t)^k · e_i for k = 1…10 for each origin state
- Sensitivity analysis over plausible misclassification rates

---

## Measurement Section Architecture

Present **before** the main findings, organized around three questions:
1. How does the RELIG16-based transition matrix compare to independently measured matrices (Add Health, NSYR)? Which cells diverge and in which direction?
2. How internally consistent is the GSS data (RELIG16 vs. MARELIG/PARELIG)?
3. How stable is RELIG16 across panel waves, and is instability correlated with religious switching?

Conclude with: direction and approximate magnitude of bias (citing Hayward et al. 2012 for mechanism), judgment about whether bias affects the *shape* of memory curves or only their *level*, and sensitivity analysis plan.

---

## Recall Bias Analysis Plan

### Add Health implementation (complete — `matrix-validation.R`)
Uses `PA22` (parent's Wave 1 religion) as ground truth and `H3RE26` (adult child's Wave 3 retrospective recall) as the recall measure. Tests implemented: mismatch rate overall and by origin category, concordant-at-W1 subsample analysis, endogeneity regression (does current Wave 3 religion predict recalling a more secular upbringing?), logistic model with and without Wave 1 child religion to decompose identity vs. trajectory effects.

### LSOG implementation (pending)

#### Relevant LSOG Variables (verify names in codebook before coding)

| Variable type | What to look for |
|---------------|-----------------|
| Linkage | `FAMID` or lineage ID linking G2 parents to G3 children across waves |
| Generation | `GEN` (G1/G2/G3/G4) |
| Wave | `WAVE` indicator |
| G2's religion, Wave 1 (1971) | `RELIG71`, `DENOM`, or equivalent — directly observed parental religion, the gold standard |
| G3's current religion (each wave) | Directly observed at Waves 2–7 |
| G3's recalled childhood religion | `RELCHILD`, `RELRAISE`, `RELGRWUP`, or equivalent — the RELIG16 analogue |
| Religiosity covariates | Attendance, importance/salience at each wave |

### Step 1 — Variable audit
Identify exact variable names for (a) G2's Wave 1 religion, (b) G3's recalled childhood religion at each available wave, and (c) G3's current religion at each wave. Confirm which waves contain each item — this determines which tests are feasible.

### Step 2 — Test A: Panel stability of retrospective recall
*What it identifies*: lower bound on measurement error.

For G3 respondents interviewed at multiple waves where the recalled childhood religion item appears, cross-tabulate response at Wave t vs. Wave t+k. Recalled childhood religion is a historical fact — agreement should be near-perfect (κ ≈ 1.0). Key diagnostic: is disagreement correlated with G3's own religious switching between waves? If yes, direct evidence of identity contamination.

### Step 3 — Test B: Direct validation against observed parental religion
*What it identifies*: aggregate mismatch between recall-based "parental religion" and actual parental religion.

Link G2's Wave 1 religion to G3's recalled childhood religion at the nearest available later wave. Cross-tabulate cell by cell. Off-diagonal entries are the composite of (i) genuine change in G2's religion after 1971, (ii) recall error, and (iii) endogenous recall bias. This gives an upper bound on bias magnitude. Compute: overall mismatch rate, direction of mismatches (do G3 who left religion recall less religious upbringings?), and which origin categories have the highest mismatch rates.

### Step 4 — Test C: Endogeneity regression
*What it identifies*: whether mismatch is endogenous to G3's current religious identity.

Construct binary outcome: does G3's recalled childhood religion match G2's observed Wave 1 religion (yes/no)? Regress on G3's current religion at the survey wave, controlling for G2's actual Wave 1 denomination, wave fixed effects, and generation/age. If G3 respondents who left religion are significantly more likely to report a mismatch — especially downgrading the religiosity of their upbringing — endogenous recall bias is confirmed.

### Step 5 — Quantify implications for transition matrices
Translate the mismatch rate into a misclassification matrix: what fraction of cells in the recall-based transition matrix are likely contaminated, and in which direction? Run sensitivity analyses on λ₂ under estimated misclassification rates (Kuha & Skinner 1997 framework or simple perturbation). This becomes the "bias bounding" section of the paper.

### Key assumptions
- G2's Wave 1 religion is treated as ground truth for G3's childhood religious environment; overstated if G2's religion changed substantially between 1971 and when G3 was school-age
- Tests A and B identify different things and cannot be summed; present as complementary, not cumulative
- Heterogamy (G2 parents with different religions) complicates origin-state assignment; flag as a separate issue

---

## Theoretical Frameworks

Three frameworks map onto specific predictions about memory-curve shapes:

1. **Boundary maintenance / strict church** (Iannaccone 1994; Kelley 1972): High-demand traditions show slower memory decay because strictness screens out the uncommitted. Orthodox Judaism is the strongest empirical case.
2. **Authority erosion** (Greeley 2004; Hout 2016): Catholic memory curves should show a historical inflection post-*Humanae Vitae* (1968): pre-1968 cohorts show slow decay, post-1968 cohorts much faster.
3. **Fuzzy fidelity** (Voas 2009): Traditions with large "fuzzy" populations show initial rapid decay as nominally affiliated children leave, followed by a slower tail as the committed core persists.

Hout (2017): ~20% of Americans are "liminal" (cycling between affiliation and nonaffiliation) — complicates the absorbing-state assumption across all three frameworks.

---

## Three-Week Action Plan

| Week | Days | Tasks |
|------|------|-------|
| 1 | 1–2 | Foundational reading |
| 1 | 3 | State-space decision with cell-count analysis |
| 1 | 4 | Cohort-specific transition matrices with bootstrapped CIs |
| 1 | 5 | Eigenvalue decomposition and preliminary memory curves |
| 2 | 6–7 | Acquire/process Add Health; construct P_AddHealth; compare to P_GSS |
| 2 | 8–9 | NSYR comparison; GSS internal consistency check |
| 2 | 10 | Synthesize measurement findings in 3–5 page memo |
| 3 | 11 | Sensitivity analysis on memory curves under alternative misclassification scenarios |
| 3 | 12 | Goodman RC model exploration |
| 3 | 13 | Formal results assembly — all figures and tables |
| 3 | 14–15 | Paper skeleton; draft introduction + measurement section |

**Current status** (as of June 2026): Weeks 1–2 complete. GSS transition matrices built, λ₂ and memory curves computed for national cohort windows (5/10-year bins). Hout (2016) replicated and extended (Figs 1–6). Add Health W1/W3 recall bias analysis complete (`matrix-validation.R`); Add Health W1/W4 religiosity-split matrices complete (`robustness-add-health.R`). Pipeline reorganized: main analysis is stages 01–06, robustness re-cuts (non-Black, GSS-period, 6-state Black-Protestant) folded in as stages 10–12 consuming `gss_clean.rds`. The 20-year cohort bins, the regional stratification, and the attitude-stratified analysis were removed. NSYR comparison, sensitivity analysis, and LSOG validation remain pending (Week 3).

---

## Key Unresolved Issues

1. **Add Health W3 recall item confirmed**: Wave 3 includes `H3RE26` (retrospective childhood religion recall), used in `matrix-validation.R`. Wave 4 does not include an equivalent item.
2. **Pew ATP**: Has the childhood religion item been asked at multiple ATP waves separated by enough time to test recall drift? Worth contacting Pew directly.
3. **State-space construction**: Requires cell-count analysis. 3–4 category scheme is leading candidate; binary affiliated/unaffiliated is the fallback.
4. **λ₂ CI width**: If CIs are too wide to detect cohort trends, emphasis shifts from eigenvalue trends to diagonal persistence and Altham indices, with memory curves becoming illustrative.
5. **Ysseldyk et al. citation error**: The Overleaf file currently has the Myers (1996) title on the Ysseldyk et al. (2010) entry — needs correction.

---

## Literature

### Methods (key papers)
- Hayward, Maselko & Meador (2012) — endogenous recall bias; only published study showing retrospective childhood religious *identity* shifts to match current adult identity
- Hout & Hastings (2016) — GSS reliability
- Brenner, LaPlante & Reed (2024) — affiliation measurement inconsistency
- Hout (2017) — liminality
- **Hout 1983 *Mobility Tables*** (Sage QASS #31) — most consequential omission; bridges class mobility and sociology of religion; not yet in Overleaf
- Steensland et al. 2000 + Woodberry et al. 2012 — RELTRAD classification; not yet in Overleaf
- Kuha & Skinner (1997) — misclassification in transition matrices

### Empirics (key papers)
- Sherkat (2001) — religious mobility tables from GSS; closest methodological template
- Voas & Crockett (2005) — asymmetric transmission rates
- Hout (2016) — Catholic persistence
- Bengtson et al. (2009) — grandparent effects / Markov violation

### Methods reading list (pure statistics — not in lit review)
Blume et al. (forthcoming), Wodtke et al. 2026, Singer & Spilerman 1976, Bartholomew 1982, Goodman 1979, Altham & Ferrie 2007, Long & Ferrie 2013, Kuha & Skinner 1997, Vermunt 2010.

---

## Coding Conventions

- Use `=` as the assignment operator, not `<-`

## Code Files

| File | Purpose |
|------|---------|
| `code/summary-stats.R` | GSS data loading (Dropbox .dta), Hout replication, Figs 1–6 |
| `code/00-run-all.R` | Sources 01–06 (main pipeline) then robustness stages 10–12 in dependency order. 08/09 (explore) and `robustness-add-health.R` are deliberately **not** sourced — run those manually |
| `code/01-prepare-data.R` | Loads `gss_all`, builds/cleans the master frame (all recodes, 5/10-year cohort bins, belief, party/polviews binaries, 5-state `_alt` + 6-state `_bp` state spaces, `oversample` flag), strips haven labels → `data/derived/gss_clean.rds`. Sample: age 30–75, cohort 1925–1994 |
| `code/02-estimate-matrices.R` | Builds every count/P/pi0/pistar list for all stratifications (national 5/10-year, binary, nativity, sex, political) → `data/derived/matrices.rds`. Raw counts (N) saved for national 5/10-year (feed homogeneity tests) |
| `code/03-diagnostics.R` | Cohort-N sample-size tables, 5/10-year windows (console) |
| `code/04-homogeneity.R` | Anderson-Goodman chi-square homogeneity tests (consumes N), figures + LaTeX table → `output/figures/homogeneity/`, `output/tables/` |
| `code/05-memory-measures.R` | National IM/memory curves (5/10-year), overall/exchange/structural mobility, MTE, transition heatmaps |
| `code/08-age-standardization-explore.R` | **Exploratory** (run manually): age-standardized λ₂/mean-diagonal sensitivity, U-vs-S matrix grids → `output/figures/explore/` |
| `code/09-period-cohort-decomp-explore.R` | **Exploratory** (run manually): age-period-cohort decomposition grids and two-way effects → `output/figures/explore/apc/` |
| `code/10-robustness-nonblack.R` | **Robustness**: re-estimates national 10-year matrices on the non-Black sample (`race != 2`) + difference heatmaps vs. full sample → `output/figures/nonblack/`. Consumes `gss_clean.rds` |
| `code/11-robustness-gss-decade.R` | **Robustness**: matrices stratified by GSS survey period (not birth cohort), plus year-by-year diagonal persistence and π*. Drops `oversample` years → `output/figures/gss-decade/`. Consumes `gss_clean.rds` |
| `code/12-robustness-bp.R` | **Robustness**: 6-state scheme (Black Protestant separate) via `reltrad_bp`; 10-year memory curves + diagonal persistence → `output/figures/bp/`. Consumes `gss_clean.rds` |
| `code/13-robustness-nativity.R` | **Robustness**: nativity-stratified 10-year matrices (US-born vs. foreign-born, 3 cohort windows 1950–1979) + diagonal persistence figure → `output/figures/nativity/`. Consumes `gss_clean.rds` |
| `code/14-robustness-sex.R` | **Robustness**: sex-stratified 10-year matrices (male/female, 5 cohort windows 1940–1989) + diagonal persistence figure → `output/figures/sex/`. Consumes `gss_clean.rds` |
| `code/15-robustness-political.R` | **Robustness**: political-stratified 10-year matrices (partyid narrow/broad, polviews narrow/broad; 5 cohort windows 1940–1989) + per-variable diagonal persistence figures → `output/figures/political/`. Consumes `gss_clean.rds` |
| `code/robustness-add-health.R` | **Robustness** (run manually — needs restricted `add-health/` files): Add Health W1→W4 religiosity-split transition matrices, console-only |
| `code/matrix-validation.R` | Add Health W1/W3 recall bias analysis (PA22 vs. H3RE26) |
| `code/utils.R` | Shared functions: matrix math, memory measures, plotting helpers |
| `code/presentation.qmd` | Quarto slide deck |

## File Locations

| Data | Location |
|------|----------|
| GSS (summary-stats.R) | Dropbox shared link (loaded via URL) |
| GSS (main pipeline, 01-prepare-data.R) | `gssr` R package (`data(gss_all)`) |
| Derived matrices/clean data (02–12) | `data/derived/*.rds` (gitignored; rebuild via `code/00-run-all.R`) |
| Add Health W1/W3 (matrix-validation.R) | `add-health/` directory (w1inhome.rds, w3inhome.rds) |
| Add Health W1/W4 (robustness-add-health.R) | `add-health/` directory (w1inhome.rds, w4inhome.rds) |
| LSOG (ICPSR 22100) | Pending — not yet loaded in any script |
