# Religious Persistence Project — Claude Context

## Project Overview

This project adapts **Markov-chain memory measures** from the class mobility literature — specifically the "memory curve" framework from Blume, Cholli, Durlauf & Lukina (forthcoming, *SMR*; NBER w33166) and its empirical application in Wodtke, Wang, Butaeva & Durlauf (NBER w34800, February 2026) — to measure **intergenerational persistence of religious affiliation** in the United States using the **General Social Survey (GSS)**.

**Core question**: How rapidly does the influence of religious origin dissipate across (synthetic) generations, and how has that rate changed across birth cohorts?

**Reframed research question**: Instead of "what is the memory of religious origin in a stationary Markov chain?" → *How has the one-step transition matrix changed across birth cohorts, and what do those changes imply for projected memory under each cohort-specific regime?* Each cohort-specific P_t yields implied memory curves "as if" that regime persisted. The comparison across cohorts is the finding.

---

## Seven Critical Weaknesses (and Status)

1. **Markov assumption violated**: Bengtson's LSOG data show grandparental religiosity has a direct independent effect on grandchildren net of parents — first-order Markov violation. Memory curves will understate true persistence.
2. **Nominal state space**: Religious affiliation has no natural ordering; "movement" between categories is treated symmetrically by the model. The scalar memory summary is not self-evidently interpretable.
3. **Non-stationarity**: The rise of the nones (5% → 25%+) is not noise; it is the phenomenon. A pooled P and its ergodic distribution have no empirical referent.
4. **RELIG16 endogenous recall bias**: Retrospective recall of childhood religion is contaminated by current adult identity. This bias is correlated with the transitions being modeled — not classical random error. Errors compound across P^k.
5. **Ergodic concept is sociologically incoherent** for contemporary U.S. religion.
6. **Ill-conditioned matrix for minority groups**: Jews, Muslims, Buddhists, Hindus will have very small N in any cohort slice.
7. **Conceptual equivalence not established**: "Religious mobility" risks being a category error — treating voluntary identity change as analogous to escaping structural disadvantage.

---

## Data Architecture

### Primary estimation sample
**GSS RELIG16 → RELIG** (full time series, 1974–present). Maximum statistical power. Transition matrices stratified by birth cohort in 15- or 20-year windows.

### Calibration/diagnostic samples

| Dataset | Purpose | Birth cohort covered |
|---------|---------|----------------------|
| **LSOG** (ICPSR 22100) | Direct validation — G2 parents observed at Wave 1 (1971); G3 adults recall childhood religion in later waves. The unique dataset for recall bias estimation. | Multi-generational |
| **Add Health** | Prospective transition matrix for birth cohort ~1976–1983; compare cell-by-cell to GSS RELIG16 matrix for same cohort | ~1976–1983 |
| **NSYR** | Wave 1 parent + youth separately → Waves 3–4 young adult religion | ~1984–1990 |
| **GSS panel (2006–2010)** | Test-retest stability of RELIG16; instability correlated with religious switching = evidence of identity contamination | — |
| **MIDUS twin subsample** | Twins who diverged in adult religiosity reporting on same childhood household | — |

**Critical data gap**: No U.S. dataset combines parent-reported affiliation at time 1 with the child's later retrospective report of childhood religion. Add Health and NSYR parent interviews provide the parent's religion but later waves do not ask the adult child "what religion were you raised in?" The paper should note this gap and recommend it for future data collection.

---

## Measures to Report

- Diagonal persistence rates by origin category across cohort windows (most robust)
- λ₂ across cohort windows with bootstrapped CIs (headline memory summary)
- Altham index d(P_t, I) and d(P_t, P_{t+1}) (does not require ordinality)
- Memory curves (P_t)^k · e_i for k = 1…10 for each origin state
- Sensitivity analysis over plausible misclassification rates

---

## State Space Decision

**Leading candidate**: 3–4 category scheme grounded in the boundary-maintenance literature (Iannaccone, Kelley): high-boundary (conservative evangelical, Mormon, Orthodox Jewish), low-boundary (mainline Protestant, Catholic), none, and possibly other/residual.

**Fallback**: Binary affiliated/unaffiliated if cell counts are too thin.

**Fine-grained (7+ categories)**: Not viable due to small-cell problems.

A Goodman RC(1) association model could estimate a latent ordering — if it recovers something like the boundary-maintenance gradient, the nominality problem is partially dissolved empirically.

---

## Measurement Section Architecture

Present **before** the main findings, organized around three questions:
1. How does the RELIG16-based transition matrix compare to independently measured matrices (Add Health, NSYR)? Which cells diverge and in which direction?
2. How internally consistent is the GSS data (RELIG16 vs. MARELIG/PARELIG)?
3. How stable is RELIG16 across panel waves, and is instability correlated with religious switching?

Conclude with: direction and approximate magnitude of bias (citing Hayward et al. 2012 for mechanism), judgment about whether bias affects the *shape* of memory curves or only their *level*, and sensitivity analysis plan.

---

## Recall Bias Analysis Plan (matrix-validation.R)

### Relevant LSOG Variables (verify names in codebook before coding)

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

**Current status** (as of April 2026): GSS downloaded; Hout (2016) replicated/extended; basic exploratory graphs complete. LSOG data downloaded to Dropbox and loaded via shared link in `matrix-validation.R`.

---

## Key Unresolved Issues

1. **Add Health codebook**: Does any later wave ask the adult child a retrospective "what religion were you raised in?" question? If yes, Add Health becomes the ideal recall-bias diagnostic.
2. **Pew ATP**: Has the childhood religion item been asked at multiple ATP waves separated by enough time to test recall drift? Worth contacting Pew directly.
3. **State-space construction**: Requires cell-count analysis. Boundary-maintenance 3–4 category scheme is leading candidate; binary affiliated/unaffiliated is the fallback.
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

## Code Files

| File | Purpose |
|------|---------|
| `Code/summary-stats.R` | GSS data loading, Hout replication, exploratory figures |
| `Code/matrix-validation.R` | LSOG-based recall bias estimation (in progress) |

## File Locations

| Data | Location |
|------|----------|
| GSS | Dropbox shared link (loaded via URL in summary-stats.R) |
| LSOG (ICPSR 22100) | Dropbox shared link (loaded via URL in matrix-validation.R) |
