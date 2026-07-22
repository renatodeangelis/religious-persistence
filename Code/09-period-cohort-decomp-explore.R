# ── 09 · PERIOD–COHORT DECOMPOSITION (EXPLORATORY) ────────────────────────────
# EXPLORATORY / NOT part of the 00–07 pipeline. Self-contained (rebuilds its own
# frame from gss_all). Follow-up to 08: 08 showed the cross-cohort decline in
# memory survives age-standardization, so age composition is not the driver. This
# script asks the remaining question — how much of the decline is COHORT vs
# PERIOD — and is explicit about what is and isn't identified.
#
# APC identity: age = period − cohort. The three LINEAR slopes are jointly
# unidentified (null space along (+1,−1,+1)); only NONLINEAR (curvature) contrasts
# and, once age is set aside, the two-way cohort/period split are identified.
#
# Layers:
#   A  descriptive cohort×period grid of transition-matrix functionals (+ the
#      structurally banded N grid that visualizes the identification problem)
#   B  two-way additive fit y_{c,p} = μ + γ_c + δ_p (identified given age set
#      aside); variance shares answer "how much reattributes to period"
#   C  individual estimable second-difference (curvature) contrasts for two focal
#      transitions (leaving religion; returning from none) — identified regardless
#   D  HAPC variance partition (lme4) as an assumption-dependent robustness cross-
#      check; skipped if lme4 is unavailable
#
# Output: console tables + output/figures/explore/apc/*.png

suppressMessages({library(dplyr); library(tidyr); library(gssr); library(ggplot2); library(splines)})
source("code/utils.R")

S = rel_level_order
lambda2 = function(P) sort(Mod(eigen(P)$values), decreasing = TRUE)[2]

# transition matrix + π₀ from a subframe
pm = function(df, levels = S) {
  N = table(factor(df$reltrad16_alt, levels), factor(df$reltrad_alt, levels))
  N = matrix(as.numeric(N), length(levels), length(levels), dimnames = list(levels, levels))
  list(P = N / rowSums(N), pi0 = rowSums(N) / sum(N), N = N)
}
retention = function(m) sum(m$pi0 * diag(m$P))                 # π₀-weighted diagonal

# ── DATA (self-contained; config B core cohorts 1925–1984) ───────────────────

reltrad_labels = c("1"="evangelical","2"="mainline","3"="black protestant",
  "4"="catholic","5"="jewish","6"="other","7"="none")
data(gss_all)
d = gss_all |>
  select(year, cohort, reltrad, reltrad16) |>
  filter(!(year %in% c(1972, 1982, 1987, 2021, 2022, 2024)), year < 2020) |>
  mutate(across(c(reltrad, reltrad16), ~ reltrad_labels[as.character(as.numeric(.))])) |>
  filter(!is.na(reltrad), !is.na(reltrad16)) |>
  mutate(across(c(reltrad, reltrad16),
    ~ case_when(. == "jewish" ~ "other", . == "black protestant" ~ "evangelical", TRUE ~ .),
    .names = "{.col}_alt")) |>
  mutate(cohort = as.numeric(cohort), age = year - cohort) |>
  filter(age >= 30, age <= 75, cohort >= 1925, cohort <= 1984)
d = as.data.frame(d)
d$coh <- (floor((d$cohort - 1905) / 10) * 10 + 1905) + 5      # 5-basis 10-yr midpoint
d$per <- floor((d$year   - 1900) / 10) * 10 + 1900            # period decade (1970..2010)

MINCELL = 150
dir.create("output/figures/explore/apc", recursive = TRUE, showWarnings = FALSE)

# ── PART A · COHORT × PERIOD GRID ────────────────────────────────────────────

grid = expand.grid(coh = sort(unique(d$coh)), per = sort(unique(d$per)))
grid = do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
  cc = grid$coh[k]; pp = grid$per[k]
  s = d[d$coh == cc & d$per == pp, ]
  n = nrow(s)
  if (n < MINCELL) return(data.frame(coh = cc, per = pp, n = n, age = NA,
                                     retention = NA, lambda2 = NA, pistar_none = NA))
  m = pm(s)
  data.frame(coh = cc, per = pp, n = n, age = round(mean(s$age), 1),
             retention = round(retention(m), 3), lambda2 = round(lambda2(m$P), 3),
             pistar_none = round(pi_star(m$P)["none"], 3))
}))

cat("==== PART A: cohort × period cell counts (blank = structurally impossible / n<150) ====\n")
print(pivot_wider(grid[, c("coh","per","n")], names_from = per, values_from = n), n = Inf)
cat("\n==== mean age per cell (note: pinned by age = period − cohort) ====\n")
print(pivot_wider(grid[, c("coh","per","age")], names_from = per, values_from = age), n = Inf)
cat("\n==== overall retention per cell ====\n")
print(pivot_wider(grid[, c("coh","per","retention")], names_from = per, values_from = retention), n = Inf)

# N heatmap (shows the band) + retention heatmap
gf = grid |> mutate(coh = factor(coh), per = factor(per))
ggsave("output/figures/explore/apc/grid_N.png",
  ggplot(gf, aes(per, coh, fill = n)) + geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(n), "", n)), size = 3) +
    scale_fill_distiller(palette = "Greens", direction = 1, na.value = "grey92") +
    labs(x = "Survey period (decade)", y = "Birth cohort (10-yr midpoint)",
         title = "Cohort × period cell N — the structurally banded APC grid") + healy_theme,
  width = 8, height = 5, dpi = 200)
ggsave("output/figures/explore/apc/grid_retention.png",
  ggplot(gf, aes(per, coh, fill = retention)) + geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(retention), "", sprintf("%.2f", retention))), size = 3) +
    scale_fill_distiller(palette = "Blues", direction = 1, na.value = "grey92", limits = c(0.5, 0.9)) +
    labs(x = "Survey period (decade)", y = "Birth cohort (10-yr midpoint)",
         title = "Overall retention by cohort × period") + healy_theme,
  width = 8, height = 5, dpi = 200)
ggsave("output/figures/explore/apc/grid_pistar_none.png",
  ggplot(gf, aes(per, coh, fill = pistar_none)) + geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(pistar_none), "", sprintf("%.2f", pistar_none))), size = 3) +
    scale_fill_distiller(palette = "Oranges", direction = 1, na.value = "grey92") +
    labs(x = "Survey period (decade)", y = "Birth cohort (10-yr midpoint)",
         title = "Steady-state share 'None' (π*) by cohort × period") + healy_theme,
  width = 8, height = 5, dpi = 200)

# ── PART B · TWO-WAY ADDITIVE DECOMPOSITION (identified given age set aside) ──
# y_{c,p} = μ + γ_c + δ_p ; weighted by cell N. Variance shares partition the
# explained variation into cohort vs period. Assumes age effects absorbed by the
# (age-homogeneous) cell structure — read RELATIVE shares, not absolute levels.

twoway = function(yvar) {
  g = grid[!is.na(grid[[yvar]]), ]
  g$cohf = factor(g$coh); g$perf = factor(g$per)
  fit = lm(reformulate(c("cohf","perf"), yvar), data = g, weights = g$n)
  aov_tab = anova(fit)
  ss = aov_tab[["Sum Sq"]]; names(ss) = rownames(aov_tab)
  share = round(ss / sum(ss), 3)
  cat(sprintf("\n---- Two-way decomposition of %s (n_cells=%d) ----\n", yvar, nrow(g)))
  cat("Variance share  cohort:", share["cohf"], " period:", share["perf"],
      " residual:", share["Residuals"], "\n")
  # effect estimates (sum-to-zero) for plotting
  eff = function(f) { e = tapply(g[[yvar]], g[[f]], function(x) mean(x)); e - mean(e) }
  list(fit = fit, share = share,
       gamma = eff("cohf"), delta = eff("perf"))
}
cat("\n==== PART B: two-way additive fits ====\n")
tw_ret <- twoway("retention")
tw_l2  <- twoway("lambda2")

eff_df = bind_rows(
  data.frame(factor = "cohort", level = as.numeric(names(tw_ret$gamma)), effect = tw_ret$gamma, y = "retention"),
  data.frame(factor = "period", level = as.numeric(names(tw_ret$delta)), effect = tw_ret$delta, y = "retention"),
  data.frame(factor = "cohort", level = as.numeric(names(tw_l2$gamma)),  effect = tw_l2$gamma,  y = "lambda2"),
  data.frame(factor = "period", level = as.numeric(names(tw_l2$delta)),  effect = tw_l2$delta,  y = "lambda2"))
ggsave("output/figures/explore/apc/twoway_effects.png",
  ggplot(eff_df, aes(level, effect, color = factor, group = factor)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_line() + geom_point(size = 2) +
    facet_wrap(~ y, scales = "free_y") +
    scale_color_manual(values = c(cohort = "#D55E00", period = "#0072B2"), name = NULL) +
    labs(x = "Cohort midpoint / period decade", y = "Additive effect (sum-to-zero)",
         title = "Two-way cohort vs period effects (age set aside)") + healy_theme,
  width = 10, height = 5, dpi = 200)

# ── PART C · ESTIMABLE SECOND-DIFFERENCE (CURVATURE) CONTRASTS ────────────────
# Identified regardless of the linear-slope indeterminacy. Two focal transitions:
#   leaving : among origin != none, P(current == none)
#   return  : among origin == none, P(current != none)
# Model: logit = α + ns(age) + factor(cohort) + factor(period). Extract Δ²γ (and
# Δ²δ) with delta-method SEs from vcov — the estimable curvature.

second_diff = function(fit, prefix) {
  b = coef(fit); V = vcov(fit)
  idx = grep(paste0("^", prefix), names(b))
  if (length(idx) < 2) return(NULL)
  # full effect vector with reference level = 0 prepended
  levs = as.numeric(sub(prefix, "", names(b)[idx]))
  ord  = order(levs); idx <- idx[ord]; levs <- levs[ord]
  e_names = names(b)[idx]
  # build contrast rows for Δ² over the ordered levels including implicit ref (0)
  full_lv = c(min(levs) - diff(levs)[1], levs)         # reference sits one step before
  k = length(idx)
  out = lapply(2:(k - 1), function(j) {
    L = numeric(length(b))
    L[idx[j - 1]] = 1; L[idx[j]] = -2; L[idx[j + 1]] = 1
    est = sum(L * b); se = sqrt(as.numeric(t(L) %*% V %*% L))
    data.frame(center = levs[j], d2 = est, se = se, z = est / se)
  })
  do.call(rbind, out)
}

fit_focal = function(sub, ylab) {
  sub$y = ylab(sub)
  sub$cohf = factor(sub$coh); sub$perf = factor(sub$per)
  glm(y ~ ns(age, 4) + cohf + perf, data = sub, family = binomial())
}
leave_sub  = d[d$reltrad16_alt != "none", ]
return_sub = d[d$reltrad16_alt == "none", ]
fit_leave  = fit_focal(leave_sub,  function(x) as.integer(x$reltrad_alt == "none"))
fit_return = fit_focal(return_sub, function(x) as.integer(x$reltrad_alt != "none"))

cat("\n==== PART C: estimable cohort curvature (Δ²) — 'leaving religion' logit ====\n")
sd_leave_c = second_diff(fit_leave, "cohf"); print(round(sd_leave_c, 3), row.names = FALSE)
cat("\n==== estimable period curvature (Δ²) — 'leaving religion' logit ====\n")
sd_leave_p = second_diff(fit_leave, "perf"); print(round(sd_leave_p, 3), row.names = FALSE)
cat("\n==== estimable cohort curvature (Δ²) — 'return from none' logit ====\n")
sd_ret_c = second_diff(fit_return, "cohf"); print(round(sd_ret_c, 3), row.names = FALSE)

curv_df = bind_rows(
  transform(sd_leave_c, contrast = "leaving · cohort"),
  transform(sd_leave_p, contrast = "leaving · period"),
  transform(sd_ret_c,   contrast = "return · cohort"))
ggsave("output/figures/explore/apc/curvature_contrasts.png",
  ggplot(curv_df, aes(center, d2)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_pointrange(aes(ymin = d2 - 1.96 * se, ymax = d2 + 1.96 * se)) +
    facet_wrap(~ contrast, scales = "free") +
    labs(x = "Center of the three-point window (cohort midpoint / period)",
         y = "Estimable second difference Δ² (log-odds)",
         title = "Identified APC curvature: where cohort/period effects bend",
         subtitle = "Nonzero Δ² = a real bend, invariant to the unidentified linear slope") + healy_theme,
  width = 11, height = 4, dpi = 200)

# ── PART D · HAPC VARIANCE PARTITION (robustness; skipped if no lme4) ─────────

if (requireNamespace("lme4", quietly = TRUE)) {
  cat("\n==== PART D: HAPC cross-classified random effects (leaving logit) ====\n")
  m = lme4::glmer(y ~ ns(age, 3) + (1 | cohf) + (1 | perf),
                  data = transform(leave_sub, y = as.integer(reltrad_alt == "none"),
                                   cohf = factor(coh), perf = factor(per)),
                  family = binomial())
  vc = as.data.frame(lme4::VarCorr(m))
  cat("Between-cohort var:", round(vc$vcov[vc$grp == "cohf"], 4),
      " between-period var:", round(vc$vcov[vc$grp == "perf"], 4), "\n")
  cat("Cohort share of group variance:",
      round(vc$vcov[vc$grp == "cohf"] / sum(vc$vcov[vc$grp %in% c("cohf","perf")]), 3), "\n")
} else {
  cat("\n[Part D skipped: lme4 not installed. install.packages('lme4') to enable HAPC.]\n")
}

cat("\nDone. Figures in output/figures/explore/apc/.\n",
    "Read Part B for the relative cohort/period split (given age set aside) and\n",
    "Part C for the identified bends (e.g., Catholic/1968 cohort, 1990s period).\n")
