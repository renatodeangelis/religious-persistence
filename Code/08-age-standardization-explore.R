# ── 08 · AGE-STANDARDIZATION (EXPLORATORY) ────────────────────────────────────
# EXPLORATORY / NOT part of the 00–07 pipeline. Self-contained: rebuilds its own
# frame from gss_all so it does not depend on 01's (still-being-decided) cohort
# floor / bin basis. Purpose: test whether the cross-cohort trend in memory (λ₂)
# and diagonal persistence survives holding AGE roughly fixed — the central
# life-cycle-bias robustness check.
#
# Design (config B core): 5-basis 10-year cohort bins, floor relaxed to 1925 so
# the leading bin (1925–34) is a full decade. Three estimates per cohort bin:
#   (U) unstandardized — full age support within the bin (what 02 currently does)
#   (R) band-restricted, unweighted — keep only the common age band
#   (S) direct age-standardization — within-band, reweight each cohort to a shared
#       reference age distribution so every cohort has the SAME age composition
# If the cohort trend in (U) is also present in (R)/(S), it is not an age artifact.
#
# The youngest bin (1985–94) has almost no age spread (SD ≈ 2.7) and cannot cover
# the band; it is reported UNSTANDARDIZED ONLY and flagged, never standardized.
#
# Output: console tables + output/figures/explore/*.png

suppressMessages({library(dplyr); library(tidyr); library(gssr); library(ggplot2)})
source("code/utils.R")

S = rel_level_order   # c("catholic","evangelical","mainline","other","none")

# ── LOCAL HELPERS (weighted matrices + second eigenvalue) ────────────────────

# Modulus of the second-largest eigenvalue of P (the memory summary λ₂).
lambda2 = function(P) sort(Mod(eigen(P)$values), decreasing = TRUE)[2]

# Weighted row-stochastic transition matrix. w defaults to 1 (unweighted).
wp_matrix = function(df, w = NULL, levels = S) {
  if (is.null(w)) w = rep(1, nrow(df))
  N = xtabs(w ~ factor(df$reltrad16_alt, levels) + factor(df$reltrad_alt, levels))
  N = matrix(as.numeric(N), length(levels), length(levels), dimnames = list(levels, levels))
  list(P = N / rowSums(N), pi0 = rowSums(N) / sum(N), N = N)
}

# π₀-weighted mean diagonal (overall retention = 1 − overall_mobility).
mean_diag = function(P, pi0) sum(pi0 * diag(P))

# ── DATA (self-contained; mirrors 01's recodes, floor = 1925) ─────────────────

reltrad_labels = c("1"="evangelical","2"="mainline","3"="black protestant",
  "4"="catholic","5"="jewish","6"="other","7"="none")
data(gss_all)
d = gss_all |>
  select(year, cohort, reltrad, reltrad16) |>
  filter(!(year %in% c(1972, 2021))) |>
  mutate(across(c(reltrad, reltrad16), ~ reltrad_labels[as.character(as.numeric(.))])) |>
  filter(!is.na(reltrad), !is.na(reltrad16)) |>
  mutate(across(c(reltrad, reltrad16),
    ~ case_when(. == "jewish" ~ "other", . == "black protestant" ~ "evangelical", TRUE ~ .),
    .names = "{.col}_alt")) |>
  mutate(cohort = as.numeric(cohort), age = year - cohort) |>
  filter(age >= 30, age <= 75, cohort >= 1925, cohort <= 1994)
d = as.data.frame(d)
d$bin = floor((d$cohort - 1905) / 10) * 10 + 1905          # 5-basis 10-yr bins
d$mid = d$bin + 5                                          # midpoint label (x-axis)

bins   = sort(unique(d$bin))
young  = 1985                                              # age-degenerate bin, U only

# ── ESTIMATION ────────────────────────────────────────────────────────────────
# Two candidate common age bands: a wide one (all core bins overlap 39–49) and a
# tight "religion at ~age 41" fixed-age band (39–43, rec. #2 in the referee memo).

run_band = function(band_lo, band_hi, band_lbl) {
  # cohorts whose age support covers the whole band can be standardized
  cov = d |> group_by(bin) |>
    summarise(lo = min(age), hi = max(age), .groups = "drop") |>
    mutate(covers = lo <= band_lo & hi >= band_hi)
  standardizable = cov$bin[cov$covers]

  # reference age distribution: pooled over standardizable cohorts within the band
  ref = d |> filter(age >= band_lo, age <= band_hi, bin %in% standardizable) |>
    count(age) |> mutate(p_ref = n / sum(n)) |> select(age, p_ref)

  mats = list()                                            # keep P/π₀ per bin × method
  rows = lapply(bins, function(b) {
    sub  = d[d$bin == b, ]                                  # (U) full support
    subb = sub[sub$age >= band_lo & sub$age <= band_hi, ]   # within band

    U = wp_matrix(sub)
    entry = list(U = U)
    out = data.frame(mid = b + 5, bin = paste0(b, "-", b + 9),
                     method = "U (unstd, full age)", n = nrow(sub),
                     lambda2 = lambda2(U$P), mean_diag = mean_diag(U$P, U$pi0),
                     cath_diag = U$P["catholic","catholic"], none_diag = U$P["none","none"])

    if (b %in% standardizable && nrow(subb) >= 100) {
      R = wp_matrix(subb)                                   # (R) band, unweighted
      entry$R = R
      out = rbind(out, data.frame(mid = b + 5, bin = paste0(b, "-", b + 9),
                     method = "R (band, unwtd)", n = nrow(subb),
                     lambda2 = lambda2(R$P), mean_diag = mean_diag(R$P, R$pi0),
                     cath_diag = R$P["catholic","catholic"], none_diag = R$P["none","none"]))

      # (S) direct standardization: weight = p_ref(age) / p_cohort(age)
      fc = subb |> count(age) |> mutate(p_coh = n / sum(n)) |> select(age, p_coh)
      wsub = subb |> left_join(ref, by = "age") |> left_join(fc, by = "age") |>
        mutate(w = ifelse(is.na(p_ref) | p_coh == 0, 0, p_ref / p_coh))
      Sm = wp_matrix(wsub, w = wsub$w)
      entry$S = Sm
      out = rbind(out, data.frame(mid = b + 5, bin = paste0(b, "-", b + 9),
                     method = "S (age-standardized)", n = nrow(subb),
                     lambda2 = lambda2(Sm$P), mean_diag = mean_diag(Sm$P, Sm$pi0),
                     cath_diag = Sm$P["catholic","catholic"], none_diag = Sm$P["none","none"]))
    }
    mats[[as.character(b)]] <<- entry
    out
  })
  res = do.call(rbind, rows)
  res[, c("lambda2","mean_diag","cath_diag","none_diag")] =
    round(res[, c("lambda2","mean_diag","cath_diag","none_diag")], 3)

  cat("\n================ COMMON BAND:", band_lbl, "(ages", band_lo, "-", band_hi, ") ================\n")
  cat("Standardizable cohorts (cover full band):", paste(standardizable, collapse = ", "), "\n")
  cat("Reference age dist: pooled within band over standardizable cohorts (N =",
      sum(d$age >= band_lo & d$age <= band_hi & d$bin %in% standardizable), ")\n\n")
  print(res[order(res$mid), ], row.names = FALSE)
  res$band = band_lbl
  list(res = res, mats = mats)
}

wide_out  = run_band(39, 49, "wide 39-49")
fixed_out = run_band(39, 43, "fixed ~41 (39-43)")
res_wide  = wide_out$res
res_fixed = fixed_out$res

# ── FIGURES ───────────────────────────────────────────────────────────────────

dir.create("output/figures/explore", recursive = TRUE, showWarnings = FALSE)

plot_metric = function(res, ycol, ylab, ttl, fname) {
  df = res
  df$method = factor(df$method,
    levels = c("U (unstd, full age)", "R (band, unwtd)", "S (age-standardized)"))
  p = ggplot(df, aes(x = mid, y = .data[[ycol]], color = method, group = method)) +
    geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
    scale_x_continuous(breaks = bins + 5) +
    scale_color_manual(values = c("U (unstd, full age)" = "#999999",
                                  "R (band, unwtd)"     = "#0072B2",
                                  "S (age-standardized)"= "#D55E00"), name = NULL) +
    labs(x = "Birth cohort (10-yr bin midpoint)", y = ylab, title = ttl) +
    healy_theme
  ggsave(fname, p, width = 8, height = 5, dpi = 200)
}

plot_metric(res_wide,  "lambda2",   "λ₂ (second eigenvalue)",
            "Memory (λ₂) by cohort: does the trend survive age-standardization? [band 39–49]",
            "output/figures/explore/age_std_lambda2_wide.png")
plot_metric(res_wide,  "mean_diag", "π₀-weighted mean diagonal (retention)",
            "Overall retention by cohort, unstandardized vs age-standardized [band 39–49]",
            "output/figures/explore/age_std_meandiag_wide.png")
plot_metric(res_fixed, "lambda2",   "λ₂ (second eigenvalue)",
            "Memory (λ₂) by cohort, fixed-age ~41 robustness [band 39–43]",
            "output/figures/explore/age_std_lambda2_fixed.png")

# ── TRANSITION-MATRIX HEATMAPS (make_combined, as in 05) ─────────────────────
# Per-cohort P + π₀ + π* panels for the WIDE band, unstandardized (U) and
# age-standardized (S); plus a faceted P-heatmap grid so the cross-cohort change
# is legible at a glance. Uses the wide band (39–49) — the reliable one.

dir.create("output/figures/explore/matrices", recursive = TRUE, showWarnings = FALSE)

for (b in names(wide_out$mats)) {
  e = wide_out$mats[[b]]; edge = as.integer(b); rng = paste0(edge, "-", edge + 9)
  ggsave(sprintf("output/figures/explore/matrices/P_U_%s.png", edge),
         make_combined(e$U$P, e$U$pi0, pi_star(e$U$P), levels = rel_level_order,
                       title_str = paste0("Unstandardized · cohort ", rng)),
         width = 10, height = 7, dpi = 200)
  if (!is.null(e$S)) {
    ggsave(sprintf("output/figures/explore/matrices/P_S_%s.png", edge),
           make_combined(e$S$P, e$S$pi0, pi_star(e$S$P), levels = rel_level_order,
                         title_str = paste0("Age-standardized 39–49 · cohort ", rng)),
           width = 10, height = 7, dpi = 200)
  }
}

# Faceted P-heatmap grid (cohort × method) — the "see the change" figure.
# Restrict to standardizable cohorts so U and S rows share the same columns.
sbins   = names(wide_out$mats)[sapply(wide_out$mats, function(e) !is.null(e$S))]
cells_df = function(mats, slot, lbl) do.call(rbind, lapply(names(mats), function(b) {
  e = mats[[b]]; if (is.null(e[[slot]])) return(NULL)
  df = as.data.frame(as.table(e[[slot]]$P)); names(df) = c("origin", "current", "p")
  df$cohort = paste0(as.integer(b), "-", as.integer(b) + 9); df$panel = lbl; df
}))
matsS   = wide_out$mats[sbins]
grid_df = rbind(cells_df(matsS, "U", "Unstandardized"),
                cells_df(matsS, "S", "Age-standardized (39–49)"))
grid_df$origin  = factor(grid_df$origin,  levels = rev(rel_level_order))
grid_df$current = factor(grid_df$current, levels = rel_level_order)
grid_df$panel   = factor(grid_df$panel, levels = c("Unstandardized", "Age-standardized (39–49)"))

p_grid = ggplot(grid_df, aes(current, origin, fill = p)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", p)), size = 2.6) +
  facet_grid(panel ~ cohort) +
  scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1), name = "P[i→j]") +
  labs(x = "Current religion (RELIG)", y = "Origin (RELIG16)",
       title = "RELIG16 → RELIG by cohort: unstandardized vs age-standardized (holding age fixed)",
       subtitle = "Watch the diagonal fade left→right; the age-standardized row shows the change is not a life-cycle artifact") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7), panel.grid = element_blank(),
        legend.position = "bottom", plot.subtitle = element_text(size = 9))
ggsave("output/figures/explore/matrices/P_grid_U_vs_S.png", p_grid,
       width = 16, height = 6, dpi = 200)

# ── HEADLINE COMPARISON: trend slope of λ₂ across cohort, by method ───────────
# A crude but legible summary: OLS slope of λ₂ on cohort-midpoint within each
# method (core standardizable cohorts only). If U's slope ≈ R's ≈ S's, the
# cohort trend is not an age artifact.

cat("\n================ λ₂ TREND SLOPE (per decade) BY METHOD ================\n")
slope_tab = do.call(rbind, lapply(list(res_wide, res_fixed), function(res) {
  do.call(rbind, lapply(split(res, res$method), function(m) {
    m = m[is.finite(m$lambda2), ]
    if (nrow(m) < 3) return(NULL)
    b = coef(lm(lambda2 ~ mid, data = m))[["mid"]] * 10   # per decade
    data.frame(band = m$band[1], method = m$method[1], n_bins = nrow(m),
               lambda2_slope_per_decade = round(b, 4))
  }))
}))
print(slope_tab, row.names = FALSE)

cat("\nDone. Figures in output/figures/explore/. Interpretation: compare the U line\n",
    "to R and S. If age-standardized (S) flattens or reverses the U trend, the\n",
    "cross-cohort memory finding is substantially a life-cycle artifact.\n")
