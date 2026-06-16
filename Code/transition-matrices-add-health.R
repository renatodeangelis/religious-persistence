library(haven)
library(dplyr)
library(tidyr)
library(survey)
library(ggplot2)
library(forcats)
library(patchwork)

# ── PATHS ─────────────────────────────────────────────────────────────────────
# Set dir_data to the folder containing the Add Health .rds and .csv files

dir_data  = "add-health"

w1inhome = readRDS(file.path(dir_data, "w1inhome.rds")) |>
  select(AID, H1RE1, H1RE3, H1RE4, H1RE5, H1RE6, H1RE7, PA22, PA23, PA24, PA25, PA26) |>
  janitor::clean_names()
w4inhome = readRDS(file.path(dir_data, "w4inhome.rds")) |>
  select(AID, H4RE1, H4RE2, H4RE3, H4RE6, H4RE7, H4RE8, H4RE9, H4RE10) |>
  janitor::clean_names()

addhealth = inner_join(w1inhome, w4inhome, by = "aid") |>
  filter(!is.na(pa22), pa22 != 96, !(h4re1 %in% c(96, 98))) |>
  mutate(reltrad_parent = case_when(
    pa22 == 7 ~ "catholic",
    pa22 %in% c(6, 11, 14, 16, 18, 17, 19, 24, 27) ~ "other",
    pa22 == 28 ~ "none",
    pa22 %in% c(8, 9, 10, 12, 13, 20, 21, 26) ~ "mainline",
    pa22 %in% c(1, 2, 3, 5, 15, 23, 25) ~ "conservative"),
  reltrad_child_w4 = case_when(
    h4re1 == 1 ~ "none",
    h4re1 == 3 ~ "catholic",
    h4re1 %in% c(5, 6, 7, 8, 9) ~ "other",
    h4re6 %in% c(2, 9, 14, 15, 16, 21, 28, 30, 34, 36, 38, 39) ~ "mainline",
    h4re6 %in% c(3, 4, 5, 6, 7, 8, 10, 12, 13, 18, 19, 20, 22, 23, 24, 25, 26, 27, 31, 32, 33, 41, 42) ~ "conservative")) |>
  filter(!is.na(reltrad_child_w4))


# ── SCALE LABELS (parent Wave I benchmark) ────────────────────────────────────

att_lbl = c("1" = ">=once/wk",  "2" = ">=once/mo",
            "3" = "<once/mo",   "4" = "never")
sal_lbl = c("1" = "very important",      "2" = "fairly important",
            "3" = "fairly unimportant",  "4" = "not important at all")
pra_lbl = c("1" = ">=once/day", "2" = ">=once/wk",
            "3" = ">=once/mo",  "4" = "<once/mo", "5" = "never")
scr_lbl = c("1" = "agree", "2" = "disagree", "3" = "no sacred scriptures")

# ── HELPER FUNCTIONS ──────────────────────────────────────────────────────────

# Heatmap of a row-normalised transition matrix
make_hm = function(pmat, title, zmax, n = NULL, cell_size = 4.2, axis_size = 11) {
  subtitle = if (!is.null(n)) paste0("N = ", format(n, big.mark = ","), " complete dyads") else ""
  df = as.data.frame(pmat)
  colnames(df) = c("parent_state", "child_state", "pct")
  df$pct = df$pct * 100
  lvls_r = rownames(pmat); lvls_c = colnames(pmat)
  df = df |>
    mutate(
      parent_state = factor(parent_state, levels = lvls_r),
      child_state  = factor(child_state,  levels = lvls_c),
      text_col     = ifelse(pct > zmax * 0.55, "white", "grey15")
    )
  ggplot(df, aes(x = child_state, y = fct_rev(parent_state), fill = pct)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", pct), color = text_col),
              fontface = "bold", size = cell_size) +
    scale_color_identity() +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b",
                        limits = c(0, zmax), name = "Row %") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid  = element_blank(),
      axis.text.x = element_text(angle = 20, hjust = 1, size = axis_size, face = "bold"),
      axis.text.y = element_text(size = axis_size, face = "bold"),
      axis.title  = element_text(size = axis_size, face = "bold"),
      plot.title  = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(color = "grey40", size = 8)
    ) +
    labs(title = title, subtitle = subtitle, x = "Child State", y = "Parent State")
}

# Shared colour ceiling from a list of row-normalised matrices
compute_zmax = function(list_of_pmat) {
  max(unlist(lapply(list_of_pmat, as.numeric)), na.rm = TRUE) * 100
}

# Replace integer dimnames with string labels from a named lookup vector
relabel_tab = function(tab, row_lbl, col_lbl = row_lbl) {
  rn = rownames(tab); cn = colnames(tab)
  if (!is.null(row_lbl)) rownames(tab) = row_lbl[rn]
  if (!is.null(col_lbl)) colnames(tab) = col_lbl[cn]
  tab
}

# Bin birth years into cohort intervals of a given width
make_cohort_bins = function(birth_years, binwidth, start_year = NULL) {
  if (is.null(start_year)) start_year = min(birth_years, na.rm = TRUE)
  bin_idx   = floor((birth_years - start_year) / binwidth)
  bin_start = start_year + bin_idx * binwidth
  bin_end   = bin_start + binwidth - 1L
  ifelse(is.na(birth_years), NA_character_, paste0(bin_start, "-", bin_end))
}

# Survey-weighted transition matrices stratified by birth-cohort bin
cohort_svytable = function(df, wt_col, fmla,
                           complete_col = "complete_dyad",
                           lbl_vec      = NULL,
                           binwidth     = 5) {
  df2  = df |> mutate(cohort_bin = make_cohort_bins(birth_year, binwidth))
  des  = svydesign(ids     = ~CLUSTER2,
                   weights = as.formula(paste0("~", wt_col)),
                   data    = subset(df2, weight_valid),
                   nest    = TRUE)
  bins = df2 |>
    filter(.data[[complete_col]] == TRUE, !is.na(cohort_bin)) |>
    pull(cohort_bin) |> unique() |> sort()
  lapply(setNames(bins, bins), function(b) {
    cond  = parse(text = sprintf("%s == TRUE & cohort_bin == '%s'", complete_col, b))
    sub_d = subset(des, eval(cond))
    tab   = svytable(fmla, design = sub_d)
    if (!is.null(lbl_vec)) tab = relabel_tab(tab, lbl_vec)
    n = sum(df2[[complete_col]] == TRUE & !is.na(df2$cohort_bin) &
              df2$cohort_bin == b, na.rm = TRUE)
    list(P = prop.table(tab, 1), tab = tab, n = n)
  })
}

# Patchwork grid of per-cohort heatmaps
hm_cohort_grid = function(mats, section_title, zmax, cell_size = 4.2, axis_size = 11) {
  plots = lapply(names(mats), function(b) {
    make_hm(mats[[b]]$P, title = b, zmax = zmax, n = mats[[b]]$n,
            cell_size = cell_size, axis_size = axis_size) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  })
  patchwork::wrap_plots(plots, ncol = 1) +
    patchwork::plot_annotation(
      title = section_title,
      theme = theme(plot.title = element_text(face = "bold", size = 11))
    )
}

# ── TRANSITION MATRICES ───────────────────────────────────────────────────────
# Add matrix estimation code below, adapting state variable names to match
# whichever state coding is implemented above.

# ── 5×5 RELTRAD MATRIX SPLIT BY PARENT RELIGIOSITY (PA23–PA26) ───────────────
# PA23–PA26 are used as splitting variables, not states. For each variable,
# dichotomize and estimate the main 5×5 reltrad transition matrix separately
# for each subgroup. Unweighted for now.

ah_split = addhealth |>
  mutate(
    # Out-of-range codes (96/97/98) become NA via case_when
    pa23_d = case_when(pa23 %in% 1:2 ~ "high", pa23 %in% 3:4 ~ "low"),   # >= monthly vs. < monthly
    pa24_d = case_when(pa24 %in% 1:2 ~ "high", pa24 %in% 3:4 ~ "low"),   # important vs. unimportant
    pa25_d = case_when(pa25 %in% 1:2 ~ "high", pa25 %in% 3:5 ~ "low"),   # >= weekly vs. less frequent
    pa26_d = case_when(pa26 == 1L    ~ "high", pa26 %in% 2:3 ~ "low")    # agree vs. disagree/no scriptures
  )

reltrad_lvls = c("conservative", "mainline", "catholic", "other", "none")

# ── HELPER: unweighted 5×5 reltrad matrix for a data subset ──────────────────

reltrad_mat = function(data, label) {
  sub = data |>
    filter(!is.na(reltrad_parent), !is.na(reltrad_child_w4)) |>
    mutate(
      parent_f = factor(reltrad_parent,   levels = reltrad_lvls),
      child_f  = factor(reltrad_child_w4, levels = reltrad_lvls)
    )
  cat(sprintf("\n─── %s  (N = %d) ───\n", label, nrow(sub)))
  tab = table(parent = sub$parent_f, child = sub$child_f)
  if (sum(tab == 0) > 0) cat("NOTE:", sum(tab == 0), "zero cell(s)\n")
  cat("Row-normalised transition matrix:\n")
  print(round(prop.table(tab, 1), 3))
  cat("Cell counts:\n")
  print(tab)
  invisible(list(P = prop.table(tab, 1), tab = tab, n = nrow(sub)))
}

# ── PA23: service attendance ──────────────────────────────────────────────────

cat("\n\n======= PA23: Service attendance =======\n")
cat("Cutpoint: >= monthly (high) vs. < monthly (low)\n")
print(table(ah_split$pa23_d, useNA = "ifany"))

res_att_high = reltrad_mat(filter(ah_split, pa23_d == "high"), "PA23 high (>= monthly)")
res_att_low  = reltrad_mat(filter(ah_split, pa23_d == "low"),  "PA23 low  (< monthly)")

# ── PA24: importance of religion ──────────────────────────────────────────────

cat("\n\n======= PA24: Importance of religion =======\n")
cat("Cutpoint: important (high) vs. unimportant (low)\n")
print(table(ah_split$pa24_d, useNA = "ifany"))

res_sal_high = reltrad_mat(filter(ah_split, pa24_d == "high"), "PA24 high (important)")
res_sal_low  = reltrad_mat(filter(ah_split, pa24_d == "low"),  "PA24 low  (unimportant)")

# ── PA25: prayer frequency ────────────────────────────────────────────────────

cat("\n\n======= PA25: Prayer frequency =======\n")
cat("Cutpoint: >= weekly (high) vs. less frequent (low)\n")
print(table(ah_split$pa25_d, useNA = "ifany"))

res_pray_high = reltrad_mat(filter(ah_split, pa25_d == "high"), "PA25 high (>= weekly)")
res_pray_low  = reltrad_mat(filter(ah_split, pa25_d == "low"),  "PA25 low  (< weekly)")

# ── PA26: scripture inerrancy ─────────────────────────────────────────────────

cat("\n\n======= PA26: Scripture inerrancy =======\n")
cat("Cutpoint: agree (high) vs. disagree or no sacred scriptures (low)\n")
print(table(ah_split$pa26_d, useNA = "ifany"))

res_scr_high = reltrad_mat(filter(ah_split, pa26_d == "high"), "PA26 high (agree)")
res_scr_low  = reltrad_mat(filter(ah_split, pa26_d == "low"),  "PA26 low  (disagree / no scriptures)")
