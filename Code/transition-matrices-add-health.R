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
  filter(!is.na(pa22), pa22 != 96) |>
  mutate(reltrad_parent = case_when(
    pa22 == 7 ~ "catholic",
    pa22 %in% c(6, 11, 14, 16, 18, 17, 19, 24, 27) ~ "other",
    pa22 == 28 ~ "none",
    pa22 %in% c(8, 9, 10, 12, 13, 20, 21, 26) ~ "mainline",
    pa22 %in% c(1, 2, 3, 5, 15, 23, 25) ~ "conservative"),
  reltrad_child_w1 = case_when(
    h1re1 = 0 ~ "none",
    
  ))

# ── STATE SPACE ───────────────────────────────────────────────────────────────
# state3 = Christian / Non-Christian / None  (pre-coded in .rds files)
# state4 = Catholic / Protestant+other Christian / Non-Christian / None
# state5 = RELTRAD-style 5-category scheme (added below):
#   Mainline Protestant    — Disciples, Congregational, Episcopal, Friends,
#                            Lutheran, Methodist, Presbyterian, UCC,
#                            "other Protestant"
#   Evangelical Protestant — Baptist, Assemblies of God, Holiness, Pentecostal,
#                            Adventist, AME/CME, National Baptist
#   Catholic               — Catholic only
#   Other                  — non-Christian religions (Jewish, Muslim, Buddhist,
#                            Hindu, Baha'i, other religion) + ambiguous Christian
#                            groups (JW, Mormon, Eastern Orthodox, Christian
#                            Science, Unitarian)
#   None                   — no religion / atheist / agnostic
#
# Wave IV/V constraint: H4RE1/H5RE1 code 2 ("Protestant") cannot be split into
# Mainline vs. Evangelical without the H4RE6 denomination follow-up. Those
# respondents are coded NA in child_state5; use denomination_wave4 for a
# fine-grained Wave IV matrix.

# ── state4 on mapping table ───────────────────────────────────────────────────

mapping_affiliation = mapping_affiliation |>
  mutate(
    state4 = case_when(
      variable == "PA22"  & raw_code == 7L  ~ "Catholic",
      variable == "H1RE1" & raw_code == 22L ~ "Catholic",
      variable == "H4RE1" & raw_code == 3L  ~ "Catholic",
      variable == "H5RE1" & raw_code == 3L  ~ "Catholic",
      state3 == "Christian"                 ~ "Protestant / other Christian",
      !is.na(state3)                        ~ state3,
      TRUE                                  ~ NA_character_
    ),
    state5 = case_when(
      # Catholic
      variable == "PA22"  & raw_code == 7L  ~ "Catholic",
      variable == "H1RE1" & raw_code == 22L ~ "Catholic",
      variable == "H4RE1" & raw_code == 3L  ~ "Catholic",
      variable == "H5RE1" & raw_code == 3L  ~ "Catholic",
      # None
      variable == "PA22"  & raw_code == 28L ~ "None",
      variable == "H1RE1" & raw_code == 0L  ~ "None",
      variable %in% c("H4RE1", "H5RE1") & raw_code == 1L ~ "None",
      # Mainline Protestant (Wave I variables only — PA22 and H1RE1)
      variable == "PA22"  & raw_code %in% c(8L, 10L, 12L, 13L, 20L, 21L, 23L, 26L) ~ "Mainline Protestant",
      variable == "H1RE1" & raw_code %in% c(5L, 7L, 8L, 9L, 13L, 14L, 17L, 18L, 19L) ~ "Mainline Protestant",
      # Evangelical Protestant (Wave I variables only)
      variable == "PA22"  & raw_code %in% c(1L, 2L, 3L, 5L, 15L, 25L) ~ "Evangelical Protestant",
      variable == "H1RE1" & raw_code %in% c(1L, 2L, 3L, 4L, 10L, 15L, 16L) ~ "Evangelical Protestant",
      # Other: non-Christian + ambiguous Christian groups
      variable == "PA22"  & raw_code %in% c(6L, 9L, 11L, 14L, 16L, 17L, 18L, 19L, 24L, 27L) ~ "Other",
      variable == "H1RE1" & raw_code %in% c(6L, 11L, 12L, 20L, 21L, 23L, 24L, 25L, 26L, 27L, 28L) ~ "Other",
      # H4RE1/H5RE1 code 2 (Protestant) = NA — can't split without H4RE6
      variable %in% c("H4RE1", "H5RE1") & raw_code == 2L ~ NA_character_,
      variable == "H4RE1" & raw_code %in% c(4L, 5L, 6L, 7L, 8L, 9L) ~ "Other",
      variable == "H5RE1" & raw_code %in% c(4L, 5L, 6L, 9L)          ~ "Other",
      TRUE ~ NA_character_
    )
  )

# ── state4 on affiliation dataframes ─────────────────────────────────────────

affiliation_wave1 = affiliation_wave1 |>
  mutate(
    parent_state4 = case_when(
      PA22  == 7L                   ~ "Catholic",
      parent_state3 == "Christian"  ~ "Protestant / other Christian",
      !is.na(parent_state3)         ~ parent_state3,
      TRUE                          ~ NA_character_
    ),
    child_state4 = case_when(
      H1RE1 == 22L                  ~ "Catholic",
      child_state3 == "Christian"   ~ "Protestant / other Christian",
      !is.na(child_state3)          ~ child_state3,
      TRUE                          ~ NA_character_
    ),
    parent_state5 = case_when(
      PA22 == 7L                                          ~ "Catholic",
      PA22 == 28L                                         ~ "None",
      PA22 %in% c(8L, 10L, 12L, 13L, 20L, 21L, 23L, 26L) ~ "Mainline Protestant",
      PA22 %in% c(1L, 2L, 3L, 5L, 15L, 25L)              ~ "Evangelical Protestant",
      PA22 %in% c(6L, 9L, 11L, 14L, 16L, 17L,
                  18L, 19L, 24L, 27L)                    ~ "Other",
      TRUE                                               ~ NA_character_
    ),
    child_state5 = case_when(
      H1RE1 == 22L                                        ~ "Catholic",
      H1RE1 == 0L                                         ~ "None",
      # "other Protestant": born-again (H1RE5 == 1) → Evangelical; else → Mainline
      H1RE1 == 19L & H1RE5 == 1L                          ~ "Evangelical Protestant",
      H1RE1 == 19L                                        ~ "Mainline Protestant",
      H1RE1 %in% c(5L, 7L, 8L, 9L, 13L, 14L, 17L, 18L)  ~ "Mainline Protestant",
      H1RE1 %in% c(1L, 2L, 3L, 4L, 10L, 15L, 16L)        ~ "Evangelical Protestant",
      H1RE1 %in% c(6L, 11L, 12L, 20L, 21L, 23L,
                   24L, 25L, 26L, 27L, 28L)              ~ "Other",
      TRUE                                               ~ NA_character_
    )
  )

affiliation_wave4 = affiliation_wave4 |>
  mutate(
    parent_state4 = case_when(
      PA22  == 7L                   ~ "Catholic",
      parent_state3 == "Christian"  ~ "Protestant / other Christian",
      !is.na(parent_state3)         ~ parent_state3,
      TRUE                          ~ NA_character_
    ),
    child_state4 = case_when(
      H4RE1 == 3L                   ~ "Catholic",
      child_state3 == "Christian"   ~ "Protestant / other Christian",
      !is.na(child_state3)          ~ child_state3,
      TRUE                          ~ NA_character_
    ),
    parent_state5 = case_when(
      PA22 == 7L                                          ~ "Catholic",
      PA22 == 28L                                         ~ "None",
      PA22 %in% c(8L, 10L, 12L, 13L, 20L, 21L, 23L, 26L) ~ "Mainline Protestant",
      PA22 %in% c(1L, 2L, 3L, 5L, 15L, 25L)              ~ "Evangelical Protestant",
      PA22 %in% c(6L, 9L, 11L, 14L, 16L, 17L,
                  18L, 19L, 24L, 27L)                    ~ "Other",
      TRUE                                               ~ NA_character_
    ),
    # H4RE1 code 2 (Protestant) → NA: cannot split Mainline vs. Evangelical
    child_state5 = case_when(
      H4RE1 == 1L                              ~ "None",
      H4RE1 == 3L                              ~ "Catholic",
      H4RE1 %in% c(4L, 5L, 6L, 7L, 8L, 9L)   ~ "Other",
      TRUE                                    ~ NA_character_
    )
  )

affiliation_wave5 = affiliation_wave5 |>
  mutate(
    parent_state4 = case_when(
      PA22  == 7L                   ~ "Catholic",
      parent_state3 == "Christian"  ~ "Protestant / other Christian",
      !is.na(parent_state3)         ~ parent_state3,
      TRUE                          ~ NA_character_
    ),
    child_state4 = case_when(
      H5RE1 == 3L                   ~ "Catholic",
      child_state3 == "Christian"   ~ "Protestant / other Christian",
      !is.na(child_state3)          ~ child_state3,
      TRUE                          ~ NA_character_
    ),
    parent_state5 = case_when(
      PA22 == 7L                                          ~ "Catholic",
      PA22 == 28L                                         ~ "None",
      PA22 %in% c(8L, 10L, 12L, 13L, 20L, 21L, 23L, 26L) ~ "Mainline Protestant",
      PA22 %in% c(1L, 2L, 3L, 5L, 15L, 25L)              ~ "Evangelical Protestant",
      PA22 %in% c(6L, 9L, 11L, 14L, 16L, 17L,
                  18L, 19L, 24L, 27L)                    ~ "Other",
      TRUE                                               ~ NA_character_
    ),
    # H5RE1 code 2 (Protestant) → NA: cannot split Mainline vs. Evangelical
    child_state5 = case_when(
      H5RE1 == 1L                        ~ "None",
      H5RE1 == 3L                        ~ "Catholic",
      H5RE1 %in% c(4L, 5L, 6L, 9L)      ~ "Other",
      TRUE                              ~ NA_character_
    )
  )

# ── COMPLETE-CASE INDICATORS: AFFILIATION ─────────────────────────────────────

affiliation_wave1 = affiliation_wave1 |>
  mutate(
    weight_valid     = !is.na(GSWGT1)   & GSWGT1   > 0,
    complete_dyad_s3 = weight_valid & !is.na(parent_state3) & !is.na(child_state3),
    complete_dyad_s4 = weight_valid & !is.na(parent_state4) & !is.na(child_state4),
    complete_dyad_s5 = weight_valid & !is.na(parent_state5) & !is.na(child_state5)
  )
affiliation_wave4 = affiliation_wave4 |>
  mutate(
    weight_valid     = !is.na(GSWGT4_2) & GSWGT4_2 > 0,
    complete_dyad_s3 = weight_valid & !is.na(parent_state3) & !is.na(child_state3),
    complete_dyad_s4 = weight_valid & !is.na(parent_state4) & !is.na(child_state4),
    complete_dyad_s5 = weight_valid & !is.na(parent_state5) & !is.na(child_state5)
  )
affiliation_wave5 = affiliation_wave5 |>
  mutate(
    weight_valid     = !is.na(GSW5)     & GSW5     > 0,
    complete_dyad_s3 = weight_valid & !is.na(parent_state3) & !is.na(child_state3),
    complete_dyad_s4 = weight_valid & !is.na(parent_state4) & !is.na(child_state4),
    complete_dyad_s5 = weight_valid & !is.na(parent_state5) & !is.na(child_state5)
  )

# ── COMPLETE-CASE INDICATORS: DENOMINATION ───────────────────────────────────

denomination_wave1 = denomination_wave1 |>
  mutate(
    weight_valid  = !is.na(GSWGT1)   & GSWGT1   > 0,
    complete_dyad = weight_valid & !is.na(parent_state7) & !is.na(child_state7)
  )
denomination_wave4 = denomination_wave4 |>
  mutate(
    weight_valid  = !is.na(GSWGT4_2) & GSWGT4_2 > 0,
    complete_dyad = weight_valid & !is.na(parent_state7) & !is.na(child_state7)
  )

# ── COMPLETE-CASE INDICATORS: PRACTICE / BELIEF ───────────────────────────────

attendance_wave1 = attendance_wave1 |>
  mutate(weight_valid  = !is.na(GSWGT1)   & GSWGT1   > 0,
         complete_dyad = weight_valid & !is.na(PA23_harm) & !is.na(H1RE3_harm))
attendance_wave4 = attendance_wave4 |>
  mutate(weight_valid  = !is.na(GSWGT4_2) & GSWGT4_2 > 0,
         complete_dyad = weight_valid & !is.na(PA23_harm) & !is.na(H4RE7_harm))
attendance_wave5 = attendance_wave5 |>
  mutate(weight_valid  = !is.na(GSW5)     & GSW5     > 0,
         complete_dyad = weight_valid & !is.na(PA23_harm) & !is.na(H5RE2_harm))

salience_wave1 = salience_wave1 |>
  mutate(weight_valid  = !is.na(GSWGT1)   & GSWGT1   > 0,
         complete_dyad = weight_valid & !is.na(PA24_harm) & !is.na(H1RE4_harm))
salience_wave4 = salience_wave4 |>
  mutate(weight_valid  = !is.na(GSWGT4_2) & GSWGT4_2 > 0,
         complete_dyad = weight_valid & !is.na(PA24_harm) & !is.na(H4RE9_harm))
salience_wave5 = salience_wave5 |>
  mutate(weight_valid  = !is.na(GSW5)     & GSW5     > 0,
         complete_dyad = weight_valid & !is.na(PA24_harm) & !is.na(H5RE3_harm))

prayer_wave1 = prayer_wave1 |>
  mutate(weight_valid  = !is.na(GSWGT1)   & GSWGT1   > 0,
         complete_dyad = weight_valid & !is.na(PA25_harm) & !is.na(H1RE6_harm))
prayer_wave4 = prayer_wave4 |>
  mutate(weight_valid  = !is.na(GSWGT4_2) & GSWGT4_2 > 0,
         complete_dyad = weight_valid & !is.na(PA25_harm) & !is.na(H4RE10_harm))
prayer_wave5 = prayer_wave5 |>
  mutate(weight_valid  = !is.na(GSW5)     & GSW5     > 0,
         complete_dyad = weight_valid & !is.na(PA25_harm) & !is.na(H5RE4_harm))

scripture_wave1 = scripture_wave1 |>
  mutate(weight_valid  = !is.na(GSWGT1)   & GSWGT1   > 0,
         complete_dyad = weight_valid & !is.na(PA26_harm) & !is.na(H1RE2_harm))

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
