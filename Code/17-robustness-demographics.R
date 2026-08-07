# ── 16 · ROBUSTNESS: DEMOGRAPHIC STRATIFICATION (SIBLINGS & CHILDREN) ─────────
# Decadal transition matrices stratified by family-size variables:
#   sibs_group  — siblings of origin: "few" (0-1), "mid" (2-4), "many" (5+)
#   childs_group — respondent's own children: "none" (0), "small" (1-2), "large" (3+)
#
# NOTE: childs_group is partly endogenous to current religion (more-religious
# respondents tend to have more children). Treat that stratification as
# descriptive rather than causal; sibs_group is the cleaner moderator.
#
# Five cohort windows: 1940-1949 through 1980-1989.
#
# Input:  data/derived/gss_clean.rds
# Output: data/derived/matrices_demographics.rds
#         output/figures/demographics/*.png

source("code/utils.R")

clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_alt = clean$states_alt

rel_level_order = c("catholic", "evangelical", "mainline", "other", "none")
reltrad_colors  = c(
  catholic    = "#0072B2",
  evangelical = "#D55E00",
  mainline    = "#009E73",
  other       = "#CC79A7",
  none        = "#999999"
)
reltrad_labels_tc = c(
  catholic = "Catholic", evangelical = "Evangelical", mainline = "Mainline",
  other = "Other", none = "None"
)

demo_vars = list(
  sibs_group  = c("few", "mid", "many"),
  childs_group = c("none", "small", "large")
)

# Readable labels for figure titles
demo_group_labels = list(
  sibs_group  = c(few = "0-1 siblings", mid = "2-4 siblings", many = "5+ siblings"),
  childs_group = c(none = "No children", small = "1-2 children", large = "3+ children")
)

demo_var_titles = c(
  sibs_group   = "Family of Origin Size (Siblings)",
  childs_group = "Own Family Size (Children)"
)

# 10-year bin midpoints (edges 1940-1980)
mids_demo = c(1940, 1950, 1960, 1970, 1980)

# ── BUILD MATRICES ────────────────────────────────────────────────────────────

P_demo = pi0_demo = pistar_demo = n_demo = list()

for (vname in names(demo_vars)) {
  for (grp in demo_vars[[vname]]) {
    for (mid in mids_demo) {
      sub = data[!is.na(data$cohort_10)     & data$cohort_10 == mid  &
                 !is.na(data[[vname]])      & data[[vname]]  == grp  &
                 !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]
      if (nrow(sub) < 30) next
      key = paste(vname, grp, mid, sep = "_")
      P_demo[[key]]      = p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
      pi0_demo[[key]]    = pi_0(sub, "reltrad16_alt")
      pistar_demo[[key]] = pi_star(P_demo[[key]])
      n_demo[[key]]      = nrow(sub)
    }
  }
}

# ── CONSOLE OUTPUT ───────────────────────────────────────────────────────────

for (key in names(P_demo)) {
  edge = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
  lbl  = gsub("_", " ", sub("_\\d{4}$", "", key))
  cat("\n──", lbl, "| Cohort", edge, "-", edge + 9,
      " (N =", n_demo[[key]], ") ──\n")
  print(round(P_demo[[key]], 3))
}

# ── FIGURES ──────────────────────────────────────────────────────────────────

dir.create("output/figures/demographics", recursive = TRUE, showWarnings = FALSE)

# Individual heatmaps
for (key in names(P_demo)) {
  edge = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
  # Extract variable name and group: key = "<vname>_<grp>_<mid>"
  # grp may contain underscores only if we built it that way — it doesn't here
  mid_str = sub(".*_(\\d{4})$", "\\1", key)
  stem    = sub(paste0("_", mid_str, "$"), "", key)
  vname   = names(demo_vars)[sapply(names(demo_vars), function(v) startsWith(stem, v))]
  grp     = sub(paste0("^", vname, "_"), "", stem)
  grp_lbl = demo_group_labels[[vname]][[grp]]
  p = make_combined(
    P_demo[[key]], pi0_demo[[key]], pistar_demo[[key]],
    levels    = rel_level_order,
    title_str = paste0(grp_lbl, " – Cohort ", edge, "-", edge + 9,
                       "  (N = ", n_demo[[key]], ")")
  )
  ggsave(paste0("output/figures/demographics/trans_", key, "_10yr.png"),
         p, width = 10, height = 7, dpi = 200)
}

# Diagonal persistence — one figure per variable (faceted by group)
for (vname in names(demo_vars)) {
  keys_v = grep(paste0("^", vname, "_"), names(P_demo), value = TRUE)
  if (length(keys_v) == 0) next

  diag_v = do.call(rbind, lapply(keys_v, function(key) {
    edge    = as.integer(sub(".*_(\\d{4})$", "\\1", key)) - 5
    mid_str = sub(".*_(\\d{4})$", "\\1", key)
    stem    = sub(paste0("_", mid_str, "$"), "", key)
    grp     = sub(paste0("^", vname, "_"), "", stem)
    grp_lbl = demo_group_labels[[vname]][[grp]]
    data.frame(cohort  = edge,
               group   = grp_lbl,
               origin  = rel_level_order,
               persist = diag(P_demo[[key]])[rel_level_order],
               row.names = NULL)
  }))
  diag_v$origin = factor(diag_v$origin, levels = rel_level_order)
  # Order groups sensibly
  diag_v$group = factor(diag_v$group,
                        levels = unname(demo_group_labels[[vname]]))

  p_diag = ggplot(diag_v, aes(x = cohort, y = persist, color = origin, group = origin)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    facet_wrap(~ group) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
    labs(x = "Birth cohort",
         y = "Probability to Stay",
         color = NULL,
         title = paste0("Retention by ", demo_var_titles[[vname]], " and Birth Cohort")) +
    healy_theme

  ggsave(paste0("output/figures/demographics/diagonal_persistence_", vname, ".png"),
         p_diag, width = 12, height = 5, dpi = 200)
}

# ── ORIGIN RELIGION BY SIBLINGS × COHORT ─────────────────────────────────────

origin_sibs = data[
  !is.na(data$cohort_10) & data$cohort_10 %in% mids_demo &
  !is.na(data$sibs_group) & !is.na(data$reltrad16_alt), ,
  drop = FALSE
]

origin_sibs$sibs_label = factor(
  origin_sibs$sibs_group,
  levels = c("few", "mid", "many"),
  labels = c("0–1 siblings", "2–4 siblings", "5+ siblings")
)

origin_counts = as.data.frame(table(
  cohort_10     = origin_sibs$cohort_10,
  sibs_label    = origin_sibs$sibs_label,
  reltrad16_alt = origin_sibs$reltrad16_alt
))
names(origin_counts)[names(origin_counts) == "Freq"] = "n"
origin_counts = origin_counts[origin_counts$n > 0, ]
origin_counts$cohort_10 = as.integer(as.character(origin_counts$cohort_10))

origin_tots = aggregate(n ~ cohort_10 + sibs_label, data = origin_counts, FUN = sum)
names(origin_tots)[3] = "total"
origin_counts = merge(origin_counts, origin_tots, by = c("cohort_10", "sibs_label"))
origin_counts$prop = origin_counts$n / origin_counts$total

origin_counts$origin = factor(origin_counts$reltrad16_alt, levels = rel_level_order)

p_origin_sibs = ggplot(
  origin_counts,
  aes(x = cohort_10, y = prop, color = origin, group = origin)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ sibs_label, nrow = 1) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_x_continuous(breaks = mids_demo,
                     labels = paste0(mids_demo, "–", mids_demo + 9)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, NA)) +
  labs(
    x     = "Birth cohort",
    y     = "Share of origin-religion sample",
    color = NULL,
    title = "Origin Religion by Siblings of Origin and Birth Cohort"
  ) +
  healy_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave("output/figures/demographics/origin_religion_by_sibs_cohort.png",
       p_origin_sibs, width = 13, height = 5, dpi = 200)

# ── CHILDREN BY BIRTH YEAR × ORIGIN RELIGION ─────────────────────────────────

childs_yr = data[
  !is.na(data$cohort) & data$cohort >= 1925 & data$cohort <= 1984 &
  !is.na(data$childs) & !is.na(data$reltrad_alt), ,
  drop = FALSE
]
childs_yr$childs  = as.numeric(childs_yr$childs)
childs_yr$origin  = factor(childs_yr$reltrad_alt, levels = rel_level_order)

p_childs = ggplot(childs_yr, aes(x = cohort, y = childs, color = origin)) +
  stat_summary(geom = "point", fun = mean, size = 1.2, alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, span = 0.4) +
  scale_color_manual(values = reltrad_colors, labels = reltrad_labels_tc) +
  scale_x_continuous(breaks = seq(1925, 1985, 10)) +
  labs(
    x     = "Birth year",
    y     = "Number of children",
    color = NULL,
    title = "Number of Children by Birth Year and Origin Religion"
  ) +
  healy_theme

ggsave("output/figures/demographics/childs_by_birthyear_religion.png",
       p_childs, width = 10, height = 6, dpi = 200)

# ── CHILDREN BY BIRTH YEAR — WOMEN WHO SWITCHED TO EVANGELICAL OR NONE ───────

sw_colors   = c(evangelical = reltrad_colors[["evangelical"]],
                none        = reltrad_colors[["none"]])
sw_linetypes = c(switcher = "solid", stayer = "dashed")

switchers = data |>
  filter(
    as.numeric(sex) == 2,
    between(cohort, 1925, 1984),
    !is.na(childs), !is.na(reltrad_alt), !is.na(reltrad16_alt),
    reltrad_alt %in% c("evangelical", "none")
  ) |>
  mutate(
    rel_dest   = reltrad_alt,
    trajectory = if_else(reltrad_alt == reltrad16_alt, "stayer", "switcher")
  )

p_switchers = ggplot(switchers,
                     aes(x = cohort, y = childs,
                         color = rel_dest, linetype = trajectory,
                         group = interaction(rel_dest, trajectory))) +
  stat_summary(geom = "point", fun = mean, size = 1.2, alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.9, span = 0.4) +
  scale_color_manual(values    = sw_colors,
                     labels    = c(evangelical = "Evangelical", none = "None")) +
  scale_linetype_manual(values = sw_linetypes,
                        labels = c(switcher = "Switcher", stayer = "Stayer")) +
  scale_x_continuous(breaks = seq(1925, 1985, 10)) +
  scale_y_continuous(limits = c(0, 4)) +
  labs(
    x        = "Birth year",
    y        = "Number of children",
    color    = "Destination",
    linetype = NULL,
    title    = "Number of Children by Birth Year — Women Switching or Staying"
  ) +
  healy_theme

ggsave("output/figures/demographics/childs_by_birthyear_switchers.png",
       p_switchers, width = 10, height = 6, dpi = 200)

saveRDS(
  list(P = P_demo, pi0 = pi0_demo, pistar = pistar_demo, n = n_demo),
  "data/derived/matrices_demographics.rds"
)
cat("Wrote output/figures/demographics/\n")
cat("Wrote data/derived/matrices_demographics.rds\n")
