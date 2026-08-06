# ── 06 · λ₂ TREND WITH BOOTSTRAPPED CIs ──────────────────────────────────────
# Computes λ₂ (modulus of second-largest eigenvalue of P_t) for each 10-year
# birth cohort window and produces 95% bootstrap CIs (250 replicates, seed 42).
#
# Input:  data/derived/gss_clean.rds
# Output: output/figures/lambda2_trend_10yr.png

library(dplyr)
library(ggplot2)
source("code/utils.R")

clean      = readRDS("data/derived/gss_clean.rds")
data       = clean$data
states_alt = clean$states_alt

lambda2 = function(P) sort(Mod(eigen(P)$values), decreasing = TRUE)[2]

mids   = c(1930, 1940, 1950, 1960, 1970, 1980)
n_boot = 250
set.seed(42)

# ── BOOTSTRAP ────────────────────────────────────────────────────────────────

boot_rows = lapply(mids, function(mid) {
  sub = data[!is.na(data$cohort_10) & data$cohort_10 == mid &
             !is.na(data$reltrad16_alt) & !is.na(data$reltrad_alt), ]

  P_obs  = suppressWarnings(
    p_matrix(sub, "reltrad16_alt", "reltrad_alt", levels = states_alt)
  )
  l2_obs = lambda2(P_obs)

  l2_boot = replicate(n_boot, {
    idx = sample(nrow(sub), replace = TRUE)
    P_b = suppressWarnings(
      p_matrix(sub[idx, ], "reltrad16_alt", "reltrad_alt", levels = states_alt)
    )
    lambda2(P_b)
  })

  data.frame(
    mid    = mid,
    cohort = paste0(mid - 5, "–", sprintf("%02d", (mid + 4) %% 100)),
    l2     = l2_obs,
    lo     = quantile(l2_boot, 0.025),
    hi     = quantile(l2_boot, 0.975),
    n      = nrow(sub),
    row.names = NULL
  )
})

boot_df = do.call(rbind, boot_rows)
boot_df$cohort = factor(boot_df$cohort, levels = boot_df$cohort)

cat("\nλ₂ estimates with 95% bootstrap CIs (n_boot = 250, seed = 42):\n")
print(boot_df[, c("cohort", "l2", "lo", "hi", "n")], row.names = FALSE, digits = 3)

# ── FIGURE ───────────────────────────────────────────────────────────────────

p_l2 = ggplot(boot_df, aes(x = cohort, y = l2, group = 1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#0072B2", alpha = 0.15) +
  geom_line(color = "#0072B2", linewidth = 0.9) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  color     = "#0072B2",
                  linewidth = 0.7,
                  size      = 0.6) +
  scale_y_continuous(limits = c(0.5, 0.9), breaks = seq(0.5, 0.9, 0.1)) +
  labs(
    x       = "Birth cohort (10-year window)",
    y       = expression(lambda[2]),
    title   = expression(lambda[2] ~ "by birth cohort with 95% bootstrap CIs")) +
  healy_theme

ggsave("output/figures/lambda2_trend_10yr.png", p_l2,
       width = 8, height = 5, dpi = 200)

cat("\nSaved: output/figures/lambda2_trend_10yr.png\n")
