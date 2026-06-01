library(dplyr)
library(tidyr)
library(ggplot2)
library(gssr)

reltrad_labels = c(
  "1" = "evangelical", "2" = "mainline",  "3" = "black protestant",
  "4" = "catholic",    "5" = "jewish",    "6" = "other", "7" = "none")

data(gss_all)
data = gss_all |>
  select(year, cohort, reltrad, reltrad16) |>
  filter(!(year %in% c(1972, 2021))) |>
  mutate(across(c(reltrad, reltrad16),
                ~ reltrad_labels[as.character(as.numeric(.))])) |>
  filter(!is.na(reltrad), !is.na(reltrad16)) |>
  mutate(across(c(reltrad, reltrad16),
                ~ case_when(
                  . == "jewish"           ~ "other",
                  . == "black protestant" ~ "evangelical",
                  TRUE                    ~ .
                ),
                .names = "{.col}_alt")) |>
  mutate(age = year - cohort) |>
  filter(age > 30, cohort >= 1900) |>
  mutate(
    cohort_5  = floor((cohort - 1900) / 5)  * 5  + 1900,
    cohort_10 = floor((cohort - 1900) / 10) * 10 + 1900
  )

# ── FUNCTIONS ────────────────────────────────────────────────────────────────────

# Matrix power operator (replaces expm::`%^%`)
d`%^%` = function(M, k) {
  if (k == 0) return(diag(nrow(M)))
  result = M
  for (i in seq_len(k - 1)) result = result %*% M
  result
}

tv_norm = function(mu, nu) {
  0.5 * sum(abs(mu - nu))
}

# Initial distribution (π₀): unweighted share of each origin state.
pi_0 = function(data, origin) {
  tab = table(data[[origin]])
  v = as.numeric(tab / sum(tab))
  names(v) = names(tab)
  v
}

# Stationary distribution (π*): left eigenvector of P corresponding to eigenvalue 1.
pi_star = function(P) {
  eig = eigen(t(as.matrix(P)))
  v = Re(eig$vectors[, which.min(abs(eig$values - 1))])
  if (any(v < 0)) v = abs(v)
  setNames(v / sum(v), rownames(P))
}

# Unweighted row-stochastic transition matrix.
p_matrix = function(data, origin, current) {
  tab = table(data[[origin]], data[[current]])
  P = tab / rowSums(tab)
  class(P) = "matrix"
  P
}

# Individual memory (IM): log TV distance from π* for each origin state at step t.
im = function(data, origin, current, t = 1) {
  P_mat = p_matrix(data, origin, current)
  pi_s  = pi_star(P_mat)
  P_t   = P_mat %^% t
  im_i  = apply(P_t, 1, function(row_i) tv_norm(row_i, pi_s))
  log(im_i)
}

