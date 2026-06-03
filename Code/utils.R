library(ggplot2)
library(patchwork)

# ── MATRIX MATH ──────────────────────────────────────────────────────────────

# Matrix power operator (replaces expm::`%^%`)
`%^%` = function(M, k) {
  if (k == 0) {
    I = diag(nrow(M))
    dimnames(I) = dimnames(M)
    return(I)
  }
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
# levels: optional character vector to fix state ordering across cohort subsets.
p_matrix = function(data, origin, current, levels = NULL) {
  f   = if (!is.null(levels)) function(x) factor(x, levels = levels) else factor
  tab = table(f(data[[origin]]), f(data[[current]]))

  # Warn on true zero cells (sampling or structural)
  zero_idx = which(tab == 0, arr.ind = TRUE)
  if (nrow(zero_idx) > 0) {
    zero_labels = paste(rownames(tab)[zero_idx[, 1]], "->",
                        colnames(tab)[zero_idx[, 2]], collapse = "; ")
    warning("Zero cells in P: ", zero_labels)
  }

  P = tab / rowSums(tab)
  class(P) = "matrix"
  P
}

# ── MEMORY MEASURES ──────────────────────────────────────────────────────────

# Individual memory (IM): log TV distance from π* for each origin state at step t.
im = function(data, origin, current, t = 1) {
  P_mat = p_matrix(data, origin, current)
  pi_s  = pi_star(P_mat)
  P_t   = P_mat %^% t
  im_i  = apply(P_t, 1, function(row_i) tv_norm(row_i, pi_s))
  log(im_i)
}

# IM from a pre-built matrix (reuses P_list_* to avoid recomputing P).
im_from_P = function(P, t = 1) {
  pi_s = pi_star(P)
  P_t  = P %^% t
  im_i = apply(P_t, 1, function(row_i) tv_norm(row_i, pi_s))
  setNames(log(im_i), rownames(P))
}

# π₀-weighted mean probability of leaving origin state; defaults to uniform π₀
overall_mobility = function(P, pi0 = NULL) {
  if (is.null(pi0)) pi0 = rep(1 / nrow(P), nrow(P))
  1 - sum(pi0 * diag(P))
}

# Marginal distribution of origin states at generation t.
mu_t = function(pi0, P, t = 0) {
  if (t == 0) return(as.numeric(pi0))
  as.numeric(pi0 %*% (P %^% t))
}

# Structural mobility: TV distance between the marginal at t and t+1.
sm = function(P, pi0, t = 0) {
  mu  = mu_t(pi0, P, t)
  mu1 = as.numeric(mu %*% P)
  tv_norm(mu, mu1)
}

# Exchange mobility: overall mobility minus structural component.
em = function(P, pi0, t = 0) {
  overall_mobility(P, mu_t(pi0, P, t)) - sm(P, pi0, t)
}

# ── PLOTTING HELPERS ─────────────────────────────────────────────────────────

plot_pmat_heatmap = function(P, levels = NULL, text_size = 5, title_str = "P") {
  df = as.data.frame(as.table(P))
  names(df) = c("origin", "current", "est")
  if (!is.null(levels)) {
    df$origin  = factor(df$origin,  levels = levels)
    df$current = factor(df$current, levels = rev(levels))
  }
  ggplot(df, aes(x = current, y = origin, fill = est)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", est)), size = text_size) +
    scale_fill_gradient(low = "lightyellow", high = "firebrick") +
    labs(x = "Current religion", y = "Origin religion",
         fill = "Transition Prob.", title = title_str) +
    theme_minimal() +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 13),
          axis.text.y     = element_text(angle = 45, hjust = 1, size = 13),
          axis.ticks.x    = element_blank(),
          axis.ticks.y    = element_blank(),
          axis.title      = element_text(size = 13),
          legend.position = "bottom",
          legend.text     = element_text(size = 11),
          legend.title    = element_text(size = 12),
          plot.title      = element_text(hjust = 0.5, size = 16),
          panel.grid      = element_blank())
}

plot_pi_column = function(vec, title_str, levels = NULL) {
  df = data.frame(origin = names(vec), value = as.numeric(vec))
  if (!is.null(levels)) {
    df$origin = factor(df$origin, levels = levels)
  } else {
    df$origin = factor(df$origin, levels = rev(unique(df$origin)))
  }
  ggplot(df, aes(x = 1, y = origin, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", value)), size = 5, color = "black") +
    scale_fill_gradient(low = "lightyellow", high = "firebrick") +
    labs(title = title_str) +
    theme_minimal() +
    theme(axis.title.x = element_blank(), axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(), axis.title.y = element_blank(),
          axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "none",
          plot.title   = element_text(hjust = 0.5, size = 16),
          panel.grid   = element_blank())
}

make_combined = function(P, pi0, pistar, levels = NULL, title_str = "P") {
  g      = plot_pmat_heatmap(P,      levels = levels, title_str = title_str)
  g0     = plot_pi_column(pi0,    title_str = "π₀", levels = levels)
  g_star = plot_pi_column(pistar, title_str = "π*",      levels = levels)
  patchwork::wrap_plots(g, g0, g_star, widths = c(6, 1, 1))
}