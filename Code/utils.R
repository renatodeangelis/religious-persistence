library(ggplot2)
library(patchwork)

# ── SHARED CONSTANTS ──────────────────────────────────────────────────────────
# Config used across the split pipeline (01–12). Data-derived constants such as
# `states_alt` live in gss_clean.rds instead, because they depend on the data.

# Column order for the 5-state scheme in every figure and homogeneity test.
rel_level_order = c("catholic", "evangelical", "mainline", "other", "none")

# Shared Healy theme (theme_bw + Okabe-Ito palette conventions)
healy_theme = theme_bw(base_size = 12) +
  theme(
    panel.grid.minor  = element_blank(),
    legend.position   = "bottom",
    legend.title      = element_text(size = 12),
    legend.text       = element_text(size = 11),
    axis.title        = element_text(size = 12),
    plot.title        = element_text(size = 13, face = "plain"),
    strip.background  = element_rect(fill = "grey92", color = NA),
    strip.text        = element_text(size = 11)
  )

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

# Mean time to exit state i: geometric sojourn with exit prob (1 - P_ii).
# Closed form: MTE_i = 1 / (1 - P_ii). Returns Inf for absorbing states (P_ii = 1).
mte = function(P) {
  d = diag(as.matrix(P))
  setNames(1 / (1 - d), rownames(P))
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
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", est)), size = text_size) +
    scale_fill_distiller(palette = "Blues", direction = 1,
                         limits = c(0, 1), name = "Prob.") +
    labs(x = "Current religion", y = "Origin religion", title = title_str) +
    theme_bw(base_size = 12) +
    theme(axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
          axis.text.y     = element_text(size = 11),
          axis.ticks      = element_blank(),
          axis.title      = element_text(size = 12),
          legend.position = "bottom",
          legend.text     = element_text(size = 10),
          legend.title    = element_text(size = 11),
          plot.title      = element_text(hjust = 0.5, size = 14, face = "plain"),
          panel.grid      = element_blank(),
          panel.border    = element_blank())
}

plot_pi_column = function(vec, title_str, levels = NULL, text_size = 5) {
  df = data.frame(origin = names(vec), value = as.numeric(vec))
  if (!is.null(levels)) {
    df$origin = factor(df$origin, levels = levels)
  } else {
    df$origin = factor(df$origin, levels = rev(unique(df$origin)))
  }
  ggplot(df, aes(x = 1, y = origin, fill = value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", value)), size = text_size, color = "black") +
    scale_fill_distiller(palette = "Blues", direction = 1, limits = c(0, 1)) +
    labs(title = title_str) +
    theme_bw(base_size = 12) +
    theme(axis.title.x   = element_blank(), axis.text.x  = element_blank(),
          axis.ticks.x   = element_blank(), axis.title.y = element_blank(),
          axis.text.y    = element_blank(), axis.ticks.y = element_blank(),
          legend.position = "none",
          plot.title     = element_text(hjust = 0.5, size = 13, face = "plain"),
          panel.grid     = element_blank(),
          panel.border   = element_blank())
}

make_combined = function(P, pi0, pistar, levels = NULL, title_str = "P", text_size = 5) {
  g      = plot_pmat_heatmap(P,      levels = levels, title_str = title_str, text_size = text_size)
  g0     = plot_pi_column(pi0,    title_str = "π₀", levels = levels, text_size = text_size)
  g_star = plot_pi_column(pistar, title_str = "π*", levels = levels, text_size = text_size)
  patchwork::wrap_plots(g, g0, g_star, widths = c(6, 1, 1))
}

count_matrix = function(data, origin, current, levels = NULL) {
  f = if (!is.null(levels)) function(x) factor(x, levels = levels) else factor
  tab = table(f(data[[origin]]), f(data[[current]]))
  class(tab) = "matrix"
  tab
}

calculate_phat_ij = function(count_list) {
  # Step 1: element-wise sum across all T matrices  
  N_pool = 
    Reduce(
      "+", 
      lapply(
        count_list, 
        function(x) 
          matrix(as.numeric(x),
                 nrow = nrow(x), 
                 ncol = ncol(x))
      )
    )
  # Step 2: pooled row totals  
  n_star_i = rowSums(N_pool)
  # Step 3: divide each row by its pooled row total
  p_hat_pool = N_pool / n_star_i   
  # Preserve state names if present
  dimnames(p_hat_pool) = dimnames(count_list[[1]])
  return(p_hat_pool)
}

calculate_phat_ijt = function(count_lst){
  phat_ijt =
    lapply(
      count_lst,
      function(mat) {
        row_sums = rowSums(mat)
        mat / row_sums
      }
    )
}

chi2_row = function(count_list, alpha = 0.05) {
  T_steps = length(count_list)
  p_hat_pool = calculate_phat_ij(count_list)
  p_hat_each = calculate_phat_ijt(count_list)
  states  = rownames(p_hat_pool)
  m       = nrow(p_hat_pool)
  chi2_i = numeric(m)
  mats = count_list
  
  for (i in seq_len(m)) {
    p0_i   = p_hat_pool[i, ]          
    active = p0_i > 0                  
    for (t in seq_len(T_steps)) {
      n_i_t  = rowSums(mats[[t]])[i]   
      if (n_i_t == 0) next
      pt_i   = p_hat_each[[t]][i, ]     
      chi2_i[i] = chi2_i[i] +
        n_i_t * sum((pt_i[active] - p0_i[active])^2 / p0_i[active])
    }
  }
  
  df_i = sapply(seq_len(m), function(i) {
    n_zero = sum(p_hat_pool[i, ] == 0)
    (m - 1 - n_zero) * (T_steps - 1)
  })
  
  p_val = mapply(function(x, df) pchisq(x, df, lower.tail = FALSE), chi2_i, df_i)
  
  data.frame(
    state       = states,
    chi2        = round(chi2_i, 4),
    df          = df_i,
    p_value     = round(p_val,  4),
    significant = p_val < alpha,
    row.names   = NULL
  )
}

chi2_joint = function(count_list, alpha = 0.05) {
  row_results = chi2_row(count_list, alpha)
  chi2_total  = sum(row_results$chi2)
  df_total    = sum(row_results$df)
  p_val       = pchisq(chi2_total, df = df_total, lower.tail = FALSE)
  
  cat("\n--- Joint test (eq. 3.8) ---\n")
  cat(sprintf("chi2 = %.4f,  df = %d,  p-value = %.4f\n",
              chi2_total, df_total, p_val))
  cat(sprintf("Decision at alpha = %.2f: %s\n", alpha,
              ifelse(p_val < alpha,
                     "REJECT H0 — transition probabilities are not stationary",
                     "Fail to reject H0 — stationarity is consistent with data")))
  
  joint_result = data.frame(
    state       = "Joint",
    chi2        = round(chi2_total, 4),
    df          = df_total,
    p_value     = round(p_val, 4),
    significant = p_val < alpha
  )
  return(joint_result)
}