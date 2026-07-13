# ================================================
# Script Name: Statistic Test of TM Homogeneity
# Purpose: Calculate the Chi-square Tests based on
#          Anderson and Goodman [1957]
# Input: complete_info_filled_0704.rds, 
#         highest_job_0704.rds, data1.dta,
#.        approach_3_ID_baseline.csv, marriage_info_uniq_final.rds
# Run order: 3 of 4
# ================================================
#----------Preliminaries----------#
rm(list = ls())
section = "Data"
subsection = "Chi_Square_Test"

dir_root = "~/Desktop/Qing_Final" # Please Change the directory to yours 
setwd(dir_root)

# Define subdirectories for logs and figures:
dir_plot = paste0(dir_root, "/" ,"Plot", "/", subsection)
# Ensure all necessary directories exist under your root folder
# if not, the function will create folders for you

create_dir_if_missing = function(dir) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message("Created directory: ", dir)
  } else {
    message("Directory already exists: ", dir)
  }
}

create_dir_if_missing(dir_root)
create_dir_if_missing(dir_plot)
options(warn = -1)

#----------Preliminaries----------#

packages =
  c("readxl",
    "tidyverse",
    "jtools",
    "MASS",
    "dplyr",
    "haven",
    "imputeTS",
    "data.table",
    "assertthat",
    "foreign",
    "philentropy",
    "ggplot2"
  )

install_and_load = function(pkg_list) {
  for (pkg in pkg_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {  
      message("Installing missing package: ", pkg)
      install.packages(pkg, dependencies = TRUE)  
    }
    library(pkg, character.only = TRUE)  
  }
}

install_and_load(packages)



#================================
# Data Preparation
# ================================

# ---- 1. 6gen Dyads ----

tm_6gen_data =
  readRDS(
    file.path(
      dir_root,
      "Data",
      "main_results",
      "contingency_6gen.rds")
    ) 

# ---- 2. Survey Wave Dyads ----
sw_file_path =
  file.path(
    dir_root,
    "Data",
    "data_process",
    "transition_matrix_6jobs_num_sw3_0717.xlsx"
  )
sheet_names =  excel_sheets(file.path(sw_file_path))

sw_list = lapply(sheet_names, function(sheet) {
  read_excel(sw_file_path, sheet = sheet)
})

df_names = sub(".*_(\\d+)$", "\\1",sheet_names)
names(sw_list) = df_names

# 3.----- Generate 6gen Contingency Table -----

contingency_6gen_list =
  lapply(
    tm_6gen_data,
    function(df){
      contingency = table(
        Parent = factor(df$stat_p, levels = 6:1),
        Child  = factor(df$stat_ch, levels = 6:1)
      )
      as.matrix(contingency)
    }
  )

transition_6gen_list =
  lapply(
    contingency_6gen_list,
    function(df){
      transition = df / rowSums(df)
    }
  )

# 4.----- Generate SW Contingency Table:

sw_6class_list = lapply(sw_list, function(df) {
  m = as.matrix(df |> dplyr::select(-rowsum, -validation))
  rownames(m) = 6:1
  colnames(m) = 6:1
  dimnames(m) = list(Parent = 6:1, Child = 6:1)
  return(m)
})


#======================================
# Calculate Chi-Square Test Statistics:
# =====================================

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


#================================
# Calculation:
# ================================

# ------ 1. All 6gen: row-wise & joint ------
gen6_test_row = chi2_row (contingency_6gen_list, alpha = 0.05)
gen6_test_all = chi2_joint(contingency_6gen_list, alpha = 0.05)

# ------ 2. Gem 2,3 & 4: row-wise & joint ------
gen6_test_234_row = chi2_row(contingency_6gen_list[c(2,3,4)], alpha = 0.05)
gen6_test_234_all = chi2_joint(contingency_6gen_list[c(2,3,4)], alpha = 0.05)


# -------3. Cross- Section ------------

years = c(1783,1792,1795,1798,1801,1804,1807,1810,1813,1816,1819,1822,
           1825,1828,1831,1837,1840,1843,1846,1849,1852,1855,1858,1861,
           1864,1867,1870,1873,1876,1879,1882,1885,1888,1903,1906,1909)
sw_6class_list_f =
  sw_6class_list[
    as.numeric(names(sw_6class_list)) %in% years
]

years_idx = seq_along(years)   # index 1 to 36

# Bin widths to compare — keep small to avoid pooling too many cohorts
bin_widths = c(2, 3, 4, 5, 6)

results = map_dfr(bin_widths, function(k) {
  
  half = floor(k / 2)
  
  # Only iterate over indices where a full window of size k is available
  valid_idx = (half + 1):(length(years) - (k - half - 1))
  
  map_dfr(valid_idx, function(i) {
    
    # k consecutive cohorts centered on index i
    window_idx   = (i - half):(i - half + k - 1)
    window_years = years[window_idx]
    center_year  = years[i]
    
    matrices = sw_6class_list_f[
      as.numeric(names(sw_6class_list_f)) %in% window_years
    ]
    
    if (length(matrices) < 2) return(NULL)
    
    test = chi2_joint(matrices, alpha = 0.05)
    
    tibble(
      k           = k,
      center_year = center_year,
      year_min    = min(window_years),
      year_max    = max(window_years),
      n_matrices  = length(matrices),
      chi2        = test$chi2,      # adjust to your chi2_joint output
      p_value     = test$p_value,   # adjust to your chi2_joint output
      significant = test$significant
    )
  })
})

# Visualize the results:

plot_data = results |>
  mutate(k_label = paste0("k = ", k))

threshold = 0.05


p = ggplot(plot_data, aes(x = center_year, y = p_value,
                           color = k_label, group = k_label,
                           shape = k_label)) +
  geom_rect(data = data.frame(xmin = 1846, xmax = 1858, ymin = -Inf, ymax = Inf),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey60", alpha = 0.5) +
  geom_hline(yintercept = threshold, linetype = "dashed",
             color = "grey30", linewidth = 0.5) +
  geom_vline(xintercept = 1850, linetype = "dotted",
             color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.65, alpha = 0.85) +
  geom_point(size = 4, alpha = 0.9) +
  scale_color_manual(
    values = c("indianred4", "goldenrod1", "darkgreen", "royalblue4", "darkslateblue"),
    name   = "Bin width"
  ) +
  scale_shape_manual(
    values = c(4, 15, 18, 17, 19),
    name   = "Bin width"
  ) +
  scale_x_continuous(breaks = years[seq(1, length(years), by = 2)],
                     guide  = guide_axis(angle = 45)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
  labs(
    x       = "Center cohort year",
    y       = "Joint test p-value",
    caption = "Each point is a joint homogeneity test (Anderson-Goodman 1957, §3.4) on k cohorts centered on that year. Shaded region: 1846–1858. Dashed line: p = 0.05."
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position    = "bottom",
    legend.key.width   = unit(1.5, "cm"),
    legend.text        = element_text(size = 11),
    legend.title       = element_text(size = 11),
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank(),
    plot.caption       = element_text(size = 8, color = "grey40", hjust = 0),
    axis.title.y       = element_text(margin = margin(r = 8))
  )

ggsave(
  file.path(
    dir_plot,
    "rolling_window_sensitivity.png"
  ), 
  p, 
  width = 12, 
  height = 7, 
  dpi = 600
  )

