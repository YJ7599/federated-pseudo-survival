## =============================================================================
## Figure 2: average model-based SE / empirical SD of the treatment log-HR
##           under the PH setting (Cox, ODAC, federated pseudo-value regression).
##
## Input : results/scenario1_ph/<scenario>/all_simulation_results.csv
## Output: figures/figure2_ph_se_calibration.png
##         results/scenario1_ph/figure2_summary.csv
##           (bias, RMSE, empirical SD, average SE, SE/SD ratio, 95% CI coverage)
##
## Run from the project root:  Rscript visualization/figure2_ph_se_calibration.R
## =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

repo_root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f), "..")) else getwd()
})
results_root <- Sys.getenv("FEDSURV_PH_RESULTS", file.path(repo_root, "results", "scenario1_ph"))

true_beta <- log(1.15)
z_975 <- qnorm(0.975)

scenario_info <- data.frame(
  scenario = c(
    "ph_er10_k10_1000x5_500x5",
    "ph_er10_k20_500x10_250x10",
    "ph_er30_k10_1000x5_500x5",
    "ph_er30_k18_1000x2_500x6_250x10"
  ),
  event_rate = c("event_rate = 10%", "event_rate = 10%", "event_rate = 30%", "event_rate = 30%"),
  site_size = c("1000*5, 500*5", "500*10, 250*10", "1000*5, 500*5", "500*10, 250*10"),
  stringsAsFactors = FALSE
)

method_info <- data.frame(
  method = c("Cox", "ODAC", "Pseudo Reg"),
  beta_col = c("cox_beta_Treatment", "odac_beta_Treatment", "if_beta_Treatment"),
  se_col = c("cox_se_Treatment", "odac_se_Treatment", "if_se_Treatment"),
  stringsAsFactors = FALSE
)

read_one_scenario <- function(info_row) {
  result_file <- file.path(results_root, info_row$scenario, "all_simulation_results.csv")
  if (!file.exists(result_file)) stop("Cannot find: ", result_file)
  raw <- read.csv(result_file, check.names = FALSE)
  bind_rows(lapply(seq_len(nrow(method_info)), function(i) {
    beta <- raw[[method_info$beta_col[i]]]
    se <- raw[[method_info$se_col[i]]]
    data.frame(
      scenario = info_row$scenario, event_rate = info_row$event_rate, site_size = info_row$site_size,
      sim_id = raw$sim_id, method = method_info$method[i], beta = beta, se = se,
      valid = is.finite(beta) & is.finite(se),
      covered = beta - z_975 * se <= true_beta & beta + z_975 * se >= true_beta,
      stringsAsFactors = FALSE
    )
  }))
}

plot_data <- bind_rows(lapply(seq_len(nrow(scenario_info)), function(i) read_one_scenario(scenario_info[i, ]))) %>%
  mutate(
    event_rate = factor(event_rate, levels = c("event_rate = 10%", "event_rate = 30%")),
    site_size = factor(site_size, levels = c("1000*5, 500*5", "500*10, 250*10")),
    method = factor(method, levels = c("Cox", "ODAC", "Pseudo Reg"))
  )

summary_data <- plot_data %>%
  filter(valid) %>%
  group_by(scenario, event_rate, site_size, method) %>%
  summarise(
    nsim_success = n(),
    bias = mean(beta - true_beta),
    rmse = sqrt(mean((beta - true_beta)^2)),
    empirical_sd = sd(beta),
    average_model_se = mean(se),
    se_to_sd_ratio = average_model_se / empirical_sd,
    coverage = mean(covered),
    .groups = "drop"
  )
write.csv(summary_data, file.path(results_root, "figure2_summary.csv"), row.names = FALSE)
print(as.data.frame(summary_data), digits = 4)

method_colors <- c("Cox" = "#EF3167", "ODAC" = "#99FF99", "Pseudo Reg" = "#87CEFA")

p_ratio <- ggplot(summary_data, aes(x = site_size, y = se_to_sd_ratio, color = method, group = method)) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "red", linewidth = 1) +
  geom_point(position = position_dodge(width = 0.45), size = 3) +
  facet_wrap(~ event_rate, ncol = 2) +
  scale_color_manual(values = method_colors, breaks = c("Cox", "ODAC", "Pseudo Reg")) +
  scale_y_continuous(limits = c(0.90, 1.10), breaks = seq(0.90, 1.10, by = 0.10)) +
  labs(x = "The Size of the Local Sites", y = "Average Model-based SE / Empirical SD") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_line(color = "grey92", linewidth = 0.35),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.55),
    strip.text = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold"),
    legend.position = "right",
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 13, face = "bold"),
    axis.text.y = element_text(size = 13, face = "bold")
  )

fig_dir <- Sys.getenv("FEDSURV_FIG_DIR", file.path(repo_root, "figures"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
png(file.path(fig_dir, "figure2_ph_se_calibration.png"), width = 2915, height = 1645, res = 300)
print(p_ratio)
invisible(dev.off())
