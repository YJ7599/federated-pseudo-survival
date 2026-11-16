## =============================================================================
## Figure 1: bias of the estimated treatment log-HR under the PH setting
##           (Cox, ODAC and federated pseudo-value regression).
##
## Input : results/scenario1_ph/<scenario>/all_simulation_results.csv
## Output: figures/Figure1_bias_boxplot.png
##
## Run from the project root:  Rscript visualization/Figure1_bias_boxplot.R
## =============================================================================
library(ggplot2)

repo_root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f), "..")) else getwd()
})
res_dir <- Sys.getenv("FEDSURV_PH_RESULTS", file.path(repo_root, "results", "scenario1_ph"))

read_estimates <- function(scenario) {
  x <- read.csv(file.path(res_dir, scenario, "all_simulation_results.csv"))
  list(cox  = x$cox_beta_Treatment[is.finite(x$cox_beta_Treatment)],
       odac = x$odac_beta_Treatment[is.finite(x$odac_beta_Treatment)],
       pr   = x$if_beta_Treatment[is.finite(x$if_beta_Treatment)])
}
s1 <- read_estimates("scenario1_er30_10sites_1000x5_500x5")
s3 <- read_estimates("scenario3_er10_10sites_1000x5_500x5")
s5 <- read_estimates("scenario5_er30_mixed18sites_1000x2_500x6_250x10")
s7 <- read_estimates("scenario7_er10_20sites_500x10_250x10")

# True treatment effect: log(1.15)
log_1_15 <- log(1.15)

# Reshape estimates and compute bias
reshape_data <- function(cox_vector, odac_vector, pr_vector, event_rate_label, dataset_label) {
  cox_df <- data.frame(Bias = cox_vector - log_1_15, Method = "Cox", Condition = event_rate_label, Dataset = dataset_label)
  odac_df <- data.frame(Bias = odac_vector - log_1_15, Method = "ODAC", Condition = event_rate_label, Dataset = dataset_label)
  pr_df <- data.frame(Bias = pr_vector - log_1_15, Method = "Pseudo Reg", Condition = event_rate_label, Dataset = dataset_label)
  return(rbind(cox_df, odac_df, pr_df))
}

data_1 <- reshape_data(s1$cox, s1$odac, s1$pr, "event_rate = 30%", "Scenario_1")
data_2 <- reshape_data(s3$cox, s3$odac, s3$pr, "event_rate = 10%", "Scenario_3")
data_3 <- reshape_data(s5$cox, s5$odac, s5$pr, "event_rate = 30%", "Scenario_5")
data_4 <- reshape_data(s7$cox, s7$odac, s7$pr, "event_rate = 10%", "Scenario_7")

combined_data <- rbind(
  data_1, data_2, data_3, data_4
)

custom_labels <- c(
  "Scenario_1" = "1000*5, 500*5",
  "Scenario_3" = "1000*5, 500*5",
  "Scenario_5" = "500*10, 250*10",
  "Scenario_7" = "500*10, 250*10"
)

combined_data$Dataset <- factor(
  combined_data$Dataset,
  levels = names(custom_labels),
  labels = custom_labels
)

p <- ggplot(combined_data, aes(x = Dataset, y = Bias, fill = Method)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5) +
  facet_wrap(~ Condition, ncol = 2) +
  scale_fill_manual(values = c("Cox" = "#EF3167", "ODAC" = "#99FF99", "Pseudo Reg" = "#87CEFA")) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red", linewidth = 1) +
  theme_minimal(base_size = 12) +
  labs(
    x = "The Size of the Local Sites",
    y = "Bias in Beta Estimates"
  ) +
  scale_y_continuous(limits = c(-0.3, 0.3)) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold")
  )

fig_dir <- Sys.getenv("FEDSURV_FIG_DIR", file.path(repo_root, "figures"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
png(file.path(fig_dir, "Figure1_bias_boxplot.png"), width = 2418, height = 1242, res = 190)
print(p)
invisible(dev.off())
