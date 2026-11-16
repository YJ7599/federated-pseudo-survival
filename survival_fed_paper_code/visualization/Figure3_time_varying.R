## =============================================================================
## Figure 3: time-varying treatment effect under the non-PH setting.
##
## Boxplots of the federated estimates beta_hat(t) at the seven landmark times,
## with the true trajectory (log cumulative-hazard ratio at X = 0) implied by the
## Weibull model of simulation/scenario2_time_varying/sim_time_varying.R:
##   beta(t) = alpha1 (log t - log 15 - b_A) - alpha0 (log t - log 15),
##   alpha0 = 0.5, alpha1 = 1.0, b_A = log(1.5).
##
## Input : results/scenario2_time_varying/landmark_times_and_truth.csv
##         results/scenario2_time_varying/time_depen_comb_peter_7_timepoints.csv
## Output: figures/Figure3_time_varying.png
##
## Run from the project root:  Rscript visualization/Figure3_time_varying.R
## =============================================================================
library(ggplot2)

repo_root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f), "..")) else getwd()
})
res_dir <- Sys.getenv("FEDSURV_TV_RESULTS", file.path(repo_root, "results", "scenario2_time_varying"))

landmark_df <- read.csv(file.path(res_dir, "landmark_times_and_truth.csv"))
beta_hat <- as.matrix(read.csv(file.path(res_dir, "time_depen_comb_peter_7_timepoints.csv"), check.names = FALSE))
beta_hat <- beta_hat[stats::complete.cases(beta_hat), , drop = FALSE]
quantile_time <- landmark_df$time

alpha0 <- 0.5
alpha1 <- 1.0
b_A <- log(1.5)
true_beta_fn <- function(t) {
  alpha1 * (log(t) - log(15)) -
    alpha1 * b_A -
    alpha0 * (log(t) - log(15))
}
stopifnot(isTRUE(all.equal(true_beta_fn(quantile_time), landmark_df$true_beta)))

df_long_comb <- data.frame(
  time = rep(quantile_time, each = nrow(beta_hat)),
  beta = as.vector(beta_hat)
)
time_cont <- seq(min(quantile_time), max(quantile_time), length.out = 200)
real_df <- data.frame(time = time_cont, estimate = true_beta_fn(time_cont))

p <- ggplot() +
  geom_boxplot(data = df_long_comb,
               aes(x = time, y = beta, group = time, fill = "Beta Estimates"),
               color = "black", linewidth = 0.6, width = 1.5, outlier.size = 0.8) +
  geom_line(data = real_df,
            aes(x = time, y = estimate, color = "True Beta trajectory"),
            linewidth = 0.8) +
  scale_fill_manual(values = c("Beta Estimates" = "#118675")) +
  scale_color_manual(values = c("True Beta trajectory" = "#fd7960")) +
  scale_x_continuous(breaks = quantile_time,
                     labels = scales::number_format(accuracy = 0.01)) +
  labs(x = "Landmark time t", y = "Estimated Treatment Effect", color = NULL, fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85"),
    legend.position = "bottom",
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, face = "bold", color = "black")
  )

fig_dir <- Sys.getenv("FEDSURV_FIG_DIR", file.path(repo_root, "figures"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
png(file.path(fig_dir, "Figure3_time_varying.png"), width = 1964, height = 1242, res = 200)
print(p)
invisible(dev.off())
