## =============================================================================
## Figure 4: simulation-average RMSE of the site-specific treatment-effect
##           estimates under sparse site heterogeneity (K = 50, n_k = 100).
##
##   Global   : federated global estimate used for every site
##   Local    : site-wise Cox estimate
##   Debiased : GSURE soft-thresholding (method "GSURE_soft")
## RMSE of a replicate = sqrt(mean_k (beta_hat_k - beta_k)^2); points are the
## averages over replicates.
##
## Input : results/scenario3_sparse_heterogeneity/**/metrics_by_sim.csv
## Output: figures/Figure4_debiasing_rmse.png
##         results/scenario3_sparse_heterogeneity/figure4_summary.csv
##
## Run from the project root:  Rscript visualization/Figure4_debiasing_rmse.R
## =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

repo_root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f), "..")) else getwd()
})
res_root <- Sys.getenv("FEDSURV_HET_RESULTS",
                       file.path(repo_root, "results", "scenario3_sparse_heterogeneity"))
debias_method <- Sys.getenv("FEDSURV_DEBIAS_METHOD", "GSURE_soft")

files <- list.files(res_root, pattern = "^metrics_by_sim\\.csv$", recursive = TRUE, full.names = TRUE)
if (!length(files)) stop("No metrics_by_sim.csv found under ", res_root)
metrics <- bind_rows(lapply(files, read.csv)) %>%
  filter(method == debias_method) %>%
  distinct(sim, M, pi, tau, method, .keep_all = TRUE)

summary_df <- metrics %>%
  group_by(M, pi, tau) %>%
  summarise(
    nsim = n(),
    Global = mean(sqrt(mse_global)),
    Local = mean(sqrt(mse_local)),
    Debiased = mean(sqrt(mse_debiased)),
    Global_sqrt_mean_mse = sqrt(mean(mse_global)),
    Local_sqrt_mean_mse = sqrt(mean(mse_local)),
    Debiased_sqrt_mean_mse = sqrt(mean(mse_debiased)),
    .groups = "drop"
  )
write.csv(summary_df, file.path(res_root, "figure4_summary.csv"), row.names = FALSE)
print(as.data.frame(summary_df), digits = 4)

plot_df <- summary_df %>%
  select(pi, tau, Global, Local, Debiased) %>%
  tidyr::pivot_longer(c(Global, Local, Debiased), names_to = "Method", values_to = "RMSE") %>%
  mutate(Method = factor(Method, levels = c("Global", "Local", "Debiased")))

p <- ggplot(plot_df, aes(x = tau, y = RMSE, color = Method, linetype = Method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ pi, nrow = 1, labeller = label_bquote(pi == .(pi))) +
  scale_x_continuous(breaks = c(0.1, 0.3, 0.7, 1.0), labels = c("0.1", "0.3", "0.7", "1.0")) +
  scale_color_manual(values = c(Global = "#E41A1C", Local = "black", Debiased = "#B04FE0")) +
  scale_linetype_manual(values = c(Global = "solid", Local = "solid", Debiased = "dashed")) +
  labs(x = expression(tau), y = "Root Mean Square Error (RMSE)", color = NULL, linetype = NULL) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

fig_dir <- Sys.getenv("FEDSURV_FIG_DIR", file.path(repo_root, "figures"))
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
png(file.path(fig_dir, "Figure4_debiasing_rmse.png"), width = 2130, height = 1300, res = 300)
print(p)
invisible(dev.off())
