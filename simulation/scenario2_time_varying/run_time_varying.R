############################################################
## Scenario 2 (Section 3.2; Figure 3): time-varying treatment effect.
##
## Weibull event times with treatment-dependent shape (non-PH):
##   shape alpha(A) = 0.5 + 0.5 A,
##   scale lambda(X, A) = 15 exp(log(0.7) X + log(1.5) A),
##   X ~ N(0, 1), A ~ Bernoulli(expit(0.5 X)),
##   censoring Exp(rate = -log(0.5)/20) truncated at 40.
## N = 5,500 subjects over 19 sites (3000, 500 x 3, 100 x 5, 50 x 10).
## Landmark times: 20th-80th percentiles of the pooled KM curve of a large
## reference sample. Truth: log cumulative-hazard ratio at X = 0,
##   beta(t) = 1.0 (log t - log 15 - log 1.5) - 0.5 (log t - log 15).
## Estimation: federated KM splines -> pseudo-observations at the landmark
## times -> renewable GEE (cloglog) with treatment-by-time interactions.
## Replicate j uses set.seed(base_seed + j); results are cached per replicate.
##
## Outputs (results/scenario2_time_varying/): per-replicate estimates, SEs and
## CIs, landmark_times_and_truth.csv, reviewer_summary_bias_se_coverage.csv,
## summary_bias_variance_mse_coverage_by_time.csv and diagnostic plots.
############################################################

library(survival)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pseudo)
library(pda)
library(jsonlite)
library(parallel)

############################################################
## User settings
############################################################

project_root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f), "..", "..")) else getwd()
})
setwd(project_root)

nsim <- 1000L
ncore <- as.integer(Sys.getenv("FEDSURV_NCORE", "50"))
base_seed <- 20260703L

output_dir <- Sys.getenv("FEDSURV_OUTPUT_DIR",
                         file.path(project_root, "results", "scenario2_time_varying"))
per_sim_dir <- file.path(output_dir, "per_sim_rds")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(per_sim_dir, recursive = TRUE, showWarnings = FALSE)

source_required <- function(file_name) {
  path <- file.path(project_root, file_name)
  if (!file.exists(path)) stop("Required source file not found: ", path)
  source(path)
}

source_optional <- function(file_name) {
  path <- file.path(project_root, file_name)
  if (file.exists(path)) source(path)
}

source_optional("datagenerator.inloop.R")
source_optional("datagenerator.outloop.R")
source_optional("RenewGLM.inloop.R")
source_optional("RenewGLM.outloop.R")
source_required(file.path("R", "renewable_gee.R"))
source_required(file.path("R", "federated_km_pseudo.R"))

############################################################
## Helper functions
############################################################

true_beta_fn <- function(t) {
  1.0 * (log(t) - log(15)) -
    1.0 * log(1.5) -
    0.5 * (log(t) - log(15))
}

save_numeric_matrix <- function(x, file) {
  write.csv(as.data.frame(x), file, row.names = FALSE)
}

read_numeric_matrix <- function(file) {
  as.matrix(read.csv(file, check.names = FALSE))
}

theme_manuscript <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "top",
      legend.title = element_blank()
    )
}

############################################################
## Global landmark times and truth
############################################################

set.seed(base_seed)

X_global <- rnorm(1000000, 0, 1)
Treatment_global <- rbinom(
  1000000, 1,
  exp(0.5 * X_global) / (1 + exp(0.5 * X_global))
)

baseline_rate_T_global <- 15 * exp(X_global * log(0.7)) * exp(Treatment_global * log(1.5))
alpha_global <- 0.5 + 0.5 * Treatment_global
Survival_Time_global <- rweibull(1000000, shape = alpha_global, scale = baseline_rate_T_global)

Censoring_Time_global <- rexp(1000000, rate = -log(0.5) / 20)
Censoring_Time_global <- pmin(Censoring_Time_global, 40)

Observed_Time_global <- pmin(Survival_Time_global, Censoring_Time_global)
Event_global <- ifelse(Survival_Time_global <= Censoring_Time_global, 1, 0)

global_curve <- survfit(Surv(Observed_Time_global, Event_global) ~ 1)
quantile_time <- quantile(global_curve, probs = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))[[1]]
true_beta <- true_beta_fn(quantile_time)
Tn <- length(quantile_time)
time_cols <- paste0("t", seq_len(Tn))

write.csv(
  data.frame(time_index = seq_along(quantile_time), time = quantile_time, true_beta = true_beta),
  file.path(output_dir, "landmark_times_and_truth.csv"),
  row.names = FALSE
)

############################################################
## One simulation replicate
############################################################

run_one_sim <- function(j) {
  per_sim_file <- file.path(per_sim_dir, sprintf("sim_%04d.rds", j))
  if (file.exists(per_sim_file)) return(readRDS(per_sim_file))
  
  set.seed(base_seed + j)
  
  fail_result <- function(msg) {
    list(
      sim = j,
      success = FALSE,
      message = msg,
      beta_ind = rep(NA_real_, Tn),
      beta_comb = rep(NA_real_, Tn),
      bias = rep(NA_real_, Tn),
      se = rep(NA_real_, Tn),
      lower = rep(NA_real_, Tn),
      upper = rep(NA_real_, Tn),
      coverage = rep(NA_real_, Tn),
      quantile_time_check = rep(NA_real_, Tn + 1L),
      p_value = NA_real_
    )
  }
  
  res <- tryCatch({
    Ns <- c(3000L, rep(500L, 3), rep(100L, 5), rep(50L, 10))
    N <- sum(Ns)
    Site <- rep(seq_along(Ns), times = Ns)
    
    X <- rnorm(N, 0, 1)
    Treatment <- rbinom(N, 1, exp(0.5 * X) / (1 + exp(0.5 * X)))
    
    baseline_rate_T <- 15 * exp(X * log(0.7)) * exp(Treatment * log(1.5))
    alpha <- 0.5 + 0.5 * Treatment
    
    Survival_Time <- rweibull(N, shape = alpha, scale = baseline_rate_T)
    Censoring_Time <- rexp(N, rate = -log(0.5) / 20)
    Censoring_Time <- pmin(Censoring_Time, 40)
    
    Observed_Time <- pmin(Survival_Time, Censoring_Time)
    Event <- ifelse(Survival_Time <= Censoring_Time, 1, 0)
    
    sim_data_ori <- data.frame(
      Site = Site,
      Treatment = Treatment,
      X = X,
      Survival_Time = Survival_Time,
      Censoring_Time = Censoring_Time,
      delta = Event,
      time = Observed_Time
    )
    
    first_site_data <- sim_data_ori[sim_data_ori$Site == 1, , drop = FALSE]
    updated_splines <- get_survival_splines_cloglog(
      first_site_data,
      degree = 2,
      knots = 1,
      bknots = c(0.1, 39.9),
      risk_set_minimum = 0.05
    )
    
    max_knots <- 7
    current_knots <- 1
    chunksize <- 5
    counter <- 0
    
    for (i in 2:length(Ns)) {
      site_data <- sim_data_ori[sim_data_ori$Site == i, , drop = FALSE]
      if (nrow(site_data) == 0) next
      
      n_batches <- nrow(site_data) / chunksize
      knots_list <- NULL
      
      if (nrow(site_data) > 100) {
        for (k in seq_len(nrow(site_data) / 100)) {
          knots_list <- c(
            knots_list,
            c(
              rep(current_knots, 100 / chunksize),
              rep(current_knots + k, n_batches - 100 / chunksize)
            )
          )
        }
      } else {
        if (counter > 200 && current_knots < max_knots) {
          knots_list <- rep(current_knots + 1, n_batches)
          counter <- 0
        } else {
          knots_list <- rep(current_knots, n_batches)
        }
        counter <- counter + nrow(site_data)
      }
      
      knots_list <- ifelse(knots_list > max_knots, max_knots, knots_list)
      updated_splines <- update_survival_splines_chunk(
        updated_splines,
        breakpoints = 10,
        time = site_data$time,
        delta = site_data$delta,
        chunksize = chunksize,
        knots_list = knots_list,
        loud = FALSE,
        restriction = 0.2,
        risk_floor = 0.005
      )
      current_knots <- knots_list[length(knots_list)]
    }
    
    final_ST <- function(t) {
      survival_function_cloglog.v(
        t,
        updated_splines$Event$Knots,
        updated_splines$Event$Degree,
        updated_splines$Event$BKnots,
        updated_splines$Event$Coefficients
      )
    }
    
    quantile_function <- function(surv_func, grid_points, q) {
      ST_values <- surv_func(grid_points)
      out <- numeric(length(q))
      for (ii in seq_along(q)) {
        candidate_times <- grid_points[ST_values > q[ii]]
        out[ii] <- if (length(candidate_times) == 0) NA_real_ else max(candidate_times)
      }
      out
    }
    
    quantile_time_check <- rev(
      quantile_function(final_ST, seq(0, max(sim_data_ori$time), 0.01), seq(0.2, 0.8, length.out = 8))
    )
    
    pseudo_list <- lapply(
      quantile_time,
      function(t0) {
        theoretical_KM_PS(
          spline_params = updated_splines,
          t_ps = t0,
          time = sim_data_ori$time,
          delta = sim_data_ori$delta,
          restriction = 2
        )
      }
    )
    
    length_ori_dat <- nrow(sim_data_ori)
    sim_data <- do.call(rbind, rep(list(sim_data_ori), Tn))
    sim_data$txx <- as.factor(rep(quantile_time, each = length_ori_dat))
    sim_data$Pseudo_Spline <- 1 - unlist(pseudo_list)
    
    ########################################################
    ## Individual time-specific fits
    ########################################################
    
    beta_ind_vec <- rep(NA_real_, Tn)
    
    for (q in seq_len(Tn)) {
      sim_dat_ind <- sim_data[sim_data$txx == quantile_time[q], , drop = FALSE]
      p_ind <- 4
      
      sim_dat_ind_site <- sim_dat_ind[sim_dat_ind$Site == 1, , drop = FALSE]
      y <- sim_dat_ind_site$Pseudo_Spline
      X_mat <- cbind(
        1,
        sim_dat_ind_site$Treatment,
        sim_dat_ind_site$X,
        sim_dat_ind_site$Treatment * sim_dat_ind_site$X
      )
      
      infomats <- diag(0, p_ind, p_ind)
      betahat <- rep(0, p_ind)
      summary_ind <- RenewGLM_out_ori(
        X_mat, y, "cloglog", betahat, infomats,
        intercept = FALSE, s = NULL, phi = NULL
      )
      
      betahat <- summary_ind[[1]]
      infomats <- summary_ind[[2]]
      
      for (i in 2:length(Ns)) {
        sim_dat_ind_site <- sim_dat_ind[sim_dat_ind$Site == i, , drop = FALSE]
        X_mat <- cbind(
          1,
          sim_dat_ind_site$Treatment,
          sim_dat_ind_site$X,
          sim_dat_ind_site$Treatment * sim_dat_ind_site$X
        )
        y <- sim_dat_ind_site$Pseudo_Spline
        
        summary_ind <- tryCatch(
          RenewGLM_out_ori(
            X_mat, y, "cloglog", betahat, infomats,
            intercept = FALSE, s = NULL, phi = NULL
          ),
          error = function(e) NULL
        )
        
        if (!is.null(summary_ind)) {
          betahat <- summary_ind[[1]]
          infomats <- summary_ind[[2]]
        }
      }
      
      beta_ind_vec[q] <- betahat[2]
    }
    
    ########################################################
    ## Combined time-varying treatment-effect fit
    ########################################################
    
    make_combined_design <- function(dat) {
      X_mat <- cbind(
        1,
        dat$Treatment,
        dat$X,
        dat$Treatment * dat$X,
        as.factor(dat$txx)
      )
      X_mat <- as.data.frame(X_mat)
      X_mat[, 5] <- as.factor(X_mat[, 5])
      encoded <- model.matrix(~ V5 - 1, data = X_mat)
      encoded <- encoded[, -1, drop = FALSE]
      as.matrix(cbind(X_mat[, c(1, 2, 3, 4)], encoded, X_mat[, 2] * encoded))
    }
    
    sim_site_1 <- sim_data[sim_data$Site == 1, , drop = FALSE]
    y <- sim_site_1$Pseudo_Spline
    X_mat <- make_combined_design(sim_site_1)
    p_comb <- ncol(X_mat)
    infomats <- diag(0, p_comb, p_comb)
    betahat <- rep(0, p_comb)
    id <- rep(seq_len(Ns[1]), Tn)
    meat_list <- list()
    
    summary_comb <- RenewGLM_out(
      X_mat, y, id, "cloglog", betahat, infomats, meat_list,
      intercept = FALSE, s = NULL, phi = NULL
    )
    
    betahat <- summary_comb[[1]]
    infomats <- summary_comb[[2]]
    meat_list[[1]] <- summary_comb[[3]]
    
    for (i in 2:length(Ns)) {
      site_data <- sim_data[sim_data$Site == i, , drop = FALSE]
      if (nrow(site_data) == 0) next
      
      y <- site_data$Pseudo_Spline
      X_mat <- make_combined_design(site_data)
      id <- rep(seq_len(Ns[i]), Tn)
      
      summary_comb <- tryCatch(
        RenewGLM_out(
          X_mat, y, id, "cloglog", betahat, infomats, meat_list,
          intercept = FALSE, s = NULL, phi = NULL
        ),
        error = function(e) NULL
      )
      
      if (!is.null(summary_comb)) {
        betahat <- summary_comb[[1]]
        infomats <- summary_comb[[2]]
        meat_list[[i]] <- summary_comb[[3]]
      }
    }
    
    valid_meat <- meat_list[!vapply(meat_list, is.null, logical(1))]
    I_combined <- infomats
    M_combined <- Reduce(`+`, valid_meat)
    vcov_beta <- solve(I_combined) %*% M_combined %*% solve(I_combined)
    
    interaction_cols <- (4 + Tn):((4 + Tn * 2) - 2)
    L <- matrix(0, nrow = Tn, ncol = p_comb)
    L[, 2] <- 1
    for (q in 2:Tn) {
      L[q, interaction_cols[q - 1]] <- 1
    }
    
    beta_comb_vec <- as.vector(L %*% betahat)
    se_comb_vec <- sqrt(pmax(diag(L %*% vcov_beta %*% t(L)), 0))
    lower_vec <- beta_comb_vec - qnorm(0.975) * se_comb_vec
    upper_vec <- beta_comb_vec + qnorm(0.975) * se_comb_vec
    bias_vec <- beta_comb_vec - true_beta
    cover_vec <- as.numeric(lower_vec <= true_beta & true_beta <= upper_vec)
    
    beta_int <- betahat[interaction_cols]
    V_int <- vcov_beta[interaction_cols, interaction_cols, drop = FALSE]
    p_value <- tryCatch({
      wald_stat <- as.numeric(t(beta_int) %*% solve(V_int) %*% beta_int)
      pchisq(wald_stat, df = length(interaction_cols), lower.tail = FALSE)
    }, error = function(e) NA_real_)
    
    list(
      sim = j,
      success = TRUE,
      message = NA_character_,
      beta_ind = beta_ind_vec,
      beta_comb = beta_comb_vec,
      bias = bias_vec,
      se = se_comb_vec,
      lower = lower_vec,
      upper = upper_vec,
      coverage = cover_vec,
      quantile_time_check = quantile_time_check,
      p_value = p_value
    )
  }, error = function(e) {
    fail_result(conditionMessage(e))
  })
  
  saveRDS(res, per_sim_file)
  res
}

############################################################
## Run simulations in parallel
############################################################

message("Running nsim = ", nsim, " with ncore = ", ncore)

sim_ids <- seq_len(nsim)
## Optional: FEDSURV_SIM_IDS="1:20" runs a subset of the replicates.
if (nzchar(Sys.getenv("FEDSURV_SIM_IDS"))) sim_ids <- as.integer(eval(parse(text = Sys.getenv("FEDSURV_SIM_IDS"))))
if (.Platform$OS.type == "unix" && ncore > 1L) {
  results_list <- mclapply(sim_ids, run_one_sim, mc.cores = ncore, mc.preschedule = FALSE)
} else {
  results_list <- lapply(sim_ids, run_one_sim)
}

extract_matrix <- function(field, ncol_expected = Tn, col_names = time_cols) {
  mat <- do.call(rbind, lapply(results_list, function(x) x[[field]]))
  mat <- as.matrix(mat)
  if (ncol(mat) == ncol_expected) colnames(mat) <- col_names
  mat
}

betahat_ind <- extract_matrix("beta_ind")
betahat_combined <- extract_matrix("beta_comb")
bias_combined <- extract_matrix("bias")
model_se_combined <- extract_matrix("se")
lower_ci_combined <- extract_matrix("lower")
upper_ci_combined <- extract_matrix("upper")
coverage_indicator <- extract_matrix("coverage")
quantile_time_dat <- extract_matrix(
  "quantile_time_check",
  ncol_expected = Tn + 1L,
  col_names = paste0("q", seq_len(Tn + 1L))
)
p_val_wald <- vapply(results_list, function(x) x$p_value, numeric(1))

sim_status <- data.frame(
  sim = vapply(results_list, function(x) x$sim, integer(1)),
  success = vapply(results_list, function(x) x$success, logical(1)),
  message = vapply(results_list, function(x) ifelse(is.na(x$message), "", x$message), character(1))
)

save_numeric_matrix(betahat_ind, file.path(output_dir, "time_depen_ind_peter_7_timepoints.csv"))
save_numeric_matrix(betahat_combined, file.path(output_dir, "time_depen_comb_peter_7_timepoints.csv"))
save_numeric_matrix(bias_combined, file.path(output_dir, "bias_combined_7_timepoints.csv"))
save_numeric_matrix(model_se_combined, file.path(output_dir, "model_se_combined_7_timepoints.csv"))
save_numeric_matrix(lower_ci_combined, file.path(output_dir, "lower_b_list_7_timepoints.csv"))
save_numeric_matrix(upper_ci_combined, file.path(output_dir, "upper_b_list_7_timepoints.csv"))
save_numeric_matrix(coverage_indicator, file.path(output_dir, "coverage_indicator_7_timepoints.csv"))
save_numeric_matrix(quantile_time_dat, file.path(output_dir, "quantile_time_mat_7_timepoints.csv"))
write.csv(data.frame(p_value = p_val_wald), file.path(output_dir, "significance_test_7_timepoints.csv"), row.names = FALSE)
write.csv(sim_status, file.path(output_dir, "simulation_status.csv"), row.names = FALSE)

############################################################
## Visualization and reviewer summaries
############################################################

landmark_df <- read.csv(file.path(output_dir, "landmark_times_and_truth.csv"))
times <- landmark_df$time
true_beta <- landmark_df$true_beta

complete_rows <- rowSums(is.na(betahat_combined)) == 0 &
  rowSums(is.na(model_se_combined)) == 0 &
  rowSums(is.na(lower_ci_combined)) == 0 &
  rowSums(is.na(upper_ci_combined)) == 0

beta_hat <- betahat_combined[complete_rows, , drop = FALSE]
bias_mat <- bias_combined[complete_rows, , drop = FALSE]
model_se <- model_se_combined[complete_rows, , drop = FALSE]
lower_ci <- lower_ci_combined[complete_rows, , drop = FALSE]
upper_ci <- upper_ci_combined[complete_rows, , drop = FALSE]

n_success <- nrow(beta_hat)

coverage <- colMeans(
  sweep(lower_ci, 2, true_beta, FUN = "<=") &
    sweep(upper_ci, 2, true_beta, FUN = ">="),
  na.rm = TRUE
)

summary_df <- data.frame(
  time_index = seq_along(times),
  time = times,
  true_beta = true_beta,
  n_success = n_success,
  mean_beta = colMeans(beta_hat, na.rm = TRUE),
  bias = colMeans(bias_mat, na.rm = TRUE),
  empirical_sd = apply(beta_hat, 2, sd, na.rm = TRUE),
  average_model_se = colMeans(model_se, na.rm = TRUE),
  coverage = coverage
) %>%
  mutate(
    rmse = sqrt(colMeans((sweep(beta_hat, 2, true_beta, FUN = "-"))^2, na.rm = TRUE)),
    se_to_sd_ratio = average_model_se / empirical_sd,
    mean_ci_lower = mean_beta - qnorm(0.975) * average_model_se,
    mean_ci_upper = mean_beta + qnorm(0.975) * average_model_se
  )

write.csv(summary_df, file.path(output_dir, "reviewer_summary_bias_se_coverage.csv"), row.names = FALSE)

bias_long <- as.data.frame(bias_mat) %>%
  mutate(sim = row_number()) %>%
  pivot_longer(-sim, names_to = "time_col", values_to = "bias") %>%
  mutate(
    time_index = as.integer(sub("^t", "", time_col)),
    time = times[time_index],
    time_label = sprintf("%.2f", time)
  )

p_bias <- ggplot(bias_long, aes(x = time_label, y = bias)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.7) +
  geom_boxplot(fill = "#56B4E9", color = "grey25", width = 0.65, outlier.size = 0.5, outlier.alpha = 0.35) +
  labs(
    x = "Landmark time",
    y = expression(hat(beta)(t) - beta(t)),
    title = "Bias in time-varying treatment-effect estimates"
  ) +
  theme_manuscript()

p_coverage <- ggplot(summary_df, aes(x = time, y = coverage)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "#D55E00", linewidth = 0.8) +
  geom_line(color = "#0072B2", linewidth = 1) +
  geom_point(color = "#0072B2", size = 2.7) +
  scale_y_continuous(limits = c(0.90, 1.00), breaks = seq(0.90, 1.00, by = 0.02)) +
  labs(
    x = "Landmark time",
    y = "Empirical 95% CI coverage",
    title = "Empirical coverage of model-based 95% confidence intervals"
  ) +
  theme_manuscript()

se_long <- summary_df %>%
  select(time, empirical_sd, average_model_se) %>%
  pivot_longer(-time, names_to = "quantity", values_to = "value") %>%
  mutate(
    quantity = recode(
      quantity,
      empirical_sd = "Empirical SD of estimates",
      average_model_se = "Average model-based SE"
    )
  )

p_se <- ggplot(se_long, aes(x = time, y = value, color = quantity, shape = quantity)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.7) +
  scale_color_manual(values = c(
    "Empirical SD of estimates" = "#D55E00",
    "Average model-based SE" = "#0072B2"
  )) +
  scale_shape_manual(values = c(
    "Empirical SD of estimates" = 16,
    "Average model-based SE" = 17
  )) +
  labs(
    x = "Landmark time",
    y = "Standard error",
    title = "Agreement between model-based SE and empirical SD"
  ) +
  theme_manuscript()

p_se_ratio <- ggplot(summary_df, aes(x = time, y = se_to_sd_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.8) +
  geom_line(color = "#009E73", linewidth = 1) +
  geom_point(color = "#009E73", size = 2.7) +
  scale_y_continuous(limits = c(0.85, 1.15), breaks = seq(0.85, 1.15, by = 0.05)) +
  labs(
    x = "Landmark time",
    y = "Average model-based SE / empirical SD",
    title = "SE calibration across landmark times"
  ) +
  theme_manuscript()

p_mean_ci <- ggplot(summary_df, aes(x = time)) +
  geom_ribbon(aes(ymin = mean_ci_lower, ymax = mean_ci_upper), fill = "#0072B2", alpha = 0.18) +
  geom_line(aes(y = mean_beta, color = "Mean estimate"), linewidth = 1) +
  geom_point(aes(y = mean_beta, color = "Mean estimate"), size = 2.5) +
  geom_line(aes(y = true_beta, color = "Truth"), linewidth = 1, linetype = "dashed") +
  scale_color_manual(values = c("Mean estimate" = "#0072B2", "Truth" = "#D55E00")) +
  labs(
    x = "Landmark time",
    y = "Treatment effect",
    title = "Mean estimate with average model-based 95% CI"
  ) +
  theme_manuscript()

ggsave(file.path(output_dir, "figure_bias_boxplot.png"), p_bias, width = 7.0, height = 4.8, dpi = 600)
ggsave(file.path(output_dir, "figure_bias_boxplot.pdf"), p_bias, width = 7.0, height = 4.8)

ggsave(file.path(output_dir, "figure_empirical_coverage.png"), p_coverage, width = 7.0, height = 4.8, dpi = 600)
ggsave(file.path(output_dir, "figure_empirical_coverage.pdf"), p_coverage, width = 7.0, height = 4.8)

ggsave(file.path(output_dir, "figure_se_vs_sd.png"), p_se, width = 7.0, height = 4.8, dpi = 600)
ggsave(file.path(output_dir, "figure_se_vs_sd.pdf"), p_se, width = 7.0, height = 4.8)

ggsave(file.path(output_dir, "figure_se_to_sd_ratio.png"), p_se_ratio, width = 7.0, height = 4.8, dpi = 600)
ggsave(file.path(output_dir, "figure_se_to_sd_ratio.pdf"), p_se_ratio, width = 7.0, height = 4.8)

ggsave(file.path(output_dir, "figure_mean_estimate_with_average_ci.png"), p_mean_ci, width = 7.0, height = 4.8, dpi = 600)
ggsave(file.path(output_dir, "figure_mean_estimate_with_average_ci.pdf"), p_mean_ci, width = 7.0, height = 4.8)

print(summary_df)
print(p_bias)
print(p_coverage)
print(p_se)
print(p_se_ratio)
print(p_mean_ci)


library(dplyr)
library(readr)

############################################################
## 1. Load saved simulation outputs
############################################################

landmark_df <- read.csv(
  file.path(output_dir, "landmark_times_and_truth.csv")
)

beta_hat <- as.matrix(read.csv(
  file.path(output_dir, "time_depen_comb_peter_7_timepoints.csv"),
  check.names = FALSE
))

lower_ci <- as.matrix(read.csv(
  file.path(output_dir, "lower_b_list_7_timepoints.csv"),
  check.names = FALSE
))

upper_ci <- as.matrix(read.csv(
  file.path(output_dir, "upper_b_list_7_timepoints.csv"),
  check.names = FALSE
))

model_se <- as.matrix(read.csv(
  file.path(output_dir, "model_se_combined_7_timepoints.csv"),
  check.names = FALSE
))

############################################################
## 2. Extract truth and time points
############################################################

times <- landmark_df$time
true_beta <- landmark_df$true_beta
Tn <- length(true_beta)

stopifnot(ncol(beta_hat) == Tn)
stopifnot(ncol(lower_ci) == Tn)
stopifnot(ncol(upper_ci) == Tn)

############################################################
## 3. Keep complete successful simulations
############################################################

complete_rows <- rowSums(is.na(beta_hat)) == 0 &
  rowSums(is.na(lower_ci)) == 0 &
  rowSums(is.na(upper_ci)) == 0

beta_hat_complete <- beta_hat[complete_rows, , drop = FALSE]
lower_ci_complete <- lower_ci[complete_rows, , drop = FALSE]
upper_ci_complete <- upper_ci[complete_rows, , drop = FALSE]
model_se_complete <- model_se[complete_rows, , drop = FALSE]

n_success <- nrow(beta_hat_complete)

############################################################
## 4. Calculate bias, variance, MSE, coverage by time point
############################################################

summary_by_time <- data.frame(
  time_index = seq_len(Tn),
  time = times,
  true_beta = true_beta,
  n_success = n_success,
  
  mean_beta = colMeans(beta_hat_complete, na.rm = TRUE),
  
  bias = colMeans(
    sweep(beta_hat_complete, 2, true_beta, FUN = "-"),
    na.rm = TRUE
  ),
  
  variance = apply(beta_hat_complete, 2, var, na.rm = TRUE),
  
  mse = colMeans(
    sweep(beta_hat_complete, 2, true_beta, FUN = "-")^2,
    na.rm = TRUE
  ),
  
  empirical_sd = apply(beta_hat_complete, 2, sd, na.rm = TRUE),
  
  average_model_se = colMeans(model_se_complete, na.rm = TRUE),
  
  coverage = colMeans(
    sweep(lower_ci_complete, 2, true_beta, FUN = "<=") &
      sweep(upper_ci_complete, 2, true_beta, FUN = ">="),
    na.rm = TRUE
  )
)

############################################################
## 5. Add RMSE and SE calibration
############################################################

summary_by_time <- summary_by_time %>%
  mutate(
    rmse = sqrt(mse),
    se_to_sd_ratio = average_model_se / empirical_sd
  )

############################################################
## 6. Save and print
############################################################

write.csv(
  summary_by_time,
  file.path(output_dir, "summary_bias_variance_mse_coverage_by_time.csv"),
  row.names = FALSE
)

print(summary_by_time)

