## =============================================================================
## Scenario 3 (Section 3.3; Figure 4): sparse site heterogeneity.
##
## K = 50 sites with n_k = 100 (N = 5,000). Site-specific treatment effects
##   beta_treat_k = log(1.5) + theta_k,  theta_k ~ (1 - pi) delta_0 + pi N(0, tau^2),
## pi in {0.1, 0.2, 0.3, 0.5}, tau in {0.1, 0.3, 0.7, 1.0}. The true site effects
## are drawn once per (pi, tau) cell (seed 1234) and held fixed over replicates;
## replicate j generates data with seed 1234 + j.
##
## Estimators of the site-specific treatment effects:
##   Global   : federated pseudo-value GEE (cloglog), same estimate for every site
##   Local    : site-wise Cox model
##   GSURE_soft (Debiased in Figure 4): Global + variance-adaptive soft-thresholded
##              deviation Delta_k = Local_k - Global, threshold tau_hat * sqrt(V(Delta_k)),
##              tau_hat selected by GSURE (R/gsure_debiasing.R; Section 2.4)
##   PWSURE_soft / PWSURE_ridge : prewhitened SURE soft-thresholding / ridge shrinkage
##              (additional comparators)
##
## Output: <results>/M_050_pi_<pi>_tau_<tau>/metrics_by_sim.csv, one row per
## replicate and method with the site-averaged squared errors (mse_*).
## =============================================================================

## =============================================================================
## Setup
## =============================================================================
repo_root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f), "..", "..")) else getwd()
})
setwd(repo_root)
getwd()

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(pseudo)
  library(ggplot2)
  library(pda)
  library(glmtrans)
  library(glmnet)
})

source(file.path("R", "renewable_gee.R"))
source(file.path("R", "federated_km_pseudo.R"))
source(file.path("R", "gsure_debiasing.R"))

## =============================================================================
## PWSURE helpers
## =============================================================================

## ---- Soft-thresholding SURE (Gaussian sequence) ----
choose_c1_pwsure_soft <- function(Z, grid = seq(0.5, 2.5, by = 0.05), m = NULL){
  d <- length(Z)
  if (is.null(m)) m <- d
  ln_m <- log(m)
  
  sure_val <- function(lambda){
    d + sum(pmin(Z^2, lambda^2) - 2 * (abs(Z) <= lambda))
  }
  
  vals  <- sapply(grid, function(c1) sure_val(c1 * sqrt(ln_m)))
  c1hat <- grid[which.min(vals)]
  list(c1 = c1hat, curve = data.frame(c1 = grid, SURE = vals))
}

## ---- Ridge SURE ----
## Model: Z ~ N(eta, I_K)
## Estimator: eta_hat = a Z, a = 1/(1+gamma)
## SURE(gamma) = ||(a-1)Z||^2 + 2K a - K
choose_gamma_ridge_sure <- function(Z, grid = 10^seq(-3, 2, length.out = 60)){
  K <- length(Z)
  sure_ridge <- function(gamma){
    a <- 1 / (1 + gamma)
    sum(((a - 1) * Z)^2) + 2 * K * a - K
  }
  vals <- sapply(grid, sure_ridge)
  ghat <- grid[which.min(vals)]
  list(gamma = ghat, curve = data.frame(gamma = grid, SURE = vals))
}

## ---- Build Sigma_Delta from single-pass sandwich summaries ----
## Σ_Δ ≈ diag(Vloc_k) + Vglob * 11' - (C 1' + 1 C')
##  where C_k ≈ [A^{-1} M_k A^{-T}]_{tt}
build_SigmaDelta_singlepass <- function(A_glob, meat_list, se_loc_treat, treat_idx_glm){
  K <- length(se_loc_treat)
  Ainv <- solve(A_glob)
  
  ## global sandwich variance (treat coordinate)
  Vglob_mat <- Ainv %*% Reduce(`+`, meat_list) %*% t(Ainv)
  Vglob_tt  <- Vglob_mat[treat_idx_glm, treat_idx_glm]
  
  ## approximate cross-cov contributions C_k
  Ck <- sapply(meat_list, function(Mk){
    tmp <- Ainv %*% Mk %*% t(Ainv)
    tmp[treat_idx_glm, treat_idx_glm]
  })
  
  ## local variances from Cox
  Vloc_tt <- se_loc_treat^2
  
  ## Sigma for Delta_k = beta_loc,k - beta_glob
  Sigma <- diag(Vloc_tt) + Vglob_tt * matrix(1, K, K) -
    ( Ck %*% t(rep(1, K)) + rep(1, K) %*% t(Ck) )
  
  ## stabilize PD (eigenvalue truncation)
  eig <- eigen(Sigma, symmetric = TRUE)
  eig$values[eig$values < 1e-10] <- 1e-10
  Sigma <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  
  list(Sigma = Sigma, Vglob_tt = Vglob_tt, Ck = Ck, Vloc_tt = Vloc_tt)
}

## =============================================================================
## Truth generators
## =============================================================================

## ---- Sparse-mixture heterogeneity ----
## theta_k ~ (1-pi) * 0 + pi * N(0, tau^2)
## beta_treat_k = log(hr_treat) + theta_k
init_taskwise_layout_sparse <- function(
    M               = 10,
    pi              = 0.10,
    tau             = 0.50,
    hr_treat        = 1.5,
    beta_conf_fixed = c(0.30, -0.20),
    seed_layout     = 1234
){
  set.seed(seed_layout)
  beta0 <- log(hr_treat)
  
  is_het <- rbinom(M, 1, pi)
  theta  <- rep(0, M)
  if (sum(is_het == 1) > 0) {
    theta[is_het == 1] <- rnorm(sum(is_het == 1), mean = 0, sd = tau)
  }
  
  beta_treat_draw <- beta0 + theta
  
  per_site <- lapply(seq_len(M), function(m) {
    list(beta_m = c(beta_treat_draw[m], beta_conf_fixed))
  })
  
  truth_df <- data.frame(
    site       = 1:M,
    beta_treat = beta_treat_draw,
    theta      = theta,
    is_het     = is_het,
    beta_X1    = beta_conf_fixed[1],
    beta_X2    = beta_conf_fixed[2]
  )
  
  list(M = M, d = 3L, pi = pi, tau = tau, hr_treat = hr_treat,
       per_site = per_site, truth_df = truth_df)
}

## ---- Equal-N-per-site survival data generator ----
make_surv_from_layout_equalN <- function(
    layout,
    N_per_site   = 500,
    lambda0      = 1/20,
    censor_rate  = -log(0.4)/15,
    admin_censor = 40,
    seed_data    = NULL
){
  if (!is.null(seed_data)) set.seed(seed_data)
  M <- layout$M
  rows <- vector("list", M)
  
  for (m in 1:M) {
    n_m <- N_per_site
    X1  <- rbinom(n_m, 1, 0.5)
    X2  <- rnorm(n_m)
    pA  <- plogis(0.4*X1 + 0.3*X2)
    A   <- rbinom(n_m, 1, pA)
    Z   <- cbind(A, X1, X2)
    
    beta_m <- layout$per_site[[m]]$beta_m
    rate    <- lambda0 * exp(as.numeric(Z %*% beta_m))
    T_event <- -log(runif(n_m)) / rate
    T_cens  <- pmin(rexp(n_m, rate = censor_rate), admin_censor)
    
    time  <- pmin(T_event, T_cens)
    delta <- as.integer(T_event <= T_cens)
    
    rows[[m]] <- data.frame(site = m, time, delta, treat = A, X1, X2)
  }
  
  dat <- dplyr::bind_rows(rows)
  list(dat = dat, truth_df = layout$truth_df)
}

## =============================================================================
## Experiment grid and I/O
## =============================================================================
nsim           <- 500
start_sim      <- 1
M_grid         <- 50       # e.g. c(20, 50, 100)
pi_grid        <- c(0.1, 0.2, 0.3, 0.5)
tau_grid       <- c(0.1, 0.3, 0.7, 1.0)
## Optional environment overrides (e.g., to run one (pi, tau) cell per job):
##   FEDSURV_PI="0.1"  FEDSURV_TAU="0.3,0.7"  FEDSURV_SIMS="1:500"
if (nzchar(Sys.getenv("FEDSURV_PI")))  pi_grid  <- as.numeric(strsplit(Sys.getenv("FEDSURV_PI"),  ",")[[1]])
if (nzchar(Sys.getenv("FEDSURV_TAU"))) tau_grid <- as.numeric(strsplit(Sys.getenv("FEDSURV_TAU"), ",")[[1]])
if (nzchar(Sys.getenv("FEDSURV_SIMS"))) { .r <- as.integer(strsplit(Sys.getenv("FEDSURV_SIMS"), ":")[[1]]); start_sim <- .r[1]; nsim <- .r[2] }
N_total_site   <- 5000
c1_grid_fixed  <- NULL  # baseline fixed c1

debug_survive  <- FALSE

## SURE grids
c1_sure_grid_soft <- seq(0.50, 2.0, by = 0.05)
gamma_sure_grid   <- 10^seq(-3, 2, length.out = 60)

base_outdir <- Sys.getenv("FEDSURV_OUTPUT_DIR", file.path("results", "scenario3_sparse_heterogeneity"))
dir.create(base_outdir, recursive = TRUE, showWarnings = FALSE)

## =============================================================================
## Landmark times for the pseudo-observations (10th-80th KM percentiles)
## =============================================================================
global_n     <- 100000
lambda0      <- 1/20
censor_rate  <- -log(0.4)/15
admin_censor <- 40

set.seed(as.integer(Sys.getenv("FEDSURV_LANDMARK_SEED", "20260301")))  # same landmark times in every run
X1  <- rbinom(global_n, 1, 0.5)
X2  <- rnorm(global_n)
pA  <- plogis(0.4*X1 + 0.3*X2)
A   <- rbinom(global_n, 1, pA)
Z   <- cbind(A, X1, X2)
beta <- c(log(1.5), 0.30, -0.20)

rate    <- lambda0 * exp(as.numeric(Z %*% beta))
T_event <- -log(runif(global_n)) / rate
T_cens  <- pmin(rexp(global_n, rate = censor_rate), admin_censor)
time    <- pmin(T_event, T_cens)
delta   <- as.integer(T_event <= T_cens)

global_curve_event <- survfit(Surv(time, delta) ~ 1)
quantile_time_new  <- quantile(global_curve_event, probs = seq(0.1, 0.8, 0.1))[[1]]

## =============================================================================
## Main simulations
## =============================================================================
for (M in M_grid) {
  for (pi in pi_grid) {
    for (tau in tau_grid) {
      
      ## sanity checks (avoid tau being overwritten to 0)
      stopifnot(is.numeric(M), is.numeric(pi), is.numeric(tau))
      stopifnot(M >= 2, pi >= 0, pi <= 1, tau >= 0)
      
      ## ----- fix true local betas ONCE per (M, pi, tau) -----
      layout <- init_taskwise_layout_sparse(
        M = M, pi = pi, tau = tau, hr_treat = 1.5, beta_conf_fixed = c(0.30, -0.20),
        seed_layout = 1234
      )
      
      out_dir  <- file.path(base_outdir, sprintf("M_%03d_pi_%.2f_tau_%.2f", M, pi, tau))
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      out_file <- file.path(out_dir, "metrics_by_sim.csv")
      
      ## write header if not exists
      if (!file.exists(out_file)) {
        write.table(
          data.frame(
            sim = integer(),
            M = integer(),
            pi = numeric(),
            tau = numeric(),
            n_per_site = integer(),
            n_het_true = integer(),
            method = character(),
            c1_soft = numeric(),     # only for PWSURE_soft
            lambda_soft = numeric(), # only for PWSURE_soft
            gamma_ridge = numeric(), # only for PWSURE_ridge
            a_ridge = numeric(),     # only for PWSURE_ridge
            c1_fixed = numeric(),    # only for fixed baseline
            survive_frac = numeric(),
            mean_abs_debiased = numeric(),
            mean_abs_global   = numeric(),
            mean_abs_local    = numeric(),
            mse_debiased      = numeric(),
            mse_global        = numeric(),
            mse_local         = numeric(),
            tau_gsure         = numeric()  # only for GSURE_soft
          ),
          out_file, sep = ",", row.names = FALSE, col.names = TRUE
        )
      }
      
      set.seed(1234)
      
      for (j in start_sim:nsim) {
        ## a replicate that fails (e.g. a numerically singular fit) is skipped; the remaining replicates continue
        rep_ok <- tryCatch({
        
          ## --- Generate equal-N data with the FIXED truth ---
          n_per_site <- as.integer(N_total_site / M)
          sim_fixed <- make_surv_from_layout_equalN(
            layout, N_per_site = n_per_site, seed_data = 1234 + j,
            lambda0 = lambda0, censor_rate = censor_rate, admin_censor = admin_censor
          )
          true_local_beta <- sim_fixed$truth_df
          sim_data        <- sim_fixed$dat
          sim_data_ori    <- sim_data
        
          truth_treat <- as.numeric(true_local_beta$beta_treat)   # beta_0^(k)
          K           <- length(truth_treat)
          n_het_true  <- sum(true_local_beta$is_het == 1)
        
          ## ====== Federated KM splines and pseudo-observations ======
          first_site_data     <- dplyr::filter(sim_data, site == 1)
          first_site_estimate <- get_survival_splines_cloglog(
            first_site_data, degree = 2, knots = 5,
            bknots = c(0.1, 39.9), risk_set_minimum = 0.05
          )
          updated_splines <- first_site_estimate
        
          max_knots <- 7
          current_knots <- 5
          chunksize <- 5
          counter <- 0
        
          for (i in 2:length(unique(sim_data$site))) {
            site_data <- dplyr::filter(sim_data, site == i)
            if (nrow(site_data) != 0) {
              n_batches <- nrow(site_data) / chunksize
              knots_list <- NULL
            
              if (nrow(site_data) > 100) {
                for (k2 in seq(1, nrow(site_data)/100)) {
                  knots_list <- c(
                    knots_list,
                    c(rep(current_knots, 100/chunksize),
                      rep(current_knots + k2, n_batches - 100/chunksize))
                  )
                }
                knots_list <- ifelse(knots_list > max_knots, max_knots, knots_list)
              } else {
                if (counter > 200 & current_knots < max_knots) {
                  knots_list <- rep(current_knots + 1, n_batches)
                  counter <- 0
                } else {
                  knots_list <- rep(current_knots, n_batches)
                }
                knots_list <- ifelse(knots_list > max_knots, max_knots, knots_list)
                counter <- counter + nrow(site_data)
              }
            
              updated_splines <- update_survival_splines_chunk(
                updated_splines, breakpoints = 10,
                time = site_data$time, delta = site_data$delta,
                chunksize = chunksize, knots_list = knots_list,
                loud = FALSE, restriction = 0.2, risk_floor = 0.005
              )
              current_knots <- tail(knots_list, 1)
            }
          }
        
          PS_fun <- function(qt) 1 - theoretical_KM_PS(
            updated_splines, qt, sim_data$time, sim_data$delta, restriction = 2
          )
          Pseudo_list <- lapply(quantile_time_new, PS_fun)
        
          len0     <- nrow(sim_data)
          sim_data <- do.call(rbind, rep(list(sim_data), length(quantile_time_new)))
          sim_data$txx <- factor(rep(quantile_time_new, each = len0))
          sim_data$Pseudo_Spline <- unlist(Pseudo_list)
        
          ## ====== Federated GLM (RenewGLM_out) ======
          build_X <- function(dd) {
            X <- cbind(1, dd$treat, dd$X1, dd$X2, factor(dd$txx))
            X <- as.data.frame(X)
            X[, 5] <- as.factor(X[, 5])
            enc <- model.matrix(~ V5 - 1, data = X)
            enc <- enc[, -1]
            as.matrix(cbind(X[, c(1, 2, 3, 4)], enc))
          }
        
          y1  <- sim_data$Pseudo_Spline[sim_data$site == 1]
          X1m <- build_X(sim_data[sim_data$site == 1, ])
          p_glm <- ncol(X1m)
          betahat  <- rep(0, p_glm)
          infomats <- diag(0, p_glm, p_glm)
        
          id <- rep(1:(nrow(X1m)/length(unique(sim_data$txx))),
                    length(unique(sim_data$txx)))
        
          summary <- RenewGLM_out(
            X1m, y1, id, "cloglog", betahat, infomats, list(),
            intercept = FALSE, s = NULL, phi = NULL
          )
          betahat   <- summary[[1]]
          infomats  <- summary[[2]]
          meat_list <- list(summary[[3]])
        
          for (i in 2:length(unique(sim_data$site))) {
            yi <- sim_data$Pseudo_Spline[sim_data$site == i]
            Xi <- build_X(sim_data[sim_data$site == i, ])
            id <- rep(1:(nrow(Xi)/length(unique(sim_data$txx))),
                      length(unique(sim_data$txx)))
          
            summary <- RenewGLM_out(
              Xi, yi, id, "cloglog", betahat, infomats, meat_list,
              intercept = FALSE, s = NULL, phi = NULL
            )
            betahat        <- summary[[1]]
            infomats       <- summary[[2]]
            meat_list[[i]] <- summary[[3]]
          }
        
          betahat_global <- betahat[2:4]  # treat, X1, X2
          treat_idx_glm  <- 2
        
          ## ====== Local site-wise Cox ======
          cox_sim_beta_ind <- NULL
          cox_sim_se_ind   <- NULL
          for (s in unique(sim_data_ori$site)) {
            dat_s <- subset(sim_data_ori, site == s)
            cox_model <- coxph(Surv(time, delta) ~ treat + X1 + X2, data = dat_s)
            cox_sim_beta_ind <- rbind(cox_sim_beta_ind, coef(cox_model))
            cox_sim_se_ind   <- rbind(cox_sim_se_ind, sqrt(diag(vcov(cox_model))))
          }
          colnames(cox_sim_beta_ind) <- c("beta_treat","beta_X1","beta_X2")
          colnames(cox_sim_se_ind)   <- c("beta_treat","beta_X1","beta_X2")
        
          beta_loc_treat  <- cox_sim_beta_ind[, "beta_treat"]
          se_loc_treat    <- cox_sim_se_ind[, "beta_treat"]
          diff_site_treat <- beta_loc_treat - betahat_global[1]  # Delta_k
        
          ## ========== Baseline: fixed c1, sitewise soft-thres using local SE only ==========
          for (c1_val in c1_grid_fixed) {
            k_fac <- sqrt(log(K))
            tau_site_treat <- c1_val * k_fac * se_loc_treat
          
            deb_treat <- betahat_global[1] + sign(diff_site_treat) *
              pmax(abs(diff_site_treat) - tau_site_treat, 0)
          
            err_deb   <- deb_treat                 - truth_treat
            err_glob  <- rep(betahat_global[1], K) - truth_treat
            err_local <- beta_loc_treat            - truth_treat
          
            survive_frac <- mean(abs(diff_site_treat) > tau_site_treat)
          
            row_out <- data.frame(
              sim = j, M = M, pi = pi, tau = tau,
              n_per_site = n_per_site,
              n_het_true = n_het_true,
              method = "fixed_soft_localSE",
              c1_soft = NA_real_,
              lambda_soft = NA_real_,
              gamma_ridge = NA_real_,
              a_ridge = NA_real_,
              c1_fixed = c1_val,
              survive_frac = survive_frac,
              mean_abs_debiased = mean(abs(err_deb)),
              mean_abs_global   = mean(abs(err_glob)),
              mean_abs_local    = mean(abs(err_local)),
              mse_debiased      = mean(err_deb^2),
              mse_global        = mean(err_glob^2),
              mse_local         = mean(err_local^2),
              tau_gsure         = NA_real_
            )
            write.table(row_out, out_file, sep = ",",
                        row.names = FALSE, col.names = FALSE, append = TRUE)
          }
        
          ## ========== Build Sigma_Delta and whiten ==========
          SDobj <- build_SigmaDelta_singlepass(
            A_glob        = infomats,
            meat_list     = meat_list,
            se_loc_treat  = se_loc_treat,
            treat_idx_glm = treat_idx_glm
          )
          Sigma_Delta <- SDobj$Sigma
        
          ## Whitening: Sigma = L %*% t(L) (L lower)
          L <- t(chol(Sigma_Delta))
          ZW <- solve(L, diff_site_treat)  # Z = L^{-1} Delta
        
          ## ========== PWSURE_soft ==========
          pw_fit <- choose_c1_pwsure_soft(ZW, grid = c1_sure_grid_soft, m = K)
          c1_pw  <- as.numeric(pw_fit$c1)
          lambda <- c1_pw * sqrt(log(K))
        
          eta_hat_soft   <- sign(ZW) * pmax(abs(ZW) - lambda, 0)
          theta_hat_soft <- as.vector(L %*% eta_hat_soft)
          deb_treat_soft <- betahat_global[1] + theta_hat_soft
        
          err_deb_soft <- deb_treat_soft           - truth_treat
          err_glob     <- rep(betahat_global[1], K) - truth_treat
          err_local    <- beta_loc_treat           - truth_treat
        
          survive_frac_soft <- mean(theta_hat_soft != 0)
        
          row_out_soft <- data.frame(
            sim = j, M = M, pi = pi, tau = tau,
            n_per_site = n_per_site,
            n_het_true = n_het_true,
            method = "PWSURE_soft",
            c1_soft = c1_pw,
            lambda_soft = lambda,
            gamma_ridge = NA_real_,
            a_ridge = NA_real_,
            c1_fixed = NA_real_,
            survive_frac = survive_frac_soft,
            mean_abs_debiased = mean(abs(err_deb_soft)),
            mean_abs_global   = mean(abs(err_glob)),
            mean_abs_local    = mean(abs(err_local)),
            mse_debiased      = mean(err_deb_soft^2),
            mse_global        = mean(err_glob^2),
            mse_local         = mean(err_local^2),
            tau_gsure         = NA_real_
          )
          write.table(row_out_soft, out_file, sep = ",",
                      row.names = FALSE, col.names = FALSE, append = TRUE)
        
          ## ========== PWSURE_ridge ==========
          ridge_fit <- choose_gamma_ridge_sure(ZW, grid = gamma_sure_grid)
          gamma_hat <- as.numeric(ridge_fit$gamma)
          a_hat     <- 1 / (1 + gamma_hat)
        
          eta_hat_ridge   <- a_hat * ZW
          theta_hat_ridge <- as.vector(L %*% eta_hat_ridge)
          deb_treat_ridge <- betahat_global[1] + theta_hat_ridge
        
          err_deb_ridge <- deb_treat_ridge - truth_treat
        
          ## ridge shrinkage gives no exact zeros: report the fraction with |theta_hat| > 1e-12
          survive_frac_ridge <- mean(abs(theta_hat_ridge) > 1e-12)
        
          row_out_ridge <- data.frame(
            sim = j, M = M, pi = pi, tau = tau,
            n_per_site = n_per_site,
            n_het_true = n_het_true,
            method = "PWSURE_ridge",
            c1_soft = NA_real_,
            lambda_soft = NA_real_,
            gamma_ridge = gamma_hat,
            a_ridge = a_hat,
            c1_fixed = NA_real_,
            survive_frac = survive_frac_ridge,
            mean_abs_debiased = mean(abs(err_deb_ridge)),
            mean_abs_global   = mean(abs(err_glob)),
            mean_abs_local    = mean(abs(err_local)),
            mse_debiased      = mean(err_deb_ridge^2),
            mse_global        = mean(err_glob^2),
            mse_local         = mean(err_local^2),
            tau_gsure         = NA_real_
          )
          write.table(row_out_ridge, out_file, sep = ",",
                      row.names = FALSE, col.names = FALSE, append = TRUE)
        
          ## ========== GSURE_soft (Section 2.4; "Debiased" in Figure 4) ==========
          SD_g       <- build_SigmaDelta_singlepass_coord_safe(infomats, meat_list, se_loc_treat, treat_idx_glm,
                                                               eig_floor = 1e-10, ridge_floor = 1e-10)
          SigmaDelta <- SD_g$Sigma
          fit_tau    <- choose_tau_euclidean_gsure_scaled(Delta = diff_site_treat, Sigma = SigmaDelta,
                                                          n_grid = 80, tau_min = 0, tau_max = NULL)
          tau_hat    <- fit_tau$tau
          se_delta   <- sqrt(diag(SigmaDelta))
          thresh_vec <- tau_hat * se_delta
          hhat       <- sign(diff_site_treat) * pmax(abs(diff_site_treat) - thresh_vec, 0)
          deb_treat_gsure <- betahat_global[1] + hhat
        
          err_deb_gsure <- deb_treat_gsure - truth_treat
          row_out_gsure <- data.frame(
            sim = j, M = M, pi = pi, tau = tau,
            n_per_site = n_per_site,
            n_het_true = n_het_true,
            method = "GSURE_soft",
            c1_soft = NA_real_,
            lambda_soft = NA_real_,
            gamma_ridge = NA_real_,
            a_ridge = NA_real_,
            c1_fixed = NA_real_,
            survive_frac = mean(abs(diff_site_treat) > thresh_vec),
            mean_abs_debiased = mean(abs(err_deb_gsure)),
            mean_abs_global   = mean(abs(err_glob)),
            mean_abs_local    = mean(abs(err_local)),
            mse_debiased      = mean(err_deb_gsure^2),
            mse_global        = mean(err_glob^2),
            mse_local         = mean(err_local^2),
            tau_gsure         = tau_hat
          )
          write.table(row_out_gsure, out_file, sep = ",",
                      row.names = FALSE, col.names = FALSE, append = TRUE)
        
          cat(sprintf(
            "[M=%d, pi=%.2f, tau=%.2f] sim %d saved: soft(c1=%.2f, lam=%.2f, surv=%.3f), ridge(gamma=%.3g, a=%.3f).\n",
            M, pi, tau, j, c1_pw, lambda, survive_frac_soft, gamma_hat, a_hat
          ))
          TRUE
        }, error = function(e) {
          message(sprintf("[M=%d, pi=%.2f, tau=%.2f] sim %d skipped: %s", M, pi, tau, j, conditionMessage(e)))
          FALSE
        })
      }
      
      cat(sprintf("Done scenario: M=%d, pi=%.2f, tau=%.2f -> %s\n", M, pi, tau, out_file))
    }
  }
}