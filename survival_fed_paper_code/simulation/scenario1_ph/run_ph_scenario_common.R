## =============================================================================
## Scenario 1 (Section 3.1; Figures 1-2): simulation engine for the PH setting.
##
## For each replicate: generate data (exponential event times, X ~ N(0,1),
## A ~ Bernoulli(expit(0.5 X))), split them into sites, and fit
##   * the pooled Cox model,
##   * ODAC (R/odac_original_0503.R),
##   * the federated pseudo-value regression: sequential federated Kaplan-Meier
##     splines, influence-function pseudo-observations at the landmark times and a
##     renewable GEE with cloglog link and sandwich variance
##     (R/PSFederated_Essential_New.R, R/new_renewglm_out.R).
## Replicate j uses set.seed(config$seed + j). Per-replicate results are cached in
## <results>/<scenario_name>/per_sim_rds/, so an interrupted run can be restarted.
## Settings are passed through `config` by the run_scenario*.R scripts.
## =============================================================================

fedsurv_src_dir <- function() getOption("fedsurv.src_dir", "R")
fedsurv_results_root <- function() getOption("fedsurv.results_root",
  Sys.getenv("FEDSURV_OUTPUT_DIR", file.path("results", "scenario1_ph")))

load_required_packages <- function() {
  suppressPackageStartupMessages({
    library(survival)
    library(dplyr)
    library(pseudo)
    library(ggplot2)
    library(jsonlite)
    library(parallel)
    library(data.table)
  })
  
  if (requireNamespace("pda", quietly = TRUE)) {
    suppressPackageStartupMessages(library(pda))
  }
  
  odac_source_file <- file.path(fedsurv_src_dir(), c("odac_original_0503.R", "odac_source_original.R"))
  odac_source_file <- odac_source_file[file.exists(odac_source_file)]
  if (length(odac_source_file) > 0) {
    source(odac_source_file[1])
  }

  source(file.path(fedsurv_src_dir(), "new_renewglm_out.R"))
  source(file.path(fedsurv_src_dir(), "PSFederated_Essential_New.R"))
  
  needed_odac <- c("pda", "getCloudConfig", "pdaGet")
  missing_odac <- needed_odac[!vapply(needed_odac, exists, logical(1), mode = "function")]
  if (length(missing_odac) > 0) {
    stop(
      "ODAC/pda functions are missing: ",
      paste(missing_odac, collapse = ", "),
      ". Install/load the pda package or place odac_original_0503.R or odac_source_original.R in the working directory."
    )
  }
}

make_sim_data <- function(config, sim_id) {
  set.seed(config$seed + sim_id)
  
  Ns <- config$Ns
  N <- sum(Ns)
  Site <- rep(seq_along(Ns), Ns)
  Site_chen <- paste0("site", Site)
  
  X <- rnorm(N, 0, 1)
  Treatment <- rbinom(N, 1, plogis(0.5 * X))
  
  rate_T <- config$base_hazard *
    exp(X * log(0.7)) *
    exp(Treatment * config$true_treatment_beta)
  
  Survival_Time <- rexp(N, rate = rate_T)
  Censoring_Time <- rexp(N, rate = config$censoring_rate)
  Censoring_Time <- pmin(Censoring_Time, config$admin_censor)
  
  Observed_Time <- pmin(Survival_Time, Censoring_Time)
  if (isTRUE(config$round_observed_time)) {
    Observed_Time <- round(Observed_Time, 4)
  }
  
  Event <- as.integer(Survival_Time <= Censoring_Time)
  
  sim_data <- data.frame(
    Site = Site,
    Treatment = Treatment,
    X = X,
    Survival_Time = Survival_Time,
    Censoring_Time = Censoring_Time,
    delta = Event,
    time = Observed_Time
  )
  
  sim_data_chen <- data.frame(
    site = Site_chen,
    time = Observed_Time,
    status = Event,
    Treatment = Treatment,
    X = X
  )
  
  list(
    sim_data = sim_data,
    sim_data_chen = sim_data_chen,
    event_rate = mean(Event),
    treatment_rate = mean(Treatment)
  )
}

run_cox <- function(sim_data) {
  fit <- coxph(Surv(time, delta) ~ Treatment + X, data = sim_data)
  list(
    beta = summary(fit)$coefficients[, "coef"],
    se = summary(fit)$coefficients[, "se(coef)"]
  )
}

run_odac <- function(sim_data_chen, sim_id, config) {
  if (isTRUE(config$round_odac_time)) {
    sim_data_chen$time <- round(sim_data_chen$time, 4)
  }
  
  data_chen_split <- split(sim_data_chen, sim_data_chen$site)
  sites <- unique(sim_data_chen$site)
  
  mydir <- file.path(config$pda_work_dir, sprintf("sim_%04d", sim_id))
  dir.create(mydir, recursive = TRUE, showWarnings = FALSE)
  
  control <- list(
    project_name = paste0("PseudoRegRevision_", config$scenario_name),
    step = "initialize",
    sites = sites,
    heterogeneity = FALSE,
    model = "ODAC",
    family = "cox",
    outcome = "Surv(time, status)",
    variables = c("Treatment", "X"),
    optim_maxit = 1000,
    lead_site = "site1",
    upload_date = as.character(Sys.time())
  )
  
  pda(site_id = "site1", control = control, dir = mydir)
  
  for (site_id in rev(sites)) {
    pda(site_id = site_id, ipdata = data_chen_split[[site_id]], dir = mydir)
    
    if (isTRUE(config$drop_zero_var_odac_sites)) {
      init_path <- file.path(mydir, paste0(site_id, "_initialize.json"))
      if (file.exists(init_path)) {
        init_json <- jsonlite::fromJSON(init_path)
        var_i <- unlist(init_json$Vhat_i)
        beta_i <- unlist(init_json$bhat_i)
        should_drop <- (any(var_i == 0, na.rm = TRUE) && !anyNA(beta_i)) || anyNA(var_i)
        if (isTRUE(should_drop)) {
          init_json$bhat_i <- NA
          init_json$Vhat_i <- NA
          jsonlite::write_json(init_json, path = init_path, pretty = TRUE, auto_unbox = TRUE)
        }
      }
    }
  }
  
  for (site_id in rev(sites)) {
    pda(site_id = site_id, ipdata = data_chen_split[[site_id]], dir = mydir)
  }
  
  for (site_id in rev(sites)) {
    pda(site_id = site_id, ipdata = data_chen_split[[site_id]], dir = mydir)
  }
  
  pda(site_id = "site1", ipdata = data_chen_split[["site1"]], dir = mydir)
  
  config_cloud <- getCloudConfig(site_id = "site1", dir = mydir)
  fit_pda <- pdaGet(name = "site1_estimate", config = config_cloud)
  
  beta <- fit_pda$btilde
  se <- sqrt(diag(solve(fit_pda$Htilde) / nrow(sim_data_chen)))
  
  list(beta = beta, se = se)
}

make_bknots <- function(sim_data, config) {
  if (config$bknots_type == "fixed_40") {
    return(c(0.1, 39.9))
  }
  if (config$bknots_type == "fixed_20") {
    return(c(0.1, 19.9))
  }
  if (config$bknots_type == "first_site_max") {
    first_site <- dplyr::filter(sim_data, Site == 1)
    upper <- max(first_site$time, na.rm = TRUE) - 0.1
    if (!is.finite(upper) || upper <= 0.1) {
      stop("Invalid first-site max time for spline boundary knots.")
    }
    return(c(0.1, upper))
  }
  stop("Unknown bknots_type: ", config$bknots_type)
}

update_splines_all_sites <- function(sim_data, Ns, config) {
  first_site_data <- dplyr::filter(sim_data, Site == 1)
  spline_args <- list(
    data = first_site_data,
    degree = 2,
    knots = 5,
    bknots = make_bknots(sim_data, config)
  )
  if (!is.null(config$risk_set_minimum)) {
    spline_args$risk_set_minimum <- config$risk_set_minimum
  }
  
  updated_splines <- do.call(get_survival_splines_cloglog, spline_args)
  
  max_knots <- 7
  current_knots <- 5
  chunksize <- 5
  counter <- 0
  
  for (site_index in seq(2, length(Ns))) {
    site_data <- dplyr::filter(sim_data, Site == site_index)
    if (nrow(site_data) == 0) next
    
    n_batches <- nrow(site_data) / chunksize
    knots_list <- NULL
    
    if (nrow(site_data) > 100) {
      for (k in seq(1, nrow(site_data) / 100)) {
        knots_list <- c(
          knots_list,
          c(
            rep(current_knots, 100 / chunksize),
            rep(current_knots + k, ceiling(n_batches) - 100 / chunksize)
          )
        )
      }
      knots_list <- ifelse(knots_list > max_knots, max_knots, knots_list)
    } else {
      if (counter > 200 && current_knots < max_knots) {
        knots_list <- rep(current_knots + 1, ceiling(n_batches))
        counter <- 0
      } else {
        knots_list <- rep(current_knots, ceiling(n_batches))
      }
      knots_list <- ifelse(knots_list > max_knots, max_knots, knots_list)
      counter <- counter + nrow(site_data)
    }
    
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
  
  updated_splines
}

get_spline_survival_function <- function(updated_splines) {
  function(t) {
    survival_function_cloglog.v(
      t,
      updated_splines$Event$Knots,
      updated_splines$Event$Degree,
      updated_splines$Event$BKnots,
      updated_splines$Event$Coefficients
    )
  }
}

quantile_from_survival <- function(surv_func, grid_points, q) {
  st_values <- surv_func(grid_points)
  times <- vapply(q, function(qi) {
    candidate <- grid_points[st_values > qi]
    if (length(candidate) == 0) return(NA_real_)
    max(candidate)
  }, numeric(1))
  
  if (any(!is.finite(times))) {
    stop("Could not determine all pseudo-observation landmark times.")
  }
  times
}

choose_landmark_times <- function(final_ST, sim_data, config) {
  if (isTRUE(config$use_event_time_grid_for_landmarks)) {
    grid_max <- max(sim_data$time[sim_data$delta == 1], na.rm = TRUE)
  } else {
    grid_max <- max(sim_data$time, na.rm = TRUE)
  }
  grid_points <- seq(0, grid_max, by = 0.1)
  
  if (config$landmark_style == "scenario3") {
    qt <- quantile_from_survival(
      final_ST,
      grid_points,
      c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    )
    qt <- unique(qt)
    qt <- rev(qt)
    qt <- qt[-length(qt)]
    return(seq(min(qt), max(qt), length.out = 8))
  }
  
  qt <- quantile_from_survival(final_ST, grid_points, c(0.2, 0.9))
  qt <- rev(qt)
  seq(qt[1], qt[2], length.out = 8)
}

make_design_matrix <- function(site_long_data) {
  xmat <- cbind(
    1,
    site_long_data$Treatment,
    site_long_data$X,
    as.factor(site_long_data$txx)
  )
  xmat <- as.data.frame(xmat)
  xmat[, 4] <- as.factor(xmat[, 4])
  
  encoded <- model.matrix(~ V4 - 1, data = xmat)
  encoded <- encoded[, -1, drop = FALSE]
  xmat <- cbind(xmat[, c(1, 2, 3)], encoded)
  as.matrix(xmat)
}

run_pseudo_reg <- function(sim_data, config) {
  Ns <- as.integer(table(sim_data$Site))
  updated_splines <- update_splines_all_sites(sim_data, Ns, config)
  final_ST <- get_spline_survival_function(updated_splines)
  quantile_time <- choose_landmark_times(final_ST, sim_data, config)
  
  pseudo_list <- lapply(
    quantile_time,
    function(t_ps) {
      theoretical_KM_PS(
        spline_params = updated_splines,
        t_ps = t_ps,
        time = sim_data$time,
        delta = sim_data$delta,
        restriction = 2
      )
    }
  )
  
  n0 <- nrow(sim_data)
  sim_data_long <- do.call(rbind, replicate(length(quantile_time), sim_data, simplify = FALSE))
  sim_data_long$txx <- as.factor(rep(quantile_time, each = n0))
  sim_data_long$Pseudo_Spline <- unlist(lapply(pseudo_list, function(x) 1 - x))
  
  site1_data <- dplyr::filter(sim_data_long, Site == 1)
  y <- site1_data$Pseudo_Spline
  xmat <- make_design_matrix(site1_data)
  
  p <- ncol(xmat)
  infomats <- diag(0, p, p)
  betahat <- rep(0, p)
  id <- rep(seq_len(Ns[1]), length(quantile_time))
  meat_list <- list()
  
  gee_summary <- RenewGLM_out(
    xmat, y, id, "cloglog",
    betahat, infomats, meat_list,
    intercept = FALSE, s = NULL, phi = NULL
  )
  
  betahat <- gee_summary[[1]]
  infomats <- gee_summary[[2]]
  meat_list[[1]] <- gee_summary[[3]]
  
  for (site_index in seq(2, length(Ns))) {
    site_data <- dplyr::filter(sim_data_long, Site == site_index)
    if (nrow(site_data) == 0) next
    
    y <- site_data$Pseudo_Spline
    xmat <- make_design_matrix(site_data)
    id <- rep(seq_len(Ns[site_index]), length(quantile_time))
    
    gee_summary <- tryCatch(
      RenewGLM_out(
        xmat, y, id, "cloglog",
        betahat, infomats, meat_list,
        intercept = FALSE, s = NULL, phi = NULL
      ),
      error = function(e) {
        warning(paste0(site_index, "th site excluded: ", conditionMessage(e)))
        NULL
      }
    )
    
    if (!is.null(gee_summary)) {
      betahat <- gee_summary[[1]]
      infomats <- gee_summary[[2]]
      meat_list[[site_index]] <- gee_summary[[3]]
    }
  }
  
  meat_nonnull <- meat_list[!vapply(meat_list, is.null, logical(1))]
  global_meat <- Reduce("+", meat_nonnull)
  se <- sqrt(diag(solve(infomats) %*% global_meat %*% solve(infomats)))
  
  list(
    beta = betahat,
    se = se,
    landmark_times = quantile_time
  )
}

run_quietly <- function(expr, log_file, quiet = TRUE) {
  if (!isTRUE(quiet)) {
    return(force(expr))
  }
  
  dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
  con <- file(log_file, open = "at")
  sink(con, type = "output")
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(con)
  }, add = TRUE)
  
  force(expr)
}

empty_result_row <- function(sim_id, status, error_message = NA_character_) {
  out <- data.frame(
    sim_id = sim_id,
    status = status,
    error_message = error_message,
    cox_error = NA_character_,
    odac_error = NA_character_,
    pseudo_reg_error = NA_character_,
    event_rate = NA_real_,
    treatment_rate = NA_real_,
    cox_beta_Treatment = NA_real_,
    cox_beta_X = NA_real_,
    cox_se_Treatment = NA_real_,
    cox_se_X = NA_real_,
    odac_beta_Treatment = NA_real_,
    odac_beta_X = NA_real_,
    odac_se_Treatment = NA_real_,
    odac_se_X = NA_real_,
    if_beta_Treatment = NA_real_,
    if_beta_X = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (nm in c("Intercept", "Treatment", "X", paste0("Time", 2:8))) {
    out[[paste0("if_se_", nm)]] <- NA_real_
  }
  for (j in seq_len(8)) {
    out[[paste0("landmark_time_", j)]] <- NA_real_
  }
  out
}

result_row <- function(sim_id, dat, cox_fit, odac_fit, if_fit) {
  out <- empty_result_row(sim_id, "ok", NA_character_)
  out$event_rate <- dat$event_rate
  out$treatment_rate <- dat$treatment_rate
  
  out$cox_beta_Treatment <- unname(cox_fit$beta[1])
  out$cox_beta_X <- unname(cox_fit$beta[2])
  out$cox_se_Treatment <- unname(cox_fit$se[1])
  out$cox_se_X <- unname(cox_fit$se[2])
  
  out$odac_beta_Treatment <- unname(odac_fit$beta[1])
  out$odac_beta_X <- unname(odac_fit$beta[2])
  out$odac_se_Treatment <- unname(odac_fit$se[1])
  out$odac_se_X <- unname(odac_fit$se[2])
  
  out$if_beta_Treatment <- unname(if_fit$beta[2])
  out$if_beta_X <- unname(if_fit$beta[3])
  
  if_se_names <- c("Intercept", "Treatment", "X", paste0("Time", 2:8))
  for (idx in seq_along(if_se_names)) {
    out[[paste0("if_se_", if_se_names[idx])]] <- unname(if_fit$se[idx])
  }
  for (idx in seq_along(if_fit$landmark_times)) {
    out[[paste0("landmark_time_", idx)]] <- unname(if_fit$landmark_times[idx])
  }
  out
}

fill_cox_result <- function(out, cox_fit) {
  out$cox_beta_Treatment <- unname(cox_fit$beta[1])
  out$cox_beta_X <- unname(cox_fit$beta[2])
  out$cox_se_Treatment <- unname(cox_fit$se[1])
  out$cox_se_X <- unname(cox_fit$se[2])
  out
}

fill_odac_result <- function(out, odac_fit) {
  out$odac_beta_Treatment <- unname(odac_fit$beta[1])
  out$odac_beta_X <- unname(odac_fit$beta[2])
  out$odac_se_Treatment <- unname(odac_fit$se[1])
  out$odac_se_X <- unname(odac_fit$se[2])
  out
}

fill_pseudo_result <- function(out, if_fit) {
  out$if_beta_Treatment <- unname(if_fit$beta[2])
  out$if_beta_X <- unname(if_fit$beta[3])
  
  if_se_names <- c("Intercept", "Treatment", "X", paste0("Time", 2:8))
  for (idx in seq_along(if_se_names)) {
    out[[paste0("if_se_", if_se_names[idx])]] <- unname(if_fit$se[idx])
  }
  for (idx in seq_along(if_fit$landmark_times)) {
    out[[paste0("landmark_time_", idx)]] <- unname(if_fit$landmark_times[idx])
  }
  out
}

run_one_sim <- function(sim_id, config) {
  per_sim_file <- file.path(config$per_sim_dir, sprintf("sim_%04d.rds", sim_id))
  if (file.exists(per_sim_file)) {
    return(readRDS(per_sim_file))
  }
  
  dat <- tryCatch(
    make_sim_data(config, sim_id),
    error = function(e) e
  )
  
  if (inherits(dat, "error")) {
    res <- empty_result_row(
      sim_id,
      "failed",
      paste0("[data_generation] ", conditionMessage(dat))
    )
    saveRDS(res, per_sim_file)
    return(res)
  }
  
  res <- empty_result_row(sim_id, "ok", NA_character_)
  res$event_rate <- dat$event_rate
  res$treatment_rate <- dat$treatment_rate
  
  log_file <- file.path(config$per_sim_dir, sprintf("sim_%04d.log", sim_id))
  
  cox_fit <- tryCatch(
    run_quietly(run_cox(dat$sim_data), log_file, quiet = config$quiet),
    error = function(e) e
  )
  if (inherits(cox_fit, "error")) {
    res$cox_error <- conditionMessage(cox_fit)
  } else {
    res <- fill_cox_result(res, cox_fit)
  }
  
  odac_fit <- tryCatch(
    run_quietly(run_odac(dat$sim_data_chen, sim_id, config), log_file, quiet = config$quiet),
    error = function(e) e
  )
  if (inherits(odac_fit, "error")) {
    res$odac_error <- conditionMessage(odac_fit)
  } else {
    res <- fill_odac_result(res, odac_fit)
  }
  
  if_fit <- tryCatch(
    run_quietly(run_pseudo_reg(dat$sim_data, config), log_file, quiet = config$quiet),
    error = function(e) e
  )
  if (inherits(if_fit, "error")) {
    res$pseudo_reg_error <- conditionMessage(if_fit)
  } else {
    res <- fill_pseudo_result(res, if_fit)
  }
  
  method_success <- c(
    is.finite(res$cox_beta_Treatment) && is.finite(res$cox_se_Treatment),
    is.finite(res$odac_beta_Treatment) && is.finite(res$odac_se_Treatment),
    is.finite(res$if_beta_Treatment) && is.finite(res$if_se_Treatment)
  )
  
  if (all(method_success)) {
    res$status <- "ok"
  } else if (any(method_success)) {
    res$status <- "partial"
    res$error_message <- paste(
      na.omit(c(
        if (!method_success[1]) paste0("[cox] ", res$cox_error),
        if (!method_success[2]) paste0("[odac] ", res$odac_error),
        if (!method_success[3]) paste0("[pseudo_reg] ", res$pseudo_reg_error)
      )),
      collapse = " | "
    )
  } else {
    res$status <- "failed"
    res$error_message <- paste(
      na.omit(c(
        paste0("[cox] ", res$cox_error),
        paste0("[odac] ", res$odac_error),
        paste0("[pseudo_reg] ", res$pseudo_reg_error)
      )),
      collapse = " | "
    )
  }
  
  saveRDS(res, per_sim_file)
  res
}

write_result_tables <- function(summary_df, config) {
  summary_df <- summary_df[order(summary_df$sim_id), ]
  cox_ok <- is.finite(summary_df$cox_beta_Treatment) & is.finite(summary_df$cox_se_Treatment)
  odac_ok <- is.finite(summary_df$odac_beta_Treatment) & is.finite(summary_df$odac_se_Treatment)
  if_ok <- is.finite(summary_df$if_beta_Treatment) & is.finite(summary_df$if_se_Treatment)
  
  write.csv(summary_df, file.path(config$results_dir, "all_simulation_results.csv"), row.names = FALSE)
  write.csv(summary_df[summary_df$status != "ok", c("sim_id", "status", "error_message")],
            file.path(config$results_dir, "failed_simulations.csv"), row.names = FALSE)
  write.csv(summary_df[summary_df$status != "ok",
                       c("sim_id", "status", "cox_error", "odac_error", "pseudo_reg_error")],
            file.path(config$results_dir, "failed_methods.csv"), row.names = FALSE)
  
  write.csv(summary_df[cox_ok, c("cox_beta_Treatment", "cox_beta_X")],
            file.path(config$results_dir, "cox_beta_list.csv"), row.names = FALSE)
  write.csv(summary_df[odac_ok, c("odac_beta_Treatment", "odac_beta_X")],
            file.path(config$results_dir, "chen_beta_list.csv"), row.names = FALSE)
  write.csv(summary_df[if_ok, c("if_beta_Treatment", "if_beta_X")],
            file.path(config$results_dir, "if_beta_list.csv"), row.names = FALSE)
  
  write.csv(summary_df[cox_ok, c("cox_se_Treatment", "cox_se_X")],
            file.path(config$results_dir, "cox_se_list.csv"), row.names = FALSE)
  write.csv(summary_df[odac_ok, c("odac_se_Treatment", "odac_se_X")],
            file.path(config$results_dir, "chen_se_list.csv"), row.names = FALSE)
  write.csv(summary_df[if_ok, grep("^if_se_", names(summary_df), value = TRUE)],
            file.path(config$results_dir, "if_se_list.csv"), row.names = FALSE)
}

make_performance_summary <- function(summary_df, config) {
  ok_df <- summary_df
  z <- qnorm(0.975)
  true_beta <- config$true_treatment_beta
  
  one_method <- function(method, beta_col, se_col) {
    beta <- ok_df[[beta_col]]
    se <- ok_df[[se_col]]
    valid <- is.finite(beta) & is.finite(se)
    beta <- beta[valid]
    se <- se[valid]
    covered <- beta - z * se <= true_beta & beta + z * se >= true_beta
    coverage <- mean(covered)
    n <- length(beta)
    
    data.frame(
      scenario = config$scenario_name,
      method = method,
      nsim_attempted = config$nsim,
      nsim_success = n,
      true_beta = true_beta,
      mean_beta = mean(beta, na.rm = TRUE),
      bias = mean(beta - true_beta, na.rm = TRUE),
      rmse = sqrt(mean((beta - true_beta)^2, na.rm = TRUE)),
      empirical_sd = sd(beta, na.rm = TRUE),
      average_model_se = mean(se, na.rm = TRUE),
      se_to_sd_ratio = mean(se, na.rm = TRUE) / sd(beta, na.rm = TRUE),
      coverage = coverage,
      coverage_mc_se = sqrt(coverage * (1 - coverage) / n),
      stringsAsFactors = FALSE
    )
  }
  
  perf <- rbind(
    one_method("Cox", "cox_beta_Treatment", "cox_se_Treatment"),
    one_method("ODAC", "odac_beta_Treatment", "odac_se_Treatment"),
    one_method("Pseudo Reg", "if_beta_Treatment", "if_se_Treatment")
  )
  
  write.csv(perf, file.path(config$results_dir, "performance_summary.csv"), row.names = FALSE)
  perf
}

run_ph_scenario <- function(config) {
  load_required_packages()
  
  if (is.null(config$quiet)) {
    config$quiet <- TRUE
  }
  
  config$results_dir <- file.path(fedsurv_results_root(), config$scenario_name)
  config$pda_work_dir <- file.path(config$results_dir, "pda_work")
  config$per_sim_dir <- file.path(config$results_dir, "per_sim_rds")
  
  dir.create(config$results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$pda_work_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$per_sim_dir, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(config, file.path(config$results_dir, "scenario_config.rds"))
  
  message("Running ", config$scenario_name)
  message("nsim = ", config$nsim, ", ncore = ", config$ncore)
  message("Results directory: ", normalizePath(config$results_dir, mustWork = FALSE))
  
  sim_ids <- if (is.null(config$sim_ids)) seq_len(config$nsim) else config$sim_ids
  if (.Platform$OS.type == "unix" && config$ncore > 1) {
    rows <- parallel::mclapply(
      sim_ids,
      run_one_sim,
      config = config,
      mc.cores = config$ncore,
      mc.preschedule = FALSE,
      mc.set.seed = TRUE
    )
  } else {
    rows <- lapply(sim_ids, run_one_sim, config = config)
  }
  
  summary_df <- do.call(rbind, rows)
  write_result_tables(summary_df, config)
  perf <- make_performance_summary(summary_df, config)
  print(perf)
  
  invisible(list(results = summary_df, performance = perf, config = config))
}
