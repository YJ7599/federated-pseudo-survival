## Scenario 5: PH setting, approximately 30% event rate,
## 18 sites with sizes 1000 x 2, 500 x 6 and 250 x 10.

repo_root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f), "..", "..")) else getwd()
})
setwd(repo_root)
source(file.path("simulation", "scenario1_ph", "run_ph_scenario_common.R"))

config <- list(
  scenario_name = "scenario5_er30_mixed18sites_1000x2_500x6_250x10",
  nsim = 1000L,
  ncore = as.integer(Sys.getenv("FEDSURV_NCORE", "50")),
  seed = 2026070205L,
  true_treatment_beta = log(1.15),
  true_x_beta = log(0.7),
  Ns = sort(c(1000L, 1000L, rep(500L, 6), rep(250L, 10)), decreasing = TRUE),
  base_hazard = 1 / 20,
  censoring_rate = -log(0.3) / 10,
  admin_censor = 40,
  round_observed_time = FALSE,
  round_odac_time = TRUE,
  bknots_type = "first_site_max",
  risk_set_minimum = NULL,
  landmark_style = "two_endpoint",
  use_event_time_grid_for_landmarks = TRUE,
  drop_zero_var_odac_sites = FALSE
)

## Optional: FEDSURV_SIM_IDS="1:20" runs a subset of the replicates.
if (nzchar(Sys.getenv("FEDSURV_SIM_IDS"))) config$sim_ids <- as.integer(eval(parse(text = Sys.getenv("FEDSURV_SIM_IDS"))))

run_ph_scenario(config)
