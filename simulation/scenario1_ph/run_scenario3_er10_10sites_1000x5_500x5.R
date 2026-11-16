## Scenario 3: PH setting, approximately 10% event rate,
## 10 sites with sizes 1000 x 5 and 500 x 5.

repo_root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f), "..", "..")) else getwd()
})
setwd(repo_root)
source(file.path("simulation", "scenario1_ph", "run_ph_scenario_common.R"))

config <- list(
  scenario_name = "scenario3_er10_10sites_1000x5_500x5",
  nsim = 1000L,
  ncore = as.integer(Sys.getenv("FEDSURV_NCORE", "50")),
  seed = 2026070203L,
  true_treatment_beta = log(1.15),
  true_x_beta = log(0.7),
  Ns = c(rep(1000L, 5), rep(500L, 5)),
  base_hazard = 1 / 50,
  censoring_rate = -log(0.2) / 10,
  admin_censor = 20,
  round_observed_time = TRUE,
  round_odac_time = FALSE,
  bknots_type = "fixed_20",
  risk_set_minimum = 0.05,
  landmark_style = "scenario3",
  use_event_time_grid_for_landmarks = FALSE,
  drop_zero_var_odac_sites = FALSE
)

## Optional: FEDSURV_SIM_IDS="1:20" runs a subset of the replicates.
if (nzchar(Sys.getenv("FEDSURV_SIM_IDS"))) config$sim_ids <- as.integer(eval(parse(text = Sys.getenv("FEDSURV_SIM_IDS"))))

run_ph_scenario(config)
