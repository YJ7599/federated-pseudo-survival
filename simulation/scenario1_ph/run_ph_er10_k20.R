## Scenario 1 (proportional hazards): approximately 10% event rate,
## 20 sites with sizes 500 x 10 and 250 x 10.

repo_root <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) normalizePath(file.path(dirname(f), "..", "..")) else getwd()
})
setwd(repo_root)
source(file.path("simulation", "scenario1_ph", "ph_simulation_engine.R"))

config <- list(
  scenario_name = "ph_er10_k20_500x10_250x10",
  nsim = 1000L,
  ncore = as.integer(Sys.getenv("FEDSURV_NCORE", "50")),
  seed = 2026070207L,
  true_treatment_beta = log(1.15),
  true_x_beta = log(0.7),
  Ns = sort(c(rep(500L, 10), rep(250L, 10)), decreasing = TRUE),
  base_hazard = 1 / 50,
  censoring_rate = -log(0.2) / 10,
  admin_censor = 20,
  round_observed_time = TRUE,
  round_odac_time = TRUE,
  bknots_type = "first_site_max",
  risk_set_minimum = NULL,
  landmark_style = "two_endpoint",
  use_event_time_grid_for_landmarks = TRUE,
  drop_zero_var_odac_sites = TRUE
)

## Optional: FEDSURV_SIM_IDS="1:20" runs a subset of the replicates.
if (nzchar(Sys.getenv("FEDSURV_SIM_IDS"))) config$sim_ids <- as.integer(eval(parse(text = Sys.getenv("FEDSURV_SIM_IDS"))))

run_ph_scenario(config)
