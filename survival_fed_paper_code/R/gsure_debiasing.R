## =============================================================================
## Variance-adaptive soft-thresholding of site-specific deviations with the
## threshold selected by GSURE (Section 2.4).
##
##   build_SigmaDelta_singlepass_coord_safe(): covariance matrix of the deviations
##       Delta_k = beta_local_k - beta_global from single-pass sandwich summaries
##   sure_euclidean_tau_scaled(): GSURE criterion for threshold tau * sqrt(V(Delta_k))
##   choose_tau_euclidean_gsure_scaled(): grid search for the GSURE-minimizing tau
## =============================================================================

safe_inv <- function(A, ridge_floor = 1e-10) {
  tryCatch(
    solve(A),
    error = function(e) {
      ridge <- ridge_floor * mean(diag(A), na.rm = TRUE)
      if (!is.finite(ridge) || ridge <= 0) ridge <- ridge_floor
      solve(A + diag(ridge, nrow(A)))
    }
  )
}

make_psd <- function(S, eig_floor = 1e-10) {
  e <- eigen(S, symmetric = TRUE)
  e$values[e$values < eig_floor] <- eig_floor
  e$vectors %*% diag(e$values) %*% t(e$vectors)
}

build_SigmaDelta_singlepass_coord_safe <- function(A_glob, meat_by_site, se_loc_vec, idx_glm, eig_floor = 1e-10, ridge_floor = 1e-10) {
  K <- length(se_loc_vec)
  Ainv <- safe_inv(A_glob, ridge_floor = ridge_floor)
  Vglob_mat <- Ainv %*% Reduce(`+`, meat_by_site) %*% t(Ainv)
  Vglob_jj  <- as.numeric(Vglob_mat[idx_glm, idx_glm])
  
  Ck <- vapply(meat_by_site, function(Mk) { as.numeric((Ainv %*% Mk %*% t(Ainv))[idx_glm, idx_glm]) }, numeric(1))
  Vloc_jj <- se_loc_vec^2
  
  Sigma <- diag(Vloc_jj) + Vglob_jj * matrix(1, K, K) - (outer(Ck, rep(1, K)) + outer(rep(1, K), Ck))
  list(Sigma = make_psd(Sigma, eig_floor = eig_floor))
}

## GSURE criterion with site-specific threshold tau * se(Delta_k)
sure_euclidean_tau_scaled <- function(Delta, Sigma, tau) {
  se_delta <- sqrt(diag(Sigma))
  thresh <- tau * se_delta
  h <- sign(Delta) * pmax(abs(Delta) - thresh, 0)
  div_term <- sum(diag(Sigma) * (abs(Delta) > thresh))
  as.numeric(sum(h^2) - 2 * sum(h * Delta) + 2 * div_term)
}

## Grid search over tau (tau_max scaled to the largest |Delta_k| / se(Delta_k))
choose_tau_euclidean_gsure_scaled <- function(Delta, Sigma, n_grid = 80, tau_min = 0, tau_max = NULL) {
  se_delta <- sqrt(diag(Sigma))
  if (is.null(tau_max)) {
    valid_idx <- which(se_delta > 1e-12)
    if (length(valid_idx) > 0) {
      tau_max <- max(abs(Delta[valid_idx]) / se_delta[valid_idx]) + 1e-4
    } else {
      tau_max <- 1
    }
  }
  tau_max <- max(tau_max, tau_min + 1e-12)
  grid <- sort(unique(c(0, seq(tau_min, tau_max, length.out = n_grid))))
  vals <- vapply(grid, function(t) sure_euclidean_tau_scaled(Delta, Sigma, t), numeric(1))
  list(tau = grid[which.min(vals)], grid = grid, sure = vals)
}

