# =============================================================================
# R/summary_stats.R
# Two responsibilities:
#
#   1. .sira_setup_summary_stats()
#      Precomputes all fixed summary statistics that are functions of the
#      data matrices (X, Y, Z, Psi_star) only.  These never change during
#      the algorithm and are stored once in env.
#
#   2. .sira_initialize_all_regions()
#      MUA-based initializer.  For each covariate of interest j, runs a
#      marginal regression, applies spatial smoothing, fits a ridge regression
#      on candidate voxels, and extracts initial connected regions.  Also
#      builds the initial sparse XTX via .expand_active_set().
#
# Summary statistics stored in env after .sira_setup_summary_stats():
#
#   env$XTY            p1 x V   = t(X) %*% Y
#   env$XTZ            p1 x p2  = t(X) %*% Z
#   env$ZTZ1ZT         p2 x n   = solve(t(Z) %*% Z) %*% t(Z)
#   env$ZTZ1ZTY        p2 x V   = ZTZ1ZT %*% Y
#   env$ZTZ1ZTX        p2 x p1  = ZTZ1ZT %*% X
#   env$Y_Psi_starT    n  x L   = Y %*% t(Psi_star)
#   env$scaling_matrix_2  L x L = solve(Psi_star %*% t(Psi_star))
#
# After .sira_initialize_all_regions():
#
#   env$XTX            V x V sparse  (active block filled by .expand_active_set)
#   env$active_voxels  integer vector
#   (env$XTXB and env$alphahat_full are built in sira() immediately after)
# =============================================================================


# =============================================================================
# 1.  FIXED SUMMARY STATISTICS
# =============================================================================

#' @keywords internal
.sira_setup_summary_stats <- function(env) {

  X        <- env$X          # n x p1
  Y        <- env$Y          # n x V
  Z        <- env$Z          # n x p2
  Psi_star <- env$Psi_star   # L x V

  if (env$verbose)
    message("  Computing fixed summary statistics...")

  # ---- Cross-products involving X and Y ----------------------------------
  env$XTX_p1 <- crossprod(X)       # p1 x p1  (fixed throughout algorithm)
  env$XTY    <- crossprod(X, Y)    # p1 x V
  env$XTZ    <- crossprod(X, Z)    # p1 x p2

  # ---- Confounder projection matrix  (ZTZ)^{-1} Z^T  --------------------
  ZTZ      <- crossprod(Z)      # p2 x p2
  ZTZ1ZT   <- solve(ZTZ, t(Z)) # p2 x n
  env$ZTZ1ZT  <- ZTZ1ZT
  env$ZTZ1ZTY <- ZTZ1ZT %*% Y  # p2 x V
  env$ZTZ1ZTX <- ZTZ1ZT %*% X  # p2 x p1

  # ---- Spatial individual effects basis ----------------------------------
  # Y_Psi_starT    =  Y %*% Psi_star^T   [n x L]
  # scaling_matrix_2 = (Psi_star %*% Psi_star^T)^{-1}   [L x L]
  # Used in thetahat update:
  #   thetahat = (Y_Psi_starT
  #               - X %*% betahat %*% Psi_star^T
  #               - Z %*% gammahat %*% Psi_star^T) %*% scaling_matrix_2
  env$Y_Psi_starT      <- Y %*% t(Psi_star)   # n x L

  # scaling_matrix_2 = (Psi_star %*% t(Psi_star))^{-1}  [L x L]
  # When Psi_star was built internally via .sira_compute_psi_star:
  #   Psi_star = diag(sqrt_Lambda) %*% Psi_orth  (Psi_orth has orthonormal rows)
  #   => Psi_star %*% t(Psi_star) = diag(Lambda)
  #   => inverse = diag(1/Lambda)  — O(L) instead of O(L^3)
  # When the user supplies their own Psi_star, env$Lambda is not available
  # so we fall back to the full matrix solve.
  env$scaling_matrix_2 <- if (!is.null(env$Lambda)) {
    diag(1 / env$Lambda)
  } else {
    solve(tcrossprod(Psi_star))
  }

  if (env$verbose)
    message("  Fixed summary statistics done.")
}


# =============================================================================
# 2.  MUA-BASED INITIALIZER
# =============================================================================

#' @keywords internal
.sira_initialize_all_regions <- function(env) {
  # Returns a list of length p1.  Element j is the initial region_list
  # for covariate j (used to seed the coordinate-descent loop).
  # Also initialises env$XTX and env$active_voxels via .expand_active_set().
  #
  # NOTE: env$alphahat_full and env$XTXB are NOT built here; they are
  # constructed in sira() immediately after this function returns, once
  # alphahat_full is available to pass to .recalibrate_XTXB().

  p1 <- env$p1
  V  <- env$V

  if (env$verbose)
    message("  Initialising regions (MUA + ridge, per covariate)...")

  initial_rl_full <- vector("list", p1)
  all_init_voxels <- integer(0)

  for (j in seq_len(p1)) {
    rl_j <- .sira_initialize_one_covariate(j, env)
    initial_rl_full[[j]] <- rl_j
    all_init_voxels <- unique(c(all_init_voxels,
                                .voxels_in_regions(rl_j)))
    if (env$verbose) {
      nv <- length(.voxels_in_regions(rl_j))
      nr <- length(rl_j)
      message(sprintf("    Covariate %d: %d initial region(s), %d voxel(s)",
                      j, nr, nv))
    }
  }

  # Note: env$XTX_p1 (p1 x p1) is already set in .sira_setup_summary_stats.
  # env$XTXB (V x p1) is built in sira() once alphahat_full is available.
  initial_rl_full
}


# =============================================================================
# Per-covariate initializer
# =============================================================================

#' @keywords internal
.sira_initialize_one_covariate <- function(j, env) {
  # Steps (mirror the paper's section 2.5):
  #   1. Residualise Y on Z  ->  Y_tilde
  #   2. Compute marginal t-statistics for covariate j
  #   3. Select top-3% candidate voxels by |t|
  #   4. Ridge regression (X as design, Y_tilde[, cand_v] as response)
  #      -> extract row j of beta_cand as ridge_beta_j on cand_v
  #   5. Spatially smooth ridge_beta_j  (3 passes)
  #   6. Extract largest positive and negative connected components
  #   7. Return as a region_list (list of list(voxels, beta_value) pairs)

  X      <- env$X               # n x p1
  Y      <- env$Y               # n x V
  Z      <- env$Z               # n x p2
  V      <- env$V
  n      <- env$n
  p1     <- env$p1
  df     <- n - p1 - env$p2    # residual df for marginal t-stats

  # ---- 1 & 2. Marginal OLS statistics for covariate j ---------------
  # Two paths depending on whether Y is in memory (standard) or has been
  # pre-accumulated into summary statistics (batched via sira_batched()).

  if (!is.null(env$XTY_tilde_0)) {
    # --- Batched path: Y is not in memory ---
    # XTY_tilde_0 = t(X) %*% Y_tilde  (p1 x V, precomputed)
    # YtY_tilde   = colSums(Y_tilde^2) (length V, precomputed)
    # Use uncentred X for t-statistics (centering is for numerical stability
    # only; minor differences don't affect region initialisation quality).
    XtY_j <- as.numeric(env$XTY_tilde_0[j, ])   # length V
    XtX_j <- max(env$XTX_p1[j, j], 1e-10)       # scalar
    naive_beta_j <- XtY_j / XtX_j

    sigma2_hat <- pmax((env$YtY_tilde - XtY_j^2 / XtX_j) / df, 1e-10)
    se_beta    <- sqrt(sigma2_hat / XtX_j)
    t_values   <- naive_beta_j / se_beta

  } else {
    # --- Standard path: Y is in memory ---
    Y_tilde <- Y - Z %*% env$ZTZ1ZTY   # n x V

    xj   <- X[, j]
    xj_c <- xj - mean(xj)              # centre for numerical stability

    XtY_j <- as.numeric(crossprod(xj_c, Y_tilde))
    XtX_j <- max(sum(xj_c^2), 1e-10)
    naive_beta_j <- XtY_j / XtX_j

    YtY_j      <- colSums(Y_tilde^2)
    sigma2_hat <- pmax((YtY_j - XtY_j^2 / XtX_j) / df, 1e-10)
    se_beta    <- sqrt(sigma2_hat / XtX_j)
    t_values   <- naive_beta_j / se_beta
  }

  # ---- 3. Candidate voxels: top 3% by |t|  ---------------------------
  n_cand  <- max(floor(V * 0.03), 10L)
  pos_ord <- order(t_values, decreasing = TRUE)
  neg_ord <- order(t_values, decreasing = FALSE)
  cand_v  <- unique(c(head(pos_ord, n_cand), head(neg_ord, n_cand)))

  # ---- 4. Ridge regression on candidates  ----------------------------
  # Design matrix: X  (n x p1)
  # Response:      Y_tilde[, cand_v]  (n x |cand_v|)
  # In the batched path Y_tilde is not available, but:
  #   t(X) %*% Y_tilde[:, cand_v] = XTY_tilde_0[, cand_v]
  # which is already precomputed. In the standard path we compute directly.
  lam_ridge <- n * 0.01
  XtX_ridge <- env$XTX_p1 + lam_ridge * diag(p1)   # p1 x p1
  XtY_ridge <- if (!is.null(env$XTY_tilde_0)) {
    env$XTY_tilde_0[, cand_v, drop = FALSE]          # p1 x |cand_v|
  } else {
    crossprod(X, Y_tilde[, cand_v, drop = FALSE])    # p1 x |cand_v|
  }
  beta_cand  <- tryCatch(
    solve(XtX_ridge, XtY_ridge),   # p1 x |cand_v|
    error = function(e) {
      warning(sprintf(
        "Ridge solve failed for covariate %d; using marginal OLS fallback.", j))
      matrix(naive_beta_j[cand_v], nrow = p1, ncol = length(cand_v),
             byrow = TRUE)
    }
  )

  # Row j of beta_cand gives the ridge estimate for covariate j at each
  # candidate voxel. Assign into a full-length zero vector.
  ridge_beta_j          <- numeric(V)
  ridge_beta_j[cand_v]  <- beta_cand[j, ]

  # ---- 5. Spatial smoothing (3 passes)  ------------------------------
  sm       <- env$smoothing_matrix   # sparse V x V
  smoothed <- ridge_beta_j
  for (pass in seq_len(3L))
    smoothed <- as.numeric(sm %*% smoothed)

  # ---- 6. Extract connected regions  ---------------------------------
  cand_pos <- cand_v[smoothed[cand_v] > 0]
  cand_neg <- cand_v[smoothed[cand_v] < 0]

  region_list <- list()

  if (length(cand_neg) > 0L) {
    neg_region <- .get_largest_connected_component(cand_neg,
                                                   env$neighbor_list)
    if (length(neg_region) >= 5L) {
      beta_val <- median(ridge_beta_j[neg_region])
      region_list[[length(region_list) + 1L]] <- list(neg_region, beta_val)
    }
  }

  if (length(cand_pos) > 0L) {
    pos_region <- .get_largest_connected_component(cand_pos,
                                                   env$neighbor_list)
    if (length(pos_region) >= 5L) {
      beta_val <- median(ridge_beta_j[pos_region])
      region_list[[length(region_list) + 1L]] <- list(pos_region, beta_val)
    }
  }

  # ---- 7. Fallback if nothing found  ---------------------------------
  if (length(region_list) == 0L) {
    warning(sprintf(
      paste0("No strong initial regions found for covariate %d. ",
             "Using single strongest voxel and its neighbors as fallback."), j))
    seed_v   <- which.max(abs(t_values))
    seed_nb  <- env$neighbor_list[[seed_v]]
    init_v   <- c(seed_v, head(seed_nb, 4L))
    beta_val <- median(ridge_beta_j[init_v])
    region_list[[1L]] <- list(init_v, beta_val)
  }

  region_list
}
