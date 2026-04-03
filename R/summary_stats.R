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
#      MUA-based initializer. For each covariate of interest j, computes the
#      MUA estimate and corresponding t-statistics, smooths the t-map, applies
#      locfdr, and constructs positive/negative initial regions.
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
  env$ZTY_raw <- crossprod(Z, Y)   # p2 x V
  env$YtY     <- colSums(Y^2)      # length V

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
  # Returns a list of length p1. Element j is the initial region_list
  # for covariate j (used to seed the coordinate-descent loop).
  #
  # NOTE: env$alphahat_full and env$XTXB are not built here; they are
  # constructed in sira() immediately after this function returns.

  p1 <- env$p1
  V  <- env$V

  if (env$verbose)
    message("  Initializing regions (MUA + locfdr, per covariate)...")

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
  # Mirror the original implementation:
  #   1. Compute the MUA solution and corresponding t-values.
  #   2. Smooth the t-values twice.
  #   3. Run locfdr with explicit settings.
  #   4. If nothing is selected, use the top and bottom 2%.
  #   5. Form one positive and one negative starting region when available.

  V   <- env$V
  df  <- env$n - env$p1 - env$p2
  sm  <- env$smoothing_matrix

  QtQ_inv <- .sira_qtq_inv(env)
  QTY     <- .sira_qty(env)
  mua_fit <- QtQ_inv %*% QTY
  beta_j  <- as.numeric(mua_fit[j, ])

  sse <- pmax(env$YtY - colSums(mua_fit * QTY), 1e-10)
  sigma2_hat <- pmax(sse / max(df, 1L), 1e-10)
  se_beta_j  <- sqrt(pmax(QtQ_inv[j, j], 1e-10) * sigma2_hat)
  t_values   <- beta_j / pmax(se_beta_j, 1e-10)

  smoothed_t <- as.numeric(sm %*% t_values)
  smoothed_t <- as.numeric(sm %*% smoothed_t)
  lfdr_keep  <- .sira_locfdr_select(smoothed_t)

  pos_vox <- which(lfdr_keep & smoothed_t > 0)
  neg_vox <- which(lfdr_keep & smoothed_t < 0)

  if (length(pos_vox) == 0L && length(neg_vox) == 0L) {
    n_tail  <- max(floor(V * 0.02), 1L)
    pos_vox <- utils::head(order(smoothed_t, decreasing = TRUE), n_tail)
    neg_vox <- utils::head(order(smoothed_t, decreasing = FALSE), n_tail)
  }

  region_list <- list()

  if (length(neg_vox) > 0L) {
    region_list[[length(region_list) + 1L]] <- list(neg_vox, mean(beta_j[neg_vox]))
  }
  if (length(pos_vox) > 0L) {
    region_list[[length(region_list) + 1L]] <- list(pos_vox, mean(beta_j[pos_vox]))
  }

  region_list
}

#' @keywords internal
.sira_qtq_inv <- function(env) {
  if (!is.null(env$QtQ_inv)) return(env$QtQ_inv)

  QtQ <- rbind(
    cbind(env$XTX_p1, env$XTZ),
    cbind(t(env$XTZ), crossprod(env$Z))
  )
  env$QtQ_inv <- solve(QtQ)
  env$QtQ_inv
}

#' @keywords internal
.sira_qty <- function(env) {
  rbind(env$XTY, env$ZTY_raw)
}

#' @keywords internal
.sira_locfdr_select <- function(smoothed_t) {
  fit <- tryCatch(
    suppressWarnings(
      suppressMessages(
        locfdr::locfdr(smoothed_t, bre = 120, df = 7, nulltype = 1, plot = 0)
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || is.null(fit$fdr)) {
    return(rep(FALSE, length(smoothed_t)))
  }

  fit$fdr < 0.05
}
