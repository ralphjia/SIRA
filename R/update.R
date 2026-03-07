# =============================================================================
# R/update.R
# Outer coordinate-descent loop for SIRA.
#
# Public entry point:
#   .sira_full_algorithm(env, max_iter)
#     Runs up to max_iter outer iterations.  Each iteration:
#       1. Sweeps through all p1 covariates, applying the five region ops.
#       2. Updates gammahat, thetahat, XTY_tilde.
#       3. Every XTXB_RECALIB_FREQ iterations, recalibrates XTXB from scratch
#          to prevent floating-point drift from incremental updates.
#       4. Checks convergence: stops when the Frobenius norm of the change in
#          betahat (the p1 rows of alphahat_full) falls below env$delta.
#
#   Returns list(region_list_full, convergence).
#   region_list_full : list of length p1, final region configuration.
#   convergence      : data.frame with columns
#                        iteration, num_regions, num_voxels, alpha_norm_diff
# =============================================================================

# Internal constant: full XTXB recomputation every this many outer iterations.
# Incremental updates in region_ops accumulate O(iter * n_ops) floating-point
# error; periodic recalibration keeps max absolute error below ~1e-8.
# Set to 1 to recompute every iteration (safe but ~20% slower).
.XTXB_RECALIB_FREQ <- 5L


# =============================================================================
# MAIN LOOP
# =============================================================================

#' @keywords internal
.sira_full_algorithm <- function(env, max_iter) {

  p1           <- env$p1
  V            <- env$V
  conv_rows    <- vector("list", max_iter)

  # Working copy of the region list, one element per covariate of interest.
  region_list_full <- env$initial_region_list_full

  for (iter in seq_len(max_iter)) {

    # Snapshot betahat before the sweep for convergence checking
    alpha_old <- env$alphahat_full[env$p1_index, , drop = FALSE]  # p1 x V

    # ------------------------------------------------------------------
    # 1.  Covariate sweep — apply all five region ops to each beta_j
    # ------------------------------------------------------------------
    for (j in seq_len(p1)) {
      region_list_full[[j]] <- .sira_region_ops_one_covariate(
        j, region_list_full[[j]], env
      )
      # Sync alphahat_full row j from the (possibly reordered / pruned)
      # region list returned by the ops.
      env$alphahat_full[j, ] <- .region_list_to_betahat(
        region_list_full[[j]], V
      )
    }

    # ------------------------------------------------------------------
    # 2.  Update nuisance parameters (gamma, theta) and working response
    # ------------------------------------------------------------------
    env$gammahat <- .sira_update_gammahat(env)             # p2 x V
    env$alphahat_full[env$p2_index, ] <- env$gammahat

    env$thetahat  <- .sira_update_thetahat(env)            # n x L
    env$XTY_tilde <- .sira_update_XTY_tilde(env)          # p1 x V

    # ------------------------------------------------------------------
    # 3.  Periodic XTXB recalibration
    # ------------------------------------------------------------------
    if (iter %% .XTXB_RECALIB_FREQ == 0L) {
      env$XTXB <- .recalibrate_XTXB(env)
    }

    # ------------------------------------------------------------------
    # 4.  Convergence check
    # ------------------------------------------------------------------
    alpha_new   <- env$alphahat_full[env$p1_index, , drop = FALSE]
    norm_diff   <- norm(alpha_new - alpha_old, type = "F")

    num_regions <- sum(vapply(region_list_full, length, integer(1L)))
    num_voxels  <- sum(vapply(region_list_full,
                              function(rl) length(.voxels_in_regions(rl)),
                              integer(1L)))

    conv_rows[[iter]] <- data.frame(
      iteration       = iter,
      num_regions     = num_regions,
      num_voxels      = num_voxels,
      alpha_norm_diff = norm_diff
    )

    if (env$verbose)
      message(sprintf(
        "  Iter %3d | regions: %d | voxels: %d | ||Δbeta||_F: %.6f",
        iter, num_regions, num_voxels, norm_diff
      ))

    if (norm_diff < env$delta) {
      if (env$verbose)
        message(sprintf(
          "  Converged at iteration %d (||Δbeta||_F = %.6f < delta = %.6f).",
          iter, norm_diff, env$delta
        ))
      break
    }
  }

  # One final recalibration to ensure XTXB is exact at exit, regardless of
  # where we stopped relative to the recalibration schedule.
  env$XTXB <- .recalibrate_XTXB(env)

  convergence <- do.call(rbind, conv_rows[!vapply(conv_rows, is.null, logical(1L))])
  rownames(convergence) <- NULL

  list(
    region_list_full = region_list_full,
    convergence      = convergence
  )
}
