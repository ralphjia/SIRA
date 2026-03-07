# =============================================================================
# R/region_ops.R
# The five region operations of the SIRA coordinate-descent algorithm:
#
#   1. revalue  — optimise the scalar beta for each existing region
#   2. expand   — grow regions by greedily absorbing zero-valued neighbours
#   3. shrink   — remove boundary voxels whose deletion improves the loss
#   4. merge    — fuse two adjacent regions into one
#   5. split    — divide a large region into two connected sub-regions
#
# Public entry point:
#   .sira_region_ops_one_covariate(j, region_list, env)
#     Runs all five operations in sequence for covariate j and returns
#     the updated region_list.
#
# Each operation reads/writes:
#   env$alphahat_full[j, ]   — the coefficient image for covariate j
#   env$XTXB[, j]            — incremental XTX*betahat column (utils.R)
#   env$XTX, env$active_voxels — expanded as new voxels are considered
#
# Loss differences are evaluated via .loss_difference() (utils.R).
# The optimal delta for revalue is found by the C++ bracket_quadratic().
# =============================================================================


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

#' @keywords internal
.sira_region_ops_one_covariate <- function(j, region_list, env) {
  # Refresh local betahat before each operation so that incremental XTXB
  # updates made inside ops are reflected in subsequent ops.
  betahat <- as.numeric(env$alphahat_full[j, ])

  region_list <- .sira_op_revalue(j, region_list, betahat, env)
  betahat     <- as.numeric(env$alphahat_full[j, ])

  region_list <- .sira_op_expand(j, region_list, betahat, env)
  betahat     <- as.numeric(env$alphahat_full[j, ])

  region_list <- .sira_op_shrink(j, region_list, betahat, env)
  betahat     <- as.numeric(env$alphahat_full[j, ])

  region_list <- .sira_op_merge(j, region_list, betahat, env)
  betahat     <- as.numeric(env$alphahat_full[j, ])

  region_list <- .sira_op_split(j, region_list, betahat, env)

  # House-keeping: drop zero regions, collapse duplicates
  region_list <- .clean_zero_regions(region_list)
  region_list <- .merge_duplicate_regions(region_list)

  region_list
}


# =============================================================================
# 1.  REVALUE
# Finds the optimal scalar shift delta for every region via bracket_quadratic,
# then accepts the move if the full loss_difference is negative.
# =============================================================================

#' @keywords internal
.sira_op_revalue <- function(j, region_list, betahat, env) {
  for (k in seq_along(region_list)) {
    vox    <- region_list[[k]][[1]]
    beta_k <- region_list[[k]][[2]]

    nb <- .get_region_neighbors(vox, env)

    bq    <- .bracket_quadratic_setup(
      region_voxels    = vox,
      region_size      = length(vox),
      region_beta      = beta_k,
      region_neighbors = nb,
      betahat          = betahat,
      which_beta       = j,
      XTXB_local       = env$XTXB,
      env              = env
    )
    # C++ bracket_quadratic expects b as a flat vector: lower bounds first,
    # then upper bounds (length 2*m). .bracket_quadratic_setup returns an
    # n_rows x 2 matrix, so we unpack the two columns.
    b_flat <- c(bq[[2]][, 1L], bq[[2]][, 2L])
    delta  <- bracket_quadratic(bq[[1]], b = b_flat)$par

    if (abs(delta) < 1e-10) next

    # Verify with the full loss difference (captures cross-covariate terms
    # that bracket_quadratic does not model).
    d <- .loss_difference(betahat, j, delta, vox, env = env)
    if (d < -1e-10) {
      .update_XTXB_incremental(j, delta, vox, env = env)
      new_beta            <- beta_k + delta
      env$alphahat_full[j, vox] <- new_beta
      betahat[vox]        <- new_beta
      region_list[[k]][[2]] <- new_beta
    }
  }
  region_list
}


# =============================================================================
# 2.  EXPAND
# For each region, tries to absorb each unoccupied face-adjacent neighbour.
# Greedy: accepts the first beneficial voxel, then continues.
# =============================================================================

#' @keywords internal
.sira_op_expand <- function(j, region_list, betahat, env) {
  occupied <- .voxels_in_regions(region_list)

  for (k in seq_along(region_list)) {
    vox    <- region_list[[k]][[1]]
    beta_k <- region_list[[k]][[2]]
    if (abs(beta_k) < 1e-10) next   # zero region — nothing to expand into

    # Unoccupied face-adjacent neighbours
    nb_all     <- unique(unlist(env$neighbor_list[vox], use.names = FALSE))
    candidates <- setdiff(nb_all, occupied)
    if (length(candidates) == 0L) next

    # Greedy single-voxel additions
    for (v in candidates) {
      if (v %in% occupied) next   # earlier iteration in this loop may have added it
      d <- .loss_difference(betahat, j, beta_k, v, env = env)
      if (d < -1e-10) {
        .update_XTXB_incremental(j, beta_k, v, env = env)
        env$alphahat_full[j, v] <- beta_k
        betahat[v]              <- beta_k
        occupied                <- c(occupied, v)
        region_list[[k]][[1]]   <- c(region_list[[k]][[1]], v)
      }
    }
  }
  region_list
}


# =============================================================================
# 3.  SHRINK
# For each region, tries to remove boundary voxels one at a time.
# A removal is only attempted if it (a) improves the loss and (b) does not
# disconnect the region.
# =============================================================================

#' @keywords internal
.sira_op_shrink <- function(j, region_list, betahat, env) {
  for (k in seq_along(region_list)) {
    vox    <- region_list[[k]][[1]]
    beta_k <- region_list[[k]][[2]]
    if (length(vox) <= 1L) next

    # Boundary voxels: inside the region and with at least one external neighbour
    is_boundary <- vapply(vox, function(v) {
      any(!(env$neighbor_list[[v]] %in% vox))
    }, logical(1L))
    boundary_vox <- vox[is_boundary]
    if (length(boundary_vox) == 0L) next

    current_vox <- vox

    for (v in boundary_vox) {
      if (!(v %in% current_vox)) next   # already removed in this loop
      if (length(current_vox) <= 1L) break

      # Fast path: if v has >1 neighbour in the region it might be an
      # articulation point — do a full connectivity check.  If it has ≤1
      # in-region neighbour, removal cannot disconnect the component.
      n_in <- sum(env$neighbor_list[[v]] %in% current_vox)
      if (n_in > 1L) {
        remaining <- setdiff(current_vox, v)
        if (length(.get_connected_components(remaining, env$neighbor_list)) > 1L)
          next   # would disconnect — skip
      }

      # Check loss diff: zero out voxel v (delta = -beta_k)
      d <- .loss_difference(betahat, j, -beta_k, v, env = env)
      if (d < -1e-10) {
        .update_XTXB_incremental(j, -beta_k, v, env = env)
        env$alphahat_full[j, v] <- 0
        betahat[v]              <- 0
        current_vox             <- setdiff(current_vox, v)
      }
    }

    region_list[[k]][[1]] <- current_vox
    if (length(current_vox) == 0L) region_list[[k]][[2]] <- 0
  }
  region_list
}


# =============================================================================
# 4.  MERGE
# Scans all pairs of adjacent regions and applies the single best merge.
# Repeats until no beneficial merge remains (greedy best-first).
# =============================================================================

#' @keywords internal
.sira_op_merge <- function(j, region_list, betahat, env) {
  if (length(region_list) < 2L) return(region_list)

  changed <- TRUE
  while (changed && length(region_list) >= 2L) {
    changed  <- FALSE
    n_reg    <- length(region_list)
    best_d   <- -1e-10   # threshold: only accept strictly beneficial merges
    best_k1  <- NA_integer_
    best_k2  <- NA_integer_
    best_new <- NA_real_

    for (k1 in seq_len(n_reg - 1L)) {
      vox1  <- region_list[[k1]][[1]]
      beta1 <- region_list[[k1]][[2]]
      nb1   <- unique(unlist(env$neighbor_list[vox1], use.names = FALSE))

      for (k2 in seq(k1 + 1L, n_reg)) {
        vox2  <- region_list[[k2]][[1]]
        beta2 <- region_list[[k2]][[2]]

        if (!any(nb1 %in% vox2)) next   # not face-adjacent

        # Candidate merged beta: unconstrained OLS over the combined region
        beta_new <- .optimal_merge_beta(vox1, vox2, beta1, beta2, j, env)
        d <- .loss_difference(betahat, j,
                              beta_new - beta1, vox1,
                              beta_new - beta2, vox2, env)
        if (d < best_d) {
          best_d   <- d
          best_k1  <- k1
          best_k2  <- k2
          best_new <- beta_new
        }
      }
    }

    if (!is.na(best_k1)) {
      vox1  <- region_list[[best_k1]][[1]]
      vox2  <- region_list[[best_k2]][[1]]
      beta1 <- region_list[[best_k1]][[2]]
      beta2 <- region_list[[best_k2]][[2]]

      .update_XTXB_incremental(j, best_new - beta1, vox1,
                               best_new - beta2, vox2, env = env)
      env$alphahat_full[j, c(vox1, vox2)] <- best_new
      betahat[c(vox1, vox2)]              <- best_new

      region_list[[best_k1]][[1]] <- c(vox1, vox2)
      region_list[[best_k1]][[2]] <- best_new
      region_list[[best_k2]]      <- NULL
      region_list <- region_list[!vapply(region_list, is.null, logical(1L))]

      changed <- TRUE
    }
  }
  region_list
}


# =============================================================================
# 5.  SPLIT
# Tries to divide each large region along each of the 3 coordinate axes,
# taking the largest connected component on each side.  Accepts the best
# axis split if it improves the loss.
#
# Leftover voxels (trimmed when taking LCCs) are zeroed immediately.
# Their contribution to the loss check is computed approximately (cross-terms
# with sub1/sub2 are ignored), which is acceptable since leftover sets are
# small and any remaining benefit is recovered by subsequent shrink ops.
# =============================================================================

#' @keywords internal
.sira_op_split <- function(j, region_list, betahat, env) {
  result <- vector("list", 0L)

  for (k in seq_along(region_list)) {
    vox    <- region_list[[k]][[1]]
    beta_k <- region_list[[k]][[2]]

    # Minimum region size to produce two viable sub-regions of ≥5 voxels
    if (length(vox) < 10L) {
      result[[length(result) + 1L]] <- region_list[[k]]
      next
    }

    coords     <- .get_voxel_coords(vox, env)
    best_d     <- -1e-10
    best_split <- NULL

    for (axis in seq_len(3L)) {
      ax_vals <- coords[, axis]
      if (length(unique(ax_vals)) < 2L) next

      med_ax <- median(ax_vals)
      sub1   <- .get_largest_connected_component(
        vox[ax_vals <= med_ax], env$neighbor_list)
      sub2   <- .get_largest_connected_component(
        vox[ax_vals >  med_ax], env$neighbor_list)
      if (length(sub1) < 5L || length(sub2) < 5L) next

      b1 <- .optimal_split_beta(sub1, beta_k, j, env)
      b2 <- .optimal_split_beta(sub2, beta_k, j, env)

      leftover <- setdiff(vox, c(sub1, sub2))

      # Loss diff for the primary split (sub1 → b1, sub2 → b2)
      d <- .loss_difference(betahat, j, b1 - beta_k, sub1,
                            b2 - beta_k, sub2, env)

      # Approximate loss diff for zeroing leftover (ignores cross-terms with
      # sub1/sub2 since those betas have already changed in the proposal)
      if (length(leftover) > 0L)
        d <- d + .loss_difference(betahat, j, -beta_k, leftover, env = env)

      if (d < best_d) {
        best_d     <- d
        best_split <- list(sub1 = sub1, sub2 = sub2,
                           b1 = b1, b2 = b2, leftover = leftover)
      }
    }

    if (!is.null(best_split)) {
      sub1 <- best_split$sub1;  sub2 <- best_split$sub2
      b1   <- best_split$b1;    b2   <- best_split$b2
      lo   <- best_split$leftover

      # Apply sub1 and sub2 updates together (captures cross-terms)
      .update_XTXB_incremental(j, b1 - beta_k, sub1,
                               b2 - beta_k, sub2, env = env)
      env$alphahat_full[j, sub1] <- b1
      env$alphahat_full[j, sub2] <- b2
      betahat[sub1] <- b1
      betahat[sub2] <- b2

      # Zero out any leftover voxels
      if (length(lo) > 0L) {
        .update_XTXB_incremental(j, -beta_k, lo, env = env)
        env$alphahat_full[j, lo] <- 0
        betahat[lo]              <- 0
      }

      result[[length(result) + 1L]] <- list(sub1, b1)
      result[[length(result) + 1L]] <- list(sub2, b2)
    } else {
      result[[length(result) + 1L]] <- region_list[[k]]
    }
  }

  result
}


# =============================================================================
# INTERNAL HELPERS
# =============================================================================


# -----------------------------------------------------------------------------
# .optimal_merge_beta
# Unconstrained OLS estimate for the common beta after merging vox1 and vox2.
# Accounts for cross-terms with the rest of the image via env$XTXB.
#
# Derivation (MS loss, setting dL/d(beta_new) = 0):
#   beta_new = (xty_merged - cross_outside) / xtx_merged
# where
#   cross_outside = XTXB[vox_merged] - beta1 * XTX[vox1, vox_merged]
#                                     - beta2 * XTX[vox2, vox_merged]
# -----------------------------------------------------------------------------

#' @keywords internal
.optimal_merge_beta <- function(vox1, vox2, beta1, beta2, j, env) {
  # Closed-form OLS optimal beta for the merged region.
  # Setting dL_MS/d(beta_new) = 0 (voxels are independent):
  #   beta_new = (xty_m - xtxb_m + XTX_jj*(|v1|*beta1 + |v2|*beta2))
  #              / (|merged| * XTX_jj)
  vox_m   <- c(vox1, vox2)
  XTX_jj  <- env$XTX_p1[j, j]
  xtx_m   <- length(vox_m) * XTX_jj
  if (xtx_m < 1e-10) return((beta1 + beta2) / 2)

  xtxb_m  <- sum(env$XTXB[vox_m, j])
  xty_m   <- sum(env$XTY_tilde[j, vox_m])
  self    <- XTX_jj * (length(vox1) * beta1 + length(vox2) * beta2)

  (xty_m - xtxb_m + self) / xtx_m
}


# -----------------------------------------------------------------------------
# .optimal_split_beta
# Unconstrained OLS estimate for a sub-region's beta after splitting, holding
# all other regions (including the complementary sub-region) at their current
# betahat values.
#
# Derivation (MS loss, setting dL/d(beta_sub) = 0):
#   beta_sub = current_beta + (xty_sub - xtxb_sub) / xtx_sub
# where xtxb_sub already encodes the self-contribution at current_beta, so
# subtracting it isolates the cross-terms with external voxels.
# -----------------------------------------------------------------------------

#' @keywords internal
.optimal_split_beta <- function(sub_vox, current_beta, j, env) {
  # Closed-form OLS optimal beta for a sub-region after splitting.
  # Setting dL_MS/d(beta_sub) = 0:
  #   beta_sub = current_beta + (xty_sub - xtxb_sub) / (|sub| * XTX_jj)
  xtx_s <- length(sub_vox) * env$XTX_p1[j, j]
  if (xtx_s < 1e-10) return(current_beta)

  xtxb_s <- sum(env$XTXB[sub_vox, j])
  xty_s  <- sum(env$XTY_tilde[j, sub_vox])

  current_beta + (xty_s - xtxb_s) / xtx_s
}
