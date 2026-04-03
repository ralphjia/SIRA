# =============================================================================
# R/region_ops.R
# Region updates for SIRA.
#
# For a single covariate j, Algorithm 2 in the paper repeatedly:
#   1. Enumerates all candidate region operations.
#   2. Computes the loss decrease for every candidate.
#   3. Accepts the single best operation.
#   4. Repeats until no candidate decreases the loss.
#
# The operations are:
#   - revalue
#   - expand
#   - shrink
#   - merge
#   - split
# =============================================================================


# =============================================================================
# MAIN ENTRY POINT
# =============================================================================

#' @keywords internal
.sira_region_ops_one_covariate <- function(j, region_list, env) {
  updated <- TRUE

  while (updated) {
    updated <- FALSE
    betahat <- as.numeric(env$alphahat_full[j, ])

    best_op <- .sira_find_best_operation(j, region_list, betahat, env)
    if (!is.null(best_op) && is.finite(best_op$loss_diff) &&
        best_op$loss_diff < -1e-10) {
      region_list <- .sira_apply_operation(j, region_list, best_op, env)
      region_list <- .clean_zero_regions(region_list)
      region_list <- .sort_regions_by_beta(region_list)
      updated <- TRUE
    }
  }

  region_list
}


# =============================================================================
# CANDIDATE SEARCH
# =============================================================================

#' @keywords internal
.sira_find_best_operation <- function(j, region_list, betahat, env) {
  best_op <- NULL

  best_op <- .better_operation(best_op,
                               .sira_best_revalue_candidate(j, region_list, betahat, env))
  best_op <- .better_operation(best_op,
                               .sira_best_expand_candidate(j, region_list, betahat, env))
  best_op <- .better_operation(best_op,
                               .sira_best_shrink_candidate(j, region_list, betahat, env))
  best_op <- .better_operation(best_op,
                               .sira_best_merge_candidate(j, region_list, betahat, env))
  best_op <- .better_operation(best_op,
                               .sira_best_split_candidate(j, region_list, betahat, env))

  best_op
}

#' @keywords internal
.better_operation <- function(current_best, candidate) {
  if (is.null(candidate)) return(current_best)
  if (is.null(current_best)) return(candidate)
  if (candidate$loss_diff < current_best$loss_diff) candidate else current_best
}


# =============================================================================
# REVALUE
# =============================================================================

#' @keywords internal
.sira_best_revalue_candidate <- function(j, region_list, betahat, env) {
  best <- NULL

  for (k in seq_along(region_list)) {
    vox    <- region_list[[k]][[1]]
    beta_k <- region_list[[k]][[2]]
    if (length(vox) == 0L) next

    nb <- .get_region_neighbors(vox, env)
    bq <- .bracket_quadratic_setup(
      region_voxels    = vox,
      region_size      = length(vox),
      region_beta      = beta_k,
      region_neighbors = nb,
      betahat          = betahat,
      which_beta       = j,
      XTXB_local       = env$XTXB,
      env              = env
    )

    b_flat <- c(bq[[2]][, 1L], bq[[2]][, 2L])
    delta  <- bracket_quadratic(bq[[1]], b = b_flat)$par
    if (abs(delta) < 1e-10) next

    d <- .loss_difference(betahat, j, delta, vox, env = env)
    cand <- list(type = "revalue", k = k, vox = vox, delta = delta,
                 beta_new = beta_k + delta, loss_diff = d)
    best <- .better_operation(best, cand)
  }

  best
}


# =============================================================================
# EXPAND
# =============================================================================

#' @keywords internal
.sira_best_expand_candidate <- function(j, region_list, betahat, env) {
  best     <- NULL
  occupied <- .voxels_in_regions(region_list)

  for (k in seq_along(region_list)) {
    vox    <- region_list[[k]][[1]]
    beta_k <- region_list[[k]][[2]]
    if (length(vox) == 0L || abs(beta_k) < 1e-10) next

    nb_all     <- unique(unlist(env$neighbor_list[vox], use.names = FALSE))
    candidates <- setdiff(nb_all, occupied)
    if (length(candidates) == 0L) next

    for (v in candidates) {
      d <- .loss_difference(betahat, j, beta_k, v, env = env)
      cand <- list(type = "expand", k = k, vox = v, beta_new = beta_k,
                   loss_diff = d)
      best <- .better_operation(best, cand)
    }
  }

  best
}


# =============================================================================
# SHRINK
# =============================================================================

#' @keywords internal
.sira_best_shrink_candidate <- function(j, region_list, betahat, env) {
  best <- NULL

  for (k in seq_along(region_list)) {
    vox    <- region_list[[k]][[1]]
    beta_k <- region_list[[k]][[2]]
    if (length(vox) <= 1L) next

    is_boundary <- vapply(vox, function(v) {
      any(!(env$neighbor_list[[v]] %in% vox))
    }, logical(1L))
    boundary_vox <- vox[is_boundary]
    if (length(boundary_vox) == 0L) next

    for (v in boundary_vox) {
      n_in <- sum(env$neighbor_list[[v]] %in% vox)
      if (n_in > 1L) {
        remaining <- setdiff(vox, v)
        if (length(.get_connected_components(remaining, env$neighbor_list)) > 1L)
          next
      }

      d <- .loss_difference(betahat, j, -beta_k, v, env = env)
      cand <- list(type = "shrink", k = k, vox = v, delta = -beta_k,
                   loss_diff = d)
      best <- .better_operation(best, cand)
    }
  }

  best
}


# =============================================================================
# MERGE
# =============================================================================

#' @keywords internal
.sira_best_merge_candidate <- function(j, region_list, betahat, env) {
  if (length(region_list) < 2L) return(NULL)

  betas <- vapply(region_list, `[[`, numeric(1L), 2L)
  ord   <- order(betas)
  best  <- NULL

  for (idx in seq_len(length(ord) - 1L)) {
    k1 <- ord[idx]
    k2 <- ord[idx + 1L]

    vox1  <- region_list[[k1]][[1]]
    vox2  <- region_list[[k2]][[1]]
    beta1 <- region_list[[k1]][[2]]
    beta2 <- region_list[[k2]][[2]]
    if (length(vox1) == 0L || length(vox2) == 0L) next

    beta_new <- .optimal_merge_beta(vox1, vox2, beta1, beta2, j, env)
    d <- .loss_difference(betahat, j,
                          beta_new - beta1, vox1,
                          beta_new - beta2, vox2, env)
    cand <- list(type = "merge", k1 = k1, k2 = k2,
                 vox1 = vox1, vox2 = vox2,
                 beta1 = beta1, beta2 = beta2,
                 beta_new = beta_new, loss_diff = d)
    best <- .better_operation(best, cand)
  }

  best
}


# =============================================================================
# SPLIT
# =============================================================================

#' @keywords internal
.sira_best_split_candidate <- function(j, region_list, betahat, env) {
  best <- NULL

  for (k in seq_along(region_list)) {
    vox    <- region_list[[k]][[1]]
    beta_k <- region_list[[k]][[2]]
    n_vox  <- length(vox)
    if (n_vox < 2L) next

    grad <- .split_gradients(j, vox, env)
    ord  <- order(grad, decreasing = TRUE)
    n1   <- ceiling(n_vox / 2)
    sub1 <- vox[ord[seq_len(n1)]]
    sub2 <- vox[ord[(n1 + 1L):n_vox]]
    if (length(sub1) == 0L || length(sub2) == 0L) next

    opt <- .optimize_split_values(j, sub1, sub2, beta_k, betahat, env)
    if (is.null(opt)) next

    d <- .loss_difference(betahat, j,
                          opt$b1 - beta_k, sub1,
                          opt$b2 - beta_k, sub2, env)
    cand <- list(type = "split", k = k,
                 sub1 = sub1, sub2 = sub2,
                 beta_old = beta_k, b1 = opt$b1, b2 = opt$b2,
                 loss_diff = d)
    best <- .better_operation(best, cand)
  }

  best
}

#' @keywords internal
.split_gradients <- function(j, vox, env) {
  2 * (env$XTXB[vox, j] - env$XTY_tilde[j, vox]) / env$n
}

#' @keywords internal
.optimize_split_values <- function(j, sub1, sub2, beta_old, betahat, env,
                                   max_cd_iter = 10L) {
  b1 <- beta_old
  b2 <- beta_old
  betahat_cur <- betahat

  for (iter in seq_len(max_cd_iter)) {
    changed <- FALSE

    d1 <- .optimal_region_delta_current(
      j = j, region_voxels = sub1, current_beta = b1,
      base_beta = beta_old, betahat_current = betahat_cur, env = env
    )
    if (is.finite(d1) && abs(d1) > 1e-10) {
      b1 <- b1 + d1
      betahat_cur[sub1] <- b1
      changed <- TRUE
    }

    d2 <- .optimal_region_delta_current(
      j = j, region_voxels = sub2, current_beta = b2,
      base_beta = beta_old, betahat_current = betahat_cur, env = env
    )
    if (is.finite(d2) && abs(d2) > 1e-10) {
      b2 <- b2 + d2
      betahat_cur[sub2] <- b2
      changed <- TRUE
    }

    if (!changed) break
  }

  list(b1 = b1, b2 = b2)
}

#' @keywords internal
.optimal_region_delta_current <- function(j, region_voxels, current_beta,
                                          base_beta, betahat_current, env) {
  if (length(region_voxels) == 0L) return(0)

  region_neighbors <- .get_region_neighbors(region_voxels, env)
  bq <- .bracket_quadratic_setup_custom(
    region_voxels    = region_voxels,
    region_size      = length(region_voxels),
    region_beta      = current_beta,
    base_beta        = base_beta,
    region_neighbors = region_neighbors,
    betahat          = betahat_current,
    which_beta       = j,
    env              = env
  )

  b_flat <- c(bq[[2]][, 1L], bq[[2]][, 2L])
  bracket_quadratic(bq[[1]], b = b_flat)$par
}


# =============================================================================
# OPERATION APPLICATION
# =============================================================================

#' @keywords internal
.sira_apply_operation <- function(j, region_list, op, env) {
  if (op$type == "revalue") {
    .update_XTXB_incremental(j, op$delta, op$vox, env = env)
    env$alphahat_full[j, op$vox] <- op$beta_new
    region_list[[op$k]][[2]] <- op$beta_new
    return(region_list)
  }

  if (op$type == "expand") {
    .update_XTXB_incremental(j, op$beta_new, op$vox, env = env)
    env$alphahat_full[j, op$vox] <- op$beta_new
    region_list[[op$k]][[1]] <- c(region_list[[op$k]][[1]], op$vox)
    return(region_list)
  }

  if (op$type == "shrink") {
    .update_XTXB_incremental(j, op$delta, op$vox, env = env)
    env$alphahat_full[j, op$vox] <- 0
    region_list[[op$k]][[1]] <- setdiff(region_list[[op$k]][[1]], op$vox)
    return(region_list)
  }

  if (op$type == "merge") {
    .update_XTXB_incremental(j, op$beta_new - op$beta1, op$vox1,
                             op$beta_new - op$beta2, op$vox2, env = env)
    merged_vox <- c(op$vox1, op$vox2)
    env$alphahat_full[j, merged_vox] <- op$beta_new
    region_list[[op$k1]] <- list(merged_vox, op$beta_new)
    region_list[[op$k2]] <- NULL
    return(region_list[!vapply(region_list, is.null, logical(1L))])
  }

  if (op$type == "split") {
    .update_XTXB_incremental(j, op$b1 - op$beta_old, op$sub1,
                             op$b2 - op$beta_old, op$sub2, env = env)
    env$alphahat_full[j, op$sub1] <- op$b1
    env$alphahat_full[j, op$sub2] <- op$b2
    region_list[[op$k]] <- list(op$sub1, op$b1)
    region_list[[length(region_list) + 1L]] <- list(op$sub2, op$b2)
    return(region_list)
  }

  stop("Unknown region operation type: ", op$type)
}


# =============================================================================
# HELPERS
# =============================================================================

#' @keywords internal
.sort_regions_by_beta <- function(region_list) {
  if (length(region_list) <= 1L) return(region_list)
  betas <- vapply(region_list, `[[`, numeric(1L), 2L)
  region_list[order(betas)]
}

#' @keywords internal
.bracket_quadratic_setup_custom <- function(region_voxels, region_size,
                                            region_beta, base_beta,
                                            region_neighbors, betahat,
                                            which_beta, env) {
  nb_size <- length(region_neighbors)
  n_rows  <- 2L * nb_size + 3L

  coefs <- matrix(0, nrow = n_rows, ncol = 3L)
  b     <- matrix(0, nrow = n_rows, ncol = 2L)

  lam <- env$lambda
  mu  <- env$mu
  n   <- env$n

  if (nb_size > 0L) {
    nb_diffs <- betahat[region_neighbors] - region_beta

    coefs[seq_len(nb_size), 2L] <- -lam
    coefs[seq_len(nb_size), 3L] <-  lam * nb_diffs

    rows2 <- nb_size + seq_len(nb_size)
    coefs[rows2, 2L] <-  lam
    coefs[rows2, 3L] <- -lam * nb_diffs

    b[seq_len(nb_size),  2L] <-  nb_diffs
    b[rows2,             1L] <-  nb_diffs
    b[seq_len(nb_size),  1L] <- -Inf
    b[rows2,             2L] <-  Inf
  }

  r1 <- 2L * nb_size + 1L
  coefs[r1, 2L] <- -mu * region_size
  coefs[r1, 3L] <- -mu * region_size * region_beta
  b[r1, 1L]     <- -Inf
  b[r1, 2L]     <- -region_beta

  r2 <- 2L * nb_size + 2L
  coefs[r2, 2L] <-  mu * region_size
  coefs[r2, 3L] <-  mu * region_size * region_beta
  b[r2, 1L]     <- -region_beta
  b[r2, 2L]     <-  Inf

  r3 <- 2L * nb_size + 3L
  xtx_jj    <- region_size * env$XTX_p1[which_beta, which_beta]
  xty_term  <- sum(env$XTY_tilde[which_beta, region_voxels])
  xtxb_term <- sum(env$XTXB[region_voxels, which_beta]) +
    region_size * (region_beta - base_beta) * env$XTX_p1[which_beta, which_beta]
  coefs[r3, 1L] <- xtx_jj / n
  coefs[r3, 2L] <- 2 * (xtxb_term - xty_term) / n
  b[r3, 1L]     <- -Inf
  b[r3, 2L]     <- Inf

  list(coefs, b)
}

#' @keywords internal
.optimal_merge_beta <- function(vox1, vox2, beta1, beta2, j, env) {
  vox_m   <- c(vox1, vox2)
  xtx_jj  <- env$XTX_p1[j, j]
  xtx_m   <- length(vox_m) * xtx_jj
  if (xtx_m < 1e-10) return((beta1 + beta2) / 2)

  xtxb_m <- sum(env$XTXB[vox_m, j])
  xty_m  <- sum(env$XTY_tilde[j, vox_m])
  self   <- xtx_jj * (length(vox1) * beta1 + length(vox2) * beta2)

  (xty_m - xtxb_m + self) / xtx_m
}
