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
  inner_iter <- 0L

  while (updated) {
    updated <- FALSE
    betahat <- as.numeric(env$alphahat_full[j, ])

    best_op <- .sira_find_best_operation(j, region_list, betahat, env)
    if (!is.null(best_op) && is.finite(best_op$loss_diff) &&
        best_op$loss_diff < -1e-10) {
      if (isTRUE(env$verbose_ops)) .sira_log_operation(j, best_op)
      region_list <- .sira_apply_operation(j, region_list, best_op, env)
      region_list <- .clean_zero_regions(region_list)
      region_list <- .sort_regions_by_beta(region_list)
      inner_iter <- inner_iter + 1L

      if (inner_iter %% 5L == 0L) {
        before_merge <- length(region_list)
        region_list <- .merge_duplicate_regions(region_list, digits = 3L)
        region_list <- .sort_regions_by_beta(region_list)
        if (isTRUE(env$verbose_ops) && length(region_list) < before_merge) {
          message(sprintf(
            "[SIRA-ops] X%d duplicate-merge: %d -> %d regions after rounding betas",
            j, before_merge, length(region_list)
          ))
        }
      }

      updated <- TRUE
    }
  }

  before_merge <- length(region_list)
  region_list <- .merge_duplicate_regions(region_list, digits = 3L)
  region_list <- .sort_regions_by_beta(region_list)
  if (isTRUE(env$verbose_ops) && length(region_list) < before_merge) {
    message(sprintf(
      "[SIRA-ops] X%d duplicate-merge: %d -> %d regions after rounding betas",
      j, before_merge, length(region_list)
    ))
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
  if (length(occupied) > 0.2 * env$V) return(NULL)

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

    grad <- .split_gradients(j, vox, betahat, env)
    halfway <- floor(n_vox / 2)
    ranks <- rank(grad)
    sub1 <- vox[which(ranks <= halfway)]
    sub2 <- vox[which(ranks > halfway)]
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
.split_gradients <- function(j, vox, betahat, env) {
  vapply(vox, function(v) {
    total <- 2 * (env$XTXB[v, j] - env$XTY_tilde[j, v]) / env$n
    total <- total + env$lambda *
      sum(sign(betahat[v] - betahat[env$neighbor_list[[v]]]))
    total <- total + env$mu * sign(betahat[v])
    total
  }, numeric(1L))
}

#' @keywords internal
.optimize_split_values <- function(j, sub1, sub2, beta_old, betahat, env,
                                   max_cd_iter = 100L) {
  b1 <- beta_old
  b2 <- beta_old
  starting_betahat <- betahat
  split_loss_diff <- Inf
  split_iter <- 0L

  while (abs(split_loss_diff) > 1e-4 && split_iter <= max_cd_iter) {
    temp_betahat <- starting_betahat

    d1 <- .optimal_region_delta_current(
      j = j,
      region_voxels = sub1,
      current_beta = b1,
      betahat_current = temp_betahat,
      env = env
    )
    b1 <- b1 + d1
    temp_betahat[sub1] <- b1

    d2 <- .optimal_region_delta_current(
      j = j,
      region_voxels = sub2,
      current_beta = b2,
      betahat_current = temp_betahat,
      env = env
    )
    b2 <- b2 + d2
    temp_betahat[sub2] <- b2

    split_loss_diff <- .loss_difference(
      starting_betahat, j, d1, sub1, d2, sub2, env
    )
    starting_betahat <- temp_betahat
    split_iter <- split_iter + 1L
  }

  if (split_iter > max_cd_iter) return(NULL)

  if (abs(b1) < 1e-4) b1 <- 0
  if (abs(b2) < 1e-4) b2 <- 0

  list(b1 = b1, b2 = b2)
}

#' @keywords internal
.optimal_region_delta_current <- function(j, region_voxels, current_beta,
                                          betahat_current, env) {
  if (length(region_voxels) == 0L) return(0)

  region_neighbors <- .get_region_neighbors(region_voxels, env)
  bq <- .bracket_quadratic_setup(
    region_voxels    = region_voxels,
    region_size      = length(region_voxels),
    region_beta      = current_beta,
    region_neighbors = region_neighbors,
    betahat          = betahat_current,
    which_beta       = j,
    XTXB_local       = .xtxb_from_betahat(j, betahat_current, env),
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

#' @keywords internal
.xtxb_from_betahat <- function(j, betahat, env) {
  xj <- env$XTX_p1[j, j]
  out <- matrix(0, nrow = env$V, ncol = env$p1)
  out[, j] <- xj * betahat
  out
}

#' @keywords internal
.sira_log_operation <- function(j, op) {
  prefix <- sprintf("[SIRA-ops] X%d", j)

  if (identical(op$type, "revalue")) {
    message(sprintf(
      "%s revalue: region=%d size=%d beta=%.6f loss_diff=%.6g",
      prefix, op$k, length(op$vox), op$beta_new, op$loss_diff
    ))
    return(invisible(NULL))
  }

  if (identical(op$type, "expand")) {
    message(sprintf(
      "%s expand: region=%d voxel=%d beta=%.6f loss_diff=%.6g",
      prefix, op$k, op$vox, op$beta_new, op$loss_diff
    ))
    return(invisible(NULL))
  }

  if (identical(op$type, "shrink")) {
    message(sprintf(
      "%s shrink: region=%d voxel=%d loss_diff=%.6g",
      prefix, op$k, op$vox, op$loss_diff
    ))
    return(invisible(NULL))
  }

  if (identical(op$type, "merge")) {
    message(sprintf(
      "%s merge: regions=(%d,%d) sizes=(%d,%d) beta=%.6f loss_diff=%.6g",
      prefix, op$k1, op$k2, length(op$vox1), length(op$vox2),
      op$beta_new, op$loss_diff
    ))
    return(invisible(NULL))
  }

  if (identical(op$type, "split")) {
    message(sprintf(
      "%s split: region=%d sizes=(%d,%d) betas=(%.6f, %.6f) loss_diff=%.6g",
      prefix, op$k, length(op$sub1), length(op$sub2),
      op$b1, op$b2, op$loss_diff
    ))
    return(invisible(NULL))
  }

  invisible(NULL)
}
