# =============================================================================
# R/utils.R
# Core utility functions used throughout the algorithm:
#
#   - Region list <-> betahat conversion
#   - XTXB management (V x p1 working matrix)
#   - Loss difference computation (the hot path of the algorithm)
#   - bracket_quadratic_setup (builds the piecewise-quadratic inputs)
#
# All functions take `env` as their first argument and read/write state via
# env$variable instead of global <<- assignments.
#
# Key matrices:
#   env$XTX_p1  : p1 x p1  = t(X) %*% X   (fixed; computed once in setup)
#   env$XTXB    : V  x p1  = t(betahat_p1) %*% XTX_p1
#                 where betahat_p1 (p1 x V) = alphahat_full[p1_index, ]
#
# Because voxels are independent in image-on-scalar regression, there are NO
# cross-voxel coupling terms in the loss.  Consequently XTX_p1 is p1 x p1
# (not V x V) and the incremental XTXB update only touches the rows of
# the voxels whose beta values actually changed.
# =============================================================================


# =============================================================================
# REGION LIST HELPERS
# =============================================================================

# A region list for covariate j is a plain R list where each element is a
# two-element list:  list(voxel_indices, beta_value)
# e.g. region_list[[k]] = list(c(1L, 2L, 5L), 3.14)

#' @keywords internal
.region_list_to_betahat <- function(region_list, V) {
  betahat <- numeric(V)
  for (k in seq_along(region_list)) {
    vox <- region_list[[k]][[1]]
    val <- region_list[[k]][[2]]
    if (!is.null(vox) && length(vox) > 0L)
      betahat[vox] <- val
  }
  betahat
}

#' @keywords internal
.voxels_in_regions <- function(region_list) {
  unlist(lapply(region_list, `[[`, 1L), use.names = FALSE)
}

#' @keywords internal
.region_list_is_null <- function(region_list) {
  if (length(region_list) == 0L) return(TRUE)
  all(sapply(region_list, function(r) is.null(r[[1L]]) || r[[2L]] == 0))
}

#' @keywords internal
.clean_zero_regions <- function(region_list) {
  keep <- vapply(region_list,
                 function(r) !is.null(r[[1L]]) && abs(r[[2L]]) > 1e-10,
                 logical(1L))
  region_list[keep]
}

#' @keywords internal
.merge_duplicate_regions <- function(region_list) {
  n <- length(region_list)
  if (n <= 1L) return(region_list)

  betas   <- vapply(region_list, `[[`, numeric(1L), 2L)
  ord     <- order(betas)
  rl      <- region_list[ord]
  betas   <- betas[ord]
  rounded <- round(betas, 3L)
  dups    <- duplicated(rounded)

  i <- n
  while (i >= 2L) {
    if (dups[i]) {
      rl[[i - 1L]][[1L]] <- c(rl[[i - 1L]][[1L]], rl[[i]][[1L]])
      rl[[i]] <- NULL
    }
    i <- i - 1L
  }
  rl
}


# =============================================================================
# XTXB MANAGEMENT
# =============================================================================

#' @keywords internal
.recalibrate_XTXB <- function(env) {
  # Recompute XTXB = t(betahat_p1) %*% XTX_p1 from scratch.
  # Returns a V x p1 matrix; the caller assigns: env$XTXB <- .recalibrate_XTXB(env)
  betahat_mat <- t(env$alphahat_full[env$p1_index, , drop = FALSE])  # V x p1
  betahat_mat %*% env$XTX_p1                                         # V x p1
}

#' @keywords internal
.update_XTXB_incremental <- function(which_beta, c1, v1, c2 = NULL, v2 = NULL, env) {
  # Update env$XTXB after betahat[which_beta, v] shifts by c1 for v in v1
  # (and optionally by c2 for v in v2).
  #
  # Derivation: XTXB[v, j'] = sum_{j''} XTX_p1[j', j''] * betahat[j'', v]
  # So when betahat[which_beta, v] changes by c1:
  #   ΔXTXB[v, j'] = XTX_p1[j', which_beta] * c1   for all j', for v in v1
  # i.e. each row v in v1 gets the same additive vector XTX_p1[, which_beta] * c1.

  if (length(v1) > 0L) {
    delta <- c1 * env$XTX_p1[, which_beta]   # length p1
    env$XTXB[v1, ] <- env$XTXB[v1, , drop = FALSE] +
      matrix(delta, nrow = length(v1), ncol = env$p1, byrow = TRUE)
  }
  if (!is.null(v2) && length(v2) > 0L) {
    delta2 <- c2 * env$XTX_p1[, which_beta]
    env$XTXB[v2, ] <- env$XTXB[v2, , drop = FALSE] +
      matrix(delta2, nrow = length(v2), ncol = env$p1, byrow = TRUE)
  }
}


# =============================================================================
# LOSS DIFFERENCE COMPUTATION
# =============================================================================
#
# The loss function is:
#   L = (1/n) * ||Y - X*beta - Z*gamma - theta*Psi_star||_F^2
#       + lambda * SAD(beta)
#       + mu    * ||beta||_1
#
# Because each voxel's MS contribution is independent, the change caused by
# shifting betahat[j, v1] by c1 (and optionally betahat[j, v2] by c2) is:
#
#   ΔL_MS = (1/n) * [|v1|*c1^2*XTX_jj + 2*c1*sum_{v1}(XTXB[v,j] - XTY_j[v])
#                  + |v2|*c2^2*XTX_jj + 2*c2*sum_{v2}(XTXB[v,j] - XTY_j[v])]
#
# There are NO cross-voxel coupling terms between v1 and v2.

#' @keywords internal
.matrix_sq_norm_diff <- function(which_beta, c1, v1, c2 = NULL, v2 = NULL, env) {
  # Change in ||X * betahat||_F^2 for the proposed shift.
  XTX_jj   <- env$XTX_p1[which_beta, which_beta]   # scalar
  XTXB_col <- env$XTXB[, which_beta]               # length-V vector

  d <- c1^2 * length(v1) * XTX_jj + 2 * c1 * sum(XTXB_col[v1])
  if (!is.null(v2))
    d <- d + c2^2 * length(v2) * XTX_jj + 2 * c2 * sum(XTXB_col[v2])
  d
}

#' @keywords internal
.compute_sad_diff <- function(betahat, c1, v1, c2 = NULL, v2 = NULL, env) {
  mod_voxels <- if (is.null(v2)) v1 else c(v1, v2)
  edge_idx   <- unique(unlist(
    lapply(mod_voxels, function(v)
      env$edge_structure$index[[as.character(v)]]),
    use.names = FALSE
  ))
  if (length(edge_idx) == 0L) return(0)

  delta <- 0
  edges <- env$edge_structure$edges

  for (idx in edge_idx) {
    ei <- edges[idx, 1L];  ej <- edges[idx, 2L]
    orig  <- abs(betahat[ei] - betahat[ej])
    new_i <- betahat[ei] +
      (if (ei %in% v1) c1 else 0) +
      (if (!is.null(v2) && ei %in% v2) c2 else 0)
    new_j <- betahat[ej] +
      (if (ej %in% v1) c1 else 0) +
      (if (!is.null(v2) && ej %in% v2) c2 else 0)
    delta <- delta + abs(new_i - new_j) - orig
  }
  delta
}

#' @keywords internal
.abs_sum_diff <- function(betahat, c1, v1, c2 = NULL, v2 = NULL) {
  d <- sum(abs(betahat[v1] + c1)) - sum(abs(betahat[v1]))
  if (!is.null(v2))
    d <- d + sum(abs(betahat[v2] + c2)) - sum(abs(betahat[v2]))
  d
}

#' @keywords internal
.loss_difference <- function(betahat, which_beta, c1, v1,
                             c2 = NULL, v2 = NULL, env) {
  chg_v <- if (is.null(v2)) v1 else c(v1, v2)
  chg_a <- c(rep(c1, length(v1)),
             if (!is.null(v2)) rep(c2, length(v2)))

  sq_diff  <- .matrix_sq_norm_diff(which_beta, c1, v1, c2, v2, env)
  xty_term <- sum(env$XTY_tilde[which_beta, chg_v] * chg_a)
  ms_diff  <- (sq_diff - 2 * xty_term) / env$n

  sad_diff <- .compute_sad_diff(betahat, c1, v1, c2, v2, env)
  l1_diff  <- .abs_sum_diff(betahat, c1, v1, c2, v2)

  ms_diff + env$lambda * sad_diff + env$mu * l1_diff
}


# =============================================================================
# BRACKET QUADRATIC SETUP
# =============================================================================

#' @keywords internal
.bracket_quadratic_setup <- function(region_voxels, region_size, region_beta,
                                     region_neighbors, betahat,
                                     which_beta, XTXB_local, env) {
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

  # Quadratic (MS) term.
  # xtx_jj = region_size * XTX_p1[j,j]:  the per-voxel quadratic coefficient
  # summed over the region (no cross-voxel coupling in image-on-scalar model).
  r3 <- 2L * nb_size + 3L
  xtx_jj    <- region_size * env$XTX_p1[which_beta, which_beta]
  xty_term  <- sum(env$XTY_tilde[which_beta, region_voxels])
  xtxb_term <- sum(XTXB_local[region_voxels, which_beta])
  coefs[r3, 1L] <- xtx_jj / n
  coefs[r3, 2L] <- 2 * (xtxb_term - xty_term) / n
  b[r3, 1L]     <- -Inf
  b[r3, 2L]     <-  Inf

  list(coefs, b)
}
