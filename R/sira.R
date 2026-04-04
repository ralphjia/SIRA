# =============================================================================
# R/sira.R
# Scalable Image-on-Scalar Regression Algorithm (SIRA)
# Main exported function and S3 methods.
# =============================================================================


#' Scalable Image-on-Scalar Regression Algorithm (SIRA)
#'
#' Fits an image-on-scalar regression model using the Scalable Image-on-Scalar
#' Regression Algorithm (SIRA). The model is:
#'
#' SIRA initializes each coefficient image with a mass-univariate estimate,
#' smooths the resulting t-statistics, uses \code{locfdr} to define positive
#' and negative starting regions, and then repeatedly applies the single
#' best loss-decreasing region operation among revalue, expand, shrink, merge,
#' and split. Nuisance coefficients and spatial random-effect scores are
#' updated in closed form between region sweeps.
#'
#' \deqn{Y_i(s_v) = \sum_{j=1}^{p_1} \beta_j(s_v) X_{ij} +
#'   \sum_{k=1}^{p_2} \gamma_k(s_v) Z_{ik} + \eta_i(s_v) + \epsilon_i(s_v)}
#'
#' SIRA enforces sparsity and spatial homogeneity on the coefficient images
#' \eqn{\beta_j} using a coordinate-descent algorithm with five region
#' operations (expand, shrink, revalue, merge, split). Merge candidates are
#' restricted to adjacent coefficient values after ranking the current regions,
#' and split candidates are constructed from voxelwise loss gradients. The
#' implementation uses precomputed summary statistics to remain scalable to
#' high-dimensional neuroimaging data.
#'
#' @param Y An \eqn{n \times V} numeric matrix of image outcomes. Each row is
#'   one subject's vectorized image (column-major / Fortran order).
#' @param X An \eqn{n \times p_1} numeric matrix of covariates of interest.
#' @param Z An \eqn{n \times p_2} numeric matrix of confounders. Should
#'   typically include a column of ones for the intercept.
#' @param d1 Integer. Number of voxels along image axis 1.
#' @param d2 Integer. Number of voxels along image axis 2.
#' @param d3 Integer. Number of voxels along image axis 3. Use \code{1} for
#'   2D images.
#' @param lambda Non-negative numeric. Spatial penalty controlling boundary
#'   smoothness. Larger values produce larger, blockier regions.
#' @param mu Non-negative numeric. L1 sparsity penalty. Larger values produce
#'   sparser solutions (fewer nonzero voxels).
#' @param coords Optional \eqn{V \times 3} numeric matrix of voxel coordinates
#'   for irregular or masked layouts. Supply \code{coords} instead of
#'   \code{d1}, \code{d2}, \code{d3} when voxels do not lie on a full
#'   rectangular grid.
#' @param Psi_star Optional \eqn{L \times V} numeric matrix. Spatial individual
#'   effects basis computed by \code{BayesGPfit}. If \code{NULL} (default),
#'   SIRA calls \code{BayesGPfit} internally with \code{poly_degree = 8}.
#' @param max_iter Integer. Maximum outer coordinate-descent iterations.
#'   Default \code{50}.
#' @param delta Numeric. Convergence threshold on the Frobenius norm of the
#'   change in \eqn{\hat\alpha} between iterations. Default \code{0.01}.
#' @param verbose Controls progress output. Use \code{FALSE} for silent mode,
#'   \code{TRUE} for iteration-level progress, or \code{"ops"} for detailed
#'   region-operation logging.
#'
#' @return An object of class \code{"sira"} containing:
#' \describe{
#'   \item{\code{betahat}}{A \eqn{p_1 \times V} matrix of coefficient image
#'     estimates. Row \eqn{j} is the vectorized image for covariate \eqn{j}.}
#'   \item{\code{gammahat}}{A \eqn{p_2 \times V} matrix of confounder
#'     coefficient estimates.}
#'   \item{\code{region_list}}{List of length \eqn{p_1}. Element \eqn{j} is
#'     the final region configuration for covariate \eqn{j}: a list of
#'     \code{list(voxel_indices, beta_value)} pairs.}
#'   \item{\code{convergence}}{Data frame with one row per outer iteration:
#'     \code{iteration}, \code{num_regions}, \code{num_voxels},
#'     \code{alpha_norm_diff}.}
#'   \item{\code{lambda}, \code{mu}}{Tuning parameters used.}
#'   \item{\code{dims}}{Named integer vector \code{c(d1, d2, d3)}.}
#'   \item{\code{n}, \code{V}}{Number of subjects and voxels.}
#'   \item{\code{runtime}}{Fitting time in seconds.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' d1 <- d2 <- d3 <- 20L
#' V  <- d1 * d2 * d3
#' n  <- 100L
#'
#' beta_true       <- numeric(V)
#' beta_true[1:64] <- 3   # a small cube near the origin
#'
#' X <- matrix(rnorm(n), ncol = 1)
#' Z <- cbind(1, rnorm(n))
#' Y <- X %*% t(beta_true) +
#'      Z %*% matrix(rnorm(2 * V), nrow = 2) +
#'      matrix(rnorm(n * V, sd = 3), nrow = n)
#'
#' fit <- sira(Y = Y, X = X, Z = Z,
#'             d1 = d1, d2 = d2, d3 = d3,
#'             lambda = 0.3, mu = 0.1)
#' print(fit)
#' }
#'
#' @seealso \code{\link{coef.sira}}, \code{\link{print.sira}},
#'   \code{\link{summary.sira}}
#'
#' @importFrom Matrix Matrix sparseMatrix
#' @importFrom locfdr locfdr
#' @export
sira <- function(Y, X, Z,
                 d1 = NULL, d2 = NULL, d3 = NULL,
                 lambda, mu,
                 coords   = NULL,
                 Psi_star = NULL,
                 max_iter = 50L,
                 delta    = 0.01,
                 verbose  = TRUE) {

  # ---- 1. Coerce and validate inputs ------------------------------------
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  storage.mode(Y) <- "double"
  storage.mode(X) <- "double"
  storage.mode(Z) <- "double"

  n  <- nrow(Y);  V  <- ncol(Y)
  p1 <- ncol(X);  p2 <- ncol(Z)

  if (nrow(X) != n)
    stop("X must have nrow(X) == nrow(Y) (", n, " subjects).")
  if (nrow(Z) != n)
    stop("Z must have nrow(Z) == nrow(Y) (", n, " subjects).")
  if (lambda < 0) stop("lambda must be >= 0.")
  if (mu     < 0) stop("mu must be >= 0.")
  if (max_iter < 1L) stop("max_iter must be a positive integer.")

  # Exactly one of (d1, d2, d3) or coords must be supplied.
  use_grid   <- !is.null(d1) && !is.null(d2) && !is.null(d3)
  use_coords <- !is.null(coords)
  if (!use_grid && !use_coords)
    stop("Provide either (d1, d2, d3) for grid data or a V x 3 coords ",
         "matrix for irregularly spaced voxels.")
  if (use_grid && use_coords)
    stop("Provide either (d1, d2, d3) or coords, not both.")
  if (use_grid && as.integer(d1) * as.integer(d2) * as.integer(d3) != V)
    stop(sprintf("d1 * d2 * d3 = %d but ncol(Y) = %d.", d1 * d2 * d3, V))
  if (use_coords) {
    if (!is.matrix(coords)) coords <- as.matrix(coords)
    if (nrow(coords) != V || ncol(coords) != 3L)
      stop(sprintf("coords must be a V x 3 matrix (V = %d).", V))
    storage.mode(coords) <- "double"
  }

  if (!is.null(Psi_star)) {
    if (!is.matrix(Psi_star)) Psi_star <- as.matrix(Psi_star)
    if (ncol(Psi_star) != V)
      stop("Psi_star must have ncol(Psi_star) == V (", V, " voxels).")
  }

  # ---- 2. Create a fresh local state environment ------------------------
  #  All internal functions receive `env` instead of using <<-.
  #  This makes each sira() call fully self-contained.
  env <- new.env(parent = emptyenv())
  env$n  <- n;   env$V  <- V
  env$p1 <- p1;  env$p2 <- p2
  env$d1 <- as.integer(d1)
  env$d2 <- as.integer(d2)
  env$d3 <- as.integer(d3)
  env$coords <- if (use_coords) coords else NULL
  env$X  <- X;   env$Y  <- Y;   env$Z <- Z
  env$lambda  <- lambda
  env$mu      <- mu
  env$delta   <- delta
  .sira_set_verbose_flags(env, verbose)

  # Partition column indices: rows 1:p1 are beta (interest),
  # rows (p1+1):(p1+p2) are gamma (confounders).
  env$p1_index <- seq_len(p1)
  env$p2_index <- p1 + seq_len(p2)

  # ---- 3. Spatial individual effects: Psi_star --------------------------
  if (env$verbose) message("[SIRA] Setting up spatial basis (Psi_star)...")
  if (is.null(Psi_star)) {
    Psi_star <- .sira_compute_psi_star(env)
  }
  env$Psi_star <- Psi_star        # L x V
  env$L        <- nrow(Psi_star)

  # ---- 4. Neighbor structures -------------------------------------------
  if (env$verbose) message("[SIRA] Building neighbor structures...")
  .sira_build_neighbors(env)      # sets env$neighbor_list, env$edge_list,
  #      env$edge_structure, env$smoothing_matrix

  # ---- 5. Summary statistics --------------------------------------------
  if (env$verbose) message("[SIRA] Computing summary statistics...")
  .sira_setup_summary_stats(env)  # sets env$XTY, env$XTZ, env$ZTZ1ZT,
  #      env$ZTZ1ZTY, env$ZTZ1ZTX,
  #      env$Y_Psi_starT, env$scaling_matrix_2

  # ---- 6. Initialization ------------------------------------------------
  if (env$verbose) message("[SIRA] Initializing regions (MUA-based)...")
  env$initial_region_list_full <- .sira_initialize_all_regions(env)

  # Build alphahat_full: (p1+p2) x V matrix, rows = covariates then confounders
  env$alphahat_full <- matrix(0, nrow = p1 + p2, ncol = V)
  for (j in seq_len(p1)) {
    env$alphahat_full[j, ] <- .region_list_to_betahat(
      env$initial_region_list_full[[j]], V
    )
  }

  # Initialize theta, gamma, XTY_tilde, XTXB
  env$thetahat  <- matrix(0, nrow = n, ncol = env$L)   # n x L
  env$gammahat  <- .sira_update_gammahat(env)           # p2 x V
  env$alphahat_full[env$p2_index, ] <- env$gammahat
  env$XTY_tilde <- .sira_update_XTY_tilde(env)         # p1 x V
  # XTXB (V x p1) = t(betahat_p1) %*% XTX_p1
  # Rows are voxels, cols are covariates of interest.
  # Kept transposed relative to alphahat_full for fast column look-ups
  # in the region operations (see region_ops.R).
  env$XTXB <- .recalibrate_XTXB(env)   # V x p1

  # ---- 7. Main algorithm ------------------------------------------------
  if (env$verbose) message("[SIRA] Running coordinate-descent algorithm...")
  t0     <- proc.time()
  result <- .sira_full_algorithm(env, max_iter = as.integer(max_iter))
  runtime <- (proc.time() - t0)[["elapsed"]]

  if (env$verbose)
    message(sprintf("[SIRA] Done. %d iteration(s), %.1f sec.",
                    nrow(result$convergence), runtime))

  # ---- 8. Build return object -------------------------------------------
  betahat  <- env$alphahat_full[env$p1_index, , drop = FALSE]  # p1 x V
  gammahat <- env$alphahat_full[env$p2_index, , drop = FALSE]  # p2 x V

  if (!is.null(colnames(X))) rownames(betahat)  <- colnames(X)
  if (!is.null(colnames(Z))) rownames(gammahat) <- colnames(Z)

  structure(
    list(
      betahat     = betahat,
      gammahat    = gammahat,
      region_list = result$region_list_full,
      convergence = result$convergence,
      lambda      = lambda,
      mu          = mu,
      dims        = if (use_grid) c(d1 = d1, d2 = d2, d3 = d3) else NULL,
      coords      = if (use_coords) coords else NULL,
      n           = n,
      V           = V,
      runtime     = runtime,
      call        = match.call()
    ),
    class = "sira"
  )
}


# =============================================================================
# Internal: compute Psi_star via BayesGPfit
# =============================================================================

#' @keywords internal
.sira_compute_psi_star <- function(env) {
  if (!requireNamespace("BayesGPfit", quietly = TRUE))
    stop("BayesGPfit is required when Psi_star = NULL. ",
         "Install it with: install.packages('BayesGPfit')")

  if (!requireNamespace("pracma", quietly = TRUE))
    stop("pracma is required for Gram-Schmidt orthogonalization. ",
         "Install it with: install.packages('pracma')")

  V <- env$V

  # Build normalized coordinate grid in [-1, 1]^3.
  # Two paths: use supplied coords (irregular layout) or build from d1/d2/d3.
  if (!is.null(env$coords)) {
    # Normalize each axis of the supplied coordinates to [-1, 1].
    v_list <- apply(env$coords, 2L, function(ax) {
      rng <- range(ax)
      if (diff(rng) < 1e-8) rep(0, length(ax))   # degenerate axis
      else (ax - rng[1L]) / diff(rng) * 2 - 1
    })
    if (!is.matrix(v_list)) v_list <- matrix(v_list, ncol = 3L)
  } else {
    d1 <- env$d1; d2 <- env$d2; d3 <- env$d3
    x1 <- rep(rep(seq(-1, 1, length.out = d1), times = d2), times = d3)
    x2 <- rep(rep(seq(-1, 1, length.out = d2), each  = d1), times = d3)
    x3 <- rep(seq( 1, -1, length.out = d3),    each  = d1 * d2)
    v_list <- cbind(x1, x2, x3)   # V x 3
  }

  # poly_degree = 8 gives L = C(11,3) = 165 basis functions for 3D images,
  # matching the paper's simulation settings (Kang, 2022).
  poly_degree <- 8L
  a <- 0.01; b <- 1

  # Step 1: evaluate basis functions at grid points  (V x L)
  Psi_raw <- BayesGPfit::GP.eigen.funcs.fast(v_list,
                                             poly_degree = poly_degree,
                                             a = a, b = b)

  # Step 2: orthogonalize columns via Gram-Schmidt  (V x L)
  Psi_orth <- pracma::gramSchmidt(Psi_raw)$Q

  # Step 3: eigenvalues  (length L)
  Lambda <- BayesGPfit::GP.eigen.value(poly_degree = poly_degree,
                                       a = a, b = b,
                                       d = ncol(v_list))
  sqrt_Lambda <- sqrt(Lambda)

  # Step 4: Psi_star = diag(sqrt_Lambda) %*% t(Psi_orth)   (L x V)
  Psi_star <- diag(sqrt_Lambda) %*% t(Psi_orth)

  # Store Lambda so summary_stats can build scaling_matrix_2 = diag(1/Lambda)
  # cheaply, without inverting the full L x L matrix.
  env$Lambda <- Lambda

  if (ncol(Psi_star) != V)
    stop("Psi_star has unexpected dimensions after construction.")

  Psi_star
}


# =============================================================================
# Internal: parameter update helpers (wrappers used in full_algorithm)
# =============================================================================

#' @keywords internal
.sira_update_gammahat <- function(env) {
  # Closed-form least-squares update for gamma, with beta and theta fixed.
  # gammahat = (ZTZ)^{-1} ZT (Y - X betahat - thetahat Psi_star)
  # Using precomputed ZTZ1ZT = (ZTZ)^{-1} ZT  [p2 x n]
  betahat_p1 <- env$alphahat_full[env$p1_index, , drop = FALSE]   # p1 x V
  env$ZTZ1ZTY -
    env$ZTZ1ZTX %*% betahat_p1 -
    (env$ZTZ1ZT %*% env$thetahat) %*% env$Psi_star
}

#' @keywords internal
.sira_update_thetahat <- function(env) {
  # Closed-form update for theta, with beta and gamma fixed.
  # thetahat = (Y - X betahat Psi_star^T - Z gammahat Psi_star^T)
  #            %*% scaling_matrix_2
  # scaling_matrix_2 = (Psi_star Psi_star^T)^{-1}  [L x L]
  betahat_p1 <- env$alphahat_full[env$p1_index, , drop = FALSE]   # p1 x V
  gammahat   <- env$alphahat_full[env$p2_index, , drop = FALSE]   # p2 x V
  residual   <- env$Y_Psi_starT -                                  # n x L
    env$X %*% tcrossprod(betahat_p1, env$Psi_star) -
    env$Z %*% tcrossprod(gammahat,   env$Psi_star)
  residual %*% env$scaling_matrix_2                                # n x L
}

#' @keywords internal
.sira_update_XTY_tilde <- function(env) {
  # XTY_tilde = X^T (Y - Z gammahat - thetahat Psi_star)
  # Used as the working response in the beta region-finding step.
  # Returns a p1 x V matrix (rows = covariates of interest).
  env$XTY -
    env$XTZ %*% env$gammahat -
    crossprod(env$X, env$thetahat) %*% env$Psi_star
}

#' @keywords internal
.sira_set_verbose_flags <- function(env, verbose) {
  if (is.logical(verbose)) {
    env$verbose <- isTRUE(verbose)
    env$verbose_ops <- FALSE
    return(invisible(NULL))
  }

  if (is.numeric(verbose) && length(verbose) == 1L) {
    env$verbose <- verbose >= 1
    env$verbose_ops <- verbose >= 2
    return(invisible(NULL))
  }

  if (is.character(verbose) && length(verbose) == 1L) {
    key <- tolower(verbose)
    env$verbose <- key %in% c("true", "yes", "all", "ops", "detail", "detailed")
    env$verbose_ops <- key %in% c("ops", "all", "detail", "detailed")
    return(invisible(NULL))
  }

  stop("verbose must be FALSE, TRUE, a numeric level, or 'ops'.")
}


# =============================================================================
# S3 methods
# =============================================================================

#' Print a SIRA fit
#'
#' @param x A \code{"sira"} object.
#' @param ... Ignored.
#' @export
print.sira <- function(x, ...) {
  cat("SIRA fit\n")
  if (!is.null(x$dims)) {
    cat(sprintf("  Image dims : %d x %d x %d  (%d voxels)\n",
                x$dims["d1"], x$dims["d2"], x$dims["d3"], x$V))
  } else {
    cat(sprintf("  Voxels     : %d  (irregular/masked layout)\n", x$V))
  }
  cat(sprintf("  Subjects   : %d\n", x$n))
  cat(sprintf("  Covariates : p1 = %d (interest),  p2 = %d (confounders)\n",
              nrow(x$betahat), nrow(x$gammahat)))
  cat(sprintf("  Tuning     : lambda = %.4g,  mu = %.4g\n",
              x$lambda, x$mu))
  cat("  Active regions per covariate:\n")
  for (j in seq_len(nrow(x$betahat))) {
    nz   <- sum(x$betahat[j, ] != 0)
    nreg <- length(x$region_list[[j]])
    lbl  <- if (!is.null(rownames(x$betahat))) rownames(x$betahat)[j] else
      paste0("X", j)
    cat(sprintf("    %-12s : %5d voxels in %d region(s)\n", lbl, nz, nreg))
  }
  cat(sprintf("  Converged  : %d iteration(s),  %.2f sec\n",
              nrow(x$convergence), x$runtime))
  invisible(x)
}

#' Summarize a SIRA fit
#'
#' Prints the call, the brief \code{print} output, and the full convergence
#' history.
#'
#' @param object A \code{"sira"} object.
#' @param ... Ignored.
#' @export
summary.sira <- function(object, ...) {
  cat("=== SIRA Summary ===\n\n")
  cat("Call:\n  ")
  print(object$call)
  cat("\n")
  print(object)
  cat("\nConvergence history:\n")
  print(object$convergence)
  invisible(object)
}

#' Extract coefficient images from a SIRA fit
#'
#' @param object A \code{"sira"} object.
#' @param type \code{"beta"} (default) for the covariate-of-interest coefficient
#'   images (\eqn{p_1 \times V}), or \code{"gamma"} for the confounder
#'   estimates (\eqn{p_2 \times V}).
#' @param ... Ignored.
#' @return A numeric matrix.
#' @export
coef.sira <- function(object, type = c("beta", "gamma"), ...) {
  switch(match.arg(type),
         beta  = object$betahat,
         gamma = object$gammahat
  )
}

#' Reshape a coefficient row into a 3D array
#'
#' Convenience function to turn row \code{j} of \code{fit$betahat} (a length-V
#' vector) back into a \eqn{d_1 \times d_2 \times d_3} array for visualisation.
#'
#' @param x A \code{"sira"} object.
#' @param j Integer. Which covariate of interest to extract (default \code{1}).
#' @param type \code{"beta"} or \code{"gamma"}.
#' @param ... Ignored.
#' @return A 3D numeric array with dimensions \code{c(d1, d2, d3)}.
#' @export
as.array.sira <- function(x, j = 1L, type = c("beta", "gamma"), ...) {
  if (is.null(x$dims))
    stop("as.array() is not supported for fits with irregular voxel layouts ",
         "(coords were supplied instead of d1/d2/d3). ",
         "Use coef() to extract the length-V coefficient vector instead.")
  mat  <- stats::coef(x, type = match.arg(type))
  vec  <- mat[j, ]
  array(vec, dim = x$dims)
}
