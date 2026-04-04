# =============================================================================
# R/sira_batched.R
# Batched preprocessing and fitting for SIRA.
# =============================================================================


#' Precompute reusable batched summaries for SIRA
#'
#' Reads the outcome matrix \eqn{Y} from disk in batches and computes the
#' summary statistics required to run SIRA repeatedly with different tuning
#' parameters without rereading the image data.
#'
#' @param Y_files Character vector of file paths. Each file must be readable by
#'   \code{Y_reader} and return a numeric matrix with \eqn{V} columns.
#' @param X An \eqn{n \times p_1} numeric matrix of covariates of interest.
#' @param Z An \eqn{n \times p_2} numeric matrix of confounders.
#' @param d1 Integer. Number of voxels along image axis 1.
#' @param d2 Integer. Number of voxels along image axis 2.
#' @param d3 Integer. Number of voxels along image axis 3. Use \code{1} for
#'   2D images.
#' @param coords Optional \eqn{V \times 3} numeric matrix of voxel coordinates
#'   for irregular or masked layouts. Supply \code{coords} instead of
#'   \code{d1}, \code{d2}, \code{d3}.
#' @param Y_reader Function with signature \code{function(path) -> matrix}.
#'   Default is \code{readRDS}.
#' @param Psi_star Optional \eqn{L \times V} numeric matrix. Spatial individual
#'   effects basis. If \code{NULL}, SIRA computes it internally.
#' @param delta Numeric. Convergence threshold carried forward to fitting
#'   helpers. Default \code{0.01}.
#' @param inner_loss_tol Optional non-negative numeric. Minimum loss decrease
#'   required to accept a region operation within a covariate update. If
#'   \code{NULL}, defaults to \code{1e-7}.
#' @param verbose Controls progress output. Use \code{FALSE} for silent mode,
#'   \code{TRUE} for iteration-level progress, or \code{"ops"} for detailed
#'   region-operation logging.
#'
#' @return An object of class \code{"sira_batched_preprocessed"} containing
#'   the reusable summaries and metadata required by
#'   \code{\link{sira_batched_fit}}.
#'
#' @examples
#' \dontrun{
#' prep <- sira_batched_preprocess(
#'   Y_files = Y_files,
#'   X = X,
#'   Z = Z,
#'   d1 = d1, d2 = d2, d3 = d3
#' )
#'
#' fit1 <- sira_batched_fit(prep, lambda = 0.1, mu = 0.01)
#' fit2 <- sira_batched_fit(prep, lambda = 0.3, mu = 0.1)
#' }
#'
#' @export
sira_batched_preprocess <- function(Y_files, X, Z,
                                    d1 = NULL, d2 = NULL, d3 = NULL,
                                    coords   = NULL,
                                    Y_reader = readRDS,
                                    Psi_star = NULL,
                                    delta    = 0.01,
                                    inner_loss_tol = NULL,
                                    verbose  = TRUE) {

  meta <- .sira_prepare_batched_inputs(
    Y_files = Y_files,
    X = X,
    Z = Z,
    d1 = d1,
    d2 = d2,
    d3 = d3,
    coords = coords,
    Y_reader = Y_reader,
    Psi_star = Psi_star,
    delta = delta,
    inner_loss_tol = inner_loss_tol,
    verbose = verbose
  )

  env <- meta$env

  if (env$verbose) message("[SIRA-batched] Setting up spatial basis (Psi_star)...")
  if (is.null(meta$Psi_star_input)) {
    Psi_star <- .sira_compute_psi_star(env)
  } else {
    Psi_star <- meta$Psi_star_input
  }
  env$Psi_star <- Psi_star
  env$L        <- nrow(Psi_star)

  if (env$verbose) message("[SIRA-batched] Building neighbor structures...")
  .sira_build_neighbors(env)

  if (env$verbose) message("[SIRA-batched] Accumulating summary statistics from batches...")
  .sira_setup_summary_stats_batched(env, meta$first_batch)

  object <- list(
    X = env$X,
    Z = env$Z,
    n = env$n,
    V = env$V,
    p1 = env$p1,
    p2 = env$p2,
    dims = if (meta$use_grid) c(d1 = env$d1, d2 = env$d2, d3 = env$d3) else NULL,
    coords = env$coords,
    delta = env$delta,
    inner_loss_tol = env$inner_loss_tol,
    verbose = env$verbose,
    verbose_ops = isTRUE(env$verbose_ops),
    Psi_star = env$Psi_star,
    Lambda = env$Lambda,
    L = env$L,
    neighbor_list = env$neighbor_list,
    edge_list = env$edge_list,
    edge_structure = env$edge_structure,
    smoothing_matrix = env$smoothing_matrix,
    XTX_p1 = env$XTX_p1,
    XTY = env$XTY,
    XTZ = env$XTZ,
    ZTZ1ZT = env$ZTZ1ZT,
    ZTZ1ZTY = env$ZTZ1ZTY,
    ZTZ1ZTX = env$ZTZ1ZTX,
    ZTY_raw = env$ZTY_raw,
    YtY = env$YtY,
    Y_Psi_starT = env$Y_Psi_starT,
    scaling_matrix_2 = env$scaling_matrix_2
  )

  class(object) <- "sira_batched_preprocessed"
  object
}


#' Fit SIRA from precomputed batched summaries
#'
#' Runs SIRA using a \code{\link{sira_batched_preprocess}} object. This avoids
#' rereading the outcome data when fitting multiple tuning parameter
#' combinations.
#'
#' @param preprocessed A \code{"sira_batched_preprocessed"} object created by
#'   \code{\link{sira_batched_preprocess}}.
#' @param lambda Non-negative numeric. Spatial penalty controlling boundary
#'   smoothness.
#' @param mu Non-negative numeric. L1 sparsity penalty.
#' @param max_iter Integer. Maximum outer coordinate-descent iterations.
#'   Default \code{50}.
#' @param delta Optional numeric. Overrides the stored convergence threshold.
#'   Default \code{NULL}, meaning use the value stored in \code{preprocessed}.
#' @param inner_loss_tol Optional non-negative numeric. Overrides the stored
#'   minimum accepted loss decrease. If \code{NULL}, uses the value stored in
#'   \code{preprocessed}.
#' @param verbose Optional verbosity setting. Overrides the stored setting from
#'   \code{preprocessed}. Use \code{FALSE}, \code{TRUE}, or \code{"ops"}.
#'
#' @return An object of class \code{"sira"}.
#'
#' @export
sira_batched_fit <- function(preprocessed, lambda, mu,
                             max_iter = 50L,
                             delta    = NULL,
                             inner_loss_tol = NULL,
                             verbose  = NULL) {
  if (!inherits(preprocessed, "sira_batched_preprocessed"))
    stop("preprocessed must be a 'sira_batched_preprocessed' object.")
  if (lambda < 0) stop("lambda must be >= 0.")
  if (mu < 0) stop("mu must be >= 0.")
  if (max_iter < 1L) stop("max_iter must be a positive integer.")

  env <- .sira_env_from_preprocessed(
    preprocessed,
    lambda = lambda,
    mu = mu,
    delta = delta,
    inner_loss_tol = inner_loss_tol,
    verbose = verbose
  )

  if (env$verbose) message("[SIRA-batched] Initializing regions (MUA-based)...")
  env$initial_region_list_full <- .sira_initialize_all_regions(env)

  env$alphahat_full <- matrix(0, nrow = env$p1 + env$p2, ncol = env$V)
  for (j in seq_len(env$p1)) {
    env$alphahat_full[j, ] <- .region_list_to_betahat(
      env$initial_region_list_full[[j]], env$V
    )
  }

  env$thetahat  <- matrix(0, nrow = env$n, ncol = env$L)
  env$gammahat  <- .sira_update_gammahat(env)
  env$alphahat_full[env$p2_index, ] <- env$gammahat
  env$XTY_tilde <- .sira_update_XTY_tilde(env)
  env$XTXB      <- .recalibrate_XTXB(env)

  if (env$verbose) message("[SIRA-batched] Running coordinate-descent algorithm...")
  t0      <- proc.time()
  result  <- .sira_full_algorithm(env, max_iter = as.integer(max_iter))
  runtime <- (proc.time() - t0)[["elapsed"]]

  if (env$verbose)
    message(sprintf("[SIRA-batched] Done. %d iteration(s), %.1f sec.",
                    nrow(result$convergence), runtime))

  betahat  <- env$alphahat_full[env$p1_index, , drop = FALSE]
  gammahat <- env$alphahat_full[env$p2_index, , drop = FALSE]
  if (!is.null(colnames(env$X))) rownames(betahat)  <- colnames(env$X)
  if (!is.null(colnames(env$Z))) rownames(gammahat) <- colnames(env$Z)

  structure(
    list(
      betahat     = betahat,
      gammahat    = gammahat,
      region_list = result$region_list_full,
      convergence = result$convergence,
      lambda      = lambda,
      mu          = mu,
      dims        = preprocessed$dims,
      coords      = preprocessed$coords,
      n           = env$n,
      V           = env$V,
      runtime     = runtime,
      call        = match.call()
    ),
    class = "sira"
  )
}


#' Scalable Image-on-Scalar Regression Algorithm (SIRA) â€” Batched variant
#'
#' Fits the same model as \code{\link{sira}} but reads the outcome matrix
#' \eqn{Y} from disk in batches, so the full \eqn{n \times V} matrix is
#' never held in memory simultaneously. For repeated fits across many
#' \code{(lambda, mu)} combinations, prefer
#' \code{\link{sira_batched_preprocess}} followed by
#' \code{\link{sira_batched_fit}}.
#'
#' @param Y_files Character vector of file paths. Each file must be readable by
#'   \code{Y_reader} and return a numeric matrix with \eqn{V} columns.
#' @param Y_reader Function with signature \code{function(path) -> matrix}.
#'   Default is \code{readRDS}. Use a custom reader for CSV, HDF5, etc.
#' @inheritParams sira
#'
#' @return An object of class \code{"sira"}.
#'
#' @examples
#' \dontrun{
#' fit <- sira_batched(
#'   Y_files = Y_files,
#'   X = X, Z = Z,
#'   d1 = d1, d2 = d2, d3 = d3,
#'   lambda = 0.3, mu = 0.1
#' )
#' }
#'
#' @seealso \code{\link{sira}}, \code{\link{sira_batched_preprocess}},
#'   \code{\link{sira_batched_fit}}
#' @export
sira_batched <- function(Y_files, X, Z,
                         d1 = NULL, d2 = NULL, d3 = NULL,
                         lambda, mu,
                         coords   = NULL,
                         Y_reader = readRDS,
                         Psi_star = NULL,
                         max_iter = 50L,
                         delta    = 0.01,
                         inner_loss_tol = NULL,
                         verbose  = TRUE) {
  prep <- sira_batched_preprocess(
    Y_files = Y_files,
    X = X,
    Z = Z,
    d1 = d1,
    d2 = d2,
    d3 = d3,
    coords = coords,
    Y_reader = Y_reader,
    Psi_star = Psi_star,
    delta = delta,
    inner_loss_tol = inner_loss_tol,
    verbose = verbose
  )

  sira_batched_fit(
    preprocessed = prep,
    lambda = lambda,
    mu = mu,
    max_iter = max_iter,
    delta = delta,
    inner_loss_tol = inner_loss_tol,
    verbose = verbose
  )
}


# =============================================================================
# HELPERS
# =============================================================================

#' @keywords internal
.sira_prepare_batched_inputs <- function(Y_files, X, Z,
                                         d1 = NULL, d2 = NULL, d3 = NULL,
                                         coords   = NULL,
                                         Y_reader = readRDS,
                                         Psi_star = NULL,
                                         delta    = 0.01,
                                         inner_loss_tol = NULL,
                                         verbose  = TRUE) {
  if (!is.character(Y_files) || length(Y_files) == 0L)
    stop("Y_files must be a non-empty character vector of file paths.")
  missing_files <- Y_files[!file.exists(Y_files)]
  if (length(missing_files) > 0L)
    stop("The following Y_files do not exist:\n  ",
         paste(missing_files, collapse = "\n  "))
  if (!is.function(Y_reader))
    stop("Y_reader must be a function with signature function(path) -> matrix.")

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Z)) Z <- as.matrix(Z)
  storage.mode(X) <- "double"
  storage.mode(Z) <- "double"

  n  <- nrow(X)
  p1 <- ncol(X)
  p2 <- ncol(Z)

  if (nrow(Z) != n)
    stop("Z must have nrow(Z) == nrow(X) (", n, " subjects).")

  use_grid   <- !is.null(d1) && !is.null(d2) && !is.null(d3)
  use_coords <- !is.null(coords)
  if (!use_grid && !use_coords)
    stop("Provide either (d1, d2, d3) for grid data or a V x 3 coords matrix for irregularly spaced voxels.")
  if (use_grid && use_coords)
    stop("Provide either (d1, d2, d3) or coords, not both.")

  if (verbose) message("[SIRA-batched] Validating batch files...")
  first_batch <- Y_reader(Y_files[1L])
  if (!is.matrix(first_batch)) first_batch <- as.matrix(first_batch)
  storage.mode(first_batch) <- "double"
  V <- ncol(first_batch)

  if (use_grid && as.integer(d1) * as.integer(d2) * as.integer(d3) != V)
    stop(sprintf("d1 * d2 * d3 = %d but ncol of batch files = %d.",
                 d1 * d2 * d3, V))
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
    storage.mode(Psi_star) <- "double"
  }

  env <- new.env(parent = emptyenv())
  env$n  <- n;   env$V  <- V
  env$p1 <- p1;  env$p2 <- p2
  env$delta   <- delta
  env$inner_max_ops <- 200L
  env$inner_loss_tol <- if (is.null(inner_loss_tol)) 1e-7 else inner_loss_tol
  if (!is.numeric(env$inner_loss_tol) || length(env$inner_loss_tol) != 1L ||
      !is.finite(env$inner_loss_tol) || env$inner_loss_tol < 0) {
    stop("inner_loss_tol must be NULL or a single non-negative finite number.")
  }
  .sira_set_verbose_flags(env, verbose)
  env$X  <- X;   env$Z  <- Z

  if (use_grid) {
    env$d1 <- as.integer(d1); env$d2 <- as.integer(d2); env$d3 <- as.integer(d3)
    env$coords <- NULL
  } else {
    env$d1 <- NULL; env$d2 <- NULL; env$d3 <- NULL
    env$coords <- coords
  }

  env$p1_index <- seq_len(p1)
  env$p2_index <- p1 + seq_len(p2)
  env$Y_files  <- Y_files
  env$Y_reader <- Y_reader

  list(
    env = env,
    first_batch = first_batch,
    use_grid = use_grid,
    Psi_star_input = Psi_star
  )
}

#' @keywords internal
.sira_env_from_preprocessed <- function(preprocessed, lambda, mu,
                                        delta = NULL,
                                        inner_loss_tol = NULL,
                                        verbose = NULL) {
  env <- new.env(parent = emptyenv())

  env$n <- preprocessed$n
  env$V <- preprocessed$V
  env$p1 <- preprocessed$p1
  env$p2 <- preprocessed$p2
  env$d1 <- if (!is.null(preprocessed$dims)) as.integer(preprocessed$dims["d1"]) else NULL
  env$d2 <- if (!is.null(preprocessed$dims)) as.integer(preprocessed$dims["d2"]) else NULL
  env$d3 <- if (!is.null(preprocessed$dims)) as.integer(preprocessed$dims["d3"]) else NULL
  env$coords <- preprocessed$coords
  env$X <- preprocessed$X
  env$Z <- preprocessed$Z
  env$lambda <- lambda
  env$mu <- mu
  env$delta <- if (is.null(delta)) preprocessed$delta else delta
  env$inner_max_ops <- 200L
  env$inner_loss_tol <- if (is.null(inner_loss_tol)) {
    preprocessed$inner_loss_tol
  } else {
    inner_loss_tol
  }
  if (!is.numeric(env$inner_loss_tol) || length(env$inner_loss_tol) != 1L ||
      !is.finite(env$inner_loss_tol) || env$inner_loss_tol < 0) {
    stop("inner_loss_tol must be NULL or a single non-negative finite number.")
  }
  .sira_set_verbose_flags(
    env,
    if (is.null(verbose)) {
      if (isTRUE(preprocessed$verbose_ops)) "ops" else preprocessed$verbose
    } else {
      verbose
    }
  )
  env$p1_index <- seq_len(env$p1)
  env$p2_index <- env$p1 + seq_len(env$p2)
  env$Psi_star <- preprocessed$Psi_star
  env$Lambda <- preprocessed$Lambda
  env$L <- preprocessed$L
  env$neighbor_list <- preprocessed$neighbor_list
  env$edge_list <- preprocessed$edge_list
  env$edge_structure <- preprocessed$edge_structure
  env$smoothing_matrix <- preprocessed$smoothing_matrix
  env$XTX_p1 <- preprocessed$XTX_p1
  env$XTY <- preprocessed$XTY
  env$XTZ <- preprocessed$XTZ
  env$ZTZ1ZT <- preprocessed$ZTZ1ZT
  env$ZTZ1ZTY <- preprocessed$ZTZ1ZTY
  env$ZTZ1ZTX <- preprocessed$ZTZ1ZTX
  env$ZTY_raw <- preprocessed$ZTY_raw
  env$YtY <- preprocessed$YtY
  env$Y_Psi_starT <- preprocessed$Y_Psi_starT
  env$scaling_matrix_2 <- preprocessed$scaling_matrix_2

  env
}


# =============================================================================
# BATCHED SUMMARY STATISTICS
# =============================================================================

#' @keywords internal
.sira_setup_summary_stats_batched <- function(env, first_batch) {
  X        <- env$X
  Z        <- env$Z
  Psi_star <- env$Psi_star
  V        <- env$V
  n        <- env$n
  p1       <- env$p1
  p2       <- env$p2
  n_files  <- length(env$Y_files)

  env$XTX_p1 <- crossprod(X)
  env$XTZ    <- crossprod(X, Z)
  ZTZ        <- crossprod(Z)
  ZTZ1ZT     <- solve(ZTZ, t(Z))
  env$ZTZ1ZT <- ZTZ1ZT

  XTY_acc     <- matrix(0, nrow = p1, ncol = V)
  ZTY_acc     <- matrix(0, nrow = p2, ncol = V)
  YtY_acc     <- numeric(V)
  Y_Psi_starT <- matrix(0, nrow = n, ncol = env$L)

  row_start <- 1L

  for (i in seq_len(n_files)) {
    if (env$verbose)
      message(sprintf("  Batch %d / %d ...", i, n_files))

    Y_b <- if (i == 1L) first_batch else {
      yb <- env$Y_reader(env$Y_files[i])
      if (!is.matrix(yb)) yb <- as.matrix(yb)
      storage.mode(yb) <- "double"
      yb
    }

    n_b <- nrow(Y_b)
    if (ncol(Y_b) != V)
      stop(sprintf("Batch %d has %d columns but V = %d.", i, ncol(Y_b), V))

    row_end <- row_start + n_b - 1L
    if (row_end > n)
      stop(sprintf(
        "Batch files contain more rows than nrow(X) = %d. Check that X/Z and Y_files correspond to the same subjects.",
        n
      ))

    X_b <- X[row_start:row_end, , drop = FALSE]
    Z_b <- Z[row_start:row_end, , drop = FALSE]

    XTY_acc <- XTY_acc + crossprod(X_b, Y_b)
    ZTY_acc <- ZTY_acc + crossprod(Z_b, Y_b)
    YtY_acc <- YtY_acc + colSums(Y_b^2)
    Y_Psi_starT[row_start:row_end, ] <- Y_b %*% t(Psi_star)

    row_start <- row_end + 1L
    rm(Y_b); gc(verbose = FALSE)
  }

  if (row_start - 1L != n)
    stop(sprintf(
      "Batch files contain %d rows total but nrow(X) = %d.",
      row_start - 1L, n
    ))

  env$ZTZ1ZTY <- solve(ZTZ, ZTY_acc)
  env$ZTZ1ZTX <- solve(ZTZ, crossprod(Z, X))
  env$ZTY_raw <- ZTY_acc
  env$YtY <- YtY_acc
  env$Y_Psi_starT <- Y_Psi_starT
  env$scaling_matrix_2 <- if (!is.null(env$Lambda)) {
    diag(1 / env$Lambda)
  } else {
    solve(tcrossprod(Psi_star))
  }
  env$XTY <- XTY_acc

  if (env$verbose)
    message("  Batched summary statistics done.")
}
