# =============================================================================
# R/sira_batched.R
# Out-of-core variant of sira() for datasets where Y is too large to hold
# in memory (e.g. UK Biobank with ~35,000 subjects x ~500,000 voxels).
#
# Instead of passing Y as a matrix, the caller passes a character vector of
# file paths.  Each file is read by Y_reader (default: readRDS) and must
# return an n_b x V numeric matrix for that batch of subjects.  The files
# are read once during setup to accumulate summary statistics, then discarded.
# The main coordinate-descent algorithm is identical to sira().
#
# Memory footprint during setup (UK Biobank example, p1=5, p2=3, L=165):
#   XTY          : p1 x V  ≈  5 x 500K  = 2.5M doubles  (~20 MB)
#   ZTY_raw      : p2 x V  ≈  3 x 500K  = 1.5M doubles  (~12 MB)
#   YtY          : V       ≈  500K       = 0.5M doubles  (~4 MB)
#   Y_Psi_starT  : n x L   ≈  35K x 165 = 5.8M doubles  (~46 MB)
#   One Y batch  : n_b x V ≈  700 x 500K = 350M doubles  (~2.8 GB peak)
# =============================================================================


#' Scalable Image-on-Scalar Regression Algorithm (SIRA) — Batched variant
#'
#' Fits the same model as \code{\link{sira}} but reads the outcome matrix
#' \eqn{Y} from disk in batches, so the full \eqn{n \times V} matrix is
#' never held in memory simultaneously.
#'
#' @param Y_files Character vector of file paths.  Each file must be readable
#'   by \code{Y_reader} and return a numeric matrix with \eqn{V} columns.
#'   Rows (subjects) are assumed to be in the same order as the rows of
#'   \code{X} and \code{Z}.  The batch sizes (number of rows per file) may
#'   differ.
#' @param Y_reader Function with signature \code{function(path) -> matrix}.
#'   Default is \code{readRDS}.  Use a custom reader for CSV, HDF5, etc.
#' @inheritParams sira
#'
#' @return An object of class \code{"sira"} — identical structure to the
#'   return value of \code{\link{sira}}.
#'
#' @examples
#' \dontrun{
#' # Split a large Y matrix into batches and save as RDS files
#' batch_size <- 700L
#' batches    <- split(seq_len(nrow(Y)),
#'                     ceiling(seq_len(nrow(Y)) / batch_size))
#' Y_files <- file.path(tempdir(), paste0("Y_batch_", seq_along(batches), ".rds"))
#' for (i in seq_along(batches))
#'   saveRDS(Y[batches[[i]], ], Y_files[i])
#'
#' fit <- sira_batched(Y_files = Y_files,
#'                     X = X, Z = Z,
#'                     d1 = d1, d2 = d2, d3 = d3,
#'                     lambda = 0.3, mu = 0.1)
#' }
#'
#' @seealso \code{\link{sira}}
#' @export
sira_batched <- function(Y_files, X, Z,
                         d1 = NULL, d2 = NULL, d3 = NULL,
                         lambda, mu,
                         coords   = NULL,
                         Y_reader = readRDS,
                         Psi_star = NULL,
                         max_iter = 50L,
                         delta    = 0.01,
                         verbose  = TRUE) {

  # ---- 1. Validate non-Y inputs -----------------------------------------
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
  if (lambda < 0) stop("lambda must be >= 0.")
  if (mu     < 0) stop("mu must be >= 0.")
  if (max_iter < 1L) stop("max_iter must be a positive integer.")

  use_grid   <- !is.null(d1) && !is.null(d2) && !is.null(d3)
  use_coords <- !is.null(coords)
  if (!use_grid && !use_coords)
    stop("Provide either (d1, d2, d3) for grid data or a V x 3 coords ",
         "matrix for irregularly spaced voxels.")
  if (use_grid && use_coords)
    stop("Provide either (d1, d2, d3) or coords, not both.")

  # ---- 2. Peek at the first file to get V --------------------------------
  if (verbose) message("[SIRA-batched] Validating batch files...")
  first_batch <- Y_reader(Y_files[1L])
  if (!is.matrix(first_batch)) first_batch <- as.matrix(first_batch)
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
  }

  # ---- 3. Build env ------------------------------------------------------
  env <- new.env(parent = emptyenv())
  env$n  <- n;   env$V  <- V
  env$p1 <- p1;  env$p2 <- p2
  env$lambda  <- lambda
  env$mu      <- mu
  env$delta   <- delta
  env$verbose <- isTRUE(verbose)
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

  # ---- 4. Spatial basis --------------------------------------------------
  if (env$verbose) message("[SIRA-batched] Setting up spatial basis (Psi_star)...")
  if (is.null(Psi_star)) {
    Psi_star <- .sira_compute_psi_star(env)
  }
  env$Psi_star <- Psi_star
  env$L        <- nrow(Psi_star)

  # ---- 5. Neighbor structures --------------------------------------------
  if (env$verbose) message("[SIRA-batched] Building neighbor structures...")
  .sira_build_neighbors(env)

  # ---- 6. Batched summary statistics ------------------------------------
  if (env$verbose) message("[SIRA-batched] Accumulating summary statistics from batches...")
  .sira_setup_summary_stats_batched(env, first_batch)
  # Sets: env$XTX_p1, env$XTY, env$XTZ, env$ZTZ1ZT, env$ZTZ1ZTY,
  #       env$Y_Psi_starT, env$scaling_matrix_2,
  #       env$XTY_tilde_0, env$YtY_tilde   (used by batched initialiser)

  # ---- 7. Initialisation (uses precomputed stats, no Y needed) ----------
  if (env$verbose) message("[SIRA-batched] Initializing regions (MUA-based)...")
  env$initial_region_list_full <- .sira_initialize_all_regions(env)

  env$alphahat_full <- matrix(0, nrow = p1 + p2, ncol = V)
  for (j in seq_len(p1)) {
    env$alphahat_full[j, ] <- .region_list_to_betahat(
      env$initial_region_list_full[[j]], V
    )
  }

  env$thetahat  <- matrix(0, nrow = n, ncol = env$L)
  env$gammahat  <- .sira_update_gammahat(env)
  env$alphahat_full[env$p2_index, ] <- env$gammahat
  env$XTY_tilde <- .sira_update_XTY_tilde(env)
  env$XTXB      <- .recalibrate_XTXB(env)

  # ---- 8. Main algorithm (identical to sira()) --------------------------
  if (env$verbose) message("[SIRA-batched] Running coordinate-descent algorithm...")
  t0      <- proc.time()
  result  <- .sira_full_algorithm(env, max_iter = as.integer(max_iter))
  runtime <- (proc.time() - t0)[["elapsed"]]

  if (env$verbose)
    message(sprintf("[SIRA-batched] Done. %d iteration(s), %.1f sec.",
                    nrow(result$convergence), runtime))

  # ---- 9. Return object (identical structure to sira()) -----------------
  betahat  <- env$alphahat_full[env$p1_index, , drop = FALSE]
  gammahat <- env$alphahat_full[env$p2_index, , drop = FALSE]
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
# BATCHED SUMMARY STATISTICS
# =============================================================================

#' @keywords internal
.sira_setup_summary_stats_batched <- function(env, first_batch) {
  # Reads Y_files one at a time, accumulates summary statistics, discards
  # each batch before reading the next.  first_batch is the already-read
  # first file (from the validation step in sira_batched), passed in to
  # avoid reading it twice.
  #
  # Quantities accumulated across batches (running sums):
  #   XTY      : p1 x V  = t(X) %*% Y
  #   ZTY_raw  : p2 x V  = t(Z) %*% Y   (→ ZTZ1ZTY after the loop)
  #   YtY      : V       = colSums(Y^2) (for t-stat sigma^2 in initialiser)
  #
  # Assembled row-by-row (manageable n x L):
  #   Y_Psi_starT : n x L = Y %*% t(Psi_star)

  X        <- env$X          # n x p1
  Z        <- env$Z          # n x p2
  Psi_star <- env$Psi_star   # L x V
  V        <- env$V
  n        <- env$n
  p1       <- env$p1
  p2       <- env$p2
  n_files  <- length(env$Y_files)

  # Fixed cross-products (don't need Y)
  env$XTX_p1 <- crossprod(X)    # p1 x p1
  env$XTZ    <- crossprod(X, Z) # p1 x p2
  ZTZ        <- crossprod(Z)    # p2 x p2
  ZTZ1ZT    <- solve(ZTZ, t(Z)) # p2 x n
  env$ZTZ1ZT <- ZTZ1ZT

  # Initialise accumulators
  XTY_acc     <- matrix(0, nrow = p1, ncol = V)
  ZTY_acc     <- matrix(0, nrow = p2, ncol = V)
  YtY_acc     <- numeric(V)
  Y_Psi_starT <- matrix(0, nrow = n, ncol = env$L)

  # Determine row ranges for X and Z slicing
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
        "Batch files contain more rows than nrow(X) = %d. ",
        "Check that X/Z and Y_files correspond to the same subjects.", n))

    X_b <- X[row_start:row_end, , drop = FALSE]   # n_b x p1
    Z_b <- Z[row_start:row_end, , drop = FALSE]   # n_b x p2

    # Accumulate
    XTY_acc <- XTY_acc + crossprod(X_b, Y_b)       # p1 x V
    ZTY_acc <- ZTY_acc + crossprod(Z_b, Y_b)       # p2 x V
    YtY_acc <- YtY_acc + colSums(Y_b^2)            # V

    # Assemble Y_Psi_starT row block  (n_b x L)
    Y_Psi_starT[row_start:row_end, ] <- Y_b %*% t(Psi_star)

    row_start <- row_end + 1L
    rm(Y_b); gc(verbose = FALSE)
  }

  if (row_start - 1L != n)
    stop(sprintf(
      "Batch files contain %d rows total but nrow(X) = %d.",
      row_start - 1L, n))

  # Finalise confounder terms
  ZTZ1ZTY        <- solve(ZTZ, ZTY_acc)   # p2 x V
  env$ZTZ1ZTY    <- ZTZ1ZTY
  env$ZTZ1ZTX    <- solve(ZTZ, crossprod(Z, X))   # p2 x p1
  env$Y_Psi_starT <- Y_Psi_starT

  # Finalise scaling matrix for thetahat update
  if (!is.null(env$Lambda)) {
    env$scaling_matrix_2 <- diag(1 / env$Lambda)
  } else {
    env$scaling_matrix_2 <- solve(tcrossprod(Psi_star))
  }

  # Store XTY for .sira_update_XTY_tilde
  env$XTY <- XTY_acc   # p1 x V

  # Precompute quantities used by the batched initialiser
  # (so .sira_initialize_one_covariate never needs Y directly):
  #
  #   XTY_tilde_0 = t(X) %*% Y_tilde  where Y_tilde = Y - Z (ZTZ)^{-1} ZT Y
  #               = XTY - XTZ %*% ZTZ1ZTY
  #   YtY_tilde   = colSums(Y_tilde^2)
  #               = YtY - colSums(ZTZ1ZTY * ZTY_raw)
  #   (uses the identity: P_Z is a projection, so ||Y_tilde||^2 = ||Y||^2 - ||P_Z Y||^2
  #    and ||P_Z Y||^2_col = colSums(ZTZ1ZTY * ZTY_raw))
  env$XTY_tilde_0 <- XTY_acc - env$XTZ %*% ZTZ1ZTY    # p1 x V
  env$YtY_tilde   <- YtY_acc - colSums(ZTZ1ZTY * ZTY_acc)  # length V

  if (env$verbose)
    message("  Batched summary statistics done.")
}
