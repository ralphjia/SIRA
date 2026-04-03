make_toy_sira_data <- function() {
  set.seed(123)

  d1 <- 3L
  d2 <- 3L
  d3 <- 2L
  V  <- d1 * d2 * d3
  n  <- 24L

  X <- cbind(rnorm(n))
  Z <- cbind(1, rnorm(n))

  beta <- numeric(V)
  beta[c(1L, 2L, 4L)] <- 1.8
  beta[c(15L, 18L)]   <- -1.4

  gamma <- rbind(rep(0.5, V), rep(-0.2, V))
  Y <- X %*% t(beta) +
    Z %*% gamma +
    matrix(rnorm(n * V, sd = 0.35), nrow = n, ncol = V)

  list(
    Y = Y,
    X = X,
    Z = Z,
    d1 = d1,
    d2 = d2,
    d3 = d3,
    V = V,
    Psi_star = matrix(1, nrow = 1L, ncol = V)
  )
}

test_that("sira returns the expected fit structure", {
  dat <- make_toy_sira_data()

  fit <- sira(
    Y = dat$Y,
    X = dat$X,
    Z = dat$Z,
    d1 = dat$d1,
    d2 = dat$d2,
    d3 = dat$d3,
    lambda = 0.05,
    mu = 0.01,
    Psi_star = dat$Psi_star,
    max_iter = 2L,
    verbose = FALSE
  )

  expect_s3_class(fit, "sira")
  expect_equal(dim(fit$betahat), c(ncol(dat$X), dat$V))
  expect_equal(dim(fit$gammahat), c(ncol(dat$Z), dat$V))
  expect_true(is.list(fit$region_list))
  expect_length(fit$region_list, ncol(dat$X))
  expect_true(is.data.frame(fit$convergence))
  expect_true(all(c("iteration", "num_regions", "num_voxels", "alpha_norm_diff") %in%
                    names(fit$convergence)))
})

test_that("initializer produces at most two signed starting regions", {
  dat <- make_toy_sira_data()

  env <- new.env(parent = emptyenv())
  env$n <- nrow(dat$Y)
  env$V <- dat$V
  env$p1 <- ncol(dat$X)
  env$p2 <- ncol(dat$Z)
  env$d1 <- dat$d1
  env$d2 <- dat$d2
  env$d3 <- dat$d3
  env$coords <- NULL
  env$X <- dat$X
  env$Y <- dat$Y
  env$Z <- dat$Z
  env$Psi_star <- dat$Psi_star
  env$Lambda <- 1
  env$verbose <- FALSE

  SIRA:::.sira_build_neighbors(env)
  SIRA:::.sira_setup_summary_stats(env)
  rl <- SIRA:::.sira_initialize_one_covariate(1L, env)

  expect_true(length(rl) >= 1L)
  expect_true(length(rl) <= 2L)
  expect_true(all(vapply(rl, function(r) length(r[[1]]) > 0L, logical(1L))))
  expect_true(all(vapply(rl, function(r) is.finite(r[[2]]), logical(1L))))
})

test_that("batched fitting matches in-memory fitting on a toy example", {
  dat <- make_toy_sira_data()

  fit_mem <- sira(
    Y = dat$Y,
    X = dat$X,
    Z = dat$Z,
    d1 = dat$d1,
    d2 = dat$d2,
    d3 = dat$d3,
    lambda = 0.05,
    mu = 0.01,
    Psi_star = dat$Psi_star,
    max_iter = 2L,
    verbose = FALSE
  )

  idx <- split(seq_len(nrow(dat$Y)), c(rep(1L, 12L), rep(2L, 12L)))
  files <- file.path(tempdir(), paste0("sira-batch-", seq_along(idx), ".rds"))
  on.exit(unlink(files), add = TRUE)
  for (i in seq_along(idx)) {
    saveRDS(dat$Y[idx[[i]], , drop = FALSE], files[i])
  }

  fit_batch <- sira_batched(
    Y_files = files,
    X = dat$X,
    Z = dat$Z,
    d1 = dat$d1,
    d2 = dat$d2,
    d3 = dat$d3,
    lambda = 0.05,
    mu = 0.01,
    Psi_star = dat$Psi_star,
    max_iter = 2L,
    verbose = FALSE
  )

  expect_equal(fit_batch$betahat, fit_mem$betahat, tolerance = 1e-8)
  expect_equal(fit_batch$gammahat, fit_mem$gammahat, tolerance = 1e-8)
  expect_equal(lengths(fit_batch$region_list), lengths(fit_mem$region_list))
})

test_that("batched preprocessing can be reused across tuning parameters", {
  dat <- make_toy_sira_data()

  idx <- split(seq_len(nrow(dat$Y)), c(rep(1L, 12L), rep(2L, 12L)))
  files <- file.path(tempdir(), paste0("sira-prep-batch-", seq_along(idx), ".rds"))
  on.exit(unlink(files), add = TRUE)
  for (i in seq_along(idx)) {
    saveRDS(dat$Y[idx[[i]], , drop = FALSE], files[i])
  }

  read_counter <- 0L
  counting_reader <- function(path) {
    read_counter <<- read_counter + 1L
    readRDS(path)
  }

  prep <- sira_batched_preprocess(
    Y_files = files,
    X = dat$X,
    Z = dat$Z,
    d1 = dat$d1,
    d2 = dat$d2,
    d3 = dat$d3,
    Y_reader = counting_reader,
    Psi_star = dat$Psi_star,
    verbose = FALSE
  )

  expect_s3_class(prep, "sira_batched_preprocessed")
  expect_equal(read_counter, length(files))

  fit1 <- sira_batched_fit(
    prep, lambda = 0.05, mu = 0.01,
    max_iter = 2L, verbose = FALSE
  )
  fit2 <- sira_batched_fit(
    prep, lambda = 0.3, mu = 0.1,
    max_iter = 2L, verbose = FALSE
  )

  expect_equal(read_counter, length(files))
  expect_s3_class(fit1, "sira")
  expect_s3_class(fit2, "sira")
  expect_equal(dim(fit1$betahat), dim(fit2$betahat))
})
