# SIRA: Scalable Image-on-Scalar Regression Algorithm

SIRA fits voxelwise image-on-scalar regression models where the coefficient
images are **sparse** (most voxels have no effect) and **spatially
homogeneous** (nonzero voxels form contiguous regions with a constant value).
It is designed to scale to high-dimensional neuroimaging data.

---

## Model

For subject $i$ at voxel $v$:

$$Y_i(s_v) = \sum_{j=1}^{p_1} \beta_j(s_v) \, X_{ij} + \sum_{k=1}^{p_2} \gamma_k(s_v) \, Z_{ik} + \eta_i(s_v) + \epsilon_i(s_v)$$

- **$X$** — covariates of interest (e.g. age, diagnosis, cognitive score)
- **$Z$** — confounders, typically including an intercept
- **$\beta_j$** — spatially sparse, piecewise-constant coefficient images (the primary output)
- **$\eta_i$** — subject-level spatial random effect, modelled via a GP basis
- **$\epsilon_i$** — independent noise

Sparsity and spatial homogeneity are jointly enforced via an L1 penalty
($\mu$) and a sum-of-absolute-differences boundary penalty ($\lambda$).

---

## Installation

```r
# Install dependencies first
install.packages(c("Matrix", "locfdr", "pracma", "BayesGPfit", "Rcpp"))

# Install SIRA from source
devtools::install_github("ralphjia/SIRA")
```

---

## Quick start

```r
library(SIRA)

set.seed(1)
d1 <- d2 <- d3 <- 10L
V  <- d1 * d2 * d3
n  <- 50L

# True signal: a 3x3x3 cube of voxels near the origin with beta = 3
beta_true        <- numeric(V)
beta_true[1:27]  <- 3

X <- matrix(rnorm(n), ncol = 1)
Z <- cbind(1, rnorm(n))
Y <- X %*% t(beta_true) +
     Z %*% matrix(rnorm(2 * V), nrow = 2) +
     matrix(rnorm(n * V, sd = 3), nrow = n)

# Y : n x V matrix of subject images (vectorised, column-major)
# X : n x p1 matrix of covariates of interest
# Z : n x p2 matrix of confounders (include a column of 1s for the intercept)
fit <- sira(Y = Y, X = X, Z = Z,
            d1 = d1, d2 = d2, d3 = d3,
            lambda = 0.05, mu = 0.01)

print(fit)
```

The two key tuning parameters are:

| Parameter | Controls | Larger value → |
|-----------|----------|----------------|
| `lambda`  | Spatial smoothness | Larger, blockier regions |
| `mu`      | Sparsity | Fewer nonzero voxels |

---

## Irregular and brain-masked data

For brain-masked or irregularly spaced data, pass a V x 3 matrix of voxel
coordinates instead of `d1/d2/d3`:

```r
# coords : V x 3 numeric matrix of (x, y, z) spatial coordinates
fit <- sira(Y = Y, X = X, Z = Z,
            coords = coords,
            lambda = 0.05, mu = 0.01)
```

SIRA infers adjacency from the coordinate spacing — two voxels are neighbours
if they are one voxel-width apart along exactly one axis.

---

## Working with results

```r
# p1 x V matrix of coefficient image estimates
beta_mat <- coef(fit, type = "beta")

# Reshape covariate j=1 back into a 3D array for visualisation
beta_array <- as.array(fit, j = 1)   # dim: c(d1, d2, d3)

# Region structure: voxel indices and beta value for each region
fit$region_list[[1]]

# Convergence history
fit$convergence
```

---

## Large datasets

When the full $n \times V$ matrix cannot be held in memory, use
`sira_batched()`. Instead of passing `Y` directly, provide a character vector
of file paths — each file is read once, used to accumulate summary statistics,
and then discarded. An example with reading NIfTI batch files is given in the 
vignette. 

```r
# Y is split into batch files, e.g. ~700 subjects per file
fit <- sira_batched(
  Y_files  = list.files("path/to/batches", pattern = "\\.rds$",
                        full.names = TRUE),
  Y_reader = readRDS,   # or a custom reader
  X = X, Z = Z,
  d1 = d1, d2 = d2, d3 = d3,
  lambda = 0.05, mu = 0.01
)
```

The return value is identical to `sira()` — all downstream code works
unchanged.

---

## Vignette

A full worked example with two covariates, spatially separated positive and
negative regions, and diagnostic output is available in the package vignette:

```r
vignette("sira-introduction", package = "SIRA")
```

---

## Citation

A manuscript describing the SIRA methodology is in preparation. In the
meantime, please cite this software repository directly:

> Ralph Jiang (2026). *SIRA: Scalable Image-on-Scalar Regression Algorithm*.
> R package. https://github.com/ralphjia/SIRA

---

## License

MIT
