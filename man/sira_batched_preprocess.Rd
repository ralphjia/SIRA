% Generated manually
\name{sira_batched_preprocess}
\alias{sira_batched_preprocess}
\title{Precompute Reusable Batched Summaries for SIRA}
\usage{
sira_batched_preprocess(
  Y_files,
  X,
  Z,
  d1 = NULL,
  d2 = NULL,
  d3 = NULL,
  coords = NULL,
  Y_reader = readRDS,
  Psi_star = NULL,
  delta = 0.01,
  verbose = TRUE
)
}
\arguments{
\item{Y_files}{Character vector of file paths. Each file must be readable by \code{Y_reader} and return a numeric matrix with \eqn{V} columns.}
\item{X}{An \eqn{n \times p_1} numeric matrix of covariates of interest.}
\item{Z}{An \eqn{n \times p_2} numeric matrix of confounders.}
\item{d1}{Integer. Number of voxels along image axis 1.}
\item{d2}{Integer. Number of voxels along image axis 2.}
\item{d3}{Integer. Number of voxels along image axis 3. Use \code{1} for 2D images.}
\item{coords}{Optional \eqn{V \times 3} numeric matrix of voxel coordinates for irregular or masked layouts. Supply \code{coords} instead of \code{d1}, \code{d2}, \code{d3}.}
\item{Y_reader}{Function with signature \code{function(path) -> matrix}. Default is \code{readRDS}.}
\item{Psi_star}{Optional \eqn{L \times V} numeric matrix. Spatial individual effects basis. If \code{NULL}, SIRA computes it internally.}
\item{delta}{Numeric. Convergence threshold carried forward to fitting helpers.}
\item{verbose}{Logical. Print progress messages.}
}
\value{
An object of class \code{"sira_batched_preprocessed"} containing reusable
summary statistics and metadata for \code{\link{sira_batched_fit}}.
}
\description{
Reads the outcome matrix \eqn{Y} from disk in batches and computes the summary
statistics required to run SIRA repeatedly with different tuning parameters
without rereading the image data.
}
