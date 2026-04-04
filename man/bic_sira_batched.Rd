\name{bic_sira_batched}
\alias{bic_sira_batched}
\title{Compute the Boundary-Voxel BIC for a Batched SIRA Fit}
\usage{
bic_sira_batched(object, Y_files, X, Z, Y_reader = readRDS, ...)
}
\arguments{
\item{object}{A \code{"sira"} object, typically returned by
\code{\link{sira_batched_fit}} or \code{\link{sira_batched}}.}
\item{Y_files}{Character vector of file paths. Each file must be readable by
\code{Y_reader} and return a numeric matrix with \eqn{V} columns.}
\item{X}{An \eqn{n \times p_1} numeric matrix of covariates of interest.}
\item{Z}{An \eqn{n \times p_2} numeric matrix of confounders.}
\item{Y_reader}{Function with signature \code{function(path) -> matrix}.
Default is \code{readRDS}.}
\item{...}{Ignored.}
}
\value{
A named numeric vector containing \code{BIC}, \code{MSE}, and
\code{boundary_voxels}.
}
\description{
Computes the Bayesian Information Criterion (BIC) for a fitted
\code{"sira"} object using the degrees-of-freedom definition from the SIRA
paper, where the effective model size is the total number of boundary voxels
across the covariates of interest. This variant reads the outcome data from
disk in batches and does not require the full \eqn{n \times V} outcome
matrix to be held in memory.
}
\details{
The criterion is
\deqn{\mathrm{BIC} = \log(nV)\,\mathrm{df} + nV\log(\mathrm{MSE}),}
where \eqn{\mathrm{df}} is the number of boundary voxels and
\eqn{\mathrm{MSE}} is computed from the fitted mean
\eqn{X\hat\beta + Z\hat\gamma + \hat\theta\Psi^\ast}.
}
