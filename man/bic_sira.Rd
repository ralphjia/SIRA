\name{bic_sira}
\alias{bic_sira}
\title{Compute the Boundary-Voxel BIC for a SIRA Fit}
\usage{
bic_sira(object, Y, X, Z, ...)
}
\arguments{
\item{object}{A \code{"sira"} object.}
\item{Y}{An \eqn{n \times V} numeric matrix of image outcomes.}
\item{X}{An \eqn{n \times p_1} numeric matrix of covariates of interest.}
\item{Z}{An \eqn{n \times p_2} numeric matrix of confounders.}
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
across the covariates of interest.
}
\details{
The criterion is
\deqn{\mathrm{BIC} = \log(nV)\,\mathrm{df} + nV\log(\mathrm{MSE}),}
where \eqn{\mathrm{df}} is the number of boundary voxels and
\eqn{\mathrm{MSE}} is computed from the fitted mean
\eqn{X\hat\beta + Z\hat\gamma + \hat\theta\Psi^\ast}.
}
