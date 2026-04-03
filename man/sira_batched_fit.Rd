% Generated manually
\name{sira_batched_fit}
\alias{sira_batched_fit}
\title{Fit SIRA from Precomputed Batched Summaries}
\usage{
sira_batched_fit(
  preprocessed,
  lambda,
  mu,
  max_iter = 50L,
  delta = NULL,
  verbose = NULL
)
}
\arguments{
\item{preprocessed}{A \code{"sira_batched_preprocessed"} object created by \code{\link{sira_batched_preprocess}}.}
\item{lambda}{Non-negative numeric. Spatial penalty controlling boundary smoothness.}
\item{mu}{Non-negative numeric. L1 sparsity penalty.}
\item{max_iter}{Integer. Maximum outer coordinate-descent iterations.}
\item{delta}{Optional numeric. Overrides the stored convergence threshold.}
\item{verbose}{Optional logical. Overrides the stored verbosity flag.}
}
\value{
An object of class \code{"sira"}.
}
\description{
Runs SIRA using a precomputed batched summary object. This avoids rereading
the outcome data when fitting multiple tuning parameter combinations.
}
