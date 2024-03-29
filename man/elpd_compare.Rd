% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elpd_compare.R
\name{elpd_compare}
\alias{elpd_compare}
\title{Compare Sitewise Differences in Log Predictive Density; Expected Difference with Standard Errors}
\usage{
elpd_compare(elpds, refname = NULL)
}
\arguments{
\item{elpds}{A data.frame or matrix of expected log posterior densities. Each row is ModelSite, each column is a Model. Typically the columns will have names.}

\item{refname}{The column name (i.e. Model) to compare all other models to. If NULL then uses the model with the highest sum of elpd.}
}
\value{
A matrix with a column for the aggregate (summed) difference of model sites between models,
and the standard error of this difference (computed as the sample standard deviation of difference, multiplied by the square root of the number of ModelSites)
}
\description{
Compare Sitewise Differences in Log Predictive Density; Expected Difference with Standard Errors
}
\details{
Very similar to loo::loo_compare(), but with fewer checks
}
\examples{
waics_l <- readRDS("./tmpdata/7_3_01_waics.rds")
mat1 <- elpd_compare(elpds = do.call(cbind, lapply(waics_l, function(x) x$loo$pointwise[, 1])))
matref <- loo::loo_compare(lapply(waics_l, function(x) x$loo))
all(mat1 == matref[, c("elpd_diff", "se_diff")])
}
