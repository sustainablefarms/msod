% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction_summaries_jsodm_lv.R
\name{poccupancy_margotherspecies.jsodm_lv}
\alias{poccupancy_margotherspecies.jsodm_lv}
\title{Prediction Summaries for jsodm_lv Model for New Sites}
\usage{
poccupancy_margotherspecies.jsodm_lv(fit)
}
\arguments{
\item{Xocc}{A matrix of unprocessed occupancy covariates. Each row is a site, each column a covariate.}

\item{occ.b}{An array of occupancy covariate loadings. Each row is a species, each column a covariate, and each layer a draw from the posterior.}
}
\description{
Interesting summaries of predictions made by a fitted jsodm_lv.
Includes
occupancy probability of species independent of all other species, for each site, with error bars;
expected species richness for single site;
}
\examples{
fit <- readRDS("../sflddata/private/data/testdata/cutfit_7_4_11_2LV.rds")
fit <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_model_2lv_e13.rds")
fit <- translatefit(fit)
Xocc <- sflddata::unstandardise.designmatprocess(fit$XoccProcess, fit$data$Xocc[1:5, , drop = FALSE])
fittmp <- supplant_new_data(fit, Xocc)
pocc <- poccupancy_margotherspecies.jsodm_lv(fittmp)
pocc <- poccupancy_margotherspeciespmaxsite.jsodm_lv(fittmp)
pocc <- poccupancy_mostfavourablesite.jsodm_lv(fittmp)
pocc <- poccupancy_randomsite.jsodm_lv(fittmp)
sprich1 <- occspecrichness.jsodm_lv(fittmp)
sprich <- occspecrichnessRV.jsodm_lv(fittmp)
system.time(sprich <- occspecrichness_avsite.jsodm_lv(fittmp))
}
