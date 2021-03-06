% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcpredictions.R
\name{Endetect_modelsite}
\alias{Endetect_modelsite}
\alias{pdetect_indvisit}
\alias{pdetect_condoccupied}
\alias{poccupy_species}
\title{Predicted Probabilities
** would be great to test these functions by monitoring z, p and mu.p in runjags (perhaps using hidden.monitor parameter)}
\usage{
Endetect_modelsite(fit, type = "median", conditionalLV = TRUE)

pdetect_indvisit(fit, type = "median", conditionalLV = TRUE)

pdetect_condoccupied(fit, type = "median")

poccupy_species(fit, type = "median", conditionalLV = TRUE, numLVsims = 10000)
}
\arguments{
\item{fit}{is a fitted runjags model}

\item{type}{is the type of point estimate to use. See get_theta().
If type is \code{marginal} then gives probability of each species seperately (marginal over other species?) using full posterior draw.}

\item{conditionalLV}{If TRUE returned probabilities are conditioned on estimated latent variable values (and species are independent due to model structure)
If FALSE returned probabilities assume no knowledge of the latent variable values and that species are independent.}

\item{Xobs}{A matrix of observation (detection) coefficients. Default is the observation coefficients saved in \code{fit}}

\item{ModelSite}{A list mapping each row in \code{Xobs} to the row in \code{Xocc} that represents the ModelSite visited.}

\item{Xocc}{A matrix of occupancy coefficient, with each row corresponding to a ModelSite (i.e. a spatial location and year).
If \code{NULL} the Xocc data saved in \code{fit} will be used.}
}
\value{
A 2 dimensional array. For each species (column) and each model site (row), the expected number of detections.

A matrix of detection probabilities. Each row is a visit, corresponding to the rows in Xobs. Each column is a species.

A matrix of detection probabilities. Each row is a visit, corresponding to the rows in Xobs. Each column is a species.

A matrix of occupany probabilities. Each row is a ModelSite, corresponding to the rows in Xocc. Each column is a species.
}
\description{
Predicted Probabilities
** would be great to test these functions by monitoring z, p and mu.p in runjags (perhaps using hidden.monitor parameter)
}
\section{Functions}{
\itemize{
\item \code{Endetect_modelsite}: For a point estimate of model parameters,
for each species and each ModelSite,
computes the expected number of detections of each species \emph{independent} of other species.

\item \code{pdetect_indvisit}: For a point estimate of model parameters,
for each species and each visit,
computes the probability of detection ignoring results from other visits (and potentially optionally species).

\item \code{pdetect_condoccupied}: For a point estimate of model parameters, the probability of detecting a species for visits,
conditional on the species occupying the ModelSite.

\item \code{poccupy_species}: For a point estimate of model parameters,
for each species,
computes the probability of the species occupying ModelSites.
@param condtionalLV Whether to condition on estimated LV values, or ignore them.
If TRUE returns probabilities of the species occupying ModelSites
given estimated latent variable values for each site.
If FALSE returns the probabilities of species occupying ModelSites
ignoring other species and assuming no knowledge of the latent variable values.
}}

\examples{
fit <- readRDS("./tmpdata/deto_wind.rds")
pDetection <- pdetect_indvisit(fit, type = "median", conditionalLV = FALSE)
pOccupancy <- poccupy_species(fit, type = "median", conditionalLV = FALSE)
}
