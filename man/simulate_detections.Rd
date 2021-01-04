% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_fit.R
\name{simulate_detections}
\alias{simulate_detections}
\alias{simulate_detections_LV}
\title{Simulate observations from parameters of a fitted object.}
\usage{
simulate_detections(fit, esttype = "median")

simulate_detections_LV(fit, esttype = "median")
}
\arguments{
\item{fit}{A runjags fitted object created by \code{\link[=run.detectionoccupancy]{run.detectionoccupancy()}}}

\item{esttype}{Specifies parameter set to extract from fit, see \code{\link[=get_theta]{get_theta()}}}

\item{UseFittedLV}{Logical. If TRUE, the simulation uses fitted latent variable values.
If FALSE, latent variable values will be simulated for each ModelSite}
}
\description{
Simulate observations from parameters of a fitted object.
}
\details{
Simulate detections from all given parameters (including LV and random effects)
}
\section{Functions}{
\itemize{
\item \code{simulate_detections_LV}: Simulate LV values and use these for simulating detections (with \code{\link[=simulate_detections]{simulate_detections()}}). For each site the simulated LV values are copied across all draws.
}}

\examples{
fit <- readRDS("./tmpdata/deto_wind.rds")
fit$data <- as_list_format(fit$data)
class(fit) <- c("jsodm", class(fit))
detected <- simulate_detections(fit, esttype = "median", UseFittedLV = TRUE)

}