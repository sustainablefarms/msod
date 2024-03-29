% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/speciesrichness.R
\name{detectednumspec}
\alias{detectednumspec}
\title{The number of observed species in a matrix of observation recordings}
\usage{
detectednumspec(y, ModelSite)
}
\arguments{
\item{y}{A matrix of species \emph{observations} with each row a visit and each column a species. Entries must be either 0 or 1.}

\item{ModelSite}{The list of ModelSite indexes corresponding to each row in y}
}
\value{
A vector of the number of species detected at each ModelSite. Names give the ModelSite index.
}
\description{
The number of observed species in a matrix of observation recordings
}
