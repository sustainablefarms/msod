#' @title Standalone Occupancy Probabilities and Biodiversity Across Sites
#' @description Operates for models without latent variables. 
#' Given occupancy parameters, computes the probability of each species being detected at each site and at *any* site.
#' @details From Xocc a matrix of the model's cannonical covariates will be computed and then centred and scaled according to [XoccProcess].
#' The chance of occupancy of each species (the model without latent variables assumes species are independent) is then computed using the supplied parameter matrix u.b
#' @examples
#' indata <- readRDS("./private/data/clean/7_2_10_input_data.rds")
#' Xocc <- indata$holdoutdata$Xocc[1:5, ]
#' fit <- readRDS("./tmpdata/7_3_02_clim_someclimate_year_woody500m_msnm_det1stO.rds")
#' XoccProcess <- fit$XoccProcess
#' u.b <- bugsvar2matrix(get_theta(fit, type = "median"), "u.b", rowidx = 1:fit$data$n, colidx = 1:ncol(fit$data$Xocc))
#' rownames(u.b) <- fit$species
#' 
#' multisiterichness_nolv(Xocc, XoccProcess, u.b)
#' poccupancy <- poccupancy_standalone_nolv(Xocc, XoccProcess, u.b)
#' pocc_any <- poccupancy_indsites_nolv(poccupancy)
#' Erichness <- sum(pocc_any) #number of species
#' Vrichness <- sum(pocc_any * (1 - pocc_any))
#' 

#' @describeIn predict_standalone_nolv Predicts occupancy of species at each site for a model without latent variables.
#' @return A matrix. Each row is a row of Xocc (a model site), each column is a species. Values are the probabilty a species occupies the model site.
#' @param Xocc a matrix of covariates in same scale as inputs to the fitted model *before* scaling and centering. Each row is a ModelSite
#' @export
poccupancy_standalone_nolv <- function(Xocc, XoccProcess, u.b){
  XoccStd <- apply.designmatprocess(XoccProcess, Xocc)

  # occupancy probability
  ModelSite.Occ.eta <- as.matrix(XoccStd) %*% t(u.b)
  ModelSite.Occ.Pred <- 1 - pnorm(-ModelSite.Occ.eta, mean = 0, sd = 1)
  return(ModelSite.Occ.Pred)
}

#' @describeIn predict_standalone_nolv The probability of occupancy in any of the ModelSites, where ModelSites are treated as independent.
#' @param poccupancy is an output of poccupancy_standalone_nolv. Each row is a ModelSite, each column is a species.
#' Values are the probability of a species occupying a ModelSite
#' @export
panyoccupancy_indsites_nolv <- function(poccupancy){
  anyoccupancy <- 1 - Rfast::colprods(1 - poccupancy)
  names(anyoccupancy) <- colnames(poccupancy)
  return(anyoccupancy)
}

#' @describeIn predict_standalone_nolv For species richness across multiple independent ModelSites 
#' @export
multisiterichness_nolv <- function(Xocc, XoccProcess, u.b){
  poccupancy <- poccupancy_standalone_nolv(Xocc, XoccProcess, u.b)
  pocc_any <- panyoccupancy_indsites_nolv(poccupancy)
  Erichness <- sum(pocc_any) #number of species
  Vrichness <- sum(pocc_any * (1 - pocc_any))
  return(c(Erichness = Erichness,
           Vrichness = Vrichness))
}
