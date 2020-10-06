# library(runjags); library(dplyr); library(tidyr); library(tibble);

#' @title Predicted Probabilities
#' ** would be great to test these functions by monitoring z, p and mu.p in runjags (perhaps using hidden.monitor parameter)

#' @examples 
#' fit <- readRDS("./tmpdata/deto_wind.rds")
#' pDetection <- pdetect_indvisit(fit, type = "median", conditionalLV = FALSE)
#' pOccupancy <- poccupy_species(fit, type = "median", conditionalLV = FALSE)


#' @describeIn predictedprobabilities  
#'  For a point estimate of model parameters,
#'  for each species and each ModelSite,
#'  computes the expected number of detections of each species *independent* of other species.
#' @param fit is a fitted runjags model
#' @param type is the type of point estimate to use. See get_theta for supported options.
#' @param Xocc A matrix of occupancy coefficient, with each row corresponding to a ModelSite (i.e. a spatial location and year).
#'  If \code{NULL} the Xocc data saved in \code{fit} will be used.
#' @param Xobs A matrix of observation (detection) coefficients. Default is the observation coefficients saved in \code{fit}
#' @param ModelSite A list mapping each row in \code{Xobs} to the row in \code{Xocc} that represents the ModelSite visited.
#' @param conditionalLV If TRUE returned probabilities are conditioned on estimated latent variable values (and species are independent due to model structure)
#' If FALSE returned probabilities assume no knowledge of the latent variable values and that species are independent.
#' @return A 2 dimensional array. For each species (column) and each model site (row), the expected number of detections.
#' @export
Endetect_modelsite <- function(fit, type = "median", Xocc = NULL, Xobs = NULL, ModelSite = NULL, conditionalLV = TRUE){
  if (!fit$summary.available){ fit <- add.summary(fit)}
  
  # Get ModelSite Occupany Predictions
  ModelSite.Occ.Pred <- poccupy_species(fit, type = type, Xocc = Xocc, conditionalLV = conditionalLV)
  
  # Get Detection Probabilities Assuming Occupied
  Visits.DetCond.Pred <- pdetect_condoccupied(fit, type = type, Xobs = Xobs)
  
  # combine with probability of occupancy 
  fitdata <- as_list_format(fit$data)
  if (is.null(ModelSite)){
    if ("ObservedSite" %in% names(fitdata)){ModelSite <- fitdata$ObservedSite} #to enable calculation on the early fitted objects with different name
    if ("ModelSite" %in% names(fitdata)){ModelSite <- fitdata$ModelSite}
  }
  
  
  
  # Expected num detections for each ModelSite conditional on occupied
  E_ndetect_condocc <- cbind(ModelSite = ModelSite, Visits.DetCond.Pred) %>%
    as_tibble() %>%
    group_by(ModelSite) %>%
    summarise_all(sum) %>%
    dplyr::select(-ModelSite)
  V_ndetect_condocc <- cbind(ModelSite = ModelSite, Visits.DetCond.Pred * (1 - Visits.DetCond.Pred)) %>%
    as_tibble() %>%
    group_by(ModelSite) %>%
    summarise_all(sum) %>%
    dplyr::select(-ModelSite)
  E_ndetect2_condocc <- V_ndetect_condocc + E_ndetect_condocc^2
  # Convert to marginal expected num detection, considering that no occupancy --> no detections
  E_ndetect <- as.matrix(E_ndetect_condocc) * ModelSite.Occ.Pred
  E_ndetect2 <- as.matrix(E_ndetect2_condocc) * ModelSite.Occ.Pred
  V_ndetect <- E_ndetect2 - E_ndetect^2
  if (!is.null(fit$species)){  # a special modification of runjags with occupation detection meta info
    colnames(E_ndetect) <- fit$species
    colnames(V_ndetect) <- fit$species}
  return(list(
    E_ndetect = E_ndetect,
    V_ndetect = V_ndetect
  ))
}



# Getting to the probability of detecting a species at a particular site
# marginal and otherwise too
#' @describeIn predictedprobabilities  
#'  For a point estimate of model parameters,
#'  for each species and each visit,
#'  computes the probability of detection ignoring results from other visits (and potentially optionally species).
#' @param fit is a fitted runjags model
#' @param type is the type of point estimate to use. See get_theta for supported options.
#' @param Xocc A matrix of occupancy coefficient, with each row corresponding to a ModelSite (i.e. a spatial location and year).
#'  If \code{NULL} the Xocc data saved in \code{fit} will be used.
#' @param Xobs A matrix of observation (detection) coefficients. Default is the observation coefficients saved in \code{fit}
#' @param ModelSite A list mapping each row in \code{Xobs} to the row in \code{Xocc} that represents the ModelSite visited.
#' @param conditionalLV If TRUE returned probabilities are conditioned on estimated latent variable values (and species are independent due to model structure)
#' If FALSE returned probabilities assume no knowledge of the latent variable values and that species are independent.
#' @return A matrix of detection probabilities. Each row is a visit, corresponding to the rows in Xobs. Each column is a species.
#' @export
pdetect_indvisit <- function(fit, type = "median", Xocc = NULL, Xobs = NULL, ModelSite = NULL, conditionalLV = TRUE){
  if (!fit$summary.available){ fit <- add.summary(fit)}
  
  # Get ModelSite Occupany Predictions
  ModelSite.Occ.Pred <- poccupy_species(fit, type = type, Xocc = Xocc, conditionalLV = conditionalLV)
  
  # Get Detection Probabilities Assuming Occupied
  Visits.DetCond.Pred <- pdetect_condoccupied(fit, type = type, Xobs = Xobs)
 
  # combine with probability of occupancy 
  fitdata <- as_list_format(fit$data)
  if (is.null(ModelSite)){
    if ("ObservedSite" %in% names(fitdata)){ModelSite <- fitdata$ObservedSite} #to enable calculation on the early fitted objects with different name
    if ("ModelSite" %in% names(fitdata)){ModelSite <- fitdata$ModelSite}
  }
  
  ## Full Probability of Detection (marginalises across occupancy, detection and latent variables)
  Detection.Pred.Marg <- ModelSite.Occ.Pred[ModelSite, ] * Visits.DetCond.Pred
  
  if (!is.null(fit$species)){colnames(Detection.Pred.Marg) <- fit$species} # a special modification of runjags with occupation detection meta info
  return(Detection.Pred.Marg)
}


#' @describeIn predictedprobabilities  For a point estimate of model parameters, the probability of detecting a species for visits,
#' conditional on the species occupying the ModelSite.
#' @param fit is a fitted runjags model
#' @param type is the type of point estimate to use. See get_theta() for available options.
#' @param Xobs A matrix of observation (detection) coefficients. Default is the observation coefficients saved in \code{fit}
#' @param ModelSite A list mapping each row in \code{Xobs} to the row in \code{Xocc} that represents the ModelSite visited.
#' @return A matrix of detection probabilities. Each row is a visit, corresponding to the rows in Xobs. Each column is a species.
#' @export
pdetect_condoccupied <- function(fit, type = "median", Xobs = NULL){
  if (!fit$summary.available){ fit <- add.summary(fit)}
  fitdata <- as_list_format(fit$data)
  
  # build list of point estimates
  theta <- get_theta(fit, type = type)
  
  ## v.b (detection coefficients)
  v.b <- bugsvar2matrix(theta, "v.b", 1:fitdata$n, 1:fitdata$Vobs) # rows are species, columns are observation (detection) covariates
 
  ## Get observation covariates 
  if (is.null(Xobs)){Xobs <- fitdata$Xobs}
  
  ## Probability of Detection, assuming occupied
  Detection.LinPred <- Xobs %*% t(v.b)
  Detection.Pred <- exp(Detection.LinPred) / (exp(Detection.LinPred) + 1)   #this is the inverse logit function
  
  if (!is.null(fit$species)){colnames(Detection.Pred) <- fit$species} # a special modification of runjags with occupation detection meta info
  return(Detection.Pred)
}

#' @describeIn predictedprobabilities  
#'  For a point estimate of model parameters,
#'  for each species,
#'  computes the probability of the species occupying ModelSites.
#'  @param condtionalLV Whether to condition on estimated LV values, or ignore them.
#'  If TRUE returns probabilities of the species occupying ModelSites
#'  given estimated latent variable values for each site.
#'  If FALSE returns the probabilities of species occupying ModelSites
#'  ignoring other species and assuming no knowledge of the latent variable values.
#' @param fit is a fitted runjags model
#' @param type is the type of point estimate to use. See get_theta().
#' @param Xocc A matrix of occupancy coefficient, with each row corresponding to a ModelSite (i.e. a spatial location and year).
#'  If \code{NULL} the Xocc data saved in \code{fit} will be used.
#' @return A matrix of occupany probabilities. Each row is a ModelSite, corresponding to the rows in Xocc. Each column is a species.
#' @export
poccupy_species <- function(fit, type = "median", Xocc = NULL, conditionalLV = TRUE){
  if (!fit$summary.available){ fit <- add.summary(fit)}
  fitdata <- as_list_format(fit$data)
  # build arrays of point estimates
  theta <- get_theta(fit, type = type)
  ## u.b (occupancy coefficients)
  u.b <- bugsvar2matrix(theta, "u.b", 1:fitdata$n, 1:fitdata$Vocc) # rows are species, columns are occupancy covariates
  
  if (conditionalLV){
    stopifnot(any(grepl("LV", names(theta)))) #means LV values not saved in model
    ## LV values
    LV <- bugsvar2matrix(theta, "LV", 1:fitdata$J, 1:fitdata$nlv) # rows are model sites, columns are latent variables
    ## LV loadings
    lv.coef <- bugsvar2matrix(theta, "lv.coef", 1:fitdata$n, 1:fitdata$nlv) # rows are species columns are occupancy covariates
    sd_u_condlv <- sqrt(1 - rowSums(lv.coef^2))
  }
  
  if (is.null(Xocc)){Xocc <- fitdata$Xocc}
  
  ## Probability of Site Occupancy
  ModelSite.Occ.eta <- Xocc %*% t(u.b) #rows are ModelSites, columns are species
  if (conditionalLV){
    ModelSite.Occ.eta <- ModelSite.Occ.eta + LV %*% t(lv.coef)
    # Make u standard deviations equal to 1 by dividing other values by sd
    # P(u < -ModelSite.Occ.eta) = P(u / sd < -ModelSite.Occ.eta / sd) = P(z < -ModelSite.Occ.eta / sd)
    ModelSite.Occ.eta <- Rfast::eachrow(ModelSite.Occ.eta, sd_u_condlv, oper = "/")     
    }
  ModelSite.Occ.Pred.CondLV <- 1 - pnorm(-ModelSite.Occ.eta, mean = 0,
                                         sd = 1)
  
  if (!is.null(fit$species)){colnames(ModelSite.Occ.Pred.CondLV) <- fit$species} # a special modification of runjags with occupation detection meta info
  return(ModelSite.Occ.Pred.CondLV)
}