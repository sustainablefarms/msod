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
#' @param conditionalLV If TRUE returned probabilities are conditioned on estimated latent variable values (and species are independent due to model structure)
#' If FALSE returned probabilities assume no knowledge of the latent variable values and that species are independent.
#' @return A 2 dimensional array. For each species (column) and each model site (row), the expected number of detections.
#' @export
Endetect_modelsite <- function(fit, type = "median", conditionalLV = TRUE){
  if (!fit$summary.available){ fit <- add.summary(fit)}
  
  # Get ModelSite Occupany Predictions
  ModelSite.Occ.Pred <- poccupy_species(fit, type = type, conditionalLV = conditionalLV)
  
  # Get Detection Probabilities Assuming Occupied
  Visits.DetCond.Pred <- pdetect_condoccupied(fit, type = type)
  
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
pdetect_indvisit <- function(fit, type = "median", conditionalLV = TRUE){
  if (!fit$summary.available){ fit <- add.summary(fit)}
  
  # Get ModelSite Occupany Predictions
  ModelSite.Occ.Pred <- poccupy_species(fit, type = type, conditionalLV = conditionalLV)
  
  # Get Detection Probabilities Assuming Occupied
  Visits.DetCond.Pred <- pdetect_condoccupied(fit, type = type)
 
  # combine with probability of occupancy 
  fit$data <- as_list_format(fit$data)
  if ("ObservedSite" %in% names(fit$data)){ModelSite <- fit$data$ObservedSite} #to enable calculation on the early fitted objects with different name
  if ("ModelSite" %in% names(fit$data)){ModelSite <- fit$data$ModelSite}
  
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
pdetect_condoccupied <- function(fit, type = "median"){
  if (!fit$summary.available){ fit <- add.summary(fit)}
  fit$data <- as_list_format(fit$data)

  if (is.character(type) && type == "marginal"){
    draws <- do.call(rbind, fit$mcmc)
  } else {
    theta <- get_theta(fit, type)
    draws <- matrix(theta, nrow = 1)
    colnames(draws) <- names(theta)
  }
  
  ## v.b (detection coefficients)
  v.b_arr <- bugsvar2array(draws, "v.b", 1:fit$data$n, 1:fit$data$Vobs) # rows are species, columns are observation (detection) covariates
  
  Detection.Pred <- pdetection_occupied_raw(fit$data$Xobs, v.b_arr)
  
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
#' If type is `marginal` then gives probability of each species seperately (marginal over other species?) using full posterior draw.
#' @param Xocc A matrix of occupancy coefficient, with each row corresponding to a ModelSite (i.e. a spatial location and year).
#'  If \code{NULL} the Xocc data saved in \code{fit} will be used.
#' @return A matrix of occupany probabilities. Each row is a ModelSite, corresponding to the rows in Xocc. Each column is a species.
#' @export
poccupy_species <- function(fit, type = "median", conditionalLV = TRUE){
  if (is.character(type) && type == "marginal"){
    draws <- do.call(rbind, fit$mcmc)
  } else {
    theta <- get_theta(fit, type)
    draws <- matrix(theta, nrow = 1)
    colnames(draws) <- names(theta)
  }
  
  Xocc <- fit$data$Xocc
  nlv <- fit$data$nlv
  nspecies <- fit$data$n
  u.b_arr <- bugsvar2array(draws, "u.b", 1:nspecies, 1:ncol(Xocc))
  
  if ((!is.null(nlv) && (nlv > 0))){
    if (conditionalLV){
      lv.coef_arr <- bugsvar2array(draws, "lv.coef", 1:nspecies, 1:nlv)
      LVvals <- bugsvar2array(draws, "LV", 1:nrow(Xocc), 1:nlv)
    } else {#simulate LV values!
      lv.coef_arr <- bugsvar2array(draws, "lv.coef", 1:nspecies, 1:nlv)
      LVvals <- matrix(rnorm(nlv * nrow(draws)),
                       nrow = nrow(draws),
                       ncol = nlv)
    }
  } else { #model has no LV
    stopifnot(!conditionalLV)
    lv.coef_arr <- NULL
    LVvals <- NULL
  }
  
  pocc <- poccupy_species_raw(fit$data$Xocc, u.b_arr, lv.coef_arr, LVvals)
  colnames(pocc) <- fit$species
  return(pocc)
}

poccupy_new_data <- function(fit, Xocc){
  Xocc <- apply.designmatprocess(fit$XoccProcess, Xocc)
  if ((!is.null(nlv) && (nlv > 0))){
    #simulate LV values!
    lv.coef_arr <- bugsvar2array(draws, "lv.coef", 1:nspecies, 1:nlv)
    LVvals <- matrix(rnorm(nlv * nrow(draws)),
                     nrow = nrow(draws),
                     ncol = nlv)
  } else { #model has no LV
    stopifnot(!conditionalLV)
    lv.coef_arr <- NULL
    LVvals <- NULL
  }
  
  pocc <- poccupy_species_raw(Xocc, u.b_arr, lv.coef_arr, LVvals)
  colnames(pocc) <- fit$species
  return(pocc)
}

poccupy_species_raw <- function(Xocc, u.b_arr, lv.coef_arr = NULL, LVvals = NULL){
  if (length(dim(LVvals)) == 3){
    stopifnot(dim(lv.coef_arr)[[3]] == dim(LVvals)[[3]])
    pocc_l <- lapply(1:nrow(Xocc), function(siteid){
      poccupy.ModelSite(Xocc[siteid, , drop = FALSE], 
                        u.b_arr, 
                        lv.coef_arr, 
                        LVvals[siteid, , , drop = FALSE])
    })
  } else {
    pocc_l <- lapply(1:nrow(Xocc), function(siteid){
      poccupy.ModelSite(Xocc[siteid, , drop = FALSE], 
                        u.b_arr, 
                        lv.coef_arr, 
                        LVvals)
    })
  }
  pocc <- do.call(rbind, pocc_l)
  
  return(pocc)
}