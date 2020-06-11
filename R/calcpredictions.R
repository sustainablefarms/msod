library(runjags); library(dplyr); library(tidyr); library(tibble);

#' @title Predicted Probabilities
#' ** would be great to test these functions by monitoring z, p and mu.p in runjags (perhaps using hidden.monitor parameter)

#' @examples 
#' fit <- readRDS("./tmpdata/deto_wind.rds")
#' pDetection <- pdetect_indvisit(fit, type = "median", conditionalLV = FALSE)
#' pOccupancy <- poccupy_species(fit, type = "median", conditionalLV = FALSE)

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
#' @value A matrix of detection probabilities. Each row is a visit, corresponding to the rows in Xobs. Each column is a species.
pdetect_indvisit <- function(fit, type = "median", Xocc = NULL, Xobs = NULL, ModelSite = NULL, conditionalLV = TRUE){
  if (!fit$summary.available){ fit <- add.summary(fit)}
  
  # Get ModelSite Occupany Predictions
  ModelSite.Occ.Pred <- poccupy_species(fit, type = type, Xocc = Xocc, conditionalLV = conditionalLV)
  
  # Get Detection Probabilities Assuming Occupied
  Visits.DetCond.Pred <- pdetect_condoccupied(fit, type = type, Xobs = Xobs)
 
  # combine with probability of occupancy 
  fitdata <- as.list.format(fit$data)
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
#' @value A matrix of detection probabilities. Each row is a visit, corresponding to the rows in Xobs. Each column is a species.
pdetect_condoccupied <- function(fit, type = "median", Xobs = NULL){
  if (!fit$summary.available){ fit <- add.summary(fit)}
  fitdata <- as.list.format(fit$data)
  
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
#' @value A matrix of occupany probabilities. Each row is a ModelSite, corresponding to the rows in Xocc. Each column is a species.
poccupy_species <- function(fit, type = "median", Xocc = NULL, conditionalLV = TRUE){
  if (!fit$summary.available){ fit <- add.summary(fit)}
  fitdata <- as.list.format(fit$data)
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

#' @describeIn predictedprobabilities For converting values for array-valued parameters from the bugs variable format to a matrix
#' @param values is a list of values named according to the bugs variables names
#' @param varname is the desired variable name (e.g. 'u.b')
#' @param rowidx is a list of rows to extract, by number
#' @param colidx is a list of columns to extract, by number
bugsvar2array <- function(values, varname, rowidx, colidx){
  if (is.vector(values)) {
    values <- matrix(values, nrow = 1, dimnames = list(row = NULL, col = names(values)))
  }
  idx <- expand.grid(row = rowidx, col = colidx)
  bugsnames <- paste0(varname, "[",idx$row, ",", idx$col, "]") #order matters, expand.grid must go through rows and then columns
  value <- array(t(values[, bugsnames]), 
                 dim = c(length(rowidx), length(colidx), nrow(values)), 
                 dimnames = list(row = rowidx, col = colidx, draw = 1:nrow(values)))
  return(value)
}
bugsvar2matrix <- function(values, varname, rowidx, colidx){
  arr <- bugsvar2array(values, varname, rowidx, colidx)
  mat <- drop_to_matrix(arr, 3)
  return(return(mat))
}

#' @param theta is a parameter arrays
#' @param name parameter
#' @return a named vector. Names are given by name and the index in the array
matrix2bugsvar <- function(theta, name){
  values <- as.vector(theta) #runs through first dimension, then second dimension, then third dimension...
  idx <- expand.grid(row = 1:nrow(theta), col = 1:ncol(theta))
  bugsnames <- paste0(name, "[",idx$row, ",", idx$col, "]") #order matters, expand.grid must go through rows and then columns
  names(values) <- bugsnames
  return(values)
}

#' @describeIn A helper function to get a vector of parameters from a fitted object
#' @param fit fitted runjags object with summary included
#' @param type An integer will select the draw from the posterior
#' (e.g. type = 5 will return the vector of parameters in the 5th sample from the posterior, mcmc chains are concatenated)
#' A charactor vector can be used to the parameters based on summarary statistics like quantiles and moments.
#' Supported values so far "median" and "mean".
get_theta <- function(fit, type){
  if (type == "median"){theta <- fit$summary$quantiles[, "50%"]}
  if (type == "mean"){theta <- fit$summary$statistics[, "Mean"]}
  if (is.numeric(type)){
    chainidx <- floor(type / 1000) + 1
    sampleinchain <- type - 1000 * (chainidx - 1)
    theta <- fit$mcmc[[chainidx]][sampleinchain, ]
  }
  return(theta)
}

#' @describeIn A quick replacement to runjags::list.format that does nothing if the argument is already a list.
#' Takes the same arguments as runjags::list.format
as.list.format <- function(data, checkvalid = TRUE){
  if ("list" %in% class(data)){return(data)}
  out <- runjags::list.format(data, checkvalid = checkvalid)
  return(out)
}

drop_to_matrix <- function(array, dimdrop = 3){
  stopifnot(dim(array)[dimdrop] == 1) #must use subsetting [,,1, drop = FALSE] first
  if ((dim(array))[[1]] == 1){
    return(matrix(array, nrow = 1))
  } else if ((dim(array))[[2]] == 1){
    return(matrix(array, ncol = 1))
  } else {
    return(drop(array))
  }
}
