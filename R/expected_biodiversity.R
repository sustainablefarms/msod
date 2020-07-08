#' @title Expected Biodiversity
#' @description The expected number of species occupying a ModelSite, or the expected number of species detected at a ModelSite.
#' 
#' @examples 
#' fit <- readRDS("./tmpdata/7_2_10_clim_AnnTemp_AnnPrec.rds")
#' theta <- get_theta(fit, type = 1)
#' Xocc <- fit$data$Xocc[2, , drop = FALSE]
#' Xobs <- fit$data$Xobs[fit$data$ModelSite == 2, , drop = FALSE]
#' numspecies <- fit$data$n
#' 
#' @param Xocc A matrix of occupancy covariates. Must have a single row. Columns correspond to covariates.
#' @param Xobs A matrix of detection covariates, each row is a visit.
#' @param theta A vector of model parameters, labelled according to the BUGS labelling convention seen in runjags
#' @param lvsim A matrix of simulated LV values. Columns correspond to latent variables, each row is a simulation
#' @export
expectedspeciesnum.ModelSite.theta <- function(Xocc, Xobs, numspecies, theta){
  u.b <- bugsvar2array(theta, "u.b", 1:numspecies, 1:ncol(Xocc))[,,1]  # rows are species, columns are occupancy covariates
  v.b <- bugsvar2array(theta, "v.b", 1:numspecies, 1:ncol(Xobs))[,,1]  # rows are species, columns are observation covariates
  Xocc <- as.matrix(Xocc)
  Xobs <- as.matrix(Xobs)
  
  ## Probability of Site Occupancy without LV
  ModelSite.Occ.eta_external <- as.matrix(Xocc) %*% t(u.b) #columns are species
  ModelSite.Occ.Pred <- 1 - pnorm(-ModelSite.Occ.eta_external, mean = 0, sd = 1)
  
  ## Expected number of species occupying modelsite, given model and theta
  Enumspecies_occ <- Rfast::rowsums(ModelSite.Occ.Pred)
  
  
  ## Probability of Detection, CONDITIONAL on occupied
  Detection.LinPred <- as.matrix(Xobs) %*% t(v.b)
  Detection.Pred.Cond <- exp(Detection.LinPred) / (exp(Detection.LinPred) + 1)   #this is the inverse logit function
  # probability of no detections, given occupied
  NoDetections.Pred.Cond <- Rfast::colprods(1 - Detection.Pred.Cond)
  
  ## probability of no detections, marginal on occupancy
  NoDetections.Pred.marg_1 <- Rfast::eachrow(ModelSite.Occ.Pred, NoDetections.Pred.Cond)  #occupied component
  NoDetections.Pred.marg_2 <- 1 - ModelSite.Occ.Pred  #plus unoccupied component
  NoDetections.Pred.marg <- NoDetections.Pred.marg_1 + NoDetections.Pred.marg_2
  Enumspecies_detected <- numspecies - sum(NoDetections.Pred.marg)
  return(Enumspecies_detected)
}
