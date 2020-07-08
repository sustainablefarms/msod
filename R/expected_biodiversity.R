#' @title Expected Biodiversity
#' @description The expected number of species occupying a ModelSite, or the expected number of species detected at a ModelSite.
#' 
#' @examples 
#' fit <- readRDS("./tmpdata/7_2_10_clim_AnnTemp_AnnPrec.rds")
#' theta <- get_theta(fit, type = 1)
#' Xocc <- fit$data$Xocc[2, , drop = FALSE]
#' Xobs <- fit$data$Xobs[fit$data$ModelSite == 2, , drop = FALSE]
#' numspecies <- fit$data$n
#' lvsim <- matrix(rnorm(2 * 1), ncol = 2, nrow = 2) #dummy lvsim vars
#' lv.coef.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$n, ncol = 2), "lv.coef")
#' theta <- c(theta, lv.coef.bugs)
#' 
#' @param Xocc A matrix of occupancy covariates. Must have a single row. Columns correspond to covariates.
#' @param Xobs A matrix of detection covariates, each row is a visit. If NULL then expected number of species in occupation is returned
#' @param theta A vector of model parameters, labelled according to the BUGS labelling convention seen in runjags
#' @export
expectedspeciesnum.ModelSite.theta <- function(Xocc, Xobs = NULL, numspecies, theta, LVvals){
  stopifnot(nrow(Xocc) == 1)
  Xocc <- as.matrix(Xocc)
  
  u.b <- bugsvar2array(theta, "u.b", 1:numspecies, 1:ncol(Xocc))[,,1]  # rows are species, columns are occupancy covariates
  lv.coef <- bugsvar2array(theta, "lv.coef", 1:numspecies, 1:ncol(LVvals))[,,1] # rows are species, columns are lv
  sd_u_condlv <- sqrt(1 - rowSums(lv.coef^2)) #for each species the standard deviation of the indicator random variable 'u', conditional on values of LV

  
  ## Probability of Site Occupancy
  # external
  ModelSite.Occ.eta_external <- as.matrix(Xocc) %*% t(u.b) #columns are species
  
  # probability of occupancy given LV
  ModelSite.Occ.eta_LV <- LVvals %*% t(lv.coef) #occupancy contribution from latent variables, performed all together
  ModelSite.Occ.eta <- Rfast::eachrow(ModelSite.Occ.eta_LV, ModelSite.Occ.eta_external, oper = "+") #add the external part to each simulation
  # Make u standard deviations equal to 1 by dividing other values by sd
  # P(u < -ModelSite.Occ.eta) = P(u / sd < -ModelSite.Occ.eta / sd) = P(z < -ModelSite.Occ.eta / sd)
  ModelSite.Occ.eta_standardised <- Rfast::eachrow(ModelSite.Occ.eta, sd_u_condlv, oper = "/") 
  ModelSite.Occ.Pred.CondLV <- 1 - pnorm(-ModelSite.Occ.eta_standardised, mean = 0, sd = 1)
  
  ## Expected number of species occupying modelsite, given model and theta, and marginal across LV
  Enumspecies_occ <- mean(Rfast::rowsums(ModelSite.Occ.Pred.CondLV))
  if (is.null(Xobs)){
    return(Enumspeciesocc = Enumspecies_occ)
  }
  
  ## Probability of Detection, CONDITIONAL on occupied
  if (!is.null(Xobs)){
    Xobs <- as.matrix(Xobs)
    v.b <- bugsvar2array(theta, "v.b", 1:numspecies, 1:ncol(Xobs))[,,1]  # rows are species, columns are observation covariates
    Detection.LinPred <- as.matrix(Xobs) %*% t(v.b)
    Detection.Pred.Cond <- exp(Detection.LinPred) / (exp(Detection.LinPred) + 1)   #this is the inverse logit function
    # probability of no detections, given occupied
    NoDetections.Pred.Cond <- Rfast::colprods(1 - Detection.Pred.Cond)
  
    ## probability of no detections, marginal on occupancy
    NoDetections.Pred.marg_1 <- Rfast::eachrow(ModelSite.Occ.Pred.CondLV, NoDetections.Pred.Cond, oper = "*")  #occupied component
    NoDetections.Pred.marg_2 <- 1 - ModelSite.Occ.Pred.CondLV  #plus unoccupied component
    NoDetections.Pred.marg <- NoDetections.Pred.marg_1 + NoDetections.Pred.marg_2
    Enumspecies_detected <- numspecies - Rfast::rowsums(NoDetections.Pred.marg)
    Enumspecies_detected_margLV <- mean(Enumspecies_detected)
    return(Enumdetected = Enumspecies_detected_margLV)
  }
}

#' @describeIn expectedspeciesnum For ModelSite information, predicts both the expected number of species in occupation, and the expected number of species detected.
#' @param LVvals A matrix of LV values. Each column corresponds to a LV. To condition on specific LV values, provide a matrix of row 1.
#' If LVvals isn't provided then expections are marginalised across possible LVvalues through simulation
expectedspeciesnum.ModelSite <- function(fit, Xocc, Xobs, chains = NULL, LVvals = NULL){
  stopifnot(nrow(Xocc) == 1)
  Xocc <- as.matrix(Xocc)
  Xobs <- as.matrix(Xobs)
  if (is.null(chains)){chains <- 1:length(fit$mcmc)}
  draws <- do.call(rbind, fit$mcmc[chains])
  numspecies <- fit$data$n
  
  if (is.null(LVvals)) {
    if ( (is.null(fit$data$nlv)) || (fit$data$nlv == 0)){ #make dummy lvsim and and 0 loadings to draws
      LVvals <- matrix(rnorm(2 * 1), ncol = 2, nrow = 2) #dummy lvsim vars
      lv.coef.bugs <- matrix2bugsvar(matrix(0, nrow = numspecies, ncol = 2), "lv.coef")
      lv.coef.draws <- Rfast::rep_row(lv.coef.bugs, nrow(draws))
      colnames(lv.coef.draws) <- names(lv.coef.bugs)
      draws <- cbind(draws, lv.coef.draws)
    } else {
      LVvals <- matrix(rnorm(fit$data$nlv * 1000), ncol = fit$data$nlv, nrow = 1000) #simulated lv values, should average over thousands
    }
  }
  
  Expectedspeciesnum.ModelSite.alltheta <- apply(draws, 1,
      function(theta) expectedspeciesnum.ModelSite.theta(Xocc, Xobs,
                                                         numspecies = numspecies,
                                                         theta = theta,
                                                         LVvals = LVvals))
  # Expectedspeciesnum.ModelSite.alltheta is a matrix with two rows, each column is an entry from 'draws'.
  # First row is expected number of species that will be detected, second row is expected number of species occupied
  out <- mean(Expectedspeciesnum.ModelSite.alltheta)
  return(out)
}
