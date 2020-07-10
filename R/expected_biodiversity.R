#' @title Expected Biodiversity
#' @description The expected number of species occupying a ModelSite, or the expected number of species detected at a ModelSite.
#' 
#' @examples 
#' fit <- readRDS("./tmpdata/7_2_9_addyear_msnm_year_time_2lv.rds")
#' theta <- get_theta(fit, type = 1)
#' Xocc <- fit$data$Xocc[2, , drop = FALSE]
#' Xobs <- fit$data$Xobs[fit$data$ModelSite == 2, , drop = FALSE]
#' numspecies <- fit$data$n
#' lvsim <- matrix(rnorm(2 * 1), ncol = 2, nrow = 2) #dummy lvsim vars
#' lv.coef.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$n, ncol = 2), "lv.coef")
#' theta <- c(theta, lv.coef.bugs)
#' 
#' expectedspeciesnum.ModelSiteIdx(fit, 2, chains = NULL, LVvals = NULL)
#' Enumspec <- predsumspecies(fit, usefittedLV = TRUE)
#' 
#' @param Xocc A matrix of occupancy covariates. Must have a single row. Columns correspond to covariates.
#' @param Xobs A matrix of detection covariates, each row is a visit. If NULL then expected number of species in occupation is returned
#' @param theta A vector of model parameters, labelled according to the BUGS labelling convention seen in runjags
#' @export
expectedspeciesnum.ModelSite.theta <- function(Xocc, Xobs = NULL, u.b, v.b = NULL, lv.coef = NULL, LVvals = NULL){
  ## Probability of Site Occupancy
  stopifnot(nrow(Xocc) == 1)
  if (is.null(Xobs)) {stopifnot(is.null(v.b))}
  Xocc <- as.matrix(Xocc)
  if (is.null(lv.coef)){ # create dummy LV  
    stopifnot(is.null(LVvals))
    lv.coef <- matrix(0, rows = nrow(u.b), ncol = 2)
    LVvals <- matrix(0, rows = 2, ncol = 2)
    } # dummy LV
  ModelSite.Occ.Pred.CondLV <- poccupy.ModelSite.theta(Xocc, u.b, lv.coef, LVvals)
  
  ## Expected number of species occupying modelsite, given model and theta, and marginal across LV
  EVnumspec_occ <- Erowsum_margrow(ModelSite.Occ.Pred.CondLV)
  names(EVnumspec_occ) <- paste0(names(EVnumspec_occ), "_occ")
  if (is.null(Xobs)){
    return(EVnumspec_occ)
  }
  
  ## Probability of Detection
  if (!is.null(Xobs)){
    Detection.Pred.Cond <- pdetection_occupied.ModelSite.theta(Xobs, v.b) # probability conditional on occupied
    # probability of no detections, given occupied
    NoDetections.Pred.Cond <- Rfast::colprods(1 - Detection.Pred.Cond)

  
    ## probability of no detections, marginal on occupancy
    NoDetections.Pred.marg_1 <- Rfast::eachrow(ModelSite.Occ.Pred.CondLV, NoDetections.Pred.Cond, oper = "*")  #occupied component
    NoDetections.Pred.marg_2 <- 1 - ModelSite.Occ.Pred.CondLV  #plus unoccupied component
    AnyDetections.Pred.marg <- 1 - (NoDetections.Pred.marg_1 + NoDetections.Pred.marg_2)
    EVnumspec_det <- Erowsum_margrow(AnyDetections.Pred.marg)
    names(EVnumspec_det) <- paste0(names(EVnumspec_det), "_det")
    return(c(EVnumspec_occ, EVnumspec_det))
  }
}

#' @title Probability of detection and occupation for a given ModelSite
#' @param v.b Covariate loadings. Each row is a species, each column a detection covariate
#' @return A matrix of detection probabilities for each row of Xobs and each species.
#' Rows in returned value correspond to rows in Xobs, columns correspond to species.
#' @export
pdetection_occupied.ModelSite.theta <- function(Xobs, v.b){
  Xobs <- as.matrix(Xobs)
  Detection.LinPred <- as.matrix(Xobs) %*% t(v.b)
  Detection.Pred.Cond <- exp(Detection.LinPred) / (exp(Detection.LinPred) + 1)   #this is the inverse logit function
  return(Detection.Pred.Cond)
}

#' @describeIn pdetection_occupied.ModelSite.theta The probability of occupation for given LV values.
#' @param u.b Covariate loadings for occupancy. Each row is a species, each column an occupancy covariate.
#' @param lv.coef Loadings for the latent variables. Each row is a species, each column corresponds to a LV.
#' @return A matrix of occupancy probabilities. Each row corresponds to a row of LVvals, each column to a species.
poccupy.ModelSite.theta <- function(Xocc, u.b, lv.coef, LVvals){
  sd_u_condlv <- sqrt(1 - rowSums(lv.coef^2)) #for each species the standard deviation of the indicator random variable 'u', conditional on values of LV
  
  # external
  ModelSite.Occ.eta_external <- as.matrix(Xocc) %*% t(u.b) #columns are species
  
  # probability of occupancy given LV
  ModelSite.Occ.eta_LV <- LVvals %*% t(lv.coef) #occupancy contribution from latent variables, performed all together
  ModelSite.Occ.eta <- Rfast::eachrow(ModelSite.Occ.eta_LV, ModelSite.Occ.eta_external, oper = "+") #add the external part to each simulation
  # Make u standard deviations equal to 1 by dividing other values by sd
  # P(u < -ModelSite.Occ.eta) = P(u / sd < -ModelSite.Occ.eta / sd) = P(z < -ModelSite.Occ.eta / sd)
  ModelSite.Occ.eta_standardised <- Rfast::eachrow(ModelSite.Occ.eta, sd_u_condlv, oper = "/") 
  ModelSite.Occ.Pred.CondLV <- 1 - pnorm(-ModelSite.Occ.eta_standardised, mean = 0, sd = 1)
  return(ModelSite.Occ.Pred.CondLV)
}

#' @param pmat A matrix of success probabilities for Bernoulli random variables.
#' Each column corresponds to an independent Bernoulli random variable, each row is a different set of parameters.
#' @return
#' The expected sum of the Bernoulli random variables, marginal across rows, E[E[N | L]]
#' The variance of the sum, marginal across rows E[N^2] - E[N]^2 = E[V[N | L] + E[N | L]^2] - E[E[N | L]]^2
Erowsum_margrow <- function(pmat){
  ## Expected sum of the Bournilli for each parameter set
  En_row <- Rfast::rowsums(pmat)
  ## Variance of sum for each parameter set (variance adds for species when independent)
  Vn_row <- Rfast::rowsums(pmat * (1 - pmat))
  ## second moment of sum for each row
  M2n_row <- Vn_row + En_row^2
  ## second moment of sum, marginal rows
  M2n <- mean(M2n_row)
  ## 1st moment of sum, marginal rows
  En <- mean(En_row)
  ## Variance of sum, marginal rows
  Vn <- M2n - En^2
  return(c(Esum = En, Vsum = Vn))
}

# use supplied LVs for modelsite
# use fitted LVs for modelsite for each theta
# use simulated LVs for each theta

#' @describeIn expectedspeciesnum For ModelSite information, predicts both the expected number of species in occupation, and the expected number of species detected.
#' @param LVvals A matrix of LV values. Each column corresponds to a LV. To condition on specific LV values, provide a matrix of row 1.
#' If LVvals isn't provided then the fitted LV values for each draw are used
expectedspeciesnum.ModelSiteIdx <- function(fit, ModelSiteIdx, chains = NULL, LVvals = NULL){
  Xocc <- as.matrix(fit$data$Xocc[ModelSiteIdx, , drop = FALSE])
  stopifnot(nrow(Xocc) == 1)
  Xobs <- as.matrix(fit$data$Xobs[fit$data$ModelSite == ModelSiteIdx, , drop = FALSE])
  if (is.null(chains)){chains <- 1:length(fit$mcmc)}
  draws <- do.call(rbind, fit$mcmc[chains])
  numspecies <- fit$data$n
  
  if (!is.null(LVvals)) {
    Moms.ModelSite.alltheta <- apply(draws, 1,
                                     function(theta) expectedspeciesnum.ModelSite.theta(Xocc, Xobs,
                                                                                        numspecies = numspecies,
                                                                                        theta = theta,
                                                                                        LVvals = LVvals))
  }
  
  if (is.null(LVvals)) {
    if ( (is.null(fit$data$nlv)) || (fit$data$nlv == 0)){ #make dummy lvsim and and 0 loadings to draws
      LVvals <- matrix(rnorm(2 * 1), ncol = 2, nrow = 2) #dummy lvsim vars
      lv.coef.bugs <- matrix2bugsvar(matrix(0, nrow = numspecies, ncol = 2), "lv.coef")
      lv.coef.draws <- Rfast::rep_row(lv.coef.bugs, nrow(draws))
      colnames(lv.coef.draws) <- names(lv.coef.bugs)
      draws <- cbind(draws, lv.coef.draws)
      Moms.ModelSite.alltheta <- apply(draws, 1,
                                       function(theta) expectedspeciesnum.ModelSite.theta(Xocc, Xobs,
                                                                                          numspecies = numspecies,
                                                                                          theta = theta,
                                                                                          LVvals = LVvals))
    } else { #LVs are part of model!
      Moms.ModelSite.alltheta <- apply(draws, 1,
                                       function(theta) {
                                         LVattheta <- bugsvar2matrix(theta, "LV", 1:nrow(fit$data$Xocc),
                                                                     1:fit$data$nlv)[ModelSiteIdx, , drop = FALSE]
                                         out <- expectedspeciesnum.ModelSite.theta(Xocc, Xobs,
                                                         numspecies = numspecies,
                                                         theta = theta,
                                                         LVvals = LVattheta)
                                         return(out)}
      )
    }
  }
  

  
  M2n_det_theta <- Moms.ModelSite.alltheta["Vsum_det", ] + Moms.ModelSite.alltheta["Esum_det", , drop = FALSE]^2
  M2n_det <- Rfast::rowmeans(M2n_det_theta)
  En_det <- Rfast::rowmeans(Moms.ModelSite.alltheta["Esum_det", , drop = FALSE])
  Vn_det <- M2n_det - En_det^2
  
  M2n_occ_theta <- Moms.ModelSite.alltheta["Vsum_occ", ] + Moms.ModelSite.alltheta["Esum_occ", , drop = FALSE]^2
  M2n_occ <- Rfast::rowmeans(M2n_occ_theta)
  En_occ <- Rfast::rowmeans(Moms.ModelSite.alltheta["Esum_occ", , drop = FALSE])
  Vn_occ <- M2n_occ - En_occ^2
  return(c(Esum_occ = En_occ,
           Vsum_occ = Vn_occ,
           Esum_det = En_det,
           Vsum_det = Vn_det))
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
  
  Moms.ModelSite.alltheta <- apply(draws, 1,
      function(theta) expectedspeciesnum.ModelSite.theta(Xocc, Xobs,
                                                         numspecies = numspecies,
                                                         theta = theta,
                                                         LVvals = LVvals))
  
  M2n_det_theta <- Moms.ModelSite.alltheta["Vsum_det", ] + Moms.ModelSite.alltheta["Esum_det", , drop = FALSE]^2
  M2n_det <- Rfast::rowmeans(M2n_det_theta)
  En_det <- Rfast::rowmeans(Moms.ModelSite.alltheta["Esum_det", , drop = FALSE])
  Vn_det <- M2n_det - En_det^2
  
  M2n_occ_theta <- Moms.ModelSite.alltheta["Vsum_occ", ] + Moms.ModelSite.alltheta["Esum_occ", , drop = FALSE]^2
  M2n_occ <- Rfast::rowmeans(M2n_occ_theta)
  En_occ <- Rfast::rowmeans(Moms.ModelSite.alltheta["Esum_occ", , drop = FALSE])
  Vn_occ <- M2n_occ - En_occ^2
  return(c(Esum_occ = En_occ,
           Vsum_occ = Vn_occ,
           Esum_det = En_det,
           Vsum_det = Vn_det))
}

#' @param usefittedLV If TRUE the fitted LV variables are used, if false then 1000 LV values are simulated.
#' @param chains The chains of MCMC to use. Default is all chains.
#' @param cl A cluster object created by parallel::makeCluster. If NULL no cluster is used.
predsumspecies <- function(fit, chains = NULL, usefittedLV = TRUE, cl = NULL){
  fit$data <- as_list_format(fit$data)
  if (is.null(chains)){chains <- 1:length(fit$mcmc)}
  draws <- do.call(rbind, fit$mcmc[chains])
  u.b <- bugsvar2array(draws, "u.b", 1:fit$data$n, 1:ncol(fit$data$Xocc))
  v.b <- bugsvar2array(draws, "v.b", 1:fit$data$n, 1:ncol(fit$data$Xobs))
  if ( (is.null(fit$data$nlv)) || (fit$data$nlv == 0)){ #LVs not in model, add dummy variables
    LVvals <- array(0, dim = c(nrow(fit$data$Xocc), 2, nrow(draws))) #dummy LVvals
    lv.coef <- array(0, dim = c(fit$data$n, 2, nrow(draws)))
    fit$data$nlv <- 2
    usefittedLV <- TRUE #calculations faster when not simulating 1000s of LV values, especially since they are all ignored here.
  } else {
    lv.coef <- bugsvar2array(draws, "lv.coef", 1:fit$data$n, 1:fit$data$nlv)
    LVvals <- bugsvar2array(draws, "LV", 1:nrow(fit$data$Xocc), 1:fit$data$nlv)
  }
  
  sitedrawidxs <- expand.grid(siteidx = 1:nrow(fit$data$Xocc), drawidx = 1:nrow(draws)) 
  
  if (!usefittedLV){ # predicting as if LVs not known, so simulate from their distribution
    lvsim <- matrix(rnorm(fit$data$nlv * 1000), ncol = fit$data$nlv, nrow = 1000)
  }
  

  # for each modelsite and each draw apply the following function:
  Enumspec <- pbapply::pbapply(sitedrawidxs, MARGIN = 1,
        function(sitedrawidx){
          Xocc <- fit$data$Xocc[sitedrawidx[["siteidx"]], , drop = FALSE]
          Xobs <- fit$data$Xobs[fit$data$ModelSite == sitedrawidx[["siteidx"]], , drop = FALSE]
          u.b_theta <- matrix(u.b[,, sitedrawidx[["drawidx"]] ], nrow = fit$data$n, ncol = ncol(Xocc))
          v.b_theta <- matrix(v.b[,, sitedrawidx[["drawidx"]] ], nrow = fit$data$n, ncol = ncol(Xobs))
          lv.coef_theta <- matrix(lv.coef[,, sitedrawidx[["drawidx"]] ], nrow = fit$data$n, ncol = fit$data$nlv)
          if (usefittedLV){
            LVvals_thetasite <- matrix(LVvals[sitedrawidx[["siteidx"]], , sitedrawidx[["drawidx"]], drop = FALSE], nrow = 1, ncol = ncol(lv.coef_theta))
          } else {
            LVvals_thetasite <- lvsim
          }
          Enumspec_sitetheta <- expectedspeciesnum.ModelSite.theta(Xocc = Xocc, Xobs = Xobs,
                                                                   u.b = u.b_theta,
                                                                   v.b = v.b_theta,
                                                                   lv.coef = lv.coef_theta,
                                                                   LVvals = LVvals_thetasite)
        },
        cl = cl)
  # each column of Enumspec is a model site
  # make each row a draw, each column a site.
  return(Enumspec)
}
