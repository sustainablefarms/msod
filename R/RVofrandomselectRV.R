#' @title For a list of random variables, return random variable that is a random selection one of the RVs
#' @param rvs A list of Random Variable (distributions) of class 'RV'
#' @examples 
#' library(discreteRV)
#' X1 <- RV(c(1,0), c(.5,.5))
#' X2 <- RV(c(1,0), c(.5,.5))
#' X3 <- RV(c(1,-1), c(.5,.5))
#' randselRV(list(X1, X2, X3))

#' full distribution of species number given fixed theta and unknown LV value (i.e. many simulated LV values). This is *within* model uncertainy to biodiversity
#' fit <- readRDS("./tmpdata/7_2_9_addyear_msnm_year_time_2lv.rds")
#' theta <- get_theta(fit, type = "median")
#' u.b <- bugsvar2matrix(theta, "u.b", 1:fit$data$n, 1:ncol(fit$data$Xocc))
#' v.b <- bugsvar2matrix(theta, "v.b", 1:fit$data$n, 1:ncol(fit$data$Xobs))
#' lv.coef <- bugsvar2matrix(theta, "lv.coef", 1:fit$data$n, 1:fit$data$nlv)
#' nLVsim <- 1000
#' lvsim <- matrix(rnorm(fit$data$nlv * nLVsim), ncol = fit$data$nlv, nrow = nLVsim)
#' LVvals <- bugsvar2matrix(theta, "LV", 1:fit$data$n, 1:fit$data$nlv)
#' sumRV_margLV <- speciesnum.ModelSite.theta(Xocc = fit$data$Xocc[1, , drop = FALSE], u.b = u.b, lv.coef = lv.coef, LVvals = lvsim)
#' quantile(sumRV_margLV, probs = c(0.025, 0.5, 0.975))
#' sumRV_fitLV <- speciesnum.ModelSite.theta(Xocc = fit$data$Xocc[1, , drop = FALSE], u.b = u.b, lv.coef = lv.coef, LVvals = LVvals[1, , drop = FALSE])
#' quantile(sumRV_fitLV, probs = c(0.025, 0.5, 0.975))



#' Uncertainty of biodiversity *incorporating* parameter uncertainty is below
predsumspeciesRV <- function(fit, chains = NULL, UseFittedLV = TRUE, nLVsim = 1000, type = "marginal", cl = NULL){
  stopifnot(type %in% c("draws", "marginal"))
  
  fit$data <- as_list_format(fit$data)
  if (is.null(chains)){chains <- 1:length(fit$mcmc)}
  draws <- do.call(rbind, fit$mcmc[chains])
  if ( (is.null(fit$data$nlv)) || (fit$data$nlv == 0)){ #LVs not in model, add dummy variables
    LVbugs <- matrix2bugsvar(matrix(0, nrow = nrow(fit$data$Xocc), ncol = 2), "LV")
    LVbugs.draws <- Rfast::rep_row(LVbugs, nrow(draws))
    colnames(LVbugs.draws) <- names(LVbugs)
    
    lv.coef.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$n, ncol = 2), "lv.coef")
    lv.coef.draws <- Rfast::rep_row(lv.coef.bugs, nrow(draws))
    colnames(lv.coef.draws) <- names(lv.coef.bugs)
    draws <- cbind(draws, lv.coef.draws, LVbugs.draws)
    fit$data$nlv <- 2
    UseFittedLV <- TRUE #calculations faster when not simulating 1000s of LV values, especially since they are all ignored here.
  } 
  
  if (UseFittedLV){nLVsim <- NULL} #don't pass number of simulations if not going to use them
  
  numspecRV_drawsites <- sumspeciesRV_raw(
    Xocc = fit$data$Xocc,
    Xobs = fit$data$Xobs,
    ModelSite = fit$data$ModelSite,
    numspecies = fit$data$n,
    nlv = fit$data$nlv,
    draws = draws,
    useLVindraws = UseFittedLV,
    nLVsim = nLVsim,
    cl = cl
  )
  # convert predictions for each site and theta into predictions for each site, marginal across theta distribution
  if (type == "draws") {out <- numspecRV_drawsites}
  if (type == "marginal") {
    numspecRVs <- numspecRV_drawsites[["numspecRVs"]]
    sitedrawidxs <- numspecRV_drawsites[["sitedrawidxs"]]
    sumRVs_margpost <- pbapply::pblapply(unique(sitedrawidxs[["siteidx"]]),
                             function(siteidx){
                               sumRV_margpost <- randselRV(numspecRVs[sitedrawidxs[["siteidx"]] == siteidx])
                               return(sumRV_margpost)
                             },
                             cl = cl)
    return(sumRVs_margpost)
  }
}


sumspeciesRV_raw <- function(Xocc, Xobs, ModelSite, numspecies, nlv, draws, useLVindraws = TRUE, nLVsim = NULL, cl = NULL){
  # prepare parameters
  nspecall <- numspecies
  ndraws <- nrow(draws)
  nsites <- nrow(Xocc)
  noccvar <- ncol(Xocc)
  nobsvar <- ncol(Xobs)
  ModelSiteIdxs <- ModelSite
  stopifnot(length(ModelSiteIdxs) == nrow(Xobs))
  stopifnot(all(ModelSiteIdxs %in% 1:nsites))
  if (!all(1:nsites %in% ModelSiteIdxs)){warning("Some ModelSite do not have observation covariate information.")}
  u.b <- bugsvar2array(draws, "u.b", 1:nspecall, 1:noccvar)
  v.b <- bugsvar2array(draws, "v.b", 1:nspecall, 1:nobsvar)
  lv.coef <- bugsvar2array(draws, "lv.coef", 1:nspecall, 1:nlv)
  
  if (useLVindraws){stopifnot(is.null(nLVsim))}
  if (!useLVindraws){stopifnot(is.numeric(nLVsim))}
  
  
  sitedrawidxs <- expand.grid(siteidx = 1:nsites, drawidx = 1:ndraws) 
  
  if (!useLVindraws){ # predicting as if LVs not known, so simulate from their distribution
    lvsim <- matrix(rnorm(nlv * nLVsim), ncol = nlv, nrow = nLVsim)
  } else {
    LVvals <- bugsvar2array(draws, "LV", 1:nsites, 1:nlv)
  }
  
  
  # for each modelsite and each draw apply the following function:
  # note that I've had to use pblapply because pbapply was doing some weird simplifications
  numspecRVs <- pbapply::pblapply(1:nrow(sitedrawidxs), 
                               function(sitedrawidx_row){
                                 sitedrawidx <- sitedrawidxs[sitedrawidx_row, , drop = FALSE]
                                 Xocc <- Xocc[sitedrawidx[["siteidx"]], , drop = FALSE]
                                 Xobs <- Xobs[ModelSiteIdxs == sitedrawidx[["siteidx"]], , drop = FALSE]
                                 u.b_theta <- matrix(u.b[,, sitedrawidx[["drawidx"]] ], nrow = nspecall, ncol = noccvar)
                                 v.b_theta <- matrix(v.b[,, sitedrawidx[["drawidx"]] ], nrow = nspecall, ncol = nobsvar)
                                 lv.coef_theta <- matrix(lv.coef[,, sitedrawidx[["drawidx"]] ], nrow = nspecall, ncol = nlv)
                                 if (useLVindraws){
                                   LVvals_thetasite <- matrix(LVvals[sitedrawidx[["siteidx"]], , sitedrawidx[["drawidx"]], drop = FALSE], nrow = 1, ncol = nlv)
                                 } else {
                                   LVvals_thetasite <- lvsim
                                 }
                                 sumRV <- speciesnum.ModelSite.theta(Xocc = Xocc, Xobs = Xobs,
                                                                                          u.b = u.b_theta,
                                                                                          v.b = v.b_theta,
                                                                                          lv.coef = lv.coef_theta,
                                                                                          LVvals = LVvals_thetasite)
                                 return(sumRV)
                               },
                               cl = cl)
  
  return(list(
    numspecRVs = numspecRVs,
    sitedrawidxs = sitedrawidxs
  ))
}

  # number of species distribution incorporating parameter uncertaintyd
  # sumRV_withmodeluncertainty <- randselRV(numspecRVs[sitedrawidxs[["siteidx"]] == 1])
  # quantile(sumRV_withmodeluncertainty, probs = c(0.025, 0.5, 0.975))

speciesnum.ModelSite.theta <- function(Xocc, Xobs = NULL, u.b, v.b = NULL, lv.coef = NULL, LVvals = NULL){
  ## Probability of Site Occupancy
  stopifnot(nrow(Xocc) == 1)
  if (is.null(Xobs)) {stopifnot(is.null(v.b))}
  Xocc <- as.matrix(Xocc)
  if (is.null(lv.coef)){ # create dummy LV  
    stopifnot(is.null(LVvals))
    lv.coef <- matrix(0, nrow = nrow(u.b), ncol = 2)
    LVvals <- matrix(0, nrow = 2, ncol = 2)
  } # dummy LV
  ModelSite.Occ.Pred.CondLV <- poccupy.ModelSite.theta(Xocc, u.b, lv.coef, LVvals)
  
  if (is.null(Xobs)){
    sum_occ_margLV <- sumRV_margrow(ModelSite.Occ.Pred.CondLV)
    return(sum_occ_margLV)
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
    sum_det_margLV <- sumRV_margrow(AnyDetections.Pred.marg)
    return(sum_det_margLV)
  }
}

sumRV_margrow <- function(pmat){
  sum_RVs <- lapply(1:nrow(pmat),
                                   function(row) {
                                     LVs <- lapply(pmat[row, ], function(p) RV(c(0, 1), probs = c(1 - p, p)))
                                     sumRV <- do.call(SofI, LVs)
                                     return(sumRV)
                                   })
  sumRV_marg <- randselRV(sum_RVs)
  return(sumRV_marg)
}

randselRV <- function(rvs, weights = rep(1, length(rvs))){
  outcomes <- unique(unlist(lapply(rvs, outcomes)))
  pmf <- vapply(outcomes, function(value) Poutcome(value, rvs, weights), FUN.VALUE = 0.2)
  outRV <- RV(outcomes,
              probs = pmf,
              fractions = FALSE)
  return(outRV)
}

# For a given value, compute the probability that a randomly selected RV will have that value
# weights are the selection weights of the RVs
Poutcome <- function(value, rvs, weights = rep(1, length(rvs))){
  Ps <- vapply(rvs, function(X) P(X == value), FUN.VALUE = 0.2)
  prob <- weighted.mean(Ps, w = weights)
  return(prob)
}
