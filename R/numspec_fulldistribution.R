#' @title Compute predicted distributions of number species detected
#' @description Uses [discreteRV] to compute distribution of number of species 
#' @examples
#' fit <- readRDS("./tmpdata/7_3_02_clim_someclimate_year_woody500m_msnm_det1stO.rds")
#' indata <- readRDS("./private/data/clean/7_2_10_input_data.rds")
#' det_dist_2060 <- predsumspeciesRV_newdata(fit,
#'                          Xocc = indata$holdoutdata$Xocc[indata$holdoutdata$Xocc$ModelSiteID == 2060, , drop = FALSE],
#'                          Xobs = indata$holdoutdata$yXobs[indata$holdoutdata$yXobs$ModelSiteID == 2060, , drop = FALSE],
#'                          ModelSiteVars = "ModelSiteID")

#' @param fit A runjags fitted object
#' @param UseFittedLV If TRUE the fitted LV variables are used, if false then 100 LV values are simulated.
#' @param chains The chains of MCMC to use. Default is all chains.
#' @param cl A cluster object created by parallel::makeCluster. If NULL no cluster is used.
#' @param type If "draws" then a list of distributions for each site and draw is returned
#' If "marginal" then a single distribution for each site is returned and this distribution is marginalised over the draws.
#' @export
predsumspeciesRV <- function(fit, chains = NULL, UseFittedLV = TRUE, nLVsim = 100, type = "marginal", cl = NULL){
  stopifnot(type %in% c("draws", "marginal"))
  
  fit$data <- as_list_format(fit$data)
  if (is.null(chains)){chains <- 1:length(fit$mcmc)}
  draws <- do.call(rbind, fit$mcmc[chains])
  if ( (is.null(fit$data$nlv)) || (fit$data$nlv == 0)){ #LVs not in model, add dummy variables
    LVbugs <- matrix2bugsvar(matrix(0, nrow = nrow(fit$data$Xocc), ncol = 2), "LV")
    LVbugs.draws <- Rfast::rep_row(LVbugs, nrow(draws))
    colnames(LVbugs.draws) <- names(LVbugs)
    
    ldet.b.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$nspecies, ncol = 2), "ldet.b")
    ldet.b.draws <- Rfast::rep_row(ldet.b.bugs, nrow(draws))
    colnames(ldet.b.draws) <- names(ldet.b.bugs)
    draws <- cbind(draws, ldet.b.draws, LVbugs.draws)
    fit$data$nlv <- 2
    UseFittedLV <- TRUE #calculations faster when not simulating 1000s of LV values, especially since they are all ignored here.
  } 
  
  if (UseFittedLV){nLVsim <- NULL} #don't pass number of simulations if not going to use them
  
  numspecRV_drawsites <- sumspeciesRV_raw(
    Xocc = fit$data$Xocc,
    Xobs = fit$data$Xobs,
    ModelSite = fit$data$ModelSite,
    numspecies = fit$data$nspecies,
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

#' @describeIn predsumspeciesRV Distribution of species numbers for new ModelSite
predsumspeciesRV_newdata <- function(fit, Xocc, Xobs = NULL, ModelSiteVars = NULL, chains = NULL, nLVsim = 100, type = "marginal", cl = NULL){
  stopifnot(type %in% c("draws", "marginal"))
  
  if (!is.null(Xobs)) {datalist <- prep_new_data(fit, Xocc, Xobs, ModelSite = ModelSiteVars)}
  else{datalist <- list(
    Xocc = apply.designmatprocess(fit$XoccProcess, Xocc),
    Xobs = NULL,
    ModelSite = NULL
  )}
  UseFittedLV <- FALSE # no LV available for new model sites
  
  fit$data <- as_list_format(fit$data)
  if (is.null(chains)){chains <- 1:length(fit$mcmc)}
  draws <- do.call(rbind, fit$mcmc[chains])
  
  if ( (is.null(fit$data$nlv)) || (fit$data$nlv == 0)){ #LVs not in model, add dummy variables
    ldet.b.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$nspecies, ncol = 2), "ldet.b")
    ldet.b.draws <- Rfast::rep_row(ldet.b.bugs, nrow(draws))
    colnames(ldet.b.draws) <- names(ldet.b.bugs)
    draws <- cbind(draws, ldet.b.draws)
    fit$data$nlv <- 2
    nLVsim = 2 #calculations faster when not simulating 1000s of LV values, especially since they are all ignored here.
  }
  
  numspecRV_drawsites <- sumspeciesRV_raw(
    Xocc = datalist$Xocc,
    Xobs = datalist$Xobs,
    ModelSite = datalist$ModelSite,
    numspecies = fit$data$nspecies,
    nlv = fit$data$nlv,
    draws = draws,
    useLVindraws = FALSE,
    nLVsim = nLVsim,
    cl = cl
  )
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
  return(out)
}


#' @param Xocc A matrix of occupancy covariates. Must have a single row. Columns correspond to covariates.
#' @param Xobs A matrix of detection covariates, each row is a visit. If NULL then expected number of species in occupation is returned
#' @param theta A vector of model parameters, labelled according to the BUGS labelling convention seen in runjags
#' @param LVvals A matrix of LV values. Each column corresponds to a LV. To condition on specific LV values, provide a matrix of row 1.
#' @describeIn predsumspeciesRV For raw model information. Returns distributions per draw.
#' @export
sumspeciesRV_raw <- function(Xocc, Xobs = NULL, ModelSite = NULL, numspecies, nlv, draws, useLVindraws = TRUE, nLVsim = NULL, cl = NULL){
  # prepare parameters
  nspecall <- numspecies
  ndraws <- nrow(draws)
  nsites <- nrow(Xocc)
  noccvar <- ncol(Xocc)
  if (!is.null(Xobs)) {
    nobsvar <- ncol(Xobs)
    ModelSiteIdxs <- ModelSite
    stopifnot(length(ModelSiteIdxs) == nrow(Xobs))
    stopifnot(all(ModelSiteIdxs %in% 1:nsites))
    if (!all(1:nsites %in% ModelSiteIdxs)){warning("Some ModelSite do not have observation covariate information.")}
  }
  occ.b <- bugsvar2array(draws, "occ.b", 1:nspecall, 1:noccvar)
  if (!is.null(Xobs)) {det.b <- bugsvar2array(draws, "det.b", 1:nspecall, 1:nobsvar)}
  ldet.b <- bugsvar2array(draws, "ldet.b", 1:nspecall, 1:nlv)
  
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
                                 if (!is.null(Xobs)) {Xobs <- Xobs[ModelSiteIdxs == sitedrawidx[["siteidx"]], , drop = FALSE]}
                                 occ.b_theta <- matrix(occ.b[,, sitedrawidx[["drawidx"]] ], nrow = nspecall, ncol = noccvar)
                                 if (!is.null(Xobs)) {det.b_theta <- matrix(det.b[,, sitedrawidx[["drawidx"]] ], nrow = nspecall, ncol = nobsvar)}
                                 else {det.b_theta = NULL}
                                 ldet.b_theta <- matrix(ldet.b[,, sitedrawidx[["drawidx"]] ], nrow = nspecall, ncol = nlv)
                                 if (useLVindraws){
                                   LVvals_thetasite <- matrix(LVvals[sitedrawidx[["siteidx"]], , sitedrawidx[["drawidx"]], drop = FALSE], nrow = 1, ncol = nlv)
                                 } else {
                                   LVvals_thetasite <- lvsim
                                 }
                                 sumRV <- speciesnum.ModelSite.theta(Xocc = Xocc, Xobs = Xobs,
                                                                                          occ.b = occ.b_theta,
                                                                                          det.b = det.b_theta,
                                                                                          ldet.b = ldet.b_theta,
                                                                                          LVvals = LVvals_thetasite)
                                 return(sumRV)
                               },
                               cl = cl)
  
  return(list(
    numspecRVs = numspecRVs,
    sitedrawidxs = sitedrawidxs
  ))
}

speciesnum.ModelSite.theta <- function(Xocc, Xobs = NULL, occ.b, det.b = NULL, ldet.b = NULL, LVvals = NULL){
  ## Probability of Site Occupancy
  stopifnot(nrow(Xocc) == 1)
  if (is.null(Xobs)) {stopifnot(is.null(det.b))}
  Xocc <- as.matrix(Xocc)
  ModelSite.Occ.Pred.CondLV <- poccupy.ModelSite.theta(Xocc, occ.b, ldet.b, LVvals)
  
  if (is.null(Xobs)){
    sum_occ_margLV <- sumRV_margrow(ModelSite.Occ.Pred.CondLV)
    return(sum_occ_margLV)
  }
  
  ## Probability of Detection
  if (!is.null(Xobs)){
    Detection.Pred.Cond <- pdetection_occupied.ModelSite.theta(Xobs, det.b) # probability conditional on occupied
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
