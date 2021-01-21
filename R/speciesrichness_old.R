# OLD species richness code

colprods_arr <- function(arr){
  arr2 <- aperm(arr, c(1, 3, 2))
  dim(arr2) <- c(dim(arr)[[1]] * dim(arr)[[3]], dim(arr)[[2]])
  colsum2 <- Rfast::rowprods(arr2)
  dim(colsum2) <- dim(arr)[c(1, 3)]
  colsum <- colsum2
  dimnames(colsum) <- list(dimnames(arr)[[1]], dimnames(arr)[[3]])
}

rowprods_arr <- function(arr){
  arr2 <- arr #aperm(arr, c(1, 3, 2))
  dim(arr2) <- c(dim(arr)[[1]],  dim(arr)[[2]] * dim(arr)[[3]])
  colsum2 <- Rfast::colprods(arr2)
  dim(colsum2) <- dim(arr)[c(2, 3)]
  colsum <- colsum2
  dimnames(colsum) <- list(dimnames(arr)[[2]], dimnames(arr)[[3]])
  return(colsum)
}

#' @describeIn predsumspecies Computes expected numbers of species for a single parameter set and single ModelSite
#' @param Xocc A matrix of occupancy covariates. Must have a single row. Columns correspond to covariates.
#' @param Xobs A matrix of detection covariates, each row is a visit. If NULL then expected number of species in occupation is returned
#' @param theta A vector of model parameters, labelled according to the BUGS labelling convention seen in runjags
#' @param lv.v A matrix of LV values. Each column corresponds to a LV. To condition on specific LV values, provide a matrix of row 1.
#' @return A named vector of the expectation and variance of the numbers of species occupying the ModelSite and given parameter set.
#' If observational covariates are supplied then the expection and variance of numbers of species detected is also returned.
#' @export
# the following is really slow - too slow to use in a webapp
expectedspeciesnum.ModelSite.theta <- function(Xocc, Xobs = NULL, occ.b, det.b = NULL, lv.b = NULL, lv.v = NULL){
  ## Probability of Site Occupancy
  stopifnot(nrow(Xocc) == 1)
  if (is.null(Xobs)) {stopifnot(is.null(det.b))}
  Xocc <- as.matrix(Xocc)
  ModelSite.Occ.Pred.CondLV <- poccupy.ModelSite.theta(Xocc, occ.b, lv.b, lv.v)
  
  ## Expected number of species occupying modelsite, given model and theta, and marginal across LV
  EVnumspec_occ <- Erowsum_margrow(ModelSite.Occ.Pred.CondLV)
  names(EVnumspec_occ) <- paste0(names(EVnumspec_occ), "_occ")
  if (is.null(Xobs)){
    return(EVnumspec_occ)
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
    EVnumspec_det <- Erowsum_margrow(AnyDetections.Pred.marg)
    names(EVnumspec_det) <- paste0(names(EVnumspec_det), "_det")
    return(c(EVnumspec_occ, EVnumspec_det))
  }
}


# @param pmat A matrix of success probabilities for Bernoulli random variables.
# Each column corresponds to an independent Bernoulli random variable, each row is a different set of parameters.
# @return
# The expected sum of the Bernoulli random variables, marginal across rows, E[E[N | L]]
# The variance of the sum, marginal across rows E[N^2] - E[N]^2 = E[V[N | L] + E[N | L]^2] - E[E[N | L]]^2
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

#' @param UseFittedLV If TRUE the fitted LV variables are used, if false then 1000 LV values are simulated.
#' @param chains The chains of MCMC to use. Default is all chains.
#' @param cl A cluster object created by parallel::makeCluster. If NULL no cluster is used.
#' @param type If "draws" then predictions are given *per draw* and a 3-array is returned with dimensions of draws, sites and statistical summaries.
#' If "marginal" then the statistical summaries *marginalise* the posterior distribution and a matrix is returned with rows of statistical summaries and columns of sites.
#' There will be four summaries: the expection and variance of the number of species occupied or detected. 
#' If "median" then the expected number of species occupying and observed is returned for the median of the posterior distribution,
#' the variance and expected number of species, marginal over the posterior is also returned as they is useful for showing variation due to model parameter uncertainty.
#' @return A matrix or 3 array with each column a ModelSite. Dimensions are labelled. The statistical summaries returned are the predicted expection and variance of the number of species occupied or detected.
#' These expectations are with respect to the full posterior distribution of the model parameters, with the exception of the LV values which depends on UseFittedLV.
#' @export
predsumspecies <- function(fit, 
                           desiredspecies = fit$species,
                           chains = NULL,
                           UseFittedLV = TRUE,
                           nLVsim = 100,
                           type = "marginal",
                           cl = NULL){
  stopifnot(type == "marginal")
  warning("Obsolete function, modified to fit with new backend. nLVsim will be divided by 100.")
  if (!is.null(chains)){fit$mcmc <- fit$mcmc[chains]}
  specrich <- speciesrichness(fit,
                              occORdetection = "detection",
                              usefittedlvv = UseFittedLV,
                              desiredspecies = desiredspecies,
                              nlvperdraw = as.integer(nLVsim/100))
  return(specrich)
}

OLDpredsumspecies <- function(fit, desiredspecies = fit$species, chains = NULL, UseFittedLV = TRUE, nLVsim = 1000, type = "median", cl = NULL){
  stopifnot(type %in% c("draws", "median", "marginal"))
  stopifnot(all(desiredspecies %in% fit$species))
  
  fit$data <- as_list_format(fit$data)
  if (is.null(chains)){chains <- 1:length(fit$mcmc)}
  draws <- do.call(rbind, fit$mcmc[chains])
  if ( (is.null(fit$data$nlv)) || (fit$data$nlv == 0)){ #LVs not in model, add dummy variables
    LVbugs <- matrix2bugsvar(matrix(0, nrow = nrow(fit$data$Xocc), ncol = 2), "lv.v")
    LVbugs.draws <- Rfast::rep_row(LVbugs, nrow(draws))
    colnames(LVbugs.draws) <- names(LVbugs)
    
    lv.b.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$nspecies, ncol = 2), "lv.b")
    lv.b.draws <- Rfast::rep_row(lv.b.bugs, nrow(draws))
    colnames(lv.b.draws) <- names(lv.b.bugs)
    draws <- cbind(draws, lv.b.draws, LVbugs.draws)
    fit$data$nlv <- 2
    UseFittedLV <- TRUE #calculations faster when not simulating 1000s of LV values, especially since they are all ignored here.
  } 
  
  if (UseFittedLV){nLVsim <- NULL} #don't pass number of simulations if not going to use them
  
  numspec_drawsitesumm <- predsumspecies_raw(
    Xocc = fit$data$Xocc,
    Xobs = fit$data$Xobs,
    ModelSite = fit$data$ModelSite,
    numspeciesinmodel = fit$data$nspecies,
    desiredspecies = (1:fit$data$nspecies)[fit$species %in% desiredspecies],
    nlv = fit$data$nlv,
    draws = draws,
    useLVindraws = UseFittedLV,
    nLVsim = nLVsim,
    cl = cl
  )
  # convert predictions for each site and theta into predictions for each site, marginal across theta distribution
  if (type == "draws") {out <- numspec_drawsitesumm}
  if (type == "marginal") {out <- EVnumspec_marginalposterior_drawssitessummaries(numspec_drawsitesumm)}
  if (type == "median"){
    posterior_numspec <- EVnumspec_marginalposterior_drawssitessummaries(numspec_drawsitesumm)
    
    thetamedian <- apply(draws, MARGIN = 2, median)
    thetamedian <- matrix(thetamedian, nrow = 1, ncol = length(thetamedian),
                          dimnames = list(row = "median", cols = names(thetamedian)))
    numspec_site_median <- predsumspecies_raw(
      Xocc = fit$data$Xocc,
      Xobs = fit$data$Xobs,
      ModelSite = fit$data$ModelSite,
      numspecies = fit$data$nspecies,
      nlv = fit$data$nlv,
      draws = thetamedian,
      useLVindraws = UseFittedLV,
      nLVsim = nLVsim,
      cl = cl
    )
    out <- rbind(
      Esum_occ_median = numspec_site_median[1, , "Esum_occ"],
      Vsum_occ_median = numspec_site_median[1, , "Vsum_occ"],
      Esum_occ_margpost = posterior_numspec["Esum_occ", ],
      Vsum_occ_margpost = posterior_numspec["Vsum_occ", ],
      
      Esum_det_median = numspec_site_median[1, , "Esum_det"],
      Vsum_det_median = numspec_site_median[1, , "Vsum_det"],
      Esum_det_margpost = posterior_numspec["Esum_det", ],
      Vsum_det_margpost = posterior_numspec["Vsum_det", ]
    )
  }
  return(out)
}

#' @describeIn predsumspecies Called by predsumspecies and predsumspecies_newdata.
#' Given ModelSite data and model information, compute expected and variance of the number of species.
#' @param ModelSite is a list of integers giving the row in Xocc corresponding to a row in Xobs
#' @param Xocc A matrix of occupancy covariates, each row is a ModelSite
#' @param Xobs A matrix of detection covariates. Each row is a visit. The visited ModelSite (row of Xocc) is given by ModelSite
#' @param numspeciesinmodel Integer. The number of species in the model, which is needed to match up with the BUGS names in [draws]
#' @param desiredspecies Integer vector. The indexes for species to include in the species sum.
#' Useful if interested in only 1 species, or a certain class of species.
#' @param draws A matrix of posterior parameter draws. Each row is a draw. Column names follow the BUGS naming convention
#' @param useLVindraws Use the LV values corresponding to each draw from within the \code{draws} object.
#' If FALSE nLVsim simulated LV values will be used for each draw.
#' @param nLVsim The number of simulated LV values if not using fitted LV values (only applies if  useLVindraws = FALSE).
#' @details No scaling or centering of Xocc or Xobs is performed by predsumspecies_raw
#' @return A 3-dimensional array. Dimensions are draws, sites and statistical summaries of the number of species random variables.
#' There will be four summaries: the expection and variance of the number of species occupied or detected.
#' @export
predsumspecies_raw <- function(Xocc, Xobs = NULL, ModelSite = NULL,
                               numspeciesinmodel, desiredspecies = 1:numspeciesinmodel,
                               nlv, draws, useLVindraws = TRUE, nLVsim = NULL, cl = NULL){
  # prepare parameters
  nspmodel <- numspeciesinmodel
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
  occ.b <- bugsvar2array(draws, "occ.b", 1:nspmodel, 1:noccvar)[desiredspecies, , , drop = FALSE]
  if (!is.null(Xobs)) {det.b <- bugsvar2array(draws, "det.b", 1:nspmodel, 1:nobsvar)[desiredspecies, , , drop = FALSE]}
  lv.b <- bugsvar2array(draws, "lv.b", 1:nspmodel, 1:nlv)[desiredspecies, , , drop = FALSE]
  
  if (useLVindraws){stopifnot(is.null(nLVsim))}
  if (!useLVindraws){stopifnot(is.numeric(nLVsim))}

  
  sitedrawidxs <- expand.grid(siteidx = 1:nsites, drawidx = 1:ndraws) 
  
  if (!useLVindraws){ # predicting as if LVs not known, so simulate from their distribution
    lvsim <- matrix(rnorm(nlv * nLVsim), ncol = nlv, nrow = nLVsim)
  } else {
    lv.v <- bugsvar2array(draws, "lv.v", 1:nsites, 1:nlv)
  }
  

  # for each modelsite and each draw apply the following function:
  Enumspec <- pbapply::pblapply(1:nrow(sitedrawidxs),
        function(rowidx){
          sitedrawidx <- sitedrawidxs[rowidx, , drop = FALSE]
          Xocc <- Xocc[sitedrawidx[["siteidx"]], , drop = FALSE]
          if (!is.null(Xobs)) {Xobs <- Xobs[ModelSiteIdxs == sitedrawidx[["siteidx"]], , drop = FALSE]}
          occ.b_theta <- drop_to_matrix(occ.b[,, sitedrawidx[["drawidx"]], drop = FALSE])
          if (!is.null(Xobs)) {det.b_theta <- drop_to_matrix(det.b[,, sitedrawidx[["drawidx"]], drop = FALSE])}
          else {det.b_theta <- NULL}
          lv.b_theta <- drop_to_matrix(lv.b[,, sitedrawidx[["drawidx"]] , drop = FALSE])
          if (useLVindraws){
            lv.v_thetasite <- matrix(lv.v[sitedrawidx[["siteidx"]], , sitedrawidx[["drawidx"]], drop = FALSE], nrow = 1, ncol = nlv)
          } else {
            lv.v_thetasite <- lvsim
          }
          Enumspec_sitetheta <- expectedspeciesnum.ModelSite.theta(Xocc = Xocc, Xobs = Xobs,
                                                                   occ.b = occ.b_theta,
                                                                   det.b = det.b_theta,
                                                                   lv.b = lv.b_theta,
                                                                   lv.v = lv.v_thetasite)
          return(Enumspec_sitetheta)
        },
        cl = cl)
  Enumspec <- simplify2array(Enumspec)
  # each column of Enumspec is row of sitedrawidxs
  Enumspec <- rbind(Enumspec, t(sitedrawidxs))
  
  ## arrange Enumspec in a 3 dimensional array using existing arrangement of Enumspec
  a <- array(Enumspec, dim = c(nrow(Enumspec), nsites, ndraws),
             dimnames = list(
               Summaries = rownames(Enumspec),
               ModelSites = 1:nsites,
               Draws = 1:ndraws
             ))
  Enumspec_drawsitesumm <- aperm(a, perm = c("Draws", "ModelSites", "Summaries")) # arrange Enumspec in a 3 dimensional array: rows = draws, cols = site, depth = summaries
  
  #following checks that conversion is correct
  stopifnot(all(Rfast::eachrow(matrix(Enumspec_drawsitesumm[ , , "siteidx", drop = FALSE],
                                      nrow = dim(Enumspec_drawsitesumm)[[1]], ncol = dim(Enumspec_drawsitesumm)[[2]]), # call to matrix keeps the input in the correct dimension
                               y = as.numeric(1:nsites), oper = "==")))
  stopifnot(sum(Rfast::eachcol.apply(matrix(Enumspec_drawsitesumm[ , , "drawidx", drop = TRUE],
                                            nrow = dim(Enumspec_drawsitesumm)[[1]], ncol = dim(Enumspec_drawsitesumm)[[2]]),
                               y = as.numeric(1:ndraws), oper = "-")) == 0)
  
  return(Enumspec_drawsitesumm)
}

EVnumspec_marginalposterior_drawssitessummaries <- function(draws_sites_summaries){
  if ("Esum_det" %in% unlist(dimnames(draws_sites_summaries))) {detavailable <- TRUE}
  else {detavailable <- FALSE}
  out_occ <- EVtheta2EVmarg(
    Vsum = matrix(draws_sites_summaries[, ,"Vsum_occ", drop = FALSE],  nrow = dim(draws_sites_summaries)[[1]], ncol = dim(draws_sites_summaries)[[2]]),
    Esum = matrix(draws_sites_summaries[, , "Esum_occ", drop = FALSE], nrow = dim(draws_sites_summaries)[[1]], ncol = dim(draws_sites_summaries)[[2]])
  )
  rownames(out_occ) <- paste0(rownames(out_occ), "_occ")
  
  if (detavailable){
    out_det <- EVtheta2EVmarg(
      Vsum = matrix(draws_sites_summaries[, ,"Vsum_det", drop = FALSE],  nrow = dim(draws_sites_summaries)[[1]], ncol = dim(draws_sites_summaries)[[2]]),
      Esum = matrix(draws_sites_summaries[, , "Esum_det", drop = FALSE], nrow = dim(draws_sites_summaries)[[1]], ncol = dim(draws_sites_summaries)[[2]])
    )
    rownames(out_det) <- paste0(rownames(out_det), "_det")
    out <- rbind(
      occ = out_occ,
      det = out_det)
  } else {
    out <- out_occ
  }
  return(out)
}

# convert predictions for each site and theta into median prediction of expected Enumspec for each site and a confidence interval for Enumspec
# @param probs are the quantile probabilities to return 
Enumspec_quantiles_drawssitessummaries <- function(draws_sites_summaries, probs = seq(0, 1, 0.25)){
  Esum_occ = matrix(draws_sites_summaries[, , "Esum_occ", drop = FALSE], nrow = dim(draws_sites_summaries)[[1]], ncol = dim(draws_sites_summaries)[[2]])
  Esum_occ_quants <- apply(Esum_occ, MARGIN = 1, function(x) quantile(x, probs = probs))
  rownames(Esum_occ_quants) <- paste0("Esum_occ_", rownames(Esum_occ_quants))
  
  Esum_det = matrix(draws_sites_summaries[, , "Esum_det", drop = FALSE], nrow = dim(draws_sites_summaries)[[1]], ncol = dim(draws_sites_summaries)[[2]])
  Esum_det_quants <- apply(Esum_det, MARGIN = 1, function(x) quantile(x, probs = probs))
  rownames(Esum_det_quants) <- paste0("Esum_det_", rownames(Esum_det_quants))
  return(rbind(Esum_occ_quants, Esum_det_quants))
}

Vnumspec_quantiles_drawssitessummaries <- function(draws_sites_summaries, probs = seq(0, 1, 0.25)){
  Vsum_occ = matrix(draws_sites_summaries[, , "Vsum_occ", drop = FALSE], nrow = dim(draws_sites_summaries)[[1]], ncol = dim(draws_sites_summaries)[[2]])
  Vsum_occ_quants <- apply(Vsum_occ, MARGIN = 1, function(x) quantile(x, probs = probs))
  rownames(Vsum_occ_quants) <- paste0("Vsum_occ_", rownames(Vsum_occ_quants))
  
  Vsum_det = matrix(draws_sites_summaries[, , "Vsum_det", drop = FALSE], nrow = dim(draws_sites_summaries)[[1]], ncol = dim(draws_sites_summaries)[[2]])
  Vsum_det_quants <- apply(Vsum_det, MARGIN = 1, function(x) quantile(x, probs = probs))
  rownames(Vsum_det_quants) <- paste0("Vsum_det_", rownames(Vsum_det_quants))
  return(rbind(Vsum_occ_quants, Vsum_det_quants))
}

# convert predictions for each site and theta into predictions for each site, marginal across theta distribution
# @param draws_sites_summaries is 3-dimensional array. First dimension is draws, second dimension is sites, third dimension is summarieis
# approximate 95% posterior density interval for sum of species detected using Gaussian approximation and variance.
numspec_posteriorinterval_Gaussian_approx <- function(draws_sites_summaries){
  moments <- EVnumspec_marginalposterior_drawssitessummaries(draws_sites_summaries)
  sum_occ_low <- moments["Esum_occ", ] - 2 * sqrt(moments["Vsum_occ", ])
  sum_occ_high <- moments["Esum_occ", ] + 2 * sqrt(moments["Vsum_occ", ])
  sum_det_low <- moments["Esum_det", ] - 2 * sqrt(moments["Vsum_det", ])
  sum_det_high <- moments["Esum_det", ] + 2 * sqrt(moments["Vsum_det", ])
  out <- rbind(
    sum_occ_low = sum_occ_low,
    sum_occ_pred = moments["Esum_occ", ],
    sum_occ_high = sum_occ_high,
    sum_det_low = sum_det_low,
    sum_det_pred = moments["Esum_det", ],
    sum_det_high = sum_det_high
  )
  return(out)
}

# @param Vsum A matrix of variance of sum of species. Each row corresponds to a different theta. Each column a ModelSite.
# @param Esum A matrix of expected sum of species. Each row corresponds to a different theta. Each column a ModelSite.
# @return The expectation and variance of the sum of species marginal across theta (across the supplied rows)
EVtheta2EVmarg <- function(Vsum, Esum){
  stopifnot(all(dim(Vsum) == dim(Esum)))
  # per draw the 2nd moments
  M2n_theta <- Vsum + Esum^2
  
  # marginal across draw moments
  M2n <- Rfast::colmeans(M2n_theta)
  En <- Rfast::colmeans(Esum)
  Vn <- M2n - En^2
  return(rbind(
    Esum = En,
    Vsum = Vn
  ))
}
