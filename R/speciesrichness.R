#' @title Expected Biodiversity
#' @description The expected number of species occupying a ModelSite, or the expected number of species detected at a ModelSite.
#' 
#' @examples 
#' fit <- readRDS("./tmpdata/7_2_9_addyear_msnm_year_time_2lv.rds")
#' theta <- get_theta(fit, type = 1)
#' Xocc <- fit$data$Xocc[2, , drop = FALSE]
#' Xobs <- fit$data$Xobs[fit$data$ModelSite == 2, , drop = FALSE]
#' numspecies <- fit$data$nspecies
#' lvsim <- matrix(rnorm(2 * 1), ncol = 2, nrow = 2) #dummy lvsim vars
#' lv.b.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$nspecies, ncol = 2), "lv.b")
#' theta <- c(theta, lv.b.bugs)
#' 
#' Enumspec <- predsumspecies(fit, UseFittedLV = TRUE, return = "median")
#' 
#' 
#' indata <- readRDS("./private/data/clean/7_2_10_input_data.rds")
#' predsumspecies_newdata(fit, Xocc = indata$holdoutdata$Xocc, Xobs = indata$holdoutdata$yXobs, ModelSiteVars = "ModelSiteID", return = "median", cl = NULL)
#' draws_sites_summaries <- predsumspecies_newdata(fit, 
#'                                                 Xocc = indata$holdoutdata$Xocc,
#'                                                 Xobs = indata$holdoutdata$yXobs,
#'                                                 ModelSiteVars = "ModelSiteID",
#'                                                 bydraw = TRUE,
#'                                                 cl = NULL)
#' # Median expected biodiversity with 95% credible intervals for expected biodiversity
#' Enumspec_quantiles_drawssitessummaries(draws_sites_summaries, probs  = c(0.025, 0.5, 0.0975))
#' 
#' # approximate 95% posterior density interval for sum of species detected using Gaussian approximation and variance.
#' numspec_interval <- numspec_posteriorinterval_Gaussian_approx(draws_sites_summaries)

# parameters the same as poccupy_raw
# returns a matrix with each column a site. Rows of expectation and variance.
OLDoccspeciesrichness_raw.jsodm_lv <- function(fixedcovar, loadfixed, randomcovar, loadrandom){
  pocc <- poccupy_raw.jsodm_lv(fixedcovar, loadfixed, randomcovar, loadrandom)
  Evals <- prednumsuccess(pocc)
  return(Evals)
}
occspeciesrichness_raw.jsodm <- function(fixedcovar, loadfixed){
  pocc <- poccupy_raw.jsodm_lv(fixedcovar, loadfixed)
  Evals <- prednumsuccess(pocc)
  return(Evals)
}

# probarr is a 3-array of probabilities of sucess Bernoulli RV. Within each layer columns (species) are independent (i.e. conditional on draw)
# computes the Expectation and Variance of the number of success per site assuming each layer (draw)
# is a sample from the true distribution of model parameters
prednumsuccess <- function(probarr){
  # En_sitedraw <- apply(probarr, MARGIN = c(1, 3), sum)
  # En_site <- apply(En_sitedraw, MARGIN = 1, mean)
  En_site <- apply(probarr, MARGIN = 1, sum) / dim(probarr)[[3]] #equivalent to two step process above
  # similarly V can be summed, then dived by number of draws squared
  Vrv <- probarr * (1 - probarr)
  # use total law of variance
  Edrawvariance <- apply(Vrv, MARGIN = 1, sum) / dim(probarr)[[3]]
  # Vdrawexpectation <- apply(apply(probarr, MARGIN = c(1, 3), sum)^2,
  #                           MARGIN = 1, mean) - En_site^2  #warning this is a biased estimate of variance: better would be to use the var function below
  Vdrawexpectation <- Rfast::rowVars(arr3_sumalong2(probarr))
  Vn_site <- Edrawvariance + Vdrawexpectation
  Evals <- rbind(E = En_site,
                 V = Vn_site)
  return(Evals)
}

#returns the variance of the "expected number of species", not the variance of the number of species
Eoccspeciesrichness_raw.jsodm_lv <- function(fixedcovar, loadfixed, randomcovar, loadrandom){
  pocc <- poccupy_raw.jsodm_lv(fixedcovar, loadfixed, randomcovar, loadrandom)
  EEn_site <- apply(pocc, MARGIN = 1, sum) / dim(pocc)[[3]] 
  En_sitedraw <- arr3_sumalong2(pocc)
  # VEn_site <- apply(En_sitedraw, MARGIN = 1, var)
  VEn_site <- Rfast::rowVars(En_sitedraw)
  vals <- rbind(E = EEn_site,
                 V = VEn_site)
  return(vals)
}

# pocc = 3 array of probability of occupancy. Within each draw occupancy of species considered independent.
# pdet_occ = 3 array of probability of detection, condition on occupancy. Within each draw species detection considered independent given occupied.
detspeciesrichness_probarr <- function(pocc, pdet_occ, ModelSite){
  ModSiteVisits <- plyr::split_indices(ModelSite)
  
  NoDetections_occ <- vapply(ModSiteVisits,
                             function(visits){
                               NoDetectProb <- rowprods_arr(1 - pdet_occ[visits, , , drop = FALSE])
                               return(NoDetectProb)
                             },
                             FUN.VALUE = matrix(1.1, nrow = ncol(pdet_occ), ncol = dim(pdet_occ)[[3]]))
  #NoDetections_occ array of species x draw x ModelSite
  NoDetections_occ <- aperm(NoDetections_occ, perm = c(3, 1, 2))

  ## probability of no detections for all visits to a site, marginal on occupancy
  NoDetections <- NoDetections_occ * pocc + 1 - pocc
  pdet <- 1 - NoDetections
  Evals <- prednumsuccess(pdet)
  return(Evals)
}

OLDdetspeciesrichness_raw.jsodm_lv <- function(Xocc, occ.b, Xobs, det.b, ModelSite, lv.v, lv.b){
  pocc <- poccupy_raw.jsodm_lv(Xocc, occ.b, lv.v, lv.b)

  ## Probability of Detection
  pdet_occ <- pdet_occ_raw.jsodm(Xobs, det.b)
  Evals <- detspeciesrichness_probarr(pocc, pdet_occ, ModelSite)
  return(Evals)
}

detspeciesrichness_raw.jsodm <- function(Xocc, occ.b, Xobs, det.b, ModelSite){
  pocc <- poccupy_raw.jsodm(Xocc, occ.b)
  
  ## Probability of Detection
  pdet_occ <- pdet_occ_raw.jsodm(Xobs, det.b)
  Evals <- detspeciesrichness_probarr(pocc, pdet_occ, ModelSite)
  return(Evals)
}

#' @export
speciesrichness <- function(fit, occORdetection, ...){
  UseMethod("speciesrichness")
}

#' @export
speciesrichness.jsodm_lv <- function(fit, 
                                     occORdetection,
                                    desiredspecies = fit$species,
                                    nlvperdraw = 1){
  stopifnot(all(desiredspecies %in% fit$species))
  occ.v <- fit$data$Xocc
  occ.b <- get_occ_b(fit)[desiredspecies, , , drop = FALSE]
  lv.b <- get_lv_b(fit)[desiredspecies, , , drop = FALSE]
  if (nlvperdraw > 1){
    occ.bs <- lapply(1:nlvperdraw, function(x){occ.b})
    occ.b <- abind::abind(occ.bs, along = 3)
    
    lv.bs <- lapply(1:nlvperdraw, function(x){lv.b})
    lv.b <- abind::abind(lv.bs, along = 3)
  }
  lv.v <- array(rnorm(dim(occ.v)[[1]] * dim(lv.b)[[2]] *  dim(lv.b)[[3]]), 
                dim = c(dim(occ.v)[[1]], dim(lv.b)[[2]],  dim(lv.b)[[3]]),
                dimnames = list(ModelSite = rownames(occ.v),
                                LV = paste0("lv", 1:dim(lv.b)[[2]], ".v"),
                                Draw = 1:dim(lv.b)[[3]]))
  
  pocc <- poccupy_raw.jsodm_lv(occ.v, occ.b, lv.v, lv.b)
  if (occORdetection == "detection"){
    det.v <- fit$data$Xobs
    det.b <- get_det_b(fit)[desiredspecies, , , drop = FALSE]
    if (nlvperdraw > 1){    
      det.bs <- lapply(1:nlvperdraw, function(x){det.b})
      det.b <- abind::abind(det.bs, along = 3) 
    }
    
    pdet_occ <- pdet_occ_raw.jsodm(det.v, det.b)
    specrich  <- detspeciesrichness_probarr(pocc, pdet_occ, fit$data$ModelSite)
  } else if (occORdetection == "occupancy") {
    specrich <- prednumsuccess(pocc)
  }
  return(specrich)
}

#' @export
speciesrichness.jsodm <- function(fit, 
                                     occORdetection,
                                     desiredspecies = fit$species){
  stopifnot(all(desiredspecies %in% fit$species))
  occ.v <- fit$data$Xocc
  occ.b <- get_occ_b(fit)[desiredspecies, , , drop = FALSE]

  pocc <- poccupy_raw.jsodm(occ.v, occ.b)
  if (occORdetection == "detection"){
    det.b <- get_det_b(fit)[desiredspecies, , , drop = FALSE]
    det.v <- fit$data$Xobs

    pdet_occ <- pdet_occ_raw.jsodm(det.v, det.b)
    specrich  <- detspeciesrichness_probarr(pocc, pdet_occ, fit$data$ModelSite)
  } else if (occORdetection == "occupancy") {
    specrich <- prednumsuccess(pocc)
  }
  return(specrich)
}


#nlvperdraw = 1 by default.
speciesrichness_newdata.jsodm_lv <- function(fit, Xocc, Xobs = NULL, ModelSite = NULL,
                                    desiredspecies = fit$species,
                                    nlvperdraw = 1){
  fit <- supplant_new_data(fit, Xocc, Xobs, ModelSite = ModelSite)
  occORdetection <- switch(
    as.character(is.null(Xobs)),
    `FALSE` = "occupancy",
    `TRUE` = "detection"
  )
  specrich <- speciesrichness.jsodm_lv(fit, 
                                       occORdetection = occORdetection,
                                       desiredspecies = desiredspecies, 
                                       nlvperdraw = nlvperdraw)
  return(specrich)
}



#' @describeIn predsumspecies For new ModelSite occupancy covariates and detection covariates, predicted number of expected species
#' @param desiredspecies List of species to sum over. Names must match names in fit$species. Default is to sum over all species.
#' @export
predsumspecies_newdata <- function(fit, Xocc, Xobs = NULL, ModelSiteVars = NULL,
                                   desiredspecies = fit$species,
                                   chains = NULL, nLVsim = 1000, type = "marginal", cl = NULL){
  if (type == "marginal" && ("jsodm_lv" %in% class(fit))){
    out <- speciesrichness_newdata.jsodm_lv(fit, Xocc, Xobs,
                                            ModelSiteVars,  desiredspecies, 
                                            nlvperdraw = as.integer(nLVsim/100))
    return(out)
  }
  
  warning("Using old code to predict species richness")
  
  stopifnot(type %in% c("draws", "marginal", "median"))
  stopifnot(all(desiredspecies %in% fit$species))
  
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
    lv.b.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$nspecies, ncol = 2), "lv.b")
    lv.b.draws <- Rfast::rep_row(lv.b.bugs, nrow(draws))
    colnames(lv.b.draws) <- names(lv.b.bugs)
    draws <- cbind(draws, lv.b.draws)
    fit$data$nlv <- 2
    nLVsim = 2 #calculations faster when not simulating 1000s of LV values, especially since they are all ignored here.
  }
  
  numspec_drawsitesumm <- predsumspecies_raw(
    Xocc = datalist$Xocc,
    Xobs = datalist$Xobs,
    ModelSite = datalist$ModelSite,
    numspeciesinmodel = fit$data$nspecies,
    desiredspecies = (1:fit$data$nspecies)[fit$species %in% desiredspecies],
    nlv = fit$data$nlv,
    draws = draws,
    useLVindraws = FALSE,
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
      Xocc = datalist$Xocc,
      Xobs = datalist$Xobs,
      ModelSite = datalist$ModelSite,
      numspecies = fit$data$nspecies,
      nlv = fit$data$nlv,
      draws = thetamedian,
      useLVindraws = FALSE,
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

#' @title The number of observed species in a matrix of observation recordings
#' @param y A matrix of species *observations* with each row a visit and each column a species. Entries must be either 0 or 1.
#' @param ModelSite The list of ModelSite indexes corresponding to each row in y
#' @return A vector of the number of species detected at each ModelSite. Names give the ModelSite index. 
#' @export
detectednumspec <- function(y, ModelSite){
  stopifnot(length(ModelSite) == nrow(y))
  stopifnot(all(as.matrix(y) %in% c(1, 0)))
  
  my <- cbind(ModelSite = ModelSite, y)
  SpDetected <- my %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise_all(~sum(.) > 0)
  NumSpecies <- as.vector(rowSums(SpDetected[, -1]))
  names(NumSpecies) <- SpDetected[, 1, drop = TRUE]
  return(NumSpecies)
}
