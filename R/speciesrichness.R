#' @title Expected Biodiversity
#' @description The expected number of species occupying a ModelSite, or the expected number of species detected at a ModelSite.
#' @examples 
#' fit <- readRDS("../sflddata/private/data/testdata/cutfit_7_4_11_2LV.rds")
#' fit <- translatefit(fit)
# cl <- NULL #parallel::makeCluster(2)
# srich <- speciesrichness(fit,  occORdetection = "detection",
#                 desiredspecies = fit$species, usefittedlvv = FALSE,
#                                      nlvperdraw = 1000, chunksize = 2, cl = cl)

# probarr is a 3-array of probabilities of sucess Bernoulli RV. Within each layer columns (species) are independent (i.e. conditional on draw)
# computes the Expectation and Variance of the number of success per site assuming each layer (draw)
# is a sample from the true distribution of model parameters
prednumsuccess <- function(probarr){
  Vrv <- probarr * (1 - probarr)
  Evals <- prednumsuccess_ErvVrv(probarr, Vrv)
  return(Evals)
}

# Erv and Vrv are 3-arrays of the mean and variance, respectively, of RV. Within each layer columns (species) are independent (i.e. conditional on draw)
# computes the Expectation and Variance of the number of success per site assuming each layer (draw)
# is a sample from the true distribution of model parameters
# returns the expected sum for each row, along with the variance of this sum.
prednumsuccess_ErvVrv <- function(Erv, Vrv){
  stopifnot(dim(Vrv)[[3]] == dim(Erv)[[3]])
  # Esum_sitedraw <- apply(probarr, MARGIN = c(1, 3), sum)
  # Esum_site <- apply(En_sitedraw, MARGIN = 1, mean)
  Esum_site <- apply(Erv, MARGIN = 1, sum) / dim(Erv)[[3]] #equivalent to two step process above
  # use total law of variance
  Edrawvariance <- apply(Vrv, MARGIN = 1, sum) / dim(Vrv)[[3]]
  # Vdrawexpectation <- apply(apply(probarr, MARGIN = c(1, 3), sum)^2,
  #                           MARGIN = 1, mean) - En_site^2  #warning this is a biased estimate of variance: better would be to use the var function below
  if (dim(Vrv)[[3]] == 1){Vdrawexpectation <- 0 * Edrawvariance} #only one draw so variance between draws is 0
  else {Vdrawexpectation <- Rfast::rowVars(arr3_sumalong2(Erv))} #this is sample variance estimate, multiply by (dim(Vrv)[[3]] - 1) / dim(Vrv)[[3]] to get old variance numbers
  Vsum_site <- Edrawvariance + Vdrawexpectation
  Evals <- rbind(E = Esum_site,
                 V = Vsum_site)
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
  if (is.null(dim(NoDetections_occ))) {
    dim(NoDetections_occ) <- dim(pocc)[c(2, 3, 1)]
    dimnames(NoDetections_occ) <- dimnames(pocc)[c(2, 3, 1)]
  }
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
#' @param occORdetection 'detection' to predict detected richness, or 'occupancy' to predict richness of occupancy species
#' @param desiredspecies A list of species to include in the calculations, defaults to all species in the model.
speciesrichness <- function(fit, occORdetection, desiredspecies, ...){
  UseMethod("speciesrichness")
}

#' @export
speciesrichness.jsodm_lv <- function(fit, 
                                     occORdetection,
                                    desiredspecies = fit$species,
                                     usefittedlvv = FALSE,
                                    nlvperdraw = 1,
                                    chunksize = ceiling(1E7 / (nrow(fit$mcmc[[1]]) * nlvperdraw)),
                                    cl = NULL){
  stopifnot(chunksize > 1)
  if ("cluster" %in% class(cl)){
    parallel::clusterExport(cl = cl,
                            varlist = c("fit", "occORdetection", "desiredspecies",
                                        "usefittedlvv", "nlvperdraw"),
                            envir = environment())
  }
  sitesplits <- split(1:nrow(fit$data$Xocc), ceiling(1:nrow(fit$data$Xocc) / chunksize))
  if ("cluster" %in% class(cl)){
    parallel::clusterExport(cl = cl,
                            varlist = c("fit", "occORdetection", "desiredspecies",
                                        "usefittedlvv", "nlvperdraw"),
                            envir = environment())
    srich.l <- parallel::parLapply(cl = cl , sitesplits, function(modelsites){
      subfit <- subsetofmodelsites.jsodm_lv(fit, modelsites)
      srich <- speciesrichness.jsodm_lv_allsites(subfit, occORdetection = occORdetection, 
                                                 desiredspecies = desiredspecies, usefittedlvv = usefittedlvv,
                                                 nlvperdraw = nlvperdraw)
    }, chunk.size = 1)
    srich <- do.call(cbind, srich.l)
  } else {
    srich <- NULL
    for (modelsites in sitesplits){
      subfit <- subsetofmodelsites.jsodm_lv(fit, modelsites)
      srich.t <- speciesrichness.jsodm_lv_allsites(subfit, occORdetection = occORdetection, 
                                                 desiredspecies = desiredspecies, usefittedlvv = usefittedlvv,
                                                 nlvperdraw = nlvperdraw)
      srich <- cbind(srich, srich.t)
    }
  }
  return(srich)
}

speciesrichness.jsodm_lv_allsites <- function(fit, 
                                     occORdetection,
                                    desiredspecies = fit$species,
                                     usefittedlvv = FALSE,
                                    nlvperdraw = 1){
  stopifnot(all(desiredspecies %in% fit$species))
  occ.v <- fit$data$Xocc
  occ.b <- get_occ_b(fit)[desiredspecies, , , drop = FALSE]
  lv.b <- get_lv_b(fit)[desiredspecies, , , drop = FALSE]
  if (sum(usefittedlvv, nlvperdraw > 1) > 1){
    stop("If using fitted lv.v then can not simulate multiple lv.v per draw.")
  }
  if (nlvperdraw > 1){
    occ.bs <- lapply(1:nlvperdraw, function(x){occ.b})
    occ.b <- abind::abind(occ.bs, along = 3)
    
    lv.bs <- lapply(1:nlvperdraw, function(x){lv.b})
    lv.b <- abind::abind(lv.bs, along = 3)
  }
  if (usefittedlvv){
    lv.v <- get_lv_v(fit)
  } else {
    lv.v <- array(rnorm(dim(occ.v)[[1]] * dim(lv.b)[[2]] *  dim(lv.b)[[3]]), 
                  dim = c(dim(occ.v)[[1]], dim(lv.b)[[2]],  dim(lv.b)[[3]]),
                  dimnames = list(ModelSite = rownames(occ.v),
                                  LV = paste0("lv", 1:dim(lv.b)[[2]], ".v"),
                                  Draw = 1:dim(lv.b)[[3]]))
  }
  
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
