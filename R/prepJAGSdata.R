#' @title Prepare JAGS-like Data for JSODM
#' @details Given the input data parameters of run.detectionoccupancy prepare the data list for JAGS
#' @param Xocc A dataframe of covariates related to occupancy. One row per ModelSite.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param yXobs A dataframe of species observations (1 or 0) and covariates related to observations. One row per visit.
#' Each column is either a covariate or a species.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param species A list of species names (corresponding to columns of yXobs) to model.
#' @param OccFmla A formula specifying the occupancy model in terms of the covariates in Xocc.
#' @param ObsFmla A formula specifying the model for detection in terms of the covariates in Xobs.
#' @param ModelSite A list of column names in y, Xocc and Xobs that uniquely specify the ModelSite. Can be simply a ModelSite index
#' @param nlv The number of latent variables
#' @param XoccProcess An object create by prep.designmatprocess for the occupancy covariates
#' @param XobsProcess An object create by prep.designmatprocess for the observation covariates
#' @export
prepJAGSdata <- function(modeltype, Xocc, yXobs, ModelSite, species, XoccProcess, XobsProcess, ...){ # nlv){
  stopifnot(modeltype %in% availmodeltypes)
  
  data.list <- do.call(paste0("prepJAGSdata.", modeltype),
                       c( list(Xocc, yXobs, ModelSite, species, XoccProcess, XobsProcess), ...))
  return(data.list)
}

#' @describeIn prepJAGSdata Data preparation for the jsodm model.
prepJAGSdata.jsodm <- function(Xocc, yXobs, ModelSite, species, XoccProcess, XobsProcess, ...){
  # check data inputs
  stopifnot(all(ModelSite %in% colnames(Xocc)))
  stopifnot(all(ModelSite %in% colnames(yXobs)))
  stopifnot(anyDuplicated(Xocc[, ModelSite]) == 0) #model site uniquely specified
  if (all(species %in% colnames(yXobs))) {stopifnot(all(as.matrix(yXobs[, species]) %in% c(1, 0)))}
  
  # create model site indexes
  # if (ModelSiteID %in% c(names(yXobs), Xocc)){warning("Overwriting ModelSiteID column in input data.")}
  ModelSiteMat <- cbind(1:nrow(Xocc), tibble::as_tibble(Xocc[, ModelSite, drop = FALSE]))
  visitedModelSiteMat <- dplyr::right_join(ModelSiteMat, tibble::as_tibble(yXobs[, ModelSite, drop = FALSE]), by = ModelSite, suffix = c("", ".in"))
  visitedModelSite <- visitedModelSiteMat[, 1, drop = TRUE]
  stopifnot(is.integer(visitedModelSite))
  stopifnot(all(visitedModelSite <= nrow(Xocc)))
  
  XoccDesign <- apply.designmatprocess(XoccProcess, Xocc)
  XobsDesign <- apply.designmatprocess(XobsProcess, yXobs)
  # rownames(XobsDesign) <- visitedModelSite
  
  n = length(species) #number of species
  J <- nrow(XoccDesign)  #number of unique sites should also be max(occ_covariates$SiteID)
  if (all(species %in% colnames(yXobs))) { #if this is true the y is part of yXobs
    y <- as.matrix(yXobs[, species])
  } else { #if not then situation of prepping data of new ModelSites
    y <- NULL
  }
  ModelSite <- visitedModelSite
  data.list = list(n=n, J=J, y=y,
                  ModelSite = ModelSite, #a list of the site visited at each visit
                  nvisits = nrow(XobsDesign), #number of visits in total - not sure what this is for
                  Xocc=XoccDesign,Xobs=XobsDesign,noccvar=ncol(XoccDesign),nobsvar=ncol(XobsDesign))
  return(data.list)
}

#' @describeIn prepJAGSdata Data preparation for the jsodm_lv model
prepJAGSdata.jsodm_lv <- function(Xocc, yXobs, ModelSite, species, XoccProcess, XobsProcess, nlv){
  data.list <- prepJAGSdata.jsodm(Xocc, yXobs, ModelSite, species, XoccProcess, XobsProcess)
  data.list <- c(data.list, nlv = nlv)
  return(data.list)
}

#' @describeIn prepJAGSdata Data preparation for the jsodm_lv_sepexp model (latent variables with separable covariance that is an exponential function)
prepJAGSdata.jsodm_lv_sepexp <- function(Xocc, yXobs, ModelSite, species, XoccProcess, XobsProcess, nlv, SpatDist, TimeDist, ...){
  data.list <- prepJAGSdata.jsodm_lv(Xocc, yXobs, ModelSite, species, XoccProcess, XobsProcess, nlv)
  spatdistmat <- SpatDist(Xocc)
  stopifnot("matrix" %in% class(spatdistmat))
  timedistmat <- TimeDist(Xocc)
  stopifnot("matrix" %in% class(timedistmat))
  warning("distances are not standardised to have mean 1, and sd of 1 (so prior of covariance scale may be inappropriate)")
  stopifnot((ncol(spatdistmat) == nrow(data.list$Xocc)) && (nrow(spatdistmat) == nrow(data.list$Xocc)))
  stopifnot((ncol(timedistmat) == nrow(data.list$Xocc)) && (nrow(timedistmat) == nrow(data.list$Xocc)))
  zero.lv.v <- rep(1, nrow(data.list$Xocc))
  warning("prior means of LV set to 1")
  data.list <- c(data.list, list(spatdistmat = spatdistmat, timedistmat = timedistmat, zero.lv.v = zero.lv.v))
  return(data.list)
}

#' @describeIn prepJAGSdata A short function that applies the prepJAGSdata function to new data, given an object created by run.detectionoccupancy
#' Xocc, yXobs, ModelSite must follow some rules as for run.detectionoccupancy, except yXobs may omit the species detections
#' @export
prep_new_data <- function(fit, Xocc, yXobs, ModelSite, ...){
  data.list <- do.call(prepJAGSdata, 
          c(list(class(fit)[[1]], Xocc, yXobs, ModelSite, fit$species, fit$XoccProcess, fit$XobsProcess), fit))
  # data.list <- prepJAGSdata(class(fit)[[1]], Xocc, yXobs, ModelSite, fit$species, fit$XoccProcess, fit$XobsProcess, fit$nlv)
  return(data.list)
}
