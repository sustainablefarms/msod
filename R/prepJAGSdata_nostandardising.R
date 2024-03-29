#' @title Prepare JAGS-like Data for JSODM
#' @details Given the input data parameters of run.detectionoccupancy prepare the data list for JAGS
#' @param Xocc A dataframe of covariates related to occupancy. One row per ModelSite.
#' @param y A dataframe of species observations (1 or 0). One row per visit. Each column is a species.
#' @param Xobs A dataframe of covariates related to observations. Each column is a covariate.
#' @param ModelSite A list of integers giving the number each row of Xocc (ModelSite) visited for each visit (row of y and Xobs)
#' @param nlv The number of latent variables
#' @export
# improve: make Xocc the MODEL matrix (DESIGN matrix), Xobs also the ModelMatrix, y the response matrix and ModelSite visited row of Xocc
prepJAGSdata2 <- function(modeltype, Xocc, Xobs, y, ModelSite, ...){ # nlv){
  stopifnot(modeltype %in% availmodeltypes)
  ModelSite <- as.integer(ModelSite)
  data.list <- do.call(paste0("prepJAGSdata2.", modeltype),
                       c( list(Xocc, Xobs, y, ModelSite), list(...)))
  return(data.list)
}

#' @describeIn prepJAGSdata Data preparation for the jsodm model.
prepJAGSdata2.jsodm <- function(Xocc, Xobs, y, ModelSite, ...){
  # do the prep that doesn't involve detections, so can exit early if detection components not included
  if (is.null(Xobs)){ #exit early if no Xobs included
    return(list(nspecies=ncol(y),
                nmodelsites = nrow(Xocc),
                Xocc=Xocc,
                noccvar=ncol(Xocc)))
  }
  
  stopifnot(setequal(ModelSite, 1:nrow(Xocc)))
  stopifnot(!anyNA(Xocc))
  stopifnot(!anyNA(Xobs))
  stopifnot("integer" %in% class(ModelSite))
  stopifnot(nrow(Xobs) == length(ModelSite))
  stopifnot(!anyNA(ModelSite))
  stopifnot(is.numeric(Xocc))
  stopifnot(is.numeric(Xobs))
  stopifnot("matrix" %in% class(Xocc))
  stopifnot("matrix" %in% class(Xobs))
  
  if (!is.null(y)){ #do fewer checks if y isn't passed - means function is altering built in data
    stopifnot(nrow(Xobs) == nrow(y))
    stopifnot(nrow(y) == length(ModelSite))
    stopifnot(!anyNA(y))
    stopifnot(is.numeric(y))
    stopifnot("matrix" %in% class(y))
    stopifnot(all(y %in% c(1, 0)))
    if (is.null(colnames(y))){
      warning("y doesn't have column names. The species names are derived from the alphabet.")
      colnames(y) <- make.names(rep(LETTERS, ceiling(ncol(y) / length(LETTERS))), unique = TRUE)[1:ncol(y)]
      }
  }
  
  
  
  nspecies = ncol(y) #number of species
  nmodelsites <- nrow(Xocc)  #number of unique sites should also be max(occ_covariates$SiteID)
  
  data.list = list(nspecies=nspecies, nmodelsites=nmodelsites, y=y,
                  ModelSite = ModelSite, #a list of the site visited at each visit
                  nvisits = nrow(Xobs), #number of visits in total - not sure what this is for
                  Xocc=Xocc,Xobs=Xobs,noccvar=ncol(Xocc),nobsvar=ncol(Xobs))
  return(data.list)
}

#' @describeIn prepJAGSdata Data preparation for the jsodm_lv model
prepJAGSdata2.jsodm_lv <- function(Xocc, Xobs, y, ModelSite, nlv, ...){
  data.list <- prepJAGSdata2.jsodm(Xocc, Xobs, y, ModelSite)
  data.list <- c(data.list, nlv = nlv)
  return(data.list)
}

#' @describeIn prepJAGSdata Data preparation for the jsodm_lv_re model
#' @param ModelSiteGroup is a vector of integers for grouping model sites for the purpose of random effects. For each row of Xocc it gives an index that defines the group. Each group has its own random effect for occupancy, and another for detection
# ModelSiteGroup <- as.integer(as.factor(Xocc$SiteCode) )
prepJAGSdata2.jsodm_lv_re <- function(Xocc, Xobs, y, ModelSite, nlv, ModelSiteGroup, ...){
  stopifnot(nrow(Xocc) == length(ModelSiteGroup))
  stopifnot(all(1:max(ModelSiteGroup) %in% ModelSiteGroup))
  data.list <- prepJAGSdata2.jsodm_lv(Xocc, Xobs, y, ModelSite, nlv)
  SiteOccInd <- ModelSiteGroup # a map from model site to group
  SiteObsInd <- ModelSiteGroup[ModelSite]# a map from visit to group 
  nsitegroups <- max(ModelSiteGroup) # number of different groups
  data.list <- c(data.list, list(nsitegroups = nsitegroups, SiteOccInd = SiteOccInd, SiteObsInd = SiteObsInd))
  return(data.list)
}

#' @describeIn prepJAGSdata Data preparation for the jsodm_lv_sepexp model (latent variables with separable covariance that is an exponential function)
prepJAGSdata2.jsodm_lv_sepexp <- function(Xocc, Xobs, y, ModelSite, nlv, spatdistmat, timedistmat, ...){
  data.list <- prepJAGSdata2.jsodm_lv(Xocc, Xobs, y, ModelSite, nlv)
  stopifnot("matrix" %in% class(spatdistmat))
  stopifnot("matrix" %in% class(timedistmat))
  stopifnot((ncol(spatdistmat) == nrow(data.list$Xocc)) && (nrow(spatdistmat) == nrow(data.list$Xocc)))
  stopifnot((ncol(timedistmat) == nrow(data.list$Xocc)) && (nrow(timedistmat) == nrow(data.list$Xocc)))
  zero.lv.v <- rep(1, nrow(data.list$Xocc))
  warning("prior means of LV set to 1")
  data.list <- c(data.list, list(spatdistmat = spatdistmat, timedistmat = timedistmat, zero.lv.v = zero.lv.v))
  return(data.list)
}

