#' @title Expected number of detections per species per site
#' @param fit A fitted jsodm
#' @param desiredspecies List of species names to include in the calculations
#' @return An array of expected number of detections per modelsite x species x draw.
#' @description For each species individually calculates the expected number of detections at each ModelSite.
#' As this is the expected value, and each species treated individually, the same method can be applied to jsodm_lv models.
#' @examples 
#' fit <- readRDS("../sflddata/private/data/testdata/cutfit_7_4_11_2LV.rds")
#' fit <- translatefit(fit)
#' Endet <- Edetperspecies.jsodm(fit)
#' @export
Edetperspecies.jsodm <- function(fit, 
                                  desiredspecies = fit$species){
  stopifnot(all(desiredspecies %in% fit$species))
  occ.v <- fit$data$Xocc
  occ.b <- get_occ_b(fit)[desiredspecies, , , drop = FALSE]
  pocc <- poccupy_raw.jsodm(occ.v, occ.b)
  det.b <- get_det_b(fit)[desiredspecies, , , drop = FALSE]
  det.v <- fit$data$Xobs
    
  pdet_occ <- pdet_occ_raw.jsodm(det.v, det.b)
  Endet <- detperspecies_probarr(pocc, pdet_occ, fit$data$ModelSite)
  return(Endet)
}
 
# pocc = 3 array of probability of occupancy. Within each draw occupancy of species considered independent.
# pdet_occ = 3 array of probability of detection, condition on occupancy. Within each draw species detection considered independent given occupied.
# returns an array of expected detection number. Dimensions of modelsite, species, draw
detperspecies_probarr <- function(pocc, pdet_occ, ModelSite){
  ModSiteVisits <- plyr::split_indices(ModelSite)
  names(ModSiteVisits) <- lapply(ModSiteVisits, function(x) ModelSite[x[[1]]]) %>% unlist()
  
  # per species expected number of detections in a site when the site is occupied. 
  # Depends on the number of visits
  Endet_occ <- vapply(ModSiteVisits,
                             function(visits){
                               Endet <- apply(pdet_occ[visits, , , drop = FALSE], MARGIN = c(2, 3), sum)
                               return(Endet)
                             },
                             FUN.VALUE = matrix(1.1, nrow = ncol(pdet_occ), ncol = dim(pdet_occ)[[3]]))
  if (is.null(dim(Endet_occ))) {
    dim(Endet_occ) <- dim(pocc)[c(2, 3, 1)]
    dimnames(Endet_occ) <- dimnames(pocc)[c(2, 3, 1)]
  }
  #Endet_occ array of species x draw x ModelSite
  Endet_occ <- aperm(Endet_occ, perm = c(3, 1, 2)) #now Model Site x Species x Draw
  
  ## probability of no detections for all visits to a site, marginal on occupancy
  Endet <- Endet_occ * pocc + 0
  return(Endet)
}

#' @describeIn Edetperspecies.jsodm Get observed detections per model site x species from the 'y' data in the fitted object.
#' @export
Odetperspecies <- function(fit, y = NULL){
  ModSiteVisits <- plyr::split_indices(fit$data$ModelSite)
  names(ModSiteVisits) <- lapply(ModSiteVisits, function(x) fit$data$ModelSite[x[[1]]]) %>% unlist()
  
  if (is.null(y)){y <- fit$data$y}
  
  # per species number of detections at each site 
  # Depends on the number of visits
  Ondet <- vapply(ModSiteVisits,
                      function(visits){
                        ndet <- colSums(y[visits, , drop = FALSE])
                        ndet <- matrix(ndet, nrow = 1, ncol = ncol(y)) 
                        colnames(ndet) <- colnames(y)
                        return(ndet)
                      },
                      FUN.VALUE = matrix(1.1, nrow = 1, ncol = ncol(y)))
  #Ondet array of 1 x species x modelsite
  Ondet <- aperm(Ondet, perm = c(3, 2, 1)) #now Model Site x Species x 1
  return(Ondet)
}
