#' @title Expected Biodiversity Across Multiple Sites
#' @description Assuming independent patches, estimates the posterior expected number of species in across all patches.

expectedspeciesnum.theta <- function(Xocc, occ.b, lv.b = NULL, LVvals = NULL){
  ## Probability of Site Occupancy
  stopifnot(nrow(Xocc) == nrow(lv.b))
  if (is.null(Xobs)) {stopifnot(is.null(det.b))}
  Xocc <- as.matrix(Xocc)
  if (is.null(lv.b)){ # create dummy LV  
    stopifnot(is.null(LVvals))
    lv.b <- matrix(0, nrow = nrow(occ.b), ncol = 2)
    LVvals <- matrix(0, nrow = 2, ncol = 2)
  } # dummy LV
  lapply()
  ModelSite.Occ.Pred.CondLV <- poccupy.ModelSite.theta(Xocc, occ.b, lv.b, LVvals)
  
  ## Expected number of species occupying modelsite, given model and theta, and marginal across LV
  EVnumspec_occ <- Erowsum_margrow(ModelSite.Occ.Pred.CondLV)
  names(EVnumspec_occ) <- paste0(names(EVnumspec_occ), "_occ")
  return(EVnumspec_occ)
}