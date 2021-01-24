#' @title Likelihood for Plain JSODM

#' @export
likelihood.jsodm <- function(fit, Xocc = NULL, yXobs = NULL, ModelSite = NULL, cl = NULL, numlvsims = NULL){
  stopifnot(class(fit)[[1]] %in% c("jsodm"))
  
  if (is.null(Xocc)){ #Extract the Xocc, yXobs etc from the fitted object, no preprocessing required
    sitedata <- fit$data
  } else {
    sitedata <- prep_new_data(fit, Xocc, yXobs, ModelSite)
  }
  
  Xocc <- sitedata$Xocc
  Xobs <- sitedata$Xobs
  y <- sitedata$y
  ModelSite <- sitedata$ModelSite
  
  occ.b_arr <- get_occ_b(fit) # rows are species, columns are occupancy covariates
  det.b_arr <- get_det_b(fit)  # rows are species, columns are observation covariates
  
  if (is.null(cl)) {
    likel.l <- lapply(1:dim(occ.b_arr)[[3]], function(drawid) {
      lkl <- likelihood_draw.jsodm(
        Xocc, Xobs, y, ModelSite,
        occ.b_arr[,,drawid, drop = FALSE],
        det.b_arr[,,drawid, drop = FALSE]
      )
      return(lkl)
    })
  }
  else {
    likel.l <- parallel::parLapply(cl = cl, 1:dim(occ.b_arr)[[3]], function(drawid) {
      lkl <- likelihood_draw.jsodm(
        Xocc, Xobs, y, ModelSite,
        occ.b_arr[,,drawid, drop = FALSE],
        det.b_arr[,,drawid, drop = FALSE]
      )
      return(lkl)
    }, FUN.VALUE = 1:nrow(Xocc) * 1.0001)
  }  
  
  likel.mat <- do.call(rbind, likel.l) # each row is a draw, each column is a modelsite (which are independent data points)
  return(likel.mat)
}


likelihood_draw.jsodm <- function(Xocc, Xobs, y, ModelSite, occ.b, det.b){
  stopifnot( dim(occ.b)[[3]] == 1 )
  stopifnot( dim(det.b)[[3]] == 1 )
  
  ## Probability of Detection, CONDITIONAL on occupied
  Detection.Pred.Cond <- pdet_occ_raw.jsodm(Xobs, det.b)[,,1]
  
  ## Likelihood (probability) of detection for each visit, given occupied (detections are independent given occupation status)
  Likl_condoccupied <- Detection.Pred.Cond * y + (1 - Detection.Pred.Cond) * (1 - y) # non-detection is complement of detection probability
  
  ModSiteVisits <- plyr::split_indices(ModelSite) #10 times faster than split(1:nrow(Likl_condoccupied), ModelSite). Gives a list of visit id (row) for each model site
  
  ## Joint likelihood (probability) of detections of each species for all visits to each model site, CONDITIONAL on occupied
  a <- vapply(ModSiteVisits,
              function(visits){
                Likl_condoccupied.JointVisit.ModelSite <- apply(Likl_condoccupied[visits, , drop = FALSE], MARGIN = 2, prod)
                return(Likl_condoccupied.JointVisit.ModelSite)
              },
              FUN.VALUE = (1:nrow(occ.b)) * 1.001)
  Likl_condoccupied.JointVisit <- t(a)
  
  ## Likelihood (probability) of y given unoccupied is either 1 or 0 for detections. Won't include that here yet.
  a <- vapply(ModSiteVisits,
              function(visits){
                NoneDetected_modelsite <- as.numeric(colSums(y[visits, , drop = FALSE]) == 0)
                return(NoneDetected_modelsite)
              },
              FUN.VALUE = (1:nrow(occ.b)) * 1.001)
  NoneDetected <- t(a)
  
  ## Probability of Site Occupancy, for each simulated LV, for the given draw. #this seems to be the SLOWEST part
  Occ.Pred  <- poccupy_raw.jsodm(Xocc, occ.b)[,,1]

  # combine with likelihoods of detections
  Likl.JointVisit <- Occ.Pred * Likl_condoccupied.JointVisit
  Likl.JointVisit <- Likl.JointVisit + 
    (1 - Occ.Pred) * NoneDetected #add probability of unoccupied for zero detections
  
  # combine likelihoods of detections between species
  Likl.JointVisitSp <- apply(Likl.JointVisit, MARGIN = 1, prod) # multiply probabilities of each species together because species are conditionally independent
  
  return(Likl.JointVisitSp)
}
