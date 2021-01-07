#' @describeIn likelihoods.fit Compute the joint-species LV-marginal likelihood for all ModelSites
#' @param occ.b_arr Occupancy covariate loadings. Each row is a species, each column an occupancy covariate, each layer (dim = 3) is a draw
#' @param det.b_arr Detection covariate loadings. Each row is a species, each column an detection covariate, each layer (dim = 3) is a draw
#' @param lv.b_arr LV loadings. Each row is a species, each column a LV, each layer (dim = 3) is a draw
#' @param data_i A row of a data frame created by \code{prep_data_by_modelsite}. Each row contains data for a single ModelSite. 
#' @param lvsim A matrix of simulated LV values. Columns correspond to latent variables, each row is a simulation
#' @param Xocc A matrix of processed occupancy covariate values for the model site. Must have 1 row.
#' @param Xobs A matrix of processed detection covariate values for each visit to the model site. 
#' @param y Matrix of species detections for each visit to the model site.
#' @return `likelihoods.fit` returns a matrix. Each row corresponds to a draw of the parameters from the posterior. Each column to a ModelSite.
#' The value in each cell is the probability density, given the parameters from the draw, evaluated at the observations for the model site.
#' @export
likelihoods.jsodm_lv <- function(fit, Xocc = NULL, yXobs = NULL, ModelSite = NULL, nlvsim = 1000, cl = NULL){
  stopifnot(class(fit)[[1]] %in% c("jsodm_lv"))

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
  lv.b_arr <- get_lv_b(fit) # rows are species, columns are lv
  
  if (is.null(cl)) {
    likel.l <- lapply(1:dim(occ.b_arr)[[3]], function(drawid) {
      lkl <- likelihood_joint_LVvmarg_draw.jsodm_lv(
        Xocc, Xobs, y, ModelSite,
        occ.b_arr[,,drawid, drop = FALSE],
        det.b_arr[,,drawid, drop = FALSE],
        lv.b_arr[,,drawid, drop = FALSE],
        nlvsim
      )
      return(lkl)
    })
  }
  else {
    likel.l <- parallel::parLapply(cl = cl, 1:dim(occ.b_arr)[[3]], function(drawid) {
      lkl <- likelihood_joint_LVvmarg_draw.jsodm_lv(
        Xocc, Xobs, y, ModelSite,
        occ.b_arr[,,drawid, drop = FALSE],
        det.b_arr[,,drawid, drop = FALSE],
        lv.b_arr[,,drawid, drop = FALSE],
        nlvsim
      )
      return(lkl)
    }, FUN.VALUE = 1:nrow(Xocc) * 1.0001)
  }  
  
  likel.mat <- do.call(rbind, likel.l) # each row is a draw, each column is a modelsite (which are independent data points)
  return(likel.mat)
}


likelihood_joint_LVvmarg_draw.jsodm_lv <- function(Xocc, Xobs, y, ModelSite, occ.b, det.b, lv.b, nlvsim){
  stopifnot( dim(occ.b)[[3]] == 1 )
  stopifnot( dim(det.b)[[3]] == 1 )
  stopifnot( dim(lv.b)[[3]] == 1 )
  
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
  Occ.Pred.CondLV <- vapply(1:nlvsim,
         function(i){
           lv.v <- array(rnorm(nrow(Xocc) * dim(lv.b)[[2]]), dim = c(nrow(Xocc), dim(lv.b)[[2]], 1))
           pocc <- poccupy_raw.jsodm_lv(Xocc, occ.b, lv.v, lv.b)
           return(pocc)
         },
         FUN.VALUE = array(0, dim = c(nrow(Xocc), dim(occ.b)[[1]], 1))
  )
  
  # combine with likelihoods of detections, use Rfast::eachrow after reshaping the arrays
  dim(Occ.Pred.CondLV) <- c(nrow(Xocc) * nrow(occ.b), nlvsim)
  dim(Likl_condoccupied.JointVisit) <- c(nrow(Xocc) * nrow(occ.b), 1)
  Likl.JointVisit.condLV <- Rfast::eachrow(Occ.Pred.CondLV, Likl_condoccupied.JointVisit, oper = "*") #per species likelihood, occupied component. Works because species conditionally independent given LV
  
  dim(NoneDetected) <- c(nrow(Xocc) * nrow(occ.b), 1)
  Likl.JointVisit.condLV <- Likl.JointVisit.condLV + 
    Rfast::eachrow((1 - Occ.Pred.CondLV), NoneDetected, oper = "*") #add probability of unoccupied for zero detections
  
  
  # combine likelihoods of detections between species
  dim(Likl.JointVisit.condLV) <- c(nrow(Xocc), dim(occ.b)[[1]], nlvsim)
  Likl.JointVisitSp.condLV <- apply(Likl.JointVisit.condLV, MARGIN = c(1, 3), prod) # multiply probabilities of each species together because species are conditionally independent
  
  # take mean of all LV sims to get likelihood marginalised across LV values, given the single draw
  Likl_margLV <- rowMeans(Likl.JointVisitSp.condLV)
  
  return(Likl_margLV)
}
