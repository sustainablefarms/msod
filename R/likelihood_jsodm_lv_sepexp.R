#' @title Likelihood Computation for JSODM_LV Model

#' @description  Computes the likelihood of observations at each ModelSite. At data in the fitted model, or on new data supplied.
#' It uses [likelihood_joint_LVvmarg_draw.jsodm_lv].
#' @param chains is a vector indicator which mcmc chains to extract draws from. If NULL then all chains used.
#' @param numlvsims the number of simulated latent variable values to use for computing likelihoods
#' @param cl a cluster created by parallel::makeCluster()
#' @return Returns a matrix. Each row corresponds to a draw of the parameters from the posterior. Each column to a ModelSite.
#' The value in each cell is the probability density, given the parameters from the draw, evaluated at the observations for the model site.
#' @export
likelihood.jsodm_lv_sepexp <- function(fit,
                                numlvsims = 1000, cl = NULL, simseed = NULL){
  stopifnot(class(fit)[[1]] %in% c("jsodm_lv"))

  Xocc <- fit$data$Xocc
  Xobs <- fit$data$Xobs
  y <- fit$data$y
  ModelSite <- fit$data$ModelSite
  
  occ.b_arr <- get_occ_b(fit) # rows are species, columns are occupancy covariates
  det.b_arr <- get_det_b(fit)  # rows are species, columns are observation covariates
  lv.b_arr <- get_lv_b(fit) # rows are species, columns are lv
  lv.v.spatscale <- get_lv_v_spatscale(fit)
  lv.v.timescale <- get_lv_v_timescale(fit)
  
  if (is.null(cl)) {
    likel.l <- lapply(1:dim(occ.b_arr)[[3]], function(drawid) {
      lkl <- likelihood_LVvmarg_draw.jsodm_lv_sepexp(
        Xocc, Xobs, y, ModelSite,
        occ.b_arr[,,drawid, drop = FALSE],
        det.b_arr[,,drawid, drop = FALSE],
        lv.b_arr[,,drawid, drop = FALSE],
        lv.v.spatscale[,,drawid, drop = FALSE],
        lv.v.timescale[,,drawid, drop = FALSE],
        fit$data$spatdistmat,
        fit$data$timedistmat,
        fit$data$zero.lv.v,
        numlvsims,
        simseed
      )
      return(lkl)
    })
  }
  else {
    likel.l <- parallel::parLapply(cl = cl, 1:dim(occ.b_arr)[[3]], function(drawid) {
      lkl <- likelihood_LVvmarg_draw.jsodm_lv_sepexp(
        Xocc, Xobs, y, ModelSite,
        occ.b_arr[,,drawid, drop = FALSE],
        det.b_arr[,,drawid, drop = FALSE],
        lv.b_arr[,,drawid, drop = FALSE],
        lv.v.spatscale[,,drawid, drop = FALSE],
        lv.v.timescale[,,drawid, drop = FALSE],
        fit$data$spatdistmat,
        fit$data$timedistmat,
        fit$data$zero.lv.v,
        numlvsims,
        simseed
      )
      return(lkl)
    })
  }  
  
  likel.mat <- do.call(rbind, likel.l) # each row is a draw, each column is a modelsite (which are independent data points)
  return(likel.mat)
}

#' @title Likelihood Computation for JSODM_LV Model for a Given Draw
#' @description Computes the joint-species LV-marginal likelihood for all ModelSites, for a given draw.
#' @param occ.b_arr Occupancy covariate loadings. Each row is a species, each column an occupancy covariate, each layer (dim = 3) is a draw
#' @param det.b_arr Detection covariate loadings. Each row is a species, each column an detection covariate, each layer (dim = 3) is a draw
#' @param lv.b_arr LV loadings. Each row is a species, each column a LV, each layer (dim = 3) is a draw
#' @param data_i A row of a data frame created by \code{prep_data_by_modelsite}. Each row contains data for a single ModelSite. 
#' @param lvsim A matrix of simulated LV values. Columns correspond to latent variables, each row is a simulation
#' @param Xocc A matrix of processed occupancy covariate values for the model site. Must have 1 row.
#' @param Xobs A matrix of processed detection covariate values for each visit to the model site. 
#' @param y Matrix of species detections for each visit to the model site.
likelihood_LVvmarg_draw.jsodm_lv_sepexp <- function(Xocc, Xobs, y, ModelSite, occ.b, det.b, lv.b, lv.v.spatscale, lv.v.timescale, spatdistmat, timedistmat, zero.lv.v, numlvsims, simseed = NULL){
  stopifnot( dim(occ.b)[[3]] == 1 )
  stopifnot( dim(det.b)[[3]] == 1 )
  stopifnot( dim(lv.b)[[3]] == 1 )
  
  ## Probability of Detection, CONDITIONAL on occupied
  Detection.Pred.Cond <- pdet_occ_raw.jsodm(Xobs, det.b)[,,1]
  
  ## Likelihood (probability) of detection for each visit, given occupied (detections are independent given occupation status)
  Likl_condoccupied <- Detection.Pred.Cond * y + (1 - Detection.Pred.Cond) * (1 - y) # non-detection is complement of detection probability
  
  ModSiteVisits <- plyr::split_indices(ModelSite) #10 times faster than split(1:nrow(Likl_condoccupied), ModelSite). Gives a list of visit id (row) for each model site
  
  ## Joint likelihood (probability) of detections of each species for all visits to each model site, CONDITIONAL on occupied
  # many of the following arrays are transposed to work with Rfast::eachrow, with each row being an LVv sim
  Likl_condoccupied.JointVisit <- vapply(ModSiteVisits,
         function(visits){
           Likl_condoccupied.JointVisit.ModelSite <- apply(Likl_condoccupied[visits, , drop = FALSE], MARGIN = 2, prod)
           return(Likl_condoccupied.JointVisit.ModelSite)
         },
         FUN.VALUE = (1:nrow(occ.b)) * 1.001)
  Likl_condoccupied.JointVisit <- t(Likl_condoccupied.JointVisit)
  # Likl_condoccupied.JointVisit is now a matrix of site x species
  
  ## Likelihood (probability) of y given unoccupied is either 1 or 0 for detections. Won't include that here yet.
  NoneDetected <- vapply(ModSiteVisits,
              function(visits){
                NoneDetected_modelsite <- as.numeric(colSums(y[visits, , drop = FALSE]) == 0)
                return(NoneDetected_modelsite)
              },
              FUN.VALUE = (1:nrow(occ.b)) * 1.001)
  NoneDetected <- t(NoneDetected)
  
  ## Probability of Site Occupancy, for each simulated LV, for the given draw. #this seems to be the SLOWEST part
  Occ.Pred.CondLV <- vapply(1:numlvsims,
         function(i){
           if (!is.null(simseed)){set.seed(simseed)}
           lv.v_site <- Rfast::rmvnorm(1, mu = zero.lv.v,
                                       sigma = exp(-spatdistmat/lv.v.spatscale[1]) * exp(-timedistmat/lv.v.timescale[1]) )
           lv.v <- matrix(lv.v_site, nrow = nrow(Xocc), ncol = dim(lv.b)[[2]], byrow = TRUE)
           dim(lv.v) <- c(dim(lv.v), 1)
           pocc <- poccupy_raw.jsodm_lv(Xocc, occ.b, lv.v, lv.b)
           return(pocc)
         },
         FUN.VALUE = array(0, dim = c(nrow(Xocc), dim(occ.b)[[1]], 1))
  )
  # Occ.Pred.CondLV is array sites x species x draw x LV sim
  
  # combine with likelihoods of detections, use Rfast::eachrow after reshaping the arrays
  dim(Occ.Pred.CondLV) <- c(nrow(Xocc) * nrow(occ.b), numlvsims) #Occ.Pred.CondLV is now rows that are (site, species; site moving first), with columns of lv sim
  dim(Likl_condoccupied.JointVisit) <- c(nrow(Xocc) * nrow(occ.b), 1) # Likl_condoccupied.JointVisit is now rows that are (site, species; site moving first), with column of 1
  Likl.JointVisit.condLV <- Rfast::eachrow(t(Occ.Pred.CondLV), t(Likl_condoccupied.JointVisit), oper = "*") #per species likelihood, occupied component. Works because species conditionally independent given LV
  
  dim(NoneDetected) <- c(nrow(Xocc) * nrow(occ.b), 1) #vector, first through sites, then through species
  Likl.JointVisit.condLV <- Likl.JointVisit.condLV + 
    Rfast::eachrow((1 - t(Occ.Pred.CondLV)), t(NoneDetected), oper = "*") #add probability of unoccupied for zero detections
  
  
  # combine likelihoods of detections between species
  Likl.JointVisit.condLV <- t(Likl.JointVisit.condLV)
  dim(Likl.JointVisit.condLV) <- c(nrow(Xocc), dim(occ.b)[[1]], numlvsims)
  Likl.JointVisitSp.condLV <- apply(Likl.JointVisit.condLV, MARGIN = c(1, 3), prod) # multiply probabilities of each species together because species are conditionally independent
  
  # take mean of all LV sims to get likelihood marginalised across LV values, given the single draw
  Likl_margLV <- Rfast::rowmeans(Likl.JointVisitSp.condLV)
  
  return(Likl_margLV)
}
