#' @title Probability of Detection Given Occupied
#' @param Xobs A matrix of detection covariates, each row is a visit.
#' @param det.b_arrA 3-array of detection covariates. Each row is a species, each column is a covariate, and layer is a draw from the posterior.
#' @return A matrix with each row a visit, each column a species. The entry is the probability, according to the posterior, of the species being detected IF it is in occupation.
#' @export
pdetection_occupied_raw <- function(Xobs, det.b_arr){
  drawids <- 1:dim(det.b_arr)[[3]]
  pdet_l <- lapply(drawids,
         function(drawid) {
          pdet <- pdetection_occupied.ModelSite.theta(Xobs,
                                              drop_to_matrix(det.b_arr[,,drawid, drop = FALSE]))
          return(pdet)
         })
  pdet <- simplify2array(pdet_l)
  pdet_post <- apply(pdet, MARGIN = c(1, 2), mean)
  return(pdet_post)
}


#' @title Probability of Detection, Given Occupied
#' @param fixedcovar An array of detection covariate values. Each row is a visit, each column is a covariate.
#' @param loadfixed An array of loadings for the covariates in 'fixedcovar'. Each row is a species,
#'  each columns is a covariate (in same order as in fixedcovar), and each layer is a draw from the distribution of loadings.
#' @return An array of detection probability values given species occupies the visited site. Each row is a visit, each column a species, 
#' each layer a draw corresponding to the loadings.
#' @export
pdet_occ_raw.jsodm <- function(fixedcovar, loadfixed){
  stopifnot(length(dim(fixedcovar)) == 2)
  stopifnot(length(dim(loadfixed)) <= 3)
  
  det_linpred <-  tensor::tensor(fixedcovar, loadfixed, alongA = 2, alongB = 2)
  
  pdet <- exp(det_linpred) / (exp(det_linpred) + 1)   #this is the inverse logit function
  return(pdet)
}

#' Raw version for random effects in the detection takes a set of offsets
#' Random effect value is the value of a random effect for each row of the fixed covariates (in this case visits)
#' Can be constructed from a fitted object by `get_det_re(fit)[fit$data$SiteObsInd, ]`
#' @export
pdet_occ_raw.jsodm_lv_re <- function(fixedcovar, loadfixed, det.re){
  stopifnot(length(dim(fixedcovar)) == 2)
  stopifnot(length(dim(loadfixed)) <= 3)
  stopifnot(nrow(det.re) == nrow(fixedcovar))
  stopifnot(dplyr::last(dim(loadfixed)) == dplyr::last(dim(det.re)))
  # convert randomeffectval to duplicated across species (as it is per site not species)
  det.re_allspec <- rep(det.re, nrow(loadfixed))
  dim(det.re_allspec) <- c(nrow(fixedcovar), ncol(det.re), nrow(loadfixed))
  det.re_allspec <- aperm(det.re_allspec, perm = c(1, 3, 2))
  
  # add to usual linear predictor
  det_linpred <-  det.re_allspec + tensor::tensor(fixedcovar, loadfixed, alongA = 2, alongB = 2)
  
  pdet <- exp(det_linpred) / (exp(det_linpred) + 1)   #this is the inverse logit function
  return(pdet)
}

#' @export
pdet_occ <- function(fit, usethetasummary = NULL, ...){
  UseMethod("pdet_occ")
}
#' @export
pdet_occ.jsodm <- function(fit, usethetasummary = NULL){
  det.v <- fit$data$Xobs
  det.b <- get_det_b(fit, usesummary = usethetasummary)
  return(pdet_occ_raw.jsodm(det.v, det.b))
}

#' @export
pdet_occ.jsodm_lv <- pdet_occ.jsodm

pdet_occ.jsodm_lv_re <- function(fit, usethetasummary = NULL, refromposterior = TRUE){
  det.v <- fit$data$Xobs
  det.b <- get_det_b(fit, usesummary = usethetasummary)
  if (refromposterior){
    det.re <- get_det_re(fit)[fit$data$SiteObsInd, ]
  } else {
    stop("Function for simulating random effect values not implemented. Please us random effects fitted in the posterior.")
  }
  return(pdet_occ_raw.jsodm_lv_re(det.v, det.b, det.re))
}  
  