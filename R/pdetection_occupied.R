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
