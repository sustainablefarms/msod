#' @title Probability of each species occupancy marginal across other species, from the full posterior distribution.
#' @param lv.v If a matrix then lv.v are considered simulated and for each possible set of parameters,
#' the probabilities are marginalised of the (prior) LV distribution. 
#' If `lv.v` is an array with 3 dimensions, then the third dimension corresponds to the posterior distribution draw (same as occ.b_arr).
#' In this situation the array must have only one row: the row corresponding the model site.
#' @param lv.b_arr A 3-array of LV loadings. If `NULL` then model is assumed to not have latent variables.
#' @param occ.b_arr A 3-array of occupancy covariates. Each row is a species, each column is a covariate, and layer is a draw from the posterior.
#' @details Compute the probability of each species individually (non-joint) using the full posterior distribution for a single ModelSite.
poccupy.ModelSite <- function(Xocc, occ.b_arr, lv.b_arr = NULL, lv.v = NULL){
  if (is.null(lv.b_arr) && is.null(lv.v)){ #model doesn't have latent variables
    pocc_l <- lapply(1:dim(occ.b_arr)[[3]],
                     function(drawid){
                       pocc <- poccupy.ModelSite.theta(Xocc, 
                                                       drop_to_matrix(occ.b_arr[,, drawid, drop = FALSE])
                       )
                       return(pocc)
                     })
  }
  
  if (length(dim(lv.v)) == 3){ #model has LV and the probabilities are conditional on the fitted LV value.
    stopifnot(dim(lv.b_arr)[[3]] == dim(lv.v)[[3]])
    stopifnot(dim(lv.v)[[1]] == 1)
    
    pocc_l <- lapply(1:dim(lv.v)[[3]],
           function(drawid){
             pocc <- poccupy.ModelSite.theta(Xocc, 
                                     drop_to_matrix(occ.b_arr[,, drawid, drop = FALSE]),
                                     drop_to_matrix(lv.b_arr[,, drawid, drop = FALSE]),
                                     drop_to_matrix(lv.v[,, drawid, drop = FALSE]))
             return(pocc)
           })
  }
  
  if (!is.null(lv.b_arr) && length(dim(lv.v)) == 2){ #this is when lv.v are meant to be marginalised for each theta
    pocc_l <- lapply(1:dim(occ.b_arr)[[3]],
                     function(drawid){
                       pocc <- poccupy.ModelSite.theta(Xocc, 
                                                       drop_to_matrix(occ.b_arr[,, drawid, drop = FALSE]),
                                                       drop_to_matrix(lv.b_arr[,, drawid, drop = FALSE]),
                                                       lv.v)
                       return(pocc)
                     })
  }
  
  pocc <- do.call(rbind, pocc_l)
  pocc_marg <- colMeans(pocc)
  return(pocc_marg)
}
