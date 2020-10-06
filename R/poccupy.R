
#' @param LVvals If a matrix then LVvals are considered simulated and for each possible set of parameters,
#' the probabilities are marginalised of the (prior) LV distribution. 
#' If `LVvals` is an array with 3 dimensions, then the third dimension corresponds to the posterior distribution draw (same as u.b_arr).
#' In this situation the array must have only one row: the row corresponding the model site.
#' @param lv.coef_arr A 3-array of LV loadings. If `NULL` then model is assumed to not have latent variables.
#' @param u.b_arr A 3-array of occupancy covariates. Each row is a species, each column is a covariate, and layer is a draw from the posterior.
#' @details Compute the probability of each species individually (non-join) using the full posterior distribution.
poccupy.ModelSite <- function(Xocc, u.b_arr, lv.coef_arr = NULL, LVvals = NULL){
  if (is.null(lv.coef_arr) && is.null(LVvals)){
    lv.coef <- matrix(0, nrow = nrow(u.b_arr), ncol = 2)
    LVvals <- matrix(0, nrow = 1, ncol = 2)
    pocc_l <- lapply(1:dim(u.b_arr)[[3]],
                     function(drawid){
                       pocc <- poccupy.ModelSite.theta(Xocc, 
                                                       u.b_arr[,, drawid],
                                                       lv.coef,
                                                       LVvals)
                       return(pocc)
                     })
  }
  
  if (length(dim(LVvals)) == 3){
    stopifnot(dim(lv.coef_arr)[[3]] == dim(LVvals)[[3]])
    stopifnot(dim(LVvals)[[1]] == 1)
    
    pocc_l <- lapply(1:dim(LVvals)[[3]],
           function(drawid){
             pocc <- poccupy.ModelSite.theta(Xocc, 
                                     u.b_arr[,, drawid],
                                     lv.coef_arr[,, drawid],
                                     LVvals[,, drawid])
             return(pocc)
           })
  }
  
  if (!is.null(lv.coef_arr) && length(dim(LVvals)) == 2){ #this is when LVvals are meant to be marginalised for each theta
    pocc_l <- lapply(1:dim(u.b_arr)[[3]],
                     function(drawid){
                       pocc <- poccupy.ModelSite.theta(Xocc, 
                                                       u.b_arr[,, drawid],
                                                       lv.coef_arr[,, drawid],
                                                       LVvals)
                       return(pocc)
                     })
  }
  
  pocc <- do.call(rbind, pocc_l)
  pocc_marg <- colMeans(pocc)
  return(pocc_marg)
}
