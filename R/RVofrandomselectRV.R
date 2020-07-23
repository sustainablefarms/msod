#' @title For a list of random variables, return random variable that is a random selection one of the RVs
#' @param rvs A list of Random Variable (distributions) of class 'RV'
#' @examples 
#' library(discreteRV)
#' X1 <- RV(c(1,0), c(.5,.5))
#' X2 <- RV(c(1,0), c(.5,.5))
#' X3 <- RV(c(1,-1), c(.5,.5))
randselRV(list(X1, X2, X3))

#' full distribution of species number given fixed theta and unknown LV value (i.e. many simulated LV values)
fit <- readRDS("./tmpdata/7_2_9_addyear_msnm_year_time_2lv.rds")
theta <- get_theta(fit, type = "median")
u.b <- bugsvar2matrix(theta, "u.b", 1:fit$data$n, 1:ncol(fit$data$Xocc))
v.b <- bugsvar2matrix(theta, "v.b", 1:fit$data$n, 1:ncol(fit$data$Xobs))
lv.coef <- bugsvar2matrix(theta, "lv.coef", 1:fit$data$n, 1:fit$data$nlv)
nLVsim <- 1000
lvsim <- matrix(rnorm(fit$data$nlv * nLVsim), ncol = fit$data$nlv, nrow = nLVsim)
sumRV_margLV <- speciesnum.ModelSite.theta(Xocc = fit$data$Xocc[1, , drop = FALSE], u.b = u.b, lv.coef = lv.coef, LVvals = lvsim)
quantile(sumRV_margLV, probs = c(0.025, 0.5, 0.975))

speciesnum.ModelSite.theta <- function(Xocc, Xobs = NULL, u.b, v.b = NULL, lv.coef = NULL, LVvals = NULL){
  ## Probability of Site Occupancy
  stopifnot(nrow(Xocc) == 1)
  if (is.null(Xobs)) {stopifnot(is.null(v.b))}
  Xocc <- as.matrix(Xocc)
  if (is.null(lv.coef)){ # create dummy LV  
    stopifnot(is.null(LVvals))
    lv.coef <- matrix(0, nrow = nrow(u.b), ncol = 2)
    LVvals <- matrix(0, nrow = 2, ncol = 2)
  } # dummy LV
  ModelSite.Occ.Pred.CondLV <- poccupy.ModelSite.theta(Xocc, u.b, lv.coef, LVvals)
  
  if (is.null(Xobs)){
    sum_occ_margLV <- sumRV_margrow(ModelSite.Occ.Pred.CondLV)
    return(sum_occ_margLV)
  }
  
  ## Probability of Detection
  if (!is.null(Xobs)){
    Detection.Pred.Cond <- pdetection_occupied.ModelSite.theta(Xobs, v.b) # probability conditional on occupied
    # probability of no detections, given occupied
    NoDetections.Pred.Cond <- Rfast::colprods(1 - Detection.Pred.Cond)
    
    
    ## probability of no detections, marginal on occupancy
    NoDetections.Pred.marg_1 <- Rfast::eachrow(ModelSite.Occ.Pred.CondLV, NoDetections.Pred.Cond, oper = "*")  #occupied component
    NoDetections.Pred.marg_2 <- 1 - ModelSite.Occ.Pred.CondLV  #plus unoccupied component
    AnyDetections.Pred.marg <- 1 - (NoDetections.Pred.marg_1 + NoDetections.Pred.marg_2)
    sum_det_margLV <- sumRV_margrow(AnyDetections.Pred.marg)
    return(sum_det_margLV)
  }
}

sumRV_margrow <- function(pmat){
  sum_RVs <- pbapply::pblapply(1:nrow(pmat),
                                   function(row) {
                                     LVs <- lapply(pmat[row, ], function(p) RV(c(0, 1), probs = c(1 - p, p)))
                                     sumRV <- do.call(SofI, LVs)
                                     return(sumRV)
                                   })
  sumRV_marg <- randselRV(sum_RVs)
  return(sumRV_marg)
}

randselRV <- function(rvs, weights = rep(1, length(rvs))){
  outcomes <- unique(unlist(lapply(rvs, outcomes)))
  pmf <- vapply(outcomes, function(value) Poutcome(value, rvs, weights), FUN.VALUE = 0.2)
  outRV <- RV(outcomes,
              probs = pmf,
              fractions = FALSE)
  return(outRV)
}

# For a given value, compute the probability that a randomly selected RV will have that value
# weights are the selection weights of the RVs
Poutcome <- function(value, rvs, weights = rep(1, length(rvs))){
  Ps <- vapply(rvs, function(X) P(X == value), FUN.VALUE = 0.2)
  prob <- weighted.mean(Ps, w = weights)
  return(prob)
}
