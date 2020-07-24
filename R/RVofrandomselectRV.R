#' @title For a list of random variables, return random variable that is a random selection one of the RVs
#' @param rvs A list of Random Variable (distributions) of class 'RV'
#' @examples 
#' library(discreteRV)
#' X1 <- RV(c(1,0), c(.5,.5))
#' X2 <- RV(c(1,0), c(.5,.5))
#' X3 <- RV(c(1,-1), c(.5,.5))
#' randselRV(list(X1, X2, X3))
sumRV_margrow <- function(pmat){
  sum_RVs <- lapply(1:nrow(pmat),
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
