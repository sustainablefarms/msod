#' @title For a list of random variables, return random variable that is a random selection one of the RVs
#' @param rvs A list of Random Variable (distributions) of class 'RV'
#' @examples 
#' library(discreteRV)
#' X1 <- RV(c(1,0), c(.5,.5))
#' X2 <- RV(c(1,0), c(.5,.5))
#' X3 <- RV(c(1,-1), c(.5,.5))
randselRV(list(X1, X2, X3))
fit <- readRDS("./tmpdata/7_2_9_addyear_msnm_year_time_2lv.rds")
theta <- get_theta(fit, type = "median")
u.b <- bugsvar2matrix(theta, "u.b", 1:fit$data$n, 1:ncol(fit$data$Xocc))
v.b <- bugsvar2matrix(theta, "v.b", 1:fit$data$n, 1:ncol(fit$data$Xobs))
lv.coef <- bugsvar2matrix(theta, "lv.coef", 1:fit$data$n, 1:fit$data$nlv)
LVvals <- bugsvar2matrix(theta, "LV", 1:fit$data$J, 1:fit$data$nlv)
poccupy <- poccupy.ModelSite.theta(fit$data$Xocc[1, , drop = FALSE], u.b, lv.coef, LVvals[1, , drop = FALSE])
oRVs <- lapply(poccupy, function(p) RV(c(0, 1), probs = c(1 - p, p)))
sumRV <- do.call(SofI, oRVs)

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
