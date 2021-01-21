#' @title Simulate many independent Bernoulli random variables
#' @param probarr An array of probability of success.
#'  Each element of the array represents the probability of success of a Bernoulli random variable that is independent of all other random variables.
#' @return 
#' An array the same dimensions of probarr with simulated values of TRUE (success) or FALSE (failure).
#' @export
rmanybern <- function(probarr){
  ind <- runif(length(probarr), min = 0, max = 1)
  dim(ind) <- dim(probarr)
  success <- (ind <= probarr) * 1 #the * 1 converts it to numeric without altering dimensions, as.integer converts to vector as well as integer values
  dimnames(success) <- dimnames(probarr)
  return(success)
}
