#' @title A helper function to get a vector of parameters from a fitted object
#' @param fit fitted runjags object with summary included
#' @param type An integer will select the draw from the posterior
#' (e.g. type = 5 will return the vector of parameters in the 5th sample from the posterior, mcmc chains are concatenated)
#' A charactor vector can be used to the parameters based on summarary statistics like quantiles and moments.
#' Supported values so far "median" and "mean".
#' If [type] is "marginal" or "all" then all draws are returned.
#' @export
get_theta <- function(fit, type){
  if (is.numeric(type) && length(type) > 1 && !is.null(names(type))){ #assumed passed 'type' is actually the desired theta
    return(type)
  }
  if (is.numeric(type) && length(type) == 1){
    chainidx <- floor(type / (fit$sample + 1 )) + 1
    sampleinchain <- type - fit$sample * (chainidx - 1)
    theta <- fit$mcmc[[chainidx]][sampleinchain, ]
    return(theta)
  }
  if (type == "median"){theta <- apply(do.call(rbind, fit$mcmc), MARGIN = 2, median)}
  if (type == "mean"){theta <- apply(do.call(rbind, fit$mcmc), MARGIN = 2, mean)}
  if (type == "marginal" | type == "all"){theta <- do.call(rbind, fit$mcmc)}
  return(theta)
}