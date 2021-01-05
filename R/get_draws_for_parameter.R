#' Get Arrays of Draws from Posterior for Each Parameter
#' @param fit Fitted model
#' @return A named array with the final dimension corresponding to draws from the posterior.
get_occ_b <- function(fit){
  draws <- do.call(rbind, fit$mcmc)
  occ.b <- bugsvar2array(draws, "occ.b", 1:fit$data$nspecies, 1:fit$data$noccvar)
  dimnames(occ.b) <- list(Species = fit$species, Covariate = colnames(fit$data$Xocc), Draw = 1:nrow(draws))
  return(occ.b)
}
