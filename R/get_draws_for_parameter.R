#' Get Arrays of Draws from Posterior for Each Parameter
#' @param fit Fitted model
#' @return A named array with the final dimension corresponding to draws from the posterior.
get_occ_b <- function(fit){
  draws <- do.call(rbind, fit$mcmc)
  occ.b <- bugsvar2array(draws, "occ.b", 1:fit$data$nspecies, 1:fit$data$noccvar)
  dimnames(occ.b) <- list(Species = fit$species, Covariate = colnames(fit$data$Xocc), Draw = 1:nrow(draws))
  return(occ.b)
}

get_det_b <- function(fit){
  draws <- do.call(rbind, fit$mcmc)
  det.b <- bugsvar2array(draws, "det.b", 1:fit$data$nspecies, 1:fit$data$nobsvar)
  dimnames(det.b) <- list(Species = fit$species, Covariate = colnames(fit$data$Xobs), Draw = 1:nrow(draws))
  return(det.b)
}

get_lv_v <- function(fit){
  draws <- do.call(rbind, fit$mcmc)
  lv.v <- bugsvar2array(draws, "lv.v", 1:fit$data$nmodelsites, 1:fit$data$nlv)
  dimnames(lv.v) <- list(ModelSite = rownames(fit$data$Xocc), LV = paste0("lv", 1:fit$data$nlv, ".v"), Draw = 1:nrow(draws))
  return(lv.v)
}

get_lv_b <- function(fit){
  draws <- do.call(rbind, fit$mcmc)
  lv.b <- bugsvar2array(draws, "lv.b", 1:fit$data$nspecies, 1:fit$data$nlv)
  dimnames(lv.b) <- list(Species = fit$species, LV = paste0("lv", 1:fit$data$nlv, ".b"), Draw = 1:nrow(draws))
  return(lv.b)
}
