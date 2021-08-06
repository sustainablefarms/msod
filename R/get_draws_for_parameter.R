#' Get Arrays of Draws from Posterior for Each Parameter
#' @param fit Fitted model
#' @param usesummary is either media, mean, an integer, or a function to apply across draws (it is applied using the `apply` function). An integer will select that draw. Chains are concatenated, so if there are 500 samples per chain then the 5th draw of the second chain can be selected by `usesummary = 5 + 500`. `usesummary = NULL` will not summarise the parameters.
#' @return A named array with the final dimension corresponding to draws from the posterior.
#' @export
get_occ_b <- function(fit, usesummary = NULL){
  draws <- do.call(rbind, fit$mcmc)
  occ.b <- bugsvar2array(draws, "occ.b", 1:fit$data$nspecies, 1:fit$data$noccvar)
  dimnames(occ.b) <- list(Species = fit$species, Covariate = colnames(fit$data$Xocc), Draw = 1:nrow(draws))
  occ.b <- picksummary(occ.b, usesummary = usesummary)
  return(occ.b)
}

#' @describeIn get_occ_b Get detection loadings for external covariates.
#' @export
get_det_b <- function(fit, usesummary = NULL){
  draws <- do.call(rbind, fit$mcmc)
  det.b <- bugsvar2array(draws, "det.b", 1:fit$data$nspecies, 1:fit$data$nobsvar)
  dimnames(det.b) <- list(Species = fit$species, Covariate = colnames(fit$data$Xobs), Draw = 1:nrow(draws))
  det.b <- picksummary(det.b, usesummary = usesummary)
  return(det.b)
}

#' @describeIn get_occ_b Get the fitted detection random effects
get_det_re <- function(fit, usesummary = NULL){
  draws <- do.call(rbind, fit$mcmc)
  det.re <- bugsvar2array_vector(draws, "det.re", 1:fit$data$nsitegroups)
  return(det.re)
}

#' @describeIn get_occ_b Get fitted latent variable values.
#' @export
get_lv_v <- function(fit, usesummary = NULL){
  draws <- do.call(rbind, fit$mcmc)
  lv.v <- bugsvar2array(draws, "lv.v", 1:fit$data$nmodelsites, 1:fit$data$nlv)
  dimnames(lv.v) <- list(ModelSite = rownames(fit$data$Xocc), LV = paste0("lv", 1:fit$data$nlv, ".v"), Draw = 1:nrow(draws))
  lv.v <- picksummary(lv.v, usesummary = usesummary)
  return(lv.v)
}

#' @describeIn get_occ_b Get fitted latent variable loadings.
#' @export
get_lv_b <- function(fit, usesummary = NULL){
  draws <- do.call(rbind, fit$mcmc)
  lv.b <- bugsvar2array(draws, "lv.b", 1:fit$data$nspecies, 1:fit$data$nlv)
  dimnames(lv.b) <- list(Species = fit$species, LV = paste0("lv", 1:fit$data$nlv, ".b"), Draw = 1:nrow(draws))
  lv.b <- picksummary(lv.b, usesummary = usesummary)
  return(lv.b)
}

#' @describeIn get_occ_b Get fitted latent variable spatial correlation scale for jsodm_lv_sepexp model
#' @export
get_lv_v_spatscale <- function(fit, usesummary = NULL){
  draws <- do.call(rbind, fit$mcmc)
  lv.v.spatscale <- draws[ , "lv.v.spatscale[1]", drop = FALSE]
  dim(lv.v.spatscale) <- c(1, 1, nrow(draws)) #convert so that third dimension is draws
  lv.v.spatscale <- picksummary(lv.v.spatscale, usesummary = usesummary)
  return(lv.v.spatscale)
}

#' @describeIn get_occ_b Get fitted latent variable temporal correlation scale for jsodm_lv_sepexp model
#' @export
get_lv_v_timescale <- function(fit, usesummary = NULL){
  draws <- do.call(rbind, fit$mcmc)
  lv.v.timescale <- draws[ , "lv.v.timescale[1]", drop = FALSE]
  dim(lv.v.timescale) <- c(1, 1, nrow(draws)) #convert so that third dimension is draws
  lv.v.timescale <- picksummary(lv.v.timescale, usesummary = usesummary)
  return(lv.v.timescale)
}


picksummary <- function(arr, usesummary = NULL){
  if (is.null(usesummary)){return(arr)}
  if (is.numeric(usesummary) && length(usesummary) == 1){
    theta <- arr[,, usesummary, drop = FALSE]
    return(theta)
  } else {
    theta <- apply(arr, MARGIN = c(1, 2), usesummary)
    if (isTRUE(all.equal(dim(theta), dim(arr)[1:2]))){ #use summary returns a single number so apply has dropped a dimension
      dim(theta) <- c(1, dim(theta))
      dimnames(theta)[2:3] <- dimnames(arr)[1:2]
    }
    theta <- aperm(theta, perm = c(2, 3, 1)) # so that first and second dimension correspond to the same thing
    return(theta)
  }
  return(theta)
}

