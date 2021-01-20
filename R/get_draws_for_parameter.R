#' Get Arrays of Draws from Posterior for Each Parameter
#' @param fit Fitted model
#' @param usesummary is either media, mean, or a function to apply across draws (it is applied using the `apply` function)
#' @return A named array with the final dimension corresponding to draws from the posterior.
#' @export
get_occ_b <- function(fit, usesummary = NULL){
  draws <- do.call(rbind, fit$mcmc)
  occ.b <- bugsvar2array(draws, "occ.b", 1:fit$data$nspecies, 1:fit$data$noccvar)
  occ.b <- picksummary(occ.b, usesummary = usesummary)
  dimnames(occ.b) <- list(Species = fit$species, Covariate = colnames(fit$data$Xocc), Draw = 1:nrow(draws))
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

#' @describeIn get_occ_b Get fitted latent variable values.
#' @export
get_lv_v <- function(fit, usesummary = NULL){
  draws <- do.call(rbind, fit$mcmc)
  lv.v <- bugsvar2array(draws, "lv.v", 1:fit$data$nmodelsites, 1:fit$data$nlv)
  lv.v <- picksummary(lv.v, usesummary = usesummary)
  dimnames(lv.v) <- list(ModelSite = rownames(fit$data$Xocc), LV = paste0("lv", 1:fit$data$nlv, ".v"), Draw = 1:nrow(draws))
  return(lv.v)
}

#' @describeIn get_occ_b Get fitted latent variable loadings.
#' @export
get_lv_b <- function(fit, usesummary = NULL){
  draws <- do.call(rbind, fit$mcmc)
  lv.b <- bugsvar2array(draws, "lv.b", 1:fit$data$nspecies, 1:fit$data$nlv)
  lv.b <- picksummary(lv.b, usesummary = usesummary)
  dimnames(lv.b) <- list(Species = fit$species, LV = paste0("lv", 1:fit$data$nlv, ".b"), Draw = 1:nrow(draws))
  return(lv.b)
}

picksummary <- function(arr, usesummary = NULL){
  if (is.null(usesummary)){return(arr)}
  if (is.numeric(type) && length(type) == 1){
    chainidx <- floor(type / (fit$sample + 1 )) + 1
    sampleinchain <- type - fit$sample * (chainidx - 1)
    theta <- fit$mcmc[[chainidx]][sampleinchain, ]
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
