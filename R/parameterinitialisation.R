#' @title Initial condition functions
#' @describe Specifies the initial conditions for the MCMC chains.
#' @param chain Integer. Index of the chain.
#' @param indata A list of data that is typically passed to [runjags::run.jags()]
#' @export
paraminits <- function(modeltype, chain, indata, ...){
  out <- do.call(paste0("paraminits.", modeltype),
    c( list(chain, indata), ...)
  )
  return(out)
}

#' @describeIn paraminits Initialise paramers for the jsodm_lv model
#' @export
paraminits.jsodm_lv <- function(chain, indata, ...){
  out <- paraminits.jsodm(chain, indata, ...)
  # lv.coef<-matrix(1, indata$n, indata$nlv)
  # lv.coef[1:indata$nlv,1:indata$nlv]<-0
  # for(l in 1:indata$nlv-1){
  #   lv.coef[l,(l+1):indata$nlv]<-NA
  # }
  # LV <- matrix(rnorm(indata$nlv * indata$J), indata$J, indata$nlv)
  # out <- c(out,
  #          lv.coef = lv.coef,
  #          LV = LV)
  return(out)
}

paraminits.jsodm_lv_sepexp <- function(chain, indata, ...){
  out <- paraminits.jsodm(chain, indata, ...)
  return(out)
}

#' @describeIn paraminits Initialise parameters for the plain jsodm model
#' @export
paraminits.jsodm <- function(chain, indata, ...) {
  v.b.proto <- lapply(colnames(indata$y),
                      function(x) {unname(coef(glm(((indata$y>0)*1)[, x] ~ . - 1, #intercept is built in
                                                   data = data.frame(indata$Xobs),
                                                   family=binomial(link=logit))))})
  v.b <- t(do.call(cbind, v.b.proto))
  v.b[v.b > 5] <- 5  #remove extremes as coefficients are assumed to be from a standard gaussian
  v.b[v.b < -5] <- -5  #remove extremes as coefficients are assumed to be from a standard gaussian

  ## this is calculated just to get initial values for occupancy covariates and occupancy estimates
  y.occ.mock <- cbind(ModelSiteID = indata$ModelSite, indata$y) %>%
    tibble::as_tibble() %>%
    dplyr::group_by(ModelSiteID) %>%
    dplyr::summarise_all(max) %>%
    dplyr::select(-ModelSiteID)
  u.b.proto <- lapply(colnames(y.occ.mock),
                      function(x) {unname(coef(glm( ((y.occ.mock>0)*1)[, x] ~ . - 1, #intercept is built in
                                                    data = data.frame(indata$Xocc),
                                                    family=binomial(link=probit))))})
  
  u.b <- t(do.call(cbind, u.b.proto))
  u.b[u.b > 5] <- 5  #remove extremes as coefficients are assumed to be from a standard gaussian
  u.b[u.b < -5] <- -5  #remove extremes as coefficients are assumed to be from a standard gaussian

  .RNG.seed <- c(1, 2, 3, 4, 5)[chain] # run jags likes to have this and the following argument defined in the initlalization functions.
  .RNG.name <- c(rep(c("base::Super-Duper", "base::Wichmann-Hill"),3))[chain]
  
  
  out <- list(
    u.b= u.b,  #initial values guestimated from u.b.proto are erroring! "u[14,1]: Node inconsistent with parents"
    v.b= v.b,
    u=(y.occ.mock>0)-runif(1,0.1,0.8),  #this looks strange -> step(u) is an indicator of whether occupied or not
    #mu.a = matrix(rbinom((n)*J, size=1, prob=1),
    #              nrow=J, ncol=(n)),
    #z = matrix(rbinom((n)*J, size=1, prob=1),
    #           nrow=J, ncol=(n))
    .RNG.seed = .RNG.seed,
    .RNG.name = .RNG.name
  )
  return(out)
}