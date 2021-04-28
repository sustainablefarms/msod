#' @title Prediction Summaries for jsodm_lv Model for New Sites
#' @description Interesting summaries of predictions made by a fitted jsodm_lv.
#'  Includes
#'   occupancy probability of species independent of all other species, for each site, with error bars;
#'   expected species richness for single site;
#' @param Xocc A matrix of unprocessed occupancy covariates. Each row is a site, each column a covariate.
#' @param occ.b An array of occupancy covariate loadings. Each row is a species, each column a covariate, and each layer a draw from the posterior.
#' @examples 
#' fit <- readRDS("../sflddata/private/data/testdata/cutfit_7_4_11_2LV.rds")
#' fit <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_model_2lv_e13.rds")
#' fit <- translatefit(fit)
#' Xocc <- unstandardise.designmatprocess(fit$XoccProcess, fit$data$Xocc[1, , drop = FALSE])
#' pocc <- poccupancy_margotherspecies.jsodm_lv(fit, Xocc)
#' pocc <- poccupancy_mostfavourablesite.jsodm_lv(fit, Xocc)
#' pocc <- poccupancy_randomsite.jsodm_lv(fit, Xocc)
#' sprich1 <- occspecrichness.jsodm_lv(fit, Xocc)
#' sprich <- occspecrichnessRV.jsodm_lv(fit, Xocc)
#' system.time(sprich <- occspecrichness_avsite.jsodm_lv(fit, Xocc))

#' @export
# returns probability of occupancy of each species ignoring other species, for each site
poccupancy_margotherspecies.jsodm_lv <- function(fit, Xocc, toXocc = fit$toXocc){
  Xocc <- sflddata::apply.designmatprocess(fit$XoccProcess, Xocc)
  stopifnot("jsodm_lv" %in% class(fit))
  occ.b <- get_occ_b(fit)
  pocc <- poccupy_raw.jsodm(Xocc, occ.b)
  limits <- hpd_narray(pocc)
  pocc_median <- apply(pocc, MARGIN = c(1, 2), median)
  out <- abind::abind(limits, median = pocc_median)
  return(out)
}

# returns probability of occupancy of a species at the most favourable (according to median) patch, along with corresponding 95% HPD limits of this most favourable patch
#' @export
poccupancy_mostfavourablesite.jsodm_lv <- function(fit, Xocc){
  pocc <- poccupancy_margotherspecies.jsodm_lv(fit, Xocc)
  bestsite <- apply(pocc[,,"median", drop = FALSE], MARGIN = 2, which.max)
  pocc_best <- pocc[1, , ] * 0
  for (i in 1:ncol(pocc)){
    pocc_best[i, ] <- pocc[bestsite[i], i, ]
  }
  out <- cbind(pocc_best, bestsite)
  return(out)
}

# returns probability of occupancy in a randomly selected site
#' @export
poccupancy_randomsite.jsodm_lv <- function(fit, Xocc){
  Xocc <- sflddata::apply.designmatprocess(fit$XoccProcess, Xocc)
  stopifnot("jsodm_lv" %in% class(fit))
  occ.b <- get_occ_b(fit)
  pocc <- poccupy_raw.jsodm(Xocc, occ.b)
  
  Epocc <- apply(pocc, MARGIN = c(1, 2), mean)
  Vpocc <- apply(pocc, MARGIN = c(1, 2), sd)^2
  randselsiteEV <- EVtheta2EVmarg(Vsum = Vpocc, Esum = Epocc)
  colnames(randselsiteEV) <- colnames(pocc)
  rownames(randselsiteEV) <- c("E", "V")
  lower <- randselsiteEV["E", ] - 2 * sqrt(randselsiteEV["V", ])
  upper <- randselsiteEV["E", ] + 2 * sqrt(randselsiteEV["V", ])
  out <- rbind(randselsiteEV, lower = lower, upper = upper)
  if ( any(lower < 0) || any(upper > 1)){warning("The Gaussian approximation for upper and lower limits is not appropriate: some limits are outside 0 and 1")}
  return(out)
}
  

# returns species richness predicted mean and variance
#' @export
occspecrichness.jsodm_lv <- function(fit, Xocc){
  stopifnot("jsodm_lv" %in% class(fit))
  specrich <- apply_to_new_data(speciesrichness, fit = fit, Xocc = Xocc,  
                                funargs = list(occORdetection = "occupancy",
                                               nlvperdraw = 5))
  warning("Using 5 simulated LV values per draw")
  return(specrich)
}

#' @export
occspecrichnessRV.jsodm_lv <- function(fit, Xocc){
  stopifnot("jsodm_lv" %in% class(fit))
  specrich <- predsumspeciesRV_newdata(fit, Xocc, nLVsim = 100, type = "marginal")
  return(specrich)
}

#' @export
# returns species richness predicted mean and variance when the site is chosen randomly with equal probability
occspecrichness_avsite.jsodm_lv <- function(fit, Xocc){
  specrich <- occspecrichness.jsodm_lv(fit, Xocc)
  Epred <- mean(specrich["E",])
  m2_site <- specrich["V", ] + specrich["E",]^2
  Em2 <- mean(m2_site)
  Vpred <- Em2 - Epred^2
  return(c(
    E = Epred,
    V = Vpred
  ))
}

# given an n array, with one dimension representing draws, compute the hpd interval. Returns an array with the dimension representing draws removed.
hpd_narray <- function(arr, prob = 0.95, drawdim = length(dim(arr))){
  inputdim <- 1:length(dim(arr))
  inputdimsizes <- dim(arr)
  arr_t <- aperm(arr, perm = c(inputdim[inputdim != drawdim], drawdim)) #reshape so that draw dimension is the final dimension
  dim(arr_t) <- c(prod(inputdimsizes[inputdim != drawdim]), inputdimsizes[drawdim]) # each column is a draw, rows go through other dimensions in order
  mcmcobj <- coda::as.mcmc(t(arr_t))
  intervals <- coda::HPDinterval(mcmcobj, prob = prob)
  out <- array(intervals, dim = c(inputdimsizes[inputdim != drawdim], 2),
        dimnames = c(dimnames(arr)[inputdim != drawdim], list(limit = c("lower", "upper"))))
  return(out)
}

