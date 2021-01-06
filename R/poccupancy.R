#' @title Occupancy Probability Calculations
#' @param fixedcovar An array of occupancy covariate values. Each row is a model site, each column is a covariate.
#' @param loadfixed An array of loadings for the covariates in 'fixedcovar'. Each row is a species,
#'  each columns is a covariate (in same order as in fixedcovar), and each layer is a draw from the distribution of loadings.
#' @param randomcovar An array of occupancy covariate values that samples the distribution of the covariate.
#' Each row is a model site, each column a covariate, and each layer is an iid draw from the (usually posterior) distribution for the covariate value.
#' @param loadrandom An array of loadings for 'randomcovar'. Each row is a species, each column a covariate, each layer a draw from the covariate distribution, and must be the same draw as loadfixed.
#' @return An array of occupancy probability values. Each row is a modelsite, each column a species, 
#' each layer a draw corresponding to the loadings, and each layer of the 4th dimension is a draw corresponding to 'randomcovar'.

#' @examples
#' fixedcovar <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5,
#'                      dimnames = list(paste0("Site", LETTERS[1:10]), paste0("C", letters[1:5])))
#' loadfixed <- array(unlist(lapply(seq(0, 12, by = 1), function(x) rnorm(7 * 5, x))), dim = c(7, 5, 12), #each layer has a larger mean
#'                    dimnames = list(paste0("Species", LETTERS[1:7]), paste0("C", letters[1:5]), paste0("LfD", 1:12)))
#' poccupy_raw.jsodm(fixedcovar, loadfixed)
#' 
#' fitold <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_allhyp_vif_logwoody500m_msnm_year_Time_Wind.rds")
#' fit <- translatefit(fitold)
#' poccupy.jsodm(fit)
#' 
#' 
#' randomcovar <- array(rnorm(10 * 2), dim = c(10, 2, 30),
#'                       dimnames = list(paste0("Site", LETTERS[1:10]), paste0("C", letters[1:2]), paste0("RcD", 1:30)))
#' loadrandom <- array(unlist(lapply(seq(0, 0.3, by = 0.01), function(x) rnorm(7 * 2, x, sd = 0.01))), dim = c(7, 2, 12), #each layer has a larger mean
#'                     dimnames = list(paste0("Species", LETTERS[1:7]), paste0("C", letters[1:2]), paste0("LrD", 1:12)))
#' pocc <- poccupy_raw.jsodm_lv(fixedcovar, loadfixed, randomcovar, loadrandom)
#' model2lv <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_model_2lv_e13.rds")
#' model2lv_new <- translatefit(model2lv)
#' pocc <- poccupy_indlvv.jsodm_lv(fixedcovar, loadfixed, randomcovar, loadrandom)
#' 
#' @export
poccupy_raw.jsodm <- function(fixedcovar, loadfixed, randomcovar = NULL, loadrandom = NULL){
  stopifnot(is.null(randomcovar))
  stopifnot(is.null(loadrandom))
  
  eta <- eta_fixed(fixedcovar, loadfixed)
  pocc <- 1 - pnorm(-eta, mean = 0, sd = 1)
  return(pocc)
}

eta_fixed <- function(fixedcovar, loadfixed){
  stopifnot(length(dim(fixedcovar)) == 2)
  stopifnot(length(dim(loadfixed)) <= 3)
  stopifnot(dim(fixedcovar)[[2]] == dim(loadfixed)[[2]])
  
  eta <- tensor::tensor(fixedcovar, loadfixed, alongA = 2, alongB = 2)

  stopifnot(all.equal(dim(eta), c( #test that final dimensions are correct
    dim(fixedcovar)[[1]],
    dim(loadfixed)[[1]],
    dim(loadfixed)[[3]]
  )))
  return(eta)
}


#' @export
poccupy.jsodm <- function(fit){
  fixedcovar <- fit$data$Xocc
  dimnames(fixedcovar) <- list(ModelSite = rownames(fixedcovar), Covariate = colnames(fixedcovar))
  
  loadfixed <- get_occ_b(fit)
  pocc <- poccupy_raw.jsodm(fixedcovar, loadfixed)
  return(pocc)
}



# if no randomness included then it is assumed no information about the LV values or loadings is known and the prediction becomes a plain jsodm prediction
# MISTAKE: loadings of fixed and lv.v must come from the same draw, that means fixedcovar and loadrandom must be of the same final dimension
# basically want to marginalise randomcovar without marginalising loadrandom
# the situation where lv.b and lv.v are tied together is fixedcovar, and loadfixed, but with a changed sd of occ_indicator
#sd_occ_indicator must be a matrix with rows that are species and columns that are draws
poccupy_posterior_lvv.jsodm_lv <- function(fixedcovar, loadfixed, sd_occ_indicator){
  eta_f <- eta_fixed(fixedcovar, loadfixed)
  #for each column and each drow, scale eta_f by sd_occ_indicator
  for (modelsite in 1:nrow(eta_f)){
    eta_f[modelsite, , ] <- eta_f[modelsite, , ] / sd_occ_indicator
  }
  pocc <- 1 - pnorm(-eta_f, mean = 0, sd = 1)
  return(pocc)
}

poccupy_indlvv.jsodm_lv <- function(fixedcovar, loadfixed, randomcovar, loadrandom){
  stopifnot(dim(loadfixed)[[3]] == dim(loadrandom)[[3]]) # stop if ndraw of loadrandom differs from loadfixed, this means the draws of each can't be tied together
  eta_f <- eta_fixed(fixedcovar, loadfixed)
  
  sd_occ_indicator <- apply(loadrandom, MARGIN = c(1, 3), function(x) sqrt(1- sum(x^2)))

  eta_rand <- tensor::tensor(randomcovar, loadrandom, alongA = 2, alongB = 2)
  # eta_rand has dimensions corresponding to site x random covariate draw x species x load draw
  eta_rand <- aperm(eta_rand, perm = c(1, 3, 4, 2)) # site x species x load draw x random covariate draw

  # add eta_f, which is the same for each random covariate draw (dimension 2)
  eta <- apply(eta_rand, MARGIN = 4, function(x) {x + eta_f})
  # 3 array is converted to vector. I'm guessing it does this by going down row (fixed column and layer), dim(eta) corrects this in the right automatic fashion!
  dim(eta) <- dim(eta_rand)
  dimnames(eta) <- dimnames(eta_rand)
  
  # check conversion
  eta_t <- apply(eta_rand, MARGIN = 4, function(x) {x + eta_f * 0})
  dim(eta_t) <- dim(eta_rand)
  stopifnot(all(eta_t == eta_rand))
  
  # for each species for each load draw, eta needs to be divided by sd_occ_indicator
  eta_s <- apply(eta, MARGIN = c(1, 4), function(x) {x / sd_occ_indicator})
  # 3 array is converted to vector. I'm guessing it does this by going down row (fixed column and layer), dim(eta) corrects this in the right automatic fashion!
  dim(eta_s) <- dim(eta)[c(2, 3, 1, 4)]
  eta_s <- aperm(eta_s, perm = c(3, 1, 2, 4))
  dimnames(eta_s) <- dimnames(eta_rand)
  
  # check conversion
  eta_t <- apply(eta, MARGIN = c(1, 4), function(x) {x / (1 + 0 * sd_occ_indicator)})
  # 3 array is converted to vector. I'm guessing it does this by going down row (fixed column and layer), dim(eta) corrects this in the right automatic fashion!
  dim(eta_t) <- dim(eta)[c(2, 3, 1, 4)]
  eta_t <- aperm(eta_t, perm = c(3, 1, 2, 4))
  all(eta_t == eta)
  
  pocc <- 1 - pnorm(-eta, mean = 0, sd = 1)
  
  stopifnot(all.equal(dim(pocc), c( #test that final dimensions are correct
    dim(fixedcovar)[[1]],
    dim(loadfixed)[[1]],
    dim(loadfixed)[[3]],
    dim(randomcovar)[[3]]
  )))
  return(pocc)
}

poccupy.jsodm_lv <- function(fit, simlvnum = 100){
  occ.v <- fit$data$Xocc
  dimnames(occ.v) <- list(ModelSite = rownames(occ.v), Covariate = colnames(occ.v))
  
  occ.b <- get_occ_b(fit)
  lv.b <- get_lv_b(fit)
  lv.v <- array(rnorm(dim(occ.v)[[1]] * dim(lv.b)[[2]] * simlvnum), 
                dim = c(dim(occ.v)[[1]], dim(lv.b)[[2]], simlvnum),
                dimnames = list(ModelSite = rownames(occ.v),
                                lv.v = paste0("lv.v", 1:dim(lv.b)[[2]]),
                                lvdraw = 1:simlvnum))
  pocc <- poccupy_indlvv.jsodm_lv(occ.v, occ.b, lv.v, lv.b)
  return(pocc)
}
