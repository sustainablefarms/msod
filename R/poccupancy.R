#' @title Occupancy Probability Calculations
#' @description Predicts occupancy probability given a draw of loadings and (random) covariate loadings.
#' The draw would typically be either (a) from the full posterior of a fitted model, which includes fitted values for the random covariate loadings
#' or (b) drawn from the posterior of the loadings with random covariates simulated.
#' @param fixedcovar An array of occupancy covariate values. Each row is a model site, each column is a covariate.
#' @param loadfixed An array of loadings for the covariates in 'fixedcovar'. Each row is a species,
#'  each columns is a covariate (in same order as in fixedcovar), and each layer is a draw from the distribution of loadings.
#' @param randomcovar An array of occupancy covariate values that samples the distribution of the covariate.
#' Each row is a model site, each column a covariate, and each layer is a draw.
#' @param loadrandom An array of loadings for 'randomcovar'. Each row is a species, each column a covariate, each layer a draw from the covariate distribution, and must be the same draw as loadfixed.
#' @return An array of occupancy probability values. Each row is a modelsite, each column a species, 
#' each layer a draw corresponding to the loadings.

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
#' randomcovar <- array(rnorm(10 * 2), dim = c(10, 2, 12),
#'                       dimnames = list(paste0("Site", LETTERS[1:10]), paste0("C", letters[1:2]), paste0("RcD", 1:12)))
#' loadrandom <- array(unlist(lapply(seq(0, 0.3, by = 0.01), function(x) rnorm(7 * 2, x, sd = 0.01))), dim = c(7, 2, 12), #each layer has a larger mean
#'                     dimnames = list(paste0("Species", LETTERS[1:7]), paste0("C", letters[1:2]), paste0("LrD", 1:12)))
#' pocc <- poccupy_raw.jsodm_lv(fixedcovar, loadfixed, randomcovar, loadrandom)
#' model2lv <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_model_2lv_e13.rds")
#' model2lv_new <- translatefit(model2lv)
#' pocc <- poccupy_raw.jsodm_lv(fixedcovar, loadfixed, randomcovar, loadrandom)
#' pocc <- poccupy.jsodm_lv(model2lv_new, fullposterior = FALSE)
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

#' @export
poccupy_raw.jsodm_lv <- function(fixedcovar, loadfixed, randomcovar, loadrandom){
  stopifnot(dim(loadfixed)[[3]] == dim(loadrandom)[[3]]) # stop if ndraw of loadrandom differs from loadfixed, this means the draws of each can't be tied together
  stopifnot(dim(loadfixed)[[3]] == dim(randomcovar)[[3]]) # all are draws are from the same joint distribution of fitted parameter x latent variable value.
  # It is either the full posterior distribution or the posterior of the loadings with independent latent variable values
  eta_f <- eta_fixed(fixedcovar, loadfixed)
  
  sd_occ_indicator <- sqrt(1 - arr3_sumalong2(loadrandom^2))

  eta_rand_l <- lapply(1:dim(loadrandom)[[3]], function(d)
           randomcovar[,,d] %*% t(loadrandom[,,d]))
  names(eta_rand_l) <- dimnames(loadrandom)[[3]]
  eta_rand <- abind_alongNp1(eta_rand_l)
  #slower version:
  #   eta_rand <- abind::abind(eta_rand_l, along = 3, force.array = TRUE, use.dnns = TRUE)
  #   dimnames(eta_rand) <- list(dimnames(randomcovar)[[1]], dimnames(loadrandom)[[1]], dimnames(loadrandom)[[3]])
  
  # for each draw combine eta_f and eta_rand 
  eta <- eta_f + eta_rand
  
  # standardise by sd_occ_indicator
  eta_s <- apply(eta, MARGIN = 1, function(x) x / sd_occ_indicator)
  dim(eta_s) <- dim(eta)[c(2, 3, 1)]
  eta_s <- aperm(eta_s, perm = c(3, 1, 2))
  dimnames(eta_s) <- dimnames(eta)
  
  # check conversion
  # eta_t <- apply(eta, MARGIN = 1, function(x) {x / (1 + 0 * sd_occ_indicator)})
  # dim(eta_t) <- dim(eta)[c(2, 3, 1)]
  # eta_t <- aperm(eta_t, perm = c(3, 1, 2))
  # stopifnot(all(eta_t == eta))
  
  pocc <- 1 - pnorm(-eta_s, mean = 0, sd = 1)
  
  stopifnot(all.equal(dim(pocc), c( #test that final dimensions are correct
    dim(fixedcovar)[[1]],
    dim(loadfixed)[[1]],
    dim(loadfixed)[[3]]
  )))
  return(pocc)
}

# equivalent, but hopefully faster than apply(arr, MARGIN = c(1, 3), sum)
# stopifnot(all.equal(apply(arr, MARGIN = c(1, 3), sum), arr3_sumalong2(arr)))
arr3_sumalong2 <- function(arr){
  stopifnot(length(dim(arr))==3)
  arr2 <- aperm(arr, perm = c(1, 3, 2))
  dim(arr2) <- c(dim(arr)[[1]] * dim(arr)[[3]], dim(arr)[[2]])
  colsum2 <- Rfast::rowsums(arr2)
  dim(colsum2) <- dim(arr)[c(1, 3)]
  colsum <- colsum2
  dimnames(colsum) <- list(dimnames(arr)[[1]], dimnames(arr)[[3]])
  return(colsum)
}

abind_alongNp1 <- function(listofarrays){
  a <- unlist(listofarrays)
  dim(a) <- c(dim(listofarrays[[1]])[c(1,2)],
              length(listofarrays))
  dimnames(a)[[1]] <- dimnames(listofarrays[[1]])[[1]]
  dimnames(a)[[2]] <- dimnames(listofarrays[[1]])[[2]]
  dimnames(a)[[3]] <- names(listofarrays)
  return(a)
}

# default treats lv.v and loadings as independent, simulating the former. 
#' lvfromposterior = TRUE draws lv.v from the posterior with all other loadings.
# margLV = TRUE returns probabilities marginal over fitted lv.b and unknown lv.v. Marginalising over lv.v (non-fitted) is equivalent to the plain jsodm occupancy probability given the occ.b loadings.
#' @export
poccupy.jsodm_lv <- function(fit, usethetasummary = NULL, lvfromposterior = TRUE, margLV = FALSE){
  occ.v <- fit$data$Xocc
  dimnames(occ.v) <- list(ModelSite = rownames(occ.v), Covariate = colnames(occ.v))
  
  if (is.null(usethetasummary)){occ.b <- get_occ_b(fit)}
  else {
    occ.b <- get_theta(fit, type = type)
  }
  
  if (margLV){
    return(poccupy_raw.jsodm(occ.v, occ.b))
  }
  
  lv.b <- get_lv_b(fit)
  if (!fullposterior){
    lv.v <- array(rnorm(dim(occ.v)[[1]] * dim(lv.b)[[2]] *  dim(lv.b)[[3]]), 
                  dim = c(dim(occ.v)[[1]], dim(lv.b)[[2]],  dim(lv.b)[[3]]),
                  dimnames = list(ModelSite = rownames(occ.v),
                                  LV = paste0("lv", 1:dim(lv.b)[[2]], ".v"),
                                  Draw = 1:dim(lv.b)[[3]]))
  } else {
    lv.v <- get_lv_v(fit)
  }
  pocc <- poccupy_raw.jsodm_lv(occ.v, occ.b, lv.v, lv.b)
  return(pocc)
}


