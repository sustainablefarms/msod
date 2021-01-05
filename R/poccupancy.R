#' @title Occupancy Probability Calculations
#' @param fixedcovar An array of occupancy covariate values. Each row is a model site, each column is a covariate.
#' @param loadfixed An array of loadings for the covariates in 'fixedcovar'. Each row is a species,
#'  each columns is a covariate (in same order as in fixedcovar), and each layer is a draw from the distribution of loadings.
#' @param randomcovar An array of occupancy covariate values that samples the distribution of the covariate.
#' Each row is a model site, each column a covariate, and each layer is an iid draw from the (usually posterior) distribution for the covariate value.
#' @param loadrandom An array of loadings for 'randomcovar'. Each row is a species, each column a covariate, each layer a draw from the covariate distribution.
#' @return An array of occupancy probability values. Each row is a modelsite, each column a species, 
#' each layer a draw corresponding to 'loadfixed', and each layer of the 4th dimension is a draw corresponding to 'loadrandom'.

#' @examples
#' fixedcovar <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5,
#'                      dimnames = list(paste0("Site", LETTERS[1:10]), paste0("C", letters[1:5])))
#' loadfixed <- array(unlist(lapply(seq(0, 12, by = 1), function(x) rnorm(7 * 5, x))), dim = c(7, 5, 12), #each layer has a larger mean
#'                    dimnames = list(paste0("Species", LETTERS[1:7]), paste0("C", letters[1:5]), paste0("D", 1:12)))
#' poccupy_raw.jsodm(fixedcovar, loadfixed)
#' 
#' fitold <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_allhyp_vif_logwoody500m_msnm_year_Time_Wind.rds")
#' fit <- translatefit(fitold)
#' poccupy.jsodm(fit)
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

randomcovar <- array(rnorm(10 * 5), dim = c(10, 5, 30),
                      dimnames = list(paste0("Site", LETTERS[1:10]), paste0("C", letters[1:5]), paste0("D", 1:30)))
loadrandom <- array(unlist(lapply(seq(0, 12, by = 1), function(x) rnorm(7 * 5, x))), dim = c(7, 5, 30), #each layer has a larger mean
                    dimnames = list(paste0("Species", LETTERS[1:7]), paste0("C", letters[1:5]), paste0("D", 1:30)))

poccupy_raw.jsodm_lv <- function(fixedcovar, loadfixed, randomcovar = NULL, loadrandom = NULL, sd_occ_indicator = NULL){
  eta_f <- eta_fixed(fixedcovar, loadfixed)
  
  ndraw <- dim(randomcovar)[[3]]
  eta_rand_l <- lapply(1:ndraw, function(d)
         randomcovar[,,d] %*% t(loadrandom[,,d]))
  eta_rand <- abind::abind(eta_rand_l, along = 3, force.array = TRUE, use.dnns = TRUE)
  dimnames(eta_rand) <- list(dimnames(randomcovar)[[1]], dimnames(loadrandom)[[1]], dimnames(loadrandom)[[3]])
  
  # need to create the 4-array eta. It needs an outer sum across the last two dimensions, but I can't see a way to do it using existing simple functions.
  eta <- array(0, dim = c(dim(eta_f), dim(eta_rand)[[3]]),
               dimnames = c(dimnames(eta_f), dimnames(eta_rand)[3]))
  sitespecindex <- expand.grid(1:nrow(eta_f), 1:ncol(eta_f))
  a <- apply(sitespecindex, MARGIN = 1,
        function(indices){
          vals <- outer(eta_f[indices[[1]], indices[[2]], ],
                eta_rand[indices[[1]], indices[[2]], ])
          return(list(vals))
        })
  for (indexrow in 1:nrow(sitespecindex)){
    eta[sitespecindex[indexrow, 1], sitespecindex[indexrow, 2], , ] <-
      a[[indexrow]][[1]]
  }
  # for (i in 1:nrow(eta_f)){
  #   for (j in 1:ncol(eta_f)){
  #     eta[i, j, , ] <- outer(eta_f[i, j, ], eta_rand[i, j, ], FUN = "+")
  #   }
  # }
  stopifnot(all(eta != 0))
  

  
  pocc <- 1 - pnorm(-eta, mean = 0, sd = 1)
  
  stopifnot(all.equal(dim(pocc), c( #test that final dimensions are correct
    dim(fixedcovar)[[1]],
    dim(loadfixed)[[1]],
    dim(loadfixed)[[3]],
    dim(loadrandom)[[3]]
  )))
  return(pocc)
}
