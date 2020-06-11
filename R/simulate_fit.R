# Simulate observations from parameters
# source("./R/calcpredictions.R")

#' @examples 
#' fit <- readRDS("./tmpdata/deto_wind.rds")
#' fit$data <- as.list.format(fit$data)
#' detected <- simulate.fit(fit, esttype = "median")
#' 
#' ## Use to check residuals
#' fit_sim <- fit
#' fit_sim$data$y <- detected
#' source("./R/DS_residuals.R")
#' source("./R/DS_residuals_plots.R")
#' resid_det <- ds_detection_residuals.fit(fit_sim, type = 100, seed = 321)
#' resid_occ <- ds_occupancy_residuals.fit(fit_sim, type = "median", seed = 123)
#' ## These residuals *should* be perfectly normally distributed for simulated data.

#' resid_occ <- ds_occupancy_residuals.fit(fit_sim, seed = 123)
#' plot_residuals_detection.fit(fit_sim, 
#'                    detectionresiduals = resid_det,
#'                    varidx = 2,
#'                    plotfunction = facet_covariate) +
#' coord_cartesian(ylim = c(-1, 1))
#' plot_residuals_detection.fit(fit_sim, 
#'                    detectionresiduals = resid_det,
#'                    varidx = 2) +
#' coord_cartesian(ylim = c(-1, 1))



simulate.fit <- function(fit, esttype = "median", conditionalLV = TRUE){
  poccupy <- poccupy_species(fit, type = esttype, conditionalLV = conditionalLV)
  pdetectcond <- pdetect_condoccupied(fit, type = esttype)
  fit$data <- as.list.format(fit$data)
  
  occupied <- apply(poccupy, c(1, 2), function(x) rbinom(1, 1, x))
  # array(NA, dim = c(replicates, dim(poccupy)[[1]], dim(poccupy)[[2]]),
  #       dimnames = list(replicates = paste0("replicate", 1:replicates),
  #                       modelsite = 1:nrow(poccupy),
  #                       species = colnames(poccupy)))
  
  pdetect <- pdetectcond * occupied[fit$data$ModelSite, ]
  detected <- apply(pdetect, c(1, 2), function(x) rbinom(1, 1, x))
  
  return(detected)
}

## function for creating a fully artificial fitted object
#' @return A list that has enough similarities to runjags objects that residual calculations are possible.
#' The true parameter set is the first (and only row) of the first MCMC chain.
#' It can be accessed using get_theta(fit, type = 1)
#' @param nspecies Number of species
#' @param nsites Number of sites
#' @param nvisitspersite Number of visits per site
#' @param nlv Number of latent variables. Must be 4 or less
#' @param OccFmla Formula for occupancy. Available variables: UpSite, Sine1 and Sine2
#' @param ObsFmla Formula for detection. Available variables: Upvisit, Step
#' @param u.b.min, u.b.max, v.b.min, v.b.max The upper and lower bouonds of the u.b and v.b parameters.
#'  May be a single number or an array with rows corresponding to species and columns to covariates.
#'  @examples 
#'  artfit <- artificial_runjags(nspecies = 2, nsites = 10, nvisitspersite = 4, nlv = 2)
artificial_runjags <- function(nspecies = 4, nsites = 100, nvisitspersite  = 2, nlv = 2,
                               OccFmla = "~ UpSite + Sine1 + Sine2",
                               ObsFmla = "~ UpVisit + Step",
                               u.b.min = -1,
                               u.b.max = 1,
                               v.b.min = -1,
                               v.b.max = 1,
                               lv.coef.min = -0.5,
                               lv.coef.max = 0.5
                               ){
  species <- LETTERS[1:nspecies]
  covardfs <- simulate_covar_data(nsites, nvisitspersite)

  LV <- scale(cbind(sites %% 2,
              sites < (nsites / 2), 
              sites > (nsites / 5),
                    ((sites %/% 5) * 5 == sites ) | (sites %/% 3) * 3 == sites))
  if (nlv == 0) {LV <- NULL}
  else {LV <- LV[, 1:nlv, drop = FALSE]}
  XoccProcess <- prep.designmatprocess(covardfs$Xocc, OccFmla)
  XobsProcess <- prep.designmatprocess(covardfs$Xobs, ObsFmla)

  data.list <- prep.data(covardfs$Xocc, yXobs = covardfs$Xobs,
                         ModelSite = "ModelSite", 
                         species = NULL,
                         nlv = nlv, 
                         XoccProcess = XoccProcess,
                         XobsProcess = XobsProcess)
  fit <- list()
  fit$data <- data.list
  fit$data$n <- length(species)
  fit$species <- species
  fit$summary.available <- TRUE

  # set parameters
  u.b <- matrix(runif( fit$data$n * fit$data$Vocc, min = u.b.min, max = u.b.max), nrow = fit$data$n, ncol = fit$data$Vocc, byrow = FALSE)
  v.b <- matrix(runif(  fit$data$n * fit$data$Vobs, min = v.b.min, max = v.b.max), nrow = fit$data$n, ncol = fit$data$Vobs, byrow = FALSE)
  theta <- c(matrix2bugsvar(u.b, name = "u.b"),
             matrix2bugsvar(v.b, name = "v.b"))
  if (nlv > 0){
    lv.coef <- matrix(runif(  fit$data$n * fit$data$nlv, min = lv.coef.min, max = lv.coef.max), nrow = fit$data$n, ncol = fit$data$nlv) #0.5 constraint makes sure rowSum(lv.coef^2) < 1
    theta <- c(theta, 
               matrix2bugsvar(lv.coef, name = "lv.coef"),
               matrix2bugsvar(LV, name = "LV"))
  }
  fit$mcmc <- list()
  fit$mcmc[[1]] <- t(as.matrix(theta))

  # simulate data
  fit$data$y <- simulate.fit(fit, esttype = 1, conditionalLV = (nlv > 0) )
  colnames(fit$data$y) <- species
  return(fit)
}

#' @param nsites Number of ModelSites to simulate
#' @param nvisitspersite Number of visits per site (the same number for each site)
#' @describeIn artificial_runjags Simulates covariate dataframes for occupancy and detection.
#' Occupancy covariate names are UpSite, Sine1 and Sine2.
#' Detection dovariate names are UpVisit and Step.
#' @return A list with elements Xocc, and Xobs for the occupancy and detection covariates respectively
simulate_covar_data <- function(nsites, nvisitspersite){
  sites <- c(1:nsites)
  XoccIn <- data.frame(ModelSite = sites,
                     UpSite = sites,
                     Sine1 = sin(2 * pi * sites / nsites),
                     Sine2 = sin(4 * pi * sites / nsites))
  XobsIn <- data.frame(ModelSite = rep(sites, nvisitspersite),
                     UpVisit = 1:(nvisitspersite*nsites),
                     Step = c(rep(0, floor(nvisitspersite*nsites / 2)),
                              rep(1, ceiling(nvisitspersite*nsites / 2)))
                     )
  return(list(Xocc = XoccIn, Xobs = XobsIn))
}

#' @param nspecies is the number of species to simulate. Species are named A, B, C... (max of 26 allowed)
#' @param nvisits Number of visits in total to simulate
#' @param p The probability of detection, constant for all species and visits.
#' @return A simulated dataframe of detections. Column names are the species names, each row is a visit.
#' Elements of the data frame are either 1 (detected) or 0 (not detected).
simulate_iid_detections <- function(nspecies, nvisits, p = 0.5){
  stopifnot(nspecies <= 26)
  species <- LETTERS[1:nspecies]
  sim_y <- as.data.frame(matrix(rbinom(length(species) * nvisits, 1, p), ncol = length(species)))
  colnames(sim_y) <- species
  return(sim_y)
}
