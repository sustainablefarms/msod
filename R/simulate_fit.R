# Simulate observations from parameters
# source("./R/calcpredictions.R")

#' @examples 
#' fit <- readRDS("./tmpdata/deto_wind.rds")
#' fit$data <- as_list_format(fit$data)
#' class(fit) <- c("jsodm", class(fit))
#' detected <- simulate_detections(fit, esttype = "median", UseFittedLV = TRUE)
#' 


#' @title Simulate observations from parameters of a fitted object.
#' @details Simulate detections from all given parameters (including LV and random effects)
#' @param fit A runjags fitted object created by [run.detectionoccupancy()]
#' @param esttype Specifies parameter set to extract from fit. See [get_occ_b()] and similar for specification.
#' @export
simulate_detections <- function(fit, esttype = "median", seed = NULL){
  fit$data <- as_list_format(fit$data)
  pocc <- poccupy(fit, usethetasummary = esttype)
  if (inherits(fit, "jsodm_lv")) {
    pocc <- poccupy(fit, usethetasummary = esttype, lvvfromposterior = TRUE)
  } else if (inherits(fit, "jsodm")) {
    pocc <- poccupy(fit, usethetasummary = esttype)}
  pdet_occ <- pdet_occ(fit, usethetasummary = esttype)
  
  # simulate occupied
  set.seed(seed)
  occupied <- rmanybern(pocc)
  # dim(occupied) <- dim(pocc)
  # array(NA, dim = c(replicates, dim(poccupy)[[1]], dim(poccupy)[[2]]),
  #       dimnames = list(replicates = paste0("replicate", 1:replicates),
  #                       modelsite = 1:nrow(poccupy),
  #                       species = colnames(poccupy)))
  
  pdetect <- pdet_occ * occupied[fit$data$ModelSite, , , drop = FALSE]
  
  if (!is.null(seed)){set.seed(seed + 10)}
  detected <- rmanybern(pdetect)
  return(drop_to_matrix(detected))
}

simulate_lv.v <- function(fit, replaceinsitu = FALSE){
  simLV <- matrix(rnorm(fit$data$nlv * nrow(fit$data$Xocc)), ncol = fit$data$nlv)
  simLVbugsname <- matrix2bugsvar(simLV, name = "lv.v")
  if (!replaceinsitu) {return(simLVbugsname)}
  
  fit$mcmc <- lapply(fit$mcmc, function(arr) {
    arr[, names(simLVbugsname)] <- Rfast::rep_row(simLVbugsname, nrow(arr))
    return(arr)
  })
  return(fit)
}

#' @describeIn simulate_detections Simulate LV values and detections. 
#' Uses [simulate_detections()] and [simulate_lv.v()].
#' For each site the simulated LV values are copied across all draws.
#' @export
simulate_detections_lv.v <- function(fit, esttype = "median", seed = NULL){
  if (!is.null(seed)){set.seed(seed + 20)}
  fit <- simulate_lv.v(fit, replaceinsitu = TRUE)
  detected <- simulate_detections(fit, esttype = esttype, seed = seed)
  return(detected)
}




#' @title Create a fully artificial fitted object
#' @return A list that has enough similarities to runjags objects that residual calculations are possible.
#' The true parameter set is the first (and only row) of the first MCMC chain.
#' It can be accessed using get_theta(fit, type = 1)
#' @param nspecies Number of species
#' @param nsites Number of sites
#' @param nvisitspersite Number of visits per site
#' @param nlv Number of latent variables. Must be 4 or less
#' @param OccFmla Formula for occupancy. Available variables: UpSite, Sine1 and Sine2
#' @param ObsFmla Formula for detection. Available variables: Upvisit, Step
#' @param occ.b.min, occ.b.max, det.b.min, det.b.max The upper and lower bouonds of the occ.b and det.b parameters.
#'  May be a single number or an array with rows corresponding to species and columns to covariates.
#' @param lv.b.min, lv.b.max Same as occ.b.min and occ.b.max for the latent variable loadings.
#' @examples 
#' artfit <- artificial_runjags(nspecies = 2, nsites = 20, nvisitspersite = 4, 
#' modeltype = "jsodm_lv", nlv = 2, seed = NULL)
#' \# with high correlation between occupancy of species
#' artfit <- artificial_runjags(nspecies = 2, nsites = 10, nvisitspersite = 4,
#'                               OccFmla = "~ 1",
#'                               occ.b.min = 0.8,
#'                               lv.b.min = 0.3,
#'                               modeltype = "jsodm_lv",
#'                               nlv = 2)
#'  cor(artfit$data$y)
#' @export
artificial_runjags <- function(nspecies = 4, nsites = 100, nvisitspersite  = 2,
                               OccFmla = "~ UpSite + Sine1 + Sine2",
                               ObsFmla = "~ UpVisit + Step",
                               occ.b.min = -1,
                               occ.b.max = 1,
                               det.b.min = -1,
                               det.b.max = 1,
                               lv.b.min = -0.5,
                               lv.b.max = 0.5,
                               modeltype = "jsodm_lv",
                               seed = NULL,
                               ...
                               ){
  stopifnot(modeltype %in% availmodeltypes)
  
  # first step is to populate predictor values (including LV) and parameter values
  species <- make.names(rep(LETTERS, ceiling(nspecies / length(LETTERS))), unique = TRUE)[1:nspecies]
  covardfs <- artificial_covar(nsites, nvisitspersite, OccFmla = OccFmla, ObsFmla = ObsFmla)
  
  data.list <- prepJAGSdata2(modeltype,
                             Xocc = covardfs$Xocc,
                             Xobs = covardfs$Xobs,
                             y = NULL,
                             ModelSite = covardfs$ModelSite,
                             ...)
  
  fit <- list()
  fit$data <- data.list
  fit$data$nspecies <- length(species)
  fit$species <- species
  fit$ModelSite <- "ModelSite"
  fit$summary.available <- TRUE
  fit$sample <- 1

  # set parameters
  if (!is.null(seed)){set.seed(seed)}
  occ.b <- matrix(runif( fit$data$nspecies * fit$data$noccvar, min = occ.b.min, max = occ.b.max), nrow = fit$data$nspecies, ncol = fit$data$noccvar, byrow = FALSE)
  if (!is.null(seed)){set.seed(seed + 1)}
  det.b <- matrix(runif(  fit$data$nspecies * fit$data$nobsvar, min = det.b.min, max = det.b.max), nrow = fit$data$nspecies, ncol = fit$data$nobsvar, byrow = FALSE)
  theta <- c(matrix2bugsvar(occ.b, name = "occ.b"),
             matrix2bugsvar(det.b, name = "det.b"))
  
  if (modeltype == "jsodm_lv"){
    nlv <- list(...)$nlv
    sites <- 1:nrow(covardfs$Xocc)
    if ((nsites <= 10) && (nlv >= 2)) {stop("More than 10 sites required to use 3rd LV.")}
    LV <- scale(cbind(sites %% 2,
                sites < (nsites / 2), 
                sites < (nsites / 10),
                      ((sites %/% 5) * 5 == sites ) | (sites %/% 3) * 3 == sites))
    if (nlv == 0) {LV <- NULL}
    else {LV <- LV[, 1:nlv, drop = FALSE]}
  if (!is.null(seed)){set.seed(seed + 2)}
    lv.b <- matrix(runif(  fit$data$nspecies * fit$data$nlv, min = lv.b.min, max = lv.b.max), nrow = fit$data$nspecies, ncol = fit$data$nlv) #0.5 constraint makes sure rowSum(lv.b^2) < 1
    theta <- c(theta, 
               matrix2bugsvar(lv.b, name = "lv.b"),
               matrix2bugsvar(LV, name = "lv.v"))
  }
  fit$mcmc <- list()
  fit$mcmc[[1]] <- t(as.matrix(theta))
  
  class(fit) <- c(modeltype, class(fit))

  # simulate data using the LV values given above
  fit$data$y <- simulate_detections(fit, esttype = 1, seed = seed)
  colnames(fit$data$y) <- species
  
  ellipsis::check_dots_used()
  return(fit)
}

#' @describeIn artificial_runjags Generate fake covariate data.
#' @param nsites Number of ModelSites to simulate
#' @param nvisitspersite Number of visits per site (the same number for each site)
#' @details 
#' Occupancy covariate names are UpSite, Sine1 and Sine2.
#' Detection dovariate names are UpVisit and Step.
#' @return A list with elements Xocc, and Xobs for the occupancy and detection covariates respectively
#' @export
artificial_covar_data <- function(nsites, nvisitspersite){
  warning("This function will be obsolete in next version.")
  out <- artificial_covar(nsites = nsites, nvisitspersite = nvisitspersite)
  Xocc <- cbind(ModelSite = 1:nrow(out$Xocc), out$Xocc)
  Xobs <- cbind(ModelSite = out$ModelSite, out$Xobs)
  return(list(Xocc = Xocc, Xobs = Xobs))
}

artificial_covar <- function(nsites, nvisitspersite, OccFmla = "~ .", ObsFmla = "~ ."){
  sites <- c(1:nsites)
  XoccIn <- data.frame(
                     UpSite = sites,
                     Sine1 = 10 * sin(2 * pi * sites / nsites) + 100,
                     Sine2 = sin(4 * pi * sites / nsites))
  XoccIn <- scale(XoccIn)
  Xocc <- model.matrix(as.formula(OccFmla), data = data.frame(XoccIn))
  XobsIn <- data.frame(ModelSite = rep(sites, nvisitspersite),
                     UpVisit = 1:(nvisitspersite*nsites),
                     Step = c(rep(0, floor(nvisitspersite*nsites / 2)),
                              rep(1, ceiling(nvisitspersite*nsites / 2)))
                     )
  ModelSite <- XobsIn$ModelSite
  XobsIn <- scale(XobsIn[, -1])
  Xobs <- model.matrix(as.formula(ObsFmla), data = data.frame(XobsIn))
  return(list(Xocc = Xocc, Xobs = Xobs, ModelSite = ModelSite))
}

#' @describeIn artificial_runjags Generate fake observation data.
#' @details Every species is equally likely to be detected at every visit.
#' @param nspecies is the number of species to simulate. Species are named A, B, C... (max of 26 allowed)
#' @param nvisits Number of visits in total to simulate
#' @param p The probability of detection, constant for all species and visits.
#' @return A simulated dataframe of detections. Column names are the species names, each row is a visit.
#' Elements of the data frame are either 1 (detected) or 0 (not detected).
#' @export
simulate_iid_detections <- function(nspecies, nvisits, p = 0.5){
  stopifnot(nspecies <= 26)
  species <- LETTERS[1:nspecies]
  sim_y <- as.data.frame(matrix(rbinom(length(species) * nvisits, 1, p), ncol = length(species)))
  colnames(sim_y) <- species
  return(sim_y)
}
