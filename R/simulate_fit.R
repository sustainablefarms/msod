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
#' @param esttype Specifies parameter set to extract from fit, see [get_theta()]
#' @param UseFittedLV Logical. If TRUE, the simulation uses fitted latent variable values.
#'  If FALSE, latent variable values will be simulated for each ModelSite
#' @export
simulate_detections <- function(fit, esttype = "median"){
  fit$data <- as_list_format(fit$data)
  if (inherits(fit, "jsodm_lv")) {poccupy <- poccupy_species(fit, type = esttype, conditionalLV = TRUE)
  } else if (inherits(fit, "jsodm")) {poccupy <- poccupy_species(fit, type = esttype, conditionalLV = FALSE)}
  pdetectcond <- pdetect_condoccupied(fit, type = esttype)
  
  occupied <- apply(poccupy, c(1, 2), function(x) rbinom(1, 1, x))
  # array(NA, dim = c(replicates, dim(poccupy)[[1]], dim(poccupy)[[2]]),
  #       dimnames = list(replicates = paste0("replicate", 1:replicates),
  #                       modelsite = 1:nrow(poccupy),
  #                       species = colnames(poccupy)))
  
  pdetect <- pdetectcond * occupied[fit$data$ModelSite, ]
  detected <- apply(pdetect, c(1, 2), function(x) rbinom(1, 1, x))
  
  return(detected)
}

simulate_LV <- function(fit, replaceinsitu = FALSE){
  simLV <- matrix(rnorm(fit$data$nlv * nrow(fit$data$Xocc)), ncol = fit$data$nlv)
  simLVbugsname <- matrix2bugsvar(simLV, name = "LV")
  if (!replaceinsitu) {return(simLVbugsname)}
  
  fit$mcmc <- lapply(fit$mcmc, function(arr) {
    arr[, names(simLVbugsname)] <- Rfast::rep_row(simLVbugsname, nrow(arr))
    return(arr)
  })
  return(fit)
}

#' @describeIn simulate_detections Simulate LV values and use these for simulating detections (with [simulate_detections()]). For each site the simulated LV values are copied across all draws.
simulate_detections_LV <- function(fit, esttype = "median"){
  fit <- simulate_LV(fit, replaceinsitu = TRUE)
  detected <- simulate_detections(fit, esttype = esttype)
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
#' @param ldet.b.min, ldet.b.max Same as occ.b.min and occ.b.max for the latent variable loadings.
#' @examples 
#' artfit <- artificial_runjags(nspecies = 2, nsites = 10, nvisitspersite = 4, modeltype = "jsodm_lv", nlv = 2)
#' \# with high correlation between occupancy of species
#' artfit <- artificial_runjags(nspecies = 2, nsites = 10, nvisitspersite = 4,
#'                               OccFmla = "~ 1",
#'                               occ.b.min = 0.8,
#'                               ldet.b.min = 0.3,
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
                               ldet.b.min = -0.5,
                               ldet.b.max = 0.5,
                               modeltype = "jsodm_lv",
                               ...
                               ){
  stopifnot(modeltype %in% availmodeltypes)
  
  # first step is to populate predictor values (including LV) and parameter values
  species <- make.names(rep(LETTERS, ceiling(nspecies / length(LETTERS))), unique = TRUE)[1:nspecies]
  covardfs <- simulate_covar_data(nsites, nvisitspersite)
  XoccProcess <- prep.designmatprocess(covardfs$Xocc, OccFmla)
  XobsProcess <- prep.designmatprocess(covardfs$Xobs, ObsFmla)

  
  data.list <- prepJAGSdata(modeltype,
                         covardfs$Xocc, yXobs = covardfs$Xobs,
                         ModelSite = "ModelSite", 
                         species = NULL,
                         XoccProcess = XoccProcess,
                         XobsProcess = XobsProcess,
                         ...) 
  fit <- list()
  fit$data <- data.list
  fit$data$n <- length(species)
  fit$species <- species
  fit$XoccProcess <- XoccProcess
  fit$XobsProcess <- XobsProcess
  fit$ModelSite <- "ModelSite"
  fit$summary.available <- TRUE
  fit$sample <- 1

  # set parameters
  occ.b <- matrix(runif( fit$data$n * fit$data$noccvar, min = occ.b.min, max = occ.b.max), nrow = fit$data$n, ncol = fit$data$noccvar, byrow = FALSE)
  det.b <- matrix(runif(  fit$data$n * fit$data$nobsvar, min = det.b.min, max = det.b.max), nrow = fit$data$n, ncol = fit$data$nobsvar, byrow = FALSE)
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
    ldet.b <- matrix(runif(  fit$data$n * fit$data$nlv, min = ldet.b.min, max = ldet.b.max), nrow = fit$data$n, ncol = fit$data$nlv) #0.5 constraint makes sure rowSum(ldet.b^2) < 1
    theta <- c(theta, 
               matrix2bugsvar(ldet.b, name = "ldet.b"),
               matrix2bugsvar(LV, name = "LV"))
  }
  fit$mcmc <- list()
  fit$mcmc[[1]] <- t(as.matrix(theta))
  
  class(fit) <- c(modeltype, class(fit))

  # simulate data using the LV values given above
  fit$data$y <- simulate_detections(fit, esttype = 1)
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
simulate_covar_data <- function(nsites, nvisitspersite){
  sites <- c(1:nsites)
  XoccIn <- data.frame(ModelSite = sites,
                     UpSite = sites,
                     Sine1 = 10 * sin(2 * pi * sites / nsites) + 100,
                     Sine2 = sin(4 * pi * sites / nsites))
  XobsIn <- data.frame(ModelSite = rep(sites, nvisitspersite),
                     UpVisit = 1:(nvisitspersite*nsites),
                     Step = c(rep(0, floor(nvisitspersite*nsites / 2)),
                              rep(1, ceiling(nvisitspersite*nsites / 2)))
                     )
  return(list(Xocc = XoccIn, Xobs = XobsIn))
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
