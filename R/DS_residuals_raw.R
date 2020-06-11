#' @title Dunn-Smyth Residuals for Detection and Occupancy for Given Precalculated Predictions of Probabilities
#' @references D. I. Warton, J. Stoklosa, G. Guillera-Arroita, D. I. MacKenzie, and A. H. Welsh, 
#' "Graphical diagnostics for occupancy models with imperfect detection," Methods in Ecology and Evolution, vol. 8, no. 4, pp. 408-419, 2017, doi: 10.1111/2041-210X.12761.
# source("./R/calcpredictions.R")

#' @param fit Is a runjags object created by fitting using package runjags.

#' @describeIn DunnSmythResidualsRaw For supplied (discrete) CDF evaluations,  compute Dunn-Syth residuals.
#' @param cdf The CDF of a discrete random variable evaluated at observations of iid copies of this random variable. Must be a list or vector.
#' @param cdfminus For observations of iid copies of a discrete random variable,
#'  the CDF of the random variable evaluated at the value just below the observed values. Must be a list or vector of same length as \code{cdf}.
cdfvals_2_dsres_discrete <- function(cdf, cdfminus, seed = NULL){
  set.seed(seed)
  u<-runif(length(cdf));  # Standard uniform value to "jitter" the cdf.
  residDet<-qnorm(cdf*u + cdfminus*(1-u));      # Dunn-Smyth residual, standard normal if cdf correct.
  return(residDet)
}

#' @describeIn DunnSmythResidualsRaw For a given ModelSite and species,
#'  computes the (discrete) probability density function (pdf) for the number of detections of the species at the ModelSite.
#' Based heavily on Rpresence code for hetpdf() initially (but not the rest this file) WARNING: do not know what the license for Rpresense source code is!! Can't use code directly at least.
#' @param pDetected Is a list of the detection probabilities for each visit to the ModelSite for a single species
#' @param x Is the value at which to evaluate the pdf
numdet_pdf<-function(x, pDetected){
  p0 <- prod(1 - pDetected) #probability of no detections
  if (x == 0){return(p0)}
  if ((x < 0) | (x > length(pDetected))) {return(0)} #x outside support
  if ((abs(as.integer(x) - x)) > 1E-8){return(0)} #x not an integer
  ind<-combn(length(pDetected), x) #index of visits that the species could be detected at
  
  tmp_pdens <- matrix(pDetected[as.vector(ind)], byrow = FALSE, nrow = x)
  pdends_ind <- apply(tmp_pdens / (1 - tmp_pdens), 2, prod) *  #the detections in numerator and denominator because p0 used for non-detections
    p0

  return(sum(pdends_ind));
}

#' @describeIn DunnSmythResidualsRaw For a given ModelSite and species,
#'  computes the (discrete) cumulative distribution function (cdf) for the number of detections of the species at the ModelSite.
#' Based heavily on Rpresence code for numdet_cdf initially (but not the rest this file) WARNING: do not know what the license for Rpresense source code is!! Can't use code directly at least.
#' @param pDetected Is a list of the detection probabilities for each visit to the ModelSite for a single species
#' @param x Is the value at which to evaluate the pdf
numdet_cdf <- function(x, pDetected){
  if (x < 0) {return(0)}
  x <- floor(x)
  x <- min(x, length(pDetected))
  pvals <- vapply(0:x, numdet_pdf, pDetected = pDetected, FUN.VALUE = 0.33)
  return(sum(pvals))
}
 
#' @describeIn DunnSmythResidualsRaw Given predictions for detection probability conditioned on the site being occupied,
#'  and corresponding and detection observations, compute Dunn-Smyth residuals for detection
#' @param preds is a dataframe with columns Species, ModelSite, and pDetected
#' @param obs is a dataframe with columns Species, ModelSite, and Detected
#' @return A dataframe with a columns for Species, ModelSite, and detection residual. 
#' The residual is only computed for species detected at least once at a site.
#' @export
ds_detection_residuals.raw <- function(preds, obs, seed = NULL){
  stopifnot(all(c("Species", "ModelSite", "pDetected") %in% names(preds)))
  stopifnot(all(c("Species", "ModelSite", "Detected") %in% names(obs)))
  stopifnot(isTRUE(all.equal(preds[, c("Species", "ModelSite")], obs[, c("Species", "ModelSite")])))
  combined <- cbind(preds, Detected = obs$Detected)
  persite <- combined %>%
    dplyr::group_by(Species, ModelSite) %>%
    dplyr::summarise(numdet = sum(Detected),
                     pDetected = list(pDetected)) %>%
    dplyr::filter(numdet > 0)  # detection residuals only use sites where a detection occured.
  # can't use mutate because numdet_cdf is not vectorised for the x and pDetected argument yet
  stopifnot(nrow(persite) > 0) #means there is no detection residuals to compute

  #the following are not yet conditional on detection greater than 0
  cdfminus <- unlist(mapply(numdet_cdf, x = persite$numdet - 1, pDetected = persite$pDetected, SIMPLIFY = FALSE))
  pdfx <- unlist(mapply(numdet_pdf, x = persite$numdet, pDetected = persite$pDetected, SIMPLIFY = FALSE))

  # condition on non-zero detection
  cdf0 <- unlist(mapply(numdet_cdf, x = 0, pDetected = persite$pDetected, SIMPLIFY = FALSE))
  cdfminus_cond <- condition_nonzero.cdf(cdf0, cdfminus)
  pdfx_cond <- condition_nonzero.pdf(cdf0, pdfx)
  
  ds_resids <- cdfvals_2_dsres_discrete(pdfx_cond + cdfminus_cond, cdfminus_cond, seed = seed)
  return(data.frame(Species = persite$Species, ModelSite = persite$ModelSite, DetectionResidual = ds_resids))
}

# @param cdfval must be cdf for 0 for higher
condition_nonzero.cdf <- function(cdf0, cdfval){
  return((cdfval - cdf0)/(1 - cdf0))
}

# @param pdfval must be pdf for values > 0
condition_nonzero.pdf <- function(cdf0, pdfval){
  return(pdfval / (1 - cdf0))
}

#' @describeIn DunnSmythResidualsRaw Given occupancy predictions and detection observations, compute Dunn-Smyth residuals for occupancy
#' @param preds is a dataframe with columns Species, ModelSite, pOccupancy, and pDetected_cond.
#'  pOccupancy is the probability of ModelSite being occupied.
#'  pDetected_cond is the probability of detecting the species, given the species occupies the ModelSite.
#' @param obs is a dataframe with columns Species, ModelSite, and Detected
#' @return A dataframe with a columns for Species, ModelSite, and occupancy residual. 
#' @export
ds_occupancy_residuals.raw <- function(preds, obs, seed = NULL){
  stopifnot(all(c("Species", "ModelSite", "pOccupancy", "pDetected_cond") %in% names(preds)))
  stopifnot(all(c("Species", "ModelSite", "Detected") %in% names(obs)))
  stopifnot(isTRUE(all.equal(preds[, c("Species", "ModelSite")], obs[, c("Species", "ModelSite")])))
  combined <- cbind(preds, Detected = obs$Detected)
  persite <- combined %>%
    dplyr::group_by(Species, ModelSite) %>%
    dplyr::summarise(anyDetected = sum(Detected) > 0,
                     pDetected_cond = list(pDetected_cond),
                     pOccupancy = first(pOccupancy))  #using first here as a shortcut --> should be all identical
                     # pOccUnique = (1 == length(unique((pOccupancy)))))  # a check that pOccupancy unique, removed for speed
  # stopifnot(all(persite$pOccUnique)) #check that pOccupancy values are unique to site x species
  
  # probability of no detection
  pNoDetect <- unlist(mapply(function(pOcc, pDetected_cond) (1 - pOcc) + pOcc * prod(1 - pDetected_cond), 
         pOcc = persite$pOccupancy,
         pDetected_cond = persite$pDetected_cond,
         SIMPLIFY = FALSE))
  
  # cdf value
  persite$cdf <- 1
  persite$cdf[persite$anyDetected == FALSE] <- pNoDetect[persite$anyDetected == FALSE]
  
  # cdf minus value
  persite$cdfminus <- 0
  persite$cdfminus[persite$anyDetected == TRUE] <- pNoDetect[persite$anyDetected == TRUE]
  
  ds_resids <- cdfvals_2_dsres_discrete(persite$cdf, persite$cdfminus, seed = seed)
  return(data.frame(Species = persite$Species, ModelSite = persite$ModelSite, OccupancyResidual = ds_resids))
}
