#' @title Dunn-Smyth Residuals for Detection and Occupancy for Fitted Models
#' @references D. I. Warton, J. Stoklosa, G. Guillera-Arroita, D. I. MacKenzie, and A. H. Welsh, 
#' "Graphical diagnostics for occupancy models with imperfect detection," Methods in Ecology and Evolution, vol. 8, no. 4, pp. 408-419, 2017, doi: 10.1111/2041-210X.12761.
# source("./R/calcpredictions.R")

#' @param fit Is a runjags object created by fitting using package runjags.
#' @examples 
#' fit <- readRDS("./tmpdata/7_2_10_clim_year.rds")
#' detection_resids <- ds_detection_residuals.fit(fit, type = "median")
#' occupancy_resids <- ds_occupancy_residuals.fit(fit, type = "median", conditionalLV = FALSE)
#' fitold <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_allhyp_vif_logwoody500m_msnm_year_Time_Wind.rds")
#' fit <- translatefit(fitold)
#' ds_detection_residuals.fit(fit)
#' ds_occupancy_residuals.fit(fit)

##### Full Dunn-Smyth Residual Functions ####
#' @describeIn DunnSmythResiduals Given a fitted occupancy detection model (variable names must match). Computes Dunn-Smyth residuals for detection, marginalising the latent variables.
#' @param fit A fitted occupancy-detection model.
#' @param seed A seed to fix randomness of Dunn-Smyth residual jitter.
#' @param type The type of point estimate to use for parameter estimates. See \code{\link{get_theta}}
#' @return A matrix, each row is a ModelSite and each column is a species.
#' Detection residuals are only computed for species and sites that have at least one detection. Other values are NA.
#' @export
ds_detection_residuals.fit <- function(fit, type = "median", seed = NULL){
  fit$data <- as_list_format(fit$data)
  det.v <- fit$data$Xobs
  det.b <- get_det_b(fit, usesummary = type)
  pDetection <- pdet_occ_raw.jsodm(det.v, det.b)  #the detection probabilities, assuming occupied
  detections <-  fit$data$y
  if (is.null(colnames(pDetection))){colnames(pDetection) <- paste0("S", 1:ncol(pDetection))} #name the species S1....Sn
  if (is.null(colnames(detections))) {#name the columns if possible
    if (!is.null(fit$species)) {colnames(detections) <- fit$species}
    else {colnames(detections) <- paste0("S", 1:ncol(detections))}
  }
  if ("ObservedSite" %in% names(fit$data)){ModelSite <- fit$data$ObservedSite} #to enable calculation on the early fitted objects with different name
  if ("ModelSite" %in% names(fit$data)){ModelSite <- fit$data$ModelSite}

  # Convert the above into format suitable for ds_detection_residuals.raw
  preds <- cbind(ModelSite = as.numeric(ModelSite), VisitId = 1:nrow(fit$data$Xobs), drop_to_matrix(pDetection)) %>%
    as_tibble() %>%
    tidyr::pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "pDetected") %>%
    arrange(VisitId, Species, ModelSite)
  obs <- cbind(ModelSite = as.numeric(ModelSite), VisitId = 1:nrow(fit$data$Xobs), detections) %>%
    as_tibble() %>%
    tidyr::pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "Detected") %>%
    arrange(VisitId, Species, ModelSite)
  
  # Compute residuals
  detection_resids <- ds_detection_residuals.raw(preds, obs, seed = seed)
  detection_resids %>%
    tidyr::pivot_wider(names_from = "Species",
                values_from = "DetectionResidual") %>%
    return()
}

#' @describeIn DunnSmythResiduals Given a fitted occupancy detection model (variable names must match).
#'  Computes Dunn-Smyth residuals for occupancy, marginalising the latent variables.
#' @param fit A fitted occupancy-detection model.
#' @param seed A seed to fix randomness of Dunn-Smyth residual jitter.
#' @param type The type of point estimate to use for parameter estimates. See \code{\link{get_theta}}
#' @return A matrix, each row is a ModelSite and each column is a species.
#' @export
ds_occupancy_residuals.fit <- function(fit, type = "median", seed = NULL, conditionalLV = TRUE){
  fit$data <- as_list_format(fit$data)
  pOccupancy <- poccupy(fit, usethetasummary = type)
  if (is.null(colnames(pOccupancy))){colnames(pOccupancy) <- paste0("S", 1:ncol(pOccupancy))} #name the species S1....Sn
  pDetected_cond <- pdet_occ(fit, usethetasummary = type)  #the detection probabilities if sites occupied
  if (is.null(colnames(pDetected_cond))){colnames(pDetected_cond) <- paste0("S", 1:ncol(pDetected_cond))} #name the species S1....Sn
  detections <- fit$data$y
  if (is.null(colnames(detections))) {#name the columns if possible
    if (!is.null(fit$species)) {colnames(detections) <- fit$species}
    else {colnames(detections) <- paste0("S", 1:ncol(detections))}
  }
  if ("ObservedSite" %in% names(fit$data)){ModelSite <- fit$data$ObservedSite} #to enable calculation on the early fitted objects with different name
  if ("ModelSite" %in% names(fit$data)){ModelSite <- fit$data$ModelSite}
  
  # convert to format for raw function
  pOccupancy <- cbind(ModelSite = 1:nrow(fit$data$Xocc), drop_to_matrix(pOccupancy)) %>%
    as_tibble() %>%
    tidyr::pivot_longer(-ModelSite, names_to = "Species", values_to = "pOccupancy")
  pDetCondOcc <- cbind(ModelSite = as.numeric(ModelSite), VisitId = 1:nrow(fit$data$Xobs), drop_to_matrix(pDetected_cond)) %>%
    as_tibble() %>%
    tidyr::pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "pDetected_cond")
  preds <- inner_join(pOccupancy, pDetCondOcc, by = c("ModelSite", "Species")) %>% arrange(VisitId, Species, ModelSite)
  obs <- cbind(ModelSite = as.numeric(ModelSite), VisitId = 1:nrow(fit$data$Xobs), detections) %>%
    as_tibble() %>%
    tidyr::pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "Detected") %>%
    arrange(VisitId, Species, ModelSite)
  
  # apply raw occupancy residuals function
  residuals <- ds_occupancy_residuals.raw(preds, obs, seed = seed)
  residuals %>%
    tidyr::pivot_wider(names_from = "Species",
                values_from = "OccupancyResidual") %>%
    return()
}


#### Extra Functions ####
# simulate number of detection (replicated n times), for visit with probability of detection given by pDetected
simDetectedDistr <- function(n, pDetected){
  numdetected <- replicate(n,
            sum(vapply(pDetected, function(x) rbinom(n = 1, size = 1, prob = x), FUN.VALUE = 1)))
  return(numdetected)
}

