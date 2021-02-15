#' Compute Model Quality Information
#' @param fit a model fitted using run.detectionoccupancy
#' @param holdoutXocc A data frame of occupancy predictors for sites.
#' @param holdoutyXobs A data frame of detection information and species detections
#' @param ModelSite Column names used to match the holdoutyXobs to sites in holdoutXocc
#' @param cl A cluster created by parallel::makeCluster
#' @examples 
#' fit <- readRDS("../sflddata/private/data/testdata/cutfit_7_4_11_2LV.rds")
#' fit <- translatefit(fit)
#' Xocc <- unstandardise.designmatprocess(fit$XoccProcess, fit$data$Xocc)
#' Xocc$ModelSite <- 1:nrow(Xocc)
#' Xobs <- unstandardise.designmatprocess(fit$XobsProcess, fit$data$Xobs)
#' Xobs$ModelSite <- fit$data$ModelSite
#' y <- fit$data$y
#' yXobs <- cbind(Xobs, y)
#' cl <- parallel::makeCluster(1)
#' quality <- modelqualstats(fit, holdoutXocc = Xocc, holdoutyXobs = yXobs, ModelSite = 'ModelSite', cl = cl)
#' @export
modelqualstats <- function(fit, holdoutXocc, holdoutXobs, holdoutModelSite, holdouty, cl, nlvperdraw = 1000, ...){
  fit_wholdoutdata <- supplant_new_data(fit, Xocc = holdoutXocc, Xobs = holdoutXobs, ModelSite = holdoutModelSite, y = holdouty, ...)
  
  holdout_quality <- modelqualstats_holdout(fit_wholdoutdata, cl = cl, nlvperdraw = nlvperdraw)
  insample_quality <- modelqualstats_insample(fit, cl = cl, nlvperdraw = nlvperdraw)
  
  quality <- list(
    holdout = holdout_quality,
    insample = insample_quality)
  return(quality)
}

modelqualstats_holdout <- function(fit_wholdoutdata, cl, nlvperdraw = 1000){
  likel.mat_holdout <- likelihood(fit_wholdoutdata, numlvsims = 1)
  lpd_holdout <- elpd(likel.mat_holdout)
  print("Estimated LPD of Holdout Sites")
  
  prednumbers_holdout <- NULL
  if ("jsodm_lv" %in% class(fit_wholdoutdata)){
    prednumbers_holdout <- speciesrichness(fit_wholdoutdata, occORdetection = "detection",
                                           usefittedlvv = FALSE,
                                           nlvperdraw = nlvperdraw)
    print("Computed: Predicted Number of Species for Holdout Data")
  } else if ("jsodm" %in% class(fit_wholdoutdata)){
    prednumbers_holdout <- speciesrichness(fit_wholdoutdata, occORdetection = "detection")
    print("Computed: Predicted Number of Species for Holdout Data")
  }
  return(list(
    lpd = lpd_holdout,
    predspecnum = prednumbers_holdout)
  )
}

modelqualstats_insample <- function(fit, cl, nlvperdraw = 1000){
  likel.mat <- likelihood(fit, cl = cl, numlvsims = nlvperdraw)
  chain_id <- lapply(1:length(fit$mcmc), function(x) rep(x, nrow(fit$mcmc[[x]])))
  chain_id <- as.integer(unlist(chain_id))
  waic <- loo::waic(log(likel.mat))
  r_eff <- loo::relative_eff(likel.mat, chain_id = chain_id, cores = length(cl))
  looest <- loo::loo(log(likel.mat), r_eff = r_eff, cores = length(cl))
  print("Computed: LOO-PSIS Estimate")
  
  
  prednumbers_insample <- prednumbers_insample_fitLV <- prednumbers_insample_margLV <- NULL
  if ("jsodm_lv" %in% class(fit)){
    prednumbers_insample_fitLV <- speciesrichness(fit, occORdetection = "detection", usefittedlvv = TRUE)
    prednumbers_insample_margLV <- speciesrichness(fit, occORdetection = "detection", usefittedlvv = FALSE, nlvperdraw = nlvperdraw)
    print("Computed: Predicted Number of Species for Insample Data")
  } else if ("jsodm" %in% class(fit)){
    prednumbers_insample <- speciesrichness(fit, occORdetection = "detection")
    print("Computed: Predicted Number of Species for Insample Data")
  }
  
  quality <- list(
      waic = waic,
      loo = looest,
      predspecnum = prednumbers_insample,
      predspecnum_fitLV = prednumbers_insample_fitLV,
      predspecnum_margLV = prednumbers_insample_margLV)
  return(quality)
}