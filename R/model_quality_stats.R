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
modelqualstats <- function(fit, holdoutXocc, holdoutyXobs, ModelSite, cl, nlvperdraw = 1000){
  likel.mat_holdout <- apply_to_new_data(likelihood, fit, 
                                   Xocc = holdoutXocc, 
                                   Xobs = holdoutyXobs[, !(colnames(holdoutyXobs) %in% fit$species), drop = FALSE],
                                   ModelSite = ModelSite,
                                   y = holdoutyXobs[, colnames(holdoutyXobs) %in% fit$species, drop = FALSE],
                                   funargs = list(cl = cl, numlvsims = nlvperdraw))
  lpd_holdout <- elpd(likel.mat_holdout)
  print("Estimated LPD of Holdout Sites")
  
  likel.mat <- likelihood(fit, cl = cl, numlvsims = nlvperdraw)
  chain_id <- lapply(1:length(fit$mcmc), function(x) rep(x, nrow(fit$mcmc[[x]])))
  chain_id <- as.integer(unlist(chain_id))
  waic <- loo::waic(log(likel.mat))
  r_eff <- loo::relative_eff(likel.mat, chain_id = chain_id, cores = length(cl))
  looest <- loo::loo(log(likel.mat), r_eff = r_eff, cores = length(cl))
  print("Computed: LOO-PSIS Estimate")
  

  prednumbers_holdout <- prednumbers_insample <- prednumbers_insample_fitLV <- prednumbers_insample_margLV <- NULL
  if ("jsodm_lv" %in% class(fit)){
    prednumbers_holdout <- 
      apply_to_new_data(speciesrichness, fit, 
                        holdoutXocc, holdoutyXobs[, !(colnames(holdoutyXobs) %in% fit$species), drop = FALSE],
                        ModelSite = ModelSite,
                        funargs = list(occORdetection = "detection",
                                       usefittedlvv = FALSE,
                                       nlvperdraw = nlvperdraw))
    print("Computed: Predicted Number of Species for Holdout Data")
    
    prednumbers_insample_fitLV <- speciesrichness(fit, occORdetection = "detection", usefittedlvv = TRUE)
    prednumbers_insample_margLV <- speciesrichness(fit, occORdetection = "detection", usefittedlvv = FALSE, nlvperdraw = nlvperdraw)
    print("Computed: Predicted Number of Species for Insample Data")
  } else if ("jsodm" %in% class(fit)){
    prednumbers_holdout <- 
      apply_to_new_data(speciesrichness,fit, 
                        holdoutXocc, holdoutyXobs[, !(colnames(holdoutyXobs) %in% fit$species), drop = FALSE],
                        ModelSite = ModelSite,
                        funargs = list(occORdetection = "detection"))
    print("Computed: Predicted Number of Species for Holdout Data")
    prednumbers_insample <- speciesrichness(fit, occORdetection = "detection")
    print("Computed: Predicted Number of Species for Insample Data")
  }
  
  quality <- list(
    holdout = list(
      lpd = lpd_holdout,
      predspecnum = prednumbers_holdout),
    insample = list(
      waic = waic,
      loo = looest,
      predspecnum = prednumbers_insample,
      predspecnum_fitLV = prednumbers_insample_fitLV,
      predspecnum_margLV = prednumbers_insample_margLV)
  )
  return(quality)
}