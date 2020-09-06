#' Compute Model Quality Information
#' @param fit a model fitted using run.detectionoccupancy
#' @param holdoutXocc A data frame of occupancy predictors for sites.
#' @param holdoutyXobs A data frame of detection information and species detections
#' @param ModelSite Column names used to match the holdoutyXobs to sites in holdoutXocc
#' @param cl A cluster created by parallel::makeCluster
#' @export
modelqualstats <- function(fit, holdoutXocc, holdoutyXobs, ModelSite, cl){
  lpd_holdout <- lppd.newdata(fit,
                       Xocc = holdoutXocc,
                       yXobs = holdoutyXobs,
                       ModelSite = ModelSite,
                       cl = cl)
  print("Computed: LPD of Holdout Sites")
  
  likel.mat <- likelihoods.fit(fit, chain = NULL,
                               cl = cl)
  chain_id <- lapply(1:length(fit$mcmc), function(x) rep(x, nrow(fit$mcmc[[x]])))
  chain_id <- as.integer(unlist(chain_id))
  waic <- loo::waic(log(likel.mat))
  r_eff <- loo::relative_eff(likel.mat, chain_id = chain_id, cores = length(cl))
  looest <- loo::loo(log(likel.mat), r_eff = r_eff, cores = length(cl))
  print("Computed: LOO-PSIS Estimate")
  
  prednumbers_holdout <- 
    predsumspecies_newdata(fit,
                           holdoutXocc,
                           holdoutyXobs,
                           ModelSiteVars = ModelSite,
                           cl = cl)
  print("Computed: Predicted Number of Species for Holdout Data")
  
  prednumbers_insample <- prednumbers_insample_fitLV <- prednumbers_insample_margLV <- NULL
  if (!is.null(fit$data$nlv) && fit$data$nlv > 0){
    prednumbers_insample_fitLV <- predsumspecies(fit, UseFittedLV = TRUE, cl = cl)
    prednumbers_insample_margLV <- predsumspecies(fit, UseFittedLV = FALSE, cl = cl)
  } else {
    prednumbers_insample <- predsumspecies(fit, UseFittedLV = FALSE, cl = cl)
  }
  print("Computed: Predicted Number of Species for Insample Data")
  
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