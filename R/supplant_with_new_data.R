#' @title Supplant data in fitted object with new data for the purposes of predicting on new data
#' @param fit Fitted object
#' @param Xocc A dataframe of covariates related to occupancy. One row per ModelSite.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param Xobs A dataframe of covariates related to observations. One row per visit.
#' Each column is a covariate.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param y A dataframe of species observations (1 or 0). Each column is a species. Rows must correspond to Xobs
#' @param ModelSite A list of column names in Xocc and yXobs that uniquely specify the ModelSite. Can be simply a ModelSite index
#' @examples 
#' fitold <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_model_2lv_e13.rds")
#' fit <- translatefit(fitold)
#' originalXocc <- unstandardise.designmatprocess(fit$XoccProcess, fit$data$Xocc)
#' originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
#' originalXobs <- unstandardise.designmatprocess(fit$XobsProcess, fit$data$Xobs)
#' originalXobs <- cbind(ModelSite = fit$data$ModelSite, originalXobs)
#' Xocc <- originalXocc[1:10, ]
#' Xobs <- originalXobs[originalXobs$ModelSite %in% Xocc$ModelSite, ]
#' y <- fit$data$y[originalXobs$ModelSite %in% Xocc$ModelSite, ]
#' fitwnewdata <- supplant_new_data.jsodm_lv(fit, Xocc, Xobs, ModelSite = "ModelSite", y = y)
#' ds_detection_residuals.fit(fitwnewdata) #detection residuals on this new data
#' ds_occupancy_residuals.fit(fitwnewdata) #detection residuals on this new data
#' 
#' @export
supplant_new_data <- function(fit, Xocc, Xobs = NULL, ModelSite = NULL, y = NULL, ...){
  UseMethod("supplant_new_data")
}

#' @export
supplant_new_data.jsodm <- function(fit, Xocc, Xobs = NULL, ModelSite = NULL, y = NULL){
  if (is.null(Xobs)){
    yXobs <- NULL
  } else {
    if (is.null(y)) {yXobs <- Xobs}
    else {yXobs <- cbind(Xobs, y)}
  }
  sitedata <- prep_new_data(fit, Xocc, yXobs, ModelSite)
  fit$data <- sitedata
  return(fit)
}

#' @export
supplant_new_data.jsodm_lv <- function(fit, Xocc, Xobs = NULL, ModelSite = NULL, y = NULL){
  fit <- supplant_new_data.jsodm(fit, Xocc, Xobs, ModelSite, y = y)
  #remove fitted lvv as no longer valid to new sites
  bugsnames_lvv <- grepl("^lv.v", colnames(fit$mcmc[[1]]))
  fit$mcmc <- lapply(fit$mcmc,
         function(x){
           out <- x[, !bugsnames_lvv]
           return(out)
         })
  return(fit)
}
