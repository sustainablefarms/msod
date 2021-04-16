#' @title Supplant data in fitted object with new data for the purposes of predicting on new data
#' @param fit Fitted object
#' @param Xocc A dataframe of covariates related to occupancy. One row per ModelSite.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param Xobs A dataframe of covariates related to observations. One row per visit.
#' Each column is a covariate.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param y A dataframe of species observations (1 or 0). Each column is a species. Rows must correspond to Xobs
#' @param ModelSite A vector of integers of equal length to the number of visits. Each entry species for the corresponding visit, the of Xocc that was visited 
#' @param toXocc NULL to use a process saved by [sflddata::save_process], otherwise a function that takes Xocc as the only input and returns a model matrix for the occupancy component.
#' @param toXobs like `toXocc`, except for Xobs and the detection component of the model.
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
#' ModelSite <- as.integer(Xobs$ModelSite)
#' fitwnewdata <- supplant_new_data(fit, Xocc, Xobs, ModelSite = ModelSite, y = y)
#' ds_detection_residuals.fit(fitwnewdata) #detection residuals on this new data
#' ds_occupancy_residuals.fit(fitwnewdata) #detection residuals on this new data
#'
#' fitwnewdata <- supplant_new_data(fit, Xocc)
#' 
#' @export
supplant_new_data <- function(fit, Xocc, Xobs = NULL, ModelSite = NULL, y = NULL, toXocc = NULL, toXobs = NULL, ...){
  UseMethod("supplant_new_data")
}

#' @export
supplant_new_data.jsodm <- function(fit, Xocc, Xobs = NULL, ModelSite = NULL, y = NULL, toXocc = NULL, toXobs = NULL){
  if (is.null(toXocc)) { toXocc <- fit$toXocc }
  if (is.null(toXobs)) { toXobs <- fit$toXobs }
  if (is.null(toXocc)){
    warning("Using saved XoccProcess. This functionality will become obsolete.")
    toXoccFun <- function(indf, mainparams = fit$XoccProcess){
      Xocc <- sflddata::apply.designmatprocess(mainparams, indf)
      return(Xocc)
    }
    toXocc <- sflddata::save_process(toXoccFun, checkwith = Xocc, params = list(mainparams = fit$XoccProcess))
  }
  if (is.null(toXobs) && !is.null(Xobs)){
    warning("Using saved XobsProcess. This functionality will become obsolete.")
    toXobsFun <- function(indf, mainparams = fit$XobsProcess){
      Xobs <- sflddata::apply.designmatprocess(mainparams, indf)
      return(Xobs)
    }
    toXobs <- sflddata::save_process(toXoccFun, checkwith = Xobs, params = list(mainparams = fit$XobsProcess))
  }
  
  if (is.function(toXocc)){
    Xocc <- toXocc(Xocc)
  } else {
    Xocc <- sflddata::apply_saved_process(toXocc, Xocc)
  }
  if (!is.null(Xobs)){
    if (is.function(toXobs)){
      Xobs <- toXobs(Xobs)
    } else {
      Xobs <- sflddata::apply_saved_process(toXobs, Xobs)
    }
  }
  
  data.list <- prepJAGSdata2("jsodm",
                Xocc = Xocc,
                Xobs = Xobs,
                y = y,
                ModelSite = ModelSite)
  fit$data$Xocc <- data.list$Xocc
  fit$data$Xobs <- data.list$Xobs
  fit$data$y <- data.list$y
  fit$data$ModelSite <- data.list$ModelSite
  fit$data$nmodelsites <- data.list$nmodelsites
  fit$data$nvisits <- data.list$nvisits
  return(fit)
}

#' @export
supplant_new_data.jsodm_lv <- function(fit, Xocc, Xobs = NULL, ModelSite = NULL, y = NULL, toXocc = NULL, toXobs = NULL){
  fit <- supplant_new_data.jsodm(fit, Xocc, Xobs = Xobs, ModelSite = ModelSite, y = y, toXocc = toXocc, toXobs = toXobs)
  #remove fitted lvv as no longer valid to new sites
  bugsnames_lvv <- grepl("^lv.v", colnames(fit$mcmc[[1]]))
  fit$mcmc <- lapply(fit$mcmc,
         function(x){
           out <- x[, !bugsnames_lvv]
           return(out)
         })
  return(fit)
}
