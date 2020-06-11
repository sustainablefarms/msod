## DS_residual_plotting
source("./R/DS_residuals.R")

# two situations:
# 1) ploting residuals against a covariate already fitted  in the model
# 2) plotting residauls against a covariate that is not yet included in the model

library(ggplot2); library(dplyr);
#' @param fit Is a runjags object created by fitting using package runjags.
#' @examples 
#' fit <- readRDS("./tmpdata/7_1_mcmcchain_20200424.rds")
#' fit <- runjags::add.summary(fit)
#' source("./R/calcpredictions.R")
#' source("./R/DS_residuals.R")
#' fit$data <- as.list.format(fit$data)
#' detection_resids <- ds_detection_residuals.fit(fit, type = "median", conditionalLV = FALSE)
#' 
#' # Plot Detection Residual 
#' plt <- plot_residuals_detection.fit(fit, detection_resids, varidx = 2, aggregatefcn = max, conditionalLV = FALSE)
#' plt
#' plt + coord_cartesian(ylim = c(-1, +1))
#' plot_residuals_detection.fit(fit, varidx = 2, esttype = "median", conditionalLV = FALSE)
#' plot_residuals_detection.fit(fit, varidx = c(2, 3), esttype = "median",
#'   plotfunction = facet_covariate, conditionalLV = FALSE) + 
#'   coord_cartesian(ylim = c(-1, +1))
#' 
#' plot_residuals_occupancy.fit(fit, varidx = 3, esttype = "median", conditionalLV = FALSE)
#' plot_residuals_occupancy.fit(fit, varidx = c(2, 3), esttype = "median",
#'  plotfunction = facet_covariate, conditionalLV = FALSE)
#' 
#' # Residuals against an unincluded covariate:
#' source("./scripts/7_1_import_site_observations.R")
#' covar <- occ_covariates[ , "ms_cover", drop = FALSE] %>%
#'   rowid_to_column(var = "ModelSite")
#' residuals <- ds_occupancy_residuals.fit(fit, type = "median", conditionalLV = FALSE)
#' pltobj <- plot_residuals.residual(residuals, covar, 
#'                  plotfunction = facet_species_covariate)
#' pltobj + scale_y_continuous(name = "Occupancy Residuals") 
#' pltobj + scale_y_continuous(name = "Occupancy Residuals") + coord_cartesian(ylim = c(-1, +1))
#' 
#' ## Example using 7_2_1 data
#' data_7_2_1 <- readRDS("./private/data/clean/7_2_1_input_data.rds")
#' fitp <- readRDS("./tmpdata/7_2_2_detectiononly_smallAmodel_run_20200529.rds")
#' fit <- results.jags("./runjagsfiles_13", read.monitor = c("LV", fitp$monitor))
#' fit <- add.summary(fit)
#' covar <- data_7_2_1$plotsmerged_detection %>%
#'   dplyr::select(-all_of(data_7_2_1$detection_data_specieslist)) %>%
#'   dplyr::select_if(is.numeric) %>%
#'   rename(ModelSite = ModelSiteID)
#' plot_residuals_detection.fit(fit, varidx = 2, esttype = "median", conditionalLV = TRUE) + coord_cartesian(ylim = c(-1, 1))
#' plt <- plot_LVvscovar.fit(fit, esttype = "median", covar = covar, aggregatefcn = min)
#' plt + coord_cartesian(ylim = c(-1, 1))

#' @describeIn plot_residuals For table of provided residuals and covariate values, makes residual plots.
#' Residual and covariate values must be provided with the ModelSite.
#' @param residuals is a table returned by \code{ds_occupancy_residuals.fit} or \code{ds_detection_residuals.fit}. Each row is a unique ModelSite.
#' It must have a column named 'ModelSite', all other columns are assumed contain residuals, and the column name corresponds to species names.
#' @param covar Is a table of covariate values, it must have a column labelled 'ModelSite' that gives the ModelSite of covariate value.
#' Rows with duplicated ModelSite values are averaged, which is in keeping with Warton et al for detection residuals.
#' @param plotfunction A plotting method to use. Default is \code{facet_species_covariate}.
#' @param ... Extra arguments to pass to plot function.
#' @value A ggplot object. Data is in the \code{data} slot.
plot_residuals.residual <- function(residuals, covar, plotfunction = facet_species_covariate,
                                    aggregatefcn = mean, ...){
  # Average covariates to ModelSite level. This is what Warton, Mackenzie et al do. In future it could be possible to present residuals per visit
  if (anyDuplicated(covar[, "ModelSite"]) > 0){
    warning("Multiple rows in 'covar' have the same ModelSite. These rows will aggregated using aggregatefcn")
    covar <- covar %>%
      as_tibble() %>%
      group_by(ModelSite) %>%
      summarise_all(aggregatefcn)
  }
  
  # prepare species names from input residuals
  speciesnames <- setdiff(colnames(residuals), "ModelSite")
  if (all(grepl("^S[0-9]+$", speciesnames))){  #order the species by index
    speciesnames <- speciesnames[order(as.integer(gsub("^S", "", speciesnames)))] 
  }
  
  # prepare covarnames
  covarnames <- setdiff(colnames(covar), "ModelSite")
  
  # prepare to print
  outframe <- left_join(residuals, covar, by = "ModelSite") %>%
    pivot_longer(any_of(speciesnames),
                 names_to = "Species",
                 values_to = "Residual",
                 values_drop_na = TRUE) %>%
    mutate(Species = ordered(Species, levels = speciesnames)) #to order species according to 'speciesnames'
  outframe <- outframe %>%
    pivot_longer(any_of(covarnames),
                 names_to = "Covariate",
                 values_to = "CovariateValue") 
  if (is.null(plotfunction)){plotfunction <- facet_species_covariate} #for when plotfunction = NULL is passed
  pltobject <- do.call(plotfunction, c(list(data = outframe), ...))
  return(pltobject)
}

#' @describeIn plot_residuals A function that prepares plot of residuals, one facet for each species and covariate
#' @param data Input tibble. Columns of Species, Residual, Covariate and CovariateValue
#' @param ... Extra arguments to pass. Currently ignored.
#' @value A ggplot object.
facet_species_covariate <- function(data, ...){
  pltobj <- data %>% 
    ggplot() +
    facet_wrap(~ Covariate + Species, as.table = TRUE, scale = "free_x") +
    geom_point(aes(x = CovariateValue, y = Residual)) +
    geom_hline(yintercept = 0, col = "blue", lty = "dashed") +
    geom_smooth(aes(x = CovariateValue, y = Residual), method = "gam", level = 0.95, formula = y ~ s(x, bs = "cs")) +
    scale_x_continuous(name = "Covariate Value")
  return(pltobj)
}

#' @describeIn plot_residuals A function that prepares plot of residuals, one facet for each covariate. Species ignored.
#' @param data Input tibble. Columns of Species, Residual, Covariate and CovariateValue
#' @param ... Extra arguments to pass. Currently ignored (no extra arguments accepted).
#' @value A ggplot object.
facet_covariate <- function(data, ...){
  data %>% 
    ggplot() +
    facet_wrap(~ Covariate, as.table = TRUE, scales = "free_x") +
    geom_point(aes(x = CovariateValue, y = Residual)) +
    geom_hline(yintercept = 0, col = "blue", lty = "dashed") +
    geom_smooth(aes(x = CovariateValue, y = Residual), method = "gam", level = 0.95, formula = y ~ s(x, bs = "cs")) +
    scale_x_continuous(name = "Covariate Value")
}

#' @describeIn plot_residuals Prepares tibbles and plots of residuals for covariates that are part of a fitted object
#' @param fit The fitted runjags object.
#' @param varidx The index of the covariate to plot. If NULL, then all covariates will be plotted
#' @param detectionresiduals Optional. A tibble of already calculated detection residuals.
#' Must have column names identical to output of \code{ds_detection_residuals.fit}.
#' If not supplied then detection residuals are computed from \code{fit} using \code{ds_detection_residuals.fit}.
#' @param esttype The point estimate extracted from fit.
#'  Passed to \code{ds_detection_residuals.fit} as argument \code{type}.
#'  See get_theta() for available options.
#' @value A ggplot object. Data is saved in the \code{data} slot.
plot_residuals_detection.fit <- function(fit, detectionresiduals = NULL, varidx = NULL, esttype = NULL, 
                                         conditionalLV = TRUE, aggregatefcn = mean, ...){
  stopifnot(is.null(detectionresiduals) | is.null(esttype))  #error if est type is supplied when detection residuals is also supplied
  fitdata <- as.list.format(fit$data)
  if (is.null(detectionresiduals)) {detectionresiduals <- ds_detection_residuals.fit(fit, type = esttype, conditionalLV = conditionalLV)}
  
  # get detection covariates
  detectioncovars <- fitdata$Xobs
  colnames(detectioncovars)[1:fitdata$Vobs] <- paste0("Xobs", 1:fitdata$Vobs)
  if (!is.null(varidx)){detectioncovars <- detectioncovars[, varidx, drop = FALSE]} 
  if ("ObservedSite" %in% names(fitdata)){ModelSite <- fitdata$ObservedSite} #to enable calculation on the early fitted objects with different name
  if ("ModelSite" %in% names(fitdata)){ModelSite <- fitdata$ModelSite}
  detectioncovars <- cbind(detectioncovars, ModelSite = ModelSite)
 
  extraargs = list(...)
  if ("plotfunction" %in% names(extraargs)){
    plotfunction <- extraargs[["plotfunction"]]
    extraargs["plotfunction"] <- NULL
  } else {plotfunction <- NULL}
  pltobject <- do.call(plot_residuals.residual, c(list(
    residual = detectionresiduals,
    covar = detectioncovars,
    plotfunction = plotfunction,
    aggregatefcn = aggregatefcn),
    extraargs)
    ) + 
      scale_y_continuous(name = "Detection Residual") 
  return(pltobject)
}

#' @describeIn ?? Prepares tibbles and plots of occupancy residuals for covariates that are part of a fitted object
#' @param fit The fitted runjags object.
#' @param varidx The index of the covariate to plot. If NULL, then all occupancy covariates will be plotted
#' @param occupancyresiduals Optional. A tibble of already calculated occupancy residuals.
#' Must have column names identical to output of \code{ds_detection_residuals.fit}.
#' If not supplied then detection residuals are computed from \code{fit} using \code{ds_detection_residuals.fit}.
#' @param esttype The point estimate extracted from fit. Passed to \code{ds_detection_residuals.fit} as argument \code{type}.
#' @value A ggplot object. Data is saved in the \code{data} slot.
plot_residuals_occupancy.fit <- function(fit, occupancyresidual = NULL, varidx = NULL,
                                         esttype = NULL, conditionalLV = TRUE, aggregatefcn = mean, ...){
  stopifnot(is.null(occupancyresidual) | is.null(esttype))  #error if est type is supplied when detection residuals is also supplied
  fitdata <- as.list.format(fit$data)
  if (is.null(occupancyresidual)) {occupancyresidual <- ds_occupancy_residuals.fit(fit, type = esttype, conditionalLV = conditionalLV)}
  
  # get occupancy covariates
  occupancycovars <- fitdata$Xocc
  colnames(occupancycovars)[1:fitdata$Vocc] <- paste0("Xocc", 1:fitdata$Vocc)
  if (!is.null(varidx)){occupancycovars <- occupancycovars[, varidx, drop = FALSE]}  
  occupancycovars <- occupancycovars %>% as_tibble() %>% rowid_to_column(var = "ModelSite")
  
  extraargs = list(...)
  if ("plotfunction" %in% names(extraargs)){
    plotfunction <- extraargs[["plotfunction"]]
    extraargs["plotfunction"] <- NULL
  } else {plotfunction <- NULL}
  pltobject <- do.call(plot_residuals.residual, c(list(
    residual = occupancyresidual,
    covar = occupancycovars,
    plot = FALSE,
    plotfunction = plotfunction,
    aggregatefcn = aggregatefcn),
    extraargs)
    ) + 
    scale_y_continuous(name = "Occupancy Residual")
  return(pltobject)
}


#' @describeIn ?? Plot estimated latent variable values against covariate values for occupancy
#' @param theta a vector of model parameters with BUGS names
#' @param fit A fitted runjags object.
#' @param esttype When fit is non-NULL, then esttype is used to extract parameter values
#' @param covar a dataframe of covariate values. It must have a column labelled 'ModelSite' that gives the ModelSite of covariate value.
#' @param aggregatefcn Rows in covar with duplicate ModelSite values are aggregrated using this function
plot_LVvscovar.fit <- function(fit, esttype = "median", theta = NULL, covar, facetvars = NULL, cuts = 3, aggregatefcn = mean){
  df <- plotdf_LVvscovar.fit(fit, esttype = esttype, theta = theta, covar = covar, aggregatefcn = aggregatefcn)
 
  if (!is.null(facetvars)){
    df <- df %>%
      mutate_at(facetvars, ~cut(., cuts))
  } 
  dflong <- df %>%
    pivot_longer(starts_with("LV"), names_to = "LV Name", values_to = "LV Value") %>%
    pivot_longer(setdiff(names(covar), c("ModelSite", facetvars)), names_to = "Covariate Name", values_to = "Covariate Value")
  if (is.null(facetvars)){
    plt <- dflong %>%
      ggplot() +
      geom_point(aes_(y = ~`LV Value`, x = ~`Covariate Value`)) +
      facet_wrap(vars(`LV Name`, `Covariate Name`), scales = "free") +
      geom_hline(yintercept = 0, col = "blue", lty = "dashed") +
      geom_smooth(aes_(y = ~`LV Value`, x = ~`Covariate Value`), method = "gam", level = 0.95, formula = y ~ s(x, bs = "cs")) 
  } else {
    plt <- dflong %>%
      ggplot() +
      geom_point(aes_(y = ~`LV Value`, x = ~`Covariate Value`, col = as.name(facetvars))) +
      facet_wrap(vars(`LV Name`, `Covariate Name`), scales = "free") +
      geom_hline(yintercept = 0, col = "blue", lty = "dashed") +
      geom_smooth(aes_(y = ~`LV Value`, x = ~`Covariate Value`, col = as.name(facetvars)), method = "gam", level = 0.95, formula = y ~ s(x, bs = "cs")) 
  }
  return(plt)
}

#' @describeIn ?? Prepare plotting data frames for latent variable values against covariate values for occupancy
#' @param theta a vector of model parameters with BUGS names
#' @param fit A fitted runjags object.
#' @param esttype When fit is non-NULL, then esttype is used to extract parameter values
#' @param covar a dataframe of covariate values. It must have a column labelled 'ModelSite' that gives the ModelSite of covariate value.
#' @param aggregatefcn Rows in covar with duplicate ModelSite values are aggregrated using this function
plotdf_LVvscovar.fit <- function(fit, esttype = "median", theta = NULL, covar, aggregatefcn = mean){
  if (is.null(theta)){theta <- get_theta(fit, type = esttype)}
  if (anyDuplicated(covar[, "ModelSite"]) > 0){
    warning("Multiple rows in 'covar' have the same ModelSite. These rows will be aggregated using aggregatefcn")
    covar <- covar %>%
      as_tibble() %>%
      group_by(ModelSite) %>%
      summarise_all(aggregatefcn)
  }
  
  fitdata <- as.list.format(fit$data)
  ## LV values
  LV <- bugsvar2matrix(theta, "LV", 1:fitdata$J, 1:fitdata$nlv) # rows are model sites, columns are latent variables
  colnames(LV) <- paste0("LV", 1:ncol(LV))
  LV <- cbind(ModelSite = 1:nrow(LV), LV) %>% as_tibble() 
  df <- dplyr::left_join(LV, covar, by = "ModelSite", suffix = c(".LV", ".X"))
  return(df)
}