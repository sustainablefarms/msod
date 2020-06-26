#' @title Plot Estimated Latent Variables Against Covariates
#' @describeIn LatentVariablePlots Plot estimated latent variable values against covariate values for occupancy
#' @param theta a vector of model parameters with BUGS names
#' @param fit A fitted runjags object.
#' @param esttype When fit is non-NULL, then esttype is used to extract parameter values
#' @param covar a dataframe of covariate values. It must have a column labelled 'ModelSite' that gives the ModelSite of covariate value.
#' @param aggregatefcn Rows in covar with duplicate ModelSite values are aggregrated using this function
#' @export
plot_LVvscovar.fit <- function(fit, esttype = "median", theta = NULL, covar, facetvars = NULL, cuts = 3, aggregatefcn = mean){
  df <- plotdf_LVvscovar.fit(fit, esttype = esttype, theta = theta, covar = covar, aggregatefcn = aggregatefcn)
 
  if (!is.null(facetvars)){
    df <- df %>%
      mutate_at(facetvars, ~cut(., cuts))
  } 
  dflong <- df %>%
    tidyr::pivot_longer(starts_with("LV"), names_to = "LV Name", values_to = "LV Value") %>%
    tidyr::pivot_longer(setdiff(names(covar), c("ModelSite", facetvars)), names_to = "Covariate Name", values_to = "Covariate Value")
  if (is.null(facetvars)){
    plt <- dflong %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(aes_(y = ~`LV Value`, x = ~`Covariate Value`)) +
      ggplot2::facet_wrap(vars(`LV Name`, `Covariate Name`), scales = "free") +
      ggplot2::geom_hline(yintercept = 0, col = "blue", lty = "dashed") +
      ggplot2::geom_smooth(aes_(y = ~`LV Value`, x = ~`Covariate Value`), method = "gam", level = 0.95, formula = y ~ s(x, bs = "cs")) 
  } else {
    plt <- dflong %>%
      ggplot2::ggplot() +
      ggplot2::geom_point(aes_(y = ~`LV Value`, x = ~`Covariate Value`, col = as.name(facetvars))) +
      ggplot2::facet_wrap(vars(`LV Name`, `Covariate Name`), scales = "free") +
      ggplot2::geom_hline(yintercept = 0, col = "blue", lty = "dashed") +
      ggplot2::geom_smooth(aes_(y = ~`LV Value`, x = ~`Covariate Value`, col = as.name(facetvars)), method = "gam", level = 0.95, formula = y ~ s(x, bs = "cs")) 
  }
  return(plt)
}

#' @describeIn LatentVariablePlots Prepare plotting data frames for latent variable values against covariate values for occupancy
#' @param theta a vector of model parameters with BUGS names
#' @param fit A fitted runjags object.
#' @param esttype When fit is non-NULL, then esttype is used to extract parameter values
#' @param covar a dataframe of covariate values. It must have a column labelled 'ModelSite' that gives the ModelSite of covariate value.
#' @param aggregatefcn Rows in covar with duplicate ModelSite values are aggregrated using this function
#' @export
plotdf_LVvscovar.fit <- function(fit, esttype = "median", theta = NULL, covar, aggregatefcn = mean){
  if (is.null(theta)){theta <- get_theta(fit, type = esttype)}
  if (anyDuplicated(covar[, "ModelSite"]) > 0){
    warning("Multiple rows in 'covar' have the same ModelSite. These rows will be aggregated using aggregatefcn")
    covar <- covar %>%
      as_tibble() %>%
      group_by(ModelSite) %>%
      summarise_all(aggregatefcn)
  }
  
  fitdata <- as_list_format(fit$data)
  ## LV values
  LV <- bugsvar2matrix(theta, "LV", 1:fitdata$J, 1:fitdata$nlv) # rows are model sites, columns are latent variables
  colnames(LV) <- paste0("LV", 1:ncol(LV))
  LV <- cbind(ModelSite = 1:nrow(LV), LV) %>% as_tibble() 
  df <- dplyr::left_join(LV, covar, by = "ModelSite", suffix = c(".LV", ".X"))
  return(df)
}