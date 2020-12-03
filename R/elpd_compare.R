#' @title Compare Sitewise Differences in Log Predictive Density; Expected Difference with Standard Errors
#' @details Very similar to loo::loo_compare(), but with fewer checks
#' @param elpds A data.frame or matrix of expected log posterior densities. Each row is ModelSite, each column is a Model. Typically the columns will have names.
#' @param refname The column name (i.e. Model) to compare all other models to. If NULL then uses the model with the highest sum of elpd.
#' @return A matrix with a column for the aggregate (summed) difference of model sites between models,
#'  and the standard error of this difference (computed as the sample standard deviation of difference, multiplied by the square root of the number of ModelSites)
#' @examples 
#' waics_l <- readRDS("./tmpdata/7_3_01_waics.rds")
#' mat1 <- elpd_compare(elpds = do.call(cbind, lapply(waics_l, function(x) x$loo$pointwise[, 1])))
#' matref <- loo::loo_compare(lapply(waics_l, function(x) x$loo))
#' all(mat1 == matref[, c("elpd_diff", "se_diff")])
#' @export
elpd_compare <- function(elpds, refname = NULL){
  # elpds <- do.call(cbind, lapply(x, function(x) x$pointwise[, 1]))
  if (!is.null(refname)){
    refidx <- which(colnames(elpds) == refname)
  } else {
    refidx <- which.max(colSums(elpds))
  }
  lpds_diff_df <- elpds %>%
    dplyr::as_tibble() %>%
    dplyr::mutate_all(~ . - elpds[, refidx, drop = TRUE])
  diff <- lpds_diff_df %>%
    dplyr::summarise_all(sum)
  sediff <- lpds_diff_df %>%
    dplyr::summarise_all(~sd(.) * sqrt(nrow(lpds_diff_df)))
  compare_summary <- bind_cols(elpd_diff = t(diff), se_diff = t(sediff)) %>%
    dplyr::mutate(model =  names(sediff)) %>%
    dplyr::arrange(-elpd_diff)
  mat <- as.matrix(compare_summary[, 1:2])
  rownames(mat) <- compare_summary$model
  return(mat)
}