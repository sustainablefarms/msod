#' @title Translate fit object into a smaller object for using in a webapp
#' @param fit An object fitted by this package
#' @return An object of the same class as `fit` with nearly all data removed. Only data used for predictions remains.
#' @examples
#' fit <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_model_2lv_e13.rds")
#' fit <- translatefit(fit)
#' fit <- minimise_fit.jsodm_lv(fit)
#' saveRDS(fit, file = "../farm_biodiversity_app/data/model_data.rds")
#' @export
minimise_fit.jsodm_lv <- function(fit){
  fit2 <- fit[c("mcmc", "XoccProcess", "XobsProcess", "species")]
  fit2$mcmc <- do.call(rbind, fit2$mcmc)
  fit2$mcmc <- fit2$mcmc[, grep("^(occ.b|det.b|lv.b)", colnames(fit2$mcmc))]
  return(fit2)
}
