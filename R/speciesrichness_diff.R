#' difference in species richness to reference environmental variables.
#' @examples
#' fit <- readRDS("../sflddata/private/data/testdata/cutfit_7_4_11_2LV.rds")
#' fit <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_model_2lv_e13.rds")
#' fit <- translatefit(fit)
#' Xocc <- unstandardise.designmatprocess(fit$XoccProcess, fit$data$Xocc[10, , drop = FALSE])
#' refXocc <- colMeans(unstandardise.designmatprocess(fit$XoccProcess, fit$data$Xocc[, , drop = FALSE]))
#' refXocc["NMdetected"] <- 1
#' refXocc <- t(refXocc)
#' refXocc <- as.data.frame(refXocc)
#' occspecrich_difference(fit, Xocc, refXocc)
#' refXocc <- Xocc
#' refXocc$woody500m <- 2
#' refXocc$NMdetected <- 1
#' refXocc$MaxTWarmMonth <- 370
#' occspecrich_difference(fit, Xocc, refXocc)

#' @export
occspecrich_difference <- function(fit, Xocc, refXocc, sameLVvalues = TRUE){
  stopifnot(sameLVvalues)
  
  set.seed(3548) #so LV are the same for reference info and pocc
  pocc <- apply_to_new_data(poccupy, fit,
                            Xocc = Xocc,
                            funargs = list(lvvfromposterior = FALSE))
  set.seed(3548) #so LV are the same for reference info and pocc
  refpocc <- apply_to_new_data(poccupy, fit,
                               Xocc = refXocc,
                               funargs = list(lvvfromposterior = FALSE))
  
  Ediffrv <- pocc - refpocc
  Vdiffrv <- pocc * (1 - pocc) + refpocc * (1 - refpocc)
  
  return(prednumsuccess_ErvVrv(Ediffrv, Vdiffrv))
}

