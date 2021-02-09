#' @title Species richness as number of species more likely than not to occur on the patch / farm.
#' @description For each draw of lvv and parameters, computes the number of species with any probability in any patch is greater than 0.5
#' @examples
#' fit <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_model_2lv_e13.rds")
#' fit <- translatefit(fit)
#' Xocc <- unstandardise.designmatprocess(fit$XoccProcess, fit$data$Xocc[c(1, 500), , drop = FALSE])
#' speciesrichness_mltn(fit, Xocc)
#' @export
speciesrichness_mltn <- function(fit, Xocc){
  pocc <- apply_to_new_data(poccupy, fit,
                            Xocc = Xocc,
                            funargs = list(lvvfromposterior = FALSE))
  pocc <- apply(pocc, MARGIN = c(2, 3), max) #maximum probability across all sites in Xocc
  nmltn <- apply(pocc > 0.5, MARGIN = 2, sum)
  
  return(quantile(nmltn, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)))
  # low_upp <- hpd_narray(nmltn, prob = prob, drawdim = 2)
  # return(cbind(low_upp, median = Mednmltn))
}
