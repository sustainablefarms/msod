#' @title Convert old fit objects to use new parameter naming convention
#' @description Change the variable names to new variables names
#' @export translatefit
#' @examples 
#' model2lv <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_model_2lv_e13.rds")
#' model2lv_new <- translatefit(model2lv)
#' bestmodel_nolv <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_allhyp_vif_logwoody500m_msnm_year_Time_Wind.rds")
#' bestmodel_nolv_new <- translatefit(bestmodel_nolv)
#' 
# places to change in fitted model object: mcmc, summary, summaries, monitor, some in $data too, model
# places to change in codebase: everywhere
parnames <- matrix(c(
  "u.b", "occ.b",
  "lv.coef", "lv.b",
  "LV", "lv.v",
  "v.b", "det.b",
  "mu.u.b",    "mu.occ.b",
  "tau.u.b",   "tau.occ.b",
  "sigma.u.b", "sigma.occ.b",
  "mu.v.b",    "mu.det.b",
  "tau.v.b",   "tau.det.b",
  "sigma.v.b", "sigma.det.b",
  "n", "nspecies",
  "J", "nmodelsites",
  "Vvisits", "nvisits",
  "ModelSite", "ModelSite",
  "Vocc", "noccvar",
  "Vobs", "nobsvar",
  "nlv", "nlv",
  "zeroLV", "zero.lv.v",
  "invSigmaLV", "invSigma.lv.v",
  "lv.covspatscale", "lv.v.spatscale",
  "lv.covtimescale", "lv.v.timescale",
  "z", "occ.v",
  "u", "occ.indicator"
  ),
  byrow = TRUE,
  ncol = 2
)
colnames(parnames) <- c("old", "new")
old2new <- tibble::deframe(data.frame(parnames))

replacebeforebracket <- function(x){
  for (i in 1:nrow(parnames)){
    x <- gsub(paste0("^", parnames[i,"old"], "\\["), paste0(parnames[i,"new"], "\\["), x)
  }
  return(x)
}

replacebetweenticks <- function(x){
  for (i in 1:nrow(parnames)){
    x <- gsub(paste0("`", parnames[i,"old"], "`"), paste0("`", parnames[i,"new"], "`"), x)
  }
  return(x)
}

translatefit <- function(fit){
  #mcmc column names
  newcolnames <- replacebeforebracket(colnames(fit$mcmc[[1]]))
  fit$mcmc <- lapply(fit$mcmc, function(x){
    colnames(x) <- newcolnames
    return(x)
  })
  
  #deviance.table, pd and more don't know structure
  
  # end.state
  fit$end.state <- replacebetweenticks(fit$end.state)
  
  #samplers
  fit$samplers$Node <- replacebeforebracket(fit$samplers$Node)
  
  # burnin, sample, thin pure integers with out names
  # model - this is impossible to do automatically. Must use another place.
  if (!is.null(fit$data$nlv) && (fit$data$nlv > 0)){
    modelfile <- system.file("modeldescriptions", "jsodm_lv.txt", package = "msod")
    class(fit) <- c("jsodm_lv", class(fit))
  } else {
    modelfile <- system.file("modeldescriptions", "jsodm.txt", package = "msod")
    class(fit) <- c("jsodm", class(fit))
  }
  modellines <- readLines(modelfile)
  model <- paste0(modellines, collapse = "\n")
  class(model) <- "runjagsmodel"
  fit$model <- model
  
  # data
  fit$data <- as_list_format(fit$data)
  newnames <- old2new[names(fit$data)]
  newnames[is.na(newnames)] <- names(fit$data)[is.na(newnames)]
  names(fit$data) <- newnames
  
  # monitor
  fit$monitor <- old2new[fit$monitor]
  
  # noread.monitor
  fit$noread.monitor
  fit$noread.monitor <- old2new[fit$noread.monitor]
  
  # modules, factories, response, residual, fitted, method, method.options, timetaken, runjags.version  need none?
  
  # not stored: stochastic, hist, dic, 
  
  # summaries
  if (!is.null(fit$summaries)) {rownames(fit$summaries) <- replacebeforebracket(rownames(fit$summaries))}
  
  # summary
  if (!is.null(fit$statistics)) {rownames(fit$summary$statistics) <- replacebeforebracket(rownames(fit$summary$statistics))}
  if (!is.null(fit$quantiles)) {rownames(fit$summary$quantiles) <- replacebeforebracket(rownames(fit$summary$quantiles))}
  
  # HPD
  if (!is.null(fit$HPD)) {rownames(fit$HPD) <- replacebeforebracket(rownames(fit$HPD))}
  if (!is.null(fit$hpd)) {rownames(fit$hpd) <- replacebeforebracket(rownames(fit$hpd))}
  
  # msce
  if (!is.null(fit$sseff)) {names(fit$mcse$sseff) <- replacebeforebracket(rownames(fit$mcse$sseff))}
  if (!is.null(fit$ssd)) {names(fit$mcse$ssd) <- replacebeforebracket(rownames(fit$mcse$ssd))}
  if (!is.null(fit$mcse)) {names(fit$mcse$mcse) <- replacebeforebracket(rownames(fit$mcse$mcse))}
  
  # psrf
  if (class(fit$psrf) == "gelmanwithtarget"){ rownames(fit$psrf$psrf) <- replacebeforebracket(rownames(fit$psrf$psrf))}
  
  # autocorr
  if (!is.null(fit$autocorr)) {colnames(fit$autocorr) <- replacebeforebracket(colnames(fit$autocorr))}
  
  # crosscorr ignored
  
  # ___stochastic
  if (!is.null(fit$truestochastic)) {names(fit$truestochastic) <- replacebeforebracket(names(fit$truestochastic))}
  if (!is.null(fit$semistochastic)) {names(fit$semistochastic) <- replacebeforebracket(names(fit$semistochastic))}
  if (!is.null(fit$nonstochastic)) {names(fit$nonstochastic) <- replacebeforebracket(names(fit$nonstochastic))}
  
  if (!is.null(fit$discrete)) {names(fit$discrete) <- replacebeforebracket(names(fit$discrete))}
  
  # trace, density, histogram, ecdfplot, key, acplot, ccplot ignored
  # ignored: "summary.available" "summary.pars"      "XoccProcess"       "XobsProcess"       "ModelSite"        
  # [50] "species"           "quality"

  return(fit)
}

