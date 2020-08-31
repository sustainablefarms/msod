#' @title Extend MCMC chains of fitted model
#' @param fit Model fitted by run.detectionoccupancy
#' @param ... arguments passed to extent.jags
#' @param filename File to save extended chains (with modelling information) to.
#' @export
extend.fit <- function(fit, ..., filename = NULL){
  if (!is.null(filename)){checkwritable(filename)} # check that file can be written before continuing
  
  # save fit$data in list form
  data.list <- fit$data
  
  fit$data <- dump.format(fit$data)
  fit2 <- extend.jags(fit, ...)

  # add summary of parameter distributions
  if (fit2$sample >= 100) {
    fit2 <- add.summary(fit2)
    fit2$crosscorr <- "Crosscorrelation removed to conserve disk size. See ?add.summary to compute it."
  }
  
  # pass meta information
  fit2$data <- list.format(fit2$data)
  colnames(fit2$data$y) <- fit$species
  colnames(fit2$data$Xocc) <- colnames(data.list$Xocc)
  rownames(fit2$data$Xocc) <- 1:nrow(data.list$Xocc)
  colnames(fit2$data$Xobs) <- colnames(data.list$Xobs)
  rownames(fit2$data$Xobs) <- fit$ModelSite
  
  
  # attach data preparation methods
  fit2$XoccProcess <- fit$XoccProcess
  fit2$XobsProcess <- fit$XobsProcess
  fit2$ModelSite <- fit$ModelSite
  fit2$species <- fit$species
  if (!is.null(filename)){try(saveRDS(fit2, filename)) }
  return(invisible(fit2))
}
