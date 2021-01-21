#' @title Apply function to fitted object after supplanting data with new data
#' 
#' @export
apply_to_new_data <- function(FUN, fit, Xocc, Xobs = NULL, ModelSite = NULL, y = NULL, dataargs = NULL, funargs = NULL){
  fitnew <- do.call(supplant_new_data, c(list(fit = fit, Xocc = Xocc, Xobs = Xobs, ModelSite = ModelSite, y = y),
                                            dataargs), quote = TRUE) #quote somehow gets the generic to notice the class of fit
  out <- do.call(FUN, c(list(fitnew), funargs), quote = TRUE)
  return(out)
}
