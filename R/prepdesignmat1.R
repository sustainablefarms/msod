prep.designmatprocess_v1 <- function(indata, fmla, stoponhighcorrelation = FALSE, ...){
  designmat1 <- model.matrix(as.formula(fmla), as.data.frame(indata))
  
  ## Check correlation between covariates
  if (sum(colnames(designmat1) != "(Intercept)") >= 2){
    cormat <- cor(designmat1[, colnames(designmat1) != "(Intercept)"])
    diag(cormat) <- NA
    if (max(abs(cormat), na.rm = TRUE) > 0.75) {
      # cormat[upper.tri(cormat)] <- NA
      # highcorr <- which(abs(cormat) > 0.2, arr.ind = TRUE)
      if (stoponhighcorrelation) {stop("Very high correlation between covariates")}
      else {warning("Very high correlation between covariates")}
    }
  }
  
  c_n_s <- get_center_n_scale(designmat1)
  out <- c(
    fmla = fmla,
    c_n_s,
    version = 1
  )
  class(out) <- c("designmatprocess", class(out))
  return(out)
}

apply.designmatprocess_v1 <- function(designmatprocess, indata){
  designmat1 <- model.matrix(as.formula(designmatprocess$fmla), as.data.frame(indata))
  designmat <- apply_center_n_scale(designmat1, center = designmatprocess$center, scale = designmatprocess$scale)
  return(designmat)
}

unstandardise.designmatprocess_v1 <- function(designmatprocess, indata){
  stopifnot(ncol(indata) == length(designmatprocess$scale))
  stopifnot(all(colnames(indata) == names(designmatprocess$center)))
  uncentered <- Rfast::eachrow(Rfast::eachrow(indata, designmatprocess$scale, oper = "*"),
                                 designmatprocess$center, oper = "+")
  colnames(uncentered) <- names(designmatprocess$center)
  return(data.frame(uncentered))
}
