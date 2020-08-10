#' @title Preprocessing Input Data
#' @describeIn apply.designmatprocess Calculates the parameters required to build a centred and scaled design matrix from input data.
#' @param indata Input dataframe to be processed.
#' @param fmla A model formula (predictor side only).
#' @param stoponhighcorrelation If TRUE the preparations for a design matrix will fail if correlations between covariates are higher than 0.75
#' @return A special list containing parameters for applying a preprocessing step to data.
#' @export
prep.designmatprocess <- function(indata, fmla, stoponhighcorrelation = FALSE){
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
  
  means <- colMeans(designmat1)
  sds <- ((nrow(designmat1) - 1) / nrow(designmat1)) * apply(designmat1, 2, sd)
  center <- means
  scale <- sds
  isconstant <- (sds < 1E-8)
  center[isconstant] <- means[isconstant] - 1 #centering of constant columns to 1
  scale[isconstant] <- 1 #no scaling of constant columns - they are already set to 1
  preprocessobj <- list(fmla = fmla, center = center, scale = scale)
  return(preprocessobj)
}
#' @describeIn apply.designmatprocess Builds a centred and scaled design matrix from input data.
#' @param designmatprocess Are instructions for preprocessing input data, created by [prep.designmatprocess()] 
#' @param indata Input dataframe to be processed.
#' @details The input data is turned into a design matrix using [stats::model.matrix()].
#' Each non-constant column is then centered and scaled.
#' @return A design matrix.
#' @export
apply.designmatprocess <- function(designmatprocess, indata){
  designmat1 <- model.matrix(as.formula(designmatprocess$fmla), as.data.frame(indata))
  designmat <- scale(designmat1, center = designmatprocess$center, scale = designmatprocess$scale)
  return(designmat)
}

#' @describeIn apply.designmatprocess Builds a centred and scaled design matrix from input data.
#' @param designmatprocess Are instructions for preprocessing input data, created by [prep.designmatprocess()] 
#' @param data Dataframe to be processed.
#' @return The columns of indata before centering and scaling 
#' @export
uncentre.designmatprocess <- function(designmatprocess, indata){
  stopifnot(ncol(indata) == length(designmatprocess$scale))
  stopifnot(all(colnames(indata) == names(designmatprocess$center)))
  uncentered <- Rfast::eachrow(Rfast::eachrow(indata, designmatprocess$scale, oper = "*"),
                                 designmatprocess$center, oper = "+")
  colnames(uncentered) <- names(designmatprocess$center)
  return(uncentered)
}
