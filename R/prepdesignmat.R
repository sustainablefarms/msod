#' @title Preprocessing Input Data
#' @describeIn apply.designmatprocess Calculates the parameters required to build a centred and scaled design matrix from input data.
#' @param indata Input dataframe to be processed.
#' @param fmla A model formula (predictor side only).
#' @param stoponhighcorrelation If TRUE the preparations for a design matrix will fail if correlations between covariates are higher than 0.75
#' @return A special list containing parameters for applying a preprocessing step to data.
#' @export
prep.designmatprocess <- function(indata, fmla, version = 2, ...){
  out <- NULL
  out <- switch(version,
                prep.designmatprocess_v1(indata, fmla, ...),
                prep.designmatprocess_v2(indata, fmla, ...))
  stopifnot(!is.null(out))
  return(out)
}

#' @describeIn apply.designmatprocess Builds a centred and scaled design matrix from input data.
#' @param designmatprocess Are instructions for preprocessing input data, created by [prep.designmatprocess()] 
#' @param indata Input dataframe to be processed.
#' @details The input data is turned into a design matrix using [stats::model.matrix()].
#' Each non-constant column is then centered and scaled.
#' @return A design matrix.
#' @export
apply.designmatprocess <- function(designmatprocess, indata){
  if (!("designmatprocess" %in% class(designmatprocess))){designmatprocess$version = 1}
  out <- switch(designmatprocess$version,
                apply.designmatprocess_v1(designmatprocess, indata),
                apply.designmatprocess_v2(designmatprocess, indata))
  return(out)
}

#' @describeIn apply.designmatprocess Uncentres and unscales already standardised data.
#' @param designmatprocess Are instructions for preprocessing input data, created by [prep.designmatprocess()] 
#' @param data Dataframe to be processed.
#' @return The columns of indata before centering and scaling 
#' @export
unstandardise.designmatprocess <- function(designmatprocess, indata){
  if (!("designmatprocess" %in% class(designmatprocess))){designmatprocess$version = 1}
  out <- switch(designmatprocess$version,
                unstandardise.designmatprocess_v1(designmatprocess, indata),
                unstandardise.designmatprocess_v2(designmatprocess, indata))
  return(out)
}

# Gets centres and scales for a matrix/data.frame. Columns that are constant are shifted to 1
# const_tol is the tolerance on the (population) SD which determines whether a column is treated as constant.
get_center_n_scale <- function(indata, const_tol = 1E-8, preserve = NULL){
  means <- colMeans(indata)
  sds <- ((nrow(indata) - 1) / nrow(indata)) * apply(indata, 2, sd)
  center <- means
  scale <- sds
  isconstant <- (sds < const_tol)
  center[isconstant] <- means[isconstant] - 1 #centering of constant columns to 1
  scale[isconstant] <- 1 #no scaling of constant columns - they are already set to 1
  center[names(center) %in% preserve] <- 0 #no centering of preserved columns
  scale[names(scale) %in% preserve] <- 1 #no scaling of preserved columns
  return(list(
    center = center,
    scale = scale
  ))
}

# centre and scale are named vectors
apply_center_n_scale <- function(indata, center, scale){
  stopifnot(names(center) == names(scale))
  stopifnot("matrix" %in% class(indata) || 
              "data.frame" %in% class(indata)
              )
  indata <- indata[, names(center), drop = FALSE]
  out <- scale(indata, center = center, scale = scale)
  return(out)
}

