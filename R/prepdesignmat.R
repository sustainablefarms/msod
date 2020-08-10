#' @title Preprocessing Input Data
#' @describeIn apply.designmatprocess Calculates the parameters required to build a centred and scaled design matrix from input data.
#' @param indata Input dataframe to be processed.
#' @param fmla A model formula (predictor side only).
#' @param stoponhighcorrelation If TRUE the preparations for a design matrix will fail if correlations between covariates are higher than 0.75
#' @return A special list containing parameters for applying a preprocessing step to data.
#' @export
prep.designmatprocess <- function(indata, fmla, version = 2, ...){
  out <- NULL
  if (version == 2) {out <- prep.designmatprocess_v2(indata, fmla, ...)}
  if (version == 1) {out <- prep.designmatprocess_v1(indata, fmla, ...)}
  stopifnot(!is.null(out))
  return(out)
}

apply.designmatprocess <- function(designmatprocess, indata){
  if (!("designmatprocess" %in% class(designmatprocess))){designmatprocess$version = 1}
  out <- switch(designmatprocess$version,
                apply.designmatprocess_v1(designmatprocess, indata),
                apply.designmatprocess_v2(designmatprocess, indata))
  return(out)
}

# keep variables to keep in the design matrix preparations
# drop variables to forced to drop in the design matrix preparations
prep.designmatprocess_v2 <- function(indata, fmla, keep = NULL, drop = NULL){
  
  fmlaNdata <- computelogsnow(fmla, indata)
  
  # get wanted columns (which be default aren't precomputed)
  ts <- terms(fmlaNdata$fmla, data = fmlaNdata$indata)
  varnames <- rownames(attr(ts, "factor"))
  
  tokens <- unlist(strsplit(varnames, "(I|\\(|\\)|,| )"))
  keep <- union(intersect(tokens, names(fmlaNdata$indata)), keep)
  keep <- setdiff(keep, drop)
  
  # extract wanted columns
  wanteddata <- fmlaNdata$indata[, keep, drop = FALSE]
  
  # check that above extraction got all required data
  tryCatch(mf <- model.frame(fmlaNdata$fmla, wanteddata),
           error = function(e) stop(paste("Didn't parse formula correctly and required columns have been removed.",
                                           "Use argument 'keep' to ensure column remains.", 
                                           e)))
  rm(mf)
  
  # center and scale before computing interactions
  c_n_s <- get_center_n_scale(wanteddata)
  out <- c(
    fmla = fmla,
    c_n_s,
    version = 2
  )
  class(out) <- c("designmatprocess", class(out))
  return(out)
}

apply.designmatprocess_v2 <- function(designmatprocess, indata){
  fmlaNdata <- computelogsnow(designmatprocess$fmla, indata)
  datastd <- apply_center_n_scale(fmlaNdata$indata, designmatprocess$center, designmatprocess$scale)
  designmat <- model.matrix(fmlaNdata$fmla, as.data.frame(datastd))
  return(designmat)
}

# function edits indata and formula so that logged variables are computed NOW
computelogsnow <- function(fmla, indata){
  indata <- as.data.frame(indata)
  fmla <- as.formula(fmla)
  rhschar <- tail(as.character(fmla), 1)
  ts <- terms(fmla, data = indata)
  varnames <- rownames(attr(ts, "factor"))

  ### remove any variables that want standardised BEFORE computing
  ### precompute some variables (like logged variables)
  computenow <- grep("^log\\(", varnames, value = TRUE)
  if (length(computenow) > 0){
    vals <- lapply(computenow, function(x) with(indata, eval(parse(text = x))))
    names(vals) <- computenow
    vals <- do.call(data.frame, c(vals, check.names = TRUE))
    for (i in 1:length(computenow)){
      rhschar <- gsub(computenow[[i]], names(vals)[[i]], rhschar, fixed = TRUE)
    }
    fmla <- reformulate(termlabels = rhschar)
    indata <- cbind(indata, vals)
  }
  return(list(
    fmla = fmla,
    indata = indata
  ))
}

# Gets centres and scales for a matrix/data.frame. Columns that are constant are shifted to 1
# const_tol is the tolerance on the (population) SD which determines whether a column is treated as constant.
get_center_n_scale <- function(indata, const_tol = 1E-8){
  means <- colMeans(indata)
  sds <- ((nrow(indata) - 1) / nrow(indata)) * apply(indata, 2, sd)
  center <- means
  scale <- sds
  isconstant <- (sds < const_tol)
  center[isconstant] <- means[isconstant] - 1 #centering of constant columns to 1
  scale[isconstant] <- 1 #no scaling of constant columns - they are already set to 1
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
#' @describeIn apply.designmatprocess Builds a centred and scaled design matrix from input data.
#' @param designmatprocess Are instructions for preprocessing input data, created by [prep.designmatprocess()] 
#' @param indata Input dataframe to be processed.
#' @details The input data is turned into a design matrix using [stats::model.matrix()].
#' Each non-constant column is then centered and scaled.
#' @return A design matrix.
#' @export
apply.designmatprocess_v1 <- function(designmatprocess, indata){
  designmat1 <- model.matrix(as.formula(designmatprocess$fmla), as.data.frame(indata))
  designmat <- apply_center_n_scale(designmat1, center = designmatprocess$center, scale = designmatprocess$scale)
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
