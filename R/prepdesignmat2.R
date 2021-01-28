# keep variables to keep in the design matrix preparations
# drop variables to forced to drop in the design matrix preparations
# if model only has intercept then scale and center parameters have length 0
# computes and stardises log(vars), but not anything else in I... (just standardises components if possible)
# Squares in interactions (e.g. I(UpSite^2) : Sine1) are not standardised, but the components are
prep.designmatprocess_v2 <- function(indata, fmla, keep = NULL, drop = NULL, preserve = NULL){
  fmlaNdata <- computelogsnow(fmla, indata)
  
  # get wanted columns (which be default aren't precomputed)
  ts <- terms(fmlaNdata$fmla, data = fmlaNdata$indata)
  varnames <- rownames(attr(ts, "factor"))
  
  tokens <- NULL
  if (length(varnames) > 0){tokens <- unlist(strsplit(varnames, "(^I\\(|\\(|\\)|,|`|\\^|\\+)"))}
  keep <- union(intersect(tokens, names(fmlaNdata$indata)), keep)
  keep <- setdiff(keep, drop)
  
  # extract wanted columns
  wanteddata <- fmlaNdata$indata[, keep, drop = FALSE]
  
  # check that above extraction got all required data
  tryCatch(mf <- model.matrix(fmlaNdata$fmla, data = wanteddata),
           error = function(e) stop(paste("Didn't parse formula correctly and required columns have been removed.",
                                           "Use argument 'keep' to ensure column remains.", 
                                           e)))
  rm(mf)
  
  # center and scale before computing interactions
  c_n_s <- get_center_n_scale(wanteddata, preserve = preserve)
  out <- c(
    fmla = fmla,
    c_n_s,
    version = 2
  )
  class(out) <- c("designmatprocess", class(out))
  return(out)
}

apply.designmatprocess_v2 <- function(designmatprocess, indata){
  stopifnot(designmatprocess$version == 2)
  fmlaNdata <- computelogsnow(designmatprocess$fmla, indata)
  datastd <- apply_center_n_scale(fmlaNdata$indata, designmatprocess$center, designmatprocess$scale)
  origNAaction <- options("na.action" = "na.pass") #setting so that model.matrix preserves the NA values
  designmat <- model.matrix(fmlaNdata$fmla, as.data.frame(datastd), na.rm = FALSE)
  options(origNAaction) # returning setting
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

undologsnow <- function(indata){
  indata <- as.data.frame(indata)
  
  ### remove any variables that want standardised BEFORE computing
  ### precompute some variables (like logged variables)
  undolog <- grep("^log\\.", names(indata), value = TRUE)
  if (length(undolog) > 0){
    origvals <- exp(indata[, undolog, drop = FALSE])
    names(origvals) <- gsub("^log\\.", "", names(origvals))
    names(origvals) <- gsub("\\.$", "", names(origvals))
    indata <- indata[, setdiff(names(indata), undolog), drop = FALSE]
    neworigvalnames <- setdiff(names(origvals), names(indata))
    indata <- cbind(indata, origvals[, neworigvalnames, drop = FALSE])
  }
  return(indata)
}


unstandardise.designmatprocess_v2 <- function(designmatprocess, indata){
  stopifnot(designmatprocess$version == 2)
  uncentered <- Rfast::eachrow(Rfast::eachrow(indata[, names(designmatprocess$scale), drop = FALSE], designmatprocess$scale, oper = "*"),
                               designmatprocess$center, oper = "+")
  colnames(uncentered) <- names(designmatprocess$center)
  unlogged <- undologsnow(uncentered)
  return(unlogged)
}