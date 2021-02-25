#' @title Removing variables using ViF and correlation
#' @details The function first removes variables based on pairwise correlation, and then based on ViF.
#' Variables are removed one at a time.
#' First a variable is removed due to having high correlation, then pairwise correlation is recomputed.
#' This is repeated until no pairwise correlations are above the threshold `corrthresh`.
#' Then generalised Variance Inflation Factors (ViF) are computed using [car::vif()]. 
#' The variable with the highest ViF is removed and ViFs are recomputed.
#' This is repeated until there are no ViFs higher than `vifthresh`.
#' 
#' @examples 
#' indata <- readRDS("./private/data/clean/7_2_10_input_data.rds")
#' remove_bycorrvif("~ AnnMeanTemp + AnnPrec + MaxTWarmMonth + PrecWarmQ + 
#'                    MinTColdMonth + PrecColdQ + PrecSeasonality + longitude * latitude",
#'                  data = indata$insampledata$Xocc,
#'                  corrthresh = 0.9,
#'                  vifthresh = 30)

#' @param fmla A model formula, specifies a possible set of main effects
#' @param data A data frame to extract a the main effects from
#' @param corrthresh A threshold.
#' The variable with the highest correlation, and appearing later in the model matrix, 
#' is removed until there are no pairwise correlations above `corrthresh`.
#' @param vifthresh A threshold. The variable with the highest ViF is removed until no variables have ViF above `vifthresh`.
#' @param centrescalemains If TRUE then [prep.designmatprocess()] and [apply.designmatprocess()] are used to centre and scale main effects (after any logarithms).
#' @export
remove_bycorrvif <- function(fmla, data, corrthresh, vifthresh, centrescalemains = FALSE){
  if (centrescalemains) {
    prepprocess <- prep.designmatprocess(data, fmla)
    mat <- apply.designmatprocess(prepprocess, data)
  } else {
    mat <- model.matrix(as.formula(fmla),
                       data = data)
  }
  mat <- mat[, colnames(mat) != "(Intercept)"] #remove the intercept column

  # remove correlations above threshold, one at a time, largest correlation to smallest,
  corr_removeinfo <- matrix(nrow = 0, ncol = 3)
  colnames(corr_removeinfo) <- c("Removed", "Correlation", "KeptPartner")
  corr_removeinfo <- as.data.frame(corr_removeinfo)
  repeat {
    cormat <- cor(mat)
    diag(cormat) <- NA
    if (max(abs(cormat), na.rm = TRUE) <= corrthresh) { break }
    max_rowcol <- arrayInd(which.max(abs(cormat)), .dim = dim(cormat)) #row and column of the maximum correlation
    corrtoremove <- cormat[max_rowcol] #logging prep
    nametoremove <- colnames(mat)[max(max_rowcol)] #logging prep
    partnerthatremains <- colnames(mat)[min(max_rowcol)] #logging prep
    corr_removeinfo <- rbind(corr_removeinfo, list(Removed = nametoremove,  #logging prep
                           Correlation = corrtoremove, 
                           KeptPartner = partnerthatremains))
    mat <- mat[, -max(max_rowcol)] #remove the variable which is later in the mat matrix
  }
  correlationremained <- cor(mat)
  # print(correlationremained)
  
  # apply ViF to the remaining variables, remove one by one
  mat <- cbind(mat, yran = rnorm(nrow(mat))) #simulate a y value, its value doesn't actually matter for ViF, just needed to create an lm
  mat <- data.frame(mat, check.names = FALSE) #need a dataframe for lm()
  ViF_removeinfo <- matrix(nrow = 0, ncol = 2) #create an empty data frame for logging
  colnames(ViF_removeinfo) <- c("Removed", "ViF")
  ViF_removeinfo <- as.data.frame(ViF_removeinfo)
  repeat {
    mod <- lm(yran ~ . , data = mat)
    
    # check for NA fitted coefficients
    coefs <- coefficients(mod) #rows are the coefficients, each model is a column
    isna <- is.na(coefs)
    if (sum(isna) > 0){
      stop(paste("Fitted loading for", names(isna[isna]), "is NA. Please modify input matrix to compute VIFs."))
    }
    
    # for each of these models compute the generalised variance inflation factors
    gvifs <- car::vif(mod)
    if (max(gvifs) <= vifthresh){ break }
    maxind <- which.max(gvifs)
    nametoremove <- names(gvifs)[[maxind]]
    nametoremove <- gsub("`", "", nametoremove)
    ViF_removeinfo <- rbind(ViF_removeinfo,
                            list(Removed = nametoremove,
                                 ViF = gvifs[[maxind]]))
    stopifnot(!any(duplicated(ViF_removeinfo$Removed))) #error if hasn't been able to remove the relevant column
    mat <- mat[, colnames(mat) != nametoremove]
    if (ncol(mat) <= 2){ break } #means there is only one covariate in the model
  }
  mat <- mat[, colnames(mat) != "yran"]
  
  return(list(
    Kept = colnames(mat),
    Corr_Removed = corr_removeinfo,
    Corr_Remained = correlationremained,
    ViF_Removed = ViF_removeinfo,
    ViF_Remained = gvifs
         ))
}
