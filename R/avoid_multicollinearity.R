#' @title Removing variables using ViF and correlation
#' @details The function first removes variables based on pairwise correlation, and then based on ViF.
#' 
#' @examples 
indata <- readRDS("./private/data/clean/7_2_10_input_data.rds")
corrthresh <- 0.5
vifthresh <- 10

remove_bycorrvif <- function(fmla, data, corrthresh, vifthresh){
  mat <- model.frame(as.formula("~ -1 + AnnMeanTemp + AnnPrec + MaxTWarmMonth + PrecWarmQ + 
                   MinTColdMonth + PrecColdQ + AnnTempRange + PrecSeasonality + longitude + latitude"),
                     data = indata$insampledata$Xocc)

  # remove correlations above threshold, one at a time, largest correlation to smallest,
  corr_removeinfo <- matrix(, nrow = 0, ncol = 3)
  colnames(corr_removeinfo) <- c("Removed", "Correlation", "KeptPartner")
  corr_removeinfo <- as.data.frame(corr_removeinfo)
  repeat {
    cormat <- cor(mat)
    diag(cormat) <- NA
    if (max(abs(cormat), na.rm = TRUE) < corrthresh) { break }
    max_rowcol <- arrayInd(which.max(abs(cormat)), .dim = dim(cormat)) #row and column of the maximum correlation
    corrtoremove <- cormat[max_rowcol] #logging prep
    nametoremove <- colnames(mat)[max(max_rowcol)] #logging prep
    partnerthatremains <- colnames(mat)[min(max_rowcol)] #logging prep
    corr_removeinfo <- rbind(corr_removeinfo, list(Removed = nametoremove,  #logging prep
                           Correlation = corrtoremove, 
                           KeptPartner = partnerthatremains))
    mat <- mat[, -max(max_rowcol)] #remove the variable which is later in the mat matrix
  }
  
  # apply ViF to the remaining variables, remove one by one
  mat$yran <- rnorm(nrow(mat)) #simulate a y value, its value doesn't actually matter for ViF, just needed to create an lm
  ViF_removeinfo <- matrix(, nrow = 0, ncol = 2) #create an empty data frame for logging
  colnames(ViF_removeinfo) <- c("Removed", "ViF")
  ViF_removeinfo <- as.data.frame(ViF_removeinfo)
  repeat {
    mod <- lm(yran ~ . , data = mat)
    
    # check for NA fitted coefficients
    coefs <- coefficients(mod) #rows are the coefficients, each model is a column
    isna <- is.na(coefs)
    if (sum(isna) > 0){
      stop(paste("Loading for ", names(isna[isna]), " is fitted as NA. Please modify input matrix to compute VIFs."))
    }
    
    # for each of these models compute the generalised variance inflation factors
    gvifs <- car::vif(mod)
    if (max(gvifs) < vifthresh){ break }
    maxind <- which.max(gvifs)
    ViF_removeinfo <- rbind(ViF_removeinfo,
                            list(Removed = names(gvifs)[[maxind]],
                                 ViF = gvifs[[maxind]]))
    mat <- mat[, colnames(mat) != names(gvifs)[[maxind]]]
  }
  
  return(list(
    Kept = colnames(mat),
    Corr_Removed = corr_removeinfo,
    ViF_Removed = ViF_removeinfo
         ))
}