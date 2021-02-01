#' @title Compare Expected Number of Species at ModelSites
#' @details Very similar to elpd_compare(). 
#' Assume each ModelSite with detection parameters are iid samples from a common distribution.
#' In this regime the expected number of species for a ModelSite and the observed number of species, are both random variables.
#' 
#' For a correct model and given covariates for a ModelSite, the difference between the observed number of species and expected number of species is a random variable
#'  with mean 0, and fixed variance.
#' For a ModelSite drawn from the common distribution, the mean is still 0, and the variance is a random variable.
#' The variance of species numbers for a ModelSite drawn at random can be estimated from
#'  the fixed variance computed for many iid ModelSite covariate draws by averaging:
#' \deqn{ V[D] = E[D^2] - E[D]^2 = E[E[D^2|X]] - E[E[D|X]]^2 = E[E[D^2|X]] = E[E[(N - E[N|X])^2|X]] = E[E[N^2|X] - E[N|X]^2] = E[V[N|X]] }
#' The variance of species numbers for a ModelSite drawn at random can also be estimated from many observed ModelSites by averaging:
#' \deqn{ V[D] = V[N - E[N|X]]}
#' The mean difference can also be estimated from observations (and should be zero):
#' \deqn{ 0 = E[D] = E[E[D|X]] ~ \frac{1}{m}\sum_{i = 1}^m d_i }
#' 
#' *No bias correction for use of insample data is included.*
#' @return A matrix with a column for the aggregate (summed) difference of model sites between models,
#'  and the standard error of this difference (computed as the sample standard deviation of difference, multiplied by the square root of the number of ModelSites)
#' @examples 
#' nsites <- 30
#' artfitA <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, modeltype = "jsodm_lv", nlv = 4)
#' artfitB <- artfitA; artfitB$mcmc[[1]] <- artfitB$mcmc[[1]]*0.01 + 10 # a second, poorer model
#' predspecrichA <- speciesrichness(artfitA, occORdetection = "occupancy",
#'                            usefittedlvv = FALSE)
#' predspecrichB <- speciesrichness(artfitA, occORdetection = "occupancy",
#'                                  usefittedlvv = FALSE)
#' 
#' ObsSpecRich <- detectednumspec(y = artfitA$data$y, ModelSite = artfitA$data$ModelSite)
#' 
#' Enum_compare_sum <- Enum_compare(ObsSpecRich,
#'                                  data.frame(A = predspecrichA["E", ], B = predspecrichB["E", ]),
#'                                  data.frame(A = predspecrichA["V", ], B = predspecrichB["V", ])
#' )

#' @param observed A list of the number of species observed for each ModelSite
#' @param predicted A dataframe or matrix with each column the expected number of species detected from a single model. Each row is a ModelSite in the same order as [observed].
#' @param predictedV Same as [predicted], but the variance of the number of species detected.
#' @export
Enum_compare <- function(observed, predicted, predictedV){
  resids <- drop(observed) - predicted
  Eresid_fromobs <- colMeans(resids)
  Eresid_frommodel <- Eresid_fromobs * 0
  Vresid_fromobs <- apply(resids, MARGIN = 2, var)
  Vresid_frommodel <- colMeans(predictedV)
  SE_Eresid_obs_frommodel <- sqrt(Vresid_frommodel/length(drop(observed)))
  SE_Eresid_obs_fromobs <- sqrt(Vresid_fromobs/length(drop(observed)))

  out <- data.frame("E[D]_model" = Eresid_frommodel,
             "E[D]_obs" = Eresid_fromobs,
             "SE(E[D]_obs)_model" = SE_Eresid_obs_frommodel,
             "SE(E[D]_obs)_obs" = SE_Eresid_obs_fromobs,
             "V[D]_model" = Vresid_frommodel,
             "V[D]_obs" = Vresid_fromobs,
             check.names = FALSE
             )
  
  stopifnot(identical(names(Eresid_frommodel), names(Eresid_fromobs)))
  stopifnot(identical(names(Eresid_frommodel), names(SE_Eresid_obs_frommodel)))
  stopifnot(identical(names(Eresid_frommodel), names(SE_Eresid_obs_fromobs)))
  stopifnot(identical(names(Eresid_frommodel), names(Vresid_frommodel)))
  stopifnot(identical(names(Eresid_frommodel), names(Vresid_fromobs)))
  rownames(out) <- make.unique(names(Eresid_frommodel))
  return(out)
}

#' @describeIn Enum_compare The coverage of approximate 95% credence interval of each model.
#' @export
Enum_coverage <- function(observed, predicted, predictedV){
  resids <- drop(observed) - predicted
  in95ci <- abs(resids) < 2 * sqrt(predictedV)
  out <- list(
    mean = colMeans(in95ci),
    in_ci = in95ci
  )
  return(out)
}
