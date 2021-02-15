#' @title Computing Likelihoods for Occupancy Detection Models


#' @details Any predictinve accuracy measure requires a choice of 
#' 1. the part of the model that is considered the 'likelihood' and 
#' 2. factorisation of the likelihood into 'data points' [Vehtari 2017]
#' 
#' On 1: New data will look like a new location or visit for a new season in our exisitng region, and observing only the species included in the model.
#' This means we have zero knowledge of the latent variable value at the new ModelSite. This means likelihood:
#'         *  conditional on the covariates occ.b and det.b (not using the fitted values of mu.occ.b, tau.occ.b etc)
#'         *  is conditional on the lv.b values of each species
#'         *  is conditional on the latent variable value for (each) new ModelSite being drawn from a standard Gaussian distribution.
#'         
#' On 2: Factoring the likelihood using the inbuilt independence properties of the model means 
#' a single 'data point' is all the data for all visits of a single ModelSite.
#' The likelihood could also be partitioned by each visit, but then data points are dependent (they have the same occupancy value).
#'         
#' The output of [likelihoods.fit()] can be easily passed to [loo::waic()] and [loo::loo()].

# For WAIC:
## function(data_i = data[i, , drop = FALSE], draws = draws)  --> returns a vector, each entry given by draw in draws.
# data: dataframe or matrix containing predictor and observed outcome data. For each observation, i, the ith row of data will be passed to the data_i argument
#       This is like the combination of Xocc joined to Xobs and y via ModelSite?
#       Except the multiple visits to the same ModelSite are *dependent*. Perhaps it is best to combine all visits to a model site!?
# draws: a posterior draws object, passed unaltered to the function
# ...  May be used too, it is passed to each call of the function (all i).
# This function can also be used to perform the PSIS-LOO estimate of PSIS. So long as the rows satisfy conditional independence in the data model.


#' @references A. Vehtari, A. Gelman, and J. Gabry, "Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC," Stat Comput, vol. 27, pp. 1413-1432, Sep. 2017, doi: 10.1007/s11222-016-9696-4.

#' @examples
#' # simulate data
#' covars <- artificial_covar_data(nsites = 50, nvisitspersite = 2)
#' y <- simulate_iid_detections(3, nrow(covars$Xocc))
#' 
#' fittedmodel <- run.detectionoccupancy(
#'   Xocc = covars$Xocc,
#'   yXobs = cbind(covars$Xobs, y),
#'   species = colnames(y),
#'   ModelSite = "ModelSite",
#'   OccFmla = "~ UpSite + Sine1",
#'   ObsFmla = "~ UpVisit + Step",
#'   nlv = 2,
#'   MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 3, thin = 1)
#' )
#' 
#' # run likelihood computations, waic, and psis-loo
#' insamplell <- likelihoods.fit(fittedmodel)
#' waic <- loo::waic(log(insamplell))
#' looest <- loo::loo(log(insamplell), cores = 2)
#' 
#' 
#' 
#' outofsample_covars <- artificial_covar_data(nsites = 10, nvisitspersite = 2)
#' outofsample_y <- simulate_iid_detections(3, nrow(outofsample_covars$Xocc))
#' outofsample_lppd <- lppd.newdata(fittedmodel,
#'              Xocc = outofsample_covars$Xocc,
#'              yXobs = cbind(outofsample_covars$Xobs, outofsample_y),
#'              ModelSite = "ModelSite")
#' 
#' # Recommend using multiple cores:
#' cl <- parallel::makeCluster(2)
#' insamplell <- likelihoods.fit(fittedmodel, cl = cl)
#' 
#' outofsample_lppd <- lppd.newdata(fittedmodel,
#'                                  Xocc = outofsample_covars$Xocc,
#'                                  yXobs = cbind(outofsample_covars$Xobs, outofsample_y),
#'                                  ModelSite = "ModelSite",
#'                                  cl = cl)
#' parallel::stopCluster(cl)

#' @describeIn likelihoods.fit Compute the log pointwise posterior density of new (out-of-sample) data
#' @return `lppd.newdata` returns a list with components
#' lpds: a list of the log predictive density for each ModelSite in the supplied data (for each model site this is the average of the predictive density conditioned on each draw from the posterior).
#' lppd: the computed log pointwise predictive density (sum of the lpds). This is equation (5) in Gelman et al 2014
#' @export
lppd_newdata <- function(fit, Xocc, yXobs, ModelSite, chains = 1, ...){
  warning("Function will become obsolete soon, apply_to_new_data is easier to use")
  likel.mat <- apply_to_new_data(likelihood, fit, Xocc, 
                    Xobs = yXobs[, !(colnames(yXobs) %in% fit$species), drop = FALSE],
                    ModelSite = ModelSite,
                    y = yXobs[, colnames(yXobs) %in% fit$species, drop = FALSE])
  # likel.mat <- likelihood(fit, Xocc = Xocc, yXobs = yXobs, ModelSite = ModelSite, ...)
  elpd_out <- elpd(likel.mat)
  return(elpd_out)
}

#' Expected log predictive density: the probability density for each model site, average across posterior draws
#' @param likelmat A matrix of predictive density for each draw (row) x modelsite (column)
#' @export
elpd <- function(likelmat){
  likel.marg <- Rfast::colmeans(likelmat) # the loglikelihood marginalised over theta (poseterior distribution)
  return(
    list(
      lppd = sum(log(likel.marg)),
      lpds = log(likel.marg) # a list of the log likelihood of the observations for each ModelSite in the supplied data
    )
  )
}

#' @describeIn likelihoods.fit Compute the likelihood of observations at each ModelSite. At data in the fitted model, or on new data supplied.
#' @param chains is a vector indicator which mcmc chains to extract draws from. If NULL then all chains used.
#' @param numlvsims the number of simulated latent variable values to use for computing likelihoods
#' @param cl a cluster created by parallel::makeCluster()
#' @return Returns a matrix. Each row corresponds to a draw of the parameters from the posterior. Each column to a ModelSite.
#' The value in each cell is the probability density, given the parameters from the draw, evaluated at the observations for the model site.
#' @examples 
#' model2lv <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_model_2lv_e13.rds")
#' model2lv_new <- translatefit(model2lv)
#' lkl <- likelihood(model2lv_new, nlvsim = 2)
#' @export
likelihood <- function(fit, ...){
  UseMethod("likelihood")
}

# #### OLD TIMING WORK ####
# library(loo)
# waic <- loo::waic(pdetect_joint_marginal.data_i,
#                   data = data[1:10, ],
#                   draws = draws[1:10, ],
#                   lvsim = lvsim)
# Above took 10 000 milliseconds on first go.
# After bugsvar2array faster, took 8000ms. Could pool bugsvar2array work to be even faster (this has been done as of Oct 6).
# After avoiding all dataframe use, dropped to 3000ms
# Can do all of JointSpVst_Liklhood.LV()'s work as matrix manipulations, dropped to 1800ms
# Down to 860ms: replaced use of "rep" with Rfast's functions eachrow, and also replaced row product with Rfast::rowprod. 
