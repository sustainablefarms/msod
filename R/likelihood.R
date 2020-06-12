#' @title Computing Likelihoods for Occupancy Detection Models


#' @details Any predictinve accuracy measure requires a choice of 
#' 1. the part of the model that is considered the 'likelihood' and 
#' 2. factorisation of the likelihood into 'data points' [Vehtari 2017]
#' 
#' On 1: New data will look like a new location or visit for a new season in our exisitng region, and observing only the species included in the model.
#' This means we have zero knowledge of the latent variable value at the new ModelSite. This means likelihood:
#'         *  conditional on the covariates u.b and v.b (not using the fitted values of mu.u.b, tau.u.b etc)
#'         *  is conditional on the lv.coef values of each species
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
#' covars <- simulate_covar_data(nsites = 50, nvisitspersite = 2)
#' y <- simulate_iid_detections(3, nrow(covars$Xocc))
#' 
#' fittedmodel <- run.detectionoccupany(
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
#' outofsample_covars <- simulate_covar_data(nsites = 10, nvisitspersite = 2)
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
#' lpds: a list of the log likelihood of the observations for each ModelSite in the supplied data
#' lppd: the computed log pointwise predictive density (sum of the lpds). This is equation (5) in Gelman et al 2014
#' @export
lppd.newdata <- function(fit, Xocc, yXobs, ModelSite, chains = 1, numlvsims = 1000, cl = NULL){
  likel.mat <- likelihoods.fit(fit, Xocc = Xocc, yXobs = yXobs, ModelSite = ModelSite,
                               chains = chains, numlvsims = numlvsims, cl = cl)
  likel.marg <- Rfast::colmeans(likel.mat) # the loglikelihood marginalised over theta (poseterior distribution)
  return(
    list(
      lppd = sum(log(likel.marg)),
      lpds = log(likel.marg) # a list of the log likelihood of the observations for each ModelSite in the supplied data
    )
  )
}

#' @describeIn likelihoods.fit Compute the likelihood of observations at each ModelSites. At data in the fitted model, or on new data supplied.
#' @param chains is a vector indicator which mcmc chains to extract draws from
#' @param numlvsims the number of simulated latent variable values to use for computing likelihoods
#' @param cl a cluster created by parallel::makeCluster()
#' @return `likelihoods.fit` returns a matrix. Each row corresponds to a draw of the parameters from the posterior. Each column to a ModelSite
#' Compute the likelihoods of each ModelSite's observations given each draw of parameters in the posterior.
#' @export
likelihoods.fit <- function(fit, Xocc = NULL, yXobs = NULL, ModelSite = NULL, chains = 1, numlvsims = 1000, cl = NULL){
  fit$data <- as.list.format(fit$data)
  draws <- do.call(rbind, fit$mcmc[chains])
  
  if ( (is.null(fit$data$nlv)) || (fit$data$nlv == 0)){ #make dummy lvsim and and 0 loadings to draws
    lvsim <- matrix(rnorm(2 * 1), ncol = 2, nrow = 2) #dummy lvsim vars
    lv.coef.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$n, ncol = 2), "lv.coef")
    lv.coef.draws <- Rfast::rep_row(lv.coef.bugs, nrow(draws))
    colnames(lv.coef.draws) <- names(lv.coef.bugs)
    draws <- cbind(draws, lv.coef.draws)
  } else {
    lvsim <- matrix(rnorm(fit$data$nlv * numlvsims), ncol = fit$data$nlv, nrow = numlvsims) #simulated lv values, should average over thousands
  }
  if (is.null(Xocc)){ #Extract the Xocc, yXobs etc from the fitted object
    Xocc <- cbind(ModelSite = 1:nrow(fit$data$Xocc), fit$data$Xocc)
    yXobs <- cbind(ModelSite = fit$data$ModelSite, fit$data$Xobs, fit$data$y)
    ModelSite <- "ModelSite"
  }
  arraydata.list <- prep_data_by_modelsite.newdata(fit, Xocc, yXobs, ModelSite)
  if (is.null(cl)) {
    likel.l <- lapply(arraydata.list, likelihood_joint_marginal.ModelSiteDataRow, draws = draws, lvsim = lvsim)
  }
  else {
    likel.l <- parallel::parLapply(cl = cl, arraydata.list, likelihood_joint_marginal.ModelSiteDataRow, draws = draws, lvsim = lvsim)
  }
  likel.mat <- do.call(cbind, likel.l) # each row is a draw, each column is a modelsite (which are independent data points)
  return(likel.mat)
}



#' @describeIn likelihoods.fit Given a fitted model, and input data, prepare the data for computing likelihood
#' @param fit is an object created by run.detectionoccupany
#' @param Xocc A dataframe of covariates related to occupancy. One row per ModelSite.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param yXobs A dataframe of species observations (1 or 0) and covariates related to observations. One row per visit.
#' Each column is either a covariate or a species.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param ModelSite A list of column names in y, Xocc and Xobs that uniquely specify the ModelSite. Can be simply a ModelSite index
#' @return An array with each row a model site and elements that are dataframes for each of Xocc, Xobs, y.
#' @export
prep_data_by_modelsite.newdata <- function(fit, Xocc, yXobs, ModelSite){
  data.list <- prep_new_data(fit, Xocc, yXobs, ModelSite)
  data <- prep_data_by_modelsite(data.list$Xocc, data.list$Xobs, data.list$y, data.list$ModelSite)
}

#' @describeIn likelihoods.fit Given data prepared by prep.data (or prep_new_data),
#'  convert to an array with each row a model site and elements that each a dataframe for 
#'  the Xocc, Xobs, y. ModelSite must be a vector that indicates the row in Xocc corresponding to the observation in Xobs 
prep_data_by_modelsite <- function(Xocc, Xobs, y, ModelSite, outformat = "list"){
  Xocc <- Xocc %>%
    as_tibble(.name_repair = "minimal") %>%
    rowid_to_column(var = "ModelSite") %>%
    tidyr::nest(Xocc = -ModelSite)
  Xobs <- Xobs %>%
    as_tibble(.name_repair = "minimal") %>%
    mutate(ModelSite = ModelSite) %>%
    tidyr::nest(Xobs = -ModelSite)
  y <- y %>%
    as_tibble(.name_repair = "minimal") %>%
    mutate(ModelSite = ModelSite) %>%
    tidyr::nest(y = -ModelSite)
  data <- inner_join(Xocc, Xobs, by = "ModelSite", suffix = c("occ", "obs")) %>%
    inner_join(y, by = "ModelSite", suffix = c("X", "y"))
  if (outformat == "list") {data <- lapply(1:nrow(data), function(i) data[i,, drop = FALSE])}
  return(data)
}

#' @describeIn likelihoods.fit Compute the joint-species LV-marginal likelihood for a ModelSite
#' @param draws A large matrix. Each column is a model parameter, with array elements named according to the BUGS naming convention.
#' Each row of \code{draws} is a simulation from the posterior.
#' @param data_i A row of a data frame created by \code{prep_data_by_modelsite}. Each row contains data for a single ModelSite. 
#' @param lvsim A matrix of simulated LV values. Columns correspond to latent variables, each row is a simulation
#' @export
likelihood_joint_marginal.ModelSiteDataRow <- function(data_i, draws, lvsim, cl = NULL){
  Xocc <- data_i[, "Xocc", drop = TRUE][[1]]
  Xobs <- data_i[, "Xobs", drop = TRUE][[1]]
  y <- data_i[, "y", drop = TRUE][[1]]

  if (is.null(cl)){
    Likl_margLV <- apply(draws, 1, function(theta) likelihood_joint_marginal.ModelSite.theta(
      Xocc, Xobs, y, theta, lvsim))
  } else {
    parallel::clusterExport(cl, list("Xocc", "Xobs", "y", "lvsim"))
    parallel::clusterEvalQ(cl, library(dplyr))
    Likl_margLV <- parallel::parApply(cl, draws, 1, function(theta) likelihood_joint_marginal.ModelSite.theta(
      Xocc, Xobs, y, theta, lvsim))
  }
  return(Likl_margLV)
}

#' @describeIn likelihoods.fit Compute the joint-species LV-marginal likelihood for a ModelSite
#' @param Xocc A matrix of occupancy covariates. Must have a single row. Columns correspond to covariates.
#' @param Xobs A matrix of detection covariates, each row is a visit.
#' @param y A matrix of detection data for a given model site. 1 corresponds to detected. Each row is visit, each column is a species.
#' @param theta A vector of model parameters, labelled according to the BUGS labelling convention seen in runjags
#' @param lvsim A matrix of simulated LV values. Columns correspond to latent variables, each row is a simulation
#' @export
likelihood_joint_marginal.ModelSite.theta <- function(Xocc, Xobs, y, theta, lvsim){
stopifnot(nrow(Xocc) == 1)
stopifnot(nrow(Xobs) == nrow(y))
y <- as.matrix(y)
Xocc <- as.matrix(Xocc)
Xobs <- as.matrix(Xobs)
u.b <- bugsvar2array(theta, "u.b", 1:ncol(y), 1:ncol(Xocc))[,,1]  # rows are species, columns are occupancy covariates
v.b <- bugsvar2array(theta, "v.b", 1:ncol(y), 1:ncol(Xobs))[,,1]  # rows are species, columns are observation covariates
lv.coef <- bugsvar2array(theta, "lv.coef", 1:ncol(y), 1:ncol(lvsim))[,,1] # rows are species, columns are lv
sd_u_condlv <- sqrt(1 - rowSums(lv.coef^2)) #for each species the standard deviation of the indicator random variable 'u', conditional on values of LV

## Probability of Detection, CONDITIONAL on occupied
Detection.LinPred <- as.matrix(Xobs) %*% t(v.b)
Detection.Pred.Cond <- exp(Detection.LinPred) / (exp(Detection.LinPred) + 1)   #this is the inverse logit function


## Likelihood (probability) of single visit given occupied
Likl_condoccupied <- Detection.Pred.Cond * y + (1 - Detection.Pred.Cond) * (1 - y) # non-detection is complement of detection probability
# Likl_condoccupied[y == 0] <- (1 - Detection.Pred.Cond)[y == 0]   # non-detection is complement of detection probability

## Joint likelihood (probability) of detections of all visits CONDITIONAL on occupied
Likl_condoccupied.JointVisit <- apply(Likl_condoccupied, 2, prod)

## Likelihood (proability) of y given unoccupied is either 1 or 0 for detections. Won't include that here yet.
NoneDetected <- as.numeric(colSums(y) == 0)

## Probability of Site Occupancy
ModelSite.Occ.eta_external <- as.matrix(Xocc) %*% t(u.b) #columns are species

# probability of occupancy given LV
ModelSite.Occ.eta_LV <- lvsim %*% t(lv.coef) #occupancy contribution from latent variables, performed all together
ModelSite.Occ.eta <- Rfast::eachrow(ModelSite.Occ.eta_LV, ModelSite.Occ.eta_external, oper = "+") #add the external part to each simulation
# Make u standard deviations equal to 1 by dividing other values by sd
# P(u < -ModelSite.Occ.eta) = P(u / sd < -ModelSite.Occ.eta / sd) = P(z < -ModelSite.Occ.eta / sd)
ModelSite.Occ.eta_standardised <- Rfast::eachrow(ModelSite.Occ.eta, sd_u_condlv, oper = "/") 
ModelSite.Occ.Pred.CondLV <- 1 - pnorm(-ModelSite.Occ.eta_standardised, mean = 0, sd = 1)

# likelihood given LV
Likl.JointVisit.condLV <- Rfast::eachrow(ModelSite.Occ.Pred.CondLV, Likl_condoccupied.JointVisit, oper = "*") #per species likelihood, occupied component. Works because species conditionally independent given LV
Likl.JointVisit.condLV <- Likl.JointVisit.condLV + 
  Rfast::eachrow((1 - ModelSite.Occ.Pred.CondLV), NoneDetected, oper = "*") #add probability of unoccupied for zero detections

# combine with likelihoods of detections
Likl.JointVisitSp.condLV <- Rfast::rowprods(Likl.JointVisit.condLV)  # multiply probabilities of each species together because species are conditionally independent

# take mean of all LV sims to get likelihood marginalised across LV values
Likl_margLV <- mean(Likl.JointVisitSp.condLV)

return(Likl_margLV)
}

# #### TIMING WORK ####
# library(loo)
# waic <- loo::waic(pdetect_joint_marginal.data_i,
#                   data = data[1:10, ],
#                   draws = draws[1:10, ],
#                   lvsim = lvsim)
# Above took 10 000 milliseconds on first go.
# After bugsvar2array faster, took 8000ms. Could pool bugsvar2array work to be even faster.
# After avoiding all dataframe use, dropped to 3000ms
# Can do all of JointSpVst_Liklhood.LV()'s work as matrix manipulations, dropped to 1800ms
# Down to 860ms: replaced use of "rep" with Rfast's functions eachrow, and also replaced row product with Rfast::rowprod. 
