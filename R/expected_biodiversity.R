#' @title Expected Biodiversity
#' @description The expected number of species occupying a ModelSite, or the expected number of species detected at a ModelSite.
#' 
#' @examples 
#' fit <- readRDS("./tmpdata/7_2_9_addyear_msnm_year_time_2lv.rds")
#' theta <- get_theta(fit, type = 1)
#' Xocc <- fit$data$Xocc[2, , drop = FALSE]
#' Xobs <- fit$data$Xobs[fit$data$ModelSite == 2, , drop = FALSE]
#' numspecies <- fit$data$n
#' lvsim <- matrix(rnorm(2 * 1), ncol = 2, nrow = 2) #dummy lvsim vars
#' lv.coef.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$n, ncol = 2), "lv.coef")
#' theta <- c(theta, lv.coef.bugs)
#' 
#' expectedspeciesnum.ModelSiteIdx(fit, 2, chains = NULL, LVvals = NULL)
#' Enumspec <- predsumspecies(fit, usefittedLV = TRUE)
#' 
#' 
#' indata <- readRDS("./private/data/clean/7_2_10_input_data.rds")
#' predsumspecies_newdata(fit, Xocc <- indata$holdoutdata$Xocc, Xobs = indata$holdoutdata$yXobs, ModelSiteVars = "ModelSiteID", draws, cl = NULL)
#' 
#' @param Xocc A matrix of occupancy covariates. Must have a single row. Columns correspond to covariates.
#' @param Xobs A matrix of detection covariates, each row is a visit. If NULL then expected number of species in occupation is returned
#' @param theta A vector of model parameters, labelled according to the BUGS labelling convention seen in runjags
#' @param LVvals A matrix of LV values. Each column corresponds to a LV. To condition on specific LV values, provide a matrix of row 1.
#' @export
expectedspeciesnum.ModelSite.theta <- function(Xocc, Xobs = NULL, u.b, v.b = NULL, lv.coef = NULL, LVvals = NULL){
  ## Probability of Site Occupancy
  stopifnot(nrow(Xocc) == 1)
  if (is.null(Xobs)) {stopifnot(is.null(v.b))}
  Xocc <- as.matrix(Xocc)
  if (is.null(lv.coef)){ # create dummy LV  
    stopifnot(is.null(LVvals))
    lv.coef <- matrix(0, rows = nrow(u.b), ncol = 2)
    LVvals <- matrix(0, rows = 2, ncol = 2)
    } # dummy LV
  ModelSite.Occ.Pred.CondLV <- poccupy.ModelSite.theta(Xocc, u.b, lv.coef, LVvals)
  
  ## Expected number of species occupying modelsite, given model and theta, and marginal across LV
  EVnumspec_occ <- Erowsum_margrow(ModelSite.Occ.Pred.CondLV)
  names(EVnumspec_occ) <- paste0(names(EVnumspec_occ), "_occ")
  if (is.null(Xobs)){
    return(EVnumspec_occ)
  }
  
  ## Probability of Detection
  if (!is.null(Xobs)){
    Detection.Pred.Cond <- pdetection_occupied.ModelSite.theta(Xobs, v.b) # probability conditional on occupied
    # probability of no detections, given occupied
    NoDetections.Pred.Cond <- Rfast::colprods(1 - Detection.Pred.Cond)

  
    ## probability of no detections, marginal on occupancy
    NoDetections.Pred.marg_1 <- Rfast::eachrow(ModelSite.Occ.Pred.CondLV, NoDetections.Pred.Cond, oper = "*")  #occupied component
    NoDetections.Pred.marg_2 <- 1 - ModelSite.Occ.Pred.CondLV  #plus unoccupied component
    AnyDetections.Pred.marg <- 1 - (NoDetections.Pred.marg_1 + NoDetections.Pred.marg_2)
    EVnumspec_det <- Erowsum_margrow(AnyDetections.Pred.marg)
    names(EVnumspec_det) <- paste0(names(EVnumspec_det), "_det")
    return(c(EVnumspec_occ, EVnumspec_det))
  }
}

#' @title Probability of detection and occupation for a given ModelSite
#' @param v.b Covariate loadings. Each row is a species, each column a detection covariate
#' @return A matrix of detection probabilities for each row of Xobs and each species.
#' Rows in returned value correspond to rows in Xobs, columns correspond to species.
#' @export
pdetection_occupied.ModelSite.theta <- function(Xobs, v.b){
  Xobs <- as.matrix(Xobs)
  Detection.LinPred <- as.matrix(Xobs) %*% t(v.b)
  Detection.Pred.Cond <- exp(Detection.LinPred) / (exp(Detection.LinPred) + 1)   #this is the inverse logit function
  return(Detection.Pred.Cond)
}

#' @describeIn pdetection_occupied.ModelSite.theta The probability of occupation for given LV values.
#' @param u.b Covariate loadings for occupancy. Each row is a species, each column an occupancy covariate.
#' @param lv.coef Loadings for the latent variables. Each row is a species, each column corresponds to a LV.
#' @return A matrix of occupancy probabilities. Each row corresponds to a row of LVvals, each column to a species.
poccupy.ModelSite.theta <- function(Xocc, u.b, lv.coef, LVvals){
  sd_u_condlv <- sqrt(1 - rowSums(lv.coef^2)) #for each species the standard deviation of the indicator random variable 'u', conditional on values of LV
  
  # external
  ModelSite.Occ.eta_external <- as.matrix(Xocc) %*% t(u.b) #columns are species
  
  # probability of occupancy given LV
  ModelSite.Occ.eta_LV <- LVvals %*% t(lv.coef) #occupancy contribution from latent variables, performed all together
  ModelSite.Occ.eta <- Rfast::eachrow(ModelSite.Occ.eta_LV, ModelSite.Occ.eta_external, oper = "+") #add the external part to each simulation
  # Make u standard deviations equal to 1 by dividing other values by sd
  # P(u < -ModelSite.Occ.eta) = P(u / sd < -ModelSite.Occ.eta / sd) = P(z < -ModelSite.Occ.eta / sd)
  ModelSite.Occ.eta_standardised <- Rfast::eachrow(ModelSite.Occ.eta, sd_u_condlv, oper = "/") 
  ModelSite.Occ.Pred.CondLV <- 1 - pnorm(-ModelSite.Occ.eta_standardised, mean = 0, sd = 1)
  return(ModelSite.Occ.Pred.CondLV)
}

#' @param pmat A matrix of success probabilities for Bernoulli random variables.
#' Each column corresponds to an independent Bernoulli random variable, each row is a different set of parameters.
#' @return
#' The expected sum of the Bernoulli random variables, marginal across rows, E[E[N | L]]
#' The variance of the sum, marginal across rows E[N^2] - E[N]^2 = E[V[N | L] + E[N | L]^2] - E[E[N | L]]^2
Erowsum_margrow <- function(pmat){
  ## Expected sum of the Bournilli for each parameter set
  En_row <- Rfast::rowsums(pmat)
  ## Variance of sum for each parameter set (variance adds for species when independent)
  Vn_row <- Rfast::rowsums(pmat * (1 - pmat))
  ## second moment of sum for each row
  M2n_row <- Vn_row + En_row^2
  ## second moment of sum, marginal rows
  M2n <- mean(M2n_row)
  ## 1st moment of sum, marginal rows
  En <- mean(En_row)
  ## Variance of sum, marginal rows
  Vn <- M2n - En^2
  return(c(Esum = En, Vsum = Vn))
}

# use supplied LVs for modelsite
# use fitted LVs for modelsite for each theta
# use simulated LVs for each theta

#' @param UseFittedLV If TRUE the fitted LV variables are used, if false then 1000 LV values are simulated.
#' @param chains The chains of MCMC to use. Default is all chains.
#' @param cl A cluster object created by parallel::makeCluster. If NULL no cluster is used.
#' @return A matrix with each column a ModelSite.
#' Each row is labelled and corresponds to the predicted expection and variance of the number of species occupied or detected.
#' These expectations are with respect to the full posterior distribution of the model parameters, with the exception of the LV values which depends on usefittedLV.
#' @export
predsumspecies <- function(fit, chains = NULL, usefittedLV = TRUE, nLVsim = 1000, cl = NULL){
  fit$data <- as_list_format(fit$data)
  if (is.null(chains)){chains <- 1:length(fit$mcmc)}
  draws <- do.call(rbind, fit$mcmc[chains])
  if ( (is.null(fit$data$nlv)) || (fit$data$nlv == 0)){ #LVs not in model, add dummy variables
    LVbugs <- matrix2bugsvar(matrix(0, nrow = nrow(fit$data$Xocc), ncol = 2), "LV")
    LVbugs.draws <- Rfast::rep_row(LVbugs, nrow(draws))
    colnames(LVbugs.draws) <- names(LVbugs)
    
    lv.coef.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$n, ncol = 2), "lv.coef")
    lv.coef.draws <- Rfast::rep_row(lv.coef.bugs, nrow(draws))
    colnames(lv.coef.draws) <- names(lv.coef.bugs)
    draws <- cbind(draws, lv.coef.draws, LVbugs.draws)
    fit$data$nlv <- 2
    usefittedLV <- TRUE #calculations faster when not simulating 1000s of LV values, especially since they are all ignored here.
  } 
  
  if (usefittedLV){nLVsim <- NULL} #don't pass number of simulations if not going to use them
  
  out <- predsumspecies_raw(
    Xocc = fit$data$Xocc,
    Xobs = fit$data$Xobs,
    ModelSite = fit$data$ModelSite,
    numspecies = fit$data$n,
    nlv = fit$data$nlv,
    draws = draws,
    useLVindraws = usefittedLV,
    nLVsim = nLVsim,
    cl = cl
  )
  return(out)
}

#' @param ModelSite is a list of integers giving the row in Xocc corresponding to a row in Xobs
#' @param Xocc A matrix of occupancy covariates, each row is a ModelSite
#' @param Xobs A matrix of detection covariates. Each row is a visit. The visited ModelSite (row of Xocc) is given by ModelSite
#' @param draws A matrix of posterior parameter draws. Each row is a draw. Column names follow the BUGS naming convention
#' @param useLVindraws Use the LV values corresponding to each draw from within the \code{draws} object.
#' If FALSE nLVsim simulated LV values will be used for each draw.
#' @param nLVsim The number of simulated LV values if not using fitted LV values (only applies if  useLVindraws = FALSE).
#' @details No scaling or centering of Xocc or Xobs is performed by predsumspecies_raw
#' @return A matrix. Each column is a model site, each row is a different summary of the number of species.
#' There will be four rows: the expection and variance of the number of species occupied or detected.
#' @export
predsumspecies_raw <- function(Xocc, Xobs, ModelSite, numspecies, nlv, draws, useLVindraws = TRUE, nLVsim = NULL, cl = NULL){
  # prepare parameters
  nspecall <- numspecies
  ndraws <- nrow(draws)
  nsites <- nrow(Xocc)
  noccvar <- ncol(Xocc)
  nobsvar <- ncol(Xobs)
  ModelSiteIdxs <- ModelSite
  stopifnot(length(ModelSiteIdxs) == nrow(Xobs))
  stopifnot(all(ModelSiteIdxs %in% 1:nsites))
  if (!all(1:nsites %in% ModelSiteIdxs)){warning("Some ModelSite do not have observation covariate information.")}
  u.b <- bugsvar2array(draws, "u.b", 1:nspecall, 1:noccvar)
  v.b <- bugsvar2array(draws, "v.b", 1:nspecall, 1:nobsvar)
  lv.coef <- bugsvar2array(draws, "lv.coef", 1:nspecall, 1:nlv)
  
  if (useLVindraws){stopifnot(is.null(nLVsim))}
  if (!useLVindraws){stopifnot(is.numeric(nLVsim))}

  
  sitedrawidxs <- expand.grid(siteidx = 1:nsites, drawidx = 1:ndraws) 
  
  if (!useLVindraws){ # predicting as if LVs not known, so simulate from their distribution
    lvsim <- matrix(rnorm(nlv * nLVsim), ncol = nlv, nrow = nLVsim)
  } else {
    LVvals <- bugsvar2array(draws, "LV", 1:nsites, 1:nlv)
  }
  

  # for each modelsite and each draw apply the following function:
  Enumspec <- pbapply::pbapply(sitedrawidxs, MARGIN = 1,
        function(sitedrawidx){
          Xocc <- Xocc[sitedrawidx[["siteidx"]], , drop = FALSE]
          Xobs <- Xobs[ModelSiteIdxs == sitedrawidx[["siteidx"]], , drop = FALSE]
          u.b_theta <- matrix(u.b[,, sitedrawidx[["drawidx"]] ], nrow = nspecall, ncol = noccvar)
          v.b_theta <- matrix(v.b[,, sitedrawidx[["drawidx"]] ], nrow = nspecall, ncol = nobsvar)
          lv.coef_theta <- matrix(lv.coef[,, sitedrawidx[["drawidx"]] ], nrow = nspecall, ncol = nlv)
          if (useLVindraws){
            LVvals_thetasite <- matrix(LVvals[sitedrawidx[["siteidx"]], , sitedrawidx[["drawidx"]], drop = FALSE], nrow = 1, ncol = nlv)
          } else {
            LVvals_thetasite <- lvsim
          }
          Enumspec_sitetheta <- expectedspeciesnum.ModelSite.theta(Xocc = Xocc, Xobs = Xobs,
                                                                   u.b = u.b_theta,
                                                                   v.b = v.b_theta,
                                                                   lv.coef = lv.coef_theta,
                                                                   LVvals = LVvals_thetasite)
        },
        cl = cl)
  # each column of Enumspec is row of sitedrawidxs
  Enumspec <- rbind(Enumspec, t(sitedrawidxs))
  
  ## arrange Enumspec in a 3 dimensional array using existing arrangement of Enumspec
  a <- array(Enumspec, dim = c(nrow(Enumspec), nsites, ndraws),
             dimnames = list(
               Summaries = rownames(Enumspec),
               ModelSites = 1:nsites,
               Draws = 1:ndraws
             ))
  Enumspec_drawsitesumm <- aperm(a, perm = c("Draws", "ModelSites", "Summaries")) # arrange Enumspec in a 3 dimensional array: rows = draws, cols = site, depth = summaries
  
  #following checks that conversion is correct
  stopifnot(all(Rfast::eachrow(matrix(Enumspec_drawsitesumm[ , , "siteidx", drop = FALSE],
                                      nrow = dim(Enumspec_drawsitesumm)[[1]], ncol = dim(Enumspec_drawsitesumm)[[2]]), # call to matrix keeps the input in the correct dimension
                               y = as.numeric(1:nsites), oper = "==")))
  stopifnot(sum(Rfast::eachcol.apply(matrix(Enumspec_drawsitesumm[ , , "drawidx", drop = TRUE],
                                            nrow = dim(Enumspec_drawsitesumm)[[1]], ncol = dim(Enumspec_drawsitesumm)[[2]]),
                               y = as.numeric(1:ndraws), oper = "-")) == 0)
  

  # convert predictions for each site and theta into predictions for each site, marginal across theta distribution
  out_det <- EVtheta2EVmarg(
    Vsum = matrix(Enumspec_drawsitesumm[, ,"Vsum_det", drop = FALSE],  nrow = ndraws, ncol = nsites),
    Esum = matrix(Enumspec_drawsitesumm[, , "Esum_det", drop = FALSE], nrow = ndraws, ncol = nsites)
  )
  rownames(out_det) <- paste0(rownames(out_det), "_det")
  out_occ <- EVtheta2EVmarg(
    Vsum = matrix(Enumspec_drawsitesumm[, ,"Vsum_occ", drop = FALSE],  nrow = ndraws, ncol = nsites),
    Esum = matrix(Enumspec_drawsitesumm[, , "Esum_occ", drop = FALSE], nrow = ndraws, ncol = nsites)
  )
  rownames(out_occ) <- paste0(rownames(out_occ), "_occ")
  
  out <- rbind(
   occ = out_occ,
   det = out_det)
  return(out)
}

#' @param Vsum A matrix of variance of sum of species. Each row corresponds to a different theta. Each column a ModelSite.
#' @param Esum A matrix of expected sum of species. Each row corresponds to a different theta. Each column a ModelSite.
#' @return The expectation and variance of the sum of species marginal across theta (across the supplied rows)
EVtheta2EVmarg <- function(Vsum, Esum){
  stopifnot(all(dim(Vsum) == dim(Esum)))
  # per draw the 2nd moments
  M2n_theta <- Vsum + Esum^2
  
  # marginal across draw moments
  M2n <- Rfast::colmeans(M2n_theta)
  En <- Rfast::colmeans(Esum)
  Vn <- M2n - En^2
  return(rbind(
    Esum = En,
    Vsum = Vn
  ))
}

predsumspecies_newdata <- function(fit, Xocc, Xobs, ModelSiteVars, chains = NULL, nLVsim = 1000, cl = NULL){
  datalist <- prep_new_data(fit, Xocc, Xobs, ModelSite = ModelSiteVars)
  usefittedLV <- FALSE # no LV available for new model sites
  
  fit$data <- as_list_format(fit$data)
  if (is.null(chains)){chains <- 1:length(fit$mcmc)}
  draws <- do.call(rbind, fit$mcmc[chains])
  
  if ( (is.null(fit$data$nlv)) || (fit$data$nlv == 0)){ #LVs not in model, add dummy variables
    lv.coef.bugs <- matrix2bugsvar(matrix(0, nrow = fit$data$n, ncol = 2), "lv.coef")
    lv.coef.draws <- Rfast::rep_row(lv.coef.bugs, nrow(draws))
    colnames(lv.coef.draws) <- names(lv.coef.bugs)
    draws <- cbind(draws, lv.coef.draws)
    fit$data$nlv <- 2
    nLVsim = 2 #calculations faster when not simulating 1000s of LV values, especially since they are all ignored here.
  }
  
  out <- predsumspecies_raw(
    Xocc = datalist$Xocc,
    Xobs = datalist$Xobs,
    ModelSite = datalist$ModelSite,
    numspecies = fit$data$n,
    nlv = fit$data$nlv,
    draws = draws,
    useLVindraws = FALSE,
    nLVsim = nLVsim,
    cl = cl
  )
  return(out)
}

#' @param y A matrix of species *observations* with each row a visit and each column a species. Entries must be either 0 or 1.
#' @param ModelSite The list of ModelSite indexes corresponding to each row in y
#' @return A vector of the number of species detected at each ModelSite. Names give the ModelSite index. 
detectednumspec <- function(y, ModelSite){
  stopifnot(length(ModelSite) == nrow(y))
  stopifnot(all(y %in% c(1, 0)))
  
  my <- cbind(ModelSite = ModelSite, y)
  SpDetected <- my %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise_all(~sum(.) > 0)
  NumSpecies <- as.vector(rowSums(SpDetected[, -1]))
  names(NumSpecies) <- SpDetected[, 1, drop = TRUE]
  return(NumSpecies)
}