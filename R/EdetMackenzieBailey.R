#' @title Compare observed detection numbers to expected
#' @description Expected detections compared to observed detections can be used to test model fit
#' (see MacKenzie and Bailey (2004) and Lemeshow, Sturdivant, Hosmer et al (2013)).
#' Note that the latter suggest that the p-value for the test is not very useful for comparing different models.
#' This function follows MacKenzie and Bailey (2004) by using cohorts given by 
#' sites with the same number of survey visits.
#' Optionally, cohorts may be further narrowed to each species.
#' @details For each draw simulates LV twice and simulates detections for each of these once.
#' The expected number of detections is based on the assumpation that the full posterior is correct (and thus averages across posterior draws).
#' The distribution of test statistic according to the full posterior is obtained by simulation of the observations.
#' @param fit a jsodm or jsodm_lv object
#' @param nperdraw The number of simulations per posterior draw (ignored if simulations passed)
#' @param seed The seed for simulations (ignored if simulations passed)
#' @param simulations a set of simulations created by mb_simulate(). If absent it will be created.
#' @param perspecies TRUE for cohorts that are species detections in sites with the same number of visits.
#' FALSE are defined only by the number of visits per site. 
#' @examples 
#' fit <- readRDS("../sflddata/private/data/testdata/cutfit_7_4_11_2LV.rds")
#' fit <- translatefit(fit)
#' X2 <- mb_stat(fit, perspecies = FALSE)
#' hist(X2$sim)
#' abline(v = X2$obs)

#' @export
mb_stat <- function(fit, nperdraw = 2, seed = 122, perspecies = FALSE, simulations = NULL){
  Endet <- Edetperspecies.jsodm(fit)
  Endet_margpost <- apply(Endet, MARGIN = c(1, 2), FUN = mean) #average across draws
  Ondet <- drop_to_matrix(Odetperspecies(fit))

  # cohorts based on number of site visits
  vistnum <- as.data.frame(ftable(fit$data$ModelSite))
  sitecohorts <- split(vistnum$Var1, vistnum$Freq)
  sitecohorts <- lapply(sitecohorts, as.numeric)

  # get value of test statistic
  X2obs <- mb_X2(sitecohorts, Endet = Endet_margpost, Ondet = Ondet, perspecies = perspecies)

  # get distribution of test statistic assuming full posterior is correct
  if (is.null(simulations)){sims <- mb_simulate(fit, nperdraw = nperdraw, seed = seed)
  } else {sims <- simulations}
  Sndet <- lapply(sims, function(y) drop_to_matrix(Odetperspecies(fit, y = y)))
  X2sim_l <- lapply(Sndet, function(simy) mb_X2(sitecohorts, Endet = Endet_margpost, Ondet = simy, perspecies = perspecies)$val)
  X2sim <- unlist(X2sim_l)
  
  # compute probability of X2 being greater than X2obs according to the posterior (simulations)
  distr <- ecdf(X2sim)
  p <- 1 - distr(X2obs$val)
  
  return(list(obs = X2obs$val, obstbl = as.data.frame(X2obs[c("Obs", "E")]), sim = X2sim, p = p))
}


mb_X2 <- function(sitecohorts, Endet, Ondet, perspecies = FALSE){
  if (perspecies){
    ## get E detections in each cohort, per species
    cohort_Edetspec <- lapply(sitecohorts, function(sites) colSums(Endet[sites, , drop = FALSE]))
    ## get O detections in each cohort, per species
    cohort_Ondet <- lapply(sitecohorts, function(sites) colSums(Ondet[sites, , drop = FALSE]))
  } else {
    ## get E detections in each cohort
    cohort_Edetspec <- lapply(sitecohorts, function(sites) sum(Endet[sites, , drop = FALSE]))
    ## get O detections in each cohort
    cohort_Ondet <- lapply(sitecohorts, function(sites) sum(Ondet[sites, , drop = FALSE]))
  }
  ## compute statistic (eq 2.7 of MacKenzie and Bailey 2004)
  cohort_Edetspec_u <- unlist(cohort_Edetspec)
  cohort_Ondet_u <- unlist(cohort_Ondet)
  stopifnot(all(names(cohort_Edetspec_u) == names(cohort_Ondet_u)))
  X2 <- sum( (cohort_Ondet_u - cohort_Edetspec_u)^2 / cohort_Edetspec_u ) 
  return(list(val = X2, Obs = cohort_Ondet_u, E = cohort_Edetspec_u))
}

mb_simulate <- function(fit, nperdraw = 2, seed = 11234){
  # get 2 simulations per draw - 1
  sims <- lapply(1:(length(fit$mcmc) * nrow(fit$mcmc[[1]])),
         function(drawidx){
           out <- lapply(1:nperdraw, function(rep){
             if ("jsodm_lv" %in% class(fit)){
              simulate_detections_lv.v(fit, esttype = drawidx, seed = (seed + nperdraw * 7)*drawidx)
             } else if ("jsodm" %in% class(fit)) {
              simulate_detections(fit, esttype = drawidx, seed = (seed + nperdraw * 7)*drawidx)
             } else {
               NULL
             }
           })
           return(out)
         })
  sims <- unlist(sims, recursive = FALSE)
  return(sims)
}
