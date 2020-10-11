# test that the following are consistent: 
# simulation data, runjags MCMC (posterior distribution matches true), predicted likelihoods, expected species detections, DS residuals

context("Wholistic tests using identical, independent, ModelSites")
skip_if(parallel::detectCores() < 10)

# Create a process with known parameters
artmodel <- artificial_runjags(nspecies = 5, nsites = 2000, nvisitspersite = 2, nlv = 0,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1")

# fit to data and simulations using runjags
originalXocc <- Rfast::eachrow(Rfast::eachrow(artmodel$data$Xocc, artmodel$XoccProcess$scale, oper = "*"),
                               artmodel$XoccProcess$center, oper = "+")
colnames(originalXocc) <- colnames(artmodel$data$Xocc)
originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
originalXobs <- Rfast::eachrow(Rfast::eachrow(artmodel$data$Xobs, artmodel$XobsProcess$scale, oper = "*"),
                               artmodel$XobsProcess$center, oper = "+")
colnames(originalXobs) <- colnames(artmodel$data$Xobs)
originalXobs <- cbind(ModelSite = artmodel$data$ModelSite, originalXobs)

fit_runjags <- run.detectionoccupancy(originalXocc, cbind(originalXobs, artmodel$data$y), 
                       species = colnames(artmodel$data$y),
                       ModelSite = "ModelSite",
                       OccFmla = artmodel$XoccProcess$fmla,
                       ObsFmla = artmodel$XobsProcess$fmla,
                       initsfunction = function(chain, indata){return(NULL)},
                       MCMCparams = list(n.chains = 2, adapt = 1000, burnin = 10000, sample = 500, thin = 40),
                       nlv = 0)

cl <- parallel::makeCluster(10)
lkl_runjags <- likelihoods.fit(fit_runjags, cl = cl)
lkl_artmodel <- likelihoods.fit(artmodel, cl = cl)
Enumspec <- predsumspecies(fit_runjags, UseFittedLV = FALSE, type = "marginal", cl = cl)
parallel::stopCluster(cl)

save(fit_runjags, artmodel, originalXocc, originalXobs,
     lkl_runjags, lkl_artmodel,  Enumspec, file = "../../tests/testthat/benchmark_identicalsitesmodel.Rdata")

load("benchmark_identicalsitesmodel.Rdata")

test_that("Posterior credible distribution overlaps true parameters", {
  var2compare <- colnames(artmodel$mcmc[[1]])
  inCI <- (fit_runjags$summaries[var2compare, "Lower95"] <= artmodel$mcmc[[1]][1, var2compare]) &
    (fit_runjags$summaries[var2compare, "Upper95"] >= artmodel$mcmc[[1]][1, var2compare])
  expect_equal(mean(inCI), 0.95, tol = 0.051)
})



test_that("Predicted likelihoods match observations", {
  my <- cbind(ModelSite = fit_runjags$data$ModelSite, fit_runjags$data$y)
  obs_per_site <- lapply(1:nrow(fit_runjags$data$Xocc), function(x) my[my[, "ModelSite"] == x, -1])
  jointoutcomes <- vapply(obs_per_site, paste0, collapse = ",", FUN.VALUE = "achar")
  
  # extra, fake observations for more accurate simulated likelihood
  makemoreobs <- function(){
    my_add <- cbind(ModelSite = fit_runjags$data$ModelSite, simulate_fit(artmodel, esttype = 1, UseFittedLV = FALSE))
    obs_per_site_add <- lapply(1:nrow(fit_runjags$data$Xocc), function(x) my_add[my_add[, "ModelSite"] == x, -1])
    jointoutcomes_add <- vapply(obs_per_site_add, paste0, collapse = ",", FUN.VALUE = "achar")
    return(jointoutcomes_add)
  }
  cl <- parallel::makeCluster(10)
  parallel::clusterExport(cl, c("makemoreobs", "fit_runjags", "artmodel"))
  parallel::clusterEvalQ(cl, library(sustfarmld))
  jointoutcomes_more <- pbapply::pbreplicate(10000, makemoreobs(), simplify = FALSE, cl = cl)
  
  # likelihood by simulation
  jointoutcomes_all <- c(jointoutcomes, unlist(jointoutcomes_more))
  sim_distr <- table(factor(jointoutcomes_all))/length(jointoutcomes_all)
  sim_distr_v <- as.vector(sim_distr)
  names(sim_distr_v) <- names(sim_distr)
  lkl_sim <- sim_distr_v[jointoutcomes]
  parallel::stopCluster(cl)
  
  # sim vs artmodel
  expect_equivalent(Rfast::colmeans(lkl_artmodel), lkl_sim, tolerance = 0.01)
  reldiff_art_sim <- abs(Rfast::colmeans(lkl_artmodel) - lkl_sim) / lkl_sim
  expect_lt(quantile(reldiff_art_sim, probs = 0.9), 0.1)
  
  # sim vs runjags
  expect_equivalent(Rfast::colmeans(lkl_runjags), lkl_sim, tolerance = 0.01)
  reldiff_jags_sim <- abs(Rfast::colmeans(lkl_runjags) - lkl_sim) / lkl_sim
  expect_lt(quantile(reldiff_art_sim, probs = 0.9), 0.1)
  
   # runjags vs artmodel
  expect_equivalent(Rfast::colmeans(lkl_runjags), Rfast::colmeans(lkl_artmodel), tolerance = 0.01)
  reldiff_jags_art <- abs(Rfast::colmeans(lkl_runjags) - Rfast::colmeans(lkl_artmodel)) / Rfast::colmeans(lkl_artmodel)
  expect_lt(quantile(reldiff_jags_art, probs = 0.9), 0.1)
  
  hist(reldiff_jags_sim)
  hist(reldiff_art_sim)
  plt <- cbind(simlkl = lkl_sim,
               lklrunjags = Rfast::colmeans(lkl_runjags),
               lklartmodel = Rfast::colmeans(lkl_artmodel)) %>%
    dplyr::as_tibble(rownames = "Observations") %>%
    dplyr::mutate(ModelSiteID = 1:nrow(fit_runjags$data$Xocc)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = ModelSiteID, y = lklrunjags)) +
    ggplot2::geom_point(ggplot2::aes(x = ModelSiteID, y = lklartmodel), col = "red", shape = "+") +
    ggplot2::geom_point(ggplot2::aes(x = ModelSiteID, y = simlkl), col = "blue", shape = "+") +
    # ggplot2::geom_point(ggplot2::aes(x = ModelSiteID, y = (lkl - simlkl) / simlkl), col = "blue", shape = "+") +
    ggplot2::scale_y_continuous(trans = "log10")
  print(plt)
})

test_that("Expected Number of Detected Species", {
  NumSpecies <- detectednumspec(y = fit_runjags$data$y, ModelSite = fit_runjags$data$ModelSite)
  
  meandiff <- dplyr::cummean(NumSpecies - Enumspec["Esum_det", ])
  meanvar <- cumsum(Enumspec["Vsum_det", ])/((1:ncol(Enumspec))^2)
  plt <- cbind(diff = meandiff, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  # check with predicted standard error once the software is computed
  sd_final <- sqrt(meanvar[ncol(Enumspec)])
  expect_equal(meandiff[ncol(Enumspec)], 0, tol = 3 * sd_final)
  
  Enum_compare_sum <- Enum_compare(NumSpecies,
                                   as.matrix(Enumspec["Esum_det", ], ncol = 1),
                                   as.matrix(Enumspec["Vsum_det", ], ncol = 1)
  )
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]])
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]])
  expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
})

