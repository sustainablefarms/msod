# test likelihood computations

library(testthat)
library(runjags)
runjags.options(silent.jags = TRUE)
runjags.options(silent.runjags = TRUE)
pbopt <- pbapply::pboptions(type = "none")

context("Tests of Likelihood Computations")

test_that("Likelihood computations run in sample data with LV", {
  covars <- simulate_covar_data(nsites = 50, nvisitspersite = 2)
  y <- simulate_iid_detections(3, nrow(covars$Xocc))
  
  suppressWarnings(fittedmodel <- run.detectionoccupancy(
    Xocc = covars$Xocc,
    yXobs = cbind(covars$Xobs, y),
    species = colnames(y),
    ModelSite = "ModelSite",
    OccFmla = "~ UpSite + Sine1",
    ObsFmla = "~ UpVisit + Step",
    nlv = 2,
    MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 2, thin = 1)
  ))
  
  lkl <- likelihoods.fit(fittedmodel)
  expect_equal(dim(lkl), c(2, 50))
})

test_that("Likelihood computations run in sample data with out LV", {
  covars <- simulate_covar_data(nsites = 20, nvisitspersite = 2)
  y <- simulate_iid_detections(3, nrow(covars$Xocc))
  
  suppressWarnings(fittedmodel <- run.detectionoccupancy(
    Xocc = covars$Xocc,
    yXobs = cbind(covars$Xobs, y),
    species = colnames(y),
    ModelSite = "ModelSite",
    OccFmla = "~ UpSite + Sine1",
    ObsFmla = "~ UpVisit + Step",
    nlv = 0,
    MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 2, thin = 1)
  ))
  
  lkl <- likelihoods.fit(fittedmodel)
  expect_equal(dim(lkl), c(2, 20))
})

test_that("lppds insample and outsample data identical when observations identical", {
  artmodel <- artificial_runjags(nlv = 0)
  
  lkl <- likelihoods.fit(artmodel)
  lkl <- Rfast::rep_row(lkl, 50)
  waic <- loo::waic(log(lkl))
  
  originalXocc <- unstandardise.designmatprocess(artmodel$XoccProcess, artmodel$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- unstandardise.designmatprocess(artmodel$XobsProcess, artmodel$data$Xobs)
  originalXobs <- cbind(ModelSite = artmodel$data$ModelSite, originalXobs)
  
  outofsample_y <- artmodel$data$y
  outofsample_lppd <- lppd.newdata(artmodel,
                                   Xocc = originalXocc,
                                   yXobs = cbind(originalXobs, outofsample_y),
                                   ModelSite = "ModelSite")
  
  # of a randomly selected NEW ModelSite
  lppd_pred_fromoutofsample <- mean(outofsample_lppd$lpds)
  lppd_pred_fromwaic <- waic$estimates["elpd_waic", "Estimate"] / nrow(artmodel$data$Xocc)
  expect_equal(lppd_pred_fromoutofsample, lppd_pred_fromwaic)
})

test_that("lppds insample and outsample data similar on very artifical simple situation", {
  artmodel <- artificial_runjags(nlv = 0)
  
  lkl <- likelihoods.fit(artmodel)
  lkl <- Rfast::rep_row(lkl, 50)
  waic <- loo::waic(log(lkl))
  
  originalXocc <- unstandardise.designmatprocess(artmodel$XoccProcess, artmodel$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- unstandardise.designmatprocess(artmodel$XobsProcess, artmodel$data$Xobs)
  originalXobs <- cbind(ModelSite = artmodel$data$ModelSite, originalXobs)
  
  outofsample_y <- simulate_fit(artmodel, esttype = 1, UseFittedLV = FALSE)
  outofsample_lppd <- lppd.newdata(artmodel,
               Xocc = originalXocc,
               yXobs = cbind(originalXobs, outofsample_y),
               ModelSite = "ModelSite")
  
  # of a randomly selected NEW ModelSite
  lppd_pred_fromoutofsample <- mean(outofsample_lppd$lpds)
  lppd_pred_fromwaic <- waic$estimates["elpd_waic", "Estimate"] / nrow(artmodel$data$Xocc)
  expect_equal(lppd_pred_fromoutofsample, lppd_pred_fromwaic, tolerance = 0.1)
})

test_that("Likelihood computations match simulations without LV for nearly certain detection", {
  artfit <- artificial_runjags(nspecies = 2, nsites = 1000, nvisitspersite = 1, nlv = 0,
                                ObsFmla = "~ 1",
                                OccFmla = "~ 1",
                                u.b.min = -1,  u.b.max = -0.9, 
                                v.b.min = 20, v.b.max = 20.1, #makes detection almost certain
                                )
  jointoutcomes <- apply(artfit$data$y, 1, paste0, collapse = ",")
  
  # theory likelihoods
  poccupancy_all <- poccupy_species(artfit, type = 1, conditionalLV = FALSE)[1, ]
  rvA <- discreteRV::RV(c(1, 0), poccupancy_all["A"], fractions = FALSE)
  rvB <- discreteRV::RV(c(1, 0), poccupancy_all["B"], fractions = FALSE)
  jointRV <- discreteRV::joint(rvA, rvB)
  lkl_th <- vapply(jointoutcomes, function(x) discreteRV::P(jointRV == x), FUN.VALUE = 3.3) 
  
  # by simulation
  sim_distr <- summary(factor(jointoutcomes))/1000
  lkl_sim <- sim_distr[jointoutcomes]
  
  # from this package
  lkl <- likelihoods.fit(artfit)
  
  # compare all
  expect_equivalent(lkl_th, lkl)
  expect_equivalent(lkl_th, lkl_sim, tolerance = 0.05)
})

test_that("Likelihood computations match simulations without LV for multiple visits", {
  artfit <- artificial_runjags(nspecies = 2, nsites = 1000, nvisitspersite = 2, nlv = 0,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1",
                               u.b.min = 0,  u.b.max = 0.001, 
                               v.b.min = 0, v.b.max = 0.01)
  my <- cbind(ModelSite = artfit$data$ModelSite, artfit$data$y)
  obs_per_site <- lapply(1:nrow(artfit$data$Xocc), function(x) my[my[, "ModelSite"] == x, -1])
  jointoutcomes <- vapply(obs_per_site, paste0, collapse = ",", FUN.VALUE = "achar")
  
  # likelihood by simulation
  sim_distr <- summary(factor(jointoutcomes))/length(jointoutcomes)
  lkl_sim <- sim_distr[jointoutcomes]
  
  # from this package
  lkl <- likelihoods.fit(artfit)
  
  # compare 
  expect_equivalent(lkl, lkl_sim, tolerance = 0.05)
})

test_that("Likelihood computations match simulations with LV, single visits", {
  # the third LV is not evenly distributed across the sites
  artfit <- artificial_runjags(nspecies = 2, nsites = 10000, nvisitspersite = 1, nlv = 3,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1",
                               u.b.min = 0,  u.b.max = 0.001, 
                               v.b.min = 1, v.b.max = 1.001,
                               lv.coef.min = matrix(c(0, 0, 0.45), nrow = 2, ncol = 3, byrow = TRUE),
                               lv.coef.max = matrix(c(0, 0, 0.65), nrow = 2, ncol = 3, byrow = TRUE))
                               
  # resimulate y as if LV are not known (which is the situation for likelihood compuations)
  artfit$data$y <- simulate_fit(artfit, esttype = 1, UseFittedLV = FALSE)
  
  # get joint outcomes
  my <- cbind(ModelSite = artfit$data$ModelSite, artfit$data$y)
  obs_per_site <- lapply(1:nrow(artfit$data$Xocc), function(x) my[my[, "ModelSite"] == x, -1])
  jointoutcomes <- vapply(obs_per_site, paste0, collapse = ",", FUN.VALUE = "achar")
  
  # likelihood by simulation
  sim_distr <- summary(factor(jointoutcomes))/length(jointoutcomes)
  lkl_sim <- sim_distr[jointoutcomes]
  
  # from this package
  lkl <- likelihoods.fit(artfit)
  
  # compare 
  expect_equivalent(lkl, lkl_sim, tolerance = 0.01)
})

test_that("Likelihood computations match simulations with LV, multiple visits", {
  artfit <- artificial_runjags(nspecies = 2, nsites = 10000, nvisitspersite = 2, nlv = 2,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1",
                               u.b.min = 0,  u.b.max = 0.001, 
                               v.b.min = 1, v.b.max = 1.001)
  # resimulate y as if LV are not known (which is the situation for likelihood compuations)
  artfit$data$y <- simulate_fit(artfit, esttype = 1, UseFittedLV = TRUE)
  
  # get joint outcomes
  my <- cbind(ModelSite = artfit$data$ModelSite, artfit$data$y)
  obs_per_site <- lapply(1:nrow(artfit$data$Xocc), function(x) my[my[, "ModelSite"] == x, -1])
  jointoutcomes <- vapply(obs_per_site, paste0, collapse = ",", FUN.VALUE = "achar")
  
  # likelihood by simulation
  sim_distr <- summary(factor(jointoutcomes))/length(jointoutcomes)
  lkl_sim <- sim_distr[jointoutcomes]
  
  # from this package
  lkl <- likelihoods.fit(artfit)
  
  # compare 
  expect_equivalent(lkl, lkl_sim, tolerance = 0.01)
})


pbapply::pboptions(pbopt)