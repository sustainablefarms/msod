# test likelihood computations

library(testthat)
library(runjags)
runjags.options(silent.jags = TRUE)
runjags.options(silent.runjags = TRUE)
pbopt <- pbapply::pboptions(type = "none")

context("Tests of Likelihood Computations")

test_that("Likelihood computations run in sample data with lv.v", {
  covars <- simulate_covar_data(nsites = 50, nvisitspersite = 2)
  y <- simulate_iid_detections(3, nrow(covars$Xocc))
  
  suppressWarnings(fittedmodel <- run.detectionoccupancy(
    Xocc = covars$Xocc,
    yXobs = cbind(covars$Xobs, y),
    species = colnames(y),
    ModelSite = "ModelSite",
    OccFmla = "~ UpSite + Sine1",
    ObsFmla = "~ UpVisit + Step",
    modeltype = "jsodm_lv",
    nlv = 2,
    MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 2, thin = 1)
  ))
  
  lkl <- likelihoods.fit(fittedmodel)
  expect_equal(dim(lkl), c(2, 50))
})

test_that("Likelihood computations run in sample data with out lv.v", {
  covars <- simulate_covar_data(nsites = 20, nvisitspersite = 2)
  y <- simulate_iid_detections(3, nrow(covars$Xocc))
  
  suppressWarnings(fittedmodel <- run.detectionoccupancy(
    Xocc = covars$Xocc,
    yXobs = cbind(covars$Xobs, y),
    species = colnames(y),
    ModelSite = "ModelSite",
    OccFmla = "~ UpSite + Sine1",
    ObsFmla = "~ UpVisit + Step",
    modeltype = "jsodm",
    MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 2, thin = 1)
  ))
  
  lkl <- likelihoods.fit(fittedmodel)
  expect_equal(dim(lkl), c(2, 20))
})

test_that("lppds insample and outsample data identical when observations identical", {
  artmodel <- artificial_runjags(modeltype = "jsodm")
  
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
  artmodel <- artificial_runjags(modeltype = "jsodm")
  
  lkl <- likelihoods.fit(artmodel)
  lkl <- Rfast::rep_row(lkl, 50)
  waic <- loo::waic(log(lkl))
  
  originalXocc <- unstandardise.designmatprocess(artmodel$XoccProcess, artmodel$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- unstandardise.designmatprocess(artmodel$XobsProcess, artmodel$data$Xobs)
  originalXobs <- cbind(ModelSite = artmodel$data$ModelSite, originalXobs)
  
  outofsample_y <- simulate_detections(artmodel, esttype = 1)
  outofsample_lppd <- lppd.newdata(artmodel,
               Xocc = originalXocc,
               yXobs = cbind(originalXobs, outofsample_y),
               ModelSite = "ModelSite")
  
  # of a randomly selected NEW ModelSite
  lppd_pred_fromoutofsample <- mean(outofsample_lppd$lpds)
  lppd_pred_fromwaic <- waic$estimates["elpd_waic", "Estimate"] / nrow(artmodel$data$Xocc)
  expect_equal(lppd_pred_fromoutofsample, lppd_pred_fromwaic, tolerance = 0.1)
})

test_that("Likelihood computations match simulations without lv.v for nearly certain detection", {
  artfit <- artificial_runjags(nspecies = 2, nsites = 1000, nvisitspersite = 1, 
                                ObsFmla = "~ 1",
                                OccFmla = "~ 1",
                                occ.b.min = -1,  occ.b.max = -0.9, 
                                det.b.min = 20, det.b.max = 20.1, #makes detection almost certain
                               modeltype = "jsodm"
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

test_that("Likelihood computations match simulations without lv.v for multiple visits", {
  artfit <- artificial_runjags(nspecies = 2, nsites = 1000, nvisitspersite = 2, 
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1",
                               occ.b.min = 0,  occ.b.max = 0.001, 
                               det.b.min = 0, det.b.max = 0.01,
                               modeltype = "jsodm")
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

test_that("Likelihood computations match simulations with lv.v, single visits", {
  # the third lv.v is not evenly distributed across the sites
  artfit <- artificial_runjags(nspecies = 2, nsites = 10000, nvisitspersite = 1,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1",
                               occ.b.min = 0,  occ.b.max = 0.001, 
                               det.b.min = 1, det.b.max = 1.001,
                               lv.coef.min = matrix(c(0, 0, 0.45), nrow = 2, ncol = 3, byrow = TRUE),
                               lv.coef.max = matrix(c(0, 0, 0.65), nrow = 2, ncol = 3, byrow = TRUE),
                               modeltype = "jsodm_lv",
                               nlv = 3)
                               
  # resimulate y as if lv.v are not known (which is the situation for likelihood compuations)
  artfit$data$y <- simulate_detections_lv.v(artfit, esttype = 1)
  
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

test_that("Likelihood computations match simulations with lv.v, multiple visits", {
  artfit <- artificial_runjags(nspecies = 2, nsites = 10000, nvisitspersite = 2,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1",
                               occ.b.min = 0,  occ.b.max = 0.001, 
                               det.b.min = 1, det.b.max = 1.001,
                               modeltype = "jsodm_lv",
                               nlv = 2)
  # resimulate y as if lv.v are not known (which is the situation for likelihood compuations)
  artfit$data$y <- simulate_detections(artfit, esttype = 1)
  
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