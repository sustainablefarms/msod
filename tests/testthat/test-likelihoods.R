# test likelihood computations

library(testthat)
runjags.options(silent.jags = TRUE)
runjags.options(silent.runjags = TRUE)

test_that("Likelihood computations run in sample data with LV", {
  covars <- simulate_covar_data(nsites = 50, nvisitspersite = 2)
  y <- simulate_iid_detections(3, nrow(covars$Xocc))
  
  suppressWarnings(fittedmodel <- run.detectionoccupany(
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
  
  suppressWarnings(fittedmodel <- run.detectionoccupany(
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

test_that("lppds insample and outsample data similar on very artifical simple situation", {
  artmodel <- artificial_runjags(nlv = 0)
  
  lkl <- likelihoods.fit(artmodel)
  lkl <- Rfast::rep_row(lkl, 50)
  waic <- loo::waic(log(lkl))
  
  
  outofsample_y <- simulate_fit(artmodel, esttype = 1, conditionalLV = FALSE)
  outofsample_lppd <- lppd.newdata(artmodel,
               Xocc = cbind(ModelSite = 1:nrow(artmodel$data$Xocc), artmodel$data$Xocc),
               yXobs = cbind(ModelSite = artmodel$data$ModelSite, artmodel$data$Xobs, outofsample_y),
               ModelSite = "ModelSite")
  
  # of a randomly selected NEW ModelSite
  lppd_pred_fromoutofsample <- mean(outofsample_lppd$lpds)
  lppd_pred_fromwaic <- waic$estimates["elpd_waic", "Estimate"] / nrow(artmodel$data$Xocc)
  expect_equal(lppd_pred_fromoutofsample, lppd_pred_fromwaic, tolerance = 0.5)
})

