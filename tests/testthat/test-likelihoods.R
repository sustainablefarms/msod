# test likelihood computations

library(testthat)

test_that("Likelihood computations run in sample data with LV", {
  covars <- simulate_covar_data(nsites = 50, nvisitspersite = 2)
  y <- simulate_iid_detections(3, nrow(covars$Xocc))
  
  fittedmodel <- run.detectionoccupany(
    Xocc = covars$Xocc,
    yXobs = cbind(covars$Xobs, y),
    species = colnames(y),
    ModelSite = "ModelSite",
    OccFmla = "~ UpSite + Sine1",
    ObsFmla = "~ UpVisit + Step",
    nlv = 2,
    MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 2, thin = 1)
  )
  
  lkl <- likelihoods.fit(fittedmodel)
  expect_equal(dim(lkl), c(2, 50))
})

test_that("Likelihood computations run in sample data with out LV", {
  covars <- simulate_covar_data(nsites = 20, nvisitspersite = 2)
  y <- simulate_iid_detections(3, nrow(covars$Xocc))
  
  fittedmodel <- run.detectionoccupany(
    Xocc = covars$Xocc,
    yXobs = cbind(covars$Xobs, y),
    species = colnames(y),
    ModelSite = "ModelSite",
    OccFmla = "~ UpSite + Sine1",
    ObsFmla = "~ UpVisit + Step",
    nlv = 0,
    MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 2, thin = 1)
  )
  
  lkl <- likelihoods.fit(fittedmodel)
  expect_equal(dim(lkl), c(2, 20))
})