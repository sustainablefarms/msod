# test running detection occupancy mcmc

library(testthat)

context("Tests of Running JAGS Preparation")

test_that("Covariate scaling scales non-constant covariates correctly", {
  covars <- artificial_covar_data(100, 2)
  
  OccFmla = "~ UpSite + Sine1 + Sine2"
  ObsFmla = "~ UpVisit + Step"
  
  XoccProcess <- prep.designmatprocess(covars$Xocc, OccFmla)
  XobsProcess <- prep.designmatprocess(covars$Xobs, ObsFmla)
  
  XoccDesign <- apply.designmatprocess(XoccProcess, covars$Xocc)
  XobsDesign <- apply.designmatprocess(XobsProcess, covars$Xobs)
  
  expect_equivalent(colMeans(XoccDesign[, -1]), rep(0, 3))
  expect_equivalent(colMeans(XobsDesign[, -1]), rep(0, 2))
 
  expect_equivalent(((100 - 1) / 100) * apply(XoccDesign[, -1], MARGIN = 2, sd), rep(1, 3))
  expect_equivalent(((200 - 1) / 200) * apply(XobsDesign[, -1], MARGIN = 2, sd), rep(1, 2))
})

