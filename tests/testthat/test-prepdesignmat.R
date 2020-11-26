
context("Prep design matrix functions")

test_that("Prep design matrix version 2 works", {
  indata <- simulate_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + Sine1 + UpSite : Sine1 + I(Sine1^2)"
  desmatproc <- prep.designmatprocess_v2(indata, fmla)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  sds <- ((nrow(indata) - 1) / nrow(indata)) * apply(desmat, 2, sd)
  
  # means
  expect_equivalent(means[c("(Intercept)", "UpSite", "Sine1")], c(1, 0, 0))
  expect_gt(abs(means["UpSite:Sine1"]), 1E-6)
  expect_gt(abs(means["I(Sine1^2)"]), 1E-6)
  
  # sds
  expect_equivalent(sds[c("(Intercept)", "UpSite", "Sine1")], c(0, 1, 1))
  expect_gt(abs(sds["UpSite:Sine1"]), 1E-6)
  expect_gt(abs(sds["I(Sine1^2)"]), 1E-6)
})

test_that("Prep design matrix version 2 works for logs", {
  indata <- simulate_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + log(UpSite)"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equivalent(means[c("(Intercept)", "UpSite", "log.UpSite.")], c(1, 0, 0))
})

test_that("Prep design matrix version 2 works for squares", {
  indata <- simulate_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + I(UpSite^2)"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equivalent(means[c("(Intercept)", "UpSite")], c(1, 0))
  expect_gt(abs(means["I(UpSite^2)"]), 1E-6)
})

test_that("Prep design matrix version 2 doesn't standardise squares with interactions", {
  indata <- simulate_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + I(UpSite^2) * Sine1"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equivalent(means[c("(Intercept)", "Sine1")], c(1, 0))
  expect_gt(abs(means["I(UpSite^2)"]), 1E-6)
  expect_gt(abs(means["I(UpSite^2):Sine1"]), 1E-6)
})

test_that("Prep design matrix version 2 works for spaces in quotes", {
  indata <- simulate_covar_data(10, 3)[[1]]
  names(indata)[3] <- "Sine 1"
  fmla <- "~ 1 + `Sine 1`"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equivalent(means[c("(Intercept)", "`Sine 1`")], c(1, 0))
})

test_that("Prep design matrix version 2 works for covariates with I in their name", {
  indata <- simulate_covar_data(10, 3)[[1]]
  names(indata)[3] <- "SIne1"
  fmla <- "~ 1 + SIne1 + Sine2 + I(Sine2^2)"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equivalent(means[c("(Intercept)", "SIne1", "Sine2")], c(1, 0, 0))
  expect_gt(abs(means["I(Sine2^2)"]), 1E-6)
})

test_that("Prep design matrix version 2 works for intercept only models", {
  indata <- simulate_covar_data(10, 3)[[1]]
  fmla <- "~ 1"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equivalent(means[c("(Intercept)")], 1)
  expect_equivalent(sd(desmat[, "(Intercept)"]), 0)
})


test_that("Prep design matrix version 1 works", {
  indata <- simulate_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + Sine1 + UpSite : Sine1 + I(Sine1^2) "
  suppressWarnings(desmatproc <- prep.designmatprocess_v1(indata, fmla))
  desmat <- apply.designmatprocess_v1(desmatproc, indata)
  means <- colMeans(desmat)
  sds <- ((nrow(indata) - 1) / nrow(indata)) * apply(desmat, 2, sd)
  
  # means
  expect_equivalent(means[c("(Intercept)", "UpSite", "Sine1", "UpSite:Sine1", "I(Sine1^2)")], c(1, 0, 0, 0, 0))

  # sds
  expect_equivalent(sds[c("(Intercept)", "UpSite", "Sine1", "UpSite:Sine1", "I(Sine1^2)")], c(0, 1, 1, 1, 1))
})

test_that("Selection of correct design matrix processing", {
  indata <- simulate_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + Sine1 + Sine2 + UpSite*Sine2 + I(Sine1^2) + log(UpSite)"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equivalent(means[c("(Intercept)", "UpSite", "Sine1", "Sine2", "log.UpSite.")], c(1, 0, 0, 0, 0))
  expect_gt(abs(means["I(Sine1^2)"]), 1E-6)
  expect_gt(abs(means["UpSite:Sine2"]), 1E-6)
  
  suppressWarnings(desmatproc <- prep.designmatprocess(indata, fmla, version = 1))
  expect_equal(desmatproc$version, 1)
  desmat <- apply.designmatprocess(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equivalent(means, c(1, 0, 0, 0, 0, 0, 0))
})

test_that("Undoing scaling and centering works", {
  indata <- simulate_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + Sine1 + Sine2 + UpSite*Sine2 + I(Sine1^2) + log(UpSite)"
  
  suppressWarnings(desmatproc <- prep.designmatprocess(indata, fmla, version = 1))
  desmat <- apply.designmatprocess(desmatproc, indata)
  origmat <- unstandardise.designmatprocess(desmatproc, desmat)
  expect_equivalent(origmat[, c("UpSite", "Sine1", "Sine2")], indata[, c("UpSite", "Sine1", "Sine2")])
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  desmat <- apply.designmatprocess(desmatproc, indata)
  origmat <- unstandardise.designmatprocess(desmatproc, desmat)
  expect_equivalent(origmat[, c("UpSite", "Sine1", "Sine2")], indata[, c("UpSite", "Sine1", "Sine2")])
})

test_that("Version 2 works inside artificial model building, with prepdata()", {
  artmodel <- artificial_runjags(nspecies = 60, nsites = 100, nvisitspersite = 2, 
                                 OccFmla = "~ 1 + UpSite + Sine1 + Sine2 + UpSite*Sine2 + I(Sine1^2) + log(UpSite)",
                                 ObsFmla = "~ 1 + UpVisit + log(UpVisit) + I(UpVisit^2)",
                                 modeltype = "jsodm_lv",
                                 nlv = 4)
  expect_equivalent(colMeans(artmodel$data$Xocc[, c("(Intercept)", "UpSite", "Sine1", "Sine2", "log.UpSite.")]), c(1, 0, 0, 0, 0))
  expect_gt(abs(mean(artmodel$data$Xocc[, "I(Sine1^2)"])), 1E-6)
  expect_gt(abs(mean(artmodel$data$Xocc[, "UpSite:Sine2"])), 1E-6)
  
  expect_equivalent(colMeans(artmodel$data$Xobs)[c("(Intercept)", "UpVisit", "log.UpVisit.")], c(1, 0, 0))
  expect_gt(abs(mean(artmodel$data$Xobs[, "I(UpVisit^2)"])), 1E-6)
})
