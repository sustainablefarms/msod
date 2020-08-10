
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

test_that("Prep design matrix version 1 works", {
  indata <- simulate_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + Sine1 + UpSite : Sine1 + I(Sine1^2) "
  desmatproc <- prep.designmatprocess_v1(indata, fmla)
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
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equivalent(means[c("(Intercept)", "UpSite", "Sine1", "Sine2", "log.UpSite.")], c(1, 0, 0, 0, 0))
  expect_gt(abs(means["I(Sine1^2)"]), 1E-6)
  expect_gt(abs(means["UpSite:Sine2"]), 1E-6)
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 1)
  expect_equal(desmatproc$version, 1)
  desmat <- apply.designmatprocess_v1(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equivalent(means, c(1, 0, 0, 0, 0, 0, 0))
})