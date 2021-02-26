local_edition(3)

test_that("Prep design matrix version 2 works", {
  indata <- artificial_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + Sine1 + UpSite : Sine1 + I(Sine1^2)"
  desmatproc <- prep.designmatprocess_v2(indata, fmla)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  sds <- ((nrow(indata) - 1) / nrow(indata)) * apply(desmat, 2, sd)
  
  # means
  expect_equal(means[c("(Intercept)", "UpSite", "Sine1")], c(1, 0, 0), ignore_attr = TRUE)
  expect_gt(abs(means["UpSite:Sine1"]), 1E-6)
  expect_gt(abs(means["I(Sine1^2)"]), 1E-6)
  
  # sds
  expect_equal(sds[c("(Intercept)", "UpSite", "Sine1")], c(0, 1, 1), ignore_attr = TRUE)
  expect_gt(abs(sds["UpSite:Sine1"]), 1E-6)
  expect_gt(abs(sds["I(Sine1^2)"]), 1E-6)
})

test_that("Prep design matrix version 2 works for logs", {
  indata <- artificial_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + log(UpSite)"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equal(means[c("(Intercept)", "UpSite", "log.UpSite.")], c(1, 0, 0), ignore_attr = TRUE)
})

test_that("Prep design matrix version 2 works for squares", {
  indata <- artificial_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + I(UpSite^2)"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equal(means[c("(Intercept)", "UpSite")], c(1, 0), ignore_attr = TRUE)
  expect_gt(abs(means["I(UpSite^2)"]), 1E-6)
})

test_that("Prep design matrix version 2 doesn't standardise squares with interactions", {
  indata <- artificial_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + I(UpSite^2) * Sine1"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equal(means[c("(Intercept)", "Sine1")], c(1, 0), ignore_attr = TRUE)
  expect_gt(abs(means["I(UpSite^2)"]), 1E-6)
  expect_gt(abs(means["I(UpSite^2):Sine1"]), 1E-6)
})

test_that("Prep design matrix version 2 works for spaces in quotes", {
  indata <- artificial_covar_data(10, 3)[[1]]
  names(indata)[3] <- "Sine 1"
  fmla <- "~ 1 + `Sine 1`"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equal(means[c("(Intercept)", "`Sine 1`")], c(1, 0), ignore_attr = TRUE)
})

test_that("Prep design matrix version 2 works for covariates with I in their name", {
  indata <- artificial_covar_data(10, 3)[[1]]
  names(indata)[3] <- "SIne1"
  fmla <- "~ 1 + SIne1 + Sine2 + I(Sine2^2)"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equal(means[c("(Intercept)", "SIne1", "Sine2")], c(1, 0, 0), ignore_attr = TRUE)
  expect_gt(abs(means["I(Sine2^2)"]), 1E-6)
})

test_that("Prep design matrix version 2 works for intercept only models", {
  indata <- artificial_covar_data(10, 3)[[1]]
  fmla <- "~ 1"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess_v2(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equal(means[c("(Intercept)")], 1, ignore_attr = TRUE)
  expect_equal(sd(desmat[, "(Intercept)"]), 0, ignore_attr = TRUE)
})


test_that("Prep design matrix version 1 works", {
  indata <- artificial_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + Sine1 + UpSite : Sine1 + I(Sine1^2) "
  suppressWarnings(desmatproc <- prep.designmatprocess_v1(indata, fmla))
  desmat <- apply.designmatprocess_v1(desmatproc, indata)
  means <- colMeans(desmat)
  sds <- ((nrow(indata) - 1) / nrow(indata)) * apply(desmat, 2, sd)
  
  # means
  expect_equal(means[c("(Intercept)", "UpSite", "Sine1", "UpSite:Sine1", "I(Sine1^2)")], c(1, 0, 0, 0, 0), ignore_attr = TRUE)

  # sds
  expect_equal(sds[c("(Intercept)", "UpSite", "Sine1", "UpSite:Sine1", "I(Sine1^2)")], c(0, 1, 1, 1, 1), ignore_attr = TRUE)
})

test_that("Selection of correct design matrix processing", {
  indata <- artificial_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + Sine1 + Sine2 + UpSite*Sine2 + I(Sine1^2) + log(UpSite)"
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  expect_equal(desmatproc$version, 2)
  desmat <- apply.designmatprocess(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equal(means[c("(Intercept)", "UpSite", "Sine1", "Sine2", "log.UpSite.")], c(1, 0, 0, 0, 0), ignore_attr = TRUE)
  expect_gt(abs(means["I(Sine1^2)"]), 1E-6)
  expect_gt(abs(means["UpSite:Sine2"]), 1E-6)
  
  suppressWarnings(desmatproc <- prep.designmatprocess(indata, fmla, version = 1))
  expect_equal(desmatproc$version, 1)
  desmat <- apply.designmatprocess(desmatproc, indata)
  means <- colMeans(desmat)
  expect_equal(means, c(1, 0, 0, 0, 0, 0, 0), ignore_attr = TRUE)
})

test_that("Undoing scaling and centering works", {
  indata <- artificial_covar_data(10, 3)[[1]]
  fmla <- "~ 1 + UpSite + Sine1 + Sine2 + UpSite*Sine2 + I(Sine1^2) + log(UpSite)"
  
  suppressWarnings(desmatproc <- prep.designmatprocess(indata, fmla, version = 1))
  desmat <- apply.designmatprocess(desmatproc, indata)
  origmat <- unstandardise.designmatprocess(desmatproc, desmat)
  expect_equal(origmat[, c("UpSite", "Sine1", "Sine2")], indata[, c("UpSite", "Sine1", "Sine2")], ignore_attr = TRUE)
  
  desmatproc <- prep.designmatprocess(indata, fmla, version = 2)
  desmat <- apply.designmatprocess(desmatproc, indata)
  origmat <- unstandardise.designmatprocess(desmatproc, desmat)
  expect_equal(origmat[, c("UpSite", "Sine1", "Sine2")], indata[, c("UpSite", "Sine1", "Sine2")], ignore_attr = TRUE)
})

test_that("Version 2 works inside artificial model building, with prepdata()", {
  artmodel <- artificial_runjags(nspecies = 60, nsites = 100, nvisitspersite = 2, 
                                 OccFmla = "~ 1 + UpSite + Sine1 + Sine2 + UpSite*Sine2 + I(Sine1^2) + log(UpSite)",
                                 ObsFmla = "~ 1 + UpVisit + log(UpVisit) + I(UpVisit^2)",
                                 modeltype = "jsodm_lv",
                                 nlv = 4)
  expect_equal(colMeans(artmodel$data$Xocc[, c("(Intercept)", "UpSite", "Sine1", "Sine2", "log.UpSite.")]), c(1, 0, 0, 0, 0), ignore_attr = TRUE)
  expect_gt(abs(mean(artmodel$data$Xocc[, "I(Sine1^2)"])), 1E-6)
  expect_gt(abs(mean(artmodel$data$Xocc[, "UpSite:Sine2"])), 1E-6)
  
  expect_equal(colMeans(artmodel$data$Xobs)[c("(Intercept)", "UpVisit", "log.UpVisit.")], c(1, 0, 0), ignore_attr = TRUE)
  expect_gt(abs(mean(artmodel$data$Xobs[, "I(UpVisit^2)"])), 1E-6)
})

test_that("Interaction with a binary variable produces zeros and non-zeros", {
  Xocc <- data.frame(ModelSite = 1:10, IsPlanting = 0, PlantingAge = NA)
  Xocc$IsPlanting[6:10] <- 1 #these sites are plantings
  Xocc$PlantingAge[6:10] <- 0.1 + (1:5) #these sites are plantings
  Xocc$PlantingAge[1:5] <- rnorm(5) #these sites are not-plantings, random filler age
  
  # check R's build in methods
  DesMat <- model.matrix(as.formula("~ IsPlanting + IsPlanting:PlantingAge"),
               Xocc)
  expect_equal(DesMat[1:5, "IsPlanting:PlantingAge"], Xocc$IsPlanting[1:5], ignore_attr = TRUE)
  
  # this packages methods
  process <- prep.designmatprocess(Xocc, "~ IsPlanting + IsPlanting:PlantingAge", preserve = c("PlantingAge", "IsPlanting"))
  DesMat2 <- apply.designmatprocess(process, Xocc)
  expect_equal(DesMat2[1:5, "IsPlanting:PlantingAge"], Xocc$IsPlanting[1:5], ignore_attr = TRUE)
  expect_equal(DesMat2[6:10, "IsPlanting:PlantingAge"], Xocc$PlantingAge[6:10], ignore_attr = TRUE)
})


test_that("Covariate scaling scales non-constant covariates correctly", {
  covars <- artificial_covar_data(100, 2)
  
  OccFmla = "~ UpSite + Sine1 + Sine2"
  ObsFmla = "~ UpVisit + Step"
  
  XoccProcess <- prep.designmatprocess(covars$Xocc, OccFmla)
  XobsProcess <- prep.designmatprocess(covars$Xobs, ObsFmla)
  
  XoccDesign <- apply.designmatprocess(XoccProcess, covars$Xocc)
  XobsDesign <- apply.designmatprocess(XobsProcess, covars$Xobs)
  
  expect_equal(colMeans(XoccDesign[, -1]), rep(0, 3), ignore_attr = TRUE)
  expect_equal(colMeans(XobsDesign[, -1]), rep(0, 2), ignore_attr = TRUE)
 
  expect_equal(((100 - 1) / 100) * apply(XoccDesign[, -1], MARGIN = 2, sd), rep(1, 3), ignore_attr = TRUE)
  expect_equal(((200 - 1) / 200) * apply(XobsDesign[, -1], MARGIN = 2, sd), rep(1, 2), ignore_attr = TRUE)
})

