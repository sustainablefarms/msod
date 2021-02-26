# test data preparation and transformation
testthat::local_edition(3)

test_that("JAGS preparations preserve the order of input data", {
  datalist <- artificial_covar_data(nsites = 20, nvisitspersite = 2)  
  datalist$Xocc$StepA <- 1
  datalist$Xocc$StepA[1:10] <- 0
  datalist$y <- data.frame(ModelSite = datalist$Xobs$ModelSite, SpeciesA = runif(nrow(datalist$Xobs)))
  datalist$y$SpeciesA[datalist$y$ModelSite <= 10] <- (datalist$y$SpeciesA[datalist$y$ModelSite <= 10] < 0.3)*1
  datalist$y$SpeciesA[datalist$y$ModelSite > 10] <- (datalist$y$SpeciesA[datalist$y$ModelSite > 10] < 0.75)*1
  # SpeciesA has about 0.5 chance of being detected at a modelsite with StepA == 0. 
  # It has a 0.93 chance of being detected at a site with StepA == 1
  
  cbnd <- datalist$y %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise(SpeciesAdetected = max(SpeciesA)) %>%
    dplyr::left_join(datalist$Xocc, by = "ModelSite")
  ftbl1 <- ftable(SpeciesAdetected = cbnd$SpeciesAdetected, StepA = cbnd$StepA)
  
  jagsdatalist <- prepJAGSdata2("jsodm", as.matrix(datalist$Xocc), as.matrix(datalist$Xobs), as.matrix(datalist$y[, -1, drop = FALSE]), as.vector(datalist$Xobs$ModelSite))
  expect_equal(jagsdatalist$y, datalist$y[, -1], ignore_attr = TRUE)
  expect_equal(jagsdatalist$Xobs, as.matrix(datalist$Xobs), ignore_attr = TRUE)
  expect_equal(jagsdatalist$Xocc, as.matrix(datalist$Xocc), ignore_attr = TRUE)
  
  runjags::runjags.options(silent.jags=TRUE, silent.runjags=TRUE)
  fit <- suppressWarnings(fitjsodm2(
    Xocc = as.matrix(datalist$Xocc),
    Xobs = as.matrix(datalist$Xobs),
    y = as.matrix(datalist$y[, -1, drop = FALSE]),
    ModelSite = as.integer(datalist$Xobs$ModelSite),
    modeltype = "jsodm",
    MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 1, thin = 1)
  ))
  
  expect_equal(fit$data$ModelSite, datalist$Xobs$ModelSite, ignore_attr = TRUE)
  expect_equal(fit$data$y, as.matrix(datalist$y[, -1]), ignore_attr = TRUE)
  expect_equal(fit$data$Xocc, as.matrix(datalist$Xocc), ignore_attr = TRUE)
  
  cbnd <- data.frame(cbind(ModelSite = fit$data$ModelSite, fit$data$y))  %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise(SpeciesAdetected = max(SpeciesA)) %>%
    dplyr::left_join(data.frame(fit$data$Xocc), by = "ModelSite")
  ftbl2 <- ftable(SpeciesAdetected = cbnd$SpeciesAdetected, StepA = cbnd$StepA)
  
  expect_equal(ftbl2, ftbl1, ignore_attr = TRUE)
})
