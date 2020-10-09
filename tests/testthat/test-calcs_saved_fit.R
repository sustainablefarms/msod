context("Likelihoods and Predictions Historically Consistent")

fit <- readRDS(system.file("extdata", "cutfit_7_4_11_2LV.rds", package="sustfarmld"))

test_that("Likelihood is historically consistent", {
  # lvsim <- matrix(rnorm(fit$data$nlv * numlvsims), ncol = fit$data$nlv, nrow = numlvsims)
  
  # lkl_site001 <- likelihood_joint_marginal.ModelSite.theta(fit$data$Xocc[1, , drop = FALSE],
  #                                                          fit$data$Xobs[fit$data$ModelSite == 1, , drop = FALSE],
  #                                                          fit$data$y[fit$data$ModelSite == 1, , drop = FALSE],
  #                                                          fit$mcmc[[1]][1, ], 
  #                                                          lvsim = lvsim)
  # expect_known_output(lkl_site001, file = "./tests/testthat/lkl_site001.txt", print = TRUE)
  
  set.seed(333)  #sets seed for simulated LV
  lkl_sites <- likelihoods.fit(fit, numlvsims = 10)
  expect_known_output(lkl_sites, file = "./tests/testthat/lkl_sites.txt", print = TRUE)
})

test_that("Occupancy of species prediction is historically consistent", {
  pocc_theta01 <- poccupy_species(fit, type = 1, conditionalLV = TRUE)
  expect_known_output(pocc_theta01, file = "./tests/testthat/pocc_theta01.txt", print = TRUE)
})
