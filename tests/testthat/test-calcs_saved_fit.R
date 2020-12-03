context("Likelihoods and Predictions Historically Consistent")

fit <- readRDS("../../../sflddata/private/data/testdata/cutfit_7_4_11_2lv.v.rds")
fit <- translatefit(fit)

test_that("Likelihood is historically consistent", {
  # lvsim <- matrix(rnorm(fit$data$nlv * numlvsims), ncol = fit$data$nlv, nrow = numlvsims)
  
  # lkl_site001 <- likelihood_joint_marginal.ModelSite.theta(fit$data$Xocc[1, , drop = FALSE],
  #                                                          fit$data$Xobs[fit$data$ModelSite == 1, , drop = FALSE],
  #                                                          fit$data$y[fit$data$ModelSite == 1, , drop = FALSE],
  #                                                          fit$mcmc[[1]][1, ], 
  #                                                          lvsim = lvsim)
  # expect_known_output(lkl_site001, file = "./tests/testthat/lkl_site001.txt", print = TRUE)
  
  set.seed(333)  #sets seed for simulated lv.v
  lkl_sites <- likelihoods.fit(fit, numlvsims = 10)
  expect_known_output(lkl_sites, file = "lkl_sites.txt", print = TRUE, update = FALSE)
})

test_that("Occupancy of species prediction is historically consistent", {
  pocc_theta01_condlv.v <- poccupy_species(fit, type = 1, conditionallv.v = TRUE)
  expect_known_output(pocc_theta01_condlv.v, file = "pocc_theta01_condlv.v.txt", print = TRUE, update = FALSE)
  
  set.seed(232413)
  pocc_theta01_marglv.v <- poccupy_species(fit, type = 1, conditionallv.v = FALSE)
  rownames(pocc_theta01_marglv.v) <- 1:10
  names(dimnames(pocc_theta01_marglv.v)) <- c("", "row")
  expect_known_output(pocc_theta01_marglv.v, file = "pocc_theta01_marglv.v.txt", print = TRUE, update = FALSE) #values saved from code with commit a0812ddad
})

test_that("Detection Probility of Species is historically consistent", { #values generate from code at commit a0812ddad
  pdet_theta01 <- pdetect_condoccupied(fit, type = 1)
  expect_known_output(pdet_theta01, file = "pdet_theta01.txt", print = TRUE, update = FALSE)
})

test_that("Expected Biodiversity is Historically Consistent", { #values generate from code at commit a0812ddad
  pbopt <- pbapply::pboptions(type = "none")
  Enspecies_condlv.v <- predsumspecies(fit, UseFittedlv.v = TRUE, type = "marginal")
  expect_known_output(Enspecies_condlv.v, file = "Especrich_condlv.v.txt", print = TRUE, update = FALSE)
  
  set.seed(1341) #for simulated lv.v
  Enspecies_marglv.v <- predsumspecies(fit, UseFittedlv.v = FALSE, type = "marginal")
  expect_known_output(Enspecies_marglv.v, file = "Especrich_marglv.v.txt", print = TRUE, update = FALSE)
  pbapply::pboptions(pbopt)
})

test_that("Expected Individual Species Detections are Historically Consistent", { #values generate from code at commit a0812ddad
  En_condlv.v <- Endetect_modelsite(fit, type = 1, conditionallv.v = TRUE)$E_ndetect
  expect_known_output(En_condlv.v, file = "En_condlv.v.txt", print = TRUE, update = FALSE)
  
  set.seed(232413)
  En_marglv.v <- Endetect_modelsite(fit, type = 1, conditionallv.v = FALSE)$E_ndetect
  expect_known_output(En_marglv.v, file = "En_marglv.v.txt", print = TRUE, update = FALSE)
})
