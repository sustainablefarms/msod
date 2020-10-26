context("Likelihoods and Predictions Historically Consistent")

fit <- readRDS("../../../sflddata/private/data/testdata/cutfit_7_4_11_2LV.rds")

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
  expect_known_output(lkl_sites, file = "lkl_sites.txt", print = TRUE, update = FALSE)
})

test_that("Occupancy of species prediction is historically consistent", {
  pocc_theta01_condLV <- poccupy_species(fit, type = 1, conditionalLV = TRUE)
  expect_known_output(pocc_theta01_condLV, file = "pocc_theta01_condLV.txt", print = TRUE, update = FALSE)
  
  set.seed(232413)
  pocc_theta01_margLV <- poccupy_species(fit, type = 1, conditionalLV = FALSE)
  rownames(pocc_theta01_margLV) <- 1:10
  names(dimnames(pocc_theta01_margLV)) <- c("", "row")
  expect_known_output(pocc_theta01_margLV, file = "pocc_theta01_margLV.txt", print = TRUE, update = FALSE) #values saved from code with commit a0812ddad
})

test_that("Detection Probility of Species is historically consistent", { #values generate from code at commit a0812ddad
  pdet_theta01 <- pdetect_condoccupied(fit, type = 1)
  expect_known_output(pdet_theta01, file = "pdet_theta01.txt", print = TRUE, update = FALSE)
})

test_that("Expected Biodiversity is Historically Consistent", { #values generate from code at commit a0812ddad
  pbopt <- pbapply::pboptions(type = "none")
  Enspecies_condLV <- predsumspecies(fit, UseFittedLV = TRUE, type = "marginal")
  expect_known_output(Enspecies_condLV, file = "Especrich_condLV.txt", print = TRUE, update = FALSE)
  
  set.seed(1341) #for simulated LV values
  Enspecies_margLV <- predsumspecies(fit, UseFittedLV = FALSE, type = "marginal")
  expect_known_output(Enspecies_margLV, file = "Especrich_margLV.txt", print = TRUE, update = FALSE)
  pbapply::pboptions(pbopt)
})

test_that("Expected Individual Species Detections are Historically Consistent", { #values generate from code at commit a0812ddad
  En_condLV <- Endetect_modelsite(fit, type = 1, conditionalLV = TRUE)$E_ndetect
  expect_known_output(En_condLV, file = "En_condLV.txt", print = TRUE, update = FALSE)
  
  set.seed(232413)
  En_margLV <- Endetect_modelsite(fit, type = 1, conditionalLV = FALSE)$E_ndetect
  expect_known_output(En_margLV, file = "En_margLV.txt", print = TRUE, update = FALSE)
})
