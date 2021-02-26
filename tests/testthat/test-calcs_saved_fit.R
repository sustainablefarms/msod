local_edition(3)

fit <- readRDS("../../../sflddata/private/data/testdata/cutfit_7_4_11_2LV.rds")
fit <- translatefit(fit)

test_that("Likelihood is historically consistent", {
  set.seed(333)  #sets seed for simulated lv.v
  lkl_sites <- likelihood(fit, numlvsims = 10)

  out <- read.delim("lkl_sites.txt", sep = " ")
  out[1:17, 7:12] <- out[18:34, 1:6]
  savedlkl <- as.matrix(out[1:17, c(2:6, 8:12)])
  expect_equal(lkl_sites, savedlkl, tolerance = 10, ignore_attr = TRUE)
  expect_snapshot_value(lkl_sites, style = "serialize", ignore_attr = TRUE)
  # the magnitude of loglikelihood varies substantially between draw and site, as the simulation methods are different,
  # it is important that at least the magnitude matches
})

test_that("Occupancy of species prediction is historically consistent", {
  pocc_theta01_condlv.v <- poccupy_species(fit, type = 1, conditionalLV = TRUE)
  
  out <- as.matrix(read.delim("pocc_theta01_condLV_part.txt", sep = " ", header = FALSE))
  expect_equal(pocc_theta01_condlv.v[, 1:3], out, ignore_attr = TRUE)
  expect_snapshot_value(pocc_theta01_condlv.v, style = "serialize", ignore_attr = TRUE)

  set.seed(232413)
  pocc_theta01_marglv.v <- poccupy_species(fit, type = 1, conditionalLV = FALSE)
  rownames(pocc_theta01_marglv.v) <- 1:10
  names(dimnames(pocc_theta01_marglv.v)) <- c("", "row")
  expect_snapshot_value(pocc_theta01_marglv.v, style = "serialize", ignore_attr = TRUE)
})

test_that("Detection Probility of Species is historically consistent", { #values generate from code at commit a0812ddad
  pdet_theta01 <- pdetect_condoccupied(fit, type = 1)
  expect_snapshot_value(pdet_theta01, style = "serialize", ignore_attr = TRUE)
})

test_that("Expected Biodiversity is Historically Consistent", { #values generate from code at commit a0812ddad
  pbopt <- pbapply::pboptions(type = "none")
  noccspecies_condlv.v <- speciesrichness(fit,
                                          occORdetection = "occupancy",
                                          usefittedlvv = TRUE)
  expect_snapshot_value(noccspecies_condlv.v, style = "serialize", ignore_attr = TRUE)

  ndetspecies_condlv.v <- speciesrichness(fit,
                                          occORdetection = "detection",
                                          usefittedlvv = TRUE)
  expect_snapshot_value(ndetspecies_condlv.v, style = "serialize", ignore_attr = TRUE)

  set.seed(1341) #for simulated lv.v
  Enoccspecies_marglv.v <- speciesrichness(fit,
                                        occORdetection = "occupancy",
                                        usefittedlvv = FALSE,
                                        nlvperdraw = 1000)
  expect_snapshot_value(Enoccspecies_marglv.v, style = "serialize", ignore_attr = TRUE)
  set.seed(1654)
  Endetspecies_marglv.v <- speciesrichness(fit,
                                        occORdetection = "detection",
                                        usefittedlvv = FALSE,
                                        nlvperdraw = 1000)
  expect_snapshot_value(Endetspecies_marglv.v, style = "serialize", ignore_attr = TRUE)

  pbapply::pboptions(pbopt)
})

test_that("Expected Individual Species Detections are Historically Consistent", { #values generate from code at commit a0812ddad
  En_condlv.v <- Endetect_modelsite(fit, type = 1, conditionalLV = TRUE)$E_ndetect
  expect_snapshot_value(En_condlv.v, style = "serialize", ignore_attr = TRUE)
  
  set.seed(232413)
  En_marglv.v <- Endetect_modelsite(fit, type = 1, conditionalLV = FALSE)$E_ndetect
  expect_snapshot_value(En_marglv.v, style = "serialize", ignore_attr = TRUE)
})
