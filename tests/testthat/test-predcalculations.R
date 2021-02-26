library(testthat);

context("Calculations of Predicted Probabilities")

test_that("poccupy marginal across LV by simulation is equal to theoretical", {
  fit <- artificial_runjags(nspecies = 5, nsites = 20, nvisitspersite = 2, 
                            modeltype = "jsodm_lv", nlv = 2)
  simulation_occ_for_LV <- replicate(1000, poccupy(fit, lvvfromposterior = FALSE))
  bysimulation <- apply(simulation_occ_for_LV, MARGIN = c(1, 2), mean)
  bytheory <- poccupy(fit,  lvvfromposterior = FALSE, margLV = TRUE)[,,1]
  expect_equal(bytheory, bysimulation, tolerance = 0.01, ignore_attr = TRUE)
})

test_that("poccupy_species is correct without lv.v", {
  OccFmla = "~ 1"
  occ.b.min <- matrix(runif(5) , ncol = 1)
  occ.b.max <- occ.b.min + 1E-5
  pOccupancyFresh <- 1 - pnorm(0, mean = occ.b.min)
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6, 
                            occ.b.min = occ.b.min,
                            occ.b.max = occ.b.max,
                            OccFmla = OccFmla,
                            modeltype = "jsodm")
  pOccupy_sp <- poccupy(fit, usethetasummary = 1)
  diff <- Rfast::eachrow(pOccupy_sp, pOccupancyFresh, oper = "-")
  expect_true(max(abs(diff)) < 1E-3)
})

test_that("poccupy_species is correct with lv.v", {
  OccFmla = "~ 1"
  occ.b.min <- matrix(runif(5) , ncol = 1)
  occ.b.max <- occ.b.min + 1E-5
  pOccupancyFreshOdd <- 1 - pnorm(0, mean = occ.b.min + 0.6*0.995, sd = sqrt(1 - 0.6^2)) #occupancy probability for odd site indexes
                       #the 0.995 value is the standarised value of the odd-even lv.v (which is the first and only lv.v used here)
  pOccupancyFreshEven <- 1 - pnorm(0, mean = occ.b.min - 0.6*0.995, sd = sqrt(1 - 0.6^2)) #occupancy probability for odd site indexes
  # 1 - pnorm((-occ.b.min - 0.6) / (1 - 0.36), mean = 0, sd = 1) #occupancy probability for even site indexes
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6, 
                            occ.b.min = occ.b.min,
                            occ.b.max = occ.b.max,
                            lv.b.min = 0.6,
                            lv.b.max = 0.6,
                            OccFmla = OccFmla,
                            modeltype = "jsodm_lv",
                            nlv = 1)
  pOccupy_sp <- poccupy(fit, usethetasummary = 1, lvvfromposterior = TRUE)[,,1]
  diffOdd <- Rfast::eachrow(pOccupy_sp[1:nrow(pOccupy_sp) %% 2, ], pOccupancyFreshOdd, oper = "-")
  expect_true(max(abs(diffOdd)) < 1E-3)
  diffEven <- Rfast::eachrow(pOccupy_sp[!(1:nrow(pOccupy_sp) %% 2), ], pOccupancyFreshEven, oper = "-")
  expect_true(max(abs(diffEven)) < 1E-3)
})

test_that("pdet_occ() is correct", {
  ObsFmla = "~ 1"
  det.b.min <- matrix(runif(5) , ncol = 1)
  det.b.max <- det.b.min + 1E-5
  pDetCondOccFresh <- boot::inv.logit(det.b.min) #detection ease is constant across sites
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6, 
                            det.b.min = det.b.min,
                            det.b.max = det.b.max,
                            ObsFmla = ObsFmla,
                            modeltype = "jsodm_lv",
                            nlv = 2 )
  pDetCondOcc<- pdet_occ(fit, usethetasummary = 1)
  diff <- Rfast::eachrow(pDetCondOcc, pDetCondOccFresh, oper = "-")
  expect_true(max(abs(diff)) < 1E-3)
})

