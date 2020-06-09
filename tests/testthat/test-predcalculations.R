library(testthat);

context("Calculations of Predicted Probabilities")

test_that("poccupy_species is correct without LV", {
  OccFmla = "~ 1"
  u.b.min <- matrix(runif(5) , ncol = 1)
  u.b.max <- u.b.min + 1E-5
  pOccupancyFresh <- 1 - pnorm(0, mean = u.b.min)
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6, nlv = 0,
                            u.b.min = u.b.min,
                            u.b.max = u.b.max,
                            OccFmla = OccFmla)
  pOccupy_sp <- poccupy_species(fit, type = 1, conditionalLV = FALSE)
  diff <- Rfast::eachrow(pOccupy_sp, pOccupancyFresh, oper = "-")
  expect_true(max(abs(diff)) < 1E-3)
})

test_that("poccupy_species is correct with LV", {
  OccFmla = "~ 1"
  u.b.min <- matrix(runif(5) , ncol = 1)
  u.b.max <- u.b.min + 1E-5
  pOccupancyFreshOdd <- 1 - pnorm(0, mean = u.b.min + 0.6*0.995, sd = sqrt(1 - 0.6^2)) #occupancy probability for odd site indexes
                       #the 0.995 value is the standarised value of the odd-even LV values (which is the first and only LV used here)
  pOccupancyFreshEven <- 1 - pnorm(0, mean = u.b.min - 0.6*0.995, sd = sqrt(1 - 0.6^2)) #occupancy probability for odd site indexes
  # 1 - pnorm((-u.b.min - 0.6) / (1 - 0.36), mean = 0, sd = 1) #occupancy probability for even site indexes
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6, nlv = 1,
                            u.b.min = u.b.min,
                            u.b.max = u.b.max,
                            lv.coef.min = 0.6,
                            lv.coef.max = 0.6,
                            OccFmla = OccFmla)
  pOccupy_sp <- poccupy_species(fit, type = 1, conditionalLV = TRUE)
  diffOdd <- Rfast::eachrow(pOccupy_sp[1:nrow(pOccupy_sp) %% 2, ], pOccupancyFreshOdd, oper = "-")
  expect_true(max(abs(diffOdd)) < 1E-3)
  diffEven <- Rfast::eachrow(pOccupy_sp[!(1:nrow(pOccupy_sp) %% 2), ], pOccupancyFreshEven, oper = "-")
  expect_true(max(abs(diffEven)) < 1E-3)
})

test_that("pdetect_condoccupied is correct", {
  ObsFmla = "~ 1"
  v.b.min <- matrix(runif(5) , ncol = 1)
  v.b.max <- v.b.min + 1E-5
  pDetCondOccFresh <- boot::inv.logit(v.b.min) #detection ease is constant across sites
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6, nlv = 2,
                            v.b.min = v.b.min,
                            v.b.max = v.b.max,
                            ObsFmla = ObsFmla)
  pDetCondOcc<- pdetect_condoccupied(fit, type = 1)
  diff <- Rfast::eachrow(pDetCondOcc, pDetCondOccFresh, oper = "-")
  expect_true(max(abs(diffOdd)) < 1E-3)
})

test_that("pdetect_indvisit is correct with LV", {
  ObsFmla = "~ 1"
  OccFmla = "~ 1"
  u.b.min <- matrix(runif(5) , ncol = 1)
  u.b.max <- u.b.min + 1E-5
  v.b.min <- matrix(runif(5) , ncol = 1)
  v.b.max <- v.b.min + 1E-5
  pOccupancyFreshOdd <- 1 - pnorm(0, mean = u.b.min + 0.6*0.995, sd = sqrt(1 - 0.6^2)) #occupancy probability for odd site indexes
  #the 0.995 value is the standarised value of the odd-even LV values (which is the first and only LV used here)
  pOccupancyFreshEven <- 1 - pnorm(0, mean = u.b.min - 0.6*0.995, sd = sqrt(1 - 0.6^2)) #occupancy probability for odd site indexes
  # 1 - pnorm((-u.b.min - 0.6) / (1 - 0.36), mean = 0, sd = 1) #occupancy probability for even site indexes
  pDetCondOccFresh <- boot::inv.logit(v.b.min) #detection ease is constant across sites
  pDetFreshOdd <- pOccupancyFreshOdd * pDetCondOccFresh
  pDetFreshEven <- pOccupancyFreshEven * pDetCondOccFresh
  
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6, nlv = 1,
                            u.b.min = u.b.min,
                            u.b.max = u.b.max,
                            v.b.min = v.b.min,
                            v.b.max = v.b.max,
                            lv.coef.min = 0.6,
                            lv.coef.max = 0.6,
                            OccFmla = OccFmla,
                            ObsFmla = ObsFmla)
  pdetect <- pdetect_indvisit(fit, type = 1)
  
  diffOdd <- Rfast::eachrow(pdetect[fit$data$ModelSite %% 2, ], pDetFreshOdd, oper = "-")
  expect_true(max(abs(diffOdd)) < 1E-3)
  diffEven <- Rfast::eachrow(pdetect[!(fit$data$ModelSite %% 2), ], pDetFreshEven, oper = "-")
  expect_true(max(abs(diffEven)) < 1E-3)
})

test_that("pdetect_condoccupied and poccupy_species keeps ordering of sites / visits", {
  ObsFmla = "~ UpVisit - 1"
  OccFmla = "~ UpSite - 1"
  u.b.min <- matrix((1:5)/8 , ncol = 1)
  u.b.max <- u.b.min + 1E-5
  v.b.min <- matrix((1:5)/8 , ncol = 1)
  v.b.max <- v.b.min + 1E-5
  
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6, nlv = 0,
                            u.b.min = u.b.min,
                            u.b.max = u.b.max,
                            v.b.min = v.b.min,
                            v.b.max = v.b.max,
                            OccFmla = OccFmla,
                            ObsFmla = ObsFmla)
  # should be a matrix that increase down the rows 
  # when scaled UpSite positive it should go up across the columns
  # when scaled UpSite negative it should go down the columns
  pOccupancy <- poccupy_species(fit, type = 1, conditionalLV = FALSE)
  expect_true(min(pOccupancy[-1, ] - pOccupancy[-nrow(pOccupancy), ]) > 0)
  expect_true(min(sign(fit$data$Xocc[, "UpSite"]) * (pOccupancy[, -1] - pOccupancy[, -ncol(pOccupancy)])) > 0)
  
  # Below should go up with visit (row)
  # when scaled UpVisit positive it should go up across the columns
  # when scaled UpVisit negative it should go down the columns
  pDetCondOcc <- pdetect_condoccupied(fit, type = 1)
  expect_true(min(pDetCondOcc[-1, ] - pDetCondOcc[-nrow(pDetCondOcc), ]) > 0)
  expect_true(min(sign(fit$data$Xobs[, "UpVisit"]) * (pDetCondOcc[, -1] - pDetCondOcc[, -ncol(pDetCondOcc)])) > 0)
  
})
