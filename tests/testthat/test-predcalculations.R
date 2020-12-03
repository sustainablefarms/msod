library(testthat);

context("Calculations of Predicted Probabilities")

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
  pOccupy_sp <- poccupy_species(fit, type = 1, conditionalLV = FALSE)
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
  pOccupy_sp <- poccupy_species(fit, type = 1, conditionalLV = TRUE)
  diffOdd <- Rfast::eachrow(pOccupy_sp[1:nrow(pOccupy_sp) %% 2, ], pOccupancyFreshOdd, oper = "-")
  expect_true(max(abs(diffOdd)) < 1E-3)
  diffEven <- Rfast::eachrow(pOccupy_sp[!(1:nrow(pOccupy_sp) %% 2), ], pOccupancyFreshEven, oper = "-")
  expect_true(max(abs(diffEven)) < 1E-3)
})

test_that("pdetect_condoccupied is correct", {
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
  pDetCondOcc<- pdetect_condoccupied(fit, type = 1)
  diff <- Rfast::eachrow(pDetCondOcc, pDetCondOccFresh, oper = "-")
  expect_true(max(abs(diff)) < 1E-3)
})

test_that("pdetect_indvisit is correct with lv.v", {
  ObsFmla = "~ 1"
  OccFmla = "~ 1"
  occ.b.min <- matrix(runif(5) , ncol = 1)
  occ.b.max <- occ.b.min + 1E-5
  det.b.min <- matrix(runif(5) , ncol = 1)
  det.b.max <- det.b.min + 1E-5
  pOccupancyFreshOdd <- 1 - pnorm(0, mean = occ.b.min + 0.6*0.995, sd = sqrt(1 - 0.6^2)) #occupancy probability for odd site indexes
  #the 0.995 value is the standarised value of the odd-even lv.v (which is the first and only lv.v used here)
  pOccupancyFreshEven <- 1 - pnorm(0, mean = occ.b.min - 0.6*0.995, sd = sqrt(1 - 0.6^2)) #occupancy probability for odd site indexes
  # 1 - pnorm((-occ.b.min - 0.6) / (1 - 0.36), mean = 0, sd = 1) #occupancy probability for even site indexes
  pDetCondOccFresh <- boot::inv.logit(det.b.min) #detection ease is constant across sites
  pDetFreshOdd <- pOccupancyFreshOdd * pDetCondOccFresh
  pDetFreshEven <- pOccupancyFreshEven * pDetCondOccFresh
  
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6,
                            occ.b.min = occ.b.min,
                            occ.b.max = occ.b.max,
                            det.b.min = det.b.min,
                            det.b.max = det.b.max,
                            lv.b.min = 0.6,
                            lv.b.max = 0.6,
                            OccFmla = OccFmla,
                            ObsFmla = ObsFmla,
                            modeltype = "jsodm_lv",
                            nlv = 1)
  pdetect <- pdetect_indvisit(fit, type = 1)
  
  diffOdd <- Rfast::eachrow(pdetect[fit$data$ModelSite %% 2, ], pDetFreshOdd, oper = "-")
  expect_true(max(abs(diffOdd)) < 1E-3)
  diffEven <- Rfast::eachrow(pdetect[!(fit$data$ModelSite %% 2), ], pDetFreshEven, oper = "-")
  expect_true(max(abs(diffEven)) < 1E-3)
})

test_that("pdetect_condoccupied and poccupy_species keeps ordering of sites / visits", {
  ObsFmla = "~ UpVisit - 1"
  OccFmla = "~ UpSite - 1"
  occ.b.min <- matrix((1:5)/3 , ncol = 1)
  occ.b.max <- occ.b.min + 1E-5
  det.b.min <- matrix((1:5)/8 , ncol = 1)
  det.b.max <- det.b.min + 1E-5
  
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6,
                            occ.b.min = occ.b.min,
                            occ.b.max = occ.b.max,
                            det.b.min = det.b.min,
                            det.b.max = det.b.max,
                            OccFmla = OccFmla,
                            ObsFmla = ObsFmla,
                            modeltype = "jsodm")
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
  
  # I expect pdetect_indvisit to have the same rows as Xobs
  pdetect <- pdetect_indvisit(fit, type = 1, conditionalLV = FALSE)
  # probability of occupancy and detection increases with Site and Visit, 
  # so full probability of detection increases except possibly whenever ModelSite returns back to 1,
  # which never happens in currently: simulated ModelSite is non-decreasing, all the ModelSite=1 visits occur in the first few rows
  rowdiff <- pdetect[-1, ] - pdetect[-nrow(pdetect), ]
  expect_true(min( rowdiff ) > 0)
  
  # across the rows (species) when UpSite < 0 and UpVisit < 0, I expect the full probability of detection to decrease
  coldiff <-  pdetect[, -1] - pdetect[, -ncol(pdetect)]
  expect_lt(max(coldiff[fit$data$Xocc[fit$data$ModelSite, "UpSite"] < 0 & fit$data$Xobs[, "UpVisit"] < 0, ]), 0)
  # across the rows (species) when UpSite > 0 and UpVisit > 0, I expect the full probability of detection to increase
  expect_gt(min(coldiff[fit$data$Xocc[fit$data$ModelSite, "UpSite"] > 0 & fit$data$Xobs[, "UpVisit"] > 0, ]), 0)
  # in other rows I expect there to be a strong pattern, it is probably predictable
})

test_that("pdetect_condoccupied and poccupy_species ordering of sites / visits with lv.v", {
  ObsFmla = "~ 1"
  OccFmla = "~ 1"
  occ.b.min <- matrix((1:5)/3 , ncol = 1)
  occ.b.max <- occ.b.min + 1E-5
  det.b.min <- matrix((1:5)/8 , ncol = 1)
  det.b.max <- det.b.min + 1E-5
  lv.b.min <- matrix((1:5)/8, ncol = 2, nrow = 5, byrow = FALSE)
  lv.b.max <- lv.b.min + 1E-5
  
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6,
                            occ.b.min = occ.b.min,
                            occ.b.max = occ.b.max,
                            det.b.min = det.b.min,
                            det.b.max = det.b.max,
                            lv.b.min = lv.b.min,
                            lv.b.max = lv.b.max,
                            OccFmla = OccFmla,
                            ObsFmla = ObsFmla,
                            modeltype = "jsodm_lv",
                            nlv = 2)
  lv.v <- bugsvar2matrix(get_theta(fit, type = 1), "lv.v", 1:100, 1:2)
  # should be a matrix that increase down the rows 
  # when scaled UpSite positive it should go up across the columns
  # when scaled UpSite negative it should go down the columns
  pOccupancy <- poccupy_species(fit, type = 1, conditionalLV = TRUE)
  expect_equivalent(apply(pOccupancy[lv.v[, 1] > 0 & lv.v[, 2] > 0, ], 2, sd), rep(0, ncol(pOccupancy)))
  expect_equivalent(apply(pOccupancy[lv.v[, 1] > 0 & lv.v[, 2] < 0, ], 2, sd), rep(0, ncol(pOccupancy)))
  expect_equivalent(apply(pOccupancy[lv.v[, 1] < 0 & lv.v[, 2] < 0, ], 2, sd), rep(0, ncol(pOccupancy)))
  expect_equivalent(apply(pOccupancy[lv.v[, 1] < 0 & lv.v[, 2] > 0, ], 2, sd), rep(0, ncol(pOccupancy)))
  expect_gt(min(apply(pOccupancy[lv.v[, 2] > 0, ], 2, sd)), 0)
  
  # expect occupancy to go up with species
  expect_true(min(pOccupancy[, -1] - pOccupancy[, -ncol(pOccupancy)]) > 0)
  
  # I expect pdetect_indvisit to have the same rows as Xobs
  pdetect <- pdetect_indvisit(fit, type = 1, conditionalLV = FALSE)
  # Similar pattern of equal probabilities as pOccupancy
  expect_equivalent(apply(pdetect[lv.v[fit$data$ModelSite, 1] > 0 & lv.v[fit$data$ModelSite, 2] > 0, ], 2, sd), rep(0, ncol(pdetect)))
  expect_equivalent(apply(pdetect[lv.v[fit$data$ModelSite, 1] > 0 & lv.v[fit$data$ModelSite, 2] < 0, ], 2, sd), rep(0, ncol(pdetect)))
  expect_equivalent(apply(pdetect[lv.v[fit$data$ModelSite, 1] < 0 & lv.v[fit$data$ModelSite, 2] < 0, ], 2, sd), rep(0, ncol(pdetect)))
  expect_equivalent(apply(pdetect[lv.v[fit$data$ModelSite, 1] < 0 & lv.v[fit$data$ModelSite, 2] > 0, ], 2, sd), rep(0, ncol(pdetect)))
})

test_that("Expected number of detections for each model site matches the simulated mean number of observations", {
  fit <- artificial_runjags(nspecies = 100, modeltype = "jsodm_lv", nlv = 2)
  Endetections <- Endetect_modelsite(fit, type = 1)$E_ndetect
  
  simulate_combinevisits <- function(fit, esttype){
    sim <- simulate_detections(fit, esttype = 1)
    ndetect_ModelSite <- cbind(ModelSite = fit$data$ModelSite, sim) %>%
      as_tibble() %>%
      group_by(ModelSite) %>%
      summarise_all(sum) %>%
      dplyr::select(-ModelSite) %>%
      as.matrix()
    return(ndetect_ModelSite)
  }
  
  # simulate multiple times
  n <- 1000
  sims <- replicate(n, simulate_combinevisits(fit, esttype = 1))
  meandetnum <- apply(sims, MARGIN = c(1, 2), mean)
  sddetnum <- apply(sims, MARGIN = c(1, 2), sd)
  skewdetnum <- apply(sims, MARGIN = c(1, 2), Rfast::skew)
  SE <- sddetnum/ sqrt(n)
  
  expect_equal(mean( abs(Endetections - meandetnum)/ SE < 2),  
                      0.95, tolerance = 0.03)
  # image(abs(Endetections - meandetnum)/ SE < 2)
})
