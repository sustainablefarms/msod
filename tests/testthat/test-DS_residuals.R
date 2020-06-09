library(testthat);

context("Tests of Dunn-Smyth Residuals")



test_that("Seed argument works for setting final residual computation", {
  a <- cdfvals_2_dsres_discrete(rep(0.5, 6), rep(0.6, 6), seed = 5346)
  b <- cdfvals_2_dsres_discrete(rep(0.5, 6), rep(0.6, 6), seed = 5346)
  c <- cdfvals_2_dsres_discrete(rep(0.5, 6), rep(0.6, 6), seed = NULL)
  d <- cdfvals_2_dsres_discrete(rep(0.5, 6), rep(0.6, 6), seed = NULL)
  expect_identical(a, b)
  expect_false(isTRUE(all.equal(a, c)))
  expect_false(isTRUE(all.equal(c, d)))
  expect_false(isTRUE(all.equal(a, d)))
})

test_that("PDF and CDF computations are correct", {
  expect_equal(numdet_pdf(0, c(0.1, 0.2)), 0.9 * 0.8)
  expect_equal(numdet_pdf(2, c(0.1, 0.2)), 0.1 * 0.2)
  expect_equal(numdet_pdf(1, c(0.1, 0.2)), 1 - 0.9 * 0.8 - 0.1 * 0.2)
  expect_equal(numdet_pdf(-1, c(0.1, 0.2)), 0)
  expect_equal(numdet_pdf(1.3, c(0.1, 0.2)), 0)
  expect_equal(numdet_pdf(3, c(0.1, 0.2)), 0)
  
  expect_equal(numdet_cdf(0, c(0.1, 0.2)),  0.9 * 0.8)
  expect_equal(numdet_cdf(1, c(0.1, 0.2)),  1 - 0.1 * 0.2)
  
  pDetected <- runif(5)
  expect_equal(numdet_cdf(0.3, pDetected),  numdet_cdf(0, pDetected))
  expect_equal(numdet_cdf(-0.1, pDetected),  0)
  expect_equal(numdet_cdf(length(pDetected), pDetected),  1)
  expect_equal(numdet_cdf(2.1, pDetected),  numdet_cdf(2, pDetected))
  expect_gt(numdet_cdf(3, pDetected), numdet_cdf(2, pDetected))
  
  # compare to simulations
  s <- summary(as.factor(simDetectedDistr(500000, pDetected)))/500000
  p <- vapply(0:length(pDetected), numdet_pdf, pDetected = pDetected, FUN.VALUE = 3)
  names(p) <- 0:length(pDetected)
  expect_equal(s, p, tolerance = 1E-3)
})

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


test_that("Occupancy and Detection Residuals match fresh conversion from fit to raw, with LV", {
  ObsFmla = "~ 1"
  OccFmla = "~ 1"
  u.b.min <- matrix(runif(5) , ncol = 1)
  u.b.max <- u.b.min + 1E-5
  v.b.min <- matrix(runif(5) , ncol = 1)
  v.b.max <- v.b.min + 1E-5
  
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6, nlv = 1,
                            u.b.min = u.b.min,
                            u.b.max = u.b.max,
                            v.b.min = v.b.min,
                            v.b.max = v.b.max,
                            lv.coef.min = 0.6,
                            lv.coef.max = 0.6,
                            OccFmla = OccFmla,
                            ObsFmla = ObsFmla)
  
  #' @param preds is a dataframe with columns Species, ModelSite, pOccupancy, and pDetected_cond.
  #'  pOccupancy is the probability of ModelSite being occupied.
  #'  pDetected_cond is the probability of detecting the species, given the species occupies the ModelSite.
  #' @param obs is a dataframe with columns Species, ModelSite, and Detected
  pOccupancy <- poccupy_species(fit, type = 1)
  pOccupancy <- cbind(ModelSite = 1:nrow(fit$data$Xocc), pOccupancy) %>%
    as_tibble() %>%
    pivot_longer(-ModelSite, names_to = "Species", values_to = "pOccupancy")
  pDetCondOcc <- cbind(ModelSite = fit$data$ModelSite, VisitId = 1:nrow(fit$data$Xobs), pdetect_condoccupied(fit, type = 1)) %>%
    as_tibble() %>%
    pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "pDetected_cond")
  preds <- inner_join(pOccupancy, pDetCondOcc, by = c("ModelSite", "Species")) %>% arrange(VisitId, Species, ModelSite)
  obs <- cbind(ModelSite = as.numeric(fit$data$ModelSite), VisitId = 1:nrow(fit$data$Xobs), fit$data$y) %>%
    as_tibble() %>%
    pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "Detected") %>%
    arrange(VisitId, Species, ModelSite)
  
  ds_residuals <- ds_occupancy_residuals.raw(preds, obs)
  shapiroresults <- shapiro.test(ds_residuals$OccupancyResidual)
  expect_gt(shapiroresults$p.value, 0.05) #if things are good this test will fail 1/20 times
  
  #' @param preds is a dataframe with columns Species, ModelSite, and pDetected
  #' @param obs is a dataframe with columns Species, ModelSite, and Detected
  preds <- cbind(ModelSite = fit$data$ModelSite, VisitId = 1:nrow(fit$data$Xobs), pdetect_indvisit(fit, type = 1)) %>%
    as_tibble() %>%
    pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "pDetected") %>%
    arrange(VisitId, Species, ModelSite)
  obs <- cbind(ModelSite = as.numeric(fit$data$ModelSite), VisitId = 1:nrow(fit$data$Xobs), fit$data$y) %>%
    as_tibble() %>%
    pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "Detected") %>%
    arrange(VisitId, Species, ModelSite)
  
  ds_residuals <- ds_detection_residuals.raw(preds, obs)
  shapiroresults <- shapiro.test(ds_residuals$DetectionResidual)
  expect_gt(shapiroresults$p.value, 0.05) #if things are good this test will fail 1/20 times
})

test_that("Detection residuals sensible and gaussian for very raw simulated data", {
  # simulate a data set
  species <- c("A", "B", "C", "D")
  sites <- c(1:1000)
  occupancy <- expand.grid(Species = species, ModelSite = sites)
  occupancy$pOccupancy <- runif(nrow(occupancy))
  
  visits <- rbind(occupancy, occupancy, occupancy[1:(2 * length(species)), ],
                  occupancy, occupancy, occupancy)
  specdet_base <- runif(length(species), 0.2, 0.8)
  names(specdet_base) <- species
  visitdetprob_offset <- runif(nrow(visits), -0.1, 0.1)
  preds <- visits %>%
    dplyr::mutate(pDetected = pOccupancy * (specdet_base[Species] + visitdetprob_offset))
  
  # simulate observations from the given predicted values (purely for testing)
  obs <- preds %>% dplyr::select(Species, ModelSite)
  obs$Detected <- mapply(rbinom, prob = preds$pDetected, MoreArgs = list(n = 1, size = 1))

  # test residuals
  expect_equal(ds_detection_residuals.raw(preds, obs, seed = 1234),
               ds_detection_residuals.raw(preds, obs, seed = 1234))
  expect_true(any(ds_detection_residuals.raw(preds, obs) != 
               ds_detection_residuals.raw(preds, obs)))
  
  # residuals should have a standard gaussian distribution when observations generated by model given by preds
  ds_residuals <- ds_detection_residuals.raw(preds, obs)
  shapiroresults <- shapiro.test(ds_residuals$DetectionResidual)
  expect_gt(shapiroresults$p.value, 0.05) #if things are good this test will fail 1/20 times
})

test_that("Occupancy residuals sensible for simulated data", {
  # simulate a data set
  species <- c("A", "B", "C", "D")
  sites <- c(1:100)
  occupancy <- expand.grid(Species = species, ModelSite = sites)
  occupancy$pOccupancy <- runif(nrow(occupancy))
  occupancy$Occupied <- mapply(rbinom, prob = occupancy$pOccupancy, MoreArgs = list(n = 1, size = 1))
  
  visits <- rbind(occupancy, occupancy, occupancy[1:(2 * length(species)), ],
                  occupancy, occupancy, occupancy)
  specdet_base <- runif(length(species), 0.2, 0.8)
  names(specdet_base) <- species
  visitdetprob_offset <- runif(nrow(visits), -0.1, 0.1)
  preds <- visits %>%
    dplyr::mutate(pDetected_cond = specdet_base[Species] + visitdetprob_offset)
  
  # simulate observations from the given predicted values (purely for testing)
  obs <- preds %>% dplyr::select(Species, ModelSite)
  obs$Detected <- mapply(rbinom, prob = preds$Occupied * preds$pDetected_cond, MoreArgs = list(n = 1, size = 1))

  preds <- preds %>% dplyr::select(Species, ModelSite, pOccupancy, pDetected_cond)
  
  # test residuals response to seed
  expect_equal(ds_occupancy_residuals.raw(preds, obs, seed = 1234),
               ds_occupancy_residuals.raw(preds, obs, seed = 1234))
  expect_true(any(ds_occupancy_residuals.raw(preds, obs) != 
                    ds_occupancy_residuals.raw(preds, obs)))
  
  # residuals should have a standard gaussian distribution when observations generated by model given by preds
  ds_residuals <- ds_occupancy_residuals.raw(preds, obs)
  shapiroresults <- shapiro.test(ds_residuals$OccupancyResidual)
  expect_gt(shapiroresults$p.value, 0.05) #if things are good this test will fail 1/20 times
})


test_that("DS Detection Residuals are Gaussian for Artificial Fitted Object made of Common Species", {
  # simulate a fitted object
  fit <- artificial_runjags(nspecies = 5, nsites = 500, nvisitspersite = 5, nlv = 2,
                            u.b.min = 0.95)
  
  # compute residuals 
  resid_det <- ds_detection_residuals.fit(fit, type = 1)
  shapiro_det_residual <- resid_det %>% dplyr::select(-ModelSite) %>%  as.matrix() %>%  as.vector() %>%
    shapiro.test()
  expect_gt(shapiro_det_residual$p.value, 0.05)
  
  # resid_det %>% dplyr::select(-ModelSite) %>% as.matrix() %>% as.vector() %>% qqnorm()
  # abline(a = 0, b = 1)
})

test_that("DS Detection Residuals are Gaussian for Artificial Fitted Object made of Rare Species", {
  # simulate a fitted object
  fit <- artificial_runjags(nspecies = 5, nsites = 500, nvisitspersite = 5, nlv = 2,
                            u.b.max = -0.9)
  
  # compute residuals 
  resid_det <- ds_detection_residuals.fit(fit, type = 1)
  shapiro_det_residual <- resid_det %>% dplyr::select(-ModelSite) %>%  as.matrix() %>%  as.vector() %>%
    shapiro.test()
  expect_gt(shapiro_det_residual$p.value, 0.05)
  
  # resid_det %>% dplyr::select(-ModelSite) %>% as.matrix() %>% as.vector() %>% qqnorm()
  # abline(a = 0, b = 1)
})

test_that("DS Detection Residuals are Gaussian for Artificial Fitted Object made of Variety of Species", {
  # simulate a fitted object
  fit <- artificial_runjags(nspecies = 5, nsites = 500, nvisitspersite = 5, nlv = 2)
  
  # compute residuals 
  resid_det <- ds_detection_residuals.fit(fit, type = 1)
  shapiro_det_residual <- resid_det %>% dplyr::select(-ModelSite) %>%  as.matrix() %>%  as.vector() %>%
    shapiro.test()
  expect_gt(shapiro_det_residual$p.value, 0.05)
  
  # resid_det %>% dplyr::select(-ModelSite) %>% as.matrix() %>% as.vector() %>% qqnorm()
  # abline(a = 0, b = 1)
})

test_that("DS Occupancy Residuals are Gaussian for Artificial Fitted Object", {
  # simulate a fitted object
  fit <- artificial_runjags(nspecies = 5, nsites = 500, nvisitspersite = 5, nlv = 2)
  
  # compute residuals 
  resid_occ <- ds_occupancy_residuals.fit(fit, type = 1)
  shapiro_occ_residual <- resid_occ %>% dplyr::select(-ModelSite) %>%  as.matrix() %>%  as.vector() %>%
    shapiro.test()
  expect_gt(shapiro_occ_residual$p.value, 0.05)
  # resid_occ %>% dplyr::select(-ModelSite) %>% as.matrix() %>% as.vector() %>% qqnorm()
  # abline(a = 0, b = 1)
})

test_that("Number of detection residuals computed for a highly common species at a single ModelSite", {
  species <- c("A")#, "B", "C", "D")
  sites <- c(1)
  occupancy <- expand.grid(Species = species, ModelSite = sites)
  occupancy$pOccupancy <- 0.9
  
  visits <- rbind(occupancy, occupancy, occupancy, occupancy, occupancy,
                  occupancy, occupancy, occupancy, occupancy, occupancy)
  specdet_base <- runif(length(species), 0.2, 0.8)
  names(specdet_base) <- species
  visitdetprob_offset <- runif(nrow(visits), -0.1, 0.1)
  preds <- visits %>%
    dplyr::mutate(pDetected = pOccupancy * (specdet_base[Species] + visitdetprob_offset))
  
  # simulate observations from the given predicted values (purely for testing)
  obs <- preds %>% dplyr::select(Species, ModelSite)
  obs$Detected <- mapply(rbinom, prob = preds$pDetected, MoreArgs = list(n = 1, size = 1))
  
  if (sum(obs$Detected) > 0){expect_equal(nrow(ds_detection_residuals.raw(preds, obs)), 1)}
  if (sum(obs$Detected) == 0){expect_error(ds_detection_residuals.raw(preds, obs))}
})

test_that("Zero detection residuals computed for a very rare species at a single ModelSite", {
  species <- c("A")#, "B", "C", "D")
  sites <- c(1)
  occupancy <- expand.grid(Species = species, ModelSite = sites)
  occupancy$pOccupancy <- 0.01
  
  visits <- rbind(occupancy, occupancy, occupancy, occupancy, occupancy)
  specdet_base <- runif(length(species), 0.2, 0.8)
  names(specdet_base) <- species
  visitdetprob_offset <- runif(nrow(visits), -0.1, 0.1)
  preds <- visits %>%
    dplyr::mutate(pDetected = pOccupancy * (specdet_base[Species] + visitdetprob_offset))
  
  # simulate observations from the given predicted values (purely for testing)
  obs <- preds %>% dplyr::select(Species, ModelSite)
  obs$Detected <- mapply(rbinom, prob = preds$pDetected, MoreArgs = list(n = 1, size = 1))
  
  if (sum(obs$Detected) == 0){expect_error(ds_detection_residuals.raw(preds, obs))}
})

