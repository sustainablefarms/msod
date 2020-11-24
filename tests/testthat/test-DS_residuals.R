library(discreteRV)

context("Tests of Dunn-Smyth Residuals")

test_that("Seed argument works for setting cumulative distribution computations", {
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
  
  pDetected <- runif(5, min = 0.3, max = 0.7)
  expect_equal(numdet_cdf(0.3, pDetected),  numdet_cdf(0, pDetected))
  expect_equal(numdet_cdf(-0.1, pDetected),  0)
  expect_equal(numdet_cdf(length(pDetected), pDetected),  1)
  expect_equal(numdet_cdf(2.1, pDetected),  numdet_cdf(2, pDetected))
  expect_gt(numdet_cdf(3, pDetected), numdet_cdf(2, pDetected))
  
  # compare to simulations
  n <- 1000
  s <- summary(as.factor(simDetectedDistr(n, pDetected)))/n
  p <- vapply(0:length(pDetected), numdet_pdf, pDetected = pDetected, FUN.VALUE = 3)
  names(p) <- 0:length(pDetected)
  
  # can approximate *poorly* (?) expected error of s by assuming that each value is the result of a bernoulli distribution
  # s will be approximate Gaussian with mean p, and variance p(1 - p)/n
  sd_s <- sqrt(p * (1 - p) / n)
  expect_true( (abs(s - p) < 2 * sd_s)[[3]])
  #above should be true 95% of the for *each* value of p.
  # As they are highly correlated, I'm not sure how often I expect the below to produce errors.
  expect_true(all(abs(s - p) < 3 * sd_s))
})

test_that("PDF and CDF computations are correct using discreteRV", {
  nvisits <- 5
  pDetected <- runif(nvisits)
  v1 <- RV(outcomes = c(0, 1), probs = c(1-pDetected[[1]], pDetected[[1]]))
  indicator_rvs <- lapply(pDetected, function(x) RV(outcomes = c(0, 1), probs = c(1-x, x)) )
  sum_rv <- Reduce("+", indicator_rvs)
  
  p_by_discreteRV <- probs(sum_rv)
  p <- vapply(0:length(pDetected), numdet_pdf, pDetected = pDetected, FUN.VALUE = 3)
  names(p) <- 0:length(pDetected)
  expect_equivalent(p, p_by_discreteRV)
  
  cdf <- vapply(0:length(pDetected), numdet_cdf, pDetected = pDetected, FUN.VALUE = 3)
  expect_equivalent(cumsum(p_by_discreteRV), cdf)
})


test_that("Conditioning on Greater Than Zero Correct", {
  pDetected <- runif(3)  #detections at each visit
  pdf <- vapply(0:length(pDetected), function(x) numdet_pdf(x, pDetected), FUN.VALUE = 0.0)
  cdf <- vapply(0:length(pDetected), function(x) numdet_cdf(x, pDetected), FUN.VALUE = 0.0)
  expect_equivalent(cdf, cumsum(pdf))
  
  expect_equal(condition_nonzero.cdf(cdf[[1]], cdf[[2]]), pdf[[2]] / (1 - pdf[[1]]))
  expect_equal(condition_nonzero.cdf(cdf[[1]], cdf[[2]]), condition_nonzero.pdf(cdf[[1]], pdf[[2]]))
  expect_equal(condition_nonzero.cdf(cdf[[1]], cdf[[1 + length(pDetected)]]), 1)
  expect_equal(condition_nonzero.cdf(cdf[[1]], cdf[[1]]), 0)
  expect_equal(condition_nonzero.pdf(cdf[[length(pDetected)]], pdf[[1 + length(pDetected)]]), 1)
})

test_that("Occupancy and Detection Residuals Gaussian for fresh conversion from fit to raw, with LV", {
  ObsFmla = "~ 1"
  OccFmla = "~ 1"
  u.b.min <- matrix(runif(5) , ncol = 1)
  u.b.max <- u.b.min + 1E-5
  v.b.min <- matrix(runif(5) , ncol = 1)
  v.b.max <- v.b.min + 1E-5
  
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 6, 
                            u.b.min = u.b.min,
                            u.b.max = u.b.max,
                            v.b.min = v.b.min,
                            v.b.max = v.b.max,
                            lv.coef.min = 0.6,
                            lv.coef.max = 0.6,
                            OccFmla = OccFmla,
                            ObsFmla = ObsFmla,
                            modeltype = "jsodm_lv",
                            nlv = 1)
  
  #' @param preds is a dataframe with columns Species, ModelSite, pOccupancy, and pDetected_cond.
  #'  pOccupancy is the probability of ModelSite being occupied.
  #'  pDetected_cond is the probability of detecting the species, given the species occupies the ModelSite.
  #' @param obs is a dataframe with columns Species, ModelSite, and Detected
  pOccupancy <- poccupy_species(fit, type = 1)
  pOccupancy <- cbind(ModelSite = 1:nrow(fit$data$Xocc), pOccupancy) %>%
    as_tibble() %>%
    tidyr::pivot_longer(-ModelSite, names_to = "Species", values_to = "pOccupancy")
  pDetCondOcc <- cbind(ModelSite = fit$data$ModelSite, VisitId = 1:nrow(fit$data$Xobs), pdetect_condoccupied(fit, type = 1)) %>%
    as_tibble() %>%
    tidyr::pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "pDetected_cond")
  preds <- inner_join(pOccupancy, pDetCondOcc, by = c("ModelSite", "Species")) %>% arrange(VisitId, Species, ModelSite)
  obs <- cbind(ModelSite = as.numeric(fit$data$ModelSite), VisitId = 1:nrow(fit$data$Xobs), fit$data$y) %>%
    as_tibble() %>%
    tidyr::pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "Detected") %>%
    arrange(VisitId, Species, ModelSite)
  
  ds_residuals <- ds_occupancy_residuals.raw(preds, obs)
  expect_gt(shapiro.test(ds_residuals$OccupancyResidual)$p.value, 0.01) #if things are good this test will fail 1/100 times
  
  #' @param preds is a dataframe with columns Species, ModelSite, and pDetected
  #' @param obs is a dataframe with columns Species, ModelSite, and Detected
  preds <- cbind(ModelSite = fit$data$ModelSite, VisitId = 1:nrow(fit$data$Xobs), pdetect_indvisit(fit, type = 1)) %>%
    as_tibble() %>%
    tidyr::pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "pDetected") %>%
    arrange(VisitId, Species, ModelSite)
  obs <- cbind(ModelSite = as.numeric(fit$data$ModelSite), VisitId = 1:nrow(fit$data$Xobs), fit$data$y) %>%
    as_tibble() %>%
    tidyr::pivot_longer(-c(ModelSite, VisitId), names_to = "Species", values_to = "Detected") %>%
    arrange(VisitId, Species, ModelSite)
  
  ds_residuals <- ds_detection_residuals.raw(preds, obs)
  expect_gt(shapiro.test(ds_residuals$DetectionResidual)$p.value, 0.01) #if things are good this test will fail 1/100 times
})

test_that("Detection residuals sensible and gaussian for very raw simulated data", {
  # simulate a data set
  species <- LETTERS
  sites <- c(1:100)
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
  expect_gt(shapiroresults$p.value, 0.01) #if things are good this test will fail 1/100 times
})

test_that("Detection residuals sensible and gaussian for very raw simulated data", {
  # simulate a data set
  species <- LETTERS
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
  expect_gt(shapiroresults$p.value, 0.01) #if things are good this test will fail 1/100 times
})

test_that("DS Detection Residuals are Gaussian for Artificial Fitted Object made of Common Species", {
  # simulate a fitted object
  fit <- artificial_runjags(nspecies = 5, nsites = 500, nvisitspersite = 5,
                            u.b.min = 0.95,
                            modeltype = "jsodm_lv",
                            nlv = 2)
  
  # compute residuals 
  resid_det <- ds_detection_residuals.fit(fit, type = 1)
  shapiro_det_residual <- resid_det %>% dplyr::select(-ModelSite) %>%  as.matrix() %>%  as.vector() %>%
    shapiro.test()
  expect_gt(shapiro_det_residual$p.value, 0.01)
  
  # resid_det %>% dplyr::select(-ModelSite) %>% as.matrix() %>% as.vector() %>% qqnorm()
  # abline(a = 0, b = 1)
})
# 1 out of 20 fails :))

test_that("DS Detection Residuals are Gaussian for Artificial Fitted Object made of Rare Species", {
  # simulate a fitted object
  fit <- artificial_runjags(nspecies = 5, nsites = 500, nvisitspersite = 5,
                            u.b.max = -0.9,
                            modeltype = "jsodm_lv",
                            nlv = 2)
  
  # compute residuals 
  resid_det <- ds_detection_residuals.fit(fit, type = 1)
  shapiro_det_residual <- resid_det %>% dplyr::select(-ModelSite) %>%  as.matrix() %>%  as.vector() %>%
    shapiro.test()
  expect_gt(shapiro_det_residual$p.value, 0.01)
  
  # resid_det %>% dplyr::select(-ModelSite) %>% as.matrix() %>% as.vector() %>% qqnorm()
  # abline(a = 0, b = 1)
})
# 1 out of 20 fails :)

test_that("DS Detection Residuals are Gaussian for Artificial Fitted Object made of Hard to Detect Species", {
  # simulate a fitted object
  fit <- artificial_runjags(nspecies = 5, nsites = 500, nvisitspersite = 2,
                            v.b.max = -1,
                            ObsFmla = "~ 1",
                            modeltype = "jsodm_lv",
                            nlv = 2)
  
  # crnings(ompute residuals 
  resid_det <- ds_detection_residuals.fit(fit, type = 1)
  shapiro_det_residual <- resid_det %>% dplyr::select(-ModelSite) %>%  as.matrix() %>%  as.vector() %>%
    shapiro.test()
  expect_gt(shapiro_det_residual$p.value, 0.01)
  
  # resid_det %>% dplyr::select(-ModelSite) %>% as.matrix() %>% as.vector() %>% qqnorm()
  # abline(a = 0, b = 1)
})
# failed 0 out of 20 times.

test_that("DS Detection Residuals are Gaussian for Artificial Fitted Object made of Variety of Species", {
  # simulate a fitted object
  fit <- artificial_runjags(nspecies = 5, nsites = 500, nvisitspersite = 5, modeltype = "jsodm_lv", nlv = 2)
  
  # compute residuals 
  resid_det <- ds_detection_residuals.fit(fit, type = 1)
  shapiro_det_residual <- resid_det %>% dplyr::select(-ModelSite) %>%  as.matrix() %>%  as.vector() %>%
    shapiro.test()
  expect_gt(shapiro_det_residual$p.value, 0.01)
  
  # resid_det %>% dplyr::select(-ModelSite) %>% as.matrix() %>% as.vector() %>% qqnorm()
  # abline(a = 0, b = 1)
  # u.b <- bugsvar2matrix(fit$mcmc[[1]][1,  ], "u.b", 1:fit$data$n, 1:fit$data$Vocc)
  # v.b <- bugsvar2matrix(fit$mcmc[[1]][1,  ], "v.b", 1:fit$data$n, 1:fit$data$Vobs)
  # lv.coef <- bugsvar2matrix(fit$mcmc[[1]][1,  ], "lv.coef", 1:fit$data$n, 1:fit$data$nlv)
})
# failed 0 out of 20 times

test_that("DS Occupancy Residuals are Gaussian for Artificial Fitted Object", {
  # simulate a fitted object
  fit <- artificial_runjags(nspecies = 5, nsites = 500, nvisitspersite = 2, modeltype = "jsodm_lv", nlv = 2)
  
  # compute residuals 
  resid_occ <- ds_occupancy_residuals.fit(fit, type = 1)
  shapiro_occ_residual <- resid_occ %>% dplyr::select(-ModelSite) %>%  as.matrix() %>%  as.vector() %>%
    shapiro.test()
  expect_gt(shapiro_occ_residual$p.value, 0.01)
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

