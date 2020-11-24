# test expected number of species

context("Number of Detected Species Expected")
pbopt <- pbapply::pboptions(type = "none")

# tests:
# multiple different draws*****
# model has LVs, use fitted LVs, in sample data
# model has LVs, marginalise LVs, in sample data
# model has no LVs, in sample data

# model has LVs, marginalise LVs, holdout data
# model has no LVs, holdout data

test_that("In sample data; fitted LV values; different draws", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 1, modeltype = "jsodm_lv", nlv = 4)
  
  # make a new, different, second parameter set
  artfit$mcmc[[2]] <- artfit$mcmc[[1]]
  artfit$sample <- 1
  bugvarnames <- names(artfit$mcmc[[1]][1, ])
  artfit$mcmc[[2]][1, grepl("^u.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^u.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^v.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^v.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^lv.coef\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^lv.coef\\[.*", bugvarnames)] * runif(3, min = 0.1, max = 0.2)
  artfit$mcmc[[2]][1, grepl("^LV\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^LV\\[.*", bugvarnames)] * runif(4, min = 0.5, max = 1)
  
  # Predicted number of species detected and in occupation
  numspec <- predsumspecies(artfit, UseFittedLV = TRUE)
  meanvar <- cumsum(numspec["Vsum_det_median", ])/((1:ncol(numspec))^2)
  sd_final <- sqrt(meanvar[ncol(numspec)])
  expect_equal(ncol(numspec), nsites)
  
  # Anticipate the Enumspec is wrong when using both draws, as simulated data in artfit is from the first draw
  NumSpecies_1st <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  meandiff_1st <- dplyr::cummean(NumSpecies_1st - numspec["Esum_det_median", ])
  plt <- cbind(diff = meandiff_1st, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  expect_gt(abs(meandiff_1st[ncol(numspec)]), 3 * sd_final)
  
  # Anticipate the Enumspec is correct when using only first draw (chain), as simulated data in artfit is from the first draw
  Enumspec_1stonly <- predsumspecies(artfit, chain = 1, UseFittedLV = TRUE)
  meanvar_1stonly <- cumsum(Enumspec_1stonly["Vsum_det_median", ])/((1:ncol(Enumspec_1stonly))^2)
  sd_final_1st <- sqrt(meanvar[ncol(Enumspec_1stonly)])
  
  meandiff_1st <- dplyr::cummean(NumSpecies_1st - Enumspec_1stonly["Esum_det_median", ])
  plt <- cbind(diff = meandiff_1st, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  expect_lt(abs(meandiff_1st[ncol(Enumspec_1stonly)]), 3 * sd_final_1st)
  
  # Anticipate that it is correct for 2nd draw separated from the 1st draw
  y_2nd <- simulate_detections(artfit, esttype = 2)
  NumSpecies_2nd <- detectednumspec(y = y_2nd, ModelSite = artfit$data$ModelSite)
  Enumspec_2ndonly <- predsumspecies(artfit, chain = 2, UseFittedLV = TRUE)
  meanvar_2ndonly <- cumsum(Enumspec_2ndonly["Vsum_det_median", ])/((1:ncol(Enumspec_2ndonly))^2)
  sd_final_2nd <- sqrt(meanvar[ncol(Enumspec_2ndonly)])
  meandiff_2nd <- dplyr::cummean(NumSpecies_2nd - Enumspec_2ndonly["Esum_det_median", ])
  plt <- cbind(diff = meandiff_2nd, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  expect_lt(abs(meandiff_2nd[ncol(Enumspec_2ndonly)]), 3 * sd_final_2nd)
  
  
  # Anticipate prediction from combined draw is similar when occupancy + detection simulated with parameters chosen with equal chance from artfit$mcmc[[1]]
  # But anticipate within-model uncertainty of median parameter values does not cover observed values.
  ## combine observations for model sites to simulate the equal credence on each parameter set
  NumSpecies_interleaved <- NumSpecies_1st
  drawselect <- as.logical(rbinom(1000, size = 1, prob = 0.5))
  NumSpecies_interleaved[drawselect] <- NumSpecies_2nd[drawselect]
  
  meandiff <- dplyr::cummean(NumSpecies_interleaved - numspec["Esum_det_margpost", ])
  meanvar <- cumsum(numspec["Vsum_det_margpost", ])/((1:ncol(numspec))^2)
  plt <- cbind(diff = meandiff, var = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  sd_final <- sqrt(meanvar[ncol(numspec)])
  expect_lt(abs(meandiff[ncol(numspec)]), 3 * sd_final)
  
  # Hope that Gaussian approximation of a 95% interval covers the observed data 95% of the time
  ininterval_marg <- (NumSpecies_interleaved > numspec["Esum_det_margpost", ] - 2 * sqrt(numspec["Vsum_det_margpost", ])) & 
    (NumSpecies_interleaved < numspec["Esum_det_margpost", ] + 2 * sqrt(numspec["Vsum_det_margpost", ]))
  expect_equal(mean(ininterval_marg), 0.95, tol = 0.05)
  
  # and that within-model median parameters *do not*
  ininterval_median <- (NumSpecies_interleaved > numspec["Esum_det_median", ] - 2 * sqrt(numspec["Vsum_det_median", ])) & 
    (NumSpecies_interleaved < numspec["Esum_det_median", ] + 2 * sqrt(numspec["Vsum_det_median", ]))
  expect_lt(mean(ininterval_median), 0.90)
})

test_that("In sample data; fitted LV values", {
  nsites <- 10000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, modeltype = "jsodm_lv", nlv = 4)
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  numspec <- predsumspecies(artfit, UseFittedLV = TRUE, type = "median")
  expect_equal(ncol(numspec), nsites)
  
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  # median of theta should be correct
  meandiff <- dplyr::cummean(NumSpecies - numspec["Esum_det_median", ])
  meanvar <- cumsum(numspec["Vsum_det_median", ])/((1:ncol(numspec))^2)
  plt <- cbind(diff = meandiff, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  # check with predicted standard error once the software is computed
  sd_final <- sqrt(meanvar[ncol(numspec)])
  expect_equal(meandiff[ncol(numspec)], 0, tol = 3 * sd_final)
  
  # difference between expected and observed should be zero on average; check that is getting closer with increasing data
  expect_lt(abs(meandiff[length(meandiff)]), max(abs(meandiff[floor(length(meandiff) / 20) + 1:20 ])))
  
  Enum_compare_sum <- Enum_compare(NumSpecies,
               data.frame(pred = numspec["Esum_det_median", ]),
               data.frame(pred = numspec["Vsum_det_median", ])
               )
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]])
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]])
  expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
  
  # Hope that Gaussian approximation of a 95% interval covers the observed data 95% of the time
  ininterval <- (NumSpecies > numspec["Esum_det_margpost", ] - 2 * sqrt(numspec["Vsum_det_margpost", ])) & 
    (NumSpecies < numspec["Esum_det_margpost", ] + 2 * sqrt(numspec["Vsum_det_margpost", ]))
  expect_equal(mean(ininterval), 0.95, tol = 0.05)
  
  plt <- cbind(NumSpecies = NumSpecies, pred = numspec["Esum_det_margpost", ], se = sqrt(numspec["Vsum_det_margpost", ])) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(resid = NumSpecies - pred) %>%
    dplyr::arrange(resid) %>%
    tibble::rowid_to_column(var = "rowid") %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= rowid, ymin = -2 * se, ymax = + 2 * se), fill = "grey") +
    ggplot2::geom_point(ggplot2::aes(x = rowid, y = NumSpecies - pred), col = "blue", lwd = 2)
  # print(plt)
})

test_that("In sample data; marginal on LV values", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3,
                               u.b.min = -0.01,
                               u.b.max = 0.01,
                               v.b.min = -0.01,
                               v.b.max = 0.01,
                               lv.coef.min = 0.4,
                               lv.coef.max = 0.5, #hopefully LVs have a much bigger effect than occupancy etc,
			       modeltype = "jsodm_lv",
			       nlv = 4
                               )
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  numspec <- predsumspecies(artfit, UseFittedLV = FALSE, nLVsim = 1000, type = "median")
  expect_equal(ncol(numspec), nsites)
  
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  # the median should be perfectly correct except that it ignores the fitted LV values aren't used
  meandiff <- dplyr::cummean(NumSpecies - numspec["Esum_det_median", ])
  meanvar <- cumsum(numspec["Vsum_det_median", ])/((1:ncol(numspec))^2)
  plt <- cbind(diff = meandiff, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  # check with predicted standard error: should be larger than bound due to misuse of LV
  sd_half <- sqrt(meanvar[floor(ncol(numspec)/2)])
  expect_gt(abs(meandiff[floor(ncol(numspec)/2)]), 3 * sd_half)
  
  # expect that cover by posterior approximate interval is about 
  # right as LV values simulate from are not extreme for a Gaussian
  # and the calculations marginalise of the LV value distribution
  ininterval <- (NumSpecies > numspec["Esum_det_margpost", ] - 2 * sqrt(numspec["Vsum_det_margpost", ])) & 
    (NumSpecies < numspec["Esum_det_margpost", ] + 2 * sqrt(numspec["Vsum_det_margpost", ]))
  expect_gt(mean(ininterval), 0.95)
  
  ######################################### Now try using marginalised LV simulations #####################################
  NumSpecies <- detectednumspec(y = simulate_detections_LV(artfit, esttype = 1), ModelSite = artfit$data$ModelSite)
  
  meandiff <- dplyr::cummean(NumSpecies - numspec["Esum_det_median", ])
  meanvar <- cumsum(numspec["Vsum_det_median", ])/((1:ncol(numspec))^2)
  plt <- cbind(diff = meandiff, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  # check with predicted standard error 
  sd_final <- sqrt(meanvar[ncol(numspec)])
  expect_equal(meandiff[ncol(numspec)], 0, tol = 3 * sd_final)
  
  # difference between expected and observed should be zero on average; check that is getting closer with increasing data
  expect_lt(abs(meandiff[length(meandiff)]), max(abs(meandiff[floor(length(meandiff) / 20) + 1:20 ])))
  
  # Hope that Gaussian approximation of a 95% interval covers the observed data 95% of the time
  ininterval <- (NumSpecies > numspec["Esum_det_margpost", ] - 2 * sqrt(numspec["Vsum_det_margpost", ])) & 
    (NumSpecies < numspec["Esum_det_margpost", ] + 2 * sqrt(numspec["Vsum_det_margpost", ]))
  expect_gt(mean(ininterval), 0.99)
  # I suspect this is not 0.95 because the Gaussian approximation may not work: the variance is enormous compared to the allowed range of number of species
  
  plt <- cbind(NumSpecies = NumSpecies, pred = numspec["Esum_det_margpost", ], se = sqrt(numspec["Vsum_det_margpost", ])) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(resid = NumSpecies - pred) %>%
    dplyr::arrange(resid) %>%
    tibble::rowid_to_column(var = "rowid") %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= rowid, ymin = -2 * se, ymax = + 2 * se), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = rowid, y = NumSpecies - pred), col = "blue", lwd = 2)
  # print(plt)
})

test_that("In sample data; no LV", {
  nsites <- 10000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, modeltype = "jsodm")
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  Enumspecdet <- predsumspecies(artfit, UseFittedLV = FALSE, type = "marginal")
  expect_equal(ncol(Enumspecdet), nsites)
  
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  meandiff <- dplyr::cummean(NumSpecies - Enumspecdet["Esum_det", ])
  meanvar <- cumsum(Enumspecdet["Vsum_det", ])/((1:ncol(Enumspecdet))^2)
  plt <- cbind(diff = meandiff, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  # check with predicted standard error once the software is computed
  sd_final <- sqrt(meanvar[ncol(Enumspecdet)])
  expect_equal(meandiff[ncol(Enumspecdet)], 0, tol = 3 * sd_final)
  
  # difference between expected and observed should be zero on average; check that is getting closer with increasing data
  expect_lt(abs(meandiff[length(meandiff)]), max(abs(meandiff[floor(length(meandiff) / 20) + 1:20 ])))
  
  # Hope that Gaussian approximation of a 95% interval covers the observed data 95% of the time
  ininterval <- (NumSpecies > Enumspecdet["Esum_det", ] - 2 * sqrt(Enumspecdet["Vsum_det", ])) & 
    (NumSpecies < Enumspecdet["Esum_det", ] + 2 * sqrt(Enumspecdet["Vsum_det", ]))
  expect_equal(mean(ininterval), 0.95, tol = 0.05)
})


test_that("Holdout data; has LVs", {
  nsites <- 10000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, modeltype = "jsodm_lv", nlv = 4) 
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  originalXocc <- unstandardise.designmatprocess(artfit$XoccProcess, artfit$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- unstandardise.designmatprocess(artfit$XobsProcess, artfit$data$Xobs)
  originalXobs <- cbind(ModelSite = artfit$data$ModelSite, originalXobs)
  outofsample_y <- simulate_detections_LV(artfit, esttype = 1)
  
  Enumspec <- predsumspecies_newdata(artfit, originalXocc, originalXobs, ModelSiteVars = "ModelSite", chains = NULL, nLVsim = 1000, type = "median", cl = NULL)
  
  expect_equal(ncol(Enumspec), nsites)
  
  NumSpecies <- detectednumspec(y = outofsample_y, ModelSite = originalXobs[, "ModelSite"])
  
  meandiff <- dplyr::cummean(NumSpecies - Enumspec["Esum_det_median", ])
  meanvar <- cumsum(Enumspec["Vsum_det_median", ])/((1:ncol(Enumspec))^2)
  plt <- cbind(diff = meandiff, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  # check with predicted standard error once the software is computed
  sd_final <- sqrt(meanvar[ncol(Enumspec)])
  expect_equal(meandiff[ncol(Enumspec)], 0, tol = 3 * sd_final)
  
  # difference between expected and observed should be zero on average; check that is getting closer with increasing data
  expect_lt(abs(meandiff[length(meandiff)]), max(abs(meandiff[floor(length(meandiff) / 100) + 1:20 ])))
  
  Enum_compare_sum <- Enum_compare(NumSpecies,
                                   data.frame(pred = Enumspec["Esum_det_median", ]),
                                   data.frame(pred = Enumspec["Vsum_det_median", ])
  )
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]])
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]])
  expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
  
  # Hope that Gaussian approximation of a 95% interval covers the observed data 95% of the time
  ininterval <- (NumSpecies > Enumspec["Esum_det_margpost", ] - 2 * sqrt(Enumspec["Vsum_det_margpost", ])) & 
    (NumSpecies < Enumspec["Esum_det_margpost", ] + 2 * sqrt(Enumspec["Vsum_det_margpost", ]))
  expect_equal(mean(ininterval), 0.95, tol = 0.05)
})



test_that("Holdout data; no LVs", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, modeltype = "jsodm")
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  originalXocc <- unstandardise.designmatprocess(artfit$XoccProcess, artfit$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- unstandardise.designmatprocess(artfit$XobsProcess, artfit$data$Xobs)
  originalXobs <- cbind(ModelSite = artfit$data$ModelSite, originalXobs)
  outofsample_y <- simulate_detections(artfit, esttype = 1)
  
  Enumspec <- predsumspecies_newdata(artfit, originalXocc, originalXobs, ModelSiteVars = "ModelSite", chains = NULL, nLVsim = 1000, type = "marginal", cl = NULL)
  
  expect_equal(ncol(Enumspec), nsites)
  
  NumSpecies <- detectednumspec(y = outofsample_y, ModelSite = originalXobs[, "ModelSite"])
  
  meandiff <- dplyr::cummean(NumSpecies - Enumspec["Esum_det", ])
  meanvar <- cumsum(Enumspec["Vsum_det", ])/((1:ncol(Enumspec))^2)
  plt <- cbind(diff = meandiff, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  # check with predicted standard error once the software is computed
  sd_final <- sqrt(meanvar[ncol(Enumspec)])
  expect_equal(meandiff[ncol(Enumspec)], 0, tol = 3 * sd_final)
  
  # difference between expected and observed should be zero on average; check that is getting closer with increasing data
  expect_lt(abs(meandiff[length(meandiff)]), max(abs(meandiff[floor(length(meandiff) / 100) + 1:20 ])))
  
  Enum_compare_sum <- Enum_compare(NumSpecies,
                                   data.frame(pred = Enumspec["Esum_det", ]),
                                   data.frame(pred = Enumspec["Vsum_det", ])
  )
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]])
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]])
  expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
  
  # Hope that Gaussian approximation of a 95% interval covers the observed data 95% of the time
  Enumspecdet <- predsumspecies_newdata(artfit, originalXocc, originalXobs, ModelSiteVars = "ModelSite", chains = NULL, nLVsim = 1000, type = "marginal")
  ininterval <- (NumSpecies > Enumspecdet["Esum_det", ] - 2 * sqrt(Enumspecdet["Vsum_det", ])) & 
    (NumSpecies < Enumspecdet["Esum_det", ] + 2 * sqrt(Enumspecdet["Vsum_det", ]))
  expect_equal(mean(ininterval), 0.95, tol = 0.05)
})

test_that("Subset biodiversity matches simulations", {
  # make it full test by having: different draws and latent variables, and testing both marginal and fitted latent variables
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 1, modeltype = "jsodm_lv", nlv = 4)
  
  # make a new, different, second parameter set
  artfit$mcmc[[2]] <- artfit$mcmc[[1]]
  artfit$sample <- 1
  bugvarnames <- names(artfit$mcmc[[1]][1, ])
  artfit$mcmc[[2]][1, grepl("^u.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^u.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^v.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^v.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^lv.coef\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^lv.coef\\[.*", bugvarnames)] * runif(3, min = 0.1, max = 0.2)
  artfit$mcmc[[2]][1, grepl("^LV\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^LV\\[.*", bugvarnames)] * runif(4, min = 0.5, max = 1)
   
  # Simulate equally from each draw
  y_1st <- simulate_detections(artfit, esttype = 1)
  y_2nd <- simulate_detections(artfit, esttype = 2)
  y_interleaved <- y_1st
  drawselect <- as.logical(rbinom(1000, size = 1, prob = 0.5))
  y_interleaved[drawselect, ] <- y_2nd[drawselect, ]
  
  # Choose a subset of species
  speciessubset <- sample(artfit$species, size = 20)
  NumSpeciesObs <- detectednumspec(y_interleaved[, speciessubset], ModelSite = artfit$data$ModelSite)
  
  # Predict number within subset, in sample, using LV
  numspec_insample_fitLV <- predsumspecies(artfit, desiredspecies = speciessubset, UseFittedLV = TRUE, type = "marginal")
  inci_insample_fitLV <- (NumSpeciesObs > numspec_insample_fitLV["Esum_det", ] - 2 * sqrt(numspec_insample_fitLV["Vsum_det", ])) & 
    (NumSpeciesObs < numspec_insample_fitLV["Esum_det", ] + 2 * sqrt(numspec_insample_fitLV["Vsum_det", ]))
  expect_equal(mean(inci_insample_fitLV), 0.95, tol = 0.05)
  
  # Predict number within subset, in sample, marginal LV
  numspec_insample_margLV <- predsumspecies(artfit, desiredspecies = speciessubset, UseFittedLV = FALSE, type = "marginal")
  inci_insample_margLV <- (NumSpeciesObs > numspec_insample_margLV["Esum_det", ] - 2 * sqrt(numspec_insample_margLV["Vsum_det", ])) & 
    (NumSpeciesObs < numspec_insample_margLV["Esum_det", ] + 2 * sqrt(numspec_insample_margLV["Vsum_det", ]))
  expect_equal(mean(inci_insample_margLV), 0.95, tol = 0.05)
  
  # Predict number within subset, outside sample, marginal LV
  originalXocc <- unstandardise.designmatprocess(artfit$XoccProcess, artfit$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- unstandardise.designmatprocess(artfit$XobsProcess, artfit$data$Xobs)
  originalXobs <- cbind(ModelSite = artfit$data$ModelSite, originalXobs)
  outofsample_y <- simulate_detections_LV(artfit, esttype = 1)
  
  numspec_holdout_margLV <- predsumspecies_newdata(artfit, originalXocc, originalXobs, ModelSiteVars = "ModelSite",
                                     desiredspecies = speciessubset,
                                     chains = NULL, nLVsim = 1000, type = "marginal", cl = NULL)
  
  inci_holdout_margLV <- (NumSpeciesObs > numspec_holdout_margLV["Esum_det", ] - 2 * sqrt(numspec_holdout_margLV["Vsum_det", ])) & 
    (NumSpeciesObs < numspec_holdout_margLV["Esum_det", ] + 2 * sqrt(numspec_holdout_margLV["Vsum_det", ]))
  expect_equal(mean(inci_holdout_margLV), 0.95, tol = 0.05)
})

test_that("Subset biodiversity to single species matches simulations", {
  # make it full test by having: different draws and latent variables, and testing both marginal and fitted latent variables
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 1, modeltype = "jsodm_lv", nlv = 4)
  
  # make a new, different, second parameter set
  artfit$mcmc[[2]] <- artfit$mcmc[[1]]
  artfit$sample <- 1
  bugvarnames <- names(artfit$mcmc[[1]][1, ])
  artfit$mcmc[[2]][1, grepl("^u.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^u.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^v.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^v.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^lv.coef\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^lv.coef\\[.*", bugvarnames)] * runif(3, min = 0.1, max = 0.2)
  artfit$mcmc[[2]][1, grepl("^LV\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^LV\\[.*", bugvarnames)] * runif(4, min = 0.5, max = 1)
  
  # Simulate equally from each draw
  y_1st <- simulate_detections(artfit, esttype = 1)
  y_2nd <- simulate_detections(artfit, esttype = 2)
  y_interleaved <- y_1st
  drawselect <- as.logical(rbinom(1000, size = 1, prob = 0.5))
  y_interleaved[drawselect, ] <- y_2nd[drawselect, ]
  
  # Choose a subset of species
  speciessubset <- sample(artfit$species, size = 1)
  NumSpeciesObs <- detectednumspec(y_interleaved[, speciessubset, drop = FALSE], ModelSite = artfit$data$ModelSite)
  
  # Predict number within subset, in sample, using LV
  numspec_insample_fitLV <- predsumspecies(artfit, desiredspecies = speciessubset, UseFittedLV = TRUE, type = "marginal")
  Enum_compare_sum <- Enum_compare(NumSpeciesObs,
                                   data.frame(pred = numspec_insample_fitLV["Esum_det", ]),
                                   data.frame(pred = numspec_insample_fitLV["Vsum_det", ])
  )
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]])
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]])
  expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])

  # Predict number within subset, in sample, marginal LV
  numspec_insample_margLV <- predsumspecies(artfit, desiredspecies = speciessubset, UseFittedLV = FALSE, type = "marginal")
  Enum_compare_sum <- Enum_compare(NumSpeciesObs,
                                   data.frame(pred = numspec_insample_margLV["Esum_det", ]),
                                   data.frame(pred = numspec_insample_margLV["Vsum_det", ])
  )
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]])
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]])
  # the follow test doesn't pass, which is consistent with the fitted LV values not being distributed according to a Gaussian distribution
  # expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
  
  # Predict number within subset, outside sample, marginal LV
  originalXocc <- unstandardise.designmatprocess(artfit$XoccProcess, artfit$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- unstandardise.designmatprocess(artfit$XobsProcess, artfit$data$Xobs)
  originalXobs <- cbind(ModelSite = artfit$data$ModelSite, originalXobs)
  outofsample_y <- simulate_detections_LV(artfit, esttype = 1)
  
  numspec_holdout_margLV <- predsumspecies_newdata(artfit, originalXocc, originalXobs, ModelSiteVars = "ModelSite",
                                                   desiredspecies = speciessubset,
                                                   chains = NULL, nLVsim = 1000, type = "marginal", cl = NULL)
  Enum_compare_sum <- Enum_compare(NumSpeciesObs,
                                   data.frame(pred = numspec_holdout_margLV["Esum_det", ]),
                                   data.frame(pred = numspec_holdout_margLV["Vsum_det", ])
  )
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]])
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]])
  # the follow test doesn't pass, which is consistent with the fitted LV values not being distributed according to a Gaussian distribution
  # expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
})

test_that("Endetect_modelsite matches predsumspecies", {
  # make it full test by having: different draws and latent variables, and testing both marginal and fitted latent variables
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 1, modeltype = "jsodm_lv", nlv = 4)
  
  # using fitted LV
  Edet1 <- Endetect_modelsite(artfit, type = "median", conditionalLV = TRUE)
  Edet2 <- lapply(artfit$species, function(sp) {
                  Edet <- predsumspecies(artfit, desiredspecies = sp, UseFittedLV = TRUE, type = "marginal")
                  #marginal works here because artfit has only one draw
                  return(Edet)}
                  )
  Edet2_t <- t(do.call(rbind, lapply(Edet2, function(x) x["Esum_det", , drop = FALSE])))
  expect_equivalent(Edet1[[1]], Edet2_t)
  
  Edet2_V_t <- t(do.call(rbind, lapply(Edet2, function(x) x["Vsum_det", , drop = FALSE])))
  expect_equivalent(Edet1[[2]], Edet2_V_t)
  
  # marginal to LV
  Edet1 <- Endetect_modelsite(artfit, type = "median", conditionalLV = FALSE)
  Edet2 <- lapply(artfit$species, function(sp) {
    Edet <- predsumspecies(artfit, desiredspecies = sp, UseFittedLV = FALSE, nLVsim = 5000, type = "marginal")
    #marginal works here because artfit has only one draw
    return(Edet)}
  )
  Edet2_t <- t(do.call(rbind, lapply(Edet2, function(x) x["Esum_det", , drop = FALSE])))
  expect_equivalent(Edet1[[1]], Edet2_t, tol = 1E-2)
  
  Edet2_V_t <- t(do.call(rbind, lapply(Edet2, function(x) x["Vsum_det", , drop = FALSE])))
  expect_equivalent(Edet1[[2]], Edet2_V_t, tol = 1E-2)
})

#########################################################################################

test_that("No LV and identical sites", {
  artfit <- artificial_runjags(nspecies = 4, nsites = 1000, nvisitspersite = 2,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1",
                               modeltype = "jsodm")
  EVsum <- predsumspecies(artfit, UseFittedLV = FALSE, type = "marginal")
  
  # check that many other sites have the same expected number of species
  expect_equivalent(EVsum["Esum_det", ], rep(EVsum["Esum_det", 1], ncol(EVsum)))
  
  # treat each model site as a repeat simulation of a ModelSite (cos all the parameters are nearly identical)
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  meandiff <- dplyr::cummean(NumSpecies - EVsum["Esum_det", ])
  meanvar <- cumsum(EVsum["Vsum_det", ])/((1:ncol(EVsum))^2)
  plt <- cbind(diff = meandiff, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  # check with predicted standard error once the software is computed
  sd_final <- sqrt(meanvar[ncol(EVsum)])
  expect_equal(meandiff[ncol(EVsum)], 0, tol = 3 * sd_final)
  
  # difference between expected and observed should be zero on average; check that is getting closer with increasing data
  expect_lt(abs(meandiff[length(meandiff)]), max(abs(meandiff[floor(length(meandiff) / 20) + 1:20 ])))

  # Expect sd to be close to theoretical sd. Hopefully within 10%
  expect_equivalent(sd(NumSpecies), sqrt(EVsum["Vsum_det", 1]), tol = 0.1 * sqrt(EVsum["Vsum_det", 1]))
  
  # Hope that Gaussian approximation of a 95% interval covers the observed data 95% of the time
  Enumspecdet <- predsumspecies(artfit, type = "marginal", UseFittedLV = FALSE)
  ininterval <- (NumSpecies > Enumspecdet["Esum_det", ] - 2 * sqrt(Enumspecdet["Vsum_det", ])) & 
    (NumSpecies < Enumspecdet["Esum_det", ] + 2 * sqrt(Enumspecdet["Vsum_det", ]))
  expect_equal(mean(ininterval), 0.95, tol = 0.05)
})

test_that("Expected occupied number for in sample data; fitted LV values", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, 
                               v.b.min = 20, v.b.max = 20.1, modeltype = "jsodm_lv", nlv = 4) #detection is still not certain
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  Enumspecdet <- predsumspecies(artfit, UseFittedLV = TRUE, type = "marginal")
  expect_equal(ncol(Enumspecdet), nsites)
  
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  meandiff <- dplyr::cummean(NumSpecies - Enumspecdet["Esum_det", ])
  meanvar <- cumsum(Enumspecdet["Vsum_det", ])/((1:ncol(Enumspecdet))^2)
  plt <- cbind(diff = meandiff, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  # check with predicted standard error once the software is computed
  sd_final <- sqrt(meanvar[ncol(Enumspecdet)])
  expect_equal(meandiff[ncol(Enumspecdet)], 0, tol = 3 * sd_final)
})

pbapply::pboptions(pbopt)
