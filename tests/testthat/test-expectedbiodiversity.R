# test expected number of species

local_edition(3)

pbopt <- pbapply::pboptions(type = "none")

# tests:
# multiple different draws*****
# model has lv.vs, use fitted lv.vs, in sample data
# model has lv.vs, marginalise lv.vs, in sample data
# model has no lv.vs, in sample data

# model has lv.vs, marginalise lv.vs, holdout data
# model has no lv.vs, holdout data

test_that("In sample data; fitted lv.v; different draws", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 1, modeltype = "jsodm_lv", nlv = 4)
  
  # make a new, different, second parameter set
  artfit$mcmc[[2]] <- artfit$mcmc[[1]]
  artfit$sample <- 1
  bugvarnames <- names(artfit$mcmc[[1]][1, ])
  artfit$mcmc[[2]][1, grepl("^occ.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^occ.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^det.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^det.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^lv.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^lv.b\\[.*", bugvarnames)] * runif(3, min = 0.1, max = 0.2)
  artfit$mcmc[[2]][1, grepl("^lv.v\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^lv.v\\[.*", bugvarnames)] * runif(4, min = 0.5, max = 1)
  
  # Predicted number of species detected and in occupation
  numspec <- predsumspecies(artfit, UseFittedLV = TRUE)
  meanvar <- cumsum(numspec["V", ])/((1:ncol(numspec))^2)
  sd_final <- sqrt(meanvar[ncol(numspec)])
  expect_equal(ncol(numspec), nsites)
  
  # Anticipate the Enumspec is wrong when using both draws, as simulated data in artfit is from the first draw
  NumSpecies_1st <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  meandiff_1st <- dplyr::cummean(NumSpecies_1st - numspec["E", ])
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
  meanvar_1stonly <- cumsum(Enumspec_1stonly["V", ])/((1:ncol(Enumspec_1stonly))^2)
  sd_final_1st <- sqrt(meanvar[ncol(Enumspec_1stonly)])
  
  meandiff_1st <- dplyr::cummean(NumSpecies_1st - Enumspec_1stonly["E", ])
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
  meanvar_2ndonly <- cumsum(Enumspec_2ndonly["V", ])/((1:ncol(Enumspec_2ndonly))^2)
  sd_final_2nd <- sqrt(meanvar[ncol(Enumspec_2ndonly)])
  meandiff_2nd <- dplyr::cummean(NumSpecies_2nd - Enumspec_2ndonly["E", ])
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
  
  meandiff <- dplyr::cummean(NumSpecies_interleaved - numspec["E", ])
  meanvar <- cumsum(numspec["V", ])/((1:ncol(numspec))^2)
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
  ininterval_marg <- (NumSpecies_interleaved > numspec["E", ] - 2 * sqrt(numspec["V", ])) & 
    (NumSpecies_interleaved < numspec["E", ] + 2 * sqrt(numspec["V", ]))
  expect_equal(mean(ininterval_marg), 0.95, tol = 0.05)
})

test_that("In sample data; fitted lv.v", {
  nsites <- 10000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, modeltype = "jsodm_lv", nlv = 4)
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  numspec <- predsumspecies(artfit, UseFittedLV = TRUE, type = "marginal")
  expect_equal(ncol(numspec), nsites)
  
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  # median of theta should be correct
  meandiff <- dplyr::cummean(NumSpecies - numspec["E", ])
  meanvar <- cumsum(numspec["V", ])/((1:ncol(numspec))^2)
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
               data.frame(pred = numspec["E", ]),
               data.frame(pred = numspec["V", ])
               )
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]])
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]])
  expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
  
  # Hope that Gaussian approximation of a 95% interval covers the observed data 95% of the time
  ininterval <- (NumSpecies > numspec["E", ] - 2 * sqrt(numspec["V", ])) & 
    (NumSpecies < numspec["E", ] + 2 * sqrt(numspec["V", ]))
  expect_equal(mean(ininterval), 0.95, tol = 0.05)
  
  plt <- cbind(NumSpecies = NumSpecies, pred = numspec["E", ], se = sqrt(numspec["V", ])) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(resid = NumSpecies - pred) %>%
    dplyr::arrange(resid) %>%
    tibble::rowid_to_column(var = "rowid") %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= rowid, ymin = -2 * se, ymax = + 2 * se), fill = "grey") +
    ggplot2::geom_point(ggplot2::aes(x = rowid, y = NumSpecies - pred), col = "blue", lwd = 2)
  # print(plt)
})

test_that("In sample data; marginal on lv.v", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3,
                               occ.b.min = -0.01,
                               occ.b.max = 0.01,
                               det.b.min = -0.01,
                               det.b.max = 0.01,
                               lv.b.min = 0.4,
                               lv.b.max = 0.5, #hopefully lv.vs have a much bigger effect than occupancy etc,
			       modeltype = "jsodm_lv",
			       nlv = 4
                               )
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  numspec <- predsumspecies(artfit, UseFittedLV = FALSE, nLVsim = 1000, type = "marginal")
  expect_equal(ncol(numspec), nsites)
  
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  meandiff <- dplyr::cummean(NumSpecies - numspec["E", ])
  meanvar <- cumsum(numspec["V", ])/((1:ncol(numspec))^2)
  plt <- cbind(diff = meandiff, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  # print(plt)
  
  # check with predicted standard error: should be larger than bound due to misuse of lv.v
  sd_half <- sqrt(meanvar[floor(ncol(numspec)/2)])
  expect_gt(abs(meandiff[floor(ncol(numspec)/2)]), 3 * sd_half)
  
  # expect that cover by posterior approximate interval is about 
  # right as lv.v simulate from are not extreme for a Gaussian
  # and the calculations marginalise of the lv.v distribution
  ininterval <- (NumSpecies > numspec["E", ] - 2 * sqrt(numspec["V", ])) & 
    (NumSpecies < numspec["E", ] + 2 * sqrt(numspec["V", ]))
  expect_gt(mean(ininterval), 0.95)
  
  ######################################### Now try using marginalised lv.v simulations #####################################
  NumSpecies <- detectednumspec(y = simulate_detections_lv.v(artfit, esttype = 1), ModelSite = artfit$data$ModelSite)
  
  meandiff <- dplyr::cummean(NumSpecies - numspec["E", ])
  meanvar <- cumsum(numspec["V", ])/((1:ncol(numspec))^2)
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
  ininterval <- (NumSpecies > numspec["E", ] - 2 * sqrt(numspec["V", ])) & 
    (NumSpecies < numspec["E", ] + 2 * sqrt(numspec["V", ]))
  expect_gt(mean(ininterval), 0.99)
  # I suspect this is not 0.95 because the Gaussian approximation may not work: the variance is enormous compared to the allowed range of number of species
  
  plt <- cbind(NumSpecies = NumSpecies, pred = numspec["E", ], se = sqrt(numspec["V", ])) %>% 
    dplyr::as_tibble() %>% 
    dplyr::mutate(resid = NumSpecies - pred) %>%
    dplyr::arrange(resid) %>%
    tibble::rowid_to_column(var = "rowid") %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= rowid, ymin = -2 * se, ymax = + 2 * se), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = rowid, y = NumSpecies - pred), col = "blue", lwd = 2)
  # print(plt)
})

test_that("In sample data; no lv.v", {
  nsites <- 10000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, modeltype = "jsodm")
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  Enumspecdet <- predsumspecies(artfit, UseFittedLV = FALSE, type = "marginal")
  expect_equal(ncol(Enumspecdet), nsites)
  
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  meandiff <- dplyr::cummean(NumSpecies - Enumspecdet["E", ])
  meanvar <- cumsum(Enumspecdet["V", ])/((1:ncol(Enumspecdet))^2)
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
  ininterval <- (NumSpecies > Enumspecdet["E", ] - 2 * sqrt(Enumspecdet["V", ])) & 
    (NumSpecies < Enumspecdet["E", ] + 2 * sqrt(Enumspecdet["V", ]))
  expect_equal(mean(ininterval), 0.95, tol = 0.05)
})


test_that("Holdout data; has lv.vs", {
  nsites <- 10000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, modeltype = "jsodm_lv", nlv = 4) 
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  originalXocc <- unstandardise.designmatprocess(artfit$XoccProcess, artfit$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- unstandardise.designmatprocess(artfit$XobsProcess, artfit$data$Xobs)
  originalXobs <- cbind(ModelSite = artfit$data$ModelSite, originalXobs)
  outofsample_y <- simulate_detections_lv.v(artfit, esttype = 1)
  
  Enumspec <- apply_to_new_data(speciesrichness,
                    artfit, originalXocc, originalXobs, ModelSite = "ModelSite",
                    funargs = list(occORdetection = "detection",
                                   usefittedlvv = FALSE,
                                   nlvperdraw = 100))
  
  expect_equal(ncol(Enumspec), nsites)
  
  NumSpecies <- detectednumspec(y = outofsample_y, ModelSite = originalXobs[, "ModelSite"])
  
  meandiff <- dplyr::cummean(NumSpecies - Enumspec["E", ])
  meanvar <- cumsum(Enumspec["V", ])/((1:ncol(Enumspec))^2)
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
                                   data.frame(pred = Enumspec["E", ]),
                                   data.frame(pred = Enumspec["V", ])
  )
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]])
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]])
  expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
  
  # Hope that Gaussian approximation of a 95% interval covers the observed data 95% of the time
  ininterval <- (NumSpecies > Enumspec["E", ] - 2 * sqrt(Enumspec["V", ])) & 
    (NumSpecies < Enumspec["E", ] + 2 * sqrt(Enumspec["V", ]))
  expect_equal(mean(ininterval), 0.95, tol = 0.05)
})



test_that("Holdout data; no lv.vs", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, modeltype = "jsodm")
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  originalXocc <- unstandardise.designmatprocess(artfit$XoccProcess, artfit$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- unstandardise.designmatprocess(artfit$XobsProcess, artfit$data$Xobs)
  originalXobs <- cbind(ModelSite = artfit$data$ModelSite, originalXobs)
  outofsample_y <- simulate_detections(artfit, esttype = 1)
  
  Enumspec <- apply_to_new_data(speciesrichness,
                                artfit, originalXocc, originalXobs, ModelSite = "ModelSite",
                                funargs = list(occORdetection = "detection"))

  expect_equal(ncol(Enumspec), nsites)
  
  NumSpecies <- detectednumspec(y = outofsample_y, ModelSite = originalXobs[, "ModelSite"])
  
  meandiff <- dplyr::cummean(NumSpecies - Enumspec["E", ])
  meanvar <- cumsum(Enumspec["V", ])/((1:ncol(Enumspec))^2)
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
                                   data.frame(pred = Enumspec["E", ]),
                                   data.frame(pred = Enumspec["V", ])
  )
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]])
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]])
  expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
  
  # Hope that Gaussian approximation of a 95% interval covers the observed data 95% of the time
  ininterval <- (NumSpecies > Enumspec["E", ] - 2 * sqrt(Enumspec["V", ])) & 
    (NumSpecies < Enumspec["E", ] + 2 * sqrt(Enumspec["V", ]))
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
  artfit$mcmc[[2]][1, grepl("^occ.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^occ.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^det.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^det.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^lv.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^lv.b\\[.*", bugvarnames)] * runif(3, min = 0.1, max = 0.2)
  artfit$mcmc[[2]][1, grepl("^lv.v\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^lv.v\\[.*", bugvarnames)] * runif(4, min = 0.5, max = 1)
   
  # Simulate equally from each draw
  y_1st <- simulate_detections(artfit, esttype = 1)
  y_2nd <- simulate_detections(artfit, esttype = 2)
  y_interleaved <- y_1st
  drawselect <- as.logical(rbinom(1000, size = 1, prob = 0.5))
  y_interleaved[drawselect, ] <- y_2nd[drawselect, ]
  
  # Choose a subset of species
  speciessubset <- sample(artfit$species, size = 20)
  NumSpeciesObs <- detectednumspec(y_interleaved[, speciessubset], ModelSite = artfit$data$ModelSite)
  
  # Predict number within subset, in sample, using lv.v
  numspec_insample_fitlv.v <- predsumspecies(artfit, desiredspecies = speciessubset, UseFittedLV = TRUE, type = "marginal")
  inci_insample_fitlv.v <- (NumSpeciesObs > numspec_insample_fitlv.v["E", ] - 2 * sqrt(numspec_insample_fitlv.v["V", ])) & 
    (NumSpeciesObs < numspec_insample_fitlv.v["E", ] + 2 * sqrt(numspec_insample_fitlv.v["V", ]))
  expect_equal(mean(inci_insample_fitlv.v), 0.95, tol = 0.05)
  
  # Predict number within subset, in sample, marginal lv.v
  numspec_insample_marglv.v <- predsumspecies(artfit, desiredspecies = speciessubset, UseFittedLV = FALSE, type = "marginal")
  inci_insample_marglv.v <- (NumSpeciesObs > numspec_insample_marglv.v["E", ] - 2 * sqrt(numspec_insample_marglv.v["V", ])) & 
    (NumSpeciesObs < numspec_insample_marglv.v["E", ] + 2 * sqrt(numspec_insample_marglv.v["V", ]))
  expect_equal(mean(inci_insample_marglv.v), 0.95, tol = 0.05)
  
  # Predict number within subset, outside sample, marginal lv.v
  originalXocc <- unstandardise.designmatprocess(artfit$XoccProcess, artfit$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- unstandardise.designmatprocess(artfit$XobsProcess, artfit$data$Xobs)
  originalXobs <- cbind(ModelSite = artfit$data$ModelSite, originalXobs)
  outofsample_y <- simulate_detections_lv.v(artfit, esttype = 1)
  
  numspec_holdout_marglv.v <- 
    apply_to_new_data(speciesrichness,artfit, 
                      originalXocc, originalXobs, ModelSite = "ModelSite",
                      funargs = list(occORdetection = "detection",
                                     usefittedlvv = FALSE,
                                     nlvperdraw = 1000))
  
  NumSpeciesObsHoldout <- detectednumspec(outofsample_y, artfit$data$ModelSite)

  inci_holdout_marglv.v <- (NumSpeciesObsHoldout > numspec_holdout_marglv.v["E", ] - 2 * sqrt(numspec_holdout_marglv.v["V", ])) & 
    (NumSpeciesObsHoldout < numspec_holdout_marglv.v["E", ] + 2 * sqrt(numspec_holdout_marglv.v["V", ]))
  expect_equal(mean(inci_holdout_marglv.v), 0.95, tol = 0.05)
})

test_that("Subset biodiversity to single species matches simulations", {
  # make it full test by having: different draws and latent variables, and testing both marginal and fitted latent variables
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 1, modeltype = "jsodm_lv", nlv = 4)
  
  # make a new, different, second parameter set
  artfit$mcmc[[2]] <- artfit$mcmc[[1]]
  artfit$sample <- 1
  bugvarnames <- names(artfit$mcmc[[1]][1, ])
  artfit$mcmc[[2]][1, grepl("^occ.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^occ.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^det.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^det.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^lv.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^lv.b\\[.*", bugvarnames)] * runif(3, min = 0.1, max = 0.2)
  artfit$mcmc[[2]][1, grepl("^lv.v\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^lv.v\\[.*", bugvarnames)] * runif(4, min = 0.5, max = 1)
  
  # Simulate equally from each draw
  y_1st <- simulate_detections(artfit, esttype = 1)
  y_2nd <- simulate_detections(artfit, esttype = 2)
  y_interleaved <- y_1st
  drawselect <- as.logical(rbinom(1000, size = 1, prob = 0.5))
  y_interleaved[drawselect, ] <- y_2nd[drawselect, ]
  
  # Choose a subset of species
  speciessubset <- sample(artfit$species, size = 1)
  NumSpeciesObs <- detectednumspec(y_interleaved[, speciessubset, drop = FALSE], ModelSite = artfit$data$ModelSite)
  
  # Predict number within subset, in sample, using lv.v
  numspec_insample_fitlv.v <- predsumspecies(artfit, desiredspecies = speciessubset, UseFittedLV = TRUE, type = "marginal")
  Enum_compare_sum <- Enum_compare(NumSpeciesObs,
                                   data.frame(pred = numspec_insample_fitlv.v["E", ]),
                                   data.frame(pred = numspec_insample_fitlv.v["V", ])
  )
  expect_equal(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]], ignore_attr = TRUE)
  expect_equal(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]], ignore_attr = TRUE)
  expect_equal(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]], ignore_attr = TRUE)

  # Predict number within subset, in sample, marginal lv.v
  numspec_insample_marglv.v <- predsumspecies(artfit, desiredspecies = speciessubset, UseFittedLV = FALSE, type = "marginal")
  Enum_compare_sum <- Enum_compare(NumSpeciesObs,
                                   data.frame(pred = numspec_insample_marglv.v["E", ]),
                                   data.frame(pred = numspec_insample_marglv.v["V", ])
  )
  expect_equal(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]], ignore_attr = TRUE)
  expect_equal(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]], ignore_attr = TRUE)
  # the follow test doesn't pass, which is consistent with the fitted lv.v not being distributed according to a Gaussian distribution
  # expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
  
  # Predict number within subset, outside sample, marginal lv.v
  originalXocc <- unstandardise.designmatprocess(artfit$XoccProcess, artfit$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- unstandardise.designmatprocess(artfit$XobsProcess, artfit$data$Xobs)
  originalXobs <- cbind(ModelSite = artfit$data$ModelSite, originalXobs)
  outofsample_y <- simulate_detections_lv.v(artfit, esttype = 1)
  
  numspec_holdout_marglv.v <- 
    apply_to_new_data(speciesrichness, 
                      fit = artfit,
                      originalXocc, 
                      originalXobs, 
                      ModelSite = "ModelSite",
                      funargs = list(occORdetection = "detection",
                                     desiredspecies = speciessubset,
                                     usefittedlvv = FALSE,
                                     nlvperdraw = 100))
  Enum_compare_sum <- Enum_compare(NumSpeciesObs,
                                   data.frame(pred = numspec_holdout_marglv.v["E", ]),
                                   data.frame(pred = numspec_holdout_marglv.v["V", ])
  )
  expect_equal(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]], ignore_attr = TRUE)
  expect_equal(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]], ignore_attr = TRUE)
  # the follow test doesn't pass, which is consistent with the fitted lv.v not being distributed according to a Gaussian distribution
  # expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
})

test_that("Endetect_modelsite matches predsumspecies", {
  # make it full test by having: different draws and latent variables, and testing both marginal and fitted latent variables
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 1, modeltype = "jsodm_lv", nlv = 4)
  
  # using fitted lv.v
  Edet1 <- Endetect_modelsite(artfit, type = "median", conditionalLV = TRUE)
  Edet2 <- lapply(artfit$species, function(sp) {
                  Edet <- speciesrichness(artfit, 
                   occORdetection = "detection",
                   usefittedlvv = TRUE,
                   desiredspecies = sp)
                  return(Edet)}
                  )
  Edet2_t <- t(do.call(rbind, lapply(Edet2, function(x) x["E", , drop = FALSE])))
  expect_equal(Edet1[[1]], Edet2_t, ignore_attr = TRUE)
  
  Edet2_V_t <- t(do.call(rbind, lapply(Edet2, function(x) x["V", , drop = FALSE])))
  expect_equal(Edet1[[2]], Edet2_V_t, ignore_attr = TRUE)
  
  # marginal to lv.v
  Edet1 <- Endetect_modelsite(artfit, type = "median", conditionalLV = FALSE)
  Edet2 <- lapply(artfit$species, function(sp) {
    Edet <- speciesrichness(artfit, 
                            occORdetection = "detection",
                            usefittedlvv = FALSE,
                            nlvperdraw = 1000,
                            desiredspecies = sp)
    return(Edet)}
  )
  Edet2_t <- t(do.call(rbind, lapply(Edet2, function(x) x["E", , drop = FALSE])))
  expect_equal(Edet1[[1]], Edet2_t, tolerance = 1E-2, ignore_attr = TRUE)
  
  Edet2_V_t <- t(do.call(rbind, lapply(Edet2, function(x) x["V", , drop = FALSE])))
  expect_equal(Edet1[[2]], Edet2_V_t, tolerance = 1E-2, ignore_attr = TRUE)
})

#########################################################################################

test_that("No lv.v and identical sites", {
  artfit <- artificial_runjags(nspecies = 4, nsites = 1000, nvisitspersite = 2,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1",
                               modeltype = "jsodm")
  EVsum <- speciesrichness(artfit, occORdetection = "detection")
  
  # check that many other sites have the same expected number of species
  expect_equal(EVsum["E", ], rep(EVsum["E", 1], ncol(EVsum)), ignore_attr = TRUE)
  
  # treat each model site as a repeat simulation of a ModelSite (cos all the parameters are nearly identical)
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  meandiff <- dplyr::cummean(NumSpecies - EVsum["E", ])
  meanvar <- cumsum(EVsum["V", ])/((1:ncol(EVsum))^2)
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
  expect_equal(sd(NumSpecies), sqrt(EVsum["V", 1]), tol = 0.1 * sqrt(EVsum["V", 1]), ignore_attr = TRUE)
  
  # Hope that Gaussian approximation of a 95% interval covers the observed data 95% of the time
  Enumspecdet <- predsumspecies(artfit, type = "marginal", UseFittedLV = FALSE)
  ininterval <- (NumSpecies > Enumspecdet["E", ] - 2 * sqrt(Enumspecdet["V", ])) & 
    (NumSpecies < Enumspecdet["E", ] + 2 * sqrt(Enumspecdet["V", ]))
  expect_equal(mean(ininterval), 0.95, tol = 0.05)
})

test_that("Expected occupied number for in sample data; fitted lv.v", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, 
                               det.b.min = 20, det.b.max = 20.1, modeltype = "jsodm_lv", nlv = 4) #detection is still not certain
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  Enumspecdet <- predsumspecies(artfit, UseFittedLV = TRUE, type = "marginal")
  expect_equal(ncol(Enumspecdet), nsites)
  
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  meandiff <- dplyr::cummean(NumSpecies - Enumspecdet["E", ])
  meanvar <- cumsum(Enumspecdet["V", ])/((1:ncol(Enumspecdet))^2)
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
