# test expected number of species

context("Number of Detected Species Expected")

# tests:
# multiple different draws*****
# model has LVs, use fitted LVs, in sample data
# model has LVs, marginalise LVs, in sample data
# model has no LVs, in sample data

# model has LVs, marginalise LVs, holdout data
# model has no LVs, holdout data

test_that("In sample data; fitted LV values; different draws", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 1, nlv = 4)
  
  # make a new, different, second parameter set
  artfit$mcmc[[2]] <- artfit$mcmc[[1]]
  artfit$sample <- 1
  bugvarnames <- names(artfit$mcmc[[1]][1, ])
  artfit$mcmc[[2]][1, grepl("^u.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^u.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^v.b\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^v.b\\[.*", bugvarnames)] * runif(3, min = 5, max = 10)
  artfit$mcmc[[2]][1, grepl("^lv.coef\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^lv.coef\\[.*", bugvarnames)] * runif(3, min = 0.1, max = 0.2)
  artfit$mcmc[[2]][1, grepl("^LV\\[.*", bugvarnames)] <- artfit$mcmc[[1]][1, grepl("^LV\\[.*", bugvarnames)] * runif(4, min = 0.5, max = 1)
  
  # Predicted number of species detected and in occupation
  Enumspec <- predsumspecies(artfit, usefittedLV = TRUE)
  meanvar <- cumsum(Enumspec["Vsum_det", ])/((1:ncol(Enumspec))^2)
  sd_final <- sqrt(meanvar[ncol(Enumspec)])
  expect_equal(ncol(Enumspec), nsites)
  
  # Anticipate the Enumspec is wrong when using both draws, as simulated data in artfit is from the first draw
  NumSpecies_1st <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  meandiff_1st <- dplyr::cummean(NumSpecies_1st - Enumspec["Esum_det", ])
  plt <- cbind(diff = meandiff_1st, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  print(plt)
  expect_gt(abs(meandiff_1st[ncol(Enumspec)]), 3 * sd_final)
  
  # Anticipate the Enumspec is correct when using only first draw (chain), as simulated data in artfit is from the first draw
  Enumspec_1stonly <- predsumspecies(artfit, chain = 1, usefittedLV = TRUE)
  meanvar_1stonly <- cumsum(Enumspec_1stonly["Vsum_det", ])/((1:ncol(Enumspec_1stonly))^2)
  sd_final_1st <- sqrt(meanvar[ncol(Enumspec_1stonly)])
  
  meandiff_1st <- dplyr::cummean(NumSpecies_1st - Enumspec_1stonly["Esum_det", ])
  plt <- cbind(diff = meandiff_1st, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  print(plt)
  
  expect_lt(abs(meandiff_1st[ncol(Enumspec)]), 3 * sd_final_1st)
  
  # Anticipate that it is correct for 2nd draw separated from the 1st draw
  y_2nd <- simulate_fit(artfit, esttype = 2, UseFittedLV = TRUE)
  NumSpecies_2nd <- detectednumspec(y = y_2nd, ModelSite = artfit$data$ModelSite)
  Enumspec_2ndonly <- predsumspecies(artfit, chain = 2, usefittedLV = TRUE)
  meanvar_2ndonly <- cumsum(Enumspec_2ndonly["Vsum_det", ])/((1:ncol(Enumspec_2ndonly))^2)
  sd_final_2nd <- sqrt(meanvar[ncol(Enumspec_2ndonly)])
  meandiff_2nd <- dplyr::cummean(NumSpecies_2nd - Enumspec_2ndonly["Esum_det", ])
  plt <- cbind(diff = meandiff_2nd, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  print(plt)
  expect_lt(abs(meandiff_2nd[ncol(Enumspec)]), 3 * sd_final_2nd)
  
  
  # Anticipate prediction from combined draw is similar when occupancy + detection simulated with parameters chosen with equal chance from artfit$mcmc[[1]]
  ## combine observations for model sites to simulate the equal credence on each parameter set
  NumSpecies_interleaved <- NumSpecies_1st
  drawselect <- as.logical(rbinom(1000, size = 1, prob = 0.5))
  NumSpecies_interleaved[drawselect] <- NumSpecies_2nd[drawselect]
  
  meandiff <- dplyr::cummean(NumSpecies_interleaved - Enumspec["Esum_det", ])
  meanvar <- cumsum(Enumspec["Vsum_det", ])/((1:ncol(Enumspec))^2)
  plt <- cbind(diff = meandiff, var = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  print(plt)
  
  sd_final <- sqrt(meanvar[ncol(Enumspec)])
  expect_lt(abs(meandiff[ncol(Enumspec)]), 3 * sd_final)
})


test_that("In sample data; fitted LV values", {
  nsites <- 10000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, nlv = 4)
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  Enumspecdet <- predsumspecies(artfit, usefittedLV = TRUE)
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
  expect_lt(abs(meandiff[length(meandiff)]), abs(mean(meandiff[floor(length(meandiff) / 4) + 1:10 ])))
})

test_that("In sample data; marginal on LV values", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, nlv = 4,
                               u.b.min = -0.01,
                               u.b.max = 0.01,
                               v.b.min = -0.01,
                               v.b.max = 0.01,
                               lv.coef.min = 0.4,
                               lv.coef.max = 0.5 #hopefully LVs have a much bigger effect than occupancy etc
                               )
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  Enumspecdet <- predsumspecies(artfit, usefittedLV = FALSE, nLVsim = 1000)
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
  
  # check with predicted standard error: should be larger than bound due to misuse of LV
  sd_half <- sqrt(meanvar[floor(ncol(Enumspecdet)/2)])
  expect_gt(abs(meandiff[floor(ncol(Enumspecdet)/2)]), 3 * sd_half)
  
  ######################################### Now try using marginalised LV simulations #####################################
  NumSpecies <- detectednumspec(y = simulate_fit(artfit, esttype = 1, UseFittedLV = FALSE), ModelSite = artfit$data$ModelSite)
  
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
  expect_lt(abs(meandiff[length(meandiff)]), abs(mean(meandiff[floor(length(meandiff) / 4) + 1:10 ])))
})

test_that("In sample data; no LV", {
  nsites <- 10000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, nlv = 0)
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  Enumspecdet <- predsumspecies(artfit, usefittedLV = FALSE)
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
  expect_lt(abs(meandiff[length(meandiff)]), abs(mean(meandiff[floor(length(meandiff) / 4) + 1:10 ])))
})


test_that("Holdout data; has LVs", {
  nsites <- 10000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, nlv = 4) #means LV nearly don't matter
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  originalXocc <- Rfast::eachrow(Rfast::eachrow(artfit$data$Xocc, artfit$XoccProcess$scale, oper = "*"),
                                 artfit$XoccProcess$center, oper = "+")
  colnames(originalXocc) <- colnames(artfit$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- Rfast::eachrow(Rfast::eachrow(artfit$data$Xobs, artfit$XobsProcess$scale, oper = "*"),
                                 artfit$XobsProcess$center, oper = "+")
  colnames(originalXobs) <- colnames(artfit$data$Xobs)
  originalXobs <- cbind(ModelSite = artfit$data$ModelSite, originalXobs)
  outofsample_y <- simulate_fit(artfit, esttype = 1, UseFittedLV = FALSE)
  
  Enumspec <- predsumspecies_newdata(artfit, originalXocc, originalXobs, ModelSiteVars = "ModelSite", chains = NULL, nLVsim = 1000, cl = NULL)
  
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
  expect_lt(abs(meandiff[length(meandiff)]), abs(mean(meandiff[floor(length(meandiff) / 4) + 1:10 ])))
})



test_that("Holdout data; no LVs", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, nlv = 0) #means LV nearly don't matter
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  originalXocc <- Rfast::eachrow(Rfast::eachrow(artfit$data$Xocc, artfit$XoccProcess$scale, oper = "*"),
                                 artfit$XoccProcess$center, oper = "+")
  colnames(originalXocc) <- colnames(artfit$data$Xocc)
  originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
  originalXobs <- Rfast::eachrow(Rfast::eachrow(artfit$data$Xobs, artfit$XobsProcess$scale, oper = "*"),
                                 artfit$XobsProcess$center, oper = "+")
  colnames(originalXobs) <- colnames(artfit$data$Xobs)
  originalXobs <- cbind(ModelSite = artfit$data$ModelSite, originalXobs)
  outofsample_y <- simulate_fit(artfit, esttype = 1, UseFittedLV = FALSE)
  
  Enumspec <- predsumspecies_newdata(artfit, originalXocc, originalXobs, ModelSiteVars = "ModelSite", chains = NULL, nLVsim = 1000, cl = NULL)
  
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
  expect_lt(abs(meandiff[length(meandiff)]), abs(mean(meandiff[floor(length(meandiff) / 4) + 1:10 ])))
})

#########################################################################################

test_that("No LV and identical sites", {
  artfit <- artificial_runjags(nspecies = 4, nsites = 1000, nvisitspersite = 2, nlv = 0,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1")
  EVsum <- predsumspecies(artfit, usefittedLV = FALSE)
  
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
  expect_lt(abs(meandiff[length(meandiff)]), abs(mean(meandiff[floor(length(meandiff) / 4) + 1:10 ])))

  # Expect sd to be close to theoretical sd. Hopefully within 10%
  expect_equivalent(sd(NumSpecies[, "numspecies"]), sqrt(EVsum["Vsum_det", 1]), tol = 0.1 *sqrt(EVsum["Vsum_det", 1]))
})

test_that("Expected occupied number for in sample data; fitted LV values", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 3, nlv = 4,
                               v.b.min = 20, v.b.max = 20.1) #makes detection almost certain)
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  
  Enumspecdet <- predsumspecies(artfit, usefittedLV = TRUE)
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

