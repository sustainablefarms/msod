# test expected number of species

context("Number of Detected Species Expected")

test_that("Correct for artifical fitted model with no LV and identical sites", {
  artfit <- artificial_runjags(nspecies = 4, nsites = 10000, nvisitspersite = 2, nlv = 0,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1")
  theta <- get_theta(artfit, type = 1)
  lvsim <- matrix(rnorm(2 * 1), ncol = 2, nrow = 2) #dummy lvsim vars
  lv.coef.bugs <- matrix2bugsvar(matrix(0, nrow = artfit$data$n, ncol = 2), "lv.coef")
  theta <- c(theta, lv.coef.bugs)
  EVsum_det_single <- expectedspeciesnum.ModelSite.theta(artfit$data$Xocc[1, , drop = FALSE],
                                                    artfit$data$Xobs[artfit$data$ModelSite == 1, , drop = FALSE],
                                                    numspecies = 4,
                                                    theta = theta,
                                                    LVvals = lvsim)
  
  # check that many other sites have the same expected number of species
  Enumspecdet_l <- lapply(1:100,
         function(idx){
           expectedspeciesnum.ModelSite.theta(artfit$data$Xocc[idx, , drop = FALSE],
                                              artfit$data$Xobs[artfit$data$ModelSite == idx, , drop = FALSE],
                                              numspecies = 4,
                                              theta = theta,
                                              LVvals = lvsim)
         })
  Enumspecdet <- simplify2array(Enumspecdet_l)
  expect_equivalent(Enumspecdet["Esum_det", ], rep(EVsum_det_single["Esum_det"], length(Enumspecdet_l)))
  
  # treat each model site as a repeat simulation of a ModelSite (cos all the parameters are nearly identical)
  my <- cbind(ModelSite = artfit$data$ModelSite, artfit$data$y)
  SpDetected <- my %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise_all(~sum(.) > 0)
  NumSpecies <- cbind(SpDetected[, 1] , numspecies = rowSums(SpDetected[, -1]))
  cmeannumspecies <- dplyr::cummean(NumSpecies[, "numspecies"])
  
  # check within standard error
  sd <- sd(NumSpecies[, "numspecies"]) / sqrt(nrow(NumSpecies))
  expect_equivalent(cmeannumspecies[length(cmeannumspecies)], EVsum_det_single["Esum_det"], tol = 3*sd)
  
  # check that is getting closer with increasing data
  expect_lt(abs(cmeannumspecies[length(cmeannumspecies)] -  EVsum_det_single["Esum_det"]),
            abs( mean(cmeannumspecies[floor(length(cmeannumspecies) / 4) + 1:10]) - EVsum_det_single["Esum_det"]))
  
  # Expect sd to be close to theoretical sd. Hopefully within 10%
  expect_equivalent(sd(NumSpecies[, "numspecies"]), sqrt(EVsum_det_single["Vsum_det"]), tol = 0.1 *sqrt(EVsum_det_single["Vsum_det"]))
})

test_that("Correct for artifical fitted model with covariates but no LV", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 5, nsites = nsites, nvisitspersite = 2, nlv = 0)
  # check that many other sites have the same expected number of species
  Enumspecdet_l <- lapply(1:nsites,
                          function(idx){
                            expectedspeciesnum.ModelSite(artfit,
                                                         artfit$data$Xocc[idx, , drop = FALSE],
                                                         artfit$data$Xobs[artfit$data$ModelSite == idx, , drop = FALSE],
                                                         LVvals = NULL)
                          })
  Enumspecdet <- simplify2array(Enumspecdet_l)

  # treat each model site as a repeat simulation of a ModelSite (cos all the parameters are nearly identical)
  my <- cbind(ModelSite = artfit$data$ModelSite, artfit$data$y)
  SpDetected <- my %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise_all(~sum(.) > 0)
  NumSpecies <- cbind(SpDetected[, 1] , numspecies = rowSums(SpDetected[, -1]))
  
  meandiff <- dplyr::cummean(NumSpecies[, "numspecies"] - Enumspecdet["Esum_det", ])
  # difference between expected and observed should be zero on average; check that is getting closer with increasing data
  expect_lt(abs(meandiff[length(meandiff)]), abs(mean(meandiff[floor(length(meandiff) / 4) + 1:10 ])))
  # plot(meandiff); abline(h = 0, col = "blue")
  
  # check with predicted standard error once the software is computed
  meanvar <- cumsum(Enumspecdet["Vsum_det", ])/((1:ncol(Enumspecdet))^2)
  sd_final <- sqrt(meanvar[ncol(Enumspecdet)])
  expect_equal(meandiff[ncol(Enumspecdet)], 0, tol = 3 * sd_final)
})

test_that("Correct for artifical fitted model with covariates and LVs", {
  nsites <- 10000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 2, nlv = 2)
  
  Enumspecdet <- predsumspecies(artfit, usefittedLV = TRUE)
  
  my <- cbind(ModelSite = artfit$data$ModelSite, artfit$data$y)
  SpDetected <- my %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise_all(~sum(.) > 0)
  NumSpecies <- cbind(SpDetected[, 1] , numspecies = rowSums(SpDetected[, -1]))
  
  meandiff <- dplyr::cummean(NumSpecies[, "numspecies"] - Enumspecdet["Esum_det", ])
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

test_that("Correct for artifical fitted model marginal on LV", {
  artfit <- artificial_runjags(nspecies = 4, nsites = 10000, nvisitspersite = 2, nlv = 4,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1")
  theta <- get_theta(artfit, type = 1)
  lvsim <- matrix(rnorm(artfit$data$nlv * 1000), ncol = artfit$data$nlv, nrow = 1000) #dummy lvsim vars
  LVvals_fixed <- bugsvar2array(get_theta(artfit, type = 1), "LV", 1:nrow(artfit$data$Xocc), 1:artfit$data$nlv)[ , , 1]
  EVsum_single <- expectedspeciesnum.ModelSite.theta(artfit$data$Xocc[1, , drop = FALSE],
                                                         artfit$data$Xobs[artfit$data$ModelSite == 1, , drop = FALSE],
                                                         numspecies = 4,
                                                         theta = theta,
                                                         LVvals = lvsim)
  
  # check that many other sites have the same expected number of species
  Enumspecdet_l <- lapply(1:10,
                          function(idx){
                            expectedspeciesnum.ModelSite.theta(artfit$data$Xocc[idx, , drop = FALSE],
                                                               artfit$data$Xobs[artfit$data$ModelSite == idx, , drop = FALSE],
                                                               numspecies = 4,
                                                               theta = theta,
                                                               LVvals = lvsim)
                          })
  Enumspecdet <- simplify2array(Enumspecdet_l)
  expect_equivalent(Enumspecdet["Esum_det", ], rep(EVsum_single["Esum_det"], length(Enumspecdet_l)))
  
  # simulate observations for each model site from iid LV values; treat each model site as a repeat simulation of a ModelSite (cos all the parameters are nearly identical)
  my <- cbind(ModelSite = artfit$data$ModelSite, simulate_fit(artfit, esttype = 1, UseFittedLV = FALSE))
  SpDetected <- my %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise_all(~sum(.) > 0)
  NumSpecies <- cbind(SpDetected[, 1] , numspecies = rowSums(SpDetected[, -1]))
  cmeannumspecies <- dplyr::cummean(NumSpecies[, "numspecies"])
  
  # check within standard error
  sd <- sqrt(EVsum_single["Vsum_det"]) / sqrt(nrow(NumSpecies))
  expect_equivalent(cmeannumspecies[length(cmeannumspecies)], EVsum_single["Esum_det"], tol = 3*sd)
  
  # Expect sd to be close to theoretical sd. Hopefully within 10%
  expect_equivalent(sd(NumSpecies[, "numspecies"]), sqrt(EVsum_single["Vsum_det"]), tol = 0.1 *sqrt(EVsum_single["Vsum_det"]))
  
  # simulate observations for each model site from prefixed non-symmatric LV values - should make the E and SD incorrect.
  my <- cbind(ModelSite = artfit$data$ModelSite, simulate_fit(artfit, esttype = 1, UseFittedLV = TRUE))
  SpDetected <- my %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise_all(~sum(.) > 0)
  NumSpecies <- cbind(SpDetected[, 1] , numspecies = rowSums(SpDetected[, -1]))
  cmeannumspecies <- dplyr::cummean(NumSpecies[, "numspecies"])
  
  expect_gt(abs(cmeannumspecies[length(cmeannumspecies)] - EVsum_single["Esum_det"]), 3 * sd)
})

context("Number of Occupied Species Expected")

test_that("Correct for artifical fitted model with covariates and LVs, certain detection.", {
  nsites <- 1000
  artfit <- artificial_runjags(nspecies = 60, nsites = nsites, nvisitspersite = 2, nlv = 2,
                                ObsFmla = "~ 1",
                                v.b.min = 20, v.b.max = 20.1) #makes detection almost certain
  LVvals <- bugsvar2array(get_theta(artfit, type = 1), "LV", 1:nrow(artfit$data$Xocc), 1:artfit$data$nlv)[ , , 1]
  # check that many other sites have the same expected number of species
  Enumspec_l <- lapply(1:nsites,
                          function(idx){
                            expectedspeciesnum.ModelSite(artfit,
                                                         artfit$data$Xocc[idx, , drop = FALSE],
                                                         artfit$data$Xobs[artfit$data$ModelSite == idx, , drop = FALSE],
                                                         LVvals = LVvals[idx, , drop = FALSE])
                          })
  Enumspec <- simplify2array(Enumspec_l)
  
  # treat each model site as a repeat simulation of a ModelSite (cos all the parameters are nearly identical)
  my <- cbind(ModelSite = artfit$data$ModelSite, artfit$data$y)
  SpDetected <- my %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise_all(~sum(.) > 0)
  NumSpecies <- cbind(SpDetected[, 1] , numspecies = rowSums(SpDetected[, -1]))
  
  meandiff <- dplyr::cummean(NumSpecies[, "numspecies"] - Enumspec["Esum_occ", ])
  meanvar <- cumsum(Enumspec["Vsum_occ", ])/((1:ncol(Enumspec))^2)
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
})

