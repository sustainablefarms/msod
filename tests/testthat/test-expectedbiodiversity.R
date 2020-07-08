# test expected number of species

context("Number of Species Expected")

test_that("Number of detected species expected on artifical fitted model with no LV and no covariates", {
  artfit <- artificial_runjags(nspecies = 4, nsites = 10000, nvisitspersite = 2, nlv = 0,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1")
  theta <- get_theta(artfit, type = 1)
  lvsim <- matrix(rnorm(2 * 1), ncol = 2, nrow = 2) #dummy lvsim vars
  lv.coef.bugs <- matrix2bugsvar(matrix(0, nrow = artfit$data$n, ncol = 2), "lv.coef")
  theta <- c(theta, lv.coef.bugs)
  Enumspecdet <- expectedspeciesnum.ModelSite.theta(artfit$data$Xocc[1, , drop = FALSE],
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
  expect_equal(unlist(Enumspecdet_l) - Enumspecdet, rep(0, length(Enumspecdet_l)))
  
  # treat each model site as a repeat simulation of a ModelSite (cos all the parameters are nearly identical)
  my <- cbind(ModelSite = artfit$data$ModelSite, artfit$data$y)
  SpDetected <- my %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise_all(~sum(.) > 0)
  NumSpecies <- cbind(SpDetected[, 1] , numspecies = rowSums(SpDetected[, -1]))
  cmeannumspecies <- dplyr::cummean(NumSpecies[, "numspecies"])
  
  # check within standard error
  se <- sd(NumSpecies[, "numspecies"]) / sqrt(nrow(NumSpecies))
  expect_equal(cmeannumspecies[length(cmeannumspecies)], Enumspecdet, tol = 3*se)
  
  # check that is getting closer with increasing data
  expect_lt(abs(cmeannumspecies[length(cmeannumspecies)] - Enumspecdet),
            abs( mean(cmeannumspecies[floor(length(cmeannumspecies) / 4) + 1:10]) - Enumspecdet))
  
})


test_that("Number of detected species expected on artifical fitted model with covariates and no LV", {
  nsites <- 5000
  artfit <- artificial_runjags(nspecies = 5, nsites = nsites, nvisitspersite = 2, nlv = 0)
  theta <- get_theta(artfit, type = 1)
  # check that many other sites have the same expected number of species
  Enumspecdet_l <- lapply(1:nsites,
                          function(idx){
                            expectedspeciesnum.ModelSite(artfit,
                                                         artfit$data$Xocc[idx, , drop = FALSE],
                                                         artfit$data$Xobs[artfit$data$ModelSite == idx, , drop = FALSE],
                                                         LVvals = NULL)
                          })

  # treat each model site as a repeat simulation of a ModelSite (cos all the parameters are nearly identical)
  my <- cbind(ModelSite = artfit$data$ModelSite, artfit$data$y)
  SpDetected <- my %>%
    dplyr::as_tibble() %>%
    dplyr::group_by(ModelSite) %>%
    dplyr::summarise_all(~sum(.) > 0)
  NumSpecies <- cbind(SpDetected[, 1] , numspecies = rowSums(SpDetected[, -1]))
  
  # difference between expected and observed should be zero on average
  meandiff <- dplyr::cummean(NumSpecies[, "numspecies"] - unlist(Enumspecdet_l))
  plot(meandiff); abline(h = 0, col = "blue")
  
  # check that is getting closer with increasing data
  expect_lt(abs(meandiff[length(meandiff)]), abs(mean(meandiff[floor(length(meandiff) / 4) + 1:10 ])))
  
  # check with predicted standard error once the software is computed
})
