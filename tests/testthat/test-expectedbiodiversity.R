# test expected number of species

context("Number of Species Expected")

test_that("Number of species expected to be detected matches on artifical fitted model", {
  artfit <- artificial_runjags(nspecies = 2, nsites = 10000, nvisitspersite = 2, nlv = 0,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1")
  theta <- get_theta(artfit, type = 1)
  Enumspecdet <- expectedspeciesnum.ModelSite.theta(artfit$data$Xocc[1, , drop = FALSE],
                                                    artfit$data$Xobs[artfit$data$ModelSite == 1, , drop = FALSE],
                                                    numspecies = 2,
                                                    theta = get_theta(artfit, type = 1))
  
  # check that many other sites have the same expected number of species
  Enumspecdet_l <- lapply(1:100,
         function(idx){
           expectedspeciesnum.ModelSite.theta(artfit$data$Xocc[idx, , drop = FALSE],
                                              artfit$data$Xobs[artfit$data$ModelSite == idx, , drop = FALSE],
                                              numspecies = 2,
                                              theta = theta)
         })
  expect_equal(unlist(Enumspecdet_l) - Enumspecdet, rep(0, length(Enumspecdet_l)))
  
  # treat each model site as a repeat simulation of a ModelSite (cos all the parameters are nearly identical)
  my <- cbind(ModelSite = artfit$data$ModelSite, artfit$data$y)
  SpDetected <- my %>%
    as_tibble() %>%
    dplyr::group_by(ModelSite) %>%
    summarise_all(~sum(.) > 0)
  NumSpecies <- cbind(SpDetected[, 1] , numspecies = rowSums(SpDetected[, -1]))
  csum_numspecies <- cumsum(NumSpecies[, "numspecies"])
  
  # check that is getting closer with increasing data
  sim_Enumspecdet_small <- csum_numspecies[floor(nrow(artfit$data$Xocc)/4)] / floor(nrow(artfit$data$Xocc)/4)
  sim_Enumspecdet_big <- csum_numspecies[nrow(artfit$data$Xocc)] / nrow(artfit$data$Xocc)
  expect_gt(abs(sim_Enumspecdet_small - Enumspecdet), abs(sim_Enumspecdet_big - Enumspecdet))
  
  # check within standard error
  se <- sd(NumSpecies[, "numspecies"]) / sqrt(nrow(NumSpecies))
  expect_equal(sim_Enumspecdet_big, Enumspecdet, tol = 3*se)
})
