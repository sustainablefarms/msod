# test of simulation methods.

context("Simulation Methods")

test_that("Opposite loadings of lv.v gives anticorrelated simulated detections for species", {
  # simulate for species with high Detection, opposite impact by the lv.v
  artfit1 <- artificial_runjags(nspecies = 2, nsites = 1000, nvisitspersite = 1,
                                ObsFmla = "~ 1",
                                OccFmla = "~ 1",
                                occ.b.min = -1,  occ.b.max = -0.9, 
                                det.b.min = 20, det.b.max = 20.1, #makes detection almost certain
                                lv.b.min = matrix(c(-0.6, 0.6), ncol = 1), lv.b.max = matrix(c(-0.55, 0.65), ncol = 1),
                                modeltype = "jsodm_lv",
                                nlv = 1)
  pocc_corr <- cor(poccupy_species(artfit1, type = 1, conditionalLV = TRUE))
  expect_equivalent(pocc_corr, matrix(c(1, -1, -1, 1), ncol = 2, nrow = 2, byrow = FALSE))
  
  expect_equivalent(colMeans(artfit1$data$y), colMeans(poccupy_species(artfit1, type = 1, conditionalLV = TRUE)), tolerance = 0.1)
  anticor_occ <- cor(artfit1$data$y)
  expect_lt(anticor_occ[1, 2], 0)
})
  
test_that("High, equal loadings of lv.v gives correlated simulated detections for species", {
  artfit2 <- artificial_runjags(nspecies = 2, nsites = 1000, nvisitspersite = 1,
                                ObsFmla = "~ 1",
                                OccFmla = "~ 1",
                                occ.b.min = -1,  occ.b.max = -0.9,
                                det.b.min = 20, det.b.max = 20.1,
                                lv.b.min = matrix(c(0.6, 0.6), ncol = 1), lv.b.max = matrix(c(0.65, 0.65), ncol = 1),
                                modeltype = "jsodm_lv",
                                nlv = 1)

  expect_equivalent(cor(poccupy_species(artfit2, type = 1, conditionalLV = TRUE)), matrix(1, ncol = 2, nrow = 2, byrow = FALSE))
  expect_equivalent(colMeans(artfit2$data$y), colMeans(poccupy_species(artfit2, type = 1, conditionalLV = TRUE)), tolerance = 0.1)
  cor_occ <- cor(artfit2$data$y)
  expect_gt(cor_occ[1, 2], 0)
})


test_that("Distributions of Model Site Values differ when using lv.vs that aren't Gaussian", {
  # the third lv.v is not evenly distributed across the sites
  artfit <- artificial_runjags(nspecies = 2, nsites = 10000, nvisitspersite = 1,
                               ObsFmla = "~ 1",
                               OccFmla = "~ 1",
                               occ.b.min = 0,  occ.b.max = 0.001, 
                               det.b.min = 1, det.b.max = 1.001,
                               lv.b.min = matrix(c(0, 0, 0.45), nrow = 2, ncol = 3, byrow = TRUE),
                               lv.b.max = matrix(c(0, 0, 0.65), nrow = 2, ncol = 3, byrow = TRUE),
                               modeltype = "jsodm_lv",
                               nlv = 3)
  
  # resimulate y as if lv.v are not known (which is the situation for likelihood compuations)
  even_y <- simulate_detections_lv.v(artfit, esttype = 1)
  
  # get joint outcomes
  my <- cbind(ModelSite = artfit$data$ModelSite, even_y)
  obs_per_site <- lapply(1:nrow(artfit$data$Xocc), function(x) my[my[, "ModelSite"] == x, -1])
  jointoutcomes <- vapply(obs_per_site, paste0, collapse = ",", FUN.VALUE = "achar")
  
  # likelihood by simulation
  even_sim_distr <- summary(factor(jointoutcomes))/length(jointoutcomes)
  
  # resimulate y as if lv.v are not known (which is the situation for likelihood compuations)
  uneven_y <- simulate_detections(artfit, esttype = 1)
  
  # get joint outcomes
  my <- cbind(ModelSite = artfit$data$ModelSite, uneven_y)
  obs_per_site <- lapply(1:nrow(artfit$data$Xocc), function(x) my[my[, "ModelSite"] == x, -1])
  jointoutcomes <- vapply(obs_per_site, paste0, collapse = ",", FUN.VALUE = "achar")
  
  # likelihood by simulation
  uneven_sim_distr <- summary(factor(jointoutcomes))/length(jointoutcomes)
  
  expect_gt(uneven_sim_distr[1] - 0.01, even_sim_distr[1]) 
  expect_lt(uneven_sim_distr[4] - 0.01, even_sim_distr[4]) 
})
 
