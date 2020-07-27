library(discreteRV)

context("Full distributions of species numbers")

ci95generous <- function(rv){
  cs <- cumsum(probs(rv))
  valuesin <- outcomes(rv)[(0.025 < cs) & (0.975 > cs)]
  low <- min(valuesin)
  high <- max(valuesin)
  return(c(low = low, high = high))
}

test_that("Full artificial model same summaries as expected_biodiversity", {
  artfit <- artificial_runjags(nspecies = 60, nsites = 20, nvisitspersite = 3, nlv = 4)
  numspec_summ <- predsumspecies(artfit, UseFittedLV = TRUE, type = "marginal")
  numspec_distr <- predsumspeciesRV(artfit, UseFittedLV = TRUE, type = "marginal")
  
  expect_equivalent(numspec_summ["Esum_det", ], vapply(numspec_distr, E, FUN.VALUE = 0.23))
  expect_equivalent(numspec_summ["Vsum_det", ], vapply(numspec_distr, V, FUN.VALUE = 0.23), tol = 1E-5)
})

test_that("Full artificial model with excellent detection, same summaries as expected_biodiversity", {
  artfit <- artificial_runjags(nspecies = 60, nsites = 20, nvisitspersite = 3, nlv = 4,
                               ObsFmla = "~ 1",
                               v.b.min = 20, v.b.max = 20.1) #makes detection almost certain
  numspec_summ <- predsumspecies(artfit, UseFittedLV = TRUE, type = "marginal")
  numspec_distr <- predsumspeciesRV(artfit, UseFittedLV = TRUE, type = "marginal")
  
  expect_equivalent(numspec_summ["Esum_occ", ], numspec_summ["Esum_det", ])
  expect_equivalent(numspec_summ["Esum_occ", ], vapply(numspec_distr, E, FUN.VALUE = 0.23))
  expect_equivalent(numspec_summ["Vsum_occ", ], vapply(numspec_distr, V, FUN.VALUE = 0.23), tol = 1E-5)
})

test_that("Uncertainty dominated by latent variables.", {
  artfit <- artificial_runjags(nspecies = 60, nsites = 100, nvisitspersite = 3, nlv = 4,
                               u.b.min = -0.01,
                               u.b.max = 0.01,
                               v.b.min = -0.01,
                               v.b.max = 0.01,
                               lv.coef.min = 0.4,
                               lv.coef.max = 0.5 #hopefully LVs have a much bigger effect than occupancy etc
  )
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)

  sumRVs_margpost_fittedLV <- predsumspeciesRV(artfit, UseFittedLV = TRUE, type = "marginal")
  # par(mfrow = c(4, 5))
  # lapply(sumRVs_margpost_fittedLV[1:20], plot)
  ci_fittedLV <- lapply(sumRVs_margpost_fittedLV, ci95generous)
  ci_fittedLVsize <- vapply(ci_fittedLV, function(x) abs(x[[2]]- x[[1]]), FUN.VALUE = 0.0)
  
  inci_fittedLV <- mapply(function(ci, obsnum) {
                    return((ci[[1]] <= obsnum) & (ci[[2]] >= obsnum))
                  },
                 ci = ci_fittedLV,
                 obsnum = NumSpecies,
                 SIMPLIFY = TRUE)
  # as draws are true representation of distribution expect the 95% credible interval to cover observation 95% of the time
  expect_equal(mean(inci_fittedLV), 0.95, tol = 0.05)

  sumRVs_margpost_margLV <- predsumspeciesRV(artfit, UseFittedLV = FALSE, nLVsim = 100, type = "marginal")
  # par(mfrow = c(4, 5))
  # lapply(sumRVs_margpost_margLV, plot)
  ci_margLV <- lapply(sumRVs_margpost_margLV, ci95generous)
  ci_margLVsize <- vapply(ci_margLV, function(x) abs(x[[2]]- x[[1]]), FUN.VALUE = 0.0)
  
  inci_margLV <- mapply(function(ci, obsnum) {
                    return((ci[[1]] <= obsnum) & (ci[[2]] >= obsnum))
                  },
                  ci = ci_margLV,
                  obsnum = NumSpecies,
                  SIMPLIFY = TRUE)
  # although LV not supplied, the LV vals are not far from Gaussian so expect the predictions marginal across the LV distribution will work
  expect_gt(mean(inci_margLV), 0.9)
  
  # expect the confidence intervals to be larger when LV isn't known
  expect_true(all(ci_margLVsize - ci_fittedLVsize >= 0))
})

test_that("Credible intervals accurate for model without restrictions.", {
  artfit <- artificial_runjags(nspecies = 60, nsites = 200, nvisitspersite = 3, nlv = 4)
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  sumRVs_margpost_fittedLV <- predsumspeciesRV(artfit, UseFittedLV = TRUE, type = "marginal")
  # par(mfrow = c(4, 5))
  # lapply(sumRVs_margpost_fittedLV, plot)
  ci_fittedLV <- lapply(sumRVs_margpost_fittedLV, ci95generous)
  ci_fittedLVsize <- vapply(ci_fittedLV, function(x) abs(x[[2]]- x[[1]]), FUN.VALUE = 0.0)
  
  inci_fittedLV <- mapply(function(ci, obsnum) {
                            return((ci[[1]] <= obsnum) & (ci[[2]] >= obsnum))
                          },
                          ci = ci_fittedLV,
                          obsnum = NumSpecies,
                          SIMPLIFY = TRUE)
  # as draws are true representation of distribution expect the 95% credible interval to cover observation 95% of the time
  expect_equal(mean(inci_fittedLV), 0.95, tol = 0.05)
  
  
  sumRVs_margpost_margLV <- predsumspeciesRV(artfit, UseFittedLV = FALSE, nLVsim = 100, type = "marginal")
  # par(mfrow = c(4, 5))
  # lapply(sumRVs_margpost_margLV, plot)
  ci_margLV <- lapply(sumRVs_margpost_margLV, ci95generous)
  ci_margLVsize <- vapply(ci_margLV, function(x) abs(x[[2]]- x[[1]]), FUN.VALUE = 0.0)
  
  inci_margLV <- mapply(function(ci, obsnum) {
                          return((ci[[1]] <= obsnum) & (ci[[2]] >= obsnum))
                        },
                        ci = ci_margLV,
                        obsnum = NumSpecies,
                        SIMPLIFY = TRUE)
  # although LV not supplied, the LV vals are not far from Gaussian so expect the predictions marginal across the LV distribution will work
  expect_equal(mean(inci_margLV), 0.95, tol = 0.05)
  
  # expect the confidence intervals to be larger when LV isn't known
  expect_true(all(ci_margLVsize - ci_fittedLVsize >= 0))
})
