library(discreteRV)

context("Full distributions of species numbers")
pbopt <- pbapply::pboptions(type = "none")

ci95generous <- function(rv){
  cs <- cumsum(probs(rv))
  valuesin <- outcomes(rv)[(0.025 < cs) & (0.975 > cs)]
  low <- min(valuesin)
  high <- max(valuesin)
  return(c(low = low, high = high))
}

test_that("Full artificial model same summaries as expected_biodiversity", {
  artfit <- artificial_runjags(nspecies = 60, nsites = 20, nvisitspersite = 3, modeltype = "jsodm_lv", nlv = 4)
  numspec_summ <- speciesrichness(artfit, occORdetection = "detection", usefittedlvv = TRUE)
  numspec_distr <- predsumspeciesRV(artfit, UseFittedLV = TRUE, type = "marginal")
  
  expect_equal(numspec_summ["E", ], vapply(numspec_distr, E, FUN.VALUE = 0.23), ignore_attr = TRUE)
  expect_equal(numspec_summ["V", ], vapply(numspec_distr, V, FUN.VALUE = 0.23), tolerance = 1E-5, ignore_attr = TRUE)
})

test_that("Full artificial model with excellent detection, same summaries as expected_biodiversity", {
  artfit <- artificial_runjags(nspecies = 60, nsites = 20, nvisitspersite = 3,
                               ObsFmla = "~ 1",
                               det.b.min = 20, det.b.max = 20.1, #makes detection almost certain
			       modeltype = "jsodm_lv", nlv = 4) 
  numspec_summ <- speciesrichness(artfit, occORdetection = "occupancy", usefittedlvv = TRUE)
  numspec_distr <- predsumspeciesRV(artfit, UseFittedLV = TRUE, type = "marginal")
  
  expect_equal(numspec_summ["E", ], vapply(numspec_distr, E, FUN.VALUE = 0.23), ignore_attr = TRUE)
  expect_equal(numspec_summ["V", ], vapply(numspec_distr, V, FUN.VALUE = 0.23), tolerance = 1E-5, ignore_attr = TRUE)
})

test_that("Uncertainty dominated by latent variables.", {
  artfit <- artificial_runjags(nspecies = 60, nsites = 100, nvisitspersite = 3,
                               occ.b.min = -0.01,
                               occ.b.max = 0.01,
                               det.b.min = -0.01,
                               det.b.max = 0.01,
                               lv.b.min = 0.4,
                               lv.b.max = 0.5, #hopefully lv.vs have a much bigger effect than occupancy etc
                               modeltype = "jsodm_lv", nlv = 4
  )
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)

  sumRVs_margpost_fittedlv.v <- predsumspeciesRV(artfit, UseFittedLV = TRUE, type = "marginal")
  # par(mfrow = c(4, 5))
  # lapply(sumRVs_margpost_fittedlv.v[1:20], plot)
  ci_fittedlv.v <- lapply(sumRVs_margpost_fittedlv.v, ci95generous)
  ci_fittedlv.vsize <- vapply(ci_fittedlv.v, function(x) abs(x[[2]]- x[[1]]), FUN.VALUE = 0.0)
  
  inci_fittedlv.v <- mapply(function(ci, obsnum) {
                    return((ci[[1]] <= obsnum) & (ci[[2]] >= obsnum))
                  },
                 ci = ci_fittedlv.v,
                 obsnum = NumSpecies,
                 SIMPLIFY = TRUE)
  # as draws are true representation of distribution expect the 95% credible interval to cover observation 95% of the time
  expect_equal(mean(inci_fittedlv.v), 0.95, tolerance = 0.05)

  sumRVs_margpost_marglv.v <- predsumspeciesRV(artfit, UseFittedLV = FALSE, nLVsim = 100, type = "marginal")
  # par(mfrow = c(4, 5))
  # lapply(sumRVs_margpost_marglv.v, plot)
  ci_marglv.v <- lapply(sumRVs_margpost_marglv.v, ci95generous)
  ci_marglv.vsize <- vapply(ci_marglv.v, function(x) abs(x[[2]]- x[[1]]), FUN.VALUE = 0.0)
  
  inci_marglv.v <- mapply(function(ci, obsnum) {
                    return((ci[[1]] <= obsnum) & (ci[[2]] >= obsnum))
                  },
                  ci = ci_marglv.v,
                  obsnum = NumSpecies,
                  SIMPLIFY = TRUE)
  # although lv.v not supplied, the lv.vs are not far from Gaussian so expect the predictions marginal across the lv.v distribution will work
  expect_gt(mean(inci_marglv.v), 0.8)
  
  # expect the confidence intervals to be larger when lv.v isn't known
  expect_true(all(ci_marglv.vsize - ci_fittedlv.vsize >= 0))
})

test_that("Credible intervals accurate for model without restrictions.", {
  artfit <- artificial_runjags(nspecies = 60, nsites = 200, nvisitspersite = 3, modeltype = "jsodm_lv", nlv = 4)
  artfit$mcmc[[1]] <- rbind(artfit$mcmc[[1]][1, ], artfit$mcmc[[1]][1, ])
  NumSpecies <- detectednumspec(y = artfit$data$y, ModelSite = artfit$data$ModelSite)
  
  sumRVs_margpost_fittedlv.v <- predsumspeciesRV(artfit, UseFittedLV = TRUE, type = "marginal")
  # par(mfrow = c(4, 5))
  # lapply(sumRVs_margpost_fittedlv.v, plot)
  ci_fittedlv.v <- lapply(sumRVs_margpost_fittedlv.v, ci95generous)
  ci_fittedlv.vsize <- vapply(ci_fittedlv.v, function(x) abs(x[[2]]- x[[1]]), FUN.VALUE = 0.0)
  
  inci_fittedlv.v <- mapply(function(ci, obsnum) {
                            return((ci[[1]] <= obsnum) & (ci[[2]] >= obsnum))
                          },
                          ci = ci_fittedlv.v,
                          obsnum = NumSpecies,
                          SIMPLIFY = TRUE)
  # as draws are true representation of distribution expect the 95% credible interval to cover observation 95% of the time
  expect_equal(mean(inci_fittedlv.v), 0.95, tolerance = 0.05)
  
  
  sumRVs_margpost_marglv.v <- predsumspeciesRV(artfit, UseFittedLV = FALSE, nLVsim = 100, type = "marginal")
  # par(mfrow = c(4, 5))
  # lapply(sumRVs_margpost_marglv.v, plot)
  ci_marglv.v <- lapply(sumRVs_margpost_marglv.v, ci95generous)
  ci_marglv.vsize <- vapply(ci_marglv.v, function(x) abs(x[[2]]- x[[1]]), FUN.VALUE = 0.0)
  
  inci_marglv.v <- mapply(function(ci, obsnum) {
                          return((ci[[1]] <= obsnum) & (ci[[2]] >= obsnum))
                        },
                        ci = ci_marglv.v,
                        obsnum = NumSpecies,
                        SIMPLIFY = TRUE)
  # although lv.v not supplied, the lv.vs are not far from Gaussian so expect the predictions marginal across the lv.v distribution will work
  expect_equal(mean(inci_marglv.v), 0.95, tolerance = 0.05)
  
  # expect the confidence intervals to be larger when lv.v isn't known
  expect_gt(mean(ci_marglv.vsize - ci_fittedlv.vsize >= 0), 0.9) #nearly all of the time exceptions could be small rounding differences
  expect_true(all(ci_marglv.vsize - ci_fittedlv.vsize >= -1)) #differences could be small rounding differences
  
  smallerci <- (ci_marglv.vsize - ci_fittedlv.vsize < 0)
  comparervplot <- function(rv1, rv2){
    plot(outcomes(rv1), probs(rv1))
    points(outcomes(rv2), probs(rv2), col = "red", pch = "+")
  }
  par(mfrow = c(4, 5))
  # mapply(comparervplot, rv1 = sumRVs_margpost_marglv.v[smallerci], rv2 = sumRVs_margpost_fittedlv.v[smallerci])
})

pbapply::pboptions(pbopt)
