context("Test strange prep design matrix twice error")

# these tests passed on:
# e1c15fe132493a64a3075533e449ac47f8f2a66e, 
# b4eace67748910756d05733322e0b82d15c867a5 on Sep 9
# October 9 commit 24cc8c17c82ce582ce45f335fc1617362f60268f
# Oct 12 commit deed04c0371d6111d48e46853c51abc790c83ab9
# Final Oct 15 commit b70b9d646cd23b46ba0b15046111a6f6fc1faf2a
# First errors on Oct 20 commit b81b7c02ca1affcb989debd0d8badf7adc08b335
# Error was due to stuffing up resetting global option for na.action
library(sustfarmld); library(testthat)

test_that("Error does not occur using prep design matrix process twice", {
  covardf <- simulate_covar_data(20, 2)
  OccFmla = "~ 1"
  XoccProcess <- prep.designmatprocess(covardf$Xocc, OccFmla)
  Xocc <- apply.designmatprocess(XoccProcess, covardf$Xocc)
  XoccProcess <- prep.designmatprocess(covardf$Xocc, OccFmla)
})

test_that("Error does not occur via artificial_runjags", {
  artmodel <- artificial_runjags(nspecies = 5, nsites = 2000, nvisitspersite = 2, nlv = 0,
                                 ObsFmla = "~ 1",
                                 OccFmla = "~ 1")
  artificial_runjags(nspecies = 5, nsites = 2000, nvisitspersite = 2, nlv = 0,
                                  ObsFmla = "~ 1",
                                  OccFmla = "~ 1")
})

test_that("Error does not occur using prep design matrix process twice", {
  covardf <- simulate_covar_data(20, 2)
  OccFmla = "~ 1"
  XoccProcess <- prep.designmatprocess(covardf$Xocc, OccFmla)
  XoccProcess <- prep.designmatprocess(covardf$Xocc, OccFmla)
})
