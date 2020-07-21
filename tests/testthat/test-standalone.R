context("Compare standalone occupancy probability to other occupancy functions")

test_that("standalone matches poccupy_species", {
  fit <- artificial_runjags(nspecies = 5, nsites = 100, nvisitspersite = 1, nlv = 0)
  a <- poccupy_species(fit, type = 1, conditionalLV = FALSE)
  XoccOrig <- uncentre.designmatprocess(fit$XoccProcess, fit$data$Xocc)
  theta <- get_theta(fit, type = 1)
  u.b <- bugsvar2matrix(theta, "u.b", 1:fit$data$n, 1:ncol(fit$data$Xocc))
  b <- poccupancy_standalone_nolv(XoccOrig, fit$XoccProcess, u.b)
  expect_equivalent(a, b)
})


test_that("Multisite richness function matches others for single sites", {
  fit <- artificial_runjags(nspecies = 60, nsites = 100, nvisitspersite = 1, nlv = 0)
  Erichnesspersite <- predsumspecies(fit, UseFittedLV = FALSE)
  
  XoccOrig <- uncentre.designmatprocess(fit$XoccProcess, fit$data$Xocc)
  theta <- get_theta(fit, type = 1)
  u.b <- bugsvar2matrix(theta, "u.b", 1:fit$data$n, 1:ncol(fit$data$Xocc))
  XoccOrig[1, , drop = FALSE]
  b <- lapply(1:nrow(XoccOrig), function(i) {
    multisiterichness_nolv(XoccOrig[i, , drop = FALSE], fit$XoccProcess, u.b)})
  b <- simplify2array(b)
  expect_equivalent( b["Erichness", ], Erichnesspersite["Esum_occ", ])
  expect_equivalent( b["Vrichness", ], Erichnesspersite["Vsum_occ", ])
})