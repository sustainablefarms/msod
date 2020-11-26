
context("Runjags with few iterations")

test_that("Prep works with complicated formula", {
  artmodel <- artificial_runjags(nspecies = 60, nsites = 100, nvisitspersite = 2,
                                 OccFmla = "~ 1 + UpSite + Sine1 + Sine2 + UpSite*Sine2 + I(Sine1^2) + log(UpSite)",
                                 modeltype = "jsodm_lv", 
                                 nlv = 4)
  
  origXocc <- unstandardise.designmatprocess(artmodel$XoccProcess, artmodel$data$Xocc)
  origXocc <- cbind(ModelSite = 1:nrow(origXocc), origXocc)
  
  origXobs <- unstandardise.designmatprocess(artmodel$XobsProcess, artmodel$data$Xobs)
  origXobs <- cbind(ModelSite = artmodel$data$ModelSite, origXobs)
  
  rjo <- runjags::runjags.options("silent.jags" = TRUE)
  
  fit_runjags <- run.detectionoccupancy(origXocc, cbind(origXobs, artmodel$data$y), 
                                        species = colnames(artmodel$data$y),
                                        ModelSite = "ModelSite",
                                        OccFmla = artmodel$XoccProcess$fmla,
                                        ObsFmla = artmodel$XobsProcess$fmla,
                                        initsfunction = function(chain, indata){return(NULL)},
                                        MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 1, thin = 1),
                                        modeltype = "jsodm_lv",
                                        nlv = 4)
  
  runjags.options(rjo)
  
  expect_equal(fit_runjags$XoccProcess$version, 2)
  expect_equal(fit_runjags$XobsProcess$version, 2)
  expect_gt(sum(grepl("sigma", colnames(fit_runjags$mcmc[[1]]))), 0)
})

