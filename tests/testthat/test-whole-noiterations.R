
context("Runjags with few iterations")

test_that("Prep with complicated formula, save preparations", {
  indata <- artificial_covar_data(100, 2)
  # get y data
  artmodel <- artificial_runjags(nspecies = 3, nsites = 100, nvisitspersite = 2,
                                 OccFmla = "~ 1 + UpSite + Sine1 + Sine2 + UpSite*Sine2 + I(Sine1^2) + log(UpSite)",
                                 modeltype = "jsodm_lv", 
                                 nlv = 1)
  y <- artmodel$data$y
  rm(artmodel)
  
  toXoccParams <- sflddata::prep.designmatprocess(indata$Xocc, "~ 1 + UpSite + Sine1 + Sine2 + UpSite*Sine2 + I(Sine1^2) + log(UpSite)")
  toXoccFun <- function(indf, params = toXoccParams){
    Xocc <- apply.designmatprocess(params, indf)
    return(Xocc)
  }
  toXocc <- save_process(toXoccFun, checkwith = indata$Xocc, params = list(params = toXoccParams))
  Xocc <- apply_saved_process(toXocc, indata$Xocc)
  cmns <- colMeans(Xocc)
  sds <- apply(Xocc, MARGIN = 2, sd)
  expect_equivalent(cmns[c("UpSite", "Sine1", "Sine2", "log.UpSite.")], rep(0, 4))
  expect_equivalent(sds[c("UpSite", "Sine1", "Sine2", "log.UpSite.")], rep(1, 4), tolerance = 0.1)
  
  toXobsParams <- sflddata::prep.designmatprocess(indata$Xobs, "~ 1")
  toXobsFun <- function(indf, params = toXobsParams){
    X <- apply.designmatprocess(params, indf)
    return(X)
  }
  toXobs <- save_process(toXobsFun, checkwith = indata$Xobs, params = list(params = toXobsParams))
  Xobs <- apply_saved_process(toXobs, indata$Xobs)
  
  rjo <- runjags::runjags.options("silent.jags" = TRUE, silent.runjags = TRUE)
  
  suppressWarnings(fit <- fitjsodm2(Xocc = Xocc,
            Xobs = Xobs,
            y = y, 
            ModelSite = indata$Xobs$ModelSite, 
            modeltype = "jsodm_lv",
            nlv = 1,
            MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 1, thin = 1)
            ))
  
  runjags::runjags.options(rjo)
  
  fit$toXocc <- toXocc
  fit$toXobs <- toXobs
})

