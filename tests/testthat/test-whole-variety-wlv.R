############### Wholistic Test on different ModelSites and lv.vs #########
# compare artmodel generated likelihood and Enum species to fit_runjags versions, and Enum species to observations

context("Wholistic tests on model with different ModelSites and lv.vs")
skip_if(parallel::detectCores() < 10)

rjo <- runjags::runjags.options("silent.jags" = TRUE)

# Create a process with known parameters
artmodel <- artificial_runjags(nspecies = 10, nsites = 2000, nvisitspersite = 2,modeltype = "jsodm_lv",  nlv = 4)

# fit to data and simulations using runjags
originalXocc <- unstandardise.designmatprocess(artmodel$XoccProcess, artmodel$data$Xocc)
originalXocc <- cbind(ModelSite = 1:nrow(originalXocc), originalXocc)
originalXobs <- unstandardise.designmatprocess(artmodel$XobsProcess, artmodel$data$Xobs)
originalXobs <- cbind(ModelSite = artmodel$data$ModelSite, originalXobs)

fit_runjags <- run.detectionoccupancy(originalXocc, cbind(originalXobs, artmodel$data$y), 
                                      species = colnames(artmodel$data$y),
                                      ModelSite = "ModelSite",
                                      OccFmla = artmodel$XoccProcess$fmla,
                                      ObsFmla = artmodel$XobsProcess$fmla,
                                      # initsfunction = function(chain, indata){return(NULL)},
                                      # MCMCparams = list(n.chains = 2, adapt = 1000, burnin = 10000, sample = 500, thin = 40),
                                      nlv = 4)
save(fit_runjags, artmodel, originalXocc, originalXobs, file = "benchmark_varietysitesmodel_initsgiven.Rdata")

runjags.options(rjo)

gwk <- tibble::enframe(coda::geweke.diag(fit_runjags, frac1=0.1, frac2=0.5)$z, name = "varname")
qqnorm(gwk$value)
qqline(gwk$value)
varname2type <- function(varnames){
  types <- dplyr::case_when(
    grepl("lv.b", varnames) ~ "lv.v Load",
    grepl("lv.v.*,1\\]", varnames) ~ "lv.v1",
    grepl("lv.v.*,2\\]", varnames) ~ "lv.v2",
    grepl("lv.v.*,3\\]", varnames) ~ "lv.v3",
    grepl("lv.v.*,4\\]", varnames) ~ "lv.v4",
    grepl("^(mu|tau)", varnames) ~ "Comm Param", #parameters of community distributions
    grepl("^occ.b", varnames) ~ "Occu Coef",
    grepl("^det.b", varnames) ~ "Detn Coef",
    TRUE ~ "other"
  )
  return(types)
}
gwk %>%
  dplyr::mutate(type = varname2type(varname)) %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(swp = shapiro.test(value)$p.value) %>%
  ggplot2::ggplot() +
  ggplot2::facet_wrap(~type) +
  ggplot2::geom_vline(ggplot2::aes(xintercept = 0.01)) +
  ggplot2::geom_point(ggplot2::aes(y = 0, x = swp)) + 
  ggplot2::scale_x_continuous(name = "Shapiro-Wilk p-value",
                              trans = "identity") +
  ggplot2::ggtitle("Geweke Convergence Statistics Normality Tests",
                   subtitle = "0.01 threshold shown")
# lv.v did not converge! --> more burnin and higher filter?
hist(fit_runjags$summaries[,"AC.300"]) #Autocorellation looks ok.


test_that("Posterior credible distribution overlaps true parameters", {
  var2compare <- colnames(artmodel$mcmc[[1]])
  inCI <- (fit_runjags$summaries[var2compare, "Lower95"] <= artmodel$mcmc[[1]][1, var2compare]) &
    (fit_runjags$summaries[var2compare, "Upper95"] >= artmodel$mcmc[[1]][1, var2compare])
  inCI[!grepl("^lv.v", names(inCI))]
  expect_equal(mean(inCI), 0.95, tol = 0.05)
})

test_that("Median of Posterior is close to true", {
  var2compare <- colnames(artmodel$mcmc[[1]])
  res <- artmodel$mcmc[[1]][1, var2compare] - fit_runjags$summaries[var2compare, "Median"] 
  relres <- res / artmodel$mcmc[[1]][1, var2compare]
  plot(res[!grepl("^lv.v", names(res))])
  plot(res[grepl("^lv.v", names(res))])
  plot(relres[!grepl("^lv.v", names(relres))])
  plot(relres[grepl("^lv.v", names(relres))])
  expect_equivalent(relres[!grepl("^lv.v", names(res))],
                    rep(0, sum(!grepl("^lv.v", names(res)))),
                    tol = 0.1)
  expect_equivalent(relres[!grepl("lv.v", names(res))],
                    rep(0, sum(!grepl("lv.v", names(res)))),
                    tol = 0.2)
})

test_that("Fitted likelihood matches true likelihood", {
  cl <- parallel::makeCluster(15)
  lkl_runjags <- likelihoods.fit(fit_runjags, cl = cl)
  lkl_artmodel <- likelihoods.fit(artmodel, cl = cl)
  parallel::stopCluster(cl)
  expect_equivalent(Rfast::colmeans(lkl_runjags), Rfast::colmeans(lkl_artmodel), tol = 0.01)
  expect_equivalent(Rfast::colmeans(lkl_runjags) - Rfast::colmeans(lkl_artmodel),
                    rep(0, ncol(lkl_artmodel)),
                    tol = 0.1 * Rfast::colmeans(lkl_artmodel))
})

test_that("Expected Number of Detected Species", {
  cl <- parallel::makeCluster(10)
  pbopt <- pbapply::pboptions(type = "none")
  Enumspec <- predsumspecies(fit_runjags, UseFittedLV = TRUE, type = "marginal", cl = cl)
  Enumspec_art <- predsumspecies(artmodel, UseFittedLV = TRUE, type = "marginal", cl = cl)
  pbapply::pboptions(pbopt)
  parallel::stopCluster(cl)
  cbind(rj = t(Enumspec), art = t(Enumspec_art)) %>%
    as_tibble() %>%
    tibble::rowid_to_column() %>%
    ggplot2::ggplot()
  expect_equivalent(Enumspec, Enumspec_art)
  
  NumSpecies <- detectednumspec(y = fit_runjags$data$y, ModelSite = fit_runjags$data$ModelSite)
  
  meandiff <- dplyr::cummean(NumSpecies - Enumspec["Esum_det", ])
  meanvar <- cumsum(Enumspec["Vsum_det", ])/((1:ncol(Enumspec))^2)
  plt <- cbind(diff = meandiff, var  = meanvar) %>% 
    dplyr::as_tibble(rownames = "CumSites") %>% 
    dplyr::mutate(CumSites = as.double(CumSites)) %>%
    ggplot2::ggplot() +
    ggplot2::geom_ribbon(ggplot2::aes(x= CumSites, ymin = -2 * sqrt(var), ymax = 2 * sqrt(var)), fill = "grey") +
    ggplot2::geom_line(ggplot2::aes(x = CumSites, y = diff), col = "blue", lwd = 2)
  #print(plt)
  
  # check 'average' model site numbers are correct within standard errors, and that standard error is close.
  Enum_compare_sum <- Enum_compare(NumSpecies,
                                   as.matrix(Enumspec["Esum_det", ], ncol = 1),
                                   as.matrix(Enumspec["Vsum_det", ], ncol = 1)
  )
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_model"]])
  expect_equivalent(Enum_compare_sum[["E[D]_obs"]], 0, tol = 3 * Enum_compare_sum[["SE(E[D]_obs)_obs"]])
  expect_equivalent(Enum_compare_sum[["V[D]_model"]], Enum_compare_sum[["V[D]_obs"]], tol = 0.05 * Enum_compare_sum[["V[D]_obs"]])
})

