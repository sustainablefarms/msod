# test saving processing function

test_that("Saving and running a data processing function", {
  set.seed(3546)
  B <- matrix(rnorm(4 * 3), nrow = 3, ncol = 4)
  center <- colMeans(B)
  scale <- apply(B, MARGIN = 2, sd)
  
  myfun <- function(indf, center, scale){
    out <- (indf - center) * 1/scale
    return(out)
  }
  Bstd <- myfun(B, center, scale)
  
  saved <- save_process(myfun, checkwith = B, params = list(center = center, scale = scale))
  rm(list = setdiff(ls(), c("saved", "Bstd")))
  
  set.seed(3546)
  B <- matrix(rnorm(4 * 3), nrow = 3, ncol = 4)
  Bstd2 <- apply_process(saved$fun, indf = B, params = saved$params)
  expect_identical(Bstd, Bstd2)
})

test_that("Saving a processing function fails when global names used inside function body", {
  B <- matrix(rnorm(4 * 3), nrow = 3, ncol = 4)
  cmns <- colMeans(B)
  sd <- apply(B, MARGIN = 2, sd)
  myfun <- function(indf, sd){
    out <- (indf - cmns) * sd
    return(out)
  }
  expect_error(save_process(myfun, checkwith = B, params = list(sd = sd)))
})

test_that("Saving process with a column that has standardised values for some rows and zeros for the other rows", {
  Xocc <- data.frame(ModelSite = 1:10, IsPlanting = 0, PlantingAge = NA)
  Xocc$IsPlanting[6:10] <- 1 #these sites are plantings
  Xocc$PlantingAge[6:10] <- 0.1 + (1:5) #these sites are plantings
  centers <- c(IsPlanting = 0, PlantingAge = mean(Xocc$PlantingAge[6:10]))
  scales <- c(IsPlanting = 1, PlantingAge = sd(Xocc$PlantingAge[6:10]))

  BldModMat <- function(indf = Xocc, scale, center){
    indfstd <- msod:::apply_center_n_scale(indf, center, scale)
    ModMat <- model.matrix(as.formula("~ IsPlanting"), as.data.frame(indfstd))
    ModMat <- cbind(ModMat, PlantingAge = indfstd[,"PlantingAge"])
    ModMat[indf$IsPlanting == 0, "PlantingAge"] <- 0
    return(ModMat)
  }
  
  saved <- save_process(BldModMat, checkwith = Xocc, params = list(center = centers, scale = scales))
  rm(list = setdiff(ls(), c("saved", "Xocc")))
  
  ModMat <- apply_process(saved$fun, Xocc, params = saved$params)
  expect_equal(ModMat, apply_saved_process(saved, Xocc))
  
  expect_equivalent(ModMat[1:5, "PlantingAge"], Xocc$IsPlanting[1:5])
  expect_equivalent(ModMat[1:5, "IsPlanting"], Xocc$IsPlanting[1:5])
})
  