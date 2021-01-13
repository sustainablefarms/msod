#' @title Functional richness, eveness and diversity from the FD package
#' @param pocc An array of occupancy probabiliies. Each row is a site (can be used as a draw), each column a species.
#' @param traits A data frame with each species a row, each column a trait. Row names must match column names of `pocc`
#' @param nsim The number of simulation of occupancy to perform for each site.
#' @return 
#' @details
#' Occupancy of each species between each site is considered independent.
#' The occupancy simulation of each species converge relatively slowly for species with very low or high probability.
#' I am not sure what effect this has on the summary metrics.
#' 
#' @examples 
#' fit <- readRDS("../Experiments/7_4_modelrefinement/fittedmodels/7_4_13_model_2lv_e13.rds")
#' fit <- translatefit(fit)
#' Xocc <- unstandardise.designmatprocess(fit$XoccProcess, fit$data$Xocc[1:2, , drop = FALSE])
#' pocc <- poccupancy_margotherspecies.jsodm_lv(fit, Xocc)[,,"median"]
#' pocc2_all <- poccupy.jsodm_lv(fit, fullposterior = FALSE)
#' pocc2 <- pocc2_all[1, , , drop = TRUE]
#' pocc2 <- t(pocc2)
#' pocc <- pocc2
#' #pocc2 <- aperm(pocc2, perm = c(3, 1, 2)) 
#' #dim(pocc2) <- c(500 * 2, 60) #going iterating through row goes through draws, then through sites
#' 
#' 
#' traits_raw <- readRDS("../sflddata/private/data/raw/bird_traits.rds")
#' traits_raw <- traits_raw[, c("CommonName", 
#'                              primary_diet = "MFScore", 
#'                              foraging_substrate = "FMScore", 
#'                              feeding_aggregation = "SSScore", 
#'                              nesting_aggregation = "BScore", 
#'                              seasonal_movement = "MScore", 
#'                              body_size = "CubeRootBodyWeight")]
#' species <- colnames(pocc)
#' stopifnot(setequal(intersect(traits_raw$CommonName, species), species))
#' traits <- traits_raw[traits_raw$CommonName %in% species, ]
#' traits <- traits[complete.cases(traits), ] #remove incomplete rows
#' stopifnot(0 == sum(duplicated(traits$CommonName))) #stop if species are duplicated (i.e. common name used for species with different functionality)
#' stopifnot(setequal(species, traits$CommonName)) #stop if after the above cleaning, some species are not present in the traits data
#' 
#' rownames(traits) <- traits$CommonName
#' traits <- traits[, colnames(traits) != "CommonName"]
#' system.time(fd_av <- fd_pocc(traits, pocc, nsim = 1)) 
#' #27 seconds for 100 simulations, 60 birds and 2 sites.
#' #132 seconds for 500 simulations, 60 birds and 2 sites.
#' #27 seconds for 500 draws, 60 birds and 1 site. At this speed all 2000 insample sites would take 20 hours.
# profvis::profvis(fd_av <- fd_pocc(traits, pocc, nsim = 100))
# profvis::profvis(fd_av <- fd_pocc(traits, pocc, nsim = 1))
#' summary(fd_av$E)
#' 
#' @export
fd_pocc <- function(traits, pocc, nsim = 100){
  distmat <- FD::gowdis(traits)
  
  # simulate occupancy
  occ <- simulateocc(pocc, nsim = nsim)

  # remove any species not at any of the sites for the set of simulations
  occslim <- occ[, colSums(occ) > 0]
  distmat2 <- as.matrix(distmat)
  distmat2 <- distmat2[colnames(occslim), colnames(occslim)]
  
  fd <- FD::dbFD(distmat2, 
             a = occslim,
             w.abun = FALSE,
             calc.CWM = FALSE) # here the a indicates that the species is present
  fd_df <- as.data.frame(fd[names(fd) != "qual.FRic"])
  Efd <- vapply(1:nrow(pocc), function(site){
    Evalsite <- colMeans(fd_df[site + seq(1, nsim, by = nrow(pocc)), ])
    return(Evalsite)
  }, FUN.VALUE = as.numeric(fd_df[1, , drop = TRUE]))
  Efd <- t(Efd)
  
  Vfd <- vapply(1:nrow(pocc), function(site){
    valsite <- apply(fd_df[site + seq(1, nsim, by = nrow(pocc)), ], MARGIN = 2, FUN = sd)
    return(valsite)
  }, FUN.VALUE = as.numeric(fd_df[1, , drop = TRUE]))^2
  Vfd <- t(Vfd)
  return(list(E = Efd, V = Vfd))
}

#returns a matrix where pocc is row-binded nsim times and simulated. 
#For large nsim the column average should equal the column averages of pocc
#Each occupancy is considered independent.
simulateocc <- function(pocc, nsim = 1E5){
  # in the following ind is a set of iid uniform random variables.
  # If each is below the occupancy probablity threshold then the species will be deemed as in-occupancy.
  # In-occupancy occurs with probability equal to the threshold.
  ind <- runif(length(pocc) * nsim, min = 0, max = 1)
  dim(ind) <- c(dim(pocc)[[1]] * nsim, dim(pocc)[[2]])
  colnames(ind) <- colnames(pocc)
  occ_t <- as.vector(t(ind)) <= as.vector(t(pocc))   #t() is used because as.vector goes through rows first. Want it to finish with each row, before moving onto next one (rather than finish with each column)
  dim(occ_t) <- c(ncol(pocc), dim(pocc)[[1]] * nsim)
  occ <- t(occ_t)
  colnames(occ) <- colnames(pocc)
  return(occ)
}

# diffs <- vapply((1:100) * 1E3, function(nsim){
#   occ <- simulateocc(pocc, nsim = nsim)
#   maxdiff <- max(abs(colMeans(occ) - colMeans(pocc)))
#   return(maxdiff)
# }, FUN.VALUE = 0.07)
# plot(diffs)
# abline(h = 0)

    