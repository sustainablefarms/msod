
library(dplyr)

#' 
#' @param Xocc A dataframe of covariates related to occupancy. One row per ModelSite.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param yXobs A dataframe of species observations (1 or 0) and covariates related to observations. One row per visit.
#' Each column is either a covariate or a species.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param species A list of species names (corresponding to columns of yXobs) to model.
#' @param OccFmla A formula specifying the occupancy model in terms of the covariates in Xocc.
#' @param ObsFmla A formula specifying the model for detection in terms of the covariates in Xobs.
#' @param ModelSite A list of column names in y, Xocc and Xobs that uniquely specify the ModelSite. Can be simply a ModelSite index
#' @param nlv The number of latent variables
#' @param MCMCparams A named list containing values for n.chains, adapt, burnin, sample and thin (to pass to run.jags). 
#' Also include keep.jags.files to specify the directory that JAGS data will be saved.
#' @param filename If non-null the runjags object (with some extra information) is saved to filename as an RDS.
#' @return A runjags object with some modifications:
#' *   the data slot is a in list.format, converted using [as.list.format()]
#' *   elements in the data slot have dimension names given by the input data frames
#' *   the slot XoccProcess is the process used to prepare the occupancy covariates (scaling and centering). 
#' It can be applied using [apply.designmatprocess()]
#' *   the slot XobsProcess is the process used to prepare the detection covariates (scaling and centering). 
#' It can be applied using [apply.designmatprocess()]
#' *   ModelSite slot is the list with values for each visit (row in detection covariates) giving the row of the ModelSite in the occupancy covariates.
#' *   species slot is the list of species names.
#' @export
run.detectionoccupany <- function(Xocc, yXobs, species, ModelSite, OccFmla = "~ 1", ObsFmla = "~ 1", nlv = 2,
                                  initsfunction = defaultinitsfunction,
                                  MCMCparams = list(n.chains = 1, adapt = 2000, burnin = 25000, sample = 1000, thin = 30),
                                  filename = NULL){
  XoccProcess <- prep.designmatprocess(Xocc, OccFmla)
  XobsProcess <- prep.designmatprocess(yXobs, ObsFmla)
  
  #Specify the data
  data.list <- prep.data(Xocc = Xocc,
            yXobs = yXobs,
            ModelSite = ModelSite,
            species = species,
            nlv = nlv,
            XoccProcess = XoccProcess,
            XobsProcess = XobsProcess)
  if (nlv > 0){
    ### Latent variable multi-species co-occurence model
    modelFile='./scripts/7_2_1_model_description.txt'
    #Specify the parameters to be monitored
    monitor.params = c('u.b','v.b','mu.u.b','tau.u.b','mu.v.b','tau.v.b','lv.coef', "LV")
  } else {
    modelFile='./scripts/7_2_3_model_description_nolv.txt'
    #Specify the parameters to be monitored
    monitor.params = c('u.b','v.b','mu.u.b','tau.u.b','mu.v.b','tau.v.b')
  }


  # set up initial values
  if (is.null(MCMCparams$n.chains)){n.chains <- 1}
  else {n.chains <- MCMCparams$n.chains}
  inits <- lapply(1:n.chains, initsfunction, indata = data.list)

#### RUNNING JAGS ######
  runjagsargs <- list(
    model = modelFile,
    n.chains = n.chains,
    data = data.list,
    inits = inits,
    method = 'parallel',
    monitor = monitor.params,
    adapt = 2000,
    burnin = 25000,
    sample = 1000,
    thin = 30,
    keep.jags.files = FALSE,
    tempdir = FALSE
  )
  runjagsargs[names(MCMCparams)] <- MCMCparams
  
  
  library(runjags)
  fit.runjags <- do.call(run.jags, args = runjagsargs)
  
  # note that for simulation studies Tobler et al says they ran 3 chains drew 15000 samples, after a burnin of 10000 samples and thinning rate of 5.
  # In the supplementary material it appears these parameters were used: n.chains=3, n.iter=20000, n.burnin=10000, n.thin=10. Experiment 7_1 suggested a higher thinning rate

  # add summary of parameter distributions
  if (runjagsargs$sample >= 100) {
    fit.runjags <- add.summary(fit.runjags)
    fit.runjags$crosscorr <- "Crosscorrelation removed to conserve disk size. See ?add.summary to compute it."
    }
 
  # convert input data of model into nice format (saves a lot of computational to avoid the list.format operation in every argument) 
  fit.runjags$data <- list.format(fit.runjags$data)
  colnames(fit.runjags$data$y) <- species
  colnames(fit.runjags$data$Xocc) <- names(XoccProcess$center)
  rownames(fit.runjags$data$Xocc) <- 1:nrow(fit.runjags$data$Xocc)
  colnames(fit.runjags$data$Xobs) <- names(XobsProcess$center)
  rownames(fit.runjags$data$Xobs) <- data.list$ModelSite
  
  # attach data preparation methods
  fit.runjags$XoccProcess <- XoccProcess
  fit.runjags$XobsProcess <- XobsProcess
  fit.runjags$ModelSite <- data.list$ModelSite
  fit.runjags$species <- species
  if (!is.null(filename)){saveRDS(fit.runjags, filename) }
  invisible(fit.runjags)
}


#' @describeIn apply.designmatprocess Used to create instructions for processing data.
#' @param indata Input dataframe to be processed.
#' @param fmla A model formula (predictor side only).
#' @return A special list containing parameters for applying a preprocessing step to data.
#' @export
prep.designmatprocess <- function(indata, fmla){
  designmat1 <- model.matrix(as.formula(fmla), as.data.frame(indata))
  means <- colMeans(designmat1)
  sds <- ((nrow(designmat1) - 1) / nrow(designmat1)) * apply(designmat1, 2, sd)
  center <- means
  scale <- sds
  center[sds * means < 1E-8] <- means[sds * means < 1E-8] - 1 #centering of constant columns to 1
  scale[sds * means < 1E-8] <- 1 #no scaling of constant columns
  preprocessobj <- list(fmla = fmla, center = center, scale = scale)
  return(preprocessobj)
}
#' @param designmatprocess Are instructions for preprocessing input data, created by [prep.designmatprocess()] 
#' @param indata Input dataframe to be processed.
#' @details The input data is turned into a design matrix using [stats::model.matrix()].
#' Each non-constant column is then centered and scaled.
#' @return A design matrix.
#' @export
apply.designmatprocess <- function(designmatprocess, indata){
  designmat1 <- model.matrix(as.formula(designmatprocess$fmla), as.data.frame(indata))
  designmat <- scale(designmat1, center = designmatprocess$center, scale = designmatprocess$scale)
  return(designmat)
}


#' Given the input data parameters of run.detectionoccupancy prepare the data list for JAGS
#' @param XoccProcess An object create by prep.designmatprocess for the occupancy covariates
#' @param XobsProcess An object create by prep.designmatprocess for the observation covariates
prep.data <- function(Xocc, yXobs, ModelSite, species, nlv, XoccProcess, XobsProcess){
  # check data inputs
  stopifnot(all(ModelSite %in% colnames(Xocc)))
  stopifnot(all(ModelSite %in% colnames(yXobs)))
  stopifnot(all(species %in% colnames(yXobs)))
  stopifnot(anyDuplicated(Xocc[, ModelSite]) == 0) #model site uniquely specified
  stopifnot(all(as.matrix(yXobs[, species]) %in% c(1, 0)))
  
  # create model site indexes
  # if (ModelSiteID %in% c(names(yXobs), Xocc)){warning("Overwriting ModelSiteID column in input data.")}
  ModelSiteMat <- cbind(1:nrow(Xocc), tibble::as_tibble(Xocc[, ModelSite, drop = FALSE]))
  visitedModelSiteMat <- dplyr::right_join(ModelSiteMat, tibble::as_tibble(yXobs[, ModelSite, drop = FALSE]), by = ModelSite, suffix = c("", ".in"))
  visitedModelSite <- visitedModelSiteMat[, 1, drop = TRUE]
  stopifnot(is.integer(visitedModelSite))
  stopifnot(all(visitedModelSite <= nrow(Xocc)))
  
  XoccDesign <- apply.designmatprocess(XoccProcess, Xocc)
  XobsDesign <- apply.designmatprocess(XobsProcess, yXobs)
  
  n = length(species) #number of species
  J <- nrow(XoccDesign)  #number of unique sites should also be max(occ_covariates$SiteID)
  y <- as.matrix(yXobs[, species])
  ModelSite <- visitedModelSite
  data.list = list(n=n, J=J, y=y,
                  ModelSite = ModelSite, #a list of the site visited at each visit
                  Vvisits = nrow(XobsDesign), #number of visits in total - not sure what this is for
                  Xocc=XoccDesign,Xobs=XobsDesign,Vocc=ncol(XoccDesign),Vobs=ncol(XobsDesign),nlv=nlv)
  return(data.list)
}

#' A short function that applies the prep.data function to new data, given an object created by run.detectionoccupancy
#' Xocc, yXobs, ModelSite must follow some rules as for run.detectionoccupancy
prep_new_data <- function(fit, Xocc, yXobs, ModelSite){
  data.list <- prep.data(Xocc, yXobs, ModelSite, fit$species, fit$nlv, fit$XoccProcess, fit$XobsProcess)
  return(data.list)
}

### Initial conditions function
#Specify the initial values using a function
defaultinitsfunction <- function(chain, indata, ...) {
  if (!is.null(indata$nlv) && (indata$nlv > 0)){
    lv.coef<-matrix(1, indata$n, indata$nlv)
    lv.coef[1:indata$nlv,1:indata$nlv]<-0
    for(l in 1:indata$nlv-1){
      lv.coef[l,(l+1):indata$nlv]<-NA
    }
    LV <- matrix(rnorm(indata$nlv * indata$J), indata$J, indata$nlv)
  } else {
    lv.coef <- NULL
    LV <- NULL
  }
  v.b.proto <- lapply(colnames(indata$y),
                      function(x) {unname(coef(glm(((indata$y>0)*1)[, x] ~ . - 1, #intercept is built in
                                                   data = data.frame(indata$Xobs),
                                                   family=binomial(link=logit))))})
  v.b <- t(do.call(cbind, v.b.proto))
  v.b[v.b > 5] <- 5  #remove extremes as coefficients are assumed to be from a standard gaussian
  v.b[v.b < -5] <- -5  #remove extremes as coefficients are assumed to be from a standard gaussian

  ## this is calculated just to get initial values for occupancy covariates and occupancy estimates
  y.occ.mock <- cbind(ModelSiteID = indata$ModelSite, indata$y) %>%
    as_tibble() %>%
    group_by(ModelSiteID) %>%
    summarise_all(max) %>%
    dplyr::select(-ModelSiteID)
  u.b.proto <- lapply(colnames(y.occ.mock),
                      function(x) {unname(coef(glm( ((y.occ.mock>0)*1)[, x] ~ . - 1, #intercept is built in
                                                    data = data.frame(indata$Xocc),
                                                    family=binomial(link=probit))))})
  
  u.b <- t(do.call(cbind, u.b.proto))
  u.b[u.b > 5] <- 5  #remove extremes as coefficients are assumed to be from a standard gaussian
  u.b[u.b < -5] <- -5  #remove extremes as coefficients are assumed to be from a standard gaussian

  .RNG.seed <- c(1, 2, 3, 4, 5)[chain] # run jags likes to have this and the following argument defined in the initlalization functions.
  .RNG.name <- c(rep(c("base::Super-Duper", "base::Wichmann-Hill"),3))[chain]
  
  
  out <- list(
    u.b= u.b,  #initial values guestimated from u.b.proto are erroring! "u[14,1]: Node inconsistent with parents"
    v.b= v.b,
    u=(y.occ.mock>0)-runif(1,0.1,0.8),  #this looks strange -> step(u) is an indicator of whether occupied or not
    #mu.a = matrix(rbinom((n)*J, size=1, prob=1),
    #              nrow=J, ncol=(n)),
    #lv.coef=matrix(runif(nlv*n,-1,1),n,nlv)*lv.coef,
    LV= LV, 
    #z = matrix(rbinom((n)*J, size=1, prob=1),
    #           nrow=J, ncol=(n))
    .RNG.seed = .RNG.seed,
    .RNG.name = .RNG.name
  )
  if (is.null(indata$nlv) || (indata$nlv == 0)){
    out[["LV"]] <- NULL
  }
}


#### Examples #####
#' 
#' @examples 
#' inputdata <- readRDS("./private/data/clean/7_2_1_input_data.rds")
#' fitjags <- run.detectionoccupany(
#'   Xocc = inputdata$occ_covariates,
#'   yXobs = inputdata$plotsmerged_detection,
#'   species = inputdata$detection_data_specieslist,
#'   ModelSite = "ModelSiteID",
#'   OccFmla = "~ 1",
#'   ObsFmla = "~ MeanWind + 1",
#'   nlv = 2,
#'   MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 1, thin = 1,
#'                     keep.jags.files = "./runjags_smallB"),
#'   filename = "./runjags_demo.rds"
#' )
#' fitjags_nolv <- run.detectionoccupany(
#'   Xocc = inputdata$occ_covariates,
#'   yXobs = inputdata$plotsmerged_detection,
#'   species = inputdata$detection_data_specieslist,
#'   ModelSite = "ModelSiteID",
#'   OccFmla = "~ 1",
#'   ObsFmla = "~ MeanWind + 1",
#'   nlv = 0,
#'   MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 1, thin = 1,
#'                     keep.jags.files = "./runjags_demo"),
#'   filename = "./runjags_demo.rds"
#' )
