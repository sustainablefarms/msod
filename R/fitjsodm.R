#' @title Fit detection occupancy models using runjags.
#' @import dplyr

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
#' *   the data slot is a in list.format, converted using [as_list_format()]
#' *   elements in the data slot have dimension names given by the input data frames
#' *   the slot XoccProcess is the process used to prepare the occupancy covariates (scaling and centering). 
#' It can be applied using [apply.designmatprocess()]
#' *   the slot XobsProcess is the process used to prepare the detection covariates (scaling and centering). 
#' It can be applied using [apply.designmatprocess()]
#' *   ModelSite slot is the list with values for each visit (row in detection covariates) giving the row of the ModelSite in the occupancy covariates.
#' *   species slot is the list of species names.
#' @export
fitjsodm              <- function(Xocc, yXobs, species, ModelSite, modeltype,  #all these arguments are mandatory
                                  ..., #ellipsis for arguments specific to each model type
                                  OccFmla = "~ 1", ObsFmla = "~ 1",
                                  initsfunction = get0(paste0("paraminits.", modeltype)),
                                  MCMCparams = list(n.chains = 1, adapt = 2000, burnin = 25000, sample = 1000, thin = 30),
                                  filename = NULL){
  stopifnot(modeltype %in% availmodeltypes)
  
  if (!is.null(filename)){checkwritable(filename)} # check that file can be written before continuing
  XoccProcess <- prep.designmatprocess(Xocc, OccFmla)
  XobsProcess <- prep.designmatprocess(yXobs, ObsFmla)
  
  #Specify the data
  data.list <- prepJAGSdata(modeltype,
            Xocc = Xocc,
            yXobs = yXobs,
            ModelSite = ModelSite,
            species = species,
            XoccProcess = XoccProcess,
            XobsProcess = XobsProcess,
            ...)
  
  # parameters to monitor 
  monitor.params = c('occ.b','det.b','mu.occ.b', 'tau.occ.b','sigma.occ.b', 'mu.det.b','tau.det.b', 'sigma.det.b')
  if (modeltype == "jsodm"){
    modelFile=system.file("modeldescriptions",
			  "jsodm.txt",
			  package = "msod")
  }
  if (modeltype == "jsodm_lv"){
    ### Latent variable multi-species co-occurence model
    modelFile=system.file("modeldescriptions",
                          "jsodm_lv.txt",
                          package = "msod")
    #Specify the parameters to be monitored
    monitor.params = c(monitor.params, 'lv.b', "LV")
  }
  if (modeltype == "jsodm_lv_sepexp"){
    ### Latent variable multi-species co-occurence model
    modelFile=system.file("modeldescriptions",
                          "jsodm_lv_sepexp.txt",
                          package = "msod")
    #Specify the parameters to be monitored
    monitor.params = c(monitor.params, 'lv.b', "LV", "lv.v.spatscale", "lv.v.timescale")
  }

  # set up initial values
  if (is.null(MCMCparams$n.chains)){
    n.chains <- 1
  } else {
    n.chains <- MCMCparams$n.chains
  }
  inits <- lapply(1:n.chains, initsfunction, indata = data.list)

#### RUNNING JAGS ######
  runjagsargs <- list(
    model = modelFile,
    n.chains = n.chains,
    data = data.list,
    inits = inits,
    method = 'parallel',
    monitor = monitor.params,
    keep.jags.files = FALSE,
    tempdir = FALSE
  )
  runjagsargs[names(MCMCparams)] <- MCMCparams
  
  fit.runjags <- do.call(runjags::run.jags, args = runjagsargs)
  
  # note that for simulation studies Tobler et al says they ran 3 chains drew 15000 samples, after a burnin of 10000 samples and thinning rate of 5.
  # In the supplementary material it appears these parameters were used: n.chains=3, n.iter=20000, n.burnin=10000, n.thin=10. Experiment 7_1 suggested a higher thinning rate

  # add summary of parameter distributions
  if (runjagsargs$sample >= 100) {
    fit.runjags <- runjags::add.summary(fit.runjags)
    fit.runjags$crosscorr <- "Crosscorrelation removed to conserve disk size. See ?add.summary to compute it."
    }
 
  # replace data form with the data.list above
  # (saves a lot of computational to avoid the list.format operation in every argument) 
  fit.runjags$data <- runjags::list.format(fit.runjags$data)
  # check data is equal to data.list
  stopifnot(isTRUE(all.equal(fit.runjags$data, data.list, check.attributes = FALSE)))
  fit.runjags$data <- data.list
  
  # attach data preparation methods
  fit.runjags$XoccProcess <- XoccProcess
  fit.runjags$XobsProcess <- XobsProcess
  fit.runjags$species <- species
  
  if (modeltype == "jsodm_lv_sepexp"){
    fit.runjags$SpatDist <- SpatDist
    fit.runjags$TimeDist <- TimeDist
  }
  
  # add more things to save here (e.g. items that create random effects)
  
  class(fit.runjags) <- c(modeltype, class(fit.runjags))
  
  if (!is.null(filename)){try(saveRDS(fit.runjags, filename)) }
  
  ellipsis::check_dots_used()
  invisible(fit.runjags)
}


#' @describeIn fitjsodm An obsolete alias to fitjsodm()
#' @export
run.detectionoccupancy <- function(Xocc, yXobs, species, ModelSite, OccFmla = "~ 1", ObsFmla = "~ 1",
                                  modeltype, ...,
                                  initsfunction = get0(paste0("paraminits.", modeltype)),
                                  MCMCparams = list(n.chains = 1, adapt = 2000, burnin = 25000, sample = 1000, thin = 30),
                                  filename = NULL){
  warning("This function is obsolete, consider using fitjsodm instead.")
  fit <- fitjsodm(Xocc, yXobs, species, ModelSite, modeltype,
                  ...,
                  OccFmla = OccFmla, ObsFmla = ObsFmla,
                  initsfunction = initsfunction,
                  MCMCparams = MCMCparams,
                  filename = filename)
  return(fit)
}



# x is filename
checkwritable <- function(x){
  if (!file.exists(x)){
    saveRDS(NA, x)
    file.remove(x)
  } else {
    stopifnot(file.access(x, mode = 2) == 0)
  }
  return(invisible(x))
}

#### Examples #####
#' 
#' @examples 
#' inputdata <- readRDS(readLines("../Experiments/link_7_2_10_input_data.txt")[[1]])
# spatdistfun <- function(df) {as.matrix(dist(data.frame(x = df$latitude, y = df$longitude)))}
# timedistfun <- function(df) {as.matrix(dist(data.frame(z = df$SurveyYear)))}
# fitjags <- fitjsodm(
#   Xocc = inputdata$insampledata$Xocc[1:100, ],
#   yXobs = inputdata$insampledata$yXobs[inputdata$insampledata$yXobs$ModelSiteID <= 100, ],
#   species = inputdata$species,
#   ModelSite = "ModelSiteID",
#   modeltype = "jsodm_lv_sepexp",
#   OccFmla = "~ 1",
#   ObsFmla = "~ 1",
#   nlv = 2,
#   SpatDist = spatdistfun,
#   TimeDist = timedistfun,
#   MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 1, thin = 1),
#   filename = "./runjags_demo.rds"
# )
#' fitjags <- fitjsodm(
#'   Xocc = inputdata$insampledata$Xocc[1:100, ],
#'   yXobs = inputdata$insampledata$yXobs[inputdata$insampledata$yXobs$ModelSiteID <= 100, ],
#'   species = inputdata$species,
#'   ModelSite = "ModelSiteID",
#'   modeltype = "jsodm_lv",
#'   OccFmla = "~ 1",
#'   ObsFmla = "~ 1",
#'   nlv = 2,
#'   MCMCparams = list(n.chains = 1, adapt = 0, burnin = 0, sample = 1, thin = 1),
#'   filename = "./runjags_demo.rds"
#' )
#' fitjags_nolv <- run.detectionoccupancy(
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
