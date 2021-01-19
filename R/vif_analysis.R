#' @examples
#' inputdata <- readRDS("./private/data/clean/7_2_10_input_data.rds")
#' 
#' OccFmla <- reformulate(paste0("\`",setdiff(names(inputdata$insample$Xocc), c("ModelSiteID", "SiteCode", "SurveySiteId", "StudyId", "holdout", "gc", "AnnTempRange")), "\`"))
#' ObsFmla <- "~ 1 + MeanWind + MeanTime + MeanClouds + MeanTemp"
#' vifs <- gvifs_lm(inputdata$insample$Xocc,
#'          inputdata$insample$yXobs,
#'          species = inputdata$species,
#'          ModelSite = "ModelSiteID", 
#'          OccFmla = OccFmla,
#'          ObsFmla = ObsFmla
#'          )
#' 
#' vifs[[1]] %>%
#'   tibble::enframe(name = "var", value = "vif") %>%
#'   arrange(vif) %>%
#'   ggplot() +
#'   geom_point(aes(y = var, x = vif)) +
#'   geom_vline(xintercept = 30)

#' @title Compute VIF using CAR package
#' @param Xocc A dataframe of covariates related to occupancy. One row per ModelSite.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param yXobs A dataframe of species observations (1 or 0) and covariates related to observations. One row per visit.
#' Each column is either a covariate or a species.
#' Must also include the ModelSiteVars columns to uniquely specify ModelSite.
#' @param species A list of species names (corresponding to columns of yXobs) to model.
#' @param OccFmla A formula specifying the occupancy model in terms of the covariates in Xocc.
#' @param ObsFmla A formula specifying the model for detection in terms of the covariates in Xobs.
#' @param ModelSite A list of column names in y, Xocc and Xobs that uniquely specify the ModelSite. Can be simply a ModelSite index
#' @export
gvifs_lm <- function(Xocc, yXobs, species, ModelSite, OccFmla = "~ 1", ObsFmla = "~ 1"){
  XoccProcess <- prep.designmatprocess(Xocc, OccFmla)
  XobsProcess <- prep.designmatprocess(yXobs, ObsFmla)

  #prepare the data
  data.list <- prepJAGSdata("jsodm",
                         Xocc = Xocc,
                         yXobs = yXobs,
                         ModelSite = ModelSite,
                         species = species,
                         XoccProcess = XoccProcess,
                         XobsProcess = XobsProcess)

  #prepare to compute gvifs using mock version of occupancy
  y.occ.mock_bysp <- cbind(ModelSiteID = data.list$ModelSite, data.list$y) %>%
    tibble::as_tibble() %>%
    dplyr::group_by(ModelSiteID) %>%
    dplyr::summarise_all(max) %>%
    dplyr::select(-ModelSiteID)
  # y.occ.mock_sum <- rowSums(y.occ.mock_bysp)

  # for each species fit a linear model
  # alldata <- data.frame(y.occ.mock_bysp, data.list$Xocc, check.names = FALSE)
  # fmlas <- paste0("`", species, "`", " ~ . -`(Intercept)`")
  # names(fmlas) <- species
  mods <- lapply(species, function(x){
    fmla <- paste0("`", x, "`", " ~ . -`(Intercept)`")
    alldata <- data.frame(y.occ.mock_bysp[, x, drop = FALSE], data.list$Xocc, check.names = FALSE)
    return(lm(fmla, alldata))
  })
  stopifnot(length(mods) > 0)
  names(mods) <- species
                 
  
  # check for NA fitted coefficients
  coefs <- simplify2array(lapply(mods, coefficients)) #rows are the coefficients, each model is a column
  varhasna <- apply(coefs, MARGIN = 1, anyNA)
  if (sum(varhasna) > 0){
    stop(paste(names(varhasna[varhasna]), " is fitted as NA. Please modify model to compute VIFs."))
  }
  
  # for each of these models compute the generalised variance inflation factors
  gvifs <- lapply(mods, car::vif)
  gvifs <- simplify2array(gvifs)
  mingvifs <- apply(gvifs, MARGIN = 1, min)
  return(list(
    min_vif = mingvifs,
    vif = gvifs
  ))
}
