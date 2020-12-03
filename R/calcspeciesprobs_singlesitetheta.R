#' @title Probability of detection and occupation for a single ModelSite and parameter set
#' @param Xocc A matrix of occupancy covariates. Must have a single row. Columns correspond to covariates.
#' @param Xobs A matrix of detection covariates, each row is a visit. 
#' @param det.b Covariate loadings. Each row is a species, each column a detection covariate
#' @return A matrix of detection probabilities for each row of Xobs and each species.
#' Rows in returned value correspond to rows in Xobs, columns correspond to species.
pdetection_occupied.ModelSite.theta <- function(Xobs, det.b){
  Xobs <- as.matrix(Xobs)
  Detection.LinPred <- as.matrix(Xobs) %*% t(det.b)
  Detection.Pred.Cond <- exp(Detection.LinPred) / (exp(Detection.LinPred) + 1)   #this is the inverse logit function
  return(Detection.Pred.Cond)
}

#' @describeIn pdetection_occupied.ModelSite.theta The probability of occupation for given LV values.
#' @param occ.b Covariate loadings for occupancy. Each row is a species, each column an occupancy covariate.
#' @param lv.b Loadings for the latent variables. Each row is a species, each column corresponds to a LV.
#' @param LVvals A matrix of LV values. Each column corresponds to a LV. To condition on specific LV values, provide a matrix of row 1.
#' @return A matrix of occupancy probabilities. Each row corresponds to a row of LVvals, each column to a species.
poccupy.ModelSite.theta <- function(Xocc, occ.b, lv.b = NULL, LVvals = NULL){
  # external
  ModelSite.Occ.eta_external <- as.matrix(Xocc) %*% t(occ.b) #columns are species
  
  if (!is.null(lv.b)){# probability of occupancy given LV
    sd_u_condlv <- sqrt(1 - rowSums(lv.b^2)) #for each species the standard deviation of the indicator random variable 'occ.indicator', conditional on values of LV
    ModelSite.Occ.eta_LV <- LVvals %*% t(lv.b) #occupancy contribution from latent variables, performed all together
    ModelSite.Occ.eta <- Rfast::eachrow(ModelSite.Occ.eta_LV, ModelSite.Occ.eta_external, oper = "+") #add the external part to each simulation
    # Make occ.indicator standard deviations equal to 1 by dividing other values by sd
    # P(occ.indicator < -ModelSite.Occ.eta) = P(occ.indicator / sd < -ModelSite.Occ.eta / sd) = P(occ.v < -ModelSite.Occ.eta / sd)
    ModelSite.Occ.eta_standardised <- Rfast::eachrow(ModelSite.Occ.eta, sd_u_condlv, oper = "/")
  } else {
    ModelSite.Occ.eta_standardised <- ModelSite.Occ.eta_external
  }
  ModelSite.Occ.Pred.CondLV <- 1 - pnorm(-ModelSite.Occ.eta_standardised, mean = 0, sd = 1)
  return(ModelSite.Occ.Pred.CondLV)
}
