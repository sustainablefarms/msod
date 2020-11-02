# the code in this file has largely been taken from boral's get.envir.cor and get.residual.cor

#' @title Compute inter-species correlations due to similarity of the occupancy linear predictor at model sites
#' @param fit An object created by [run.detectionoccupancy()]
#' @param est Stastical summary of draws to return.
#' @param prob The target probability coverage of highest probability density intervals.
#' @details 
#' For each draw from the posterior, 
#' computes the correlation and covariance between the species-specific occupancy linear predictions without latent variables across all sites.
#' The `est`, for each species, is returned.
#' Significantly non-zero covariance or correlation is computed by [coda::HPDinterval()] according to the signficance level `prob`.
#' @value
#' A list of matrices. Column and row names are species names.
#' `cor` is the correlation estimate computed according to `est`
#' `cor.lower` is the lower bound of the highest posterior density interval for correlation
#' `cor.upper` is the upper bound of the highest posterior density interval for correlation
#' `sig.cor` is `cor` where all correlations with hpd that includes 0 are set to 0.
#' `cov` is the covariance estimate computed accordting to `est`
#' @examples 
fit <- readRDS("../sflddata/private/data/testdata/cutfit_7_4_11_2LV.rds")
cor_envir <- get.enviro.cor(fit, est = "median", prob = 0.95)

#' @export
get.enviro.cor <- function(fit, est = "median", prob = 0.95) 
{
  nsites <- nrow(fit$data$Xocc)
  nsp <- ncol(fit$data$y)
  spnames <- colnames(fit$data$y)
  draws <- do.call(rbind, fit$mcmc)
  
  enviro_cor_mat <- enviro_cor_mat_cilower <- enviro_cor_mat_ciupper <- enviro_cov_mat <- matrix(0,nsp,nsp)
  sig_enviro_cor_mat <- matrix(0,nsp,nsp)

  rownames(enviro_cor_mat) <- rownames(enviro_cor_mat_cilower) <- rownames(enviro_cor_mat_ciupper) <- rownames(enviro_cov_mat) <- rownames(sig_enviro_cor_mat) <- spnames
  colnames(enviro_cor_mat) <- colnames(enviro_cor_mat_cilower) <- colnames(enviro_cor_mat_ciupper) <- colnames(enviro_cov_mat) <- colnames(sig_enviro_cor_mat) <- spnames
  all_enviro_cov_mat <- all_enviro_cor_mat <- array(0, dim = c(nrow(draws),nsp,nsp))
  
  
  for(k0 in 1:nrow(draws)) 
  { ## compute correlations per draw
    cw_X_coefs <- bugsvar2matrix(draws[k0, ], "u.b", 1:nsp, 1:ncol(fit$data$Xocc))
    enviro.linpreds <- tcrossprod(fit$data$Xocc,as.matrix(cw_X_coefs))
    #enviro.linpreds is the occupancy linear predictor without LV for posterior draw k0, for each species and ModelSite.
    #each row is a ModelSite, each column is a species
    all_enviro_cov_mat[k0,,] <- cov(enviro.linpreds)
    all_enviro_cor_mat[k0,,] <- cor(enviro.linpreds) 
  }
  
  for(j in 1:nsp) { for(j2 in 1:nsp) 
  { ## Average/Median over the MCMC samples
    if(est == "median") 
    { 
      enviro_cov_mat[j,j2] <- median(all_enviro_cov_mat[,j,j2])
      enviro_cor_mat[j,j2] <- median(all_enviro_cor_mat[,j,j2]) 
    }
    if(est == "mean") 
    {
      enviro_cov_mat[j,j2] <- mean(all_enviro_cov_mat[,j,j2])
      enviro_cor_mat[j,j2] <- mean(all_enviro_cor_mat[,j,j2]) 
    } 
    
    sig_enviro_cor_mat[j,j2] <- enviro_cor_mat[j,j2]
    get.hpd.cors <- coda::HPDinterval(coda::as.mcmc(all_enviro_cor_mat[,j,j2]), prob = prob)
    enviro_cor_mat_cilower[j,j2] <- get.hpd.cors[1] 
    enviro_cor_mat_ciupper[j,j2] <- get.hpd.cors[2]
    if(0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) 
      sig_enviro_cor_mat[j,j2] <- 0
  } }
  
  #return(list(residual.correlations = enviro_cor_mat))
  #corrplot(enviro_cor_mat, title = "Environmental correlations", type = "lower")
  return(list(cor = enviro_cor_mat, cor.lower = enviro_cor_mat_cilower, cor.upper = enviro_cor_mat_ciupper, sig.cor = sig_enviro_cor_mat, cov = enviro_cov_mat))
}

## Produce the residual correlation based on latent variables
get.residual.cor <- function(object, est = "median", prob = 0.95) 
{
  if(is.null(object$jags.model)) 
    stop("MCMC samples not found")
  
  fit.mcmc <- get.mcmcsamples(object)
  y <- object$y
  X <- object$X
  num.lv <- object$num.lv
  
  if(length(grep("lvs", colnames(fit.mcmc))) == 0) 
    stop("Cannot find MCMC samples corresponding to latent variables")
  
  n <- nrow(y)
  p <- ncol(y)
  sig_rescor_mat <- rescor_mat <- rescor_mat_cilower <- rescor_mat_ciupper <- rescov_mat <- matrix(0, nrow=p, ncol=p)
  sig_respres_mat <- respres_mat <- respres_mat_cilower <- respres_mat_ciupper <- matrix(0, nrow=p, ncol=p)
  if(is.null(colnames(y))) 
    colnames(y) <- 1:ncol(y); 
  rownames(rescor_mat) <- rownames(rescor_mat_cilower) <- rownames(rescor_mat_ciupper) <- rownames(sig_rescor_mat) <- colnames(rescor_mat) <- colnames(rescor_mat_cilower) <- colnames(rescor_mat_ciupper) <- colnames(sig_rescor_mat) <- colnames(y)
  rownames(respres_mat) <- rownames(respres_mat_cilower) <- rownames(respres_mat_ciupper) <- rownames(sig_respres_mat) <- colnames(respres_mat) <- colnames(respres_mat_cilower) <- colnames(respres_mat_ciupper) <- colnames(sig_respres_mat) <- colnames(y)
  rownames(rescov_mat) <- colnames(rescov_mat) <- colnames(y)
  all_rescor_mat <- all.rescov_mat <- all_respres_mat <- array(0, dim = c(nrow(fit.mcmc),p,p))
  all_trace_rescor <- numeric(nrow(fit.mcmc))
  
  for(k0 in 1:nrow(fit.mcmc)) 
  {
    lv.coefs <- matrix(fit.mcmc[k0,grep("lv.coefs", colnames(fit.mcmc))],nrow=p)
    # 		if(all(object$family == "binomial") & all(object$trial.size == 1)) 
    # 			lv.coefs[,2:(num.lv+1)] <- lv.coefs[,2:(num.lv+1)]/matrix(sqrt(1-rowSums(lv.coefs[,2:(num.lv+1)]^2)),nrow=p,ncol=num.lv,byrow=FALSE) ## If data is Bernoulli, then scale the coefficients to acocunt for constraints (see Knott and Bartholomew, Chapter 4)
    
    lambdalambdaT <- tcrossprod(as.matrix(lv.coefs[,2:(num.lv+1)]))
    all.rescov_mat[k0,,] <- (lambdalambdaT) 
    all_trace_rescor[k0] <- sum(diag(lambdalambdaT))
    
    #  		if(all(object$family == "negative.binomial")) {
    #    			get.var.phis <- numeric(p); 
    #    			## Multiplicative Poisson gamma model implies a log gamma random effect on the linear predictors
    #    			for(j in 1:p) 
    # 				get.var.phis[j] <- var(log(rgamma(2000,shape=1/lv.coefs[j,ncol(lv.coefs)],rate=1/lv.coefs[j,ncol(lv.coefs)])))
    # 			all.rescov_mat[k0,,] <- lambdalambdaT + diag(x=get.var.phis,nrow=p)
    # 			}
    all_rescor_mat[k0,,] <- cov2cor(all.rescov_mat[k0,,]) 
    all_respres_mat[k0,,] <- cor2pcor(lambdalambdaT)
  }
  
  for(j in 1:p) { for(j2 in 1:p) 
  { ## Average/Median over the MCMC samples
    if(est == "median") 
    { 
      rescov_mat[j,j2] <- median(all.rescov_mat[,j,j2]) 
      rescor_mat[j,j2] <- median(all_rescor_mat[,j,j2]); 
      respres_mat[j,j2] <- median(all_respres_mat[,j,j2]); 
    }
    
    if(est == "mean") 
    { 
      rescov_mat[j,j2] <- mean(all.rescov_mat[,j,j2]) 
      rescor_mat[j,j2] <- mean(all_rescor_mat[,j,j2]); 
      respres_mat[j,j2] <- mean(all_respres_mat[,j,j2]); 
    }
    
    sig_rescor_mat[j,j2] <- rescor_mat[j,j2]
    get.hpd.cors <- HPDinterval(as.mcmc(all_rescor_mat[,j,j2]), prob = prob)
    rescor_mat_cilower[j,j2] <- get.hpd.cors[1]
    rescor_mat_ciupper[j,j2] <- get.hpd.cors[2]
    if(0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) 
      sig_rescor_mat[j,j2] <- 0
    
    sig_respres_mat[j,j2] <- respres_mat[j,j2]
    get.hpd.cors <- HPDinterval(as.mcmc(all_respres_mat[,j,j2]), prob = prob)
    respres_mat_cilower[j,j2] <- get.hpd.cors[1]
    respres_mat_ciupper[j,j2] <- get.hpd.cors[2]
    if(0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) 
      sig_respres_mat[j,j2] <- 0
  } }
  
  if(est == "median") 
    final_trace <- median(all_trace_rescor)
  if(est == "mean") 
    final_trace <- mean(all_trace_rescor) 	
  
  return(list(cor = rescor_mat, cor.lower = rescor_mat_cilower, cor.upper = rescor_mat_ciupper, sig.cor = sig_rescor_mat, cov = rescov_mat, prec = respres_mat, prec.lower = respres_mat_cilower, prec.upper = respres_mat_ciupper, sig.prec = sig_respres_mat, trace = final_trace))
}

