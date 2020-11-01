# the code in this file has largely been taken from boral's [https://cran.r-project.org/package=boral] get.envir.cor and get.residual.cor
## Produce the correlation due to similarity of responses to X
get.enviro.cor <- function(object, est = "median", prob = 0.95) 
{
  
  if(is.null(object$jags.model)) 
    stop("MCMC samples not found")
  fit.mcmc <- get.mcmcsamples(object)
  y <- object$y
  X <- object$X
  
  if(length(grep("X.coefs", colnames(fit.mcmc))) == 0) 
    stop("Cannot find MCMC sample corresponding to coefficients for X")
  
  n <- nrow(y); p <- ncol(y)
  enviro_cor_mat <- enviro_cor_mat_cilower <- enviro_cor_mat_ciupper <- enviro_cov_mat <- matrix(0,p,p)
  sig_enviro_cor_mat <- matrix(0,p,p)
  if(is.null(colnames(y))) 
    colnames(y) <- 1:ncol(y)
  rownames(enviro_cor_mat) <- rownames(enviro_cor_mat_cilower) <- rownames(enviro_cor_mat_ciupper) <- rownames(enviro_cov_mat) <- rownames(sig_enviro_cor_mat) <- colnames(y)
  colnames(enviro_cor_mat) <- colnames(enviro_cor_mat_cilower) <- colnames(enviro_cor_mat_ciupper) <- colnames(enviro_cov_mat) <- colnames(sig_enviro_cor_mat) <- colnames(y)
  all_enviro_cov_mat <- all_enviro_cor_mat <- array(0, dim = c(nrow(fit.mcmc),p,p))
  
  
  for(k0 in 1:nrow(fit.mcmc)) 
  {
    cw_X_coefs <- matrix(fit.mcmc[k0,grep("X.coefs", colnames(fit.mcmc))], nrow=p)
    enviro.linpreds <- tcrossprod(X,as.matrix(cw_X_coefs))
    all_enviro_cov_mat[k0,,] <- cov(enviro.linpreds)
    all_enviro_cor_mat[k0,,] <- cor(enviro.linpreds) 
  }
  
  for(j in 1:p) { for(j2 in 1:p) 
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
    get.hpd.cors <- HPDinterval(as.mcmc(all_enviro_cor_mat[,j,j2]), prob = prob)
    enviro_cor_mat_cilower[j,j2] <- get.hpd.cors[1] 
    enviro_cor_mat_ciupper[j,j2] <- get.hpd.cors[2]
    if(0 > get.hpd.cors[1] & 0 < get.hpd.cors[2]) 
      sig_enviro_cor_mat[j,j2] <- 0
  } }
  
  #return(list(residual.correlations = enviro_cor_mat))
  #corrplot(enviro_cor_mat, title = "Environmental correlations", type = "lower")
  return(list(cor = enviro_cor_mat, cor.lower = enviro_cor_mat_cilower, cor.upper = enviro_cor_mat_ciupper, sig.cor = sig_enviro_cor_mat, cov = enviro_cov_mat))
}
