#' Cross-sectional FoSR using GLS
#' 
#' Fitting function for function-on-scalar regression for cross-sectional data.
#' This function estimates model parameters using GLS: first, an OLS estimate of 
#' spline coefficients is estimated; second, the residual covariance is estimated
#' using an FPC decomposition of the OLS residual curves; finally, a GLS estimate
#' of spline coefficients is estimated. Although this is in the `BayesFoSR` package,
#' there is nothing Bayesian about this FoSR.
#' 
#' @param formula a formula indicating the structure of the proposed model. 
#' @param Kt number of spline basis functions used to estimate coefficient functions
#' @param data an optional data frame, list or environment containing the 
#' variables in the model. If not found in data, the variables are taken from 
#' environment(formula), typically the environment from which the function is 
#' called.
#' @param basis basis type; options are "bs" for b-splines and "pbs" for periodic
#' b-splines
#' @param verbose logical defaulting to \code{TRUE} -- should updates on progress be printed?
#' @param sigma optional covariance matrix used in GLS; if \code{NULL}, OLS will be
#' used to estimated fixed effects, and the covariance matrix will be estimated from
#' the residuals.
#' @param CI.type Indicates CI type for coefficient functions; options are "pointwise" and
#' "simultaneous"
#'  
#' @references
#' Goldsmith, J., Kitago, T. (Under Review).
#' Assessing Systematic Effects of Stroke on Motor Control using Hierarchical 
#' Function-on-Scalar Regression.
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @importFrom splines bs
#' @importFrom pbs pbs
#' @export
#' 
#' @examples
#' 
#' library(refund)
#' library(dplyr)
#' 
#' data(DTI)
#' DTI = subset(DTI, select = c(cca, case, pasat))
#' DTI = DTI[complete.cases(DTI),]
#' DTI$gender = factor(sample(c("male","female"), dim(DTI)[1], replace = TRUE))
#' DTI$status = factor(sample(c("RRMS", "SPMS", "PPMS"), dim(DTI)[1], replace = TRUE))
#' 
#' fosr.dti = fosr_gls(cca ~ pasat * gender + status, data = DTI)
#' 
gls_cs = function(formula, data=NULL, Kt=5, basis = "bs", sigma = NULL, verbose = TRUE, CI.type = "pointwise"){
  
  # not used now but may need this later
  call <- match.call()
  
  tf <- terms.formula(formula, specials = "re")
  trmstrings <- attr(tf, "term.labels")
  specials <- attr(tf, "specials")    # if there are no random effects this will be NULL
  where.re <-specials$re - 1
  
  # gets matrix of fixed and random effects
  if(length(where.re)!=0){
    mf_fixed <- model.frame(tf[-where.re], data = data)
    formula = tf[-where.re]
    
    # get random effects matrix
    responsename <- attr(tf, "variables")[2][[1]]
    REs = eval(parse(text=attr(tf[where.re], "term.labels")))
    
    # set up dataframe if data = NULL
    formula2 <- paste(responsename, "~", REs[[1]],sep = "")
    newfrml <- paste(responsename, "~", REs[[2]],sep = "")
    newtrmstrings <- attr(tf[-where.re], "term.labels")
    
    formula2 <- formula(paste(c(formula2, newtrmstrings), collapse = "+"))
    newfrml <- formula(paste(c(newfrml, newtrmstrings), collapse = "+"))
    mf <- model.frame(formula2, data = data)
    
    # creates the Z matrix. get rid of $zt if you want a list with more stuff.
    if(length(data)==0){Z = lme4::mkReTrms(lme4::findbars(newfrml),fr=mf)$Zt
    }else
    {Z = lme4::mkReTrms(lme4::findbars(newfrml),fr=data)$Zt}
    
    
  } else {
    mf_fixed <- model.frame(tf, data = data)
  }
  mt_fixed <- attr(mf_fixed, "terms")
  
  # get response (Y)
  Y <- model.response(mf_fixed, "numeric")
  
  # x is a matrix of fixed effects
  # automatically adds in intercept
  X <- model.matrix(mt_fixed, mf_fixed, contrasts)
  
  ### model organization ###
  D = dim(Y)[2]
  I = dim(X)[1]
  p = dim(X)[2]
  
  if(basis == "bs"){
    Theta = bs(1:D, df = Kt, intercept=TRUE, degree=3)
  } else if(basis == "pbs"){
    Theta = pbs(1:D, df = Kt, intercept=TRUE, degree=3)
  }

  X.des = X
  Y.vec = as.vector(t(Y)) 
  X = kronecker(X.des, Theta)
  n.coef = dim(X.des)[2]
  
  if(is.null(sigma)){
    ## OLS model fitting and processing results
    if(verbose) { cat("Using OLS to estimate residual covariance \n") }
    model.ols = lm(Y.vec ~ -1 + X)
    Bx.ols = matrix(model.ols$coef, nrow = Kt, ncol = n.coef)  
    beta.hat.ols = t(Bx.ols) %*% t(Theta)
    
    resid.mat = matrix(resid(model.ols), I, D, byrow = TRUE)
    
    ## Get Residual Structure using FPCA
    ## note: this is commented out because, in simulations based on the headstart data, 
    ## using FPCA lead to higher-than-nominal sizes for tests of nested models. 
    ## using the raw covariance worked better. using FPCA is possible, but relies
    ## on some case-specific choices.
    # raw.resid.cov = cov(resid.mat)
    # fpca.resid = fpca.sc(resid.mat, pve = .9995, nbasis = 20)
    # resid.cov = with(fpca.resid, efunctions %*% diag(evalues) %*% t(efunctions))
    
    ## account for (possibly non-constant) ME nugget effect
    # sm.diag = Theta %*% solve(crossprod(Theta)) %*% t(Theta) %*% (diag(raw.resid.cov) - diag(resid.cov))
    # if(sum( sm.diag < 0 ) >0) { sm.diag[ sm.diag < 0] = min((diag(raw.resid.cov) - diag(resid.cov))[ sm.diag < 0])}
    # diag(resid.cov) = diag(resid.cov) + sm.diag
    
    sigma = cov(resid.mat) * (I - 1) / (I - p)

  }
  
  ## GLS fit through prewhitening
  if(verbose) { cat("GLS \n") }
  
  S = chol(solve(sigma))
  Y.t = t(Y)
  Z = as.vector(S %*% Y.t)
  T = S %*% Theta
  M = kronecker (X.des,  T)
  model = lm(Z ~ -1 + M)

  ## process results
  Bx = matrix(model$coef, nrow = Kt, ncol = n.coef)  
  beta.hat = t(Bx) %*% t(Theta)
  re = model$residuals
  Re = matrix(re, I, D, byrow = TRUE)
  cov<-vcov(model)
  
  ## get confidence intervals
  beta.UB = beta.LB = matrix(NA, p, D)
  wald.val = rep(NA, p)
  for(p.cur in 1:p){
    a = Kt*p.cur-(Kt-1)
    b = Kt*p.cur
    cov.cur = Theta %*% cov[a:b,a:b] %*%t(Theta)
    if(CI.type == "pointwise"){
      beta.UB[p.cur,] = beta.hat[p.cur,] + 1.96 * sqrt(diag(cov.cur))
      beta.LB[p.cur,] = beta.hat[p.cur,] - 1.96 * sqrt(diag(cov.cur))
    } else if(CI.type == "simultaneous") {
      norm.samp = mvrnorm(2500, mu = rep(0, D), Sigma = cov.cur)/
        matrix(sqrt(diag(cov.cur)), nrow = 2500, ncol = D, byrow = TRUE)
      crit.val = quantile(apply(abs(norm.samp), 1, max), .95)    
      beta.UB[p.cur,] = beta.hat[p.cur,] + crit.val * sqrt(diag(cov.cur))
      beta.LB[p.cur,] = beta.hat[p.cur,] - crit.val * sqrt(diag(cov.cur))
    }
    wald.val[p.cur] = Bx[,p.cur] %*% solve(cov[a:b,a:b]) %*% Bx[,p.cur]
  }
  
  Yhat = X.des %*% beta.hat
  
  ret = list(beta.hat, beta.UB, beta.LB, Yhat, mt_fixed, data, model, sigma, wald.val)
  names(ret) = c("beta.hat", "beta.UB", "beta.LB", "Yhat", "terms", "data", "model.gls", "sigma", "wald.val")
  class(ret) = "fosr"
  ret
    
}

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################