#' Multilevel FoSR using a Gibbs sampler and Wishart prior
#' 
#' Fitting function for function-on-scalar regression for multilevel data.
#' This function estimates model parameters using a Gibbs sampler and estimates
#' the residual covariance surface using a Wishart prior.
#' 
#' @param formula a formula indicating the structure of the proposed model. 
#' @param Kt number of spline basis functions used to estimate coefficient functions
#' @param data an optional data frame, list or environment containing the 
#' variables in the model. If not found in data, the variables are taken from 
#' environment(formula), typically the environment from which the function is 
#' called.
#' @param N.iter number of iterations used in the Gibbs sampler
#' @param N.burn number of iterations discarded as burn-in
#' @param alpha tuning parameter balancing second-derivative penalty and
#' zeroth-derivative penalty (alpha = 0 is all second-derivative penalty)
#' 
#' @references
#' Goldsmith, J., Kitago, T. (Under Review).
#' Assessing Systematic Effects of Stroke on Motor Control using Hierarchical 
#' Function-on-Scalar Regression.
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @importFrom splines bs
#' @importFrom MCMCpack riwish
#' @export
#' 
gibbs_mult_wish = function(formula, Kt=5, data=NULL, N.iter = 5000, N.burn = 1000, alpha = .1){

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
  
  set.seed(1)
  
  ## fixed and random effect design matrices
  W.des = X
  Z.des = t(as.matrix(Z))
  
  I = dim(Z.des)[2]
  D = dim(Y)[2]
  Ji = as.numeric(apply(Z.des, 2, sum))
  IJ = sum(Ji)
  p = dim(W.des)[2]

  ## bspline basis and penalty matrix
  Theta = bs(1:D, df=Kt, intercept=TRUE, degree=3)

  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2

  SUBJ = factor(apply(Z.des %*% 1:dim(Z.des)[2], 1, sum))

  ## find first observation
  firstobs = rep(NA, length(unique(SUBJ)))
  for(i in 1:length(unique(SUBJ))){
    firstobs[i] = which(SUBJ %in% unique(SUBJ)[i])[1]
  }
  Wi = W.des[firstobs,]

  ## data organization; these computations only need to be done once
  Y.vec = as.vector(t(Y))
  IIP = kronecker(diag(1, I, I), P.mat)
  WIk = kronecker(Wi, diag(1, Kt, Kt))
  tWIW = t(WIk) %*% IIP %*% WIk
  tWI = t(WIk) %*% IIP


  ## initial estimation and hyperparameter choice
  vec.bz = solve(kronecker(t(Z.des)%*% Z.des, t(Theta) %*% Theta)) %*% t(kronecker(Z.des, Theta)) %*% Y.vec
  bz = matrix(vec.bz, nrow = Kt, ncol = I)
  
  w.temp = kronecker(t(Wi), diag(1, Kt, Kt))
  vec.bw = solve(tWIW) %*% tWI %*% vec.bz
  bw = matrix(vec.bw, Kt, p)

  Yhat = as.matrix(Z.des %*% t(bz) %*% t(Theta))
  varhat = var(as.vector(Y - Yhat))

  Psi = diag(varhat*IJ, D, D)
  v = IJ
  inv.sig = solve(Psi/v)

  Az = I*Kt / 2
  Bz = sum(diag((t(bz) - Wi %*% t(bw)) %*% P.mat %*% t(t(bz) - Wi %*% t(bw))))
  
  Aw = Kt / 2
  Bw = sapply(1:p, function(u) max(1, sum(diag( t(bw[,u]) %*% P.mat %*% (bw[,u])))))


  ## matrices to store within-iteration estimates 
  BW = array(NA, c(Kt, p, N.iter))
      BW[,,1] = bw
  BZ = array(NA, c(Kt, I, N.iter))
      BZ[,,1] = bz
  INV.SIG = array(NA, c(D, D, N.iter))
      INV.SIG[,,1] = inv.sig 
  LAMBDA.BW = matrix(NA, nrow = N.iter, ncol = p)
      LAMBDA.BW[1,] = lambda.bw = Aw/Bw
  LAMBDA.BZ = rep(NA, N.iter)
      LAMBDA.BZ[1] = lambda.ranef = Az/Bz
  
  y.post = array(NA, dim = c(IJ, D, (N.iter - N.burn)))

  cat("Beginning Sampler \n")

  for(i in 1:N.iter){

    ###############################################################
    ## update b-spline parameters for subject random effects
    ###############################################################

    for(subj in 1:length(unique(SUBJ))){

      t.designmat.Z = t(kronecker(rep(1, Ji[subj]), Theta))

      sigma = solve(t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% t(t.designmat.Z) + 
                    (lambda.ranef) * P.mat)
      mu = sigma %*% (t.designmat.Z %*% kronecker(diag(1, Ji[subj], Ji[subj]), inv.sig) %*% (as.vector(t(Y[which(SUBJ == unique(SUBJ)[subj]),]))) +
                       ((lambda.ranef) * P.mat) %*% bw %*% (Wi[subj,]))
      
      bz[,subj] = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = Kt, ncol = 1)
    }
    ranef.cur = Z.des %*% t(bz) %*% t(Theta)

    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################
    
    sigma = solve( kronecker(diag(lambda.bw), P.mat) + ((lambda.ranef) * tWIW) )
    mu = sigma %*% ((lambda.ranef) * tWI %*% as.vector(bz))
      
    bw = matrix(mvrnorm(1, mu = mu, Sigma = sigma), nrow = Kt, ncol = p)

    beta.cur = t(bw) %*% t(Theta)

    ###############################################################
    ## update inverse covariance matrix
    ###############################################################

    resid.cur = Y - ranef.cur
    inv.sig = solve(riwish(v + IJ, Psi + t(resid.cur) %*% resid.cur))

    ###############################################################
    ## update variance components
    ###############################################################

    ## lambda for beta's
    for(term in 1:p){
      a.post = Aw + Kt/2
      b.post = Bw[term] + 1/2 * bw[,term] %*% P.mat %*% bw[,term]
      lambda.bw[term] = rgamma(1, a.post, b.post)
    }
      
    ## lambda for random effects
    a.post = Az + I*Kt/2
    b.post = Bz + .5 * sum(sapply(1:I, function(u) (t(bz[,u]) - Wi[u,] %*% t(bw)) %*% P.mat %*% t(t(bz[,u]) - Wi[u,] %*% t(bw)) ))
    lambda.ranef = rgamma(1, a.post, b.post)

    ###############################################################
    ## save this iteration's parameters
    ###############################################################

    BW[,,i] = as.matrix(bw)
    BZ[,,i] = as.matrix(bz)
    
    INV.SIG[,,i] = inv.sig
    LAMBDA.BW[i,] = lambda.bw
    LAMBDA.BZ[i] = lambda.ranef

    if(i > N.burn){
      y.post[,,i - N.burn] = ranef.cur
    }

  }

  ###############################################################
  ## compute posteriors for this dataset
  ###############################################################

  ## main effects
  beta.post = array(NA, dim = c(p, D, (N.iter - N.burn)))
  for(n in 1:(N.iter - N.burn)){
    beta.post[,,n] = t(BW[,, n + N.burn]) %*% t(Theta)
  }
  beta.pm = apply(beta.post, c(1,2), mean)
  beta.LB = apply(beta.post, c(1,2), quantile, c(.025))
  beta.UB = apply(beta.post, c(1,2), quantile, c(.975))

  ## random effects
  b.pm = matrix(NA, nrow = I, ncol = D)
  for(i in 1:I){
  	b.post = matrix(NA, nrow = (N.iter - N.burn), ncol = D)
  	for(n in 1:(N.iter - N.burn)){
  	  b.post[n,] = BZ[,i, n + N.burn] %*% t(Theta)
  	}
  	b.pm[i,] = apply(b.post, 2, mean)
  }

  ## covariance matrix
  sig.pm = solve(apply(INV.SIG, c(1,2), mean))
  
  ## export fitted values
  ranef.pm = Z.des %*% b.pm
  Yhat = apply(y.post, c(1,2), mean)
  
  ret = list(beta.pm, beta.LB, beta.UB, ranef.pm, sig.pm, Yhat)
  names(ret) = c("beta.pm", "beta.LB", "beta.UB", "ranef.pm", "sig.pm", "Yhat")

  ret
  
}


###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################