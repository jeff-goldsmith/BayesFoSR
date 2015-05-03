#' Cross-sectional FoSR using Variational Bayes and FPCA
#' 
#' Fitting function for function-on-scalar regression for cross-sectional data.
#' This function estimates model parameters using a VB and estimates
#' the residual covariance surface using FPCA.
#' 
#' @param formula a formula indicating the structure of the proposed model. 
#' @param Kt number of spline basis functions used to estimate coefficient functions
#' @param Kp number of FPCA basis functions to be estimated
#' @param data an optional data frame, list or environment containing the 
#' variables in the model. If not found in data, the variables are taken from 
#' environment(formula), typically the environment from which the function is 
#' called.
#' 
#' @references
#' Goldsmith, J., Kitago, T. (Under Review).
#' Assessing Systematic Effects of Stroke on Motor Control using Hierarchical 
#' Function-on-Scalar Regression.
#' 
#' @author Jeff Goldsmith \email{ajg2202@@cumc.columbia.edu}
#' @importFrom splines bs
#' @export
#' 
vb_cs_fpca = function(formula, data=NULL, Kt=5, Kp=2){
  
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
  
  ## subject covariates
  I = dim(X)[1]
  D = dim(Y)[2]
  p = dim(X)[2]
  
  ## bspline basis and penalty matrix
  Theta = bs(1:D, df=Kt, intercept=TRUE, degree=3)
  
  ## hyper parameters for inverse gaussians, bernoulli
  v0 = .001
  v1 = 100
  A = .5
  B = .5
  theta = .1
  
  ## matrices to to approximate paramater values
  sigma.q.BW = vector("list", p)
  for(k in 1:p){
    sigma.q.BW[[k]] = diag(1, Kt)
  }
  mu.q.BW = matrix(0, nrow = Kt, ncol = p)  
  
  mu.q.gamma = rep(1, p)  
  mu.q.dinv = kronecker(diag((1-mu.q.gamma)/v0 + mu.q.gamma/v1), diag(1, Kt, Kt))
  
  sigma.q.Bpsi = vector("list", Kp)
  for(k in 1:Kp){
    sigma.q.Bpsi[[k]] = diag(1, Kt)
  }
  mu.q.Bpsi = matrix(0, nrow = Kt, ncol = Kp)  
  
  sigma.q.C = vector("list", I)
  for(k in 1:I){
    sigma.q.C[[k]] = diag(1, Kp)
  }
  mu.q.C = matrix(rnorm(I*Kp, 0, .01), I, Kp)
  
  b.q.lambda.Bpsi = rep(1, Kp)
  b.q.sigma.me = 1
  
  ## data organization; these computations only need to be done once
  Y.vec = as.vector(t(Y))
  obspts.vec = !is.na(Y.vec)
  Y.vec = Y.vec[obspts.vec]
  J = sum(obspts.vec)
  
  t.designmat.X = t(kronecker(X, Theta)[obspts.vec,])
  Xstar = vector("list", length = I)
  XtX = matrix(0, Kt*p, Kt*p)
  sumXtX = matrix(0, Kt*p, Kt*p)
  for(i in 1:I){
    obs.points = which(!is.na(Y[i, ]))
    X.cur = Xstar[[i]] = kronecker(matrix(X[i,], nrow = 1, ncol = p), Theta)[obs.points,]
    XtX = XtX + crossprod(X.cur)
    sumXtX = sumXtX + t(X.cur)%*% X.cur
  }  

  ## initialize estimates of fixed, random and pca effects
  fixef.cur = matrix(0, nrow = I, ncol = D)
  pcaef.cur = matrix(0, I, D)

  for(i in 1:20){
    
    ###############################################################
    ## update regression coefficients
    ###############################################################
  
    mean.cur = as.vector(t(pcaef.cur))[obspts.vec]
    
    sigma.q.beta = solve(as.numeric((A + I*D/2)/(b.q.sigma.me)) * XtX + mu.q.dinv)
    mu.q.beta = matrix(sigma.q.beta %*% (as.numeric((A + I*D/2)/(b.q.sigma.me)) * t.designmat.X %*% (Y.vec - mean.cur)), nrow = Kt, ncol = p)
  
    beta.cur = t(mu.q.beta) %*% t(Theta)
    fixef.cur = as.matrix(X %*% beta.cur)
    
    ###############################################################
    ## update gammas -- not in this version!
    ###############################################################
  

    ###############################################################
    ## update b-spline parameters for PC basis functions
    ###############################################################
    
    mean.cur = as.vector(t(fixef.cur))[obspts.vec]
    designmat = kronecker(mu.q.C, Theta)[obspts.vec,]
    
    sigma.q.Bpsi = solve( 
      kronecker(diag(1, Kt, Kt), diag((A+Kt/2)/b.q.lambda.Bpsi)) + 
        as.numeric((A + J/2)/(b.q.sigma.me)) * f_sum(mu.q.c = mu.q.C, sig.q.c = sigma.q.C, theta = t(Theta), obspts.mat = !is.na(Y))
    )
    mu.q.Bpsi = matrix(((A + J/2)/(b.q.sigma.me)) * sigma.q.Bpsi %*% f_sum2(y = Y, fixef = fixef.cur, mu.q.c = mu.q.C, kt = Kt, theta = t(Theta)), nrow = Kt, ncol = Kp)
        
    psi.cur = t(mu.q.Bpsi) %*% t(Theta)
    ppT = (psi.cur) %*% t(psi.cur)
    
    ###############################################################
    ## scores for each individual
    ###############################################################
    
    for(subj in 1:I){
      obs.points = which(!is.na(Y[subj, ]))
      Theta_i = t(Theta)[,obs.points]
      sigma.q.C[[subj]] = solve( 
        diag(1, Kp, Kp ) +
          ((A + J/2)/(b.q.sigma.me)) * (f_trace(Theta_i = Theta_i, Sig_q_Bpsi = sigma.q.Bpsi, Kp = Kp, Kt = Kt) + 
                                          t(mu.q.Bpsi) %*% Theta_i %*% t(Theta_i) %*% mu.q.Bpsi)
      )
      
      mu.q.C[subj,] = ((A + J/2)/(b.q.sigma.me)) * sigma.q.C[[subj]] %*% as.matrix(psi.cur[,obs.points]) %*%  (Y[subj,obs.points] - fixef.cur[subj,obs.points] )
    }
    
    pcaef.cur =  as.matrix(mu.q.C %*% psi.cur)
    
    ###############################################################
    ## update variance components
    ###############################################################
  
    ## measurement error variance
    resid = as.vector(Y - fixef.cur - pcaef.cur)
    b.q.sigma.me = as.numeric(B + .5 * (crossprod(resid[!is.na(resid)]) + 
                                        sum(diag(sumXtX %*% sigma.q.beta)) + 
                                        f_sum4(mu.q.c= mu.q.C, sig.q.c = sigma.q.C, mu.q.bpsi = mu.q.Bpsi, sig.q.bphi = sigma.q.Bpsi, theta= Theta, obspts.mat = !is.na(Y))) )
        
    ## lambda for FPCA basis functions
    for(K in 1:Kp){
      b.q.lambda.Bpsi[K] = B + .5 * (t(mu.q.Bpsi[,K]) %*% diag(1, Kt, Kt) %*% mu.q.Bpsi[,K] + 
                                       sum(diag(diag(1, Kt, Kt) %*% sigma.q.Bpsi[(Kt*(K-1)+1):(Kt*K),(Kt*(K-1)+1):(Kt*K)])))
    }
  }

  ## export fitted values
  Yhat = X %*% beta.cur

  ## export variance components
  sigeps.pm = 1 / as.numeric((A + J/2)/(b.q.sigma.me))

  ## do svd to get rotated fpca basis
  temp = svd(t(psi.cur))
  psi.cur = t(temp$u)
  lambda.pm = temp$d
  
  ret = list(beta.cur, psi.cur, Yhat, sigeps.pm, lambda.pm)
  names(ret) = c("beta.pm", "psi.pm", "Yhat", "sigeps.pm", "lambda.pm")

  ret

}

###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################