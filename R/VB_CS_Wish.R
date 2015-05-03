#' Cross-sectional FoSR using Variational Bayes and Wishart prior
#' 
#' Fitting function for function-on-scalar regression for cross-sectional data.
#' This function estimates model parameters using VB and estimates
#' the residual covariance surface using a Wishart prior.
#' 
#' @param formula a formula indicating the structure of the proposed model. 
#' @param Kt number of spline basis functions used to estimate coefficient functions
#' @param data an optional data frame, list or environment containing the 
#' variables in the model. If not found in data, the variables are taken from 
#' environment(formula), typically the environment from which the function is 
#' called.
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
#' @export
#' 
vb_cs_wish = function(formula, data=NULL, Kt=5, alpha = .1){
  
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
  
  
  I = dim(Y)[1]
  
  ## fixed effect design matrix
  W.des = X
  
  ## bspline basis and penalty matrix
  D = dim(Y)[2]
  Theta = bs(1:D, df = Kt, intercept=TRUE, degree=3)
  
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(Theta) %*% t(diff0) %*% diff0 %*% Theta
  P2 = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  P.mat = alpha * P0 + (1-alpha) * P2
  
  ## hyper parameters for inverse gamma
  A = .01
  B = .01
  
  ## number of fixed effects
  p = dim(W.des)[2]
  
  ## data organization; these computations only need to be done once
  Y.vec = as.vector(t(Y))
  t.designmat.X = t(kronecker(W.des, Theta))
  sig.X = kronecker(t(W.des) %*% W.des, t(Theta)%*% Theta)
  
  ## matrices to to approximate paramater values
  sigma.q.BW = diag(1, p * Kt)
  mu.q.BW = matrix(0, nrow = Kt, ncol = p)  
  
  b.q.lambda.BW = rep(1, p)
  
  Psi = diag(1*I, D, D)
  v = I
  inv.sig = solve(Psi/v)
  
  lpxq=c(0,1)
  j=2
  
  cat("Beginning Algorithm \n")
  
  #  while(j<4 | (lpxq[j]-lpxq[j-1])>1.0E-1){
  while(j<11){
    
    ###############################################################
    ## update b-spline parameters for fixed effects
    ###############################################################
    
    sigma.q.BW = solve(t.designmat.X %*% kronecker(diag(1, I, I), inv.sig) %*% t(t.designmat.X) +
                         kronecker(diag((A+Kt/2)/b.q.lambda.BW), P.mat ))
    mu.q.BW = matrix( sigma.q.BW %*% (t.designmat.X %*% kronecker(diag(1, I, I), inv.sig) %*%  Y.vec), nrow = Kt, ncol = p)
    
    beta.cur = t(mu.q.BW) %*% t(Theta)
    
    ###############################################################
    ## update inverse covariance matrix
    ###############################################################
    
    T.BW.W = Theta %*% mu.q.BW %*% t(W.des)
    
    mu.q.v = v + I
    mu.q.Psi = Psi + t(Y) %*% Y - 
      t(Y) %*% t(T.BW.W) -
      T.BW.W %*% Y +
      T.BW.W %*% t(T.BW.W) + 
      matrix(apply(sapply(1:p, function(u) I * Theta %*% sigma.q.BW[(Kt * (u-1) +1):(Kt * u),(Kt * (u-1) +1):(Kt * u)] %*% t(Theta)), 1, sum), D, D)
    
    inv.sig = solve(mu.q.Psi/mu.q.v)
    
    ###############################################################
    ## update variance components
    ###############################################################
    
    ## lambda for fixed effects
    for(term in 1:dim(W.des)[2]){
      b.q.lambda.BW[term] = B + .5 * (t(mu.q.BW[,term]) %*% P.mat %*% mu.q.BW[,term] + 
                                        sum(diag(P.mat %*% sigma.q.BW[(Kt*(term-1)+1):(Kt*term),(Kt*(term-1)+1):(Kt*term)])))
    }
    
    ###############################################################
    ## lower bound
    ###############################################################
    
    curlpxq = 10
    lpxq = c(lpxq, curlpxq)
    j=j+1
    
    cat(".")
    
  }
  
  lpxq=lpxq[-(1:2)]
  
  ## compute CI for fixed effects
  beta.sd = beta.LB = beta.UB = matrix(NA, nrow = p, ncol = D)
  for(i in 1:p){
    beta.sd[i,] = sqrt(diag((Theta) %*% sigma.q.BW[(Kt*(i-1)+1):(Kt*i),(Kt*(i-1)+1):(Kt*i)] %*% t(Theta)))
    beta.LB[i,] = beta.cur[i,]-1.96*beta.sd[i,]
    beta.UB[i,] = beta.cur[i,]+1.96*beta.sd[i,]
  }
  
  ## convert objects from spam to matrix
  beta.cur = as.matrix(beta.cur)
  
  ## export fitted values
  Yhat.fixed = W.des %*% beta.cur
  Yhat = Yhat.fixed
  
  ## export various r2 values
  #r2.f = 1 - (sum((Y - Yhat.fixed)^2)/(IJ*D)) / (sum((Y)^2)/(IJ*D))
  #r2.fr = 1 - (sum((Y - Yhat)^2)/(IJ*D)) / (sum((Y)^2)/(IJ*D))
  
  ret = list(beta.cur, beta.UB, beta.LB, Yhat)
  names(ret) = c("beta.pm", "beta.UB", "beta.LB", "Yhat")
  
  ret
  
}


###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################