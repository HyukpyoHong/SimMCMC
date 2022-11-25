library(ramcmc) 
library(mvtnorm)
library(MASS)


TimeDelayGillespieforXY <- function(A.X, B.X, alpha.X, beta.X, A.Y, B.Y, alpha.Y, beta.Y, K.M, repnum = 300000, maxT = 100, Volume = 1, Hc = 1){
  X <- 0
  XList <- rep(NA, repnum)
  Y <- 0
  YList <- rep(NA, repnum)
  currentTime <- 0
  TList <- rep(NA, repnum)
  Xbirth <- rep(0, maxT)
  Xdeath <- rep(0, maxT)
  Ybirth <- rep(0, maxT)
  Ydeath <- rep(0, maxT)
  n <- Hc
  k <- 1
  stackTimeX <- c()
  stackTimeY <- c()
  
  K.M <- K.M * Volume
  
  for (i in 1:repnum){
    a1 <- Volume * A.X
    a2 <- B.X * X
    # a3 <- lambda2 * X # Linear 
    a3 <- Volume * A.Y * (X^n / (K.M^n + X^n)) # Michaelis-Menten or Hill-Type
    a4 <- B.Y * Y
    a0 <- sum(a1,a2,a3,a4)
    # r2 <- runif(1)
    currentTime <- currentTime + rexp(1, rate = a0)
    
    stackTimeX <- sort(stackTimeX)  
    stackTimeY <- sort(stackTimeY)
    if(!(is.null(stackTimeX) & is.null(stackTimeY))){
      minStack <- min(stackTimeX, stackTimeY)
    } else {
      minStack <- Inf
    }
    if (currentTime < minStack){                                    
      r1 <- runif(1)
      if (r1 < a1/a0){
        XList[i] <-X
        YList[i] <- Y
        TList[i] <- currentTime
        stackTimeX <- c(stackTimeX, currentTime + k*rgamma(n=1, shape = alpha.X, rate = beta.X))
        # stackTimeX <- c(stackTimeX, currentTime) # without the delay of births of X
      } else if (r1 < (a1+a2)/a0){
        X <- X-1;
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- currentTime
        Xdeath[ceiling(currentTime)] = Xdeath[ceiling(currentTime)] + 1
      } else if (r1 < (a1+a2+a3)/a0){
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- currentTime
        stackTimeY <- c(stackTimeY, currentTime + k*rgamma(n=1, shape = alpha.Y, rate = beta.Y))
        # stackTimeY <- c(stackTimeY, currentTime) # without the delay of births of Y
      } else {
        Y <- Y-1;
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- currentTime
        Ydeath[ceiling(currentTime)] = Ydeath[ceiling(currentTime)] + 1
      }
    } else{
      if (min(stackTimeX) < min(stackTimeY)){
        X <- X+1;
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- minStack
        currentTime <- minStack
        stackTimeX <- stackTimeX[-1]
        Xbirth[ceiling(currentTime)] = Xbirth[ceiling(currentTime)] + 1
      } else {
        Y <- Y+1;
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- minStack
        currentTime <- minStack
        stackTimeY <- stackTimeY[-1]
        Ybirth[ceiling(currentTime)] = Ybirth[ceiling(currentTime)] + 1
      }
      if (currentTime > maxT){
        break
      }
    }
  }
  XList <- XList/Volume
  YList <- YList/Volume
  Xbirth <- Xbirth/Volume
  Xdeath <- Xdeath/Volume
  Ybirth <- Ybirth/Volume
  Ydeath <- Ydeath/Volume
    
  my_list <- list("XList" = XList, "YList" = YList, "TList" = TList, "Xbirth" = Xbirth, "Xdeath" = Xdeath, "Ybirth" = Ybirth, "Ydeath" = Ydeath)
  return(my_list)
}

tauleapingforXY <- function(A.X, B.X, alpha.X, beta.X, A.Y, B.Y, alpha.Y, beta.Y, K.M, maxT = 100, Volume = 1, Hc = 1){
  x.birth.propensity = A.X * Volume * KI(P = c(alpha.X, beta.X), maxt = maxT)
  x.birth.pois = rpois(maxT,lambda = x.birth.propensity)
  x.death.pois = rep(NA, maxT)
  x.trj.pois = rep(NA, maxT+1)
  x.trj.pois[1] = 0
  for(tt in 1:maxT){
    x.death.pois[tt] = rpois(1, lambda = 2/(2+B.X) * (B.X * x.trj.pois[tt] + B.X/2 * x.birth.pois[tt]))
    x.trj.pois[tt+1] = x.trj.pois[tt] + (x.birth.pois[tt] - x.death.pois[tt])
    if(x.trj.pois[tt+1] < 0){
      x.trj.pois[tt+1] = 0
      x.death.pois[tt] = x.trj.pois[tt] + x.birth.pois[tt]
    }
  }
  
  y.birth.propensity = A.Y * Volume * KI.Y(P = c(alpha.Y, beta.Y), in.X = x.trj.pois, K.M = K.M, Hc = Hc)
  
  y.birth.pois = rpois(maxT,lambda = y.birth.propensity)
  y.death.pois = rep(NA, maxT)
  y.trj.pois = rep(NA, maxT+1)
  y.trj.pois[1] = 0
  for(tt in 1:maxT){
    y.death.pois[tt] = rpois(1, lambda = 2/(2+B.Y) * (B.Y * y.trj.pois[tt] + B.Y/2 * y.birth.pois[tt]))
    y.trj.pois[tt+1] = y.trj.pois[tt] + (y.birth.pois[tt] - y.death.pois[tt])
    if(y.trj.pois[tt+1] < 0){
      y.trj.pois[tt+1] = 0
      y.death.pois[tt] = y.trj.pois[tt] + y.birth.pois[tt]
    }
  }
  
  my_list <- list("Xbirth" = x.birth.pois, "Xdeath" = x.death.pois, "Xsample" = x.trj.pois, 
                  "Ybirth" = y.birth.pois, "Ydeath" = y.death.pois, "Ysample" = y.trj.pois)
  return(my_list)
}


# This function calculates kappa Eq.(S9) in Sup. and gamma_k(m, Delta) 5th Eq. in page 2.
# Returning the cumulative sum of kappa: sum_{m=0}^{i}kappa(delta,m) in Eq. (S8) in Sup. 
# P: shape parameter alpha and rate parameter of beta in gamma delay distribution. 
KI <- function(P, maxt){
  a <- P[1]         #shape parameter alpha of gamma distribution in Eq.(S9) in Sup.
  b <- P[2]         #rate parameter beta of gamma distribution in Eq.(S9) in Sup.
  f <-function(x) pgamma(x,a,rate=b)
  k.j = rep(1,maxt)
  for (i in 1:maxt){
    k.j[i] = integrate(f,i-1,i)$value #equation in line 11 on page 3 in Sup.
  }
  return(k.j)
}   

# this function is an approximate propensity for birth process Y in eq (S9) & (S10) on page XY_DDE docu. 
KI.Y <-function(P,in.X, K.M, Hc = 1){
  a <- P[1]             #shape parameter alpha of gamma distribution in Eq.(S9) in Sup.
  b <- P[2]             #rate parameter beta of gamma distribution in Eq.(S9) in Sup.
  X <- in.X
  maxt <- length(X) - 1;
  A.m = function(t){
    # (1-t)*(pgamma(t,a,b)-pgamma(t-1,a,b)) + (gamma(a+1)/(b*gamma(a)))*(pgamma(t,a+1,b)-pgamma(t-1,a+1,b))
    (1-t)*(pgamma(t,a,b)-pgamma(t-1,a,b)) + a/b * (pgamma(t,a+1,b) - pgamma(t-1,a+1,b))
  }
  B.m = function(t){
    # t*(pgamma(t,a,b)-pgamma(t-1,a,b)) - (gamma(a+1)/(b*gamma(a)))*(pgamma(t,a+1,b)-pgamma(t-1,a+1,b)) 
    t*(pgamma(t,a,b)-pgamma(t-1,a,b)) - a/b * (pgamma(t,a+1,b) - pgamma(t-1,a+1,b)) 
  } 
  
  k.j =rep(1,maxt); A = rep(0,maxt); B=rep(0,maxt)
  for (m in 0:(maxt-1)){
    A[m+1] = max(integrate(A.m,m,m+1)$value, 0)
    B[m+1] = max(integrate(B.m,m,m+1)$value, 0)
  }
  for (i in 0:(maxt-1)){
    msum = 0
    for(m in 0:i){
      msum = msum + X[i-m+1]^Hc/(K.M^Hc + X[i-m+1]^Hc)*A[m+1] +X[i-m+2]^Hc/(K.M^Hc+X[i-m+2]^Hc)*B[m+1]
      #msum = msum + A[m+1] +B[m+1]
    }
    k.j[i+1]  = msum
  }
  k.j = cbind(k.j, A, B)
  return(k.j)
}

MH.P.X.all <- function(P,S,rep, r.X.birth, Ax, tun, pri.alpha.X, pri.beta.X, maxt, flatpri = FALSE, lowbnds = c(0,0)) {
  count = 0
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>lowbnds[1] && P.star[2]>lowbnds[2]){
      break
    }
  }
  KI.star <- KI(P.star, maxt = maxt)      # calculating kappa using current alpha & candidate of beta
  KI.m <- KI(P, maxt = maxt)              # calculating kappa using current alpha & beta
  l.lik.st <- 0
  l.lik <- 0
  
  for (i in 1:ncol(r.X.birth)){
    Rii <- r.X.birth[, i]
    l.lik.st <- l.lik.st + sum(Rii * log(KI.star + 1e-300), na.rm = T) - Ax * (sum(KI.star))
    l.lik <- l.lik + sum(Rii * log(KI.m + 1e-300), na.rm = T) - Ax * (sum(KI.m))
  }
  l.prior1.st <- dgamma(P.star[1], pri.alpha.X[1], pri.alpha.X[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , pri.alpha.X[1], pri.alpha.X[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], pri.beta.X[1], pri.beta.X[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , pri.beta.X[1], pri.beta.X[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l.lik.st - l.lik 
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l.lik.st - l.lik + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] - lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star;count = 1;
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-2/2)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count))
}



MH.P.Y.all <- function(P,S,rep, r.Y.birth, in.X.all, Ay, K.M, tun, pri.alpha.Y, pri.beta.Y, maxt, Hc = 1, flatpri = FALSE, lowbnds = c(0,0), noise_add = F){
  count = 0
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>lowbnds[1] && P.star[2]>lowbnds[2]){
      break
    }
  }
  noise_prob = c(0.001, 0.007, 0.027, 0.09, 0.21, 0.33, 0.21, 0.09, 0.027, 0.007, 0.001)
  noise_val = (-5):5

  l.lik.st <- 0
  l.lik <- 0

  for (i in 1:ncol(r.Y.birth)){
    Rii <- r.Y.birth[, i]

    lambda = Ay * KI.Y(P = P, in.X = in.X.all[,i], K.M = K.M, Hc = Hc)
    lambda.st = Ay * KI.Y(P = P.star, in.X = in.X.all[,i], K.M = K.M, Hc = Hc)

    if(noise_add){
      # log_lik_val = 0
      L = length(Rii)
      for(ii in 1:L){
        p_val = 0
        p_val.st = 0
        for(jj in 1:length(noise_val)){
          p_val= p_val + dpois(Rii[ii] - noise_val[jj], lambda[ii, 1]) * noise_prob[jj]
          p_val.st = p_val.st + dpois(Rii[ii] - noise_val[jj], lambda.st[ii, 1]) * noise_prob[jj]
        }
        l.lik = l.lik + sum(log(p_val + 1e-300), na.rm = T)
        l.lik.st = l.lik.st + sum(log(p_val.st + 1e-300), na.rm = T)
      }
    }else{
      l.lik.st = l.lik.st + sum(log(dpois(Rii,lambda.st[,1])+1e-300))
      l.lik    = l.lik    + sum(log(dpois(Rii,lambda[,1]   )+1e-300))
    }


  }
  l.prior1.st <- dgamma(P.star[1], pri.alpha.Y[1], pri.alpha.Y[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , pri.alpha.Y[1], pri.alpha.Y[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], pri.beta.Y[1], pri.beta.Y[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , pri.beta.Y[1], pri.beta.Y[2], log = TRUE)

  if(flatpri){
    logMH <- (l.lik.st - l.lik
              + log(pmvnorm(upper = c(P[1] - lowbnds[1] ,P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l.lik.st - l.lik + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] -lowbnds[1],P[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] - lowbnds[1],P.star[2] - lowbnds[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star;count = 1;
  }
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-2/2)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count))
}


##################################################################################################
#Functions for MH with independence chain for X trajectory
##################################################################################################

# m is the mean vector of the proposal distributuion of X trajectory. m is log transformed. 
# bi is the sd vector of the proposal distributuion of X trajectory. bi can be tuned for optimal acceptance ratio.  

##################################################################################################
#The rest of the functions were not used in this study where delay only occurs in X.
##################################################################################################


MH.KM.all <- function(km, s, rep, r.all, x.all, b, pri.KM, Delta.Y, ay, Hc =1, flatpri = FALSE, noise_add = F){
  # km: Michaelis-Menten constant in the current iteration.
  # s: scaling factor for tunning
  # rep: iteration number
  # x.all: X trajectories to construct the likelihood for the birth of Y. Dim: (maxT + 1) * (number of experiments)
  # r.all: numbers of the birth reactions of Y. Dim: maxT * (number of experiments)
  # b: tunning parameters for KM
  # pri.KM : prior distribution parameters for KM
  noise_prob = c(0.001, 0.007, 0.027, 0.09, 0.21, 0.33, 0.21, 0.09, 0.027, 0.007, 0.001)
  noise_val = (-5):5
  count = 0
  repeat{
    u = rnorm(1,0,b)
    km.star = km + s * u
    if(km.star > 0) break
  }

  l.lik.st = 0;
  l.lik = 0;

  for(j in 1:ncol(r.all)){
    x = x.all[,j]
    r = r.all[,j]

    lambda = ay * KI.Y(P = Delta.Y, in.X = x , K.M = km, Hc = Hc)
    lambda.st = ay * KI.Y(P = Delta.Y, in.X = x , K.M = km.star, Hc = Hc)

    if(noise_add){
      # log_lik_val = 0
      L = length(r)
      for(ii in 1:L){
        p_val = 0
        p_val.st = 0
        for(jj in 1:length(noise_val)){
          p_val= p_val + dpois(r[ii] - noise_val[jj], lambda[ii, 1]) * noise_prob[jj]
          p_val.st = p_val.st + dpois(r[ii] - noise_val[jj], lambda.st[ii, 1]) * noise_prob[jj]
        }
        l.lik = l.lik + sum(log(p_val + 1e-300), na.rm = T)
        l.lik.st = l.lik.st + sum(log(p_val.st + 1e-300), na.rm = T)
      }
    }else{
      l.lik.st = l.lik.st + sum(log(dpois(r,lambda.st[,1])+1e-300))
      l.lik    = l.lik    + sum(log(dpois(r,lambda[,1]   )+1e-300))
    }
  }

  if(flatpri){
    logMH = (l.lik.st - l.lik #+ dgamma(km.star, pri.KM[1], pri.KM[2],log=T) - dgamma(km, pri.KM[1], pri.KM[2],log=T)
             + pnorm(km, 0, s*b, log.p = T) - pnorm(km.star, 0, s*b, log.p = T))
  }else{
    logMH = (l.lik.st - l.lik + dgamma(km.star, pri.KM[1], pri.KM[2],log=T) - dgamma(km, pri.KM[1], pri.KM[2],log=T)
             + pnorm(km, 0, s*b, log.p = T) - pnorm(km.star, 0, s*b, log.p = T))
  }
  if(!is.nan(logMH) && runif(1)<exp(logMH)) {
    km=km.star; count = 1;
  }
  alpha = min(exp(logMH),1)
  s=ramcmc::adapt_S(s,u,alpha,rep,gamma = min(1,(1*rep)^(-2/3)))
  return(list(km=km,s=s, count=count))
}





