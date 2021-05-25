library(ramcmc) 
library(mvtnorm)
library(MASS)


impute_r.X <- function(x, B.X){
  X.diff = diff(x)
  r4 = rpois(length(X.diff),B.X*(x[-length(x)]+x[-1])/2)
  r3 = round(X.diff + r4) # input x could be a non-integer, so r3 should rounded to be an integer.
  if(min(r3) >= 0){ # success to find non-negative pair of the reaction numbers.
    return(cbind(r3,r4)) # return normally
  }else{
    return(matrix(-1, nrow=2,ncol=2)) # return abnormally: (at least) one of the reaction numbers is still negative.
  }
}

#sampling r_1i and r_2i
# Step 1 & 2
# generating only 1 set of number of reaction of Y; 
impute_r.Y <- function(y, B.Y){
  Y.diff = diff(y)
  r2 = rpois(length(Y.diff),B.Y*(y[-length(y)]+y[-1])/2)
  r1 = Y.diff + r2
  for(jj in 1:length(r1)){
    if (r1[jj] < 0){
      r1[jj] = 0
      r2[jj] = -Y.diff[jj]
    }
  }
  return(cbind(r1,r2))
}

#generating 1000 sets of number of reaction Y and calculating average of number of reaction of Y; 
impute_r.Y.mean <- function(y, B.Y){
  Y.diff = diff(y)
  r2.all = rep(0,length(Y.diff)); r1.all =rep(0,length(Y.diff));
  for(i in 1:1000){
    r2 = rpois(length(Y.diff),B.Y*(y[-length(y)]+y[-1])/2)
    r1 = Y.diff + r2
    r2.all = r2.all + r2
    r1.all = r1.all + r1
  }
  r1 = round(r1.all/1000); r2 = round(r2.all/1000)
  return(cbind(r1,r2))
}

dpois_conti <- function(x, lambda){

   return(lambda^x/gamma(x+1) * exp(-lambda))

  }


#propensity of X 
# birth & death process. birth process has delay. 
# Same to bioinformatics 
H.X <-function(x, ki, th){
  G = c(ki,x)
  return(th * G)
}

#propensity of Y: enzyme kinetics 
H.Y <-function(x,y,th){
  G = c(th[1]*x/(th[3]+x  ), th[2]*y )
  return(G)
}

#propensity of Y: linear in X 
H.Y.linear <-function(x,y,th){
  G = c(x, y)
  return(th * G)
}


# function for sampling # of reaction 
# Boys et al.'s method
# draw sample from Skellman distribution
cand <-function(r,tun){
  repeat{
    z = 1+r^2/tun;# z=0+1/b;
    u = rpois(1,z)-rpois(1,z)
    r.star = round(r) + u
    if(r.star>=0) break
  }
  return(r.star)
}


# function of f(r_3i|tilde{r}_3i) in step 3.3
# Proposal distribution using Skellman distribution 
prop.R <-function(r.star,r,tun){
  lambda = 1+r^2/tun
  prop.st = log(besselI(2*lambda   ,abs(r.star-r),expon.scaled = T)+1e-300)
  return(prop.st)
}

# function of p() in step 3.3 
# log likelihood function of Poisson distribution. 
q.R <- function(r,H0){
  q.st = sum(dpois(r,H0, log = T))
  return(q.st)
}

# 2020.03.25 = Hyukpyo. Generate delayed X for the proposal mean to implement the independent chain MCMC.
# 2020.09.01. (Hyukpyo) This code 'TimeDelayGillespieforX' is not used anymore.
TimeDelayGillespieforX <- function(A.X, B.X, alpha.X, beta.X, repnum = 3000, maxT = 100){
  X <- 0
  XList <- rep(NA, repnum)
  currentTime <- 0
  TList <- rep(NA, repnum)
  n <- 1
  k <- 1
  stackTimeX <- c()
  
  for (i in 1:repnum){
    a1 <- A.X
    a2 <- B.X * X
    a0 <- sum(a1,a2)
    # r2 <- runif(1)
    currentTime <- currentTime + rexp(1, rate = a0)
    
    stackTimeX <- sort(stackTimeX)  
    if(!(is.null(stackTimeX))){
      minStack <- min(stackTimeX)
    } else {
      minStack <- Inf
    }
    if (currentTime < minStack){                                    
      r1 <- runif(1)
      if (r1 < a1/a0){
        XList[i] <- X
        TList[i] <- currentTime
        stackTimeX <- c(stackTimeX, currentTime + k*rgamma(n=1, shape = alpha.X, rate = beta.X))
        #stackTimeX <- c(stackTimeX, currentTime)
      } else{
        X <- X-1;
        XList[i] <- X
        TList[i] <- currentTime
      } 
    } else{
      X <- X+1;
      XList[i] <- X
      TList[i] <- minStack
      currentTime <- minStack
      stackTimeX <- stackTimeX[-1]
    }
    if (currentTime > maxT){
      break
    }
  }
  my_list <- list("XList" = XList, "TList" = TList)
  return(my_list)
}

TimeDelayGillespieforXR <- function(A.X, B.X, alpha.X, beta.X, repnum = 3000, maxT = 100, Volume = 1){
  X <- 0
  XList <- rep(NA, repnum)
  Xbirth <- rep(0, maxT)
  Xdeath <- rep(0, maxT)
  currentTime <- 0
  TList <- rep(NA, repnum)
  n <- 1
  k <- 1
  stackTimeX <- c()
  
  for (i in 1:repnum){
    a1 <- Volume * A.X
    a2 <- B.X * X
    a0 <- sum(a1,a2)
    # r2 <- runif(1)
    currentTime <- currentTime + rexp(1, rate = a0)
    
    stackTimeX <- sort(stackTimeX)  
    if(!(is.null(stackTimeX))){
      minStack <- min(stackTimeX)
    } else {
      minStack <- Inf
    }
    if (currentTime < minStack){                                    
      r1 <- runif(1)
      if (r1 < a1/a0){
        XList[i] <- X
        TList[i] <- currentTime
        stackTimeX <- c(stackTimeX, currentTime + k*rgamma(n=1, shape = alpha.X, rate = beta.X))
        #stackTimeX <- c(stackTimeX, currentTime)
      } else{
        # death of X occurs
        X <- X-1;
        XList[i] <- X
        TList[i] <- currentTime
        Xdeath[ceiling(currentTime)] = Xdeath[ceiling(currentTime)] + 1
      } 
    } else{
      X <- X+1;
      XList[i] <- X
      TList[i] <- minStack
      currentTime <- minStack
      stackTimeX <- stackTimeX[-1]
      Xbirth[ceiling(currentTime)] = Xbirth[ceiling(currentTime)] + 1
    }
    if (currentTime > maxT){
      break
    }
  }
  
  XList <- XList/Volume
  Xbirth <- Xbirth/Volume
  Xdeath <- Xdeath/Volume
  
  my_list <- list("XList" = XList, "TList" = TList, "Xbirth" = Xbirth, "Xdeath" = Xdeath)
  return(my_list)
}

TimeDelayGillespieforXY <- function(A.X, B.X, alpha.X, beta.X, A.Y, B.Y, alpha.Y, beta.Y, K.M, repnum = 300000, maxT = 100, Volume = 1){
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
  n <- 1
  k <- 1
  stackTimeX <- c()
  stackTimeY <- c()
  
  K.M <- K.M * Volume^n 
  
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


GillespieforXYonlyXdelay <- function(A.X, B.X, alpha.X, beta.X, A.Y, B.Y, K.M, repnum = 300000, maxT = 100, Volume = 1){
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
  n <- 1
  k <- 1
  stackTimeX <- c()

  K.M <- K.M * Volume^n 
  
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
    if(!(is.null(stackTimeX))){
      minStack <- min(stackTimeX)
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
        Y <- Y+1
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- currentTime
        Ybirth[ceiling(currentTime)] = Ybirth[ceiling(currentTime)] + 1
      } else {
        Y <- Y-1;
        XList[i] <- X
        YList[i] <- Y
        TList[i] <- currentTime
        Ydeath[ceiling(currentTime)] = Ydeath[ceiling(currentTime)] + 1
      }
    } else{
      X <- X+1;
      XList[i] <- X
      YList[i] <- Y
      TList[i] <- minStack
      currentTime <- minStack
      stackTimeX <- stackTimeX[-1]
      Xbirth[ceiling(currentTime)] = Xbirth[ceiling(currentTime)] + 1
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

GillespieforXYonlyYdelay <- function(A.X, B.X, A.Y, B.Y, alpha.Y, beta.Y, K.M, repnum = 300000, maxT = 100, Volume = 1){
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
  n <- 1
  k <- 1
  stackTimeY <- c()
  
  K.M <- K.M * Volume^n 
  
  for (i in 1:repnum){
    a1 <- Volume * A.X
    a2 <- B.X * X
    # a3 <- lambda2 * X # Linear 
    a3 <- Volume * A.Y * (X^n / (K.M^n + X^n)) # Michaelis-Menten or Hill-Type
    a4 <- B.Y * Y
    a0 <- sum(a1,a2,a3,a4)
    # r2 <- runif(1)
    currentTime <- currentTime + rexp(1, rate = a0)
    
    stackTimeY <- sort(stackTimeY)
    if(!(is.null(stackTimeY))){
      minStack <- min(stackTimeY)
    } else {
      minStack <- Inf
    }
    if (currentTime < minStack){                                    
      r1 <- runif(1)
      if (r1 < a1/a0){
        X <- X+1;
        XList[i] <-X
        YList[i] <- Y
        TList[i] <- currentTime
        Xbirth[ceiling(currentTime)] = Xbirth[ceiling(currentTime)] + 1
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
      Y <- Y+1;
      XList[i] <- X
      YList[i] <- Y
      TList[i] <- minStack
      currentTime <- minStack
      stackTimeY <- stackTimeY[-1]
      Ybirth[ceiling(currentTime)] = Ybirth[ceiling(currentTime)] + 1
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




# Metropolis-Hastings algorithm for updating trajectory of X and number of reaction of X 
MH.X.R <-function(r,x,ki,bi=tun.B){
  maxt <- length(x) - 1;
  count = matrix(0, nrow = maxt, ncol = 1)
  for (i in 1:(length(X)-1)) {
    r.X = c(r[i,3],r[i,4]); 
    r.Y = c(r[i,1],r[i,2])  
    xi=x[i];xi1=x[i+1]
    k <- ki[i]
    #bi=b[i,]
    r3.star = cand(r.X[1],bi[3])
    r4.star = cand(r.X[2],bi[4])
    r.X.star = c(r3.star,r4.star)
    xi1.star = r.X.star[1]-r.X.star[2] + xi
    if(xi1.star>=0){
      prop.r3.st = prop.R(r.X.star[1], r.X[1]     , bi[3])
      prop.r4.st = prop.R(r.X.star[2], r.X[2]     , bi[4])
      prop.r3    = prop.R(r.X[1]     , r.X.star[1], bi[3])
      prop.r4    = prop.R(r.X[2]     , r.X.star[2], bi[4])
      
      H0.X.st = (H.X(xi, k, theta.X)+H.X(xi1.star,k, theta.X))/2
      H0.X    = (H.X(xi, k, theta.X)+H.X(xi1     ,k, theta.X))/2
      H0.Y.st = (H.Y(xi, Y[i], theta.Y)+H.Y(xi1.star, Y[i+1], theta.Y))/2
      H0.Y    = (H.Y(xi, Y[i], theta.Y)+H.Y(xi1,      Y[i+1], theta.Y))/2
      if(min(H0.X.st,H0.X,H0.Y.st,H0.Y) >= 0){
        q.X.st = sum(log(dpois(r.X.star,H0.X.st)+1e-300), na.rm = T)
        q.X    = sum(log(dpois(r.X     ,H0.X   )+1e-300), na.rm = T)
        q.Y.st = sum(log(dpois(r.Y,H0.Y.st)+1e-300), na.rm = T)
        q.Y    = sum(log(dpois(r.Y,H0.Y   )+1e-300), na.rm = T)
        
        logMH = prop.r3 + prop.r4 - prop.r3.st - prop.r4.st + q.X.st + q.Y.st - q.X - q.Y
        
        if(!is.nan(logMH) && runif(1)<exp(logMH)) {
          r[i,3:4]=r.X.star; x[i+1] = xi1.star; count[i] = 1;
        }
      }
    }  
  }
  return(list(r=r, x=x, count=count))
}


MH.R <-function(r,x,ki,thetax, bi=tun.B){# if b=tun.B then remove $ from bi=b[i,]
  maxt <- length(x) - 1;
  count = matrix(0, nrow = maxt, ncol = 1)
  for (i in 1:maxt) {
    r.X = c(r[i,3],r[i,4]); 
    xi=x[i];xi1=x[i+1]; d = xi1 - xi;
    k <- ki[i]
    #bi=b[i,]
    r3.star = cand(r.X[1],bi[3])
    r.X.star = cbind(r3.star, r3.star - d)
    # cat(i, ", ",rep,", ",r.X.star, ", ",r.X, ", ","\n")    
    if(min(r.X.star)>=0){
      prop.r3.st = prop.R(r.X.star[1], r.X[1]     , bi[3])
      prop.r3    = prop.R(r.X[1]     , r.X.star[1], bi[3])
      
      H0.X    = (H.X(xi, k, thetax)+H.X(xi1 ,k, thetax))/2
      q.X.st = sum(log(dpois(r.X.star,H0.X)+1e-300), na.rm = T)
      q.X    = sum(log(dpois(r.X     ,H0.X)+1e-300), na.rm = T)
      logMH = prop.r3 - prop.r3.st + q.X.st  - q.X 
      # cat(i, ", ",rep,", ",r.X.star, ", ",r.X, ", ", prop.r3.st, ", ", prop.r3, ", ", x[i+1],", ",
      #     q.X.st,", ", q.X, ", ", exp(logMH),"\n ")
      if(!is.nan(logMH) && runif(1)<exp(logMH)) {
        r[i,3:4]=r.X.star; count[i] = 1;
      }
    }  
  }
  return(list(r=r, count=count))
}

# ouput 100th iteration value
MH.R.repeat <-function(r,x,ki,thetax, bi=tun.B){# if b=tun.B then remove $ from bi=b[i,]
  maxt <- length(x) - 1;
  count = matrix(0, nrow = maxt, ncol = 1)
  X.diff = diff(X) #x(i+1) - x(i)    
  r = matrix(0,ncol = 4, nrow = maxt) #saving number of reaction 
  for (i in 1:maxt) {
    r[i,3] = max(X.diff[i],0)  # # of birth reaction of Y
    r[i,4] = max(-X.diff[i],0) # # of death reaction of Y
  }
  for (rep in 1:rrepeat){
    for (i in 1:maxt) {
      r.X = c(r[i,3],r[i,4]); 
      xi=x[i];xi1=x[i+1]; d = X.diff[i];
      k <- ki[i]
      #bi=b[i,]
      r3.star = cand(r.X[1],bi[3])
      r.X.star = cbind(r3.star, r3.star - d)
      # cat(i, ", ",rep,", ",r.X.star, ", ",r.X, ", ","\n")    
      if(min(r.X.star)>=0){
        prop.r3.st = prop.R(r.X.star[1], r.X[1]     , bi[3])
        prop.r3    = prop.R(r.X[1]     , r.X.star[1], bi[3])
        
        H0.X    = (H.X(xi, k, thetax)+H.X(xi1 ,k, thetax))/2
        q.X.st = sum(log(dpois(r.X.star,H0.X)+1e-300), na.rm = T)
        q.X    = sum(log(dpois(r.X     ,H0.X)+1e-300), na.rm = T)
        logMH = prop.r3 - prop.r3.st + q.X.st  - q.X 
        if(!is.nan(logMH) && runif(1)<exp(logMH)) {
          r[i,3:4]=r.X.star; count[i] = count[i]+ 1;
        }
      }  
    }
  }
  return(list(r=r, count=count))
}

# ouput average of last 50 iteration's values
MH.R.repeat.ave <-function(r,x,ki,thetax, bi=tun.B){# if b=tun.B then remove $ from bi=b[i,]
  maxt <- length(x) - 1;
  count = matrix(0, nrow = maxt, ncol = 1)
  X.diff = diff(X) #x(i+1) - x(i)    
  r = matrix(0,ncol = 4, nrow = maxt) #saving number of reaction 
  temp = matrix(0,ncol = 2, nrow = maxt) #saving number of reaction 
  for (i in 1:maxt) {
    r[i,3] = max(X.diff[i],0)  # # of birth reaction of Y
    r[i,4] = max(-X.diff[i],0) # # of death reaction of Y
  }
  for (rep in 1:rrepeat){
    for (i in 1:maxt) {
      r.X = c(r[i,3],r[i,4]); 
      xi=x[i];xi1=x[i+1]; d = X.diff[i];
      k <- ki[i]
      #bi=b[i,]
      r3.star = cand(r.X[1],bi[3])
      r.X.star = cbind(r3.star, r3.star - d)
      if(min(r.X.star)>=0){
        prop.r3.st = prop.R(r.X.star[1], r.X[1]     , bi[3])
        prop.r3    = prop.R(r.X[1]     , r.X.star[1], bi[3])
        H0.X    = (H.X(xi, k, thetax)+H.X(xi1 ,k, thetax))/2
        q.X.st = sum(log(dpois(r.X.star,H0.X)+1e-300), na.rm = T)
        q.X    = sum(log(dpois(r.X     ,H0.X)+1e-300), na.rm = T)
        logMH = prop.r3 - prop.r3.st + q.X.st  - q.X 
        if(!is.nan(logMH) && runif(1)<exp(logMH)) {
          r[i,3:4]=r.X.star; count[i] = count[i] + 1;
        }
      }
      temp[i,] = temp[i,] + r[i,3:4]
    }
  }
  r[,3:4] = round(temp/(rrepeat))
  return(list(r=r, count=count))
}

MH.X <-function(r,x,kix,km, thetax, bi=tun.X, ptnum = 4, useall = FALSE){
  count=0
  x.star = x
  prop.X.st = rep(0,length(x))
  prop.X    = rep(0,length(x))
  H0.X.st = matrix(rep(0,2*(length(x)-1)),ncol=2,nrow = (length(x)-1))
  H0.X    = matrix(rep(0,2*(length(x)-1)),ncol=2,nrow = (length(x)-1))
  
  for (i in 2:length(x)){
    x.star[i] = cand(x[i],bi[i]) 
    prop.X.st[i] = prop.R(x.star[i], x[i]     , bi[i])
    prop.X[i]    = prop.R(x[i]     , x.star[i], bi[i])
    H0.X.st[i-1,] = (H.X(x.star[i-1], kix[i-1], thetax)+H.X(x.star[i], kix[i-1], thetax))/2
    H0.X[i-1,]    = (H.X(x[i-1]     , kix[i-1], thetax)+H.X(x[i]     , kix[i-1], thetax))/2
  }
  if (useall == TRUE){
    fy.st = A.Y * KI.Y(Delta.Y,x.star, K.M=km)
    fy    = A.Y * KI.Y(Delta.Y,x     , K.M=km)
  }else{
    fy.st = A.Y * KI.Ynt(Delta.Y,x.star, N = ptnum, K.M=km)
    fy    = A.Y * KI.Ynt(Delta.Y,x     , N = ptnum, K.M=km)
  }
  
  q.X.st = sum(log(dpois(r[,3:4], H0.X.st)+1e-300), na.rm = T)
  q.X    = sum(log(dpois(r[,3:4], H0.X   )+1e-300), na.rm = T)
  q.Y.st = sum(log(dpois(r[,1],fy.st[,1])+1e-300), na.rm = T)
  q.Y    = sum(log(dpois(r[,1],fy[,1]   )+1e-300), na.rm = T)
  prop.X.st = sum(prop.X.st, na.rm = T)
  prop.X    = sum(prop.X   , na.rm = T)
  
  logMH = prop.X - prop.X.st  + q.Y.st - q.Y  + q.X.st - q.X 
  if(!is.nan(logMH) && runif(1)<exp(logMH)) {
    x=x.star; count = 1;
  }    
  return(list(x=x, count=count))
}

MH.X.R_XY <-function(r,x,kix,kiy,km, thetax, bi=tun.B){# if b=tun.B then remove $ from bi=b[i,]
  maxt <- length(x) - 1;
  count = matrix(0, nrow = maxt, ncol = 1)
  for (i in 1:maxt) {
    r.X = c(r[i,3],r[i,4]); 
    r.Y = c(r[i,1],r[i,2])  
    xi=x[i];xi1=x[i+1]
    k.x <- kix[i]
    #bi=b[i,]
    r3.star = cand(r.X[1],bi[3])
    r4.star = cand(r.X[2],bi[4])
    r.X.star = c(r3.star,r4.star)
    xi1.star = r.X.star[1]-r.X.star[2] + xi
    if(xi1.star>=0){
      x.st=x; x.st[i+1]=xi1.star;
      fyi = 0; fyi.st=0;
      for(m in 0:(i-1)){
        fyi.st = fyi.st + x.st[i-m]/(km + x.st[i-m])*kiy[m+1,2] + x.st[i-m+1]/(km + x.st[i-m+1])*kiy[m+1,3]
        fyi    = fyi    +    x[i-m]/(km +    x[i-m])*kiy[m+1,2] +    x[i-m+1]/(km +    x[i-m+1])*kiy[m+1,3]
      }
      fyi.st = A.Y * fyi.st
      fyi    = A.Y * fyi 
      
      prop.r3.st = prop.R(r.X.star[1], r.X[1]     , bi[3])
      prop.r3    = prop.R(r.X[1]     , r.X.star[1], bi[3])
      prop.r4.st = prop.R(r.X.star[2], r.X[2]     , bi[4])
      prop.r4    = prop.R(r.X[2]     , r.X.star[2], bi[4])
      
      H0.X.st = (H.X(xi, k.x, thetax)+H.X(xi1.star,k.x, thetax))/2
      H0.X    = (H.X(xi, k.x, thetax)+H.X(xi1     ,k.x, thetax))/2
      #if(min(H0.X.st,H0.X,fyi.st,fyi) < 0) break
      q.X.st = sum(dpois(r.X.star,H0.X.st, log = T ))
      q.X    = sum(dpois(r.X     ,H0.X, log = T   ))
      q.Y.st = dpois(r.Y[1],fyi.st, log = T)
      q.Y    = dpois(r.Y[1],fyi,    log = T)
      logMH = prop.r3 + prop.r4 - prop.r3.st - prop.r4.st + q.X.st + q.Y.st - q.X - q.Y
      #cat(i, ", ",r.Y[1], ", ",fyi.st, ", ", fyi, ", ", xi1.star, ", ", x[i+1],", ", 
      #                 q.X.st,", ", q.X, ", ", q.Y.st, ", ", q.Y,", ", exp(logMH),"\n ")
      if(!is.nan(logMH) && runif(1)<exp(logMH)) {
        r[i,3:4]=r.X.star; x[i+1] = xi1.star; count[i] = 1;
      }    
    } 
  }
  return(list(r=r, x=x, count=count))
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
KI.Y <-function(P,in.X, K.M){
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
      msum = msum + X[i-m+1]/(K.M + X[i-m+1])*A[m+1] +X[i-m+2]/(K.M+X[i-m+2])*B[m+1]
      #msum = msum + A[m+1] +B[m+1]
    }
    k.j[i+1]  = msum
  }
  k.j = cbind(k.j, A, B)
  return(k.j)
}

# this function is an approximate propensity for birth process Y in eq (S9) & (S10) on page XY_DDE docu. 
# this function assume that Y(t) is affected by only X(t-5) and X(t-6).
# if t<6 then Y(t) is affected by X(1)
KI.Ynt <-function(P,in.X, N = 4, K.M){
  a <- P[1]             #shape parameter alpha of gamma distribution in Eq.(S9) in Sup.
  b <- P[2]             #rate parameter beta of gamma distribution in Eq.(S9) in Sup.
  X <- in.X
  maxt <- length(X) - 1;
  A.m = function(t){
    (1-t)*(pgamma(t,a,b)-pgamma(t-1,a,b)) + (gamma(a+1)/(b*gamma(a)))*(pgamma(t,a+1,b)-pgamma(t-1,a+1,b))
  }
  B.m = function(t){
    t*(pgamma(t,a,b)-pgamma(t-1,a,b)) - (gamma(a+1)/(b*gamma(a)))*(pgamma(t,a+1,b)-pgamma(t-1,a+1,b)) 
  } 
  
  k.j =rep(1,maxt); A = rep(0,maxt); B=rep(0,maxt)
  for (m in 0:(maxt-1)){
    A[m+1] = integrate(A.m,m,m+1)$value
    B[m+1] = integrate(B.m,m,m+1)$value
  }
  
  highIdx = sort(x=A, decreasing = TRUE, index.return=TRUE)$ix[1:N]
  #highIdxB = sort(x=B, decreasing = TRUE, index.return=TRUE)$ix[1:N]
  for(i in 0:(maxt-1)){
    msum = 0;
    if(i >= min(highIdx)){
      for(m in highIdx){
        if(i-(m-1)+1 > 0){
          msum = msum + 0.5 / sum(A[highIdx]) * X[i-(m-1)+1]/(K.M + X[i-(m-1)+1])*A[m] + 0.5 / sum(B[highIdx]) * X[i-(m-1)+2]/(K.M+X[i-(m-1)+2])*B[m]
        }
        else if(i-(m-1)+2 > 0){
          msum = msum + 0.5 / sum(B[highIdx]) * X[i-(m-1)+2]/(K.M+X[i-(m-1)+2])*B[m]
        }
        #msum = msum + A[m+1] +B[m+1]
      }
    }else{
      msum = msum + 0.5 / sum(A[highIdx]) * X[1]/(K.M + X[1])*A[1] + 0.5 / sum(B[highIdx]) * X[2]/(K.M+X[2])*B[1]
    }
    k.j[i+1]  = msum
  }
  k.j = cbind(k.j, A, B)
  return(k.j)
}

# RAM (robust adaptive Metropolis) method for rate parameter alpha & beta of gamma delay distribution
# p =(alpha, beta): parameter value at the current step in MCMC
# Ri: T by 1 vertor (r_{11}, r_{12}, ..., r_{i,T-1} )
# A ; A in Eq. (5) 
# This function returns upadated beta & alpha 
# refering Vihola(2012), Statistics and Computing, 22(5):997-1008. 
# using informative gamma prior for beta

MH.P.X <- function(P,S,rep, r.X.birth, Ax, tun, pri.alpha.X, pri.beta.X, maxt, flatpri = FALSE){
  count = 0
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>0 && P.star[2]>0){
      break
    }
  }
  
  KI.star <- KI(P.star, maxt = maxt)      # calculating kappa using current alpha & candidate of beta
  KI.m <- KI(P, maxt = maxt)              # calculating kappa using current alpha & beta
  
  l.lik.st <- sum(r.X.birth * log(KI.star), na.rm = T) - Ax * (sum(KI.star))
  l.lik <- sum(r.X.birth * log(KI.m), na.rm = T) - Ax * (sum(KI.m))
  
  l.prior1.st <- dgamma(P.star[1], pri.alpha.X[1], pri.alpha.X[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , pri.alpha.X[1], pri.alpha.X[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], pri.beta.X[1], pri.beta.X[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , pri.beta.X[1], pri.beta.X[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l.lik.st - l.lik 
              + log(pmvnorm(upper = c(P[1] ,P[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] ,P.star[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l.lik.st - l.lik + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] ,P[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] ,P.star[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  # print(S)
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star;count = 1;
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-2/2)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count))
}

MH.P.Y <- function(P,S,rep, r.Y.birth, in.X, Ay, K.M, tun, pri.alpha.Y, pri.beta.Y, maxt, flatpri = FALSE){
  count = 0
  repeat{
    u = mvrnorm(1,c(0,0),diag(c(tun[1],tun[2])))
    P.star = P + S%*%u
    # if(P.star[1]>1 && P.star[2]>0){
    #   break
    # }
    if(P.star[1]>0 && P.star[2]>0){
      break
    }
  }
  
  ## Under editing. - 2020.06.30
  
  KI.Y.star <- KI.Y(P = P.star, in.X = in.X,  K.M = K.M)      # calculating kappa using current alpha & candidate of beta
  KI.Y.m <- KI.Y(P = P, in.X = in.X,  K.M = K.M)              # calculating kappa using current alpha & beta
  
  l.lik.st <- sum(r.Y.birth * log(KI.Y.star[,1]), na.rm = T) - Ay * (sum(KI.Y.star[,1])) 
  l.lik <- sum(r.Y.birth * log(KI.Y.m[,1]), na.rm = T) - Ay * (sum(KI.Y.m[,1]))
  
  l.prior1.st <- dgamma(P.star[1], pri.alpha.Y[1], pri.alpha.Y[2], log = TRUE)
  l.prior1    <- dgamma(P[1]     , pri.alpha.Y[1], pri.alpha.Y[2], log = TRUE)
  l.prior2.st <- dgamma(P.star[2], pri.beta.Y[1], pri.beta.Y[2], log = TRUE)
  l.prior2    <- dgamma(P[2]     , pri.beta.Y[1], pri.beta.Y[2], log = TRUE)
  
  if(flatpri){
    logMH <- (l.lik.st - l.lik 
              + log(pmvnorm(upper = c(P[1] ,P[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] ,P.star[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }else{
    logMH <- (l.lik.st - l.lik + l.prior1.st - l.prior1 + l.prior2.st - l.prior2
              + log(pmvnorm(upper = c(P[1] ,P[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1])
              - log(pmvnorm(upper = c(P.star[1] ,P.star[2]), sigma = S%*%t(S)%*%diag(c(tun[1],tun[2])))[1]))
  }
  
  # print(S)
  alpha = min (exp(logMH),1)
  if(!is.nan(alpha) && runif(1)<alpha){
    P <- P.star;count = 1;
  } 
  S.update=ramcmc::adapt_S(S,u,alpha,rep,gamma = min(1,(2*rep)^(-2/2)))
  if(is.nan(min(S.update))) S.update=S
  S.update=(S.update+S)/2;
  return(list(P=P, S=S.update, count=count))
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
    l.lik.st <- l.lik.st + sum(Rii * log(KI.star), na.rm = T) - Ax * (sum(KI.star))
    l.lik <- l.lik + sum(Rii * log(KI.m), na.rm = T) - Ax * (sum(KI.m))
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





MH.P.X.all.share3 <- function(P,S,rep, r.X.birth1,r.X.birth2,r.X.birth3, Ax, tun, pri.alpha.X, pri.beta.X, maxt, flatpri = FALSE, lowbnds = c(0,0)) {
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
  
  for (i in 1:ncol(r.X.birth1)){
    Rii <- r.X.birth1[, i]
    l.lik.st <- l.lik.st + sum(Rii * log(KI.star), na.rm = T) - Ax[1] * (sum(KI.star))
    l.lik <- l.lik + sum(Rii * log(KI.m), na.rm = T) - Ax[1] * (sum(KI.m))
  }
  
  for (i in 1:ncol(r.X.birth2)){
    Rii <- r.X.birth2[, i]
    l.lik.st <- l.lik.st + sum(Rii * log(KI.star), na.rm = T) - Ax[2] * (sum(KI.star))
    l.lik <- l.lik + sum(Rii * log(KI.m), na.rm = T) - Ax[2] * (sum(KI.m))
  }
  
  for (i in 1:ncol(r.X.birth3)){
    Rii <- r.X.birth3[, i]
    l.lik.st <- l.lik.st + sum(Rii * log(KI.star), na.rm = T) - Ax[3] * (sum(KI.star))
    l.lik <- l.lik + sum(Rii * log(KI.m), na.rm = T) - Ax[3] * (sum(KI.m))
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






MH.P.X.all.share <- function(P,S,rep, r.X.birth, nsample, Ax, tun, pri.alpha.X, pri.beta.X, maxt, flatpri = FALSE, lowbnds = c(0,0)){
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
  
  l.lik.st <- 0
  l.lik <- 0
  
  for(gg in 1:length(maxt)){
    KI.star <- KI(P.star, maxt = maxt[gg])      # calculating kappa using current alpha & candidate of beta
    KI.m <- KI(P, maxt = maxt[gg])              # calculating kappa using current alpha & beta
    for (i in 1:nsample[gg]){
      Rii <- r.X.birth[1:maxt[gg],i,gg]
      l.lik.st <- l.lik.st + sum(Rii * log(KI.star), na.rm = T) - Ax[gg] * (sum(KI.star))
      l.lik <- l.lik + sum(Rii * log(KI.m), na.rm = T) - Ax[gg] * (sum(KI.m))
    }
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




MH.P.Y.all <- function(P,S,rep, r.Y.birth, in.X.all, Ay, K.M, tun, pri.alpha.Y, pri.beta.Y, maxt, flatpri = FALSE, lowbnds = c(0,0)) {
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
  
  l.lik.st <- 0
  l.lik <- 0
  
  for (i in 1:ncol(r.Y.birth)){
    Rii <- r.Y.birth[, i]
    KI.Y.star <- KI.Y(P.star, in.X = in.X.all[,i], K.M = K.M)      # calculating kappa using current alpha & candidate of beta
    KI.Y.m <- KI.Y(P, in.X = in.X.all[,i], K.M = K.M)              # calculating kappa using current alpha & beta
    l.lik.st <- l.lik.st + sum(Rii * log(KI.Y.star[,1]), na.rm = T) - Ay * (sum(KI.Y.star[,1]))
    l.lik <- l.lik + sum(Rii * log(KI.Y.m[,1]), na.rm = T) - Ay * (sum(KI.Y.m[,1]))
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

MH.X.ind <-function(r,x,kix,km, thetax, m = mean.x, bi=tun.X, ptnum = 4, useall = FALSE){
  count=0
  x.star = x
  prop.X.st = rep(0,length(x))
  prop.X    = rep(0,length(x))
  H0.X.st = matrix(rep(0,2*(length(x)-1)),ncol=2,nrow = (length(x)-1))
  H0.X    = matrix(rep(0,2*(length(x)-1)),ncol=2,nrow = (length(x)-1))
  
  for (i in 2:length(x)){
    x.star[i] = exp(rnorm(1, mean = m[i], sd = bi[i])) # mean of the proposal distribution should be tranformed by logarithm
    # prop.X.st[i] = dnorm(log(x.star[i]+1e-300), mean = m[i], sd = bi[i], log = T)
    # prop.X[i]    = dnorm(log(x[i]+1e-300)     , mean = m[i], sd = bi[i], log = T)
    prop.X.st[i] = dnorm(log(x.star[i]+1e-1), mean = m[i], sd = bi[i], log = T)
    prop.X[i]    = dnorm(log(x[i]+1e-1)     , mean = m[i], sd = bi[i], log = T)
    H0.X.st[i-1,] = (H.X(x.star[i-1], kix[i-1], thetax)+H.X(x.star[i], kix[i-1], thetax))/2
    H0.X[i-1,]    = (H.X(x[i-1]     , kix[i-1], thetax)+H.X(x[i]     , kix[i-1], thetax))/2
  }
  if (useall == TRUE){
    fy.st = A.Y * KI.Y(Delta.Y,in.X = x.star, K.M=km)
    fy    = A.Y * KI.Y(Delta.Y,in.X = x     , K.M=km)
  }else{
    fy.st = A.Y * KI.Ynt(Delta.Y,in.X = x.star, N = ptnum, K.M=km)
    fy    = A.Y * KI.Ynt(Delta.Y,in.X = x     , N = ptnum, K.M=km)
  }
  GlobH0.X.st <<- H0.X.st
  
  q.X.st = sum(log(dpois(r[,3:4], H0.X.st)+1e-300), na.rm = T)
  q.X    = sum(log(dpois(r[,3:4], H0.X   )+1e-300), na.rm = T)
  q.Y.st = sum(log(dpois(r[,1],fy.st[,1])+1e-300), na.rm = T)
  q.Y    = sum(log(dpois(r[,1],fy[,1]   )+1e-300), na.rm = T)
  
  prop.X.st = sum(prop.X.st, na.rm = T)
  prop.X    = sum(prop.X   , na.rm = T)
  prior.X.st = sum(dnorm(log(x.star[i]+1e-1), mean = 0, sd = 10*3, log = T)) # non-informative normal prior
  prior.X   = sum(dnorm(log(x[i]+1e-1)     , mean = 0, sd = 10*3, log = T)) # non-informative normal prior
  
  logMH <- (q.Y.st - q.Y  + q.X.st - q.X + prior.X.st - prior.X  # posterior Dist. 
            - prop.X.st + prop.X                   # proposal Dist. 
            + sum(log(x.star+1e-1)) - sum(log(x+1e-1)))       # jacobian 
  
  # 2020.03.26 Hyukypo 
  GlobX <<- x; GlobX.st <<- x.star;
  print(paste0("x.star:", x.star)); 
  print(paste0("q.Y.st:", q.Y.st))
  print(paste0("-q.Y:", -q.Y))
  print(paste0("q.X.st:", q.X.st))
  print(paste0("-q.X:", -q.X))
  print(paste0("prior.X.st:", prior.X.st))
  print(paste0("-prior.X:", -prior.X))
  print(paste0("-prop.X.st:", -prop.X.st))
  print(paste0("prop.X:", prop.X))
  print(paste0("sum(log(x.star)):", sum(log(x.star+1e-1))))
  print(paste0("-sum(log(x)):", -sum(log(x+1e-1))))
  print(paste0("logMH:", logMH))
  print("===================")
  
  if(!is.nan(logMH) && runif(1)<exp(logMH)) {
    x=x.star; count = 1;
  }    
  return(list(x=x, count=count))
}

MH.X.ind.gamma <-function(r,x,kix,km, thetax, m = mean.x, bi=tun.X, ptnum = 4, useall = FALSE){
  count=0
  x.star = x
  prop.X.st = rep(0,length(x))
  prop.X    = rep(0,length(x))
  H0.X.st = matrix(rep(0,2*(length(x)-1)),ncol=2,nrow = (length(x)-1))
  H0.X    = matrix(rep(0,2*(length(x)-1)),ncol=2,nrow = (length(x)-1))
  
  for (i in 2:length(x)){
    x.star[i] = rgamma(1, shape = m[i]^2/bi[i]^2, rate = m[i]/bi[i]^2 + 1e-1) # mean of the proposal distribution should be tranformed by logarithm
    # prop.X.st[i] = dnorm(log(x.star[i]+1e-300), mean = m[i], sd = bi[i], log = T)
    # prop.X[i]    = dnorm(log(x[i]+1e-300)     , mean = m[i], sd = bi[i], log = T)
    
    prop.X.st[i] = dgamma(x.star[i] + 1e-1, shape = m[i]^2/bi[i]^2, rate = m[i]/bi[i]^2 + 1e-1)
    prop.X[i]    = dgamma(x[i] + 1e-1, shape = m[i]^2/bi[i]^2, rate = m[i]/bi[i]^2 + 1e-1)
    
    H0.X.st[i-1,] = (H.X(x.star[i-1], kix[i-1], thetax)+H.X(x.star[i], kix[i-1], thetax))/2
    H0.X[i-1,]    = (H.X(x[i-1]     , kix[i-1], thetax)+H.X(x[i]     , kix[i-1], thetax))/2
  }
  if (useall == TRUE){
    fy.st = A.Y * KI.Y(Delta.Y,in.X = x.star, K.M=km)
    fy    = A.Y * KI.Y(Delta.Y,in.X = x     , K.M=km)
  }else{
    fy.st = A.Y * KI.Ynt(Delta.Y,in.X = x.star, N = ptnum, K.M=km)
    fy    = A.Y * KI.Ynt(Delta.Y,in.X = x     , N = ptnum, K.M=km)
  }
  # GlobH0.X.st <<- H0.X.st
  Globprop.X <<- prop.X
  
  q.X.st = sum(log(dpois(r[,3:4], H0.X.st)+1e-300), na.rm = T)
  q.X    = sum(log(dpois(r[,3:4], H0.X   )+1e-300), na.rm = T)
  q.Y.st = sum(log(dpois(r[,1],fy.st[,1])+1e-300), na.rm = T)
  q.Y    = sum(log(dpois(r[,1],fy[,1]   )+1e-300), na.rm = T)
  
  prop.X.st = sum(prop.X.st, na.rm = T)
  prop.X    = sum(prop.X   , na.rm = T)
  prior.X.st = sum(dgamma(x.star + 1e-1 , shape = 1, rate = 1e-2)) # non-informative gamma prior
  prior.X   = sum(dgamma(x + 1e-1 , shape = 1, rate = 1e-2)) # non-informative gamma prior
  
  logMH <- (q.Y.st - q.Y  + q.X.st - q.X + prior.X.st - prior.X  # posterior Dist. 
            - prop.X.st + prop.X)                   # proposal Dist. 
  
  # 2020.03.26 Hyukypo 
  GlobX <<- x; GlobX.st <<- x.star;
  print(paste0("x.star:", x.star)); 
  print(paste0("q.Y.st:", q.Y.st))
  print(paste0("-q.Y:", -q.Y))
  print(paste0("q.X.st:", q.X.st))
  print(paste0("-q.X:", -q.X))
  print(paste0("prior.X.st:", prior.X.st))
  print(paste0("-prior.X:", -prior.X))
  print(paste0("-prop.X.st:", -prop.X.st))
  print(paste0("prop.X:", prop.X))
  print(paste0("logMH:", logMH))
  print("===================")
  
  if(!is.nan(logMH) && runif(1)<exp(logMH)) {
    x=x.star; count = 1;
  }    
  return(list(x=x, count=count))
}

MH.XR.ind.gamma <-function(r,x,kix,km, thetax, m = mean.x, bi, ptnum = 4, useall = TRUE, A.Y){
  x.star <- x
  rx.st <- matrix(0, ncol=2,nrow = (length(x)-1))
  # rx <- matrix(0, ncol=2,nrow = (length(x)-1))
  
  prop.X.st = rep(0,length(x))
  prop.X    = rep(0,length(x))
  H0.X.st = matrix(0, ncol=2,nrow = (length(x)-1))
  H0.X    = matrix(0, ncol=2,nrow = (length(x)-1))
  
  for (i in 2:length(x)){
    x.star[i] = rgamma(1, shape = m[i]^2/bi[i]^2, rate = m[i]/bi[i]^2 + 1e-1) # x.star: a candidate for the next X.
    prop.X.st[i] = dgamma(x.star[i] + 1e-2, shape = m[i]^2/bi[i]^2, rate = m[i]/bi[i]^2 + 1e-1)
    prop.X[i]    = dgamma(x[i] + 1e-2, shape = m[i]^2/bi[i]^2, rate = m[i]/bi[i]^2 + 1e-1)
    
    H0.X.st[i-1,] = (H.X(x.star[i-1], kix[i-1], thetax)+H.X(x.star[i], kix[i-1], thetax))/2
    H0.X[i-1,]    = (H.X(x[i-1]     , kix[i-1], thetax)+H.X(x[i]     , kix[i-1], thetax))/2
    #H0.X = (birth propensity of X, death propensity of X) at time i.
  }
  if (useall == TRUE){
    fy.st = A.Y * KI.Y(Delta.Y,in.X = x.star, K.M=km)
    fy    = A.Y * KI.Y(Delta.Y,in.X = x     , K.M=km)
  }else{
    fy.st = A.Y * KI.Ynt(Delta.Y,in.X = x.star, N = ptnum, K.M=km)
    fy    = A.Y * KI.Ynt(Delta.Y,in.X = x     , N = ptnum, K.M=km)
  }
  gg<<-x.star
  rx.st <- impute_r.X(x.star, B.X = thetax[2])
  rx <- r[,3:4]
  if(rx.st[1,1] == -1){
    return(list(x=x, count = 0, rx = rx, errflg = 1))
  }
  
  q.X.bir.st = sum(log(dpois(rx.st[,1], H0.X.st[,1]) + 1e-300), na.rm = T)
  q.X.bir    = sum(log(dpois(rx[,1], H0.X[,1]   ) + 1e-300), na.rm = T)
  
  # q.X.dea.st = sum(log(dpois(rx.st[,2], H0.X.st[,2])+1e-300), na.rm = T)
  # q.X.dea    = sum(log(dpois(rx[,2], H0.X[,2]   )+1e-300), na.rm = T)
  
  q.X.dea.st <- 0
  q.X.dea    <- 0
  
  q.Y.st = sum(log(dpois(r[,1],fy.st[,1])+1e-300), na.rm = T)
  q.Y    = sum(log(dpois(r[,1],fy[,1]   )+1e-300), na.rm = T)
  
  prop.X.st = sum(log(prop.X.st + 1e-300), na.rm = T);# prop.X.st <- 0;
  prop.X    = sum(log(prop.X + 1e-300)   , na.rm = T);# prop.X <- 0;
  
  # prior.X.st = sum(dgamma(x.star[i] + 1e-1 , shape = 1, rate = 1e-2)) # non-informative gamma prior
  # prior.X   = sum(dgamma(x[i] + 1e-1 , shape = 1, rate = 1e-2)) # non-informative gamma prior
  
  prior.X.st = sum(log(dgamma(x.star , shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
  prior.X   = sum(log(dgamma(x, shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
  
  # logMH <- (q.Y.st - q.Y  + q.X.bir.st + q.X.dea.st - q.X.bir - q.X.dea + prior.X.st - prior.X  # posterior Dist. 
  #           - prop.X.st + prop.X - prop.rx.st + prop.rx)                   # proposal Dist. 
  logMH <- (q.Y.st - q.Y  + q.X.bir.st + q.X.dea.st - q.X.bir - q.X.dea + prior.X.st - prior.X  # posterior Dist. 
            - prop.X.st + prop.X)
  count <- 0;
  if(!is.nan(logMH) && runif(1)<exp(logMH)) {
    x=x.star; count = 1; rx <- rx.st
  }    
  
  # return(list(x=x, count=count, rx=rx,errflg = 0))
  return(list(x=x, count=count, rx=rx, q.X.bir.st = q.X.bir.st, q.X.dea.st = q.X.dea.st, q.X.bir = q.X.bir, q.X.dea = q.X.dea,
              q.Y.st=q.Y.st, q.Y=q.Y, prop.X.st=prop.X.st, prop.X = prop.X,
              prior.X.st = prior.X.st, prior.X = prior.X, errflg = 0))
}



MH.XR.ind.gamma.all <-function(r.all , x.all ,kix,km, thetax, m = mean.x, bi=tun.X, ptnum = 4, useall = TRUE, A.Y){
  # r.all: (data.num) * (max.T) * 2 matrix. The birth and death numbers of Y's;
  # x.all: (data.num) * (max.T+1) matrix. X trajectory for each Y.
  data.num <- dim(r.all)[1];
  maxt <- dim(r.all)[2];
  
  x.all.star = x.all
  # rx.st <- matrix(0, ncol=2, nrow = maxt)
  # rx <- matrix(0, ncol=2, nrow = maxt)
  rx.st <- array(0, dim = c(data.num, maxt, 2))
  rx <- array(0, dim = c(data.num, maxt, 2))
  
  prop.X.st = rep(0,length(x))
  prop.X    = rep(0,length(x))
  H0.X.st = array(0, dim = c(data.num, maxt, 2))
  H0.X    = array(0, dim = c(data.num, maxt, 2))
  
  fy <- 0; fy.st <- 0;
  q.X.bir.st <- 0; q.X.bir <- 0;
  q.X.dea.st <- 0; q.X.dea <- 0;
  q.Y.st <- 0; q.Y <- 0;
  prior.X.st <- 0; prior.X <- 0;   
  
  
  for(j in 1:data.num){
    for(i in 2:(maxt+1)){
      x.all.star[j,i] = rgamma(1, shape = m[i]^2/bi[i]^2, rate = m[i]/bi[i]^2 + 1e-1) # x.star: a candidate for the next X.
      # prop.X.st[i] = dnorm(log(x.star[i]+1e-300), mean = m[i], sd = bi[i], log = T)
      # prop.X[i]    = dnorm(log(x[i]+1e-300)     , mean = m[i], sd = bi[i], log = T)
      
      prop.X.st[i] = prop.X.st[i] + dgamma(x.all.star[j,i] + 1e-2, shape = m[i]^2/bi[i]^2, rate = m[i]/bi[i]^2 + 1e-1)
      prop.X[i]    = prop.X[i] + dgamma(x.all[j,i] + 1e-2, shape = m[i]^2/bi[i]^2, rate = m[i]/bi[i]^2 + 1e-1)
      
      H0.X.st[j,i-1,] = H0.X.st[j,i-1,] + (H.X(x.all.star[j,i-1], kix[i-1], thetax)+H.X(x.all.star[j,i], kix[i-1], thetax))/2
      H0.X[j,i-1,]    = H0.X[j,i-1,] + (H.X(x.all[j,i-1]     , kix[i-1], thetax)+H.X(x.all[j,i]     , kix[i-1], thetax))/2
      #H0.X = (birth propensity of X, death propensity of X) at time i.
    }
    if (useall == TRUE){
      fy.st = fy.st + A.Y * KI.Y(Delta.Y,in.X = x.all.star[j,], K.M=km)
      fy    = fy + A.Y * KI.Y(Delta.Y,in.X = x.all[j,]     , K.M=km)
    }else{
      fy.st = fy.st + A.Y * KI.Ynt(Delta.Y,in.X = x.all.star[j,], N = ptnum, K.M=km)
      fy    = fy + A.Y * KI.Ynt(Delta.Y,in.X = x.all[j,]     , N = ptnum, K.M=km)
    }

    rx.st[j,,] <- impute_r.X(x.all.star[j,], B.X = thetax[2])
    rx[j,,] <- impute_r.X(x.all[j,], B.X = thetax[2])
    if(rx.st[j,1,1] == -1){
      return(list(x=x.all, count = 0, rx = r[,,3:4], errflg = 1))
    }else if(rx[j,1,1] == -1 & rx.st[1,1] != -1){
      return(list(x=x.all.star, count= 1 , rx = rx.st, errflg = 2))
    }
    
    q.X.bir.st = q.X.bir.st + sum(log(dpois(rx.st[j,,1], H0.X.st[j,,1])+1e-300), na.rm = T)
    q.X.bir    = q.X.bir    + sum(log(dpois(rx[j,,1], H0.X[j,,1]   )+1e-300), na.rm = T)
    # q.X.dea.st = q.X.dea.st + sum(log(dpois(rx.st[j,,1], H0.X.st[j,,2])+1e-300), na.rm = T)
    # q.X.dea    = q.X.dea    + sum(log(dpois(rx[j,,1], H0.X[j,,2]   )+1e-300), na.rm = T)
    q.X.dea.st <- 0;
    q.X.dea <- 0;
    
    q.Y.st = q.Y.st + sum(log(dpois(r[,1],fy.st[,1])+1e-300), na.rm = T)
    q.Y    = q.Y    + sum(log(dpois(r[,1],fy[,1]   )+1e-300), na.rm = T)
    

  }
  
  prior.X.st = sum(log(dgamma(x.all.star, shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
  prior.X   = sum(log(dgamma(x.all, shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
  prop.X.st = sum(log(prop.X.st + 1e-300), na.rm = T)
  prop.X    = sum(log(prop.X + 1e-300)   , na.rm = T)
  
  logMH <- (q.Y.st - q.Y  + q.X.bir.st + q.X.dea.st - q.X.bir - q.X.dea + prior.X.st - prior.X  # posterior Dist. 
            - prop.X.st + prop.X)                   # proposal Dist. 
  count <- 0;
  if(!is.nan(logMH) && runif(1)<exp(logMH)) {
    x=x.all.star; count = 1; rx <- rx.st
  }    
  # return(list(x=x, count=count, rx=rx,errflg = 0))
  return(list(x=x, count=count, rx=rx, q.X.bir.st = q.X.bir.st, q.X.dea.st = q.X.dea.st, q.X.bir = q.X.bir, q.X.dea = q.X.dea,
              q.Y.st=q.Y.st, q.Y=q.Y, prop.X.st=prop.X.st, prop.X = prop.X,
              prior.X.st = prior.X.st, prior.X = prior.X, errflg = 0))
}

MH.XR.ind <-function(r,x,kix,km, thetax, m = mean.x, bi=tun.X, ptnum = 4, useall = FALSE){
  count=0
  x.star = x
  
  rx.st <- matrix(rep(0,2*(length(x)-1)),ncol=2,nrow = (length(x)-1))
  rx <- matrix(rep(0,2*(length(x)-1)),ncol=2,nrow = (length(x)-1))
  
  prop.X.st = rep(0,length(x))
  prop.X    = rep(0,length(x))
  H0.X.st = matrix(rep(0,2*(length(x)-1)),ncol=2,nrow = (length(x)-1))
  H0.X    = matrix(rep(0,2*(length(x)-1)),ncol=2,nrow = (length(x)-1))
  
  for (i in 2:length(x)){
    x.star[i] = exp(rnorm(1, mean = m[i], sd = bi[i])) # mean of the proposal distribution should be tranformed by logarithm
    # prop.X.st[i] = dnorm(log(x.star[i]+1e-300), mean = m[i], sd = bi[i], log = T)
    # prop.X[i]    = dnorm(log(x[i]+1e-300)     , mean = m[i], sd = bi[i], log = T)
    prop.X.st[i] = dnorm(log(x.star[i]+1e-1), mean = m[i], sd = bi[i], log = T)
    prop.X[i]    = dnorm(log(x[i]+1e-1)     , mean = m[i], sd = bi[i], log = T)
    H0.X.st[i-1,] = (H.X(x.star[i-1], kix[i-1], thetax)+H.X(x.star[i], kix[i-1], thetax))/2
    H0.X[i-1,]    = (H.X(x[i-1]     , kix[i-1], thetax)+H.X(x[i]     , kix[i-1], thetax))/2
  }
  if (useall == TRUE){
    fy.st = A.Y * KI.Y(Delta.Y,in.X = x.star, K.M=km)
    fy    = A.Y * KI.Y(Delta.Y,in.X = x     , K.M=km)
  }else{
    fy.st = A.Y * KI.Ynt(Delta.Y,in.X = x.star, N = ptnum, K.M=km)
    fy    = A.Y * KI.Ynt(Delta.Y,in.X = x     , N = ptnum, K.M=km)
  }
  rx.st <- impute_r.X(x.star)
  rx <- impute_r.X(x)
  
  GlobH0.X.st <<- H0.X.st
  
  q.X.st = sum(log(dpois(rx.st[,1:2], H0.X.st)+1e-300), na.rm = T)
  q.X    = sum(log(dpois(rx[,1:2], H0.X   )+1e-300), na.rm = T)
  q.Y.st = sum(log(dpois(r[,1],fy.st[,1])+1e-300), na.rm = T)
  q.Y    = sum(log(dpois(r[,1],fy[,1]   )+1e-300), na.rm = T)
  prop.X.st = sum(prop.X.st, na.rm = T)
  prop.X    = sum(prop.X   , na.rm = T)
  
  prior.X.st = sum(dnorm(log(x.star[i]+1e-1), mean = 0, sd = 10*3, log = T)) # non-informative normal prior
  prior.X   = sum(dnorm(log(x[i]+1e-1)     , mean = 0, sd = 10*3, log = T)) # non-informative normal prior
  
  logMH <- (q.Y.st - q.Y  + q.X.st - q.X + prior.X.st - prior.X  # posterior Dist. 
            - prop.X.st + prop.X                   # proposal Dist. 
            + sum(log(x.star+1e-1)) - sum(log(x+1e-1)))       # jacobian 
  
  # 2020.03.26 Hyukypo 
  GlobX <<- x; GlobX.st <<- x.star;
  print(paste0("x.star:", x.star)); 
  print(paste0("q.Y.st:", q.Y.st))
  print(paste0("-q.Y:", -q.Y))
  print(paste0("q.X.st:", q.X.st))
  print(paste0("-q.X:", -q.X))
  print(paste0("prior.X.st:", prior.X.st))
  print(paste0("-prior.X:", -prior.X))
  print(paste0("-prop.X.st:", -prop.X.st))
  print(paste0("prop.X:", prop.X))
  print(paste0("sum(log(x.star)):", sum(log(x.star+1e-1))))
  print(paste0("-sum(log(x)):", -sum(log(x+1e-1))))
  print(paste0("logMH:", logMH))
  print("===================")
  
  if(!is.nan(logMH) && runif(1)<exp(logMH)) {
    x=x.star; count = 1;
  }    
  return(list(x=x, count=count))
}

##################################################################################################
#The rest of the functions were not used in this study where delay only occurs in X.
##################################################################################################


# Metropolis-Hastings algorithm for updating number of reaction of X given the trajectory of X.
# It is used for checking only 

# MH.R <-function(r,x,ki, b=tun.B){
#   count = matrix(0, nrow = maxt, ncol = 1)
#   for (i in 1:(length(X)-1)) {
#     r.X = c(r[i,3],r[i,4])
#     r.Y = c(r[i,1],r[i,2])  
#     xi=x[i];xi1=x[i+1]
#     k <- ki[i]
#     #bi=b[i,]
#     
#     r3.star = cand(r.X[1],b[3])
#     r4.star = r3.star - (xi1-xi)
#     r.X.star = c(r3.star,r4.star)
#     #if(min(r.X.star)>=0){
#       prop.r3.st = prop.R(r.X.star[1], r.X[1], b[3])
#       prop.r3    = prop.R(r.X[1]     , r.X[1], b[3])
#       H0.X    = (H.X(xi, k, theta.X)+H.X(xi1, k, theta.X))/2
#       q.X.st = sum(dpois(r.X.star,H0.X, log = T))
#       q.X    = sum(dpois(r.X     ,H0.X, log = T))
#       logMH = prop.r3 - prop.r3.st + q.X.st - q.X 
#       if(!is.nan(logMH) && runif(1)<exp(logMH)) {
#         r[i,3:4]=r.X.star; count[i] = 1;
#       }    
#     #}
#   }
#   return(list(r=r, count=count)) 
# }

# Metropolis-Hastings algorithm for updating the trajectory of X given the number of reaction of X.
# This function does not work. Need to correction. 

# MH.X <- function(x, R=RR[,c(1,4)],tun=tun.X, B = B.X, A=A.Y, KM=K.M){
#   count = rep(0,maxt)
#   for(i in 2:(maxt)){
#     x.m = x[i]
#     lambda.i = rgamma(1,R[i,1]*tun,rate = tun)
#     x.star = max(round(lambda.i*KM/(A-lambda.i)),0)
#     l.lik.st = R[i,2] * log(B*x.star) - B*x.star + R[i,1] * log(A*x.star/(KM+x.star)) - A*x.star/(KM+x.star)
#     l.lik    = R[i,2] * log(B*x.m) - B*x.m + R[i,1] * log(A*x.star/(KM+x.m)) - A*x.m/(KM+x.m)
#     logMH = l.lik.st - l.lik - dgamma(lambda.i,R[i,1]*tun, rate = tun, log = T) + dgamma( A*x.m/(KM+x.m),R[i,1]*tun,rate = tun, log = T)
#     if(!is.nan(logMH) && runif(1)<exp(logMH)) {
#       x[i] = x.star; count[i] = 1;
#     }    
#   }
#   return(list(x=x,count=count))
# } 



# RAM (robust adaptive Metropolis) method for KM constant
# refering Vihola(2012), Statistics and Computing, 22(5):997-1008. 
# using informative gamma prior for beta

MH.KM <- function(km, s, rep, r, x, b, pri.KM, Delta.Y, flatpri = FALSE){
  # km: Michaelis-Menten constant in the current iteration.
  # s: scaling factor for tunning
  # rep: iteration number
  # x: X trajectory to construct the likelihood for the birth of Y 
  # r: numbers of the birth reactions of Y.
  # b: tunning parameters for KM
  # pri.KM : prior distribution parameters for KM
  
  count = 0
  repeat{
    u = rnorm(1,0,b)
    km.star = km + s * u
    if(km.star > 0) break
  }
  # lambda.st = ((A.Y*x[-length(x)])/(as.vector(km.star) + x[-length(x)]) 
  #              + (A.Y*x[-1])/(as.vector(km.star) + x[-1]))/2
  # lambda    = ((A.Y*x[-length(x)])/(km      + x[-length(x)]) 
  #              + (A.Y*x[-1])/(km      + x[-1]))/2
  lambda.st = A.Y * KI.Y(Delta.Y,x, K.M=km.star)
  lambda    = A.Y * KI.Y(Delta.Y,x     , K.M=km)
  l.lik.st = sum(log(dpois(r,lambda.st[,1])+1e-300))
  l.lik    = sum(log(dpois(r,lambda[,1])+1e-300))
  
  if(flatpri){
    logMH = (l.lik.st - l.lik #+ dgamma(km.star, pri.KM[1], pri.KM[2],log=T) - dgamma(km, pri.KM[1], pri.KM[2],log=T)
             + pnorm(km, 0, s*b, log.p = T) - pnorm(km.star, 0, s*b, log.p = T))  
  }else{
    logMH = (l.lik.st - l.lik + dgamma(km.star, pri.KM[1], pri.KM[2],log=T) - dgamma(km, pri.KM[1], pri.KM[2],log=T)
             + pnorm(km, 0, s*b, log.p = T) - pnorm(km.star, 0, s*b, log.p = T))  
  }
  
  
  if(!is.nan(logMH) && runif(1)<exp(logMH)){
    km=km.star; count = 1;
  }
  alpha = min(exp(logMH),1)
  s = ramcmc::adapt_S(s,u,alpha,rep,gamma = min(1,(1*rep)^(-2/3)))
  return(list(km=km, s=s, count=count, l.lik = l.lik, l.lik.st = l.lik.st))
}


MH.KM.all <- function(km, s, rep, r.all, x.all, b, pri.KM, Delta.Y, flatpri = FALSE){
  # km: Michaelis-Menten constant in the current iteration.
  # s: scaling factor for tunning
  # rep: iteration number
  # x.all: X trajectories to construct the likelihood for the birth of Y. Dim: (maxT + 1) * (number of experiments)
  # r.all: numbers of the birth reactions of Y. Dim: maxT * (number of experiments)
  # b: tunning parameters for KM
  # pri.KM : prior distribution parameters for KM
  
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
    lambda.st = A.Y * KI.Y(Delta.Y,x, K.M=km.star)
    lambda    = A.Y * KI.Y(Delta.Y,x,      K.M=km)
    l.lik.st = l.lik.st + sum(log(dpois(r,lambda.st[,1])+1e-300))
    l.lik    = l.lik    + sum(log(dpois(r,lambda[,1]   )+1e-300))
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

MH.KM.singleX <- function(km, s, rep, r.all, x, b, pri.KM, Delta.Y, flatpri = FALSE){
  # km: Michaelis-Menten constant in the current iteration.
  # s: scaling factor for tunning
  # rep: iteration number
  # x: X trajectory to construct the likelihood for the birth of Y. 
  # r.all: numbers of the birth reactions of Y. Dim: maxT * (number of experiments)
  # b: tunning parameters for KM
  # pri.KM : prior distribution parameters for KM
  
  count = 0
  repeat{
    u = rnorm(1,0,b)
    km.star = km + s * u
    if(km.star > 0) break
  }
  
  l.lik.st = 0;
  l.lik = 0;
  
  for(j in 1:ncol(r.all)){
    r = r.all[,j]
    lambda.st = A.Y * KI.Y(Delta.Y,x, K.M=km.star)
    lambda    = A.Y * KI.Y(Delta.Y,x,      K.M=km)
    l.lik.st = l.lik.st + sum(log(dpois(r,lambda.st[,1])+1e-300))
    l.lik    = l.lik    + sum(log(dpois(r,lambda[,1]   )+1e-300))
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



MH.A.X <- function(A,s,rep, r=RR[,3], g=g_11, b=1){
  count=0
  repeat{
    u = rnorm(1,0,b)
    A.star = A + s * u
    if(A.star > 0) break
  }
  l.lik.st = dgamma(A.star, shape = sum(r), rate = g_11, log = T)
  l.lik    = dgamma(A     , shape = sum(r), rate = g_11, log = T)
  
  logMH = (l.lik.st - l.lik + dgamma(A.star, a_A, b_A, log = T) -dgamma(A, a_A, b_A, log = T) 
           + pnorm(A, 0, s*b, log.p = T) - pnorm(A.star, 0, s*b, log.p = T)  ) 
  if(!is.nan(logMH) && runif(1)<exp(logMH)) {
    A=A.star; count = 1;
  }
  alpha = min(exp(logMH),1)
  s=ramcmc::adapt_S(s,u,alpha,rep,gamma = min(1,(1*rep)^(-2/3)))
  return(list(A=A,s=s, count=count))
}

log_factorial <- function(x){
  val = log(gamma(x+2 - ceiling(x)));
  if(length(val)==1){
    if (x >=1){
      for (jj in 0:(ceiling(x)-1)){
        val = val + log(x-jj)
      }
    }
  }else{
    for (ii in 1:length(val)){
      if (x[ii] >=1){
        for (jj in 0:(ceiling(x[ii])-1)){
          val[ii] = val[ii] + log(x[ii]-jj)
        }
      } 
    }
  }
  return(val)
}

