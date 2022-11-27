setwd("~/Dropbox/Twostep_delay/TwoStepDelay/")
source('2stepDDE_functions_public.R')

# library("invgamma")
# rndseed <- round((as.numeric(Sys.time())*1000)%% 10000)

effrepeat = 110 # the number of desired posterior samples.
max.T = 50
nsample = 5
t_int = 1 # time interval between measurements. 

lowbnds = c(1,0.001) # lower bounds for the shape parameter and the rate parameter of gamma delay distribution

param_est <- c(1,0,1,1,1,0,0,0,0) # AX KM alphaX betaX AY alphaY betaY B sigma
# if the entry is 1 then the corresponding parameter is estimated
# otherwise, it is fixed.

# true parameter value setting 
B<- 0.05*t_int;  
A.Y <- 60*t_int; 
alpha.Y <- 3.6; beta.Y <- 0.6*t_int;  
A.X <- 10*t_int; alpha.X <- 3.6; beta.X <- 0.6*t_int; 
K.M <- 100; 
var.noise <- 10;
max.T <- max.T/t_int; # rescale the time according to the time interval between measurements.

tspan <- 0:max.T

birthX.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true birth number of X
deathX.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true death number of X
birthY.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true birth number of Y
deathY.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true death number of Y
X.all <- matrix(0, nrow = max.T+1, ncol = nsample)
Y.all <- matrix(0, nrow = max.T+1, ncol = nsample)

obs.Y.all <- matrix(0, nrow = max.T+1, ncol = nsample)

obs.Y.all <- read.csv("Example_data.csv", header = FALSE) # User input data.
#Each column of 'obs.Y.all' is observed data at time 0, 1, 2, ..., max.T

pri.A.X <- c(0.001, 0.001); # non-informative prior for A.X
pri.alpha.X <- c(0.001, 0.001); # non-inormative prior for alpha.X
pri.beta.X <- c(0.001, 0.001); # non-inormative prior for beta.X
pri.KM <- c(0.001, 0.001); # non-informative prior for KM
pri.A.Y <- c(0.001, 0.001); # non-informative prior for A.X
pri.alpha.Y <- c(0.001, 0.001); # non-inormative prior for alpha.X
pri.beta.Y <- c(0.001, 0.001); # non-inormative prior for beta.X
pri.B <- c(0.001, 0.001); # non-informative prior for KM
pri.var.noise <- c(0.001, 0.001); # non-informative prior for var.noise

tun.KM <- 1; # tuning parameter of MH for K_M
tun.Delta.X <- c(1.0, 1); # tuning parameter of MH for (alpha_X, beta_X)
tun.Delta.Y <- c(1.0, 1); # tuning parameter of MH for (alpha_X, beta_X)
tun.sigma <- 40; # tuning parameter for the acceptance of proposed trajectories. Higher values give a larger acceptance rate.

burn <- 0; # the length of burn-in period for MCMC iteration.
thin <- 1; # the reciprocal of the thining rate for MCMC iteration.
nrepeat <- burn + thin*effrepeat; # the total number of MCMC iteration.

selrow <- seq(from = burn + thin, by = thin, length.out = effrepeat) # remaining iteration after burning and thining.

# count variables represent the number of acceptance during the MCMC iteration.
count_KM <- 0; 
count_X <- rep(0, nsample); 
count_Delta.X <- 0; 
count_Delta.Y <- 0;  

theta <- matrix(0,nrow = nrepeat, ncol = 9) # matrix saving the sampled parameters from MCMC iteration.

#initial scales of KM, delta.X(=alpha.X and beta.X), and delta.Y(=alpha.Y and beta.Y) used in RAM method 
KM.S <- 10
Delta.X.S <- diag(2)
Delta.Y.S <- diag(2)

# initial parameter value for MCMC algorithm 
theta[1,] = c(A.X, K.M, alpha.X, beta.X, A.Y, alpha.Y, beta.Y, B, var.noise)


RR.all = array(0, dim = c(max.T, 4, nsample)) #saving number of reaction 

for(jj in 1:nsample){
  myListXY <- TimeDelayGillespieforXY(A.X = theta[1,1], B.X = theta[1,8], alpha.X = theta[1,3], beta.X = theta[1,4], 
                                      A.Y = theta[1,5], B.Y = theta[1,8], alpha.Y = theta[1,6], beta.Y = theta[1,7], 
                                      K.M = theta[1,2], repnum = max.T*10000, maxT = max.T+3)
  
  X.bir.st <- myListXY$Xbirth[1:max.T]
  X.dea.st <- myListXY$Xdeath[1:max.T]
  Y.bir.st <- myListXY$Ybirth[1:max.T]
  Y.dea.st <- myListXY$Ydeath[1:max.T]
  
  X.star <- c(0, cumsum(X.bir.st - X.dea.st));
  Y.star <- c(0, cumsum(Y.bir.st - Y.dea.st));
  
  X.all[(0:max.T)+1,jj] <- X.star
  Y.all[(0:max.T)+1,jj] <- Y.star
  RR.all[1:max.T,1,jj] <- Y.bir.st
  RR.all[1:max.T,2,jj] <- Y.dea.st
  RR.all[1:max.T,3,jj] <- X.bir.st
  RR.all[1:max.T,4,jj] <- X.dea.st
}

for(rep in 2:nrepeat){
  # step 1 & 2: sampling  r2 and r1 (death and birth of Y)
  # step 3: sampling X, Y & r1, r2, r3, and r4
  # updating X and Y using independent chain MH
  
  # generate a proposal mean trajectory using the current parameter set.
  for(jj in 1:nsample){
    
    myListXY <- tauleapingforXY(A.X = theta[rep-1,1], B.X = theta[rep-1,8], alpha.X = theta[rep-1,3], beta.X = theta[rep-1,4], 
                                A.Y = theta[rep-1,5], B.Y = theta[rep-1,8], alpha.Y = theta[rep-1,6], beta.Y = theta[rep-1,7], 
                                K.M = theta[rep-1,2], maxT = max.T)
    
    X.bir.st <- myListXY$Xbirth
    X.dea.st <- myListXY$Xdeath
    Y.bir.st <- myListXY$Ybirth
    Y.dea.st <- myListXY$Ydeath
    
    X.star <- c(0, cumsum(X.bir.st - X.dea.st));
    Y.star <- c(0, cumsum(Y.bir.st - Y.dea.st));
    
    
    q.Y = sum(dnorm(x = obs.Y.all[,jj], mean = Y.all[,jj], sd = sqrt(theta[rep-1,9]) + tun.sigma, log = TRUE))
    q.Y.st = sum(dnorm(x = obs.Y.all[,jj], mean = Y.star, sd = sqrt(theta[rep-1,9]) + tun.sigma, log = TRUE))
    
    
    # logMH <- q.Y.st - q.Y + prior.X.st - prior.X; # considering prior.
    logMH <- q.Y.st - q.Y; # Completely non-informative, i.e., always prior.X.st == prior.X  // improper prior
    
    if(!is.nan(logMH) && runif(1)<exp(logMH)){
      X.all[,jj] <- X.star
      Y.all[,jj] <- Y.star
      RR.all[,1,jj] <- Y.bir.st
      RR.all[,2,jj] <- Y.dea.st
      RR.all[,3,jj] <- X.bir.st
      RR.all[,4,jj] <- X.dea.st
      count_X[jj] = count_X[jj] + 1
    }
  }
  # cat(mean(count_X))
  
  
  # step  4: samping A.X 
  if (param_est[1] == 0){
    theta[rep,1] = theta[rep-1,1]
  }else{
    K.i <- KI(P = theta[rep-1,3:4], maxt = max.T);
    g_11 <- sum(K.i);
    theta[rep,1] = rgamma(1,shape = sum(RR.all[,3,]) + pri.A.X[1], rate = nsample * g_11 + pri.A.X[2]);
  }
  
  # step 5 & 6: sampling alpha.X and beta.X: the delay parameters for the birth reaction of X.
  if (param_est[3] == 0){
    count_Delta.X = count_Delta.X
    theta[rep, 3:4] = theta[rep-1, 3:4]
  }else{
    p.update <- MH.P.X.all(P = theta[rep-1,3:4], Delta.X.S, rep, RR.all[,3,], Ax = theta[rep,1],  tun = tun.Delta.X,
                           pri.alpha.X = pri.alpha.X, pri.beta.X = pri.beta.X, maxt = max.T, lowbnds = lowbnds)
    theta[rep,3:4] = p.update$P
    Delta.X.S = p.update$S
    count_Delta.X = count_Delta.X + p.update$count
  }
  
  # step 7: sampling the Michaelis-Menten constant K.M
  if (param_est[2] == 0){
    theta[rep,2] = theta[rep-1,2]
  }else{
    KM.update = MH.KM.all(theta[rep-1,2] , KM.S, rep, RR.all[,1,], X.all, b = tun.KM, pri.KM = pri.KM, Delta.Y = c(theta[rep-1,6], theta[rep-1,7]), ay = theta[rep-1,5], noise_add = FALSE)
    theta[rep,2] = KM.update$km
    KM.S = KM.update$s
    count_KM = count_KM + KM.update$count
  }
  
  # step 8: sampling alpha.Y and beta.Y: the delay parameters for the birth reaction of Y.
  if (param_est[6] == 0){
    theta[rep,6:7] = theta[rep-1, 6:7]
    count_Delta.Y = count_Delta.Y
  }else{
    p.update <- MH.P.Y.all(P = theta[rep-1,6:7], S = Delta.Y.S, rep = rep, r.Y.birth = RR.all[,1,], in.X.all = X.all, 
                           Ay = theta[rep-1,5], K.M = theta[rep, 2], tun = tun.Delta.Y, pri.alpha.Y = pri.alpha.Y, pri.beta.Y = pri.beta.Y, maxt = max.T, lowbnds = lowbnds, noise_add = FALSE)
    theta[rep,6:7] = p.update$P
    Delta.Y.S = p.update$S
    count_Delta.Y = count_Delta.Y + p.update$count
  }
  
  # 
  # step 9: sampling A.Y
  if(param_est[5] == 0){
    theta[rep,5] = theta[rep-1,5]
  }else{
    KY.sum <- 0
    for(ii in 1:nsample){
      KY.i <- KI.Y(P = theta[rep,6:7], in.X = X.all[,ii], K.M = theta[rep,2])
      KY.sum <- KY.sum + sum(KY.i[,1])
    }
    theta[rep,5] = rgamma(1,shape = sum(RR.all[,1,]) + pri.A.Y[1], rate = KY.sum + pri.A.Y[2]);
  }
  
  # step 10: sampling B(=B.X=B.Y); the common dilution rate 
  if(param_est[8] == 0){
    theta[rep,8] = theta[rep-1,8]
  }else{
    theta[rep,8] = rgamma(1,shape = sum(RR.all[,2,]) + sum(RR.all[,4,]) + pri.B[1],
                          rate = sum(X.all) + sum(Y.all) - 0.5*sum(X.all[max.T+1,]) - 0.5*sum(Y.all[max.T+1,])+ pri.B[2]);
  }
  
  # step 11: sampling var.noise ; the variance of observation noise.
  if(param_est[9] == 0){
    theta[rep,9] = theta[rep-1,9]
  }else{
    theta[rep,9] = rinvgamma(n=1, shape = pri.var.noise[1] + nsample*max.T/2, rate = pri.var.noise[2] + sum((obs.Y.all-Y.all)^2)/2)
  }
  
  if(rep%%100 ==0){
    cat(rep)
    cat(" ")
  } # check if the code is running properly.
  
}

save.image(file = paste("Estimation_results2.RData", sep = ""))
