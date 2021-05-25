rm(list=ls())
cat("\014")  

setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode")
source('2stepDDE_functions.R')

rndseed <- round((as.numeric(Sys.time())*1000)%% 10000)
set.seed(rndseed)


int <- 1; 
# known (fixed) parameters
B.X <- 0.05*int;  
B.Y <- 0.05*int; 
A.Y <- 60*int; 
alpha.Y <- 3.6; beta.Y <- 0.6*int;  

# unknowns parameters
A.X <- 10*int; alpha.X <- 3.6; beta.X <- 0.6*int; 
K.M <- 200; 

max.T <- 150 # simulated data will be given from t = 0, ..., max.T
tspan <- 0:max.T
nsample <- 10
vol.gen <- 1
vol.est <- 1
scaling.rate <- 1

birthX.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true birth number of X
deathX.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true death number of X
birthY.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true birth number of Y
deathY.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true death number of Y
sim.X.all <- matrix(0, nrow = max.T+1, ncol = nsample)
sim.Y.all <- matrix(0, nrow = max.T+1, ncol = nsample)

for(jj in 1:nsample){
  # myList is raw simulated data. 
  myList <- TimeDelayGillespieforXY(A.X = A.X, B.X = B.X, alpha.X = alpha.X, beta.X = beta.X, A.Y = A.Y, B.Y = B.Y, alpha.Y = alpha.Y, beta.Y = beta.Y, K.M = K.M, repnum = max.T*10000, maxT = max.T+3, Volume = vol.gen)
  # sim.X is true X data, and sim.Y is true Y data.
  birthX.sim[,jj] <- myList$Xbirth[1:max.T]
  deathX.sim[,jj] <- myList$Xdeath[1:max.T]
  birthY.sim[,jj] <- myList$Ybirth[1:max.T]
  deathY.sim[,jj] <- myList$Ydeath[1:max.T]
  
  sim.X.all[,jj] <- c(0, cumsum(birthX.sim[,jj] - deathX.sim[,jj]))
  sim.Y.all[,jj] <- c(0, cumsum(birthY.sim[,jj] - deathY.sim[,jj]))
}

Y.all <- sim.Y.all
X.all <- sim.X.all

pri.A.X <- c(0.001, 0.001); # non-informative prior for A.X
pri.alpha.X <- c(0.001, 0.001); # inormative prior for alpha.X
pri.beta.X <- c(0.001, 0.001); # inormative prior for beta.X
pri.KM <- c(0.001, 0.001); # non-informative prior for KM
pri.A.Y <- c(0.001, 0.001); # non-informative prior for A.X
pri.alpha.Y <- c(0.001, 0.001); # inormative prior for alpha.X
pri.beta.Y <- c(0.001, 0.001); # inormative prior for beta.X
pri.B <- c(0.001, 0.001); # non-informative prior for KM

tun.KM <- 1; 
tun.Delta.X <- c(1.0, 1);
tun.Delta.Y <- c(1.0, 1);

effrepeat <- 200
burn <- 0; thin <- 1;
nrepeat <- burn + thin*effrepeat;

selrow <- seq(from = burn + thin, by = thin, length.out = effrepeat)

#initial value setting 
theta.X <- c(A.X, B.X)
theta.Y <- c(A.Y, B.Y, K.M)

Delta.X <- c(alpha.X, beta.X) #initial & true values of delay parameter of X 
Delta.Y <- c(alpha.Y, beta.Y) #initial & true values of delay parameter of Y 

RR.all = array(0, dim = c(max.T, 4, nsample)) #saving number of reaction 

# for(jj in 1:nsample){
#   Y.diff <- diff(Y.all[,jj]) #y(i+1) - y(i)
#   X.diff <- diff(X.all[,jj]) #x(i+1) - x(i)    
#   for (i in 1:max.T) {
#     RR.all[i,1,jj] <- max(Y.diff[i],0)  # # of birth reaction of Y
#     RR.all[i,2,jj] <- max(-Y.diff[i],0) # # of death reaction of Y
#     RR.all[i,3,jj] <- max(X.diff[i],0)  # # of birth reaction of X
#     RR.all[i,4,jj] <- max(-X.diff[i],0) # # of death reaction of X
#   }
# }

################################################
# iteration start!!!
################################################

# matrix & vector for saving MCMC results
count_KM <- 0; count_X <- rep(0, nsample); count_Delta.X <- 0; count_Delta.Y <- 0;  

theta <- matrix(0,nrow = nrepeat, ncol=8)
X.fit <- array(0, dim = c(nrepeat, max.T+1, nsample))
R.fit <- array(0, dim = c(4*nrepeat, max.T, nsample))

#initial scales of KM, delata.X, A.X used in RAM method 
KM.S <- 10
Delta.X.S <- diag(2)
Delta.Y.S <- diag(2)

K.i <- KI(Delta.X, maxt = max.T); 

ptnum <- 4;
useall <- TRUE;
theta[1,] = c(theta.X[1], theta.Y[3], Delta.X[1], Delta.X[2], theta.Y[1], Delta.Y[1], Delta.Y[2], theta.X[2])


myList <- TimeDelayGillespieforXY(A.X = A.X, B.X = B.X, alpha.X = alpha.X, beta.X = beta.X, A.Y = A.Y, B.Y = B.Y, alpha.Y = alpha.Y, beta.Y = beta.Y, K.M = K.M, repnum = max.T*10000, maxT = max.T+3, Volume = vol.est)
for(jj in 1:nsample){
  RR.all[,1,jj] <- myList$Ybirth[1:max.T]
  RR.all[,2,jj] <- myList$Ydeath[1:max.T]
  RR.all[,3,jj] <- myList$Xbirth[1:max.T]
  RR.all[,4,jj] <- myList$Xdeath[1:max.T]
  
  X.all[,jj] <- c(0, cumsum(RR.all[,3,jj] - RR.all[,4,jj]))
  
}

# RR.all[,1,] <- birthY.sim
# RR.all[,2,] <- deathY.sim
# RR.all[,3,] <- birthX.sim
# RR.all[,4,] <- deathX.sim

# fix KM
# theta[,2] <- scaling.rate * theta.Y[3];

# fix B
theta[,8] <- theta.X[2];
# fix AY
# theta[,5] <- theta.Y[1];
# fix alphaY and beta Y
# theta[,6] <- Delta.Y[1];
# theta[,7] <- Delta.Y[2];
theta[,6] <- 9;
theta[,7] <- 1;


X.all <- X.all/2
RR.all <- RR.all/2

sim.X.all <- sim.X.all * scaling.rate
sim.Y.all <- sim.Y.all * scaling.rate
X.all <- X.all * scaling.rate
Y.all <- Y.all * scaling.rate
RR.all <- RR.all * scaling.rate

for(rep in 2:nrepeat){
  # step 1 & 2: sampling  r2 and r1 (death and birth of Y)
  K.i <- KI(P = theta[rep-1,3:4], maxt = max.T);
  
  for(jj in 1:nsample){
    RR.all[,1:2,jj] <- impute_r.Y(Y.all[,jj], B.Y = theta[rep-1,8])
  }
  
  # step 3: sampling X &   r3, r4
  # updating X using independent chain MH
  
  # generate a proposal mean trajectory using the current parameter set.
  for(jj in 1:nsample){
    myListX <- TimeDelayGillespieforXR(A.X = theta[rep-1,1], B.X = theta[rep-1,8], alpha.X = theta[rep-1,3], beta.X = theta[rep-1,4], repnum = round(max.T*10000), maxT = max.T+5, Volume = vol.est)
    # X.bir.st <- floor(myListX$Xbirth[1:max.T])
    # X.dea.st <- floor(myListX$Xdeath[1:max.T])
    X.bir.st <- myListX$Xbirth[1:max.T]
    X.dea.st <- myListX$Xdeath[1:max.T]
    X.star <- c(0, cumsum(X.bir.st - X.dea.st));
    # print(X.update$errflg)
    if (useall == TRUE){
      fy.st = A.Y * KI.Y(P = theta[rep-1,6:7],in.X = X.star, K.M=theta[rep-1,2])
      fy    = A.Y * KI.Y(P = theta[rep-1,6:7],in.X = X.all[,jj]     , K.M=theta[rep-1,2])
    }else{
      fy.st = A.Y * KI.Ynt(P = theta[rep-1,6:7],in.X = X.star, N = ptnum, K.M=theta[rep-1,2])
      fy    = A.Y * KI.Ynt(P = theta[rep-1,6:7],in.X = X.all[,jj]     , N = ptnum, K.M=theta[rep-1,2])
    }
    # l.lik.st <- sum(r.X.birth * log(KI.star), na.rm = T) - Ax * (sum(KI.star))
    
    # lik.X.bir <- sum(log(dpois(RR.all[,3,jj], theta[rep-1,1] * K.i)+1e-300), na.rm = T)
    # lik.X.bir.st <- sum(log(dpois(X.bir.st, theta[rep-1,1] * K.i)+1e-300), na.rm = T)
    
    
    
    lik.X.bir <- sum(log(dpois_conti(vol.gen * RR.all[,3,jj], vol.gen * theta[rep-1,1] * K.i)+1e-300), na.rm = T)
    lik.X.bir.st <- sum(log(dpois_conti(vol.gen * X.bir.st, vol.gen * theta[rep-1,1] * K.i)+1e-300), na.rm = T)
    
    # lik.X.dea <- sum(log(dpois(RR.all[,4,jj], theta[rep-1, 8] * (X.all[-1,jj] + X.all[-(max.T+1),jj])/2)+1e-300), na.rm = T)
    # lik.X.dea.st <- sum(log(dpois(X.dea.st, theta[rep-1, 8] * (X.star[-1] + X.star[-(max.T+1)])/2)+1e-300), na.rm = T)
    
    
    lik.X.dea <- sum(log(dpois_conti(vol.gen * RR.all[,4,jj], vol.gen * theta[rep-1, 8] * (X.all[-1,jj] + X.all[-(max.T+1),jj])/2) + 1e-300), na.rm = T)
    lik.X.dea.st <- sum(log(dpois_conti(vol.gen * X.dea.st, vol.gen * theta[rep-1, 8] * (X.star[-1] + X.star[-(max.T+1)])/2) + 1e-300), na.rm = T)
    
    lik.Y.bir.st = sum(log(dpois_conti(vol.gen * RR.all[,1,jj], vol.gen * fy.st[,1])+1e-300), na.rm = T)
    lik.Y.bir    = sum(log(dpois_conti(vol.gen * RR.all[,1,jj], vol.gen * fy[,1]   )+1e-300), na.rm = T)
    
    prop.X.bir <- sum(log(dpois_conti(vol.est * RR.all[,3,jj], vol.est * theta[rep-1,1] * K.i)+1e-300), na.rm = T)
    prop.X.bir.st <- sum(log(dpois_conti(vol.est * X.bir.st, vol.est * theta[rep-1,1] * K.i)+1e-300), na.rm = T)

    prop.X.dea <- sum(log(dpois_conti(vol.est * RR.all[,4,jj], theta[rep-1, 8] * vol.est * (X.all[-1,jj] + X.all[-(max.T+1),jj])/2)+1e-300), na.rm = T)
    prop.X.dea.st <- sum(log(dpois_conti(vol.est * X.dea.st, theta[rep-1, 8] * vol.est * (X.star[-1] + X.star[-(max.T+1)])/2)+1e-300), na.rm = T)
    
    
    # lik.Y.bir.st = sum(log(dpois(RR.all[,1,jj],fy.st[,1])+1e-300), na.rm = T)
    # lik.Y.bir    = sum(log(dpois(RR.all[,1,jj],fy[,1]   )+1e-300), na.rm = T)
    
    
    
    # prior.X.st = sum(log(dgamma(X.star , shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
    # prior.X   = sum(log(dgamma(X.all[,jj], shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
    
    # logMH <- q.Y.st - q.Y + prior.X.st - prior.X; # considering prior.
    # logMH <- q.Y.st - q.Y; # Completely non-informative, i.e., always prior.X.st == prior.X
    
    logMH <- (lik.Y.bir.st - lik.Y.bir +
                lik.X.bir.st - lik.X.bir +
                lik.X.dea.st - lik.X.dea +
                prop.X.bir - prop.X.bir.st +
                prop.X.dea - prop.X.dea.st) # Completely non-informative & considering tuning of proposal distribution (cell volume scaling).
    
    # print(logMH);
    if(!is.nan(logMH) && runif(1)<exp(logMH)){
      X.all[,jj] <- X.star; RR.all[,3,jj] <- X.bir.st; RR.all[,4,jj] <- X.dea.st;
      count_X[jj] = count_X[jj] + 1;
    }
  }
  
  # step  4: samping A.X 
  g_11 <- sum(K.i);
  theta[rep,1] = rgamma(1,shape = sum(RR.all[,3,]) + nsample * pri.A.X[1], rate = nsample * (g_11 + pri.A.X[2]));
  
  
  # step 5 & 6: sampling alpha.X and beta.X: the delay parameters for the birth reaction of X.
  p.update <- MH.P.X.all(P = theta[rep-1,3:4], Delta.X.S, rep, RR.all[,3,], Ax = theta[rep,1],  tun = tun.Delta.X, pri.alpha.X = pri.alpha.X, pri.beta.X = pri.beta.X, maxt = max.T)
  theta[rep,3:4] = p.update$P
  Delta.X.S = p.update$S
  count_Delta.X = count_Delta.X + p.update$count
  
  
  # step 7: sampling the Michaelis-Menten constant K.M
  KM.update = MH.KM.all(theta[rep-1,2] , KM.S, rep, RR.all[,1,], X.all, b = tun.KM, pri.KM = pri.KM, Delta.Y = Delta.Y, flatpri = TRUE)
  theta[rep,2] = KM.update$km;
  KM.S = KM.update$s
  count_KM = count_KM + KM.update$count
  
  # step 8: sampling alpha.Y and beta.Y: the delay parameters for the birth reaction of Y.
  # p.update <- MH.P.Y.all(P = theta[rep-1,6:7], S = Delta.Y.S, rep = rep, r.Y.birth = RR.all[,1,], in.X.all = X.all, Ay = theta[rep,5], K.M = theta[rep, 2], tun = tun.Delta.Y, pri.alpha.Y = pri.alpha.Y, pri.beta.Y = pri.beta.Y, maxt = max.T)
  # theta[rep,6:7] = p.update$P
  # Delta.Y.S = p.update$S
  # count_Delta.Y = count_Delta.Y + p.update$count
  
  # step 9: sampling A.Y
  KY.sum <- 0
  for(ii in 1:nsample){
    KY.i <- KI.Y(P = theta[rep,6:7], in.X = X.all[,ii], K.M = theta[rep,2])
    KY.sum <- KY.sum + sum(KY.i[,1])
  }
  theta[rep,5] = rgamma(1,shape = sum(RR.all[,1,]) + nsample * pri.A.Y[1], rate = KY.sum + nsample*pri.A.Y[2]);

  
  # step 10: sampling B(=B.X=B.Y); the common dilution rate 
  # theta[rep,8] = rgamma(1,shape = sum(RR.all[,2,]) + sum(RR.all[,4,]) + nsample*pri.B[1],
  #                       rate = sum(X.all) + sum(Y.all) - 0.5*sum(X.all[max.T+1,]) - 0.5*sum(Y.all[max.T+1,])+ nsample*pri.B[2]);
  
  
  X.fit[rep,,] = X.all
  R.fit[4*rep-3,,] = RR.all[,1,] # birth number of Y
  R.fit[4*rep-2,,] = RR.all[,2,] # death number of Y
  R.fit[4*rep-1,,] = RR.all[,3,] # birth number of X
  R.fit[4*rep-0,,] = RR.all[,4,] # death number of X
  
  
  if(rep%%3 ==0){
    cat(rep)
    cat(" ")
  }
  if(theta[rep,1] > 300){
    print("Estimated Ax > 300")
    break
  } 
}

# the estimated reaction numbers from MCMC algorithm.
birthY<- R.fit[4*(selrow-1) + 1,,]  
deathY<- R.fit[4*(selrow-1) + 2,,]
birthX<- R.fit[4*(selrow-1) + 3,,]
deathX<- R.fit[4*(selrow-1) + 4,,]

print(paste0("Acceptance ratio for X: ", count_X / nrepeat))  # Acceptance ratio for X
print(paste0("Acceptance ratio for alpha.X and beta.X: ", count_Delta.X / nrepeat))  # Acceptance ratio for alpha.X and beta.X
print(paste0("Acceptance ratio for K.M: ", count_KM / nrepeat))  # Acceptance ratio for alpha.X and beta.X
print(paste0("Acceptance ratio for alpha.Y and beta.Y: ", count_Delta.Y / nrepeat))  # Acceptance ratio for alpha.X and beta.X

# generate a Y trajectory from the mean of the estimated parameters from the effective iteration indexed by 'selrow'.

# gen.num <- 10
summary_row = selrow[seq(from = 100, to = length(selrow), by = 5)[-1]]
gen.y2 <- matrix(0, nrow = length(summary_row), ncol = max.T+1)
tmp_idx <- 0
for(jj in summary_row){
  tmp_idx = tmp_idx + 1
  myList2 <- TimeDelayGillespieforXY(A.X = theta[jj,1], B.X = theta[jj,8], alpha.X = theta[jj,3], beta.X = theta[jj,4], 
                                     A.Y = theta[jj,5], B.Y = theta[jj,8], alpha.Y = theta[jj,6], beta.Y = theta[jj,7], 
                                     K.M = theta[jj,2], repnum = max.T*500, maxT = max.T+3, Volume = vol.est)
  gen.y2[tmp_idx,] <- c(0, cumsum(myList2$Ybirth[1:max.T] - myList2$Ydeath[1:max.T]))
}
mean.y <- colMeans(gen.y2)
plot(rowMeans(sim.Y.all), ylim = c(0,1000))
lines(mean.y)

# Y(\infty) fitting check code

# plot(theta[,1]*theta[,5] / (theta[,8] *(theta[,1] + theta[,8] * theta[,2])))


jj <- 16

myList2 <- TimeDelayGillespieforXY(A.X = theta[jj,1], B.X = theta[jj,8], alpha.X = theta[jj,3], beta.X = theta[jj,4], 
                                   A.Y = theta[jj,5], B.Y = theta[jj,8], alpha.Y = theta[jj,6], beta.Y = theta[jj,7], 
                                   K.M = theta[jj,2], repnum = max.T*500, maxT = max.T+3, Volume = vol.est)
plot(c(0, cumsum(myList2$Ybirth[1:max.T] - myList2$Ydeath[1:max.T])))

plot(rowMeans(sim.Y.all))
lines(c(0, cumsum(myList2$Ybirth[1:max.T] - myList2$Ydeath[1:max.T])))

max.T <- 100

gen.y3 <- matrix(0, nrow = 10, ncol = max.T+1)
gen.y4 <- matrix(0, nrow = 10, ncol = max.T+1)
for(jj in 1:10){
  myList3 <- TimeDelayGillespieforXY(A.X = 10, B.X = 0.05, alpha.X = 3.6, beta.X = 0.6, 
                                     A.Y = 60, B.Y = 0.05, alpha.Y = 3.6, beta.Y = 0.6, 
                                     K.M = 200, repnum = max.T*500, maxT = max.T+3, Volume = 1)
  myList4 <- TimeDelayGillespieforXY(A.X = 10, B.X = 0.05, alpha.X = 3.6, beta.X = 0.6,
                                     A.Y = 60, B.Y = 0.05, alpha.Y = 3.6, beta.Y = 0.6,
                                     K.M = 200, repnum = max.T*500, maxT = max.T+3, Volume = 5)
  gen.y3[jj,] <- c(0, cumsum(myList3$Xbirth[1:max.T] - myList3$Xdeath[1:max.T]))
  gen.y4[jj,] <- c(0, cumsum(myList4$Xbirth[1:max.T] - myList4$Xdeath[1:max.T]))
}

matplot(t(gen.y3), type ="l", ylim = c(0,260), xlab = "Time", ylab = "X(t), Volume = 1")
matplot(t(gen.y4), type ="l", ylim = c(0,260), xlab = "Time", ylab = "X(t), Volume = 5")


plot()

plot(rowMeans(sim.Y.all))
lines(c(0, cumsum(myList2$Ybirth[1:max.T] - myList2$Ydeath[1:max.T])))
