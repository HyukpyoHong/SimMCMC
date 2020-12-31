setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode")
# setwd("C:/Users/HongHyukpyo/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/ProfChoi_Code_20200326")
# setwd("C:/Users/Hyukpyo Hong/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/ProfChoi_Code_20200326")
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

max.T <- 100 # simulated data will be given from t = 0, ..., max.T
tspan <- 0:max.T
nsample <- 20;
vol <- 3;

birthX.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true birth number of X
deathX.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true death number of X
birthY.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true birth number of Y
deathY.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true death number of Y
sim.X.all <- matrix(0, nrow = max.T+1, ncol = nsample)
sim.Y.all <- matrix(0, nrow = max.T+1, ncol = nsample)

for(jj in 1:nsample){
  # myList is raw simulated data. 
  myList <- TimeDelayGillespieforXY(A.X = A.X, B.X = B.X, alpha.X = alpha.X, beta.X = beta.X, A.Y = A.Y, B.Y = B.Y, alpha.Y = alpha.Y, beta.Y = beta.Y, K.M = K.M, repnum = max.T*10000, maxT = max.T+3, Volume = vol)
  # sim.X is true X data, and sim.Y is true Y data.
  birthX.sim[,jj] <- floor(myList$Xbirth[1:max.T])
  deathX.sim[,jj] <- floor(myList$Xdeath[1:max.T])
  birthY.sim[,jj] <- floor(myList$Ybirth[1:max.T])
  deathY.sim[,jj] <- floor(myList$Ydeath[1:max.T])
  
  sim.X.all[,jj] <- c(0, cumsum(birthX.sim[,jj] - deathX.sim[,jj]))
  sim.Y.all[,jj] <- c(0, cumsum(birthY.sim[,jj] - deathY.sim[,jj]))
}
Y.all <- sim.Y.all
X.all <- sim.X.all

pri.A.X <- c(0.001, 0.001); # non-informative prior for A.X
pri.alpha.X <- c(0.001, 0.001); # inormative prior for alpha.X
pri.beta.X <- c(0.001, 0.001); # inormative prior for beta.X
pri.KM <- c(0.001, 0.001); # non-informative prior for KM

tun.KM <- 1; 
tun.Delta.X <- c(1.0, 1);

effrepeat <- 100
burn <- 0; thin <- 1;
nrepeat <- burn + thin*effrepeat;

selrow <- seq(from = burn + thin, by = thin, length.out = effrepeat)

#initial value setting 
theta.X <- c(A.X , B.X)
theta.Y <- c(A.Y, B.Y, K.M)

Delta.X <- c(alpha.X, beta.X) #initial & true values of delay parameter of X 
Delta.Y <- c(alpha.Y, beta.Y) #initial & true values of delay parameter of Y 

RR.all = array(0, dim = c(max.T, 4, nsample)) #saving number of reaction 

for(jj in 1:nsample){
  Y.diff <- diff(Y.all[,jj]) #y(i+1) - y(i)
  X.diff <- diff(X.all[,jj]) #x(i+1) - x(i)    
  
  for (i in 1:max.T) {
    RR.all[i,1,jj] <- max(Y.diff[i],0)  # # of birth reaction of Y
    RR.all[i,2,jj] <- max(-Y.diff[i],0) # # of death reaction of Y
    RR.all[i,3,jj] <- max(X.diff[i],0)  # # of birth reaction of X
    RR.all[i,4,jj] <- max(-X.diff[i],0) # # of death reaction of X
  }
}

################################################
# iteration start!!!
################################################

# matrix & vector for saving MCMC results
count_KM <- 0; count_X <- rep(0, nsample); count_Delta.X <- 0;

theta <- matrix(0,nrow = nrepeat, ncol=4)
X.fit <- array(0, dim = c(nrepeat, max.T+1, nsample))
R.fit <- array(0, dim = c(4*nrepeat, max.T, nsample))

#initial scales of KM, delata.X, A.X used in RAM method 
KM.S <- 10
Delta.X.S <- diag(2)

K.i <- KI(Delta.X, maxt = max.T); 

KM.lik.fit <- rep(0, nrepeat)
KM.star.lik.fit <- rep(0, nrepeat)

ptnum <- 4;
useall <- TRUE;
theta[1,] = c(theta.X[1], theta.Y[3], Delta.X[1], Delta.X[2])
RR.all[,3,] <- birthX.sim

for(rep in 2:nrepeat) {
  # step 1 & 2: sampling  r2 and r1 (death and birth of Y)
  K.i <- KI(P = theta[rep-1,3:4], maxt = max.T);
  
  for(jj in 1:nsample){
    RR.all[,1:2,jj] <- impute_r.Y(Y.all[,jj], B.Y = B.Y)
  }
  
  # step 3: sampling X &   r3, r4 
  # updating X using independent chain MH

  # generate a proposal mean trajectory using the current parameter set.
  for(jj in 1:nsample){
    myListX <- TimeDelayGillespieforXR(A.X = theta[rep-1,1], B.X = B.X, alpha.X = theta[rep-1,3], beta.X = theta[rep-1,4], repnum = round(max.T*10000), maxT = max.T+5)
    X.bir.st <- myListX$Xbirth[1:max.T]
    X.dea.st <- myListX$Xdeath[1:max.T]
    X.star <- c(0, cumsum(X.bir.st - X.dea.st));
    # print(X.update$errflg)
    if (useall == TRUE){
      fy.st = A.Y * KI.Y(P = Delta.Y,in.X = X.star, K.M=theta[rep-1,2])
      fy    = A.Y * KI.Y(P = Delta.Y,in.X = X.all[,jj]     , K.M=theta[rep-1,2])
    }else{
      fy.st = A.Y * KI.Ynt(P = Delta.Y,in.X = X.star, N = ptnum, K.M=theta[rep-1,2])
      fy    = A.Y * KI.Ynt(P = Delta.Y,in.X = X.all[,jj]     , N = ptnum, K.M=theta[rep-1,2])
    }
    q.Y.st = sum(log(dpois(RR.all[,1,jj],fy.st[,1])+1e-300), na.rm = T)
    q.Y    = sum(log(dpois(RR.all[,1,jj],fy[,1]   )+1e-300), na.rm = T)
    
    # prior.X.st = sum(log(dgamma(X.star , shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
    # prior.X   = sum(log(dgamma(X.all[,jj], shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
    
    # logMH <- q.Y.st - q.Y + prior.X.st - prior.X; # considering prior.
    logMH <- q.Y.st - q.Y; # Completely non-informative, i.e., always prior.X.st == prior.X 
    
    # print(logMH);
    if(!is.nan(logMH) && runif(1)<exp(logMH)){
      X.all[,jj] <- X.star; RR.all[,3,jj] <- X.bir.st; RR.all[,4,jj] <- X.dea.st;
      count_X[jj] = count_X[jj] + 1;
    }
  }
  
  # step  4: samping A.X 
  g_11 <- sum(K.i);
  theta[rep,1] = rgamma(1,shape = sum(RR.all[,3,]) + nsample * pri.A.X[1], rate = nsample * (g_11 + pri.A.X[2]));
  
  # theta[rep,1] = rgamma(1,shape = sum(RR.all[,3,]), rate = nsample * g_11); # Completely non-informative, i.e., always prior.X.st == prior.X 
  
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
  
  X.fit[rep,,] = X.all
  R.fit[4*rep-3,,] = RR.all[,1,] # birth number of Y
  R.fit[4*rep-2,,] = RR.all[,2,] # death number of Y
  R.fit[4*rep-1,,] = RR.all[,3,] # birth number of X
  R.fit[4*rep-0,,] = RR.all[,4,] # death number of X
  
  
  if(rep%%10 ==0 ){
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

# generate a Y trajectory from the mean of the estimated parameters from the effective iteration indexed by 'selrow'.

gen.num <- 10
gen.y2 <- matrix(0, nrow = gen.num, ncol = max.T+1)
for(jj in 1:gen.num){
  myList2 <- TimeDelayGillespieforXY(A.X = colMeans(theta[selrow,])[1], B.X = 0.05, alpha.X = colMeans(theta[selrow,])[3], beta.X = colMeans(theta[selrow,])[4], A.Y = 60, B.Y = 0.05, alpha.Y = 3.6, beta.Y = 0.6, K.M = colMeans(theta[selrow,])[2], repnum = max.T*500, maxT = max.T+3)
  sim.X2 <- approx(myList2$TList[!is.na(myList2$TList)], myList2$XList[!is.na(myList2$XList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList2$XList[!is.na(myList2$XList)]))$y
  gen.y2[jj,] <- approx(myList2$TList[!is.na(myList2$TList)], myList2$YList[!is.na(myList2$YList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList2$YList[!is.na(myList2$YList)]))$y
}
mean.y <- colMeans(gen.y2)

matplot(0:max.T, X.all1, type = "l", ylim = c(0,250))
matplot(0:max.T, X.all3, type = "l", ylim = c(0,250))



matplot(0:max.T, Y.all, type = "l", ylim = c(0,600))
matplot(0:max.T, Y.all, type = "l", ylim = c(0,150), xlim = c(0,40))

matplot(0:max.T, Y.all1, type = "l", ylim = c(0,250))



