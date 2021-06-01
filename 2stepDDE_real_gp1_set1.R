setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode/")
source('2stepDDE_functions.R')

rndseed <- round((as.numeric(Sys.time())*1000)%% 10000)
set.seed(rndseed)


effrepeat <- 1000 # number of MCMC iterations

int <- 1; 

# read data
raw.data = read.csv("2mM_IPTG.csv",header = F)
raw.data = read.csv("0_2mM_IPTG.csv",header = F)
raw.data = read.csv("0_05mM_IPTG.csv",header = F)

# rawdata <- round(all.data[1:(max.T+1),1:data.num]);     
diffdata <- matrix(0, nrow = nrow(raw.data), ncol = ncol(raw.data)); 
for(i in 1:ncol(diffdata)){
  diffdata[,i]=raw.data[,i] - raw.data[1,i]
}

for (i in 1:nrow(raw.data)) {
  for (j in 1:ncol(raw.data)) {
    if(diffdata[i,j]<0) diffdata[i,j] =0
  }
}

nmlzdata1 = round(diffdata/0.089)[1:27,1:49]; # data from experiment 1
nmlzdata2 = round(diffdata/0.089)[1:32,50:88]; # data from experiment 2
nmlzdata1 = round(diffdata/0.089)[,1:33]; # data from experiment 1
nmlzdata2 = round(diffdata/0.089)[,34:73]; # data from experiment 2


rownum_without_outlier1 <- which(nmlzdata1[47,] > 200 & nmlzdata1[27,] < 600)
rownum_without_outlier2 <- which(nmlzdata2[32,] > 200 & nmlzdata2[32,] < 500)

# nmlzdata1[, rownum_without_outlier1]
# nmlzdata2[, rownum_without_outlier2]

Y.all <- nmlzdata1[, rownum_without_outlier1]
matplot(0:26, Y.all, type = "l", xlab = "Time (min)", ylab = "Molecular count (Y(t))")
Y.all <- nmlzdata2[, rownum_without_outlier2]
matplot(0:31, Y.all, type = "l", xlab = "Time (min)", ylab = "Molecular count (Y(t))")


X.all <- Y.all # X is hidden, but we set it as Y. It is almost meaningless.

max.T <- nrow(Y.all) - 1
nsample <- ncol(Y.all)

tun.B <- c(50,50, 100, 100);

# known (fixed) parameters
alpha.Y = 5.89; beta.Y = 0.89*int;
B.X <- 0.022;  B.Y <- 0.022; A.Y <- 35.4;

# unknowns parameters
A.X = 35.4*int; alpha.X = 5.89; beta.X = 0.89*int;
K.M = 100; 

tspan <- 0:max.T

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

burn <- 0; thin <- 1;
nrepeat <- burn + thin*effrepeat;

selrow <- seq(from = burn + thin, by = thin, length.out = effrepeat)

#initial value setting 
theta.X <- c(A.X, B.X)
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


birthX.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true birth number of X
deathX.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true death number of X
birthY.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true birth number of Y
deathY.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true death number of Y
sim.X.all <- matrix(0, nrow = max.T+1, ncol = nsample)
sim.Y.all <- matrix(0, nrow = max.T+1, ncol = nsample)

for(jj in 1:nsample){
  # myList is raw simulated data. 
  myList <- TimeDelayGillespieforXY(A.X = A.X, B.X = B.X, alpha.X = alpha.X, beta.X = beta.X, A.Y = A.Y, B.Y = B.Y, alpha.Y = alpha.Y, beta.Y = beta.Y, K.M = K.M, repnum = max.T*10000, maxT = max.T+3)
  # sim.X is true X data, and sim.Y is true Y data.
  birthX.sim[,jj] <- myList$Xbirth[1:max.T]
  deathX.sim[,jj] <- myList$Xdeath[1:max.T]
  birthY.sim[,jj] <- myList$Ybirth[1:max.T]
  deathY.sim[,jj] <- myList$Ydeath[1:max.T]
  
  sim.X.all[,jj] <- c(0, cumsum(birthX.sim[,jj] - deathX.sim[,jj]))
  sim.Y.all[,jj] <- c(0, cumsum(birthY.sim[,jj] - deathY.sim[,jj]))
}


RR.all[,1,] <- birthY.sim
RR.all[,2,] <- deathY.sim
RR.all[,3,] <- birthX.sim
RR.all[,4,] <- deathX.sim

for(rep in 2:nrepeat){
  # step 1 & 2: sampling  r2 and r1 (death and birth of Y)

  for(jj in 1:nsample){
    RR.all[,1:2,jj] <- impute_r.Y(Y.all[,jj], B.Y = theta[rep-1,8])
  }
  
  # step 3: sampling X & r3, r4
  # updating X using independent chain MH
  
  # generate a proposal mean trajectory using the current parameter set.
  for(jj in 1:nsample){
    myListX <- TimeDelayGillespieforXR(A.X = theta[rep-1,1], B.X = theta[rep-1,8], alpha.X = theta[rep-1,3], beta.X = theta[rep-1,4], repnum = round(max.T*10000), maxT = max.T+5)
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
    
    q.Y.st = sum(log(dpois(RR.all[,1,jj],fy.st[,1])+1e-300), na.rm = T)
    q.Y    = sum(log(dpois(RR.all[,1,jj],fy[,1]   )+1e-300), na.rm = T)
    
    # prior.X.st = sum(log(dgamma(X.star , shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
    # prior.X   = sum(log(dgamma(X.all[,jj], shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
    
    # logMH <- q.Y.st - q.Y + prior.X.st - prior.X; # considering prior.
    logMH <- q.Y.st - q.Y; # Completely non-informative, i.e., always prior.X.st == prior.X 
    
    # print(logMH);
    if(!is.nan(logMH) && runif(1)<exp(logMH)){
      X.all[,jj] <- X.star
      RR.all[,3,jj] <- X.bir.st
      RR.all[,4,jj] <- X.dea.st
      count_X[jj] = count_X[jj] + 1
    }
  }
  
  # step  4: samping A.X 
  K.i <- KI(P = theta[rep-1,3:4], maxt = max.T);
  g_11 <- sum(K.i);
  theta[rep,1] = rgamma(1,shape = sum(RR.all[,3,]) + nsample * pri.A.X[1], rate = nsample * (g_11 + pri.A.X[2]));
  # theta[rep,1] = theta[rep-1,1]
  
  
  # step 5 & 6: sampling alpha.X and beta.X: the delay parameters for the birth reaction of X.
  p.update <- MH.P.X.all(P = theta[rep-1,3:4], Delta.X.S, rep, RR.all[,3,], Ax = theta[rep,1],  tun = tun.Delta.X, pri.alpha.X = pri.alpha.X, pri.beta.X = pri.beta.X, maxt = max.T)
  theta[rep,3:4] = p.update$P
  Delta.X.S = p.update$S
  count_Delta.X = count_Delta.X + p.update$count
  
  
  # step 7: sampling the Michaelis-Menten constant K.M
  KM.update = MH.KM.all(theta[rep-1,2] , KM.S, rep, RR.all[,1,], X.all, b = tun.KM, pri.KM = pri.KM, Delta.Y = c(theta[rep-1,6], theta[rep-1,7]))
  theta[rep,2] = KM.update$km;
  KM.S = KM.update$s
  count_KM = count_KM + KM.update$count
  # theta[rep,2] = theta[rep-1,2]
  
  # step 8: sampling alpha.Y and beta.Y: the delay parameters for the birth reaction of Y.
  # p.update <- MH.P.Y.all(P = theta[rep-1,6:7], S = Delta.Y.S, rep = rep, r.Y.birth = RR.all[,1,], in.X.all = X.all, 
  #                        Ay = theta[rep-1,5], K.M = theta[rep, 2], tun = tun.Delta.Y, pri.alpha.Y = pri.alpha.Y, pri.beta.Y = pri.beta.Y, maxt = max.T)
  # theta[rep,6:7] = p.update$P
  # Delta.Y.S = p.update$S
  # count_Delta.Y = count_Delta.Y + p.update$count
  theta[rep,6:7] = theta[rep-1,6:7]
  
  # 
  # step 9: sampling A.Y
  # KY.sum <- 0
  # for(ii in 1:nsample){
  #   KY.i <- KI.Y(P = theta[rep,6:7], in.X = X.all[,ii], K.M = theta[rep,2])
  #   KY.sum <- KY.sum + sum(KY.i[,1])
  # }
  # theta[rep,5] = rgamma(1,shape = sum(RR.all[,1,]) + nsample * pri.A.Y[1], rate = KY.sum + nsample*pri.A.Y[2]);
  theta[rep,5] = theta[rep-1,5]
  
  
  # step 10: sampling B(=B.X=B.Y); the common dilution rate 
  # theta[rep,8] = rgamma(1,shape = sum(RR.all[,2,]) + sum(RR.all[,4,]) + nsample*pri.B[1], 
  #                       rate = sum(X.all) + sum(Y.all) - 0.5*sum(X.all[max.T+1,]) - 0.5*sum(Y.all[max.T+1,])+ nsample*pri.B[2]);
  
  theta[rep,8] = theta[rep-1,8]
  

  X.fit[rep,,] = X.all
  R.fit[4*rep-3,,] = RR.all[,1,] # birth number of Y
  R.fit[4*rep-2,,] = RR.all[,2,] # death number of Y
  R.fit[4*rep-1,,] = RR.all[,3,] # birth number of X
  R.fit[4*rep-0,,] = RR.all[,4,] # death number of X
  
  
  if(rep%%100 ==0){
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

currentT <- Sys.time()
timestamp0 <- paste(substr(currentT, 1,4), substr(currentT, 6,7), substr(currentT, 9,10),substr(currentT, 12,13), substr(currentT, 15,16), substr(currentT, 18,19), sep = "")

compath <- "/home/users/hphong/TimeDelay/real_group1//";

par.filename <- paste("Parameters_real_gp1_set1_maxT", toString(max.T), "_", toString(timestamp0), ".csv", sep = "")

result.par<- matrix(0, nrow = max(max.T, nrepeat, 4*nsample+10), ncol = 4*nsample+11)
result.par[1:nrepeat, 1:8] <- theta
result.par[1:max.T, 9:(4*nsample+8)] <- cbind(birthX.sim, deathX.sim, birthY.sim, deathY.sim)
result.par[1:(nsample+4),4*nsample+9] <- c(count_X, count_Delta.X, count_Delta.Y, count_KM, rndseed)
result.par[1:7, (4*nsample+10):(4*nsample+11)] <- rbind(pri.A.X, pri.alpha.X, pri.beta.X, pri.A.Y, pri.alpha.Y, pri.beta.Y, pri.B)

write.table(result.par, paste(compath, par.filename, sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")

for(jj in 1:nsample){
  Reaction.filename <- paste("ReactionNumbers_real_gp1_set1_maxT", toString(max.T), "_",  toString(timestamp0), "_", toString(jj), ".csv", sep = "")
  result.Reactions <- cbind(birthX[,,jj], deathX[,,jj], birthY[,,jj], deathY[,,jj])
  write.table(result.Reactions, paste(compath, Reaction.filename, sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
}


