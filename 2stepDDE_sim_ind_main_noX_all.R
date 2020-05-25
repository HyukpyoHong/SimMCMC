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

max.T <- 25 # simulated data will be given from t = 0, ..., max.T
tspan <- 0:max.T
nsample <- 2;

birthX.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true birth number of X
deathX.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true death number of X
birthY.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true birth number of Y
deathY.sim <- matrix(0, nrow = max.T, ncol = nsample) # a list for the true death number of Y
sim.X.all <- matrix(0, nrow = max.T+1, ncol = nsample)
sim.Y.all <- matrix(0, nrow = max.T+1, ncol = nsample)

for(jj in 1:nsample){
  # myList is raw simulated data. 
  myList <- TimeDelayGillespieforXY(A.X = A.X, B.X = B.X, alpha.X = alpha.X, beta.X = beta.X, A.Y = A.Y, B.Y = B.Y, alpha.Y = alpha.Y, beta.Y = beta.Y, K.M = K.M, repnum = max.T*500, maxT = max.T+3)
  
  birthX.sim.tmp <- as.numeric(diff(c(0, myList$XList)) == 1) # binary for the birth reaction of X
  deathX.sim.tmp <- as.numeric(diff(c(0, myList$XList)) == -1) # binary for the death reaction of X
  birthY.sim.tmp <- as.numeric(diff(c(0, myList$YList)) == 1) # binary for the birth reaction of Y
  deathY.sim.tmp <- as.numeric(diff(c(0, myList$YList)) == -1) # binary for the death reaction of Y
  
  # sim.X is true X data, and sim.Y is true Y data.
  sim.X.all[,jj] <- approx(myList$TList[!is.na(myList$TList)], myList$XList[!is.na(myList$XList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList$XList[!is.na(myList$XList)]))$y
  sim.Y.all[,jj] <- approx(myList$TList[!is.na(myList$TList)], myList$YList[!is.na(myList$YList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList$YList[!is.na(myList$YList)]))$y
  
  # stacking the birth and death reaction numbers of X and Y.
  for(i in 1:length(ceiling(myList$TList))){
    t <- ceiling(myList$TList)[i]
    if(is.na(t)){
      warning("Current time t is NA.")
      warning("repnum in the gillespie algorithm might be not enough.")
      break;
    }
    if(t > max.T){
      break;
    }
    birthX.sim[t,jj] <- birthX.sim[t,jj] + birthX.sim.tmp[i]
    deathX.sim[t,jj] <- deathX.sim[t,jj] + deathX.sim.tmp[i]
    birthY.sim[t,jj] <- birthY.sim[t,jj] + birthY.sim.tmp[i]
    deathY.sim[t,jj] <- deathY.sim[t,jj] + deathY.sim.tmp[i]
  }
}
Y.all <- sim.Y.all
X.all <- sim.X.all

tun.B <- c(50,50, 100, 100);
tmp <- seq(from=0.1, by=1, length.out = max.T+1) # tunning parameter for the setting 1.
# tun.X <- 8 * tmp^2 / (600 + tmp^2) + 0.3
tun.X <- seq(from=0.1, by=0.03, length.out = max.T+1)
# plot(tun.X)

pri.A.X <- c(10 * 0.01, 0.01); # non-informative prior for A.X
pri.alpha.X <- c(3.6 * 0.1, 0.1); # inormative prior for alpha.X
pri.beta.X <- c(0.01, 0.01); # inormative prior for beta.X
pri.KM <- c(200* 0.01, 0.01); # non-informative prior for KM

tun.KM <- 1; 
tun.Delta.X <- c(1.0, 1);

effrepeat <- 500
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
count_KM <- 0; count_X <- 0; count_Delta.X <- 0;

theta <- matrix(0,nrow = nrepeat, ncol=4)
X.fit <- array(0, dim = c(nrepeat, max.T+1, nsample))
R.fit <- array(0, dim = c(4*nrepeat, max.T, nsample))

#initial scales of KM, delata.X, A.X used in RAM method 
KM.S <- 10
Delta.X.S <- diag(2)

# edited by Hyukpyo. X is fixed as the true X data.

K.i <- KI(Delta.X, maxt = max.T); 

AcceptR.fit <- matrix(0,nrow = nrepeat, ncol = 5)
AcceptR.star.fit <- matrix(0,nrow = nrepeat, ncol = 5)

KM.lik.fit <- rep(0, nrepeat)
KM.star.lik.fit <- rep(0, nrepeat)

ptnum <- 4;
useall <- TRUE;
theta[1,] = c(theta.X[1], theta.Y[3], Delta.X[1], Delta.X[2])

X.all.star <- X.all
X.bir.st <- matrix(0, nrow = max.T, ncol = nsample)
X.dea.st <- matrix(0, nrow = max.T, ncol = nsample)

q.Y.st <- 0; q.Y <- 0;

for(rep in 2:nrepeat) {
  # step 1 & 2: sampling  r2 and r1 (death and birth of Y)
  for(jj in 1:nsample){
    RR.all[,1:2,jj] <- impute_r.Y(Y.all[,jj], B.Y = B.Y)
  }
  # step 3: sampling X &   r3, r4 
  # updating X using independent chain MH
  
  # generate a proposal mean trajectory using the current parameter set.
  for(jj in 1:nsample){
    myListX <- TimeDelayGillespieforXR(A.X = theta[rep-1,1], B.X = B.X, alpha.X = theta[rep-1,3], beta.X = theta[rep-1,4], repnum = round(max.T*500), maxT = max.T+5)
    X.bir.st[,jj] <- myListX$Xbirth[1:max.T]
    X.dea.st[,jj] <- myListX$Xdeath[1:max.T]
    X.all.star[,jj] <- c(0, cumsum(X.bir.st[,jj] - X.dea.st[,jj]));
    # print(X.update$errflg)
    if (useall == TRUE){
      fy.st = A.Y * KI.Y(Delta.Y,in.X = X.all.star[,jj], K.M=theta[rep-1,2])
      fy    = A.Y * KI.Y(Delta.Y,in.X = X.all[,jj]     , K.M=theta[rep-1,2])
    }else{
      fy.st = A.Y * KI.Ynt(Delta.Y,in.X = X.all.star[,jj], N = ptnum, K.M=theta[rep-1,2])
      fy    = A.Y * KI.Ynt(Delta.Y,in.X = X.all[,jj]     , N = ptnum, K.M=theta[rep-1,2])
    }
    q.Y.st = q.Y.st + sum(log(dpois(RR.all[,1,jj],fy.st[,1])+1e-300), na.rm = T)
    q.Y    = q.Y + sum(log(dpois(RR.all[,1,jj],fy[,1]   )+1e-300), na.rm = T)
    
  }
  
  prior.X.st = sum(log(dgamma(X.all.star , shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
  prior.X   = sum(log(dgamma(X.all, shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
  
  logMH <- q.Y.st - q.Y + prior.X.st - prior.X;
  # print(logMH);
  if(!is.nan(logMH) && runif(1)<exp(logMH)){
    X.all=X.all.star; RR.all[,3,] <- X.bir.st; RR.all[,4,] <- X.dea.st;
    count_X = count_X + 1; 
  }
  
  # step  4: samping A.X 
  g_11 <- sum(K.i);
  theta[rep,1] = rgamma(1,shape = sum(RR.all[,3,]) + pri.A.X[1], rate = nsample * (g_11 + pri.A.X[2]));
  
  # step 5 & 6: sampling alpha.X and beta.X: the delay parameters for the birth reaction of X.
  p.update <- MH.P.X.all(P = theta[rep-1,3:4], Delta.X.S, rep, RR.all[,3,], Ax = theta[rep-1,1],  tun = tun.Delta.X, pri.alpha.X = pri.alpha.X, pri.beta.X = pri.beta.X, maxt = max.T)
  theta[rep,3:4] = p.update$P
  Delta.X.S = p.update$S
  count_Delta.X = count_Delta.X + p.update$count
  
  K.i <- KI(P = theta[rep-1,3:4], maxt = max.T);
  
  # step 7: sampling the Michaelis-Menten constant K.M
  KM.update = MH.KM.all(theta[rep-1,2] , KM.S, rep, RR.all[,1,], X.all, b = tun.KM, pri.KM = pri.KM)
  theta[rep,2] = KM.update$km;
  KM.S = KM.update$s
  count_KM = count_KM + KM.update$count
  
  X.fit[rep,,] = X.all
  R.fit[4*rep-3,,] = RR.all[,1,] # birth number of Y
  R.fit[4*rep-2,,] = RR.all[,2,] # death number of Y
  R.fit[4*rep-1,,] = RR.all[,3,] # birth number of X
  R.fit[4*rep-0,,] = RR.all[,4,] # death number of X
  
  # S.fit[rep,] = as.vector(Delta.X.S)
  if(rep%%100 ==0 ) cat("0")
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


colMeans(theta)
# plot(colMeans(birthY)); lines(birthY.sim);

result.par <- theta[selrow,]
result.X.trj <- X.fit[selrow,,]
# result.etc <- matrix(data= 0, nrow = max.T + 1, ncol = 4)
# result.etc[,1] <- sim.X; result.etc[,2] <- sim.Y; result.etc[1,3] <- rndseed; 
# result.etc[1,4] <- count_X / nrepeat; result.etc[2,4] <- count_Delta.X / nrepeat; result.etc[3,4] <- count_KM / nrepeat;  

currentT <- Sys.time()
timestamp <- paste(substr(currentT, 1,4), substr(currentT, 6,7), substr(currentT, 9,10),substr(currentT, 12,13), substr(currentT, 15,16), substr(currentT, 18,19), sep = "")

# compath <- "D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/ProfChoi_Code_20200326//";
compath <- "/home/users/hphong/TimeDelay//";

par.X.filename <- paste("Parameters_maxT", toString(max.T), "_", toString(timestamp), ".csv", sep = "")
X.trj.filename <- paste("Xtrj_maxT", toString(max.T), "_",  toString(timestamp), ".csv", sep = "")
etc.filename <- paste("Etc_maxT", toString(max.T), "_",  toString(timestamp), ".csv", sep = "")

# write.table(result.par, paste(compath, par.X.filename, sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
# write.table(result.X.trj, paste(compath, X.trj.filename, sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")
# write.table(result.etc, paste(compath, etc.filename, sep = ""), row.names = FALSE, col.names = FALSE, sep = ",")

### print result ###

# plot(A.Y * KI.Y(Delta.Y, in.X = colMeans(result.X.trj), K.M = mean(theta[,2]))[,1])
# lines(birthY.sim)
plot(tspan, rowMeans(sim.Y.all))
lines(tspan, mean.y, col = "red")

# gen.y3 <- matrix(0, nrow = effrepeat, ncol = max.T+1)
# for(jj in 1:effrepeat){
#   myList2 <- TimeDelayGillespieforXY(A.X = theta[selrow[jj],1], B.X = 0.05, alpha.X = theta[selrow[jj],3], beta.X = theta[selrow[jj],4], A.Y = 60, B.Y = 0.05, alpha.Y = 3.6, beta.Y = 0.6, K.M = theta[selrow[jj],2], repnum = max.T*500, maxT = max.T+3)
#   sim.X2 <- approx(myList2$TList[!is.na(myList2$TList)], myList2$XList[!is.na(myList2$XList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList2$XList[!is.na(myList2$XList)]))$y
#   gen.y3[jj,] <- approx(myList2$TList[!is.na(myList2$TList)], myList2$YList[!is.na(myList2$YList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList2$YList[!is.na(myList2$YList)]))$y
#   if(rep%%1000 ==0 ) cat("0")
# }
# mean.y2 <- colMeans(gen.y3)

# tmp <- 1177
# lines(tspan, gen.y3[tmp,], col = "blue")
# theta[tmp,]

colMeans(theta)

plot(theta[,1], type = "l")
plot(theta[,2], type = "l")
plot(theta[,1]/theta[,2], type = "l")
plot(theta[,3], type = "l")
plot(theta[,3]/theta[,4], type = "l")
mean(theta[1:2500,3]/theta[1:2500,4])
mean(theta[1:2500,3])
mean(theta[1:2500,4])
plot(theta[,3]/theta[,4]^2, type = "l")
plot(theta[,4], type = "l")

