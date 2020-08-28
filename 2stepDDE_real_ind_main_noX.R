setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode")
# setwd("C:/Users/HongHyukpyo/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/ProfChoi_Code_20200326")
# setwd("C:/Users/Hyukpyo Hong/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/ProfChoi_Code_20200326")
source('2stepDDE_functions.R')

rndseed <- round((as.numeric(Sys.time())*1000)%% 10000)
set.seed(rndseed)

int = 1 ;

# read data
raw.data = read.csv("2mM_IPTG.csv",header = F)

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

Y = nmlzdata1[,23] # choose one data.
max.T <- length(Y) - 1;
X <- Y # X is hidden, but we set it as Y. It is almost meaningless.

# known (fixed) parameters
alpha.Y = 5.89; beta.Y = 0.89*int;
B.X <- 0.015;  B.Y <- 0.015; A.Y <- 35.4;

# unknowns parameters
A.X = 35.4*int; alpha.X = 5.89; beta.X = 0.89*int;
K.M = 50; 

tspan <- 0:max.T

tun.B=c(50,50, 100, 100);
tun.X = seq(from=1, by=0.2, length.out=length(Y)) # tunning parameter for the setting 1.


pri.A.X <- c(10 * 1, 1); # non-informative prior for A.X
pri.alpha.X <- c(3.6 , 1); # inormative prior for alpha.X
pri.beta.X <- c(0.01, 0.01); # inormative prior for beta.X
pri.KM <- c(200* 0.01, 0.01); # non-informative prior for KM

tun.KM =1;
tun.Delta.X = c(1.0, 1);

effrepeat <- 3000
burn <- 0; thin = 1;
nrepeat <- burn + thin*effrepeat;

selrow <- seq(from = burn + thin, by = thin, length.out = effrepeat)

#initial value setting 
theta.X = c(A.X , B.X)
theta.Y = c(A.Y, B.Y, K.M)

Delta.X = c(alpha.X, beta.X) #initial & true values of delay parameter of X 
Delta.Y = c(alpha.Y, beta.Y) #initial & true values of delay parameter of Y 

Y.diff = diff(Y) #y(i+1) - y(i)
X.diff = diff(X) #x(i+1) - x(i)    
RR = matrix(0,ncol = 4, nrow = max.T) #saving number of reaction 

for (i in 1:max.T) {
  RR[i,1] = max(Y.diff[i],0)  # # of birth reaction of Y
  RR[i,2] = max(-Y.diff[i],0) # # of death reaction of Y
  RR[i,3] = max(X.diff[i],0)  # # of birth reaction of X
  RR[i,4] = max(-X.diff[i],0) # # of death reaction of X
}


################################################
# iteration start!!!
################################################

# matrix & vector for saving MCMC results
count_KM <- 0; count_X <- 0; count_Delta.X <- 0; 

theta <- matrix(0,nrow = nrepeat, ncol=4)
X.fit <- matrix(0,nrow = nrepeat, ncol = max.T+1)
R.fit <- matrix(0,nrow = 4*nrepeat, ncol = max.T)
# S.fit = matrix(0,nrow = nrepeat, ncol = 4)

#initial scales of KM, delata.X, A.X used in RAM method 
KM.S <- 10
Delta.X.S <- diag(2)

# edited by Hyukpyo. X is fixed as the true X data.

K.i <- KI(Delta.X, maxt = max.T)

AcceptR.fit <- matrix(0,nrow = nrepeat, ncol = 5)
AcceptR.star.fit <- matrix(0,nrow = nrepeat, ncol = 5)

KM.lik.fit <- rep(0, nrepeat)
KM.star.lik.fit <- rep(0, nrepeat)

ptnum <- 4
useall <- TRUE

theta[1,] = c(theta.X[1], theta.Y[3], Delta.X[1], Delta.X[2])
for(rep in 2:nrepeat) {
  # step 1 & 2: sampling  r2 and r1 (death and birth of Y)
  RR[,1:2] <- impute_r.Y(Y, B.Y = B.Y)
  
  # step 3: sampling X &   r3, r4 
  # updating X using independent chain MH
  
  # generate a proposal mean trajectory using the current parameter set.
  myListX <- TimeDelayGillespieforXR(A.X = theta[rep-1,1], B.X = B.X, alpha.X = theta[rep-1,3], beta.X = theta[rep-1,4], repnum = round(max.T*500), maxT = max.T+5)
  X.bir.st <- myListX$Xbirth[1:max.T]
  X.dea.st <- myListX$Xdeath[1:max.T]
  X.star <- c(0, cumsum(X.bir.st - X.dea.st));
  # print(X.update$errflg)
  if (useall == TRUE){
    fy.st = A.Y * KI.Y(Delta.Y,in.X = X.star, K.M=theta[rep-1,2])
    fy    = A.Y * KI.Y(Delta.Y,in.X = X     , K.M=theta[rep-1,2])
  }else{
    fy.st = A.Y * KI.Ynt(Delta.Y,in.X = X.star, N = ptnum, K.M=theta[rep-1,2])
    fy    = A.Y * KI.Ynt(Delta.Y,in.X = X     , N = ptnum, K.M=theta[rep-1,2])
  }
  q.Y.st = sum(log(dpois(RR[,1],fy.st[,1])+1e-300), na.rm = T)
  q.Y    = sum(log(dpois(RR[,1],fy[,1]   )+1e-300), na.rm = T)
  prior.X.st = sum(log(dgamma(X.star , shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
  prior.X   = sum(log(dgamma(X, shape = 1, rate = 1e-2) + 1e-300)) # non-informative gamma prior
  
  # logMH <- q.Y.st - q.Y + prior.X.st - prior.X;
  logMH <- q.Y.st - q.Y; # Completely non-informative, i.e., always prior.X.st == prior.X  
  # print(logMH);
  if(!is.nan(logMH) && runif(1)<exp(logMH)) {
    X=X.star; RR[,3] <- X.bir.st; RR[,4] <- X.dea.st;
    count_X = count_X + 1; 
  }
  
  # step  4: samping A.X 
  g_11 <- sum(K.i);
  theta[rep,1] = rgamma(1,shape = sum(RR[,3]) + pri.A.X[1], rate = g_11 + pri.A.X[2]);
  
  # step 5 & 6: sampling alpha.X and beta.X: the delay parameters for the birth reaction of X.
  p.update <- MH.P.X(P = theta[rep-1,3:4], Delta.X.S, rep, RR[,3], Ax = theta[rep,1], tun = tun.Delta.X, pri.alpha.X = pri.alpha.X, pri.beta.X = pri.beta.X, maxt = max.T)
  theta[rep,3:4] = p.update$P
  Delta.X.S = p.update$S
  count_Delta.X = count_Delta.X + p.update$count
  
  K.i <- KI(P = theta[rep,3:4], maxt = max.T);
  
  # step 7: sampling the Michaelis-Menten constant K.M
  KM.update = MH.KM(theta[rep-1,2] , KM.S, rep, RR[,1], X, b = tun.KM, pri.KM = pri.KM, Delta.Y = Delta.Y, flatpri = TRUE)
  theta[rep,2] = KM.update$km;
  KM.S = KM.update$s
  count_KM = count_KM + KM.update$count
  
  X.fit[rep,] = X
  R.fit[4*rep-3,] = RR[,1] # birth number of Y
  R.fit[4*rep-2,] = RR[,2] # death number of Y
  R.fit[4*rep-1,] = RR[,3] # birth number of X
  R.fit[4*rep-0,] = RR[,4] # death number of X
  
  # S.fit[rep,] = as.vector(Delta.X.S)
  if(rep%%1000 ==0 ) cat("0")
  if(theta[rep,1] > 300){
    print("Estimated Ax > 300")
    break
  } 
}

# the estimated reaction numbers from MCMC algorithm.
birthY<- R.fit[4*(selrow-1) + 1,]  
deathY<- R.fit[4*(selrow-1) + 2,]
birthX<- R.fit[4*(selrow-1) + 3,]
deathX<- R.fit[4*(selrow-1) + 4,]

print(paste0("Acceptance ratio for X: ", count_X / nrepeat))  # Acceptance ratio for X
print(paste0("Acceptance ratio for alpha.X and beta.X: ", count_Delta.X / nrepeat))  # Acceptance ratio for alpha.X and beta.X
print(paste0("Acceptance ratio for K.M: ", count_KM / nrepeat))  # Acceptance ratio for alpha.X and beta.X

# generate a Y trajectory from the mean of the estimated parameters from the effective iteration indexed by 'selrow'.
gen.num <- 10
gen.y2 <- matrix(0, nrow = gen.num, ncol = max.T+1)
for(jj in 1:gen.num){
  myList2 <- TimeDelayGillespieforXY(A.X = colMeans(theta[selrow,])[1], B.X = 0.015, alpha.X = colMeans(theta[selrow,])[3], beta.X = colMeans(theta[selrow,])[4], A.Y = 35.4, B.Y = 0.015, alpha.Y = 5.89, beta.Y = 0.89, K.M = colMeans(theta[selrow,])[2], repnum = max.T*500, maxT = max.T+3)
  sim.X2 <- approx(myList2$TList[!is.na(myList2$TList)], myList2$XList[!is.na(myList2$XList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList2$XList[!is.na(myList2$XList)]))$y
  gen.y2[jj,] <- approx(myList2$TList[!is.na(myList2$TList)], myList2$YList[!is.na(myList2$YList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList2$YList[!is.na(myList2$YList)]))$y
}
mean.y <- colMeans(gen.y2)


colMeans(theta)
# plot(colMeans(birthY)); lines(birthY.sim);

result.par <- theta[selrow,]
result.X.trj <- X.fit[selrow,]
result.etc <- matrix(data= 0, nrow = max.T + 1, ncol = 4)
result.etc[,1] <- sim.X; result.etc[,2] <- sim.Y; result.etc[1,3] <- rndseed; 
result.etc[1,4] <- count_X / nrepeat; result.etc[2,4] <- count_Delta.X / nrepeat; result.etc[3,4] <- count_KM / nrepeat;  
plot(tspan, sim.X/200);   


### print result ###

plot(A.Y * KI.Y(Delta.Y, in.X = colMeans(result.X.trj), K.M = mean(theta[,2]))[,1])
lines(birthY.sim)
plot(tspan, sim.Y)
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


plot(theta[,1], type = "l", xlab = "iteration", ylab = "AX")
plot(theta[,2], type = "l", xlab = "iteration", ylab = "KM")
plot(theta[,1]/theta[,2], type = "l", xlab = "iteration", ylab = "AX/KM")
plot(theta[,3], type = "l", xlab = "iteration", ylab = "alphaX")
plot(theta[,3]/theta[,4], type = "l", xlab = "iteration", ylab = "Mean Delay")
plot(theta[,3]/theta[,4]^2, type = "l", xlab = "iteration", ylab = "Var. Delay")
plot(theta[,4], type = "l", xlab = "iteration", ylab = "betaX") 


hist(theta[,1], breaks = 30)
hist(theta[,2], breaks = 30)
hist(theta[,1]/theta[,2], breaks = 30)
hist(theta[,3]/theta[,4], breaks = 30)


