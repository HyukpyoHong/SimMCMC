setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/ProfChoi_Code_20200326")
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

Y = nmlzdata1[,1] # choose one data.
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

pri.A.X = c(0.01, 0.01); # non-informative prior for A.X
pri.alpha.X=c(alpha.X*1, 1); # inormative prior for alpha.X
pri.beta.X=c(beta.X*1,1); # inormative prior for beta.X
pri.KM = c(1*0.01, 0.01 ); # non-informative prior for KM

tun.KM =1;
tun.Delta.X = c(1.0, 1);

effrepeat <- 3000
burn <- 1000; thin = 3;
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
count_R = rep(0,max.T); count_KM=0;count_X = 0;
count_al <- 0; count_be <- 0; count_Delta.X <- 0;count_A.X <- 0;

theta = matrix(0,nrow = nrepeat, ncol=4)
X.fit = matrix(0,nrow = nrepeat, ncol = max.T+1)
R.fit = matrix(0,nrow = 4*nrepeat, ncol = max.T)
S.fit = matrix(0,nrow = nrepeat, ncol = 4)

#initial scales of KM, delata.X, A.X used in RAM method 
KM.S = 10; #A.X.S=1;
Delta.X.S = diag(2);

# edited by Hyukpyo. X is fixed as the true X data.

K.i <- KI(Delta.X, maxt = max.T); 

gen.num <- 3
gen.x <- matrix(0, nrow = gen.num, ncol = length(X))

for (rep in 1:nrepeat) {
  # step 1 & 2: sampling  r2 and r1 (death and birth of Y)
  RR[,1:2] = impute_r.Y(Y)
  
  # step 3: sampling X &   r3, r4 
  # updating X using independent chain MH
  
  #step 3: sampling X &   r3, r4 
  # updating X
  X.update = MH.X(RR,X,kix=K.i,km=K.M, thetax=theta.X, ptnum = ptnum, useall = TRUE)
  X = X.update$x
  count_X = count_X + X.update$count
  
  # updating number of birth and death reaction of X
  R.update = MH.R(RR,X,ki=K.i,thetax=theta.X);# if b=tun.B then remove $ from bi=b[i,]
  # R.update = MH.R.repeat.ave(RR,X,ki=K.i,thetax=theta.X);# if b=tun.B then remove $ from bi=b[i,]
  # R.update = MH.R.repeat.ave(RR,X,ki=K.i,thetax=theta.X);# if b=tun.B then remove $ from bi=b[i,]
  RR = R.update$r
  count_R = count_R + R.update$count
  
  # step  4: samping A.X 
  g_11 <- sum(K.i);
  theta.X[1] = rgamma(1,shape = sum(RR[,3]) + pri.A.X[1], rate = g_11 + pri.A.X[2]);
  theta[rep,1] =theta.X[1]
  
  # step 5 & 6: sampling alpha.X and beta.X: the delay parameters for the birth reaction of X.
  p.update <- MH.P.X(Delta.X, Delta.X.S, rep, RR[,3:4], theta.X[1], maxt = max.T)
  Delta.X = p.update$P
  Delta.X.S = p.update$S
  count_Delta.X = count_Delta.X + p.update$count
  
  K.i <- KI(Delta.X, maxt = max.T);
  theta[rep,3:4]=Delta.X
  
  # step 7: sampling the Michaelis-Menten constant K.M
  KM.update = MH.KM(theta.Y[3] , KM.S, rep, RR[,1], X)
  theta.Y[3] = KM.update$km;
  K.M=theta.Y[3];
  KM.S = KM.update$s
  count_KM = count_KM + KM.update$count
  
  theta[rep,2] = theta.Y[3]
  X.fit[rep,] = X
  R.fit[4*rep-3,] = RR[,1] # birth number of Y
  R.fit[4*rep-2,] = RR[,2] # death number of Y
  R.fit[4*rep-1,] = RR[,3] # birth number of X
  R.fit[4*rep-0,] = RR[,4] # death number of X
  
  # S.fit[rep,] = as.vector(Delta.X.S)
  if(rep%%1000 ==0 ) cat("0")
}

# the estimated reaction numbers from MCMC algorithm.
birthY<- R.fit[4*(selrow-1) + 1,]  
deathY<- R.fit[4*(selrow-1) + 2,]
birthX<- R.fit[4*(selrow-1) + 3,]
deathX<- R.fit[4*(selrow-1) + 4,]

print(paste0("Acceptance ratio for X: ", count_X / nrepeat))  # Acceptance ratio for X
print(paste0("Acceptance ratio for K.M: ", count_KM / nrepeat))  # Acceptance ratio for alpha.X and beta.X
print(paste0("Acceptance ratio for alpha.X and beta.X: ", count_Delta.X / nrepeat))  # Acceptance ratio for alpha.X and beta.X

plot(theta[,1])
plot(theta[,2])

plot(theta[,2]/theta[,1])

# generate a Y trajectory from the mean of the estimated parameters from the effective iteration indexed by 'selrow'.
gen.num <- 10
gen.y2 <- matrix(0, nrow = gen.num, ncol = max.T+1)
for(jj in 1:gen.num){
  myList2 <- TimeDelayGillespieforXY(A.X = colMeans(theta[selrow,])[1], B.X = 0.05, alpha.X = colMeans(theta[selrow,])[3], beta.X = colMeans(theta[selrow,])[4], A.Y = 60, B.Y = 0.05, alpha.Y = 3.6, beta.Y = 0.6, K.M = colMeans(theta[selrow,])[2], repnum = max.T*150, maxT = max.T+3)
  sim.X2 <- approx(myList2$TList[!is.na(myList2$TList)], myList2$XList[!is.na(myList2$XList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList2$XList[!is.na(myList2$XList)]))$y
  gen.y2[jj,] <- approx(myList2$TList[!is.na(myList2$TList)], myList2$YList[!is.na(myList2$YList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList2$YList[!is.na(myList2$YList)]))$y
}
mean.y <- colMeans(gen.y2)

plot(tspan, mean.y); lines(tspan, Y, col = "red");

result.par <- theta[selrow,]
result.X.trj <- X.fit[selrow,]
result.etc <- matrix(data= 0, nrow = max.T + 1, ncol = 4)
result.etc[1,4] <- count_X / nrepeat; result.etc[2,4] <- count_Delta.X / nrepeat; result.etc[3,4] <- count_KM / nrepeat;  

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


