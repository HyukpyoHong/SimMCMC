setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode/")
source("2stepDDE_functions.R")
# setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode/")
# setwd("/Users/hyukpyohong/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode/")
# setwd("/Users/hyukpyohong/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/CodeForSobolev/Data_from_simulation/")


args1 = 1

flg1 = 0
set_comb = matrix(NA, nrow = 60, ncol = 5)

for(effrepeat in c(1100000)){
  for(max.T in c(30, 50, 150)){
    for(nsample in c(5, 10, 20, 40)){
      for(param_num in c(3)){
        for (km_ratio in c(0.2, 0.5, 1, 2, 5)){
          flg1 = flg1 +1
          set_comb[flg1, ] = c(effrepeat, max.T, nsample, param_num, km_ratio)
        }
      }
    }
  }
}

effrepeat = set_comb[args1, 1]
max.T = set_comb[args1, 2]
nsample = set_comb[args1, 3]
param_num = set_comb[args1, 4]
km_ratio = set_comb[args1, 5]


if (param_num == 1){
  param_est <- c(1,1,1,1,1,1,1,0) # AX KM alphaX betaX AY alphaY betaY B
}else if(param_num == 2){
  param_est <- c(1,1,1,1,1,0,0,0) # AX KM alphaX betaX AY alphaY betaY B
}else if(param_num == 3){
  param_est <- c(1,0,1,1,1,0,0,0) # AX KM alphaX betaX AY alphaY betaY B
}else if(param_num == 4){
  param_est <- c(1,1,1,1,0,0,0,0) # AX KM alphaX betaX AY alphaY betaY B
}else if(param_num == 5){
  param_est <- c(0,1,1,1,0,0,0,0) # AX KM alphaX betaX AY alphaY betaY B
}else{
  param_est <- c(0,0,1,1,0,0,0,0) # AX KM alphaX betaX AY alphaY betaY B
}

int <- 1; 

# known (fixed) parameters
B.X <- 0.05*int;  
B.Y <- 0.05*int; 
A.Y <- 60*int; 
alpha.Y <- 3.6; beta.Y <- 0.6*int;  

# unknowns parameters
A.X <- 10*int; alpha.X <- 3.6; beta.X <- 0.6*int; 
K.M <- 200; 

tspan <- 0:max.T

# install.packages("stringr")
library(stringr)


for(fileidx in 1:dim(filename.list)[1]){
  filename.current = filename.list[fileidx, 1]
  param_num = as.numeric(substr(filename.current, 22,22))
  theta.est = read.table(paste("Theta", substr(filename.current,7,str_length(filename.current)), sep = ""), sep = ",")
  repnum = sum(theta.est[,1] > 0)
  reaction.numbers.true = read.table(paste("Reaction_numbers", substr(filename.current,7,str_length(filename.current)), sep = ""), sep = ",")
  max.T = dim(reaction.numbers.true)[1]
  nsample = dim(reaction.numbers.true)[2]/4
  update.numbers.est = read.table(paste("Update_numbers", substr(filename.current,7,str_length(filename.current)), sep = ""), sep = ",")
  priors.true = read.table(paste("Priors", substr(filename.current,7,str_length(filename.current)), sep = ""), sep = ",")
  
  Xbirth.true = reaction.numbers.true[,1:nsample]
  Xdeath.true = reaction.numbers.true[,(nsample+1):(2*nsample)]
  Ybirth.true = reaction.numbers.true[,(2*nsample+1):(3*nsample)]
  Ydeath.true = reaction.numbers.true[,(3*nsample+1):(4*nsample)]
  
  if (param_num == 1 & max.T = 150){
    selrow = seq(from = 510, by = 10, to = repnum)
    theta.est
  }
}


fileidx = 1

filename.current = filename.list[fileidx, 1]
param_num = as.numeric(substr(filename.current, 22,22))
theta.est = read.table(paste("Theta", substr(filename.current,7,str_length(filename.current)), sep = ""), sep = ",")
repnum = sum(theta.est[,1] > 0)
selrow = seq(from = 1010, by = 10, to = repnum)
theta.est = theta.est[selrow,]
reaction.numbers.true = read.table(paste("Reaction_numbers", substr(filename.current,7,str_length(filename.current)), sep = ""), sep = ",")
max.T = dim(reaction.numbers.true)[1]
nsample = dim(reaction.numbers.true)[2]/4
update.numbers.est = read.table(paste("Update_numbers", substr(filename.current,7,str_length(filename.current)), sep = ""), sep = ",")
priors.true = read.table(paste("Priors", substr(filename.current,7,str_length(filename.current)), sep = ""), sep = ",")

Xbirth.true = reaction.numbers.true[,1:nsample]
Xdeath.true = reaction.numbers.true[,(nsample+1):(2*nsample)]
Ybirth.true = reaction.numbers.true[,(2*nsample+1):(3*nsample)]
Ydeath.true = reaction.numbers.true[,(3*nsample+1):(4*nsample)]

param.name.list = c("AX", "KM", "alphaX", "betaX", "AY", "alphaY", "betaY", "B")
par(mfrow = c(3, 3))
for (param.id in 1:8){
  plot(theta.est[,param.id], type = "l" , col = rgb(0/255,0/255,0/255, 0.6),
       xlab = "", ylab = "", bty = "n")
  # xlim = c(0,15), ylim = c(0,15)
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  mtext(side=1, line=2, "Iteration", font=1.6,cex=1.2)
  mtext(side=2, line=2, param.name.list[param.id], font=1.6,cex=1.2)
}

plot(theta.est[,3]/theta.est[,4], type = "l" , col = rgb(0/255,0/255,0/255, 0.6),
     xlab = "Iteration", ylab = "X delay", bty = "n")
# xlim = c(0,15), ylim = c(0,15)
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
mtext(side=1, line=2, "Iteration", font=1.6,cex=1.2)
mtext(side=2, line=2, "", font=1.6,cex=1.2)

par(mfrow = c(1,1))
plot(theta.est[,3]/theta.est[,4],theta.est[,6]/theta.est[,7], col = rgb(0/255,0/255,0/255, 0.6),
     xlab = "", ylab = "", bty = "n", xlim = c(0,15), ylim = c(0,15))
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
lines(c(0,12), c(12,0), lwd = 2, col = rgb(216/255,75/255,36/255, 0.8))
mtext(side=1, line=2, "X delay", font=1.6,cex=1.2)
mtext(side=2, line=2, "Y delay", font=1.6,cex=1.2)


plot(theta.est[,1]/theta.est[,2],col = rgb(0/255,0/255,0/255, 0.6),
     xlab = "", ylab = "", bty = "n")

plot(theta.est[,4],col = rgb(0/255,0/255,0/255, 0.6),
     xlab = "", ylab = "", bty = "n")
theta.est[,3]/theta.est[,4]

plot(theta.est[,2], type = "l")
plot(theta.est[,3], type = "l")
plot(theta.est[,4], type = "l")
plot(theta.est[,5], type = "l")
plot(theta.est[,6], type = "l")
plot(theta.est[,7], type = "l")
plot(theta.est[,8], type = "l")

totalData <- as.matrix(cbind(A.X, KM, KM/A.X, alpha.X, beta.X, delayMean, delayVar))
NmlzData <- totalData / (matrix(rep(c(10,200,20,3.6,0.6,6,10), length(A.X)),ncol = 7,byrow=TRUE))
dataVar <- (totalData - matrix(rep(colMeans(totalData), length(A.X)),ncol = 7,byrow=TRUE))^2
resultVar <- colSums(dataVar)/(length(A.X)-1)
resultMean <- colMeans(totalData)

# Print the mean and variance of the estimated parameters
print(resultMean)
print(resultVar)
print(cbind(AR.A.X, AR.KM, AR.delay))
# print(mean(KM/A.X))
# print(var(KM/A.X))

# Comparison between mean of generated Y from the estimated parameters and true Y data

gen.num <- 200
gen.y <- matrix(0, nrow = gen.num, ncol = max.T+1)
for(jj in 1:gen.num){
  myList <- TimeDelayGillespieforXY(A.X = mean(A.X), B.X = 0.05, alpha.X = mean(alpha.X), beta.X = mean(beta.X), A.Y = 60, B.Y = 0.05, alpha.Y = 3.6, beta.Y = 0.6, K.M = mean(KM), repnum = max.T*150, maxT = max.T+3)
  sim.X <- approx(myList$TList[!is.na(myList$TList)], myList$XList[!is.na(myList$XList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList$XList[!is.na(myList$XList)]))$y
  gen.y[jj,] <- approx(myList$TList[!is.na(myList$TList)], myList$YList[!is.na(myList$YList)], xout = seq(from = 0, to = max.T, by=1), method = "constant", yleft = 0, yright = max(myList$YList[!is.na(myList$YList)]))$y
}
mean.y <- colMeans(gen.y)
times <- 0:max.T

uQ.y <- rep(0,max.T+1)
lQ.y <- rep(0,max.T+1)
for(jj in 1:(max.T+1)){
  lQ.y[jj] <- quantile(gen.y[,jj], probs = c(0.025, 0.975))[1]
  uQ.y[jj] <- quantile(gen.y[,jj], probs = c(0.025, 0.975))[2]
}

par(mfrow = c(1,1))
plot(times, true.Y, ylim = c(0,max(true.Y)), xlab = "time", ylab = "Y trajectory", cex.lab = 1.5)
polygon(c(times, rev(times)), c(lQ.y ,rev(uQ.y)), col = rgb(0, 0, 0, 0.3) )
lines(times, mean.y, col = "red", lwd=2)

plot(times, rowMeans(fit.X/matrix(rep(KM, max.T+1), ncol = length(KM), byrow = TRUE)), ylim = c(0,max(true.X/200)), xlab = "time", ylab = "X trajectory", cex.lab = 1.5)
# lines(times, rowMeans(fit.X)/mean(KM), ylim = c(0,max(true.X/200)), xlab = "time", ylab = "X(t) / K.M", cex.lab = 1.5)
lines(times, true.X/200, col = "red", lwd = 2)

# 
plot(A.X/10, KM/200, xlab = "Norm. A.X", ylab = "Norm. K.M", log = "xy", cex.lab = 1.5, col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.5))
points(x = 1,y = 1, cex=1.5, pch=17, col = "red")
plot(alpha.X/3.6, beta.X/0.6, xlab = "Norm. alpha.X", ylab = "Norm. beta.X", log = "xy", cex.lab = 1.5, col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.5))
points(x = 1,y = 1, cex=1.5, pch=17, col = "red")
plot(delayMean/6, delayVar/10, xlab = "Norm. delay Mean", ylab = "Norm. delay Var", log ="xy", cex.lab = 1.5, col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.5))
points(x = 1,y = 1, cex=1.5, pch=17, col = "red")

plot(A.X/KM * 20, delayMean/6, xlab = "Norm. A.X/K.M", ylab = "Norm. delay Mean", log = "xy", cex.lab = 1.5, col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.5))
points(x = 1,y = 1, cex=1.5, pch=17, col = "red")


# boxplot of the parameters normalized by the true parameter values
boxplot(x = NmlzData, las = 1, names = c("A.X", "KM", "KM/A.X", "alpha.X", "beta.X", "delay Mean" ,"delay Var"),log = "y", cex.lab = 1, ylab = "Normalized & log scale", cex.axis = 1)

# auto-corrleation function plots of the parameters from MCMC
par(mfrow=c(2, 3))
acf(A.X, cex.lab = 1.5)
acf(KM, cex.lab = 1.5)
acf(alpha.X, cex.lab = 1.5)
acf(beta.X, cex.lab = 1.5)
acf(delayMean, cex.lab = 1.5)
acf(delayVar, cex.lab = 1.5)

colMeans(theta)
par(mfrow=c(1,1))

plot(theta[,1], type = "l", xlab = "iteration", ylab = "AX")
plot(theta[,2], type = "l", xlab = "iteration", ylab = "KM")
plot(theta[,1]/theta[,2], type = "l", xlab = "iteration", ylab = "AX/KM")
plot(theta[,3], type = "l", xlab = "iteration", ylab = "alphaX")
plot(theta[,3]/theta[,4], type = "l", xlab = "iteration", ylab = "Mean Delay")
plot(theta[,3]/theta[,4]^2, type = "l", xlab = "iteration", ylab = "Var. Delay")
plot(theta[,4], type = "l", xlab = "iteration", ylab = "betaX") 

ggpairs(data = as.data.frame(totalData))

## ===========2021.04.14. sum of time delays are preserved.

plot(theta[,3]/theta[,4],theta[,6]/theta[,7], 
     xlim = c(0, 15), ylim = c(0,15), col = rgb(100/255,100/255,100/255, 0.3),
     xlab = "", ylab = "", bty = "n", main = NULL)
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
mtext(side=1, line=2, "alpha_X", font=1.6,cex=1.2)
mtext(side=2, line=2, "beta_X", font=1.6,cex=1.2)
# points(mean(param_now[seq(from = 1000, to = 10000, by = 10),3]), mean(param_now[seq(from = 1000, to = 10000, by = 10),4]), col=rgb(10/255,10/255,10/255, 0.8), lwd = 3, pch = 16)
# points(3.6,0.6, col=rgb(243/255,135/255,47/255, 0.8), lwd = 2, pch = 1)

plot(theta[,3]/theta[,4], type = "l")
plot(theta[,3]/theta[,4] + theta[,6]/theta[,7], type = "l")

resultMean
resultVar

