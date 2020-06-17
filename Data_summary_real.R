# install.packages("Rfast")
# library(Rfast)

setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/CodeForSobolev/Data_from_simulation")

# For max.T == 25
filenameDigit <- c("20200531000829", "20200601235510", "20200601232603")[3]
max.T <- 25
# For max.T == 50
filenameDigit <- c("20200406234838", "20200406235258", "20200406235329", "20200406235656", "20200407001113")[1]
max.T <- 50
# For max.T == 100
filenameDigit <- c("20200407034525", "20200407041011", "20200407045155", "20200407062228", "20200407083413")[5]
max.T <- 100

raw.Table <- read.csv(file = paste0("Parameters_commonX_maxT", max.T ,"_", filenameDigit, ".csv", sep = ""), header = F)
raw.Table <- read.csv(file = paste0("Parameters_multiX_maxT", max.T ,"_", filenameDigit, ".csv", sep = ""), header = F)

# raw.Table <- read.csv(file = "Parameters_maxT25_20200406220145.csv", header = F)
# raw.Table <- read.csv(file = "Parameters_maxT25_20200406220555.csv", header = F)
# raw.Table <- read.csv(file = "Parameters_maxT25_20200406221108.csv", header = F)
# raw.Table <- read.csv(file = "Parameters_maxT25_20200406221149.csv", header = F)
true.X <- read.csv(file = paste0("Etc_maxT", max.T ,"_", filenameDigit, ".csv", sep = ""), header = F)[,1]
true.Y <- read.csv(file = paste0("Etc_maxT", max.T ,"_", filenameDigit, ".csv", sep = ""), header = F)[,2]
fit.X <- t(read.csv(file = paste0("Xtrj_maxT", max.T ,"_", filenameDigit, ".csv", sep = ""), header = F))
AR.A.X <- read.csv(file = paste0("Etc_maxT", max.T ,"_", filenameDigit, ".csv", sep = ""), header = F)[1,4]
AR.KM <- read.csv(file = paste0("Etc_maxT", max.T ,"_", filenameDigit, ".csv", sep = ""), header = F)[2,4]
AR.delay <- read.csv(file = paste0("Etc_maxT", max.T ,"_", filenameDigit, ".csv", sep = ""), header = F)[3,4]

raw.Table <- theta
true.Y <- Y

A.X <- raw.Table[,1]
KM <- raw.Table[,2]
alpha.X <- raw.Table[,3]
beta.X <- raw.Table[,4]
delayMean <- alpha.X/beta.X
delayVar <- alpha.X/beta.X^2

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
  myList <- TimeDelayGillespieforXY(A.X = mean(A.X), B.X = 0.015, alpha.X = mean(alpha.X), beta.X = mean(beta.X), A.Y = 35.4, B.Y = 0.015, alpha.Y = 5.89, beta.Y = 0.89, K.M = mean(KM), repnum = max.T*150, maxT = max.T+3)
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


plot(A.X, KM, xlab = "A.X", ylab = "K.M", log = "xy", cex.lab = 1.5, col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.5))
# points(x = 1,y = 1, cex=1.5, pch=17, col = "red")
plot(alpha.X, beta.X, xlab = "alpha.X", ylab = "beta.X", log = "xy", cex.lab = 1.5, col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.5))
# points(x = 1,y = 1, cex=1.5, pch=17, col = "red")
plot(delayMean, delayVar, xlab = "Delay Mean", ylab = "Delay Var", log ="xy", cex.lab = 1.5, col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.5))
# points(x = 1,y = 1, cex=1.5, pch=17, col = "red")
plot(A.X/KM , delayMean , xlab = "A.X/K.M", ylab = "Delay Mean", log = "xy", cex.lab = 1.5, col=rgb(red=0.2, green=0.2, blue=1.0, alpha=0.5))


# boxplot of the parameters normalized by the true parameter values
# boxplot(x = NmlzData, las = 1, names = c("A.X", "KM", "KM/A.X", "alpha.X", "beta.X", "delay Mean" ,"delay Var"),log = "y", cex.lab = 1, ylab = "Normalized & log scale", cex.axis = 1)

# auto-corrleation function plots of the parameters from MCMC
par(mfrow=c(2, 3))
acf(A.X, cex.lab = 1.5)
acf(KM, cex.lab = 1.5)
acf(alpha.X, cex.lab = 1.5)
acf(beta.X, cex.lab = 1.5)
acf(delayMean, cex.lab = 1.5)
acf(delayVar, cex.lab = 1.5)

# trace plot of the parameters from MCMC
colMeans(theta)
par(mfrow=c(1,1))

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


ggpairs(data = as.data.frame(totalData))
