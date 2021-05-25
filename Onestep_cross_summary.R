setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode")
setwd("C:/Users/HongHyukpyo/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode")
source('2stepDDE_functions.R')

linearTF <- as.numeric(args[1]) 
tmppar <- as.numeric(args[2]) 

rndseed <- round((as.numeric(Sys.time())*1000)%% 10000)
set.seed(rndseed)
ptm <- proc.time()

int <- 1; 

# known (fixed) parameters
B.X <- 0.05*int;  
B.Y <- 0.05*int; 
A.Y <- 60*int; 
alpha.Y <- 3.6; beta.Y <- 0.6*int;  

# unknowns parameters
A.X <- 10*int; alpha.X <- 3.6; beta.X <- 0.6*int; 
K.M <- tmppar; 

effrepeat <- 11000

if (linearTF == 1){
  K.M = 10^8
  A.Y = K.M * tmppar
}
max.T <- 150 # simulated data will be given from t = 0, ..., max.T
tspan <- 0:max.T
nsample <- 10;

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

Y.all <- sim.Y.all
X.all <- sim.X.all




#example of YFP1_short data set. using 40 frajectories 
short1 <- BDD(indata=Y.all,scale=1/1, B=0.00,nrepeat = 11000, tun.B=300, tun.al = 0.5, tun.be = 70000)


# for (kmval in c(1,10,100,200,1000)){
for(kmval in c(1,10,100,200,1000)){
  load(paste("UsingOnestep_linearTF0tmppar",kmval,".RData", sep = ""))
  selrow = seq(from = 1000, by = 10, to = 11000)
  int <- 1; 
  # known (fixed) parameters
  B.X <- 0.05*int;  
  B.Y <- 0.05*int; 
  A.Y <- 60*int; 
  alpha.Y <- 3.6; beta.Y <- 0.6*int;  
  
  # unknowns parameters
  A.X <- 10*int; alpha.X <- 3.6; beta.X <- 0.6*int; 
  # K.M <- tmppar; 
  
  K.M = kmval
  A.Y = 60
  max.T <- 150 # simulated data will be given from t = 0, ..., max.T
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
    myList <- TimeDelayGillespieforXY(A.X = A.X, B.X = B.X, alpha.X = alpha.X, beta.X = beta.X, A.Y = A.Y, B.Y = B.Y, alpha.Y = alpha.Y, beta.Y = beta.Y, K.M = K.M, repnum = max.T*10000, maxT = max.T+3)
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
  
  errX = TimeDelayGillespieforX(A.X = mean(short1[selrow,1]), B.X = mean(short1[selrow,2]), alpha.X = mean(short1[selrow,3]), beta.X = mean(short1[selrow,4]), repnum =  10000000, maxT = 150)
  plot(0:max.T, rowMeans(Y.all), type = "l", xlab = "", ylab = "", lwd = 2, 
       col = rgb(35/255,110/255,150/255, 0.8), lty = 1, 
       bty="n", ylim = c(0,round(max(rowMeans(Y.all))/100)*100+100))
  
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  mtext(side=1, line=2, "Time (min)", font=1.6,cex=1.2)
  mtext(side=2, line=2, "Y(t)", font=1.6,cex=1.2)
  # title(main = "Simulated trajectories")
  lines(errX$TList, errX$XList, lwd = 2, 
        col = rgb(243/255,135/255,47/255, 0.8))
  
  cat(kmval); cat("\n")
  cat(colMeans(short1[selrow,])); cat("\n")
  cat(mean(short1[selrow,3]/short1[selrow,4])); cat("\n")
}


est_delay_tot = c(9.14, 11.63, 18.24, 25.49)
kmval = c(1,10,100,1000)
plot(kmval, est_delay_tot, log = "x", axes = FALSE, ylim = c(5,30), lwd =3, pch=16)
axis(1, at=kmval, labels=c(10^0, 10^1, 10^2, 10^3), lwd = 2)
axis(2, lwd = 2)


load("UsingOnestep_linearTF1tmppar0.1.RData")
load("UsingOnestep_linearTF1tmppar0.2.RData")
load("UsingOnestep_linearTF1tmppar0.3.RData")
load("UsingOnestep_linearTF1tmppar0.5.RData")
load("UsingOnestep_linearTF1tmppar1.RData")
short1
par(mfrow=c(1,1))
matplot(short1, type ="l")
selrow = seq(from = 1000, by = 5, to = 11000)
matplot(short1[selrow,1], type ="l")
matplot(short1[selrow,2], type ="l")
matplot(short1[selrow,3], type ="l")
matplot(short1[selrow,4], type ="l")
matplot(short1[selrow,3]/short1[selrow,4], type ="l")
mean(short1[selrow,3]/short1[selrow,4])

colMeans(short1[selrow,])
errX = TimeDelayGillespieforX(A.X = mean(short1[selrow,1]), B.X = mean(short1[selrow,2]), alpha.X = mean(short1[selrow,3]), beta.X = mean(short1[selrow,4]), repnum =  10000000, maxT = 150)

int <- 1; 

# known (fixed) parameters
B.X <- 0.05*int;  
B.Y <- 0.05*int; 
A.Y <- 60*int; 
alpha.Y <- 6; beta.Y <- 1*int;  

# unknowns parameters
A.X <- 10*int; alpha.X <- 6; beta.X <- 1*int; 
# K.M <- tmppar; 

effrepeat <- 11000

K.M = 10^8
A.Y = 10^8
max.T <- 150 # simulated data will be given from t = 0, ..., max.T
tspan <- 0:max.T
nsample <- 10;

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

Y.all <- sim.Y.all
X.all <- sim.X.all

lines(rowMeans(Y.all))
errX2 = TimeDelayGillespieforX(A.X = 10, B.X = 0.05, alpha.X = 13, beta.X = 1, repnum =  10000000, maxT = 150)

plot(errX2$TList, errX2$XList)










######  === 2021.04.09 ======================
matplot(0:max.T, X.all, lwd = 3, col = rgb(243/255,135/255,47/255, 0.8), type = "l",)

matplot(0:max.T, Y.all, lwd = 3, 
        col = rgb(243/255,135/255,47/255, 0.4), lty = 1, type = "l",
        xlab = "", ylab = "",
        bty="n")
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
mtext(side=1, line=2, "Time (min)", font=1.6,cex=1.2)
mtext(side=2, line=2, "Y(t)", font=1.6,cex=1.2)
# title(main = "Simulated trajectories")

setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/CodeForSobolev/Data_from_simulation")

paramset4 = as.matrix(read.table(file = "Parameters_set4_maxT150_20210119043630.csv", sep =","))[,1:8]
paramset1 = as.matrix(read.table(file = "Parameters_set1_maxT150_20210119045047.csv", sep =","))[,1:8]
param_now = paramset1
plot(param_now[seq(from = 1000, to = 10000, by = 10),3]/param_now[seq(from = 1000, to = 10000, by = 10),4], 
     param_now[seq(from = 1000, to = 10000, by = 10),6]/param_now[seq(from = 1000, to = 10000, by = 10),7],
     col = rgb(0/255,0/255,0/255, 0.6),
     xlab = "", ylab = "", xlim = c(0,15), ylim = c(0,15), bty = "n")
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
mtext(side=1, line=2, "Time delay of X", font=1.6,cex=1.2)
mtext(side=2, line=2, "Time delay of Y", font=1.6,cex=1.2)

hist(param_now[seq(from = 1000, to = 10000, by = 10),3]/param_now[seq(from = 1000, to = 10000, by = 10),4], col="grey", border="white",
     prob=T, ylim=c(0,0.6),
     xlab = "", ylab = "", bty = "n", main = NULL)
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
mtext(side=1, line=2, "Time delay of X", font=1.6,cex=1.2)
mtext(side=2, line=2, "Posterior density", font=1.6,cex=1.2)
lines(c(6,6), c(0,0.6), col=rgb(243/255,135/255,47/255, 0.8), lwd = 3)
set_mtd = mean(param_now[seq(from = 1000, to = 10000, by = 10),3]/param_now[seq(from = 1000, to = 10000, by = 10),4])
lines(c(set_mtd,set_mtd), c(0,0.6), col=rgb(10/255,10/255,10/255, 0.8), lwd = 3)


plot(param_now[seq(from = 1000, to = 10000, by = 10),1],param_now[seq(from = 1000, to = 10000, by = 10),2], 
     xlim = c(6, 18), ylim = c(120,360), col = rgb(100/255,100/255,100/255, 0.3),
     xlab = "", ylab = "", bty = "n", main = NULL)
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
mtext(side=1, line=2, "Production rate of X", font=1.6,cex=1.2)
mtext(side=2, line=2, "MM constant", font=1.6,cex=1.2)
points(mean(param_now[seq(from = 1000, to = 10000, by = 10),1]), mean(param_now[seq(from = 1000, to = 10000, by = 10),2]), col=rgb(10/255,10/255,10/255, 0.8), lwd = 3, pch = 16)
points(10,200, col=rgb(243/255,135/255,47/255, 0.8), lwd = 2, pch = 1)

plot(param_now[seq(from = 1000, to = 10000, by = 10),3],param_now[seq(from = 1000, to = 10000, by = 10),4], 
     xlim = c(0, 12), ylim = c(0,3), col = rgb(100/255,100/255,100/255, 0.3),
     xlab = "", ylab = "", bty = "n", main = NULL)
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
mtext(side=1, line=2, "alpha_X", font=1.6,cex=1.2)
mtext(side=2, line=2, "beta_X", font=1.6,cex=1.2)
points(mean(param_now[seq(from = 1000, to = 10000, by = 10),3]), mean(param_now[seq(from = 1000, to = 10000, by = 10),4]), col=rgb(10/255,10/255,10/255, 0.8), lwd = 3, pch = 16)
points(3.6,0.6, col=rgb(243/255,135/255,47/255, 0.8), lwd = 2, pch = 1)




matplot(0:max.T, X.all, lwd = 1.5, col = rgb(243/255,135/255,47/255, 0.3), type = "l")
# 
# install.packages("ggplot2")
# library(ggplot2)
# ggplot(data = xx, aes(x = time, y= trjs.1)) + geom_line()
# 
# xx = data.frame("time" = 0:max.T, "trjs" = Y.all)



### ============================


setwd("C:/Users/HongHyukpyo/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode/")
IPTG2 = as.matrix(read.table("2mM_IPTG.csv", sep =","))


paramset4 = as.matrix(read.table(file = "Parameters_set4_maxT150_20210119043630.csv", sep =","))[,1:8]
paramset1 = as.matrix(read.table(file = "Parameters_set1_maxT150_20210119045047.csv", sep =","))[,1:8]
param_now = paramset1
matplot(0:46, IPTG2 - min(IPTG2), col = rgb(243/255,135/255,47/255, 0.6), type = "l",
     xlab = "", ylab = "", lwd = 2.0, lty = 1,
     xlim = c(0,46), ylim = c(0,320), 
     bty = "n")
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
mtext(side=1, line=2, "Time (min)", font=1.6,cex=1.2)
mtext(side=2, line=2, "Time delay of Y", font=1.6,cex=1.2)





