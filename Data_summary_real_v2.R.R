setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode/")
source("2stepDDE_functions.R")
setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/CodeForSobolev/realdata")



# setwd("/Users/hyukpyohong/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/timedelayGitCode/")
# setwd("/Users/hyukpyohong/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/CodeForSobolev/Data_from_simulation/")


args1 = 1

flg1 = 0
set_comb = matrix(NA, nrow = 72, ncol = 4)

for(effrepeat in c(1100000)){
  for(max.T in c(30, 50, 150)){
    for(nsample in c(5, 10, 20, 40)){
      for(param_num in c(1,2,3,4,5,6)){
        flg1 = flg1 +1
        set_comb[flg1, ] = c(effrepeat, max.T, nsample, param_num)
      }
    }
  }
}

effrepeat = set_comb[args1, 1]
max.T = set_comb[args1, 2]
nsample = set_comb[args1, 3]
param_num = set_comb[args1, 4]

if (param_num == 1){
  param_est <- c(1,1,1,1,1,1,1,0) # AX KM alphaX betaX AY alphaY betaY B
}else if(param_num == 2){
  param_est <- c(1,1,1,1,1,0,0,0) # AX KM alphaX betaX AY alphaY betaY B
}else if(param_num == 3){
  param_est <- c(0,0,1,1,1,0,0,0) # AX KM alphaX betaX AY alphaY betaY B
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

setwd("D:/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/CodeForSobolev/realdata")
filename.list = read.table("filenamelist_realdata_lowbnds.csv", sep =",")

for (fileidx in 1:15){
  filename.current = filename.list[fileidx, 1]
  concentration_idx = as.numeric(substr(filename.current, 15,15))

  group_idx = as.numeric(substr(filename.current, 19,19))
  param_num = as.numeric(substr(filename.current, 35,35))
  if (concentration_idx == 1){
    raw.data = as.matrix(read.table("005mM_IPTG.csv", sep = ","))
  }else if (concentration_idx == 2){
    raw.data = as.matrix(read.table("02mM_IPTG.csv", sep = ","))
  }else{
    raw.data = as.matrix(read.table("2mM_IPTG.csv", sep = ","))
  }
  
  
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
  
  if (concentration_idx == 1){
    nmlzdata1 = round(diffdata/0.089)[,1:24]; # data from experiment 1
    nmlzdata2 = round(diffdata/0.089)[,25:68]; # data from experiment 2
    
    rownum_without_outlier1 <- which(nmlzdata1[74,] > 20 & nmlzdata1[74,] < 80)
    rownum_without_outlier2 <- which(nmlzdata2[74,] > 20 & nmlzdata2[74,] < 80)
  }else if (concentration_idx == 2){
    nmlzdata1 = round(diffdata/0.089)[,1:33]; # data from experiment 1
    nmlzdata2 = round(diffdata/0.089)[1:45,34:73]; # data from experiment 2
    
    rownum_without_outlier1 <- which(nmlzdata1[51,] > 250 & nmlzdata1[51,] < 800)
    rownum_without_outlier2 <- which(nmlzdata2[45,] > 270 & nmlzdata2[45,] < 700)
  }else{
    nmlzdata1 = round(diffdata/0.089)[1:27,1:49]; # data from experiment 1
    nmlzdata2 = round(diffdata/0.089)[1:32,50:88]; # data from experiment 2
    
    rownum_without_outlier1 <- which(nmlzdata1[27,] > 200 & nmlzdata1[27,] < 600)
    rownum_without_outlier2 <- which(nmlzdata2[32,] > 200 & nmlzdata2[32,] < 500)
  }
  
  
  
  if (group_idx == 1){
    Y.all <- nmlzdata1[, rownum_without_outlier1]
  }else{
    Y.all <- nmlzdata2[, rownum_without_outlier2]
  }
  
  
  theta.est.raw = as.matrix(read.table(paste("Theta", substr(filename.current,7,str_length(filename.current)), sep = ""), sep = ","))
  repnum = sum(theta.est.raw[,1] > 0)
  selrow = seq(from = 1010, by = 10, to = repnum)
  theta.est = theta.est.raw[selrow,]
  max.T = dim(Y.all)[1] - 1
  nsample = dim(Y.all)[2]
  priors.true = read.table(paste("Priors", substr(filename.current,7,str_length(filename.current)), sep = ""), sep = ",")
  concent_name = c("0.05mM_IPTG","0.2mM_IPTG","2mM_IPTG")
  title_text = paste("Real data lowbnds ", concent_name[concentration_idx], " group", group_idx," Param set", param_num, " maxT", max.T, " nsample", nsample, " repnum", repnum, sep ="")
  save.file.name = paste("Realdata_lowbnds_", concent_name[concentration_idx], "_group", group_idx,"_Paramset", param_num, "_maxT", max.T, "_nsample", nsample, " repnum", repnum, sep ="")
  
  param.name.list = c("AX", "KM", "alphaX", "betaX", "AY", "alphaY", "betaY", "B")
  png(paste(save.file.name,"summary.png", sep = ""), width = 800, height = 600)
  par(mfrow = c(3, 3))
  par(mar = c(5.1,4.1,4.1,2.1))
  for (param.id in 1:2){
    plot(theta.est[,param.id], type = "l" , col = rgb(0/255,0/255,0/255, 0.6),
         xlab = "", ylab = "", bty = "n")
    # xlim = c(0,15), ylim = c(0,15)
    axis(side = 1, lwd = 2)
    axis(side = 2, lwd = 2)
    mtext(side=1, line=2, "Iteration", font=1.6,cex=1.2)
    mtext(side=2, line=2, param.name.list[param.id], font=1.6,cex=1.2)
  }
  plot(theta.est[,1]/theta.est[,2], type = "l" , col = rgb(0/255,0/255,0/255, 0.6),
       xlab = "", ylab = "", bty = "n")
  # xlim = c(0,15), ylim = c(0,15)
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  mtext(side=1, line=2, "Iteration", font=1.6,cex=1.2)
  mtext(side=2, line=2, "AX/KM", font=1.6,cex=1.2)
  
  for (param.id in 3:4){
    plot(theta.est[,param.id], type = "l" , col = rgb(0/255,0/255,0/255, 0.6),
         xlab = "", ylab = "", bty = "n")
    # xlim = c(0,15), ylim = c(0,15)
    axis(side = 1, lwd = 2)
    axis(side = 2, lwd = 2)
    mtext(side=1, line=2, "Iteration", font=1.6,cex=1.2)
    mtext(side=2, line=2, param.name.list[param.id], font=1.6,cex=1.2)
  }
  plot(theta.est[,3]/theta.est[,4], type = "l" , col = rgb(0/255,0/255,0/255, 0.6),
       xlab = "", ylab = "", bty = "n")
  # xlim = c(0,15), ylim = c(0,15)
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  mtext(side=1, line=2, "Iteration", font=1.6,cex=1.2)
  mtext(side=2, line=2, "X delay", font=1.6,cex=1.2)
  
  for (param.id in 5:7){
    plot(theta.est[,param.id], type = "l" , col = rgb(0/255,0/255,0/255, 0.6),
         xlab = "", ylab = "", bty = "n")
    # xlim = c(0,15), ylim = c(0,15)
    axis(side = 1, lwd = 2)
    axis(side = 2, lwd = 2)
    mtext(side=1, line=2, "Iteration", font=1.6,cex=1.2)
    mtext(side=2, line=2, param.name.list[param.id], font=1.6,cex=1.2)
  }
  mtext(text= title_text, side = 3, line = -2, outer = TRUE, cex = 1.5)
  
  dev.off()
  
  
  png(paste(save.file.name,"_fitting.png", sep = ""), width = 800, height = 600)
  par(mfrow = c(2,2))
  par(mar = c(5.1,4.1,4.1,2.1))
  
  plot(theta.est[,3]/theta.est[,4],theta.est[,6]/theta.est[,7], col = rgb(0/255,0/255,0/255, 0.6),
       xlab = "", ylab = "", bty = "n")
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  lines(c(0,12), c(12,0), lwd = 2, col = rgb(216/255,75/255,36/255, 0.8))
  mtext(side=1, line=2, "X delay", font=1.6,cex=1.2)
  mtext(side=2, line=2, "Y delay", font=1.6,cex=1.2)
  
  # install.packages("dplyr")
  # library(dplyr)
  
  # Xtrj.true = rbind(rep(0, nsample),apply(Xbirth.true - Xdeath.true, 2, cumsum))
  Ytrj.true = t(Y.all)
  
  Xtrj.est.mat = c()
  Ytrj.est.mat = c()
  for (jj in 1:10){
    simul0 = TimeDelayGillespieforXY(A.X = colMeans(theta.est)[1], B.X =  colMeans(theta.est)[8],
                                     alpha.X = colMeans(theta.est)[3], beta.X = colMeans(theta.est)[4],
                                     A.Y = colMeans(theta.est)[5], K.M = colMeans(theta.est)[2],
                                     alpha.Y = colMeans(theta.est)[6], beta.Y = colMeans(theta.est)[7],
                                     repnum = 10000000, maxT = max.T, Volume = 1, B.Y = colMeans(theta.est)[8])
    
    Xtrj.est.mat = rbind(Xtrj.est.mat, c(0, cumsum(simul0$Xbirth[1:max.T] - simul0$Xdeath[1:max.T])))
    Ytrj.est.mat = rbind(Ytrj.est.mat, c(0, cumsum(simul0$Ybirth[1:max.T] - simul0$Ydeath[1:max.T])))
  }
  Xtrj.est = colMeans(Xtrj.est.mat)
  Ytrj.est = colMeans(Ytrj.est.mat)
  
  
  
  dim(theta.est)[1]
  selrow.ind = seq(from = dim(theta.est)[1]/25, by = dim(theta.est)[1]/25, to = dim(theta.est)[1])
  Xtrj.est.ind.mat = c()
  Ytrj.est.ind.mat = c()
  for (jj in 1:25){
    simul.ind = TimeDelayGillespieforXY(A.X = theta.est[selrow.ind[jj],1], B.X =  theta.est[selrow.ind[jj],8],
                                        alpha.X = theta.est[selrow.ind[jj],3], beta.X = theta.est[selrow.ind[jj],4],
                                        A.Y = theta.est[selrow.ind[jj],5], K.M = theta.est[selrow.ind[jj],2],
                                        alpha.Y = theta.est[selrow.ind[jj],6], beta.Y = theta.est[selrow.ind[jj],7],
                                        repnum = 10000000, maxT = max.T, Volume = 1, B.Y = theta.est[selrow.ind[jj],8])
    
    Xtrj.est.ind.mat = rbind(Xtrj.est.ind.mat, c(0, cumsum(simul.ind$Xbirth[1:max.T] - simul.ind$Xdeath[1:max.T])))
    Ytrj.est.ind.mat = rbind(Ytrj.est.ind.mat, c(0, cumsum(simul.ind$Ybirth[1:max.T] - simul.ind$Ydeath[1:max.T])))
  }
  

  for (cellid in 1:25){
    if (cellid == 1){
      plot(0:max.T, Xtrj.est.ind.mat[cellid, ]/theta.est[selrow.ind[cellid],2], xlab ="", ylab = "", col = rgb(250/255,157/255,0/255, 0.4), type = "l")
    }else{
      lines(0:max.T, Xtrj.est.ind.mat[cellid, ]/theta.est[selrow.ind[cellid],2], xlab ="", ylab = "", col = rgb(250/255,157/255,0/255, 0.4))
    }
  }
  
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  mtext(side=1, line=2, "Time", font=1.6,cex=1.2)
  mtext(side=2, line=2, "X(t)/KM", font=1.6,cex=1.2)
  
  
  for (cellid in 1:min(16, nsample)){
    if (cellid == 1){
      plot(0:max.T, Ytrj.true[cellid,], xlab ="", ylab = "", bty = "n", type = "l", ylim = c(0,1.3*max(Ytrj.true[cellid,])))
    }else{
      lines(0:max.T, Ytrj.true[cellid,], xlab ="", ylab = "", col = rgb(35/255,35/255,35/255, 0.8))
    }
  }
  
  for (cellid in 1:25){
    lines(0:max.T, Ytrj.est.ind.mat[cellid, ], xlab ="", ylab = "", col = rgb(250/255,157/255,0/255, 0.4))
  }
  
  lines(0:max.T, Ytrj.est, xlab ="", ylab = "", lwd = 2, col = "red")
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  mtext(side=1, line=2, "Time", font=1.6,cex=1.2)
  mtext(side=2, line=2, "Y(t)", font=1.6,cex=1.2)
  mtext(text= title_text, side = 3, line = -2, outer = TRUE, cex = 1.5)
  dev.off()
  
}


# par(mar = c(5.1,4.1,4.1,2.1))






####### FIGURE 1 - Misspecification -- ONE step model --------

setwd("/Users/hyukpyohong/OneDrive - kaist.ac.kr/Research/ResearchMaterial_HHP/TimeDelayEstimation/CodeForSobolev/")
true.twostep.param = c(10, 200, 3.6, 0.6, 60, 3.6, 0.6, 0.05)
AX.list = c(2, 5, 10, 20, 50)
AX.list.est= rep(NA, 5)
B.list.est= rep(NA, 5)
alpha.list.est= rep(NA, 5)
beta.list.est= rep(NA, 5)
delay.list.est= rep(NA, 5)


max.T = 150


for(ii in 1:5){
  # load(paste("UsingOnestep_fixBTF1_AX_",AX.list[ii],".RData", sep = ""))
  load(paste("UsingOnestep_fixBTF0_AX_",AX.list[ii],".RData", sep = ""))
  save.file.name = paste("UsingOnestep_fixBTF0_AX_",AX.list[ii], sep = "")
  title_text = paste("Onestep fix B AX:",AX.list[ii],sep ="")
  # title_text = paste("Onestep not fix B AX:",AX.list[ii],sep ="")
  
  
  png(paste(save.file.name,"_summary.png", sep = ""), width = 800, height = 600)
  
  par(mfrow = c(2, 3))
  par(mar = c(5.1,4.1,4.1,2.1))
  selrow = seq(from = 5010, by = 10, to = 55000)
  param.name.list = c("AX", "B", "alphaX", "betaX")
  
  theta.est = short1[selrow, ]
  
  for (param.id in 1:2){
    plot(theta.est[,param.id], type = "l" , col = rgb(0/255,0/255,0/255, 0.6),
         xlab = "", ylab = "", bty = "n")
    # xlim = c(0,15), ylim = c(0,15)
    axis(side = 1, lwd = 2)
    axis(side = 2, lwd = 2)
    mtext(side=1, line=2, "Iteration", font=1.6,cex=1.2)
    mtext(side=2, line=2, param.name.list[param.id], font=1.6,cex=1.2)
  }
  
  plot(theta.est[,1]/theta.est[,2], type = "l" , col = rgb(0/255,0/255,0/255, 0.6),
       xlab = "", ylab = "", bty = "n")
  # xlim = c(0,15), ylim = c(0,15)
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  mtext(side=1, line=2, "Iteration", font=1.6,cex=1.2)
  mtext(side=2, line=2, "AX/B", font=1.6,cex=1.2)
  
  for (param.id in 3:4){
    plot(theta.est[,param.id], type = "l" , col = rgb(0/255,0/255,0/255, 0.6),
         xlab = "", ylab = "", bty = "n")
    # xlim = c(0,15), ylim = c(0,15)
    axis(side = 1, lwd = 2)
    axis(side = 2, lwd = 2)
    mtext(side=1, line=2, "Iteration", font=1.6,cex=1.2)
    mtext(side=2, line=2, param.name.list[param.id], font=1.6,cex=1.2)
  }
  
  plot(theta.est[,3]/theta.est[,4], type = "l" , col = rgb(0/255,0/255,0/255, 0.6),
       xlab = "", ylab = "", bty = "n")
  # xlim = c(0,15), ylim = c(0,15)
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  mtext(side=1, line=2, "Iteration", font=1.6,cex=1.2)
  mtext(side=2, line=2, "Delay", font=1.6,cex=1.2)
  mtext(text= title_text, side = 3, line = -2, outer = TRUE, cex = 1.5)
  
  dev.off()
  
  
  AX.list.est[ii] = mean(theta.est[, 1])
  B.list.est[ii] = mean(theta.est[, 2])
  alpha.list.est[ii] = mean(theta.est[, 3])
  beta.list.est[ii] = mean(theta.est[, 4])
  delay.list.est[ii] = mean(theta.est[, 3]/theta.est[, 4])
  
  
  
  par(mfrow = c(2, 3))
  par(mar = c(5.1,4.1,4.1,2.1))
  png(paste(save.file.name,"_fitting.png", sep = ""), width = 800, height = 600)
  
  
  
  selrow.ind = seq(from = dim(theta.est)[1]/25, by = dim(theta.est)[1]/25, to = dim(theta.est)[1])
  Xtrj.est.ind.mat = c()
  
  Ytrj.est.ind.mat = c()
  for (jj in 1:25){
    simul.ind.X = TimeDelayGillespieforXR(A.X = theta.est[selrow.ind[jj],1], B.X =  theta.est[selrow.ind[jj],2],
                                          alpha.X = theta.est[selrow.ind[jj],3], beta.X = theta.est[selrow.ind[jj],4],
                                          maxT = 150, repnum = 1000000, Volume = 1)
    
    Xtrj.est.ind.mat = rbind(Xtrj.est.ind.mat, c(0, cumsum(simul.ind.X$Xbirth[1:max.T] - simul.ind.X$Xdeath[1:max.T])))
    
    simul.ind.Y = TimeDelayGillespieforXY(A.X = AX.list[ii], B.X =  true.twostep.param[8],
                                          alpha.X = true.twostep.param[3], beta.X = true.twostep.param[4],
                                          A.Y = true.twostep.param[5], K.M = true.twostep.param[2],
                                          alpha.Y = true.twostep.param[6], beta.Y = true.twostep.param[7],
                                          repnum = 10000000, maxT = 150, Volume = 1, B.Y = true.twostep.param[8])
    
    Ytrj.est.ind.mat = rbind(Ytrj.est.ind.mat, c(0, cumsum(simul.ind.Y$Ybirth[1:max.T] - simul.ind.Y$Ydeath[1:max.T])))
  }
  
  Xtrj.est.mat = c()
  
  for (jj in 1:10){
    simul0 = TimeDelayGillespieforXR(A.X = colMeans(theta.est)[1], B.X =  colMeans(theta.est)[2],
                                    alpha.X = colMeans(theta.est)[3], beta.X = colMeans(theta.est)[4],
                                    repnum = 10000000, maxT = max.T, Volume = 1)
    
    Xtrj.est.mat = rbind(Xtrj.est.mat, c(0, cumsum(simul0$Xbirth[1:max.T] - simul0$Xdeath[1:max.T])))
  }
  Xtrj.est = colMeans(Xtrj.est.mat)
  
  # par(mfrow = c(1,1))
  # par(mar = c(1,3,3,1))
  # par(mar = c(5.1,4.1,4.1,2.1))
  for (cellid in 1:25){
    if (cellid == 1){
      plot(0:max.T, Ytrj.est.ind.mat[cellid,], xlab ="", ylab = "", bty = "n", type = "l")
    }else{
      lines(0:max.T, Ytrj.est.ind.mat[cellid,], xlab ="", ylab = "", col = rgb(35/255,35/255,35/255, 0.8))
    }
  }
  
  for (cellid in 1:25){
    lines(0:max.T, Xtrj.est.ind.mat[cellid, ], xlab ="", ylab = "", col = rgb(250/255,175/255,0/255, 0.4))
  }
  
  lines(0:max.T, Xtrj.est, xlab ="", ylab = "", lwd = 2, col = "red")
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  mtext(side=1, line=2, "Time", font=1.6,cex=1.2)
  mtext(side=2, line=2, "Y(t)", font=1.6,cex=1.2)
  dev.off()
}

par(mfrow = c(1,1))
axis.lim.val = max(c(AX.list/AX.list[1], AX.list.est/AX.list.est[1]))
plot(AX.list/AX.list[1], AX.list.est/AX.list.est[1], xlim = c(0, axis.lim.val), ylim = c(0, axis.lim.val), cex = 2, xlab = "", ylab = "", bty = "n")
lines(c(0, axis.lim.val),c(0, axis.lim.val), col = "red")
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
mtext(side=1, line=2, "Norm. two step AX", font=1.6,cex=1.2)
mtext(side=2, line=2, "Norm. one step estimated AX", font=1.6,cex=1.2)





par(mfrow = c(2, 3))
par(mar = c(5.1,4.1,4.1,2.1))
png(paste("UsingOnestep_fixBTF0_fitting.png", sep = ""), width = 800, height = 600)
for(ii in 1:5){
  # load(paste("UsingOnestep_fixBTF1_AX_",AX.list[ii],".RData", sep = ""))
  load(paste("UsingOnestep_fixBTF0_AX_",AX.list[ii],".RData", sep = ""))
  save.file.name = paste("UsingOnestep_fixBTF0_AX_",AX.list[ii], sep = "")
  # title_text = paste("Onestep not fix B AX:",AX.list[ii],sep ="")
  
  
  AX.list.est[ii] = mean(theta.est[, 1])
  B.list.est[ii] = mean(theta.est[, 2])
  alpha.list.est[ii] = mean(theta.est[, 3])
  beta.list.est[ii] = mean(theta.est[, 4])
  delay.list.est[ii] = mean(theta.est[, 3]/theta.est[, 4])
  selrow = seq(from = 5010, by = 10, to = 55000)
  theta.est = short1[selrow, ]
  
  selrow.ind = seq(from = dim(theta.est)[1]/25, by = dim(theta.est)[1]/25, to = dim(theta.est)[1])
  Xtrj.est.ind.mat = c()
  Ytrj.est.ind.mat = c()
  for (jj in 1:25){
    simul.ind.X = TimeDelayGillespieforXR(A.X = theta.est[selrow.ind[jj],1], B.X =  theta.est[selrow.ind[jj],2],
                                          alpha.X = theta.est[selrow.ind[jj],3], beta.X = theta.est[selrow.ind[jj],4],
                                          maxT = 150, repnum = 1000000, Volume = 1)
    
    Xtrj.est.ind.mat = rbind(Xtrj.est.ind.mat, c(0, cumsum(simul.ind.X$Xbirth[1:max.T] - simul.ind.X$Xdeath[1:max.T])))
    
    simul.ind.Y = TimeDelayGillespieforXY(A.X = AX.list[ii], B.X =  true.twostep.param[8],
                                          alpha.X = true.twostep.param[3], beta.X = true.twostep.param[4],
                                          A.Y = true.twostep.param[5], K.M = true.twostep.param[2],
                                          alpha.Y = true.twostep.param[6], beta.Y = true.twostep.param[7],
                                          repnum = 10000000, maxT = 150, Volume = 1, B.Y = true.twostep.param[8])
    
    Ytrj.est.ind.mat = rbind(Ytrj.est.ind.mat, c(0, cumsum(simul.ind.Y$Ybirth[1:max.T] - simul.ind.Y$Ydeath[1:max.T])))
  }
  
  Xtrj.est.mat = c()
  
  for (jj in 1:10){
    simul0 = TimeDelayGillespieforXR(A.X = colMeans(theta.est)[1], B.X =  colMeans(theta.est)[2],
                                     alpha.X = colMeans(theta.est)[3], beta.X = colMeans(theta.est)[4],
                                     repnum = 10000000, maxT = max.T, Volume = 1)
    
    Xtrj.est.mat = rbind(Xtrj.est.mat, c(0, cumsum(simul0$Xbirth[1:max.T] - simul0$Xdeath[1:max.T])))
  }
  Xtrj.est = colMeans(Xtrj.est.mat)
  
  # par(mfrow = c(1,1))
  # par(mar = c(1,3,3,1))
  # par(mar = c(5.1,4.1,4.1,2.1))
  for (cellid in 1:25){
    if (cellid == 1){
      plot(0:max.T, Ytrj.est.ind.mat[cellid,], xlab ="", ylab = "", bty = "n", type = "l")
    }else{
      lines(0:max.T, Ytrj.est.ind.mat[cellid,], xlab ="", ylab = "", col = rgb(35/255,35/255,35/255, 0.8))
    }
  }
  
  for (cellid in 1:25){
    lines(0:max.T, Xtrj.est.ind.mat[cellid, ], xlab ="", ylab = "", col = rgb(250/255,175/255,0/255, 0.4))
  }
  
  lines(0:max.T, Xtrj.est, xlab ="", ylab = "", lwd = 2, col = "red")
  axis(side = 1, lwd = 2)
  axis(side = 2, lwd = 2)
  mtext(side=1, line=2, "Time", font=1.6,cex=1.2)
  mtext(side=2, line=2, "Y(t)", font=1.6,cex=1.2)
}

dev.off()


