library(pdist)
library(readxl)
library(FNN)
library(caTools) 
library(class)
library(gstat)  
library(sp)
library(geoR)
library("reticulate")
library(reticulate)
library(dplyr)
require("reticulate")
library("corrplot")
library(factoextra)
import("scipy")
pd <- import("pandas")
np <- import("numpy")
source_python("BFuns_for_Bio.py")

full_data = read_excel("C:/Users/qihan/Desktop/CTG_Qihan.xls", sheet = 2)
useful_vars = c("LB", "AC...11", "FM...12", "UC...13", "DL...14", "DS...15", "DP...16", "ASTV",
                "MSTV", "ALTV", "MLTV", "Width","Min", "Max", "Nmax", "Nzeros", "Mode", "Mean",
                "Median", "Variance", "Tendency", "CLASS", "NSP")

useful_vars = c("LB","AC...11", "FM...12", "UC...13", "DP...16", "ASTV",
                "MSTV", "ALTV", "MLTV", "Min", "Mean", "Mode", "Median", "CLASS", "NSP")

reduce_data = full_data[useful_vars]
#colnames(reduce_data)[2] <- "AC"
#colnames(reduce_data)[3] <- "FM"
#colnames(reduce_data)[4] <- "UC"
#colnames(reduce_data)[5] <- "DL"
#colnames(reduce_data)[6] <- "DS"
#colnames(reduce_data)[7] <- "DP"
##################################################################################################
# Splitting data into train and test data 
set.seed(2024) 
split <- sample.split(reduce_data$LB, SplitRatio = 0.7) 
train_data <- subset(reduce_data, split == "TRUE") 
test_data <- subset(reduce_data, split == "FALSE")


x.comb = data.frame(train_data[, 1:(length(useful_vars)-2)])
x.test = data.frame(test_data[, 1:(length(useful_vars)-2)])
y.comb = data.frame(train_data$NSP)
y.test = data.frame(test_data$NSP)


cv.k <- 20
n.comb <- dim(y.comb)[1]
if (n.comb %% cv.k == 0) {
  split.ind <- rep(1:cv.k, each = floor(n.comb/cv.k))
} else {
  split.ind <- c(rep(1:cv.k, each = floor(n.comb/cv.k)),
                 1:(n.comb %% cv.k))
}
split.ind <- sample(split.ind) 



A_set = c()
K1_set = c()
K2_set = c()
for (j in 1:cv.k) {
  print(j)
  x.tune <- data.frame(x.comb[split.ind == j,])
  y.tune <- data.frame(y.comb[split.ind == j,])
  x.train <- data.frame(x.comb[split.ind != j,])
  y.train <- data.frame(y.comb[split.ind != j,])
  dist_tt1 = as.matrix(dist(x.train))
  dist_tt2 = as.matrix(pdist(x.tune, x.train))
  mds <- cmdscale(dist_tt1, k = 2)
  x_tt <- mds[,1]
  y_tt <- mds[,2]
  X.train = as.matrix(x.train)
  X.tune = as.matrix(x.tune)
  Y.train = as.matrix(y.train)
  Y.tune = as.matrix(y.tune)
  
  res_tt = MLE_fit(y = Y.train, X = as.matrix(rep(1, dim(Y.train)[1])), D = dist_tt1, cov_model = "Exp")
  thetahat_tt = res_tt[[1]] 
  betahat_tt  = res_tt[[2]] 
  Sigma_tt = thetahat_tt[2]*exp(-dist_tt1/thetahat_tt[1]) 
  diag(Sigma_tt) = thetahat_tt[2] + thetahat_tt[3]
  
  cmat_tt = thetahat_tt[2]*exp(-dist_tt2/thetahat_tt[1]) 
  pl_tt = as.matrix(rep(1, dim(Y.tune)[1]))%*%betahat_tt
  pe_tt = cmat_tt%*% solve(Sigma_tt)%*%(Y.train - as.matrix(rep(1, dim(Y.train)[1]))%*%betahat_tt) 
  Ypred_tt = pl_tt + pe_tt
  
  N = 1000
  M = 1000
  max_A = 0
  max_K1 = 0
  max_K2 = 0
  idx=0
  for (a in 1:N) {
    K1 = 1+0.001*(a-1)
    for (b in 1:M) {
      idx = idx+1
      K2 = 2+0.001*(b-1)
      Ypred1 = Ypred_tt
      idx1 = Ypred1 < K1
      Ypred1[idx1]=1
      idx2 = Ypred1 > K2
      Ypred1[idx2]=3
      idx3 = (K1 <= Ypred1) & (Ypred1 <= K2)
      Ypred1[idx3]=2
      misClassError <- mean(Ypred1 != y.tune[,1, drop = TRUE]) 
      Accuracy = 1-misClassError
      
      if (Accuracy >= max_A){
        max_A = Accuracy
        max_K1 = K1
        max_K2 = K2
      }
    }
  }
  mylist = c(max_A, max_K1, max_K2)
  A_set[j] = max_A
  K1_set[j] = max_K1
  K2_set[j] = max_K2
}

total = cbind(K1_set, K2_set, A_set)




