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
import("scipy")
pd <- import("pandas")
np <- import("numpy")
source_python("BFuns_for_Bio.py")
#source_python("C:/Users/qihan/Desktop/Geo/BFuns_co21_in_R.py")

full_data = read_excel("C:/Users/qihan/Desktop/CTG_Qihan.xls", sheet = 2)
#useful_vars = c("LB", "AC...11", "FM...12", "UC...13", "DL...14", "DS...15", "DP...16", "ASTV",
#                "MSTV", "ALTV", "MLTV", "Width","Min", "Max", "Nmax", "Nzeros", "Mode", "Mean",
#                "Median", "Variance", "Tendency", "CLASS", "NSP")


useful_vars = c("LB","AC...11", "FM...12", "UC...13", "DP...16", "ASTV",
                "MSTV", "ALTV", "MLTV", "Min", "Mean", "Mode", "Median", "CLASS", "NSP")


reduce_data = full_data[useful_vars]
reduce_data$inter <- rep(1, dim(reduce_data)[1])





set.seed(2024) 
split <- sample.split(reduce_data$LB, SplitRatio = 0.7) 
train_data <- subset(reduce_data, split == "TRUE") 
test_data <- subset(reduce_data, split == "FALSE")


Xtrain = train_data[, 1:(length(useful_vars)-2)]
Xtest = test_data[, 1:(length(useful_vars)-2)]


Ytrain = train_data$NSP
Ytest = test_data$NSP

#Ytrain = train_data$CLASS
#Ytest = test_data$CLASS


dist1 = as.matrix(dist(Xtrain))
dist2 = as.matrix(pdist(Xtest, Xtrain))


mds <- cmdscale(dist1, k = 2)
x <- mds[,1]
y <- mds[,2]


Xtrain = as.matrix(Xtrain)
Xtest = as.matrix(Xtest)
Ytrain = as.matrix(Ytrain)
Ytest = as.matrix(Ytest)



data1 = cbind(x, y, Ytrain)
geodata1 = as.geodata(data1, coords.col = 1:2, data.col = 3)
vario = variog(geodata1, option="bin", max.dist=50, trend = ~1)
plot(vario,xlab = "h", ylab = "gamma(h)")
intvalues = eyefit(vario, silent = FALSE)
int_sigmasq = intvalues[[1]][2][[1]][1]
int_phi = intvalues[[1]][2][[1]][2]
#res1 = MLE_fit(y = Ytrain, X = as.matrix(train_data$inter), D = dist1, cov_model = "Exp", opt = "NM", i_sigma = int_sigmasq, i_phi = int_phi)



res1 = MLE_fit(y = Ytrain, X = as.matrix(train_data$inter), D = dist1, cov_model = "Exp")

thetahat1 = res1[[1]] 
betahat1  = res1[[2]] 
Sigma1 = thetahat1[2]*exp(-dist1/thetahat1[1]) 
diag(Sigma1) = thetahat1[2] + thetahat1[3]

cmat1 = thetahat1[2]*exp(-dist2/thetahat1[1]) 
pl1 = as.matrix(test_data$inter)%*%betahat1
pe1 = cmat1%*% solve(Sigma1)%*%(Ytrain - as.matrix(train_data$inter)%*%betahat1) 


Ypred = pl1 + pe1

#Ypred1 = round(Ypred1)
#idx1 = Ypred1>3
#Ypred1[idx1]=3
#idx2 = Ypred1<1
#Ypred1[idx2]=1

Kset1 = c(1.454, 1.491, 1.519, 1.386, 1.8, 1.852, 1.6, 1.672, 1.778, 1.843)
Kset2 = c(2.494, 2.415, 2.91, 2.89, 2.341, 2.865, 2.263, 2.122, 2.407, 2.405)

AA1 = c()
for (i in 1:length(Kset1)) {
  K1 = Kset1[i]
  K2 = Kset2[i]
  Ypred1 = Ypred
  
  idx1 = Ypred1 < K1
  Ypred1[idx1]=1
  idx2 = Ypred1 > K2
  Ypred1[idx2]=3
  idx3 = (K1 <= Ypred1) & (Ypred1 <= K2)
  Ypred1[idx3]=2
  
  misClassError <- mean(Ypred1 != Ytest) 
  Accuracy =1-misClassError
  AA1[i] = Accuracy
  
}

K1 = mean(Kset1)
K2 = mean(Kset2)

Ypred1 = Ypred

idx1 = Ypred1 < K1
Ypred1[idx1]=1
idx2 = Ypred1 > K2
Ypred1[idx2]=3
idx3 = (K1 <= Ypred1) & (Ypred1 <= K2)
Ypred1[idx3]=2

misClassError <- mean(Ypred1 != Ytest) 
print(paste('Accuracy =', 1-misClassError))
(cmatrix = table(Ytest, Ypred1))

plot(Ytest, pl1 + pe1, xlab = 'Y true', ylab = 'Y Predicted')
plot(pl1 + pe1, xlab = 'Index', ylab = 'Y Predicted')
plot(Ypred1, xlab = 'Index', ylab = 'projection Y Predicted')
plot(Ytest, xlab = 'Index', ylab = 'Y true')

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
    Ypred1 = Ypred
    idx1 = Ypred1 < K1
    Ypred1[idx1]=1
    idx2 = Ypred1 > K2
    Ypred1[idx2]=3
    idx3 = (K1 <= Ypred1) & (Ypred1 <= K2)
    Ypred1[idx3]=2
    misClassError <- mean(Ypred1 != Ytest) 
    Accuracy = 1-misClassError
    
    if (Accuracy > max_A){
      max_A = Accuracy
      max_K1 = K1
      max_K2 = K2
    }
  }
}


