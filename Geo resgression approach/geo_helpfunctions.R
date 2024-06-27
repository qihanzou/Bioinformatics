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



full_data = read_excel("C:/Users/qihan/Desktop/CTG_Qihan.xls", sheet = 2)
#useful_vars = c("LB", "AC...11", "FM...12", "UC...13", "DL...14", "DS...15", "DP...16", "ASTV",
#                "MSTV", "ALTV", "MLTV", "Width","Min", "Max", "Nmax", "Nzeros", "Mode", "Mean",
#                "Median", "Variance", "Tendency", "CLASS", "NSP")
useful_vars = c("LB","AC...11", "FM...12", "UC...13", "DP...16", "ASTV",
                "MSTV", "ALTV", "MLTV", "Min", "Mean", "Mode", "Median", "CLASS", "NSP")

reduce_data = full_data[useful_vars]
#colnames(reduce_data)[2] <- "AC"
#colnames(reduce_data)[3] <- "FM"
##colnames(reduce_data)[4] <- "UC"
#colnames(reduce_data)[5] <- "DL"
#colnames(reduce_data)[6] <- "DS"
#colnames(reduce_data)[7] <- "DP"

M = cor(reduce_data)
corrplot(M, method="color")


#useful_vars = c("LB", "AC", "FM", "UC", "DL", "DS", "DP", "ASTV",
#                "MSTV", "ALTV", "MLTV", "Width","Min", "Max", "Nmax", "Nzeros", "Mean",
#                "Variance", "Tendency", "CLASS", "NSP")
#reduce_data = reduce_data0[useful_vars]

M = cor(reduce_data)
corrplot(M, method="color")

set.seed(2024) 
split <- sample.split(reduce_data$LB, SplitRatio = 0.7) 
train_data <- subset(reduce_data, split == "TRUE") 
test_data <- subset(reduce_data, split == "FALSE")


Xtrain = train_data[, 1:13]
Xtest = test_data[, 1:13]
#Ytrain = train_data$NSP
#Ytest = test_data$NSP
Ytrain = train_data$NSP
Ytest = test_data$NSP

#Xtrain = Xtrain[-c(241, 794, 243, 168),]
#Ytrain = Ytrain[-c(241, 794, 243, 168)]


#dist1 = as.matrix(pdist(Xtrain, Xtrain))
dist1 = as.matrix(dist(Xtrain))
dist2 = as.matrix(pdist(Xtest, Xtrain))




mds <- cmdscale(dist1, k = 2)
x <- mds[,1]
y <- mds[,2]

#data1 = cbind(x, y, Ytrain, Xtrain)
#geodata = as.geodata(data1, coords.col = 1:2, data.col = 3, covar.col = 4:24, 
#                                          covar.names = c("LB", "AC", "FM", "UC", "DL", "DS", "DP", "ASTV",
#                                                          "MSTV", "ALTV", "MLTV", "Width","Min", "Max", "Nmax", "Nzeros", "Mode", "Mean",
#                                                          "Median", "Variance", "Tendency")) 


#duplicated(geodata, incomparables)






#vario = variog(geodata, option="bin", max.dist=500, trend = ~LB+AC+FM+UC+DL+DS+DP+ASTV+MSTV+ALTV+MLTV+Width+Min+Max+Nmax+Nzeros+Mode+Mean+Median+Variance+Tendency)
#plot(vario,xlab = "h", ylab = "gamma(h)")
#intvalues = eyefit(vario, silent = FALSE)

#cov.model sigmasq    phi tausq kappa kappa2 practicalRange
#1 exponential     6.8 104.09  0.85    NA     NA       311.8258


Xtrain = as.matrix(Xtrain)
Xtest = as.matrix(Xtest)
Ytrain = as.matrix(Ytrain)
Ytest = as.matrix(Ytest)






res1 = MLE_fit(y = Ytrain, X = Xtrain, D = dist1, cov_model = "Exp")
res1


thetahat1 = res1[[1]] 
betahat1  = res1[[2]] 
Sigma1 = thetahat1[2]*exp(-dist1/thetahat1[1]) 
diag(Sigma1) = thetahat1[2] + thetahat1[3]

cmat1 = thetahat1[2]*exp(-dist2/thetahat1[1]) 
pl1 = Xtest%*%betahat1
pe1 = cmat1%*% solve(Sigma1)%*%(Ytrain - Xtrain%*%betahat1) 
Ypred1 = pl1 + pe1

Ypred1 = round(Ypred1)
idx1 = Ypred1>3
Ypred1[idx1]=3
idx2 = Ypred1<1
Ypred1[idx2]=1

#idx1 = Ypred1>10
#Ypred1[idx1]=10
#idx2 = Ypred1<1
#Ypred1[idx2]=1

misClassError <- mean(Ypred1 != Ytest) 
print(paste('Accuracy =', 1-misClassError))




(cmatrix = table(Ytest, Ypred1))


plot(Ytest, pl1 + pe1, xlab = 'Y true', ylab = 'Y Predicted')


plot(pl1 + pe1, xlab = 'Index', ylab = 'Y Predicted')

plot(Ypred1, xlab = 'Index', ylab = 'projection Y Predicted')

plot(Ytest, xlab = 'Index', ylab = 'Y true')




