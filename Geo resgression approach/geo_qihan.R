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



full_data = read_excel("C:/Users/qihan/Desktop/CTG_Qihan.xls", sheet = 2)
useful_vars = c("LB", "AC...11", "FM...12", "UC...13", "DL...14", "DS...15", "DP...16", "ASTV",
                "MSTV", "ALTV", "MLTV", "Width","Min", "Max", "Nmax", "Nzeros", "Mode", "Mean",
                "Median", "Variance", "Tendency", "CLASS", "NSP")
reduce_data0 = full_data[useful_vars]
colnames(reduce_data0)[2] <- "AC"
colnames(reduce_data0)[3] <- "FM"
colnames(reduce_data0)[4] <- "UC"
colnames(reduce_data0)[5] <- "DL"
colnames(reduce_data0)[6] <- "DS"
colnames(reduce_data0)[7] <- "DP"


useful_vars = c("LB", "AC", "FM", "UC", "DL", "DS", "DP", "ASTV",
                "MSTV", "ALTV", "MLTV", "Width","Min", "Max", "Nmax", "Nzeros", "Mean",
                "Variance", "Tendency", "CLASS", "NSP")
reduce_data = reduce_data0[useful_vars]

M = cor(reduce_data)
corrplot(M, method="color")




set.seed(2024) 
split <- sample.split(reduce_data$LB, SplitRatio = 0.7) 
train_data <- subset(reduce_data, split == "TRUE") 
test_data <- subset(reduce_data, split == "FALSE")


Xtrain = train_data[, 1:19]
Xtest = test_data[, 1:19]
Ytrain = train_data$NSP
Ytest = test_data$NSP


dist1 = as.matrix(pdist(Xtrain, Xtrain))
dist2 = as.matrix(pdist(Xtest, Xtrain))
mds <- cmdscale(dist1, k = 2)
x <- mds[,1]
y <- mds[,2]



data1 = cbind(x, y, Ytrain, Xtrain)
geodata = as.geodata(data1, coords.col = 1:2, data.col = 3, covar.col = 4:22, 
                                          covar.names = c("LB", "AC", "FM", "UC", "DL", "DS", "DP", "ASTV",
                                                          "MSTV", "ALTV", "MLTV", "Width","Min", "Max", "Nmax", "Nzeros", "Mean",
                                                          "Variance", "Tendency")) 


duplicated(geodata, incomparables)







Xtrain[241,]
Xtrain[242,]
Ytrain[241]
Ytrain[242]

Xtrain[794,]
Xtrain[795,]
Ytrain[794]
Ytrain[795]

Xtrain[243,]
Xtrain[249,]
Ytrain[243]
Ytrain[249]


Xtrain[168,]
Xtrain[573,] # after reduce, index become 571

r_Xtrain = Xtrain[-c(241, 794, 243, 168),]
r_Ytrain = Ytrain[-c(241, 794, 243, 168)]

r_dist1 = as.matrix(pdist(r_Xtrain, r_Xtrain))
r_dist2 = as.matrix(pdist(Xtest, r_Xtrain))
r_mds <- cmdscale(r_dist1, k = 2)
r_x <- r_mds[,1]
r_y <- r_mds[,2]


r_data = cbind(r_x, r_y, r_Ytrain, r_Xtrain)
r_geodata = as.geodata(r_data, coords.col = 1:2, data.col = 3, covar.col = 4:22, 
                     covar.names = c("LB", "AC", "FM", "UC", "DL", "DS", "DP", "ASTV",
                                     "MSTV", "ALTV", "MLTV", "Width","Min", "Max", "Nmax", "Nzeros", "Mean",
                                     "Variance", "Tendency")) 


duplicated(r_geodata, incomparables)


vario = variog(r_geodata, option="bin", max.dist=500, trend = ~1+LB+AC+FM+UC+DL+DS+DP+ASTV+MSTV+ALTV+MLTV+Width+Min+Max+Nmax+Nzeros+Mean+Variance+Tendency)
plot(vario,xlab = "h", ylab = "gamma(h)")
#intvalues = eyefit(vario, silent = FALSE)


Xtrain = as.matrix(r_Xtrain)
Xtest = as.matrix(Xtest)
Ytrain = as.matrix(r_Ytrain)
Ytest = as.matrix(Ytest)

r_dist1 = as.matrix(r_dist1)

theta.ini = c(104.34, 0.17)
lo_bound = c(10, 1e-10)
up_bound = c(200,5)





source_python("BFuns_for_Bio.py")

cov_model = "Exp"
res1 = MLE_fit(y = r_Ytrain, X = r_Xtrain, D = r_dist1, cov_model = "Exp")
res1


thetahat1 = res1[[1]] 
betahat1  = res1[[4]] 
Sigma1 = thetahat1[2]*exp(-r_dist1/thetahat1[1]) 
diag(Sigma1) = thetahat1[2] + thetahat1[3]

cmat1 = thetahat1[2]*exp(-r_dist2/thetahat1[1]) 
pl1 = Xtest%*%betahat1
pe1 = cmat1%*% solve(Sigma1)%*%(Ytrain - Xtrain%*%betahat1) 
Ypred1 = pl1 + pe1

Ypred1 = round(Ypred1)
idx = Ypred1>3
Ypred1[idx]=3

misClassError <- mean(Ypred1 != Ytest) 
print(paste('Accuracy =', 1-misClassError))

plot(Ytest, Ypred1)
error = Ytest - Ypred1
plot(Ypred1, error)



(cmatrix = table(Ytest, Ypred1))

(Accuracy = (cmatrix[1] + cmatrix[5] + cmatrix[9])/sum(cmatrix))


