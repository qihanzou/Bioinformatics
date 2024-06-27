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

reduce_data = full_data[useful_vars]
colnames(reduce_data)[2] <- "AC"
colnames(reduce_data)[3] <- "FM"
colnames(reduce_data)[4] <- "UC"
colnames(reduce_data)[5] <- "DL"
colnames(reduce_data)[6] <- "DS"
colnames(reduce_data)[7] <- "DP"

pca_vars = c("LB", "AC", "FM", "UC", "DL", "DS", "DP", "ASTV",
                "MSTV", "ALTV", "MLTV", "Width","Min", "Max", "Nmax", "Nzeros", "Mode", "Mean",
                "Median", "Variance", "Tendency")


y_fullset = reduce_data$NSP
X_fullset = reduce_data[pca_vars]



res_pca = prcomp(X_fullset, scale = TRUE)

summary(res_pca)
print(res_pca$rotat, digits = 3)
print(res_pca$x, digits = 3)

lambda.v2 <- res_pca$sdev^2
# proportion of variance explained 
pve.v2 <- lambda.v2 / sum(lambda.v2)
#variance  of 1:
pve.v2[1] #0.2887278
#2 0.1670926
plot(pve.v2)


plot(res_pca, main = "", las = 1)



# Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
fviz_eig(res_pca)



fviz_pca_var(res_pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)





data = as.data.frame(res_pca$x[,1:21])
data$NSP = y_fullset
data$inter <- rep(1, dim(reduce_data)[1])


set.seed(2024) 
split <- sample.split(reduce_data$LB, SplitRatio = 0.7) 
train_data <- subset(data, split == "TRUE") 
test_data <- subset(data, split == "FALSE")


Xtrain = train_data[, 1:(length(data)-2)]
Xtest = test_data[, 1:(length(data)-2)]


Ytrain = train_data$NSP
Ytest = test_data$NSP



dist1 = as.matrix(dist(Xtrain))
dist2 = as.matrix(pdist(Xtest, Xtrain))


mds <- cmdscale(dist1, k = 2)
x <- mds[,1]
y <- mds[,2]


Xtrain = as.matrix(Xtrain)
Xtest = as.matrix(Xtest)
Ytrain = as.matrix(Ytrain)
Ytest = as.matrix(Ytest)



#data1 = cbind(x, y, Ytrain, train_data$inter)
#geodata = as.geodata(data1, coords.col = 1:2, data.col = 3, covar.col = 4, 
#                     covar.names = c("inter")) 
#vario = variog(geodata, option="bin", max.dist=500, trend = ~ inter)
#plot(vario,xlab = "h", ylab = "gamma(h)")
#intvalues = eyefit(vario, silent = FALSE)







res1 = MLE_fit(y = Ytrain, X = as.matrix(train_data$inter), D = dist1, cov_model = "Exp")

thetahat1 = res1[[1]] 
betahat1  = res1[[2]] 
Sigma1 = thetahat1[2]*exp(-dist1/thetahat1[1]) 
diag(Sigma1) = thetahat1[2] + thetahat1[3]

cmat1 = thetahat1[2]*exp(-dist2/thetahat1[1]) 
pl1 = as.matrix(test_data$inter)%*%betahat1
pe1 = cmat1%*% solve(Sigma1)%*%(Ytrain - as.matrix(train_data$inter)%*%betahat1) 
Ypred1 = pl1 + pe1

Ypred1 = round(Ypred1)
idx1 = Ypred1>3
Ypred1[idx1]=3
idx2 = Ypred1<1
Ypred1[idx2]=1



misClassError <- mean(Ypred1 != Ytest) 
print(paste('Accuracy =', 1-misClassError))
(cmatrix = table(Ytest, Ypred1))



x.comb = data.frame(train_data[, 1:(length(data)-2)])
x.test = data.frame(test_data[, 1:(length(data)-2)])

y.comb = data.frame(train_data$NSP)
y.test = data.frame(test_data$NSP)

# Partition the combined set into k parts
cv.k <- 10
n.comb <- dim(y.comb)[1]
if (n.comb %% cv.k == 0) {
  # Equal groups possible
  split.ind <- rep(1:cv.k, each = floor(n.comb/cv.k))
} else {
  # If groups are unequal, add the remaining terms
  split.ind <- c(rep(1:cv.k, each = floor(n.comb/cv.k)),
                 1:(n.comb %% cv.k))
}
split.ind <- sample(split.ind) # randomising



CV.error <- function(x.comb, y.comb, split.ind, nn.k, cv.k){
  misClassError = 0
  for (j in 1:cv.k) {
    x.tune <- data.frame(x.comb[split.ind == j,])
    y.tune <- data.frame(y.comb[split.ind == j,])
    x.train <- data.frame(x.comb[split.ind != j,])
    y.train <- data.frame(y.comb[split.ind != j,])
    
    classifier_knn <- knn(train = x.train, 
                          test = x.tune, 
                          cl = y.train[,1, drop = TRUE], 
                          k = nn.k)
    
    
    cm = table(y.tune[,1, drop = TRUE], classifier_knn)
    Accuracy = (cm[1] + cm[5] + cm[9])/sum(cm)
    misClassError = misClassError + mean(classifier_knn != y.tune[,1, drop = TRUE]) 
  }
  return(misClassError/cv.k)
}


K <- 50 # Largest choice of k
cv.err.k <- rep(0,K) # Vector of errors
for (nn.k in 1:K) {
  cv.err.k[nn.k] <- CV.error(x.comb, y.comb, split.ind,nn.k, cv.k)
}
cv.min.k <- which.min(cv.err.k) # k with smallest error
c(cv.min.k, min(cv.err.k))
plot(1:50, cv.err.k, xlab = 'k values', ylab = 'mis-Class Error', main = 'mis-Class Error and k-values plot')


# Fitting KNN Model to training dataset 
classifier_knn <- knn(train = x.comb, 
                      test = x.test, 
                      cl = y.comb[,1, drop = TRUE], 
                      k = cv.min.k) 
# Confusiin Matrix
(cm = table(y.test[,1, drop = TRUE], classifier_knn))
(Accuracy = (cm[1] + cm[5] + cm[9])/sum(cm))
misClassError <- mean(classifier_knn != y.test[,1, drop = TRUE]) 
print(paste('Accuracy =', 1-misClassError))




A_geo = c()
A_knn = c()
K_knn = c()
idx_s = c()

for (j in 1:21) {
  idx_s[j] = j
  print(j)
  data = as.data.frame(res_pca$x[,1:j])
  data$NSP = y_fullset
  data$inter <- rep(1, dim(reduce_data)[1])
  set.seed(2024) 
  split <- sample.split(reduce_data$LB, SplitRatio = 0.7) 
  train_data <- subset(data, split == "TRUE") 
  test_data <- subset(data, split == "FALSE")
  Xtrain = train_data[, 1:(length(data)-2)]
  Xtest = test_data[, 1:(length(data)-2)]
  Ytrain = train_data$NSP
  Ytest = test_data$NSP
  dist1 = as.matrix(dist(Xtrain))
  dist2 = as.matrix(pdist(as.matrix(Xtest), as.matrix(Xtrain)))
  
  mds <- cmdscale(dist1, k = 2)
  x <- mds[,1]
  y <- mds[,2]
  Xtrain = as.matrix(Xtrain)
  Xtest = as.matrix(Xtest)
  Ytrain = as.matrix(Ytrain)
  Ytest = as.matrix(Ytest)
  res1 = MLE_fit(y = Ytrain, X = as.matrix(train_data$inter), D = dist1, cov_model = "Exp")
  thetahat1 = res1[[1]] 
  betahat1  = res1[[2]] 
  
  Sigma1 = thetahat1[2]*exp(-dist1/thetahat1[1]) 
  diag(Sigma1) = thetahat1[2] + thetahat1[3]
  cmat1 = thetahat1[2]*exp(-dist2/thetahat1[1]) 
  pl1 = as.matrix(test_data$inter)%*%betahat1
  pe1 = cmat1%*% solve(Sigma1)%*%(Ytrain - as.matrix(train_data$inter)%*%betahat1) 
  Ypred1 = pl1 + pe1
  Ypred1 = round(Ypred1)
  idx1 = Ypred1>3
  Ypred1[idx1]=3
  idx2 = Ypred1<1
  Ypred1[idx2]=1
  misClassError <- mean(Ypred1 != Ytest) 
  print(paste('Accuracy =', 1-misClassError))
  cmatrix = table(Ytest, Ypred1)
  A_geo_j = 1-misClassError
  A_geo[j] = A_geo_j
  
  
  x.comb = data.frame(train_data[, 1:(length(data)-2)])
  x.test = data.frame(test_data[, 1:(length(data)-2)])
  y.comb = data.frame(train_data$NSP)
  y.test = data.frame(test_data$NSP)
  cv.k <- 10
  n.comb <- dim(y.comb)[1]
  if (n.comb %% cv.k == 0) {
    split.ind <- rep(1:cv.k, each = floor(n.comb/cv.k))
  } else {
    split.ind <- c(rep(1:cv.k, each = floor(n.comb/cv.k)),
                   1:(n.comb %% cv.k))
  }
  split.ind <- sample(split.ind) 
  CV.error <- function(x.comb, y.comb, split.ind, nn.k, cv.k){
    misClassError = 0
    for (j in 1:cv.k) {
      x.tune <- data.frame(x.comb[split.ind == j,])
      y.tune <- data.frame(y.comb[split.ind == j,])
      x.train <- data.frame(x.comb[split.ind != j,])
      y.train <- data.frame(y.comb[split.ind != j,])
      classifier_knn <- knn(train = x.train, 
                            test = x.tune, 
                            cl = y.train[,1, drop = TRUE], 
                            k = nn.k)
      cm = table(y.tune[,1, drop = TRUE], classifier_knn)
      Accuracy = (cm[1] + cm[5] + cm[9])/sum(cm)
      misClassError = misClassError + mean(classifier_knn != y.tune[,1, drop = TRUE]) 
    }
    return(misClassError/cv.k)
  }
  K <- 50 # Largest choice of k
  cv.err.k <- rep(0,K) # Vector of errors
  for (nn.k in 1:K) {
    cv.err.k[nn.k] <- CV.error(x.comb, y.comb, split.ind,nn.k, cv.k)
  }
  cv.min.k <- which.min(cv.err.k) # k with smallest error
  c(cv.min.k, min(cv.err.k))
  plot(1:50, cv.err.k, xlab = 'k values', ylab = 'mis-Class Error', main = 'mis-Class Error and k-values plot')
  classifier_knn <- knn(train = x.comb, 
                        test = x.test, 
                        cl = y.comb[,1, drop = TRUE], 
                        k = cv.min.k) 
  cm = table(y.test[,1, drop = TRUE], classifier_knn)
  Accuracy = (cm[1] + cm[5] + cm[9])/sum(cm)
  misClassError <- mean(classifier_knn != y.test[,1, drop = TRUE]) 
  print(paste('Accuracy =', 1-misClassError))
  A_knn[j] = 1-misClassError
  K_knn[j] = cv.min.k
}

