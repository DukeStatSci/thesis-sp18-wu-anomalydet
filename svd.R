setwd("~/Desktop/Stats Thesis/Dataset/")
library(plyr)
library(data.table)
library(factoextra)
library(ggbiplot)
library(fa)
#read in data
argus = read.csv("argus-anon-20170201.csv")[1:50,] 
summary(argus)
#clean data
sapply(argus, class)
argus = transform(argus, 
                  Sport = as.factor(Sport),
                  Dport = as.factor(Dport))
argus = subset(argus, select = c("Flgs", "SrcAddr", "Sport", "DstAddr", "Dport", 
                                 "SrcPkts", "DstPkts", "SrcBytes", "DstBytes", "State"))
#attach(argus)
categorical = c("Flgs", "SrcAddr", "Sport", "DstAddr", "Dport", "State")
continuous = c("SrcPkts", "DstPkts", "SrcBytes", "DstBytes")

#create matrix 
uniq_Sport = sort(unique(Sport))
uniq_Dport = sort(unique(Dport))
n_Sport = length(uniq_Sport)
n_Dport = length(uniq_Dport)

all_ports_matrix = matrix(NA, nrow = n_Sport, ncol = n_Dport)
dimnames(all_ports_matrix) = list(uniq_Sport, uniq_Dport)
for (s in 1:n_Sport){
  for (d in 1:n_Dport){
    combination = argus[is.element(Sport, uniq_Sport[s]) & is.element(Dport, uniq_Dport[d]),]
    if (length(combination$SrcPkts) > 0){
      all_ports_matrix[s,d] = combination$SrcPkts
    }
  }
}
all_ports_matrix

#softImpute SVD
fit = softImpute(all_ports_matrix,rank.max=5,lambda=3,trace=TRUE)
#check cross validation for matrices, elementwise: randomly divide nonmissing elements into k groups, 
#leave out 1 group, fit the lowrank svd and try to predict the remaining
#modulo convergence criterion
fit$d
filled = complete(all_ports_matrix, fit)
#plot nonmissing vs true

plot(a[!is.na(a)], b[!is.na(a)])
plot(a[!is.na(a)], b[!is.na(a)], xlim = c(0,25), ylim = c(0,25))

#y bar is treated like its the same 

#write out the sparse matrix optimization problem into a statistical model
# y = UDV^T + sigma * E
# mxn mxr rxr rx n (r = rank)
# least squares problem: find U,D,V to minimize ||Y-UDV^T||^2. 
# The answer is take the SVD and the first r vectors
# Use alternating least squares to fill in the missing y's

#Write out the problem where the variances are different; it is basically finding an estimator in a statistical model


