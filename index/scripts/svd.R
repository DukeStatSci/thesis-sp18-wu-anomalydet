# setwd("~/Desktop/Stats Thesis/")
# argus = readRDS("Dataset/argus_complete.rds")
setwd("C:/Users/Admin/Desktop/thesis-sp18-wu-anomalydet/index")
argus = readRDS("data/argus_complete.rds")

#matrix parameters
n_Sport = 100
n_Dport = 100

#get freqs
Sport_table = as.data.frame(table(argus$Sport))
Sport_table = Sport_table[order(-Sport_table$Freq),]
top_Sport = (head(Sport_table$Var1, n_Sport))

#get freqs
Dport_table = as.data.frame(table(argus$Dport))
Dport_table = Dport_table[order(-Dport_table$Freq),]
top_Dport = (head(Dport_table$Var1, n_Dport))

#create starting matrices
ports_combo_matrix = matrix(list(), nrow = n_Sport, ncol = n_Dport)
dimnames(ports_combo_matrix) = list(top_Sport, top_Dport)

ports_freq_matrix = matrix(0, nrow = n_Sport, ncol = n_Dport)
dimnames(ports_freq_matrix) = list(top_Sport, top_Dport)

nscore = function(x) {
  nscore = qqnorm(x, plot.it = FALSE)$x  # normal score 
  trn.table = data.frame(x=sort(x),nscore=sort(nscore))
  return (list(nscore=nscore, trn.table=trn.table))
}

#fill the ports_combo_matrix and ports_freq_matrix
for (s in 1:n_Sport){
  for (d in 1:n_Dport){
    combination = argus[is.element(argus$Sport, top_Sport[s])
                        & is.element(argus$Dport, top_Dport[d]),]
    obs = log(combination$DstPkts)
    n_obs = length(obs) #ignores NA values
    if (n_obs > 0){
      #obs = nscore(obs)$nscore #normal transformation
      for (i in 1:n_obs){
        ports_combo_matrix[[s,d]] = c(ports_combo_matrix[[s,d]],obs[i]) 
        #O(1) time to append values to a list?
        ports_freq_matrix[s,d] = ports_freq_matrix[s,d] + 1
      }
    }
  }
}

#create mean and variance matrix
ports_mean_matrix = matrix(NA, nrow = n_Sport, ncol = n_Dport)
dimnames(ports_mean_matrix) = list(top_Sport, top_Dport)

ports_variance_matrix = matrix(NA, nrow = n_Sport, ncol = n_Dport)
dimnames(ports_variance_matrix) = list(top_Sport, top_Dport)

#fill mean and variance matrix
for (s in 1:n_Sport){
  for (d in 1:n_Dport){
    if (ports_freq_matrix[s,d] == 1){
      ports_mean_matrix[s,d] = ports_combo_matrix[[s,d]]
      ports_variance_matrix[s,d] = 0
    }
    else if (ports_freq_matrix[s,d] > 1){
      ports_mean_matrix[s,d] = mean(ports_combo_matrix[[s,d]])
      ports_variance_matrix[s,d] = var(ports_combo_matrix[[s,d]])
    }
  }
}


### Filtering NA's
Y = Y[, colSums(is.na(Y)) != nrow(Y)] #remove NA cols
Y = Y[rowSums(is.na(Y)) != ncol(Y),] #remove NA rows
## REMOVING FULL 0 COLUMNS FROM MATRICES
M = M[ , !apply(M==0,2,all)]
M = M[ !apply(M==0,1,all) , ]
V = V[, colSums(is.na(V)) != nrow(V)] #remove NA cols
V = V[rowSums(is.na(V)) != ncol(V),]









































# #untuned ALS using softimpute
# library(softImpute)
# fit = softImpute(ports_mean_matrix,rank.max=3,lambda=0.9,trace=TRUE,type="als")
# fit$d
# filled = complete(ports_mean_matrix, fit)
# plot(ports_mean_matrix[!is.na(ports_mean_matrix)], filled[!is.na(ports_mean_matrix)])
# plot(ports_mean_matrix[!is.na(ports_mean_matrix)], filled[!is.na(ports_mean_matrix)],
#      xlim = c(0,5000), ylim = c(0,5000))

####Eckhart Young Theorem Implementation, Best Rank k Approximation####
matrix_complete = function(S = 1000, k = 2, n_Sport, n_Dport, Y, M){
  S = 1000
  Y_imputed = Y
  #overall mean
  n = sum(M)
  mu = sum(Y, na.rm = TRUE)/n
  #calculate row means and col means
  a_i = rowMeans(Y, na.rm = TRUE)
  b_j = colMeans(Y, na.rm = TRUE)
  #set NaN to 0 in means to fix anova fill in
  a_i = sapply(a_i, function(x) if (!is.finite(x)) {0} else {x})
  b_j = sapply(b_j, function(x) if (!is.finite(x)) {0} else {x})
  
  #Fill in missing values in Y_imputed with ANOVA
  for (s in 1:n_Sport){
    for (d in 1:n_Dport){
      if (M[s,d] == 0){
        Y_imputed[s,d] = a_i[s] + b_j[d] - mu
      }
    }
  }
  for (i in 1:S){
    #extract SVD
    svd_Y = svd(Y_imputed)
    D = diag((svd_Y$d)[1:k])
    U = svd_Y$u
    V = svd_Y$v
    #EYM theorem
    EYM = U[,1:k] %*% D %*% t(V[,1:k])
    #Replacing only missing means we cant assess fitted error
    for (s in 1:n_Sport){
      for (d in 1:n_Dport){
        if (M[s,d] == 0){
          Y_imputed[s,d] = EYM[s,d]
        }
      }
    }
  }
  return (Y_imputed)
}

# matrix_complete2 = function(Y, k = 2, lambda = 1.0){
#   fit = softImpute(Y, rank.max=k, lambda = lambda, trace=TRUE, type="als")
#   # fit$d
#   filled = complete(Y, fit)
#   return(filled)
# }

#ports_mean_matrix_imputed = matrix_complete(1000, 2, n_Sport, n_Dport, ports_mean_matrix, ports_freq_matrix)

#Relative distance using Frobenius Norm
relative_distance = function(Y, Y_imputed){
  return (frobenius.norm(Y - Y_imputed) / frobenius.norm(Y))
}

#Leave One Out Cross Validation
loocv = function (S = 1000, k = 2, nrows = n_Sport, ncols = n_Dport, Y, M){
  error = 0
  rmse = 0
  n = 0
  for (s in 1:nrows){
    for (d in 1:ncols){
      if (M[s,d] != 0){
        n = n + 1
        M_imputed = M
        true_sd = Y[s,d]
        M_imputed[s,d] = 0
        Y_imputed = matrix_complete(S, k, nrows, ncols, Y, M_imputed)
        # true_sd = Y[s,d]
        # Y[s,d] <- NA
        # Y_imputed = matrix_complete2(Y, k, 0.9)
        error = error + abs((Y_imputed[s,d] - true_sd))
        rmse = rmse + (Y_imputed[s,d] - true_sd)^2
        Y[s,d] = true_sd
      }
    }
  }
  rmse = sqrt(rmse/n)
  return (list(Error = error, RMSE = rmse, Observations = n))
}

RMSEs = lapply(seq(2,5,1), function(k) loocv(250,k,n_Sport, n_Dport, ports_mean_matrix, ports_freq_matrix)$RMSE)
# RMSEs = c()
# for (k in ranks){
#   print (k)
#   cv = loocv(250,k,n_Sport, n_Dport, ports_mean_matrix, ports_freq_matrix)
#   RMSEs = c(RMSEs, cv$RMSE)
# }
# plot(ranks, RMSEs)

# cv1 = loocv(250,1,n_Sport, n_Dport, ports_mean_matrix, ports_freq_matrix)
# cv2 = loocv(250,2,n_Sport, n_Dport, ports_mean_matrix, ports_freq_matrix)
# cv3 = loocv(250,3,n_Sport, n_Dport, ports_mean_matrix, ports_freq_matrix)
# cv5 = loocv(250,5,n_Sport, n_Dport, ports_mean_matrix, ports_freq_matrix)


##### LOW RANK MATRIX Approximation
library (softImpute)
x = matrix(sample(c(1,0),25, replace=TRUE),5,5)
m = x 
# m = matrix(list(), nrow = 5, ncol = 5)
# for (i in 1:5){
#   for(j in 1:5){
#     if(x[i,j] == 1){
#       m[i,j] = 1
#     }
#     else{
#       m[i,j] = 0
#     }
#   }
# }

cvs = lapply(seq(2,5,1), function(k) loocv(100,k,5,5,x,m)$RMSE)



