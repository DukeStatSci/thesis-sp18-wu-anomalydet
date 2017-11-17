setwd("~/Desktop/Stats Thesis/")
argus = readRDS("Dataset/argus_complete.rds")

library(softImpute)

###Full Matrix
# top_Sport = sort(unique(argus$Sport))
# top_Dport = sort(unique(argus$Dport))
# n_Sport = length(uniq_Sport)
# n_Dport = length(uniq_Dport)

# ###Subsetting Code
#matrix parameters
n_Sport = 20
n_Dport = 20

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
    obs = combination$SrcBytes
    n_obs = length(obs) #ignores NA values
    if (n_obs > 0){
      #obs = nscore(obs)$nscore #normal transformation
      for (i in 1:n_obs){
        ports_combo_matrix[[s,d]] = c(ports_combo_matrix[[s,d]],obs[i]) #O(1) time to append values to a list?
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

#untuned ALS using softimpute
fit = softImpute(ports_mean_matrix,rank.max=2,lambda=0.9,trace=TRUE,type="als")
fit$d
filled = complete(ports_mean_matrix, fit)
plot(ports_mean_matrix[!is.na(ports_mean_matrix)], filled[!is.na(ports_mean_matrix)])
plot(ports_mean_matrix[!is.na(ports_mean_matrix)], filled[!is.na(ports_mean_matrix)], 
     xlim = c(0,5000), ylim = c(0,5000))


####Eckhart Young Theorem Implementation, Best Rank k Approximation####
S = 1000
k = 2
Y = ports_mean_matrix

#calculate overall mean
n = 0
sum = 0
for (s in 1:n_Sport){
  for (d in 1:n_Dport){
    if (ports_freq_matrix[s,d] != 0){
      sum = sum + ports_mean_matrix[s,d]
      n = n + ports_freq_matrix[s,d]
    }
  }
}
overall_mean = sum/n
#calculate row means and col means
row_means = rowMeans(ports_mean_matrix, na.rm = TRUE)
col_means = colMeans(ports_mean_matrix, na.rm = TRUE)

#Fill in missing values in Y with ANOVA
for (s in 1:n_Sport){
  for (d in 1:n_Dport){
    if (ports_freq_matrix[s,d] == 0){
      Y[s,d] = row_means[s] + col_means[d] - overall_mean
    }
  }
}

for (i in 1:S){
  #extract SVD
  svd_Y = svd(Y)
  D = diag(svd_Y$d)
  U = svd_Y$u
  V = svd_Y$v
  #EYM theorem
  EYM = U[,1:k] %*% D %*% V[,1:k]
  #replace imputed values with new values 
  for (s in 1:n_Sport){
    for (d in 1:n_Dport){
      if (ports_freq_matrix[s,d] == 0){
        Y[s,d] = EYM[s,d]
      }
    }
  }
}


