setwd("~/Desktop/Stats Thesis/thesis-sp18-wu-anomalydet/index")
Y = readRDS("data/means.RDS")
M = readRDS("data/freqs.RDS")

####Eckhart Young Theorem Implementation, Best Rank k Approximation####
matrix_complete = function(S = 1000, k = 2, nrows, ncols, Y, M){
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
  #Y_imputed = outer(1:nrow(Y), 1:ncol(Y), function(r,c) ifelse(M[r,c] == 0, a_i[r] + b_j[c] - mu, Y[r,c]))
  for (i in 1:nrows){
    for (j in 1:ncols){
      if (M[i,j] == 0){
        Y_imputed[i,j] = a_i[i] + b_j[j] - mu
      }
    }
  }
  for (s in 1:S){
    #extract SVD
    svd_Y = svd(Y_imputed)
    U = svd_Y$u
    V = svd_Y$v
    #EYM theorem
    if (k == 1){
      EYM = (matrix(U[,1:k]) * (svd_Y$d)[1:k]) %*% t(matrix(V[,1:k]))
    }
    else {
      EYM = U[,1:k] %*% diag((svd_Y$d)[1:k]) %*% t(V[,1:k])
    }
    for (i in 1:nrows){
      for (j in 1:ncols){
        if (M[i,j] == 0){
          Y_imputed[i,j] = EYM[i,j]
        }
      }
    }
  }
  return (Y_imputed)
}

#Leave One Out Cross Validation
loocv = function (S = 1000, k = 2, nrows = nrows, ncols = ncols, Y, M){
  error = 0
  rmse = 0
  n = 0
  for (i in 1:nrows){
    for (j in 1:ncols){
      if (M[i,j] != 0){
        n = n + 1
        M_imputed = M
        true_sd = Y[i,j]
        M_imputed[i,j] = 0
        Y_imputed = matrix_complete(S, k, nrows, ncols, Y, M_imputed)
        error = error + abs((Y_imputed[i,j] - Y[i,j]))
        rmse = rmse + (Y_imputed[i,j] - Y[i,j])^2
      }
    }
  }
  rmse = sqrt(rmse/n)
  return (list(Error = error, RMSE = rmse, Observations = n))
}

RMSEs = lapply(seq(1,10,1), function(k) loocv(250,k, 20, 20, Y, M)$RMSE)

# generate m x n matrix with rank r, add noise
generate_low_rank_matrix = function(m, n, r, noise = FALSE){
  A = matrix(rnorm(m * r, mean=0, sd=1), m, r) 
  B = matrix(rnorm(r * n, mean=0, sd=1), r, n)
  AB = A %*% B
  if (noise){
    E = matrix(rnorm(m * n, mean=0, sd=1), m, n)
    AB = AB + E
  }
  return (AB)
}

# sets each elements to missing one at a time
full_loocv = function (S = 1000, k, nrows, ncols, Y){
  error = 0
  rmse = 0
  n = 0
  for (i in 1:nrows){
    for (j in 1:ncols){
      n = n + 1
      M = Y
      M[i,j] = 0
      Y_imputed = matrix_complete(S, k, nrows, ncols, Y, M)
      error = error + abs((Y_imputed[i,j] - Y[i,j]))
      rmse = rmse + (Y_imputed[i,j] - Y[i,j])^2
    }
  }
  rmse = sqrt(rmse/n)
  return (list(Error = error, RMSE = rmse, Observations = n))
}

approximate_rank = function(Y, M, S = 50, simulated = TRUE){
  nrows = nrow(Y)
  ncols = ncol(Y)
  ranks = c()
  all_errors = list()
  for (i in 1:S){
    cat("Approximation Round: ", i, '\n')
    if (simulated){
      cv_errors = lapply(seq(1,10,1), function(k) full_loocv(200,k,nrows,ncols,Y)$RMSE)
    }
    else{
      cv_errors = lapply(seq(1,10,1), function(k) loocv(200,k,nrows,ncols,Y, M)$RMSE)
    }
    all_errors = c()
    low_rank = which.min(cv_errors)
    ranks = c(ranks, low_rank)
  }
  return (ranks)
}

ranks = approximate_rank(Y, M, 10, simulated = FALSE)

# TESTING WITH SIMULATED DATA
Y1 = generate_low_rank_matrix(10,10,1)
ranks1 = approximate_rank(Y1, M, 20, simulated = TRUE)

Y2 = generate_low_rank_matrix(10,10,2)
ranks2 = approximate_rank(Y2, M, 20, simulated = TRUE)

Y3 = generate_low_rank_matrix(10,10,3)
ranks3 = approximate_rank(Y3, M, 20, simulated = TRUE)

Y4 = generate_low_rank_matrix(10,10,4)
ranks4 = approximate_rank(Y4, M, 20, simulated = TRUE)

Y5 = generate_low_rank_matrix(10,10,5)
ranks5 = approximate_rank(Y5, M, 20, simulated = TRUE)




