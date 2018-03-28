library(MASS)

gibbs_sampler = function (X, y, W, sigma2_0 = 1, S = 100){
  n = nrow(X)
  r = ncol(X)
  ### prior values 
  nu_0 = 2
  beta_0 = numeric(n)  
  gamma2 = 100
  W_0 = diag(n) * gamma2 #W is n x n
  
  ### starting values
  set.seed(1)
  BETAs = matrix(nrow = S, ncol = n)
  INV_SIGMA2s = c()
  beta = beta_0
  BETAs[1,] = beta
  inv_sigma2 = 1 / sigma2_0
  INV_SIGMA2s = c(INV_SIGMA2s, inv_sigma2)
  ### Gibbs sampling
  for(s in 2:S) {
    
    # generate a new Beta value from its full conditional
    Sigma_n = ginv ( ( t(X) %*% ginv(W) %*% X ) / ( 1 / inv_sigma2 ) + ginv(W_0) ) 
    beta_n = Sigma_n %*% ( ( t(X) %*% ginv(W) %*% y ) / ( 1 / inv_sigma2 ) + ginv(W_0) %*% beta_0)
    beta = mvrnorm( 1, beta_n, Sigma_n ) 
    
    # generate a new 1/sigma2 value from its full conditional
    SSR_W = ( t( y - X %*% beta ) %*% ginv(W)  %*% ( y - X %*% beta ))
    inv_sigma2 = rgamma(1, ( nu_0 + n )/2, ( nu_0 * sigma2_0 + SSR_W) / 2)
    
    BETAs[s,] = beta
    INV_SIGMA2s = c(INV_SIGMA2s, inv_sigma2)
  }
  return (list(BETAs = BETAs, INV_SIGMA2s = INV_SIGMA2s))
}

#testing gibbs sampler
n = 100
X = matrix(rnorm(n * n, mean=0, sd=1), n, n)
W = diag(x = rexp(n), nrow = n, ncol = n)
y = matrix(rnorm(n * 1, mean=0, sd=1), n, 1)
S = 1000
sigma2_0 = 1/2
PHI = gibbs_sampler(X, y, W, sigma2_0, S)

#run the algorithm more than 10 times, check that the beta posterior mean (column means of the beta matrix object)
#get a distribution of betas, calculate the posterior mean of those, get a distribution of sigma2 (make sure it is, dont use variance of 1)

#Evaluating Performance for BETA
BETAs = PHI$BETAs
beta = BETAs[S,]
y = X%*%beta  + matrix(rnorm(n * 1, mean=0, sd=1), n, 1)

#beta_samp = rnorm(n)
plot(beta, lm(y~ -1+ X)$coef)
abline(0,2, col = "red")

# Evaluating Performance for sigma2
INV_SIGMA2s = PHI$INV_SIGMA2s
1/INV_SIGMA2s[1000] #last sigma2, should be close to 1/2
1/mean(INV_SIGMA2s) #mean sigma2, should be close to 1/2
1/median(INV_SIGMA2s) #median sigma2, should be close to 1/2

##FULL procedure

Y = readRDS("data/means_SB.rds")
M = readRDS("data/freqs.rds")
m = nrow(Y)
n = ncol(Y)

S = 100 #number of iterations

##Set initial Values
#initialize U and V w/ anova and svd
mu = sum(Y, na.rm = TRUE)/n
a_i = rowMeans(Y, na.rm = TRUE)
b_j = colMeans(Y, na.rm = TRUE)
for (i in 1:m){
  for (j in 1:n){
    if (M[i,j] == 0){
      Y[i,j] = a_i[i] + b_j[j] - mu
    }
  }
}
svd_Y = svd(Y)
U = svd_Y$u
V = svd_Y$v

#initialize sigmas  taus as the overall sd 
SIGMA = rep(sd(Y, na.rm = TRUE), m)
TAU = rep(sd(Y, na.rm = TRUE), n)


##Repeat Imputation
for (s in 1:S){
  #impute U
  for (i in 1:m){
    W = diag(1, nrow = n, ncol = n) 
    #WHAT TO DO WHEN M IS NOT SQUARE AND NEED TO GET OBSERVATIONS??
    ##w_jj = tau_j^2/s[i,j]
    for (j in 1:n){
      W[j,j] = W[j,j] * TAU[j]^2
      if (M[i,j] != 0){ #WHAT HAPPENS WHEN IT IS EQUAL TO 0???
        W[j,j] = W[j,j] / M[i,j]
      }
      # for (j in 1:n){
      #   if (M[i,j] != 0 && i == j){
      #     W[i,j] = W[i,j] / M[i,j]
      #   }
      # }
    }
    y = Y[i,]
    S_U = 100
    PHI  = gibbs_sampler(V, y, W, SIGMA[s], S_U)
    U[i,] = PHI$BETAs[S_U,]
    SIGMA[i] = mean(sqrt(PHI$INV_SIGMA2s))
  }
  #impute V
  for (j in 1:m){
    W = diag(1, nrow = m, ncol = m) 
    #WHAT TO DO WHEN M IS NOT SQUARE AND NEED TO GET OBSERVATIONS??
    ##w_jj = tau_j^2/s[i,j]
    for (i in 1:m){ 
      W[i,i] = W[i,i] * SIGMA[i]^2
      if (M[i,j] != 0){ #WHAT HAPPENS WHEN IT IS EQUAL TO 0???
        W[i,i] = W[i,i] / M[i,j] # IS THIS RIGHT?
      }
      # for (j in 1:n){
      #   if (M[i,j] != 0 && i == j){
      #     W[i,j] = W[i,j] / M[i,j]
      #   }
      # }
    }
    y = Y[,j]
    S_V = 100
    PHI  = gibbs_sampler(U, y, W, SIGMA[s], S_V)
    V[j,] = PHI$BETAs[S_V,]
    SIGMA[i] = mean(PHI$INV_SIGMAs)
  }
  #impute y
  for (i in 1:m){
    for (j in 1:n){
      if (M[i,j] == 0){
        u_i = U[i,]
        v_j = V[j,]
        sigma_i = SIGMA[i]
        tau_j = TAU[j]
        Y[i,j] = rnorm(t(u_i) * v_j, sigma_i^2 * tau_j^2) ##how does u_i and v_j become scalar????
      }
    }
  }
}





# 1. Initialize $\sigma_i$ and $\tau_j$ as the overall standard deviation of the $Y^{(k)}$ matrix. 
# 2. Simulate $u_i$ and $\sigma_i$ using the generalized Gibbs Sampler. Set random values 
#     for the starting value of $X$, the algorithm will naturally converge to the true values of $v_j$.
# 3. Simulate $v_j$ and $\tau_j$ using the generalized Gibbs Sampler. 
#     Set random values for the starting value of $X$, the algorithm will naturally converge to the true values of $u_i$.
# 4. Fill in the missing values $y_{ij}$ in $Y^{(k)}$ by sampling from the normal 
#     distribution $$y_{ij} \sim N(u_i^Tv_j, \frac{\sigma_i^2\tau_j^2}{\sqrt{(n_ij)}})$$