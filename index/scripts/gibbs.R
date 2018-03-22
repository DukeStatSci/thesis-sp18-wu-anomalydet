library(MASS)

gibbs_sampler = function (X, y, W, S = 1000){
  m = nrow(X)
  ### prior values 
  nu_0 = 2
  sigma2_0 = 1
  beta_0 = numeric(m)  
  gamma2 = 100
  W_0 = diag(m) * gamma2 #W is m x m
  
  ### starting values
  set.seed(1)
  BETAs = list()
  INV_SIGMAs = list()
  beta = beta_0
  BETAs[[1]] = beta
  inv_sigma2 = 1 / sigma2_0
  INV_SIGMAs[[1]] = inv_sigma2
  ### Gibbs sampling
  for(s in 2:S) {
    
    # generate a new Beta value from its full conditional
    Sigma_n = ginv ( ( t(X) %*% ginv(W) %*% X ) / ( 1 / inv_sigma2 ) + ginv(W_0) ) 
    beta_n = Sigma_n %*% ( ( t(X) %*% ginv(W) %*% y ) / ( 1 / inv_sigma2 ) + ginv(W_0) %*% beta_0)
    beta = mvrnorm( 1, beta_n, Sigma_n ) 
    
    # generate a new 1/sigma2 value from its full conditional
    SSR_W = ( ginv( y - X %*% beta ) %*% W  %*% ( y - X %*% beta ))
    inv_sigma2 = rgamma(1, ( nu_0 + m )/2, ( nu_0 * sigma2_0 + SSR_W) / 2)
    
    BETAs[[s]] = beta
    INV_SIGMAs[[s]] = inv_sigma2
  }
  return (list(BETAs = BETAs, INV_SIGMAs = INV_SIGMAs))
}

#testing gibbs sampler
# X = matrix(rnorm(m * m, mean=0, sd=1), m, m) 
# W = diag(x = rnorm(m, mean=0, sd=1), nrow = m, ncol = m)
# y = matrix(rnorm(m * 1, mean=0, sd=1), m, 1) 
# 
# 
# X = matrix(1, m, m) 
# W = diag(1, nrow = m, ncol = m)
# y = matrix(1, m, 1) 
# 
# 
# gibbs_sampler(X, y, W, 100)

##FULL procedure

#initialize sigma_i  taU_j as the overall sd 
Y = readRDS("data/means_SB.rds")
M = readRDS("data/freqs.rds")
m = nrow(Y)
n = ncol(Y)


#calculate W matrix
sigma_i = sd(Y, na.rm = TRUE)
tau_j = sd(Y, na.rm = TRUE)

W = diag(sigma_i, nrow = m, ncol = m) # m x m
for (i in 1:m){
  for (j in 1:n){
    #how to divide by n_ij when the dimensions are different? whats the iteration for j
  }
}

# for (s in 1:S){
#   X =  matrix(nrow = m, ncol = r) ##CODE FOR SETTING X, X IS m x r but how do i get R
#   for (i in 1:m){
#     W = 
#   } 
# }





# 1. Initialize $\sigma_i$ and $\tau_j$ as the overall standard deviation of the $Y^{(k)}$ matrix. 
# 2. Simulate $u_i$ and $\sigma_i$ using the generalized Gibbs Sampler. Set random values 
#     for the starting value of $X$, the algorithm will naturally converge to the true values of $v_j$.
# 3. Simulate $v_j$ and $\tau_j$ using the generalized Gibbs Sampler. 
#     Set random values for the starting value of $X$, the algorithm will naturally converge to the true values of $u_i$.
# 4. Fill in the missing values $y_{ij}$ in $Y^{(k)}$ by sampling from the normal 
#     distribution $$y_{ij} \sim N(u_i^Tv_j, \frac{\sigma_i^2\tau_j^2}{\sqrt{(n_ij)}})$$