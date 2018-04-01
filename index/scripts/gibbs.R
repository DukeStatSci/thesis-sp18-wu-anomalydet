library(MASS)

gibbs_sampler = function (X, y, W, sigma2_0 = 1, S = 1000){
  n = nrow(X)
  r = ncol(X)
  ### prior values 
  nu_0 = 2
  beta_0 = numeric(r)  
  gamma2 = 100
  S_0 = diag(r) * gamma2 #S_0 is r x r
  
  ### starting values
  set.seed(1)
  BETAs = matrix(nrow = S, ncol = r)
  INV_SIGMA2s = c()
  beta = beta_0
  BETAs[1,] = beta
  inv_sigma2 = 1 / sigma2_0
  INV_SIGMA2s = c(INV_SIGMA2s, inv_sigma2)
  ### Gibbs sampling
  for(s in 2:S) {
    # generate a new Beta value from its full conditional
    Sigma_n = ginv ( ( t(X) %*% ginv(W) %*% X ) / ( 1 / inv_sigma2 ) + ginv(S_0) ) 
    beta_n = Sigma_n %*% ( ( t(X) %*% ginv(W) %*% y ) / ( 1 / inv_sigma2 ) + ginv(S_0) %*% beta_0)
    beta = mvrnorm( 1, beta_n, Sigma_n ) 
    
    # generate a new 1/sigma2 value from its full conditional
    SSR_W = ( t( y - X %*% beta ) %*% ginv(W)  %*% ( y - X %*% beta ))
    inv_sigma2 = rgamma( 1, ( nu_0 + n )/2, ( nu_0 * sigma2_0 + SSR_W) / 2)
    
    BETAs[s,] = beta
    INV_SIGMA2s = c(INV_SIGMA2s, inv_sigma2)
  }
  return (list(BETAs = BETAs, INV_SIGMA2s = INV_SIGMA2s))
}

#testing gibbs sampler
n = 100
r = 5
beta = rnorm(r)
X = matrix(rnorm(n * r, mean=0, sd=1), n, r)
W = diag(x = rexp(n), nrow = n, ncol = n)
e = rnorm(n,0,1)
sigma = 1/2
y = X %*% beta  + sigma * sqrt(W) %*% e

PHI = gibbs_sampler(X, y, W)

#evaluating posterior means of beta
BETAs = PHI$BETAs
beta_post = colMeans(BETAs)

#evaluating posterior mean of sigma
INV_SIGMA2s = PHI$INV_SIGMA2s
sigma_post = mean(sqrt(1/INV_SIGMA2s))

plot(beta, lm(y~ -1+ X)$coef)
abline(0,2, col = "red")

# > sigma
# [1] 0.5
# > sigma_post
# [1] 0.5267057
# > beta
# [1]  0.8960932 -0.7614885 -1.1147130 -0.7686936 -0.3393229
# > beta_post
# [1]  0.8956598 -0.7162321 -1.0892373 -0.8033861 -0.3284342

#as you increase n you should get closer to the truth, as 
#you decrease the effects of sigma*W*e you should get closer to the truth

##FULL procedure
Y = readRDS("data/means_SB.rds")
M = readRDS("data/freqs.rds")
m = nrow(Y)
n = ncol(Y)

S = 10 #number of iterations
##Set initial Values
mu = mean(Y, na.rm = TRUE) #overall mean
psi = sd(Y, na.rm = TRUE) #overall sd
a_i = rowMeans(Y, na.rm = TRUE)
b_j = colMeans(Y, na.rm = TRUE)
for (i in 1:m){
  for (j in 1:n){
    if (M[i,j] == 0){
      Y[i,j] = a_i[i] + b_j[j] - mu
    }
  }
}

# #initialize U and V w/ anova and svd
# svd_Y = svd(Y)
# U = svd_Y$u
# V = svd_Y$v
# initialize U and V w/ latent factors, r = 5 how r is selected
r = 5
U = matrix(rnorm(m * r, mean=mu, sd=psi), m, r)
V = matrix(rnorm(n * r, mean=mu, sd=psi), n, r)

#initialize sigmas  taus as the overall sd 
SIGMA = sqrt(rep(sd(Y, na.rm = TRUE), m))
TAU = sqrt(rep(sd(Y, na.rm = TRUE), n))

##Repeat Imputation
for (s in 1:S){
  #simulate U
  for (i in 1:m){
    W = diag(1, nrow = n, ncol = n) 
    ##w_jj = tau_j^2/s[i,j]
    for (j in 1:n){
      W[j,j] = W[j,j] * TAU[j]^2
      if (M[i,j] != 0){ 
        W[j,j] = W[j,j] / M[i,j]
      }
    }
    y = Y[i,]
    S_U = 3
    PHI  = gibbs_sampler(V, y, W, SIGMA[i], S_U)
    U[i,] = PHI$BETAs[S_U,]
    SIGMA[i] = sqrt(PHI$INV_SIGMA2s[S_U])
  }
  #impute V
  for (j in 1:n){
    W = diag(1, nrow = m, ncol = m) 
    for (i in 1:m){ 
      W[i,i] = W[i,i] * SIGMA[i]^2
      if (M[i,j] != 0){
        W[i,i] = W[i,i] / M[i,j]
      }
    }
    y = Y[,j]
    S_V = 3
    PHI  = gibbs_sampler(U, y, W, SIGMA[j], S_V)
    V[j,] = PHI$BETAs[S_V,]
    TAU[j] = sqrt(PHI$INV_SIGMA2s[S_V])
  }
  #impute y
  for (i in 1:m){
    for (j in 1:n){
      if (M[i,j] == 0){
        u_i = U[i,]
        v_j = V[j,]
        sigma_i = SIGMA[i]
        tau_j = TAU[j]
        Y[i,j] = rnorm(1, t(u_i) * v_j, sigma_i * tau_j) ##how does u_i and v_j become scalar????
      }
    }
  }
}

#How to select the number of latent factors?
#How to evaluate the quality of the output



# 1. Initialize $\sigma_i$ and $\tau_j$ as the overall standard deviation of the $Y^{(k)}$ matrix. 
# 2. Simulate $u_i$ and $\sigma_i$ using the generalized Gibbs Sampler. Set random values 
#     for the starting value of $X$, the algorithm will naturally converge to the true values of $v_j$.
# 3. Simulate $v_j$ and $\tau_j$ using the generalized Gibbs Sampler. 
#     Set random values for the starting value of $X$, the algorithm will naturally converge to the true values of $u_i$.
# 4. Fill in the missing values $y_{ij}$ in $Y^{(k)}$ by sampling from the normal 
#     distribution $$y_{ij} \sim N(u_i^Tv_j, \frac{\sigma_i^2\tau_j^2}{\sqrt{(n_ij)}})$$