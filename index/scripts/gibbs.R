library(MASS)

gibbs_sampler = function (X, y, W, sigma2_0 = 1, S = 1000){
  n = nrow(X)
  p = ncol(X)
  ### prior values 
  nu_0 = 2
  beta_0 = numeric(p)  
  gamma2 = 100
  S_0 = diag(p) * gamma2 #S_0 is p x p
  
  ### starting values
  set.seed(1)
  BETAs = matrix(nrow = S, ncol = p)
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
p = 5
beta = rnorm(p)
X = matrix(rnorm(n * p, mean=0, sd=1), n, p)
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

# The Bayes estimator (posterior mean) should be close to the GLS
# estimator. There are theoretical results that say how close
# the GLS estimator should be to the true value. In particular,
# the variance matrix of the GLS estimator around the true value is
# (X' V^{-1} X )^{-1}. So if you were to simulate many data sets,
# and get the posterior mean estimator for each, the variance of
# these simulated posterior means should be about (X' V^{-1} X )^{-1}.

##testing with GLS estimator:
#simulating multiple datasets and getting the posterior mean estimator for each , WHAT IS V in this case
#should I have multiple beta_post? if i simulate multiple X, what is the inputs to the 
#(X' V^{-1} X )^{-1}. whats the difference of using this to check versus using (X' V^{-1} X )^{-1}X'V^{-1}y on the beta post 
beta_post



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
# initialize U and V w/ latent factors, p = 5 how p is selected
p = 5
U = matrix(rnorm(m * p, mean=mu, sd=psi), m, p)
V = matrix(rnorm(n * p, mean=mu, sd=psi), n, p)

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
# This is a model selection choice - it might be good to try a few values
# and see how the results change. Alternatively, one could do cross
# validation, but this might be to computationally expensive.


#How to evaluate the quality of the output
#Simulate some data from the model, where the true mean
#matrix Theta is truly low rank, that is, Theta=UV' for
#some tall skinny matrices U and V

U = matrix(rnorm(m * p, mean=mu, sd=psi), m, p) 
V = matrix(rnorm(n * p, mean=0, sd=psi), n, p)
THETA = U %*% t(V)

# (a) Try to recover UV' from a full data set using your
# Gibbs sampler. This is the idealized case.
# 
# (b) Now pretend you don't have data for some cells.
# Use your Gibbs sampler to obtain an estimate of
# UV' in this case. How much worse do you do than
# in (a)? Do you do better for the cells for which
# you have data, than for the cells where you don't have data?


