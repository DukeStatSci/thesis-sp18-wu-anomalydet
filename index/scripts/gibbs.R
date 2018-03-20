library(MASS)

### data
X = matrix() #input data, what exactly is X
y = vector() #what is y
W = matrix() #what exactly is w
n = 0 #what exacltyis n

gibbs_sampler = function (X, y, W, n_ij, S = 1000){
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
  SIGMAs = list()
  BETAs[1] = beta = X # what initial starting value for X
  INV_SIGMAs[1] = inv_sigma2 = 1 / sigma2_0 # what initial starting value for inv_sigma/use inv_sigma or sigma 
  ### Gibbs sampling
  for(s in 2:S) {
    
    # generate a new Beta value from its full conditional
    Sigma_n = inv ( ( t(X) %*% W %*% X ) / ( 1 / inv_sigma2 ) + inv(W_0) ) 
    beta_n = Sigma_n %*% ( ( t(X) %*% W %*% X ) %*% y * ( 1 / inv_sigma2 ) + inv(W_0) %*% beta_0)
    beta = mvrnorm( 1, beta_n, Sigma_n ) # how to sample MVN, use MASS?
    
    # generate a new 1/sigma2 value from its full conditional
    SSR_W = ( inv( y - X %*% beta ) %*% W  %*% ( y - X %*% beta ))
    inv_sigma2 = rgamma(1, ( nu_0 + n_ij )/2, ( nu_0 * sigma2_0 + SSR_W) / 2)
    
    BETAs[s] = beta
    INV_SIGMAs[s] = inv_sigma2
  }
  return (BETAs = BETAs, INV_SIGMAs = INV_SIGMAs)
}

###
