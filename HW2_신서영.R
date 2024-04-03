library(MASS)
library(ggplot2)
library(LaplacesDemon)
library(mvtnorm)
library(viridisLite)

vir30 = viridis(n = 30)   # color set
################## Q1 ####################
# Setting for (d)
m = 30
n = 10
d = 2
N = n*m

X = array(rep(0, n*m*2), dim = c(n, d, m))
set.seed(503)
X[,1,] = rep(1, n*m)
X[,2,] = runif(n*m)

sigma2 = 0.5
mu_beta = c(0.5, 1)
cov_beta = matrix(c(2, 1.9, 1.9, 2), nrow=2)

# Generate beta_j and y_j
set.seed(506)
beta = t(mvrnorm(m, mu_beta, cov_beta))     # dxm(2x30) matrix
y = matrix(rep(NA, n*m), ncol=m)
for(j in 1:m){
  y[,j] = mvrnorm(1, X[,,j]%*%beta[,j], diag(rep(sigma2, n)))
}


# 1)
par(mfrow=c(1,1))
t=seq(0,1,0.1)
for(j in 1:m){
  plot(t, beta[1,j] + beta[2,j]*t, type="l", col=vir30[j], ylim=c(-6, 9))
  axis(side=2, at=seq(-6, 9, 2))
  par(new=T)
}



# 2) 

rho = d + 1
Psi = diag(rep(1, d))
# Use hyperparameters that make the priors for mu_beta, sigma^2 vague
xi = c(0,0)
Omega = diag(rep(1, d))
nu = 1
tau2 = 0.1


# Program and run the sampler to obtain the posterior, posterior predictive
# Assign data type
niter = 50000  # the number of iteration
nburn = 10000
count = 1

new_beta = matrix(rep(NA, d*m), ncol=m)   # each column represents beta_j(subject-specific random effects)
new_sigma2 = NA
new_mu_beta = c(rep(NA, d))
new_cov_beta = matrix(rep(NA, d^2), nrow=d)


beta_samp = array(rep(NA, d*m*(niter-nburn)), dim=c(d,m,niter-nburn))
sigma2_samp = rep(NA, niter-nburn)
mu_beta_samp = matrix(rep(NA, (niter-nburn)*d), nrow=niter-nburn)
cov_beta_samp = array(rep(NA, d*d*(niter-nburn)), dim=c(d,d,niter-nburn))
y_rep1 = array(rep(NA, n*m*(niter-nburn)), dim=c(n,m,niter-nburn))


## Set initial value
new_beta = matrix(rep(0, d*m), ncol=m)
new_sigma2 = 0.1
new_mu_beta = c(0, 0)
new_cov_beta = matrix(c(1, 0.9, 0.9, 1), nrow=2)
y = y

for(i in 1:niter){
  
  # (1) Update beta_j
  for(j in 1:m){
    post_mean_beta = solve(t(X[,,j])%*%X[,,j]/new_sigma2 + solve(new_cov_beta)) %*% t(t(y[,j])%*%X[,,j]/new_sigma2 + t(new_mu_beta)%*%solve(new_cov_beta))
    post_var_beta = solve(t(X[,,j])%*%X[,,j]/new_sigma2 + solve(new_cov_beta))
    new_beta[,j] = mvrnorm(1, post_mean_beta, post_var_beta)
  }
  
  # (2) Update sigma^2
  post_a_sigma2 = N + nu
  post_b_sigma2 = 0
  for(j in 1:m){
    post_b_sigma2 = post_b_sigma2 + t(y[,j]-X[,,j]%*%new_beta[,j])%*%(y[,j]-X[,,j]%*%new_beta[,j])
  }
  post_b_sigma2 = (post_b_sigma2 + nu*tau2)/(N + nu)
  new_sigma2 = rinvchisq(1, post_a_sigma2, post_b_sigma2) 
  
  # (3) Update mu_beta
  post_mean_mu_beta = solve(m*solve(new_cov_beta) + solve(Omega))%*%(apply(solve(new_cov_beta)%*%new_beta, 1, sum) + solve(Omega)%*%xi)
  post_var_mu_beta = solve(m*solve(new_cov_beta) + solve(Omega))
  new_mu_beta = mvrnorm(1, post_mean_mu_beta, post_var_mu_beta)
  
  # (4) Update cov_beta
  post_nu_cov_beta = m + rho
  post_S_cov_beta = matrix(rep(0, d^2), nrow=2)
  for(j in 1:m){
    post_S_cov_beta = post_S_cov_beta + (new_beta[,j] - new_mu_beta)%*%t(new_beta[,j] - new_mu_beta)
  }
  post_S_cov_beta = post_S_cov_beta + Psi
  new_cov_beta = rinvwishart(post_nu_cov_beta, post_S_cov_beta)
  
  
  
  # Burning
  if(i > nburn){
  
    # Store samples
    beta_samp[,,count] = new_beta
    sigma2_samp[count] = new_sigma2
    mu_beta_samp[count,] = new_mu_beta
    cov_beta_samp[,,count] = new_cov_beta
    
    # Sampling y_rep
    for(j in 1:m){
      y_rep1[,j,count] = mvrnorm(1, X[,,j]%*%new_beta[,j], new_sigma2*diag(rep(1,n)))
    }

    count = count+1
  }
  
}


## Trace plot
# mu_beta_samp
par(mfrow=c(1, 2))
for(i in 1:2){
  plot(mu_beta_samp[, i], type="l", main=paste0("mu_beta_samp","[", ",", i, "]"))
}

# cov_beta_samp
par(mfrow=c(2, 2))
for(i in 1:2){
  for(j in 1:2){
    plot(cov_beta_samp[i, j,], type="l", main=paste0("cov_beta_samp","[", i,",",j, "]"))
  }
}

# sigma2_samp
par(mfrow=c(1,1))
plot(sigma2_samp, type="l", main="sigma2_samp")


## Histogram

# marginal posterior for mu_beta
par(mfrow=c(1,2))
for(j in 1:2){
  hist(mu_beta_samp[,j], axes=F, main=paste0("mu_beta_samp","[",",", j, "]"))
  abline(v=c(mean(mu_beta_samp[,j]), mu_beta[j]), col = c("red", "blue"))
  axis(side=1, at=round(c(min(mu_beta_samp[,j]), mean(mu_beta_samp[,j]), max(mu_beta_samp[,j])), 3))
}

# marginal posterior for cov_beta
par(mfrow=c(2,2))
for(i in 1:2){
  for(j in 1:2){
    hist(cov_beta_samp[i,j,], axes=F, main=paste0("cov_beta_samp","[", i,",", j ,"]"), xlab="", breaks=15)
    abline(v=c(mean(cov_beta_samp[i,j,]), cov_beta[i,j]), col = c("red", "blue"))
    axis(side=1, at=round(c(min(cov_beta_samp[i,j,]), mean(cov_beta_samp[i,j,]), max(cov_beta_samp[i,j,])), 2))
  }
}

# marginal posterior for sigma^2
par(mfrow=c(1,1))

hist(sigma2_samp, axes=F, main="sigma2_samp", xlab="")
abline(v=c(mean(sigma2_samp), sigma2), col = c("red", "blue"))
axis(side=1, at=c(min(sigma2_samp), mean(sigma2_samp), sigma2, max(sigma2_samp)))



# 3) posterior predictive check
test_T1 = rep(0, niter-nburn)
test_T = rep(0, niter-nburn)

for(t in 1:(niter-nburn)){
  for(j in 1:m){
    test_T1[t] = test_T1[t] + sum((y_rep1[,j,t] - X[,,j]%*%beta_samp[,j,t])^2)
    test_T[t] = test_T[t] + sum((y[,j] - X[,,j]%*%beta_samp[,j,t])^2)
  }
}
mean(test_T1-test_T>0)


# 4) Repeat 2), 3) using the prior in (c)

# 2) using the prior in (c)
rho = d + 1

# Setting hyperparameters
xi = c(0,0)
Omega = diag(rep(3, d))
nu = 1
tau2 = 1 
psi2 = 10


# Program and run the sampler to obtain the posterior, posterior predictive
# Assign data type and initial values
niter = 50000  # the number of iteration
nburn = 10000
count = 1

new_beta = matrix(rep(NA, d*m), ncol=m)   # each column represents beta_j(subject-specific random effects)
new_sigma2 = NA
new_mu_beta = c(rep(NA, d))
new_sigma2_beta = c(rep(NA, d))
new_cov_beta = matrix(rep(NA, d^2), nrow=d)

beta_samp = array(rep(NA, d*m*(niter-nburn)), dim=c(d,m,(niter-nburn)))
sigma2_samp = rep(NA, (niter-nburn))
mu_beta_samp = matrix(rep(NA, (niter-nburn)*d), nrow=(niter-nburn))
cov_beta_samp = array(rep(NA, d*d*(niter-nburn)), dim=c(d,d,(niter-nburn)))
y_rep2 = array(rep(NA, n*m*(niter-nburn)), dim=c(n,m,niter-nburn))

## Set initial value
new_beta = matrix(rep(0, d*m), ncol=m)
new_sigma2 = 0.5
new_mu_beta = c(0.5, 1)
new_sigma2_beta = 1/rchisq(2, rho, psi2)
new_cov_beta = diag(new_sigma2_beta)
y = y


for(i in 1:niter){
  
  # (1) Update beta_j
  ### new_beta matrix is a 
  for(j in 1:m){
    post_mean_beta = solve(t(X[,,j])%*%X[,,j]/new_sigma2 + solve(new_cov_beta))%*%t(t(y[,j])%*%X[,,j]/new_sigma2 + t(new_mu_beta)%*%solve(new_cov_beta))
    post_var_beta = solve(t(X[,,j])%*%X[,,j]/new_sigma2 + solve(new_cov_beta))
    new_beta[,j] = mvrnorm(1, post_mean_beta, post_var_beta)
  }
  
  
  # (2) Update sigma^2
  post_a_sigma2 = N + nu
  post_b_sigma2 = 0
  for(j in 1:m){
    post_b_sigma2 = post_b_sigma2 + t(y[,j]-X[,,j]%*%new_beta[,j])%*%(y[,j]-X[,,j]%*%new_beta[,j])
  }
  post_b_sigma2 = (post_b_sigma2 + nu*tau2)/(N + nu)
  new_sigma2 = rinvchisq(1, post_a_sigma2, post_b_sigma2) 
  
  # (3) Update mu_beta
  post_mean_mu_beta = solve(m*solve(new_cov_beta) + solve(Omega))%*%(apply(solve(new_cov_beta)%*%new_beta, 1, sum) + solve(Omega)%*%xi)
  post_var_mu_beta = solve(m*solve(new_cov_beta) + solve(Omega))
  new_mu_beta = mvrnorm(1, post_mean_mu_beta, post_var_mu_beta)
  
  
  # (4) Update cov_beta
  ## different from (b)
  post_a_sigma2_beta = m + rho
  post_S_cov_beta = matrix(rep(0, d^2), nrow=2)
  for(k in 1:d){
    for(j in 1:m){
      post_S_cov_beta[k] = post_S_cov_beta[k] + (new_beta[k,j] - new_mu_beta[k])^2
    }
    post_S_cov_beta[k] = (post_S_cov_beta[k] + psi2)/(m+rho)
    new_sigma2_beta[k] = rinvchisq(1, post_a_sigma2_beta, post_S_cov_beta[k])
  }
  new_cov_beta = diag(new_sigma2_beta)
  
  
  # Burning
  if(i > nburn){
  
    # Store samples
    beta_samp[,,count] = new_beta
    sigma2_samp[count] = new_sigma2
    mu_beta_samp[count,] = new_mu_beta
    cov_beta_samp[,,count] = new_cov_beta
    
    # Sampling y_rep
    for(j in 1:m){
      y_rep2[,j,count] = mvrnorm(1, X[,,j]%*%new_beta[,j], new_sigma2*diag(rep(1,n)))
    }
  
    count = count+1
  
  }
}

## Trace plot
par(mfrow=c(1,2))

for(i in 1:2){
  plot(beta_samp[i,1,], type="l", main=paste0("beta_samp","[", i,",", "]"), col=vir30[j])
  for(j in 2:m){
    par(new=T)
    plot(beta_samp[i,j,], type="l", main=paste0("beta_samp","[", i,",", "]"), col=vir30[j])
  }
}

# mu_beta_samp
par(mfrow=c(1, 2))
for(i in 1:2){
  plot(mu_beta_samp[, i], type="l", main=paste0("mu_beta_samp","[", ",", i, "]"))
}

# cov_beta_samp
par(mfrow=c(1, 2))
for(i in 1:2){
  plot(cov_beta_samp[i, i,], type="l", main=paste0("cov_beta_samp","[", i,",",i, "]"))
}

# sigma2_samp
par(mfrow=c(1,1))
plot(sigma2_samp, type="l", main="sigma2_samp")


## Histogram
# marginal posterior for mu_beta
par(mfrow=c(1,2))
for(j in 1:2){
  hist(mu_beta_samp[,j], axes=F, main=paste0("mu_beta_samp","[",",", j, "]"))
  abline(v=c(mean(mu_beta_samp[,j]), mu_beta[j]), col = c("red", "blue"))
  axis(side=1, at=c(min(mu_beta_samp[,j]), mean(mu_beta_samp[,j]), max(mu_beta_samp[,j])))
}


# marginal posterior for cov_beta
par(mfrow=c(1,2))
for(i in 1:2){
  hist(cov_beta_samp[i,i,], axes=F, main=paste0("cov_beta_samp","[", i,",", i ,"]"), xlab="")
  abline(v=c(mean(cov_beta_samp[i,i,]), cov_beta[i,i]), col = c("red", "blue"))
  axis(side=1, at=round(c(min(cov_beta_samp[i,i,]), mean(cov_beta_samp[i,i,]), max(cov_beta_samp[i,i,])),4))
}

# marginal posterior for sigma^2
par(mfrow=c(1,1))

hist(sigma2_samp, axes=F, main="sigma2_samp", xlab="")
abline(v=c(mean(sigma2_samp), sigma2), col = c("red", "blue"))
axis(side=1, at=c(min(sigma2_samp), mean(sigma2_samp), sigma2, max(sigma2_samp)))




# 3) posterior predictive check using the prior in (c)
test_T2 = rep(0, niter-nburn)
test_T = rep(0, niter-nburn)

for(t in 1:(niter-nburn)){
  for(j in 1:m){
    test_T2[t] = test_T2[t] + sum((y_rep2[,j,t] - X[,,j]%*%beta_samp[,j,t])^2)
    test_T[t] = test_T[t] + sum((y[,j] - X[,,j]%*%beta_samp[,j,t])^2)
  }
}
mean(test_T2-test_T>0)


################## Q2 ####################
# Setting for (c)
n = 500
beta = c(1, 2)
set.seed(506)
X = matrix(rnorm(500*2), nrow = 500)
y = rbinom(500, 1, pnorm(X%*%beta))

# Get beta_hat and V_pos from glm fit
m_fit = glm(y ~ X[,1] + X[,2] - 1, family=binomial(link="probit"))
summary(m_fit)
beta_hat = as.vector(m_fit$coeff)
eta_hat = X%*%beta_hat
V_pos = -y/pnorm(eta_hat)^2*(eta_hat*dnorm(eta_hat)*pnorm(eta_hat) + dnorm(eta_hat)^2) + 
          (1-y)/(1-pnorm(eta_hat))^2*(eta_hat*dnorm(eta_hat)*(1-pnorm(eta_hat)) - dnorm(eta_hat)^2)
V_pos = solve(-t(X)%*%diag(c(V_pos))%*%X)


### Independent Metropolis-Hastings
count = 0
nthin = 5
nburn = 5000
niter = 20000
beta_samp1 = matrix(rep(NA, (niter-nburn)/nthin*2), ncol=2)
old_beta = mvrnorm(1, beta_hat, V_pos) # initial setting
accept_ratio1 = 0

for(t in 1:niter){
    new_beta = mvrnorm(1, beta_hat, V_pos)
    new_log_like = sum(y*pnorm(X%*%new_beta, log=T) + (1-y)*pnorm(-X%*%new_beta, log=T))
    old_log_like = sum(y*pnorm(X%*%old_beta, log=T) + (1-y)*pnorm(-X%*%old_beta, log=T))
    ratio = new_log_like - old_log_like 
    ratio = ratio + dmvnorm(old_beta, beta_hat, V_pos, log=T) - dmvnorm(new_beta, beta_hat, V_pos, log=T)
    
    if(ratio > 0.0) {accept = 1}
    else{
      u = log(runif(1))
      if(u < ratio) {accept = 1}
      else {accept = 0}
    }
    
    if(accept==1){
      old_beta = new_beta
      accept_ratio1 = accept_ratio1 + 1/niter
    }
    else {new_beta = old_beta}
  
  # Burning and Thinning
  if(t >= nburn){
    if(t %% nthin == 0){
      beta_samp1[count, ] = old_beta
      count = count+1
    }
  }
}


### Random walk Metropolis
count = 0
nthin = 5
nburn = 5000
niter = 20000
beta_samp2 = matrix(rep(NA, (niter-nburn)/nthin*2), ncol=2)
old_beta = c(0, 0) # initial setting
accept_ratio2 = 0
V_pos2 = 2.38^2/2*V_pos

for(t in 1:niter){
  new_beta = mvrnorm(1, old_beta, V_pos2)
  
  new_log_like = sum(y*pnorm(X%*%new_beta, log=T) + (1-y)*pnorm(-X%*%new_beta, log=T))
  old_log_like = sum(y*pnorm(X%*%old_beta, log=T) + (1-y)*pnorm(-X%*%old_beta, log=T))
  ratio = new_log_like - old_log_like + dmvnorm(old_beta, beta_hat, V_pos2, log=T) - dmvnorm(new_beta, beta_hat, V_pos2, log=T)
  
  if(ratio > 0.0) {accept=1}
  else{
    u = log(runif(1, 0, 1))
    if(u < ratio) {accept = 1}
    else {accept = 0}
  }
  
  if(accept==1){
    old_beta = new_beta
    accept_ratio2 = accept_ratio2 + 1/niter
  }
  else {new_beta = old_beta}
  
  # Burning and Thinning
  if(t >= nburn){
    if(t %% nthin == 0){
      beta_samp2[count, ] = old_beta
      count = count+1
    }
  }
  
}

### Data augmentation
niter = 10000
beta_samp3 = matrix(rep(NA, niter*2), ncol=2)
new_beta = beta_hat
u = mvrnorm(1, X%*%beta_hat, diag(rep(1, n)))

for(t in 1:niter){
  
  new_beta = mvrnorm(1, solve(t(X)%*%X)%*%t(X)%*%u, solve(t(X)%*%X))
  for(i in 1:n){
    tmp = rnorm(1, X[i,]%*%new_beta, 1)
    if(y[i]==0){
      while(tmp >= 0) tmp = rnorm(1, X[i,]%*%new_beta, 1)
    }
    else{
      while(tmp < 0) tmp = rnorm(1, X[i,]%*%new_beta, 1)
    }
    u[i] = tmp
  }
  
  beta_samp3[t,] = new_beta
  
}

# Trace plot
par(mfrow=c(2,3))
for(j in 1:2){
  plot(beta_samp1[,j], type="l", main=paste0("beta_samp1","[", ",", j, "]"))
  plot(beta_samp2[,j], type="l", main=paste0("beta_samp2","[", ",", j, "]"))
  plot(beta_samp3[,j], type="l", main=paste0("beta_samp3","[", ",", j, "]"))
}


# marginal posterior of beta
par(mfrow=c(2,3))

for(j in 1:2){
  hist(beta_samp1[,j], axes=F, main=paste0("beta_samp1","[", ",", j, "]"), breaks=15)
  abline(v=c(mean(beta_samp1[,j]), j), col = c("red", "blue"))
  axis(side=1, at=round(c(min(beta_samp1[,j]), mean(beta_samp1[,j]), max(beta_samp1[,j])),3))
  
  hist(beta_samp2[,j], axes=F, main=paste0("beta_samp1","[", ",", j, "]"), breaks=15)
  abline(v=c(mean(beta_samp2[,j]), j), col = c("red", "blue"))
  axis(side=1, at=round(c(min(beta_samp2[,j]), mean(beta_samp2[,j]), max(beta_samp2[,j])),3))
  
  hist(beta_samp3[,j], axes=F, main=paste0("beta_samp1","[", ",", j, "]"), breaks=15)
  abline(v=c(mean(beta_samp3[,j]), j), col = c("red", "blue"))
  axis(side=1, at=round(c(min(beta_samp3[,j]), mean(beta_samp3[,j]), max(beta_samp3[,j])),3))
}

