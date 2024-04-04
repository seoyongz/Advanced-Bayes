## Implement Generalized Linear Model with
## 1. Independent Metropolis-Hastings
## 2. Random walk Metropolis
## 3. Data augmentation

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
