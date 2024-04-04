library(gtools)
library(stats)
library(LaplacesDemon)
library(ggplot2)
library(MASS)
library(latex2exp)
library(mvtnorm)

# Load the data
setwd("/Users/seoyoung/Desktop/my_git/Advanced_bayes")
galaxy = read.table("galaxy.txt", header=T)

head(galaxy, 10)
summary(galaxy)

y = galaxy$speed
hist(y, breaks=50, freq=F)
n = length(y)


#### Gibbs Procedure
niter = 100000
nburn = 50000

## Specify the values of Hyper parameters and H
H = 20
H = 4
H = 2
H = 1

# hyper-parameters
a = rep(1/H, H)
mu0 = mean(y)
a_tau = 1
b_tau = 1
kappa0 = 10

# Initial value
newpi_h = rdirichlet(1, a)
newtau_h = rinvgamma(H, a_tau, b_tau)
newta_h = var(y)
newmu_h = rnorm(H, mu0, sqrt(kappa0*newtau_h))

mu_samp = matrix(nrow=(niter-nburn), ncol=H)
pi_samp = matrix(nrow=(niter-nburn), ncol=H)
tau_samp = matrix(nrow=(niter-nburn), ncol=H)
zi_samp = matrix(nrow=(niter-nburn), ncol=n)
kappa_samp = matrix(nrow=(niter-nburn), ncol=H)


## Sampling procedure

probs = matrix(rep(NA, n*H), ncol=H)
newzi = rep(NA, n)
nh = rep(NA, H)
ybar_h = rep(NA, H)
newb_tau = rep(NA, H)

for(i in 1:niter){
  
  # 1. Update zi from its conditional posterior
  for(h in 1:H) probs[, h] = exp(log(newpi_h[h]) + dnorm(y, newmu_h[h], newtau_h[h], log=T))
  for(j in 1:n) newzi[j] = sample(1:H, 1, prob = probs[j,])
  for(h in 1:H){
    nh[h] = sum(newzi == h)
    ybar_h[h] = ifelse(nh[h] > 0, mean(y[newzi==h]), 0)
  }
  
  # 2. Update (mu_h, tau_h) from its conditional posterior
  newkappa_h = 1/(1/kappa0 + nh)
  newa_tau = a_tau + nh/2
  for(h in 1:H){
    newb_tau[h] = ifelse(nh[h] > 0, b_tau + 0.5*(sum((y[newzi==h] - ybar_h[h])^2) + nh[h]/(1 + kappa0*nh[h])*(ybar_h[h] - mu0)^2), b_tau)
  }
  newtau_h = rinvgamma(H, newa_tau, newb_tau)
  newmu_h = rnorm(H, newkappa_h*(1/kappa0*mu0 + nh*ybar_h), sqrt(newkappa_h*newtau_h))
  
  # 3. Update (pi_1, ..., pi_H) from its conditional posterior
  newpi_h = rdirichlet(1, a + nh)
  
  # Store the samples after burn-in period
  if(i > nburn){
    mu_samp[i-nburn,] = newmu_h
    pi_samp[i-nburn,] = newpi_h
    tau_samp[i-nburn,] = newtau_h
    zi_samp[i-nburn,] = newzi
    kappa_samp[i-nburn,] = newkappa_h
  }
}



#### pointwise mean curve
y_grid = seq(5, 40, 0.1)
m=length(y_grid)
post_y = matrix(rep(0, (niter-nburn)*m), ncol = m)
for(i in 1:m){
  for(t in 1:(niter-nburn)){
    post_y[t, i] = sum(pi_samp[t, ] * dnorm(y_grid[i], mu_samp[t, ], sqrt(tau_samp[t, ])))
  }
}


post_mean = apply(post_y, 2, mean)
post_LB = apply(post_y, 2, quantile, 0.025)
post_UB = apply(post_y, 2, quantile, 0.975)



## Plot
ggplot() + 
  geom_histogram(mapping=aes(x=y, y=..density..), bins=100, fill="orange", alpha=0.6) +
  geom_line(mapping=aes(x=y_grid, y=post_mean), size=0.5, color="black") + 
  geom_line(mapping=aes(x=y_grid, y=post_LB), size=0.3, color="blue", lty = "dashed") +
  geom_line(mapping=aes(x=y_grid, y=post_UB), size=0.3, color="red",alpha=1.2 ,lty = "dashed")+
  ggtitle(paste0('Pointwise mean curve and 95% credible interval when H = ', H))+
  labs(x="speed")


