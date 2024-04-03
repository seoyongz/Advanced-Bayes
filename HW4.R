# HW4

library(gtools)
library(stats)
library(LaplacesDemon)
library(ggplot2)
library(MASS)
library(latex2exp)
library(mvtnorm)

## Q1.(c)

# Initial Setting
sigma = 0.1
n = 200
# Generate x_i
set.seed(512)
xi = runif(n, 0, 1)
mux = sin(2*pi*xi^3)^3

# Generate y_i
set.seed(512)
yi = rnorm(n, mux, sigma)
plot(xi, yi, cex=0.5)

# x tilde
x_tilde = seq(1:999)/1000

GP_fit = function(yi, xi, x_tilde, tau, l, niter){
  n=200   # initial setting
  
  # covariance function
  cov_k = function(x1, x2, tau, l){
    mat = matrix(nrow = length(x1), ncol=length(x2))
    for(i in 1:length(x2)){
      mat[,i] = tau^2*exp(-((x1-x2[i])^2)/(l^2))
    }
    return(mat)
  }
  
  # Calculate covariance function matrix
  Kxxt = cov_k(xi, x_tilde, tau, l)
  Kxtxt = cov_k(x_tilde, x_tilde, tau, l)
  Kxx = cov_k(xi, xi, tau, l)
  
  sigma_grid = seq(0.01, 2, length.out=200)
  Kxx_inv = array(dim=c(n,n,niter))
  for(i in 1:niter){
    Kxx_inv[,,i] = solve(Kxx + diag(rep(sigma_grid[i], n)))
  }
  
  # Sampling sigma with grid sampling
  sigma_prob = c()
  for(i in 1:n){
    sigma_prob[i] = exp(dmvnorm(yi, rep(0, n), Kxx + diag(n)*sigma_grid[i]^2, log=T) - 2*log(sigma_grid[i]))
  }
  sigma_samp_n = sample(seq(1, length(sigma_grid)), niter, replace=T ,prob=sigma_prob)
  sigma_samp = sigma_grid[sigma_samp_n]
  
  # Sampling mu_tilde
  mu_tilde = matrix(nrow=niter, ncol=length(x_tilde))
  for(t in 1:niter){
    
    tmp_cov = Kxtxt - t(Kxxt) %*%Kxx_inv[,,sigma_samp_n[t]]%*%Kxxt
    tmp_mu = t(Kxxt)%*%Kxx_inv[,,sigma_samp_n[t]]%*%yi
    
    mu_tilde[t, ] = mvrnorm(1, tmp_mu, tmp_cov)
  }
  
  # Compute posterior mean and 95% credible interval
  post_mean = apply(mu_tilde, 2, mean)
  post_LB = apply(mu_tilde, 2, quantile, 0.025)
  post_UB = apply(mu_tilde, 2, quantile, 0.975)
  
  
  
  # Define output
  output=list()
  output$mu_tilde = mu_tilde
  output$post_mean = post_mean
  output$post_LB = post_LB
  output$post_UB = post_UB
  
  return(output)
  
}


# Fitting with various value of tau and l
# various l and fixed tau 
fit1 = GP_fit(yi, xi, x_tilde, tau=2, l=1, niter=1000)
fit2 = GP_fit(yi, xi, x_tilde, tau=2, l=0.5, niter=1000)
fit3 = GP_fit(yi, xi, x_tilde, tau=2, l=0.2, niter=1000)
fit4 = GP_fit(yi, xi, x_tilde, tau=2, l=0.1, niter=1000)

# various tau and fixed l 
fit5 = GP_fit(yi, xi, x_tilde, tau=0.05, l=0.1, niter=1000)
fit6 = GP_fit(yi, xi, x_tilde, tau=0.1, l=0.1, niter=1000)
fit7 = GP_fit(yi, xi, x_tilde, tau=0.5, l=0.1, niter=1000)
fit8 = GP_fit(yi, xi, x_tilde, tau=1, l=0.1, niter=1000)


# Draw the pointwise posterior mean curve and 95% credible band on $(0,1)$
ggplot(mapping = aes(x=x_tilde, y=fit1$post_mean)) +
  geom_line(mapping = aes(x=xi, y=sin(2*pi*xi^3)^3), color="red", size=1.0, alpha=0.7) +
  geom_line(color="blue", size=0.5) +
  geom_ribbon(aes(ymin=fit1$post_LB, ymax=fit1$post_UB), alpha=0.4) +
  geom_point(mapping=aes(x=xi, y=yi), alpha=0.3, size=0.2)+
  ggtitle(TeX('Posterior mean and the 95% credible interval when $\\tau = 2$, $l = 1$')) 

ggplot(mapping = aes(x=x_tilde, y=fit2$post_mean)) +
  geom_line(mapping = aes(x=xi, y=sin(2*pi*xi^3)^3), color="red", size=1.0, alpha=0.7) +
  geom_line(color="blue", size=0.5) +
  geom_ribbon(aes(ymin=fit2$post_LB, ymax=fit2$post_UB), alpha=0.4) +
  geom_point(mapping=aes(x=xi, y=yi), alpha=0.3, size=0.2)+
  ggtitle(TeX('Posterior mean and the 95% credible interval when $\\tau = 2$, $l = 0.5$')) 

ggplot(mapping = aes(x=x_tilde, y=fit3$post_mean)) +
  geom_line(mapping = aes(x=xi, y=sin(2*pi*xi^3)^3), color="red", size=1.0, alpha=0.7) +
  geom_line(color="blue", size=0.5) +
  geom_ribbon(aes(ymin=fit3$post_LB, ymax=fit3$post_UB), alpha=0.4) +
  geom_point(mapping=aes(x=xi, y=yi), alpha=0.3, size=0.2)+
  ggtitle(TeX('Posterior mean and the 95% credible interval when $\\tau = 2$, $l = 0.2$')) 

ggplot(mapping = aes(x=x_tilde, y=fit4$post_mean)) +
  geom_line(mapping = aes(x=xi, y=sin(2*pi*xi^3)^3), color="red", size=1.0, alpha=0.7) +
  geom_line(color="blue", size=0.5) +
  geom_ribbon(aes(ymin=fit4$post_LB, ymax=fit4$post_UB), alpha=0.4) +
  geom_point(mapping=aes(x=xi, y=yi), alpha=0.5, size=0.2)+
  ggtitle(TeX('Posterior mean and the 95% credible interval when $\\tau = 2$, $l = 0.1$')) 



ggplot(mapping = aes(x=x_tilde, y=fit5$post_mean)) +
  geom_line(mapping = aes(x=xi, y=sin(2*pi*xi^3)^3), color="red", size=1.0, alpha=0.7) +
  geom_line(color="blue", size=0.5) +
  geom_ribbon(aes(ymin=fit5$post_LB, ymax=fit5$post_UB), alpha=0.4) +
  geom_point(mapping=aes(x=xi, y=yi), alpha=0.3, size=0.2)+
  ggtitle(TeX('Posterior mean and the 95% credible interval when $\\tau = 0.05$, $l = 0.1$')) 

ggplot(mapping = aes(x=x_tilde, y=fit6$post_mean)) +
  geom_line(mapping = aes(x=xi, y=sin(2*pi*xi^3)^3), color="red", size=1.0, alpha=0.7) +
  geom_line(color="blue", size=0.5) +
  geom_ribbon(aes(ymin=fit6$post_LB, ymax=fit6$post_UB), alpha=0.4) +
  geom_point(mapping=aes(x=xi, y=yi), alpha=0.3, size=0.2)+
  ggtitle(TeX('Posterior mean and the 95% credible interval when $\\tau = 0.1$, $l = 0.1$')) 

ggplot(mapping = aes(x=x_tilde, y=fit7$post_mean)) +
  geom_line(mapping = aes(x=xi, y=sin(2*pi*xi^3)^3), color="red", size=1.0, alpha=0.7) +
  geom_line(color="blue", size=0.5) +
  geom_ribbon(aes(ymin=fit7$post_LB, ymax=fit7$post_UB), alpha=0.4) +
  geom_point(mapping=aes(x=xi, y=yi), alpha=0.3, size=0.2)+
  ggtitle(TeX('Posterior mean and the 95% credible interval when $\\tau = 0.5$, $l = 0.1$')) 

ggplot(mapping = aes(x=x_tilde, y=fit8$post_mean)) +
  geom_line(mapping = aes(x=xi, y=sin(2*pi*xi^3)^3), color="red", size=1.0, alpha=0.7) +
  geom_line(color="blue", size=0.5) +
  geom_ribbon(aes(ymin=fit8$post_LB, ymax=fit8$post_UB), alpha=0.4) +
  geom_point(mapping=aes(x=xi, y=yi), alpha=0.5, size=0.2)+
  ggtitle(TeX('Posterior mean and the 95% credible interval when $\\tau = 1$, $l = 0.1$')) 

















########################################################################################################################
## Q2.(b)

# Load the data
setwd("/Users/seoyoung/Desktop/2/고급베이즈")
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


