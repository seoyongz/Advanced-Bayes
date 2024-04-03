
n = 500
j = 1:999
x0 = j/1000
k = 3
g = n

prob.H = c()
news2 = c()

L_temp = seq(1,40)       # prior of L
H_temp = L_temp + k + 1

niter = 10000
betaH_samp = matrix(rep(NA, niter*(max(L_temp)+k+1)), nrow=niter)
mu_x = matrix(rep(NA, niter*length(x0)), nrow=niter)
for(t in 1:niter){
  
  # Sample H from the posterior distribution
  for(h in L_temp){
    xh = seq(min(xi), max(xi), length.out = h)  # knots
    WH = bs(xi, knots=xh, degree=k, intercept=T)
    prob.H[h] = dpois(h, 1, log=T) + n/2*log(0.5 * t(y) %*% (diag(rep(1, n)) + (g/(1+g))*WH%*%ginv(t(WH)%*%WH)%*%t(WH)) %*%y)
  }
  
  log_prob_sum = max(prob.H) + log(sum(exp(prob.H - max(prob.H))))
  prob.H = exp(prob.H - log_prob_sum)
  
  # Given H, Construct B-spline basis
  newL = sample(L_temp, 1, prob=prob.H)
  newH = newL + k + 1
  xh = seq(min(xi), max(xi), length.out = newL)  # knots
  newWH = bs(xi, knots=xh, degree=k, intercept=T)
  
  # Sample sigma2 from the posterior distribution given H
  news2[t] = rinvgamma(1, n/2, 0.5 * t(y) %*% (diag(rep(1, n)) + (g/(1+g))*WH%*%ginv(t(WH)%*%WH)%*%t(WH)) %*%y)
  
  # Sample beta_H from the posterior distribution given H and sigma2
  betaH = mvrnorm(1, g/(1+g)*(ginv(t(newWH)%*%newWH)%*%t(newWH)%*%y), g*news2[t]/(1+g)*ginv(t(newWH)%*%newWH))
  
  # mu
  xh = seq(min(x0), max(x0), length.out = newL)
  WH = bs(x0, knots=xh, degree=k, intercept=T)
  for(i in 1:length(x0)){
    mu_x[t, i] = sum(betaH*WH[i, ])
  }
  
  # Store the parameters
  for(h in 1:H){
    betaH_samp[t, h] =  betaH[h]
  }
  
}

plot(news2, type="l")
plot(betaH_samp[,3], type="l")

post_mean_mux = apply(mu_x, 2, mean)
post_LB_mux = apply(mu_x, 2, quantile, 0.025)
post_UB_mux = apply(mu_x, 2, quantile, 0.975)

ggplot(mapping = aes(x=x0, y=post_mean_mux)) +
  geom_line(color="red") +
  geom_ribbon(aes(ymin=post_LB_mux, ymax=post_UB_mux), alpha=0.2) +
  geom_point(mapping=aes(x=xi, y=y), alpha=0.5)





# Stan model
stanmodel2 = "
data {
  int<lower=1> L;
  int<lower=0> N;
  int<lower=1, upper=L> ll[N];
  vector[N] x;
  vector[N] y;
}
parameters {
  real beta1;
  real beta2;
  real beta3;
  vector[L] u;
  real<lower=0> tau;
  real<lower=0> sigma2[L];
}

model {
  // priors
  for(l in 1:L){
    sigma2[l] ~ invgamma(0.001, 0.001);
  }
  beta1 ~ normal(192, 1000);
  beta2 ~ normal(728, 1000);
  beta3 ~ normal(353, 1000);
  tau ~ uniform(0, 1000);

  // likelihood
  for (l in 1:L){
    u[l] ~ normal(0, tau);
  }
  for(n in 1:N){
    y[n] ~ normal((beta1 + u[ll[n]]) / (1 + exp(-(x[n] - beta2)/beta3)), sqrt(sigma2[ll[n]]));
  }
}

generated quantities{
  real Y_mean[N];
  real Y_pred[N];
  for(n in 1:N){
    Y_mean[n] = (beta1 + u[ll[n]]) / (1 + exp(-(x[n] - beta2)/beta3));
    Y_pred[n] = normal_rng(Y_mean[n], sqrt(sigma2[ll[n]]));
  }
}
"

niter = 100000
nwarm = 50000
fit2 = stan(
  model_code = stanmodel2,
  data = list(y=Orange$circumference, x=Orange$age, N=nrow(Orange), L=5, ll=as.numeric(Orange$Tree)),
  chains = 4,
  warmup = nwarm,
  iter = niter,
  cores=4,
)

param2 = extract(fit2)
print(fit2)
