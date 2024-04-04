library(ggplot2)
library(rstan)
library(latex2exp)
library(bayestestR)
library(viridis)
library(MASS)

data(Orange)
head(Orange)
summary(Orange)

#########################################################
## Fit the oragne data with the nonlinear mixed model 

## (a) Overlay line plots for the trees to see how they grow over time.
ggplot(mapping=aes(x=age, y=circumference, group=Tree), data=Orange) + 
  geom_line(color=Orange$Tree) 

## (b) Use Stan, find the posterior distribution of each parameter
beta1 = c(100, 200, 300)
beta3 = c(100, 200, 300)
x = seq(0, 1599, length.out=400)

# with beta2 = 500
par(mfrow=c(3,3))
for(k in 1:3){
  for(i in 1:3){
    y = beta1[i]/(1 + exp(-(x-500)/beta3[k]))
    plot(x, y, type="l")
    title(main = paste0("(beta1, beta2, beta3) = (", beta1[i],",",500,",", beta3[k],")"))
  }
}

# Fit the nonlinear least squares
nlm = nls(circumference~ beta1/(1 + exp(-(age - beta2)/beta3)), start=list(beta1 = 100, beta2 = 500, beta3 = 200), data=Orange)
summary(nlm)
beta_hat = coef(nlm)


# Fitted model
y = beta_hat[1]/(1 + exp(-(x-beta_hat[2])/beta_hat[3]))
par(mfrow=c(1,1))
plot(x, y, type="l", ylim=c(min(Orange$circumference), max(Orange$circumference)))
points(x=Orange$age, y=Orange$circumference)
title(main = "Fitted Logistic Growth Model")


# bootstrap distribution for each parameter to assign priors
nboot = 1000
beta1_boot = c()
beta2_boot = c()
beta3_boot = c()
control = nls.control(maxiter = 1000)
for(i in 1:nboot){
  index.boot = sample(1:nrow(Orange), nrow(Orange) , replace=T)
  sample.boot = Orange[index.boot, ]
  beta.boot = coef(nls(circumference ~ beta1/(1 + exp(-(age - beta2)/beta3)), start=list(beta1 = 200, beta2 = 700, beta3 = 350), data=sample.boot))
  beta1_boot[i] = beta.boot[1]
  beta2_boot[i] = beta.boot[2]
  beta3_boot[i] = beta.boot[3]
}

par(mfrow=c(1,3))
hist(beta1_boot, breaks=30, main="bootstrap dist of beta1", freq=F)
abline(v=beta_hat[1], col="red", lwd=1.5)
abline(v=mean(beta1_boot), col="blue", lwd=1.5)

hist(beta2_boot, breaks=30, main="bootstrap dist of beta2", freq=F)
abline(v=beta_hat[2], col="red", lwd=1.5)
abline(v=mean(beta2_boot), col="blue", lwd=1.5)

hist(beta3_boot, breaks=30, main="bootstrap dist of beta3", freq=F)
abline(v=beta_hat[3], col="red", lwd=1.5)
abline(v=mean(beta3_boot), col="blue", lwd=1.5)
legend("topright",legend=c("beta_hat","mean of bootstrap samp"),fill=c("red","blue"),border="white",box.lty=0,cex=1.2)

# log of bootstrap distribution
par(mfrow=c(1,3))
hist(log(beta1_boot), breaks=30, main="bootstrap dist of log(beta1)", freq=F)
hist(log(beta2_boot), breaks=30, main="bootstrap dist of log(beta2)", freq=F)
hist(log(beta3_boot), breaks=30, main="bootstrap dist of log(beta3)", freq=F)


# Stan model

# beta ~ normal prior 
# tau ~ uniform prior
stanmodel = "
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
  real<lower=0> sigma2;
}

model {
  // priors
  sigma2 ~ inv_gamma(0.001, 0.001);
  beta1 ~ normal(192, 60);
  beta2 ~ normal(728, 210);
  beta3 ~ normal(353, 90);
  tau ~ uniform(0, 100);

  // likelihood
  for (l in 1:L){
    u[l] ~ normal(0, tau);
  }
  for(n in 1:N){
    y[n] ~ normal((beta1 + u[ll[n]]) / (1 + exp(-(x[n] - beta2)/beta3)), sqrt(sigma2));
  }
}

generated quantities{
  real Y_mean[N];
  real Y_pred[N];
  for(n in 1:N){
    Y_mean[n] = (beta1 + u[ll[n]]) / (1 + exp(-(x[n] - beta2)/beta3));
    Y_pred[n] = normal_rng(Y_mean[n], sqrt(sigma2));
  }
}
"

# stan fit
fit2 = stan(
  model_code = stanmodel,
  data = list(y=Orange$circumference, x=Orange$age, N=nrow(Orange), L=5, ll=as.numeric(Orange$Tree)),
  chains = 4,
  warmup = 20000,
  iter = 50000,
  cores=4,
)
print(fit2)


##### priors 1 ***
sigma2 ~ inv_gamma(0.001, 0.001);
beta1 ~ normal(192, 20);
beta2 ~ normal(728, 70);
beta3 ~ normal(353, 30);
tau ~ uniform(0, 50);


##### priors 2 ***
sigma2 ~ inv_gamma(0.001, 0.001);
beta1 ~ normal(192, 40);
beta2 ~ normal(728, 140);
beta3 ~ normal(353, 60);
tau ~ uniform(0, 100);



# trace plot of parameters
traceplot(fit, pars=c("beta1", "beta2", "beta3"))
traceplot(fit, pars=c("u[1]", "u[2]", "u[3]", "u[4]", "u[5]"))


# histogram of parameters
param = extract(fit)
par(mfrow=c(2,3))
for(k in 1:5){
  hist(param$u[,k], breaks=30)
}
hist(param$beta1)
hist(param$beta1)
hist(param$beta3)

Y_mean <- extract(fit, "Y_mean")
Y_mean_cred <- apply(Y_mean$Y_mean, 2, quantile, c(0.05, 0.95))
Y_mean_mean <- apply(Y_mean$Y_mean, 2, mean)

Y_pred <- extract(fit, "Y_pred")
Y_pred_cred <- apply(Y_pred$Y_pred, 2, quantile, c(0.05, 0.95))
Y_pred_mean <- apply(Y_pred$Y_pred, 2, mean)


res.obs=sqrt(colMeans((Orange$circumference - t(Y_mean$Y_mean))^2))
res.rep=sqrt(rowMeans((Y_pred$Y_pred-Y_mean$Y_mean)^2))	# just same as the posterior dist of sigma; hist(1/sqrt(param$tau))
plot(res.obs,res.rep); abline(a=0,b=1,col="red")
mean(res.obs>res.rep)


# fitting check with posterior mean
u_hat = c(mean(param$u[,1]), mean(param$u[,2]), mean(param$u[,3]), mean(param$u[,4]), mean(param$u[,5]))
beta_hat2 = c(mean(param$beta1), mean(param$beta2), mean(param$beta3))
y = matrix(rep(NA, 5*length(x)), ncol=5)

my_col = viridis(n = 5) 
par(mfrow=c(1,1))
u_index = unique(as.integer(Orange$Tree))
for(i in 1:5){
  y[,i] = (beta_hat2[1] + u_hat[u_index[i]])/(1 + exp(-(x - beta_hat2[2])/beta_hat2[3]))
  plot(x, y[,i], type="l", ylim=c(min(Orange$circumference), max(Orange$circumference)), col=my_col[i])
  points(x=Orange[Orange$Tree==i,]$age, y=Orange[Orange$Tree==i,]$circumference, col=my_col[i])
  par(new=T)
}





# Refit the data using a more suitable modeling framework and discuss the results.
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
  real<lower=0> sigma2[L];
  real<lower=0> tau;
}

model {
  // priors
  for(l in 1:L){
    sigma2[l] ~ normal(69, 20);
  }
  beta1 ~ normal(192, 40);
  beta2 ~ normal(728, 140);
  beta3 ~ normal(353, 60);
  tau ~ uniform(0, 100);

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

fit2 = stan(
  model_code = stanmodel2,
  data = list(y=Orange$circumference, x=Orange$age, N=nrow(Orange), L=5, ll=as.numeric(Orange$Tree)),
  chains = 4,
  warmup = 50000,
  iter = 100000,
  cores=4,
)
print(fit2)