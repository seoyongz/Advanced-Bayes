library(ggplot2)
library(lubridate)
library(GGally)
library(gtools)
library(stats)
library(LaplacesDemon)
library(MASS)
library(latex2exp)
library(mvtnorm)
library(rstan)

# Load data set
house = read.csv("/Users/seoyoung/Desktop/2/고급베이즈/final/kc_house_data.csv")
house$good_view = (ifelse(house$view>0, 1, 0))
head(house)


# house[house$grade %in% c(1,3,4), "grade2"] = 1
# house[house$grade==5, "grade2"] = 2
# house[house$grade==6, "grade2"] = 3
# house[house$grade==7, "grade2"] = 4
# house[house$grade==8, "grade2"] = 5
# house[house$grade==9, "grade2"] = 6
# house[house$grade==10, "grade2"] = 7
# house[house$grade==11, "grade2"] = 8
# house[house$grade==12, "grade2"] = 9
# house[house$grade==13, "grade2"] = 10


house_log = house[, c("long", "lat", "sqft_living", "price", "good_view")]
house_log$price = log(house_log$price)
house_log$sqft_living = log(house_log$sqft_living)


# Stan model 4 : 4 variables
GP_stan = "
data {  
  int<lower=1> N1;   // # of obs
  int<lower=1> N2;   // # of mu_tilde (we want to evaluate)
  int<lower=1> D;   // # of dimension of x
  real<lower=0> dist_tol;  // covariance matrix tolerance // if distance > tol then cov=0
  
  matrix[N1, D] x_bar;      // long, lat, log(sqft_living), grade indicator
  matrix[N2, D] x_tilde;
  vector[N1] y1;
  
  int<lower=0, upper=1> x1_train[N1];
  vector[N1] x2_train;
  int<lower=0, upper=1> x1_test[N2];
  vector[N2] x2_test;
}


parameters {
  real<lower=0> l2;
  real<lower=0> tau2;
  real<lower=0> sigma2;
  real b0;
  real b1;
  real b2;
}

transformed parameters{
  
  vector[N1] mu_bar;
  vector[N2] mu_tilde;
  real<lower=0> norm_tmp;
  real lat1;
  real lon1;
  real lat2;
  real lon2;
  real a;
  real c;
  matrix[N1, N1] Kxx;
  matrix[N2, N2] Kxtxt;
  
  // Compute covariance matrix
    for(i in 1:N1){
      mu_bar[i] = b0 + b1*x1_train[i] + b2*x2_train[i];
      for(j in i:N1){
        
        lat1 = x_bar[i, 1];
        lon1 = x_bar[i, 2];
        lat2 = x_bar[j, 1];
        lon2 = x_bar[j, 2];
        a = square(sin((lat2-lat1)/2)) + cos(lat1)*cos(lat2)*square(sin((lon2-lon1)/2));
        c = 2*atan2(sqrt(a), sqrt(1-a));
        
        norm_tmp = 6371.0 * c;
        
        // 가까우면 cov function 넣어주고
        if(norm_tmp < dist_tol){
          Kxx[i, j] = tau2 * exp(-norm_tmp/l2);
          Kxx[j, i] = tau2 * exp(-norm_tmp/l2);
        }
        // 멀면 두 점은 위치정보에 따른 가격이 independent
        else{
          Kxx[i, j] = 0;
          Kxx[j, i] = 0;
        }
      }
      Kxx[i, i] += sigma2;
    }
    
    
    for(i in 1:N2){
      mu_tilde[i] = b0 + b1*x1_test[i] + b2*x2_test[i];
      for(j in i:N2){
        norm_tmp = 0;
        
        lat1 = x_tilde[i, 1];
        lon1 = x_tilde[i, 2];
        lat2 = x_tilde[j, 1];
        lon2 = x_tilde[j, 2];
        a = square(sin((lat2-lat1)/2)) + cos(lat1)*cos(lat2)*square(sin((lon2-lon1)/2));
        c = 2*atan2(sqrt(a), sqrt(1-a));
        
        norm_tmp = 6371.0 * c;
        
        // 가까우면 cov function 넣어주고
        if(norm_tmp < dist_tol){
          Kxtxt[i, j] = tau2 * exp(-norm_tmp/l2);
          Kxtxt[j, i] = tau2 * exp(-norm_tmp/l2);
        }
        // 멀면 두 점은 위치정보에 따른 가격이 independent
        else{
          Kxtxt[i, j] = 0;
          Kxtxt[j, i] = 0;
        }
        
      }
      Kxtxt[i, i] += sigma2;
    }
}

model {
  
  // priors
  sigma2 ~ inv_gamma(3, 0.5);
  tau2 ~ normal(0, 1) T[0,];
  l2 ~ normal(0, 1) T[0,];
  
  // likelihood
  y1 ~ multi_normal(mu_bar, Kxx);
  
}

generated quantities{
  vector[N2] y2_pred;
  
  y2_pred = multi_normal_rng(mu_tilde, Kxtxt);
}
"


niter=600; nwarm=300; n_samp=600; m_pred=400
dist_tol = 2

set.seed(603)
house_samp_train = house_log[sample(1:nrow(house_log), n_samp), ]
set.seed(603)
house_samp_test = house_log[sample(1:nrow(house_log), m_pred), ]

stan_fit1 = stan(model_code = GP_stan,
                 data = list(y1=log(house_samp_train$price), 
                             x_bar = house_samp_train[,c("long", "lat")], 
                             x_tilde = house_samp_test[,c("long", "lat")],
                             x1_train = house_samp_train$good_view,
                             x1_test = house_samp_test$good_view,
                             x2_train = house_samp_train$sqft_living,
                             x2_test = house_samp_test$sqft_living,
                             N1 = n_samp, N2 = m_pred, D=2, dist_tol=dist_tol),
                 chains = 4,
                 warmup = nwarm,
                 iter = niter,
                 cores = 4
)
save(stan_fit1, file = "/Users/seoyoung/Desktop/자소서/고베방 플젝/stan_fit1.rda")

param = extract(stan_fit1)

Y_train_mean = extract(stan_fit1, "y1")
Y_mean_cred <- apply(Y_train_mean$y1_mean, 2, quantile, c(0.05, 0.95))
Y_mean_mean <- apply(Y_train_mean$y1_mean, 2, mean)

Y_pred <- extract(stan_fit1, "y2_pred")
Y_pred_cred <- apply(Y_pred$y2_pred, 2, quantile, c(0.05, 0.95))
Y_pred_mean <- apply(Y_pred$y2_pred, 2, mean)




#### Result
fit_df = house_samp_test
fit_df$log_resid = house_samp_test$price - Y_pred_mean
fit_df$resid = exp(house_samp_test$price) - exp(Y_pred_mean)
fit_df$y_pred = Y_pred_mean

ggplot(fit_df, aes(x=long, y=lat,  color=log_resid))+
  geom_point(size=0.7, alpha=0.7)+
  theme_dark() + 
  scale_color_gradientn(colours = c("blue", "white", "red"))
# ggtitle(paste0("Fitted with fixed tau2 = ", tau2, " and l2 = ", l2, "with niter = ", niter))


ggplot(house_samp_train, aes(x=long, y=lat,  color=Y_pred_mean))+
  geom_point(size=0.7, alpha=0.7)+
  theme_dark() + 
  scale_color_gradient(low="white", high="red")

ggplot(fit_df, aes(x=long, y=lat,  color=Y_pred_mean))+
  geom_point(size=0.7, alpha=0.7)+
  theme_dark() + 
  scale_color_gradient(low="white", high="red")

ggplot(fit_df, aes(x=long, y=lat,  color=price))+
  geom_point(size=0.7, alpha=0.7)+
  theme_dark() + 
  scale_color_gradient(low="white", high="red")


ggplot(fit_df, aes(x=long, y=lat,  color=Y_pred_mean))+
  geom_point(size=0.7, alpha=0.7)+
  theme_dark() + 
  scale_color_gradientn(colours = c("blue", "white", "red"))
  # ggtitle(paste0("Fitted with fixed tau2 = ", tau2, " and l2 = ", l2, "with niter = ", niter))



