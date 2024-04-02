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

library(Rcpp)
library(RcppArmadillo)
sourceCpp("/Users/seoyoung/Desktop/2/cov_mat.cpp")

# Load data set
house = read.csv("/Users/seoyoung/Desktop/2/고급베이즈/final/kc_house_data.csv")


# Haversine 공식을 사용하여 두 지점 사이의 거리를 계산하는 함수
haversine <- function(z1, z2) {
  # 지구의 반지름 (단위: 킬로미터)
  R <- 6371.0
  
  # 경도와 위도를 라디안으로 변환
  lat1 = z1[1] * (pi/180)
  lon1 = z1[2] * (pi/180)
  lat2 = z2[1] * (pi/180)
  lon2 = z2[2] * (pi/180)
  
  # Haversine 공식 적용
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  distance <- R * c
  
  return(distance)
}

# 예시: 서울과 뉴욕 사이의 거리를 계산
z1 = c(47.67, -122.392)
z2 = c(47.67, -122.394)

z1 = c(37.5665,126.9780)
z2 = c(40.7128,-74.0060)
seoul_lon <- 126.9780
seoul_lat <- 37.5665
new_york_lon <- -74.0060
new_york_lat <- 40.7128

distance <- haversine(z1, z2)



GP_fit_final1 = function(house, n_samp, m_pred, tau2, l2, niter, nburn, n_sigma, sigma2_max, file_name){
  # house input should have only "price", "long", "lat", sqft_living, "grade column
  # sqft_living and price are log transformed already
  n = n_samp   
  m = m_pred
  x_var_grid = (data.frame(long=house$long, lat=house$lat,  price=house$price))
  
  set.seed(603)
  train_idx = sample(1:nrow(house), n_samp)
  x_train = x_var_grid[train_idx,]
  y = x_train$price
  
  ## x_test : we are interested in
  set.seed(603)
  test_idx = sample(1:nrow(house), m_pred)
  x_test = x_var_grid[test_idx, c("long", "lat")]
  
  # lm_coef = lm(price ~ sqft_living, data=x_train)$coef
  
  
  # Calculate covariance function matrix
  
  Kxxt = haversine()
  Kxxt = cov_k(as.matrix(x_train), as.matrix(x_test), tau2, l2)
  Kxtxt = cov_k(as.matrix(x_test), as.matrix(x_test), tau2, l2)
  Kxx = cov_k(as.matrix(x_train), as.matrix(x_train), tau2, l2)
  #
  sigma2_grid = seq(0, sigma2_max, length.out=n_sigma + 1)[-1]
  Kxx_inv = array(dim=c(n_samp, n_samp, n_sigma))

  cat("=============== Compute inverse matrix ===============\n")
  for(i in 1:n_sigma){
    Kxx_inv[,,i] = chol2inv(Kxx + diag(rep(sigma2_grid[i], n_samp)))
    if(i %in% seq(0, n_sigma, n_sigma/10)[-c(1)]){
      cat(i/n_sigma*100, "% | ")
    }
  }
  
  # mu_bar = y[,3]
  # mu_bar = lm_coef[1] + x_train$sqft_living * lm_coef[2]
  # Sampling sigma2 with grid sampling
  sigma2_prob = c()
  cat("\n\n=============== Compute sigma2 probability ===============\n")
  for(i in 1:n_sigma){
    sigma2_prob[i] = exp(dmvnorm(y, mu_bar, Kxx + diag(n_samp)*sigma2_grid[i], log=T) - log(sigma2_grid[i]))
    if(i %in% seq(0, n_sigma, n_sigma/10)[-c(1)]){
      cat(i/n_sigma*100, "% | ")
    }
  }
  cat("\n>>>>> Max sigma2 prob :", max(sigma2_prob), "when sigma2 =", sigma2_grid[which.max(sigma2_prob)],"<<<<\n")
  
  sigma2_samp_idx = sample(seq(1, n_sigma), niter, replace=T, prob=sigma2_prob)
  sigma2_samp = sigma2_grid[sigma2_samp_idx]
  
  
  # Sampling mu_tilde
  mu_tilde = matrix(nrow=niter, ncol=m_pred)
  
  cat("\n=============== Sample mu_tilde =============== \n")
  for(t in 1:niter){
    tmp_cov =  Kxtxt - t(Kxxt) %*%Kxx_inv[,,sigma2_samp_idx[t]]%*%Kxxt
    tmp_mu = lm_coef[1] + lm_coef[2]*x_test$sqft_living + t(Kxxt)%*% Kxx_inv[,,sigma2_samp_idx[t]] %*% (y - mu_bar)
    mu_tilde[t, ] = mvrnorm(1, tmp_mu, tmp_cov)
    
    if(t %in% seq(0, niter, niter/10)[-c(1)]){
      cat(t/niter*100, "% | ")
    }
  }
  
  # Compute posterior mean and 95% credible interval
  post_mean = apply(mu_tilde, 2, mean)
  post_LB = apply(mu_tilde, 2, quantile, 0.025)
  post_UB = apply(mu_tilde, 2, quantile, 0.975)
  
  # Merge dataframe
  x_test$post_LB = post_LB
  x_test$post_UB = post_UB
  x_test$post_mean = post_mean
  
  fit_df = merge(house, x_test, by = c("long", "lat", "sqft_living"))
  fit_df$log_resid = fit_df$price - fit_df$post_mean
  fit_df$resid = exp(fit_df$price) - exp(fit_df$post_mean)
  
  # Define output
  output=list()
  output$mu_tilde = mu_tilde
  output$fit_df = fit_df
  output$tau2 = tau2
  output$l2 = l2
  output$niter = niter
  output$sigma2_prob = sigma2_prob
  output$sigma2_range = c(min(sigma2_grid), sigma2_max)
  output$coef = lm_coef
  # output$Kxx = Kxx
  # output$Kxxt = Kxxt
  # output$Kxtxt = Kxtxt
  # output$Kxx_inv = Kxx_inv
  
  
  save(output, file = paste0("/Users/seoyoung/Desktop/2/고급베이즈/final/fit_data/", file_name,".rda"))
  
  return(output)
}






# GP_fit5_mu fitting 
house_log4 = data.frame(price = log(house$price), long=house$long, lat=house$lat, sqft_living=log(house$sqft_living))


# same l2 tau2 (1)
GP5_final1_fit1 = GP_fit_final1(house_log4, n_samp=500, m_pred=700, tau2 = 1, l2=c(1, 1), 
                         niter=100, nburn=0, n_sigma = 200, sigma2_max=0.5, file_name = "GP_final_fit1")

# same l2 tau2 (2)
GP5_final1_fit2 = GP_fit_final1(house_log4, n_samp=500, m_pred=700, tau2 = 1, l2=c(0.5, 0.5), 
                                niter=100, nburn=0, n_sigma = 200, sigma2_max=0.5, file_name = "GP_final_fit2")


GP5_final1_fit7 = GP_fit_final1(house_log4, n_samp=2500, m_pred=1500, tau2 = 2, l2=c(1.5, 0.5), 
                                niter=300, nburn=0, n_sigma = 200, sigma2_max=0.5, file_name = "GP_final_fit7")


GP5_final1_fit8 = GP_fit_final1(house_log4, n_samp=2500, m_pred=1500, tau2 = 2, l2=c(1, 0.1), 
                                niter=300, nburn=0, n_sigma = 200, sigma2_max=0.5, file_name = "GP_final_fit8",
                                sigma2_prob = sigma2_prob, Kxx=Kxx, Kxxt=Kxxt, Kxtxt=Kxtxt, Kxx_inv=Kxx_inv)



fit_tmp = output

ggplot(fit_tmp$fit_df, aes(x=long, y=lat,  color=post_mean))+
  geom_point(size=0.7, alpha=0.7)+
  theme_dark() + 
  scale_color_gradient(high = "red", low =  "white", limits=c(11,15)) +
  ggtitle(paste0("Fitted with fixed tau2 = ", fit_tmp$tau2, " and l2 = ","(", fit_tmp$l2[1], ",", fit_tmp$l2[2],")", " with niter = ", fit_tmp$niter))



ggplot(fit_tmp$fit_df, aes(x=long, y=lat, color=log_resid))+
  geom_point(size=0.7, alpha=0.7)+
  theme_dark() + 
  scale_colour_gradient2(midpoint=0,high = "red", low =  "blue", mid="white", limits=c(-1,1)) +
  ggtitle(paste0("Fitted with fixed tau2 = ", fit_tmp$tau2, " and l2 = ","(", fit_tmp$l2[1], ",", fit_tmp$l2[2],")", " with niter = ", fit_tmp$niter))


ggplot(fit_tmp$fit_df, aes(x=long, y=lat, color=resid))+
  geom_point(size=0.7, alpha=0.7)+
  theme_dark() + 
  scale_colour_gradient2(midpoint=0,high = "red", low =  "blue", mid="white") +
  ggtitle(paste0("Fitted with fixed tau2 = ", fit_tmp$tau2, " and l2 = ","(", fit_tmp$l2[1], ",", fit_tmp$l2[2],")", " with niter = ", fit_tmp$niter))



ggplot(fit_tmp$fit_df, aes(x=long, y=lat,  color=price))+
  geom_point(size=0.7, alpha=0.7)+
  theme_dark() + 
  scale_color_gradient(high = "red", low =  "white", limits=c(11,15)) 


## 5p. scatter plot
ggplot(house, aes(x=sqft_living, y=price))+
  geom_point(size=0.7, alpha=0.7)
  
ggplot(house, aes(x=log(sqft_living), y=log(price)))+
  geom_point(size=0.7, alpha=0.7)


## 15p. Residual plot
ggplot(fit_tmp$fit_df, aes(x=post_mean, y=log_resid))+
  geom_point(size=0.7, alpha=0.7)+
  geom_hline(yintercept=0, color="red")

ggplot(fit_tmp$fit_df, aes(x=long, y=log_resid))+
  geom_point(size=0.7, alpha=0.7)+
  geom_hline(yintercept=0, color="red")

ggplot(fit_tmp$fit_df, aes(x=lat, y=log_resid))+
  geom_point(size=0.7, alpha=0.7)+
  geom_hline(yintercept=0, color="red")

ggplot(fit_tmp$fit_df, aes(x=sqft_living, y=log_resid))+
  geom_point(size=0.7, alpha=0.7)+
  geom_hline(yintercept=0, color="red")




ggplot(house, aes(x=long, y=lat,  color=grade))+
  geom_point(size=0.7, alpha=0.7)+
  theme_dark() + 
  scale_color_gradient(high = "red", low =  "white") +
  ggtitle("True distribution")







