library(ggplot2)
library(lubridate)
library(GGally)

house = read.csv("/Users/seoyoung/Desktop/2/고급베이즈/final/kc_house_data.csv")


# EDA 1
ggplot()+
  geom_point(aes(x=long, y=lat, color=log(price)), size=0.3, data=house)+
  scale_color_gradient(low="white", high="red")

summary(factor(house$view))
ggplot()+
  geom_point(aes(x=long, y=lat), alpha=0.1, size=0.3, data=house[house$view==0,]) +
  geom_point(aes(x=long, y=lat, color=factor(view)), size=0.3, data=house[house$view>0,])
  
ggplot()+
  geom_point(aes(x=long, y=lat), alpha=0.1, size=0.3, data=house[house$view==0,]) +
  geom_point(aes(x=long, y=lat), color="red",size=0.3, data=house[house$view>0,])
house$good_view = ifelse(house$view>0, 1, 0)



ggplot(data=house)+
  geom_boxplot(aes(x=good_view, y=log(price), group=good_view))


# EDA 2
## log(sqft_living) vs. log(price) 
ggplot()+
  geom_point(house, mapping=aes(x=log(sqft_living), y=log(price)), size=0.1) +
  theme_minimal()

ggplot()+
  geom_point(house, mapping=aes(x=sqft_living, y=log(price)), size=0.1) +
  theme_minimal()


## 그룹별 log(sqft_living) vs. log(price) 
### grade2별 log(sqft_living) vs. log(price) 

house$grade2 = house$grade
house[house$grade %in% c(1,3,4), "grade2"] = 4
house$grade2 = as.factor(house$grade2)
table(factor(house$grade2))

ggplot()+
  geom_point(house, mapping=aes(x=log(sqft_living), y=log(price)), size=0.1) +
  facet_wrap(~grade2)+
  theme_minimal()
  # stat_smooth(data=house, mapping=aes(x=log(sqft_living), y=log(price)), method = lm, colour="red", lwd=0.4)

ggplot(house, mapping=aes(x=log(sqft_living), y=log(price), group=grade2, color=grade2) )+
  geom_point(size=0.3, alpha=0.5) +
  stat_smooth(method = lm, se=F, lwd=0.7)+
  theme_minimal()+
  ggtitle("Linear regression for each grade2 group")
  


### View별 log(sqft_living) vs. log(price) 
house$view = factor(house$view)
ggplot()+
  geom_point(house, mapping=aes(x=log(sqft_living), y=log(price)), size=0.1) +
  facet_wrap(~view) + 
  theme_minimal()

ggplot(house, mapping=aes(x=log(sqft_living), y=log(price), group=view, color=view))+
  geom_point(size=0.3, alpha=0.5) +
  stat_smooth(method = lm, se=F, lwd=0.7)+
  theme_minimal()



### condition별 log(sqft_living) vs. log(price) 
house$condition = factor(house$condition)
ggplot()+
  geom_point(house, mapping=aes(x=log(sqft_living), y=log(price)), size=0.1) +
  facet_wrap(~condition)+
  theme_minimal()

ggplot(house, mapping=aes(x=log(sqft_living), y=log(price), group=condition, color=condition))+
  geom_point(size=0.3, alpha=0.5) +
  stat_smooth(method = lm, se=F, lwd=0.7) +
  theme_minimal()


### lat별 log(sqft_living) vs. log(price) 
ggplot(house, mapping=aes(x=log(sqft_living), y=log(price), color=lat))+
  geom_point(size=0.3, alpha=0.5) +
  theme_minimal()


### lat 그룹별 log(sqft_living) vs. log(price) 

#### lat histogram
summary(house$lat)
boxplot(house$lat)

lat_grid0 = seq(round(min(house$lat),2)-0.01, round(max(house$lat), 2), 0.01)
ggplot(house, aes(x=lat))+
  geom_density(alpha=0.1, size=0.05, fill="red", color="red") + 
  geom_histogram(aes(y=..density..), 
                 bins=200, color="black", fill="gray", alpha=0.3, size=0.25)+
  scale_x_continuous(breaks = lat_grid0) +
  theme_minimal()+
  theme(axis.text.x=element_text(angle=70, hjust=1)) + 
  ggtitle("Histogram of latitude")


#### lat 그룹화 : lat_group1
vintercept1 = c(47.23, 47.415, 47.605)
lat_grid = seq(round(min(house$lat),2)-0.01, round(max(house$lat), 2), 0.05)
ggplot(house, aes(x=lat))+
  geom_density(alpha=0.1, size=0.05, fill="red", color="red") + 
  geom_histogram(aes(y=..density..), 
                 bins=200, color="black", fill="gray", alpha=0.3, size=0.25)+
  geom_vline(xintercept = vintercept1, color="red", size=0.4)+
  theme_minimal()+
  scale_x_continuous(breaks = lat_grid) +
  annotate("text", x=vintercept1, y=-0.1, label=vintercept1,fontface=2, size=3.2) + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ggtitle("Histogram of latitude")


#### lat 그룹화 : lat_group2
vintercept2 = sort(c(47.245, 47.415, 47.557, 47.605, 47.28, 47.34, 47.465, 47.52, 47.655, 47.71, 47.755, 47.58))

lat_grid = seq(round(min(house$lat),2)-0.01, round(max(house$lat), 2), 0.05)
ggplot(house, aes(x=lat))+
  geom_density(alpha=0.1, size=0.05, fill="red", color="red") + 
  geom_histogram(aes(y=..density..), 
                 bins=200, color="black", fill="gray", alpha=0.3, size=0.25)+
  geom_vline(xintercept = vintercept2, color="red", size=0.4)+
  theme_minimal()+
  scale_x_continuous(breaks = lat_grid) +
  annotate("text", x=vintercept2, y=-0.1, label=vintercept2, fontface=2, size=3.2) + 
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  ggtitle("Histogram of latitude")


ggplot()+
  geom_point(aes(x=long, y=lat, color=log(price)), size=0.3, data=house)+
  scale_color_gradient(low="white", high="red")+
  geom_hline(yintercept = vintercept2)




#### lat_group1별 log(sqft_living) vs. log(price) 
house$lat_group1 = 1
for(i in 1:length(vintercept1)){
  house[house$lat > vintercept1[i], "lat_group1"] = i + 1
}
house$lat_group1 = factor(house$lat_group1)

ggplot()+
  geom_point(house, mapping=aes(x=log(sqft_living), y=log(price)), size=0.1) +
  facet_wrap(~lat_group1)+
  theme_minimal()

ggplot(house, mapping=aes(x=log(sqft_living), y=log(price), group=lat_group1, color=lat_group1))+
  geom_point(size=0.3, alpha=0.5) +
  stat_smooth(method = lm, se=F, lwd=0.7) +
  theme_minimal()


#### lat_group2별 log(sqft_living) vs. log(price) 
house$lat_group2 = 1
for(i in 1:length(vintercept2)){
  house[house$lat > vintercept2[i], "lat_group2"] = i + 1
}
house$lat_group2 = factor(house$lat_group2)

ggplot()+
  geom_point(house, mapping=aes(x=log(sqft_living), y=log(price)), size=0.1) +
  facet_wrap(~lat_group2)+
  theme_minimal()

ggplot(house, mapping=aes(x=log(sqft_living), y=log(price), group=lat_group2, color=lat_group2))+
  geom_point(size=0.3, alpha=0.5) +
  stat_smooth(method = lm, se=F, lwd=0.7) +
  theme_minimal()




### Regression summary 비교
lm_whole = lm(log(price) ~ log(sqft_living), data=house)
lm_lat1 = lm(log(price) ~ log(sqft_living) + lat_group1, data=house)
lm_lat2 = lm(log(price) ~ log(sqft_living) + lat_group2, data=house)
 
summary(lm_whole)
summary(lm_lat1)
summary(lm_lat2)




## Mixed model fitting
library(lme4)
model1 = lmer(log(price) ~ log(sqft_living) + (1 + log(sqft_living) | lat_group1), data = house)
model2 = lmer(log(price) ~ log(sqft_living) + (1 + log(sqft_living) | lat_group2), data = house)

summary(model1)
summary(model2)

# EDA 3
ggplot() +
  geom_point(mapping=aes(x=yr_built, y=log(price)), data = house)+
  geom_line(mapping=aes(x=yr_built, y=log(price)), data = house[house$yr_renovated>0,], color="red")
  
  
dim(house)
length(unique(house$id))

## 그룹별 도수가 2 이상인 집




