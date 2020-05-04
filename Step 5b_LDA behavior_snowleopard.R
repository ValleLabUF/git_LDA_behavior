set.seed(1)

library('MCMCpack')
library('Rcpp')
library(progress)
library(tidyverse)
library(lubridate)
library(viridis)
library(ggnewscale)
library(raster)


source('LDA_behavior_function.R')
source('gibbs sampler.R')
source('helper functions.R')
sourceCpp('aux1.cpp')




############################
#### Load and Prep Data ####
############################

#get data
dat<- read.csv('Snow Leopard Data_behav.csv', header = T, sep = ',')
dat$date<- dat$date %>% as_datetime()
dat.list<- df.to.list(dat)  #for later behavioral assignment

nbins<- c(5,8)  #number of bins per param (in order)
dat_red<- dat %>% dplyr::select(c(id, tseg, SL, TA))  #only keep necessary cols
obs<- get.summary.stats_behav(dat = dat_red, nbins = nbins)  #to run Gibbs sampler on


#prepare for Gibbs sampler
ngibbs=1000
nburn=ngibbs/2
nmaxclust=max(nbins) - 1  #one fewer than max number of bins used for params
ndata.types=length(nbins)

#prior
gamma1=0.1
alpha=0.1 

#####################################################
#### Run Gibbs Sampler on All IDs Simultaneously ####
#####################################################

res=LDA_behavior_gibbs(dat=obs, gamma1=gamma1, alpha=alpha,
                       ngibbs=ngibbs, nmaxclust=nmaxclust,
                       nburn=nburn, ndata.types=ndata.types)

#Check traceplot of log marginal likelihood
plot(res$loglikel, type='l')

#Extract and plot proportions of behaviors per time segment
theta.post<- res$theta[(nburn+1):ngibbs,]  #extract samples from posterior
theta.estim<- theta.post %>% apply(2, mean) %>% matrix(nrow(obs), nmaxclust) #calc mean of posterior
boxplot(theta.estim, xlab="Behavior", ylab="Proportion of Total Behavior")

#Determine proportion of behaviors (across all time segments)
#Possibly set threshold below which behaviors are excluded
round(apply(theta.estim, 2, sum)/nrow(theta.estim), digits = 3)



#################################################################
#### Visualize Histograms of Movement Parameters by Behavior ####
#################################################################

behav.res<- get_behav_hist(res = res, dat_red = dat_red)
behav.res<- behav.res[behav.res$behav <=3,]  #only select the top 3 behaviors

#Plot histograms of frequency data; order color scale from slow to fast
ggplot(behav.res, aes(x = bin, y = count, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Frequency\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = viridis(n=3)[3:1], guide = F) +
  facet_grid(param ~ behav, scales = "free_y")



#Plot histograms of proportion data; order color scale from slow to fast
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = viridis(n=3)[3:1], guide = F) +
  facet_grid(param ~ behav, scales = "fixed")



################################################
#### Visualize Behavior Estimates Over Time ####
################################################

#Assign behaviors (via theta) to each time segment
theta.estim<- apply(theta.estim[,1:3], 1, function(x) x/sum(x)) %>% t()  #normalize probs for only first 3 behaviors being used
theta.estim<- data.frame(id = obs$id, tseg = obs$tseg, theta.estim)
theta.estim$id<- as.character(theta.estim$id)
names(theta.estim)<- c("id", "tseg", "Exploratory", "Encamped", "Kill_Site")  #define behaviors
nobs<- data.frame(id = obs$id, tseg = obs$tseg, n = apply(obs[,3:7], 1, sum)) #calc obs per tseg using SL bins (more reliable than TA)

#Create augmented matrix by replicating rows (tsegs) according to obs per tseg
theta.estim2<- aug_behav_df(dat = dat, theta.estim = theta.estim, nobs = nobs)
ind1<- which(names(theta.estim) != "id")
theta.estim2<- theta.estim2 %>% mutate_at(names(theta.estim)[ind1], as.numeric)

#Change into long format
theta.estim.long<- theta.estim2 %>% gather(key, value, -id, -tseg, -time1, -date)
theta.estim.long$date<- theta.estim.long$date %>% as_datetime()
names(theta.estim.long)[5:6]<- c("behavior","prop")
theta.estim.long$behavior<- factor(theta.estim.long$behavior,
                                   levels = c("Encamped", "Kill_Site", "Exploratory"))



### Aligned by first observation

#stacked area
ggplot(theta.estim.long) +
  geom_area(aes(x=time1, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill") +
  labs(x = "\nObservation", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, scales = "free_x")



### Aligned by date


#stacked area
ggplot(theta.estim.long) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill") +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)



#bar plot (includes gaps in time)
ggplot(data = theta.estim.long) +
  geom_bar(aes(x=time1, y=prop, fill = behavior), size = 0.25, stat = "identity",
           position = "fill") +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, scales = "free_x")




########################################
#### Map Dominant Behavioral States ####
########################################


#Add cluster assignments to original data; one column for dominant behavior and another for prop/prob to use for alpha of points
dat2<- assign_behav(dat.list = dat.list, theta.estim2 = theta.estim2)
dat2$behav<- factor(dat2$behav, levels = c("Encamped", "Kill_Site", "Exploratory"))


#load DEM
setwd("~/Documents/Snail Kite Project/Data/snowleopards")

dem<- raster("dem/w001001.adf") %>%
  crop(extent(min(dat$x), max(dat$x), min(dat$y), max(dat$y)))  #crop to extent of points
dem.df<- raster::aggregate(dem, fact = 8, fun = mean) %>%
  as.data.frame(xy=TRUE)


# Facet plot of maps
ggplot() +
  geom_tile(data = dem.df, aes(x, y, fill = w001001)) +
  scale_fill_gradientn(name = "Elevation (m)", colours = grey(1:100/100),
                       na.value = "transparent") +
  new_scale_fill() +
  geom_path(data = dat2, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = dat2, aes(x, y, fill=behav), size=2.5, pch = 21, alpha=dat2$prop) +
  scale_fill_viridis_d("Behavior") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)










## Plot 'Kill Site' behavior only

dat2.kill<- dat2 %>% filter(behav == "Kill_Site")

setwd("~/Documents/Snail Kite Project/Data/snowleopards")


### Read in kill site shapefile
shp <- dir(getwd(), "*.shp$")
kill_sf<- st_read(shp[2])


ggplot() +
  geom_tile(data = dem.df, aes(x, y, fill = w001001)) +
  scale_fill_gradientn(name = "Elevation (m)", colours = grey(1:100/100),
                       na.value = "transparent") +
  new_scale_fill() +
  geom_path(data = dat2.kill, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = dat2.kill, aes(x, y, fill=behav), size=2.5, pch = 21, alpha=dat2.kill$prop) +
  geom_sf(data = kill_sf, color = "goldenrod", alpha = 0.7) +
  scale_fill_manual("Behavior", values = viridis(n=3)[2]) +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)


## There's some spatial overlap between "Kill Site" behavior and actual locations of kill sites in space for Khani and Pahlawan (and maybe Pari). Now need to verify with temporal component

kill_sf$Date<- as.POSIXct(as.character(kill_sf$Date), format = "%d-%b-%y", tz = "UTC")


ggplot(theta.estim.long) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill") +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  geom_vline(xintercept = kill_sf$Date, color = viridis(n=9)[7]) +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, scales = "free_x")





### Seems like Khani and Pahlawan overlap temporally with kill sites; time to bring it together for a few key examples


# 2nd reported kill site
kill_2<- dat2 %>%
  filter(date >= (kill_sf$Date[2] - ddays(7)) & date <= (kill_sf$Date[2] + ddays(1))) %>%
  filter(id == "Khani1")

ggplot() +
  geom_path(data = kill_2, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = kill_2, aes(x, y, fill=behav), size=2.5, pch = 21) +
  geom_sf(data = kill_sf[2,], color = "red", alpha = 0.7) +
  scale_fill_viridis_d("Behavior") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)

#dominant behavior is 'Encamped', but based on time series behavior plot, time segment is nearly evenly split between 'encamped' and 'kill site' behaviors



# 3rd reported kill site
kill_3<- dat2 %>%
  filter(date >= (kill_sf$Date[3] - ddays(4)) & date <= kill_sf$Date[3]) %>%
  filter(id == "Khani1")

ggplot() +
  geom_path(data = kill_3, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = kill_3, aes(x, y, fill=behav), size=2.5, pch = 21) +
  geom_sf(data = kill_sf[3,], color = "red", alpha = 0.7) +
  scale_fill_viridis_d("Behavior") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)


# dominant behavior is "Kill Site", which clearly fits with known time and location of carcass for Khani



# 4th reported kill site
kill_4<- dat2 %>%
  filter(date >= (kill_sf$Date[4] - ddays(5)) & date <= kill_sf$Date[4]) %>%
  filter(id == "Khani1" | id == "Pahlawan")

ggplot() +
  geom_path(data = kill_4, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = kill_4, aes(x, y, fill=behav), size=2.5, pch = 21) +
  geom_sf(data = kill_sf[4,], color = "red", alpha = 0.7) +
  scale_fill_viridis_d("Behavior") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)

# Khani, but not Pahlawan, associated with this carcass; "Encamped" is dominant behavior





# 5th reported kill site
kill_5<- dat2 %>%
  filter(date >= (kill_sf$Date[5] - ddays(7)) & date <= (kill_sf$Date[5] + ddays(1))) %>%
  filter(id == "Khani1")

ggplot() +
  geom_path(data = kill_5, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = kill_5, aes(x, y, fill=behav), size=2.5, pch = 21) +
  geom_sf(data = kill_sf[5,], color = "red", alpha = 0.7) +
  scale_fill_viridis_d("Behavior") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)

# Khani also associated with this kill site, but also showing "Encamped" as dominant behavior




# 6th reported kill site
kill_6<- dat2 %>%
  filter(date >= (kill_sf$Date[6] - ddays(7)) & date <= kill_sf$Date[6]) %>%
  filter(id == "Pahlawan")

ggplot() +
  geom_path(data = kill_6, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = kill_6, aes(x, y, fill=behav), size=2.5, pch = 21) +
  geom_sf(data = kill_sf[6,], color = "red", alpha = 0.7) +
  scale_fill_viridis_d("Behavior") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)

#Pahlawan associated with kill site, but "Encamped" is dominant behavior





# 7th reported kill site
kill_7<- dat2 %>%
  filter(date >= (kill_sf$Date[7] - ddays(6)) & date <= kill_sf$Date[7])

ggplot() +
  geom_path(data = kill_7 %>% filter(id=="Khani1"), aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = kill_7 %>% filter(id=="Khani1"), aes(x, y, fill=behav), size=2.5, pch = 21) +
  geom_sf(data = kill_sf[7,], color = "red", alpha = 0.7) +
  scale_fill_viridis_d("Behavior") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)


# Kill associated with Khani, although dominant behavior is "Exploratory"; this appears to be result of very large time segment where most movements are "Exploratory"




# 10th reported kill site
kill_10<- dat2 %>%
  filter(date >= (kill_sf$Date[10] - ddays(6)) & date <= kill_sf$Date[10])

ggplot() +
  geom_path(data = kill_10 %>% filter(id=="Pahlawan"), aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = kill_10 %>% filter(id=="Pahlawan"), aes(x, y, fill=behav), size=2.5, pch = 21) +
  geom_sf(data = kill_sf[10,], color = "red", alpha = 0.7) +
  scale_fill_viridis_d("Behavior") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)

#Pahlawan associated with kill site, but in "Encamped" behavior as dominant state


### Overall, it appears that most kill sites are associated with "Encamped" as dominant behavior, but this may be due to long time series from the segmentation model and/or the "Kill Site" behavior as a minor subset of "Encamped" behavior that never really dominates a time segment, but obviously is still present during certain periods
