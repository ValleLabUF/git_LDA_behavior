set.seed(1)

library('MCMCpack')
library('Rcpp')
library(progress)
library(tidyverse)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(viridis)


source('LDA_behavior_function.R')
source('gibbs sampler.R')
source('helper functions.R')
sourceCpp('aux1.cpp')




############################
#### Load and Prep Data ####
############################

#get data
dat<- read.csv('Armadillo Data_behav.csv', header = T, sep = ',')
dat$date<- dat$date %>% as_datetime()
dat$TA<- ifelse(!is.na(dat$TA), dat$TA, 9)  #need to assign bin to NAs
dat.list<- df.to.list(dat, ind = "id")  #for later behavioral assignment

nbins<- c(5,9)  #number of bins per param (in order)
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
theta.estim_df<- theta.estim %>% 
  as.data.frame() %>% 
  pivot_longer(., cols = 1:8, names_to = "behavior", values_to = "prop") %>% 
  modify_at("behavior", factor)
levels(theta.estim_df$behavior)<- 1:8

ggplot(theta.estim_df, aes(behavior, prop)) +
  geom_boxplot(aes(fill = behavior), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = behavior), position = position_jitter(0.1), 
              alpha = 0.3) +
  scale_color_viridis_d("", guide = F) +
  scale_fill_viridis_d("", guide = F) +
  labs(x="\nBehavior", y="Proportion of Total Behavior\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))


#Determine proportion of behaviors (across all time segments)
#Possibly set threshold below which behaviors are excluded
round(apply(theta.estim, 2, sum)/nrow(theta.estim), digits = 3)



#################################################################
#### Visualize Histograms of Movement Parameters by Behavior ####
#################################################################

behav.res<- get_behav_hist(res = res, dat_red = dat_red)
behav.res<- behav.res[behav.res$behav <=3,]  #only select the top 3 behaviors
behav.res$param<- str_replace_all(behav.res$param, "SL", "Step Length") %>% 
  str_replace_all(., "TA", "Turning Angle")
behav.res$behav<- behav.res$behav %>% 
  as.character() %>% 
  str_replace_all(., "1", "Foraging") %>% 
  str_replace_all(., "2", "Burrow") %>% 
  str_replace_all(., "3", "Transit")

#Plot histograms of frequency data; order color scale from slow to fast
ggplot(behav.res, aes(x = bin, y = count, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Frequency\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14), strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  facet_grid(behav ~ param, scales = "free")



#Plot histograms of proportion data; order color scale from slow to fast
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold")) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00)) +
  scale_x_continuous(breaks = 1:9) +
  facet_grid(behav ~ param, scales = "free_x")



################################################
#### Visualize Behavior Estimates Over Time ####
################################################

#Assign behaviors (via theta) to each time segment
theta.estim<- apply(theta.estim[,1:3], 1, function(x) x/sum(x)) %>% t()  #normalize probs for only first 3 behaviors being used
theta.estim<- data.frame(id = obs$id, tseg = obs$tseg, theta.estim)
theta.estim$id<- as.character(theta.estim$id)
names(theta.estim)<- c("id", "tseg", "Foraging", "Burrow", "Transit")  #define behaviors
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
                                   levels = c("Burrow", "Foraging", "Transit"))



### Aligned by first observation

#stacked area
ggplot(theta.estim.long %>% filter(id == "tm15" | id == "tm30")) +
  geom_area(aes(x=time1, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill") +
  labs(x = "\nObservation", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id, scales = "free_x", nrow = 2)



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
ggplot(theta.estim.long) +
  geom_bar(aes(x=time1, y=prop, fill = behavior), size = 0.25, stat = "identity",
           position = "stack", width = 300) +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("Behavior") +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  facet_wrap(~id)




########################################
#### Map Dominant Behavioral States ####
########################################


#Add cluster assignments to original data; one column for dominant behavior and another for prop/prob to use for alpha of points
dat2<- assign_behav(dat.list = dat.list, theta.estim2 = theta.estim2)
dat2$behav<- factor(dat2$behav, levels = c("Burrow", "Foraging", "Transit"))



# Facet plot of maps
ggplot() +
  geom_path(data = dat2 %>% filter(id == "tm15" | id == "tm30"), aes(x=x, y=y),
            color="gray60", size=0.25) +
  geom_point(data = dat2 %>% filter(id == "tm15" | id == "tm30"), aes(x, y, fill=behav),
             size=1.5, pch=21, alpha=0.7) +
  scale_fill_viridis_d("Behavior") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank()) +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(~id, scales = "free")



#Export DF

# write.csv(dat2, "Armadillo Data_behavior.csv", row.names = F)
