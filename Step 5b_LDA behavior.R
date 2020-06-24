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
library(wesanderson)


source('LDA_behavior_function.R')
source('gibbs sampler.R')
source('helper functions.R')
sourceCpp('aux1.cpp')




############################
#### Load and Prep Data ####
############################

#get data
dat<- read.csv('Snail Kite Gridded Data_TOHO_behav2.csv', header = T, sep = ',')
dat$date<- dat$date %>% as_datetime()
dat.list<- df.to.list(dat, ind = "id")  #for later behavioral assignment

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
theta.estim_df<- theta.estim %>% 
  as.data.frame() %>% 
  pivot_longer(., cols = 1:7, names_to = "behavior", values_to = "prop") %>% 
  modify_at("behavior", factor)
levels(theta.estim_df$behavior)<- 1:7

ggplot(theta.estim_df, aes(behavior, prop)) +
  geom_boxplot(aes(fill = behavior), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = behavior), position = position_jitter(0.1), 
              alpha = 0.3) +
  scale_color_viridis_d("", guide = F) +
  scale_fill_viridis_d("", guide = F) +
  labs(x="\nBehavior", y="Proportion of Time Segment\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))

# ggsave("Figure 5b (behavior boxplot).png", width = 6, height = 4, units = "in",
#        dpi = 330)


#Determine proportion of behaviors (across all time segments)
#Possibly set threshold below which behaviors are excluded
round(apply(theta.estim, 2, sum)/nrow(theta.estim), digits = 3)



#################################################################
#### Visualize Histograms of Movement Parameters by Behavior ####
#################################################################

behav.res<- get_behav_hist(dat = res, nburn = nburn, ngibbs = ngibbs, nmaxclust = nmaxclust,
                           var.names = c("Step Length","Turning Angle"))
behav.res<- behav.res[behav.res$behav <=3,]  #only select the top 3 behaviors
behav.res$param<- factor(behav.res$param)
behav.res$behav<- factor(behav.res$behav)
levels(behav.res$behav)<- c("Encamped","ARS","Transit")

#Plot histograms of proportion data; order color scale from slow to fast
ggplot(behav.res, aes(x = bin, y = prop, fill = as.factor(behav))) +
  geom_bar(stat = 'identity') +
  labs(x = "\nBin", y = "Proportion\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.text.x = element_text(face = "bold"),
        panel.grid = element_blank()) +
  scale_fill_manual(values = viridis(n=3), guide = F) +
  scale_y_continuous(breaks = c(0.00, 0.50, 1.00), limits = c(0,1.1)) +
  scale_x_continuous(breaks = 1:8) +
  facet_grid(behav ~ param, scales = "free_x")

# ggsave("Figure 5c (behavior histograms).png", width = 7, height = 5, units = "in",
#        dpi = 330)


################################################
#### Visualize Behavior Estimates Over Time ####
################################################

#Assign behaviors (via theta) to each time segment
theta.estim<- apply(theta.estim[,1:3], 1, function(x) x/sum(x)) %>% t()  #normalize probs for only first 3 behaviors being used
theta.estim<- data.frame(id = obs$id, tseg = obs$tseg, theta.estim)
theta.estim$id<- as.character(theta.estim$id)
names(theta.estim)<- c("id", "tseg", "Encamped", "ARS", "Transit")  #define behaviors
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
                                   levels = c("Encamped", "ARS", "Transit"))



### Aligned by first observation

#stacked area
ggplot(theta.estim.long) +
  geom_area(aes(x=time1, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill") +
  labs(x = "\nTime", y = "Proportion of Behavior\n") +
  scale_fill_viridis_d("") +
  scale_x_continuous(guide = guide_axis(check.overlap = TRUE)) +
  scale_y_continuous(breaks = c(0,0.5,1)) +
  theme_bw() +
  theme(axis.title = element_text(size = 16), axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top") +
  facet_wrap(~id, scales = "free_x")

# ggsave("Figure S1 (behavior prop time series_all).png", width = 8, height = 6, units = "in",
#        dpi = 330)


### Aligned by date

#Window of peak breeding (March 1 - June 30)
# breed<- data.frame(xmin = as_datetime(c("2016-03-01 00:00:00","2017-03-01 00:00:00",
#                                         "2018-03-01 00:00:00","2019-03-01 00:00:00")),
#                    xmax = as_datetime(c("2016-06-30 23:59:59","2017-06-30 23:59:59",
#                                         "2018-06-30 23:59:59","2019-06-30 23:59:59")),
#                    ymin = -Inf, ymax = Inf)

breed<- data.frame(xmin = as_datetime(c("2018-03-01 00:00:00","2019-03-01 00:00:00")),
                   xmax = as_datetime(c("2018-06-30 23:59:59","2019-06-30 23:59:59")),
                   ymin = -Inf, ymax = Inf)



#stacked area
ggplot(theta.estim.long %>% filter(id == "SNIK 12" | id == "SNIK 14" | id == "SNIK 15")) +
  geom_rect(data = breed, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.5) +
  geom_area(aes(x=date, y=prop, fill = behavior), color = "black", size = 0.25,
            position = "fill", alpha = 0.7) +
  labs(x = "Time", y = "Proportion of Time Segment") +
  scale_y_continuous(breaks = c(0,0.5,1), expand = c(0.05,0.05)) +
  scale_fill_viridis_d("") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x.bottom = element_text(size = 12),
        strip.text = element_text(size = 12, face = "bold"),
        panel.grid = element_blank(),
        legend.position = "n") +
  facet_wrap(~id, nrow = 3)

# ggsave("Figure 6b (behavior prop time series).png", width = 7, height = 5, units = "in",
#        dpi = 330)



#bar plot (includes gaps in time)
ggplot(theta.estim.long) +
  geom_rect(data = breed, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey", alpha = 0.5) +
  geom_bar(aes(x=date, y=prop, fill = behavior), size = 0.25, stat = "identity",
           position = "stack", width = 3600) +
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
dat2<- assign_behav(dat.list = dat.list, theta.estim.long = theta.estim.long,
                    behav.names = c("Encamped","ARS","Transit"))
dat2$behav<- factor(dat2$behav, levels = c("Encamped", "ARS", "Transit"))


#load map data
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida") %>% st_transform(fl, crs = "+init=epsg:32617")

# lakes
lakes10 <- ne_download(scale = 10, type = 'lakes', category = 'physical', returnclass = "sf")
lakes10<- sf::st_transform(lakes10, crs = "+init=epsg:32617") %>%
  sf::st_crop(xmin = min(dat$x-20000), xmax = max(dat$x+20000), ymin = min(dat$y-20000),
              ymax = max(dat$y+20000))

# Facet plot of maps
dat2_focal<- dat2 %>% filter(id == "SNIK 12" | id == "SNIK 14" | id == "SNIK 15")


## Load hydro data
library(lwgeom)
library(nhdplusTools)
setwd("~/Documents/NHD Seamless Geodatabase/NHDPlusNationalData")

download_dir<- download_nhdplushr(nhd_dir = getwd(), 
                                  hu_list = c("0307","0308","0309","0310","0311","0312","0313",
                                              "0314"), download_files = TRUE)
waterbody<- read_sf(file.path(download_dir, "nhd_hr_FL.gpkg"), "NHDWaterbody") %>%
  st_cast("MULTIPOLYGON")
waterbody2<- waterbody %>%
  filter(AreaSqKM > 50) %>%
  st_simplify(preserveTopology = TRUE, dTolerance = 250)

#mask flowline and waterbody to FL layer
waterbody_fl<- st_intersection(st_make_valid(waterbody2), fl) %>% 
  filter(GNIS_Name != "The Everglades")


ggplot() +
  geom_sf(data = fl, fill = "grey45", color = NA) +
  geom_sf(data = waterbody_fl, fill = "white", color = NA) +
  coord_sf(xlim = c(min(dat$x-60000), max(dat$x+60000)),
           ylim = c(min(dat$y-20000), max(dat$y+10000)), expand = FALSE) +
  geom_path(data = dat2_focal, aes(x=x, y=y), color="gray60", size=0.25) +
  geom_point(data = dat2_focal, aes(x, y, fill=behav), size=2.5, pch=21, alpha=0.7) +
  geom_point(data = dat2_focal %>%
               group_by(id) %>%
               slice(which(row_number() == 1)) %>%
               ungroup(), aes(x, y), color = "green", pch = 21, size = 3,
             stroke = 1.25) +
  geom_point(data = dat2_focal %>%
               group_by(id) %>%
               slice(which(row_number() == n())) %>%
               ungroup(), aes(x, y), color = "red", pch = 24, size = 3,
             stroke = 1.25) +
  scale_fill_viridis_d("") +
  scale_x_continuous(breaks = c(-83:-80)) +
  scale_y_continuous(breaks = c(26:29)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid = element_blank(),
        legend.position = "top") +
  guides(fill = guide_legend(label.theme = element_text(size = 12),
                             title.theme = element_text(size = 14))) +
  facet_wrap(~id)

# ggsave("Figure 6a (maps pf focal IDs).png", width = 7, height = 5, units = "in",
#        dpi = 330)



#Calculate activity budget (by obs proportions)
# activ.budg<- theta.estim.long %>% 
#   group_by(id, behavior) %>% 
#   summarise(mean=mean(prop), min=min(date, na.rm = T), max=max(date, na.rm = T),
#             duration=difftime(max,min, units = "days")) %>% 
#   ungroup() %>% 
#   dplyr::select(-c(min, max)) %>% 
#   modify_at("duration", as.numeric)

#Separate by stages for longest tracks
activ.budg<- theta.estim.long %>% 
  filter(id == "SNIK 12" | id == "SNIK 14" | id == "SNIK 15") %>% 
  group_by(id, behavior) %>%
  mutate(stage = case_when(date <= (date[1] + months(6)) ~ "Juvenile",
                           date > (date[1] + months(6)) ~ "Adult")) %>% 
  group_by(id, stage, behavior) %>% 
  summarise(mean=mean(prop)) %>% 
  modify_at("stage", ~factor(., levels = c("Juvenile","Adult"))) %>% 
  ungroup()


ggplot(activ.budg) +
  geom_path(aes(stage, mean, group = id, color = id)) +
  geom_point(aes(stage, mean, color = id), size = 3, alpha = 0.7) +
  facet_wrap(~ behavior) +
  scale_color_manual("", values = wes_palette("Zissou1", n=5, type = "discrete")[c(1,3,5)]) +
  labs(x=NULL, y="Mean Proportion of Behavior\n") +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "bottom",
        panel.grid = element_blank(),
        strip.text = element_text(size = 12, face = "bold"))

# ggsave("Figure 7 (activity budget line plot).png", width = 6, height = 4, units = "in",
#        dpi = 330)





#activity budget by month
activ.budg2<- dat2_focal %>% 
  pivot_longer(., cols = c(Encamped, ARS, Transit), names_to = "behavior",values_to = "prop") %>%
  dplyr::select(id, date, behavior, prop) %>% 
  mutate(yr_mo = paste(year(date), month(date, label = TRUE, abbr = F), 1, sep = "-")) %>% 
  group_by(id, yr_mo, behavior) %>%
  summarise(mean = mean(prop)) %>% 
  ungroup() %>% 
  mutate_at("yr_mo", ~as.POSIXct(., format = "%Y-%B-%d")) %>% 
  arrange(id, yr_mo) %>% 
  mutate_at("behavior", ~factor(., levels = c("Encamped", "ARS", "Transit")))


ggplot(activ.budg2, aes(as.Date(yr_mo), mean, fill = behavior)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis_d("") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(date_labels = "%b %Y") +
  # geom_smooth(aes(color=behavior), se=F, method="loess", show.legend = FALSE,
  #             lwd=0.9) +
  scale_color_viridis_d("") +
  labs(x = "\nDate", y = "Proportion of Time\n") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.position = "top",
        plot.margin = margin(t=0.5, r=1, b=0.5, l=0.5, unit = "cm")) +
  facet_wrap(~ id, scales = "fixed")

# ggsave("Figure 7 (activity budget bar plot).png", width = 12, height = 8, units = "in",
#        dpi = 330)



#Export DF

# write.csv(dat2, "Snail Kite Gridded Data_TOHO_behavior.csv", row.names = F)