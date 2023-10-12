### EGD Collective Behavior Code (Clean-up; streamline) ###

# packages needed for analysis #
library(tidyverse)
library(ggplot2)
library(readxl)
library(RColorBrewer)
library(ggpmisc)
library(ggpubr)
library(dplyr)
library(ggExtra)
library(MASS)
library(ggpointdensity)
library(tidyr)
library(patchwork)
library(png)
library(spatstat)
library(broom)
library(stringr)
library(ggbreak)
library(report)
library(sjPlot)
library(performance)
library(jtools)
library(emmeans)
library(multcomp)
library(multcompView)
library(scales)
library(GGally)

## Main Paper ##
set.seed(12345)

##### datasets #####
## annotated tank dataset for distance analyses
old_dat <- read.csv("~/Desktop/ESM_Hensley_etal_2023/data/EGD_XY_26apr23.csv",header = TRUE)
old_dat2 <- old_dat %>% dplyr::select(-c(stim_start,stim_end))
new_dat <- read.csv("~/Desktop/ESM_Hensley_etal_2023/data/EGD_XY_24apr23.csv",header = TRUE)
alldat <- rbind(old_dat2,new_dat)
alldat$eu_dist <- sqrt((alldat$x_pos)^2 + (alldat$y_pos)^2) #distance from origin
xydat <- alldat %>% mutate(Observation = interaction(date,tank,stim_duration,order)) %>% group_by(Observation) %>%
  dplyr::mutate(new_time = (time - min(time)))

##attempting to do a nested dist matrix analysis on the ex situ observations (duration = 0) only
ex_situ <- alldat %>% filter(stim_duration == 0) %>% dplyr::mutate(date_tank = interaction(date,tank)) %>% #ungroup() %>%
  dplyr::select(c(x_pos,y_pos,date_tank)) %>%
  nest_by(date_tank)

nearest_neighbor_prep <- function(data_set,index){
  input <- as.character(data_set$date_tank[[index]])
  data_date <- str_split_i(input,pattern = "\\.",i = 1)
  data_tank <- str_split_i(input,pattern = "\\.",i = 2)

  d <- data_set$data[[index]]
  nn <- nndist(d,k = 1)
  np <- xydat %>% filter(date == data_date,tank == data_tank,stim_duration == 0) %>% ungroup() %>% dplyr::select(num_pulse,new_time)
  nrdp <- data.frame(np,nearest_neighbor = nn)

  return(nrdp)
}

nn1 <- nearest_neighbor_prep(ex_situ,1)
nn2 <- nearest_neighbor_prep(ex_situ,2)
nn3 <- nearest_neighbor_prep(ex_situ,3)
nn4 <- nearest_neighbor_prep(ex_situ,4)

nn_np_all <- rbind(nn1,nn2,nn3,nn4)

ex_situ_t <- alldat %>% filter(stim_duration == 0) %>% dplyr::mutate(date_tank = interaction(date,tank)) %>% #ungroup() %>%
  dplyr::select(c(time,date_tank)) %>%
  nest_by(date_tank)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

nearest_time_prep <- function(data_set1,data_set2,index){
  input <- as.character(data_set1$date_tank[[index]])
  data_date <- str_split_i(input,pattern = "\\.",i = 1)
  data_tank <- str_split_i(input,pattern = "\\.",i = 2)

  d <- data_set1$data[[index]]
  tD <- dist(d,method = "euclidean",diag = FALSE,upper = FALSE)
  m <- as.matrix(tD)
  xy <- t(combn(colnames(m), 2))
  dm <- data.frame(xy, dist=m[xy])

  d_t <- data_set2$data[[index]]
  tDt <- dist(d_t,method = "maximum",diag = FALSE,upper = FALSE)
  mt <- as.matrix(tDt)
  xyt <- t(combn(colnames(mt), 2))
  dmm <- data.frame(xyt, time=mt[xyt])

  xxyy <- left_join(dm,dmm)
  xxyy <- xxyy %>% filter(time <= 32)

  xxyy$date_tank <- rep(input,nrow(xxyy))

  return(xxyy)
}

np1 <- nearest_time_prep(ex_situ,ex_situ_t,1)
np2 <- nearest_time_prep(ex_situ,ex_situ_t,2)
np3 <- nearest_time_prep(ex_situ,ex_situ_t,3)
np4 <- nearest_time_prep(ex_situ,ex_situ_t,4)

nn_nt_all <- rbind(np1,np2,np3,np4)
nn_nt_all$density <- get_density(nn_nt_all$time, nn_nt_all$dist, n = 100)

## manual observation datasets from 2017
egd_dat = read.csv("~/Desktop/ESM_Hensley_etal_2023/data/EGD_duration_intensity_trials.csv",header=TRUE)
egd_dat <- egd_dat %>% filter(include == 1)
dur = subset.data.frame(egd_dat,egd_dat$Experiment == "Duration")
int = subset.data.frame(egd_dat,egd_dat$Experiment == "Intensity")

## light stimulus post-hoc measurements
stim_dat <- readxl::read_xlsx("~/Desktop/ESM_Hensley_etal_2023/data/tod_results.xlsx")

## field survey data
dat1 <- read_xls("~/Desktop/ESM_Hensley_etal_2023/data/EGD_Field_surverys_data.xls",
                 sheet = 1,col_names = TRUE)
dat2 <- read_xls("~/Desktop/ESM_Hensley_etal_2023/data/EGD_2022_field_survey.xls",
                 sheet = 1,col_names = TRUE)
d1 <- dat1 %>% group_by(Date) %>% dplyr::mutate(year = "2018",dt = Time - lag(Time,n = 1))
d2 <- dat2 %>% group_by(Date) %>% dplyr::mutate(year = "2022",dt = Time - lag(Time,n = 1))
d_survey <- rbind(d1,d2)

## EGD display trait data
disp_dat <- read.csv("~/Desktop/ESM_Hensley_etal_2023/data/egd_display_data_raw.csv",header=TRUE)
bright_dat <- read.csv("~/Desktop/ESM_Hensley_etal_2023/data/brightness_analysis.csv",header=TRUE)

bright_dat <- bright_dat %>% arrange(display_num,pulse_num,frame) %>%
  mutate(brightness = (max_pix_value/roi_area) - (background_max/background_area)) %>%
  group_by(display_num) %>%
  mutate(p1_bright = max(brightness[pulse_num == 1]),
         p1_min = min(brightness[pulse_num == 1]),
         relative_brightness = brightness / p1_bright,
         display_num = display_num + max(disp_dat$display_num)) %>% ungroup() %>%
  group_by(display_num,pulse_num) %>%
  dplyr::summarise(display_num = mean(display_num),
                   pulse_num = mean(pulse_num),
                   relative_brightness = max(relative_brightness)) %>% ungroup()

sub_bright <- subset.data.frame(bright_dat,select = c("display_num","pulse_num","relative_brightness"))
sub_bright$time <- as.numeric(NA)
sub_bright$ipi_num <- as.integer(NA)
sub_bright$ipi_time_s <- as.numeric(NA)
sub_bright$ipd_distance_mm <- as.numeric(NA)
sub_bright$vert_dist_num <- as.numeric(NA)
sub_bright$source <- "NMH"

sub_bright <- sub_bright %>% dplyr::select(display_num,pulse_num,time,ipi_num,ipi_time_s,ipd_distance_mm,vert_dist_num,relative_brightness,source)
names(sub_bright) <- colnames(disp_dat)

disp_dat <- rbind(disp_dat,sub_bright)

## video time series data
timeser <- read.csv("~/Desktop/ESM_Hensley_etal_2023/data/09_27_18_right_2_beautified.csv",header=FALSE,col.names = c("frame","brightness","adj_brightness"))

t1 <- timeser %>% dplyr::summarise(median = median(adj_brightness),
                                   min = min(adj_brightness),
                                   max = max(adj_brightness),
                                   sd = sd(adj_brightness),
                                   q1 = median - sd, q3 = median + 3*sd) %>% gather(key = "measurement",value = "value")

fr_med <- median(timeser$frame[timeser$adj_brightness == t1[1,2]])
fr_low <- median(which(abs(timeser$adj_brightness - t1[6,2]) == min(abs(timeser$adj_brightness - t1[6,2]))))
fr_hi <- median(timeser$frame[timeser$adj_brightness == t1[3,2]])

dab <- data.frame(frame = c(fr_med,fr_low,fr_hi),
                  adj_b = c(timeser$adj_brightness[timeser$frame == fr_med],
                            timeser$adj_brightness[timeser$frame == fr_low],
                            timeser$adj_brightness[timeser$frame == fr_hi]))

p_med <- readPNG("~/Desktop/ESM_Hensley_etal_2023/data/23783_ed2.png",native = TRUE)
p_low <- readPNG("~/Desktop/ESM_Hensley_etal_2023/data/40297_ed2.png",native = TRUE)
p_hi <- readPNG("~/Desktop/ESM_Hensley_etal_2023/data/8819_ed2.png",native = TRUE)

## Photeros sp. display data
p_disp <- readPNG("~/Desktop/ESM_Hensley_etal_2023/data/output_238.png",native = TRUE)

phodat <- read.csv("~/Desktop/ESM_Hensley_etal_2023/data/photeros_data.csv",header=TRUE)
phodat2 <- phodat %>% group_by(species,pulse_num) %>% mutate(duration = lead(x), Direction = case_when(direction == "down" ~ "Downwards",
                                                                                                 direction == "up" ~ "Upwards"),
                                                       Species = case_when(species == "egd" ~ "P. sp. 'EGD'",
                                                                           species == "annecohenae" ~ "P. annecohenae",
                                                                           species == "graminicola" ~ "P. graminicola",
                                                                           species == "jamescasei" ~ "P. jamescasei",
                                                                           species == "johnbucki" ~ "P. johnbucki",
                                                                           species == "mcelroyi" ~ "P. mcelroyi",
                                                                           species == "morini" ~ "P. morini",
                                                                           species == "shulmanae" ~ "P. shulmanae"
                                                       ),
                                                       Species = factor(Species,levels = c("P. annecohenae","P. graminicola","P. jamescasei","P. mcelroyi","P. johnbucki","P. morini","P. shulmanae","P. sp. 'EGD'")),
) %>% ungroup()

phodat1 <- phodat2 %>% filter(row_number() %% 2 == 1) %>% mutate(y = case_when(Direction == "Downwards" ~ y*-1,
                                                                           Direction == "Upwards" ~ y))

### figure 5B/C data prep
old_dat3 <- old_dat %>% mutate(signal_start_delay = (frame - stim_start)/30, signal_end_delay = (frame - stim_end)/30,
                   eu_dist = sqrt((x_pos)^2 + (y_pos)^2),
                   signal_onset = case_when(signal_start_delay < 0 & signal_end_delay < 0 ~ "Before",
                                            signal_start_delay >= 0 & signal_end_delay <= 0 ~ "During",
                                            signal_start_delay >= 0 & signal_end_delay > 0 ~ "After"),
                   Observation = interaction(date,tank,order,stim_duration)) %>%
  pivot_longer(c(signal_start_delay,signal_end_delay),names_to = "signal_timing",values_to = "delay") %>%
  dplyr::group_by(Observation) %>%
  mutate(new_time = time - min(time)) %>% ungroup()

##### Analysis and Figures #####
##New Figure 1
fig1new <- (fig1a + fig2) + fig3a + fig3b + plot_layout(nrow=1,widths = c(0.25,1,0.65,0.65))#+ plot_annotation(tag_levels = c("A")) & theme(plot.tag = element_text(face = 'bold',size = 50))
ggsave(filename = "fig1new.jpeg",plot = fig1new,device = "jpeg",path = "~/Desktop/ESM_Hensley_etal_2023/figures/",dpi = 300,units = "px",height = 2900,width = 13500)

#Figure S1
fig1a <- wrap_elements(p_disp)
fig1b <- phodat2 %>% mutate(y = case_when(Direction == "Downwards" ~ y*-1,
                                 Direction == "Upwards" ~ y)) %>%
  ggplot(aes(x=x,y=y,color=Species)) +
  geom_point(shape=19,size=2.5,stroke=1.5) +
  geom_segment(aes(x=x,y=y,xend=duration,yend=lead(y)),linewidth=1.25) +
  geom_line(data=phodat1,aes(x=x,y=y),linetype=2,linewidth = 1.25) +
  xlab("Time (s)") + ylab("Distance (cm)") + theme_minimal(base_size = 35) +
  scale_color_manual(values = c("P. sp. 'EGD'" = "#0072B2",
                                "P. annecohenae" = "#D55E00",
                                "P. graminicola" = "#F57825",
                                "P. jamescasei" = "#FF9340",
                                "P. johnbucki" = "#7AC1FF",
                                "P. mcelroyi" = "#FFAE5A",
                                "P. morini" = "#5BA6EA",
                                "P. shulmanae" = "#3A8BCE"
                                )) +
  theme(legend.text = element_text(face="italic"),
        legend.title = element_text(face="bold"),
        axis.title = element_text(face="bold"),
        strip.text = element_text(face = "bold"),
        axis.line.x.bottom = element_line(linewidth = 1),
        axis.line.y = element_line(linewidth = 1),
        legend.position = "right") + ylim(c(-150,200))

fig1 <- fig1a + fig1b + plot_layout(nrow=1,widths = c(1.5,2))
ggsave(filename = "fig1.jpeg",plot = fig1,device = "jpeg",path = "~/Desktop/ESM_Hensley_etal_2023/figures/",dpi = 300,units = "px",height = 2500,width = 6000)

#Figure 1
ppp <- wrap_elements(panel = p_hi) + wrap_elements(panel = p_med) + wrap_elements(panel = p_low) +
  plot_layout(ncol = 3,widths = 1)

pd <- timeser %>% ggplot(aes(x=frame/(30),y=adj_brightness)) + geom_line(linewidth=1.5,color="#56B4E9") +
  geom_point(data=dab,aes(x=frame/(30),y=adj_b),pch=21,size=4,stroke=2.5,color="#000000") +
  theme_bw() + xlab("Time (s)") + ylab("Adjusted frame brightness") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text = element_text(size=20),
        axis.title = element_text(size=15,face="bold"))

fig2 <- ggarrange(ppp,pd,nrow=2)
ggsave(filename = "fig2.jpeg",plot = fig2,device = "jpeg",path = "~/Desktop/ESM_Hensley_etal_2023/figures/",dpi = 300,units = "px",height = 3300,width = 5000)

# Table 1
disp_dat %>% dplyr::group_by(display_num) %>%
  tally() %>% dplyr::summarise(mean_num = mean(n),
                               se_num = sd(n)/sqrt(n()))

disp_dat %>% group_by(display_num) %>% dplyr::summarise(max_dist = max(vert_dist_mm)) %>%
  filter(!is.na(max_dist)) %>%
  dplyr::summarise(mean_dist = mean(max_dist,na.rm=TRUE),se_dist=sd(max_dist,na.rm=TRUE)/sqrt(n()),max_dist = max(max_dist),n=n())

disp_dat %>% dplyr::mutate(duration = time*1000, ipi = ipi_time_s * 1000) %>% group_by(display_num) %>%
  dplyr::summarise(tot_dur = sum(duration,na.rm = TRUE) - sum(ipi,na.rm = TRUE)) %>%
  filter(tot_dur > 0) %>%
  dplyr::summarise(mean_td = mean(tot_dur),max_td = max(tot_dur),se_td = sd(tot_dur)/sqrt(n()),n=n())

disp_dat %>% dplyr::mutate(duration = time*1000, ipi = ipi_time_s * 1000) %>%
  dplyr::select(display_num,pulse_num,duration,ipi,ipd_distance_mm,vert_dist_mm,relative_brightness) %>%
  gather(key = "trait",value = "measurement",duration,ipi,ipd_distance_mm,vert_dist_mm,relative_brightness) %>%
  filter(!is.na(measurement)) %>%
  group_by(trait,pulse_num) %>%
  dplyr::summarise(mean = mean(measurement,na.rm=TRUE),
                   sd = sd(measurement,na.rm=TRUE),
                   se = sd(measurement,na.rm=TRUE)/sqrt(n()),
                   n = n()) %>% print(n = 29)

# Figure 1
gsum <- d_survey %>% group_by(year,Date) %>% dplyr::summarise(ndisp = n(),dafter = sum(Days_after_fullmoon)/n()) %>%
  dplyr::mutate(n = as.numeric(ndisp), Days_after_fullmoon = as.numeric(dafter), year = as.factor(year))

dd <- d_survey %>% dplyr::mutate(min_after_sunset = Time - Sunset,dt = as.numeric(dt))
m1 <- lm(dt ~ Days_after_fullmoon + year + min_after_sunset,data=dd)
m2 <- lm(dt ~ Days_after_fullmoon + year*min_after_sunset,data=dd)
anova(m1,m2)
summary(m2); qqnorm(resid(m2))
report(m2)
check_model(m2)
tab_model(m2,show.stat = TRUE,show.df = TRUE,dv.labels = c("Time (s) between collective events"),
          pred.labels = c("(Intercept)","# days after full moon","Year [2022]", "# minutes post sunset","Year x # min. post sunset"))

fig3a <- d_survey %>% mutate(min_after_sunset = Time - Sunset) %>% ggplot(aes(x=min_after_sunset)) +
  geom_histogram(aes(fill=year),alpha=0.75,position = "identity") +
  xlab("Minutes after sunset") + ylab("# of collective display events") + theme_minimal(base_size = 20) +
  theme(axis.title = element_text(face = "bold")) +
  scale_fill_manual(values=c("#0072B2","#56B4E9")) +
  guides(fill=guide_legend(title="Year"))

fig3b <- d_survey %>% ggplot(aes(x=Days_after_fullmoon,y=dt,color=year)) +
  geom_jitter(size=2.5,aes(color=year)) +
  stat_smooth(fullrange = TRUE,method = "lm",alpha=0.75,se = FALSE) +
  xlab("# days after full moon") + ylab("Time between collective display events (s)") +
  theme_minimal(base_size = 20) +
  theme(axis.title = element_text(face = "bold")) +
  scale_color_manual(values=c("#0072B2","#56B4E9")) +
  guides(color=guide_legend(title="Year"))

fig3b <- ggMarginal(p = fig3b, data = gsum,x = Days_after_full_moon,y=n+1,type = "histogram",margins = "x",bins=14,size = 2.5,groupColour = TRUE,fill=NA,linewidth=2)
fig3 <- ggarrange(fig3a,fig3b,ncol = 2,labels = "AUTO",font.label = list(size=60,face="bold"))
ggsave(filename = "fig3.jpeg",plot = fig3,device = "jpeg",path = "~/Desktop/ESM_Hensley_etal_2023/figures/",dpi = 300,units = "px",height = 2500,width = 6000)

# Figure 2
fig4a <- xydat %>% filter(stim_duration == 0) %>% mutate(date_tank = interaction(date,tank),
                                                         Trial = case_when(date_tank == "5/19/17.A" ~ "Trial 4",
                                                                           date_tank == "5/19/17.C" ~ "Trial 3",
                                                                           date_tank == "5/18/17.B" ~ "Trial 1",
                                                                           date_tank == "5/18/17.C" ~ "Trial 2")) %>%
  arrange(date,tank,time) %>%
  ggplot(aes(x=time)) +
  geom_histogram(fill=NA,binwidth = 30) + facet_wrap(.~Trial,ncol = 4)

binwidth = layer_data(fig4a) %>% mutate(w=xmax-xmin) %>% pull(w) %>% median

fig4a <- fig4a + stat_bin(geom="step",linewidth=1.75,closed = "right",pad=TRUE,
                          aes(color=Trial),binwidth=30) + ##can change the BINWIDTH here is wanted #, position=position_nudge(x=-0.25*binwidth)) +
  xlab("Time (s)") + ylab("# of displays") +
  theme_minimal(base_size = 20) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        axis.text = element_text(size=30),
        axis.title = element_text(size=50,face="bold"),
        strip.text = element_text(size=50,face="bold"),
        panel.spacing.x = unit(1,"lines")) +
  #scale_x_continuous(limits=c(0,960),n.breaks = 3) +
  scale_color_manual(values = c(
    "Trial 1" = "#56B4E9",
    "Trial 2" = "#E69F00",
    "Trial 3" = "#009E73",
    "Trial 4" = "#CC79A7"
  ))

#Fig 4B Left Column

'''
> ex_situ_t
# A tibble: 4 × 2
# Rowwise:  date_tank
  date_tank               data
  <fct>     <list<tibble[,1]>>
1 5/19/17.A          [101 × 1]
2 5/18/17.B           [44 × 1]
3 5/18/17.C           [60 × 1]
4 5/19/17.C          [137 × 1]

"5/19/17.A" ~ "Trial 4",
"5/19/17.C" ~ "Trial 3",
"5/18/17.B" ~ "Trial 1",
"5/18/17.C" ~ "Trial 2"))

index 1 is Trial 4
index 2 is Trial 1
index 3 is Trial 2
index 4 is Trial 3
'''

synchrony_poisson_exp <- function(data_set1,min_obs,index){
  input <- as.character(data_set1$date_tank[[index]])
  data_date <- str_split_i(input,pattern = "\\.",i = 1)
  data_tank <- str_split_i(input,pattern = "\\.",i = 2)

  d <- data_set1$data[[index]]
  n_d <- nrow(d)

  tot_obs_time <- min_obs*60
  rate_expected <- n_d / tot_obs_time
  h_pois_expected <- rpois(tot_obs_time,rate_expected)
  h_data_expected <- data.frame(count = h_pois_expected,time = round(seq(from = 0,to = tot_obs_time, length.out = length(h_pois_expected))))

  h_data_expected$ind <- rep(1:100,each=30)[1:nrow(h_data_expected)] #each = window size

  #h_summary_exp <- h_data_expected %>% group_by(count) %>% dplyr::summarize(n = n()) %>% mutate(percent = (n/sum(n)*100))
  h_summary_exp <- h_data_expected %>% group_by(ind) %>% dplyr::summarize(n = sum(count)) %>% mutate(percent = (n/sum(n))) %>%
    mutate(group_bin = case_when(
      n < 1 ~ "0",
      n >= 1 & n <= 5 ~ "1 - 5",
      n > 5 & n <= 10 ~ "6 - 10",
      n > 10 ~ "11+"),
      group_bin = as.factor(group_bin)) %>% group_by(group_bin) %>%
    dplyr::summarise(sum_count = sum(n)) %>% mutate(percent = (sum_count/sum(sum_count)))

  return(h_summary_exp)
}

synchrony_poisson_obs <- function(data_set1,min_obs,index){
  input <- as.character(data_set1$date_tank[[index]])
  data_date <- str_split_i(input,pattern = "\\.",i = 1)
  data_tank <- str_split_i(input,pattern = "\\.",i = 2)

  d <- data_set1$data[[index]]
  n_d <- nrow(d)

  tot_obs_time <- min_obs*60

  nt <- data.frame(time = seq(0,tot_obs_time,1))
  dd <- d %>% mutate(disp_rate = 1/(lead(time,6) - time))
  dd_sim <- data.frame(time = round(dd$time,digits = 0),count = 1)
  reg <- merge(nt,dd_sim,all.x = TRUE)
  reg[is.na(reg)] <- 0
  h_data_obs <- reg %>% group_by(time) %>% dplyr::summarize(count = sum(count)) %>% ungroup()

  h_data_obs$ind <- rep(1:100,each=30)[1:nrow(h_data_obs)]

  #h_summary_obs <- h_data_obs %>% group_by(count) %>% dplyr::summarise( n = n()) %>% mutate(percent = (n/sum(n)*100))
  h_summary_obs <- h_data_obs %>% group_by(ind) %>% dplyr::summarize(n = sum(count)) %>% mutate(percent = (n/sum(n))) %>%
    mutate(group_bin = case_when(
      n < 1 ~ "0",
      n >= 1 & n <= 5 ~ "1 - 5",
      n > 5 & n <= 10 ~ "6 - 10",
      n > 10 ~ "11+"),
      group_bin = as.factor(group_bin)) %>% group_by(group_bin) %>%
    dplyr::summarise(sum_count = sum(n)) %>% mutate(percent = (sum_count/sum(sum_count)))

  return(h_summary_obs)
}

#h_size <- seq(0,11)
#h_lab <- as.character(seq(0,11))

### INDEX 1 ### TRIAL 4
ne1.1 <- synchrony_poisson_exp(ex_situ_t,15,1)
ne1.2 <- synchrony_poisson_obs(ex_situ_t,15,1)

extra_r1 <- data.frame(group_bin = factor("11+"),sum_count = 0,percent = 0)
extra_ne1.1 <- rbind(ne1.1,extra_r1)

#extra_ne1.2 <- data.frame(count = seq(max(ne1.2$count)+1,11))
#extra_r1 <- data.frame(n = rep(0,nrow(extra_ne1.2)), percent = rep(0,nrow(extra_ne1.2)))
#ne1.2 <- rbind(ne1.2,extra_rows1)

ne1 <- ne1.2 %>% mutate(group_bin = fct_relevel(group_bin,"0","1 - 5","6 - 10","11+")) %>%
  ggplot(aes(x=group_bin)) +
  geom_bar(aes(y=percent),stat = "identity",fill="#CC79A7") +
  geom_line(data=extra_ne1.1,aes(y=percent,group=1),color="#000000",linewidth=1.5) +
  geom_point(data=extra_ne1.1,aes(y=percent),color="#000000",size=3) +
  ylab("P( Observe # of signals in 30 s )") +
  xlab("# observed signals in 30 s") +
  theme_minimal(base_size = 20) +
  ylim(c(0,1)) +
  ggtitle(label = "Trial 4") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_text(size=50,face = "bold"),
        axis.title.y = element_blank(),
        axis.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5,size=50,face = "bold"))

# ne1 <- ne1.2 %>% ggplot(aes(x=count)) +
#   geom_bar(aes(y=percent),stat = "identity",fill="#56B4E9") +
#   geom_line(data=ne1.1,aes(y=percent),color="#000000",linewidth=1.5) +
#   geom_point(data=ne1.1,aes(y=percent),color="#000000",size=3) +
#   ylab("% of time with observed # signals") + xlab("# observed signals") +
#   theme_minimal(base_size = 20) +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major.x = element_blank(),
#         #axis.title = element_text(size=50,face = "bold"),
#         axis.title = element_blank(),
#         axis.text = element_text(size=30),
#         title = element_text(hjust = 0.5,size=50,face = "bold")) +
#   scale_x_continuous(breaks = h_size,labels = h_lab) + ylim(c(0,100)) +
#   ggtitle(label = "Trial 1") +
#   scale_y_break(ticklabels = c(5,10,15,90,95,100),breaks = c(15,85),scales = "free",expand = FALSE) +
#   theme(
#     plot.title = element_text(hjust = 0.5,face="bold"),
#     axis.text.y.right = element_blank(),
#     axis.line.y.right = element_blank(),
#     axis.ticks.y.right = element_blank()
#   )


### INDEX 2 ### TRIAL 1
ne2.1 <- synchrony_poisson_exp(ex_situ_t,15,2)
ne2.2 <- synchrony_poisson_obs(ex_situ_t,15,2)

extra_r2 <- data.frame(group_bin = c(factor("6 - 10"),factor("11+")),sum_count = c(0,0),percent = c(0,0))
extra_ne2.1 <- rbind(ne2.1,extra_r2)

extra_r2.2 <- data.frame(group_bin = factor("11+"),sum_count = 0,percent = 0)
extra_ne2.2 <- rbind(ne2.2,extra_r2.2)

# extra_ne2.2 <- data.frame(count = seq(max(ne2.2$count)+1,11))
# extra_r2 <- data.frame(n = rep(0,nrow(extra_ne2.2)), percent = rep(0,nrow(extra_ne2.2)))
# extra_rows2 <- cbind(extra_ne2.2,extra_r2)
# ne2.2 <- rbind(ne2.2,extra_rows2)

ne2 <- extra_ne2.2 %>% mutate(group_bin = fct_relevel(group_bin,"0","1 - 5","6 - 10","11+")) %>%
  ggplot(aes(x=group_bin)) +
  geom_bar(aes(y=percent),stat = "identity",fill="#56B4E9") +
  geom_line(data=extra_ne2.1,aes(y=percent,group=1),color="#000000",linewidth=1.5) +
  geom_point(data=extra_ne2.1,aes(y=percent),color="#000000",size=3) +
  ylab("P( Observe # of signals in 30 s )") + xlab("# observed signals in 30 s") +
  theme_minimal(base_size = 20) +
  ylim(c(0,1)) +
  ggtitle(label = "Trial 1") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        #axis.title = element_text(size=50,face = "bold"),
        axis.title = element_blank(),
        axis.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5,size=50,face = "bold"))


### INDEX 3 ### TRIAL 2
ne3.1 <- synchrony_poisson_exp(ex_situ_t,15,3)
ne3.2 <- synchrony_poisson_obs(ex_situ_t,15,3)

extra_r3 <- data.frame(group_bin = c(factor("11+"),factor("6 - 10")),sum_count = c(0,0),percent = c(0,0))
extra_ne3.1 <- rbind(ne3.1,extra_r3)

# extra_ne3.2 <- data.frame(count = seq(max(ne3.2$count)+1,11))
# extra_r3 <- data.frame(n = rep(0,nrow(extra_ne3.2)), percent = rep(0,nrow(extra_ne3.2)))
# extra_rows3 <- cbind(extra_ne3.2,extra_r3)
# ne3.2 <- rbind(ne3.2,extra_rows3)

ne3 <- ne3.2 %>% mutate(group_bin = fct_relevel(group_bin,"0","1 - 5","6 - 10","11+")) %>%
  ggplot(aes(x=group_bin)) +
  geom_bar(aes(y=percent),stat = "identity",fill="#E69F00") +
  geom_line(data=extra_ne3.1,aes(y=percent,group=1),color="#000000",linewidth=1.5) +
  geom_point(data=extra_ne3.1,aes(y=percent),color="#000000",size=3) +
  #ylab("P( Observe # of signals in 30 s )") +
  #xlab("# observed signals in 30 s") +
  theme_minimal(base_size = 20) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        #axis.title = element_text(size=50,face = "bold"),
        axis.title = element_blank(),
        axis.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5,size=50,face = "bold")) +
  ylim(c(0,1)) +
  ggtitle(label = "Trial 2")

#### INDEX 4 ### TRIAL 3
ne4.1 <- synchrony_poisson_exp(ex_situ_t,15,4)
ne4.2 <- synchrony_poisson_obs(ex_situ_t,15,4)

extra_r4 <- data.frame(group_bin = c(factor("11+"),factor("0")),sum_count = c(0,0),percent = c(0,0))
extra_ne4.1 <- rbind(ne4.1,extra_r4)
# extra_ne4.2 <- data.frame(count = seq(max(ne4.2$count)+1,11))
# extra_r4 <- data.frame(n = rep(0,nrow(extra_ne4.2)), percent = rep(0,nrow(extra_ne4.2)))
# extra_rows4 <- cbind(extra_ne4.2,extra_r4)
# ne4.2 <- rbind(ne4.2,extra_rows4)

ne4 <- ne4.2 %>% mutate(group_bin = fct_relevel(group_bin,"0","1 - 5","6 - 10","11+")) %>%
  ggplot(aes(x=group_bin)) +
  geom_bar(aes(y=percent),stat = "identity",fill="#009E73") +
  geom_line(data=extra_ne4.1,aes(y=percent,group=1),color="#000000",linewidth=1.5) +
  geom_point(data=extra_ne4.1,aes(y=percent),color="#000000",size=3) +
  #ylab("P( Observe # of signals in 30 s )") +
  #xlab("# signals observed in 30 s") +
  theme_minimal(base_size = 20) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        #axis.title = element_text(size=50,face = "bold"),
        axis.title = element_blank(),
        axis.text = element_text(size=30),
        plot.title = element_text(hjust = 0.5,size=50,face = "bold")) +
  ylim(c(0,1)) +
  ggtitle(label = "Trial 3")

## test for statistical deviation from Poisson distribution ##
# from Safarti et al. 2021
synchrony_prep <- function(data_set1,min_obs,index){
  input <- as.character(data_set1$date_tank[[index]])
  data_date <- str_split_i(input,pattern = "\\.",i = 1)
  data_tank <- str_split_i(input,pattern = "\\.",i = 2)

  d <- data_set1$data[[index]]
  n_d <- nrow(d)

  tot_obs_time <- min_obs*60
  nt <- data.frame(time = seq(0,tot_obs_time,1))
  dt <- data.frame(time = round(d$time,digits = 0),count = 1)
  reg <- merge(nt,dt,all.x = TRUE)
  reg[is.na(reg)] <- 0
  r <- reg %>% group_by(time) %>% dplyr::summarize(count = sum(count))

  #rate_nd <- n_d / tot_obs_time
  #h0 <- rpois(tot_obs_time,rate_nd) #in seconds

  return(reg)
}
nt1 <- synchrony_prep(ex_situ_t,16,1)
nt2 <- synchrony_prep(ex_situ_t,16,2)
nt3 <- synchrony_prep(ex_situ_t,16,3)
nt4 <- synchrony_prep(ex_situ_t,16,4)

z <- function(dat){
  colnames(dat) <- c("time","count")
  v <- var(dat$count)
  u <- mean(dat$count)
  n <- length(dat$count)

  z <- ((v - u) / sqrt(2*n*u*(n*u - 1)))*n*sqrt(n-1)
  return(z)
}
t <- function(dat){
  library(e1071)
  a <- skewness(dat$count)
  b <- kurtosis(dat$count) - 3
  n <- length(dat$count)

  t <- 0.5 * sqrt( (n) / (b^2 + 24*b + 6) ) * (a^2 - b)
  return(t)
}
p_erc <- function(x){
  #x = abs(x)/(2*sqrt(5)) #only use if the variance needs to be increased due to increased flashes (i.e. pulses per display being detected)
  x = abs(x)
  e = (-1 * ( (x^2) / (log(10) ) + log10(x*sqrt(pi) )))
  p = 10^(-1 * ( (x^2) / (log(10) ) + log10(x*sqrt(pi) )))
  #return(c(paste0("Exponent is ",e),paste0("p = ",p)))
  return(c(e,p))
}

z_score1 <- p_erc(z(nt1))
t_score1 <- p_erc(t(nt1))
n1grb1 <- paste0("z = ",round(z_score1[1],digits = 2),", p[z] = ",scales::scientific(z_score1[2], digits = 3))
n1grb2 <- paste0("t = ",round(t_score1[1],digits = 2),", p[t] = ",scales::scientific(t_score1[2], digits = 3))
ne1_fin <- ne1 + annotate("text",x=3.5, y=.95,label = n1grb1,size=20) + annotate("text",x=3.5, y=.75,label = n1grb2,size=20)

z_score2 <- p_erc(z(nt2))
t_score2 <- p_erc(t(nt2))
n2grb1 <- paste0("z = ",round(z_score2[1],digits = 2),", p[z] = ",scales::scientific(z_score2[2], digits = 3))
n2grb2 <- paste0("t = ",round(t_score2[1],digits = 2),", p[t] = ",scales::scientific(t_score2[2], digits = 3))
ne2_fin <- ne2 + annotate("text",x=3.5, y=.95,label = n2grb1,size=20) + annotate("text",x=3.5, y=.75,label = n2grb2,size=20)

z_score3 <- p_erc(z(nt3))
t_score3 <- p_erc(t(nt3))
n3grb1 <- paste0("z = ",round(z_score3[1],digits = 2),", p[z] = ",scales::scientific(z_score3[2], digits = 3))
n3grb2 <- paste0("t = ",round(t_score3[1],digits = 2),", p[t] = ",scales::scientific(t_score3[2], digits = 3))
ne3_fin <- ne3 + annotate("text",x=3.5, y=.95,label = n3grb1,size=20) + annotate("text",x=3.5, y=.75,label = n3grb2,size=20)

z_score4 <- p_erc(z(nt4))
t_score4 <- p_erc(t(nt4))
n4grb1 <- paste0("z = ",round(z_score4[1],digits = 2),", p[z] = ",scales::scientific(z_score4[2], digits = 3))
n4grb2 <- paste0("t = ",round(t_score4[1],digits = 2),", p[t] = ",scales::scientific(t_score4[2], digits = 3))
ne4_fin <- ne4 + annotate("text",x=3.5, y=.95,label = n4grb1,size=20) + annotate("text",x=3.5, y=.75,label = n4grb2,size=20)

fig4b <- ggarrange(print(ne2_fin),print(ne3_fin),print(ne4_fin),print(ne1_fin),nrow = 4,align = "h")
fig4b <- annotate_figure(fig4b,left = text_grob("P( Observe # of signals in 30 s )", color = "#000000", rot = 90,size=50,face="bold"))

### attempting to compare the distribution of nearest neighbor distances from the data compared to a distribution of nndist from uniformly
### signalling individuals

nn_nt_all2 <- nn_nt_all %>% nest_by(date_tank)

ks_test_sim_dat <- function(nested_dataset,index){
  dataset <- nested_dataset$data[[index]]

  dat_length <- length(dataset$time)
  samp_rate_median <- median(dataset$time)

  sim_dat <- rnormTrunc(dat_length,samp_rate_median,samp_rate_sd,min = 0)
  real_dat <- dataset$time

  p_value_ks <- ks.test(real_dat,sim_dat,exact = TRUE)
  return(p_value_ks)
}

ks_p1 <- ks_test_sim_dat(nn_nt_all2,1) #Trial 4
ks_p2 <- ks_test_sim_dat(nn_nt_all2,2) #Trial 1
ks_p3 <- ks_test_sim_dat(nn_nt_all2,3) #Trial 2
ks_p4 <- ks_test_sim_dat(nn_nt_all2,4) #Trial 3

ks_test_sim_dist <- function(nested_dataset,index){
  dataset <- nested_dataset$data[[index]]

  dat_length <- length(dataset$dist)
  samp_rate_median <- median(dataset$dist)

  sim_dat <- rnormTrunc(dat_length,samp_rate_median,samp_rate_sd,min = 0)
  real_dat <- dataset$dist

  p_value_ks <- ks.test(real_dat,sim_dat,exact = TRUE)
  return(p_value_ks)
}

ks_pdist1 <- ks_test_sim_dist(nn_nt_all2,1) #Trial 4
ks_pdist2 <- ks_test_sim_dist(nn_nt_all2,2) #Trial 1
ks_pdist3 <- ks_test_sim_dist(nn_nt_all2,3) #Trial 2
ks_pdist4 <- ks_test_sim_dist(nn_nt_all2,4) #Trial 3

fig4c <- nn_nt_all %>% mutate(Trial = case_when(date_tank == "5/19/17.A" ~ "Trial 4",
                                                date_tank == "5/19/17.C" ~ "Trial 3",
                                                date_tank == "5/18/17.B" ~ "Trial 1",
                                                date_tank == "5/18/17.C" ~ "Trial 2")) %>%
  ggplot(aes(y=dist,x=time)) +
  geom_point(aes(fill=density,color=Trial),size=5.5,shape=21,stroke=1.25) +
  scale_fill_gradient(name="Density",high="black",low="white",guide="colorbar") +
  theme_minimal(base_size = 20) +
  theme(legend.position = "left") +
  xlab("Intersignal timing (s)") + ylab("Intersignal distance (cm)") +
  theme(axis.title = element_text(size=50,face = "bold"),
        axis.text = element_text(size=30),
        legend.title = element_text(face="bold"),
        strip.text = element_text(size=50,face="bold"),
        legend.position = "none") +
  xlim(0,32) + ylim(0,40) +
  scale_x_continuous(expand = c(0.01,0.01)) +
  facet_grid(Trial~.) +
  scale_color_manual(values = c(
    "Trial 1" = "#56B4E9",
    "Trial 2" = "#E69F00",
    "Trial 3" = "#009E73",
    "Trial 4" = "#CC79A7"
  ))

fig4c_ks_text_time <- data.frame(p_values = c(paste0("p[D_time] = ",scales::scientific(ks_p1$p.value,digits = 3)),
                                         paste0("p[D_time] = ",scales::scientific(ks_p2$p.value,digits = 3)),
                                         paste0("p[D_time] = ",scales::scientific(ks_p3$p.value,digits = 3)),
                                         paste0("p[D_time] = ",scales::scientific(ks_p4$p.value,digits = 3))),
                            Trial = c("Trial 4","Trial 1","Trial 2","Trial 3"))

fig4c_ks_text_dist <- data.frame(p_values = c(paste0("p[D_dist] = ",scales::scientific(ks_pdist1$p.value,digits = 3)),
                                         paste0("p[D_dist] = ",scales::scientific(ks_pdist2$p.value,digits = 3)),
                                         paste0("p[D_dist] = ",scales::scientific(ks_pdist3$p.value,digits = 3)),
                                         paste0("p[D_dist] = ",scales::scientific(ks_pdist4$p.value,digits = 3))),
                            Trial = c("Trial 4","Trial 1","Trial 2","Trial 3"))

fig4c_fin <- fig4c + geom_text(data=fig4c_ks_text_dist,aes(label = p_values),x = 25, y = 35,size=20) +
  geom_text(data=fig4c_ks_text_time,aes(label = p_values),x = 25, y = 30,size=20)


fig4.1 <- ggarrange(fig4b,fig4c_fin,ncol=2,widths = c(0.9,1))
fig4 <- ggarrange(fig4a,fig4.1,nrow=2,heights = c(1,2))
ggsave(filename = "fig4.2.jpeg",plot = fig4,device = "jpeg",path = "~/Desktop/ESM_Hensley_etal_2023/figures/",dpi = 300,units = "px",height = 12000,width = 12000,limitsize = FALSE)

# Figure 3
disp_db <- disp_dat %>% dplyr::mutate(duration = time, ipi = ipi_time_s) %>%
  dplyr::select(display_num,pulse_num,duration,ipi,ipd_distance_mm,vert_dist_mm,relative_brightness) %>%
  gather(key = "trait",value = "measurement",duration,ipi,ipd_distance_mm,vert_dist_mm,relative_brightness) %>%
  filter(!is.na(measurement)) %>%
  group_by(trait,pulse_num) %>%
  dplyr::summarise(mean = mean(measurement,na.rm=TRUE),
                   se = sd(measurement,na.rm=TRUE)/sqrt(n())) %>%
  dplyr::filter(pulse_num == 1)

disp_db1 <- as.data.frame(disp_db[1,]) %>% rename(Stimulus_duration = mean) %>% mutate(Response = 0)

pm = lm(Response ~ poly(Stimulus_duration,2) + Tank/Tank_order,data=dur)
summary(pm)
anova(pm) #nested interaction is significant; trial, trial_order are not significant
report(pm)
check_model(pm)
effect_plot(pm,Stimulus_duration)
tab_model(pm,dv.labels = c("# of displays"))

pm2 <- lm(Response ~ as.factor(Stimulus_duration) + Tank/Tank_order,data=dur)
pmmeans <- emmeans(pm2,pairwise ~ Stimulus_duration)
model_means_cld <- cld(object = pmmeans,
                       adjust = "sidak",
                       Letters = letters,
                       alpha = 0.05)

fig5a <- ggplot(data = dur, aes(x=Stimulus_duration,y=Response)) +
  geom_jitter(color="#56B4E9",alpha=0.75,size=8) +
  geom_boxplot(aes(group=Stimulus_duration),outlier.shape = NA,fill=NA,color="#0072B2",size=2) +
  stat_smooth(method = "lm", formula = y ~ poly(x,2),se=FALSE,fullrange = TRUE) +
  theme_minimal(base_size = 20) +
  xlab('Duration of stimulus (s)') +
  ylab("# of displays") +
  geom_text(data=model_means_cld,aes(label=.group,y=max(emmean + 10*SE)),size=10,color="#0072B2") +
  theme(axis.title = element_text(size=50,face = "bold"),
        axis.text = element_text(size=30)) +
  EnvStats::stat_n_text(color = "#0072B2",size=10) +
  geom_pointrange(data = disp_db1,stat="identity",
                  aes(xmin=Stimulus_duration-1.96*se,xmax=Stimulus_duration+1.96*se),fatten = 10,
                  color="black")

chi_dat <- old_dat3 %>% filter(signal_timing == "signal_start_delay")
chisq.test(table(chi_dat$signal_onset))

delay_dat_during <- old_dat3 %>% filter(signal_onset == "During",signal_timing == "signal_start_delay")
m3 <- lm(delay~stim_duration + order + eu_dist,data=delay_dat_during) #distance, order is statistically significantly more explanatory
summary(m3); anova(m3_alt,m3)
check_model(m3)
effect_plot(m3,stim_duration)
report(m3)
tab_model(m3,dv.labels = "Latency to respond (s)",show.df = TRUE,pred.labels = c("(Intercept)","Stimulus duration (s)","Testing order","Distance from stimulus (cm)"))

disp_db2 <- as.data.frame(disp_db[1,]) %>% rename(stim_duration = mean) %>% mutate(delay=0,xmin = -Inf,xmax = Inf)

fig5c <- delay_dat_during %>% rename(Distance = eu_dist) %>%
  ggplot(aes(x=stim_duration,y=delay)) +
  geom_rect(data=disp_db2,
            aes(ymin=stim_duration-1.96*se,
                ymax=stim_duration+1.96*se,
                xmin=xmin, xmax=xmax),
            fill="darkgrey",alpha=0.2) +
  geom_abline(slope=0, intercept=5.770833,color="black",linetype=2,linewidth = 1.25) +
  geom_jitter(aes(color=Distance),size=8,shape=21,stroke=2) +
  scale_color_gradient(name = "Distance from stimulus (cm)",low="#F0E442",high="#009E73",
                       guide="colourbar",aesthetics = ("colour")) +
  geom_boxplot(aes(group=as.factor(stim_duration)),color="#0072B2",outlier.shape = NA,fill=NA,size=2) +
  stat_smooth(method = "lm",formula = y ~ x,se = FALSE,fullrange = TRUE) +
  theme_minimal(base_size = 20) +
  theme(axis.title = element_text(size=50,face = "bold"),
        axis.text = element_text(size=30),
        legend.title = element_text(face="bold",size=30),
        legend.text = element_text(size=30),
        legend.position = "top") +
  ylim(c(0,12)) +
  xlab("Duration of stimulus (s)") +
  ylab("Latency to display\nduring stimulus (s)") +
  guides(colour = guide_colorbar(label.theme=element_text(size=15)))

fig5c <- ggMarginal(fig5c,margins = "y",type = "histogram",color="#0072B2",fill=NA,size = 10,linewidth=2)

delay_dat_after <- old_dat3 %>% filter(signal_onset == "After",signal_timing == "signal_end_delay",delay <= 2.5)
m4 <- lm(delay~poly(stim_duration,2) + order,data=delay_dat_after) #distance is not significant
summary(m4); anova(m4)
check_model(m4)
effect_plot(m4,pred = stim_duration)
report(m4)
tab_model(m4,dv.labels = "Latency to respond (s)",show.df = TRUE,pred.labels = c("(Intercept)","Stimulus duration [1st degree]","Stimulus duration [2nd degree]","Testing order"))

delay_m2 <- lm(delay ~ as.factor(stim_duration) + order,data=delay_dat_after)
dmmeans <- emmeans(delay_m2,pairwise ~ stim_duration)
dm_means_cld <- cld(object = dmmeans,
                       adjust = "sidak",
                       Letters = letters,
                       alpha = 0.05)

fig5d <- delay_dat_after %>%
  ggplot(aes(y=delay,x=stim_duration)) +
  geom_jitter(color="#56B4E9",alpha=0.75,size=8) +
  geom_boxplot(aes(group=stim_duration),outlier.shape = NA,fill=NA,color="#0072B2",size=2) +
  geom_text(data=dm_means_cld,aes(label=.group,y=max(emmean + 7*SE)),size=10,color="#0072B2") +
  stat_smooth(method = "lm", formula = y ~ poly(x,2),se=FALSE,fullrange = TRUE) +
  theme_minimal(base_size = 20) +
  theme(axis.title = element_text(size=50,face = "bold"),
        axis.text = element_text(size=30)) +
  EnvStats::stat_n_text(color = "#0072B2",size=8) +
  xlab("Duration of stimulus (s)") + ylab("Latency to display\nafter stimulus (s)") +
  geom_pointrange(data = disp_db2,stat="identity",
                  aes(xmin=stim_duration-1.96*se,xmax=stim_duration+1.96*se),fatten = 10,
                  color="black")

fig5 <- ggarrange(fig5a,fig5d,fig5c,nrow=1,ncol = 3,labels = "AUTO",font.label = list(size=60,face="bold"))
ggsave(filename = "fig5.2.jpeg",plot = fig5,device = "jpeg",path = "~/Desktop/ESM_Hensley_etal_2023/figures/",dpi = 300,units = "px",height = 4000,width = 14000,limitsize = FALSE)

# Figure 4
## WARNING: THE PERMUTED DATASETS IN THE NEXT SECTION TAKE A VERY LONG TIME TO COMPUTE. DO NOT RE-PERMUTE THE DATASETS ONCE YOU HAVE MADE THEM ALREADY

time_neighbors <- function(data_set2,index){
  d_t <- data_set2$data[[index]]
  tDt <- dist(d_t,method = "maximum",diag = FALSE,upper = FALSE)
  mt <- as.matrix(tDt)
  diag(mt) <- Inf
  smallest_distances <- apply(mt, 2, min)

  return(smallest_distances)
}
nnt1 <- time_neighbors(ex_situ_t,1)
nnt2 <- time_neighbors(ex_situ_t,2)
nnt3 <- time_neighbors(ex_situ_t,3)
nnt4 <- time_neighbors(ex_situ_t,4)

temp_nb <- c(nnt1,nnt2,nnt3,nnt4)

nn_np_all <- cbind(nn_np_all,temp_nb)

nn_np_all <- nn_np_all %>%
  dplyr::mutate(date = c(rep("5/19/17",101),rep("5/18/17",44),rep("5/18/17",60),rep("5/19/17",137)),
                tank = c(rep("A",101),rep("B",44),rep("C",60),rep("C",137)),
                date_tank = interaction(date,tank),
                Trial = case_when(date_tank == "5/19/17.A" ~ "Trial 4",
                                  date_tank == "5/19/17.C" ~ "Trial 3",
                                  date_tank == "5/18/17.B" ~ "Trial 1",
                                  date_tank == "5/18/17.C" ~ "Trial 2")) %>%
  filter(num_pulse > 0)

fig6a <- nn_np_all %>%
  ggplot(aes(y=nearest_neighbor,x=num_pulse)) +
  geom_jitter(size=8,aes(color=temp_nb)) +
  geom_boxplot(color="#0072B2",aes(group=as.factor(num_pulse)),outlier.shape = NA,fill=NA,size=2) +
  geom_smooth(method="lm",se=FALSE) +
  xlab("Number of pulses per display") +
  ylab("Distance to nearest neighbor (cm)") +
  theme_minimal(base_size = 20,) +
  scale_color_gradient(name = "Time to next soonest display (s)",low="lightgrey",high="black",
                       guide="colourbar",aesthetics = ("colour")) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9)) +
  theme(axis.title = element_text(size=50,face = "bold"),
        axis.text = element_text(size=30),
        legend.title = element_text(size=30,face="bold"),
        legend.position = "top") +
  EnvStats::stat_n_text(color = "#0072B2",size=8) +
  guides(colour = guide_colorbar(label.theme=element_text(size=10)))

##trying a permutation for each trial's nearest neighbor distances
nearest_matrix_shuffle_dist <- function(data_set,index,n){
  input <- as.character(data_set$date_tank[[index]])

  d <- data_set$data[[index]]
  tD <- dist(d,method = "euclidean",diag = FALSE,upper = FALSE)
  m <- as.matrix(tD)
  diag(m) <- Inf
  size_m <- nrow(m)

  nn1 <- list()
  for (i in 1:n) {
    nn_perm <- numeric(size_m)
    for (j in 1:n){
      z <- sample(1:size_m,size_m,replace=FALSE)
      tmp <- m[z,]
      tmp <- tmp[,z]

      smallest_distances <- apply(tmp, 2, min)
    }
    nn1[[i]] <- smallest_distances
  }
  df <- do.call("rbind",nn1)
  df <- t(df)
  #colnames(df) <- paste0(input,"_perm_",seq(1:n))
  colnames(df) <- paste0("perm_",seq(1:n))
  return(as.data.frame(df))
} #n is num. of shuffles #only use this one with distance data
nearest_matrix_shuffle_time <- function(data_set,index,n){
  input <- as.character(data_set$date_tank[[index]])

  d <- data_set$data[[index]]
  tD <- dist(d,method = "maximum",diag = FALSE,upper = FALSE)
  m <- as.matrix(tD)
  diag(m) <- Inf
  size_m <- nrow(m)

  nn1 <- list()
  for (i in 1:n) {
    nn_perm <- numeric(size_m)
    for (j in 1:n){
      z <- sample(1:size_m,size_m,replace=FALSE)
      tmp <- m[z,]
      tmp <- tmp[,z]

      smallest_distances <- apply(tmp, 2, min)
    }
    nn1[[i]] <- smallest_distances
  }
  df <- do.call("rbind",nn1)
  df <- t(df)
  #colnames(df) <- paste0(input,"_perm_",seq(1:n))
  colnames(df) <- paste0("perm_",seq(1:n))
  return(as.data.frame(df))
} #only use this one to look at nearest time neighbors

#doing a distance permutation on observed trials to test for nearest neighbor effect
nn_trial1 <- nearest_matrix_shuffle_dist(ex_situ,1,1000)
nn_trial2 <- nearest_matrix_shuffle_dist(ex_situ,2,1000)
nn_trial3 <- nearest_matrix_shuffle_dist(ex_situ,3,1000)
nn_trial4 <- nearest_matrix_shuffle_dist(ex_situ,4,1000)

nn1_alldat_perm <- cbind(nn1,time = nnt1,nn_trial1)
nn2_alldat_perm <- cbind(nn2,time = nnt2,nn_trial2)
nn3_alldat_perm <- cbind(nn3,time = nnt3,nn_trial3)
nn4_alldat_perm <- cbind(nn4,time = nnt4,nn_trial4)

nn1_tidy <- nn1_alldat_perm %>% pivot_longer(cols = c(starts_with("perm_"),"nearest_neighbor"),names_to = "dataset",values_to = "distance")
nn2_tidy <- nn2_alldat_perm %>% pivot_longer(cols = c(starts_with("perm_"),"nearest_neighbor"),names_to = "dataset",values_to = "distance")
nn3_tidy <- nn3_alldat_perm %>% pivot_longer(cols = c(starts_with("perm_"),"nearest_neighbor"),names_to = "dataset",values_to = "distance")
nn4_tidy <- nn4_alldat_perm %>% pivot_longer(cols = c(starts_with("perm_"),"nearest_neighbor"),names_to = "dataset",values_to = "distance")

nn_all_tidy <- rbind(nn1_tidy,nn2_tidy,nn3_tidy,nn4_tidy)

nn_map_graph <- nn_all_tidy %>% tidyr::nest(.by = dataset) %>%
  mutate(model = map(data, ~lm(num_pulse ~ distance + time, data = .)),
         cleaned_fit = map(model,tidy),
         results = map(model,glance)) %>%
  dplyr::select(dataset,cleaned_fit) %>% unnest(cleaned_fit) %>%
  dplyr::filter(term == "distance") %>%
  mutate(color = ifelse(dataset == "nearest_neighbor","blue","black"))

nn_observed <- subset.data.frame(nn_map_graph,dataset == 'nearest_neighbor')

p.table <- table(abs(nn_map_graph$statistic) >= abs(nn_observed$statistic))
p.test <- (p.table[2] - 1)/1000

fig6b <- nn_map_graph %>% ggplot(aes(x=statistic)) +
  geom_histogram(fill=NA,colour="black") +
  geom_vline(xintercept = nn_observed$statistic,colour="blue",linewidth=1.5,linetype=2) +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab(expression(paste(italic(t)," statistic"))) + ylab("Count") +
  ggtitle("Permutation test:\nDistance to nearest neighbour")

#doing a time permutation on observed trials to test for nearest time effect
nn_trial1_t <- nearest_matrix_shuffle_time(ex_situ_t,1,1000)
nn_trial2_t <- nearest_matrix_shuffle_time(ex_situ_t,2,1000)
nn_trial3_t <- nearest_matrix_shuffle_time(ex_situ_t,3,1000)
nn_trial4_t <- nearest_matrix_shuffle_time(ex_situ_t,4,1000)

nn1_alldat_perm_t <- cbind(num_pulse = nn1$num_pulse,distance = nn1$nearest_neighbor,temp_near = nnt1,nn_trial1_t)
nn2_alldat_perm_t <- cbind(num_pulse = nn2$num_pulse,distance = nn2$nearest_neighbor,temp_near = nnt2,nn_trial2_t)
nn3_alldat_perm_t <- cbind(num_pulse = nn3$num_pulse,distance = nn3$nearest_neighbor,temp_near = nnt3,nn_trial3_t)
nn4_alldat_perm_t <- cbind(num_pulse = nn4$num_pulse,distance = nn4$nearest_neighbor,temp_near = nnt4,nn_trial4_t)

nn1_tidy_t <- nn1_alldat_perm_t %>% pivot_longer(cols = c(starts_with("perm_"),"temp_near"),names_to = "dataset",values_to = "time")
nn2_tidy_t <- nn2_alldat_perm_t %>% pivot_longer(cols = c(starts_with("perm_"),"temp_near"),names_to = "dataset",values_to = "time")
nn3_tidy_t <- nn3_alldat_perm_t %>% pivot_longer(cols = c(starts_with("perm_"),"temp_near"),names_to = "dataset",values_to = "time")
nn4_tidy_t <- nn4_alldat_perm_t %>% pivot_longer(cols = c(starts_with("perm_"),"temp_near"),names_to = "dataset",values_to = "time")

nn_all_tidy_t <- rbind(nn1_tidy_t,nn2_tidy_t,nn3_tidy_t,nn4_tidy_t)

nn_map_all_t <- nn_all_tidy_t %>% tidyr::nest(.by = dataset) %>%
  mutate(model = map(data, ~lm(num_pulse ~ distance + time, data = .)),
         cleaned_fit = map(model,tidy),
         results = map(model,broom::glance))

nn_map_graph_t <- nn_map_all_t %>%
  dplyr::select(dataset,cleaned_fit) %>% unnest(cleaned_fit) %>%
  dplyr::filter(term == "time") %>%
  mutate(color = ifelse(dataset == "temp_near","blue","black"))

nn_observed_t <- subset.data.frame(nn_map_graph_t,dataset == 'temp_near')

p.table_t <- table(abs(nn_map_graph_t$statistic) >= abs(nn_observed_t$statistic))
p.test_t <- (p.table_t[2] - 1)/1000

fig6c <- nn_map_graph_t %>% ggplot(aes(x=statistic)) +
  geom_histogram(fill=NA,colour="black") +
  geom_vline(xintercept = nn_observed_t$statistic,colour="blue",linewidth=1.5,linetype=2) +
  xlab(expression(paste(italic(t)," statistic"))) + ylab("Count") +
  ggtitle("Permutation test:\nTime to next soonest display") +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5))

fig6.1 <- ggarrange(fig6b,fig6c,ncol = 1,labels = c("B","C"),font.label = list(size=50,face="bold"))
fig6.2 <- ggarrange(fig6a,fig6.1,ncol = 2,labels = c("A",NA),font.label = list(size=50,face="bold"),widths = c(2,1))

### ATTEMPT TO USE NEAREST NEIGHBOR AND DISTANCE TO STIMULUS ###
##attempting to do a nested dist matrix analysis on the experimental data (duration != 0) only
ex_experimental <- xydat %>% filter(stim_duration != 0,new_time <= stim_duration + 90) %>% dplyr::mutate(date_tank = interaction(date,tank,stim_duration,order)) %>% ungroup() %>%
  dplyr::select(c(x_pos,y_pos,date_tank)) %>% group_by(date_tank) %>%
  filter(n() > 2) %>% ungroup() %>%
  nest_by(date_tank)

ex_experimental_t <- xydat %>% filter(stim_duration != 0,new_time <= stim_duration + 90) %>% dplyr::mutate(date_tank = interaction(date,tank,stim_duration,order)) %>% ungroup() %>%
  dplyr::select(c(time,date_tank)) %>% group_by(date_tank) %>%
  filter(n() > 2) %>% ungroup() %>%
  nest_by(date_tank)

ex_experimental_dstim <- xydat %>% filter(stim_duration != 0,new_time <= stim_duration + 90) %>% dplyr::mutate(date_tank = interaction(date,tank,stim_duration,order)) %>% ungroup() %>%
  dplyr::select(c(eu_dist,date_tank)) %>% group_by(date_tank) %>%
  filter(n() > 2) %>% ungroup() %>%
  nest_by(date_tank)

nearest_neighbor_stim <- function(data_set,index){
  input <- as.character(data_set$date_tank[[index]])
  data_date <- str_split_i(input,pattern = "\\.",i = 1)
  data_tank <- str_split_i(input,pattern = "\\.",i = 2)

  d <- data_set$data[[index]]
  nn <- nndist(d,k = 1)
  np <- xydat %>% filter(stim_duration != 0,new_time <= stim_duration + 90) %>% dplyr::mutate(date_tank = interaction(date,tank,stim_duration,order)) %>% ungroup() %>%
    filter(date_tank == input) %>%
    dplyr::select(num_pulse,eu_dist,date_tank,new_time,stim_duration,order)
  nrdp <- data.frame(np,nearest_neighbor = nn)

  return(nrdp)
}

all_np_dat_nn <- data.frame()
for (i in 1:nrow(ex_experimental)){
  nnprep <- nearest_neighbor_stim(ex_experimental,i)
  all_np_dat_nn <- rbind(all_np_dat_nn,nnprep)
}

all_np_dat_nt <- c()
for (i in 1:nrow(ex_experimental_t)){
  nnprep <- time_neighbors(ex_experimental_t,i)
  all_np_dat_nt <- c(all_np_dat_nt,nnprep)
}

all_np_dat_nn <- cbind(all_np_dat_nn,temp_nb = all_np_dat_nt) %>% filter(new_time <= stim_duration + 90)

fig6d <- all_np_dat_nn %>%
  ggplot(aes(x=num_pulse,y=nearest_neighbor)) +
  geom_jitter(aes(color=temp_nb),size=8)+#shape=21,stroke=2) +
  geom_boxplot(color="#0072B2",aes(group=as.factor(num_pulse)),outlier.shape = NA,size=2,fill=NA) +
  geom_smooth(method="lm",se=FALSE) +
  xlab("Number of pulses per display") +
  ylab("Distance to nearest neighbor (cm)") + theme_minimal(base_size = 20) +
  scale_color_gradient(name = "Time to next soonest display (s)",low="lightgrey",high="black",
                        guide="colourbar",aesthetics = ("colour")) +
  #scale_color_gradient(name = "Distance from stimulus (cm)",low="#F0E442",high="#009E73",
  #                     guide="colourbar",aesthetics = ("colour")) +
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8,9)) +
  theme(axis.title = element_text(size=50, face = "bold"),
        axis.text = element_text(size=30),
        legend.title = element_text(size=30,face="bold"),legend.position = "top") +
  EnvStats::stat_n_text(color = "#0072B2",size=8) +
  guides(colour = guide_colorbar(label.theme=element_text(size=10)))

#permutation tests for nearest neighbor
n_shuff <- 1000
nnexp_perm_dist <- data.frame()
for (i in 1:nrow(ex_experimental)){
  dist_shuff <- nearest_matrix_shuffle_dist(ex_experimental,i,n_shuff)
  nnexp_perm_dist <- rbind(nnexp_perm_dist,dist_shuff)
}

nnexp_perm_dist_all <- subset.data.frame(all_np_dat_nn,select = c('num_pulse','eu_dist','nearest_neighbor','temp_nb'))
nnexp_perm_dist_all <- cbind(nnexp_perm_dist_all,nnexp_perm_dist)

nnexp_tidy <- nnexp_perm_dist_all %>%
  pivot_longer(cols = c(starts_with("perm_"),"nearest_neighbor"),names_to = "dataset",values_to = "distance")

nnexp_map <- nnexp_tidy %>% tidyr::nest(.by = dataset) %>%
  mutate(model = map(data, ~lm(num_pulse ~ distance + temp_nb + eu_dist, data = .)),
         cleaned_fit = map(model,tidy),
         results = map(model,broom::glance))

nnexp_map_graph <- nnexp_map %>%
  dplyr::select(dataset,cleaned_fit) %>% unnest(cleaned_fit) %>%
  dplyr::filter(term == "distance") %>%
  mutate(color = ifelse(dataset == "nearest_neighbor","blue","black"))

nnexp_observed <- subset.data.frame(nnexp_map_graph,dataset == 'nearest_neighbor')

p.table_exp <- table(abs(nnexp_map_graph$statistic) >= abs(nnexp_observed$statistic))
p.test_exp <- (p.table_exp[2] - 1)/1000

fig6e <- nnexp_map_graph %>% ggplot(aes(x=statistic)) +
  geom_histogram(fill=NA,colour="black") +
  geom_vline(xintercept = nnexp_observed$statistic,colour="blue",linewidth=1.5,linetype=2) +
  xlab(expression(paste(italic(t)," statistic"))) + ylab("Count") +
  ggtitle("Permutation test:\nDistance to nearest neighbour") +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5))

#permutation tests for next soonest display
nnexp_perm_t <- data.frame()
for (i in 1:nrow(ex_experimental_t)){
  time_shuff <- nearest_matrix_shuffle_time(ex_experimental_t,i,n_shuff)
  nnexp_perm_t <- rbind(nnexp_perm_t,time_shuff)
}

nnexp_perm_time_all <- subset.data.frame(all_np_dat_nn,select = c('num_pulse','eu_dist','nearest_neighbor','temp_nb'))
nnexp_perm_time_all <- cbind(nnexp_perm_time_all,nnexp_perm_t)

nnexp_tidy_t <- nnexp_perm_time_all %>%
  pivot_longer(cols = c(starts_with("perm_"),"temp_nb"),names_to = "dataset",values_to = "time")

nnexp_map_t <- nnexp_tidy_t %>% tidyr::nest(.by = dataset) %>%
  mutate(model = map(data, ~lm(num_pulse ~ nearest_neighbor + time + eu_dist, data = .)),
         cleaned_fit = map(model,tidy),
         results = map(model,broom::glance))

nnexp_map_graph_t <- nnexp_map_t %>%
  dplyr::select(dataset,cleaned_fit) %>% unnest(cleaned_fit) %>%
  dplyr::filter(term == "time") %>%
  mutate(color = ifelse(dataset == "temp_near","blue","black"))

nnexp_observed_t <- subset.data.frame(nnexp_map_graph_t,dataset == 'temp_nb')

p.table_exp_t <- table(abs(nnexp_map_graph_t$statistic) >= abs(nnexp_observed_t$statistic))
p.test_exp_t <- (p.table_exp_t[2] - 1)/1000

fig6f <- nnexp_map_graph_t %>% ggplot(aes(x=statistic)) +
  geom_histogram(fill=NA,colour="black") +
  geom_vline(xintercept = nnexp_observed_t$statistic,colour="blue",linewidth=1.5,linetype=2) +
  xlab(expression(paste(italic(t)," statistic"))) + ylab("Count") +
  ggtitle("Permutation test:\nTime to next soonest display") +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5))

#permutation tests for distance to stimulus
nnexp_perm_dstim <- data.frame()
for (i in 1:nrow(ex_experimental_dstim)){
  dstim_shuff <- nearest_matrix_shuffle_time(ex_experimental_dstim,i,n_shuff)
  nnexp_perm_dstim <- rbind(nnexp_perm_dstim,dstim_shuff)
}

nnexp_perm_dstim_all <- subset.data.frame(all_np_dat_nn,select = c('num_pulse','eu_dist','nearest_neighbor','temp_nb'))
nnexp_perm_dstim_all <- cbind(nnexp_perm_dstim_all,nnexp_perm_dstim)

nnexp_tidy_dstim <- nnexp_perm_dstim_all %>%
  pivot_longer(cols = c(starts_with("perm_"),"eu_dist"),names_to = "dataset",values_to = "dstim")

nnexp_map_dstim <- nnexp_tidy_dstim %>% tidyr::nest(.by = dataset) %>%
  mutate(model = map(data, ~lm(num_pulse ~ nearest_neighbor + temp_nb + dstim, data = .)),
         cleaned_fit = map(model,tidy),
         results = map(model,broom::glance))

nnexp_map_graph_dstim <- nnexp_map_dstim %>%
  dplyr::select(dataset,cleaned_fit) %>% unnest(cleaned_fit) %>%
  dplyr::filter(term == "dstim") %>%
  mutate(color = ifelse(dataset == "eu_dist","blue","black"))

nnexp_observed_dstim <- subset.data.frame(nnexp_map_graph_dstim,dataset == 'eu_dist')

p.table_exp_dstim <- table(abs(nnexp_map_graph_dstim$statistic) >= abs(nnexp_observed_dstim$statistic))
p.test_exp_dstim <- (p.table_exp_dstim[2] - 1)/1000

fig6g <- nnexp_map_graph_dstim %>% ggplot(aes(x=statistic)) +
  geom_histogram(fill=NA,colour="black") +
  geom_vline(xintercept = nnexp_observed_dstim$statistic,colour="blue",linewidth=1.5,linetype=2) +
  xlab(expression(paste(italic(t)," statistic"))) + ylab("Count") +
  ggtitle("Permutation test:\nDistance to stimulus") +
  theme_minimal(base_size = 30) +
  theme(plot.title = element_text(hjust = 0.5))

fig6.3 <- ggarrange(fig6e,fig6f,fig6g,ncol = 1,labels = c("E","F","G"),font.label = list(size=50,face="bold"))
fig6.4 <- ggarrange(fig6d,fig6.3,ncol = 2,labels = c("D",NA),font.label = list(size=50,face="bold"),widths = c(2,1))

fig6 <- ggarrange(fig6.2,fig6.4,ncol = 1)
ggsave(filename = "fig6.jpeg",plot = fig6,device = "jpeg",path = "~/Desktop/ESM_Hensley_etal_2023/figures/",dpi = 300,units = "px",height = 8000,width = 8000,limitsize = FALSE)


fig6new <- ggarrange(fig6a,fig6d,nrow = 1,labels = "AUTO",common.legend = TRUE,font.label = list(size=50,face="bold"))
ggsave(filename = "fig6new.jpeg",plot = fig6new,device = "jpeg",path = "~/Desktop/ESM_Hensley_etal_2023/figures/",dpi = 300,units = "px",height = 4800,width = 11000,limitsize = FALSE)

fig6supp_top <-ggarrange(fig6b,fig6c,nrow = 1)
fig6supp_btm <-ggarrange(fig6e,fig6f,fig6g,nrow = 1)
fig6supp <-ggarrange(fig6supp_top,fig6supp_btm,nrow = 2,labels = "AUTO",font.label = list(size=50,face="bold"))
ggsave(filename = "fig6supp.jpeg",plot = fig6supp,device = "jpeg",path = "~/Desktop/ESM_Hensley_etal_2023/figures/",dpi = 300,units = "px",height = 5000,width = 8000,limitsize = FALSE)

## Supplementary Results ##

# Figure S1
library(WaveletComp)
library(TSA)

t <- data.frame(time = seq(0,960))
d_fft <- layer_data(fig4a) %>% dplyr::select(count,x,PANEL) %>% rename(time = x)

dt1 <- d_fft %>% filter(PANEL==1) %>% dplyr::select(count,time) %>% mutate(time=round(time,digits = 0))
dt1.1 <- left_join(t,dt1,"time") %>% mutate(count = case_when(is.na(count) ~ 0,
                                                              !is.na(count) ~ count),
                                            count_diff = count - lag(count,1)) %>% filter(!is.na(count_diff))

dt2 <- d_fft %>% filter(PANEL==2) %>% dplyr::select(count,time) %>% mutate(time=round(time,digits = 0))
dt2.1 <- left_join(t,dt2,"time") %>% mutate(count = case_when(is.na(count) ~ 0,
                                                              !is.na(count) ~ count),
                                            count_diff = count - lag(count,1)) %>% filter(!is.na(count_diff))

dt3 <- d_fft %>% filter(PANEL==3) %>% dplyr::select(count,time) %>% mutate(time=round(time,digits = 0))
dt3.1 <- left_join(t,dt3,"time") %>% mutate(count = case_when(is.na(count) ~ 0,
                                                              !is.na(count) ~ count),
                                            count_diff = count - lag(count,1)) %>% filter(!is.na(count_diff))

dt4 <- d_fft %>% filter(PANEL==4) %>% dplyr::select(count,time) %>% mutate(time=round(time,digits = 0))
dt4.1 <- left_join(t,dt4,"time") %>% mutate(count = case_when(is.na(count) ~ 0,
                                                              !is.na(count) ~ count),
                                            count_diff = count - lag(count,1)) %>% filter(!is.na(count_diff))

par(mfrow = c(2, 4))

my.w1 <- analyze.wavelet(dt1.1, "count",
                        loess.span = 0,
                        dt = 1, dj = 1/250,
                        lowerPeriod = 4,
                        upperPeriod = 1024,
                        make.pval = TRUE, n.sim = 10)

wt.image(my.w1, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),main="Trial 1",graphics.reset = FALSE)

my.w2 <- analyze.wavelet(dt2.1, "count",
                         loess.span = 0,
                         dt = 1, dj = 1/250,
                         lowerPeriod = 4,
                         upperPeriod = 1024,
                         make.pval = TRUE, n.sim = 10)

wt.image(my.w2, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),main="Trial 2",graphics.reset = FALSE)

my.w3 <- analyze.wavelet(dt3.1, "count",
                         loess.span = 0,
                         dt = 1, dj = 1/250,
                         lowerPeriod = 4,
                         upperPeriod = 1024,
                         make.pval = TRUE, n.sim = 100)

wt.image(my.w3, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),main="Trial 3",graphics.reset = FALSE)

my.w4 <- analyze.wavelet(dt4.1, "count",
                         loess.span = 0,
                         dt = 1, dj = 1/250,
                         lowerPeriod = 4,
                         upperPeriod = 1024,
                         make.pval = TRUE, n.sim = 100)

wt.image(my.w4, color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels", mar = 4.7),main="Trial 4",graphics.reset = FALSE)

maximum.level = 1.001*max(my.w1$Power.avg, my.w2$Power.avg,my.w3$Power.avg,my.w4$Power.avg)
wt.avg(my.w1, maximum.level = maximum.level)
wt.avg(my.w2, maximum.level = maximum.level)
wt.avg(my.w3, maximum.level = maximum.level)
wt.avg(my.w4, maximum.level = maximum.level)

# Figure S2
x_size <- seq(0,9)
x_lab <- as.character(seq(0,9))

xydat %>% dplyr::filter(num_pulse != 0,stim_duration == 0) %>%
  ggplot(aes(x=num_pulse,fill=as.factor(num_pulse))) +
  geom_bar(aes(y=after_stat(count))) +
  theme_minimal(base_size = 20) +
  scale_x_continuous(breaks = x_size,labels = x_lab) +
  xlab("Number of pulses in a display") + ylab("Count") + theme(axis.title = element_text(face = "bold"),legend.position = "none")

# Figure S3
xydat %>% dplyr::filter(num_pulse != 0,stim_duration != 0) %>%
  ggplot(aes(x=num_pulse,fill=as.factor(num_pulse))) +
  geom_bar(aes(y=after_stat(count))) +
  theme_minimal(base_size = 20) +
  scale_x_continuous(breaks = x_size,labels = x_lab) +
  xlab("Number of pulses in a display") + ylab("Count") + theme(axis.title = element_text(face = "bold"),legend.position = "none") +
  facet_wrap(stim_duration~.)


# Figure S4
d3 <- d2 %>% dplyr::group_by(stimulus) %>% dplyr::summarise(max_radiance = max(radiance),
                                                            lambda_max = Wavelength[which(radiance == max(radiance))],
                                                            stimulus_power = stimulus_power[which(radiance == max(radiance))],
                                                            stimulus_distance = stimulus_distance[which(radiance == max(radiance))],
                                                            stimulus_total = sum(radiance)
) %>%
  filter(stimulus_distance == "close")

##custom function for plotting
nm_to_RGB <- function(wavelengths){
  sapply(wavelengths, function(wavelength) {
    red <- green <- blue <- 0
    if((wavelength >= 380) & (wavelength < 440)){
      red <- -(wavelength - 440) / (440 - 380)
      blue <- 1
    }else if((wavelength >= 440) & (wavelength<490)){
      green <- (wavelength - 440) / (490 - 440)
      blue <- 1
    }else if((wavelength >= 490) && (wavelength<510)){
      green <- 1
      blue = -(wavelength - 510) / (510 - 490)
    }else if((wavelength >= 510) && (wavelength<580)){
      red = (wavelength - 510) / (580 - 510)
      green <- 1
    }else if((wavelength >= 580) && (wavelength<645)){
      red = 1
      green <- -(wavelength - 645) / (645 - 580)
    }else if((wavelength >= 645) && (wavelength<781)){
      red = 1
    }
    if((wavelength >= 380) && (wavelength<420)){
      fac <- 0.3 + 0.7*(wavelength - 380) / (420 - 380)
    }else if((wavelength >= 420) && (wavelength<701)){
      fac <- 1
    }else if((wavelength >= 701) && (wavelength<781)){
      fac <- 0.3 + 0.7*(780 - wavelength) / (780 - 700)
    }else{
      fac <- 0
    }
    do.call(rgb, as.list((c(red, green, blue) * fac)^0.8))
  })
}

d3

# plot stimuli by differing wavelengths
d2 %>% filter(stimulus_power %in% c(5,25,50,100,200,255), stimulus_distance == "close") %>%
  ggplot(aes(x=Wavelength,y=radiance)) + geom_line(aes(linetype=as.factor(stimulus_power))) +
  geom_point(aes(color=nm_to_RGB(Wavelength),shape=as.factor(stimulus_power)),size=1.25) + scale_color_identity() +
  xlab("Wavelength (nm)") + ylab("Radiance (W sr^-1 m^-2 s)") +
  theme_minimal(base_size = 20) + theme(legend.position = "none") #+

# plot stimuli by differing power levels
d2 %>% filter(stimulus_power %in% c(5,25,50,100,200,255), stimulus_distance == "close") %>%
  group_by(stimulus_power) %>%
  ggplot(aes(x=stimulus_power,y=radiance)) +  geom_line(aes(color=nm_to_RGB(Wavelength)),alpha=0.2) +
  #geom_point(aes(color=nm_to_RGB(Wavelength)),alpha=0.2) +
  geom_point(data=d3,aes(x=stimulus_power,y=max_radiance,color=nm_to_RGB(lambda_max)),size = 4) +
  scale_color_identity() +
  xlab("Stimulus power") + ylab("Radiance (W sr^-1 m^-2 s) per wavelength") +
  theme_minimal(base_size = 20)


## Figure S5
d2 <- stim_dat %>% gather(stimulus,radiance,tod_below_50_units:tod_middle_50_units,factor_key = TRUE) %>%
  mutate(stimulus_power =
           case_when(
             stimulus == "tod_below_50_units" ~ 50,
             stimulus == "tod_255_units" ~ 255,
             stimulus == "tod_200_units" ~ 200,
             stimulus == "tod_100_units" ~ 100,
             stimulus == "tod_50_units" ~ 50,
             stimulus == "tod_25_units" ~ 25,
             stimulus == "tod_5_units" ~ 5,
             stimulus == "tod_middle_255_units" ~ 255,
             stimulus == "tod_middle_50_units" ~ 50,
           ),
         stimulus_distance =
           case_when(
             stimulus == "tod_below_50_units" ~ "far",
             stimulus == "tod_middle_255_units" ~ "middle",
             stimulus == "tod_middle_50_units" ~ "middle",
             TRUE ~ "close"
           )
  )

f0 <- d2 %>% filter(stimulus_power <= 50) %>% do(tidy(lm(radiance~0 + stimulus_power*Wavelength,data=.)))
f50 <- d2 %>% filter(stimulus_power >=50 & stimulus_power <= 100) %>% do(tidy(lm(radiance~stimulus_power*Wavelength,data=.)))
f100 <- d2 %>% filter(stimulus_power >=100 & stimulus_power <= 200) %>% do(tidy(lm(radiance~stimulus_power*Wavelength,data=.)))
f200 <- d2 %>% filter(stimulus_power >=200 & stimulus_power <= 255) %>% do(tidy(lm(radiance~stimulus_power*Wavelength,data=.)))
keyz <- sort(unique(int$Brightness))
f0_list <- keyz[1:5]
f50_list <- keyz[6:8]
f100_list <- keyz[9:10]
f200_list <- keyz[11]

newval <- f0$estimate[1]*f0_list
newval1 <- f50$estimate[1] + f50$estimate[2]*f50_list
newval2 <- f100$estimate[1] + f100$estimate[2]*f100_list
newval3 <- f200$estimate[1] + f200$estimate[2]*f200_list

int2 <- int %>% mutate(radiance =
                         case_when(Brightness == 5 ~ newval[1],
                                   Brightness == 10 ~ newval[2],
                                   Brightness == 20 ~ newval[3],
                                   Brightness == 30 ~ newval[4],
                                   Brightness == 50 ~ newval[5],
                                   Brightness == 60 ~ newval1[1],
                                   Brightness == 70 ~ newval1[2],
                                   Brightness == 90 ~ newval1[3],
                                   Brightness == 120 ~ newval2[1],
                                   Brightness == 180 ~ newval2[2],
                                   Brightness == 255 ~ newval3[1]
                         )
)

im = lm(Response ~ log10(radiance) + Tank/Tank_order + Trial,data=int2) #trial order is not significant
summary(im); anova(im)
check_model(im)
effect_plot(im,radiance,data=int2)
report(im)
tab_model(im,dv.labels = c("# of displays"))

figS5 <- ggplot(data = int2, aes(x=radiance,y=Response)) +
  geom_jitter(color="#56B4E9",alpha=0.75,width = 1/11,size=8) +
  geom_boxplot(aes(group=as.factor(log10(radiance))),outlier.shape = NA,fill=NA,width=1/11,color="#0072B2",size=2) +
  stat_smooth(method = "lm",se=FALSE,fullrange = TRUE) +
  theme_minimal(base_size = 20) +
  xlab('Estimated radiance (W sr^-1 m^-2 s)') +
  ylab('# of displays') +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "b",scaled = TRUE,short = grid::unit(4,"mm"),mid = unit(6,"mm"),long = unit(8,"mm"),colour = "#E2E2E2",size=1) +
  theme(axis.title = element_text(size=50,face = "bold"),
        axis.text = element_text(size=30),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  EnvStats::stat_n_text(color = "#0072B2",label.padding = 0.5,size=8)


## Figure S6
delay_dat_during %>% rename(Distance = eu_dist) %>%
  ggplot(aes(x=stim_duration,y=delay)) +
  geom_abline(slope=1, intercept=0,color="darkgrey",linetype=2,linewidth = 1.25) +
  geom_jitter(aes(color=Distance),size=8,shape=21,stroke=2) +
  scale_color_gradient(low="#F0E442",high="#009E73",guide="colourbar") +
  geom_boxplot(aes(group=as.factor(stim_duration)),color="#0072B2",outlier.shape = NA,fill=NA,size=2) +
  stat_smooth(method = "lm",formula = y ~ x,se = FALSE,fullrange = TRUE) +
  theme_minimal(base_size = 20) +
  theme(axis.title = element_text(size=50,face = "bold"),
        axis.text = element_text(size=30),
        legend.title = element_text(face="bold",size=30),
        legend.text = element_text(size=30),
        legend.position = "right") +
  ylim(c(0,12)) +
  xlab("Duration of stimulus (s)") +
  ylab("Latency to display:\nsignal start - stimulus start (s)") +
  guides(color=guide_legend(title = "Distance\nfrom\nstimulus\n(cm)")) + facet_wrap(.~order)

## Figure S7
old_dat3 %>% dplyr::filter(signal_onset == "After",signal_timing == "signal_end_delay") %>%
  ggplot(aes(y=log10(delay),x=stim_duration)) +
  geom_jitter(aes(color=eu_dist,shape = delay > 2.5),alpha=0.75,size=6,stroke=2) +
  scale_color_gradient(low="#F0E442",high="#009E73",guide="colourbar") +
  scale_shape_manual(values = c(19,21)) +
  geom_boxplot(aes(group=stim_duration),outlier.shape = NA,fill=NA,color="#0072B2",size=1) +
  stat_smooth(method = "lm", formula = y ~ poly(x,2),se=FALSE,fullrange = TRUE) +
  theme_minimal(base_size = 20) +
  theme(axis.title = element_text(face = "bold")) +
  EnvStats::stat_n_text(color = "#0072B2",size=8) +
  xlab("Duration of stimulus (s)") + ylab("Latency to display:\nlog10( signal start - stimulus end ) (s)") +
  facet_wrap(.~order)

old_dat3 %>% dplyr::filter(signal_onset == "After",signal_timing == "signal_end_delay") %>%
  do(tidy(lm(log10(delay)~poly(stim_duration,2) + order + eu_dist,data=.)))

## Figure S8
ggpairs(nn_np_all[,c(1,3,4)])

## Figure S9
ggpairs(all_np_dat_nn[,c(1,2,7,8)])

