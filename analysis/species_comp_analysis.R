# species_compa_analysis.R
## script to analyze the NutNet LBK site species comp data
## data spans from 2018 to 2023 
## author : Hannah German 

## load packages
library(tidyverse)
library(ggplot2)  ## do not install again queen 
library(dplyr)

## load data
spcomp_data <- read.csv('../../NutNet_LBK/Data/SPP_comp/species_comp_withinfo.csv')

## remove litter and bareground estimates
spcomp_data_plants <- subset(spcomp_data, binomial != 'litter' & binomial != 'bareground')
head(spcomp_data_plants)

## square percent covers for diversity estimate
spcomp_data_plants$Percent.Cover.squared <- (spcomp_data_plants$Percent.Cover/100) * (spcomp_data_plants$Percent.Cover/100)

## calculate diversity indices

### group by plot per year per day of year
spcomp_data_groupby_plot_year <- group_by(spcomp_data_plants, Plot, Year, DOY)

### richness per plot per year per day of year
spcomp_data_plot_year_richness <- summarise(spcomp_data_groupby_plot_year, richness = n_distinct(binomial))
head(spcomp_data_plot_year_richness)

### diversity per year per day of year
spcomp_data_plot_year_diversity <- summarise(spcomp_data_groupby_plot_year, diversity = 1/sum(Percent.Cover.squared))
head(spcomp_data_plot_year_diversity)

### combine richness and diversity datasets and calculate evenness per year per day of year
spcomp_diversity <- left_join(spcomp_data_plot_year_richness, spcomp_data_plot_year_diversity)
head(spcomp_diversity)
spcomp_diversity$evenness <- spcomp_diversity$diversity/spcomp_diversity$richness
head(spcomp_diversity)

### add in plot treatment information
plot_information <- read.csv('../../NutNet_LBK/plot_types/plot_types.csv')[,2:3]
head(plot_information)
spcomp_diversity_plottype <- left_join(spcomp_diversity, plot_information)
head(spcomp_diversity_plottype)









### first attempt at analysis; was not good bud
what <-spcomp_diversity_plottype%>%
  group_by(Plot)
ggplot (spcomp_diversity_plottype, aes(trt, diversity, color=factor(Year)))+
  geom_point()
Summ_spcomp_diversity_plottype <- spcomp_diversity_plottype %>% group_by(Plot)%>%
  dplyr::select(Plot, Year, diversity) 
AAH <- Summ_spcomp_diversity_plottype %>% group_by(Plot)%>% 
  summarise(diversity=mean(diversity,na.rm = T)) %>% ungroup() ## did not work
