# species_compa_analysis.R
## script to analyze the NutNet LBK site species comp data
## data spans from 2018 to 2023 
## author : Hannah German 

## load packages
library(tidyverse)
library(ggplot2)  ## do not install again queen 
library(dplyr)
library(lme4)
library(car)
library(emmeans)
library(vegan)
library(multcompView)
library(multcomp)
library(readxl)

## load data
spcomp_data <- read.csv('../../NutNet_LBK/Data/SPP_comp/species_comp_withinfo.csv')

## remove litter and bareground estimates
spcomp_data_plants <- subset(spcomp_data, binomial != 'litter' & binomial != 'bareground')
head(spcomp_data_plants)
### Assessment of % cover within this datasheet due to outliers
#### P15 2022: D=1666, E=555, R=3; %Cover=1, 2, 1, 70
#### P23 2022: D=1428, E=357, R=4; %Cover=2, 7
#### P24 2020: D=1111, E=1111, R=1; %Cover=3
#### P33 2022: D=1000, E=250, R=4; %Cover=1, 1, 2, 2, 15
#### P42 2022: D=909, E=303, R=3; %Cover= 1, 1, 3

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

### Mean entities in which the same plot was sampled twice in the same year. ex) plot 1, 2018, DOY 134 and plot 1, 2018, DOY 289
Summ_spcomp_diversity_plottype <- spcomp_diversity_plottype %>% group_by(trt, Year) %>%
         summarise(diversity=mean(diversity,na.rm = T), evenness=mean(evenness, na.rm =T),richness = mean(richness, na.rm=T))

Summ_spcomp_diversity_ptype_wPlots <- spcomp_diversity_plottype %>% group_by(Plot, Year, trt) %>%
         summarise(diversity=mean(diversity,na.rm = T), evenness=mean(evenness, na.rm =T),richness = mean(richness, na.rm=T)) %>%
         mutate(n = ifelse(trt == "N" | trt == "NP" | trt == "NK" | trt == "NPK", 1, 0),
                p = ifelse(trt == "P" | trt == "NP" | trt == "PK" | trt == "NPK", 1, 0),
                k = ifelse(trt == "K" | trt == "NK" | trt == "PK" | trt == "NPK", 1, 0))
plot(Summ_spcomp_diversity_plottype)

### initial plots to asses basic trends before doing stats
ggplot (Summ_spcomp_diversity_ptype_wPlots, aes(trt, diversity, color=factor(Year))) + 
        geom_point()

ggplot (subset(Summ_spcomp_diversity_ptype_wPlots, diversity<100 & trt!= 'Fence'& trt != 'NPK+Fence'& trt != 'xControl'), 
        aes(Year, diversity, fill=factor(trt))) + geom_point() + geom_bar(stat = "identity",position = "dodge") +  facet_wrap (~trt)

ggplot (subset(Summ_spcomp_diversity_ptype_wPlots, diversity<100 & trt!= 'Fence'& trt != 'NPK+Fence'& trt != 'xControl'), 
        aes(Year, evenness, fill=factor(trt))) + geom_point() + geom_bar(stat = "identity",position = "dodge") + facet_wrap (~trt)

ggplot (subset(Summ_spcomp_diversity_ptype_wPlots, diversity<100 & trt!= 'Fence'& trt != 'NPK+Fence'& trt != 'xControl'), 
        aes(Year, richness, fill=factor(trt))) + geom_point() + geom_bar(stat = "identity",position = "dodge") + facet_wrap (~trt) + theme_bw()

### model making time 
#### add in treatment binary factors
Summ_spcomp_diversity_ptype_wPlots$nfac <- as.factor(Summ_spcomp_diversity_ptype_wPlots$n)
Summ_spcomp_diversity_ptype_wPlots$pfac <- as.factor(Summ_spcomp_diversity_ptype_wPlots$p)
Summ_spcomp_diversity_ptype_wPlots$kfac <- as.factor(Summ_spcomp_diversity_ptype_wPlots$k)
Summ_spcomp_diversity_ptype_wPlots$plotfac <- as.factor(Summ_spcomp_diversity_ptype_wPlots$Plot)
Summ_spcomp_diversity_ptype_wPlots$yearfac <- as.factor(Summ_spcomp_diversity_ptype_wPlots$Year)
#### add in blocks
Summ_spcomp_diversity_ptype_wPlots$block <- 'block2'
Summ_spcomp_diversity_ptype_wPlots$block[Summ_spcomp_diversity_ptype_wPlots$Plot <15] <- 'block1'
Summ_spcomp_diversity_ptype_wPlots$block[Summ_spcomp_diversity_ptype_wPlots$Plot >28] <- 'block3'

#### remove certain plot types
spcomp_data_4lmer <- subset(Summ_spcomp_diversity_ptype_wPlots, trt!= 'Fence'& trt != 'NPK+Fence'& trt != 'xControl')

#### model time part 2
mod_div.year.trt <- lmer(log(diversity) ~ yearfac * nfac * pfac * kfac + (1| plotfac) + (1|block), data = (spcomp_data_4lmer))  
plot(resid(mod_div.year.trt) ~ fitted(mod_div.year.trt))
Anova(mod_div.year.trt) 
cld(emmeans(mod_div.year.trt, ~yearfac))
cld(emmeans(mod_div.year.trt, ~nfac*kfac))

mod_rich.year.trt <- lmer((richness) ~ yearfac * nfac * pfac * kfac + (1| plotfac) + (1|block), data = (spcomp_data_4lmer))  
plot(resid(mod_rich.year.trt) ~ fitted(mod_rich.year.trt))
Anova(mod_rich.year.trt)  
cld(emmeans(mod_rich.year.trt, ~yearfac*nfac))
cld(emmeans(mod_rich.year.trt, ~nfac*pfac))

mod_evenness.year.trt <- lmer(log(evenness) ~ yearfac * nfac * pfac * kfac + (1| plotfac) + (1|block), data = (spcomp_data_4lmer))  
plot(resid(mod_evenness.year.trt) ~ fitted(mod_evenness.year.trt))
Anova(mod_evenness.year.trt)  
cld(emmeans(mod_evenness.year.trt, ~yearfac))

## climate 
KLBB_weather <- read_excel("~/Documents/Git/HG_NutNEt_lbk_spcomp/data/KLBB_weather.xlsx")
View(KLBB_weather)
weather <- lmer(Year ~ precip_mm * temp_C + (1| Month), data = (KLBB_weather))
plot(weather)
Anova(weather)

ggplot (KLBB_weather, aes(Year, precip_mm, fill=factor(Year))) + geom_point() + geom_bar(stat = "identity",position = "dodge") + facet_wrap (~Year) + theme_bw()

## adding a comment
