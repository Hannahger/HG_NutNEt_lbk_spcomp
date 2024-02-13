# species_comp_analysis.R
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
library(lubridate)
library(ggpubr)
library(patchwork)

## load data
spcomp_data <- read.csv("../data/species_comp.csv")

## remove litter and bareground estimates
spcomp_data_plants <- subset(spcomp_data, binomial != 'litter' & binomial != 'bareground')
head(spcomp_data_plants)
### Assessment of % cover within this datasheet due to outliers   ## No longer a concern (1/23) ##
#### P15 2022: D=1666, E=555, R=3; %Cover=1, 2, 1, 70
#### P23 2022: D=1428, E=357, R=4; %Cover=2, 7
#### P24 2020: D=1111, E=1111, R=1; %Cover=3
#### P33 2022: D=1000, E=250, R=4; %Cover=1, 1, 2, 2, 15
#### P42 2022: D=909, E=303, R=3; %Cover= 1, 1, 3

## square percent covers for diversity estimate
spcomp_data_plants$Percent.Cover.squared <- (spcomp_data_plants$Percent.Cover/100) * (spcomp_data_plants$Percent.Cover/100)

## calculate diversity indices

### group by plot per year per day of year
spcomp_data_groupby_plot_year <- group_by(spcomp_data_plants, Plot, Year, DOY) # NGS: don't add species information here as we are summarising by plot

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
# Summ_spcomp_diversity_plottype <- spcomp_diversity_plottype %>% group_by(trt, Year) %>%
#          summarise(diversity=mean(diversity,na.rm = T), evenness=mean(evenness, na.rm =T),richness = mean(richness, na.rm=T))

### remove DOY 159 from year 1 (i.e., 2018)
spcomp_diversity_plottype_fall <- subset(spcomp_diversity_plottype, DOY > 160)

Summ_spcomp_diversity_ptype_wPlots <- spcomp_diversity_plottype_fall %>%
         mutate(n = ifelse(trt == "N" | trt == "NP" | trt == "NK" | trt == "NPK" | trt == "NPK+Fence", 1, 0),
                p = ifelse(trt == "P" | trt == "NP" | trt == "PK" | trt == "NPK" | trt == "NPK+Fence", 1, 0),
                k = ifelse(trt == "K" | trt == "NK" | trt == "PK" | trt == "NPK" | trt == "NPK+Fence", 1, 0))
# plot(Summ_spcomp_diversity_ptype_wPlots)

### initial plots to asses basic trends before doing stats
# ggplot (Summ_spcomp_diversity_ptype_wPlots, aes(trt, diversity, color=factor(Year))) + geom_point()
# 
# ggplot (subset(Summ_spcomp_diversity_ptype_wPlots, diversity<100 & trt!= 'Fence'& trt != 'NPK+Fence'& trt != 'xControl'), 
#         aes(Year, diversity, fill=factor(trt))) + geom_point() + geom_bar(stat = "identity",position = "dodge") +  facet_wrap (~trt)
# 
# ggplot (subset(Summ_spcomp_diversity_ptype_wPlots, diversity<100 & trt!= 'Fence'& trt != 'NPK+Fence'& trt != 'xControl'), 
#         aes(Year, evenness, fill=factor(trt))) + geom_point() + geom_bar(stat = "identity",position = "dodge") + facet_wrap (~trt)
# 
# ggplot (subset(Summ_spcomp_diversity_ptype_wPlots, diversity<100 & trt!= 'Fence'& trt != 'NPK+Fence'& trt != 'xControl'), 
#         aes(Year, richness, fill=factor(trt))) + geom_point() + geom_bar(stat = "identity",position = "dodge") + facet_wrap (~trt) + theme_bw()


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

#### statistical models of diversity across years
# testing the hypotheses about treatment impacts on community diversity, accounting for year-to-year differences
mod_div.year.trt <- lmer(log(diversity) ~ yearfac * nfac * pfac * kfac + (1| plotfac) + (1|block), 
                         data = (spcomp_data_4lmer))  
plot(resid(mod_div.year.trt) ~ fitted(mod_div.year.trt)) # might need to deal with normality issues
Anova(mod_div.year.trt) 
cld(emmeans(mod_div.year.trt, ~yearfac))
cld(emmeans(mod_div.year.trt, ~nfac*kfac))

mod_rich.year.trt <- lmer(richness ~ yearfac * nfac * pfac * kfac + (1| plotfac) + (1|block), 
                          data = (spcomp_data_4lmer))  
plot(resid(mod_rich.year.trt) ~ fitted(mod_rich.year.trt))
Anova(mod_rich.year.trt)  
cld(emmeans(mod_rich.year.trt, ~yearfac))
cld(emmeans(mod_rich.year.trt, ~yearfac*nfac))
cld(emmeans(mod_rich.year.trt, ~nfac*pfac))

mod_evenness.year.trt <- lmer(log(evenness) ~ yearfac * nfac * pfac * kfac + (1| plotfac) + (1|block), 
                              data = (spcomp_data_4lmer))  
plot(resid(mod_evenness.year.trt) ~ fitted(mod_evenness.year.trt)) # check this
Anova(mod_evenness.year.trt)  
cld(emmeans(mod_evenness.year.trt, ~yearfac))

## look at treatment impacts on each plant type
head(spcomp_data_plants)

### make new column with plant type information
### what plant types do we have?
unique(spcomp_data_plants$lifeform)
unique(spcomp_data_plants$lifespan)
unique(spcomp_data_plants$ps_path)

### figure out the unknown and NAs >> these are all unknown species
subset(spcomp_data_plants, lifeform == 'Unknown')
subset(spcomp_data_plants, lifespan == 'Unknown')
spcomp_data_plants[770:774,]

### make a subset with the unknowns removed
spcomp_data_plants_known <- subset(spcomp_data_plants, Code != 'Unknown')
head(spcomp_data_plants_known)
spcomp_data_plants_known$Code

### make a new column for plant type
#### what plant types do we have? c4_perennial_grass, c4_annual_forb, c3_perennial_forb, c3_annual_forb, c3_perennial_shrub
spcomp_data_plants_known$pft[spcomp_data_plants_known$ps_path == 'C4' & 
                               spcomp_data_plants_known$lifespan == 'annual' &
                               spcomp_data_plants_known$lifeform == 'Grass'] <- 'c4_annual_grass'
spcomp_data_plants_known$pft[spcomp_data_plants_known$ps_path == 'C4' & 
                               spcomp_data_plants_known$lifespan == 'perennial' &
                               spcomp_data_plants_known$lifeform == 'Grass'] <- 'c4_perennial_grass'
spcomp_data_plants_known$pft[spcomp_data_plants_known$ps_path == 'C4' & 
                               spcomp_data_plants_known$lifespan == 'perennial' &
                               spcomp_data_plants_known$lifeform == 'Forb'] <- 'c4_perennial_forb'
spcomp_data_plants_known$pft[spcomp_data_plants_known$ps_path == 'C4' & 
                               spcomp_data_plants_known$lifespan == 'annual' &
                               spcomp_data_plants_known$lifeform == 'Forb'] <- 'c4_annual_forb'
spcomp_data_plants_known$pft[spcomp_data_plants_known$ps_path == 'C3' & 
                               spcomp_data_plants_known$lifespan == 'annual' &
                               spcomp_data_plants_known$lifeform == 'Grass'] <- 'c3_annual_grass'
spcomp_data_plants_known$pft[spcomp_data_plants_known$ps_path == 'C3' & 
                               spcomp_data_plants_known$lifespan == 'perennial' &
                               spcomp_data_plants_known$lifeform == 'Grass'] <- 'c3_perennial_grass'
spcomp_data_plants_known$pft[spcomp_data_plants_known$ps_path == 'C3' & 
                               spcomp_data_plants_known$lifespan == 'perennial' &
                               spcomp_data_plants_known$lifeform == 'Forb'] <- 'c3_perennial_forb'
spcomp_data_plants_known$pft[spcomp_data_plants_known$ps_path == 'C3' & 
                               spcomp_data_plants_known$lifespan == 'annual' &
                               spcomp_data_plants_known$lifeform == 'Forb'] <- 'c3_annual_forb'
spcomp_data_plants_known$pft[spcomp_data_plants_known$ps_path == 'C3' & 
                               spcomp_data_plants_known$lifespan == 'perennial' &
                               spcomp_data_plants_known$lifeform == 'Shrub/tree'] <- 'c3_perennial_woody'

### summarise by pft for each plot in each year

# TAKE HOME: treatments have no impact on any metric of diversity in any year
# but there is significant year-to-year variation, with even (i.e., dry) years having
# the greatest diversity, lowest richness, but highest evenness

### NEXT STEP: summarise each pft by plot within year and fit lmer models (similar to diversity models)

#### statistical models regarding the diversity of plant types across years 
# testing hypothesis that treatment impacts on plant type (annual to outperform perennial, C3 to ____ C4, forbs to outperform grasses)

# make binomial, lifeform, lifespan, and ps_path factors
## leave to go over with NGS; I know what is going on but don't have time at the moment to fix (don't adjust things above)
#head(Summ_spcomp_diversity_ptype_wPlots)
#Summ_spcomp_diversity_ptype_wPlots$binomial_fac <- as.factor(Summ_spcomp_diversity_ptype_wPlots$binomial)
#Summ_spcomp_diversity_ptype_wPlots$lifeform_fac <- as.factor(Summ_spcomp_diversity_ptype_wPlots$lifeform)
#Summ_spcomp_diversity_ptype_wPlots$lifespan_fac <- as.factor(Summ_spcomp_diversity_ptype_wPlots$lifespan)
#Summ_spcomp_diversity_ptype_wPlots$ps_path_fac <- as.factor(Summ_spcomp_diversity_ptype_wPlots$ps_path)


## climate 
KLBB_weather <- read_excel("../data/KLBB_weather.xlsx")
head(KLBB_weather)

annual_precip <- KLBB_weather %>% mutate(Date_Time = ymd_hms(Date_Time), Year = year(Date_Time)) %>% group_by(Year) %>% 
  summarise(annual_precip = sum(precip_mm, na.rm = TRUE))

plot_df <- Summ_spcomp_diversity_ptype_wPlots %>% full_join(annual_precip)

## making figures for TTABSS

# fig theme; Thank you Evan for letting me steal this :)

figtheme <- theme_bw(base_size = 18) +
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.border = element_rect(size = 1.5, fill = NA),
        legend.box.background = element_blank(),
        legend.key = element_rect(fill = NA),
        legend.background=element_blank(),
        legend.title = element_text(face = "bold"),
        axis.ticks.length = unit(0.25, "cm"),
        panel.grid.minor.y = element_blank(),
        legend.text.align = 0)

# fig.2 Significant year-to-year variation between wet and dry years 
fig.2.D <- ggplot() + 
  stat_boxplot(data = subset(plot_df, diversity < 100 & trt != 'Fence' & trt != 'NPK+Fence' & trt != 'xControl'), 
               aes(x = as.factor(Year), y = diversity), size = 0.75, geom = "errorbar", width = 0.2)  +
  geom_boxplot(data = subset(plot_df, diversity < 100 & trt!= 'Fence' & trt != 'NPK+Fence' & trt != 'xControl'),
               aes(as.factor(Year), diversity, fill = annual_precip), outlier.shape = NA) + 
  scale_fill_gradient(low = "#d60404", high = "#0047ab") + 
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 45, 15), name = "Simpson's Diversity") + 
  labs(x = "Year") +   
  figtheme +
  theme(legend.position = "none")

fig.2.R <- ggplot() + 
  stat_boxplot(data = subset(plot_df, richness & trt!= 'Fence' & trt!= 'NPK+Fence' & trt!= 'xControl'),
               aes(as.factor(Year), richness), size = 0.75, geom = "errorbar", width = 0.2) + 
  geom_boxplot(data = subset(plot_df, richness & trt!= 'Fence' & trt != 'NPK+Fence' & trt != 'xControl'),
               aes(as.factor(Year), richness, fill = annual_precip), outlier.shape = NA) + 
  scale_fill_gradient(low = "#d60404", high = "#0047ab") +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, 2.5),  name = "Species Richness") + 
  labs(x = "Year") + 
  labs(fill = "MAP (mm)") +
  figtheme + theme(legend.position = "bottom") +

  theme(legend.key.size = unit(1.5, "cm"), legend.key.height = unit (0.5, "cm"))


fig.2.E <- ggplot() + 
  stat_boxplot(data = subset(plot_df, evenness & trt!= 'Fence' &  trt != 'NPK+Fence' & trt != 'xControl'),
               aes(as.factor(Year), evenness), size = 0.75, geom = "errorbar", width = 0.2) +  
  geom_boxplot(data = subset(plot_df, evenness & trt!= 'Fence' & trt != 'NPK+Fence' & trt != 'xControl'),
               aes(as.factor(Year), evenness, fill = annual_precip), outlier.shape = NA) + 
  scale_fill_gradient(low = "#d60404", high = "#0047ab") +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5), name = "Species Evenness") + 
  labs(x = "Year") + 
  figtheme +
  theme(legend.position = "none")


# fig.1 Treatment had no effect on any diversity metric 

### EAP note for Hannah: changed trt factor levels to reflect order of NPK 
### treatments in a bit more of an intuitive order. 
### Previous order: Control, K, N, NK, NP, NPK, P, PK
### New order:  Control, N, P, K, NP, NK, PK, NPK)
plot_df$trt <- factor(plot_df$trt, levels = c("Control", "N", "P", "K",
                                              "NP", "NK", "PK", "NPK",
                                              "xControl", "Fence", "NPK+Fence"))

fig.1.D <- ggplot() + 
  stat_boxplot(data = subset(plot_df, diversity < 60 & trt!= 'Fence' & trt != 'NPK+Fence' & trt != 'xControl'),
              aes(as.factor(trt), diversity), size = 0.75, geom = "errorbar", width = 0.2) +
  geom_boxplot(data = subset(plot_df, diversity < 60 & trt!= 'Fence' & trt != 'NPK+Fence' & trt != 'xControl'),
              aes(as.factor(trt), diversity, fill = as.factor(trt)), outlier.shape = NA) + 
  scale_fill_manual(values = c("gray", "#bb5566", "#bb5566", "#bb5566", "#bb5566",
                               "#bb5566", "#bb5566", "#bb5566", "#bb5566")) + 
  theme(legend.position = "none") + 
  labs(x = "Treatment",y = "Simpson's Diversity") + 
  scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  figtheme + 
  theme(legend.position = "none")

fig.1.R <- ggplot() +    #getting message "Can't add `scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5))` to a theme object."
  stat_boxplot(data = subset(plot_df, richness & trt!= 'Fence' &  trt != 'NPK+Fence' & trt != 'xControl'),
              aes(as.factor(trt), richness), size = 0.75, geom = "errorbar", width = 0.2) +
  geom_boxplot(data = subset(plot_df, richness & trt!= 'Fence' & trt != 'NPK+Fence' & trt != 'xControl'),
              aes(as.factor(trt), richness, fill = as.factor(trt)), outlier.shape = NA,) + 
  scale_fill_manual(values = c("gray", "#bb5566", "#bb5566", "#bb5566", "#bb5566",
                               "#bb5566", "#bb5566", "#bb5566", "#bb5566"))
  theme(legend.position = "none") + 
  labs(x = "Treatment", y = "Species Richness") + 
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  figtheme +
  theme(legend.position = "none")

fig.1.E <- ggplot() +     #getting message "Can't add `scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5))` to a theme object."
  stat_boxplot(data = subset(plot_df, evenness < 10 & trt!= 'Fence' & trt != 'NPK+Fence' & trt != 'xControl'),
              aes(as.factor(trt), evenness), size = 0.75, geom = "errorbar", width = 0.2) +
  geom_boxplot(data = subset(plot_df, evenness < 10 & trt!= 'Fence' & trt != 'NPK+Fence' & trt != 'xControl'),
              aes(as.factor(trt), evenness, fill = as.factor(trt)), outlier.shape = NA) +
  scale_fill_manual(values = c("gray", "#bb5566", "#bb5566", "#bb5566", "#bb5566",
                             "#bb5566", "#bb5566", "#bb5566", "#bb5566"))
  theme(legend.position = "none") + 
  labs(x = "Treatment", y = "Species Evenness") + 
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2.5)) +
  figtheme + 
  theme(legend.position = "none")

  
master.fig <- ggarrange(fig.1.D, fig.2.D, fig.2.R, fig.2.E)
  
## download png's of figures 

png("../plots/master.fig.png", 
    width = 9, height = 8, units = 'in', res = 600)
master.fig
dev.off()



