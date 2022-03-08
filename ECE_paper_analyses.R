## ___________________________________________________________________________##
## Author: Berengere Husson
## Year: 2021
## ___________________________________________________________________________##


## *********************************** ##
##  Successive ECE in Barents Sea      ##
## *********************************** ##

## Libraries -------------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggthemes)
library(tidyverse)
library(patchwork)
library(FactoMineR)
library(factoextra)
library(readxl)
library(ggcorrplot)
library(yarrr)
library(marmap)

## Load data -------------------------------------------------------------------
# list of species with at least 5% occurrences in all samples
species_pnp <- read.csv("input_data/species_list.csv", header = T)

# Barents Sea ecosystem survey aggregated data
load("input_data/pnp_data.RData")

# Traits from  "A trait collection of marine fish species from North Atlantic and 
# Northeast Pacific continental shelf seas." PANGAEA, 
# https://doi.org/10.1594/PANGAEA.900866, Beukhof et al., 2019
traits_Beukhof <- read_excel(
  "input_data/TraitCollectionFishNAtlanticNEPacificContShelf_Beukhofetal2019.xlsx")
traits_Beukhof

# Potential niche descriptors (from Husson et al. 2020)
# https://onlinelibrary.wiley.com/doi/epdf/10.1111/fog.12493 
load("input_data/spe_pref_dml.RData")

# Kernel of cluster's distribution
load("input_data/kernel_density_cluster_distrib.rds")
# Ice, heat and freshwater data 
ice_june <- read.csv("input_data/SeaIceExtent_june.csv")
names(ice_june) <- c("Year", "Arctic Ocean-Sea ice extent (km²)",
                     "(a) Barents Sea\nSIE (km²)")
ice_sep <- read.csv("input_data/SeaIceExtent_sep.csv")
names(ice_sep) <- c("Year", "Arctic Ocean-Sea ice extent (km²)",
                    "(a) Barents Sea\nSIE (km²)")
heat <- read.table("input_data/EnvMat.dat", sep=",", header = T)

# Format data ------------------------------------------------------------------
# Select traits in Barents Sea
traits <- traits_Beukhof %>% 
  # We select species with at least 5% occurrence
  # Icelus is not available for the Barents Sea. 
  # we take from the East Bering Sea instead
  filter( (taxon =="Icelus"& LME==2) | 
            (taxon %in% species_pnp[,2] & LME %in% c(20) ))
length(unique(traits$taxon))# number of species
names(traits)

# Trait slection based on species ability to react to short term perturbations
traits <- traits %>% 
  dplyr::select(!contains(c("level","reference", "family","genus","species",
                            "taxonomic.rank","LME","FAO", "infinity", 
                            "maturity","AR","growth"))) %>%
  distinct() %>% filter(complete.cases(.))

# Niche descriptors selection based on most limiting factors, according to 
# Husson et al. 2020
niche_data2 <- spe_niche_desc %>%  
  dplyr::select(-contains(c("SML","slope","min","max","tot_range",
                            "q025","q975"))) %>% 
  rename_with(~gsub("main_", "", .x, fixed = TRUE))

# Change ice data to long format
ice_june_lg <- ice_june %>% 
  pivot_longer(cols=-Year, names_to = "part",
               values_to = "km2")%>%
  mutate("Month"="June")
ice_sep_lg <- ice_sep %>% 
  pivot_longer(cols=-Year, names_to = "part",
               values_to = "km2")%>%
  mutate("Month"="September")

ice <- rbind.data.frame(ice_june_lg,ice_sep_lg) %>% na.omit()
ice$part <- factor(ice$part)

# Change heat and freshwater format
heat <- heat %>% 
  mutate(Heat_content_million = Heat_content_u100m/10^6)
heat_lg <- heat[, -c(2,4)]%>%
  pivot_longer(cols=-Year, names_to="variable",values_to="value")
heat_lg <- heat_lg %>% group_by(variable)%>%
  mutate(Variable_2=factor(variable))%>%
  mutate(Variable_2=fct_recode(Variable_2,
                               "(c) Heat content (MJ/m², 0-100m)"=
                                 "Heat_content_million",
                               "(b) Freshwater content (0-100m)"=
                                 "Freshwater_content_u100m"))
heat_lg$Variable_2 <- factor(heat_lg$Variable_2, 
                             levels=c("(b) Freshwater content (0-100m)",
                                      "(c) Heat content (MJ/m², 0-100m)"))

## ---------
## Figure 1 -------------------------------------------------------------------#
## ---------
## See matlab script by Sigrid Lind on the same github
## ---------
## Figure 2 -------------------------------------------------------------------#
## ---------

# Figure 2a : Sea ice time series, and ECE detection ---------------------------
# 1. Fit a loess to the time series to cature the non-linear trend
ice.lo <- loess(km2/10^6  ~ Year, 
                data= ice[ice$part %in% "(a) Barents Sea\nSIE (km²)"  & 
                            ice$Month == "June", ], span = 0.75)
ice.lo.pred <- predict(
  ice.lo,
  data.frame(Year = ice$Year[ice$part %in% "(a) Barents Sea\nSIE (km²)"  & 
                               ice$Month == "June" ]), 
  se = TRUE)

# 2. Extract trend residuals to detect ECEs (<5th percentile or>95th percentile)
resid.ice <- ice$km2[ice$part %in% "(a) Barents Sea\nSIE (km²)"   & 
                       ice$Month == "June"]/10^6 - ice.lo.pred$fit
ece.ice <- (resid.ice < quantile(resid.ice, 0.05, na.rm = T))|
  (resid.ice > quantile(resid.ice, 0.95, na.rm = T))

# 3. Join in one dataframe
ice.lo.df <- cbind.data.frame(Year= ice$Year[ice$part %in% 
                                               "(a) Barents Sea\nSIE (km²)"  & 
                                               ice$Month == "June"],
                              part = "(a) Barents Sea\nSIE (km²)",
                              Month ="June",
                              km2 = ice.lo.pred$fit,
                              SE = ice.lo.pred$se.fit,
                              fifth= quantile(resid.ice, 0.05, na.rm = T),
                              nientyfifth= quantile(resid.ice, 0.95, na.rm = T),
                              ext_val = ice$km2[ice$part %in% 
                                                  "(a) Barents Sea\nSIE (km²)"& 
                                                  ice$Month == "June"]/10^6,
                              extreme =ece.ice)

# 4. Plot 
ann_text_ice <- cbind.data.frame(Year=c(-Inf), km2=c(Inf), 
                                 part=c("(a) Barents Sea\nSIE (km²)"),
                                 Month=rep("September",1),
                                 lab=c("\n(a) Sea ice extent (million km²)"))
ice_BS <- ggplot()+ 
  # study period of the paper
  geom_rect(data= ice %>% filter(part %in% "(a) Barents Sea\nSIE (km²)"),
            aes(xmin=2004, xmax=2017, ymin=-Inf, ymax=Inf), 
            fill="grey90", alpha=0.05)+
  # SIE June data
  geom_line(data= ice %>% filter(part %in% "(a) Barents Sea\nSIE (km²)"), 
            aes(x=Year, y=km2/10^6, linetype=Month))+
  theme_pander(base_size=9)+
  # the 3 bottom temperature events
  geom_vline(xintercept = c(2006,2012,2016), color="grey60", linetype=2)+
  # non linear trend and confidence interval
  geom_line(data= ice.lo.df,aes(x=Year, y=km2), color="red")+
  geom_ribbon(data= ice.lo.df,aes(x=Year, ymin=km2 -(2*SE), ymax=km2 +(2*SE)), 
              fill="red", alpha=0.2)+
  # Upper and lower 5th percentile and detected ECEs
  geom_line(data= ice.lo.df,aes(x=Year, y=km2+fifth), 
            color="red", linetype=2)+
  geom_line(data= ice.lo.df,aes(x=Year, y=km2+nientyfifth), 
            color="red", linetype=2)+
  geom_point(data= ice.lo.df,aes(x=Year, y=ext_val,shape=extreme), 
             size=2, color="red")+
  # Esthetics
  geom_text(data=ann_text_ice,
            aes(x=Year, y=km2/10^6,label=lab),hjust = -0.1, size=3)+
  labs(x="", y="Sea ice extent (million km²)", linetype="")+ 
  scale_shape_manual(values = c(NA,16), guide=F)+
  scale_linetype_manual(values=c(1,2),guide=F)+
  xlim(1970,2020)+
  theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_line(colour = "grey60",linetype = 3),
        panel.grid.major.x = element_line(colour = "grey60",linetype = 3),
        panel.border = element_rect(colour = "black"),
        legend.position = c(0.5, 0.525), legend.direction = "horizontal",
        legend.key = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size=9),
        axis.title.x = element_blank(),axis.text.x = element_blank(),
        axis.title.y = element_blank())



# Figure 2b and c : Heat and freshwater time series, and ECE detection ---------
# 1. Fit a loess to the time series to cature the non-linear trend
heat.lo <- loess(value ~ Year, 
                 data= heat_lg[heat_lg$Variable_2 %in% 
                                 "(c) Heat content (MJ/m², 0-100m)", ], 
                 span = 0.75)
heat.lo.pred <- predict(
  heat.lo, 
  data.frame(Year = heat_lg$Year[heat_lg$Variable_2 %in% 
                                   "(c) Heat content (MJ/m², 0-100m)" ]), 
  se = TRUE)

# 2. Extract trend residuals to detect ECEs (<5th percentile or>95th percentile)
resid.heat <- heat_lg$value[heat_lg$Variable_2 %in% 
                              "(c) Heat content (MJ/m², 0-100m)" ] - 
  heat.lo.pred$fit
ece.heat <- (resid.heat > quantile(resid.heat, 0.95, na.rm = T))|
  (resid.heat < quantile(resid.heat, 0.05, na.rm = T))

# 3. Join in one dataframe
heat.lo.df <- cbind.data.frame(Year= heat_lg$Year[heat_lg$Variable_2 %in% "(c) Heat content (MJ/m², 0-100m)"],
                               Variable_2 = "(c) Heat content (MJ/m², 0-100m)",
                               value = heat.lo.pred$fit,
                               SE = heat.lo.pred$se.fit,
                               fifth= quantile(resid.heat, 0.05, na.rm = T),
                               nientyfifth= quantile(resid.heat, 0.95, na.rm = T),
                               ext_val = heat_lg$value[heat_lg$Variable_2 %in% "(c) Heat content (MJ/m², 0-100m)"],
                               extreme =ece.heat)

fresh.lo <- loess(value ~ Year, data= heat_lg[heat_lg$Variable_2 %in% "(b) Freshwater content (0-100m)", ], span = 0.75)
fresh.lo.pred <- predict(fresh.lo, data.frame(Year = heat_lg$Year[heat_lg$Variable_2 %in% "(b) Freshwater content (0-100m)" ]), se = TRUE)
resid.fresh <- heat_lg$value[heat_lg$Variable_2 %in% "(b) Freshwater content (0-100m)" ] - fresh.lo.pred$fit
ece.fresh <- (resid.fresh < quantile(resid.fresh, 0.05, na.rm = T))|(resid.fresh > quantile(resid.fresh, 0.95, na.rm = T))
fresh.lo.df <- cbind.data.frame(Year= heat_lg$Year[heat_lg$Variable_2 %in% "(b) Freshwater content (0-100m)"],
                                Variable_2 = "(b) Freshwater content (0-100m)",
                                value = fresh.lo.pred$fit,
                                SE = fresh.lo.pred$se.fit,
                                fifth= quantile(resid.fresh, 0.05, na.rm = T),
                                nientyfifth= quantile(resid.fresh, 0.95, na.rm = T),
                                ext_val = heat_lg$value[heat_lg$Variable_2 %in% "(b) Freshwater content (0-100m)"],
                                extreme =ece.fresh)

# 4. Plot
ann_text_heat <- cbind.data.frame(Year=rep(-Inf,2), value=c(Inf, -Inf), 
                                  Variable_2=c("(c) Heat content (MJ/m², 0-100m)",
                                               "(b) Freshwater content (0-100m)"),
                                  lab=c("\n(c) Heat content MJ/m² (0-100m)",
                                        "(b) Freshwater content m/m² (0-100m)\n"))

heat_facet <- ggplot()+
  # study period of the paper
  geom_rect(data= heat_lg, aes(xmin=2004, xmax=2017, ymin=-Inf, ymax=Inf),
            fill="grey90", alpha=0.1)+
  # Heat and freshwater data
  geom_line(data= heat_lg,aes(x=Year, y=value))+ 
  facet_wrap(~Variable_2, scales="free_y",ncol=1)+
  theme_pander(base_size = 9)+
  # non linear trend and confidence interval
  geom_line(data= heat.lo.df,aes(x=Year, y=value), color="red")+
  geom_line(data= fresh.lo.df,aes(x=Year, y=value), color="red")+
  geom_ribbon(data= heat.lo.df,aes(x=Year, ymin=value -(2*SE), ymax=value +(2*SE)), fill="red", alpha=0.2)+
  geom_ribbon(data= fresh.lo.df,aes(x=Year, ymin=value -(2*SE), ymax=value +(2*SE)), fill="red", alpha=0.2)+
  # Upper and lower 5th percentile and detected ECEs
  geom_line(data= heat.lo.df,aes(x=Year, y=value+fifth), color="red", linetype=2)+
  geom_line(data= heat.lo.df,aes(x=Year, y=value+nientyfifth), color="red", linetype=2)+
  geom_line(data= fresh.lo.df,aes(x=Year, y=value+fifth), color="red", linetype=2)+
  geom_line(data= fresh.lo.df,aes(x=Year, y=value+nientyfifth), color="red", linetype=2)+
  geom_point(data= heat.lo.df,aes(x=Year, y=ext_val,shape=extreme), size=2, color="red")+
  geom_point(data= fresh.lo.df,aes(x=Year, y=ext_val,shape=extreme), size=2, color="red")+
  # the 3 bottom temperature events
  geom_vline(xintercept = c(2006,2012,2016), color="grey60", linetype=2)+
  # Esthetics
  labs(x="", y="")+   
  geom_text(data=ann_text_heat, aes(x=Year, y=value,label=lab),hjust = -0.1, size=3)+
  xlim(1970,2020)+
  scale_shape_manual(values = c(NA,16), guide=F)+
  theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(colour = "grey60",linetype = 3),
        panel.border = element_rect(colour = "black"),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        legend.title = element_text(size=9),
        axis.title.x = element_blank(),axis.text.x = element_blank(),
        axis.title.y = element_blank())

# Figure 2d : Bottom temperature time series, and ECE detection ---------
env_fig <- env_longer%>%
  dplyr::select(contains(c("Year","BS_T.bottom"))) %>% 
  pivot_longer(cols = -Year) %>% 
  filter(!name %in% "BS_T.50m")
ann_text_temp <- cbind.data.frame(Year=c(-Inf), value=c(Inf), 
                                  name = "BS_T.bottom",
                                  lab=c("\n(d) Bottom temperature (°C)"))

temp_BS <- ggplot(env_fig ,aes(x=Year, y=value, linetype=name))+ 
  # study period of the paper
  geom_rect(aes(xmin=2004, xmax=2017, ymin=-Inf, ymax=Inf), fill="grey80", alpha=0.025)+
  # bottom temperature data
  geom_line()+theme_pander(base_size=9)+
  # the 3 bottom temperature events
  geom_vline(xintercept = c(2006,2012,2016), color="grey60", linetype=2)+
  # Esthetics
  geom_text(data=ann_text_temp,
            aes(label=lab),hjust = -0.1, size=3)+
  labs(x="", y="Bottom temperature (°C)", linetype="")+ 
  scale_linetype_manual(values=c(1,2),guide=F)+
  xlim(1970,2020)+
  theme(panel.grid.minor.y = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_line(colour = "grey60",linetype = 3),
        panel.border = element_rect(colour = "black"),
        legend.position = c(0.5, 0.525), legend.direction = "horizontal",
        legend.key = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(size=9),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Figure 2e-h : Species responses ----------------------------------------------

# A. Individual species responses ----------------------------------------------
# 1. scale response variables (already done for densities)
# All species will have the same weight in the plot
geog_extent_std <- geog_extent%>% group_by(species) %>% 
  mutate(std_ncell = scale(ncell))
spdat_mode_lat_std <- spdat_mode_lat %>% group_by( species)%>%
  mutate(std_lat_mode = scale(mode_yr_lat), # mode
         std_lat_gcenter = scale(gcenter_lat)) # gravity center
spdat_mode_lon_std <- spdat_mode_lon %>% group_by( species)%>%
  mutate(std_lon_mode = scale(mode_yr_lon),  # mode
         std_lon_gcenter = scale(gcenter_lon)) # gravity center

# 2. select standardized columns and rename
lat_sp <- spdat_mode_lat_std[, c("Year", "species", "std_lat_mode")]
lon_sp <- spdat_mode_lon_std[, c("Year", "species", "std_lon_mode")]
geo_ext_sp <- geog_extent_std[, c("Year", "species", "std_ncell")]
dens_sp <- main_ts_sp_ct[, c("Year", "species", "std_catch")]
names(lat_sp)[3] <- names(lon_sp)[3] <- 
  names(geo_ext_sp)[3] <- names(dens_sp)[3] <- "variable"

# 3. create variable for plot title and facetting
lat_sp$plot <- "(h) Mode of latitude"
lon_sp$plot <- "(g) Mode of longitude"
geo_ext_sp$plot <- "(f) Geographical extent" 
dens_sp$plot <- "(e) Mean density"

# 4. combine species response and order factor levels
sp_ts <- rbind.data.frame(dens_sp, geo_ext_sp,lon_sp, lat_sp)
sp_ts$plot <- factor(sp_ts$plot, levels= c("(e) Mean density",
                                           "(f) Geographical extent",
                                           "(g) Mode of longitude",
                                           "(h) Mode of latitude"))

# B. Mean species responses ----------------------------------------------------
# 1. Summarize mean per year.
mean_geog_extent <- geog_extent_std %>% group_by(Year)%>%
  summarise(mean = mean(std_ncell, na.rm=T))
mean_mode_lat <- spdat_mode_lat_std %>% group_by(Year)%>%
  summarise(mean_lat_mode = mean(std_lat_mode, na.rm=T),
            mean_lat_gcenter= mean(std_lat_gcenter, na.rm=T))
mean_mode_lon <- spdat_mode_lon_std %>% group_by(Year)%>%
  summarise(mean_lon_mode = mean(std_lon_mode, na.rm=T),
            mean_lon_gcenter= mean(std_lon_gcenter, na.rm=T))
main_ts_ct = main_ts_sp_ct %>%  group_by(Year) %>%
  summarise(meancatch=mean(std_catch, na.rm=T),
            confintlow= quantile(std_catch, 0.025, na.rm=T),
            confinthigh= quantile(std_catch, 0.975, na.rm=T),
            conf25 = quantile(std_catch, 0.25, na.rm=T),
            conf75 = quantile(std_catch, 0.75, na.rm=T))

# 2. Select only mode of long lat, nod gravity center
lat_mean <- mean_mode_lat[, c("Year","mean_lat_mode")]
lon_mean <- mean_mode_lon[, c("Year","mean_lon_mode")]
geo_ext_mean <- mean_geog_extent[, c("Year","mean")] # rename, just for consistency
dens_mean <- main_ts_ct[, c("Year","meancatch")]
names(lat_mean)[2] <- names(lon_mean)[2] <- 
  names(geo_ext_mean)[2] <- names(dens_mean)[2] <- "variable"

# 3. create variable for plot title and facetting
lat_mean$plot <- "(h) Mode of latitude"
lon_mean$plot <- "(g) Mode of longitude"
geo_ext_mean$plot <-"(f) Geographical extent" 
dens_mean$plot <- "(e) Mean density"

# 4. combine species response and order factor levels
mean_ts <- rbind.data.frame(dens_mean, geo_ext_mean,lon_mean, lat_mean)
mean_ts$plot <- factor(mean_ts$plot, 
                       levels= c("(e) Mean density",
                                 "(f) Geographical extent",
                                 "(g) Mode of longitude", 
                                 "(h) Mode of latitude"))

# C. PLOT Figure 2e-h ----------------------------------------------------------
ann_text <- cbind.data.frame(Year=-Inf, variable=Inf, 
                             plot=c("(e) Mean density",
                                    "(f) Geographical extent",
                                    "(g) Mode of longitude", 
                                    "(h) Mode of latitude"),
                             lab=c("(e) Mean density",
                                   "(f) Geographical extent",
                                   "(g) Mode of longitude", 
                                   "(h) Mode of latitude"))
fig1_biol <-ggplot()+
  # individual species lines
  geom_line(data=sp_ts, aes(x=Year, y=variable, group=species), color="grey")+
  # mean responses
  geom_line(data=mean_ts, aes(x=Year, y=variable), color="black")+
  geom_point(data=mean_ts, aes(x=Year, y=variable), color="black")+
  # the 3 events in bottom temperature
  geom_vline(xintercept = c(2006,2012,2016), color="grey60", linetype=2)+
  # Esthetics
  geom_text(data=ann_text, aes(x=Year, y=variable, label=lab),
            hjust = -0.1, vjust = 1.2, size=3)+
  geom_hline(yintercept = 0, color="grey50", linetype=2)+
  theme_pander(base_size = 9)+facet_wrap(plot~., scales="free_y", ncol=1)+
  labs(y="",x="")+
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.text.x= element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.title.x = element_blank(),axis.title.y = element_blank()
  )

# FIGURE 2: COMBINE ------------------------------------------------------------
((ice_BS / heat_facet / temp_BS +
    plot_layout(heights = c(1,2,1))) | fig1_biol) + 
  plot_layout(widths = c(3,1))



## ---------
## Figure 3 -------------------------------------------------------------------#
## ---------

# Figure 3a: clustering analysis -----------------------------------------------
# 1. Selecting uncorrelated traits
corr <- cor(traits[,-c(1:3,5,6,8)])
pval_mat <- cor_pmat(traits[,-c(1:3,5,6,8)])
ggcorrplot(corr, type = "lower",p.mat = pval_mat, insig = "blank",
           lab = TRUE)
sel_traits <- traits
names(sel_traits)[1] <- "species"

# 2. Join traits and potential niche descriptors
niche_trait_data <- niche_data2 %>% inner_join(sel_traits) %>% as.data.frame()

id_categorical <- which(names(niche_trait_data) %in% 
                          c("species","habitat","feeding.mode",
                            "body.shape","fin.shape","spawning.type"))

# 3. Check for correlations (Supplementary figure S6)
corr <- cor(niche_trait_data[, - c(id_categorical)])
pval_mat <- cor_pmat(niche_trait_data[, - c(id_categorical)])
ggcorrplot(corr, p.mat = pval_mat, insig = "blank",
           lab = TRUE, lab_size = 2)

id_keep <- which(names(niche_trait_data) %in% 
                   c("length.max","fecundity","age.max","tl","feeding.mode",
                     "offspring.size","fin.shape","body.shape",
                     "mode_T.bottom","range_T.bottom",
                     # "mode_S.bottom" is removed, too correlated with depth
                     "mode_depth","range_depth",
                     #"mode_chla", "range_chla", low importance and too correlated with surface Temp
                     "range_T.10m"#,
                     #"range_S.10m" too correlated with ice and surface temp.
                   ))

# 4. Short names for the species (to make plot more readable)
list_spe <- paste(substring(niche_trait_data$species, 1, 1),
                  list2DF(sapply(strsplit(niche_trait_data$species," "),
                                 function(x) x[-1])),
                  sep=".")
list_spe[grepl("NA",list_spe)] <- c("Icelus", "Liparidae")
rownames(niche_trait_data) <- list_spe#niche_trait_data$species

# 5. Clustering analysis and result exploration
res <- PCA(niche_trait_data[,id_keep], quali.sup = c(6,8,9), quanti.sup =10:11)

fviz_screeplot(res) # variance explained by the axes
summary(res) 
dimdesc(res)[1:3] # description of dimension by variables
fviz_pca_biplot(res, repel = T, col.var = "grey") # simple biplot
fviz_contrib(res,"var",axes = 1:2)# variables contribution to 2 first dimensions
fviz_pca_ind(res, repel = T, col.ind = "cos2") # individuals representation on the plan

clus_niche_traits <- HCPC(res) # three groups are proposed, click on the line
clus_niche_traits$desc.var # clusters' description by variables
clus_niche_traits$desc.ind # clusters' description by individuals
clus_niche_traits$desc.axes# clusters' description by axes
clus_niche_traits$desc.var$call

# 6. Plot Figure 3a 
xmen.cols <- unname(piratepal(palette = "appletv"))[c(1,3,5)]

ggclass_niche_traits <- fviz_cluster(clus_niche_traits, repel=T,  
                                     show.clust.cent = F, 
                                     labelsize = 8, 
                                     ggtheme = theme_bw())+
  scale_colour_manual(values = xmen.cols, 
                      labels=c("Artic-like","Boreal-like",
                               "Widespread predators"))+
  scale_fill_manual(values = xmen.cols, 
                    labels=c("Artic-like","Boreal-like",
                             "Widespread predators"))+
  scale_shape_manual(values = c(19,17,15), 
                     labels=c("Artic-like","Boreal-like",
                              "Widespread predators"))

ggclass_niche_traits$layers[[3]]$aes_params$fontface <- "italic"
ggclass_niche_traits <- fviz_add(ggclass_niche_traits, 
                                 df=2.7*res$var$coord[,1:2], geom = "arrow",
                                 color = "black", repel=T, labelsize = 3.5)
ggclass_niche_traits <- fviz_add(ggclass_niche_traits, 
                                 df=2.7*res$quanti.sup$coord[,1:2], 
                                 geom = "arrow",
                                 color = "grey40", repel=T, labelsize = 3)
ggclass_niche_traits <-fviz_add(ggclass_niche_traits, 
                                df=res$quali.sup$coord[,1:2], 
                                geom = "point",
                                color = "grey40",repel=T,  labelsize =3)+
  theme(plot.title = element_blank(), legend.position = "bottom")

# Figure 3b,c,d: clusters distribution -----------------------------------------
world.longlat <- map_data("world")
continent.longlat <- world.longlat[world.longlat$region %in% 
                                     c("Norway","Russia", "Sweden","Finland"),]
b = getNOAA.bathy(lon1 = 15, lon2 = 65, lat1 = 68.5, lat2 = 82, 
                  resolution = 10)
bf = fortify.bathy(b)

map_back <- ggplot()+theme_void()+
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-350),
               size=c(0.1),
               colour=c("grey40"))+
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-250),
               size=c(0.1),
               colour=c("grey60"))+
  geom_contour(data = bf, 
               aes(x=x, y=y, z=z),
               breaks=c(-150),
               size=c(0.1),
               colour=c("grey80"))


ggclust1 <- map_back +
  geom_contour(data= weight_dist %>% filter(cluster==1), 
               aes(x=x, y=y, z=z, 
                   group=cluster,color=after_stat(level)), 
               bins = 25)+
  geom_polygon(data=continent.longlat, 
               aes(x=long, y=lat, group=group),
               fill="grey10")+
  coord_cartesian(xlim=c(15,65),ylim=c(68.5,82))+
  theme(legend.position = "none")+
  scale_colour_gradient(high = xmen.cols[1], low="#c7e9c0")
ggclust2 <-  map_back+
  geom_contour(data= weight_dist %>% filter(cluster==2), 
               aes(x=x, y=y, z=z, 
                   group=cluster,color=after_stat(level)), 
               bins = 25)+
  geom_polygon(data=continent.longlat, 
               aes(x=long, y=lat, group=group),
               fill="grey10")+
  coord_cartesian(xlim=c(15,65),ylim=c(68.5,82))+
  theme(legend.position = "none")+
  scale_colour_gradient(high = xmen.cols[2], low="#fdd0a2")
ggclust3 <- map_back+
  geom_contour(data= weight_dist %>% filter(cluster==3), 
               aes(x=x, y=y, z=z, 
                   group=cluster,color=after_stat(level)), 
               bins = 25)+
  geom_polygon(data=continent.longlat, 
               aes(x=long, y=lat, group=group),
               fill="grey10")+
  coord_cartesian(xlim=c(15,65),ylim=c(68.5,82))+
  theme(legend.position = "none")+
  scale_colour_gradient(high = xmen.cols[3], low="#dadaeb")

ggclass_niche_traits/(ggclust1 + ggclust2 +ggclust3)+
  plot_layout(heights = c(5/7,2/7))+
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 10, hjust = 0, vjust = 0))

## ---------
## Figure 4 -------------------------------------------------------------------#
## ---------

# 1. Rename clusters
data_clus_niche_traits = clus_niche_traits$data.clust %>% 
  mutate(species=(niche_trait_data$species))
data_clus_niche_traits$clust <- recode(data_clus_niche_traits$clust,
                                       "1"="Arctic-like",
                                       "2"="Boreal-like",
                                       "3"="Widespread predators")

# 2. Summarise responses among clusters
sp_ct_clus <- main_ts_sp_ct %>% ungroup() %>% 
  inner_join(data_clus_niche_traits[,c("species", "clust")]) %>%
  group_by(Year, clust) %>% 
  mutate(mean_clus = mean(std_catch),
         max_clus= max(std_catch), min_clus=min(std_catch)) %>%
  ungroup()
ncell_clus <- geog_extent %>% group_by(species) %>% 
  mutate(std_ncell=scale(ncell)) %>% 
  inner_join(data_clus_niche_traits[,c("species", "clust")]) %>%
  group_by(Year, clust) %>% 
  mutate(mean_ncell = mean(std_ncell),
         max_clus= max(std_ncell), min_clus=min(std_ncell)) %>% 
  ungroup()
modelon_clus <- spdat_mode_lon %>% group_by(species) %>% 
  mutate(std_lon=scale(mode_yr_lon)) %>% 
  inner_join(data_clus_niche_traits[,c("species", "clust")]) %>%
  group_by(clust, species) %>% 
  mutate(std_lon=scale(residuals(lm(mode_yr_lon ~Year))))%>% 
  group_by(Year, clust) %>% 
  mutate(mean_lon = mean(std_lon),
         max_clus= max(std_lon), min_clus=min(std_lon)) %>% ungroup()
modelat_clus <- spdat_mode_lat %>% group_by(species) %>% 
  mutate(std_lat=scale(mode_yr_lat)) %>% 
  inner_join(data_clus_niche_traits[,c("species", "clust")]) %>%
  group_by(clust, species) %>% 
  mutate(std_lat=scale(residuals(lm(mode_yr_lat ~Year))))%>% 
  group_by(Year, clust) %>% 
  mutate(mean_lat = mean(std_lat),
         max_clus= max(std_lat), min_clus=min(std_lat)) %>% ungroup()

# 3. Plots
ggdensity <- ggplot(data= sp_ct_clus)+ theme_bw()+
  geom_vline(xintercept = c(2006,2012,2016), color="grey60",linetype=2)+
  geom_line(aes(x=Year, y=std_catch, color=clust, group=species), alpha=0.3)+
  geom_line(aes(x=Year, y=mean_clus, color=clust), size=1)+
  geom_point(aes(x=Year, y=mean_clus, color=clust))+
  geom_hline(yintercept = 0,linetype=2)+
  facet_wrap(~clust, scales="free_y", ncol=1)+
  labs(title=" (a) Mean density",color="",fill="", y="scaled responses", x="")+
  scale_colour_manual(values = xmen.cols)+
  guides(colour=guide_legend(title="Year", title.position = "top")) +
  theme( plot.title = element_text(size=9,hjust=0.5),  
         legend.position = "bottom",legend.title.align = 0.5,
         strip.background = element_blank(), 
         strip.text.x = element_blank())

ggncell <- ggplot(data= ncell_clus)+ theme_bw()+
  geom_vline(xintercept = c(2006,2012,2016), color="grey60",linetype=2)+
  geom_line(aes(x=Year, y=std_ncell, color=clust, group=species), alpha=0.3)+
  geom_line(aes(x=Year, y=mean_ncell, color=clust), size=1)+
  geom_point(aes(x=Year, y=mean_ncell, color=clust))+
  geom_hline(yintercept = 0,linetype=2)+
  facet_wrap(~clust, scales="free_y", ncol=1)+
  labs(title=" (b) Geographical\nextent",color="",fill="", y="", x="")+
  scale_colour_manual(values = xmen.cols)+ 
  guides(colour=guide_legend(title="Year", title.position = "top")) +
  theme(axis.title.y = element_blank(), 
        plot.title = element_text(size=9,hjust=0.5),  
        legend.position = "bottom",legend.title.align = 0.5,
        strip.background = element_blank(),
        strip.text.x = element_blank())

gglon <- ggplot(data= modelon_clus)+ theme_bw()+
  geom_vline(xintercept = c(2006,2012,2016), color="grey60",linetype=2)+
  geom_line(aes(x=Year, y=std_lon, color=clust, group=species), alpha=0.3)+
  geom_line(aes(x=Year, y=mean_lon, color=clust), size=1)+
  geom_point(aes(x=Year, y=mean_lon, color=clust))+
  geom_hline(yintercept = 0,linetype=2)+
  facet_wrap(~clust, scales="free_y", ncol=1)+
  labs(title=" (c) Mode of\nlongitude",color="",fill="", y="", x="")+
  scale_colour_manual(values = xmen.cols)+ 
  guides(colour=guide_legend(title="Year", title.position = "top")) +
  theme(axis.title.y = element_blank(), 
        plot.title = element_text(size=9,hjust=0.5),  
        legend.position = "bottom", legend.title.align = 0.5,
        strip.background = element_blank(),
        strip.text.x = element_blank())

gglat <- ggplot(data= modelat_clus)+ theme_bw()+
  geom_vline(xintercept = c(2006,2012,2016), color="grey60",linetype=2)+
  geom_line(aes(x=Year, y=std_lat, color=clust, group=species), alpha=0.3)+
  geom_line(aes(x=Year, y=mean_lat, color=clust), size=1)+
  geom_point(aes(x=Year, y=mean_lat, color=clust))+
  geom_hline(yintercept = 0,linetype=2)+
  facet_wrap(~clust, scales="free_y", ncol=1)+
  labs(title=" (d) Mode of\nlatitude",color="",fill="", y="", x="")+
  scale_colour_manual(values = xmen.cols)+ 
  guides(colour=guide_legend(title="Year", title.position = "top")) +
  theme(axis.title.y = element_blank(), 
        plot.title = element_text(size=9,hjust=0.5), 
        legend.position = "bottom",
        legend.title.align = 0.5,
        strip.background = element_blank(),
        strip.text.x = element_blank())

(ggdensity|ggncell|gglon|gglat) + plot_layout(guides = "collect")&
  theme(legend.position = "bottom",legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

