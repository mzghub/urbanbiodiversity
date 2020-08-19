# install and load devtools before using install_github

install.packages("devtools")
library("devtools")
install.packages("usethis")
library("devtools")

install_github("timnewbold/predicts-demo/predictsFunctions")
install_github("timnewbold/StatisticalModels")

# Install some missing packages
install.packages("DHARMa")
install.packages("MuMIn")

# Load some packages for manipulating and modelling the statistical and geospatial data.

library(predictsFunctions) # useful functions for dealing with PREDICTS data
library(StatisticalModels) # useful functions for plotting PREDICTS models and testing spatial autocorrelation
library(raster) # for dealing with spatial data
library(dplyr) # for handy data manipulation functions
library(tidyr) # ditto
library(lme4) # for mixed effects models
library(car) # for getting anova tables with significance values
library(DHARMa) # for model criticism plots
library(MuMIn) # for checking explanatory power of mixed effects models

# Set work directory
setwd("/Library/MengRui")

# Read in the diversity data from the PREDICTS database
diversity <- readRDS("/Library/MengRui/database.rds")

# Overwrite existing variables of Measurement and Sampling Efforts into the same variable
diversity <- mutate(diversity, Measurement = Effort_corrected_measurement, Sampling_effort = Rescaled_sampling_effort)

diversity <- MergeSites(diversity, silent = TRUE)

sites <- diversity %>%
  
# add Diversity_metric_is_valid column
mutate(Diversity_metric_is_valid = TRUE) %>%
  
# calculate SiteMetrics  
SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use")) %>%
  
# calculate the total abundance within each study
group_by(SS) %>%
mutate(MaxAbundance = ifelse(Diversity_metric_type == "Abundance", max(Total_abundance), NA)) %>%
  
ungroup() %>%
  
# now calculate the rescaled abundance (abundance divided by the maximum within each study)
mutate(RescaledAbundance = ifelse(Diversity_metric_type == "Abundance", Total_abundance/MaxAbundance, NA))


# Create Account at URS.Earthdata.Nasa.Gov with ZhangMR and Password London2020!!
hpd <- raster("/Library/MengRui/gpw-v4-population-density-2000/gpw-v4-population-density_2000.tif")
point_coords <- dplyr::select(sites, Longitude, Latitude)

# extract the HPD values for these coordinates
hpd_value <- raster::extract(hpd, point_coords)

# take a look at the data
head(hpd_value)
# place these data into the sites dataframe, adding a column called hpd and a log-transformed hpd column
sites <- mutate(sites, hpd = hpd_value, loghpd = log(hpd + 1))


sites <- sites %>%
mutate(
    
# collapse primary forest and non-forest together into primary vegetation as these aren't well distinguished
    
Predominant_land_use = recode_factor(Predominant_land_use, "Primary forest" = "Primary vegetation", "Primary non-forest" = "Primary vegetation"),
    
# indeterminate secondary veg and cannot decide get NA
Predominant_land_use = na_if(Predominant_land_use, "Secondary vegetation (indeterminate age)"), Predominant_land_use = na_if(Predominant_land_use, "Cannot decide"), Use_intensity = na_if(Use_intensity, "Cannot decide"),
    
    
# set reference levels
Predominant_land_use = factor(Predominant_land_use), Predominant_land_use = relevel(Predominant_land_use, ref = "Primary vegetation"), Use_intensity = factor(Use_intensity),  Use_intensity = relevel(Use_intensity, ref = "Minimal use"))


# take a look at the LandUse/Use intensity split
table(sites$Predominant_land_use, sites$Use_intensity)

source("https://highstat.com/Books/Book2/HighstatLibV10.R")
corvif(sites[ , c("Predominant_land_use", "Use_intensity", "loghpd")])

model_data <- drop_na(sites, Total_abundance, Predominant_land_use, Use_intensity, loghpd)

# Calculate the total Abundance
abundance_data <- diversity %>%
  
# pull out just the abundance measures
filter(Diversity_metric_type == "Abundance") %>%
  
# group by SSBS (each unique value corresponds to a unique site)
group_by(SSBS) %>%
  
# now add up all the abundance measurements within each site
mutate(TotalAbundance = sum(Effort_corrected_measurement)) %>%
  
# ungroup
ungroup() %>%

# pull out unique sites
distinct(SSBS, .keep_all = TRUE) %>%
  
# now group by Study ID
group_by(SS) %>%
  
# pull out the maximum abundance for each study
mutate(MaxAbundance = max(TotalAbundance)) %>%
  
# ungroup
ungroup() %>%
  
# now rescale total abundance, so that within each study, abundance varies from 0 to 1.
mutate(RescaledAbundance = TotalAbundance/MaxAbundance)


# To transform the RescaledAbundance.

model_data <- mutate(model_data, logAbundance = log(RescaledAbundance + 1), sqrtAbundance = sqrt(RescaledAbundance))

m1 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + loghpd +
             
             Predominant_land_use:Use_intensity + 
             
             Predominant_land_use:loghpd + 
             
             (1|SS) + (1|SSB), data = model_data)

# Choose the random effects
m2 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + loghpd +
             
             Predominant_land_use:Use_intensity + 
             
             Predominant_land_use:loghpd + 
             
             (1+Predominant_land_use|SS) + (1|SSB), data = model_data)
m3 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + loghpd +
             
             Predominant_land_use:Use_intensity + 
             
             Predominant_land_use:loghpd + (1+Use_intensity|SS) + (1|SSB), data = model_data)



m4 <- lmer(sqrtAbundance ~ Predominant_land_use + Use_intensity + loghpd +
             
             Predominant_land_use:Use_intensity + 
             
             Predominant_land_use:loghpd + (1+loghpd|SS) + (1|SSB), data = model_data)

# compare the models using Akaike's Information Criterion (AIC)

AIC(m1, m2, m3, m4)

# The lowest AIC value is the best model of the bunch. In this case, m2.
## have a look at the significance of the terms
Anova(m3)
# The Predominant_land_use:Use_intensity interaction is significant, 
# Start with the most complicated, least significant effect and try to remove it from the model.


m3.1 <- update(m3,~.- Predominant_land_use:loghpd)

# refitting model(s) with ML (instead of REML)
anova(m3.1, m3)

# Now let’s look at the model estimates of our mimumum adequate model
summary(m3)

PlotGLMERFactor(model = m3,
                
                data = model_data,
                
                responseVar = "sqrt(Rescaled Abundance)",
                
                xtext.srt = 20,
                
                seMultiplier = 1.96,
                
                logLink="n",
                
                catEffects = c("Predominant_land_use"))

# To see how urban population density affects the land-use change
install.packages("StatisticalModels")
library("StatisticalModels")

PlotGLMERContinuous(m1,diversity,effects,byContEffect=NULL,
                    otherContEffects=character(0),
                    otherFactors=character(0),xlab,ylab,
                    byFactor=NULL,zlab=NULL,outDir=NULL,
                    logLink="n",plotRug=FALSE,seMultiplier=1.96,
                    params=list(),add=FALSE,xlim=NULL,ylim=NULL,
                    zlim=NULL,line.cols=NULL,line.types=NULL,
                    plotUncertainty=TRUE,nPanels=1,main=NULL,
                    yDiv=1,transformX=FALSE)


SiteMetrics(extra.cols = c("SSB", "SSBS", "Predominant_land_use")) %>%
help(SiteMetrics)

vegan::diversity(diversity, index = "simpson", MARGIN = 1, base = exp(1))


