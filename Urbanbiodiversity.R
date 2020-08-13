# Load in the two rds files of PREDICTS, dtabase and sites
diversity <- file.choose("database.rds")
sites <- file.choose("sites.rds")

# Set the work directory 
work_dir <- ("/Library/MengRui")
setwd(work_dir)

# Install the packages needed for this project if missing
install.packages("robustlmm")
install.packages("lme4")

# Load the packages needed for this project
library("dplyr") # for easy data manipulation
library("robustlmm") # for Linear Mixed Effects Models in site analysis
library("lme4") # for mixed effects models
library("car") # for logit transformation with adjustment
library("raster") # for working with raster data
library("foreach") # running loops
library("doParallel") # running loops in parallel
library("tidyr") # To replace duplicated values
library("magrittr") # for piping
library("geosphere")  # calculating geographic distance between sites


# Read in the loaded data of Biodiversity
diversity <- readRDS("database.rds") %>%
  
# filter out the data for the Americas and have a glimpse Of data 
filter(UN_region == "Americas") 
glimpse(diversity)

table(diversity$Predominant_land_use, diversity$Use_intensity, diversity$Habitat_patch_area_square_metres)



# Simplify the data accordantly
diversity <- diversity %>%
  # make a level of Primary minimal. Everything else gets the coarse land use
  mutate(
    LandUse = ifelse(Predominant_land_use == "Primary vegetation" & Use_intensity == "Minimal use",
                     "Primary minimal",
                     paste(Predominant_land_use)),
    
    # collapse the secondary vegetation classes together
    LandUse = ifelse(grepl("secondary", tolower(LandUse)),
                     "Secondary vegetation",
                     paste(LandUse)),
    
    # change cannot decide into NA
    LandUse = ifelse(Predominant_land_use == "Cannot decide",
                     NA, 
                     paste(LandUse)),
    
    # relevel the factor so that Primary minimal is the first level (so that it is the intercept term in models)
    LandUse = factor(LandUse),
    LandUse = relevel(LandUse, ref = "Primary minimal")
  )
 


   

# Calculating the total Abundance
abundance_data <- diversity %>%
  
# Filter out the abundance measures
filter(Diversity_metric_type == "Abundance") %>%
  
# group by SSBS (each unique value corresponds to a unique site)
group_by(SSBS) %>%
  
# add up all the abundance measurements for each site, the value the abundance measurement here is a relative value, implying abundance per unit effort
mutate(TotalAbundance = sum(Effort_corrected_measurement)) %>%
  
# ungroup to avoid potential unintended errors after grouping
ungroup() %>%
  
# pull out unique sites
distinct(SSBS, .keep_all = TRUE) %>%
  
# now group by Study ID
group_by(SS) %>%
  
# pull out the maximum abundance for each study
mutate(MaxAbundance = max(TotalAbundance)) %>%
  
# ungroup as above
ungroup() %>%
  
# rescale the total abundance from 0 to 1 for each study.
mutate(RescaledAbundance = TotalAbundance/MaxAbundance)



# run a simple model
ab_m <- lmer(sqrt(RescaledAbundance) ~ LandUse + (1|SS) + (1|SSB), data = abundance_data)
summary("ab_m")


# Install package vegan for Simpson's diversity calculation
install.packages("vegan")

# Read more about the functions in Package vegan
browseVignettes("vegan")

# Using the newly developed diversity function in vegan to process the Simpson's index, avoid the overlapping names of diversity function and diversity.rds 
data("diversity")

H <- vegan::diversity(diversity)
simp <- diversity(diversity, "simpson")
invsimp <- vegan::diversity(diversity, "inv")
r.2 <- rarefy(diversity, 2)
alpha <- fisher.alpha(diversity)
pairs(cbind(H, simp, invsimp, r.2, alpha), pch="+", col="blue")

# Species richness (S) and Pielou's evenness (J):
S <- specnumber(diversity) 
# rowSums(diversity > 0) does the same...
J <- H/log(S)
