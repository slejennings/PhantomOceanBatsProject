#### PHANTOM OCEAN BAT ANALYSES ####
#### Author(s): Sarah L Jennings ####
#### Most Recent Update: August 28, 2024 ####
#### Step 1: Data Preparation and Organization ####


# load packages
library(here)
library(tidyverse)
library(ggplot2)
library(stringr)
library(lunar)
library(lubridate)
library(suncalc)
library(psych)
library(geodist)
library(sf)
library(geosphere)
library(osmdata)

#### Import Data ####
batswide <-read.csv(here::here("Data", "bats_freq_amp_11_24_20.csv"), header=T)
head(batswide)
nrow(batswide) # 5361 rows
str(batswide)

#### Add Moon Variables ####
# recreate calendar dates for extracting moon phase
batswide$dateint <- paste(batswide$year, batswide$julian,sep="/")
batswide$date <-as.Date(batswide$dateint, format="%Y/%j")
batswide <-batswide %>% dplyr::select(-dateint)

# add in moon illumination variable 
batswide$millum <- lunar.illumination(as.Date(batswide$date, format = "%Y-%m-%d"), shift=-8) # gives the proportion of illuminated surface on specified date (%). Continuous variable 

#### Add in Vegetation Variables ####

# import veg data
veg <- read.csv(here::here("Data", "veg_lvl_point_11.11.19.csv"), header =T)
head(veg)
veg <- veg %>% dplyr::select(-X)

# select variables for PCA 
vegvars <- matrix(c(veg$AVG.HT, veg$prop_woodyshrub, veg$prop_forbe, veg$prop_grass, veg$prop_nonplant),ncol=5)

# perform PCA on vegetation variables
veg.pca <- principal(vegvars, nfactors=2, rotate="varimax", scores=TRUE) # using psych package
veg.pca
# 83% of variance explained by first two PCs

# extract scores for PC1 and PC2
veg$pc1 <- veg.pca$scores[,1]
veg$pc2 <- veg.pca$scores[,2]

# remove locations with no vegetation data (PC scores are NA)
veg <- veg %>% filter(!is.na(pc1))

# rename sites for merging
veg$site <- as.character(veg$site)
veg$site[veg$cluster=="1"]<- paste("NB",veg$site[veg$cluster=="1"],sep="")
veg$site[veg$cluster=="2"]<- paste("SB",veg$site[veg$cluster=="2"],sep="")
veg$site[veg$cluster=="3"]<- paste("NB",veg$site[veg$cluster=="3"],sep="")
veg$site[veg$cluster=="4"]<- paste("NB",veg$site[veg$cluster=="4"],sep="")
veg$site[veg$cluster=="5"]<- paste("NB",veg$site[veg$cluster=="5"],sep="")

veg$site <- as.factor(veg$site) # change site back to factor

#### Combine Bat and Vegetation Data ####

# simplify veg to retain only the columns that are needed
veg <- veg %>% dplyr::select(site, point, richness, pc1, pc2)

batswide <- inner_join(batswide, veg, by=c( "site", "point")) %>% # combine bat data with veg pc scores
   filter(!is.na(median))  # remove rows with NAs for frequency metrics
nrow(batswide) # 5305 rows

# examine structure and summary
str(batswide)
summary(batswide)

# convert site, point, cluster, treatment and year into factors
batswide <- batswide %>% mutate_at(c("site", "point", "cluster", "year", "treatment"), as.factor)

#### Remove Bat Species with Only 1-2 Detections Total (Eupe, Myvo, Myth) ####
batswide <- batswide %>% dplyr::select(-Eupe, -Myvo, -Myth)

#### Make Adjustments to Pahe and Labl Detections ####

# labl and pahe had a few nights with a really high number of detections
# all high values were converted to 12 to help with model fit

sort(batswide$Pahe, decreasing=T) # need to change 3 entries
sort(batswide$Labl, decreasing=T) # need to change 6 entries

batswide$Pahe <- ifelse(batswide$Pahe>12, 12, batswide$Pahe)
batswide$Labl <- ifelse(batswide$Labl>12, 12, batswide$Labl)

#### Add Geographic Coordinates for Sampling Locations #####

latlong <- read.csv(here::here("Data", "point_coordinates.csv"), header=T)
head(latlong)

# convert coordinates to sf object
latlong_sf <- latlong %>% st_as_sf(coords=c("Long", "Lat")) %>%
  st_set_crs(4326) # WSG 84

# get bounding box for osm data download (US) and 
# download coastline data for this area
osm_box <- getbb (place_name = "California") %>% # set bounding box
  opq () %>% # convert to overpass query object
  add_osm_feature("natural", "coastline") %>% # define the features of interest
  osmdata_sf() 
# if this gives an error message, restart R and try again

# use dist2Line from geosphere package - only works for WGS84 
# this function will find the distance in meters between the bat locations and coastline
dist_to_coast <- geosphere::dist2Line(p = st_coordinates(latlong_sf), # coordinates for the points
                                      line = st_coordinates(osm_box$osm_lines)[,1:2]) # matrix that contains longitude and latitude of the line

# combine initial location list with distance to coastline measurements
class(latlong) # data frame
class(dist_to_coast) # matrix
dist_to_coast <- as.data.frame(dist_to_coast) # convert to data frame
 
latlong_dist  <- bind_cols(latlong, dist_to_coast) %>% 
  dplyr::select(-lon, -lat)  %>% rename(dist_to_coast = distance)
head(latlong_dist)

# add to the bat data
batswide <- inner_join(batswide, latlong_dist)
head(batswide)

#### Add Unique Identifier for Sampling Locations ####
batswide <- batswide %>% mutate(pointID = as.factor(paste(batswide$site, batswide$point, sep="_")))

#### Convert to Long Format ####
batslong <- batswide %>% 
  pivot_longer(cols = Anpa:Tabr, names_to = "species", values_to = "dets")

head(batslong)

# batslong will be used for the individual species models
# export as RDS
saveRDS(batslong, here::here("Output", "batslong.rds"))

##################### Continue with Additional Organization Steps for Community Analyses ############################

# Each location was sampled for a variable number of nights in each year
# This is potentially problematic for certain community analyses as sampling effort can impact the number of species and individuals detected
# We also have clusters with that contain S, C and T treatment, and other clusters that contain only O (ocean) treatment -
# Need to be careful with this structure for permutation based community analyses

# We will construct 3 datasets that summarize the bat community in each location in each year and use them for the following analyses
# 1) for mvGLM
# 2) for dbRDA
# 3) for RAC and species-treatment associations

# Make a list of locations that have 10 or more nights of detections 
keep_ocean <- batswide %>% mutate(location = paste(site,"_",point,"_",year)) %>% 
  dplyr::count(cluster,site,point,location,year,treatment) %>% filter(n>=10)
nrow(keep_ocean) # 224 locations remain

# use locations in the object "keep_ocean" to reduce batswide
batsdets_w2016 <- batswide %>% 
  mutate(location = paste(site,"_",point,"_",year)) %>% # add location column
  inner_join(., keep_ocean) # this will remove locations with less than 10 detection nights 
head(batsdets_w2016)

# Now, reduce the data to obtain one row per location/year that gives the total number of each bat species detected 
# the resulting data frame will be appropriate for multivariate analyses that do not rely on distance-matrices and/or permutations
# The column "n" which contains number of sampling nights will be used as model offset to account for variable sampling effort
bats.w2016 <- batsdets_w2016 %>% 
  group_by(location) %>% 
  mutate(leq=mean(leq),  median=mean(median),  # mean values for sound variables
         richness=mean(richness), pc1=mean(pc1), pc2=mean(pc2), # mean values for veg
         Anpa=sum(Anpa),Coto=sum(Coto),Epfu=sum(Epfu), Labl=sum(Labl), 
         Laci=sum(Laci), Lano=sum(Lano), Myca=sum(Myca), Mylu=sum(Mylu), 
         Myyu=sum(Myyu), Pahe=sum(Pahe), Tabr=sum(Tabr)) %>% # sum to get the total number of detections for each species
  dplyr::select(location, cluster:treatment, n, leq, median, richness, pc1, pc2, dist_to_coast, Anpa:Tabr)%>%
  distinct() %>% # keep only one row per location/year
  rowwise() %>%
  dplyr::mutate(bat_abundance = sum(c_across(Anpa:Tabr)), # add column for bat abundance
    bat_richness = sum(c_across(Anpa:Tabr)!=0)) %>% # add column for bat species richness
  filter(bat_richness!=0) # remove rows where no bats were detected

nrow(bats.w2016) # 223
head(bats.w2016)
length(unique(bats.w2016$site))
# to be used for mvGLM

#### We need a more conservative approach for other community analyses ####

# averaging the detections for each bat species at a particular location can help account for number of individuals being higher at locations with more sampling effort
# BUT... it does not account for the fact that more species are likely detected with more effort

# look at plot of bat species richness for samples with differing numbers of detection nights (n)
ggplot(data=bats.w2016, aes(x=n, y=bat_richness, color=leq))+
  geom_point()+
  theme_bw()
# we can see that higher richness is observed in sites with more sampling effort
cor(bats.w2016$n, bats.w2016$bat_richness) # 0.571

# as a more conservative approach, we will:
# keep only locations with 15 days or more
# if a location has more than 30 days of effort, we will randomly select 30 nights of data to keep and drop excess sampling dates

# create data frame of locations with 15 to 30 nights of sampling
locations_15to30 <- batswide %>% 
  mutate(location = paste(site,"_",point,"_",year)) %>% 
  count(cluster,site,point,location,year,treatment) %>% 
  filter(n>=15 & n<=30)
nrow(locations_15to30) # 136

# create a data frame of locations with more than 30 nights of sampling
locations_morethan30 <- batswide %>% mutate(location = paste(site,"_",point,"_",year)) %>% 
  count(cluster,site,point,location,year,treatment) %>%
  filter(n>30) 
nrow(locations_morethan30) # 39

# use the locations in object "locations_15to30" to reduce batswide
bats.c1 <- batswide %>% 
  mutate(location = paste(site,"_",point,"_",year)) %>% # add location column
  inner_join(., locations_15to30) # limit the data to locations with 15 or more nights
head(bats.c1)

# use the locations in object "locations_morethan30" to reduce batswide
# randomly select 30 nights of data from each location using slice_sample()
set.seed(5678)
bats.c2 <- batswide %>%
  mutate(location = paste(site,"_",point,"_",year)) %>% # add location column
  inner_join(., locations_morethan30) %>% # limit the data to locations with 30 or more nights
  group_by(location) %>%
  slice_sample(n=30) %>% # within each location, randomly select 30 dates (rows) to keep
  mutate(n = 30) # make n equal to 30 for all locations

# if the slice_sample worked correctly, there should be 30 * 39 rows in bats.c2
length(unique(bats.c2$location))*30 == nrow(bats.c2) # this should be TRUE

# combine bats.c1 and bats.c2 and make 2 data frames:
# 1) With ocean sites and 2016 data. Bat abundance for each species averaged by number of sampling days -> for RAC
# 3) Without ocean sites and only 2017/2018 data. Bat abundance for each species averaged by number of sampling days -> for dbRDA

# df1: contains ocean sites and 2016 data
# averaged bat detections
# will be used for Rank Abundance Curve and Species-Treatment Association Analyses
bats_consv.w2016.ave <- bind_rows(bats.c1, bats.c2) %>%
  group_by(location) %>% 
  mutate(leq=mean(leq),  median=mean(median),  # mean values for sound variables
         richness=mean(richness),pc1=mean(pc1),pc2=mean(pc2), # mean values for veg
         Anpa=(sum(Anpa)/n),Coto=(sum(Coto)/n),Epfu=(sum(Epfu)/n), Labl=(sum(Labl)/n), 
         Laci=(sum(Laci)/n), Lano=(sum(Lano)/n), Myca=(sum(Myca)/n), Mylu=(sum(Mylu)/n), 
         Myyu=(sum(Myyu)/n), Pahe=(sum(Pahe)/n), Tabr=(sum(Tabr)/n)) %>%
  dplyr::select(location, cluster:treatment, n, leq, median, richness, pc1, pc2, dist_to_coast, Anpa:Tabr) %>%
  distinct() %>%
  rowwise() %>%
  mutate(bat_abundance = sum(c_across(Anpa:Tabr)), # add column for bat abundance
         bat_richness = sum(c_across(Anpa:Tabr)!=0)) %>% # add column for bat species richness
  filter(bat_richness!=0) # remove rows where no bats were detected

# df2: remove ocean sites and 2016 data (only contains 2017, 2018 and treatments C, S, and P)
# calculate averaged bat detections
# to be used for dbRDA 
bats_consv.ave <- bind_rows(bats.c1, bats.c2) %>%
  filter(year != "2016") %>% filter(treatment != "O") %>%
  group_by(location) %>% 
  mutate(leq=mean(leq),  median=mean(median),  # mean values for sound variables
         richness=mean(richness),pc1=mean(pc1),pc2=mean(pc2), # mean values for veg
         Anpa=(sum(Anpa)/n),Coto=(sum(Coto)/n),Epfu=(sum(Epfu)/n), Labl=(sum(Labl)/n), 
         Laci=(sum(Laci)/n), Lano=(sum(Lano)/n), Myca=(sum(Myca)/n), Mylu=(sum(Mylu)/n), 
         Myyu=(sum(Myyu)/n), Pahe=(sum(Pahe)/n), Tabr=(sum(Tabr)/n)) %>%
  dplyr::select(location, cluster:treatment, n, leq, median, richness, pc1, pc2, dist_to_coast, Anpa:Tabr) %>%
  distinct() %>%
  rowwise() %>%
  mutate(bat_abundance = sum(c_across(Anpa:Tabr)), # add column for bat abundance
         bat_richness = sum(c_across(Anpa:Tabr)!=0)) %>% # add column for bat species richness
  filter(bat_richness!=0) # remove rows where no bats were detected


# Examine new plot of bat species richness vs detection nights

# gather data to do this comparison. We need all 3 years, non-averaged (aka total) bat detections
bats_consv.w2016 <- bind_rows(bats.c1, bats.c2) %>%
  group_by(location) %>% 
  mutate(leq=mean(leq),  median=mean(median),  # mean values for sound variables
         richness=mean(richness),pc1=mean(pc1),pc2=mean(pc2), # mean values for veg
         Anpa=sum(Anpa),Coto=sum(Coto),Epfu=sum(Epfu), Labl=sum(Labl), 
         Laci=sum(Laci), Lano=sum(Lano), Myca=sum(Myca), Mylu=sum(Mylu), 
         Myyu=sum(Myyu), Pahe=sum(Pahe), Tabr=sum(Tabr)) %>% # sum of detections for each species
  dplyr::select(location, cluster:treatment, n, leq, median, richness, pc1, pc2, dist_to_coast, Anpa:Tabr) %>%
  distinct() %>%
  rowwise() %>%
  mutate(bat_abundance = sum(c_across(Anpa:Tabr)), # add column for bat abundance
         bat_richness = sum(c_across(Anpa:Tabr)!=0)) %>% # add column for bat species richness
  filter(bat_richness!=0) # remove rows where no bats were detected

# old plot
ggplot(data=bats.w2016, aes(x=n, y=bat_richness, color=leq))+
  geom_point()+
  theme_bw()

# new plot using more conservative approach
ggplot(data=bats_consv.w2016, aes(x=n, y=bat_richness, color=leq))+
  geom_jitter()+
  theme_bw()
# lower richness values still more common in locations with fewer detection dates (lower n), but much less extreme trend
# use the three data frames above to perform the same analyses as the less conservative approach to see if it changes the results
cor(bats.w2016$n, bats.w2016$bat_richness) # Original correlation = 0.571
cor(bats_consv.w2016$n, bats_consv.w2016$bat_richness) # new correlation = 0.342. This is a reduction

#### Export data frames as rds ####
saveRDS(bats.w2016, here::here("Output", "bats_w2016.rds"))
saveRDS(bats_consv.ave, here::here("Output", "bats_consv_ave.rds"))
saveRDS(bats_consv.w2016.ave, here::here("Output", "bats_consv_w2016_ave.rds"))
