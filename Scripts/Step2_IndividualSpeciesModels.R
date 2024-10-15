#### PHANTOM OCEAN BAT ANALYSES ####
#### Author(s): Sarah L Jennings ####
#### Most Recent Update: August 28, 2024 ####
#### Step 2: Individual Species Models ####

#### Load Packages and Import Data ####
library(here)
library(DHARMa)
library(spaMM)
library(tidyverse)
library(ggplot2)
library(performance)
library(gt)
library(showtext)
#devtools::install_github("YuLab-SMU/ggtree")
library(ggtree)
library(ape)
library(geiger)
library(treeio)
library(ggnewscale)
library(colorspace)

# import data
batslong <- readRDS(here("Output", "batslong.rds"))

spp <- as.character(unique(batslong$species)) # make list of species

# load font to use for tables and figures
font_add_google(name="Open Sans", family="opensans")
showtext_auto()

########################################################################################
###### Step 1: Determine Best Random Effects Structure ########
########################################################################################

#### Run Model for Each Species with Simple Random Effects using Restricted Maximum Likelihood Estimation #####

# run models (time consuming - alternatively, skip to step that reads in models saved as RDS)
for(i in 1:length(spp)){
  assign(paste0(spp[i], "_m1_REML"), # give each model a name
         fitme(dets ~ scale(leq)*scale(median) + scale(millum) + year + julian + pc1 + pc2 + richness + (1|cluster/site/point),
               data = batslong[batslong$species==spp[i], ], 
               family = negbin2(), method="REML")
  )
}

# save models as RDS files
for(i in 1:length(spp)){
  saveRDS(get(paste0(spp[i], "_m1_REML")), # get each model so it can be saved
          file=paste0(here("Model_Objects"), "/", spp[i], "_m1_REML", ".rds") # naming scheme for files
  )
}

# read in models that were saved as RDS files
filepath_REML <- list.files(path = here::here("Model_Objects"), pattern = "_m1_REML.rds", full.names = TRUE)
speciesREML <- gsub(pattern = "_m1_REML.rds", replacement = "", x = basename(filepath_REML))

for(i in 1:length(speciesREML)){
  assign(paste0(speciesREML[i],"_m1_REML"), readRDS(filepath_REML[i]))
}


# get AIC all of the models
m1_list <- lapply(filepath_REML, readRDS) # create a list of models
AIC_m1 <- as.data.frame(sapply(m1_list, AIC)) # extract AIC for the list of models and put into a data frame
speciesm1 <- paste0(speciesREML, "_REML_m1") # create list to use as column names
colnames(AIC_m1) <- speciesm1 # assign the species name and model number as the column name
head(AIC_m1)

#### Run Model with Spatio-Temporal Random Effects for Each Species using REML ####

# run models (time consuming - alternatively, skip to step that reads in models saved as RDS)
for(i in 1:length(spp)){
  assign(paste0(spp[i], "_m2_REML"), # give each model a name
         fitme(dets ~ scale(leq)*scale(median) + scale(millum) + year + pc1 + pc2 + richness + 
                 (1|pointID) + Matern(1 | Long + Lat) + AR1(1|julian %in% year),
               data = batslong[batslong$species==spp[i], ], 
               fixed=list(nu=0.5),
               family = negbin2(), method="REML"))
}

# save models as RDS files
for(i in 1:length(spp)){
  saveRDS(get(paste0(spp[i], "_m2_REML")), # get each model so it can be saved
          file= paste0(here("Model_Objects"), "/", spp[i], "_m2_REML", ".rds") # naming scheme for files
  )
}


# read in models that were saved as RDS files
filepath2_REML <- list.files(path = here("Model_Objects"), pattern = "_m2_REML.rds", full.names = TRUE)
species2REML <- gsub(pattern = "_m2_REML.rds", replacement = "", x = basename(filepath2_REML))

for(i in 1:length(species2REML)){
  assign(paste0(species2REML[i],"_m2_REML"), readRDS(filepath2_REML[i]))
}

# get AIC for all of the models
m2_list <- lapply(filepath2_REML, readRDS) # create a list of models
AIC_m2 <- as.data.frame(sapply(m2_list, AIC)) # extract AIC for the list of models and put into a data frame
speciesm2 <- paste0(species2REML, "_REML_m2") # create list to use as column names
colnames(AIC_m2) <- speciesm2 # assign the species name and model number as the column name
head(AIC_m2)


#### Compare AIC from simple and spatio-temporal models #####
AIC_Compare <- bind_cols(AIC_m1, AIC_m2) %>% 
  rownames_to_column(., var="AIC_Type") %>% 
  filter(str_detect(AIC_Type, "conditional AIC")) %>%
  pivot_longer(., cols = Anpa_REML_m1:Tabr_REML_m2, names_to= "Species", values_to= "Conditional_AIC") %>%
  mutate(Model = str_sub(Species, -2), Species = str_sub(Species, end=4)) %>%
  dplyr::select(-AIC_Type) %>%
  arrange(Species, Model)

print(AIC_Compare, n=Inf)
# AIC lower for spatio-temporal model for all species

#### make table of AIC

# get scientific names for species
spp_scinames <- read.csv(here("Data", "PhantomSurfBatTraitsConfirmed.csv"), header=T) %>% 
  mutate(sciname = str_replace(Scientific_Name, "_", " ")) %>%
  rename(Species=Species_Code) %>%
  dplyr::select(Species, sciname) %>%
  arrange(Species)

AICtable <- AIC_Compare %>%
  pivot_wider(., names_from=Model, values_from=Conditional_AIC) %>%
  left_join(., spp_scinames) %>% 
  dplyr::select(sciname, m1, m2) %>%
  gt() %>%
  fmt_number(columns = c(m1, m2), decimals=0) %>% # report numbers in these columns to 3 decimal places
  cols_label(sciname = "Species", m1 = "Simple REs", m2 = "Spatiotemporal REs") %>% # change text for column labels
  tab_spanner(label = "Conditional AIC for Model", columns = c(m1, m2)) %>%
  tab_options(column_labels.font.size = 11, column_labels.font.weight = "bold", # make font of column headings bold
              table.font.size= 10, # set font size for remaining aspects of table
              table_body.border.bottom.color = "gray20", # set color for borders
              column_labels.border.top.color = "gray20", 
              column_labels.border.bottom.color = "gray20") %>%
  tab_style(style = cell_text(style = "italic"), # make species column use italic font
            locations = cells_body(columns = sciname)) %>%
  tab_style(style = cell_text(align = "center"), # set alignment of column labels 
          locations = cells_body(columns = c(m1, m2))) %>%
  opt_table_font(google_font(name="Open Sans"))

AICtable

gtsave(AICtable, "AICtable.pdf", path=here("Tables/Individual Species Models"))

########################################################################################
###### Step 2: Run Spatio-Temporal Models with ML ########
########################################################################################

#### Spatio-temporal models with Maximum Likelihood Estimation ####
### spaMM package creators recommend ML estimation with fitme() function for best estimates of fixed effects

# run models (time consuming - alternatively, skip to step that reads in models saved as RDS)
for(i in 1:length(spp)){
  assign(paste0(spp[i], "_m2_ML"), # give each model a name
         fitme(dets ~ scale(leq) * scale(median) + scale(millum) + year + pc1 + pc2 + richness + 
                 (1|pointID) + Matern(1 | Long + Lat)+ AR1(1|julian %in% year),
               data = batslong[batslong$species==spp[i], ], 
               fixed=list(nu=0.5),
               family = negbin2(), method="ML"))
}

# save models as RDS files
for(i in 1:length(spp)){
  saveRDS(get(paste0(spp[i], "_m2_ML")), # get each model so it can be saved
          file= paste0(here("Model_Objects"), "/", spp[i], "_m2_ML", ".rds") # naming scheme for files
  )
}

# read in models that were saved as RDS files
filepath2_ML <- list.files(path = here("Model_Objects"), pattern = "_m2_ML.rds", full.names = TRUE)
species2ML <- gsub(pattern = "_m2_ML.rds", replacement = "", x = basename(filepath2_ML))

for(i in 1:length(species2ML)){
  assign(paste0(species2ML[i],"_m2_ML"), readRDS(filepath2_ML[i]))
}

# examine outputs
summary(Anpa_m2_ML)
summary(Coto_m2_ML)
summary(Epfu_m2_ML)
summary(Labl_m2_ML)
summary(Laci_m2_ML)
summary(Lano_m2_ML)
summary(Myca_m2_ML)
summary(Mylu_m2_ML)
summary(Myyu_m2_ML)
summary(Pahe_m2_ML)
summary(Tabr_m2_ML)

#### Extract Estimates for Parameters for Each Species from Spatio-temporal Model (ML) ####

species <- as.data.frame(spp) %>% rename(Species=spp) # get species names

estimate <- list()
m2_ML_list <- lapply(filepath2_ML, readRDS) # create a list of models

for(i in 1:length(m2_ML_list)){
estimate[[i]] <- m2_ML_list[[i]]$fixef # pull fixed effects values from each model
}

# organize into data frame
estimatespecies <- bind_rows(estimate) %>% bind_cols(species) %>%
  pivot_longer(!Species, names_to = "Variable", values_to = "Estimate")

# extract standard error for the parameters for each model
SE <- list()

for(i in 1:length(m2_ML_list)){
  SE[[i]] <- sqrt(diag(vcov(m2_ML_list[[i]]))) # get SE for parameters from each model
}

# organize into data frame
SEspecies <- bind_rows(SE) %>% bind_cols(species) %>%
 pivot_longer(!Species, names_to = "Variable", values_to = "Standard_Error")

# combine estimates and SE into a data frame that can be exported as results table
AllModelsFixEff <- inner_join(estimatespecies, SEspecies) %>% 
  rowwise() %>%
  mutate(T_Statistic = Estimate/Standard_Error,
         Upper95CI = Estimate + (1.96*Standard_Error),
         Lower95CI = Estimate -(1.96*Standard_Error),
         Upper =ifelse(Upper95CI > 0, "G", "L"),
         Lower = ifelse(Lower95CI<0, "L", "G"),
         Overlaps_Zero=ifelse(Upper == Lower, "No", "Yes")) %>% 
  dplyr::select(Species, Variable, Estimate, Standard_Error, T_Statistic, Upper95CI, Lower95CI, Overlaps_Zero)

print(AllModelsFixEff)
saveRDS(AllModelsFixEff, here("Output", "AllModelsFixEff.rds"))

########################################################################################
#### Step 3: Estimate % Change in Activity in Response to Acoustic Conditions ########
########################################################################################

# grab coefficients and standard errors for calculations:

# ANPA
summary(Anpa_m2_ML)
# response for leq
Anpa_Beta_leq <- Anpa_m2_ML$fixef[2]
Anpa_SE_leq <- sqrt(diag(vcov(Anpa_m2_ML)))[2]

-(1-exp(Anpa_Beta_leq/sd(batslong$leq))^10) # # This is the difference in bat activity over 10 dB
-(1-exp((Anpa_Beta_leq-1.96*Anpa_SE_leq)/sd(batslong$leq))^10)
-(1-exp((Anpa_Beta_leq+1.96*Anpa_SE_leq)/sd(batslong$leq))^10)

# COTO
summary(Coto_m2_ML)
# response for frequency

Coto_Beta_median <- Coto_m2_ML$fixef[3]
Coto_SE_median <- sqrt(diag(vcov(Coto_m2_ML)))[3]

-(1-exp(Coto_Beta_median/sd(batslong$median))^1000) # difference in bat activity over 1kHz frequency 
-(1-exp((Coto_Beta_median-1.96*Coto_SE_median)/sd(batslong$median))^1000) # Lower 95% bound
-(1-exp((Coto_Beta_median+1.96*Coto_SE_median)/sd(batslong$median))^1000) # Upper 95% bound

# EPFU
summary(Epfu_m2_ML) 
# response for sound level and frequency 

Epfu_Beta_leq <- Epfu_m2_ML$fixef[2]
Epfu_SE_leq <- sqrt(diag(vcov(Epfu_m2_ML)))[2]

-(1-exp(Epfu_Beta_leq/sd(batslong$leq))^10) # # This is the difference in bat activity over 10 dB
-(1-exp((Epfu_Beta_leq-1.96*Epfu_SE_leq)/sd(batslong$leq))^10)
-(1-exp((Epfu_Beta_leq+1.96*Epfu_SE_leq)/sd(batslong$leq))^10)

Epfu_Beta_median <- Epfu_m2_ML$fixef[3]
Epfu_SE_median <- sqrt(diag(vcov(Epfu_m2_ML)))[3]

-(1-exp(Epfu_Beta_median/sd(batslong$median))^1000) # difference in bat activity over 1kHz frequency
-(1-exp((Epfu_Beta_median-1.96*Epfu_SE_median)/sd(batslong$median))^1000) # Lower 95% bound
-(1-exp((Epfu_Beta_median+1.96*Epfu_SE_median)/sd(batslong$median))^1000) # Upper 95% bound

# LABL

summary(Labl_m2_ML) 
# response for sound level  

Labl_Beta_leq <- Labl_m2_ML$fixef[2]
Labl_SE_leq <- sqrt(diag(vcov(Labl_m2_ML)))[2]

-(1-exp(Labl_Beta_leq/sd(batslong$leq))^10) # # This is the difference in bat activity over 10 dB
-(1-exp((Labl_Beta_leq-1.96*Labl_SE_leq)/sd(batslong$leq))^10)
-(1-exp((Labl_Beta_leq+1.96*Labl_SE_leq)/sd(batslong$leq))^10)

# LANO

summary(Lano_m2_ML) 
# response for sound level and frequency 

Lano_Beta_leq <- Lano_m2_ML$fixef[2]
Lano_SE_leq <- sqrt(diag(vcov(Lano_m2_ML)))[2]

-(1-exp(Lano_Beta_leq/sd(batslong$leq))^10) # # This is the difference in bat activity over 10 dB
-(1-exp((Lano_Beta_leq-1.96*Lano_SE_leq)/sd(batslong$leq))^10)
-(1-exp((Lano_Beta_leq+1.96*Lano_SE_leq)/sd(batslong$leq))^10)

Lano_Beta_median <- Lano_m2_ML$fixef[3]
Lano_SE_median <- sqrt(diag(vcov(Lano_m2_ML)))[3]

-(1-exp(Lano_Beta_median/sd(batslong$median))^1000) # difference in bat activity over 1kHz frequency
-(1-exp((Lano_Beta_median-1.96*Lano_SE_median)/sd(batslong$median))^1000) # Lower 95% bound
-(1-exp((Lano_Beta_median+1.96*Lano_SE_median)/sd(batslong$median))^1000) # Upper 95% bound

# MYCA
summary(Myca_m2_ML) 
# response for sound level and frequency 

Myca_Beta_leq <- Myca_m2_ML$fixef[2]
Myca_SE_leq <- sqrt(diag(vcov(Myca_m2_ML)))[2]

-(1-exp(Myca_Beta_leq/sd(batslong$leq))^10) # # This is the difference in bat activity over 10 dB
-(1-exp((Myca_Beta_leq-1.96*Myca_SE_leq)/sd(batslong$leq))^10)
-(1-exp((Myca_Beta_leq+1.96*Myca_SE_leq)/sd(batslong$leq))^10)

Myca_Beta_median <- Myca_m2_ML$fixef[3]
Myca_SE_median <- sqrt(diag(vcov(Myca_m2_ML)))[3]

-(1-exp(Myca_Beta_median/sd(batslong$median))^1000) # difference in bat activity over 1kHz frequency
-(1-exp((Myca_Beta_median-1.96*Myca_SE_median)/sd(batslong$median))^1000) # Lower 95% bound
-(1-exp((Myca_Beta_median+1.96*Myca_SE_median)/sd(batslong$median))^1000) # Upper 95% bound


# MYLU
summary(Mylu_m2_ML) 
# response for sound level  

Mylu_Beta_leq <- Mylu_m2_ML$fixef[2]
Mylu_SE_leq <- sqrt(diag(vcov(Mylu_m2_ML)))[2]

-(1-exp(Mylu_Beta_leq/sd(batslong$leq))^10) # # This is the difference in bat activity over 10 dB
-(1-exp((Mylu_Beta_leq-1.96*Mylu_SE_leq)/sd(batslong$leq))^10)
-(1-exp((Mylu_Beta_leq+1.96*Mylu_SE_leq)/sd(batslong$leq))^10)

# MYYU

summary(Myyu_m2_ML) 
# response for sound level  

Myyu_Beta_leq <- Myyu_m2_ML$fixef[2]
Myyu_SE_leq <- sqrt(diag(vcov(Myyu_m2_ML)))[2]

-(1-exp(Myyu_Beta_leq/sd(batslong$leq))^10) # # This is the difference in bat activity over 10 dB
-(1-exp((Myyu_Beta_leq-1.96*Myyu_SE_leq)/sd(batslong$leq))^10)
-(1-exp((Myyu_Beta_leq+1.96*Myyu_SE_leq)/sd(batslong$leq))^10)

# PAHE
summary(Pahe_m2_ML) 
# response for sound level  

Pahe_Beta_leq <- Pahe_m2_ML$fixef[2]
Pahe_SE_leq <- sqrt(diag(vcov(Pahe_m2_ML)))[2]

-(1-exp(Pahe_Beta_leq/sd(batslong$leq))^10) # # This is the difference in bat activity over 10 dB
-(1-exp((Pahe_Beta_leq-1.96*Pahe_SE_leq)/sd(batslong$leq))^10)
-(1-exp((Pahe_Beta_leq+1.96*Pahe_SE_leq)/sd(batslong$leq))^10)

########################################################################################
#### Step 4: Make Plot to Visualize Model Results ########
########################################################################################

# import phylogenetic tree for bats
battree <- read.nexus(here("Data", "S17613.nex")) # ape package

# trim the tree to include relevant species
battraits <- read.csv(here("Data", "PhantomSurfBatTraitsConfirmed.csv")) # import bat traits that contains scientific names for species
scinames <- battraits %>% column_to_rownames(var="Scientific_Name")

tree_scinames <- geiger::treedata(battree, scinames, sort=T) # match tree with data

# store reduced tree as a new object
battree_trim <- tree_scinames$phy

# plot simplified phylogenetic tree
p <- ggtree(battree_trim) + geom_tiplab(as_ylab=TRUE) # plot names as y-axis to avoid them being truncated when plotting

# create a data frame with different potential species labels
namescodes <- battraits %>% dplyr::select(Scientific_Name, Species_Code) %>% 
  rename(Species=Species_Code) %>%
  mutate(SciName_Labels = str_replace(Scientific_Name, "_", " ")) # create version of scientific name without the underscore for labeling plot

# make plot of phylogenetic tree
treeplot <- ggtree(battree_trim) %<+% 
  namescodes + geom_tiplab(aes(label=SciName_Labels), # use SciName_Labels column from namescode to name tree tips
                           size = 3, # font size for tip labels
                           color="black", # font color for tip labels
                           offset=5, # amount of space between tips and labels
                           family="opensans",
                           fontface="italic"
  ) +
  hexpand(0.9) # specify a fraction of x range to expand the x-axis limit by. direction = 1 is for right side

treeplot
# organize data frame of species effects for plotting           
heatmap_dat <- AllModelsFixEff %>% filter(! Variable %in% c("(Intercept)", "year2017", "year2018")) %>%
  left_join(., namescodes) %>%
  dplyr::select(Scientific_Name, Variable, T_Statistic) %>%
  pivot_wider(., names_from = Variable, values_from=T_Statistic ) %>%
  rename(SoundLevel = 'scale(leq)', SoundFrequency = 'scale(median)', 'Level*Frequency' = 'scale(leq):scale(median)', MoonIllumination = 'scale(millum)', 
         VegetationPC1 = pc1, VegetationPC2 = pc2, PlantRichness = richness) %>%
  column_to_rownames(., var="Scientific_Name")

# plot species effects as heatmap with phylogenetic tree on left

# set up divergent palette
colorspace::diverging_hcl(n = 11, h = c(294, 146), c = c(32, 68), l = c(20, 97), power = c(1, 1.4), 
                          register = "custom_greenpurple" )

treeplusheatmap <- gheatmap(
  treeplot, heatmap_dat, offset=50, width=1.2,
  colnames_angle = 75, 
  colnames_offset_y = -0.2, 
  colnames_position="top",
  font.size=3,
  family="opensans",
  hjust=0,
  custom_column_labels = c("Sound Level", "Sound Frequency", "Moon Illumination", "Vegetation PC1",
                           "Vegetation PC2", "Plant Richness", "Level: Frequency")) +
  scale_fill_continuous_diverging(palette = "custom_greenpurple", n_interp=11, mid = 0, limits=c(-4.25,4.25), name="t-statistic", rev=T)+
  vexpand(.2,1) +
  theme(axis.title.x = element_text(size = 8, face = "bold",family="opensans"),
        legend.title = element_text(size= 8, family="opensans"),
        legend.text = element_text(size= 8, family="opensans"))

treeplusheatmap

# save
ggsave(treeplusheatmap, filename = "Fig3_TreePlusHeatMap.pdf", path = here("Figures"), width=16, height=11, units = "cm")

########################################################################################
#### Step 5: Make Model Results Tables ########
########################################################################################

# change labels for variables to improve ease of interpretation
AllModelsFixEff$Variable <- paste0(rep(c("Intercept", "Sound Level", "Sound Frequency", "Moonlight Illumination", "Year (2017)", "Year (2018)", 
                    "Vegetation PC1", "Vegetation PC2", "Vegetation Richness", "Sound Level:Frequency"),11))

# get scientific names for species
spp_scinames <- read.csv(here("Data", "PhantomSurfBatTraitsConfirmed.csv"), header=T) %>% 
  mutate(sciname = str_replace(Scientific_Name, "_", " ")) %>%
  rename(Species=Species_Code) %>%
  dplyr::select(Species, sciname) %>%
  arrange(Species)

# add sci names as column
FixEff_Spp <- left_join(AllModelsFixEff, spp_scinames)
  
# split model results by species
species_split <- FixEff_Spp %>%
  group_by(Species) %>%
  group_split(.) %>%
  setNames(unique(FixEff_Spp$sciname))

# make list of tables
spp_table_list <- map2(
  species_split,
  names(species_split),
  ~gt(.) %>%
    tab_header(title=.y) %>%
    cols_hide(columns=c(Species, sciname, Overlaps_Zero)) %>% # hide several columns
    fmt_number(columns=c(Estimate, Standard_Error, T_Statistic, Lower95CI, Upper95CI), decimals=3) %>% # print numbers for these columns with 3 decimals
    cols_merge(columns=c(Lower95CI, Upper95CI), pattern="{1}, {2}") %>% # merge upper and lower boundaries of 95% CI into one column
    cols_align('left', columns = Variable) %>% # align Variable column to the left. this aligns column label and all text in the column
    cols_label(Standard_Error = "SE", T_Statistic = "t", Lower95CI = "95% CI") %>% # change text for column labels
    tab_footnote(footnote = "estimates based on scaled and centered predictor variables",
                 locations = cells_column_labels(columns = Estimate)) %>%
    tab_footnote(footnote = "estimate reflects difference in intercept between listed year and 2016",
                 locations = cells_body(columns = Variable, rows = c(5,6))) %>%
    tab_options(column_labels.font.size = 11, column_labels.font.weight = "bold", # make font of column headings bold
                heading.title.font.size = 13, # set heading font size
                table.font.size = 10, # set font size for remaining aspects of table
                heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
                table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
                table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
                table.border.bottom.color = "gray20", # set color for very bottom horizontal line
                column_labels.border.bottom.color = "gray20") %>%
    tab_style(style = cell_text(style="italic"), # change plot titles to be italics for species names
              locations = cells_title(groups= 'title')) %>% 
    tab_style(style = cell_text(align = "center"), # align column labels to center for most columns (excluding Variable)
              locations = cells_column_labels(columns = c(Estimate, Standard_Error, T_Statistic, Lower95CI))) %>%
    tab_style(style = cell_text(align = "center"), # align column text that is not the label to center for most columns (excluding Variable)
              locations = cells_body(columns = c(Estimate, Standard_Error, T_Statistic, Lower95CI))) %>%
    cols_align_decimal(columns = c(Estimate, Standard_Error, T_Statistic)) %>% # align text within columns using decimal
    opt_table_font(google_font(name="Open Sans")) # set font for entire table
)

# get vector of species names
scinames <- spp_scinames[,2]

# save tables as pdf
for(i in 1:length(scinames)){
  gtsave(spp_table_list[[i]], # cycle through each table in list to be saved
         filename = paste0(scinames[i], "_table", ".pdf"), # naming scheme
         path = here("Tables/Individual Species Models") # location to save tables
  )}

# alternatively, save as png
for(i in 1:length(scinames)){
  gtsave(spp_table_list[[i]], # cycle through each table in list to be saved
         filename = paste0(scinames[i], "_table", ".png"), # naming scheme
         path = here("Tables/Individual Species Models") # location to save tables
  )}

########################################################################################
#### Step 6: Examine Model Diagnostics for Spatio-temporal Models using ML ########
########################################################################################

# the process outlined below was repeated for each species using the spatio-temporal model with maximum likelihood estimation (m2_ML)
# just one species (ANPA) shown here as an example
# this process was performed for all species

# Get data for the species
anpa <- batslong %>% filter(species=="Anpa")

# Simulate scaled quantile residuals
anpa.res <- simulateResiduals(fittedModel = Anpa_m2_ML, plot=F, n=1000) # simulate scaled residuals

#### Examine residual plots 

# QQ-Plot
plotQQunif(anpa.res)

# Residuals vs Predicted values
plotResiduals(anpa.res) 

# Plots of residuals vs each predictor
plotResiduals(anpa.res, form = anpa$leq, xlab = "leq") 
plotResiduals(anpa.res, form = anpa$median, xlab = "median") 
plotResiduals(anpa.res, form = anpa$millum, xlab = "millum") 
plotResiduals(anpa.res, form = anpa$pc1, xlab = "pc1") 
plotResiduals(anpa.res, form = anpa$pc2, xlab = "pc2") 
plotResiduals(anpa.res, form = anpa$richness, xlab = "richness") 
plotResiduals(anpa.res, form = anpa$year, xlab = "year")

#### Perform goodness of fit tests on residuals

# quantile regression
testQuantiles(anpa.res)

# bootstrap for outliers
testOutliers(anpa.res, alternative="two.sided", margin="both", type="bootstrap", nBoot=1000, plot=T)

# test over/under dispersion and zero inflation
testDispersion(anpa.res)
testZeroInflation(anpa.res)
testUniformity(anpa.res) 

# check for multicollinearity (using the performance package)
check_collinearity(Anpa_m2_ML)
