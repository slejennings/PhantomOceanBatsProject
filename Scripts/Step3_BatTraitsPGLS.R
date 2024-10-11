#### PHANTOM OCEAN BAT ANALYSES ####
#### Author(s): Sarah L Jennings, Clinton D Francis ####
#### Most Recent Update: August 28, 2024 ####
#### Step 3: Phylogenetic Generalized Least Squares (PGLS) Bat Trait Models ####

#### Response variables include:
# t-statistics from individual species models that describe each species' change in activity in response to a) sound level, b) frequency

##### Predictors variables (aka traits) include:
# PCA of wing loading and aspect ratio = reflects flight speed and style
# Ear height:Forearm length = reflects an index from strictly echolocation to gleaning
# Call bandwidth = reflects clutter tolerance
# Call frequency =	reflects frequency of echolocation call

#### Phylogenetic tree from: Shi, J. J., and D. L. Rabosky. 2015. Speciation dynamics during the global radiation of extant bats. Evolution 69:1528–1545.

#########################################################
#### Step 1: set up  ####
#########################################################

# load some needed libraries
library(ape)
library(geiger)
library(nlme)
library(effects)
library(tidyverse)
library(ggeffects)
library(performance)
library(psych)
library(here)
library(gt)
library(showtext)
library(viridis)
library(colorspace)
library(patchwork)

# load fonts
font_add_google(name="Open Sans", family="opensans")
font_add_google(name="Roboto", family="regular")
showtext_auto()

### load phylogenetic tree
tree <- read.nexus(here("Data", "S17613.nex")) # ape package

plot(tree, type ="fan", cex=0.2) # big tree!

### load trait data
traits <- read.csv(here("Data","PhantomSurfBatTraitsConfirmed.csv"))
head(traits)
str(traits)
colnames(traits)

### load individual species responses to sound
responses<- readRDS(here("Output", "AllModelsFixEff.rds"))
head(responses)

#########################################################
#### Step 2: Pair traits with the responses  ####
#########################################################

# filter species' responses for sound variables, select t-statistics, pivot to wide format
resp_wide <- responses %>% 
  rename(Species_Code = Species) %>% # rename to facilitate join with trait info
  filter(Variable %in% c("scale(leq)", "scale(median)")) %>%
  dplyr::select(Species_Code, Variable, T_Statistic) %>% 
  pivot_wider(names_from = Variable, values_from = c(T_Statistic)) %>%
  rename(t_freq = 'scale(median)', t_leq = 'scale(leq)')

# join traits and responses to sound variables 
batTR <- left_join(resp_wide, traits, by="Species_Code")


#########################################################
#### Step 3: link tree with data ####
#########################################################
batTR <- column_to_rownames(batTR, var="Scientific_Name") # set species names as row names

et <- geiger::treedata(tree, batTR, sort=T) # match tree with data

# store reduced tree as a new object
batTree <- et$phy

# plot the reduced tree
par(mar=c(0,0,0,0))
# fan plot
plot(batTree, type ="fan", cex = 0.75) 
# regular plot
plot(batTree)

# put the trait data and sound responses into a data frame
batDat <- data.frame(et$data) 
class(batDat)

#########################################################
#### Step 4: Make relevant predictor variables ####
#########################################################
colnames(batDat)
str(batDat)

# convert several variables into numeric and create variable for ear height/forearm ratio
batDat <- batDat %>% 
  dplyr::select(t_freq, t_leq, Call_Frequency_kHz:Aspect_ratio) %>%
  mutate(across(where(is.character), as.numeric)) %>%
  mutate(E.F = Ear_height_mm/Forearm_length_mm) %>%
  rownames_to_column(., var="Species")

head(batDat)

# perform a PCA of wing loading and aspect ratio
pca.dat <- batDat %>% column_to_rownames(., var="Species") %>% dplyr::select(Wing_loading_N.m2, Aspect_ratio)
PCFlight <- principal(pca.dat, nfactors =1, rotate ="varimax")
PCFlight ## PC1 explains 91% of variance

# add the PC scores to the larger data frame
PC1 <- as.data.frame(PCFlight$scores) %>% rownames_to_column(., var="Species")
batDat <-inner_join(batDat, PC1)

####################################################################
#### Step 5: Model the responses to frequency for the 4 traits ####
####################################################################

### ear to forearm length                
medianE.F <- gls(t_freq ~ log(E.F), correlation = corPagel(0, phy=batTree,fixed=F, form = ~Species), data = batDat)
summary(medianE.F)
confint(medianE.F)

check_model(medianE.F)

saveRDS(medianE.F, here("Model_Objects", "PGLS_EF_Freq.rds"))

# medianE.F <- readRDS(here("Model_Objects", "PGLS_EF_Freq.rds"))

####################################################################                  
### wing PCA
medianWing <- gls(t_freq ~ PC1, correlation = corPagel(0, phy=batTree,fixed=F, form = ~Species), data = batDat)
summary(medianWing)
confint(medianWing)

check_model(medianWing)

saveRDS(medianWing, here("Model_Objects", "PGLS_WingPC_Freq.rds"))

#medianWing <- readRDS(here("Model_Objects", "PGLS_WingPC_Freq.rds"))

####################################################################
### call bandwidth
medianBand <- gls(t_freq ~ Call_Bandwidth_kHz, correlation = corPagel(0, phy=batTree,fixed=F, form = ~Species), data = batDat)
summary(medianBand)
confint(medianBand)

check_model(medianBand)

saveRDS(medianBand, here("Model_Objects", "PGLS_Bandwidth_Freq.rds"))

#medianBand <- readRDS(here("Model_Objects", "PGLS_Bandwidth_Freq.rds"))

####################################################################
### call freq
medianFreq <- gls(t_freq ~ Call_Frequency_kHz, correlation = corPagel(0, phy=batTree,fixed=F, form = ~Species), data = batDat)
summary(medianFreq)
confint(medianFreq)

check_model(medianFreq)

saveRDS(medianFreq, here("Model_Objects", "PGLS_CallFreq_Freq.rds"))

#medianFreq <- readRDS(here("Model_Objects", "PGLS_CallFreq_Freq.rds"))

######################################################################
#### Step 5: Model the responses to sound level for the 4 traits ####
######################################################################

## ear to forearm length                
leqE.F <- gls(t_leq ~ log(E.F), correlation = corPagel(0, phy=batTree,fixed=F, form = ~Species), data = batDat)
summary(leqE.F)
confint(leqE.F)

check_model(leqE.F)

saveRDS(leqE.F, here("Model_Objects", "PGLS_EF_SoundLevel.rds"))

#leqE.F <- readRDS(here("Model_Objects", "PGLS_EF_SoundLevel.rds"))

####################################################################                  
### wing PCA
leqWing <- gls(t_leq ~ PC1, correlation = corPagel(0, phy=batTree,fixed=F, form = ~Species), data =batDat)
summary(leqWing)
confint(leqWing)

check_model(leqWing) 

saveRDS(leqWing, here("Model_Objects", "PGLS_WingPC_SoundLevel.rds"))

#leqWing <- readRDS(here("Model_Objects", "PGLS_WingPC_SoundLevel.rds"))

####################################################################
### call bandwidth
leqBand <- gls(t_leq ~ Call_Bandwidth_kHz, correlation = corPagel(0, phy=batTree,fixed=F, form = ~Species), data =batDat)
summary(leqBand)
confint(leqBand)

check_model(leqBand)

saveRDS(leqBand, here("Model_Objects", "PGLS_Bandwidth_SoundLevel.rds"))

#leqBand <- readRDS(here("Model_Objects", "PGLS_Bandwidth_SoundLevel.rds"))

####################################################################
### call freq
leqFreq <- gls(t_leq ~ Call_Frequency_kHz, correlation = corPagel(0, phy=batTree,fixed=F, form = ~Species), data =batDat)
summary(leqFreq)
confint(leqFreq)

check_model(leqFreq) 

saveRDS(leqFreq, here("Model_Objects", "PGLS_CallFreq_SoundLevel.rds"))

#leqFreq <- readRDS(here("Model_Objects", "PGLS_CallFreq_SoundLevel.rds"))

######################################################################
#### Step 6: Model response to db by response to frequency ####
######################################################################

######################################################################
#### response to db by response to frequency
leqFreqR <- gls(t_leq ~ t_freq, correlation = corPagel(0, phy=batTree,fixed=T, form = ~Species), data = batDat) # lambda converged <0, so fixed at zero
summary(leqFreqR)
confint(leqFreqR)

check_model(leqFreqR)

saveRDS(leqFreqR, here("Model_Objects", "PGLS_LevelFreqResponse.rds"))

#leqFreqR <- readRDS(here("Model_Objects", "PGLS_LevelFreqResponse.rds"))

######################################################################
#### Step 7: Make Plots ###
#######################################################################

#### for frequency and ear:forearm length
E.Fmed_pred <- predictorEffect("E.F", medianE.F)
E.Fmed_df<-data.frame(E.Fmed_pred)

font_add_google(name="Open Sans", family="opensans")
showtext_auto()

D1 <- ggplot(data=E.Fmed_df, aes(x=log(E.F), y=fit)) +
  geom_line(color ="#00BE7D",lwd=1.2, lty=1) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#00BE7D",alpha=.2, lwd=.1)

D2 <- D1 +
  geom_point(data = batDat,aes(x=log(E.F),y=t_freq),color="#00BE7D", size=2, pch=19)+
  geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6) +
  labs(x="Natural Log(Ear:Forearm)", 
       y = "Effect of frequency on activity",
       subtitle = "(b)") +
  theme_classic() +
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=1), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9, family="opensans"),
        axis.text.y = element_text(color = "black", size = 9, family="opensans"),
        axis.title.x = element_text(size = 9, color = "black",family="opensans"),
        axis.title.y = element_text(size = 9, color = "black", family="opensans"),
        plot.subtitle = element_text(face="bold", size = 11, family="opensans")) 

D2

######################################################################
#### for sound level and ear:forearm length

E.Fleq_pred <- predictorEffect("E.F", leqE.F)
E.Fleq_df <-data.frame(E.Fleq_pred)

E1 <- ggplot(data=E.Fleq_df, aes(x=log(E.F), y=fit)) +
  geom_line(color ="#00BE7D",lwd=1.2, lty=1) +
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#00BE7D",alpha=.2, lwd=.1)

E2 <- E1 +
  geom_point(data = batDat,aes(x=log(E.F),y=t_leq),color="#00BE7D", size=2, pch=19)+
  geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6) +
  labs(x="Natural Log(Ear:Forearm)", 
       y = "Effect of sound level on activity",
       subtitle = "(a)") +
  theme_classic() +
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=1), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9, family="opensans"),
        axis.text.y = element_text(color = "black", size = 9, family="opensans"),
        axis.title.x = element_text(size = 9, color = "black",family="opensans"),
        axis.title.y = element_text(size = 9, color = "black", family="opensans"),
        plot.subtitle = element_text(face="bold", size = 11, family="opensans"))

E2

######################################################################
### for leq and call bandwidth

leqBand_pred <- predictorEffect("Call_Bandwidth_kHz", leqBand)
leqBand_df <- data.frame(leqBand_pred)

F1 <- ggplot(data = leqBand_df, aes(x=Call_Bandwidth_kHz, y=fit))+
  geom_line(color="#00BE7D",lwd=1.2, lty=1)+
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#00BE7D", alpha=.2, lwd=.1)

F2 <- F1 + 
  geom_point(data=batDat, aes(x=Call_Bandwidth_kHz, y=t_leq), color="#00BE7D", size=2, pch=19)+
  geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6) +
  labs(x="Call bandwidth (kHz)", 
       y = "Effect of sound level on activity",
       subtitle = "(d)")+
  theme_classic() +
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9, family="opensans"),
        axis.text.y = element_text(color = "black", size = 9, family="opensans"),
        axis.title.x = element_text(size = 9, color = "black",family="opensans"),
        axis.title.y = element_text(size = 9, color = "black", family="opensans"),
        plot.subtitle = element_text(face="bold", size = 11, family="opensans")) 

F2

######################################################################
#### for response to freq and response to sound level

leqFreqR_pred <- predictorEffect("t_freq", leqFreqR)
leqFreqR_df <-data.frame(leqFreqR_pred)

G1 <- ggplot(data=leqFreqR_df, aes(x=t_freq, y=fit))+
  geom_line(color="#00BE7D", lwd=1.2, lty=1)+
  geom_ribbon(aes(ymin=lower, ymax=upper), fill="#00BE7D",alpha=.2, lwd=.1) 

G2 <- G1 + geom_point(data=batDat,aes(x=t_freq, y=t_leq), color="#00BE7D", size=2, pch=19)+
  geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6) +
  labs(x="Effect of frequency on activity", 
       y = "Effect of sound level on activity",
       subtitle = "(c)") +
  theme_classic() + 
  theme(#panel.border = element_rect(colour = "black", fill=NA, size=1), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 9, family="opensans"),
        axis.text.y = element_text(color = "black", size = 9, family="opensans"),
        axis.title.x = element_text(size = 9, color = "black",family="opensans"),
        axis.title.y = element_text(size = 9, color = "black", family="opensans"),
        plot.subtitle = element_text(face="bold", size = 11, family="opensans")) 

G2

#### save plots to file ####
g <- E2 + D2 + G2 + F2
ggsave(here("Figures","Fig4_PGLS_FourPanelTraits.pdf"), g, width = 12, height = 12, units = "cm")


 ###############################################################################
#### Step 8: Make tables of model results ####
###############################################################################

# results for responses to frequency models
# 4 models - Ear:Forearm, Wing PCA, Call Frequency, Call Bandwidth

# ear:forearm and frequency
medianEF_CI <- data.frame(confint(medianE.F)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

medianEF_res <- data.frame(summary(medianE.F)$tTable) %>% # put the model results into data frame
  rownames_to_column(., var="Variable") %>%
  left_join(., medianEF_CI) # combine with 95% CI

summary(medianE.F)$modelStruct # look at lambda
r2(medianE.F) # look at R2

EF_freq <- medianEF_res %>% gt() %>%
  tab_header("Ear to Forearm Length") %>%
  tab_source_note(html("&#955; = 0.764, R<sup>2</sup> = 0.284")) %>% # put values for lambda and R2 in source note (very bottom of table)
  fmt_number(columns=c(Variable, Value, Std.Error, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Value = "Estimate", Std.Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("log(E.F)" ~ "log(Ear:Forearm)", # change the way the model variables are labeled
                  "(Intercept)"~"Intercept") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # align column labels to center for most columns (excluding Variable)
            locations = cells_column_labels(columns = c(Value, Std.Error, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "right"), # align column text that is not the label to right for most columns (excluding Variable)
            locations = cells_body(columns = c(Value, Std.Error, t.value, lowerCI))) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  opt_table_font(google_font(name="Open Sans")) # set font for entire table


# wing PCA and frequency
medianWing_CI <- data.frame(confint(medianWing)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

medianWing_res <- data.frame(summary(medianWing)$tTable) %>% # put the model results into data frame
  rownames_to_column(., var="Variable") %>%
  left_join(., medianWing_CI)

summary(medianWing)$modelStruct # look at lambda
r2(medianWing) # look at R2

WingPC_freq <- medianWing_res %>% gt() %>%
  tab_header("Wing Morphology (PC1)") %>%
  tab_source_note(html("&#955; = 0.607, R<sup>2</sup> = 0.027")) %>%
  fmt_number(columns=c(Variable, Value, Std.Error, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Value = "Estimate", Std.Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("PC1" ~ "Wing PC1", # change the way the model variables are labeled
                  "(Intercept)"~"Intercept") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # align column labels to center for most columns (excluding Variable)
            locations = cells_column_labels(columns = c(Value, Std.Error, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "right"), # align column text that is not the label to right for most columns (excluding Variable)
            locations = cells_body(columns = c(Value, Std.Error, t.value, lowerCI))) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  opt_table_font(google_font(name="Open Sans")) # set font for entire table


# Call Bandwidth
medianBand_CI <- data.frame(confint(medianBand)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

medianBand_res <- data.frame(summary(medianBand)$tTable) %>% # put the model results into data frame
  rownames_to_column(., var="Variable") %>%
  left_join(., medianBand_CI) # combine with 95% CI
  
summary(medianBand)$modelStruct # look at lambda
r2(medianBand) # look at R2
# need to add these values to the code for the table

BandW_freq <- medianBand_res %>% gt() %>%
  tab_header("Call Bandwidth") %>%
  tab_source_note(html("&#955; = 0.734, R<sup>2</sup> = 0.076")) %>% # put values for lambda and R2 in source note (very bottom of table)
  fmt_number(columns=c(Variable, Value, Std.Error, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Value = "Estimate", Std.Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("Call_Bandwidth_kHz" ~ "Call bandwidth", # change the way the model variables are labeled
                  "(Intercept)"~"Intercept") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # align column labels to center for most columns (excluding Variable)
            locations = cells_column_labels(columns = c(Value, Std.Error, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "right"), # align column text that is not the label to right for most columns (excluding Variable)
            locations = cells_body(columns = c(Value, Std.Error, t.value, lowerCI))) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  opt_table_font(google_font(name="Open Sans")) # set font for entire table


# Call Frequency
medianFreq_CI <- data.frame(confint(medianFreq)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

medianFreq_res <- data.frame(summary(medianFreq)$tTable) %>% # put the model results into data frame
  rownames_to_column(., var="Variable") %>%
  left_join(., medianFreq_CI)

summary(medianFreq)$modelStruct # look at lambda
r2(medianFreq) # look at R2

CFreq_freq <- medianFreq_res %>% gt() %>%
  tab_header("Call Frequency") %>%
  tab_source_note(html("&#955; = 0.237, R<sup>2</sup> = 0.042")) %>% # put values for lambda and R2 in source note (very bottom of table)
  fmt_number(columns=c(Variable, Value, Std.Error, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Value = "Estimate", Std.Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("Call_Frequency_kHz" ~ "Call frequency", # change the way the model variables are labeled
                  "(Intercept)"~"Intercept") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # align column labels to center for most columns (excluding Variable)
            locations = cells_column_labels(columns = c(Value, Std.Error, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "right"), # align column text that is not the label to right for most columns (excluding Variable)
            locations = cells_body(columns = c(Value, Std.Error, t.value, lowerCI))) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  opt_table_font(google_font(name="Open Sans")) # set font for entire table

# put model tables together as a group
freq_group <- gt_group(EF_freq,WingPC_freq,BandW_freq,CFreq_freq) 
freq_group

# can only save entire group at once as html
gtsave(freq_group, filename="PGLS_FreqGroup_Table.html", path=here("Tables/PGLS Trait Analyses"))

# save each individual table as pdf:
gtsave(EF_freq, filename="PGLS_EF_Freq_Table.pdf", path=here("Tables/PGLS Trait Analyses"))
gtsave(WingPC_freq, filename="PGLS_WingPC_Freq_Table.pdf", path=here("Tables/PGLS Trait Analyses"))
gtsave(BandW_freq, filename="PGLS_BandWidth_Freq_Table.pdf", path=here("Tables/PGLS Trait Analyses"))
gtsave(CFreq_freq, filename="PGLS_CallFreq_Freq_Table.pdf", path=here("Tables/PGLS Trait Analyses"))

####### results for responses to sound level (leq) models
# 5 models - Ear:Forearm, Wing PCA, Call Frequency, Call Bandwidth, and additional model for frequency response

# ear:forearm and frequency
leqEF_CI <- data.frame(confint(leqE.F)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

leqEF_res <- data.frame(summary(leqE.F)$tTable) %>% # put the model results into data frame
  rownames_to_column(., var="Variable") %>%
  left_join(., leqEF_CI) # combine with 95% CI

summary(leqE.F)$modelStruct # look at lambda
r2(leqE.F) # look at R2

EF_leq <- leqEF_res %>% gt() %>%
  tab_header("Ear to Forearm Length") %>%
  tab_source_note(html("&#955; = 0.898, R<sup>2</sup> = 0.198")) %>% # put values for lambda and R2 in source note (very bottom of table)
  fmt_number(columns=c(Variable, Value, Std.Error, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Value = "Estimate", Std.Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("log(E.F)" ~ "log(Ear:Forearm)", # change the way the model variables are labeled
                  "(Intercept)"~"Intercept") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # align column labels to center for most columns (excluding Variable)
            locations = cells_column_labels(columns = c(Value, Std.Error, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "right"), # align column text that is not the label to right for most columns (excluding Variable)
            locations = cells_body(columns = c(Value, Std.Error, t.value, lowerCI))) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  opt_table_font(google_font(name="Open Sans")) # set font for entire table


# wing PCA and frequency
leqWing_CI <- data.frame(confint(leqWing)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

leqWing_res <- data.frame(summary(leqWing)$tTable) %>% # put the model results into data frame
  rownames_to_column(., var="Variable") %>%
  left_join(., leqWing_CI)

summary(leqWing)$modelStruct # look at lambda
r2(leqWing) # look at R2

WingPC_leq <- leqWing_res %>% gt() %>%
  tab_header("Wing Morphology (PC1)") %>%
  tab_source_note(html("&#955; = 0.388, R<sup>2</sup> = 0.038")) %>% # put values for lambda and R2 in source note (very bottom of table)
  fmt_number(columns=c(Variable, Value, Std.Error, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Value = "Estimate", Std.Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("PC1" ~ "Wing PC1", # change the way the model variables are labeled
                  "(Intercept)"~"Intercept") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # align column labels to center for most columns (excluding Variable)
            locations = cells_column_labels(columns = c(Value, Std.Error, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "right"), # align column text that is not the label to right for most columns (excluding Variable)
            locations = cells_body(columns = c(Value, Std.Error, t.value, lowerCI))) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  opt_table_font(google_font(name="Open Sans")) # set font for entire table


# Call Bandwidth
leqBand_CI <- data.frame(confint(leqBand)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

leqBand_res <- data.frame(summary(leqBand)$tTable) %>% # put the model results into data frame
  rownames_to_column(., var="Variable") %>%
  left_join(., leqBand_CI) # combine with 95% CI

summary(leqBand)$modelStruct # look at lambda
r2(leqBand) # look at R2
# need to add these values to the code for the table

BandW_leq <- leqBand_res %>% gt() %>%
  tab_header("Call Bandwidth") %>%
  tab_source_note(html("&#955; = 0.623, R<sup>2</sup> = 0.256")) %>% # put values for lambda and R2 in source note (very bottom of table)
  fmt_number(columns=c(Variable, Value, Std.Error, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Value = "Estimate", Std.Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("Call_Bandwidth_kHz" ~ "Call bandwidth", # change the way the model variables are labeled
                  "(Intercept)"~"Intercept") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # align column labels to center for most columns (excluding Variable)
            locations = cells_column_labels(columns = c(Value, Std.Error, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "right"), # align column text that is not the label to right for most columns (excluding Variable)
            locations = cells_body(columns = c(Value, Std.Error, t.value, lowerCI))) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  opt_table_font(google_font(name="Open Sans")) # set font for entire table


# Call Frequency
leqFreq_CI <- data.frame(confint(leqFreq)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

leqFreq_res <- data.frame(summary(leqFreq)$tTable) %>% # put the model results into data frame
  rownames_to_column(., var="Variable") %>%
  left_join(., leqFreq_CI)

summary(leqFreq)$modelStruct # look at lambda
r2(leqFreq) # look at R2

CFreq_leq <- leqFreq_res %>% gt() %>%
  tab_header("Call Frequency") %>%
  tab_source_note(html("&#955; = 0.325, R<sup>2</sup> = 0.022")) %>% # put values for lambda and R2 in source note (very bottom of table)
  fmt_number(columns=c(Variable, Value, Std.Error, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Value = "Estimate", Std.Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("Call_Frequency_kHz" ~ "Call frequency", # change the way the model variables are labeled
                  "(Intercept)"~"Intercept") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # align column labels to center for most columns (excluding Variable)
            locations = cells_column_labels(columns = c(Value, Std.Error, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "right"), # align column text that is not the label to right for most columns (excluding Variable)
            locations = cells_body(columns = c(Value, Std.Error, t.value, lowerCI))) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  opt_table_font(google_font(name="Open Sans")) # set font for entire table

# frequency and sound level responses

leqFreqR_CI <- data.frame(confint(leqFreqR)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

leqFreqR_res <- data.frame(summary(leqFreqR)$tTable) %>% # put the model results into data frame
  rownames_to_column(., var="Variable") %>%
  left_join(., leqFreqR_CI)

summary(leqFreqR) # look at lambda. Fixed at zero for this model
r2(leqFreqR) # look at R2

leqFreqResponse <- leqFreqR_res %>% gt() %>%
  tab_header("Effect of Sound Frequency on Activity") %>%
  tab_source_note(html("&#955; = 0, R<sup>2</sup> = 0.424")) %>% # put values for lambda and R2 in source note (very bottom of table)
  fmt_number(columns=c(Variable, Value, Std.Error, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Value = "Estimate", Std.Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("t_freq" ~ "Response to Freq", # change the way the model variables are labeled
                  "(Intercept)"~"Intercept") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # align column labels to center for most columns (excluding Variable)
            locations = cells_column_labels(columns = c(Value, Std.Error, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "right"), # align column text that is not the label to right for most columns (excluding Variable)
            locations = cells_body(columns = c(Value, Std.Error, t.value, lowerCI))) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  opt_table_font(google_font(name="Open Sans")) # set font for entire table

# put all tables together into a group
leq_group <- gt_group(EF_leq,WingPC_leq,BandW_leq,CFreq_leq) 
leq_group

# can only save entire group at once as html
gtsave(leq_group, filename="PGLS_SoundLevelGroup_Table.html", path=here("Tables/PGLS Trait Analyses"))

# or save each individual table as pdf:
gtsave(EF_leq, filename="PGLS_EF_SoundLevel_Table.pdf", path=here("Tables/PGLS Trait Analyses"))
gtsave(WingPC_leq, filename="PGLS_WingPC_SoundLevel_Table.pdf", path=here("Tables/PGLS Trait Analyses"))
gtsave(BandW_leq, filename="PGLS_BandWidth_SoundLevel_Table.pdf", path=here("Tables/PGLS Trait Analyses"))
gtsave(CFreq_leq, filename="PGLS_CallFreq_SoundLevel_Table.pdf", path=here("Tables/PGLS Trait Analyses"))

# save table for final model comparing responses to leq and median
gtsave(leqFreqResponse, filename="PGLS_LevelFreqResp_Table.pdf", path = here("Tables/PGLS Trait Analyses"))

