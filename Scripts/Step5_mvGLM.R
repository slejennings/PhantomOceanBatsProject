#### PHANTOM OCEAN BAT ANALYSES ####
#### Author(s): Sarah L Jennings ####
#### Most Recent Update: August 28, 2024 ####
#### Step 5: Multivariate GLM ####

# we will use the manyglm function which fits a separate GLM to each species, using a common set of explanatory variables
# anova and summary on a manyglm object uses resampling-based hypothesis testing to make community-level and taxon-specific inferences about which env variables are associated w/ the multivariate abundances
# this approach takes into account the correlation between species (which standard glms do not)

#### Load packages
library(mvabund)
library(here)
library(tidyverse)
library(gt)
library(showtext)

#### Import data stored in rds files
# use data frame with 2016 and ocean treatment with total abundance of each species (not averaged)
bats.w2016 <- readRDS(here("Output", "bats_w2016.rds"))

#### Set up data and examine
bat_spp <- mvabund::mvabund(bats.w2016[,14:24])

# look at spread of data
boxplot(bat_spp, horizontal=T, las=2, main="abundance") # defaults to log transformation of abundances
# we can see that some species are more abundant and variable than others (as expected)

# plot the mean-variance relationship directly
meanvar.plot(bat_spp) # this has log transformation of axes
meanvar.plot(bat_spp, log="") # same plot without the log transformation
# high means (x-axis) have high variances (y-axis)
# we need to use negative binomial distribution in model


#### Set up the model 
# scale most of the continuous variables
# add n (number of detection nights) as an offset to account for varying sampling effort
# use negative binomial distribution to as we have count data with high variance relative to mean

detnight <- bats.w2016$n # to be used as an offset in the model

# run model and anova table for additive model
# alternatively, use readRDS below to import model objects
mvbat_mod1 <- manyglm(bat_spp ~ scale(median) + scale(leq) + year + pc1 + pc2 + scale(richness) + scale(dist_to_coast), 
                      data = bats.w2016, family = "negative_binomial", offset=log(detnight), cor.type="R")

anova.mod1 <-anova(mvbat_mod1, 
                        test="wald", # Wald test statistic
                        cor.type="shrink", # correlation matrix to account for relationships between species
                        nBoot=9999) # number of resamples 
# cor.type = "shrink" uses ridge regularization to shrink a matrix of correlations between each pair of species towards independence among species 
# two test statistics can be used with this type of correlation structure -> Wald and Rao's score

# save model and anova table as RDS
saveRDS(mvbat_mod1, here("Model_Objects", "mvGLM_additivemodel.rds"))
saveRDS(anova.mod1, here("Model_Objects", "mvGLM_additivemodel.rds"))

# read in model and anova table
#mvbat_mod1 <- readRDS(here("Model_Objects", "mvGLM_additivemodel.rds"))
#anova.mod1 <- readRDS(here("Model_Objects", "mvGLM_additivemodel.rds"))

# run second model that includes the interaction of acoustic variables
# alternatively, use readRDS below to import model objects
mvbat_mod2 <- manyglm(bat_spp ~ scale(median) * scale(leq) + year + pc1 + pc2 + scale(richness) + scale(dist_to_coast), 
                      data = bats.w2016, 
                      family = "negative_binomial", 
                      offset=log(detnight))

anova.mod2 <- anova(mvbat_mod2, 
                    test="wald", 
                    cor.type = "shrink",
                    nBoot=9999)
anova.mod2

# save model and anova table as RDS
saveRDS(mvbat_mod2, here("Model_Objects", "mvGLM_interactionmodel.rds"))
saveRDS(anova.mod2, here("Model_Objects", "mvGLM_interactionmodel.rds"))

# read in model and anova table
#mvbat_mod2 <- readRDS(here("Model_Objects", "mvGLM_interactionemodel.rds"))
#anova.mod2 <- readRDS(here("Model_Objects", "mvGLM_interactionmodel.rds"))

# compare models with and without interaction
compare_1and2 <- anova(mvbat_mod1, mvbat_mod2, nBoot=9999, test="LR")
compare_1and2

# save anova table as RDS
saveRDS(compare_1and2, here("Model_Objects", "mvGLM_comparemodels_ANOVA.rds"))
# read in model saved as RDS
#compare_1and2 <- readRDS(here("Model_Objects", "mvGLM_comparemodels_ANOVA.rds"))

# check model assumptions for model with the interaction
par(mfrow=c(1,3))
plotenvelope(mvbat_mod2, which=1:3, sim.method="stand.norm")
dev.off()

#### export results:
str(anova.mod2) # look at structure of manyglm anova objects
anova.mod2$table # these are the results we want to export

# load font to use for table and figure
library(showtext)
font_add_google(name="Open Sans", family="opensans")

# set up order for model variables
order_mvGLM <- c("Sound Level", "Sound Frequency", "Vegetation PC1", "Vegetation PC2", 
                "Vegetation Richness", "Year", "Distance to Coast", "Frequency: Level")

# organize data frame for table
mvGLM.df <- data.frame(
  anova.mod2$table,
  Variable = c("Intercept", "Sound Frequency", "Sound Level", "Year", "Vegetation PC1",
               "Vegetation PC2", "Vegetation Richness", "Distance to Coast", "Frequency: Level")) %>% 
  filter(Variable != "Intercept") %>%
  dplyr::select(Variable, Df.diff, Res.Df, wald, Pr..wald.) %>%
  rename(df = Df.diff, Residual_df = Res.Df, Wald = wald, p_value = Pr..wald.) %>%
  remove_rownames() %>%
  arrange(factor(Variable, levels = order_mvGLM)) # specify order for Variables column
  
# make table
mvGLM.table <- mvGLM.df %>%
  gt() %>%
  fmt_number(columns=c("Wald", "p_value"),decimals=3) %>% # print numbers for these columns with 3 decimals
  cols_align('left', columns = Variable) %>% # align Variable column to the left
  cols_label(Residual_df = "Residual df", p_value = "p") %>% # change printed column names
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              table.font.size = 14, # set font size for remaining aspects of table
              table_body.border.bottom.color = "gray20", # set color for bottom horizontal border of table
              column_labels.border.top.color = "gray20", column_labels.border.bottom.color = "gray20") %>% # set color for upper/lower border for column labels
  tab_style(style = cell_text(align = "center"), # align column labels to center
            locations = cells_column_labels(columns = c(df, Residual_df, Wald, p_value))) %>%
  tab_style(style = cell_text(align = "center"), # align column text to center
            locations = cells_body(columns = c(df, Residual_df, Wald))) %>%
  tab_style(style = cell_text(align = "right"), # align column text to right
            locations = cells_body(columns = p_value)) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = p_value)) %>% 
  opt_table_font(google_font(name="Open Sans")) %>% # set font for entire table
  sub_small_vals(columns="p_value", threshold=0.001) %>% # change low p-values to be labeled as < 0.001
  cols_width(Residual_df ~ px(105)) # give Residual df column a little extra width

mvGLM.table

gtsave(mvGLM.table, filename="mvGLM_Results_Table.pdf", path=here("Tables/dbRDA_mvGLM"))

