#### PHANTOM OCEAN BAT ANALYSES ####
#### Author(s): Sarah L Jennings ####
#### Most Recent Update: August 28, 2024 ####
#### Step 4: Distance-Based Redundancy Analysis ####

#### Load packages
library(here)
library(tidyverse)
library(vegan)
library(permute)
# devtools::install_github("gavinsimpson/ggvegan")
library(ggvegan)
library(viridis)
library(colorspace)
library(gt)
library(showtext)

#### Import data stored in rds file ####

# use data where 2016 and ocean sites have been removed. Will allow for appropriate permutation structure
bats_consv.ave <- readRDS(here("Output", "bats_consv_ave.rds")) # also import conservative version

#### Use rank index to identify best measure for distance/dissimilarity matrix ####

# what is the best distance measure? 
head(bats_consv.ave)

# for leq (sound pressure level)
rankindex(bats_consv.ave$leq, bats_consv.ave[14:24], indices = c("jaccard", "horn","gow", "bray", "kulczynski", "hellinger"))

# for median (sound frequency)
rankindex(bats_consv.ave$median, bats_consv.ave[14:24], indices = c("jaccard", "horn","gow", "bray", "kulczynski", "hellinger"))

# gower dissimilarity is the highest for both leq and median. But it is used for mixed numeric and non-numeric data, which does not apply here.
# As bray curtis is also high, we will use that 

##### Create dissimilarity matrices ####
betaBray <-vegdist(bats_consv.ave[14:24], method="bray") 

#### Perform partial distance-based redundancy analysis

# two models: with and without interaction between the two sound variables
# for each model, perform analysis with both standard and conservative approach dissimilarity matrices to see if results differ

# we are using a partial dbRDA to allow for the remove the effect of cluster (spatial grouping) prior to evaluating the effects of the other variables
# we will assess significance of the model using permutations that are performed within each cluster, meaning samples are not moved between clusters during the permutation process

# Model 1: without interaction between sound variables
dbrda.simple <- dbrda(formula = betaBray ~ richness + pc1 + pc2 + year + leq + median + Condition(cluster), # removes the effect of cluster on data prior to evaluating other variables
                         data=bats_consv.ave)
print(dbrda.simple)
vif.cca(dbrda.simple) # look at variance inflation factors for model variables
print(dbrda.simple)
# 30.11% of variation explained by spatial Cluster
# 12.71% of variation explained by predictor variables

dbrda.simple_summary <-summary(dbrda.simple) # save summary as an object
dbrda.simple_summary$concont # constrained axes eigenvalues
# dbRDA1 explains 62.14% of 12.71%
0.6214 * 0.1271 # 7.898% of variation explained by dbRDA1

# which variables are highly correlated with dbRDA1?
dbrda.simple_summary$biplot
# both sound variables. in particular frequency (median) with 0.828

# get p-values and significance for model terms (marginal effects)
set.seed(321)
perm.cluster <- with(bats_consv.ave, permute::how(nperm = 9999, blocks = cluster))
dbrda.simple_anova <- anova(dbrda.simple, by="margin", permutations= perm.cluster)
dbrda.simple_anova # examine anova table for model

# save model and anova table
saveRDS(dbrda.simple, here("Model_Objects", "dbRDA_additive_model.rds"))
saveRDS(dbrda.simple_anova, here("Model_Objects", "dbRDA_additive_ANOVA.rds"))

# read in model and anova table
#dbrda.simple <- readRDS(here("Model_Objects", "dbRDA_additive_model.rds"))
#dbrda.simple_anova <- readRDS(here("Model_Objects", "dbRDA_additive_ANOVA.rds"))


# Model 2: including interaction between sound variables
dbrda.int <- dbrda(formula = betaBray ~ richness + pc1 + pc2 + year + leq * median + Condition(cluster), # removes the effect of cluster on data prior to evaluating other variables
                      data=bats_consv.ave)
print(dbrda.int)
summary(dbrda.int)
summ.int <- summary(dbrda.int)
summ.int$concont
summ.int$biplot
vif.cca(dbrda.int) # look at variance inflation factors for model variables

# get p-values and significance for model terms (marginal effects)
set.seed(456)
perm.cluster <- with(bats_consv.ave, permute::how(nperm = 9999, blocks = cluster))
dbrda.int_anova <- anova(dbrda.int, by="margin", permutations= perm.cluster) 
dbrda.int_anova # examine anova table for model

#### Compare two models with and without interaction of sound variables
set.seed(901)
dbrda.compare_anova <- anova(dbrda.simple, dbrda.int, permutations=how(nperm=9999))
dbrda.compare_anova
# model without interaction is equivalent to model with interaction

# save models and anova table
saveRDS(dbrda.int, here("Model_Objects", "dbRDA_interaction_model.rds"))
saveRDS(dbrda.int_anova, here("Model_Objects", "dbRDA_interaction_ANOVA.rds"))
saveRDS(dbrda.compare_anova, here("Model_Objects", "dbRDA_compare_ANOVA.rds"))

# read in models and anova table
#dbrda.int <- readRDS(here("Model_Objects", "dbRDA_interaction_model.rds"))
#dbrda.int_anova <- readRDS(here("Model_Objects", "dbRDA_interaction_ANOVA.rds"))
#dbrda.compare_anova <- readRDS(here("Model_Objects", "dbRDA_compare_ANOVA.rds"))

#### Make plot for model without interaction ####

# use fortify function in ggvegan package to extract scores from dbrda to data frame, which allows for more plotting options
dbrda_fort <- fortify(dbrda.simple, axes=1:2) # this function extracts the scores of the locations from the dbrda so we can use them with ggplot

# select the environmental variables that go with the bat community data
dbrda_env <- bats_consv.ave %>% dplyr::select(location:pc2)

# separate arrows and scores
# for arrows, keeping only sound variables
arrows <- dbrda_fort %>% filter(score=="biplot") %>% filter(label %in% c("leq", "median")) # take only biplot arrows
scores <- dbrda_fort %>% filter(score=="sites") %>% bind_cols(dbrda_env, .)

# load font to use for table and figure
font_add_google(name="Open Sans", family="opensans")
showtext_auto()

# use the data frame in ggplot
dbrdaplot <- ggplot(scores, aes(x = dbRDA1, y = dbRDA2)) +
  geom_point(aes(color=median), size= 1.5) +
  colorspace::scale_color_continuous_sequential(palette="Viridis", name="Frequency (Hz)", rev=T) + # set colors
  theme_classic() +
  coord_equal()+ # important for properly scaling the plot
  xlim(-2.1, 5) + # give x-axis a little more room to allow for arrows and labels
  ylim(-4,4)

# set locations for arrow labels
label.df <- data.frame(dbRDA1 = c(3.5, 3.9), dbRDA2 = c(0.8, -0.3))

# add arrows with labels for sound variables to the plot
dbrda_ordplot <- dbrdaplot +
  geom_segment(data = arrows,
               aes(x = 0, y = 0, xend = dbRDA1*3.1, yend = dbRDA2*3.1),
               arrow = arrow(angle=25, type="open", length=unit(0.1, "inches")),
               colour = "black") +
  geom_text(data=label.df, aes(x=dbRDA1, y=dbRDA2), label=c("Sound Level", "Frequency"), family="opensans", size=3)+
  theme(legend.title = element_text(size=9, family="opensans"),
  legend.text = element_text(size=7, family="opensans"),
  axis.title.x = element_text(size=10, family="opensans"),
  axis.title.y = element_text(size=10, family="opensans"),
  axis.text.x = element_text(size=9, family="opensans"),
  axis.text.y = element_text(size=9, family="opensans"),
  legend.position = "right",
  legend.title.position = 'top')
dbrda_ordplot

ggsave(here("Figures","Fig5_dbRDA_Ordination.pdf"), dbrda_ordplot, width=12, height=7.5, units="cm")

#### Put together data for table

# create data frame for table
dbrda.plotdf <- data.frame(
  Variable = c("Vegetation Richness", "Vegetation PC1", "Vegetation PC2", "Year", "Sound Level", "Sound Frequency", "Residuals"),
  dbrda.simple_anova$Df,
  dbrda.simple_anova$SumOfSqs,
  dbrda.simple_anova$'F',
  dbrda.simple_anova$`Pr(>F)`)

# edit column names
colnames(dbrda.plotdf) <- c("Variable", "DF", "SS" , "F_Statistic", "p_value")

# set up order for model variables
order_vars <- c("Sound Level", "Sound Frequency", "Vegetation PC1", "Vegetation PC2", 
                "Vegetation Richness", "Year", "Residuals")
# make table
dbrda.table <- dbrda.plotdf %>%
  arrange(factor(Variable, levels=order_vars)) %>% # specify order for Variable column
  gt() %>%
  fmt_number(columns=c(SS, F_Statistic, p_value),decimals=3) %>% # print numbers for these columns with 3 decimals
  sub_missing(columns=everything(), rows=everything(), missing_text = "") %>% # remove NAs for F and p for Residuals
  cols_align('left', columns = Variable) %>% # align Variable column to the left
  cols_label(DF = "df", F_Statistic = "F", p_value = "p") %>% # change printed column names
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              table.font.size= 14, # set font size for remaining aspects of table
              table_body.border.bottom.color = "gray20", # set color for bottom horizontal border of table
              column_labels.border.top.color = "gray20", column_labels.border.bottom.color = "gray20") %>% # set color for upper/lower border for column labels
  tab_style(style = cell_text(align = "center"), # align column labels to center
            locations = cells_column_labels(columns = c(DF, SS, F_Statistic, p_value))) %>%
  tab_style(style = cell_text(align = "right"), # align column text to right
            locations = cells_body(columns = c(DF, SS, F_Statistic, p_value))) %>%
  tab_style(style = cell_text(style="italic"), # make font for column labels italic
            locations = cells_column_labels(columns = c(F_Statistic, p_value))) %>% 
  opt_table_font(google_font(name="Open Sans")) # set font for entire table
 
gtsave(dbrda.table, filename="dbRDA_Results_Table.pdf", path=here("Tables/dbRDA_mvGLM"))

  