#### PHANTOM OCEAN BAT ANALYSES ####
#### Author(s): Sarah L Jennings ####
#### Most Recent Update: August 28, 2024 ####
#### Step 6: Rank Abundance Curves and Species-Treatment Association Analyses ####

#### Load packages
library(here)
library(tidyverse)
library(codyn)
library(performance)
library(viridis)
library(colorspace)
library(lme4)
library(lmerTest)
library(ggeffects)
library(gt)
library(patchwork)
library(showtext)
library(indicspecies)

#### Import data
# we need data with 2016 and ocean treatments but where bat abundance has been averaged based on the sampling days
# using the more conservative approach to reduce the possibility that differences in richness are due to variation in sampling effort
bats_consv.w2016.ave <- readRDS(here("Output", "bats_consv_w2016_ave.rds"))

# identify the years with data for the various locations
years <- bats_consv.w2016.ave %>% ungroup() %>%
  arrange(site, point, year) %>%
  dplyr::select(site, point, year) %>%
  group_by(site, point) %>%
  mutate(yrs_avail = str_flatten(year, ", ")) %>%
  dplyr::select(site, point, yrs_avail) %>%
  distinct() %>%
  filter(yrs_avail !="2016" & yrs_avail != "2017, 2018" & yrs_avail != "2017" & yrs_avail != "2018") %>%
  mutate(analysis = if_else(yrs_avail=="2016, 2017, 2018", "all",
                            if_else(yrs_avail=="2016, 2017", "2017", "2018"))) 

# reduce to keep only locations with all three years
years_all <- years %>% filter(analysis=="all")
nrow(years_all) # 32 locations with sampling in 3 years

# use list of locations in years_all to reduce bats_consv.w2016.ave
keep <- left_join(years_all, bats_consv.w2016.ave) %>% 
  arrange(site, point, year) %>% 
  mutate(site_point=paste(site, point, sep="_")) %>%
  ungroup() %>%
  dplyr::select(location, site_point, year, treatment, Anpa:Tabr)
head(keep)
# confirm correct number of rows are present
nrow(keep) == nrow(years_all)*3 # this will be TRUE

# make data frame long
keep_long <- keep %>% pivot_longer(cols=Anpa:Tabr, names_to="species", values_to="abund") 
nrow(keep_long) == nrow(keep) * 11 # this should be TRUE if there is now a row for each species

#### Calculate Changes in Rank Abundance Curves #####
# look at differences in community for 2016 vs 2017 (1st year of treatment) and 2016 vs 2018 (2nd year of treatment)

# use codyn package to calculate rank abundance curve changes
# gives five measures of change: richness change, evenness change, rank change, gains and losses
rac <- codyn::RAC_change(df = keep_long,
                           time.var="year",
                           species.var = "species",
                           abundance.var = "abund",
                           replicate.var= "site_point",
                           reference.time = "2016") # compare all with 2016


head(rac)
# we need to add treatment labels based on year2 (aka what treatment did the location have in the experimental years of 2017/2018?)

keep_treat <- keep_long %>%
  dplyr::select(site_point, year, treatment) %>%
  distinct() %>%
  filter(year %in% c("2017", "2018")) %>%
  rename(year2 = year)

# combine treatments with RAC measures of change
RAC_treat <- left_join(rac, keep_treat) %>%
  mutate(treatment = str_replace_all(treatment, c(C="Control", P = "Phantom", S="Shifted") ))
nrow(RAC_treat) # 64 rows because each location is represented twice. Once for change between 2016-2017, and again for change between 2016-2018

#### Run Models #####
# One model for each measure of community change (5 in total)
# use site_point as a random intercept to account for multiple values (years of data) at each location

# make sure year2 is factor
class(RAC_treat$year2)

####  Richness ####

# examine the response variable
range(RAC_treat$richness_change)
hist(RAC_treat$richness_change)

# run model 
mrichc <-lmer(data = RAC_treat, formula = richness_change ~ treatment + year2 + (1|site_point))
summary(mrichc)
anova(mrichc)
confint(mrichc)
r2(mrichc)

# save model as RDS
saveRDS(mrichc, here("Model_Objects", "RAC_Richness.rds"))
# mrichc <- readRDS(here("Model_Objects", "RAC_Richness.rds"))

# load font to use for tables and figures
library(showtext)
font_add_google(name="Open Sans", family="opensans")
theme_set(theme_classic()) # set theme for plots


# make richness plot 
pred.richc <- ggeffects::predict_response(mrichc, c("treatment", "year2"), margin = "marginalmeans")

richness.plot <-plot(pred.richc, show_data=T, jitter=0.15, dot_size = 2, line_size=1.25, colors=c("#4B0055", "#00BE7D")) + 
  theme_classic(base_size = 14, base_family = "") +
  labs(x="", y="Change in Richness", color="Experimental Year", title="",  subtitle = "(a)") +
  theme(
    legend.title = element_text(size = 12, family="opensans"),
    legend.text = element_text(size = 12, family="opensans"),
    axis.title.x = element_text(size = 12,family="opensans", margin = margin(t=10, unit="pt")),
    axis.title.y = element_text(size = 12, family="opensans"),
    axis.text.x = element_text(size=12, family="opensans"),
    axis.text.y = element_text(size=12, family="opensans"),
    plot.subtitle = element_text(size=14,family="opensans", face="bold"),
    legend.position="bottom"
  )
  
# compare treatments
predict_response(mrichc, "treatment", margin = "marginalmeans") # get CIs for treatments

# look at model diagnostics for final model
performance::check_model(mrichc)

#### Evenness ####

# examine response variable
hist(RAC_treat$evenness_change)
range(RAC_treat$evenness_change)

# run model
mevenc <-lmer(data = RAC_treat, formula = evenness_change ~ treatment + year2  + (1|site_point))
summary(mevenc)
anova(mevenc)
r2(mevenc)
confint(mevenc)

# save model as RDS
saveRDS(mevenc, here("Model_Objects", "RAC_Evenness.rds"))
#mevenc <- readRDS(here("Model_Objects", "RAC_Evenness.rds"))

# make plot of predicted effects for additive model 
pred.evenc <- ggeffects::predict_response(mevenc, c("treatment", "year2"), margin = "marginalmeans")
plot(pred.evenc, show_data=T) # make a plot

# compare treatments
predict_response(mevenc, "treatment", margin = "marginalmeans")

# look at model diagnostics for final model
performance::check_model(mevenc)

#### Rank Change (Species Reordering) ####

# examine response variable
hist(RAC_treat$rank_change)

# run the model
mrankc <-lmer(data = RAC_treat, formula = rank_change ~ treatment + year2 + (1|site_point))
summary(mrankc)
anova(mrankc)
r2(mrankc)
confint(mrankc)

# save model as RDS
saveRDS(mrankc, here("Model_Objects", "RAC_RankAbundance.rds"))
#mrankc <- readRDS(here("Model_Objects", "RAC_RankAbundance.rds"))

# make a plot
pred.rankc <- ggeffects::predict_response(mrankc, c("treatment", "year2"), margin = "marginalmeans")
plot(pred.rankc, show_data=T) # make a plot

# look at model diagnostics for final model
performance::check_model(mrankc)

#### Gains ####

# examine response variable
hist(RAC_treat$gains) # skewed
range(RAC_treat$gains) # contains 0, so we can't easily use log tranformation

# run model 
mgainc <-lmer(data = RAC_treat, formula = gains ~ treatment + year2 + (1|site_point))
summary(mgainc)
anova(mgainc)
r2(mgainc)
confint(mgainc)

# save model as RDS
saveRDS(mgainc, here("Model_Objects", "RAC_Gains.rds"))
#mgainc <- readRDS(here("Model_Objects", "RAC_Gains.rds"))

# make a plot
pred.gainc <- ggeffects::predict_response(mgainc, c("treatment", "year2"), margin = "marginalmeans")

gains.plot <-plot(pred.gainc, show_data=T, jitter=0.15, dot_size = 2, line_size=1.25, colors=c("#4B0055", "#00BE7D")) + 
  theme_classic(base_size = 14, base_family = "") +
  labs(x="",y="Species Gains", color="Experimental Year", title="", subtitle="(b)") +
  theme(
    legend.title = element_text(size=12, family="opensans"),
    legend.text = element_text(size=12, family="opensans"),
    axis.title.x = element_text(size = 12, family="opensans", margin = margin(t=10, unit="pt")),
    axis.title.y = element_text(size = 12, family="opensans"),
    axis.text.x = element_text(size=12, family="opensans"),
    axis.text.y = element_text(size=12, family="opensans"),
    plot.subtitle = element_text(size=14, family="opensans", face="bold"),
    legend.position="bottom"
  )

# compare treatments
predict_response(mgainc, "treatment", margin = "marginalmeans") # CIs for treatments. Gives marginal means

# look at model diagnostics for final model
performance::check_model(mgainc)

### Losses ###

# examine response variable
hist(RAC_treat$losses) # slightly skewed
range(RAC_treat$losses) # contains zeros

# run model 
mlossc <-lmer(data = RAC_treat, formula = losses ~ treatment + year2 + (1|site_point))
summary(mlossc)
anova(mlossc)
r2(mlossc)

# save model as RDS
saveRDS(mlossc, here("Model_Objects", "RAC_Losses.rds"))
#mlossc <- readRDS(here("Model_Objects", "RAC_Losses.rds"))

# make a plot
pred.lossc <- ggeffects::predict_response(mlossc, c("treatment", "year2"), marin="marginalmeans")

losses.plot <- plot(pred.lossc, show_data=T, jitter=0.15, dot_size = 2, line_size=1.25, colors=c("#4B0055", "#00BE7D")) + 
  theme_classic(base_size = 14, base_family = "") +
  labs(x="", y="Species Losses", color="Experimental Year", title="", subtitle="(c)") +
  theme(
    legend.title = element_text(size=12, family="opensans"),
    legend.text = element_text(size=12, family="opensans"),
    axis.title.x = element_text(size = 12,family="opensans", margin = margin(t=10, unit="pt")),
    axis.title.y = element_text(size = 12, family="opensans"),
    axis.text.x = element_text(size=12, family="opensans"),
    axis.text.y = element_text(size=12, family="opensans"),
    plot.subtitle = element_text(size=14, family="opensans", face="bold"),
    legend.position="bottom",
    plot.margin = margin(t=0, r=15, b=0, l=0) # add a smidge of white space on right side so the 3 plots look okay when combined
  )
losses.plot
# compare treatments
predict_response(mlossc, "treatment", margin = "marginalmeans") # get CIs for treatments. Calls ggemmeans for marginal means

# look at model diagnostics for final model
performance::check_model(mlossc)

# combine plots with common legend at bottom
RACplots <- richness.plot + gains.plot + losses.plot + 
  plot_layout(guides = "collect")  & # use a single legend 
  theme(legend.position = 'bottom', # put legend at the bottom
        legend.box.spacing=unit(-15, "pt"), # move legend up closer to plots
        legend.margin = margin(1, 1, 1, 1)) # set margin around legend
RACplots

ggsave(here("Figures", "RAC_ThreePanel_Figure.pdf"), width=8.5, height=3.5)

######### Make Tables of Model Results ################

#### Table for Change in Species Richness Model ####
rich_CI <- data.frame(confint(mrichc)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>% slice(., -(1:2)) %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

rich_coef <- data.frame(summary(mrichc)$coefficients) %>% # put the model fixed effects into data frame
  rownames_to_column(., var="Variable") %>%
  rename(p.value = Pr...t..) %>%
  left_join(., rich_CI)

r2(mrichc) # look at conditional and marginal R2
# add values into table below

rich_table <- rich_coef %>% gt() %>%
  tab_header("Change in Species Richness") %>%
  tab_source_note(html("Conditional R<sup>2</sup> = 0.725, Marginal R<sup>2</sup> = 0.220")) %>%
  fmt_number(columns=c(Variable, Estimate, Std..Error, df, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Std..Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("(Intercept)"~"Intercept", # change the way the model variables are labeled
                  "treatmentPhantom" ~ "Treatment (Phantom)",
                  "treatmentShifted" ~ "Treatment (Shifted)",
                  "year22018" ~ "Year (2018)") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_borders(sides = c("top"),  weight = px(1.5), color='gray20'),
            locations = cells_source_notes()) %>%
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # set alignment of column labels for most columns 
            locations = cells_column_labels(columns = c(Estimate, Std..Error, df, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "center"), # set alignment of column text for most columns 
            locations = cells_body(columns = c(Estimate, Std..Error, df, t.value))) %>%
  tab_style(style = cell_text(align = "right"), # set alignments for CI column
            locations = cells_body(columns = lowerCI)) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  tab_footnote(footnote = "estimate reflects the difference in intercept between the specified treatment and control",
               locations = cells_body(columns = Variable, rows = c(2,3))) %>%
  tab_footnote(footnote = "estimate reflects the difference in intercept between 2018 and 2017",
               locations = cells_body(columns = Variable, rows = 4)) %>%
  cols_align_decimal(columns = c(Estimate, Std..Error, df, t.value)) %>% # align text within columns using decimal
  opt_table_font(google_font(name="Open Sans")) # set font for entire table

rich_table

gtsave(rich_table, "RAC_Richness_Table.pdf", path=here("Tables/RAC_SppTreatAssoc"))

#### Table for Change in Evenness Model ####
even_CI <- data.frame(confint(mevenc)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>% slice(., -(1:2)) %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

even_coef <- data.frame(summary(mevenc)$coefficients) %>% # put the model fixed effects into data frame
  rownames_to_column(., var="Variable") %>%
  rename(p.value = Pr...t..) %>%
  left_join(., even_CI)

r2(mevenc) # look at conditional and marginal R2
# add values into table below

even_table <- even_coef %>% gt() %>%
  tab_header("Change in Species Evenness") %>%
  tab_source_note(html("Conditional R<sup>2</sup> = 0.603, Marginal R<sup>2</sup> = 0.120")) %>%
  fmt_number(columns=c(Variable, Estimate, Std..Error, df, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Std..Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("(Intercept)"~"Intercept", # change the way the model variables are labeled
                  "treatmentPhantom" ~ "Treatment (Phantom)",
                  "treatmentShifted" ~ "Treatment (Shifted)",
                  "year22018" ~ "Year (2018)") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>%
  tab_style(style = cell_borders(sides = c("top"),  weight = px(1.5), color='gray20'),
            locations = cells_source_notes()) %>%
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # set alignment of column labels for most columns 
            locations = cells_column_labels(columns = c(Estimate, Std..Error, df, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "center"), # set alignment of column text for most columns 
            locations = cells_body(columns = c(Estimate, Std..Error, df, t.value))) %>%
  tab_style(style = cell_text(align = "right"), # set alignments for CI column
            locations = cells_body(columns = lowerCI)) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  tab_footnote(footnote = "estimate reflects the difference in intercept between the specified treatment and control",
               locations = cells_body(columns = Variable, rows = c(2,3))) %>%
  tab_footnote(footnote = "estimate reflects the difference in intercept between 2018 and 2017",
               locations = cells_body(columns = Variable, rows = 4))%>%
  cols_align_decimal(columns = c(Estimate, Std..Error, df, t.value)) %>% # align text within columns using decimal
  opt_table_font(google_font(name="Open Sans")) # set font for entire table

even_table

gtsave(even_table, "RAC_Evenness_Table.pdf", path=here("Tables/RAC_SppTreatAssoc"))

#### Table for Change in Rank Abundance Model ####
rank_CI <- data.frame(confint(mrankc)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>% slice(., -(1:2)) %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

rank_coef <- data.frame(summary(mrankc)$coefficients) %>% # put the model fixed effects into data frame
  rownames_to_column(., var="Variable") %>%
  rename(p.value = Pr...t..) %>%
  left_join(., rank_CI)

r2(mrankc) # look at conditional and marginal R2
# add values into table below

rank_table <- rank_coef %>% gt() %>%
  tab_header("Change in Rank Abundance") %>%
  tab_source_note(html("Conditional R<sup>2</sup> = 0.086, Marginal R<sup>2</sup> = 0.062")) %>%
  fmt_number(columns=c(Variable, Estimate, Std..Error, df, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Std..Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("(Intercept)"~"Intercept", # change the way the model variables are labeled
                  "treatmentPhantom" ~ "Treatment (Phantom)",
                  "treatmentShifted" ~ "Treatment (Shifted)",
                  "year22018" ~ "Year (2018)") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>%
  tab_style(style = cell_borders(sides = c("top"),  weight = px(1.5), color='gray20'),
            locations = cells_source_notes()) %>%
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # set alignment of column labels for most columns 
            locations = cells_column_labels(columns = c(Estimate, Std..Error, df, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "center"), # set alignment of column text for most columns 
            locations = cells_body(columns = c(Estimate, Std..Error, df, t.value))) %>%
  tab_style(style = cell_text(align = "right"), # set alignments for CI column
            locations = cells_body(columns = lowerCI)) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  tab_footnote(footnote = "estimate reflects the difference in intercept between the specified treatment and control",
               locations = cells_body(columns = Variable, rows = c(2,3))) %>%
  tab_footnote(footnote = "estimate reflects the difference in intercept between 2018 and 2017",
               locations = cells_body(columns = Variable, rows = 4)) %>%
  cols_align_decimal(columns = c(Estimate, Std..Error, df, t.value)) %>% # align text within columns using decimal
  opt_table_font(google_font(name="Open Sans")) # set font for entire table

rank_table

gtsave(rank_table, "RAC_RankAbund_Table.pdf", path=here("Tables/RAC_SppTreatAssoc"))

#### Table for Gains Model ####
gains_CI <- data.frame(confint(mgainc)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>% slice(., -(1:2)) %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

gains_coef <- data.frame(summary(mgainc)$coefficients) %>% # put the model fixed effects into data frame
  rownames_to_column(., var="Variable") %>%
  rename(p.value = Pr...t..) %>%
  left_join(., gains_CI)

r2(mgainc) # look at conditional and marginal R2
# add values into table below

gains_table <- gains_coef %>% gt() %>%
  tab_header("Turnover: Species Gains") %>%
  tab_source_note(html("Conditional R<sup>2</sup> = 0.706, Marginal R<sup>2</sup> = 0.157")) %>%
  fmt_number(columns=c(Variable, Estimate, Std..Error, df, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Std..Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("(Intercept)"~"Intercept", # change the way the model variables are labeled
                  "treatmentPhantom" ~ "Treatment (Phantom)",
                  "treatmentShifted" ~ "Treatment (Shifted)",
                  "year22018" ~ "Year (2018)") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>%
  tab_style(style = cell_borders(sides = c("top"),  weight = px(1.5), color='gray20'),
            locations = cells_source_notes()) %>%
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # set alignment of column labels for most columns 
            locations = cells_column_labels(columns = c(Estimate, Std..Error, df, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "center"), # set alignment of column text for most columns 
            locations = cells_body(columns = c(Estimate, Std..Error, df, t.value))) %>%
  tab_style(style = cell_text(align = "right"), # set alignments for CI column
            locations = cells_body(columns = lowerCI)) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  tab_footnote(footnote = "estimate reflects the difference in intercept between the specified treatment and control",
               locations = cells_body(columns = Variable, rows = c(2,3))) %>%
  tab_footnote(footnote = "estimate reflects the difference in intercept between 2018 and 2017",
               locations = cells_body(columns = Variable, rows = 4)) %>%
  cols_align_decimal(columns = c(Estimate, Std..Error, df, t.value)) %>% # align text within columns using decimal
  opt_table_font(google_font(name="Open Sans")) # set font for entire table

gains_table

gtsave(gains_table, "RAC_Gains_Table.pdf", path=here("Tables/RAC_SppTreatAssoc"))

#### Table for Losses Model ####
losses_CI <- data.frame(confint(mlossc)) %>% # put 95% CI for the model parameters into data frame
  rownames_to_column(., var="Variable") %>% slice(., -(1:2)) %>%
  rename(lowerCI = X2.5.., upperCI=X97.5..)

losses_coef <- data.frame(summary(mlossc)$coefficients) %>% # put the model fixed effects into data frame
  rownames_to_column(., var="Variable") %>%
  rename(p.value = Pr...t..) %>%
  left_join(., losses_CI)

r2(mlossc) # look at conditional and marginal R2
# add values into table below

losses_table <- losses_coef %>% gt() %>%
  tab_header("Turnover: Species Losses") %>%
  tab_source_note(html("Conditional R<sup>2</sup> = 0.598, Marginal R<sup>2</sup> = 0.173")) %>%
  fmt_number(columns=c(Variable, Estimate, Std..Error, df, t.value, p.value, lowerCI, upperCI), decimals=3) %>% # report numbers in these columns to 3 decimal places
  cols_merge(c(lowerCI, upperCI), pattern="{1}, {2}") %>% # put upper and lower CI into one column
  cols_label(Std..Error = "SE", t.value = "t", lowerCI = "95% CI") %>% # change column labels
  cols_hide(p.value) %>%
  text_case_match("(Intercept)"~"Intercept", # change the way the model variables are labeled
                  "treatmentPhantom" ~ "Treatment (Phantom)",
                  "treatmentShifted" ~ "Treatment (Shifted)",
                  "year22018" ~ "Year (2018)") %>%
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for horizontal border under model variables
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>%
  tab_style(style = cell_borders(sides = c("top"),  weight = px(1.5), color='gray20'),
            locations = cells_source_notes()) %>%
  tab_style(style = cell_text(style="italic"), # change plot heading to be italics
            locations = cells_title(groups= 'title')) %>% 
  tab_style(style = cell_text(align = "center"), # set alignment of column labels for most columns 
            locations = cells_column_labels(columns = c(Estimate, Std..Error, df, t.value, lowerCI ))) %>%
  tab_style(style = cell_text(align = "center"), # set alignment of column text for most columns 
            locations = cells_body(columns = c(Estimate, Std..Error, df, t.value))) %>%
  tab_style(style = cell_text(align = "right"), # set alignments for CI column
            locations = cells_body(columns = lowerCI)) %>%
  tab_style(style = cell_text(style="italic"), # make font for column label italic
            locations = cells_column_labels(columns = t.value)) %>% 
  tab_footnote(footnote = "estimate reflects the difference in intercept between the specified treatment and control",
               locations = cells_body(columns = Variable, rows = c(2,3))) %>%
  tab_footnote(footnote = "estimate reflects the difference in intercept between 2018 and 2017",
               locations = cells_body(columns = Variable, rows = 4)) %>%
  cols_align_decimal(columns = c(Estimate, Std..Error, df, t.value)) %>% # align text within columns using decimal
  opt_table_font(google_font(name="Open Sans")) # set font for entire table

losses_table

gtsave(losses_table, "RAC_Losses_Table.pdf", path=here("Tables/RAC_SppTreatAssoc"))

#####################################################################################################
#### Species-Treatment Associations #####
#####################################################################################################
head(keep)

# separate out data by treatment year
# 2017
keep2017 <- keep %>% filter(year=="2017")

# 2018
keep2018 <- keep %>% filter(year=="2018")

# 2016
keep2016 <- keep %>% filter(year=="2016")

# put abundance data into separate data frame
abund2017 <- keep2017 %>% ungroup() %>% dplyr::select(Anpa:Tabr)
abund2018 <- keep2018 %>% ungroup() %>% dplyr::select(Anpa:Tabr)
abund2016 <- keep2016 %>% ungroup() %>% dplyr::select(Anpa:Tabr)

# put treatment info for each location into vector
treat2017 <- keep2017 %>% ungroup() %>% 
  pull(treatment)
treat2017 # still has ocean as a level, need to remove
t2017 <- fct_drop(treat2017)

treat2018 <- keep2018 %>% ungroup() %>% 
  pull(treatment)
treat2018 # still has ocean as a level, need to remove
t2018 <- fct_drop(treat2018)

#### Run Association Analysis ####

# 4 association indices are possible in indispecies package
# two association indices IndVal and IndValg
# two correlation indices r and r.g
# IndValg and r.g are group-equalized, so they account for different number of sites within each treatment -> we want this!
# correlation indices indicate the degree of preference for the treatment compared to the other treatments
# whereas indicator value indices assess how much the treatment matches the set of treatments where the species is found
# see: https://esajournals.onlinelibrary.wiley.com/doi/10.1890/08-1823.1
# from this paper:
# "For determining the ecological preference of a given species among a set of alternative site groups, 
# the correlation approach is probably more useful than the indicator value approach,
# because the former naturally allows the detection of negative preferences."

# decision: using r.g, which is group-equalized correlation index

?strassoc
# this function tests the strength of species site-group associations
# it returns a matrix of association values
# a bootstrap confidence interval can be calculated for each association value

# for 2017
rg17 <- strassoc(abund2017, cluster=t2017, func="r.g", nboot.ci=9999)

# for 2018
rg18 <- strassoc(abund2018, cluster=t2018, func="r.g", nboot.ci=9999)


# save models as RDS
#saveRDS(rg17, here("Model_Objects", "SppTreatAssoc_2017.rds"))
#saveRDS(rg18, here("Model_Objects", "SppTreatAssoc_2018.rds"))

# read in saved outputs
rg17 <- readRDS(here("Model_Objects", "SppTreatAssoc_2017.rds"))
rg18 <- readRDS(here("Model_Objects", "SppTreatAssoc_2018.rds"))


# get scientific names for species
spp_scinames <- read.csv(here("Data", "PhantomSurfBatTraits.csv"), header=T) %>% 
  mutate(sciname = str_replace(Scientific_Name, "_", " ")) %>%
  rename(Species=Species_Code) %>%
  dplyr::select(Species, sciname)

# put results into data frame to be made into a table

# for 2017
corr17 <- round(rg17$stat, 3) %>% as.data.frame() %>% # get correlation scores
  rename(C_corr = C, P_corr = P, S_corr = S)
lowCI17 <- round(rg17$lowerCI, 3) %>% as.data.frame() %>% # get lower 95% CI boundary
  rename(C_lowCI = C, P_lowCI = P, S_lowCI = S)
upperCI17 <- round(rg17$upperCI, 3) %>% as.data.frame() %>% # get upper 95% CI boundary
  rename(C_upperCI = C, P_upperCI = P, S_upperCI = S)
sppassoc.results17 <- bind_cols(corr17, lowCI17, upperCI17) %>% # put correlations and 95% CI together into one data frame
  rownames_to_column(., var="Species") %>% left_join(., spp_scinames) %>%
  dplyr::select(sciname, C_corr:S_upperCI) %>%
  mutate(Cupper= if_else(C_upperCI > 0, "G", "L"), # add columns that state whether the 95% CI overlaps zero or not
         Clower = ifelse(C_lowCI<0, "L", "G"),
         C_OverlapZero=ifelse(Cupper == Clower, "No", "Yes"),
         Pupper= if_else(P_upperCI > 0, "G", "L"),
         Plower = ifelse(P_lowCI<0, "L", "G"),
         P_OverlapZero=ifelse(Pupper == Plower, "No", "Yes"),
         Supper= if_else(S_upperCI > 0, "G", "L"),
         Slower = ifelse(S_lowCI<0, "L", "G"),
         S_OverlapZero=ifelse(Supper == Slower, "No", "Yes")
  ) %>%
  dplyr::select(-Cupper, -Clower, -Supper, -Slower, -Pupper, -Plower) # remove some columns that are not needed

head(sppassoc.results17)

# implement same steps as above for 2018
corr18 <- round(rg18$stat, 3) %>% as.data.frame() %>%
  rename(C_corr = C, P_corr = P, S_corr = S)
lowCI18 <- round(rg18$lowerCI, 3) %>% as.data.frame() %>%
  rename(C_lowCI = C, P_lowCI = P, S_lowCI = S)
upperCI18 <- round(rg18$upperCI, 3) %>% as.data.frame() %>%
  rename(C_upperCI = C, P_upperCI = P, S_upperCI = S)
sppassoc.results18 <- bind_cols(corr18, lowCI18, upperCI18) %>%
  rownames_to_column(., var="Species") %>% left_join(., spp_scinames) %>%
  dplyr::select(sciname, C_corr:S_upperCI) %>%
  mutate(Cupper= if_else(C_upperCI > 0, "G", "L"),
         Clower = ifelse(C_lowCI<0, "L", "G"),
        C_OverlapZero=ifelse(Cupper == Clower, "No", "Yes"),
        Pupper= if_else(P_upperCI > 0, "G", "L"),
        Plower = ifelse(P_lowCI<0, "L", "G"),
        P_OverlapZero=ifelse(Pupper == Plower, "No", "Yes"),
        Supper= if_else(S_upperCI > 0, "G", "L"),
        Slower = ifelse(S_lowCI<0, "L", "G"),
        S_OverlapZero=ifelse(Supper == Slower, "No", "Yes")
        ) %>%
  dplyr::select(-Cupper, -Clower, -Supper, -Slower, -Pupper, -Plower)

head(sppassoc.results18)

# Make tables of results for species-treatment associations for 2017 and 2018

# results table for 2017
sppassoc.table17 <-
  sppassoc.results17 %>%
  gt() %>%
  cols_merge(c(C_corr, C_lowCI, C_upperCI), pattern="{1} ({2}, {3})") %>% # merge correlation score with 95%CI
  cols_merge(c(P_corr, P_lowCI, P_upperCI), pattern="{1} ({2}, {3})") %>%
  cols_merge(c(S_corr, S_lowCI, S_upperCI), pattern="{1} ({2}, {3})") %>%
  cols_label(sciname = "Species", C_corr = "Control", P_corr = "Phantom", S_corr = "Shifted") %>% # change column labels
  data_color(columns = C_OverlapZero, rows = C_OverlapZero == "No" & C_corr > 0, # conditionally color text of cells. purple for positive correlations, green for negative
             target_columns = C_corr, palette = "#793C8E", apply_to="text") %>%
  data_color(columns = C_OverlapZero, rows = C_OverlapZero == "No" & C_corr < 0,
             target_columns = C_corr, palette = "#009B5F", apply_to="text") %>%
  data_color(columns = P_OverlapZero, rows = P_corr <0 & P_OverlapZero == "No" ,
             target_columns = P_corr, palette = "#009B5F", apply_to="text") %>%
  data_color(columns = S_OverlapZero, rows = S_OverlapZero == "No" & S_corr > 0,
             target_columns = S_corr,  palette = "#793C8E", apply_to="text") %>%
  data_color(columns = S_OverlapZero, rows = S_OverlapZero == "No" & S_corr < 0,
             target_columns = S_corr, palette = "#009B5F", apply_to="text") %>%
  tab_style( style = cell_text(weight = "bolder"), # make colored correlations have bold text
             locations = cells_body(columns = C_corr, rows = C_OverlapZero == "No")) %>%
  tab_style( style = cell_text(weight = "bolder"), 
             locations = cells_body(columns = P_corr, rows = P_OverlapZero == "No")) %>%
  tab_style( style = cell_text(weight = "bolder"), 
             locations = cells_body(columns = S_corr, rows = S_OverlapZero == "No")) %>%
  cols_hide(columns = c(C_OverlapZero, P_OverlapZero, S_OverlapZero)) %>% # hide several columns
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
    heading.title.font.size = 15, # set heading font size
    table.font.size=14, # set font size for remaining aspects of table
    heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
    table_body.border.bottom.color = "gray20", # set color for bottom horizontal border of table
    table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
    table.border.bottom.color = "gray20", # set color for very bottom horizontal line
    column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_text(style="italic"), # make font for species column italic
    locations = cells_body(columns = sciname)) %>% 
  tab_style(style = cell_text(align = "center"), # align column labels for Control, Phantom and Shifted to center
            locations = cells_column_labels(columns = c(C_corr, P_corr, S_corr))) %>%
  tab_style(style = cell_text(align = "right"), # align column text for C, S and P to right
            locations = cells_body(columns = c(C_corr, P_corr, S_corr))) %>%
  tab_header(title = "2017") %>%
  opt_table_font(google_font(name="Open Sans")) # set font for entire table

  
# results table for 2018
sppassoc.table18 <-
  sppassoc.results18 %>%
  gt() %>%
  cols_merge(c(C_corr, C_lowCI, C_upperCI), pattern="{1} ({2}, {3})") %>%
  cols_merge(c(P_corr, P_lowCI, P_upperCI), pattern="{1} ({2}, {3})") %>%
  cols_merge(c(S_corr, S_lowCI, S_upperCI), pattern="{1} ({2}, {3})") %>%
  cols_label(sciname = "Species", C_corr = "Control", P_corr = "Phantom", S_corr = "Shifted") %>%
  data_color(columns = C_OverlapZero, rows = C_OverlapZero == "No" & C_corr > 0,
             target_columns = C_corr, palette = "#793C8E", apply_to="text") %>%
  data_color(columns = C_OverlapZero, rows = C_OverlapZero == "No" & C_corr < 0,
             target_columns = C_corr, palette = "#009B5F", apply_to="text") %>%
  data_color(columns = P_OverlapZero, rows = P_corr <0 & P_OverlapZero == "No" ,
             target_columns = P_corr, palette = "#009B5F", apply_to="text") %>%
  data_color(columns = S_OverlapZero, rows = S_OverlapZero == "No" & S_corr > 0,
             target_columns = S_corr,  palette = "#793C8E", apply_to="text") %>%
  data_color(columns = S_OverlapZero, rows = S_OverlapZero == "No" & S_corr < 0,
             target_columns = S_corr, palette = "#009B5F", apply_to="text") %>%
  tab_style( style = cell_text(weight = "bolder"), 
             locations = cells_body(columns = C_corr, rows = C_OverlapZero == "No")) %>%
  tab_style( style = cell_text(weight = "bolder"), 
             locations = cells_body(columns = P_corr, rows = P_OverlapZero == "No")) %>%
  tab_style( style = cell_text(weight = "bolder"), 
             locations = cells_body(columns = S_corr, rows = S_OverlapZero == "No")) %>%
  cols_hide(columns = c(C_OverlapZero, P_OverlapZero, S_OverlapZero)) %>% # hide several columns
  tab_options(column_labels.font.size = 13, column_labels.font.weight = "bold", # make font of column headings bold
              heading.title.font.size = 15, # set heading font size
              table.font.size = 14, # set font size for remaining aspects of table
              heading.border.bottom.color = "gray20", # set color for horizontal border under table heading
              table_body.border.bottom.color = "gray20", # set color for bottom horizontal border of table
              table.border.top.color = "white", # set color for very top horizontal line (above plot heading). using white to make invisible
              table.border.bottom.color = "gray20", # set color for very bottom horizontal line
              column_labels.border.bottom.color = "gray20") %>% 
  tab_style(style = cell_text(style="italic"), # make font for species column italic
            locations = cells_body(columns = sciname)) %>% 
  tab_style(style = cell_text(align = "center"), # align column labels for Control, Phantom and Shifted to center
            locations = cells_column_labels(columns = c(C_corr, P_corr, S_corr))) %>%
  tab_style(style = cell_text(align = "right"), # align column text for C, S and P to right
            locations = cells_body(columns = c(C_corr, P_corr, S_corr))) %>%
  tab_header(title = "2018") %>%
  opt_table_font(google_font(name="Open Sans")) # set font for entire table

# view tables
sppassoc.table17
sppassoc.table18 

# save tables as pdf and png
gtsave(data = sppassoc.table17, filename = "SpeciesTreatment_Assoc_2017_Table.pdf", path = here("Tables/RAC_SppTreatAssoc"))
gtsave(data = sppassoc.table17, filename = "SpeciesTreatment_Assoc_2017_Table.png", path = here("Tables/RAC_SppTreatAssoc"))

gtsave(data = sppassoc.table18, filename = "SpeciesTreatment_Assoc_2018_Table.pdf", path = here("Tables/RAC_SppTreatAssoc"))
gtsave(data = sppassoc.table18, filename = "SpeciesTreatment_Assoc_2018_Table.png", path = here("Tables/RAC_SppTreatAssoc"))
  