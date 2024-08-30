#### PHANTOM OCEAN BAT ANALYSES ####
#### Author(s): Sarah L Jennings ####
#### Most Recent Update: August 28, 2024 ####
#### Make Interaction Plots for Certain Individual Species Models ####

#### Load Packages and Import Data ####
library(here)
library(spaMM)
library(tidyverse)
library(ggplot2)
library(performance)
library(showtext)
library(RColorBrewer)
library(colorspace)
library(pdp)
library(patchwork)
library(viridis)


# import data
batslong <- readRDS(here("Output", "batslong.rds"))

# subset three species
myca <- batslong %>% filter(., species == "Myca")
myyu <- batslong %>% filter(., species == "Myyu")
labl <- batslong %>% filter(., species == "Labl")

# run models
# reading in already saved models often doesn't work for plotting
myca.m <- fitme(dets ~ scale(leq)*scale(median) + scale(millum) + year + pc1 + pc2 + richness +
                                 (1|pointID) + Matern(1 | Long + Lat)+ AR1(1|julian %in% year), data = myca, fixed=list(nu=0.5), family = negbin2(), method="ML")

myyu.m <- fitme(dets ~ scale(leq)*scale(median) + scale(millum) + year + pc1 + pc2 + richness +
                  (1|pointID) + Matern(1 | Long + Lat)+ AR1(1|julian %in% year), data = myyu, fixed=list(nu=0.5), family = negbin2(), method="ML")

labl.m <- fitme(dets ~ scale(leq)*scale(median) + scale(millum) + year + pc1 + pc2 + richness +
                  (1|pointID) + Matern(1 | Long + Lat)+ AR1(1|julian %in% year), data = labl, fixed=list(nu=0.5), family = negbin2(), method="ML")

par.myca <- partial(myca.m, pred.var = c("leq", "median"), chull=T)
saveRDS(par.myca, "~/Desktop/Bat_Results/Archived Files/Species Models/partialMyca.rds")
par.myca <- readRDS(here("Figures", "partialMyca.rds"))

plot.myca <- plotPartial(par.myca, contour=T)

par.labl <- partial(labl.m, pred.var = c("leq", "median"), chull=T)
saveRDS(par.labl, "~/Desktop/Bat_Results/Archived Files/Species Models/partialLabl.rds")
plot.labl <- plotPartial(par.labl, contour=T)
par.labl <- readRDS(here("Figures", "partialLabl.rds"))

par.myyu <- partial(myyu.m, pred.var = c("leq", "median"), chull=T)
saveRDS(par.myyu, "~/Desktop/Bat_Results/Archived Files/Species Models/partialMyyu.rds")
par.myyu <- readRDS(here("Figures", "partialMyyu.rds"))

plot.myyu <- plotPartial(par.myyu, contour=T)
class(plot.myyu)

# make nicer plots in ggplot

# load font to use for table and figure
font_add_google(name="Open Sans", family="opensans")
showtext_auto()

lablplot <- ggplot(par.labl, aes(x = leq, y = median, z=yhat, fill = yhat)) +
  geom_tile() +
  geom_contour(color = "white", alpha = 0.5) +
  scale_fill_viridis(option="D" , name="Predicted activity") +
  labs(x = "Sound pressure level (dBA)", y = "Median frequency (Hz)", subtitle = "Lasiurus blossevilli") +
  theme_classic(base_size = 14, base_family = "") +
  theme(
    legend.title = element_text(size=10, family="opensans", angle=90, hjust=0.5),
    legend.title.position = "left",
    legend.text = element_text(size=10, family="opensans"),
    axis.title.x = element_text(size = 12,family="opensans", margin = margin(t=10, unit="pt")),
    axis.title.y = element_text(size = 12, family="opensans"),
    axis.text.x = element_text(size=12, family="opensans"),
    axis.text.y = element_text(size=12, family="opensans"),
    plot.subtitle = element_text(size=12, family="opensans", face="italic"),
    legend.box.spacing = unit(0.01, "pt")
  )

mycaplot <- ggplot(par.myca, aes(x = leq, y = median, z=yhat, fill = yhat)) +
  geom_tile() +
  geom_contour(color = "white", alpha = 0.5) +
  scale_fill_viridis(option="D" , name="Predicted activity") +
  labs(x = "Sound pressure level (dBA)", y = "Median frequency (Hz)", subtitle = "Myotis californicus") +
  theme_classic(base_size = 14, base_family = "") +
  theme(
    legend.title = element_text(size=10, family="opensans", angle=90, hjust=0.5),
    legend.title.position = "left",
    legend.text = element_text(size=10, family="opensans"),
    axis.title.x = element_text(size = 12,family="opensans", margin = margin(t=10, unit="pt")),
    axis.title.y = element_text(size = 12, family="opensans"),
    axis.text.x = element_text(size=12, family="opensans"),
    axis.text.y = element_text(size=12, family="opensans"),
    plot.subtitle = element_text(size=12, family="opensans", face="italic"),
    legend.box.spacing = unit(0.01, "pt")
  )

myyuplot <- ggplot(par.myyu, aes(x = leq, y = median, z=yhat, fill = yhat)) +
  geom_tile() +
  geom_contour(color = "white", alpha = 0.5) +
  scale_fill_viridis(option="D" , name="Predicted activity") +
  labs(x = "Sound pressure level (dBA)", y = "Median frequency (Hz)", subtitle = "Myotis yumanensis") +
  theme_classic(base_size = 14, base_family = "") +
  theme(
    legend.title = element_text(size=10, family="opensans", angle=90, hjust=0.5),
    legend.title.position = "left",
    legend.text = element_text(size=10, family="opensans"),
    axis.title.x = element_text(size = 12,family="opensans", margin = margin(t=10, unit="pt")),
    axis.title.y = element_text(size = 12, family="opensans"),
    axis.text.x = element_text(size=12, family="opensans"),
    axis.text.y = element_text(size=12, family="opensans"),
    plot.subtitle = element_text(size=12, family="opensans", face="italic"),
    legend.box.spacing = unit(0.01, "pt")
  )
 

interaction.plots <- lablplot + mycaplot + myyuplot + 
  plot_layout(ncol=1) + 
  plot_annotation(tag_levels ='a', tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(size = 12, family="opensans", face="bold", vjust=2)) 
          
interaction.plots
ggsave("Figures/InteractionPlots.pdf", width = 6, height = 4.5, units = 'in')
