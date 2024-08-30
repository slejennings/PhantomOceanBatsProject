#### PHANTOM OCEAN BAT ANALYSES ####
#### Author(s): Sarah L Jennings ####
#### Most Recent Update: August 28, 2024 ####
#### Create Fig 2 in Manuscript (Parts A and B) of Acoustic Conditions in Different Treatments ####

# load packages
library(tidyverse)
library(showtext)
library(colorspace)
library(here)
library(ggdensity)
library(patchwork)

# import files
C <- read.delim(here("Figures", "Control512.txt"))
O <- read.delim(here("Figures", "Ocean512.txt"))
P <- read.delim(here("Figures", "Phantom512.txt"))
S <- read.delim(here("Figures", "Shifted512.txt"))

# Note about files: these text files were created by opening a sound file from each treatment in Audacity
# In Audacity, I used the Analyze > Plot Spectrum option on each file
# I exported the FFT data as a txt file using the following settings Algorithm: Spectrum, Size: 512, Function: Hann Window, Axis: Linear Frequency
# I then had to manually modify the files a little.
# This was mainly because plotting them in R was tricky
# I used geom_polygon in R for plotting each and this needed equal y-axis points on both ends of the data to draw a straight line bottom across the polygon
# I suspect there is a better way to do this, but I got this to work and decided not to spend more time figuring out another approach
# The edits I made to the txt files did not change the appearance of the spectra other than forcing R to create nice looking polygons for each treatment

# add label for treatment name
C2 <- C %>% mutate(Type = "Control")
O2 <- O %>% mutate(Type = "Ocean")
P2 <- P %>% mutate(Type = "Phantom")
S2 <- S %>% mutate(Type = "Shifted")

# combine into single data frame
sounds <- bind_rows(C2, O2, P2, S2)

# we want to adjust the peak amplitude of all files to the max of Control
# find max peak amplitudes of each treatment
maxamp <- sounds %>% group_by(Type) %>%
  summarize(max = max(Level..dB.))

# Control treatment has the lowest max
# determine amount to subtract from all other treatments
maxampsub <- maxamp %>%
  rowwise() %>%
  mutate(ctrlmax = -52,
         sub = abs(ctrlmax - max)) %>%
  dplyr::select(Type, sub)

sounds_sub <- sounds %>% left_join(., maxampsub) %>%
  rowwise() %>%
  mutate(Level..dB.adj = Level..dB. - sub)

# load font and color scheme for plots
font_add_google(name="Open Sans", family="opensans")
showtext_auto()

hcl_palettes(palette = "Viridis") # load viridis palette in colorspace package
viridis4 <- specplot(sequential_hcl(4, "Viridis"), plot=T)
print(viridis4) # get HCL, RGB and hex for 4 color viridis palette
viridis4$hex # just look at HEX

# specify the order of the colors as they relate to the levels of the treatment (Control, Ocean, Phantom, Shifted)
col4 <- c( "#4B0055", "#00BE7D", "#007094", "#FDE333") # to match with supplementary plot Clint made

# make plot with legend at bottom of the plot
# convert Hz to kHz to make x-axis less cluttered

powerspectra <- ggplot(sounds_sub, aes(x=(Frequency..Hz./1000), y=Level..dB.adj, color = factor(Type), fill=factor(Type))) + 
  geom_polygon(linewidth=1, alpha = 0.15) +
  scale_color_manual(values = col4, name="") +
  scale_fill_manual(values = col4, name="") +
  theme_classic()+
  labs(x = "Frequency (kHz)", y = "Amplitude (dBA)") +
  scale_x_continuous(limits=c(0,14), breaks = seq(0, 14, 2), expand=c(-0.005,0))+
  scale_y_continuous(limits=c(-90,-50), breaks = seq(-80, 0, 10), expand=c(-0.02,0)) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = c("center", "top"),
    legend.box.spacing = unit(5, "pt"),
    legend.margin = margin(1, 1, 1, 1),
    legend.text = element_text(size= 9, family="opensans", margin = margin(l= 1)),
    axis.title.x = element_text(size = 11,family="opensans", margin = margin(t=3, unit="pt")),
    axis.title.y = element_text(size = 11, family="opensans", margin = margin(r=2, unit="pt")),
    axis.text.x = element_text(size=10, family="opensans"),
    axis.text.y = element_blank(),
    axis.ticks.length=unit(.25, "cm")
  )

# with custom labels on y-axis tick marks
powerspectra <- ggplot(sounds_sub, aes(x=(Frequency..Hz./1000), y=Level..dB.adj, color = factor(Type), fill=factor(Type))) + 
  geom_polygon(linewidth=1, alpha = 0.15) +
  scale_color_manual(values = col4, name="") +
  scale_fill_manual(values = col4, name="") +
  theme_classic()+
  labs(x = "Frequency (kHz)", y = "Amplitude (dBA)") +
  scale_x_continuous(limits=c(0,14), breaks = seq(0, 14, 2), expand=c(-0.005,0))+
  scale_y_continuous(limits=c(-90,-50), breaks = seq(-80, -50, 10), expand=c(-0.02,0), labels = c("y1", "y2", "y3", "y4")) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.justification = c("center", "top"),
    legend.box.spacing = unit(5, "pt"),
    legend.margin = margin(1, 1, 1, 1),
    legend.text = element_text(size= 9, family="opensans", margin = margin(l= 1)),
    axis.title.x = element_text(size = 11,family="opensans", margin = margin(t=3, unit="pt")),
    axis.title.y = element_text(size = 11, family="opensans", margin = margin(r=2, unit="pt")),
    axis.text.x = element_text(size=10, family="opensans"),
    axis.text.y = element_text(size=10, family="opensans"),
    axis.ticks.length=unit(.25, "cm")
  )
powerspectra


###########################
### Plots of Average Acoustic Conditions within Each Treatment
#### Import data
# we need data with 2016 and ocean treatments but where bat abundance has been averaged based on the sampling days
bats.w2016.ave <- readRDS(here("Output", "bats_w2016_ave.rds")
  
# plot using ggdensity
library(ggdensity)

# put each treatment in a separate facet
sounddensity_facet <- ggplot(data = bats.w2016.ave, aes(x = (median/1000), y = leq, fill=treatment)) +
  geom_hdr(probs=c(0.5,0.75,0.95), method=method_mvnorm()) +
  geom_point(shape=21, show.legend = F) +
  theme_minimal()+
  scale_fill_manual(guide="none", values = col4) +
  facet_wrap(~treatment, labeller=labeller(treatment=c("C"="Control", "O"="Ocean", "P"="Phantom", "S"="Shifted"))) +
  labs(y="Amplitude (dBA)", x ="Median Frequency (kHz)",
       fill="Treatment", alpha="") +
  theme(
    panel.border = element_rect(linewidth = 0.5, color="black", fill=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box.spacing = unit(5, "pt"),
    legend.justification = c("center", "top"),
    legend.margin = margin(1, 1, 1, 1),
    strip.text.x = element_text(size = 11, family="opensans" ), # facet label font
    legend.text = element_text(size = 9, family="opensans"),
    axis.title.x = element_text(size = 11, family="opensans", margin = margin(t=3, unit="pt")),
    axis.title.y = element_text(size = 11, family="opensans", margin = margin(r=2, unit="pt")),
    axis.text.x = element_text(size = 10, family="opensans"),
    axis.text.y = element_text(size = 10, family="opensans")
  )
sounddensity_facet


#### Combine the two plots into a single figure using patchwork package ####

# figure out the layout for the combined plot
# going to use patchwork package to do this step
# I want to arrange the plots vertically
# the first plot should be shorter and less wide
# second plot should be taller and slightly wider

layout <- c(
  patchwork::area(t=1, l=1, b=2, r=3), # t = top, l = left, b = bottom, r = right 
  patchwork::area(t=3, l=1, b=6, r=4)
)

# look at this layout
plot(layout)

# combine two plots using the specified layout
powerspectra / sounddensity_facet + # put powerspectra above sounddensity_face
  plot_layout(design = layout) + # use layout determined above
  plot_annotation(tag_levels ='a', tag_prefix = "(", tag_suffix = ")") & # add labels for (a) and (b)
  theme(plot.tag = element_text(size = 12, family="opensans", face="bold")) # set font and size for annotations

ggsave(here("Figures", "Fig2_AcousticTreatments.pdf"), width = 6, height = 8.5)
