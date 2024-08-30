#### PHANTOM OCEAN BAT ANALYSES ####
#### Author(s): Sarah L Jennings ####
#### Most Recent Update: August 28, 2024 ####
#### Make Map of Study Area (Figure 1 in the manuscript) ####

#### Load Packages and Import Data ####
library(ggmap)
library(ggplot2)
library(sf)
library(sp)
library(maps)
library(mapdata)
library(showtext)
library(ggspatial)
library(here)
library(maps)
library(mapdata)
library(tidyverse)
library(colorspace)

# import locations
locs <- read.csv(here("Data", "point_coordinates.csv"))
head(locs)

# import bat data
batswide <-read.csv(here::here("Data", "bats_freq_amp_11_24_20.csv"), header=T)
# not all sites/cluster have processed bat data (1 cluster/3 sites were dropped)
# need to drop the cluster and 3 sites that were dropped

# get treatments for sites
bats_treat <- batswide %>%
  filter(year=="2017") %>% # we want treatments in experimental years
  dplyr::select(cluster, site, point, treatment) %>%
  distinct()

# join treatments with coordinates. This adds treatment and simultaneously removes extra cluster/sites
locs <- bats_treat %>% 
  left_join(., locs) %>%
  rename(Treatment = treatment)
  
# perform a couple of checks  
length(unique(locs$cluster)) == length(unique(batswide$cluster)) # for cluster
# also check sites
length(unique(locs$site)) == length(unique(batswide$site))

# keep only middle location and use that to mark sites on the map
length(unique(locs$site))
locs <- locs %>% filter(point == "PM1")
nrow(locs) # should be equal to 19

locs_sf <- locs %>% st_as_sf(coords=c("Long", "Lat")) %>%
  st_set_crs(4326) %>% # WSG 84
  st_transform(3857) # transform to 3857 which is what google maps uses
  
st_crs(locs_sf) # confirm the crs is 3857

# register google key for pulling map layer
register_google(key = "AIzaSyCnG7wMso6s0OS5PZsCoA7af8wVqJBHg8o")

# set font
showtext_auto()
font_add_google(name="Open Sans", family="opensans")

# first, get map of CA
states <- map_data("state")
california <- subset(states, region %in% "california")

# make data frame containing the location of the study site to be placed on the map of CA
vandenburg_loc <- data.frame(long=-120.5724, lat=34.6820)

# make map of CA
CA_map <- ggplot(california) + 
  geom_polygon(aes(x=long, y = lat, group = group), fill="#DAEECA", color="gray30") + # use same green as google maps to fill state
  coord_fixed(1.3) + 
  theme_nothing() + # remove all axis and gridlines
  geom_point(data=vandenburg_loc, aes(x=long, y = lat), size = 10, pch=0, color="black", stroke = 2) + # add square for study site location
  annotate("text", x = -119.5, y= 37, label = "California", size=7, family="opensans", angle=311) + # put name of the state down the middle
  theme(panel.background = element_rect(fill="#92CCE8"))
CA_map

ggsave(CA_map, filename="CaliforniaMap.pdf", path=here("Figures"))

###### make Vandenburg Base Map with Sites Marked

# pull google maps layer
vandenburg <- get_googlemap(center = c(lon=-120.57,lat=34.71),
                                      zoom=11,
                           maptype = "terrain", 
                           style = 'feature:administrative.country%7Celement:geometry.stroke%7Ccolor:0x7f7f7f&style=feature:administrative.country%7Celement:labels%7Cvisibility:off&style=feature:administrative.land_parcel%7Cvisibility:off&style=feature:administrative.locality%7Cvisibility:off&style=feature:administrative.neighborhood%7Cvisibility:off&style=feature:administrative.province%7Celement:geometry.stroke%7Ccolor:0x165413&style=feature:administrative.province%7Celement:labels.text%7Ccolor:0x015286%7Cvisibility:simplified&style=feature:landscape%7Ccolor:0xd5efc7&style=feature:poi%7Cvisibility:off&style=feature:road%7Cvisibility:off&style=feature:transit%7Cvisibility:off&style=feature:water%7Celement:geometry.fill%7Ccolor:0x92cce8&style=feature:water%7Celement:labels%7Cvisibility:off')


# now we have to do some weird stuff
# ggsn package used to work with ggmap to plot scalebars and north arrows
# this no longer works because of all the spatial packages that were retired recently
# following a hack I found here: https://stackoverflow.com/questions/47749078/how-to-put-a-geom-sf-produced-map-on-top-of-a-ggmap-produced-raster
# this will make the ggmap work with the ggspatial package that can add scalebars and north arrows

# define a function to fix the bbox of the ggmap to be in EPSG:3857
ggmap_bbox <- function(vandenburg) {
  if (!inherits(vandenburg, "ggmap")) stop("map must be a ggmap object")
  # Extract the bounding box (in lat/lon) from the ggmap to a numeric vector, 
  # and set the names to what sf::st_bbox expects:
  map_bbox <- setNames(unlist(attr(vandenburg, "bb")), 
                       c("ymin", "xmin", "ymax", "xmax"))
 
   # Convert the bbox to an sf polygon, transform it to 3857, 
  # and convert back to a bbox (convoluted, but it works)
  bbox_3857 <- st_bbox(st_transform(st_as_sfc(st_bbox(map_bbox, crs = 4326)), 3857))
  
  # Overwrite the bbox of the ggmap object with the transformed coordinates 
  attr(vandenburg, "bb")$ll.lat <- bbox_3857["ymin"]
  attr(vandenburg, "bb")$ll.lon <- bbox_3857["xmin"]
  attr(vandenburg, "bb")$ur.lat <- bbox_3857["ymax"]
  attr(vandenburg, "bb")$ur.lon <- bbox_3857["xmax"]
  vandenburg
}

# Use the function on the downloaded google map:
van_map <- ggmap_bbox(vandenburg)

# look at what we've done so far
ggmap(van_map) # doesn't work unless we set the crs with coord_sf
ggmap(van_map) + coord_sf(crs = st_crs(3857))
# this looks good but we need to crop it because it is wider that we need it to be

# to crop the map, we need to set the xlim and ylim argument in coord_sf 
# in order to do this, we need to find out the coordinates in the 3857 crs
# I visually selected x and y limits in WSG 84 and am converting them
lim_coords <- data.frame(Lat=c(34.54,34.82), Long=c(-120.72,-120.48)) %>% # put the limits I selected into a data frame
  st_as_sf(coords=c("Long", "Lat")) %>% # convert them to sf
  st_set_crs(4326)  %>% # tell it that they are WSG84
  st_transform(3857) # transform them to EPSG:3857
lim_coords

# we also need to do another stupid hack
# if I try to add the site coordinates in locs_sf using geom_sf, it overrides the cropping of the map
# as a work around, I am converting the geometry of points back into a data frame with long and lat columns
# then plotting them as geom_point
head(locs_sf) # this what we are starting with
locs_df <- locs_sf %>%  dplyr::mutate(long = sf::st_coordinates(.)[,1],
                                             lat = sf::st_coordinates(.)[,2])
# locs_df has geometry and the columns long and lat that can be provided to geom_point

# set colors for sites
colors <- c("#D55E00", "#353E7C", "#FDE333", "#E69F00", "#009B95", "#007094", "#4B0055")

vanbg_map <- ggmap(van_map) + 
  coord_sf(crs = st_crs(3857), # set the crs for the map
           xlim=c(-13411772, -13438489), ylim=c(4101543,4139447)) + # set the boundaries for the map using lim_coords
 geom_point(data=locs_df, aes(x=long, y=lat, fill = cluster, shape = Treatment), 
            size = 3.1, color="gray20",
            show.legend = c(fill=F, shape=T)) +
  annotation_north_arrow(location = "bl", which_north = "true", # add north arrow to bottom left from ggspatial package
                         width= unit(1, "cm"), height=unit(1.2, "cm"), # set size of arrow
                         pad_x = unit(0.5, "cm"), # move arrow away from map edges slightly
                         pad_y = unit(0.5, "cm")) + 
  annotation_scale(location = "tl", text_family="opensans") + # add a scalebar
  labs(x="Longitude", y="Latitude", fill="Treatment") +
  theme( panel.border = element_rect(color="black", fill=NA, linewidth = 1),
    axis.title.x = element_text(size = 12, family="opensans"),
         axis.title.y = element_text(size = 12, family="opensans", margin = margin(r=-6, unit="pt")),
         axis.text.x = element_text(size = 11, family="opensans"),
         axis.text.y = element_text(size = 11, family="opensans"),
         legend.title = element_text(size = 11, family="opensans"),
         legend.text = element_text(size = 11, family="opensans"),
         legend.position = "top",
         legend.key = element_rect(fill = NA)) +
  scale_x_continuous(breaks = c(-120.7, -120.6, -120.5)) + # make the map have fewer x-axis tick mark labels
  scale_fill_manual(values=colors) +
  scale_shape_manual(values = c(21, 22, 23, 24)) 
  
vanbg_map

ggsave(vanbg_map, filename="VandenburgSiteMap.pdf", path=here("Figures"), width=6, height=6, units="in")

#### used powerpoint to put together CA map, the Site Map and to create a zoomed in schematic of the sampling locations/speakers within each site