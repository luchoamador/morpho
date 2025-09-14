####Map of Chile
###Luis Amador

library(ggmap)
library(ggplot2)
library(raster)
library(maptools)

# Map of Chile -----------------------------------------------------

mapa <- borders("world", regions = c("Chile" ), fill = "grey70", colour = "black")

ggplot() + mapa + theme_bw() + xlab("Longitude (decimals)") + ylab("Latitude (decimals)") + 
  theme(panel.border = element_blank(), panel.grid.major = element_line(colour = "grey80"), panel.grid.minor = element_blank())
