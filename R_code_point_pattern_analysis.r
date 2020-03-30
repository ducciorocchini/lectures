# point oattern analysis frst part

################ covid data
# https://services.arcgis.com/5T5nSi527N4F7luB/arcgis/rest/services/COVID_19_HistoricCasesByCountry(pt)View/FeatureServer
library(ggplot2)

setwd("~/lab/")

# setwd("C:/RStudio/data_book/raster_data/final")

covid <- read.table("covid_agg.csv", head=T)

plot(covid$country,covid$cases)
plot(covid$country,covid$cases,las=1)

plot(covid$country,covid$cases)
plot(covid$country,covid$cases,las=1) # y axis horizontal
plot(covid$country,covid$cases,las=2) # x and y horizontal
plot(covid$country,covid$cases,las=3) # x axis horizontal

# Plotting data by ggplot2
# Every ggplot2 plot has three key components:
# 1. data,
# 2. A set of aesthetic mappings between variables in the data and visual
# properties, and
# 3. At least one layer which describes how to render each observation. Layers
# are usually created with a geom function.
# Here’s a simple example:
ggplot(mpg, aes(x = displ, y = hwy)) +
geom_point()


ggplot(covid, aes(x=lon, y=lat, size = cases)) +
    geom_point(alpha=0.3)
    


download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip", 
              destfile = 'coastlines.zip')

unzip(zipfile = "coastlines.zip", 
      exdir = 'ne-coastlines-10m')

coastlines <- readOGR("ne-coastlines-10m/ne_10m_coastline.shp")
      
plot(dcov)
points(covid,pch=19,cex=0.5,col="red")
plot(coastlines, add = TRUE)

