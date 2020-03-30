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

# Most basic bubble plot
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

