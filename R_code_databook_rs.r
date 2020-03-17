
##-- packages:
install.packages(“RStoolbox”)

##-- functions:
https://github.com/bleutner/RStoolbox/tree/master/R

##-- data:
http://book.ecosens.org/software/rstoolbox/

sudo apt-get install p7zip-full
7z x data_book.7z 

## landsat bands
https://landsat.usgs.gov/what-are-band-designations-landsat-satellites

##-- Code

setwd("~/duccio/Documents/lectures_and_seminars/monitoring_ecosystems_unibo/data_book/raster_data/final
")

# setwd("C:/RStudio/data_book/raster_data/final")

install.packages("raster") #----- installazione pacchetti
install.packages("rgdal") #----- installazione pacchetti
install.packages("ggplot2")

# install.packages("RStoolbox")

library(raster)
library(rgdal)
# library(RStoolbox)
library(ggplot2)
# library(hexbin)

# https://landsat.gsfc.nasa.gov/the-worldwide-reference-system/

# 2011 image
p224r063_2011 <- brick("p224r63_2011.grd")

p224r63_2011

 summary(p224r63_2011)
            B1_sre     B2_sre     B3_sre    B4_sre    B5_sre B6_bt     B7_sre
Min.    0.00000000 0.01154835 0.00730000 0.0000000 0.0000000 295.1 0.00000000
1st Qu. 0.01371136 0.03108290 0.02000000 0.2516721 0.1137448 296.4 0.04183780
Median  0.01654248 0.03526897 0.02386578 0.2842704 0.1290826 296.9 0.04887195
3rd Qu. 0.02157205 0.04421049 0.03364492 0.3121567 0.1567999 298.1 0.06530904
Max.    0.10590097 0.14880101 0.24380812 0.5219069 0.3955572 304.4 0.31459978
NAs    0.00000000 0.00000000 0.00000000 0.0000000 0.0000000   0.0 0.00000000

p224r63_2011m <- brick("p224r63_2011_masked.grd")
p224r63_1988m <- brick("p224r63_1988_masked.grd")

# plotRGB(...)
pairs(p224r63_2011m) # time consuming

# 1   2   3   4
# b   g   r
#     b   g   r

# plotRGB(im,3,2,1)
plotRGB(p224r63_2011m,r=3,g=2,b=1)
plotRGB(p224r63_2011m,r=3,g=2,b=1,scale=1000,stretch="Lin")
plotRGB(p224r63_2011m,r=4,g=3,b=2,scale=1000,stretch="Lin")

par(mfrow=c(2,1))
plotRGB(p224r63_2011m,r=3,g=2,b=1,scale=1000,stretch="Lin")
plotRGB(p224r63_2011m,r=4,g=3,b=2,scale=1000,stretch="Lin")

# stretch
par(mfrow=c(2,1))
plotRGB(p224r63_2011m,r=4,g=3,b=2,scale=1000,stretch="Lin")
plotRGB(p224r63_2011m,r=4,g=3,b=2,scale=1000,stretch="hist")

# output
pdf("image.pdf")
plotRGB(p224r63_2011m,r=4,g=3,b=2,scale=1000,stretch="hist")
dev.off()

## resampling
# plot(p224r63_2011m)
agg <- aggregate(p224r63_2011m, fact=20)
p224r63_2011m_res <- resample(p224r63_2011m, agg)

par(mfrow=c(2,1))
plotRGB(p224r63_2011m,r=4,g=3,b=2,scale=1000,stretch="Lin")
plotRGB(p224r63_2011m_res,r=4,g=3,b=2,scale=1000,stretch="Lin")




par(mfrow=c(2,1))
plotRGB(p224r63_2011m,r=3,g=2,b=1,scale=1000,stretch="Lin",main="natural colours")
plotRGB(p224r63_2011m,r=4,g=3,b=2,scale=1000,stretch="Lin", main="infrared")

multitemp <-  p224r63_2011m - p224r63_1988m 

plot(multitemp)

##- spectral indices

ndvi2011 <- (p224r63_2011m$B4_sre-p224r63_2011m$B3_sre) / (p224r63_2011m$B4_sre+p224r63_2011m$B3_sre)
# ~time: 1min 

ndvi1988 <- (p224r63_1988m$B4_sre-p224r63_1988m$B3_sre) / (p224r63_1988m$B4_sre+p224r63_1988m$B3_sre)
# ~time: 1min
#--- potential function: spectralIndices

par(mfrow=c(2,1))
plot(ndvi1988,main="NDVI 1988")
plot(ndvi2011,main="NDVI 2011")

dif <- ndvi2011-ndvi1988

######!!!!!!PITFALL!!!!!!!
Always standardise the legends!

### comparing NDVi values over time

ndvi1988
 -0.001903766, 0.9689609 
ndvi2011
 -0.007789201, 1  

# par
par(mfrow=c(1,2))
hist(ndvi1988,ylim=c(0,2000000))
hist(ndvi2011,ylim=c(0,2000000))
 
# io:
# Histogram Colored (blue and red)
hist(ndvi1988, col=rgb(1,0,0,0.5),main="NDVI frequencies")
hist(ndvi2011, col=rgb(0,0,1,0.5), add=T)
box()
legend("topleft", c("NDVI 1988", "NDVI 2011"), fill=c("red", "blue"))

# library(hexbin)
# hbin <- hexbin(ndvi1988,ndvi2011, xbins = 40)
# plot(hbin)

## Plot RGB space
im <- brick("~/Downloads/ndvi-radar.png")

 im
class       : RasterBrick 
dimensions  : 3232, 4807, 15536224, 3  (nrow, ncol, ncell, nlayers)
resolution  : 1, 1  (x, y)
extent      : 0, 4807, 0, 3232  (xmin, xmax, ymin, ymax)
coord. ref. : NA 
data source : /home/duccio/Downloads/ndvi-radar.png 
names       : ndvi.radar.1, ndvi.radar.2, ndvi.radar.3 
min values  :            0,            0,            0 
max values  :          255,          255,          255 

4807 * 3232
[1] 15536224


plotRGB(im,3,2,1)
plotRGB(p224r63_2011m,r=4,g=3,b=2,scale=1000,stretch="Lin")
plotRGB(p224r63_2011m,r=3,g=2,b=1,scale=1000,stretch="Lin")

library(ggplot2)
ggRGB(p224r63_2011m,4,3,2)

par(mfrow=c(2,1))
plotRGB(p224r63_2011m,r=3,g=2,b=1,scale=1000,stretch="Lin",main="natural colours")
plotRGB(p224r63_2011m,r=4,g=3,b=2,scale=1000,stretch="Lin", main="infrared")

#hexbin
#### pca

pairs(p224r63_2011m)

library(RStoolbox)
plot(p224r63_2011m$B1_sre,p224r63_2011m$B1_sre)
  
#reample data
p224r63_2011m_res <- aggregate(p224r63_2011m, fact=3, fun=mean)
# original res=30, newres=90
# aggregate(x, fact=2, fun=mean, expand=TRUE, na.rm=TRUE, filename=


> plot(r)
> manner <- aggregate(r, fact=20)
# > resares <- resample(r, resa, method="ngb") 
# in case "fact=" does not work in Wondows use "fact<-"
> r_new <- resample(r, manner)



p224r63_2011m_pca <- rasterPCA(p224r63_2011m_res)
# 3 minutes needed

summary(p224r63_2011m_pca$model) 

Importance of components:
                          Comp.1      Comp.2       Comp.3       Comp.4
Standard deviation     1.2950291 0.052987610 0.0213916820 5.551811e-03
Proportion of Variance 0.9980317 0.001670837 0.0002723173 1.834234e-05
Cumulative Proportion  0.9980317 0.999702523 0.9999748401 9.999932e-01
                             Comp.5       Comp.6       Comp.7
Standard deviation     2.621003e-03 1.710617e-03 1.288550e-03
Proportion of Variance 4.088090e-06 1.741370e-06 9.880696e-07
Cumulative Proportion  9.999973e-01 9.999990e-01 1.000000e+00


loadings(p224r63_2011m_pca$model)

Loadings:
       Comp.1 Comp.2 Comp.3 Comp.4 Comp.5 Comp.6 Comp.7
B1_sre                0.116  0.317         0.604 -0.721
B2_sre                0.161  0.684 -0.276  0.294  0.581
B3_sre                0.274  0.559  0.194 -0.697 -0.296
B4_sre        -0.892 -0.432         0.115              
B5_sre        -0.412  0.663 -0.318 -0.526              
B6_bt  -0.999                                          
B7_sre        -0.162  0.508 -0.117  0.772  0.238  0.223

               Comp.1 Comp.2 Comp.3 Comp.4 Comp.5 Comp.6 Comp.7
SS loadings     1.000  1.000  1.000  1.000  1.000  1.000  1.000
Proportion Var  0.143  0.143  0.143  0.143  0.143  0.143  0.143
Cumulative Var  0.143  0.286  0.429  0.571  0.714  0.857  1.000

#--------

p224r63_2011m_pca

$call
rasterPCA(img = p224r63_2011m_res)

$model
Call:
princomp(cor = spca, covmat = covMat[[1]])

Standard deviations:
     Comp.1      Comp.2      Comp.3      Comp.4      Comp.5      Comp.6 
1.295029050 0.052987610 0.021391682 0.005551811 0.002621003 0.001710617 
     Comp.7 
0.001288550 

 7  variables and  494500 observations.

$map
class       : RasterBrick 
dimensions  : 500, 989, 494500, 7  (nrow, ncol, ncell, nlayers)
resolution  : 90, 90  (x, y)
extent      : 579765, 668775, -522735, -477735  (xmin, xmax, ymin, ymax)
coord. ref. : +proj=utm +zone=22 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 
data source : in memory
names       :         PC1,         PC2,         PC3,         PC4,         PC5,         PC6,         PC7 
min values  : -7.29596861, -0.21213869, -0.12206057, -0.07208760, -0.02024087, -0.03419134, -0.01445724 
max values  :  2.40711155,  0.33669578,  0.16446461,  0.09348111,  0.05662598,  0.02284201,  0.04294400 


attr(,"class")
[1] "rasterPCA" "RStoolbox"

#-------- Plot the map
plot(p224r63_2011m_pca$map)

plotRGB(p224r63_2011m_pca$map,r=1,g=2,b=3,scale=1000,stretch="Lin")

library(ggplot2)
ggRGB(p224r63_2011m_pca$map,1,2,3, stretch="lin")

#--- plot components
plot(p224r63_2011m_pca$map$PC1,p224r63_2011m_pca$map$PC2)

################### 
### USING EXTERNAL FUNCTIONS
################### 
source("spectralrao.r")

ndvi2011 <- (p224r63_2011m$B4_sre-p224r63_2011m$B3_sre) / (p224r63_2011m$B4_sre+p224r63_2011m$B3_sre)

plot(ndvi2011)

# resample NDVI
agg <- aggregate(ndvi2011, fact=20)
ndvi.r <- resample(ndvi2011, agg)

par(mfrow=c(2,1))
plot(ndvi2011,main="Original NDVI")
plot(ndvi.r,main="Resampled NDVI")

raomat <- spectralrao(ndvi.r, mode="classic", distance_m="euclidean",window=3, shannon=T) # 3min

par(mfrow=c(1,3))
plot(ndvi.r,main="Resampled NDVI")
plot(raster(raomat[[1]]), main="Rao's Q")
plot(raster(raomat[[2]]), main="Shannon's H'")

################### 
### MAPPING BIODIVERSITY
################### 

plot(NDVI_07_20085km)

sample <- spsample(pol_clip, n = 1000, "random")

# extract NDVI
extract(NDVI_07_20085km,sample)

# variance
sd <- focal(NDVI_07_20085km, w=matrix(1/9,nrow=3,ncol=3), fun=var)

# extract variance
extract(sd,sample)

# species
plot(species$sd,species$species_richness)
# plot(NDVI,species_richness)

lm(species$species_richness~species$sd)

Call:
lm(formula = species$species_richness ~ species$sd)

Coefficients:
(Intercept)   species$sd  
     -51.36      3752.94  

species_div = -51.36 + 3752.94*sd



 lm(species$species_richness~species$sd+species$NDVI)

Call:
lm(formula = species$species_richness ~ species$sd + species$NDVI)

Coefficients:
 (Intercept)    species$sd  species$NDVI  
      -428.1       27217.9        -779.5  

species_div = -428.1 + 27217.9*sd -779.5*NDVI_07_20085km

