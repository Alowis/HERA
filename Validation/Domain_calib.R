# ==============================================================================
# Title: Plot of the calibrated domain of the HERA dataset, linked to the article:  
# HERA: a high-resolution pan-European hydrological reanalysis (1951-2020)
# Author: Alois Tilloy - Joint Research Centre - Unit C6 
# Date: 2024 -02 -01 
# Description:
#   This script generates a plot of the calibrated domain of HERA (See article supplement)
# ==============================================================================

#script to identify how much of HERA has been calibrated
library(rgdal)
library(raster)
library(rgdal)
library(raster)
library(ncdf4)
library(lubridate)
library(ggplot2)
library(rasterVis)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rgeos)


workDir<-("//ies.jrc.it/H07/nahaUsers/tilloal/ERA5l_x_lisflood/EFAS1rcmin_GitLFS_efas5_24032023/")
dataDir<-("//ies.jrc.it/H07/nahaUsers/tilloal/ERA5l_x_lisflood/EFAS1rcmin_GitLFS_efas5_24032023/parameters")
hydroDir<-("D:/tilloal/Documents/06_Floodrivers/DataPaper/")
#Import reservoir locations from netcdf
resOpen=function(dir,outletname){
  ncbassin=paste0(dir,outletname)
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncb,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat)
  lla=length(latdat)
  start=c(1,1)
  count=c(llo,lla)
  
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$res=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outfinal=outll[which(!is.na(outll$res)),]
  return (outfinal)
}
#Step 1: load catchment list and coordinates

data <- read.table(paste0(dataDir,"/efas_all_inter_catchments_coords.txt"))
datacl=read.csv(paste0(dataDir,"/ID_closest_geographyclimate_EFASnext_regionalization.csv"),header=TRUE)
#Step 2: Identify catchments that do not have a calibration points

datano=inner_join(data,datacl,by=c("V1"="UncalibratedID"))

datamin=data[which(data$V1=="1111111"),]
#Step 3: overlap these catchments with my HERA domain

#load the netcdf 
calname="/efas_all_inter_catchments.nc"
nccal=paste0(dataDir,calname)
ncc=nc_open(nccal)
name.vb=names(ncc[['var']])
namev=name.vb[2]


name.lon="lon"
name.lat="lat"
start=c(1,1)
londat = ncvar_get(ncc,name.lon) 
llo=length(londat)
latdat = ncvar_get(ncc,name.lat)
lla=length(latdat)
count=c(llo,lla)
idcal = ncvar_get(ncc,namev,start = start, count= count) 
idcal=as.vector(idcal)
outll=expand.grid(londat,latdat)
outll$V1=idcal

nocalpix=inner_join(outll,datano,by="V1")
nocalpix$latlong=paste(round(nocalpix$Var1,2),round(nocalpix$Var2,2))
ouptix=outll[which(outll$V1=="1111111"),]

#Step 4: Check plot

### Plot parameters ----
cord.dec=outll[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
Impdates=seq(1950,2020,by=10)
valuenames=paste0("Y",Impdates)
basemap=w2

calP <- st_as_sf(nocalpix, coords = c("Var1", "Var2"), crs = 4326)
calP <- st_transform(calP, crs = 3035)
## Mean flood date plot
paletx=c(hcl.colors(8, palette = "PurpOr", alpha = NULL, rev = TRUE, fixup = TRUE))
tsize=16
osize=16
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=calP,aes(geometry=geometry),col="orange",alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_colourbar(barwidth = 2, barheight = 15,reverse=T))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  ggtitle("ungauged catchments")

ggsave( filename="D:/tilloal/Documents/06_Floodrivers/DataPaper/Plots/Ungauged_catchments.jpg", width=16.3, height=15, units="cm",dpi=500)

#Step 5: join these pixels with my HERA domain

#load the area of HERA
GridHERA=raster( paste0(hydroDir,"/HERA_domain_01min.tif"))
GHERA=as.data.frame(GridHERA,xy=T)
GHERA=GHERA[which(!is.na(GHERA$HERA_domain_01min)),]
#load the netcdf of pixel area
aname="/pixarea_European_01min.nc"
nccal=paste0(workDir,aname)
nca=nc_open(nccal)
name.vb=names(nca[['var']])
namev=name.vb[2]
#time <- ncvar_get(ncb,"time")

#timestamp corretion
name.lon="lon"
name.lat="lat"
start=c(1,1)
count=c(llo,lla)
londat = ncvar_get(nca,name.lon,start=start[1],count=count[1]) 
llo=length(londat)
latdat = ncvar_get(nca,name.lat,start=start[2],count=count[2])
lla=length(latdat)
idcal = ncvar_get(nca,namev,start = start, count= count) 
idcal=as.vector(idcal)
max(idcal,na.rm=T)
area=expand.grid(londat,latdat)
#lonlatloop=expand.grid(c(1:llo),c(1:lla))
area$V1=idcal/1000000
area$latlong=paste(round(area$Var1,4),round(area$Var2,4),sep=" ")
#Joining area with HERA domain
GHERA$latlong=paste(round(GHERA$x,4),round(GHERA$y,4),sep=" ")

mag=match(GHERA$latlong,area$latlong)
GHERA$area=area$V1[mag]
areaH=sum(GHERA$area)

nocalpix$latlong=paste(round(nocalpix$Var1,4),round(nocalpix$Var2,4),sep=" ")
ouptix$latlong=paste(round(ouptix$Var1,4),round(ouptix$Var2,4),sep=" ")


anog=inner_join(GHERA,nocalpix,by="latlong")
anout=inner_join(GHERA,ouptix,by="latlong")

domcal=match(anog$latlong,GHERA$latlong)
ancal1=GHERA[-domcal,]
  
domcal2=match(anout$latlong,ancal1$latlong)
ancal2=ancal1[-domcal2,]

areaC=sum(ancal2$area)
areaN=sum(anog$area)
areaO=sum(anout$area)
areaC+areaO+areaN
areaN/areaH

### Plot parameters ----
outletname="efas_rnet_100km_01min"
outll=outletopen(hydroDir,outletname)
cord.dec=outll[,c(2,3)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
basemap=w2

calP1 <- st_as_sf(anout, coords = c("Var1", "Var2"), crs = 4326)
calP1 <- st_transform(calP1, crs = 3035)

calP2 <- st_as_sf(anog, coords = c("Var1", "Var2"), crs = 4326)
calP2 <- st_transform(calP2, crs = 3035)

calP3 <- st_as_sf(ancal2, coords = c("x", "y"), crs = 4326)
calP3 <- st_transform(calP3, crs = 3035)
## Mean flood date plot
paletx=c(hcl.colors(8, palette = "PurpOr", alpha = NULL, rev = TRUE, fixup = TRUE))
tsize=16
osize=16
ggplot(w2) +
  geom_sf(fill="white")+
  geom_sf(data=calP1,aes(geometry=geometry),col="red",alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(data=calP2,aes(geometry=geometry),col="orange",alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(data=calP3,aes(geometry=geometry),col="royalblue",alpha=1,size=0.25,stroke=0,shape=15)+ 
  geom_sf(fill=NA, color="grey") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude")+
  guides(color = guide_colourbar(barwidth = 2, barheight = 15,reverse=T))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        panel.grid.major = element_line(colour = "grey70"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave( filename="D:/tilloal/Documents/06_Floodrivers/DataPaper/Plots/Ungauged_catchments3.jpg", width=16.3, height=15, units="cm",dpi=500)





