# ==============================================================================
# Title: Creation of the plots of the Usage note section of the article:  
# HERA: a high-resolution pan-European hydrological reanalysis (1950-2020)
# Author: Alois Tilloy - Joint Research Centre - Unit C6 
# Date: 2024 -02 -01 
# Description:
#   This script allows to generate plots that compare hydrological regimes 
#   of selected rivers at different points in time
# Note:
# Please note that the script requires transformed inputs from the HERA dataset
# ==============================================================================


# Import functions from function file ------------------------------------------
source("~/06_Floodrivers/DataPaper/Code/HERA/functions.R")

main_path = 'D:/tilloal/Documents/06_Floodrivers/' ### CHANGE THIS PATH
valid_path = paste0(main_path,'DataPaper/')

hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]

outletname="efas_rnet_100km_01min"
#Create the upArea df that will be used to match with station data
nca=nc_open(paste0(valid_path, 'GIS/upArea_European_01min.nc'))
name.lon="lon"
name.lat="lat"
nav=names(nca[['var']])
#Band1 is the second variable
t=nca$var[[2]]
name.var=names(nca$var)[2]
tsize<-t$varsize
tdims<-t$ndims
nt1<-tsize[tdims]
lon=ncvar_get(nca,name.lon)
lat=ncvar_get(nca,name.lat)
start <- rep(1,tdims) # begin with start=(1,1,1,...,1)
count <- tsize 
#Here I need to extract only some values as this is the full EFAS domain
upArea   = ncvar_get(nca,name.var,start = start, count= count) 
#convert to sqkm
upArea=upArea/1000000


# Load data ---------------------------------------------------------------------
#To get all time steps for discharge data, 
#I post processed HERA yearly file to obtain 88 netcdf files which are spatial chunks of the HERA domain
#See the subspace_efas.csv file to see the extent of the chunks


# Selection of river pixels to be checked

#The Ebro is Spain's longest river, with low and high water levels alternating throughout the year, 
#influenced by winter snowmelt and summer evaporation/human usage. 
#The river is vital for agriculture.
Ebro=list(Nsq=42,
river=data.frame("lon"=-0.825,"lat"=41.608),
catch="Ebro @ Zaragoza")

#The Ardèche is Mediterranean river mostly known for tourism due 
#to its scenic gorges, but floods and droughts can impact the local economy and environment.
Ardeche=list(Nsq=42,
river=data.frame("lon"=4.658,"lat"=44.258),
catch="Ardeche @ Saint-Martin")

#The Rhône is France's most powerful river,characterized by a significant seasonal variation in flow rates.
#The Rhône River is crucial for transportation, hydropower generation, and irrigation in the region.
Rhone=list(Nsq=52,
river=data.frame("lon"=4.891,"lat"=45.772),
catch="Rhone @ Lyon")

#The Po is the longest river in Italy
#The Po river is crucial for transportation and irrigation in northern Italy.
Po=list(Nsq=52,
river=data.frame("lon"=11.60,"lat"=44.89),
catch="Po @ Ferrara")

#The Danube is the longest river in the EU
#the Danube River is an important waterway for international trade, connecting several countries in Central and Eastern Europe.
Danube=list(Nsq=63,
river=data.frame("lon"=16.64,"lat"=48.13),
catch="Danube @ Vienna")

#The Warta River, which flows through Poznan in western Poland, is a major waterway with a total length of 1,004 kilometers.
#It is an important source of freshwater for agriculture, industry, and urban water supply in the region.
Warta=list(river=data.frame("lon"=16.941,"lat"=52.408),
Nsq=61,
catch="Warta @ Poznan")

#The Vistula River, which flows through Warsaw is the longest and largest river in Poland.
#It is an important source of freshwater for agriculture, industry, and urban water supply in the region.
Vistula=list(river=data.frame("lon"=21.03,"lat"=52.24),
Nsq=72,
catch="Vistula @ Warsaw")

# The Schelde rivers spans approximately 350 kilometers, it is influenced by a maritime climate
# It Experiences flow rate fluctuations and tidal influences due to being an estuary
Schelde=list(river=data.frame("lon"=3.774,"lat"=51.058),
             Nsq=40,
             catch="Schelde @ Gent")

Allriver=list(Ebro,Ardeche,Rhone,Po,Danube,Warta,Vistula,Schelde)


# Creating plots for all the rivers----------------------------------------------
for (ir in 1:length(Allriver)){
  Riviere=Allriver[[ir]]
  Nsq=Riviere$Nsq
  river=Riviere$river
  catch=Riviere$catch
  Cairo::Cairo(
    20, #length
    15, #width
    file = paste("Regimes/",catch, ".png", sep = ""),
    type = "png", #tiff
    bg = "transparent", #white or transparent depending on your requirement 
    dpi = 300,
    units = "cm" #you can change to pixels etc 
  )
  Regime_comp(rspace,Nsq,river,hydroDir,catch,outletname)
  dev.off()
}