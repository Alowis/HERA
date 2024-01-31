##########################################################################################
############   THIS SCRPT IS FOR RUNNING 1 SQUARE OF THE DOMAIN ON THE HPC  ##############
##########################################################################################
#Import functions from TSEVA
source("~/LFRuns_utils/TS-EVA/functions.R")

#Library importation
suppressWarnings(suppressMessages(library(ncdf4)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(rnaturalearth)))
suppressWarnings(suppressMessages(library(rnaturalearthdata)))
suppressWarnings(suppressMessages(library(rgeos)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(lubridate)))


#Functions

main_path = 'D:/tilloal/Documents/06_Floodrivers/' ### CHANGE THIS PATH
valid_path = paste0(main_path,'DataPaper/')
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
count <- tsize # begin w/count=(nx,ny,nz,...,nt), reads entire var
#Here I need to extract only some values as this is the full EFAS domain
upArea   = ncvar_get(nca,name.var,start = start, count= count) 

#convert to sqkm
upArea=upArea/1000000



#Step 1: load data
#Squares that I have on this machine
Nsq = 42
outlets="RNetwork"

# haz="flood"
var = "dis"

# workDir = "/BGFS/CLIMEX/tilloal/HydroMeteo/"
# setwd(workDir)
# hydroDir<-"/BGFS/CLIMEX/tilloal/HydroMeteo/Timeseries/dis6_Uncal2"

hydroDir<-("D:/tilloal/Documents/LFRuns_utils/data")
#workDir<-("D:/tilloal/Documents/06_Floodrivers/dis")
rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]



nrspace=rspace[Nsq,]
outletname="efas_rnet_100km_01min"
nameout="UCRnet"
outEFAS=outletopen(hydroDir,outletname)


unikout=outEFAS$outlets
outEFAS$latlong=paste(round(outEFAS$Var1,4),round(outEFAS$Var2,4),sep=" ")

lagaronne=outEFAS[which(outEFAS$Var1<=-0.7 & outEFAS$Var2<45.4),]
lagaronne=lagaronne[which(lagaronne$Var1>=-0.8),]
lagaronne=lagaronne[which(lagaronne$Var2>45.32),]
outlag=lagaronne$outlets[1]

laloire=outEFAS[which(outEFAS$Var1<=4.1 & outEFAS$Var1>4.09),]
laloire=laloire[which(laloire$Var2>=46.1 & laloire$Var2<46.15),]


#Load the file
#loading the files as netcdf (needs to be checked offline)
Nsq=42
#load specific outlet file for this square:
outletname="efas_rnet_100km_01min"
nameout="UCRnet"
rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]
nrspace=rspace[Nsq,]
outsq=outletopen(hydroDir,outletname,nrspace)
Idstart=as.numeric(Nsq)*10000
if (length(outsq$outlets)>0){
  outsq$outlets=seq((Idstart+1),(Idstart+length(outsq$outlets)))
}


filename=paste0("timeseries/validations/dis_",Nsq,"_1950_2020_cf")

#Extract two times eries: 1951-1981 & 1990-2020

timebound1=c(as.POSIXct(("1951-01-01 00:00:00")))
timebound2=c(as.POSIXct(("1981-12-31 00:00:00")))

tb1=match(timebound1,df.dis$timeStamps)
tb2=match(timebound2,df.dis$timeStamps)

# Selection of pixels to be checked


Ebro=list(Nsq=42,
river=data.frame("lon"=-0.825,"lat"=41.608),
catch="Ebro @ Zaragoza")

Ardeche=list(Nsq=42,
river=data.frame("lon"=4.658,"lat"=44.258),
catch="Ardeche @ Saint-Martin")

Rhone=list(Nsq=52,
river=data.frame("lon"=4.891,"lat"=45.772),
catch="Rhone @ Lyon")

Po=list(Nsq=52,
river=data.frame("lon"=11.60,"lat"=44.89),
catch="Po @ Ferrara")

Danube=list(Nsq=63,
river=data.frame("lon"=16.64,"lat"=48.13),
catch="Danube @ Vienna")

Warta=list(river=data.frame("lon"=16.941,"lat"=52.408),
Nsq=61,
catch="Warta @ Poznan")

Vistula=list(river=data.frame("lon"=21.03,"lat"=52.24),
Nsq=72,
catch="Vistula @ Warsaw")


bug=list(river=data.frame("lon"=25.42,"lat"=57.54),
Nsq=71,
catch="bug @ bugland")
Nsq=bug$Nsq
river=bug$river
nrspace=rspace[Nsq,]
outsq=outletopen(hydroDir,outletname,nrspace)
Idstart=as.numeric(Nsq)*10000
if (length(outsq$outlets)>0){
  outsq$outlets=seq((Idstart+1),(Idstart+length(outsq$outlets)))
}

filename=paste0("timeseries/validations/dis_",Nsq,"_1950_2020_cf")

riv=st_as_sf(river, coords = c("lon", "lat"), crs = 4326)
outloc=st_as_sf(outsq, coords = c("Var1", "Var2"), crs = 4326)
oula=st_distance(riv,outloc)
oula=oula/1000
md=min(oula)
if (as.numeric(md)<1){
  rloc=which(oula==md)
}else{
  print("No river found")
  break
}

riv=as.matrix(round(outsq[rloc,c(2,3)],3))
outsq$V1r=round(outsq$Var1,3)
outsq$V2r=round(outsq$Var2,3)

dists=disNcopenloc(filename,hydroDir,outsq,rloc)
min(dists$outlets)
timeStamps=unique(as.Date(dists$time,origin="1979-01-01"))
timeStamps=as.POSIXct(timeStamps-1/24)
txx=timeStamps
dists$timeStamps=timeStamps
plot(dists$timeStamps[53000:54500],dists$outlets[53000:54500])

# river=data.frame("lon"=1.825,"lat"=50.09)
# Nsq=40
# catch="Somme @ Abbeville"


# river=data.frame("lon"=2.009,"lat"=50.374)
# Nsq=40
# catch="Canche @ Hesdin"

# river=data.frame("lon"=13.405,"lat"=52.5229)
# Nsq=61
# catch="Spree @ Berlin"

Schelde=list(river=data.frame("lon"=3.774,"lat"=51.058),
Nsq=40,
catch="Schelde @ Gent")

Allriver=list(Ebro,Ardeche,Rhone,Po,Danube,Warta,Vistula,Schelde)

Regime_comp(rspace,Nsq,river,hydroDir,catch)

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
  Regime_comp(rspace,Nsq,river,hydroDir,catch)
  dev.off()
}



#Sicily catchments
#Load model
#Also strange catchment in france

Imera=list(river=data.frame("lon"=13.945,"lat"=37.174),
         Nsq=66,
         catch="Imera @ Drasi")

Platani=list(river=data.frame("lon"=13.657,"lat"=37.476),
           Nsq=66,
           catch="Platani @ Passofonduto",
           sname="PASSOFONDUTO")

Simeto=list(river=data.frame("lon"=14.7933,"lat"=37.6585),
           Nsq=66,
           catch="Simeto @ Ponte Maccarrone",
           sname="PONTE MACCARRONE")


Francheville=list(river=data.frame("lon"=4.57,"lat"=47.758),
            Nsq=41,
            catch="Seine @ Francheville",
            sname="H0100010")


VEP=list(river=data.frame("lon"=4.62,"lat"=48.74),
                  Nsq=40,
                  catch="Saulx @ Vitry",
                  sname="H5172010")


Loranca=list(river=data.frame("lon"=-3.09,"lat"=40.457),
         Nsq=32,
         catch="loranca",
         sname="2003003")


#new technique
Roya=list(   catch="loranca",
             sname="2003003")

#Find Nsq with my station
sID=6139140
Roya=kgeout[which(kgeout$V1==sID),]
#load nrspace

rspace= read.csv(paste0(hydroDir,"/subspace_efas.csv"))
rspace=rspace[,-1]
#ido
rs1=rspace[which(rspace$stlon<=Roya$idlo & rspace$endlon>=Roya$idlo),]
rs2=rs1[which(rs1$stlat<=Roya$idla & rs1$endlat>=Roya$idla),]
Nsq=rs2$spID
nrspace=rs2



outall=outletopen(hydroDir,outletname)
outall$idlalo=paste(outall$idlo,outall$idla,sep=" ")
matcho=match(Roya$idlalo,outall$idlalo)
pic=outall[matcho,]
outsq=outletopen(hydroDir,outletname,nrspace)
Idstart=as.numeric(Nsq)*10000
if (length(outsq$outlets)>0){
  outsq$outlets=seq((Idstart+1),(Idstart+length(outsq$outlets)))
}

filename=paste0("timeseries/validations/dis_",Nsq,"_1950_2020_cf")

riv=st_as_sf(river, coords = c("lon", "lat"), crs = 4326)
outloc=st_as_sf(outsq, coords = c("Var1", "Var2"), crs = 4326)
oula=st_distance(riv,outloc)
oula=oula/1000
md=min(oula)
if (as.numeric(md)<1){
  rloc=which(oula==md)
}else{
  print("No river found")
  break
}

riv=as.matrix(round(outsq[rloc,c(2,3)],3))
outsq$V1r=round(outsq$Var1,3)
outsq$V2r=round(outsq$Var2,3)

pic$vv=paste(pic$Var1,pic$Var2,sep=" ")
outsq$vv=paste(outsq$Var1,outsq$Var2,sep=" ")
rloc=match(pic$vv,outsq$vv)
rloc


dists=disNcopenloc(filename,hydroDir,outsq,rloc)
min(dists$outlets)
timeStamps=unique(as.Date(dists$time,origin="1979-01-01"))
timeStamps=as.POSIXct(timeStamps-1/24)
txx=timeStamps
dists$timeStamps=timeStamps
EFAS_flow=dists
EFAS_flow$dQ=tsEvaNanRunningMean(EFAS_flow$outlets,4)
selectix=EFAS_flow$time-floor(EFAS_flow$time)
EFAS_flday=EFAS_flow[which(selectix==0.5),]

#Now load the observations
#What is the ID of the station?

#load efas matches
efas_stations=read.csv("efas_flooddriver_match.csv")
mystat=efas_stations[which(efas_stations$National_Station_Identifier==Loranca$sname),]

#load domis stations
Q_stations <- st_read(paste0(valid_path,'Europe_daily_combined_2022_v2_WGS84_NUTS.shp'))  # SHP with gauge locations
mystat=Q_stations[which(Q_stations$StationID==sID),]

#loop on all years to load the timeserie of that station:
y=1950
s=mystat$StationID
Q_data <- read.csv(paste0('Q_', y, '.csv'), header = F)  # CSVs with observations
Station_data_IDs <- as.vector(t(Q_data[1, ]))
ix <- which(!is.na(match(Station_data_IDs,s)))
Q_obs <- Q_data[-1, ix]
print(Q_obs)
yrlist=c(1951:2020)
for (y in yrlist){
  print(y)
  Q_data <- read.csv(paste0('Q_', y, '.csv'), header = F)  # CSVs with observations
  Station_data_IDs <- as.vector(t(Q_data[1, ]))
  ix <- which(!is.na(match(Station_data_IDs,s)))
  Q_s <- Q_data[-1, ix]
  Q_obs=c(Q_obs,Q_s)
}
plot(Q_obs)

Q_obs=Q_obs[-c(1,2)]

EFAS_flday$Qobs=Q_obs

Q_comp=EFAS_flday[which(!is.na(EFAS_flday$Qobs)),]

b1=2000
b2=3200
plot(Q_comp$timeStamps[b1:b2],Q_comp$Qobs[b1:b2],type="l",col="blue",xlab="Time",ylab="Q (m3/s)",main = paste0("hydrograph @ ",mystat$StationName,". River: ",mystat$River))
points(Q_comp$timeStamps[b1:b2],Q_comp$dQ[b1:b2], type="l",col="red")
legend(Q_comp$timeStamps[900],30, legend=c("simulated", "observed"),
       col=c("red", "blue"),lty=1, cex=0.8)


Cairo::Cairo(
  20, #length
  15, #width
  file = paste("Regimes/",Francheville$catch, "flows.png", sep = ""),
  type = "png", #tiff
  bg = "transparent", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
plot(Q_comp$timeStamps,Q_comp$Qobs,type="l",col="blue",xlab="Time",ylab="Q (m3/s)",main = paste0("hydrograph @ ",mystat$StationName,". River: ",mystat$River))
points(Q_comp$timeStamps,Q_comp$dQ, type="l",col="red")
legend(Q_comp$timeStamps[900],30, legend=c("simulated", "observed"),
       col=c("red", "blue"),lty=1, cex=0.8)
dev.off()


data=Q_comp[,c(6,8,7)]
data[,2]=tsEvaNanRunningMean(data[,2],120)
data[,3]=tsEvaNanRunningMean(data[,3],120)
Cairo::Cairo(
  20, #length
  15, #width
  file = paste("Regimes/",Simeto$catch, ".png", sep = ""),
  type = "png", #tiff
  bg = "transparent", #white or transparent depending on your requirement 
  dpi = 300,
  units = "cm" #you can change to pixels etc 
)
test1=plotQjm (data,catch="Roya",UpAloc=mystat$DrainingArea.km2.Provider) 
dev.off()

Q_out=Q_comp[,c(6,8,7,3,4)]
names(Q_out)[3] = "Qsim"
write.csv(Q_out,file=paste0("out/Q",mystat$StationName,".csv"))


#Q loranca from Spanish platforms
loranca=read.csv("3003_Loranca_de_Tajuña.csv",sep=";")

Q_comp$timeD=as.Date(Q_comp$timeStamps)
loranca$timeD=as.Date(loranca$fecha,format="%d/%m/%Y")

Qcomp=inner_join(Q_comp,loranca, by="timeD")


plot(Qcomp$timeD,Qcomp$Qobs,type="l",col="blue",xlab="Time",ylab="Q (m3/s)",main = paste0("hydrograph @ ",mystat$StationName,". River: ",mystat$River))
points(Qcomp$timeD,Qcomp$dQ, type="l",col="red")
points(Qcomp$timeD,Qcomp$caudal, type="l",col="green")
legend(Qcomp$timeStamps[900],30, legend=c("simulated", "observed"),
       col=c("red", "blue"),lty=1, cex=0.8)


