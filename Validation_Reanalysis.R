# ==============================================================================
# Title: Creation of the validation metrics and plots of the HERA dataset, linked to the article:  
# HERA: a high-resolution pan-European hydrological reanalysis (1950-2020)
# Author: Alois Tilloy - Joint Research Centre - Unit C6 
# Date: 2024 -02 -01 
# Description:
#   This script allows to generate plots  and assess the skills of the HERA dataset against
#   observed discharge data from 2901 river gauges across Europe
# ==============================================================================

source("~/06_Floodrivers/DataPaper/Code/HERA/functions.R")

main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
dis_path<-paste0(main_path,'dis/calibrated/filtered/Histo/')
setwd(valid_path)


# Data generation -----------------------------
#Load all Q_EFAS (HERA simulated discharge) and create a single file
#create one file with all years
y=1950
Q_d <- read.csv(paste0('out/EFAS_', y, '.csv'), header = F)  # CSVs with observations
Station_data_IDs <- as.vector(t(Q_d[1, ]))
Q_sim <- Q_d
yrlist=c(1951:2020)
for (y in yrlist){
  print(y)
  Q_d <- read.csv(paste0('out/EFAS_', y, '.csv'), header = F)  # CSVs with observations
  Q_s <- Q_d[-1, ]
  Q_sim=rbind(Q_sim,Q_s)
}
write.csv(Q_sim,file="out/HERA_Val_19502020.csv")

# Part 1: Create the Valid station file -------------------------------------------

#Loading all txt files with matching locations
#load matching coordinates
Sloc=read.table(paste0(valid_path,"out/Stations_locations_EFAS_grid_1950.txt"),sep=",")
yrlist=c(1951:2020)
for (yi in yrlist){
  print(yi)
  Slocy=read.table(paste0(valid_path,"out/Stations_locations_EFAS_grid_",yi,".txt"),sep=",")
  Slocrep=Slocy[which(Sloc$V2==0),]
  Sloc[which(Sloc$V2==0),]=Slocrep
}
Sloc_final=Sloc[which(Sloc$V2!=0),]
Sloc_final$csource="SpatialQMatch"

#transalting to R indexing after python
Sloc_final$V2=Sloc_final$V2+1
Sloc_final$V3=Sloc_final$V3+1

#load file of stations from EFAS
flefas=read.csv(paste0(valid_path,"Stations/efas_flooddriver_match.csv"))
efasmatch=match(flefas$StationID,Sloc_final$V1)
faichier=inner_join(flefas,Sloc_final,by=c("StationID"="V1"))
Sloc_final$csource[efasmatch]="EFAS"
Sloc_final$idlalo=paste(Sloc_final$V3,Sloc_final$V2, sep=" ")

#Now the upstream area
outletname="GIS/upArea_European_01min.nc"
dir=valid_path
UpArea=UpAopen(valid_path,outletname,Sloc_final)
ValidSta=UpArea[which(UpArea$upa>100),]

## Flag EFAS stations that were used for calibration --------
efas_stations_orig=read.csv("Stations/stations_efas_meta.csv",sep=";")
efas_stations_orig$latlong=paste(efas_stations_orig$LisfloodX, efas_stations_orig$LisfloodY, sep=" ")
efas_stations=flefas
efas_stations$latlong=paste(efas_stations$LisfloodX, efas_stations$LisfloodY, sep=" ")
efas_stations2=inner_join(efas_stations, efas_stations_orig,by="latlong")
efas_stations_cal=efas_stations2[which(efas_stations2$EC_calib==1),]
ValidSta$calib=FALSE
matchcal=match(efas_stations_cal$StationID,ValidSta$V1)
ValidSta$calib[matchcal]=TRUE

## Flag stations that have very different mean discharges -------------
fileobs="out/obs_meanAY.csv"
filesim="out/EFAS_meanAY.csv"
obs=read.csv(paste0(valid_path,fileobs))
names(obs)=c("Station_ID","mean")
sim=read.csv(paste0(valid_path,filesim))
names(sim)=c("Station_ID","mean")

years=c(1950:2020)
obstest=obs[which(!is.infinite(obs$mean)),]
simtest=sim[which(!is.infinite(sim$mean)),]
obs_sim=inner_join(obstest,simtest, by=c("Station_ID"))
names(obs_sim)[c(2,3)]=c("obs","sim")

#Remove stations that were removed in the first step
rmv0=which(!is.na(match(obs_sim$Station_ID,ValidSta$V1)))
obs_sim=obs_sim[rmv0,]
rmv1=match(obs_sim$Station_ID,ValidSta$V1)
obs_sim$UpA=ValidSta$upa[rmv1]


## Number of years with data ------------
Q_data <- read.csv(paste0('Q_19502020.csv'), header = F)  # CSVs with observations
Q_data=Q_data[-1,-1]
Station_data_IDs <- as.vector(t(Q_data[1, -1]))
lR=c()
for (id in 1:length(Station_data_IDs)){
  Qs=Q_data[-1,id+1]
  xl=length(Qs)
  l=length(which(!is.na(Qs)))
  lR=c(lR,l)
}
plot(lR[order(lR)])
RecordLen=data.frame(Station_data_IDs,lR)
rmv2=match(obs_sim$Station_ID,RecordLen[,1])
obs_sim$Rlen=RecordLen$lR[rmv2]
dat=obs_sim
rat1=abs((dat$sim-dat$obs))/dat$obs
rat2=abs((dat$obs-dat$sim))/dat$sim
rmv=which(rat1>3 | rat2>3)
dat$rat1=rat1
dat$rat2=rat2
dat$flag1=NA
dat$flag1[which(dat$rat1>3 | dat$rat2>3)]=1
dat$flag1[which(dat$obs>10 & dat$flag1==1)]=2
dat$flag1[which(dat$rat1>6 | dat$rat2>6)]=3

#loading stations' initial database
Q_stations <- st_read(paste0(valid_path,'Europe_daily_combined_2022_v2_WGS84_NUTS.shp')) 

#Writing result: merging dat and validSta
ValidS=inner_join(ValidSta,dat,by=c("V1"="Station_ID"))
ValidS=inner_join(ValidS,Q_stations,by = c("V1"="StationID"))

#reorganise the data
ValidSf=ValidS[,c(1,2,7,19,4,5,6,14,10,11,12,13,15,18,22)]

## Identify stations with low KGE ---------------
kgefile="out/EFAS_obs_kgeAY.csv"
kge=SpatialSkillPlot(ValidSf,"kge",kgefile)
ValidSf=kge[[1]][,-16]
#Isolate problematic stations:
ValidSf$flag2=NA
ValidSf$flag2[which(ValidSf$skill<=(-0.41))]=1
ValidSf$removal = ""
ValidSf$comment = ""

Mancheck_in=ValidSf[which(ValidSf$flag2==1),]
#write.csv(Mancheck_in,file="stations_manual_check.csv")

## Individual check of the remaining stations ------------------
#The manual check is done in excel

#load remaining stations checked
Mancheck=read.csv(file="Stations/stations_manual_check.csv")
Mancheck=Mancheck[,-1]
manmatch=match(Mancheck$V1,ValidSf$V1)
Vcheck=Mancheck[c(25,26)]
colnames(Vcheck)
colnames(ValidSf)[c(19,20)]
ValidSf[manmatch,c(19,20)]=Vcheck

#Other manual Checks
#Gagnieres 6139070
ValidSf$comment[which(ValidSf$V1==6139070)]="wrong river"
ValidSf$removal[which(ValidSf$V1==6139070)]="YES"
#6233203
ValidSf$comment[which(ValidSf$V1==6233203)]="dubious observations"
ValidSf$removal[which(ValidSf$V1==6233203)]="YES"
#6118175
ValidSf$comment[which(ValidSf$V1==6118175)]="wrong location"
ValidSf$removal[which(ValidSf$V1==6118175)]="YES"
#6124440
ValidSf$comment[which(ValidSf$V1==6124440)]="wrong location"
ValidSf$removal[which(ValidSf$V1==6124440)]="YES"
#6340215
ValidSf$comment[which(ValidSf$V1==6340215)]="wrong location"
ValidSf$removal[which(ValidSf$V1==6340215)]="YES"

ValidSf$comment[which(ValidSf$flag1==2 & ValidSf$removal=="")]="mismatch between obs and sim Qmean"
ValidSf$removal[which(ValidSf$flag1==2 & ValidSf$removal=="")]="YES"

ValidSf$comment[which(ValidSf$flag1==3 & ValidSf$removal=="")]="mismatch between obs and sim Qmean"
ValidSf$removal[which(ValidSf$flag1==3 & ValidSf$removal=="")]="YES"

## Distance between "official gauges" and efas points -----------------------
points=ValidSf[,c(1,2,3)]
Vsfloc=st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
Statloc=Q_stations
Sfloc=inner_join(points,Statloc,by=c("V1"="StationID"))
dist=c()
for (r in Vsfloc$V1){
  cat(paste0(r,"\n"))
  v1=Vsfloc[which(Vsfloc$V1==r),]
  v2=Statloc[which(Statloc$StationID==r),]
  oula=st_distance(v1,v2)
  oula=oula/1000
  dist=c(dist,oula)
}

ValidSf$distance=dist
ValidSf$flag2[which(ValidSf$flag2==1 & ValidSf$removal=="")]=2
ValidSf$removal[which(ValidSf$flag2==2)]="YES"
ValidSf$comment[which(ValidSf$flag2==2)]="mismatch between Qobs and Qsim"

length(ValidSf$comment[which(ValidSf$comment=="too far from original station")])
length(ValidSf$comment[which(ValidSf$comment=="mismatch between obs and sim Qmean")])
konar=(ValidSf[which(ValidSf$removal=="YES"),])
length(which(ValidSf$removal=="YES"))
length(ValidSf$comment[which(ValidSf$csource=="EFAS")])

#Writing of the final file
#write.csv(ValidSf,file="Stations/Stations_Validation.csv")
StatCheck=ValidSf[which(ValidSf$flag2==1),]
#write.csv(StatCheck,file="Stations/Stations_lowKGE.csv")







# Part 2: Diagnostic plots------------------------------------

ValidSf=read.csv(file="Stations/Stations_Validation.csv")[,-1]

ValidSY=ValidSf[which(ValidSf$removal!="YES"),]
#ValidSY is the final set of stations used

vEFAS=ValidSY[which(ValidSY$csource=="EFAS"),]


## Figure S2 Upstream area--------------------------------------

### Figure S2.a Map of upstream area--------------------------------------
#Plot parameters
cord.dec=ValidSY[,c(1,2)]
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
palet=c(hcl.colors(9, palette = "viridis", alpha = NULL, rev = T, fixup = TRUE))

ppl <- st_as_sf(ValidSY, coords = c("Var1", "Var2"), crs = 4326)
ppl <- st_transform(ppl, crs = 3035)
ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,fill=Rlen,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                 sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
    breaks=c(365,1825,3650,7300,14600,21900), labels=c(1,5,10,20,40,60)) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5,reverse=F),
         size= guide_legend(override.aes = list(fill = "grey50")))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "bottom",
        legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


ggsave("Plots/ValidStations.jpg", width=20, height=24, units=c("cm"),dpi=1500)

### Figure S2.b Histogram of Upstream area--------------------------------------------------
p<-ggplot(ppl, aes(x=UpA)) + 
  geom_histogram(color="steelblue", fill="slategray1",bins=15,alpha=0.9,lwd=1)+
  scale_y_continuous(breaks=seq(0,600, by=100),name="Number of stations")+
  scale_x_log10(name=expression(paste("Upstream area ", (km^2),sep = " ")),
                breaks=c(100,1000,10000,100000), minor_breaks = log10_minor_break(),
                labels=c("100","1 000","10 000","100 000")) +
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=16),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))+
  annotate("label", x=450000, y=500, label= paste0("n = ",length(ppl$UpA)),size=6)

ggsave("histo_stations_f.jpg", p, width=20, height=15, units=c("cm"),dpi=1500)


## Figure 8: Scatter plot of quantiles -------------
fileobs="out/obs_percentilesAY.csv"
filesim="out/EFAS_percentilesAY.csv"
obs=read.csv(paste0(valid_path,fileobs))
sim=read.csv(paste0(valid_path,filesim))
names(obs)[1] = names(sim)[1]="Station_ID"
names(obs)[c(2:100)] = names(sim)[c(2:100)]=paste('q',1:99,sep="_")

#quantile selection:
qtil=05
myq=paste0("q_",qtil)
if (qtil<10) qtil=paste0(0,qtil)
nobs=paste0("Q",qtil,"_obs")
nsim=paste0("Q",qtil,"_sim")
scplot=paste0("scatterplot_Q",qtil)
qloc=match(myq,colnames(obs))
obsx=data.frame(obs[,c(1,qloc)])
#Only keep valid stations
valid=match(ValidSY$V1, obsx$Station_ID)
obsv=obsx[valid,]
obstest<-melt(obsv,id.vars ="Station_ID")
obsin=obstest[which(obstest$value>0),]

simx=data.frame(sim[,c(1,qloc)])
valid=match(ValidSY$V1, simx$Station_ID)
simv=simx[valid,]
simtest<-melt(simv,id.vars ="Station_ID")
simin=simtest[which(simtest$value>0),]
Myplot2=ScatterValid(obsin,nobs,simin,nsim,scplot,valid_path,ValidSY)
Myplot2
ggsave(paste0("Plots/",scplot,".jpg"), width=20, height=15, units=c("cm"),dpi=1500)

## Figure 6: Spatial skill plots -------------------------

kgefile="out/EFAS_obs_kgeAY.csv"
corrfile="out/EFAS_obs_corr_PearsonAY.csv"
biasfile="out/EFAS_obs_biasAY.csv"
varfile="out/EFAS_obs_variabilityAY.csv"


##Add the skills for logged data
sqval=T
kgefile="out/EFAS_obs_sqrtkgeAY.csv"
corrfile="out/EFAS_obs_sqrtcorr_PearsonAY.csv"
biasfile="out/EFAS_obs_sqrtbiasAY.csv"
varfile="out/EFAS_obs_sqrtvariabilityAY.csv"




kge=SpatialSkillPlot(ValidSY,"kge",kgefile)[[1]]
median(kge$skill)
quantile(kge$skill,c(0.25,0.5,0.75))
corr=SpatialSkillPlot(ValidSY,"r",corrfile)[[1]]
bias=SpatialSkillPlot(ValidSY,"b",biasfile)[[1]]
median(bias$skill)
bg=kge[which(kge$skill<=-0.41),]



length(bg)/2901
plot(bg)
var=SpatialSkillPlot(ValidSY,"v",varfile)[[1]]
mean(var$skill)
length((var$value[which(var$value<1)]))/length(var$value)
kgeout=kge[[1]]
kgeout <- st_as_sf(kgeout, coords = c("Var1", "Var2"), crs = 4326)
kgeout <- st_transform(kgeout, crs = 3035)


Plots=FALSE
if (Plots=TRUE){
  kgep=SpatialSkillPlot(ValidSY,"kge",kgefile)
  kgep[[1]]
  corr=SpatialSkillPlot(ValidSY,"r",corrfile)
  corr[[2]]
  biasp=SpatialSkillPlot(ValidSY,"b",biasfile)
  biasp[[2]]
  varp=SpatialSkillPlot(ValidSY,"v",varfile)
  varp[[2]]

  if (sqval==F){
    ggsave("Plots/Validation_KGE.jpg", kgep.width=20, height=20, units=c("cm"),dpi=1500)
    ggsave("Plots/Validation_corr.jpg",  corr[[2]], width=20, height=20, units=c("cm"),dpi=1500)
    ggsave("Plots/Validation_bias.jpg",  biasp[[2]], width=20, height=20, units=c("cm"),dpi=1500)
    ggsave("Plots/Validation_var.jpg",  varp[[2]], width=20, height=20, units=c("cm"),dpi=1500)
  }
  if (sqval==T){
    ggsave("Plots/Validation_KGE_sqrt.jpg", kgep[[2]],width=20, height=20, units=c("cm"),dpi=1500)
    ggsave("Plots/Validation_corr_sqrt.jpg",  corr[[2]], width=20, height=20, units=c("cm"),dpi=1500)
    ggsave("Plots/Validation_bias_sqrt.jpg",  biasp[[2]], width=20, height=20, units=c("cm"),dpi=1500)
    ggsave("Plots/Validation_var_sqrt.jpg",  varp[[2]],width=20, height=20, units=c("cm"),dpi=1500)
  }


  
  #other values histogram
  length(which(bias$skill>0.8 & bias$skill<1.2))/length(bias$skill)
  length(which(corr$skill>0.5))/length(corr$skill)
  ValidSZ=biasp[[1]]
  name="bias"
  limi=c(0,2)
  
  hist(ValidSZ$skill)
  
  p<-ggplot(ValidSZ, aes(x=skill)) + 
    geom_histogram(fill="aliceblue", color="grey20",alpha=0.9,lwd=1,bins=30) +
    geom_point(x=mean(ValidSZ$skill,na.rm=T),y=0,pch=21,size=4,stroke=2,fill="royalblue")+
    geom_vline(xintercept=1,col=2,lwd=2)+
    scale_y_continuous(name="Number of stations")+
    scale_x_continuous(name=name,limit=limi, breaks=seq(-1,2,by=0.2)) +
    # scale_x_discrete(labels= c("< -0.41","-0.41-0.2", "0.2-0.5","0.5-0.7","0.7-0.8", ">0.8"), name="KGE")+
    theme(axis.title=element_text(size=16, face="bold"),
          axis.text = element_text(size=13),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          legend.position = "none",
          panel.grid.major = element_line(colour = "grey80"),
          panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  p
  ggsave("histo_var_AY.jpg", p, width=20, height=15, units=c("cm"),dpi=1500)
  
  bias=SpatialSkillPlot(ValidSY,"b",biasfile)
  ValidSZ=kgep[[1]]
  min(ValidSZ$skill)
  ptn=empdis(ValidSY$skill,1)
  ptn$y=ptn$emp.f*400
  p<-ggplot(ValidSZ, aes(x=skill)) + 
    geom_histogram(fill="aliceblue", color="grey20",alpha=0.9,lwd=1,bins=30) +
    geom_point(x=mean(ValidSZ$skill,na.rm=T),y=0,pch=21,size=4,stroke=2,fill="royalblue")+
    geom_vline(xintercept=1,col=2,lwd=2)+
    geom_vline(xintercept=-0.41,col="darkgreen",lwd=1.3,lty=2)+
    geom_line(data=ptn, aes(x=Q, y=y),col="red",lwd=1.2)+
    scale_y_continuous(name="Number of stations",    # Add a second axis and specify its features
                       sec.axis = sec_axis( trans=~./400, name="ECDF"))+
    scale_x_continuous(name=name,limit=limi, breaks=seq(-1,1,by=0.2)) +
    # scale_x_discrete(labels= c("< -0.41","-0.41-0.2", "0.2-0.5","0.5-0.7","0.7-0.8", ">0.8"), name="KGE")+
    theme(axis.title=element_text(size=16, face="bold"),
          axis.text = element_text(size=13),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          legend.position = "none",
          panel.grid.major = element_line(colour = "grey80"),
          panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  p
  ggsave("histo_kge_AY.jpg", p, width=20, height=15, units=c("cm"),dpi=1500)
  parpl=kgep[[1]]
  
  length(parpl$skill[which(parpl$skill>=0.5)])/length(parpl$skill)
  median(parpl$skill)
  colorz = c('#d73027','orange','#fee090','lightblue','royalblue',"darkblue")
  kgelabs=c("< -0.41","-0.41 - 0","0 - 0.2", "0.2 - 0.5","0.5 - 0.75", ">0.75")
  #KGE values histogram
  p<-ggplot(parpl, aes(x=factor(kgecode))) + 
    geom_histogram(aes(fill=factor(kgecode)), color="grey20",alpha=0.9,lwd=1,stat="count")+
    scale_fill_manual(values = colorz, labels= kgelabs, name="KGE")   +
    scale_y_continuous(name="Number of stations", breaks = seq(0,1500,by=200),minor_breaks = seq(0,1500,50))+
    # scale_x_continuous(name="KGE",
    #               breaks=seq(-0.4,1,by=0.1)) +
    scale_x_discrete(labels= kgelabs, name="KGE")+
    theme(axis.title=element_text(size=16, face="bold"),
          axis.text = element_text(size=13),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          legend.position = "none",
          panel.grid.major = element_line(colour = "grey80"),
          panel.grid.minor.y = element_line(colour = "grey90",linetype="dashed"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  p
  ggsave("histo_KGE_AYc.jpg", p, width=20, height=15, units=c("cm"),dpi=1500)
  
}

## Figure 7: Performaces by station categories ---------------

### Boxplot of performance by catchment area groups -----------------
#groups= 100-1000 1000-10000 10000-100000 100000-10000000
ValidSY=kgep[[1]]
ValidSY$UpAgroup=1
ValidSY$UpAgroup[which(ValidSY$UpA>200 & ValidSY$UpA<=500)]=2
ValidSY$UpAgroup[which(ValidSY$UpA>500 & ValidSY$UpA<=1000)]=3
ValidSY$UpAgroup[which(ValidSY$UpA>1000 & ValidSY$UpA<=10000)]=4
ValidSY$UpAgroup[which(ValidSY$UpA>10000 & ValidSY$UpA<=100000)]=5
ValidSY$UpAgroup[which(ValidSY$UpA>100000)]=6


agUpA=aggregate(list(upav=ValidSY$UpA),
                by = list(upa=ValidSY$UpAgroup),
                FUN = function(x) c(length(x)))

median(ValidSY$UpA)
meds <- c(by(ValidSY$skill, ValidSY$UpAgroup, median))
q <- c(by(ValidSY$skill, ValidSY$UpAgroup, quantile))
merdecol=match(ValidSY$UpAgroup,names(meds))
ValidSY$col=meds[merdecol]
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))

p1<-ggplot(ValidSY, aes(x=factor(UpAgroup), y=skill)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-0.5,1),name="KGE",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  scale_x_discrete(labels=c("1" = "100-200", "2" = "200-500",
                            "3" = "500-1000","4" = "1000-10 000",
                            "5" = "10 000-100 000","6" = ">100 000"),name="Catchment Area (km2)")+
  scale_fill_gradientn(
    colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


p1
ggsave("Plots/boxplot_KGE.jpg", p1, width=20, height=15, units=c("cm"),dpi=1500)

### Boxplot of performance by decade -----
#groups= "1950-1960","1960-1970","1970-1980","1980-1990","1990-2000","2000-2010","2010-2020"

kgefile="out/EFAS_obs_kge.csv"
skill=read.csv(paste0(valid_path,kgefile))
if (length(skill[1,])>2){
  years=c(1950:2020)
  skill=skill[-73]
  names(skill)[c(2:72)]=years
}
names(skill)[1]="Station_ID"
#Only keep valid stations
skillm<-melt(skill,id.vars ="Station_ID")
skillm$year=as.numeric(as.character(skillm$variable))
skillm$decade=1950
skillm$decade[which(skillm$year>1959 & skillm$year<1970)]=1960
skillm$decade[which(skillm$year>1969 & skillm$year<1980)]=1970
skillm$decade[which(skillm$year>1979 & skillm$year<1990)]=1980
skillm$decade[which(skillm$year>1989 & skillm$year<2000)]=1990
skillm$decade[which(skillm$year>1999 & skillm$year<2010)]=2000
skillm$decade[which(skillm$year>2009)]=2010

#averaging by decade
skillm$value[which(is.infinite(skillm$value))]=NA
dat_skill=aggregate(list(skill=skillm$value),
                      by = list(Station_ID=skillm$Station_ID,decade=skillm$decade),
                      FUN = function(x) c(mean= mean(x,na.rm=T),median=median(x,na.rm=T),max=max(x,na.rm=T)))
dat_skill <- do.call(data.frame, dat_skill)
onlyV=which(!is.na(match(dat_skill$Station_ID,ValidSY$V1)))
dat_skill=dat_skill[onlyV,]


palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))
meds <- c(by(dat_skill$skill.median, dat_skill$decade, median,na.rm=T))
mcol=match(dat_skill$decade,names(meds))
dat_skill$col=meds[mcol]

labelsX=c("1950-1960","1960-1970","1970-1980","1980-1990","1990-2000","2000-2010","2010-2020")

p3<-ggplot(dat_skill, aes(x=factor(decade), y=skill.median)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.1)+
  scale_fill_gradientn(
    colors=palet, n.breaks=6,limits=c(0.4,0.6)) +
  scale_y_continuous(limits = c(-0.5,1),name="KGE",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  scale_x_discrete(name="Decade",labels=labelsX)+
  geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=round(lslo,2)), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position="none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/boxp_KGEtime.jpg", p3, width=20, height=15, units=c("cm"),dpi=1500)



### Boxplot of performances according to the presence or not of reservoirs ------------------------

#load reservoir ratio map
res_path<-("D:/tilloal/Documents/LFRuns_utils/data/reservoirs/")
outletname="res_ratio_European_01min.nc"
dir=res_path
ResData=resOpen(res_path,outletname,ValidSY)  
rsx=ResData[[2]]
length(rsx$res.ratio[which(rsx$res.ratio<0.5)])/length(rsx$res.ratio)
ValidSYp=ResData
ValidSYp$res.group=0
ValidSYp$res.group[which(ValidSYp$res.ratio>=0.5 & ValidSYp$res.ratio<=1)]=1
ValidSYp$res.group[which(ValidSYp$res.ratio>1)]=2

meds <- c(by(ValidSYp$skill, ValidSYp$res.group, median))
len <- c(by(ValidSYp$skill, ValidSYp$res.group, length))
mcol=match(ValidSYp$res.group,names(meds))
ValidSYp$col=meds[mcol]
palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))

p4<-ggplot(ValidSYp, aes(x=factor(res.group), y=skill)) +
  geom_boxplot(notch=F,position=position_dodge(.9),alpha=.8,aes(fill=col),linewidth=0.8,outlier.alpha = 0.4)+
  scale_y_continuous(limits = c(-0.5,1),name="KGE",breaks = seq(-1,1,by=0.5),minor_breaks = seq(-1,1,0.1))+
  scale_x_discrete(labels=c("0" = "c < 0.5 \n (low reservoir impact)", "1" = "0.5 <= c <= 1  \n (medium reservoir impact)",
                            "2" = "c > 1 \n (high reservoir impact)"),name="c ratio (reservoir volume to mean annual streamflow) ")+
  scale_fill_gradientn(
    colors=palet, n.breaks=6,limits=c(0.2,0.8)) +
  geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=len), col='black', size=4,fontface="bold")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=12),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor.y = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

ggsave("Plots/boxp_KGE_reservoirs.jpg", p4, width=20, height=15, units=c("cm"),dpi=1500)

## Figure 9: Performance on annual maxima and minima timings-----------------------
Q_data <- read.csv(paste0('Q_19502020.csv'), header = F)  # CSVs with observations
Q_sim <- read.csv(file="out/EFAS_19502020.csv")
Station_data_IDo <- as.numeric(as.vector(t(Q_data[1, -c(1)])))
Station_data_IDs <- as.vector(t(Q_sim[1, ]))

### 1- Anmax for every station-----------------------------
hit_out=c()
rs_obs=c()
rs_sim=c()
for (id in 1:length(Station_data_IDo)){
  print(id)
  s=Station_data_IDo[id]

  ms=Q_data[-c(1,2,3,4),c(1,id+1)]
  ms$V2=as.numeric(ms$V2)
  ms$time=as.Date(ms$V2-1,origin="0000-01-01")
  ms[,2]=as.numeric(ms[,2])
  ms=data.frame(time=ms$time,data=ms[,2])
  
  #extract peak moments
  sl=which(Station_data_IDs==s)
  mss=data.frame(time=ms$time,data=Q_sim[-1,sl])
  iy=which(is.na(ms$data))
  hit=NA
  if (length(iy)>0){
    ms=ms[-iy,]
    mss=mss[-iy,]
  }
  if (length(ms$data>0)){
    AMAX_obs <- computeAnnualMaxima(ms)
    AMAX_sim <- computeAnnualMaxima(mss)
    if (length(AMAX_sim$annualMaxIndx)==length(AMAX_obs$annualMaxIndx)){
      AMAX_diff=AMAX_obs$annualMaxIndx-AMAX_sim$annualMaxIndx
      hit=length(which(abs(AMAX_diff)<8))/length(AMAX_diff)
    }
    data=data.frame(ms,sim=mss$data)
    names(data)=c("date","Q","Qs")
    res_obs=data.frame(Station=rep(s,length(AMAX_obs$annualMax)),AnMax=AMAX_obs$annualMax,AnMaxDate=AMAX_obs$annualMaxDate)
    res_sim=data.frame(Station=rep(s,length(AMAX_sim$annualMax)),AnMax=AMAX_sim$annualMax,AnMaxDate=AMAX_sim$annualMaxDate)
    # for (ev in 1:length(AMAX_obs$annualMaxIndx)){
    #   #print(ev)
    #   test=StatEpisode(data, date.deb = data$date[AMAX_obs$annualMaxIndx[ev]], win=30, qx.min = 100, mode = "FL",obs=TRUE)
    #   res=rbind(res,test)
    # }
    # save=c(Station_ID=s,rmean=mean(res[,2],na.rm=T),errmean=mean(res[,3],na.rm=T))
    rs_obs=rbind(rs_obs,res_obs)
    rs_sim=rbind(rs_sim,res_sim)
  }else{print("no observations")}
  hit_out=c(hit_out,hit)
}

hit_out=data.frame(hit_out,Station_data_IDo)
hitm=match(ValidSY$V1,hit_out$Station_data_IDo)
hit_out=hit_out[hitm,]
boxplot(hit_out$hit_out)

catagg_obs=aggregate(list(Q=rs_obs$AnMax),
                   by = list(Station=rs_obs$Station),
                   FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
catagg_obs <- do.call(data.frame, catagg_obs)

catagg_sim=aggregate(list(Q=rs_sim$AnMax),
                     by = list(Station=rs_sim$Station),
                     FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
catagg_sim <- do.call(data.frame, catagg_sim)

catagg_obssim=inner_join(catagg_obs,catagg_sim, by="Station",suffix = c("obs","sim"))


#reataining only good matches
matV=match(ValidSY$V1,catagg_obssim$Station)
Catagg_f=catagg_obssim[matV,]
Catagg_f$dPeaks=(Catagg_f$Q.meansim-Catagg_f$Q.meanobs)/Catagg_f$Q.meanobs*100
matO=which(!is.na(match(rs_obs$Station,ValidSY$V1)))
matS=which(!is.na(match(rs_sim$Station,ValidSY$V1)))
rs_obsF=rs_obs[matO,]
rs_simF=rs_sim[matS,]
rs_obsF$yday=yday(rs_obsF$AnMaxDate)

seasonMax_obs=aggregate(list(Q=rs_obsF$yday),
                        by = list(Station=rs_obsF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMax_obs <- do.call(data.frame, seasonMax_obs)

rs_simF$yday=yday(rs_simF$AnMaxDate)
seasonMax_sim=aggregate(list(Q=rs_simF$yday),
                        by = list(Station=rs_simF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMax_sim <- do.call(data.frame, seasonMax_sim)

#difference in seasonality
seasonMax_obs$dseason=seasonMax_obs$Q.season-seasonMax_sim$Q.season
seasonMax_obs$dseason[which(seasonMax_obs$dseason>182)]=seasonMax_obs$dseason[which(seasonMax_obs$dseason>182)]-365.25
seasonMax_obs$dseason[which(seasonMax_obs$dseason<=-182)]=365.25+seasonMax_obs$dseason[which(seasonMax_obs$dseason<=-182)]

#Bonus: Spatial plot, seasonality
palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
points=inner_join(ValidSY,seasonMax_obs,by=c("V1"="Station"))
parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
parpl <- st_transform(parpl, crs = 3035)
cord.dec=points[,c(1,2)]
cord.dec = SpatialPoints(cord.dec, proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:3035"))
nco=cord.UTM@coords
world <- ne_countries(scale = "medium", returnclass = "sf")
Europe <- world[which(world$continent == "Europe"),]
e2=st_transform(Europe,  crs=3035)
w2=st_transform(world,  crs=3035)
tsize=12
osize=12
limi=c(-100,100)
basemap=w2

ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=parpl,aes(geometry=geometry,fill=dseason,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(fill=NA, color="grey20") +
  geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 3), trans="sqrt")+
  scale_fill_gradientn(
    colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name="Qpeak season difference") +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 12,reverse=F),size="none")+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
ggsave("Plots/Spatial_Speakdiff.jpg", width=20, height=15, units=c("cm"),dpi=1500)


### 2- Anmin for every station-----------------------------
hit_outl=c()
rsl_obs=c()
rsl_sim=c()
for (id in 1:length(Station_data_IDo)){
  print(id)
  s=Station_data_IDo[id]
  ms=Q_data[-c(1,2,3,4),c(1,id+1)]
  ms$V2=as.numeric(ms$V2)
  ms$time=as.Date(ms$V2-1,origin="0000-01-01")
  ms[,2]=as.numeric(ms[,2])
  ms=data.frame(time=ms$time,data=ms[,2])
  #extract peak moments
  sl=which(Station_data_IDs==s)
  mss=data.frame(time=ms$time,data=Q_sim[-1,sl])
  iy=which(is.na(ms$data))
  hit=NA
  if (length(iy)>0){
    ms=ms[-iy,]
    mss=mss[-iy,]
  }
  if (length(ms$data>0)){
    AMIN_obs <- computeAnnualMinima(ms)
    AMIN_sim <- computeAnnualMinima(mss)
    if (length(AMIN_sim$annualMinIndx)==length(AMIN_obs$annualMinIndx)){
      AMIN_diff=AMIN_obs$annualMinIndx-AMIN_sim$annualMinIndx
      hit=length(which(abs(AMIN_diff)<15))/length(AMIN_diff)
    }
    data=data.frame(ms,sim=mss$data)
    names(data)=c("date","Q","Qs")
    res_obs=data.frame(Station=rep(s,length(AMIN_obs$annualMin)),AnMin=AMIN_obs$annualMin,AnMinDate=AMIN_obs$annualMinDate)
    res_sim=data.frame(Station=rep(s,length(AMIN_sim$annualMin)),AnMin=AMIN_sim$annualMin,AnMinDate=AMIN_sim$annualMinDate)
    rsl_obs=rbind(rsl_obs,res_obs)
    rsl_sim=rbind(rsl_sim,res_sim)
  }else{print("no observations")}
  hit_outl=c(hit_outl,hit)
}
hit_outl=data.frame(hit_outl,Station_data_IDo)
hitm=match(ValidSY$V1,hit_outl$Station_data_IDo)
hit_outl=hit_outl[hitm,]

cataggl_obs=aggregate(list(Q=rsl_obs$AnMin),
                     by = list(Station=rsl_obs$Station),
                     FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
cataggl_obs <- do.call(data.frame, cataggl_obs)

cataggl_sim=aggregate(list(Q=rsl_sim$AnMin),
                     by = list(Station=rsl_sim$Station),
                     FUN = function(x) c(med=median(x,na.rm=T),mean=mean(x,na.rm=T),dev=sd(x,na.rm=T),len=length(x)))
cataggl_sim <- do.call(data.frame, cataggl_sim)

cataggl_obssim=inner_join(cataggl_obs,cataggl_sim, by="Station",suffix = c("obs","sim"))

#reataining only good matches
matV=match(ValidSY$V1,cataggl_obssim$Station)
Cataggl_f=cataggl_obssim[matV,]
Cataggl_f$dPeaks=(Cataggl_f$Q.meansim-Cataggl_f$Q.meanobs)/Cataggl_f$Q.meanobs*100

matO=which(!is.na(match(rsl_obs$Station,ValidSY$V1)))
matS=which(!is.na(match(rsl_sim$Station,ValidSY$V1)))

rsl_obsF=rsl_obs[matO,]
rsl_simF=rsl_sim[matS,]

rsl_obsF$yday=yday(rsl_obsF$AnMinDate)
seasonMin_obs=aggregate(list(Q=rsl_obsF$yday),
                        by = list(Station=rsl_obsF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMin_obs <- do.call(data.frame, seasonMin_obs)

rsl_simF$yday=yday(rsl_simF$AnMinDate)
seasonMin_sim=aggregate(list(Q=rsl_simF$yday),
                        by = list(Station=rsl_simF$Station),
                        FUN = function(x) c(season=season1(x),len=length(x)))
seasonMin_sim <- do.call(data.frame, seasonMin_sim)

seasonMin_obs$dseason=seasonMin_obs$Q.season-seasonMin_sim$Q.season
seasonMin_obs$dseason[which(seasonMin_obs$dseason>182)]=seasonMin_obs$dseason[which(seasonMin_obs$dseason>182)]-365.25
seasonMin_obs$dseason[which(seasonMin_obs$dseason<=-182)]=365.25+seasonMin_obs$dseason[which(seasonMin_obs$dseason<=-182)]


Perfres=rbind(seasonMax_obs[,c(1,4)],seasonMin_obs[,c(1,4)])
Perfres$group=1
Perfres$group[c(2902:5802)]=2

#Violin plot of timin errors in AnnualMin and AnnualMax
ggplot(Perfres, aes(x=factor(group), y=dseason)) +
  geom_violin(aes(fill=factor(group)),position=position_dodge(.9),alpha=0.8,size=0.9)+
  geom_boxplot(width=0.05,notch=F,size=0.8)+
  coord_flip()+
  scale_y_continuous(limits = c(-60,60),name="Error (days)",breaks = seq(-60,60,by=10),minor_breaks = seq(-100,100,5))+
  scale_x_discrete(labels=c("1" = "Annual maxs", "2" = "Annual mins"),name="")+
  scale_fill_manual(values = c("1"="#9B98F6","2" ="#FF8888"), name = "")+
  theme(axis.title=element_text(size=16),
        axis.text = element_text(size=14),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "none",
        panel.grid.major = element_line(colour = "grey60"),
        panel.grid.minor = element_line(colour = "grey80",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm")) 

ggsave("Plots/Bxplot_seasonF.jpg", width=20, height=7, units=c("cm"),dpi=1500)
