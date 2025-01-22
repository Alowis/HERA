# ==============================================================================
# Title: Creation of the validation metrics and plots of the HERA dataset, linked to the article:  
# HERA: a high-resolution pan-European hydrological reanalysis (1951-2020)
# Author: Alois Tilloy - Joint Research Centre - Unit C6 
# Date: 2024 -02 -01 
# Description:
#   This script allows to generate plots  and assess the skills of the HERA dataset
#   and the mHM run against observed discharge data from 2448 river gauges across Europe
# ==============================================================================

source("~/06_Floodrivers/DataPaper/Code/HERA/functions.R")

main_path = 'D:/tilloal/Documents/06_Floodrivers/'
valid_path = paste0(main_path,'DataPaper/')
dis_path<-paste0(main_path,'dis/calibrated/filtered/Histo/')
setwd(valid_path)

hydroDir<-("D:/tilloal/Documents/06_Floodrivers/DataPaper/")

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
#write.csv(Q_sim,file="out/HERA_Val_19502020.csv")

# Part 1: Create the Valid station file -------------------------------------------

#Loading all txt files with matching locations
#load matching coordinates
Sloc=read.table(paste0(valid_path,"out/yearly/Stations_locations_EFAS_gridv2_1950.txt"),sep=",")
yrlist=c(1951:2020)
for (yi in yrlist){
  print(yi)
  Slocy=read.table(paste0(valid_path,"out/yearly/Stations_locations_EFAS_gridv2_",yi,".txt"),sep=",")
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
ValidSta=UpArea

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

## Load locations of mHM model to see how many matching stations we have -------
mHM_loc=read.table(paste0(valid_path,"Revisions/mHM_EU/mHM_locations_EFAS_grid.txt"),sep=",")

#transalting to R indexing after python
mHM_loc$V2=as.numeric(mHM_loc$V2)+1
mHM_loc$V3=as.numeric(mHM_loc$V3)+1

#load the xls file for catchment area
# Load the readxl package
library(readxl)

# Load the Excel file into a tibble
data <- read_excel(paste0(valid_path,"Revisions/mHM_EU/european_catchments_with_filename.xlsx"))


mHM_sta=full_join(data,mHM_loc, by=c("filename"="V1"))



#load the mHM discharge
mHM_dis=read.table(paste0(valid_path,"Revisions/mHM_EU/Q_mHM_filename.txt"),sep=" ")
mHM_dis=as.data.frame(mHM_dis)
colnames(mHM_dis)=mHM_dis[1,]
mHM_dis=mHM_dis[-1,]
mHM_dis[] <- lapply(mHM_dis, as.numeric)

name_vector=colnames(mHM_dis)
#find the date
date_vector=seq(as.Date("1960-01-01"),as.Date("2010-12-31"),by="days")

#extract mhm_station that are in the dataset

matm=na.omit(match(name_vector,mHM_sta$filename))
mHM_sta=mHM_sta[matm,]

match(mHM_sta$filename,name_vector)


#Plot my stations and mHM stations
### Figure S2.a Map of upstream area--------------------------------------
#Plot parameters

ValidSf=read.csv(file="Stations/Stations_ValidationF.csv")[,-1]
ValidSY=ValidSf[which(ValidSf$removal!="YES"),]

ValidSY$Rlenyr=ValidSY$Rlen/365
ValidSY=ValidSY[-which(ValidSY$Rlenyr<30),]
ValidSY$UpA=ValidSY$upa
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

mHM <- st_as_sf(mHM_sta, coords = c("LON", "LAT"), crs = 4326)
mHM <- st_transform(mHM, crs = 3035)

ggplot(basemap) +
  geom_sf(fill="gray95", color=NA) +
  geom_sf(data=ppl,aes(geometry=geometry,size=UpA),color="transparent",fill="blue",alpha=.9,shape=21,stroke=0)+ 
  geom_sf(data=mHM,aes(geometry=geometry,size=Area_given_km2),color="red",fill="transparent",alpha=.9,shape=1,stroke=0)+ 
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

#match with distance between points and upstream area
## Distance between "official gauges" and efas points -----------------------
hera=ValidSY[,c(1,2,3)]
heraloc=st_as_sf(hera, coords = c("Var1", "Var2"), crs = 4326)
mhml=mHM_sta[,c(3,2,4)]
mhmloc=st_as_sf(mhml, coords = c("LON", "LAT"), crs = 4326)
dist=c()
mlm=c()
for (r in 1:length(heraloc$upa)){
  cat(paste0(r,"\n"))
  v1=heraloc[r,]
  v2=mhmloc
  oula=st_distance(v1,v2)
  oula=oula/1000
  ziz=which.min(oula)
  guez=oula[ziz]
  dist=c(dist,guez)
  mlm=c(mlm,ziz)
}

hist(dist)

ValidSY$dist2mhm=dist
mhm_candidate=mHM_sta[mlm,]

finalcom=data.frame(ValidSY,mhm_candidate)

plot(finalcom$UpA,finalcom$Area_given_km2)
difupa=(finalcom$UpA)/(finalcom$Area_given_km2)

plot(finalcom$dist2mhm[order(finalcom$dist2mhm)], ylim=c(0,10))
hist(difupa,breaks=20000,xlim=c(0,4))
finalcom_upaclean=finalcom[which(difupa<=1.25 & difupa>=0.85 & dist<5),]

un_st=unique(finalcom_upaclean$filename)

#loop to remove double stations
rmf=c()
for (s in 1:length(un_st)){
  st=un_st[s]
  matsta=which(!is.na(match(finalcom_upaclean$filename,st)))
  if (length(matsta)>1){
    print("multi match")
    print(length(matsta))
    print(st)
    torm=which.min(finalcom_upaclean$dist2mhm[matsta])
    rmt=matsta[-torm]
    rmf=c(rmf,rmt)
  }
}

finalcom_upaclean=finalcom_upaclean[-rmf,]

finalcom_check=finalcom[which(difupa<=1.2 & difupa>=1.1 & dist<5),]
plot(finalcom_upaclean$UpA,finalcom_upaclean$Area_given_km2)

#Plot with the matching points
pp_match <- st_as_sf(finalcom_upaclean, coords = c("Var1", "Var2"), crs = 4326)
pp_match <- st_transform(pp_match, crs = 3035)

m1=ggplot(basemap) +
  geom_sf(fill="gray95", color="black") +
  geom_sf(data=ppl,aes(geometry=geometry,size=UpA,shape="H"),color="royalblue",fill="transparent",alpha=.8,stroke=0)+ 
  geom_sf(data=mHM,aes(geometry=geometry,size=Area_given_km2,shape="m"),color="red",fill="transparent",alpha=.9,stroke=0)+ 
  geom_sf(data=pp_match,aes(geometry=geometry,size=Area_given_km2,shape="ma"),color="black",fill="transparent",alpha=.9,stroke=0)+ 
  #geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
    breaks=c(30,40,50,60,70)) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  scale_color_manual(name="Legend:",
                     values = c("royalblue", "red","black"),
                     labels = c("HERA", "mHM", "match")) +
  scale_shape_manual(name="River stations:",
                     values = c("H" = 19, "m" = 1, "ma" = 8),
                     labels = c(paste0("HERA (n=",length(ValidSY$V1),")"), 
                                paste0("mHM (n=",length(mHM_sta$filename),")"), 
                                paste0("match (n=",length(finalcom_upaclean$V1),")"))) +
  # guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5,reverse=F),
         # size= guide_legend(override.aes = list(fill = "grey50")))+
  guides(fill = "none",
         size= "none",
         shape = guide_legend(override.aes = list(size = 3)))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


m1
ggsave("Revisions/map4.pdf", m1,width=30, height=20, units=c("cm"),dpi=1500)

plot(finalcom_upaclean$idlo,finalcom_upaclean$V3)


#Clean and redo that part: base also for the other analysis

#Compare the performances in terms of KGE and others

#load the mHM discharge
mHM_dis=read.table(paste0(valid_path,"Revisions/mHM_EU/Q_mHM_filename.txt"),sep=" ")
mHM_dis=as.data.frame(mHM_dis)
colnames(mHM_dis)=mHM_dis[1,]
mHM_dis=mHM_dis[-1,]
mHM_dis[] <- lapply(mHM_dis, as.numeric)

name_vector=colnames(mHM_dis)
#find the date
date_vector=seq(as.Date("1960-01-01"),as.Date("2010-12-31"),by="days")

length(date_vector)
length(mHM_dis$AT_200089_DMF)

#Reduce the amount of matched table
# vecto=which(!is.na((match(name_vector,finalcom_upaclean$filename))))
# 
# chier=name_vector[match(finalcom_upaclean$filename,name_vector)]
# mhm_dis2=mHM_dis[,vecto]
# name_vector2=name_vector[vecto]
# 
# ass=match(chier, name_vector2)
# 
# match2=match(name_vector2,finalcom_upaclean$filename)

final_dsetmatch=finalcom_upaclean
vecto=which(!is.na((match(name_vector,finalcom_upaclean$filename))))
name_vector2=name_vector[vecto]
mhm_dis2=mHM_dis[,vecto]
#load lisflood discharge
## HERA loading ---------------------
HERA_data<-read.csv(file="out/HERA_Val2_19502020.csv")
HERA_cordata<-read.csv(file="out/HERA_CorStat_19502020.csv")

replax=match(as.numeric(HERA_cordata[1,]),as.numeric(HERA_data[1,]))[-1]
HERA_data[,replax]=HERA_cordata[,-1]

date2=seq(as.Date("1950-01-04"),as.Date("2020-12-31"),by="days")

obs_Q=read.csv(file=paste0(valid_path,"out/Q_19502020.csv"))

Q_data=obs_Q[,-1]
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

lr1=match(final_dsetmatch$V1,RecordLen[,1])
final_dsetmatch$Rlen=RecordLen$lR[lr1]/365.25
plot(final_dsetmatch$Rlen[order(final_dsetmatch$Rlen)])


lr2=match(ValidSta$V1,RecordLen[,1])
ValidSta$Rlen=RecordLen$lR[lr2]/365.25
plot(ValidSta$Rlen[order(ValidSta$Rlen)])

date_vector2=seq(as.Date("1950-01-04"),as.Date("2020-12-31"),by="days")
HERA_Q=HERA_data
HERA_ID=HERA_Q[1,]
Q_ID=obs_Q[1,]
id_match=final_dsetmatch$V1


skills_hera=c()
skills_mhm=c()
skills_heraval=c()
#loop over the matching IDs
for (i in 1:length(id_match)){
  # i=44
  print(i) 
  id_check=id_match[i]
  # id_check=6246160
  id_HERA=match(id_check,HERA_ID)
  id_obs=match(id_check,Q_ID)
  HERA_loc=HERA_Q[-1,id_HERA]
  Q_loc=obs_Q[-c(1:4),id_obs]
  
  #remove year 1950
  # rmv=which(as.integer(format(date_vector2, "%Y"))==1950)
  # HERA_loc=HERA_loc[-rmv]
  # Q_loc=Q_loc[-rmv]
  #find corresponding mHM point
  Selected_river=finalcom_upaclean[match(id_check,finalcom_upaclean$V1),]
  id_mhm=Selected_river$filename
  mhmcol=match(id_mhm,name_vector2)
  mhm_qloc=mhm_dis2[,mhmcol]
  
  #put all q on same temporal interval
  
  mtime=match(date_vector,date_vector2)
  
  HERA_mloc=HERA_loc[mtime]
  Q_mloc=Q_loc[mtime]
  
  mhm_mloc=mhm_qloc
  mhm_mloc[which(is.na(Q_mloc))]=NA
  # plot(Q_mloc,type="l")
  # lines(HERA_mloc,col=6)
  # lines(mhm_mloc,col=3)
  
  #save mhm discharge for the considered locations
  #Now do a loop and compute all relevant statistics
  kge_hera=KGE(HERA_mloc,Q_mloc, na.rm=TRUE, method="2012",out.type="full")
  kge_mhm=KGE(mhm_mloc,Q_mloc, na.rm=TRUE, method="2012",out.type="full")
  
  #manual computation of skill
  # Load necessary libraries
  # library(psych)
  # iz=which(is.na(Q_loc))
  # iy=which(is.na(HERA_loc))
  # # Calculate Pearson correlation
  # rp <- cor.test(HERA_loc[-iz], Q_loc[-iz], method = "pearson")$estimate
  # # Calculate b
  # b <- mean(HERA_loc[-iz]) / mean(Q_loc[-iz])
  # 
  # # Calculate g
  # g <- (sd(HERA_loc[-iz]) / mean(HERA_loc[-iz])) / (sd(Q_loc[-iz]) / mean(Q_loc[-iz]))
  
  kge_hera_val=KGE(HERA_loc,Q_loc, na.rm=TRUE, method="2012",out.type="full")
  
  skills_heraval=rbind(skills_heraval,c(id_check,kge_hera_val$KGE.value,kge_hera_val$KGE.elements))
  
  skills_hera=rbind(skills_hera,c(id_check,kge_hera$KGE.value,kge_hera$KGE.elements))
  skills_mhm=rbind(skills_mhm,c(id_check,kge_mhm$KGE.value,kge_mhm$KGE.elements))
}

median(skills_heraval[,2],na.rm=T)
median(skills_hera[,2],na.rm=T)
median(skills_mhm[,2],na.rm=T)

#load performance estimates from python
kgefile="out/EFAS_obs_kgeAY.csv"
corrfile="out/EFAS_obs_corr_PearsonAY.csv"
biasfile="out/EFAS_obs_biasAY.csv"
varfile="out/EFAS_obs_variabilityAY.csv"


#remove the locations where there must be an issue
#bias greater than 6
mhm_dub=which(skills_mhm[,4]>6 | skills_mhm[,4]<0.17)
hera_dub=which(skills_hera[,4]>6 | skills_hera[,4]<0.17)

all_dub=unique(c(mhm_dub,hera_dub))

skills_mhm_c=skills_mhm[-all_dub,]
skills_hera_c=skills_hera[-all_dub,]

plot(skills_hera_c[,2],skills_mhm_c[,2],xlim=c(-.5,1), ylim=c(-.5,1))
abline(a=0,b=1)

plot(skills_hera_c[,3],skills_mhm_c[,3],xlim=c(0,1), ylim=c(0,1))
abline(a=0,b=1)

plot(skills_hera_c[,4],skills_mhm_c[,4],xlim=c(0,2), ylim=c(0,2))
abline(a=0,b=1)

plot(skills_hera_c[,5],skills_mhm_c[,5],xlim=c(0,2), ylim=c(0,2))
abline(a=0,b=1)


#Now some boxplots


#Only keep valid stations

valid=which(!is.na(match(ValidSY$V1,skills_hera_c[,1])))
skillv=ValidSY[valid,]
skillv=data.frame(skillv,skills_hera_c)
skillv$class="HERA"
#I want as many stations from mhm run
skillm=ValidSY[valid,]
valid2=na.omit((match(skillm$V1,skills_mhm_c[,1])))
skillm=data.frame(skillm,skills_mhm_c[valid2,])

skillm$class="mHM"


#redo the matching map


#add river network
outletname="efas_rnet_100km_01min"

outriver=outletopen(hydroDir,outletname)
outR <- st_as_sf(outriver, coords = c("Var1", "Var2"), crs = 4326)
outR <- st_transform(outR, crs = 3035)


#Plot with the matching points
pp_match <- st_as_sf(skillm, coords = c("Var1", "Var2"), crs = 4326)
pp_match <- st_transform(pp_match, crs = 3035)
tsize=20
osize=18
m1=ggplot(basemap) +
  geom_sf(fill="gray95", color="black") +
  geom_sf(data=outR,aes(geometry=geometry),color="gray20",alpha=.9,shape=15,stroke=0, size=0.1)+ 
  geom_sf(data=ppl,aes(geometry=geometry,size=UpA,shape="H"),color="royalblue",fill="transparent",alpha=.6,stroke=0.1)+ 
  geom_sf(data=mHM,aes(geometry=geometry,size=Area_given_km2,shape="m"),color="red",fill="transparent",alpha=.9,stroke=0.2)+ 
  geom_sf(data=pp_match,aes(geometry=geometry,size=upa,shape="ma"),color="black",fill="transparent",alpha=.9,stroke=0.2)+ 
  #geom_sf(fill=NA, color="grey20") +
  scale_x_continuous(breaks=seq(-30,40, by=5)) +
  scale_size(range = c(1, 4), trans="sqrt",name= expression(paste("Upstream area ", (km^2),
                                                                  sep = " ")),
             breaks=c(101,1000,10000,100000,500000), labels=c("100","1000", "10 000", "100 000", "500 000"))+
  scale_fill_gradientn(
    colors=palet,oob = scales::squish, name="record length (years)", trans="sqrt",
    breaks=c(30,40,50,60,70)) +
  coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
  labs(x="Longitude", y = "Latitude") +
  scale_color_manual(name="Legend:",
                     values = c("royalblue", "red","black"),
                     labels = c("HERA", "mHM", "match")) +
  scale_shape_manual(name="River stations:",
                     values = c("H" = 19, "m" = 1, "ma" = 8),
                     labels = c(paste0("HERA (n=",length(ValidSY$V1),")"), 
                                paste0("mHM (n=",length(mHM_sta$filename),")"), 
                                paste0("match (n=",length(skillm$V1),")"))) +
  # guides(fill = guide_colourbar(barwidth = 20, barheight = 0.5,reverse=F),
  # size= guide_legend(override.aes = list(fill = "grey50")))+
  guides(fill = "none",
         size= "none",
         shape = guide_legend(override.aes = list(size = 6, stroke=0.8)))+
  theme(axis.title=element_text(size=tsize),
        panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=tsize),
        legend.text = element_text(size=osize),
        legend.position = "right",
        legend.box = "vertical",
        panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
        panel.grid.minor = element_line(colour = "grey90"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))


m1
ggsave("Revisions/mapCF.jpg", m1,width=30, height=20, units=c("cm"),dpi=500)


median(skillv$V2)
median(skillm$V2)
#Now 2 boxplots and goodbye

points2=rbind(skillv,skillm)
# points2$skill[c(1:length(points$Var1))]=points$KGE
# points2$Rungroup=2
# points2$Rungroup[c(1:length(points$Var1))]=1


if (length(which(is.na(points2$V2)))>0){
  points2=points2[-which(is.na(points2$V2)),]
}
meds <- c(by(points2$V2, points2$class, median))
means <- c(by(points2$V2, points2$class, mean))
sds <- c(by(points2$V2, points2$class, sd))



merdecol=match(points2$class,names(meds))
points2$col=meds[merdecol]

palet=c(hcl.colors(9, palette = "Blues", alpha = NULL, rev = T, fixup = TRUE))
colorz=c("skyblue","purple")
p1<-ggplot(points2, aes(x=factor(class), y=V2)) +
  geom_violin(width=1,position=position_dodge(.9),alpha=.8,aes(fill=factor(col)),linewidth=0.8)+
  geom_boxplot(width=0.1,notch=F,color="black",outlier.alpha = 0.6,alpha=0.6,fill=NA,linewidth=1,coef=1.5)+
  #stat_boxplot(geom ='errorbar') +
  scale_y_continuous(limits = c(-0.5,1),name="KGE'",breaks = seq(-1,1,by=0.2),minor_breaks = seq(-1,1,0.1))+
  scale_x_discrete(labels=c("1" = "EFAS", "2" = "HERA"),name="Dataset")+
  scale_fill_manual(values = colorz) +
  ggtitle (paste0("Comparison of HERA and mHM run \nfor the period 1960-2010 over n = ",length(points2$Var1)/2," European stations"))+
  # scale_fill_gradientn(
  #   colors=palet, n.breaks=6,limits=c(0.4,0.8)) +
  #geom_text(data=data.frame(), aes(x=names(meds), y=meds-0.05, label=agUpA$upav), col='black', size=4,fontface="bold")+
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


ggsave("Revisions/boxplot_KGE_compF.jpg", p1, width=15, height=20, units=c("cm"),dpi=500)




#retain only stations that were used in EFAS calibration

ValidFeat=skillv[,c(1:8,15,25,28,29,30,31)]

smhm=match(skillm$V1,ValidFeat$V1)

ValidFeat=cbind(ValidFeat,skillm[smhm,c(28,29,30,31)])
names(ValidFeat)[c(11:14)]=c("KGEHERA","rHERA","bHERA","vHERA")
names(ValidFeat)[c(15:18)]=c("KGEmHM","rmHM","bmHM","vmHM")

ValidCalib=ValidFeat[which(ValidFeat$calib==TRUE),]

write.csv(ValidCalib,file="out/HERAQvsmHM_CalT.csv", row.names = FALSE)


#now extract these stations from HERAVAL and mhmQ

id_HERA=match(ValidCalib$V1,HERA_ID)
id_obs=match(ValidCalib$V1,Q_ID)


Selected_river=finalcom_upaclean[match(ValidCalib$V1,finalcom_upaclean$V1),]
id_mhm=Selected_river$filename
mhmcol=match(id_mhm,name_vector2)

HERA_St=HERA_Q[-1,id_HERA]
Q_St=obs_Q[-c(1:4),id_obs]
mhm_St=mhm_dis2[,mhmcol]

#Save the files 



#add time dimension
HERA_St=data.frame(date_vector2,HERA_St)
names(HERA_St)[1]="date"
names(HERA_St)[c(2:172)]=ValidCalib$V1

Q_St=data.frame(date_vector2,Q_St)
names(Q_St)[1]="date"
names(Q_St)[c(2:172)]=ValidCalib$V1

write.csv(HERA_St,file="out/HERAQ_calTRUE.csv", row.names = FALSE)
write.csv(Q_St,file="out/obsQ_calTRUE.csv", row.names = FALSE)

mhm_St=data.frame(date_vector,mhm_St)
names(mhm_St)[1]="date"
names(mhm_St)[c(2:172)]=ValidCalib$V1
write.csv(mhm_St,file="out/mHMQ_calTRUE.csv", row.names = FALSE)


hvh=match(skillm[,3],finalcom_upaclean$V1)
hera_vs_mhm=finalcom_upaclean[hvh,]
hera_vs_mhm$RlenYr=hera_vs_mhm$Rlen/366



points=skillv
skillm=as.data.frame(skillm)
metric="KGE"
#For correlation
if (metric=="r"){
  palet=c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
  limi=c(-1,1)
  var="r"
  
  parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
  parpl <- st_transform(parpl, crs = 3035)
  parpef=parpl[which(parpl$calib==TRUE),]
  tsize=size=12
  palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
  limi=c(-0.25,0.25)
  metric="r (HERA) - r (mHM)"
  
  points$rdiff=points$r-skillm$r
  parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
  parpl <- st_transform(parpl, crs = 3035)
  map3=ggplot(basemap) +
    geom_sf(fill="gray95", color=NA) +
    geom_sf(data=parpl,aes(geometry=geometry,fill=rdiff,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
    geom_sf(fill=NA, color="grey20") +
    geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
    scale_x_continuous(breaks=seq(-30,40, by=5)) +
    scale_size(range = c(1, 3), trans="sqrt",guide = 'none')+
    scale_fill_gradientn(
      colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name=metric) +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    ggtitle (paste0("Comparison of HERA and mHM run \nfor the period 1960-2010 over n = ",length(points2$Var1)/2," European stations"))+
    labs(x="Longitude", y = "Latitude") +
    guides(fill = guide_colourbar(barwidth = 24, barheight = 0.5,reverse=F, title.position = "top"))+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize,hjust = 0.5, vjust = 1),
          legend.text = element_text(size=osize),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  map3
  
  ggsave("Revisions/r_mHMvsHERA_final.jpg", map3, width=15, height=20, units=c("cm"),dpi=500)
  
}
#For Bias
if (metric=="b"){
  palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  limi=c(0,2)
  var="b"
}
#For Variability
if (metric=="v"){
  palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  #palet = c('#f61f0f','#f4cccc','#f6ebeb','#ffffff','#e4edf5','#91bfdb','#1c87d8')
  limi=c(0,2)
  var="v"
}
if (metric=="kge"){
  colorz = c('#d73027','orange','#fee090','lightblue','royalblue',"darkblue")
  kgelabs=c("< -0.41","-0.41 - 0","0 - 0.2", "0.2 - 0.5","0.5 - 0.75", ">0.75")
  points$skill=as.numeric(points$skill)
  points$kgecode=0
  points$kgecode[which(is.na(points$skill))]=NA
  points$kgecode[which(points$skill>(-0.41) & points$skill<=0)]=1
  points$kgecode[which(points$skill>0 & points$skill<=0.2)]=2
  points$kgecode[which(points$skill>0.2 & points$skill<=0.5)]=3
  points$kgecode[which(points$skill>0.5 & points$skill<=0.75)]=4
  points$kgecode[which(points$skill>0.75)]=5
  
  points$KGE=skillm$V2
  points$kgecode2=0
  points$kgecode2[which(is.na(points$KGE))]=NA
  points$kgecode2[which(points$KGE>(-0.41) & points$KGE<=0)]=1
  points$kgecode2[which(points$KGE>0 & points$KGE<=0.2)]=2
  points$kgecode2[which(points$KGE>0.2 & points$KGE<=0.5)]=3
  points$kgecode2[which(points$KGE>0.5 & points$KGE<=0.75)]=4
  points$kgecode2[which(points$KGE>0.75)]=5
  var="kge"
  
  parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
  parpl <- st_transform(parpl, crs = 3035)
  parpef=parpl[which(parpl$calib==TRUE),]
  tsize=size=12
  map1=ggplot(basemap) +
    geom_sf(fill="gray95", color=NA) +
    geom_sf(data=parpl,aes(geometry=geometry,colour=factor(kgecode),size=UpA),alpha=.9,stroke=0,shape=16)+ 
    geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.05,shape=1)+ 
    #geom_sf(data=parpef,aes(geometry=geometry),col="black",alpha=1,size=1.4,stroke=0.1,shape=8)+ 
    geom_sf(fill=NA, color="grey20") +
    scale_x_continuous(breaks=seq(-30,40, by=5)) +
    scale_size(range = c(1, 3), trans="sqrt")+
    labs(x="Longitude", y = "Latitude")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_manual(values = colorz, labels= kgelabs, name="KGE")   +
    labs(x="Longitude", y = "Latitude") +
    guides(colour = guide_legend(override.aes = list(size = 4)),size="none")+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  map1
  
  map2=ggplot(basemap) +
    geom_sf(fill="gray95", color=NA) +
    geom_sf(data=parpl,aes(geometry=geometry,colour=factor(kgecode2),size=UpA),alpha=.9,stroke=0,shape=16)+ 
    geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.05,shape=1)+ 
    geom_sf(fill=NA, color="grey20") +
    scale_x_continuous(breaks=seq(-30,40, by=5)) +
    scale_size(range = c(1, 3), trans="sqrt")+
    labs(x="Longitude", y = "Latitude")+
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    scale_color_manual(values = colorz, labels= kgelabs, name="KGE")   +
    labs(x="Longitude", y = "Latitude") +
    guides(colour = guide_legend(override.aes = list(size = 4)),size="none")+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize),
          legend.text = element_text(size=osize),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  map2
  
  
  palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
  limi=c(-0.2,0.2)
  metric="KGE'(HERA) - KGE'(mHM)"
  
  points$kgediff=points$skill-points$KGE
  parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
  parpl <- st_transform(parpl, crs = 3035)
  map3=ggplot(basemap) +
    geom_sf(fill="gray95", color=NA) +
    geom_sf(data=parpl,aes(geometry=geometry,fill=kgediff,size=UpA),color="transparent",alpha=.9,shape=21,stroke=0)+ 
    geom_sf(fill=NA, color="grey20") +
    geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="black",alpha=1,stroke=0.1,shape=1)+ 
    scale_x_continuous(breaks=seq(-30,40, by=5)) +
    scale_size(range = c(1, 3), trans="sqrt",guide = 'none')+
    scale_fill_gradientn(
      colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name=metric) +
    coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
    ggtitle (paste0("Comparison of HERA and mHM run \nfor the period 1960-2010 over n = ",length(points2$Var1)/2," European stations"))+
    labs(x="Longitude", y = "Latitude") +
    guides(fill = guide_colourbar(barwidth = 24, barheight = 0.5,reverse=F, title.position = "top"))+
    theme(axis.title=element_text(size=tsize),
          panel.background = element_rect(fill = "aliceblue", colour = "grey1"),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=tsize,hjust = 0.5, vjust = 1),
          legend.text = element_text(size=osize),
          legend.position = "bottom",
          panel.grid.major = element_line(colour = "grey85",linetype="dashed"),
          panel.grid.minor = element_line(colour = "grey90"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  map3
  
  ggsave("Revisions/KGE_mHMvsHERA_finalF.jpg", map3, width=15, height=20, units=c("cm"),dpi=1500)
  
}

#Do a histogram with orderer difference
#points$skill=points$r
points$skill=as.numeric(points$skill)
points$KGE=skillm$V2
points$kgediff=points$skill-points$KGE
Recap=points

Recap=Recap[order(Recap$kgediff),]
Recap$ord=c(1:length(Recap$Var1))
plot(Recap$skill[Recap$ord])
length(which(Recap$kgediff>0))/length(Recap$Var2)
# Recap <- Recap %>% 
#   mutate(mycolor = ifelse(kgediff>0, "type1", "type2"))
# 
colors <- c(hcl.colors(20, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))
#   
px=ggplot(Recap) +
  geom_segment( aes(y=ord, yend=ord, x=0, xend=kgediff, color=skill), size=1.3, alpha=0.9) +
  scale_color_gradientn(
    colors=colors,n.breaks=5,oob = scales::squish, limits=c(-0.5,1),name="KGE'(HERA)")+
  scale_x_continuous(breaks=seq(-3,3, by=.2),name="KGE'(HERA) - KGE'(mHM)") +
  scale_y_continuous(name="stations") +
  coord_cartesian(xlim=c(-1,1))+
  guides(color = guide_colourbar(barwidth = 12, barheight = 0.5,reverse=F, title.position = "top")) +
  theme(axis.title=element_text(size=14, face="bold"),
        axis.text = element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
px
ggsave("Revisions/diff_KGE_hera-mhmF.jpg", px, width=10, height=20, units=c("cm"),dpi=1500)
dev.off()



metric="b"
#same plot for correlation. bias variability
if (metric=="r"){
  palet=c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
  limi=c(-0.5,0.5)
  var="r"
  points$vmhm = skillm$r
  points$vdiff = points$r - points$vmhm
  Recap = points
  Recap$vcol=r
  Recap=Recap[order(Recap$vdiff),]
  Recap$ord=c(1:length(Recap$Var1))
  nplot="r(HERA) - r(mHM)"
}
#For Variability
if (metric=="v"){
  palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  #palet = c('#f61f0f','#f4cccc','#f6ebeb','#ffffff','#e4edf5','#91bfdb','#1c87d8')
  limi=c(-1,1)
  points$vmhm = abs(skillm$Gamma-1)
  points$vdiff = abs(points$Gamma-1) - points$vmhm
  Recap = points
  Recap$vcol=r
  Recap=Recap[order(-Recap$vdiff),]
  Recap$ord=c(1:length(Recap$Var1))
  var="v"
  nplot="abs(v(HERA)) - abs(v(mHM))"
}
if (metric=="b"){
  palet=c(hcl.colors(9, palette = "RdBu", alpha = NULL, rev = F, fixup = TRUE))
  limi=c(-1,1)
  hist(points$Beta-1)
  points$vmhm = abs(skillm$Beta-1)
  points$vdiff = abs(points$Beta-1) - points$vmhm
  Recap = points
  Recap$vcol=r
  Recap=Recap[order(-Recap$vdiff),]
  Recap$ord=c(1:length(Recap$Var1))
  var="b"
  nplot="abs(b(HERA)) - abs(b(mHM))"
}


#palet=c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
palet=c(hcl.colors(9, palette = "RdYlBu", alpha = NULL, rev = F, fixup = TRUE))

pl=ggplot(Recap) +
  geom_segment( aes(y=ord, yend=ord, x=0, xend=vdiff, color=skill), size=1.3, alpha=0.9) +
  scale_color_gradientn(
    colors=palet,n.breaks=5,oob = scales::squish, limits=c(-0.5,1),name="KGE'(HERA)")+
  scale_x_continuous(breaks=seq(-1,1, by=.2),name=nplot) +
  scale_y_continuous(name="stations") +
  coord_cartesian(xlim=limi)+
  guides(color = guide_colourbar(barwidth = 12, barheight = 0.5,reverse=F, title.position = "top")) +
  theme(axis.title=element_text(size=14, face="bold"),
        axis.text = element_text(size=10),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))
pl

ggsave(paste0("Revisions/diff_",var,"_hera-mhmF.jpg"), pl, width=10, height=20, units=c("cm"),dpi=1500)

plot(points$UpA,points$kgediff, log="x")
#deeper analysis about where most performance difference is found

plot(hera_vs_mhm$upa,hera_vs_mhm$Area_given_km2)
hera_vs_mhm=inner_join(hera_vs_mhm,points, by=c("V1"))

length(which(hera_vs_mhm$kgediff<0))
boxplot(hera_vs_mhm$kgediff, ylim=c(-1,1))
#plot ecdf of kge and other variables for both
ecdf_hera=empdis(points$skill,1)
plot(ecdf_hera$Q,ecdf_hera$emp.f)

ecdf_mhm=empdis(points$KGE,1)
points(ecdf_mhm$Q,ecdf_mhm$emp.f, col=2)

ecdf_hera$d="HERA"
ecdf_mhm$d="mHM"
ecdf_all=rbind(ecdf_hera,ecdf_mhm)
#ecdf_all=data.frame(kge_h=ecdf_hera$Q,kge_m=ecdf_mhm$Q,f_h=ecdf_hera$emp.f,f_m=ecdf_mhm$emp.f)
p<-ggplot(ecdf_all) + 
  geom_point(aes(x=Q, y=emp.f, group=d,col=d),size=2) +
  # geom_point(aes(x=kge_m, y=f_m),col="darkgrey",size=2) +
  # geom_point(x=median(points$skill,na.rm=T),y=0,pch=21,size=4,stroke=2,fill="red")+
  # geom_point(x=median(points$KGE,na.rm=T),y=0,pch=21,size=4,stroke=2,fill="darkgrey")+
  geom_vline(xintercept=1,col=2,lwd=2)+
  geom_vline(xintercept=-0.41,col="darkgreen",lwd=2)+
  geom_vline(xintercept=median(points$skill,na.rm=T),col="red",lwd=1.4, lty=2)+
  geom_vline(xintercept=median(points$KGE,na.rm=T),col="darkgrey",lwd=1.4, lty=2)+
  scale_color_manual(name= "Run",values=c("red","darkgrey"))+
  scale_y_continuous(name="ECDF",breaks = c(0,0.2,.4,.6,.8,1))+
  scale_x_continuous(name="KGE'",limit=c(-0.5,1), breaks=seq(-1,1,by=0.2)) +
  guides(colour = guide_legend(override.aes = list(size = 6)))+
  # scale_x_discrete(labels= c("< -0.41","-0.41-0.2", "0.2-0.5","0.5-0.7","0.7-0.8", ">0.8"), name="KGE")+
  theme(axis.title=element_text(size=18, face="bold"),
        axis.text = element_text(size=16),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

p
ggsave("Revisions/KGE_comp_HERA-mHMF.jpg", p, width=20, height=15, units=c("cm"),dpi=500)
mean(points$KGE)
mean(points$skill)

#same plot for small catchments
pointsmall=points[which(points$UpA<500),]
ecdf_hera=empdis(pointsmall$skill,1)
plot(ecdf_hera$Q,ecdf_hera$emp.f)

ecdf_mhm=empdis(pointsmall$KGE,1)
points(ecdf_mhm$Q,ecdf_mhm$emp.f, col=2)

ecdf_hera$d="HERA"
ecdf_mhm$d="mHM"
ecdf_all=rbind(ecdf_hera,ecdf_mhm)
#ecdf_all=data.frame(kge_h=ecdf_hera$Q,kge_m=ecdf_mhm$Q,f_h=ecdf_hera$emp.f,f_m=ecdf_mhm$emp.f)
ggplot(ecdf_all) + 
  geom_point(aes(x=Q, y=emp.f, group=d,col=d),size=2) +
  # geom_point(aes(x=kge_m, y=f_m),col="darkgrey",size=2) +
  # geom_point(x=median(points$skill,na.rm=T),y=0,pch=21,size=4,stroke=2,fill="red")+
  # geom_point(x=median(points$KGE,na.rm=T),y=0,pch=21,size=4,stroke=2,fill="darkgrey")+
  geom_vline(xintercept=1,col=2,lwd=2)+
  geom_vline(xintercept=-0.41,col="darkgreen",lwd=2)+
  geom_vline(xintercept=median(pointsmall$skill,na.rm=T),col="red",lwd=1.4, lty=2)+
  geom_vline(xintercept=median(pointsmall$KGE,na.rm=T),col="darkgrey",lwd=1.4, lty=2)+
  scale_color_manual(name= "Model",values=c("red","darkgrey"))+
  scale_y_continuous(name="ECDF",breaks = c(0,0.2,.4,.6,.8,1))+
  scale_x_continuous(name="KGE'",limit=c(-0.5,1), breaks=seq(-1,1,by=0.2)) +
  # scale_x_discrete(labels= c("< -0.41","-0.41-0.2", "0.2-0.5","0.5-0.7","0.7-0.8", ">0.8"), name="KGE")+
  theme(axis.title=element_text(size=16, face="bold"),
        axis.text = element_text(size=13),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid = element_blank(),
        panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey80"),
        panel.grid.minor.x = element_line(colour = "grey90",linetype="dashed"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"),
        legend.key.size = unit(.8, "cm"))

