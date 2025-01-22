# ==============================================================================
# Title: Functions to run the validation and usage note linked to the article:  
# HERA: a high-resolution pan-European hydrological reanalysis (1951-2020)
# Author: Alois Tilloy - Joint Research Centre - Unit C6 
# Date: 2024 -02 -01 
# Description:
#   This script countains the functions used in Validation_Reanalysis.R and Usage_note.R
# ==============================================================================


# Library importation ----------------------------------------------------------
suppressWarnings(suppressMessages(library(ncdf4)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(rnaturalearth)))
# suppressWarnings(suppressMessages(library(rnaturalearthdata)))
#suppressWarnings(suppressMessages(library(rgeos)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(lubridate)))
#suppressWarnings(suppressMessages(library(rgdal)))
suppressWarnings(suppressMessages(library(sp)))
suppressWarnings(suppressMessages(library(reshape2)))
suppressWarnings(suppressMessages(library(MASS)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(viridis)))
suppressWarnings(suppressMessages(library(hydroGOF)))



# Input opening ----------------------------------------------------------------

#' Extract Upstream Area Data from NetCDF
#'
#' @description
#' The `UpAopen` function opens a NetCDF file containing upstream area data and
#' extracts the data for a specified outlet name. The function retrieves data
#' for the entire grid defined by the longitude and latitude variables and
#' then joins this data with a final location data frame using an identifier.
#'
#' @param dir A character string specifying the directory path where the NetCDF
#'   file is located.
#' @param outletname A character string specifying the name of the outlet for which
#'   the upstream area data is to be extracted.
#' @param Sloc_final A data frame with the final location information, including
#'   an 'idlalo' column that contains identifiers used to join with the extracted
#'   upstream area data.
#'
#' @return A data frame resulting from an inner join between the extracted upstream
#'   area data and the final location data frame. The resulting data frame includes
#'   longitude, latitude, upstream area (upa), and additional location identifiers.
#'
#' @examples
#' # Assuming 'dir' is the directory path, 'outletname' is the outlet name, and
#' # 'Sloc_final' is a data frame with location information:
#' result <- UpAopen(dir, outletname, Sloc_final)
#' # The result is a data frame with upstream area data joined with final location data.
#'
#' @export
#' @importFrom ncdf4 nc_open ncvar_get
#' @importFrom dplyr inner_join
#'
#' @seealso
#' `ncdf4::nc_open`, `ncdf4::ncvar_get`, and `dplyr::inner_join` for more
#' information on opening NetCDF files, extracting data, and joining data frames.
UpAopen=function(dir,outletname,Sloc_final){
  ncbassin=paste0(dir,outletname)
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[2]
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
  outlets=as.vector(outlets)/1000000
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$upa=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outfinal=inner_join(outll, Sloc_final, by="idlalo")
  return (outfinal)
}

#' Extract Resource Ratio Data from NetCDF
#'
#' @description
#' The `resOpen` function opens a NetCDF file and extracts resource ratio data
#' for a specified outlet name. It retrieves the data for the entire grid
#' defined by longitude and latitude variables, and then processes this data
#' to match and update a validation dataset with resource ratio values.
#'
#' @param dir A character string specifying the directory path where the NetCDF
#'   file is located.
#' @param outletname A character string specifying the name of the outlet for which
#'   the resource ratio data is to be extracted.
#' @param ValidSY A data frame containing validation identifiers and other information
#'   with an 'idlalo' column to be used for matching and updating with the extracted data.
#'
#' @return A data frame that is an updated version of the 'ValidSY' input data frame,
#'   now including 'res.ratio' values matched to each validation identifier.
#'
#' @examples
#' # Assuming 'dir' is the directory path, 'outletname' is the outlet name, and
#' # 'ValidSY' is a data frame with validation information:
#' result <- resOpen(dir, outletname, ValidSY)
#' # The result is the 'ValidSY' data frame updated with resource ratio data.
#'
#' @export
#' @importFrom ncdf4 nc_open ncvar_get
#' @importFrom base match
#'
#' @seealso
#' `ncdf4::nc_open`, `ncdf4::ncvar_get` for opening NetCDF files and extracting data,
#' `base::match` for matching elements in data frames.
#'
resOpen=function(dir,outletname,ValidSY){
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
  lonlatloop=expand.grid(c(1:llo),c(lla:1))
  outll$res.ratio=outlets
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  outll$idlalo=paste(outll$idlo,outll$idla,sep=" ")
  outlp=outll[which(!is.na(outll$res.ratio)),]
  mf=match(ValidSY$idlalo,outlp$idlalo)
  ValidSY$res.ratio=outlp$res.ratio[mf]
  return (ValidSY)
}

#' Extract Outlet Data from NetCDF File
#'
#' @description
#' The `outletopen` function opens a NetCDF file and extracts outlet data based on
#' specified spatial indices. Optionally, you can define a spatial subset of the
#' data to extract using grid indices. If no grid indices are provided, the function
#' extracts all data.
#'
#' @param dir A character string specifying the directory path where the NetCDF
#'   file is located.
#' @param outletname A character string specifying the name of the outlet (without
#'   the file extension) for which data is to be extracted from the NetCDF file.
#' @param nrspace A numeric vector of length 5 providing grid indices to define
#'   a spatial subset for extraction, consisting of start and count indices for
#'   longitude and latitude (in the order of start longitude, count longitude,
#'   start latitude, count latitude). If not provided, the entire spatial extent
#'   is used.
#'
#' @return A data frame containing the outlet data ('outlets') and the corresponding
#'   longitude ('idlo') and latitude ('idla') grid indices.
#'
#' @examples
#' # Assuming 'directory' is the path and 'outlet_name' is the name of the NetCDF file:
#' outlet_data <- outletopen(directory, outlet_name, nrspace = c(NA, 1, 10, 1, 10))
#' # The result is a data frame with outlet data for the specified spatial subset.
#'
#' @export
#' @importFrom ncdf4 nc_open ncvar_get
#'
#' @seealso
#' `ncdf4::nc_open` and `ncdf4::ncvar_get` for functions used to access and
#' extract data from NetCDF files.
#'
outletopen=function(dir,outletname,nrspace=rep(NA,5)){
  ncbassin=paste0(dir,"/",outletname,".nc")
  ncb=nc_open(ncbassin)
  name.vb=names(ncb[['var']])
  namev=name.vb[1]
  #time <- ncvar_get(ncb,"time")
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  if (!is.na(nrspace[1])){
    start=as.numeric(nrspace[c(2,4)])
    count=as.numeric(nrspace[c(3,5)])-start+1
  }else{
    londat = ncvar_get(ncb,name.lon) 
    llo=length(londat)
    latdat = ncvar_get(ncb,name.lat)
    lla=length(latdat)
    start=c(1,1)
    count=c(llo,lla)
  }
  
  londat = ncvar_get(ncb,name.lon,start=start[1],count=count[1]) 
  llo=length(londat)
  latdat = ncvar_get(ncb,name.lat,start=start[2],count=count[2])
  lla=length(latdat)
  outlets = ncvar_get(ncb,namev,start = start, count= count) 
  outlets=as.vector(outlets)
  outll=expand.grid(londat,latdat)
  lonlatloop=expand.grid(c(1:llo),c(1:lla))
  outll$idlo=lonlatloop$Var1
  outll$idla=lonlatloop$Var2
  
  mierda=which(!is.na(outlets))
  outll=outll[mierda,]
  outlets=outlets[which(!is.na(outlets))]
  outll=data.frame(outlets,outll)
  return (outll)
}

#' Extract Discharge Data for Multiple Locations from a NetCDF File
#'
#' @description
#' The `disNcopen` function opens a NetCDF file and extracts discharge data
#' for multiple specified locations. The function loops through each location,
#' retrieves the discharge time series, and compiles the data into a single
#' data frame with identifiers for each location.
#'
#' @param fname A character string specifying the file name (without the extension)
#'   of the NetCDF file containing the discharge data.
#' @param dir A character string specifying the directory path where the NetCDF
#'   file is located.
#' @param outloc A data frame containing location information, including grid indices
#'   for latitude ('idla') and longitude ('idlo'), and an identifier for each location.
#'
#' @return A data frame with columns for discharge values, location identifiers,
#'   longitude, latitude, and time. Each row corresponds to a time series entry
#'   for a specific location.
#'
#' @examples
#' # Assuming the NetCDF file 'discharge_data.nc' is in the directory 'data/',
#' # and 'locations' is a data frame with grid indices and location identifiers:
#' all_discharges <- disNcopen("discharge_data", "data/", locations)
#' # The result is a data frame with discharge time series for all specified locations.
#'
#' @export
#' @importFrom ncdf4 nc_open ncvar_get
#'
#' @seealso
#' `ncdf4::nc_open` and `ncdf4::ncvar_get` for functions used to access and
#' extract data from NetCDF files.
#'
disNcopen=function(fname,dir,outloc){
  ncdis=paste0(dir,"/",fname,".nc")
  ncd=nc_open(ncdis)
  name.vb=names(ncd[['var']])
  namev=name.vb[1]
  time <- ncvar_get(ncd,"time")
  lt=length(time)
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncd,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncd,name.lat) 
  lla=length(latdat)
  
  outllplus=matrix(-9999, nrow = lt*length(outloc[,1]), ncol = 5)
  outllplus=as.data.frame(outllplus)
  for ( idc in 1:length(outloc[,1])){
    print(idc)
    idm=1+(idc-1)*lt
    start=c(outloc$idlo[idc],outloc$idla[idc],1)
    count=c(1,1,lt)
    outlets = ncvar_get(ncd,namev,start = start, count= count)
    outlets=as.vector(outlets)
    outid=rep(outloc[idc,1],length(time))
    #lonlatloop=expand.grid(c(1:lla),c(1:lt))
    lon=rep(londat[idc],length(time))
    lat=rep(latdat[idc],length(time))
    outll=data.frame(outlets,outid,lon,lat,time)
    names(outllplus)=names(outll)
    fck=length(c(idm:(idm+lt-1)))
    outllplus[c(idm:(idm+lt-1)),]= outll
  }
  
  
  # outllplus=c()
  # for ( idc in 1:llo){
  #   print(idc)
  #   start=c(idc,1,1)
  #   count=c(1,lla,1)
  #   outlets = ncvar_get(ncd,namev,start = start, count= count)
  #   lna=length(which(!is.na(outlets)))
  #   if (lna>0){
  #     print("data")
  #     exact=(which(!is.na(outlets)))
  #     for (st in exact){
  #       start=c(idc,st,1)
  #       count=c(1,1,lt)
  #       outlets = ncvar_get(ncd,namev,start = start, count= count)
  #       outlets=as.vector(outlets)
  #       #lonlatloop=expand.grid(c(1:lla),c(1:lt))
  #       lon=rep(londat[idc],length(time))
  #       lat=rep(latdat[st],length(time))
  #       outll=data.frame(outlets,lon,lat,time)
  #     }
  #     outllplus=rbind(outllplus,outll)
  #   }
  
  #}
  
  return (outllplus)
}

#' Extract Discharge Data from NetCDF File at Specific Location
#'
#' @description
#' The `disNcopenloc` function opens a NetCDF file and extracts discharge data
#' for a specified location based on its identifier. It retrieves a time series
#' for the given location and constructs a data frame with discharge values along
#' with corresponding longitude, latitude, and time.
#'
#' @param fname A character string specifying the file name (without the extension)
#'   of the NetCDF file containing the discharge data.
#' @param dir A character string specifying the directory path where the NetCDF
#'   file is located.
#' @param outloc A data frame containing location information, including IDs and
#'   corresponding grid indices for latitude ('idla') and longitude ('idlo').
#' @param idc An integer representing the location identifier index for which the
#'   discharge data will be extracted.
#'
#' @return A data frame with columns for discharge values ('outlets'), location
#'   identifier ('outid'), longitude ('lon'), latitude ('lat'), and time ('time').
#'
#' @examples
#' # Assuming the NetCDF file 'discharge_data.nc' is in the directory 'data/',
#' # 'locations' is a data frame with grid indices, and 'location_id' is the
#' # identifier index:
#' discharge_series <- disNcopenloc("discharge_data", "data/", locations, location_id)
#' # The result is a data frame with the discharge time series for the specified location.
#'
#' @export
#' @importFrom ncdf4 nc_open ncvar_get
#'
#' @seealso
#' `ncdf4::nc_open` and `ncdf4::ncvar_get` for functions used to access and
#' extract data from NetCDF files.
#'
disNcopenloc=function(fname,dir,outloc,idc){
  ncdis=paste0(dir,"/",fname,".nc")
  ncd=nc_open(ncdis)
  name.vb=names(ncd[['var']])
  namev="dis"
  time <- ncvar_get(ncd,"time")
  lt=length(time)
  
  #timestamp corretion
  name.lon="lon"
  name.lat="lat"
  londat = ncvar_get(ncd,name.lon) 
  llo=length(londat)
  latdat = ncvar_get(ncd,name.lat) 
  lla=length(latdat)
  
  idm=1+(idc-1)*lt
  start=c(outloc$idlo[idc],outloc$idla[idc],1)
  count=c(1,1,lt)
  outlets = ncvar_get(ncd,namev,start = start, count= count)
  outlets=as.vector(outlets)
  outid=rep(outloc[idc,1],length(time))
  #lonlatloop=expand.grid(c(1:lla),c(1:lt))
  lon=rep(londat[start[1]],length(time))
  lat=rep(latdat[start[2]],length(time))
  outll=data.frame(outlets,outid,lon,lat,time)
  
  return (outll)
}


# Plotting and results ---------------------------------------------------------


#' Create Spatial Skill Assessment Plots
#'
#' @description
#' The `SpatialSkillPlot` function generates spatial plots for skill assessment metrics
#' using point data and corresponding NetCDF input files. It supports various metrics
#' like correlation, bias, variability, and the Kling-Gupta Efficiency (KGE). It
#' transforms the coordinates to UTM, matches the skill assessment data with the
#' provided points, and creates a ggplot object for visualization.
#'
#' @param points A data frame or matrix with longitude and latitude coordinates
#'   of points to be plotted. The first two columns should be longitude and latitude,
#'   respectively.
#' @param metric A character string indicating the skill assessment metric to be used
#'   for the plots. Supported values are "r" for correlation, "b" for bias, "v" for
#'   variability, and "kge" for Kling-Gupta Efficiency.
#' @param inputfile A character string specifying the name of the input file containing
#'   the skill assessment data to be matched with points. The file should be a CSV
#'   located in a directory specified by `valid_path`.
#'
#' @return A list containing two elements: `points`, which is a data frame with updated
#'   skill metrics, and `map`, which is a ggplot object representing the spatial
#'   visualization of the skill assessment.
#'
#' @examples
#' # Assuming 'points' is a data frame with coordinates, 'metric' is "r", and
#' # 'inputfile' is the name of the CSV file:
#' result <- SpatialSkillPlot(points, metric, inputfile)
#' # The result is a list with the points data frame and the ggplot object.
#'
#' @export
#' @importFrom sp SpatialPoints CRS spTransform
#' @importFrom rgdal readOGR
#' @importFrom rnaturalearth ne_countries
#' @importFrom sf st_as_sf st_transform
#' @importFrom ggplot2 ggplot geom_sf scale_x_continuous labs coord_sf scale_color_manual theme
#' @importFrom dplyr inner_join
#' @importFrom reshape2 melt
#'
#' @seealso
#' `sp::SpatialPoints`, `rgdal::readOGR`, `rnaturalearth::ne_countries`,
#' `sf::st_as_sf`, `ggplot2::ggplot`, `dplyr::inner_join`, `reshape2::melt`
#' for functions used to perform spatial transformations, plot creation, and data manipulation.
#'
SpatialSkillPlot<-function(points,metric,inputfile){
  #Plot parameters
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
  basemap=w2
  
  myfile=inputfile
  skill=read.csv(paste0(valid_path,myfile))
  if (length(skill[1,])==3) skill=skill[,-1]
  if (length(skill[1,])>2){
    years=c(1950:2020)
    skill=skill[-73]
    names(skill)[c(2:72)]=years
  }
  names(skill)[1]="Station_ID"
  
  #Only keep valid stations
  valid=match(points$V1, skill$Station_ID)
  skillv=skill[valid,]
  skillm<-melt(skillv,id.vars ="Station_ID")
  skillm2=skillm[which(skillm$value!=0),]
  skillm2$value[which(is.infinite(skillm2$value))]=NA
  
  #average
  if (length(skill[1,])>2){
    dat_skill=aggregate(list(skill=skillm2$value),
                        by = list(Station_ID=skillm2$Station_ID),
                        FUN = function(x) c(mean= mean(x,na.rm=T),median=median(x,na.rm=T),max=max(x,na.rm=T)))
    
    dat_skill <- do.call(data.frame, dat_skill)
  }else{
    dat_skill=skillm2[,c(1,3)]
  }
  
  
  points=inner_join(points,dat_skill,by=c("V1"="Station_ID"))
  if (length(skill[1,])>2){
    points$skill=points$skill.mean
    
  }else{
    points$skill=points$value
  }
  
  #For correlation
  if (metric=="r"){
    palet=c(hcl.colors(9, palette = "YlGnBu", alpha = NULL, rev = T, fixup = TRUE))
    limi=c(0,1)
    var="r"
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
    var="kge"
  }
  parpl <- st_as_sf(points, coords = c("Var1", "Var2"), crs = 4326)
  parpl <- st_transform(parpl, crs = 3035)
  parpef=parpl[which(parpl$calib==TRUE),]
  if (metric =="kge"){
    tsize=size=12
    map=ggplot(basemap) +
      geom_sf(fill="gray95", color=NA) +
      geom_sf(data=parpl,aes(geometry=geometry,colour=factor(kgecode),size=UpA),alpha=.7,stroke=0,shape=16)+ 
      geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="grey20",alpha=1,stroke=0.05,shape=1)+ 
      #geom_sf(data=parpef,aes(geometry=geometry),col="black",alpha=1,size=1.4,stroke=0.1,shape=8)+ 
      geom_sf(fill=NA, color="grey20") +
      scale_x_continuous(breaks=seq(-30,40, by=5)) +
      scale_size(range = c(1, 3), trans="sqrt")+
      labs(x="Longitude", y = "Latitude")+
      coord_sf(xlim = c(min(nco[,1]),max(nco[,1])), ylim = c(min(nco[,2]),max(nco[,2])))+
      scale_color_manual(values = colorz, labels= kgelabs, name="KGE'")   +
      labs(x="Longitude", y = "Latitude") +
      guides(colour = guide_legend(override.aes = list(size = 8,shape=15), reverse=T), size = "none")+
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
  }else{
    print(metric)
    #plot for correlation even if called skill
    map=ggplot(basemap) +
      geom_sf(fill="gray95", color=NA) +
      geom_sf(data=parpl,aes(geometry=geometry,fill=skill,size=UpA),color="transparent",alpha=.7,shape=21,stroke=0)+ 
      geom_sf(data=parpl,aes(geometry=geometry,size=UpA),col="grey20",alpha=1,stroke=0.05,shape=1)+ 
      geom_sf(fill=NA, color="grey20") +
      scale_x_continuous(breaks=seq(-30,40, by=5)) +
      scale_size(range = c(1, 3), trans="sqrt")+
      scale_fill_gradientn(
        colors=palet,n.breaks=5,oob = scales::squish, limits=limi, name=metric) +
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
  }
  return(list(points,map))
}
#' Create Scatter Plot for Validation of Simulation against Observations
#'
#' @description
#' The `ScatterValid` function generates a scatter plot for the purpose of validating
#' simulation data against observational data. It performs a log transformation for
#' both datasets, calculates correlation and R-squared values, and plots the data
#' on a log-log scale with a color density representation. The function assumes that
#' both input datasets have been pre-processed to align observation and simulation
#' data by station ID.
#'
#' @param obsin A data frame containing observational data with columns including
#'   'Station_ID' and observed values.
#' @param nobs A character string representing the label for the y-axis (observed values).
#' @param simin A data frame containing simulation data with columns including
#'   'Station_ID' and simulated values.
#' @param nsim A character string representing the label for the x-axis (simulated values).
#' @param scplot A character string specifying the name of the scatter plot file to save.
#' @param valid_path A character string specifying the directory path where validation
#'   data files are located.
#' @param ValidSta A data frame containing valid station IDs.
#'
#' @return A ggplot object representing the scatter plot of the validation results.
#'
#' @examples
#' # Assuming 'obsin' and 'simin' are data frames with aligned validation data, and
#' # 'ValidSta' contains the valid station IDs:
#' scatter_plot <- ScatterValid(obsin, "Observed", simin, "Simulated", "scatter_plot", "path/to/validation/", ValidSta)
#' # The result is a ggplot object of the scatter plot.
#'
#' @export
#' @importFrom ggplot2 ggplot geom_jitter geom_abline scale_x_log10 scale_y_log10
#'   annotate scale_color_viridis theme
#' @importFrom dplyr inner_join
#' @importFrom scales log10_minor_break margin
#' @importFrom viridis viridis_pal
#'
#' @seealso
#' `ggplot2::ggplot` for creating the plot,
#' `dplyr::inner_join` for merging datasets,
#' `scales::log10_minor_break` for minor break calculation on log scales.
#'
ScatterValid=function(obsin,nobs,simin,nsim,scplot,valid_path,ValidSta){
  
  # obs=read.csv(paste0(valid_path,fileobs))
  # sim=read.csv(paste0(valid_path,filesim))
  # names(obs)[1] = names(sim)[1]="Station_ID"
  # 
  # 
  # if (length(obs[1,])>2){
  #   years=c(1950:2020)
  #   obs=obs[-73]
  #   names(obs)[c(2:72)]=years
  #   sim=sim[-73]
  #   names(sim)[c(2:72)]=years
  # }
  # 
  # #Only keep valid stations
  # valid=match(ValidSta$V1, obs$Station_ID)
  # obsv=obs[valid,]
  # obstest<-melt(obsv,id.vars ="Station_ID")
  # obstest2=obstest[which(obstest$value>0),]
  # 
  # length(unique(obstest2$Station_ID))
  # 
  # valid=match(ValidSta$V1, sim$Station_ID)
  # simv=sim[valid,]
  # simtest<-melt(simv,id.vars ="Station_ID")
  # simtest2=simtest[which(simtest$value>0),]
  # 
  # length(unique(simtest2$Station_ID))
  # 
  
  #Output plots
  
  if (length(obsin[1,])>3){
    obs_sim=inner_join(obsin,simin, by=c("Station_ID","variable"))
    obs_sim$value.x[which(is.infinite(obs_sim$value.x))]=NA
    obs_sim$value.y[which(is.infinite(obs_sim$value.y))]=NA
    length(unique(obs_sim$Station_ID))
    names(obs_sim)[c(3,4)]=c("obs","sim")
    dat=as.data.frame(obs_sim[c(1,2,3,4)])
    
    #average
    dat_av=aggregate(list(obs=dat$obs,sim=dat$sim),
                     by = list(station=dat$Station_ID),
                     FUN = function(x) c(mean=mean(x,na.rm=T)))
    
    dat <- do.call(data.frame, dat_av)
    dat=dat[which(!is.na(dat$obs)),]
  }else{
    obsin=obsin[,c(1,3)]
    simin=simin[,c(1,3)]
    obs_sim=inner_join(obsin,simin, by=c("Station_ID"))
    obs_sim$value.x[which(is.infinite(obs_sim$value.x))]=NA
    obs_sim$value.y[which(is.infinite(obs_sim$value.y))]=NA
    length(unique(obs_sim$Station_ID))
    names(obs_sim)[c(2,3)]=c("obs","sim")
    dat=obs_sim
  }
  r=cor(dat$sim,dat$obs)
  r2=round(r^2,3)
  
  mmax=max(max(dat$obs,dat$sim))
  mmin=min(min(dat$obs,dat$sim))
  tmax=round(mmax/1000)*1000
  tmin=max(mmin,1e-3)
  
  bi=2000
  if (tmax<10000) bi=1000
  if (tmax<4000) bi=500
  #dat$year=as.numeric(as.character(dat$variable))
  # dat$obsj=jitter(dat$obs,0.001)+0.0001
  # dat$simj=jitter(dat$sim,0.001)+0.0001
  dat$density <- get_density(log(dat$obs), log(dat$sim), n = 100)
  #dat$density <- get_density(log(dat$obsj), log(dat$simj), n = 100)
  ggplot(dat,environment=environment()) + 
    # geom_point(aes(x=obs, y=sim,col=density),stroke=0,size=3,alpha=0.2,shape=16) +
    # geom_point(aes(x=obs, y=sim,col=density),stroke=0,size=3,alpha=0.2,shape=16) +
    geom_jitter(aes(x=obs, y=sim,col=density),stroke=0,size=3,alpha=0.5,shape=16, height=.1, width=.0)+
    geom_abline(slope=1, intercept=0, lwd =1, alpha=1,col="darkgreen",linetype="dashed")+
    scale_x_log10(name=nobs,
                  breaks=c(0.1,1,10,100,1000,10000), minor_breaks = log10_minor_break(),
                  labels=c("0.1","1","10","100","1000","10000"), limits=c(tmin,tmax)) +
    scale_y_log10(name=nsim,
                  breaks=c(0.1,1,10,100,1000,10000), minor_breaks = log10_minor_break(),
                  labels=c("0.1","1","10","100","1000","10000"), limits=c(tmin,tmax)) +
    # scale_x_continuous(name=nobs,breaks=seq(0,tmax, by=bi))+
    # scale_y_continuous(name=nsim,breaks=seq(0,tmax, by=bi))+
    annotate("label", x=3*tmin, y=max(tmax), label= paste0("R2 = ",r2),size=5)+
    scale_color_viridis(option="A")+
    theme(axis.title=element_text(size=16, face="bold"),
          axis.text = element_text(size=16),
          panel.background = element_rect(fill = "white", colour = "white"),
          panel.grid = element_blank(),
          panel.border = element_rect(linetype = "solid", fill = NA, colour="black"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          legend.position = "none",
          plot.margin = margin(1,1,1,1, "cm"),
          panel.grid.major = element_line(colour = "grey80"),
          panel.grid.minor = element_line(colour = "grey90",linetype="dashed"),
          legend.key = element_rect(fill = "transparent", colour = "transparent"),
          legend.key.size = unit(.8, "cm"))
  
  #dat$rat=dat$sim/dat$obs
  #return(dat)
  #ggsave(paste0(scplot,".jpg"), width=20, height=15, units=c("cm"),dpi=1500)
}

#' Analyze and Plot Hydrological Episodes
#'
#' @description
#' The `StatEpisode` function analyzes hydrological episodes based on observed and
#' simulated discharge data within a specified window around a given date. It
#' calculates key statistics such as the Kling-Gupta Efficiency (KGE) and error metrics.
#' Optionally, it can generate a plot showing the discharge over time for both
#' observed and simulated data.
#'
#' @param data A data frame containing the hydrological data, including columns for
#'   date, observed discharge (Q), and simulated discharge (Qs).
#' @param date.deb A date object or character string indicating the central date around
#'   which the episode is analyzed. If NA, the function uses the full data range.
#' @param win An integer specifying the window size (in days) used for the episode
#'   analysis. Defaults to 10 days.
#' @param qx.min A numeric value indicating the minimum discharge to consider in the
#'   analysis. If NULL, the function uses the maximum observed discharge.
#' @param mode A character string indicating the type of episode to analyze: "FL" for
#'   flood or "DR" for drought. Defaults to "FL" (flood).
#' @param obs A logical value indicating whether to include observed discharge data
#'   in the analysis and plot. Defaults to TRUE.
#' @param plot A logical value indicating whether to generate a plot of the episode.
#'   Defaults to FALSE.
#'
#' @return If plotting is disabled, the function returns a numeric vector containing
#'   the KGE, error metric, maximum observed discharge (QmObs), and the year of the
#'   episode. If plotting is enabled, no explicit return value is provided, but a plot
#'   is generated in the current graphics device.
#'
#' @examples
#' # Assuming 'data' is a data frame with appropriate structure:
#' result <- StatEpisode(data, date.deb = "2000-05-15", win = 10,
#'                       mode = "FL", obs = TRUE, plot = TRUE)
#' # The result is a numeric vector with statistics, and a plot is generated.
#'
#' @export
#' @importFrom graphics par plot new axis box points legend mtext
#' @importFrom stats na.omit cor
#'
#' @seealso
#' `graphics::plot`, `stats::cor` for plotting and statistical functions used within the analysis.
#'
StatEpisode <- function (data, date.deb = NA, win = 10,
                         qx.min = NULL, mode = "FL",obs=TRUE,plot=FALSE){
  x <- data
  # if (class(x$date)[1] == "POSIXt") {
  # if (inherits(x$date,"POSIXt")){
  #   date.deb <- as.POSIXct(date.deb, tz = "utc")
  # }
  if (mode=="FL"){
    ind <- which(x$date >= date.deb-win/2)
    ip <- which(ind >= min(ind) & ind <= (min(ind) + win))
    ind <- ind[ip]
  }
  if (mode=="DR"){
    ind <- which(x$date >= date.deb-1)
    ip <- which(ind >= min(ind) & ind <= (min(ind) + win+1))
    ind <- ind[ip]
  }
  
  if (length(ind) < 3) {
    cat("Episodes du ", format(date.deb, "%d-%b-%y"), "ignore (n<3)\n")
    return(NULL)
  }
  if (obs==T) xep <- x[ind, c("date", "Q" , "Qs")]
  if (obs==F) {xep <- x[ind, c("date", "Qs")]
  xep$Q=xep$Qs}
  # if (!is.null(fp)) {
  #   tmp <- abat.niv(xep, fp)
  #   xep <- data.frame(xep, Peff = tmp$P)
  # }
  if (is.null(qx.min))
    qx.min <- max(x$Q, na.rm = TRUE)
  #     if (is.null(px.min))
  #         px.min <- max(x$P, na.rm = TRUE)
  if (mode=="FL"){
    q.max <- max(c(xep$Q, xep$Qs, qx.min), na.rm = TRUE)
    q.maxsim<-max(na.omit(xep$Qs))
    if (obs==T) {q.maxobs<-max(na.omit(xep$Q))
    ec.pic = round((q.maxobs-q.maxsim)/q.maxobs*100,1)}
    # p.max <- max(c(xep$P, px.min), na.rm = TRUE)
    q.range <- c(0, 1.5 * q.max)
    # p.range <- c(0, 3 * p.max)
    ec=ec.pic
  }
  if (mode=="DR"){
    q.max <- max(c(xep$Q, xep$Qs, qx.min), na.rm = TRUE)
    q.minsim<-min(na.omit(xep$Qs))
    if (obs==T){ q.minobs<-min(na.omit(xep$Q))
    ec.vol = round((sum(na.omit(xep$Qs)-na.omit(xep$Q)))/sum(na.omit(xep$Q))*100,1)}
    # p.max <- max(c(xep$P, px.min), na.rm = TRUE)
    q.range <- c(0, 1.2 * q.max)
    # p.range <- c(0, 3 * p.max)
    ec=ec.vol
  }
  if (obs==T) kge<-c(kge=KGE(xep$Q,xep$Qs),r=cor(xep$Q,xep$Qs))
  if (plot==T){
    par(mar = c(2, 3, 2, 3), mgp = c(1.7, 0.5, 0))
    # plot.window(range(xep$date))#, p.range)
    plot.new()
    
    #     plot.hyeto(xep$date, xep$P, curve = TRUE, ylim = rev(p.range),
    #         grid = FALSE, col = "lightblue", border = "azure", lwd = 2)
    #     abline(v = axis.POSIXct(1, xep$date), col = "lightgrey",
    #         lty = 3)
    #     if (!is.null(fp))
    #         plot.hyeto(xep$date, xep$Peff, curve = TRUE, ylim = rev(p.range),
    #             grid = FALSE, border = NA)
    #     axis(4, pretty(c(0, p.max), n = 5))
    #     mtext("PrÃ©cipitations (mm)", 4, 1.7, cex = par()$cex)
    
    par(new = TRUE)
    plot(xep$date, xep$Qs, ylim = q.range, type = "n", axes = FALSE,
         xlab = "", ylab = expression(paste("Discharge ", (m^3/s),
                                            sep = " ")))
    h = axis(2, pretty(c(0, q.max*1.5), n = 5))
    axis(1,xep$date, label=format(xep$date,"%d %b"),cex.axis=.8)
    box()
    points(xep$date, xep$Qs, type = "l", col = "red", pch = 2,
           lwd = 2, cex = 0.6)
    
    if (obs==T) {points(xep$date, xep$Q, type = "l", col = "royalblue", lwd = 2)
      legend("left", legend = c("Histo.", "Recal."), 
             col = c("royalblue","red"), pch = c(NA, NA), lwd = 2, bty = "n", cex = 0.7)
    }else{
      
    }
    if (mode =="FL"){
      #abline(h=quantile(data$Q,0.90,na.rm=T),col="darkgrey",lwd=2) 
      if (obs==T){
        txt.lab <- paste("KGE =", round(kge[1],2), "\nEc-pic =",ec.pic," %","\nr =",round(kge[2],3))
        text(mean(xep$date), 1.2*q.max-0.1*q.max, pos=4, lab=txt.lab,cex=0.7)
      }
      titre <- paste("Episode du ", format(xep$date[which.max(xep$Q)],
                                           "%a %d %B %Y"))
      ec=ec.pic
    }
    if (mode =="DR"){
      abline(h=quantile(data$Qs,0.05,na.rm=T),col="darkgrey",lwd=2) 
      if (obs==T){
        txt.lab <- paste("KGE =", round(kge[1],2), "\nEc-vol =",ec.vol," %","\nr =",round(kge[2],3))
        text(mean(xep$date), 1.2*q.max-0.1*q.max, pos=4, lab=txt.lab,cex=0.7)
      }
      titre <- paste("Episode du ", format(xep$date[1],"%a %d %B %Y") ,"au ",
                     format(xep$date[win], "%a %d %B %Y"))
      ec=ec.vol
    }
    mtext(titre,3,font = 2,adj = 0,line = 0.5,cex = .75)
    # legend("right", legend = c("Pluie", "Neige"), pch = 22, lwd = NA,
    # lty = NA, col = "lightblue", pt.bg = c("lightblue", "azure"),
    # pt.lwd = 1, pt.cex = 2, bty = "n", cex = 0.7)
    # legend("right", legend = "Pluie", pch = 22, lwd = NA,
    #         lty = NA, col = "azure", pt.bg = "lightblue",
    #         pt.lwd = 1, pt.cex = 2, bty = "n", cex = 0.7)
  }
  if (obs==T) return(c(kge,error=ec,QmObs=q.maxobs,QmObs=q.maxobs,year=year(xep$date)[1]))
}


# Usage -----------------------------------------------------------------------
#' Plot Hydrological Regime
#'
#' @description
#' The `plotRegime` function generates a plot of the hydrological regime over a period,
#' showing the mean, 10th percentile, and 90th percentile daily flows. Optionally,
#' a threshold can be set to filter the data.
#'
#' @param data A data frame with two columns: 'date' (Date objects) and 'Q' (numeric
#'   values of discharge).
#' @param catch A character string representing the name of the catchment.
#' @param seuil An optional numeric threshold for the discharge data; only values above
#'   this threshold will be included in the plot. If NULL, no threshold is applied.
#'
#' @return A data frame containing the calculated mean, 10th percentile (q10), and 90th
#'   percentile (q90) discharge values for each day of the year.
#'
#' @examples
#' # Assuming 'discharge_data' is a data frame with dates and discharge values:
#' Qc <- plotRegime(data = discharge_data, catch = "Catchment Name", seuil = 10)
#' # The result is a plot of the hydrological regime and a data frame 'Qc'.
#'
#' @export
#' @importFrom graphics plot lines axis polygon legend box mtext
#' @importFrom grDevices colorRampPalette
#'
#' @seealso
#' Graphics functions like `plot`, `lines`, `axis`, `polygon`, and `legend` for
#' creating the plot. The `grDevices::colorRampPalette` function for creating
#' color palettes.
#'
plotRegime <- function(data, catch, seuil=NULL){
  #~ data <- as.simu(data)
  names(data)=c("date","Q")
  if( is.null(seuil) ) seuil <- -Inf
  ## Suppression des lignes sans couples Qobs/Qsim
  data <- subset( data, !is.na(Q))
  ## Preparation de la sequence de mois
  deb=data$date[4]
  mois.deb <-  seq(as.Date(deb), by="month", length=12)
  
  period=paste0(format(range(data$date)[1],"%b %Y"),"-",format(range(data$date)[2],"%b %Y"))
  main = bquote(.(catch)~ " Regime | Period: "~ .(period))
  jours <- as.numeric(format(data$date,"%j"))
  ## Iddinces des lignes de meme categorie
  ind.j <- tapply(seq(length(jours)), jours, c)
  ind.j <- ind.j[-366]
  ## Calcul
  Qc <- data.frame(date=as.numeric(names(ind.j)),
                   mean=sapply(ind.j, function(x) mean(data$Q[x], na.rm=TRUE)),
                   q10=sapply(ind.j, function(x) quantile(data$Q[x],0.1, na.rm=TRUE)),
                   q90=sapply(ind.j, function(x) quantile(data$Q[x],.9, na.rm=TRUE)))
  
  qlim=c(0,1.1*max(Qc[,4]))
  plot(Qc$date, Qc$mean, type="n", axes=FALSE, ylim=qlim,
       xlab = NA, ylab = expression(paste("Debit (",m^3/s,")")))
  mtext(main,3,font = 2,line = 0.5,cex = .75)
  lines(Qc$date, Qc$mean, col="darkblue",lwd=2)
  lines(Qc$date, Qc$q10, col="blue",lwd=1,lty=2)
  lines(Qc$date, Qc$q90, col="blue",lwd=1,lty=2)
  abline(v = format(mois.deb,"%j"), col="lightgrey", lty=3)
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q10,rev(Qc$q90)),
          col=alpha("lightblue",.4),border="transparent")
  axis(2)
  axis(1, format(mois.deb,"%j"), label=format(mois.deb,"%b"),cex.axis=.8)
  box()
  legend("topleft", leg=c("mean","Q90","Q10"),
         lwd=c(2,2,1,2), col=c("darkblue","blue","blue"),
         cex=.8, lty=c(1,2,2),bty="n")
  return(Qc)
}

#' Plot Comparative Hydrological Regimes
#'
#' @description
#' The `plotQj` function plots comparative hydrological regimes using observed and
#' simulated discharge data. It visualizes the mean daily flow, the interquartile range
#' for two different periods, and the absolute deviation between observed and simulated
#' flows.
#'
#' @param data A data frame with four columns: 'date1' and 'date2' (Date objects for
#'   two periods), 'Q1' (numeric values of observed discharge), and 'Q2' (numeric values
#'   of simulated discharge).
#' @param catch A character string representing the name of the catchment area.
#' @param UpAloc A numeric value representing the upstream catchment area in square
#'   kilometers.
#'
#' @return The function does not explicitly return a value but generates a plot showing
#'   the hydrological regimes for the two periods, the interquartile range, and
#'   absolute deviation between observed and simulated flows.
#'
#' @examples
#' # Assuming 'discharge_data' is a data frame with dates and discharge values for two
#' # periods:
#' plotQj(data = discharge_data, catch = "Catchment Name", UpAloc = 1000)
#' # This generates a plot of comparative hydrological regimes.
#'
#' @export
#' @importFrom graphics plot lines axis polygon legend box mtext
#' @importFrom grDevices colorRampPalette
#'
#' @seealso
#' Graphics functions like `plot`, `lines`, `axis`, `polygon`, and `legend` for
#' creating the plot. The `grDevices::colorRampPalette` function for creating
#' color palettes.
#'
plotQj <- function(data,catch,UpAloc){
  #~ data <- as.simu(data)
  names(data)=c("date1","Q1","date2","Q2")
  ## Suppression des lignes sans couples Qobs/Qsim
  data <- subset( data, !is.na(Q1) & !is.na(Q2))
  ## Preparation de la sequence de mois
  mois.deb <-  seq(as.Date("1950-01-15"), by="month", length=12)
  ## Vecteur de catÃ©gories
  period1=paste0(format(range(data$date1)[1],"%Y"),"-",format(range(data$date1)[2],"%Y"))
  period2=paste0(format(range(data$date2)[1],"%Y"),"-",format(range(data$date2)[2],"%Y"))
  print(period2)
  n.area=UpAloc
  name.riv=catch
  #main = bquote(.(name.riv)~ " Regime (Area="~.(n.area)~ km^2~") | Periods: "~ .(period1) ~" vs" ~ .(period2))
  main = bquote(.(name.riv)~ " Regime (Area="~.(n.area)~ km^2~")")
  ## Suppression des lignes sans couples Qobs/Qsim
  ## Preparation de la sequence de mois
  jours <- as.numeric(format(data$date1,"%j"))
  ## Iddinces des lignes de mÃ©me catÃ©gorie
  ind.j <- tapply(seq(length(jours)), jours, c)
  ind.j <- ind.j[-366]
  
  Qc <- data.frame(date=as.numeric(names(ind.j)),
                   obs=sapply(ind.j, function(x) mean(data$Q1[x], na.rm=TRUE)),
                   sim=sapply(ind.j, function(x) mean(data$Q2[x], na.rm=TRUE)),
                   q25o=sapply(ind.j, function(x) quantile(data$Q1[x],0.25, na.rm=TRUE)),
                   q75o=sapply(ind.j, function(x) quantile(data$Q1[x],.75, na.rm=TRUE)),
                   q25s=sapply(ind.j, function(x) quantile(data$Q2[x],0.25, na.rm=TRUE)),
                   q75s=sapply(ind.j, function(x) quantile(data$Q2[x],.75, na.rm=TRUE)))
  
  md=mean(data$Q1)
  qlim=c(0, max(unlist(Qc[,c(5,7)]),2*md))

  
  ## ParamÃ©tres graphiques
  plot(Qc$date, Qc$obs, type="n", axes=FALSE, ylim=qlim,xaxs="i",yaxs="i",
       xlab = NA, ylab="")
  mtext(main,3,font = 2,line = 0.5,cex = 1.5)
  abline(v = format(mois.deb,"%j"), col="lightgrey", lty=2)
  lines(Qc$date, Qc$obs, col="blue",lwd=2)
  lines(Qc$date, Qc$sim, col="red",lwd=2)
  axis(2,cex.axis=1)
  title(ylab = expression(paste("Q (",m^3/s,")")),cex.lab=1.5,line=2)
  axis(1, format(mois.deb,"%j"), label=format(mois.deb,"%b"),cex.axis=1)
  box()
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q25o,rev(Qc$q75o)),
          col=alpha("lightblue",.3),border="transparent")
  
  polygon(c(Qc$date,rev(Qc$date)),c(Qc$q25s,rev(Qc$q75s)),
          col=alpha("indianred",.3),border="transparent")
  ## Trace des erreurs absolues
  polygon(c(0,Qc$date,365),c(0,abs(Qc$sim-Qc$obs),0),
          col=alpha("lightgrey",.4),border="grey")
  ## Trace de la grille mensuelle
  
  ## Trace des courbes de dÃ©bits journaliers simules et observes
  
  ## Definition de la legende
  legend("topleft", leg=c("1990-2020","1951-1981","deviation"),
         lwd=c(2,2,1,2), col=c("red","blue","grey"),
         cex=1.2, lty=c(1,1,2,1), bg=alpha("white",.6))
  #return(Qc)
} 

#' Compare River Discharge Regimes
#'
#' @description
#' The `Regime_comp` function compares the river discharge regimes for a specified
#' catchment by analyzing and plotting discharge data from different periods. It loads
#' the outlet data for a specific grid square, calculates distances to the river,
#' extracts relevant discharge data, and generates a comparative regime plot.
#'
#' @param rspace A data frame or matrix with predefined grid indices and spatial
#'   information for extraction.
#' @param Nsq An integer representing the specific grid square to be analyzed.
#' @param river Spatial data for the river, typically as an `sf` object or a matrix of
#'   longitude and latitude coordinates.
#' @param hydroDir A character string specifying the directory path where the hydrological
#'   data files are located.
#' @param catch A character string representing the name of the catchment area.
#'
#' @return The function does not return a value but generates a plot showing the comparative
#'   hydrological regimes for the selected periods.
#'
#' @examples
#' # Assuming 'grid_space' is a data frame with grid indices, 'grid_square' is the index
#' # for the square to analyze, 'river_data' is the spatial data for the river, 'hydro_dir'
#' # is the path to the hydrological data, and 'catchment' is the name of the catchment area:
#' Regime_comp(grid_space, grid_square, river_data, hydro_dir, catchment)
#' # This generates a plot comparing the discharge regimes.
#'
#' @export
#' @importFrom graphics plot lines axis polygon legend box mtext
#' @importFrom sf st_as_sf st_distance
#' @importFrom utils read.csv
#'
#' @seealso
#' `sf::st_as_sf` and `sf::st_distance` for spatial operations on river and outlet data,
#' `graphics` functions for creating the comparative regime plot.
#' `plotQj` function for the creation of regime plots
#'
Regime_comp=function(rspace,Nsq,river,hydroDir,catch,outletname){
  #load specific outlet file for this square:
  nrspace=rspace[Nsq,]
  
  outEFAS=outletopen(hydroDir,outletname)
  outEFAS$latlong=paste(round(outEFAS$Var1,4),round(outEFAS$Var2,4),sep=" ")
  
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
  
  rloc2=which(round(outEFAS$Var1,3)==riv[1] & round(outEFAS$Var2,3)==riv[2])
  outR=outEFAS[rloc2,]
  UpAloc=round(upArea[outR$idlo,outR$idla])
  print(paste0("Upstream area: ", UpAloc))
  df.dis=dists 
  timeStamps=unique(as.Date(df.dis$time,origin="1979-01-01"))
  timeStamps=as.POSIXct(timeStamps-1/24)
  txx=timeStamps
  df.dis$timeStamps=timeStamps
  
  timebound1=c(as.POSIXct(("1951-01-01 00:00:00")))
  timebound2=c(as.POSIXct(("1981-12-31 00:00:00")))
  tb1=match(timebound1,df.dis$timeStamps)
  tb2=match(timebound2,df.dis$timeStamps)
  
  timebound3=c(as.POSIXct(("1990-01-01 00:00:00")))
  timebound4=c(as.POSIXct(("2020-12-31 00:00:00")))
  tb3=match(timebound3,df.dis$timeStamps)
  tb4=match(timebound4,df.dis$timeStamps)
  
  names(df.dis)[c(1,2)]=c("dis","outlets")
  ts1=df.dis[c(tb1:tb2),]
  
  #Now the regime
  dts1=ts1[,c(6,1)]
  
  #plot regime allyears:
  dts=df.dis[,c(6,1)]
  
  ts2=df.dis[c(tb3:tb4),]
  
  #Now the regime
  dts2=ts2[,c(6,1)]
  
  data=data.frame(dts1,dts2)
  #30 days running mean for a smoother regime
  data[,2]=tsEvaNanRunningMean(data[,2],120)
  data[,4]=tsEvaNanRunningMean(data[,4],120)
  
  plotQj(data,catch,UpAloc)
  
}

# Circular statistics ----------------------------------------------------------
#' Calculate Seasonal Index and Strength
#'
#' @description
#' The `seasony` function computes the seasonal index and strength of a time series
#' based on circular statistics, often used in phenology or annual cycle studies. It
#' takes a vector of dates (day of the year) and calculates the mean direction and
#' length of the resultant vector, providing an index for the central tendency of
#' seasonal events and its consistency.
#'
#' @param x A numeric vector representing the day of the year for each event in the
#'   time series.
#'
#' @return A numeric vector containing two elements: the first element (`Di`) is the
#'   seasonal index (day of the year) representing the mean date of occurrence, and
#'   the second element (`R`) is the length of the mean resultant vector indicating
#'   the strength or consistency of the seasonal pattern. A higher `R` value indicates
#'   a more consistent seasonal pattern.
#'
#' @examples
#' # Example vector with day of the year for each event
#' days_of_year <- c(100, 105, 110, 95, 108)
#' 
#' # Calculate the seasonal index and strength
#' result <- seasony(days_of_year)
#' 
#' # The result contains the seasonal index and strength
#' # result[1] is the seasonal index, and result[2] is the strength
#'
#' @export
#'
#' @seealso
#' The concept of circular statistics as applied in phenology and the study of
#' seasonal cycles.
#'
seasony=function(x){
  theta=x*(2*pi/365.25)
  #plot(cos(theta),sin(theta),xlim=c(-1,1),ylim=c(-1,1))
  xi=1/(length(theta))*sum(cos(na.omit(theta)))
  yi=1/(length(theta))*sum(sin(na.omit(theta)))
  if (xi<=0){
    Di=(atan(yi/xi)+pi)*(365.25/(2*pi))
  }else if(xi>0 & yi>=0){
    Di=(atan(yi/xi))*(365.25/(2*pi))
  }else if(xi>0 & yi<0){
    Di=(atan(yi/xi)+2*pi)*(365.25/(2*pi))
  }
  R=sqrt(xi^2+yi^2)
  return(c(Di,R))
}

#' Calculate Seasonal Index from Date Vector
#'
#' @description
#' The `season1` function calculates the seasonal index based on circular statistics
#' from a vector of dates (day of the year). It is used to determine the central
#' tendency of seasonal events (e.g., peak timing of a recurring annual event).
#'
#' @param x A numeric vector representing the day of the year for each event in the
#'   time series. Missing values (NA) are allowed and will be ignored in the calculation.
#'
#' @return A single numeric value representing the seasonal index, which is the
#'   mean date of occurrence calculated as the day of the year. If the input vector
#'   has no valid dates (all NA), the function returns NA.
#'
#' @examples
#' # Example vector with day of the year for each event, including NA values
#' days_of_year <- c(100, NA, 110, 95, 108, NA)
#'
#' # Calculate the seasonal index
#' seasonal_index <- season1(days_of_year)
#'
#' # The result is the mean date of occurrence as the day of the year
#'
#' @export
#'
#' @seealso
#' `seasony` for a related function that also calculates the strength of the
#' seasonal pattern.
#'
season1=function(x){
  l1=length(which(!is.na(x)))
  if(l1>0){
    x=x[which(!is.na(x))]
    theta=x*(2*pi/365.25)
    # plot(theta)
    
    xi=1/(length(theta))*sum(cos(theta))
    yi=1/(length(theta))*sum(sin(theta))
    if (xi<=0){
      Di=(atan(yi/xi)+pi)*(365.25/(2*pi))
    }else if(xi>0 & yi>=0){
      Di=(atan(yi/xi))*(365.25/(2*pi))
    }else if(xi>0 & yi<0){
      Di=(atan(yi/xi)+2*pi)*(365.25/(2*pi))
    }
    R=sqrt(xi^2+yi^2)
  }else{Di=NA}
  return(Di)
}

#' Calculate Seasonality strength from Date Vector
#'
#' @description
#' The `season2` function calculates the seasonality strength based on circular statistics
#' from a vector of dates (day of the year). It is used to determine the central
#' tendency of seasonal events (e.g., peak timing of a recurring annual event).
#'
#' @param x A numeric vector representing the day of the year for each event in the
#'   time series. Missing values (NA) are allowed and will be ignored in the calculation.
#'
#' @return A single numeric value representing the seasonality strength. If the input vector
#'   has no valid dates (all NA), the function returns NA.
#'
#' @examples
#' # Example vector with day of the year for each event, including NA values
#' days_of_year <- c(100, NA, 110, 95, 108, NA)
#'
#' # Calculate the seasonal index
#' seasonal_strength <- season2(days_of_year)
#'
#' @export
#'
#' @seealso
#' `seasony` for a related function.
#'
season2=function(x){
  l1=length(which(!is.na(x)))
  if(l1>0){
    x=x[which(!is.na(x))]
    theta=x*(2*pi/365.25)
    # plot(theta)
    
    xi=1/(length(theta))*sum(cos(theta))
    yi=1/(length(theta))*sum(sin(theta))
    if (xi<=0){
      Di=(atan(yi/xi)+pi)*(365.25/(2*pi))
    }else if(xi>0 & yi>=0){
      Di=(atan(yi/xi))*(365.25/(2*pi))
    }else if(xi>0 & yi<0){
      Di=(atan(yi/xi)+2*pi)*(365.25/(2*pi))
    }
    R=sqrt(xi^2+yi^2)
  }else{Di=NA
  R=NA}
  return(R)
}


# Others -----------------------------------------------------------------------

#' Compute Density Estimates for Points on a 2D Plane
#'
#' @description
#' The `get_density` function calculates the density estimates for a set of points
#' on a two-dimensional plane using kernel density estimation. It is useful for
#' assessing the density of points in a scatter plot or other two-dimensional
#' visualization.
#'
#' @param x A numeric vector representing the x-coordinates of the points.
#' @param y A numeric vector representing the y-coordinates of the points.
#' @param ... Additional arguments passed to `MASS::kde2d`.
#'
#' @return A numeric vector of the same length as `x` and `y` containing the
#'   density estimates for each point.
#'
#' @examples
#' # Example vectors of x and y coordinates
#' x_coords <- rnorm(100)
#' y_coords <- rnorm(100)
#'
#' # Calculate the density estimates for these points
#' density_estimates <- get_density(x_coords, y_coords)
#'
#' # The result is a vector of density estimates for each (x,y) point
#'
#' @export
#' @importFrom MASS kde2d
#' @importFrom stats findInterval
#'
#' @seealso
#' `MASS::kde2d` for the underlying kernel density estimation function.
#'
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#' Generate Minor Logarithmic Breaks for Plotting
#'
#' @description
#' The `log10_minor_break` function generates a sequence of minor logarithmic breaks
#' for use in plotting functions, particularly useful when customizing the minor
#' breaks on a logarithmic scale.
#'
#' @param ... Additional arguments (currently not used).
#'
#' @return A function that takes a numeric vector `x` and returns a numeric vector
#'   of minor breaks in the original data scale (not log-transformed). These breaks
#'   are calculated between the major logarithmic breaks of `x`.
#'
#' @examples
#' # Generate a set of logarithmic minor breaks for plotting
#' x_data <- c(0.1, 1, 10, 100, 1000)
#' minor_breaks_func <- log10_minor_break()
#' minor_breaks <- minor_breaks_func(x_data)
#'
#' # Use the minor breaks in a plot (for example, with ggplot2)
#' # ggplot(data, aes(x, y)) +
#' #   geom_point() +
#' #   scale_x_log10(minor_breaks = minor_breaks)
#'
#' @export
#'
#' @seealso
#' Plotting functions in packages like `ggplot2`, which allow for customization
#' of minor breaks on a logarithmic scale.
#'
log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}

#' Compute Annual Maximum Values from Time Series Data
#'
#' @description
#' The `computeAnnualMaxima` function takes a two-column matrix or data frame 
#' with time stamps and corresponding time series data. It computes the annual 
#' maximum values for each year present in the time series.
#'
#' @param timeAndSeries A two-column matrix or data frame where the first column 
#'   contains time stamps and the second column contains the time series data.
#'
#' @return A list containing three elements: 'annualMax' which is a vector of 
#'   the annual maximum values, 'annualMaxDate' which is a vector of the dates 
#'   on which these maxima occur, and 'annualMaxIndx' which is a vector of the 
#'   indices in the original time series that correspond to the annual maxima.
#'
#' @examples
#' # Assume 'timeAndSeries' is a matrix with POSIXct time stamps in the first column
#' # and corresponding data values in the second column:
#' result <- computeAnnualMaxima(timeAndSeries)
#' # The result is a list containing the annual maximum values and their dates.
#'
#' @export
#'
#' @seealso
#' `tapply` for the function used to apply a function over subsets of data,
#' `which.max` for finding the index of the maximum value in a vector.
#'
computeAnnualMaxima <- function(timeAndSeries) {
  timeStamps <- timeAndSeries[,1]
  srs <- timeAndSeries[,2]
  
  tmvec <- as.Date(timeStamps)
  years <- format(tmvec, "%Y")
  
  findMax <- function(subIndxs) {
    subIndxMaxIndx <- which.max(srs[subIndxs])[1]
    subIndxs[subIndxMaxIndx]
  }
  annualMaxIndx <- tapply(1:length(srs), years, findMax)
  annualMaxIndx <- annualMaxIndx[!is.na(annualMaxIndx)]
  annualMax <- srs[annualMaxIndx]
  annualMaxDate <- timeStamps[annualMaxIndx]
  
  return(list(annualMax = annualMax, annualMaxDate = annualMaxDate, annualMaxIndx = annualMaxIndx))
}

#' Compute Annual Minimum Values from Time Series Data
#'
#' @description
#' The `computeAnnualMinima` function takes a two-column matrix or data frame 
#' with time stamps and corresponding time series data. It computes the annual 
#' maximum values for each year present in the time series.
#'
#' @param timeAndSeries A two-column matrix or data frame where the first column 
#'   contains time stamps and the second column contains the time series data.
#'
#' @return A list containing three elements: 'annualMin' which is a vector of 
#'   the annual minimum values, 'annualMinDate' which is a vector of the dates 
#'   on which these minima occur, and 'annualMinIndx' which is a vector of the 
#'   indices in the original time series that correspond to the annual minima.
#'
#' @examples
#' # Assume 'timeAndSeries' is a matrix with POSIXct time stamps in the first column
#' # and corresponding data values in the second column:
#' result <- computeAnnualMinima(timeAndSeries)
#' # The result is a list containing the annual minimum values and their dates.
#'
#' @export
#'
#' @seealso
#' `tapply` for the function used to apply a function over subsets of data,
#' `which.min` for finding the index of the minimum value in a vector.
#'
computeAnnualMinima <- function(timeAndSeries) {
  timeStamps <- timeAndSeries[,1]
  srs <- timeAndSeries[,2]
  
  tmvec <- as.Date(timeStamps)
  years <- format(tmvec, "%Y")
  
  findMin <- function(subIndxs) {
    subIndxMinIndx <- which.min(srs[subIndxs])[1]
    subIndxs[subIndxMinIndx]
  }
  annualMinIndx <- tapply(1:length(srs), years, findMin)
  annualMinIndx <- annualMinIndx[!is.na(annualMinIndx)]
  annualMin <- srs[annualMinIndx]
  annualMinDate <- timeStamps[annualMinIndx]
  
  return(list(annualMin = annualMin, annualMinDate = annualMinDate, annualMinIndx = annualMinIndx))
}

#' Calculate Empirical Discharge Distribution and Return Periods
#'
#' @description
#' The `empdis` function computes empirical discharge distributions and various
#' return period estimates for a given hydrological time series. It applies
#' different ranking formulas (Hazen, Cunnane, Gumbel) to estimate the return
#' periods associated with observed discharges.
#'
#' @param x A numeric vector representing the discharge time series.
#' @param nyr The number of years represented in the time series data.
#'
#' @return A data frame containing the empirical return periods (emp.RP),
#'   Hazen return periods (haz.RP), Cunnane return periods (cun.RP), Gumbel
#'   transformed values (gumbel), empirical distribution function values (emp.f),
#'   Hazen distribution function values (emp.hazen), Cunnane distribution function
#'   values (emp.cunnane), sorted discharges (Q), and associated timestamps
#'   (timestamp) for the input discharges.
#'
#' @examples
#' # Assuming 'discharges' is a numeric vector of discharge values and
#' # 'num_years' represents the number of years in the time series:
#' discharge_distributions <- empdis(discharges, num_years)
#' # The result is a data frame with empirical distribution and return period estimates.
#'
#' @export
#'
#' @seealso
#' Distribution and ranking techniques used in hydrology for estimating
#' return periods and analyzing frequency.
#'
empdis<-function(x,nyr){
  ts=seq(1:length(x))
  ts=as.data.frame(ts[order(x)])
  x=as.data.frame(x[order(x)])
  epyp=length(x[,1])/nyr
  rank <- apply(x, 2, rank, na.last = "keep")
  hazen<-(rank-0.5)/length(x[,1])
  cunnane<-(rank-0.4)/(length(x[,1])+0.2)
  gumbel=-log(-log(hazen))
  n=length(x$`x[order(x)]`)
  rpc=(1 / (1 - (1:n) / (n + 1)))/12
  nasmp <- apply(x, 2, function(x) sum(!is.na(x)))
  epdp <-  rank/rep(nasmp + 1, each = nrow(rank))
  #plot(rpc,1/(epyp*(1-epdp)))
  empip=data.frame(1/(epyp*(1-epdp)),1/(epyp*(1-hazen)),1/(epyp*(1-cunnane)),gumbel,epdp,hazen,cunnane,x,ts)
  names(empip)=c("emp.RP","haz.RP","cun.RP","gumbel","emp.f","emp.hazen","emp.cunnane",'Q',"timestamp")
  return(empip)
}
empdisl<-function(x,nyr){
  x=as.data.frame(x[order(x,decreasing = T)])
  epyp=length(x[,1])/nyr
  rank <- apply(-x, 2, rank, na.last = "keep")
  hazen<-(rank-0.5)/length(x[,1])
  gumbel=-log(-log(hazen))
  nasmp <- apply(x, 2, function(x) sum(!is.na(x)))
  epdp <-  rank/rep(nasmp + 1, each = nrow(rank))
  empip=data.frame(1/(epyp*(1-epdp)),1/(epyp*(1-hazen)),gumbel,epdp,hazen,x)
  names(empip)=c("emp.RP","haz.RP","gumbel","emp.f","emp.hazen",'Q')
  return(empip)
}
