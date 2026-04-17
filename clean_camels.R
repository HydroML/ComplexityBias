rm(list = ls())
gc()

library(Evapotranspiration)
library(zoo)
library(tibble)
library(lubridate)

#get basin attributes
df <- read.table("C:/Users/janssenj/Desktop/research_postdoc/Data/basin_dataset_public_v1p2/basin_metadata/gauge_information.txt",
                 header = TRUE, sep = "\t", quote = "", fill = TRUE,stringsAsFactors = FALSE,strip.white=TRUE, 
                 colClasses = c("character","character","character","numeric","numeric","numeric"))

df$DRAINAGE.AREA..KM.2.<-NULL
df$LAT<-NULL
df$LONG<-NULL
colnames(df)<-c("HUC","ID","Name")


camels_clim <- read.csv("C:/Users/janssenj/Desktop/research_postdoc/Data/basin_dataset_public_v1p2/camels_clim.txt",
                        sep=";", colClasses = c("character", rep(NA, 100000)))

camels_geol <- read.csv("C:/Users/janssenj/Desktop/research_postdoc/Data/basin_dataset_public_v1p2/camels_geol.txt",
                        sep=";", colClasses = c("character", rep(NA, 100000)))

camels_soil <- read.csv("C:/Users/janssenj/Desktop/research_postdoc/Data/basin_dataset_public_v1p2/camels_soil.txt",
                        sep=";", colClasses = c("character", rep(NA, 100000)))

camels_topo <- read.csv("C:/Users/janssenj/Desktop/research_postdoc/Data/basin_dataset_public_v1p2/camels_topo.txt",
                        sep=";", colClasses = c("character", rep(NA, 100000)))

camels_vege <- read.csv("C:/Users/janssenj/Desktop/research_postdoc/Data/basin_dataset_public_v1p2/camels_vege.txt",
                        sep=";", colClasses = c("character", rep(NA, 100000)))

catchment_attr<-merge(df,camels_clim,by.x="ID",by.y = "gauge_id",all.x = T)
catchment_attr<-merge(catchment_attr,camels_geol,by.x="ID",by.y = "gauge_id",all.x = T)
catchment_attr<-merge(catchment_attr,camels_soil,by.x="ID",by.y = "gauge_id",all.x = T)
catchment_attr<-merge(catchment_attr,camels_topo,by.x="ID",by.y = "gauge_id",all.x = T)
catchment_attr<-merge(catchment_attr,camels_vege,by.x="ID",by.y = "gauge_id",all.x = T)

rm(df,camels_clim,camels_geol,camels_soil,camels_topo,camels_vege)

catchment_attr$area_errors<-abs(catchment_attr$area_gages2 - catchment_attr$area_geospa_fabric)/max(catchment_attr$area_gages2,catchment_attr$area_geospa_fabric)



PERFECT_df<-data.frame(date=seq(from=as.Date("1980-01-01"),to=as.Date("2014-12-31"),by="day"))
PERFECT_df$year<-year(PERFECT_df$date)
PERFECT_df$month<-month(PERFECT_df$date)
PERFECT_df$day<-day(PERFECT_df$date)

catchment_attr$missing_prec<-NA
# clean time series, add missing data % to attribute data
for(i in 1:nrow(catchment_attr)){
  dat_force <- read.table(paste0("C:/Users/janssenj/Desktop/research_postdoc/Data/basin_dataset_public_v1p2/basin_mean_forcing/daymet/",
                                 catchment_attr$HUC[i],"/",catchment_attr$ID[i]
                                 ,"_lump_cida_forcing_leap.txt"), skip = 3, header = TRUE)
  
  
  flow_dat <- read.table(paste0("C:/Users/janssenj/Desktop/research_postdoc/Data/basin_dataset_public_v1p2/usgs_streamflow/"
                                ,catchment_attr$HUC[i],"/",catchment_attr$ID[i],"_streamflow_qc.txt"), quote="\"", comment.char="")
  print(paste(nrow(dat_force), nrow(flow_dat)))
  flow_dat$V5[flow_dat$V5<0]<-NA
  flow_dat$V1<-NULL
  colnames(flow_dat)<-c("year","month","day","flow","flag")
  #sum(is.na(flow_dat$flow))
  
  dat_flow<-merge(PERFECT_df,flow_dat,all.x = T,by=c("year","month","day"),sort = F)
  #sum(is.na(flow_dat$flow))
  dat_flow<-dat_flow[order(dat_flow$year, dat_flow$month, dat_flow$day),]
  
  #change streamflow cubic feet per second to mm/day
  dat_flow$flow<-dat_flow$flow*86400 #first change to cubic feet per day
  dat_flow$flow<-dat_flow$flow*(0.0003048)^3 # change cubic feet per day to cubic km per day
  dat_flow$flow<-dat_flow$flow/catchment_attr$area_gages2[i] # change to km per day
  dat_flow$flow<-dat_flow$flow*1000*1000 # change to mm per day
  
  
  
  
  
  
  dates <- as.Date(sprintf("%04d-%02d-%02d", dat_force$Year, dat_force$Mnth, dat_force$Day))
  
  data_for_pet<-list(Date.daily=as.Date(paste0(dat_force$Year,"-",dat_force$Mnth,"-",dat_force$Day)),
                     Tmax=zoo(dat_force$tmax.C., order.by = dates),
                     Tmin= zoo(dat_force$tmin.C., order.by = dates),
                     n=zoo(dat_force$dayl.s./(60*60), order.by = dates))
  
  results <- ET.Hamon(data_for_pet, ts="daily", message="yes", AdditionalStats="yes", save.csv="no")
  
  dat_force$PET<-as.numeric(results$ET.Daily)
  dat_force$tmean<-(dat_force$tmax.C.+ dat_force$tmin.C.)/2

  dat_force$flow<-dat_flow$flow
  dat_force$flow_flag<-dat_flow$flag
  
  write.csv(dat_force,paste0("C:/Users/janssenj/Desktop/research_postdoc/Data/cleaned_camels_ts/",
                              catchment_attr$ID[i],"_",catchment_attr$HUC[i],".csv"),row.names = F)
  
  catchment_attr$missing_prec[i]<-mean(is.na(dat_force$flow))
  
}

write.csv(catchment_attr, "C:/Users/janssenj/Desktop/research_postdoc/Data/catchmet_attr.csv",row.names = F)




