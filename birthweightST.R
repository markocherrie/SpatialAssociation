############### read in birthweight data

# read in the data
library(readr)
library(sf)
dz<-read_sf("boundaries/DZ/SG_DataZoneBdry_2011/SG_DataZone_Bdry_2011.shp") 
simd<-read_csv("data/SIMD2020indicators.csv")
# read in the BW data
BW<-read.csv("data/low-birthweight.csv", skip=7)

#merge 
BWdz<-merge(dz, BW, by.x="DataZone", by.y="Reference.Area", all.x=T)
