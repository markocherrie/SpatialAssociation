# Read in the data
coviddeaths<-read.csv("data/ScotlandCovidDeaths.csv")
library(sf)
IZ<-read_sf("boundaries/IZ/SG_IntermediateZone_Bdry_2011.shp")
IZcoviddeaths<-merge(IZ, coviddeaths, by.x="InterZone", by.y="Intermediate.Zone.code")
st_crs(IZcoviddeaths)
lookup<-read.csv("data/IZtoLAlookup.csv")
IZcoviddeaths<-merge(IZcoviddeaths, lookup, by.x="InterZone", by.y="IntZone")
IZcoviddeaths<-subset(IZcoviddeaths, CAName%in%c("Glasgow City"))

# read in the data
library(readr)
library(sf)
iz<-read_sf("boundaries/IZ/SG_IntermediateZone_Bdry_2011.shp") 
covid<-read_csv("data/ScotlandCovidDeaths.csv")
names(covid)[1]<-"InterZone"


