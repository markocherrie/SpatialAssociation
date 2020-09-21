############### STEP 1: PRE_PROCESSING

# load all the packages required
load.lib<-c("sf", "dplyr", "pscl", "spdep", "lctools")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

### Read in the data
# Old data
#coviddeaths<-read.csv("data/ScotlandCovidDeaths.csv")
# New data
coviddeaths<-read.csv("data/covid-deaths-data-week-37_Table 11.csv", skip=2)[-1,1:6]

# spatial data
IZ<-read_sf("boundaries/IZ/SG_IntermediateZone_Bdry_2011.shp")
IZcoviddeaths<-merge(IZ, coviddeaths, by.x="InterZone", by.y="Intermediate.Zone.code")
st_crs(IZcoviddeaths)

### Read in the air pollution data

# air pollution from https://uk-air.defra.gov.uk/data/pcm-data
no2<-read.csv("data/mapno22018.csv", skip = 5)
no2$no22018[no2$no22018=="MISSING"]<-NA
no2$no22018<-as.numeric(as.character(no2$no22018))
coords <- no2[ , c("x", "y")]   # coordinates
data <- as.data.frame(no2[ , 4]) # data
crs <- sp::CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs ") # proj4string of coords

# make the SpatialPointsDataFrame object
no2_sf <- SpatialPointsDataFrame(coords = coords, data = data, proj4string = crs) %>%
          st_as_sf(no2_spdf) %>%
          st_transform(27700) 

# attach NO2 points to IZ
no2_in_IZ <- st_join(no2_sf, IZcoviddeaths, join = st_intersects) %>%
             tidyr::drop_na(InterZone) %>%
             as.data.frame(.) %>%
             mutate(no2=as.numeric(as.character(`no2[, 4]`)))

# calculate the mean of NO2 for each IZ
no2_in_IZclean<-no2_in_IZ %>% 
  # group by IZ
  dplyr::group_by(InterZone) %>%
  # take mean of NO2 by IZ
  dplyr::summarize(no2_clean=mean(no2, na.rm=T))%>%
  # join it to the covid deaths dataset
  left_join(IZcoviddeaths, ., by="InterZone")

# income simd for IZ's
dztoiz<-read.csv("data/lookupdatazonetointermediatezone.csv")
simd<-read.csv("data/SIMD2020indicators.csv")

# define the cleaner function
removespecial<-function(x){stringr::str_replace_all(x, "[^[:alnum:]]", " ")}

# cleaning pipe
dztoiz_simd<-simd %>%
  # apply the cleaning function to the variables
  mutate_at(c("Income_rate", "SMR"), removespecial) %>%
  # select only the data we need
  select(Data_Zone, SMR, Income_rate)%>%
  # join with the lookup table
  left_join(., dztoiz, by=c("Data_Zone" = "DZ")) %>%
  # group by IZ
  dplyr::group_by(IZcode) %>%
  # make the variables numeric so we can do the mean calc
  mutate_at(c("Income_rate", "SMR"), funs(as.numeric)) %>%
  # now take the mean of the DZ by IZß
  dplyr::mutate_at(c("Income_rate", "SMR"), mean, na.rm = TRUE) 

# merge the NO2 covid data to the simd data
IZcoviddeaths_ap_simd<-left_join(no2_in_IZclean, iz_simd, by=c("InterZone"="IZcode"))
# lets restrict to complete data
IZcoviddeaths_ap_simd<-IZcoviddeaths_ap_simd[!is.na(IZcoviddeaths_ap_simd$no2_clean),]

############### STEP 2: MODELLING

library(pscl)
# negative binomial not scaled
nb <- zeroinfl(Number.of.Deaths ~ no2_clean + Income_rate, data = IZcoviddeaths_ap_simd, dist = "negbin")
summary(nb)

# negative binomial scaled
scalednb <- zeroinfl(Number.of.Deaths ~ scale(no2_clean) + scale(Income_rate), data = IZcoviddeaths_ap_simd, dist = "negbin")
summary(scalednb)

# poisson model scaled
scaledpoisson<-glm(formula=Number.of.Deaths ~ scale(no2_clean) + scale(Income_rate), family="poisson", data=IZcoviddeaths_ap_simd)
summary(scaledpoisson)

# get the residuals from the first model so we can test the spatial autocorrelation
IZcoviddeaths_ap_simd$residuals<- residuals(nb)
library(spdep)
#W.nb<-poly2nb(IZcoviddeaths_ap_simd, row.names =IZcoviddeaths_ap_simd$InterZone)
W.nb<-readRDS("data/weightmatrix3.rds")
W<- nb2mat(W.nb, style="B", zero.policy = T)
W.list<-nb2listw(W.nb, style="B", zero.policy = T)
set.seed(1234) 
moranresults<-moran.mc(x=IZcoviddeaths_ap_simd$residuals, listw=W.list, nsim=10000, zero.policy = T)
moran.plot(IZcoviddeaths_ap_simd$Number.of.Deaths, W.list, zero.policy=T)
plot(moranresults, main="", las=1)

# spaMM 
library(spaMM)
row.names(W) <- NULL
test<-IZcoviddeaths_ap_simd%>%
  select(InterZone,Number.of.Deaths, no2_clean, Income_rate, SMR)
test$geometry<-NULL

m3 <- fitme(Number.of.Deaths ~ scale(no2_clean) + scale(Income_rate) + adjacency(1|InterZone),
                     adjMatrix = W,
                     data = test, family = 'poisson')

# run the spatial model
remotes::install_github("gregmacfarlane/sppois")
library(sppois)

# https://linkinghub.elsevier.com/retrieve/pii/S0166046210000207
spmodel<-sarpoisson(Number.of.Deaths ~ scale(no2_clean) + scale(Income_rate), data = IZcoviddeaths_ap_simd,
           listw = W.list, method = "fiml")

# NO2 0.17082 in spatial model
# NO2 0.19638 in non-spatial model

moranresults<-moran.mc(x=spmodel$residuals, listw=W.list, nsim=10000, zero.policy = T)



library(lctools)
IZcoviddeaths_centroids<-st_coordinates(st_centroid(IZcoviddeaths_ap_simd))
IZcoviddeaths_ap_simd$Number.of.Deaths<-as.numeric(IZcoviddeaths_ap_simd$Number.of.Deaths)
IZcoviddeaths_ap_simd$X<-IZcoviddeaths_centroids[,1]
IZcoviddeaths_ap_simd$Y<-IZcoviddeaths_centroids[,2]

IZcoviddeaths_ap_simd<-IZcoviddeaths_ap_simd%>%
                      select(Number.of.Deaths, no2_clean, Income_rate, SMR, X, Y)
IZcoviddeaths_ap_simd$geometry<-NULL
IZcoviddeaths_ap_simd$Number.of.Deaths<-as.integer(IZcoviddeaths_ap_simd$Number.of.Deaths)


test<-IZcoviddeaths_ap_simd[101:200,]
test$geometry<-NULL
test$NOD<-test$Number.of.Deaths

gw.zip.bw <- gw.zi.bw(Number.of.Deaths ~ scale(no2_clean) + scale(Income_rate), 
                      "poisson", IZcoviddeaths_ap_simd, cbind(IZcoviddeaths_ap_simd$X,IZcoviddeaths_ap_simd$Y),
                      kernel = 'adaptive', b.min = 5, b.max=30)

gw.zip <- gw.zi(NOD ~ scale(no2_clean) + scale(Income_rate), 
                "poisson", test, 30, 
                kernel = 'adaptive',
                cbind(test$X,test$Y))



# NO2 0.16929 in alternative spatial model


# estimate for scaled air pollution in GW regression = 0.22528
# estimate for scaled air pollution in naive regression = 0.22528


# The main effect of such violations is that the Error SS (Sum of Squares) is underestimated (Davis, 1986 ) thus inflating the value of test statistic. 
# An inflated test statistic increases the chance of a Type I error (Incorrect rejection of a Null Hypothesis). 
# Most GIS provide tools to measure the level of spatial autocorrelation (e.g. Moran's I).











