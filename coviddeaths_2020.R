############### STEP 1: PRE_PROCESSING

# load most the packages required
load.lib<-c("sf", "dplyr", "pscl", "ggplot2", 
            "spdep", "lctools", "stars", "mapview", "leafsync")
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
crs <- sp::CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 
               +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 
               +units=m +no_defs ") # proj4string of coords

# make the SpatialPointsDataFrame object
no2_sf <- SpatialPointsDataFrame(coords = coords, 
                                 data = data, 
                                 proj4string = crs) %>%
          st_as_sf(no2_spdf) %>%
          st_transform(27700) 

# subset to IZ boundaries
no2_sf<-no2_sf[IZ,]
library(stars)
no2_rst<-st_rasterize(no2_sf, dx = 1000, dy = 1000) 
plot(no2_sf)
#devtools::install_github("michaeldorman/geobgu")
library(geobgu)
add<-as.data.frame(raster_extract(no2_rst, IZcoviddeaths, fun = mean, na.rm = TRUE))
IZcoviddeaths<-cbind(IZcoviddeaths, add)
IZcoviddeaths<-IZcoviddeaths%>%dplyr::mutate(no2=as.numeric(as.character(V1)))
library(mapview)
mapview(IZcoviddeaths, zcol="no2")

# income simd for IZ's
dztoiz<-read.csv("data/lookupdatazonetointermediatezone.csv")
simd<-read.csv("data/SIMD2020indicators.csv")

# define the cleaner function
removespecial<-function(x){stringr::str_replace_all(x, "[^[:alnum:]]", " ")}

# cleaning pipe
dztoiz_simd<-simd %>%
  # apply the cleaning function to the variables
  mutate_at(c("Income_rate", "SMR", "Total_population"), removespecial) %>%
  # select only the data we need
  select(Data_Zone, SMR, Income_rate, Total_population)%>%
  # join with the lookup table
  left_join(., dztoiz, by=c("Data_Zone" = "DZ")) %>%
  # group by IZ
  dplyr::group_by(IZcode) %>%
  # make the variables numeric so we can do the mean calc
  mutate_at(c("Income_rate", "SMR", "Total_population"), funs(as.numeric)) %>%
  # now take the mean of the DZ by IZ
  dplyr::summarise_at(c("Income_rate", "SMR", "Total_population"), mean, na.rm = TRUE) 

# merge the NO2 covid data to the simd data
IZcoviddeaths_ap_simd<-left_join(IZcoviddeaths, dztoiz_simd, by=c("InterZone"="IZcode"))
# lets restrict to complete data (lose 16)
IZcoviddeaths_ap_simd<-IZcoviddeaths_ap_simd[!is.na(IZcoviddeaths_ap_simd$no2),]

############### STEP 2: BUILD NAIVE MODELS

library(pscl)
# negative binomial not scaled
nb <- zeroinfl(Number.of.Deaths ~ no2 + Income_rate + Total_population, 
               data = IZcoviddeaths_ap_simd, 
               dist = "negbin")
summary(nb)

# t<-MASS::glm.nb(formula =  Number.of.Deaths ~ no2 + Income_rate + Total_population, 
# data = IZcoviddeaths_ap_simd, 
# link = log)

# negative binomial scaled
scalednb <- zeroinfl(Number.of.Deaths ~ scale(no2) + scale(Income_rate) + scale(Total_population), 
                     data = IZcoviddeaths_ap_simd, 
                     dist = "negbin")
summary(scalednb)

# results
terms <- c('scale(no2)', 'scale(Income_rate)', 'scale(Total_population)')
coefs <- coef(summary(scalednb))$count
outputcoefs<-function(temrs, coefs, sename, modelname){
row <- row.names(coefs) %in% terms
lower <- coefs[row,'Estimate'] - 1.96*coefs[row, sename]
upper <- coefs[row,'Estimate'] + 1.96*coefs[row, sename]
df<-data.frame(
           modelname = modelname,
           terms = terms,
           IRR = exp(coefs[row,'Estimate']),
           lower = exp(lower),
           upper = exp(upper))
return(df)
}
naivecoefs<-outputcoefs(terms, coefs, "Std. Error", "scalednb")

############### STEP 3: TEST SPATIAL AUTOCORRELATION

# get the residuals from the first model so we can IZcoviddeaths_ap_simd the spatial autocorrelation
IZcoviddeaths_ap_simd$residuals<- residuals(nb)
library(spdep)
W.nb<-poly2nb(IZcoviddeaths_ap_simd, row.names =IZcoviddeaths_ap_simd$InterZone)
#saveRDS(W.nb,"data/weightmatrix3.rds")
W.nb<-readRDS("data/weightmatrix3.rds")
W<- nb2mat(W.nb, style="B", zero.policy = T)
W.list<-nb2listw(W.nb, style="B", zero.policy = T)
set.seed(1234) 
moranresults<-moran.mc(x=IZcoviddeaths_ap_simd$residuals, listw=W.list, nsim=10000, zero.policy = T)
moranresults
# IZcoviddeaths_ap_simd using different method
#lm.moranIZcoviddeaths_ap_simd(t, listw=W.list, zero.policy=T)
moran.plot(IZcoviddeaths_ap_simd$Number.of.Deaths, W.list, zero.policy=T)
plot(moranresults, main="", las=1)

############### STEP 4: BUILD SPATIAL MODELS

# spaMM 
library(spaMM)
row.names(W) <- NULL

# Is there an effect of long term air pollution on risk of covid death?
spatialmodel <- fitme(Number.of.Deaths ~ scale(no2) + scale(Income_rate)+ scale(Total_population)
            + adjacency(1|InterZone),
            adjMatrix = W,
            data = IZcoviddeaths_ap_simd, family = 'negbin')
summary(spatialmodel)

# results
terms <- c('scale(no2)', 'scale(Income_rate)', 'scale(Total_population)')
coefs <- as.data.frame(summary(spatialmodel)$beta_table)

# Calculate the IRR (incidence rate ratio)
spatialcoefs<-outputcoefs(terms, coefs, "Cond. SE", "scalednbCAR")

############### STEP 5: COMPARE MODELS

# fitted 
IZcoviddeaths_ap_simd$fitted_car <- fitted(spatialmodel)
IZcoviddeaths_ap_simd$fitted_naive <- fitted(scalednb)

# scatterplot fitted
library(ggplot2)
require(gridExtra)
p1<-ggplot() + geom_point(aes(IZcoviddeaths_ap_simd$fitted_naive, 
                              IZcoviddeaths_ap_simd$Number.of.Deaths))
p2<-ggplot() + geom_point(aes(IZcoviddeaths_ap_simd$fitted_car, 
                              IZcoviddeaths_ap_simd$Number.of.Deaths))
grid.arrange(p1, p2, ncol=2)

# spatial fitted
GLA <-IZcoviddeaths_ap_simd %>% 
  filter(Local.Authority=="Glasgow City")

library(mapview)
# IZcoviddeaths_ap_simd
m1 <- mapview(GLA, zcol = "Number.of.Deaths", legend=F)
m2 <- mapview(GLA, zcol = "fitted_naive", legend=F)
m3 <- mapview(GLA, zcol = "fitted_car", legend=F)
library(leafsync)
sync(m1,m2,m3)# compare coefficients
allcoefs<-rbind(naivecoefs, spatialcoefs)
allcoefs$terms<-as.factor(allcoefs$terms)

require(ggplot2)
ggplot(allcoefs, aes(x = modelname, y = IRR, color=terms)) +
  geom_point(position=position_dodge(width=0.9),size = 4) +
  geom_errorbar(position=position_dodge(width=0.9),
                aes(ymax = upper, ymin = lower))+
                ylim(0.9, 1.5) + 
                geom_hline(yintercept = 1,
                           linetype="dashed")

############### STEP 6: RESULTS INTERPRETATION

# RESEARCH Q1
# Is there an effect of existing inequalities on risk of covid death?
airpolresult <- fitme(Number.of.Deaths ~ no2 + scale(Income_rate)+ scale(Total_population)
                        + adjacency(1|InterZone),
                        adjMatrix = W,
                        data = IZcoviddeaths_ap_simd, family = 'negbin')
summary(airpolresult)
terms <- c('no2')
coefs <- as.data.frame(summary(airpolresult)$beta_table)
airpolcoefs<-outputcoefs(terms, coefs, "Cond. SE", "airpol")

# RESEARCH Q2
# Is there an effect of existing inequalities on risk of covid death?

IZcoviddeaths_ap_simd$SMR10<-IZcoviddeaths_ap_simd$SMR/10
syndemicresult <- fitme(Number.of.Deaths ~ SMR10 + scale(Income_rate)+ scale(Total_population)
            + adjacency(1|InterZone),
            adjMatrix = W,
            data = IZcoviddeaths_ap_simd, family = 'negbin')
summary(syndemicresult)
terms <- c('SMR10')
coefs <- as.data.frame(summary(syndemicresult)$beta_table)
syndemiccoefs<-outputcoefs(terms, coefs, "Cond. SE", "syndemic")

# RESEARCH Q3
# Is there a multiplicative effect of existing inequalities and long term air pollution on risk of covid death?
syndemicairpolresult <- fitme(Number.of.Deaths ~ scale(SMR)*scale(no2) + scale(Income_rate)+ scale(Total_population)
                        + adjacency(1|InterZone),
                        adjMatrix = W,
                        data = IZcoviddeaths_ap_simd, family = 'negbin')
summary(syndemicairpolresult)

terms <- c('scale(SMR):scale(no2)')
coefs <- as.data.frame(summary(syndemicairpolresult)$beta_table)
syndemicairpolcoefs<-outputcoefs(terms, coefs, "Cond. SE", "syndemicairpol")




