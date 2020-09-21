
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
