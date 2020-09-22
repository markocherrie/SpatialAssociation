# run below once
library(devtools)
devtools::install_github("michaeldorman/geobgu")

# load all the packages required
load.lib<-c("sf", "dplyr", "pscl", "ggplot2", "MASS","geobgu",
            "spdep", "lctools", "stars", "mapview", "leafsync")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

