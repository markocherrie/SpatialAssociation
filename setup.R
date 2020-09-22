# load most the packages required
load.lib<-c("sf", "dplyr", "pscl", "ggplot2", "MASS","devtools",
            "spdep", "lctools", "stars", "mapview", "leafsync")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)

devtools::install_github("michaeldorman/geobgu")