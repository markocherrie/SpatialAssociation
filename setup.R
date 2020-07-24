#setup

load.lib<-c("readr", "sf", "dplyr", "tmap", "sp", "spdep")
install.lib<-load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib,dependencies=TRUE)
sapply(load.lib,require,character=TRUE)