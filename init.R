


## list_of_packages <- c("openxlsx",
##                       "matrixStats",
##                       "jsonlite")

## loadallpackages <- function(packagelist) {
##     loaded_yn <- lapply(packagelist, require,
##      character.only = TRUE, quietly = TRUE)
##     return(loaded_yn)
## }

## loadallpackages(list_of_packages)
library(matrixStats)
library("beepr")

sorry_dave <- function(){
   system("afplay /Users/keiranmaskell/Dropbox/tripartite-symbiosis/sorry_dave.mp3", 
   intern=FALSE)  
}

source('population_function.R')
source('base_prms.R')
source('meta_functions.R')
source('plot_code.R')
 
make.pop <- function(prms){
  pop.init <- array(rep(prms$init.pop.size, each=prms$npatch), dim = c(prms$npatch, 3))
  colnames(pop.init) <- c('f', 'm', 'r')
  pop.init <- sample(1:10, size=length(pop.init), replace=TRUE) +
    pop.init
}

