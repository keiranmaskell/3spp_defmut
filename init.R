
#default_wd <- "~/Dropbox/tripartite-symbiosis/new_code"
default_wd <- "/Users/keiranmaskell/Desktop/3spp_defmut_code/3spp_defmut"

#define working direcory
if(readline("Is this your first time running this? (y/n)")=='y'){
  wd <- readline("Specify working directory filepath")
  setwd(wd)
}else{
  message("Using default working directory")
  setwd(default_wd)
}

#audio message to play upon error
sorry_dave <- function(){
   system("afplay ./resources/sorry_dave.mp3
", 
   intern=FALSE)  
}

#the list of R packages needed for this project
list_of_packages <- c("matrixStats",
                       "beepr","splines2")

#check to see if packages are installed, if not, install them
lapply(list_of_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
})

#check to see if packages successfully loaded (and installed)
#and return the names of any packages that failed to load
loaded_packages <- lapply(list_of_packages, require,
      character.only = TRUE, quietly = TRUE)

if (all(unlist(loaded_packages))) {
  message("Packages successfully loaded")
} else {
  message(sprintf("The following failed to load: %s",list_of_packages[loaded_packages==FALSE]))
  sorry_dave()
  stop("One or more packages failed to load")
}


#source all the R files needed for this project
source('population_function.R')
source('base_prms.R')
source('meta_functions.R')
source('plot_code.R')
 
#initialize a population
make.pop <- function(prms){
  pop.init <- array(rep(prms$init.pop.size, each=prms$npatch), dim = c(prms$npatch, 3))
  colnames(pop.init) <- c('f', 'm', 'r')
  pop.init <- sample(1:10, size=length(pop.init), replace=TRUE) +
    pop.init
}

#make directories for the output files and figures
if(!dir.exists("./output")){
  system("mkdir output")
}
if(!dir.exists("./figs")){
  system("mkdir figs")
}
