
#define working direcory
if(readline("Use default wd (search for '3spp_defmut-main' from Downloads folder)? (y/n)")=='y'){
  
  message("Using default working directory")
  start_path <- "~/Downloads"
  target_dir_name <- "3spp_defmut-main"
  dirs <- list.dirs(path=start_path, recursive=TRUE)
  found_dir <- names(unlist(sapply(dirs, grep, pattern=target_dir_name)))[1]
  setwd(found_dir)

}else{
  wd <- readline("Specify working directory filepath")
  setwd(wd)
}

#audio message to play upon error
sorry_dave <- function(){
   system("afplay ./resources/sorry_dave.mp3
", 
   intern=FALSE)  
}

#the list of R packages needed for this project
list_of_packages <- c("matrixStats",
                       "beepr","splines2",
                      "readxl","writexl")

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
source('population_functions.R')
source('base_prms.R')
source('meta_functions.R')
source('plotting.R')
 
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
