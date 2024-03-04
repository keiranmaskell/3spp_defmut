rm(list=ls())

#sets working directory or throws an error if no path is provided in the call
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
    setwd(args[1])
} else {
    stop("SPECIFY MAIN WORKING DIRECTORY (the 3spp_defmut-main dir)")
}

#initialize
source("src/init.R")


#how many parameter randomizations to try
nsample <- 10000

#generate randomized parameter values and store in index
#see base_prms.R for the ranges over which individual parameters are sampled
prms<-make.prms()
index <- lapply(seq_len(nsample), function(i) prms.update(prms, rand = TRUE))

                
#name the experiment using the current date, appended to a random number between 1000 and 9999
test_id <- sprintf("test_%s_rID_%s",Sys.Date(),sample(1000:9999,1))


                
#nreps is the number of replicate runs using a given set of parameter values 
run.model(index=index,nreps=5,expt_id=test_id)



