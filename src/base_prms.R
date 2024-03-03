
make.prms <- function(

                      #intrinsic growth rates
                      r.f= 0.6,
                      r.forage=0,

                      #intrinsic death rates
                      g.f = 0.01,
                      g.m = 0.01,
                      g.r = 0.01,

                      #density-dependent death rates
                      d.f = 0.001,
                      d.m = 0.01,
                      d.r = 0.01,

                      #sigmoid parameter A
                      a.fm = 0,
                      a.mf = 0,
                      a.rm = 0,
                      a.mr = 0,
                      v.v = 0,
                      b.fr = 0.1,

                      #sigmoid parameter B
                      theta.mf = 1,
                      theta.fm = 1,
                      theta.v = 1,
                      theta.r = 1,
                    
                      #sigmoid parameter C
                      #mu.r = 0.0015,
                      mu.mf = 0.08,
                      mu.fm = 0.08,
                      mu.v = 0.05,
                      mu.r = -0.05,

                      #sigmoid parameter D
                      sigma.r = 0,
                      sigma.mf = 0,
                      sigma.fm = 0,
                      sigma.v = 0,
            
                      #predator efficiency
                      h.fr = 0.8,

                      #conversion rates
                      conv.ftom = 0.3,
                      conv.ftor = 0.3,
                      
                      #dispersal rates
                      frac.f.disperse = 0.01,
                      frac.m.disperse = 0.001,
                      frac.r.disperse = 0.001,

                      #other parameters
                      init.pop.size = c(f = 1000, m = 10, r = 10),
                      #Holling type of predators
                      raider.type=3,
                      #number of patches
                      npatch = 1,
                      #dispersal function type
                      dispersal.case = 'cloud',
                      #mercenary invasion timepoint
                      intrdate = 1,
                      #switch for whether to print the population matrix over time as the simulation runs
                      #print.pop=T,
                      print.pop=F,
                      #total number of generations (timepoints)
                      ngens = 4000){
  

  inputs <- as.list(environment())
  prms <- prms.update(inputs)
  prms[order(names(prms))] ## order for consistency
}

prms.update <- function(prms, rand=FALSE){

  if(rand==TRUE){
    
    
   #intrinsic growth rates
   prms <- make.prms(
                      r.f= 0.6,
                      #expt 1 ~3
                      #r.forage=0,
                      #expt 4
                      r.forage=runif(1,0,0.6),

                      #intrinsic death rates
                      g.f = 0.01,
                      g.m = 0.01,
                      g.r = 0.01,

                      #density-dependent death rates
                      d.f = 0.01,
                      d.m = 0.01,
                      d.r = 0.01,

                      #sigmoid parameter A
                      a.fm = runif(1,0,1),
                      a.mf = runif(1,-1,1),
                      a.rm = 0,
                      a.mr = 0,
                      v.v = runif(1,0,1),
                      b.fr = runif(1,0,1),

                      #sigmoid parameter B
                      theta.mf = runif(1,0.5,1),
                      theta.fm = runif(1,0.5,1),
                      theta.v = runif(1,0.5,1),
                      theta.r = runif(1,0.5,1),
                    
                      #sigmoid parameter C
                      #mu.r = 0.0015,
                      mu.mf = runif(1,-1,1),
                      mu.fm = runif(1,-1,1),
                      mu.v = runif(1,-1,1),
                      mu.r = runif(1,-1,0),

                      #sigmoid parameter D
                      sigma.r = runif(1,0,0.3),
                      sigma.mf = runif(1,-0.3,0.3),
                      sigma.fm = runif(1,-0.3,0.3),
                      sigma.v = runif(1,0,0.3),
            
                      #predator efficiency
                      h.fr = runif(1,0,1),

                      #conversion rates
                      conv.ftom = 0.3,
                      conv.ftor = 0.3,
                      
                      #dispersal rates
                      frac.f.disperse = 0.01,
                      frac.m.disperse = 0.001,
                      frac.r.disperse = 0.001,

                      #other parameters
                      init.pop.size = c(f = 1000, m = 10, r = 10),
                      #Holling type of predators
                      raider.type=3,
                      #number of patches
                      npatch = 1,
                      #dispersal function type
                      dispersal.case = 'cloud',
                      #mercenary invasion timepoint
                      intrdate = 1,
                      #switch for whether to print the population matrix over time as the simulation runs
                      #print.pop=T,
                      print.pop=F,
                      #total number of generations (timepoints)
                      ngens = 4000
   )
  
    

  }
    
return(prms)
}

