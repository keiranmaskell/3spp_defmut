
make.prms <- function(k.f = 1000,
                      k.m = 100,
                      k.r = 100,
                      a.fm = -0.5,
                      a.mf = -0.5,
                      a.rm = 0,
                      a.mr = 0,
                      b.fr = -0.3,
                      b.fm = 0.3,
                      mu.m = 0.05,
                      #mu.r = 0.0015,
                      mu.r = 0.05,
                      mu.mf = 0.08,
                      mu.fm = 0.08,
                      mu.v = 0.05,
                      theta.m = 5,
                      theta.r = 5,
                      theta.mf = 5,
                      theta.fm = 5,
                      theta.v = 5,
                      r.f = 0.3,
                      v.v = 0.1,
                      h.fr = 4,
                      h.fm = 6,
                      conv.ftom = 1,
                      conv.ftor = 1,
                      g.f = 0.1,
                      g.m = 0.1,
                      g.r = 0.1,
                      frac.f.disperse = 0.01,
                      frac.m.disperse = 0.001,
                      frac.r.disperse = 0.001,
                      sigma.m = 0.1,
                      sigma.r = 0.3,
                      sigma.mf = 0.2,
                      sigma.fm = 0.2,
                      sigma.v = 0.1,
                      init.pop.size = c(f = 3000, m = 30, r = 10),
                      raider.type=3,
                      npatch = 12,
                      forageswitchon = TRUE,
                      r.forage = 0,
                      dispersal.case = 'cloud',
                      depredation.case = '3',
                      intrdate = 1,
                      print.pop=T,
                      ngens = 2000,
                      seed = 6){
  

  inputs <- as.list(environment())
  prms <- prms.update(inputs)
  prms[order(names(prms))] ## order for consistency
}

prms.update <- function(prms, rand=FALSE){

  if(rand==T){
    mu.m <- runif(1,0.5,4)
    mu.r <- runif(1,0.05,4)
    mu.f <- runif(1,0.5,4)
    sigma.f <- runif(1,0.5,4)
    sigma.m <- runif(1,0.5,4)
    sigma.r <- runif(1,0.5,4)
    h.fr <- runif(1,0.5,4)
    h.fm <- runif(1,0.5,4)
    r.f <- runif(1,0.1,1)
    b.fm <- runif(1,0.1,1)
    b.fr<- runif(1,0.1,1)
    a.fm <- runif(1,0.1,4)
    a.mf <- runif(1,0.1,4)
    r.forage <- runif(1,0.1,1)
    #r.forage <- 0
    #r.forage <- 0.1
    conv.ftom <- runif(1,0.01,0.5)
    conv.ftor <- runif(1,0.01,0.5)
    g.f <- runif(1,0.01,0.5)
    g.m <- runif(1,0.01,0.5)
    g.r <- runif(1,0.01,0.5)
    #g.r <- 0.9
    frac.f.disperse <- runif(1,0.01,0.5)
    frac.m.disperse <- runif(1,0.01,0.5)
    frac.r.disperse <- runif(1,0.01,0.5)
    #omega.m - intrinsic attack rate level coefficient of mercs
    omega.m <- runif(1,0.01,1.0)
    psi.m <- runif(1,0.00001, 0.1)

  }
    
return(prms)
}



#before 0403
# make.prms <- function(k.f = 1000,
#                       k.m = 100,
#                       k.r = 100,
#                       a.fm = 1.6,
#                       a.mf = 0.5,
#                       a.rm = 0,
#                       a.mr = 0,
#                       b.fr = 0.5,
#                       b.fm = 0.001,
#                       mu.m = 0.5,
#                       mu.r = 0.0015,
#                       mu.f = 0.08,
#                       mu.fm = 0.08,
#                       sigma.m = 3,
#                       sigma.r = 0.01,
#                       sigma.f = 1,
#                       sigma.fm = 1,
#                       r.f = 0.3, 
#                       h.fr = 2,
#                       h.fm = 6,
#                       conv.ftom = 1,
#                       conv.ftor = 1,
#                       g.f = 0.1,
#                       g.m = 0.1,
#                       g.r = 0.1,
#                       frac.f.disperse = 0.01,
#                       frac.m.disperse = 0.001,
#                       frac.r.disperse = 0.001,
#                       omega.m = 0.1,
#                       psi.m = 0.005,
#                       sig.f = 1.1,
#                       sig.fm = 1.1,
#                       init.pop.size = c(f = 3000, m = 30, r = 10),
#                       raider.type=3,
#                       npatch = 12,
#                       forageswitchon = TRUE,
#                       r.forage = 0,
#                       dispersal.case = 'cloud',
#                       depredation.case = '3',
#                       intrdate = 1,
#                       print.pop=T,
#                       ngens = 2000,
#                       seed = 6,
#                       sig.m = 0.0001,
#                       v.m = 0.001){

