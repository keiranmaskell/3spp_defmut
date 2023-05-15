
run.seq.sim.new <- function(pop, prms, fec, set_seedlist) {
  
  if (sum(prms$v.v,prms$sigma.v) > prms$r.f){
    stop(message="V < R ; choose v.v & sigma.v such that sum(v.v,sigma.v) < r.f")
  }
  
  #set.seed(prms$seed)
  ## res <- vector(mode = 'list', length=2)
  ## names(res) <- c()
  
  res <- array(NA, dim=c(prms$ngens, 6, prms$npatch))
  dimnames(res) <- list(1:prms$ngens,
                        c('freq.f', 'freq.m', 'freq.r',
                          'num.m', 'num.f', 'num.r'),
                        1:prms$npatch)
  names(dimnames(res)) <- c('gen', 'output', 'patch')
  
  
  #03/19 add fec.res
  #should these be just fecundities or fecundities -intrinsic deaths?
  if(fec == TRUE){
    fec.res <- array(NA, dim=c((prms$ngens-1),4,prms$npatch))
    dimnames(fec.res) <- list(1:(prms$ngens-1),
                              c('fec.m', 'fec.f', 'fec.r','fec.rx'),
                              1:prms$npatch)
    names(dimnames(fec.res)) <- c('gen', 'output', 'patch')
  }
  ##
  
  gencounter <- 0
  
  #instead of using a simple seed_list and repeating the seed for each time the
  #seed needs to be reset during a given generation, use nested lists to 
  #hold sets of unique seeds to use for each of the resets, then index this 
  #later to reproduce the operations during the calculation of costs
  #old: simple seed_list
  #seed_list <- round(runif(prms$ngens,1,prms$ngens))
  #new: nested lists of seeds
  #each seed_block is the set of seeds used in a given generation
  #length(seed_metalist) = ngens
  if(typeof(set_seedlist) =='list'){
    seed_metalist <- set_seedlist
  }else{
    seed_metalist <- list()
    for(i in 1:prms$ngens){ seed_block <- list(round(runif(19,1,prms$ngens)))
    seed_metalist <- append(seed_metalist, seed_block)
    }
  }
  
  #print(seed_metalist)
  
  for(ii in 1:prms$ngens) {
    
    gencounter <- gencounter +1
    set.seed(seed_metalist[[ii]][1])
    
    res[ii, 'freq.f',] <- pop[ ,'f']/sum(pop[,'f'])
    res[ii, 'freq.m',] <- pop[ ,'m']/sum(pop[,'m'])
    res[ii, 'freq.r',] <- pop[ ,'r']/sum(pop[,'r'])
    
    res[ii,'num.f',] <- pop[,'f']
    res[ii,'num.m',] <- pop[,'m']
    res[ii,'num.r',] <- pop[,'r']
    
    ## farmer reproduction
    #
    
    #print(seed_metalist[[ii]])
    replist <- repp(pp=pop, prms=prms,sid=seed_metalist[[ii]])
    
    #03/19 fec.res
    if(fec == TRUE){
      if(gencounter < prms$ngens){
        fec.res[ii, 'fec.f',] <- pmax(replist[[1]],0)
      }
    }
    #
    pop[,'f'] <- pmax(pop[,'f'] + replist[[1]],0)
    
    
    ## predation and raider foraging; predation() output is vector of farmer death and others' fec
    pred <- pred(pp=pop, prms=prms, vm=replist[[2]], sid=seed_metalist[[ii]])
    #
    #03/19 fec.res
    if(fec == TRUE){
      if(gencounter < prms$ngens){
        fec.res[ii, 'fec.m',] <- pmax(pred[,'m.fec'],0)
        fec.res[ii, 'fec.r',] <- pmax(pred[,'r.fec'],0)
        fec.res[ii, 'fec.rx',] <- pmax(pred[,'r.fec.x'],0)
      }
    }
    #
    pop[,'f'] <- pmax(pop[,'f'] - pred[,'f.eaten'],0)
    pop[,'m'] <- pmax(pop[,'m'] + pred[,'m.fec'],0)
    pop[,'r'] <- pmax(pop[,'r'] + pred[,'r.fec'] + pred[,'r.fec.x'],0)
    ## deathrates
    dead <- rip(pp=pop, prms, sid=seed_metalist[[ii]]) 
    pop <- pmax(pop - dead,0)
    
    #pop[,'f'] <- pmax(pop[,'f'] - dead[,'f.rip'], 0)
    #pop[,'m'] <- pmax(pop[,'m'] - dead[,'r.rip'], 0)
    #pop[,'r'] <- pmax(pop[,'r'] - dead[,'r.rip'], 0)
    
    ## if(ii == 181) browser()
    ##dispersal
    if(prms$dispersal.case=='cloud')
      pop.after.disperse <- disperse.cloud.each(pp=pop, prms=prms, sid=seed_metalist[[ii]])
    
    if(prms$dispersal.case=='adjacent')
      pop.after.disperse <- disperse.adjacent(pp=pop)
    
    pop <- pop.after.disperse
    
    if(prms$print.pop==T){
      print(ii)
      print(pop)
    }    
  }
  
  if(prms$print.pop==T){
    print("final pop")
    print(pop)
  }
  
  if(fec==TRUE){
    reslist <- list(res, fec.res, seed_metalist)
  }else{
    reslist <- list(res, seed_metalist)
  }  
  return(reslist)
  
}





quicksave_prms <- function(prms){
  get_prms_names <- c(as.character(names(prms)))
  get_prms <-c(as.character(prms))
  prmsthing <- data.frame(get_prms,row.names = get_prms_names)
  write.table(prmsthing,sprintf("/Users/keiranmaskell/Desktop/3spp_defmut/vecs/paramvec_%s.txt",Sys.time()))
  
}

data.grab <- function(res.data, prms, filepath, timetag){
  
  gens <- length(dimnames(res.data)$gen)
  
  total.num.m <- sum(res.data[prms$ngens,'num.m',])
  total.num.f <- sum(res.data[prms$ngens,'num.f',])
  total.num.r <- sum(res.data[prms$ngens,'num.r',])
  
  if( total.num.m == 0){
    categorize.m <- "invasion failure"
  }else{
    categorize.m <- "invasion success"
  }
  if(total.num.f == 0){
    categorize.f <- "farmers extinct"
  }else{
    categorize.f <- "farmers survive"
  }
  if( total.num.r == 0){
    categorize.r <- "raiders extinct"
  }else{
    categorize.r <- "raiders survive"
  }
  #could make more outcome codes based on combinations or other scheme
  if(all(total.num.m!=0,total.num.f!=0,total.num.r!=0)){
    outcome_code <- 'coex'
  }else{
    if(all(total.num.m!=0,total.num.f!=0)){
      outcome_code <- 'mercs_win'
    }
    if(all(total.num.r!=0,total.num.f!=0)){
      outcome_code <- 'raiders_win'
    }
    if(all(total.num.m==0,total.num.f==0,total.num.r==0))
      outcome_code <- 'total_collapse'
    else{
      outcome_code <- 'failure'
    }
  }
  outcomes <- c(categorize.m=categorize.m,
                categorize.f=categorize.f,
                categorize.r=categorize.r)
  
  #run.data <- list(res=res, outcomes=outcomes, prms=prms)
  run.data <- list(res=res.data, prms=prms)
  
  
  #system(paste0('cd output'))
  #system(paste0('pwd'))
  #system(paste0('ls'))
  
  #outcome_code <- 1
  #scenario <- 14
  #timetag <- gsub(' ','_', gsub(':','',Sys.time()))
  
  #fp <- sprintf('output/%s_scenario%s_%s_%s',Sys.Date(),scenario,outcome_code,timetag)
  #system(paste0('mkdir ',fp,''))
  save(run.data, file=sprintf('%s/ejr_run_%s%s_%s.RData', filepath, outcome_code, gens, timetag))
}



#cutoff should either be false or an integer (row == generation index)
calculate_fitness <- function(res_data,species,cutoff){
  
  if(cutoff == FALSE){
    run_array <- as.array(res_data)
  }else{
    run_array <- as.array(res_data)
    run_array <- run_array[cutoff:nrow(res_data), ,]
  }
  #run_array <- as.array(run_env$run.data$res)
  
  
  #run_array[,'num.m',]
  metapop <- apply(run_array[,species,], MARGIN = 1, FUN = sum)
  
  rl <- rle(c(metapop)
              == 0)
  num_zero_runs <- sum(rl$values & rl$lengths >= 1)
  metapop <- replace(metapop,metapop==0,NA)
  metapop <- na.omit(metapop)

  #Calculate mean metapopulation size, counting stretches of zeroes as single zeroes
  ifelse(length(rep(0, num_zero_runs))==0,meanval<-mean(metapop),meanval<-mean(metapop,rep(0, num_zero_runs)))
  
  #Calculate stdev of metapopulation size, counting stretches of zeroes as single zeroes
  ifelse(length(rep(0, num_zero_runs))==0,stdev<-sd(metapop),stdev<-sd(metapop,rep(0, num_zero_runs)))
  
  fitness_list <- list(meanval,stdev)
  return(fitness_list)
}

fitness_calculators <- function(res_data,cutoff_list){
  
  farmer_fitness <- calculate_fitness(res_data,'num.f',cutoff=cutoff_list[1])
  merc_fitness <- calculate_fitness(res_data,'num.m',cutoff=cutoff_list[2])
  raider_fitness <- calculate_fitness(res_data,'num.r',cutoff=cutoff_list[3])
  
  dfout <- data.frame(c(farmer_fitness,merc_fitness,raider_fitness))
  colnames(dfout) <- c('farmer mean metapop','farmer metapop stdev','merc mean metapop','merc metapop stdev','raider mean metapop','raider metapop stdev')
  return(dfout)
  
}





#calculate costs


#04/26/2023 index the populations in by generation in res, and index the seed for each operation
#using the seedlist
total_cost_matrix <- function(res.data, prms, seedlist, filepath, timetag, testing){
  
  #define the arrays which will hold growth matrices, modified growth matrices, 
  #and total cost matrices
  
  #define the mercenary total cost array
  #tcostarr_m <- array(NA, dim=c(prms$ngens-1, 3, prms$npatch))
  #dimnames(tcostarr_m) <- list(1:(prms$ngens-1),
  #                      c('benefit_m_to_f','benefit_m_to_m','benefit_m_to_r'),
  #                      1:prms$npatch)
  #names(dimnames(tcostarr_m)) <- c('gen', 'costs_of_m', 'patch')
  
  #define the raider total cost array
  #tcost_m_r <- array(NA, dim=c(prms$ngens-1, 3, prms$npatch))
  #dimnames(tcost_m_r) <- list(1:(prms$ngens-1),
  #                      c('benefit_r_to_f','benefit_r_to_m','benefit_r_to_r'),
  #                      1:prms$npatch)
  #names(dimnames(tcost_m_r)) <- c('gen', 'costs_of_r', 'patch')
  
  #define the growth array
  growtharr <- array(NA, dim=c(prms$ngens-1, 3, prms$npatch))
  dimnames(growtharr) <- list(1:(prms$ngens-1),
                              c('farmer_growth',
                                'merc_growth',
                                'raider_growth'),
                              1:prms$npatch)
  names(dimnames(growtharr)) <- c('gen', 'growth', 'patch')
  
  #define the modified growth array where mercs are nerfed
  growtharr_nerfmercs <- array(NA, dim=c(prms$ngens-1, 3, prms$npatch))
  dimnames(growtharr_nerfmercs) <- list(1:(prms$ngens-1),
                                        c('farmer_mod_growth',
                                          'merc_mod_growth',
                                          'raider_mod_growth'),
                                        1:prms$npatch)
  names(dimnames(growtharr_nerfmercs)) <- c('gen', 'modified growth', 'patch')
  
  #define the modified growth array where raiders are nerfed
  growtharr_nerfraiders <- array(NA, dim=c(prms$ngens-1, 3, prms$npatch))
  dimnames(growtharr_nerfraiders) <- list(1:(prms$ngens-1),
                                          c('farmer_mod_growth','merc_mod_growth','raider_mod_growth'),
                                          1:prms$npatch)
  names(dimnames(growtharr_nerfraiders)) <- c('gen', 'modified_growth', 'patch')
  
  
  #populate the growth array with values from res
  #count <- 0
  for(i in seq(1,(prms$ngens-1),1)){
    #count <- count +1
    #print(count)
    growtharr[i, 'farmer_growth',] <- res.data[(i+1),'num.f',] - res.data[i,'num.f',]
    growtharr[i, 'merc_growth',] <- res.data[(i+1),'num.m',] - res.data[i,'num.m',]
    growtharr[i, 'raider_growth',] <- res.data[(i+1),'num.r',] - res.data[i,'num.r',]
  }
  
  #populate the modified growth arrays
  #populate the modified growth array where mercs are nerfed
  
  #test
  #test_res to compare populations with focal species
  if(testing == TRUE){
    test_res <- array(rep(0, each=prms$npatch), dim = c(prms$npatch, 3, 2))
    dimnames(test_res) <- list(1:prms$npatch,c("f","m","r"),c("n(t)","n(t+1)"))
  }
  
  
  temp_res <- array(rep(0, each=prms$npatch), dim = c(prms$npatch, 3, 2))
  dimnames(temp_res) <- list(1:prms$npatch,c("f","m","r"),c("n(t)","n(t+1)"))
  for(i in seq(1,(prms$ngens-1),1)){
    
    #test
    #for(i in seq(1,(100),1)){
    if(testing==TRUE){
      print(i)
    }
    
    #print(res.data[i,'num.r',])
    #print(res.data[(i+1),'num.r',])
    
    
    # #count <- count +1
    # #print(count)
    temp_res[,"f",1] <- res.data[i,'num.f',]
    # #temp_res[,"m",1] <- res.data[i,'num.m',]
    temp_res[,"m",1] <- 0
    temp_res[,"r",1] <- res.data[i,'num.r',]
    
    temp_res[,"f",2] <- res.data[i,'num.f',]
    # #temp_res[,"m",1] <- res.data[i,'num.m',]
    temp_res[,"m",2] <- 0
    temp_res[,"r",2] <- res.data[i,'num.r',]
    # #dim(temp_res)
    #print(temp_res)
    
    #test
    if(testing == TRUE){
      test_res[,"f",1] <- res.data[i,'num.f',]
      test_res[,"m",1] <- res.data[i,'num.m',]
      test_res[,"r",1] <- res.data[i,'num.r',]
      
      test_res[,"f",2] <- res.data[(i+1),'num.f',]
      test_res[,"m",2] <- res.data[(i+1),'num.m',]
      test_res[,"r",2] <- res.data[(i+1),'num.r',]
    }
    
    # 
    #print(test_res)
    # #match seed
    set.seed(seedlist[[i]][1])
    # 
    # ## farmer reproduction
    reprod.f.temp <- repp(pp=temp_res[,,2], prms=prms,sid=seedlist[[i]])
    # update farmers
    temp_res[,'f',2] <- pmax(temp_res[,'f',2] + reprod.f.temp[[1]],0)
    # 
    # #load mercs and raiders
    # temp_res[,'m',2] <- temp_res[,'m',1]
    # temp_res[,'r',2] <- temp_res[,'r',1]
    # 
    ## predation; output is vector of farmer death and others' fec
    pred.temp <- pred(pp=temp_res[,,2], prms=prms, vm=reprod.f.temp[[2]],sid=seedlist[[i]])
    #update farmers
    temp_res[,'f',2] <- pmax(temp_res[,'f',2] - pred.temp[,'f.eaten'], 0)
    #update mercs
    temp_res[,'m',2] <- pmax(temp_res[,'m',2] + pred.temp[,'m.fec'], 0)
    # 
    #update raiders
    temp_res[,'r',2] <- pmax(temp_res[,'r',2] + pred.temp[,'r.fec'] + pred.temp[,'r.fec.x'], 0)
    # 
    # 
    ## deathrates
    dead.temp <- rip(pp=temp_res[,,2], prms,sid=seedlist[[i]]) 
    # 
    ## Updated population sizes
    temp_res[,,2] <- pmax(temp_res[,,2] - dead.temp,0)
    # 
    if(prms$dispersal.case=='cloud'){
      pop.after.disperse <- disperse.cloud.each(pp=temp_res[,,2], prms=prms,sid=seedlist[[i]])}
    # 
    if(prms$dispersal.case=='adjacent'){
      pop.after.disperse <- disperse.adjacent(pp=temp_res[,,2])}
    # 
    temp_res[,,2] <- pop.after.disperse
    # 
    
    if(testing == TRUE){
      print("temp_res")
      print(temp_res)
      print("test_res")
      print(test_res)
      print("difference in growth: mercs - no mercs")
      print(temp_res - test_res)
    }
    
    
    
    growtharr_nerfmercs[i, 'farmer_mod_growth',] <- temp_res[,"f",2] - temp_res[,"f",1]
    growtharr_nerfmercs[i, 'merc_mod_growth',] <- temp_res[,"m",2] - temp_res[,"m",1]
    growtharr_nerfmercs[i, 'raider_mod_growth',] <- temp_res[,"r",2] - temp_res[,"r",1]
    
    if(testing == TRUE){
      print(growtharr[i,,])
      print(growtharr_nerfmercs[i,,])
      print(growtharr[i,,] - growtharr_nerfmercs[i,,])
    }
  }
  
  
  #populate the modified growth array where raiders are nerfed
  temp_res <- array(rep(0, each=prms$npatch), dim = c(prms$npatch, 3, 2))
  dimnames(temp_res) <- list(1:prms$npatch,c("f","m","r"),c('n(t)','n(t+1)'))
  for(i in seq(1,(prms$ngens-1),1)){
    #count <- count +1
    #print(count)
    temp_res[,"f",1] <- res.data[i,'num.f',]
    temp_res[,"m",1] <- res.data[i,'num.m',]
    #temp_res[,"r",1] <- res.data[i,'num.r',]
    temp_res[,"r",1] <- 0
    #dim(temp_res)
    temp_res[,"f",2] <- res.data[i,'num.f',]
    temp_res[,"m",2] <- res.data[i,'num.m',]
    #temp_res[,"r",1] <- res.data[i,'num.r',]
    temp_res[,"r",2] <- 0
    
    #match seed
    set.seed(seedlist[[i]])
    
    ## farmer reproduction
    reprod.f.temp <- repp(pp=temp_res[,,2], prms=prms,sid=seedlist[[i]])
    # update farmers
    temp_res[,'f',2] <- pmax(temp_res[,'f',2] + reprod.f.temp[[1]],0)
    
    #load mercs and raiders
    #temp_res[,'m',2] <- temp_res[,'m',1]
    #temp_res[,'r',2] <- temp_res[,'r',1]
    
    ## predation; output is vector of farmer death and others' fec
    pred.temp <- pred(pp=temp_res[,,2], prms=prms, vm=reprod.f.temp[[2]],sid=seedlist[[i]])
    
    #update farmers
    temp_res[,'f',2] <- pmax(temp_res[,'f',2] - pred.temp[,'f.eaten'], 0)
    #update mercs
    temp_res[,'m',2] <- pmax(temp_res[,'m',2] + pred.temp[,'m.fec'], 0)
    
    #update raiders
    temp_res[,'r',2] <- pmax(temp_res[,'r',2] + pred.temp[,'r.fec'] + pred.temp[,'r.fec.x'], 0)
    
    
    ## deathrates
    dead.temp <- rip(pp=temp_res[,,2], prms,sid=seedlist[[i]]) 
    
    ## Updated population sizes
    # temp_res[,'f',2] <- pmax(temp_res[,'f',2] -
    #                            dead.temp[,'f.rip'], 0)
    # 
    # temp_res[,'m',2] <- pmax(temp_res[,'m',2] - dead.temp[,'m.rip'], 0)
    # 
    # temp_res[,'r',2] <- pmax(temp_res[,'r',1] -
    #                            dead.temp[,'r.rip'], 0)
    
    
    ## Updated population sizes
    temp_res[,,2] <- pmax(temp_res[,,2] - dead.temp,0)
    
    if(prms$dispersal.case=='cloud'){
      pop.after.disperse <- disperse.cloud.each(pp=temp_res[,,2], prms=prms,sid=seedlist[[i]])}
    
    if(prms$dispersal.case=='adjacent'){
      pop.after.disperse <- disperse.adjacent(pp=temp_res[,,2])}
    
    temp_res[,,2] <- pop.after.disperse
    
    
    growtharr_nerfraiders[i, 'farmer_mod_growth',] <- temp_res[,"f",2] - temp_res[,"f",1]
    growtharr_nerfraiders[i, 'merc_mod_growth',] <- temp_res[,"m",2] - temp_res[,"m",1]
    growtharr_nerfraiders[i, 'raider_mod_growth',] <- temp_res[,"r",2] - temp_res[,"r",1]
    
  }
  #define & populate the merc total cost array
  
  tcostarr_m <- growtharr - growtharr_nerfmercs
  
  #growtharr[1999,,]
  #growtharr_nerfmercs[1999,,]
  #growtharr[1999,,] - growtharr_nerfmercs[1999,,]
  
  dimnames(tcostarr_m) <- list(1:(prms$ngens-1),
                               c('benefit_m_to_f','benefit_m_to_m','benefit_m_to_r'),
                               1:prms$npatch)
  names(dimnames(tcostarr_m)) <- c('gen', 'costs_of_m', 'patch')
  
  
  #define & populate the raider total cost array
  tcostarr_r <- growtharr - growtharr_nerfraiders
  #tcost_m_r <- array(NA, dim=c(prms$ngens-1, 3, prms$npatch))
  dimnames(tcostarr_r) <- list(1:(prms$ngens-1),
                               c('benefit_r_to_f','benefit_r_to_m','benefit_r_to_r'),
                               1:prms$npatch)
  names(dimnames(tcostarr_r)) <- c('gen', 'costs_of_r', 'patch')
  
  #save data
  total_cost.data <- list(res=res.data, growth_array = growtharr, m_totalcost_array = tcostarr_m,
                          r_totalcost_array = tcostarr_r, prms=prms)
  save(total_cost.data, file=sprintf('%s/kdm_total_costs_%s.RData',filepath, timetag))
}

#04/26/2023 index the populations in by generation in res, and index the seed for each operation
#using the seedlist
component_cost_matrix <- function(res.data, prms, curr_prms, seedlist, filepath, timetag, testing, ...){
  
  #get parameters for measuring effect in a list
  get_comps <-  list(...)
  
  if("r.f" %in% get_comps){
    print("if r.f is zero, so will be v.m")
    get_comps <- append(get_comps,'v.v')
  }
  #print(1)
  #print(get_comps)
  #print(get_comps[[1]])
  #print(c(unlist(get_comps)))
  #print(typeof(unlist(get_comps)))
  #setup string for filename 
  filnam <- paste(get_comps,collapse='_')
  #print(filnam)
  
  #get parameters which are already changed from default settings
  
  ##testing only
  #get_comps <- list("mu.r","sigma.fm")
  #curr_prms <- make.prms(mu.r=0.5,sigma.fm=0.9)
  #prms <- make.prms()
  ##
  
  #just names
  curr_diff_prms_names <- names(curr_prms[as.character(curr_prms)!=as.character(prms)])
  #print(c(curr_diff_prms_names))
  
  #names and values
  #curr_diff_prms <- as.list(c(names(Sys.getenv(curr_prms[as.character(curr_prms)!=as.character(prms)]))))
  curr_diff_prms <- as.list(c(curr_prms[as.character(curr_prms)!=as.character(prms)]))
  #print(curr_diff_prms)
  
  #reformat list of parameters for effect measurement
  input_for_make.prms <- paste(get_comps,collapse='=0,')
  
  #see if there are any doubles and remove the versions belonging to the non-default parameter list
  #print(2)
  #print(typeof(curr_diff_prms_names))
  #print(curr_diff_prms_names[curr_diff_prms_names %in% unlist(get_comps)])
  #print(3)
  #print(c(curr_diff_prms_names) %in% c(unlist(get_comps)))
  
  curr_diff_prms_names <- curr_diff_prms_names[!(curr_diff_prms_names %in% unlist(get_comps))]
  #print(curr_diff_prms_names)
  curr_diff_prms <- curr_diff_prms[!(names(curr_diff_prms) %in% unlist(get_comps))]
  #print(curr_diff_prms)
  
  #reformat list of parameters which are already changed from default settings for passing to 
  #the make.prms() function
  #curr_input_for_make.prms <- sprintf('%s=%s',paste(as.list(curr_diff_prms_names),collapse=sprintf('=%s,',curr_diff_prms[1:length(curr_diff_prms)])),curr_diff_prms[length(curr_diff_prms)])
  curr_input_for_make.prms <- sprintf('%s=%s',paste(as.list(curr_diff_prms_names),collapse='=TOBEREPLACED,'),curr_diff_prms[length(curr_diff_prms)])
  #lapply(curr_diff_prms[1:length(curr_diff_prms)-1],sprintf, fmt='=%s')
  #teststr2 <- sub('TOBEREPLACED',curr_diff_prms,curr_input_for_make.prms)
  #teststr2 <- sapply(curr_diff_prms,sub('TOBEREPLACED',curr_diff_prms,curr_input_for_make.prms))
  #library(stringr)
  #teststr2 <- str_replace_all(curr_input_for_make.prms,'TOBEREPLACED',c(curr_diff_prms))
  #curr_input_for_make.prms <- paste(as.list(curr_diff_prms_names),collapse='=0,')
  
  #slow method of replacing the separator with the non-default prm values corresponding to their names 
  for(i in curr_diff_prms[1:length(curr_diff_prms)-1]){
    #print(curr_diff_prms[i])
    #print(i)
    
    curr_input_for_make.prms <- sub('TOBEREPLACED',i,curr_input_for_make.prms)
    
    #print(input_for_make.prms)
  }
  
  
  #print('fed in')
  #tester
  #paste(list(testlist),collapse=sprintf('=%s,',list(empty2)))
  
  #instantiate make.prms() with all the parameters that are non-default or for which cost is to be measured
  #mod_prms <- (eval(parse(text = paste((sprintf('make.prms(%s=0,%s)',input_for_make.prms,curr_input_for_make.prms))))))
  
  if(length(curr_input_for_make.prms) == 0){
    print("These parameters will be evaluated for cost/benefit")
    print(paste(get_comps,collapse=','))
    print("No other parameters are modified from the default")
    mod_prms <- (eval(parse(text = paste((sprintf('make.prms(%s=0)',input_for_make.prms))))))
    
  }else{
    print("These parameters will be evaluated for cost/benefit")
    print(paste(get_comps,collapse=','))
    print("These parameters are also modified from the default")
    print(curr_input_for_make.prms)
    mod_prms <- (eval(parse(text = paste((sprintf('make.prms(%s=0,%s)',input_for_make.prms,curr_input_for_make.prms))))))
  }
  
  
  #paste((sprintf('make.prms(%s=0,%s)',input_for_make.prms,curr_input_for_make.prms)))
  #mod_prms <- (eval(parse(text = paste((sprintf('make.prms(%s=0)',input_for_make.prms))))))
  
  #print(mod_prms)
  
  
  
  
  #print(mod_prms)
  #print(paste((sprintf('make.prms(%s=0,%s)',input_for_make.prms,curr_input_for_make.prms))))
  
  
  
  #A - B = target C
  #define arrays 
  #define the growth array (A)
  growtharr <- array(NA, dim=c(curr_prms$ngens-1, 3, curr_prms$npatch))
  dimnames(growtharr) <- list(1:(curr_prms$ngens-1),
                              c('farmer_growth',
                                'merc_growth',
                                'raider_growth'),
                              1:curr_prms$npatch)
  names(dimnames(growtharr)) <- c('gen', 'growth', 'patch')
  
  #define the modified growth array where comps are nerfed (B)
  growtharr_nerfcomps <- array(NA, dim=c(curr_prms$ngens-1, 3, curr_prms$npatch))
  dimnames(growtharr_nerfcomps) <- list(1:(curr_prms$ngens-1),
                                        c('farmer_mod_growth',
                                          'merc_mod_growth',
                                          'raider_mod_growth'),
                                        1:curr_prms$npatch)
  names(dimnames(growtharr_nerfcomps)) <- c('gen', 'modified growth', 'patch')
  
  
  #populate growth array (A)
  for(i in seq(1,(curr_prms$ngens-1),1)){
    #count <- count +1
    #print(count)
    growtharr[i, 'farmer_growth',] <- res.data[(i+1),'num.f',] - res.data[i,'num.f',]
    growtharr[i, 'merc_growth',] <- res.data[(i+1),'num.m',] - res.data[i,'num.m',]
    growtharr[i, 'raider_growth',] <- res.data[(i+1),'num.r',] - res.data[i,'num.r',]
  }
  
  #test
  #test_res to compare populations with focal species
  if(testing == TRUE){
    test_res <- array(rep(0, each=prms$npatch), dim = c(prms$npatch, 3, 2))
    dimnames(test_res) <- list(1:prms$npatch,c("f","m","r"),c("n(t)","n(t+1)"))
  }
  
  #populate modified growth array (B)
  temp_res <- array(rep(0, each=mod_prms$npatch), dim = c(mod_prms$npatch, 3, 2))
  dimnames(temp_res) <- list(1:mod_prms$npatch,c("f","m","r"),1:2)
  for(i in seq(1,(mod_prms$ngens-1),1)){
    
    #test
    #for(i in seq(1,(100),1)){
    if(testing==TRUE){
      print(i)
    }
    
    #count <- count +1
    #print(count)
    temp_res[,"f",1] <- res.data[i,'num.f',]
    temp_res[,"m",1] <- res.data[i,'num.m',]
    temp_res[,"r",1] <- res.data[i,'num.r',]
    #dim(temp_res)
    
    temp_res[,"f",2] <- res.data[i,'num.f',]
    temp_res[,"m",2] <- res.data[i,'num.m',]
    temp_res[,"r",2] <- res.data[i,'num.r',]
    # #dim(temp_res)
    #print(temp_res)
    
    #test
    if(testing == TRUE){
      test_res[,"f",1] <- res.data[i,'num.f',]
      test_res[,"m",1] <- res.data[i,'num.m',]
      test_res[,"r",1] <- res.data[i,'num.r',]
      
      test_res[,"f",2] <- res.data[(i+1),'num.f',]
      test_res[,"m",2] <- res.data[(i+1),'num.m',]
      test_res[,"r",2] <- res.data[(i+1),'num.r',]
    }
    
    
    
    
    #match seed
    #set.seed(as.integer(seedlist[i]))
    set.seed(seedlist[[i]][1])
    
    ## farmer reproduction
    reprod.f.temp <- repp(pp=temp_res[,,1], prms=mod_prms,sid=seedlist[[i]])
    # update farmers
    temp_res[,'f',2] <- pmax(temp_res[,'f',1] + reprod.f.temp[[1]],0)
    
    #load mercs and raiders
    #temp_res[,'m',2] <- temp_res[,'m',1]
    #temp_res[,'r',2] <- temp_res[,'r',1]
    
    ## predation; output is vector of farmer death and others' fec
    pred.temp <- pred(pp=temp_res[,,2], prms=mod_prms, vm=reprod.f.temp[[2]],sid=seedlist[[i]])
    
    #update farmers
    temp_res[,'f',2] <- pmax(temp_res[,'f',2] - pred.temp[,'f.eaten'], 0)
    #update mercs
    temp_res[,'m',2] <- pmax(temp_res[,'m',2] + pred.temp[,'m.fec'], 0)
    
    #update raiders
    temp_res[,'r',2] <- pmax(temp_res[,'r',2] + pred.temp[,'r.fec'] + pred.temp[,'r.fec.x'], 0)
    
    
    ## deathrates
    dead.temp <- rip(pp=temp_res[,,2], mod_prms,sid=seedlist[[i]]) 
    
    # ## Updated population sizes
    temp_res[,,2] <- pmax(temp_res[,,2] - dead.temp,0)
    # temp_res[,'f',2] <- pmax(temp_res[,'f',2] -
    #                            dead.temp[,'f.rip'], 0)
    # 
    # temp_res[,'m',2] <- pmax(temp_res[,'m',2] - dead.temp[,'m.rip'], 0)
    # 
    # temp_res[,'r',2] <- pmax(temp_res[,'r',1] -
    #                            dead.temp[,'r.rip'], 0)
    
    
    if(prms$dispersal.case=='cloud'){
      pop.after.disperse <- disperse.cloud.each(pp=temp_res[,,2], prms=mod_prms,sid=seedlist[[i]])}
    
    if(prms$dispersal.case=='adjacent'){
      pop.after.disperse <- disperse.adjacent(pp=temp_res[,,2])}
    
    temp_res[,,2] <- pop.after.disperse
    
    
    if(testing == TRUE){
      print("temp_res")
      print(temp_res)
      print("test_res")
      print(test_res)
      print("difference")
      print(temp_res - test_res)
    }
    
    growtharr_nerfcomps[i, 'farmer_mod_growth',] <- temp_res[,"f",2] - temp_res[,"f",1]
    growtharr_nerfcomps[i, 'merc_mod_growth',] <- temp_res[,"m",2] - temp_res[,"m",1]
    growtharr_nerfcomps[i, 'raider_mod_growth',] <- temp_res[,"r",2] - temp_res[,"r",1]
    
    if(testing == TRUE){
      print(growtharr_nerfcomps[i,,])
      print(growtharr[i,,] - growtharr_nerfcomps[i,,])
    }
    
  }
  if(testing == TRUE){
    print("MOD PRMS")
    print(mod_prms)
    #print("")
    #print()
  }
  
  #####
  #define & populate the component cost array
  #target C = A - B
  c_costarr <- growtharr - growtharr_nerfcomps
  
  dimnames(c_costarr) <- list(1:(curr_prms$ngens-1),
                              c('benefit_to_f','benefit_to_m','benefit_to_r'),
                              1:curr_prms$npatch)
  names(dimnames(c_costarr)) <- c('gen', 'component_costs', 'patch')
  
  #save data
  component_cost.data <- list(res=res.data, growth_array = growtharr, component_cost_array = c_costarr, original_prms=curr_prms, modified_prms=mod_prms)
  save(component_cost.data, file=sprintf('%s/kdm_component_costs_%s%s.RData', filepath,filnam,timetag))
}





collect_data <- function(fp,master_df,exp,rep,scenario,prms,datenow,plot_yn){
  

experiment <- exp
repnum <- as.character(rep)
scenario <-scenario
#prms <- make.prms(ngens=4000,mu.r=0.05,sigma.m=0,b.fm=0,sigma.v=0.1,sigma.mf=-0.2,sigma.fm=-0.2,r.forage=0.05)

#collect data about prms
#Invasion gen
invasion_date <- prms$intrdate
#Parasitism inc or dec
if(sign(prms$v.v) == 1){
  parasitism_incdec <- '+'
}else{if(sign(prms$v.v == 0)){
  parasitism_incdec <- 'c'
}else{
  parasitism_incdec <- '-'
}
}
#Parasitoid effects y/n
if(prms$b.fm !=0 | prms$sigma.m != 0){
  parasitoid_yn <- 'Y'
}else{
  parasitoid_yn <-'N'
}
#Parasitoid inc or dec
if(parasitoid_yn == 'N'){
  parasitoid_incdec <- NA
}else{
  if(prms$b.fm ==0 | prms$mu.m ==0){
    parasitoid_incdec <- 'c'
  }else{
    if(prms$b.fm >0){
      parasitoid_incdec <- '+'
    }else{
      parasitoid_indec <- '-'
    }
  }
}
#Competition F to M y/n
if(prms$a.mf ==0 & prms$sigma.mf ==0){
  competition_ftom <- 'n'
}else{
  if(prms$a.mf > 0){
    if(prms$a.mf + prms$sigma.mf > 0){
      competition_ftom <- 'y'
    }else{
      competition_ftom <- 'n'
    }
  }else{
    if(prms$a.mf <0){
      if(prms$sigma.mf >0){
        competition_ftom <- 'y'
      }else{
        competition_ftom <- 'n'
      }
    }
  }
}

#Competition M to F y/n
if(prms$a.fm ==0 & prms$sigma.fm ==0){
  competition_mtof <- 'n'
}else{
  if(prms$a.fm > 0){
    if(prms$a.fm + prms$sigma.fm > 0){
      competition_mtof <- 'y'
    }else{
      competition_mtof <- 'n'
    }
  }else{
    if(prms$a.fm <0){
      if(prms$sigma.fm >0){
        competition_mtof <- 'y'
      }else{
        competition_mtof <- 'n'
      }
    }
  }
}

#Competition F to M inc / dec
if(competition_ftom == 'n'){
  competition_ftom_incdec <- NA
}else{
  if(prms$a.mf ==0 | prms$mu.mf ==0){
    competition_ftom_incdec <- 'c'
  }else{
    if(prms$a.mf >0){
      competition_ftom_incdec <- '+'
    }else{
      competition_ftom_incdec <- '-'
    }
  }
}

#Competition M to F inc / dec
if(competition_mtof == 'n'){
  competition_mtof_incdec <- NA
}else{
  if(prms$a.fm ==0 | prms$mu.fm ==0){
    competition_mtof_incdec <- 'c'
  }else{
    if(prms$a.fm >0){
      competition_mtof_incdec <- '+'
    }else{
      competition_mtof_incdec <- '-'
    }
  }
}

#Cooperation F to M y/n
if(prms$a.mf ==0 & prms$sigma.mf ==0){
  cooperation_ftom <- 'n'
}else{
  if(prms$a.mf > 0){
    if(prms$sigma.mf < 0){
      cooperation_ftom <- 'y'
    }else{
      cooperation_ftom <- 'n'
    }
  }else{
    if(prms$a.mf <0){
      if(prms$a.mf + prms$sigma.mf <0){
        cooperation_ftom <- 'y'
      }else{
        cooperation_ftom <- 'n'
      }
    }
  }
}

#Cooperation M to F y/n
if(prms$a.fm ==0 & prms$sigma.fm ==0){
  cooperation_mtof <- 'n'
}else{
  if(prms$a.fm > 0){
    if(prms$sigma.fm < 0){
      cooperation_mtof <- 'y'
    }else{
      cooperation_mtof <- 'n'
    }
  }else{
    if(prms$a.fm <0){
      if(prms$a.fm + prms$sigma.fm <0){
        cooperation_mtof <- 'y'
      }else{
        cooperation_mtof <- 'n'
      }
    }
  }
}

#Cooperation F to M inc / dec
if(cooperation_ftom == 'n'){
  cooperation_ftom_incdec <- NA
}else{
  if(prms$a.mf ==0 | prms$mu.mf ==0){
    cooperation_ftom_incdec <- 'c'
  }else{
    if(prms$a.mf >0){
      cooperation_ftom_incdec <- '-'
    }else{
      cooperation_ftom_incdec <- '+'
    }
  }
}
#Cooperation M to F inc / dec
if(cooperation_mtof == 'n'){
  cooperation_mtof_incdec <- NA
}else{
  if(prms$a.fm ==0 | prms$mu.fm ==0){
    cooperation_mtof_incdec <- 'c'
  }else{
    if(prms$a.fm >0){
      cooperation_mtof_incdec <- '-'
    }else{
      cooperation_mtof_incdec <- '+'
    }
  }
}
#R foraging
if(prms$r.forage !=0){
  r_foraging <- 'y'
}else{
  r_foraging <- 'n'
}


#run the model
pop <- make.pop(prms)

#test pairwise F&R
popfr <- pop
popfr[,'m'] <- 0
#new seedlist created in response to FALSE argument
res <- run.seq.sim.new(popfr, prms, fec=FALSE,set_seedlist=FALSE)
seedlist <- res[2][[1]]
res_fr <- res[1][[1]]
#save data
tt <- gsub(' ','_', gsub(':','',Sys.time()))
data.grab(res_fr,prms,fp,tt)
#get average metapops
rundata <- list.files(fp,pattern=tt)
rundf <- load(sprintf('%s/%s/%s',getwd(),fp,rundata))
rundf_data <- get(rundf)
res_fr <- rundf_data[1]
res_fr <- res_fr[[1]]
fr_metapop_data <- fitness_calculators(res_fr,c(FALSE,FALSE,FALSE))
fr_f_avg_metapop <- fr_metapop_data['farmer mean metapop'][1,1]
fr_f_metapop_stdev <- fr_metapop_data['farmer metapop stdev'][1,1]
fr_r_avg_metapop <- fr_metapop_data['raider mean metapop'][1,1]
fr_r_metapop_stdev <- fr_metapop_data['raider metapop stdev'][1,1]
if(sum(res_fr[prms$ngens,'num.f',])!=0 & sum(res_fr[prms$ngens,'num.r',])!=0){
  pairwise_fr_outcome <- 'coex'
}else{
  if(sum(res_fr[prms$ngens,'num.f',])==0){
    pairwise_fr_outcome <- 'collapse'
  }else{
    pairwise_fr_outcome <- 'failure'
  }
}


#test pairwise F&M but keep the seeds
popfm <- pop
popfm[,'r'] <- 0
res <- run.seq.sim.new(popfm, prms, fec=FALSE,set_seedlist=seedlist)
res_fm <- res[1][[1]]
tt <- gsub(' ','_', gsub(':','',Sys.time()))
data.grab(res_fm,prms,fp,tt)
rundata <- list.files(fp,pattern=tt)
rundf <- load(sprintf('%s/%s/%s',getwd(),fp,rundata))
rundf_data <- get(rundf)
res_fm <- rundf_data[1]
res_fm <- res_fm[[1]]
fm_metapop_data <- fitness_calculators(res_fm,c(FALSE,FALSE,FALSE))
fm_f_avg_metapop <- fm_metapop_data['farmer mean metapop'][1,1]
fm_f_metapop_stdev <- fm_metapop_data['farmer metapop stdev'][1,1]
fm_m_avg_metapop <- fm_metapop_data['merc mean metapop'][1,1]
fm_m_metapop_stdev <- fm_metapop_data['merc metapop stdev'][1,1]
if(sum(res_fm[prms$ngens,'num.f',])!=0 & sum(res_fm[prms$ngens,'num.m',])!=0){
  pairwise_fm_outcome <- 'coex'
}else{
  if(sum(res_fm[prms$ngens,'num.f',])==0){
    pairwise_fm_outcome <- 'collapse'
  }else{
    pairwise_fm_outcome <- 'failure'
  }
}

#test all 3; still same seeds
res <- run.seq.sim.new(pop, prms, fec=FALSE,set_seedlist=seedlist)
res <- res[1][[1]]
tt <- gsub(' ','_', gsub(':','',Sys.time()))
data.grab(res,prms,fp,tt)
rundata <- list.files(fp,pattern=tt)
rundf <- load(sprintf('%s/%s/%s',getwd(),fp,rundata))
rundf_data <- get(rundf)
res_all3 <- rundf_data[1]
res_all3 <- res_all3[[1]]
#res_all3 <- res
all3_metapop_data <- fitness_calculators(res_all3,c(FALSE,FALSE,FALSE))
all3_f_avg_metapop <- all3_metapop_data['farmer mean metapop'][1,1]
all3_f_metapop_stdev <- all3_metapop_data['farmer metapop stdev'][1,1]
all3_m_avg_metapop <- all3_metapop_data['merc mean metapop'][1,1]
all3_m_metapop_stdev <- all3_metapop_data['merc metapop stdev'][1,1]
all3_r_avg_metapop <- all3_metapop_data['raider mean metapop'][1,1]
all3_r_metapop_stdev <- all3_metapop_data['raider metapop stdev'][1,1]


tt <- gsub(' ','_', gsub(':','',Sys.time()))
total_cost_matrix(res_all3, prms, seedlist, fp, tt, testing =FALSE)

component_cost_matrix(res_all3, make.prms(),prms, seedlist, fp, tt, testing =FALSE,"mu.r")
component_cost_matrix(res_all3, make.prms(),prms, seedlist, fp, tt, testing =FALSE,"a.fm","a.mf")
component_cost_matrix(res_all3, make.prms(),prms, seedlist, fp, tt, testing =FALSE,"b.fm","sigma.m")
component_cost_matrix(res_all3, make.prms(),prms, seedlist, fp, tt, testing =FALSE,"v.v","sigma.v")
component_cost_matrix(res_all3, make.prms(),prms, seedlist, fp, tt, testing =FALSE,"b.fr","sigma.r")
component_cost_matrix(res_all3, make.prms(),prms, seedlist, fp, tt, testing =FALSE,"r.forage")



if(sum(res_all3[prms$ngens,'num.f',])!=0 & sum(res_all3[prms$ngens,'num.m',])!=0 & sum(res_all3[prms$ngens,'num.m',]!=0)){
  all3_outcome <- 'coex'
}else{
  if(sum(res_all3[prms$ngens,'num.f',])==0){
    all3_outcome <- 'collapse'
  }else{
    if(sum(res_all3[prms$ngens,'num.m',] & sum(res_all3[prms$ngens,'num.r',])==0)){
      all3_outcome <- 'both failure'
    }else{
      if(sum(res_all3[prms$ngens,'num.m',]==0)){
        all3_outcome <- 'merc failure'
      }else{
        all3_outcome <- 'raider failure'
      }
    }
  }
}



#also dispersal 

#also fitness landscapes 

#also tsne on other params


totaldata_list <- list.files(fp,pattern='total_costs')
compdata_list <- list.files(fp,pattern='component')

#Mercs true mutualists
for(j in totaldata_list){
  filepath_tcost <- sprintf('%s/%s',fp,j)
  dftotal_costs <- load(file.path(sprintf('%s',filepath_tcost)), tcost_env <- new.env() )
  dftotal_costs <- load(file.path(sprintf('%s',filepath_tcost)))
  df_total_costs_data <- get(dftotal_costs)
  
  #growth_array <- df_total_costs_data[2][[1]]
  m_totalcosts <- df_total_costs_data[3][[1]]
  #r_totalcosts <- df_total_costs_data[4][[1]]
  raw_avg_cost_ben_mercs_to_farmers <- mean(m_totalcosts[,'benefit_m_to_f',])
  stdev_cost_ben <- sd(m_totalcosts[,'benefit_m_to_f',])
  
  if(raw_avg_cost_ben_mercs_to_farmers>0){
    mercs_mutualists <- 'y'
  }else{
    mercs_mutualists <- 'n'
  }

}


#Notes
notes <- NA
#Datadir
datadir <- fp

#plotting
#make plot directory
plotfp <- sprintf('figs/%s_scenario%s_%s',datenow,scenario,tt)
system(paste0('mkdir ',plotfp,''))

#Figdir
figdir <- plotfp


data_vec <- c(experiment,repnum,scenario,pairwise_fr_outcome,pairwise_fm_outcome,all3_outcome,fr_f_avg_metapop,fr_r_avg_metapop,fr_f_metapop_stdev,fr_r_metapop_stdev,fm_f_avg_metapop,fm_m_avg_metapop,fm_f_metapop_stdev,fm_m_metapop_stdev,all3_f_avg_metapop,all3_m_avg_metapop,all3_r_avg_metapop,all3_f_metapop_stdev,all3_m_metapop_stdev,all3_r_metapop_stdev,invasion_date,parasitism_incdec,parasitoid_yn,parasitoid_incdec,competition_ftom,competition_mtof,cooperation_ftom,cooperation_mtof,competition_ftom_incdec,competition_mtof_incdec,cooperation_ftom_incdec,cooperation_mtof_incdec,r_foraging,mercs_mutualists,raw_avg_cost_ben_mercs_to_farmers,stdev_cost_ben,notes,datadir,figdir,unlist(prms))

master_df <- readRDS(master_df)

master_transpose <- cbind(t(master_df),data_vec)
#View(t(master_transpose))
master_df <- t(master_transpose)
#works
saveRDS(master_df,sprintf('%s/output/master_df.RDS',getwd()))





#plot populations
for(i in 1:prms$npatch){
  plot.pop.patch(patch=i,prms=prms,res=res,plotfp)
}
#plot metapopulations
plot.pop.all(prms,res,plotfp)




#plot total costs
for(j in totaldata_list){
  filepath_tcost <- sprintf('%s/%s',fp,j)
  #filepath_tcost <- "/Users/keiranmaskell/Dropbox/tripartite-symbiosis/new_code/output/03_14_scenario14_bfr_coexistence_1/kdm_total_costs_2023-03-14 21:15:14.RData"
  dftotal_costs <- load(file.path(sprintf('%s',filepath_tcost)), tcost_env <- new.env() )
  
  #Merc costs
  #plot.pop.tcosts(patch,`prms,tcostarr,species)
  for(i in 1:prms$npatch){
    plot.pop.tcosts(patch=i,prms=prms,tcostarr=tcost_env$total_cost.data$m_totalcost_array,
                    species='Mercenaries',sep='n', plotfp)
  }
  for(i in 1:prms$npatch){
    plot.pop.tcosts(patch=i,prms=prms,tcostarr=tcost_env$total_cost.data$m_totalcost_array,
                    species='Mercenaries',sep='y', plotfp)
  }
  
  #Raider costs
  for(i in 1:prms$npatch){
    plot.pop.tcosts(patch=i,prms=prms,tcostarr=tcost_env$total_cost.data$r_totalcost_array,
                    species='Raiders',sep='n', plotfp)
  }
  for(i in 1:prms$npatch){
    plot.pop.tcosts(patch=i,prms=prms,tcostarr=tcost_env$total_cost.data$r_totalcost_array,
                    species='Raiders',sep='y', plotfp)
  }
  
  
#   
#   #plot growth over gens
#   #plot.pop.growth(patch, prms, growtharr)
#   for(i in 1:prms$npatch){
#     plot.pop.growth(patch=i,prms=prms,growtharr=tcost_env$total_cost.data$growth_array,sep='n', plotfp)
#   }
#   
#   for(i in 1:prms$npatch){
#     plot.pop.growth(patch=i,prms=prms,growtharr=tcost_env$total_cost.data$growth_array,sep='y', plotfp)
#   }
  
  }




#component costs
for(j in compdata_list){
  filepath_ccost <- sprintf('%s/%s',fp,j)
  #print(filepath_ccost)
  dfcomp_costs <- load(file.path(sprintf('%s',filepath_ccost)), ccost_env <- new.env() )
  
  
  # #all together
  # for(i in 1:prms$npatch){
  #   plot.pop.comp.costs(patch=i,prms=prms,ccostarr=ccost_env$component_cost.data$component_cost_array,
  #                       orig.prms=ccost_env$component_cost.data$original_prms,
  #                       mod.prms=ccost_env$component_cost.data$modified_prms,
  #                       sep='n',
  #                       movavg=20,plotfp)
  # }
  #separated
  for(i in 1:prms$npatch){
    plot.pop.comp.costs(patch=i,prms=prms,ccostarr=ccost_env$component_cost.data$component_cost_array,
                        orig.prms=ccost_env$component_cost.data$original_prms,
                        mod.prms=ccost_env$component_cost.data$modified_prms,
                        sep='y',
                        movavg=20,plotfp)
  }
  
  for(i in 1:prms$npatch){
    plot.pop.comp.costs(patch=i,prms=prms,ccostarr=ccost_env$component_cost.data$component_cost_array,
                        orig.prms=ccost_env$component_cost.data$original_prms,
                        mod.prms=ccost_env$component_cost.data$modified_prms,
                        sep='y',
                        movavg=250,plotfp)
  }
  
  # #plot psplines of component cost
  # for(i in 1:prms$npatch){
  #   comp_costplot.splines(cost_fp = filepath_ccost, patch=i,plotfp)
  # }
  
  
}


}

#num_reps should make it repeat the same runs multiple times, recording what rep it is in repnum
run.model <- function(index, num_reps, scen){
  
  if(length(list.files('output',pattern='master_df.RDS'))==0){
    #create master file for data collection
    master_df <- data.frame(matrix(ncol=88,nrow=0))
    #length(c('Experiment','rep num','Scenario','Pairwise F&R','Pairwise F&M','All 3','F&R F avg metapop','F&R avg metapop','F&R F metapop stdev','F&R R metapop stdev','F&M F avg metapop','F&M M avg metapop','F&M F metapop stdev','F&M M metapop stdev','All 3 F avg metapop','All 3 M avg metapop','All 3 R avg metapop','All 3 F metapop stdev','All 3 M metapop stdev','All 3 R metapop stdev','Invasion gen','Parasitism inc or dec','Parasitoid effects','Parasitoid effects inc or dec','Competition F to M','Competition M to F','Cooperation F to M','Cooperation M to F','Competition F to M inc or dec','Competition M to F inc or dec','Cooperation F to M inc or dec','Cooperation M to F inc or dec','R Foraging','Mercs true mutualists','Raw avg cost/ben','Std dev cost/ben','Notes','Datadir','Figdir','a.fm','a.mf','a.mr','a.rm','b.fm','b.fr','conv.ftom','conv.ftor','depredation.case','dispersal.case','forageswitchon','frac.f.disperse','frac.m.disperse','frac.r.disperse','g.f','g.m','g.r','h.fm','h.fr','f_init','m_init','r_init','intrdate','k.f','k.m','k.r','mu.fm','mu.m','mu.mf','mu.r','mu.v','ngens','npatch','print.pop','r.f','r.forage','raider.type','seed','sigma.fm','sigma.m','sigma.mf','sigma.r','sigma.v','theta.fm','theta.m','theta.mf','theta.r','theta.v','v.v'))
    colnames(master_df) <- c('Experiment','rep num','Scenario','Pairwise F&R','Pairwise F&M','All 3','F&R F avg metapop','F&R avg metapop','F&R F metapop stdev','F&R R metapop stdev','F&M F avg metapop','F&M M avg metapop','F&M F metapop stdev','F&M M metapop stdev','All 3 F avg metapop','All 3 M avg metapop','All 3 R avg metapop','All 3 F metapop stdev','All 3 M metapop stdev','All 3 R metapop stdev','Invasion gen','Parasitism inc or dec','Parasitoid effects','Parasitoid effects inc or dec','Competition F to M','Competition M to F','Cooperation F to M','Cooperation M to F','Competition F to M inc or dec','Competition M to F inc or dec','Cooperation F to M inc or dec','Cooperation M to F inc or dec','R Foraging','Mercs true mutualists','Raw avg cost/ben','Std dev cost/ben','Notes','Datadir','Figdir','a.fm','a.mf','a.mr','a.rm','b.fm','b.fr','conv.ftom','conv.ftor','depredation.case','dispersal.case','forageswitchon','frac.f.disperse','frac.m.disperse','frac.r.disperse','g.f','g.m','g.r','h.fm','h.fr','f_init','m_init','r_init','intrdate','k.f','k.m','k.r','mu.fm','mu.m','mu.mf','mu.r','mu.v','ngens','npatch','print.pop','r.f','r.forage','raider.type','seed','sigma.fm','sigma.m','sigma.mf','sigma.r','sigma.v','theta.fm','theta.m','theta.mf','theta.r','theta.v','v.v')
    saveRDS(master_df,sprintf('%s/output/master_df.RDS',getwd()))
  }
  
  #read_master <- readRDS('output/master_df.RDS')
  read_master <- file.path("output",pattern='master_df.RDS')
  #typeof(read_master)

  scenario <- scen
  
  for(xx in 1:nrow(index)){
    #print(xx)
    
    
    prms <- make.prms(ngens=4000,
                      r.forage = index[xx,'r.forage'],
                      mu.r= index[xx,'mu.r'],
                      mu.mf= index[xx,'mu.f'],
                      mu.fm= index[xx,'mu.f'],
                      h.fr= index[xx,'h.fr'],
                      b.fr = 1.0,
                      a.mf = 0.1,
                      a.fm = 0.1,
                      sigma.fm = 0.2,
                      sigma.mf = 0.2,
                      frac.m.disperse=0.01,
                      frac.r.disperse=0.05,
                      b.fm=0.00017,
                      h.fm = 6.0,
                      seed = 6,
                      mu.m = 0.5,
    )
    
    tt <- gsub(' ','_', gsub(':','',Sys.time()))
    datenow <- Sys.Date()
    fp <- sprintf('output/%s_scenario%s_%s',datenow,scenario,tt)
    system(paste0('mkdir ',fp,''))
    
    for(jj in seq(1,num_reps,1)){
      #print(jj)
      # datalist <- run.model(xx)
      # res <- datalist[1][[1]]
      # b <- datalist[2][[1]]
      
      collect_data(fp=fp,master_df=read_master,exp=xx,rep=jj,scenario=NA,prms=prms,datenow=datenow,plot_yn=TRUE)
      
    }
    
  }
  read_master <- readRDS(sprintf('%s/output/master_df.RDS',getwd()))
  View(read_master)
}

#run.model(index=index,num_reps=5,scen=NA)








