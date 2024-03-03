

#this function will run the model with parameters determined by the "index" df
#num_reps should make it repeat the same runs multiple times, recording what rep it is in repnum
#expt_id is a unique tag for the experiment
run.model <- function(index, num_reps, expt_id){
  
  #make a new output df if none exists
  #the alternative is that the current experiment gets appended rowwise to an existing record
  if(length(list.files('output',pattern='master_df.RDS'))==0){
    #create master file for data collection
    
    master_df <- data.frame(matrix(ncol=68,nrow=0))
    
    #replace seed with seedid
    #colnames(master_df) <- c('Experiment','index','rep num','Pairwise F&R','Pairwise F&M','All 3','F&R F avg metapop','F&R avg metapop','F&R F metapop stdev','F&R R metapop stdev','F&M F avg metapop','F&M M avg metapop','F&M F metapop stdev','F&M M metapop stdev','All 3 F avg metapop','All 3 M avg metapop','All 3 R avg metapop','All 3 F metapop stdev','All 3 M metapop stdev','All 3 R metapop stdev','Invasion gen','R Foraging','Mercs true mutualists','Raw avg cost/ben','Std dev cost/ben','Notes','Datadir','Figdir','a.fm','a.mf','a.mr','a.rm','b.fr','conv.ftom','conv.ftor','dispersal.case','frac.f.disperse','frac.m.disperse','frac.r.disperse','g.f','g.m','g.r','h.fr','f_init','m_init','r_init','mu.fm','mu.m','mu.mf','mu.r','mu.v','ngens','npatch','r.f','r.forage','raider.type','seed','sigma.fm','sigma.m','sigma.mf','sigma.r','sigma.v','theta.fm','theta.m','theta.mf','theta.r','theta.v','v.v')
    
    
    colnames(master_df) <- c('Experiment','index','rep num','Pairwise F&R','Pairwise F&M','All 3',
'F&R F avg metapop','F&R R avg metapop','F&R F metapop stdev','F&R R metapop stdev',
'F&M F avg metapop','F&M M avg metapop','F&M F metapop stdev','F&M M metapop stdev',
'All 3 F avg metapop','All 3 M avg metapop','All 3 R avg metapop','All 3 F metapop stdev','All 3 M metapop stdev','All 3 R metapop stdev',
'R Foraging','Mercs true mutualists',
'Raw avg cost/ben','Std dev cost/ben',
'Notes','Datadir','Figdir',

'a.fm','a.mf','a.mr','a.rm','b.fr','conv.ftom','conv.ftor','d.f','d.m','d.r','dispersal.case',
'frac.f.disperse','frac.m.disperse','frac.r.disperse','g.f','g.m','g.r','h.fr','f_init','m_init','r_init', 'invasion gen',
'mu.fm','mu.mf','mu.r','mu.v','ngens','npatch','print.pop',
'r.f','r.forage','raider.type','sigma.fm','sigma.mf','sigma.r','sigma.v',
'theta.fm','theta.mf','theta.r','theta.v','v.v')
    
    
    saveRDS(master_df,sprintf('%s/output/master_df.RDS',getwd()))
  }
  
  read_master <- file.path("output",pattern='master_df.RDS')
  
  
  for(xx in 1:length(index)){
    print(xx)
    #writeLines(as.character(xx),sprintf("%s/%scurrent_index.txt",getwd(),Sys.Date()))
    ind <- xx 
    
    prms <- make.prms(ngens=4000,


                      #intrinsic growth rates
                      r.f= index[[xx]]$r.f,
                      r.forage= index[[xx]]$r.forage,

                      a.fm = index[[xx]]$a.fm,
                      a.mf = index[[xx]]$a.mf,
                      v.v = index[[xx]]$v.v,
                      b.fr = index[[xx]]$b.fr,

                      #a.mr
                      #a.rm

                      theta.fm = index[[xx]]$theta.fm,
                      theta.mf = index[[xx]]$theta.mf,
                      theta.r = index[[xx]]$theta.r,
                      theta.v = index[[xx]]$theta.v,

                      mu.mf= index[[xx]]$mu.mf,
                      mu.fm= index[[xx]]$mu.fm,
                      mu.v= index[[xx]]$mu.v,
                      mu.r= index[[xx]]$mu.r,

                      sigma.fm = index[[xx]]$sigma.fm,
                      sigma.mf = index[[xx]]$sigma.mf,
                      sigma.r = index[[xx]]$sigma.r,
                      sigma.v = index[[xx]]$sigma.v,

                      h.fr= index[[xx]]$h.fr,

                      
                    
                      npatch = index[[xx]]$npatch,

                      #dispersal rates
                      frac.f.disperse = index[[xx]]$frac.f.disperse,
                      frac.m.disperse = index[[xx]]$frac.m.disperse,
                      frac.r.disperse = index[[xx]]$frac.r.disperse,
                      
               
    )

  
    
    tt <- gsub(' ','_', gsub(':','',Sys.time()))
    #datenow <- Sys.Date()
    #fp <- sprintf('output/%s_expt_id%s_%s',datenow,expt_id,tt)
    fp <- sprintf('output/id%s_%s',expt_id,tt)
    system(paste0('mkdir ',fp,''))
    
    for(jj in seq(1,num_reps,1)){
      #print(jj)
      # datalist <- run.model(xx)
      # res <- datalist[1][[1]]
      # b <- datalist[2][[1]]
      


      ###now do collect_data
      collect_data(fp=fp,master_df=read_master,ind=xx,repn=jj,expt_id=expt_id,prms=prms,plot_yn=TRUE)
      
    }
    
  }
  read_master <- readRDS(sprintf('%s/output/master_df.RDS',getwd()))
  View(read_master)
}

#run.model(index=index,num_reps=5,scen=NA)









collect_data <- function(fp,master_df,ind,repn,expt_id,prms,plot_yn){
  
datenow <- Sys.Date()

experiment <- ind
repnum <- as.character(repn)
expt_id <-expt_id
#prms <- make.prms(ngens=4000,mu.r=0.05,sigma.m=0,b.fm=0,sigma.v=0.1,sigma.mf=-0.2,sigma.fm=-0.2,r.forage=0.05)

#collect data about prms
#Invasion gen
invasion_date <- prms$intrdate

 
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
res <- run.sim.inst(popfr, prms, fec=FALSE,set_seedlist=FALSE)
#use this seedlist for other comparisons
seedlist <- res[2][[1]]
res_fr <- res[1][[1]]



#save data
tt <- gsub(' ','_', gsub(':','',Sys.time()))
data.grab(res_fr,prms,fp,tt)
#load it back in and calculate average metapops
rundata <- list.files(fp,pattern=sprintf("%s.RData",tt))
rundf <- load(sprintf('%s/%s/%s',getwd(),fp,rundata))
rundf_data <- get(rundf)
res_fr <- rundf_data[1]
res_fr <- res_fr[[1]]

fr_metapop_data <- fitness_calculators(res_fr,list(FALSE,FALSE,FALSE))
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
res <- run.sim.inst(popfm, prms, fec=FALSE,set_seedlist=seedlist)
res_fm <- res[1][[1]]
tt <- gsub(' ','_', gsub(':','',Sys.time()))

#save outcome
data.grab(res_fm,prms,fp,tt)
rundata <- list.files(fp,pattern=sprintf("%s.RData",tt))
rundf <- load(sprintf('%s/%s/%s',getwd(),fp,rundata))
rundf_data <- get(rundf)
res_fm <- rundf_data[1]
res_fm <- res_fm[[1]]
fm_metapop_data <- fitness_calculators(res_fm,list(FALSE,FALSE,FALSE))
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
res <- run.sim.inst(pop, prms, fec=FALSE,set_seedlist=seedlist)
res <- res[1][[1]]
tt <- gsub(' ','_', gsub(':','',Sys.time()))
data.grab(res,prms,fp,tt)
rundata <- list.files(fp,pattern=sprintf("%s.RData",tt))
rundf <- load(sprintf('%s/%s/%s',getwd(),fp,rundata))
rundf_data <- get(rundf)
res_all3 <- rundf_data[1]
res_all3 <- res_all3[[1]]
#res_all3 <- res
all3_metapop_data <- fitness_calculators(res_all3,list(FALSE,FALSE,FALSE))
all3_f_avg_metapop <- all3_metapop_data['farmer mean metapop'][1,1]
all3_f_metapop_stdev <- all3_metapop_data['farmer metapop stdev'][1,1]
all3_m_avg_metapop <- all3_metapop_data['merc mean metapop'][1,1]
all3_m_metapop_stdev <- all3_metapop_data['merc metapop stdev'][1,1]
all3_r_avg_metapop <- all3_metapop_data['raider mean metapop'][1,1]
all3_r_metapop_stdev <- all3_metapop_data['raider metapop stdev'][1,1]


tt <- gsub(' ','_', gsub(':','',Sys.time()))

total_cost_matrix(res_all3, prms, seedlist, fp, tt, testing =FALSE)

component_cost_matrix(res_all3, make.prms(npatch=prms$npatch,ngens=prms$ngens),prms, seedlist, fp, tt, testing =FALSE,"mu.r")
component_cost_matrix(res_all3, make.prms(npatch=prms$npatch,ngens=prms$ngens),prms, seedlist, fp, tt, testing =FALSE,"a.fm","a.mf")
component_cost_matrix(res_all3, make.prms(npatch=prms$npatch,ngens=prms$ngens),prms, seedlist, fp, tt, testing =FALSE,"v.v","sigma.v")
component_cost_matrix(res_all3, make.prms(npatch=prms$npatch,ngens=prms$ngens),prms, seedlist, fp, tt, testing =FALSE,"b.fr","sigma.r")
component_cost_matrix(res_all3, make.prms(npatch=prms$npatch,ngens=prms$ngens),prms, seedlist, fp, tt, testing =FALSE,"r.forage")



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
Notes <- "NA"
#Datadir
datadir <- fp

#plotting
#make plot directory
plotfp <- sprintf('%s/%s_expt_id%s_%s_%s_%s',fp,datenow,expt_id,ind,repnum,tt)
system(paste0('mkdir ',getwd(),plotfp,''))


# length(c('Experiment','index','rep num','Pairwise F&R','Pairwise F&M','All 3',
# 'F&R F avg metapop','F&R R avg metapop','F&R F metapop stdev','F&R R metapop stdev',
# 'F&M F avg metapop','F&M M avg metapop','F&M F metapop stdev','F&M M metapop stdev',
# 'All 3 F avg metapop','All 3 M avg metapop','All 3 R avg metapop','All 3 F metapop stdev','All 3 M metapop stdev','All 3 R metapop stdev',
# 'R Foraging','Mercs true mutualists',
# 'Raw avg cost/ben','Std dev cost/ben',
# 'Notes','Datadir','Figdir',

# 'a.fm','a.mf','a.mr','a.rm','b.fr','conv.ftom','conv.ftor','d.f','d.m','d.r','dispersal.case',
# 'frac.f.disperse','frac.m.disperse','frac.r.disperse','g.f','g.m','g.r','h.fr','f_init','m_init','r_init', 'invasion gen',
# 'mu.fm','mu.mf','mu.r','mu.v','ngens','npatch','r.f','r.forage','raider.type','sigma.fm','sigma.mf','sigma.r','sigma.v',
# 'theta.fm','theta.mf','theta.r','theta.v','v.v'))


data_vec <- c(
    expt_id,
    ind,
    repnum,
    pairwise_fr_outcome,
    pairwise_fm_outcome,
    all3_outcome,

    fr_f_avg_metapop,
    fr_r_avg_metapop,
    fr_f_metapop_stdev,
    fr_r_metapop_stdev,

    fm_f_avg_metapop,
    fm_m_avg_metapop,
    fm_f_metapop_stdev,
    fm_m_metapop_stdev,

    all3_f_avg_metapop,
    all3_m_avg_metapop,
    all3_r_avg_metapop,
    all3_f_metapop_stdev,
    all3_m_metapop_stdev,
    all3_r_metapop_stdev,

    r_foraging,
    mercs_mutualists,

    raw_avg_cost_ben_mercs_to_farmers,
    stdev_cost_ben,
    Notes,
    datadir,
    plotfp,

    unlist(prms))

master_df <- readRDS(master_df)

master_transpose <- cbind(t(master_df),data_vec)
#View(t(master_transpose))
master_df <- t(master_transpose)
#works
saveRDS(master_df,sprintf('%s/output/master_df.RDS',getwd()))


}









run.sim.inst <- function(pop, prms, fec, set_seedlist){

    res <- array(NA, dim=c(prms$ngens, 3, prms$npatch))
  dimnames(res) <- list(1:prms$ngens,
                        c('num.m', 'num.f', 'num.r'),
                        1:prms$npatch)
  names(dimnames(res)) <- c('gen', 'output', 'patch')
  #update_list <- list(ppf = nf_next, ppm = nm_next, ppr= nr_next, interaction_monf = a.fm.mod, interaction_fonm = a.mf.mod, virulence_coef = vir.mod, attk_rate = attk.rate, consumed_prey = cons)

  

  if(fec == TRUE){
    fec.res <- array(NA, dim=c((prms$ngens-1),8,prms$npatch))
    dimnames(fec.res) <- list(1:(prms$ngens-1),
                              c(
                                'fec.f',
                               'fec.m',
                                'fec.r',
                                'interaction_monf',
                                'interaction_fonm',
                                'vir',
                                'attk_rate',
                                'consumed_prey'
                                ),
                              1:prms$npatch)
    names(dimnames(fec.res)) <- c('gen', 'output', 'patch')
  }
  


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
    for(i in 1:prms$ngens){ seed_block <- list(round(runif(15,1,prms$ngens)))
    seed_metalist <- append(seed_metalist, seed_block)
    }
  }
  
  #print(seed_metalist)
  
  for(ii in 1:prms$ngens) {
    
    gencounter <- gencounter +1
    set.seed(seed_metalist[[ii]][1])
    
    
    res[ii,'num.f',] <- pop[,'f']
    res[ii,'num.m',] <- pop[,'m']
    res[ii,'num.r',] <- pop[,'r']

    update_list <- update_pop(pp=pop, prms=prms, sid=seed_metalist[[ii]])

#update_pop() output is:
#update_list <- list(ppf = nf_next, ppm = nm_next, ppr= nr_next, interaction_monf = a.fm.mod, interaction_fonm = a.mf.mod, virulence_coef = vir.mod, attk_rate = attk.rate, consumed_prey = cons)
    if(fec == TRUE){
      if(gencounter < prms$ngens){
        fec.res[ii, 'fec.f',] <- pmax(update_list$ppf,0)
        fec.res[ii, 'fec.m',] <- pmax(update_list$ppm,0)
        fec.res[ii, 'fec.r',] <- pmax(update_list$ppr,0)
        fec.res[ii, 'interaction_monf',] <- update_list$interaction_monf
        fec.res[ii, 'interaction_fonm',] <- update_list$interaction_fonm
        fec.res[ii, 'vir',] <- update_list$virulence_coef
        fec.res[ii, 'attk_rate',] <- update_list$attk_rate
        fec.res[ii, 'consumed_prey',] <- update_list$consumed_prey
        
      }
    }
    #

    pop[,'f'] <- pmin(pmax(update_list$ppf,0),10000000)
    pop[,'m'] <- pmin(pmax(update_list$ppm,0),10000000)
    pop[,'r'] <- pmin(pmax(update_list$ppr,0),10000000)


    if(prms$npatch>1){

    if(prms$dispersal.case=='cloud')
      pop.after.disperse <- disperse.cloud.each(pp=pop, prms=prms, sid=seed_metalist[[ii]])
    
    if(prms$dispersal.case=='adjacent')
      pop.after.disperse <- disperse.adjacent(pp=pop)
    
    pop <- pop.after.disperse
    }
    
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

  







run.sim.seq <- function(pop, prms, fec, set_seedlist){

  #reproduction of farmers, 

}



quicksave_prms <- function(prms,filepath){
  if(dir.exists(sprintf("%s%s",getwd(),filepath))==FALSE){
    system(paste0('mkdir ',sprintf("%s%s",getwd(),filepath),''))
  }
  get_prms_names <- c(as.character(names(prms)))
  get_prms <-c(as.character(prms))
  prmsthing <- data.frame(get_prms,row.names = get_prms_names)
  tt <- gsub(' ','_', gsub(':','',Sys.time()))
  write.table(prmsthing,sprintf("%s%s/paramvec_%s.txt",getwd(),filepath,tt))
  #write.table(prmsthing,sprintf("%s/paramvec_%s.txt",getwd(),Sys.time()))

  
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
#   outcomes <- c(categorize.m=categorize.m,
#                 categorize.f=categorize.f,
#                 categorize.r=categorize.r)
  
  #run.data <- list(res=res, outcomes=outcomes, prms=prms)
  run.data <- list(res=res.data, prms=prms)
  
  save(run.data, file=sprintf('%s/ejr_run_%s%s_%s.RData', filepath, outcome_code, gens, timetag))
 
 #fix - doesn't seem to be writing
  write_xlsx(data.frame(prms[names(prms)!="init.pop.size"]),path=sprintf('%s/ejr_run_%s%s_%s.xlsx', filepath, outcome_code, gens, timetag))
}



#cutoff should either be false or an integer (row == generation index)
#cutoff should either be false or a vector of two numbers; lower and upper cutoffs to use to slice the generations
#this function uses runlength encoding to turn stretches of zeroes into single zeroes
calculate_fitness <- function(res_data,species,cutoff){
  
  if(cutoff == FALSE){
    run_array <- as.array(res_data)
  }else{
    run_array <- as.array(res_data)
    run_array <- run_array[cutoff[1]:cutoff[2]]
    #run_array <- run_array[cutoff:nrow(res_data), ,]
  }
  #run_array <- as.array(run_env$run.data$res)
  
  
  #run_array[,'num.m',]
  if(length(run_array[1,species,])<2){
    metapop <- run_array[,species,]
  }else{
    metapop <- apply(run_array[,species,], MARGIN = 1, FUN = sum)
  }
  
  
  rl <- rle(c(metapop)
              == 0)
  num_zero_runs <- sum(rl$values & rl$lengths >= 1)
  metapop <- replace(metapop,metapop==0,NA)
  metapop <- na.omit(metapop)

  if(length(metapop) > 0){
    #Calculate mean metapopulation size, counting stretches of zeroes as single zeroes
    ifelse(length(rep(0, num_zero_runs))==0,meanval<-mean(metapop),meanval<-mean(c(metapop,rep(0, num_zero_runs)),trim=0))
    
    #Calculate stdev of metapopulation size, counting stretches of zeroes as single zeroes
    ifelse(length(rep(0, num_zero_runs))==0,stdev<-sd(metapop),stdev<-sd(c(metapop,rep(0, num_zero_runs))))
    
    
  }else{
    meanval <- 0
    stdev <- 0
  }

  fitness_list <- list(meanval,stdev)
  return(fitness_list)
}

fitness_calculators <- function(res_data,cutoff_list){
  
  farmer_fitness <- calculate_fitness(res_data,'num.f',cutoff=cutoff_list[[1]])
  merc_fitness <- calculate_fitness(res_data,'num.m',cutoff=cutoff_list[[2]])
  raider_fitness <- calculate_fitness(res_data,'num.r',cutoff=cutoff_list[[3]])
  
  dfout <- data.frame(c(farmer_fitness,merc_fitness,raider_fitness))
  colnames(dfout) <- c('farmer mean metapop','farmer metapop stdev','merc mean metapop','merc metapop stdev','raider mean metapop','raider metapop stdev')
  return(dfout)
  
}













#09/10 modification to allow for npatch=1

#04/26/2023 index the populations in by generation in res, and index the seed for each operation
#using the seedlist
total_cost_matrix <- function(res.data, prms, seedlist, filepath, timetag, testing){
  
  
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
    test_res <- array(rep(0, each=prms$npatch*3*2), dim = c(prms$npatch, 3, 2))
    dimnames(test_res) <- list(1:prms$npatch,c("f","m","r"),c("n(t)","n(t+1)"))
  }
  
  
  temp_res <- array(rep(0, each=prms$npatch*3*2), dim = c(prms$npatch, 3, 2))
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

    temp_res_nt <- t(as.matrix(temp_res[,,'n(t)']))
    temp_res_nt_1 <- t(as.matrix(temp_res[,,'n(t+1)']))
    

    # 
    #print(test_res)
    # #match seed
    set.seed(seedlist[[i]][1])

    
    update_list <- update_pop(temp_res_nt,prms=prms,seedlist[[i]])

    temp_res_nt_1[,'f'] <- pmin(pmax(update_list$ppf,0),10000000)
    temp_res_nt_1[,'m'] <- pmin(pmax(update_list$ppm,0),10000000)
    temp_res_nt_1[,'r'] <- pmin(pmax(update_list$ppr,0),10000000)

    if(prms$npatch > 1){

        if(prms$dispersal.case=='cloud'){
          pop.after.disperse <- disperse.cloud.each(pp=temp_res_nt_1, prms=prms,sid=seedlist[[i]])}

        if(prms$dispersal.case=='adjacent'){
          pop.after.disperse <- disperse.adjacent(pp=temp_res_nt_1)}

        temp_res_nt_1 <- pop.after.disperse
    
    }

    
    temp_res[,,2] <- temp_res_nt_1


    growtharr_nerfmercs[i, 'farmer_mod_growth',] <- temp_res[,"f",2] - temp_res[,"f",1]
    growtharr_nerfmercs[i, 'merc_mod_growth',] <- temp_res[,"m",2] - temp_res[,"m",1]
    growtharr_nerfmercs[i, 'raider_mod_growth',] <- temp_res[,"r",2] - temp_res[,"r",1]
    


    if(testing == TRUE){
      print("temp_res")
      print(temp_res)
      print("test_res")
      print(test_res)
      print("difference in growth: mercs - no mercs")
      print(temp_res - test_res)
    }
    
    
    if(testing == TRUE){
      print(growtharr[i,,])
      print(growtharr_nerfmercs[i,,])
      print(growtharr[i,,] - growtharr_nerfmercs[i,,])
    }
  }
  


   
 #define & populate the merc total cost array
  
  tcostarr_m <- growtharr - growtharr_nerfmercs

  dimnames(tcostarr_m) <- list(1:(prms$ngens-1),
                               c('benefit_m_to_f','benefit_m_to_m','benefit_m_to_r'),
                               1:prms$npatch)
  names(dimnames(tcostarr_m)) <- c('gen', 'costs_of_m', 'patch')


#raiders

  
  #populate the modified growth array where raiders are nerfed
  temp_res <- array(rep(0, each=prms$npatch), dim = c(prms$npatch, 3, 2))
  dimnames(temp_res) <- list(1:prms$npatch,c("f","m","r"),c('n(t)','n(t+1)'))
  for(i in seq(1,(prms$ngens-1),1)){


    #test
    #for(i in seq(1,(100),1)){
    if(testing==TRUE){
      print(i)
    }
    
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

    # temp_res_nt <- t(as.matrix(temp_res[,,'n(t)']))
    # temp_res_nt_1 <- t(as.matrix(temp_res[,,'n(t+1)']))
    

    #test
    if(testing == TRUE){
      test_res[,"f",1] <- res.data[i,'num.f',]
      test_res[,"m",1] <- res.data[i,'num.m',]
      test_res[,"r",1] <- res.data[i,'num.r',]
      
      test_res[,"f",2] <- res.data[(i+1),'num.f',]
      test_res[,"m",2] <- res.data[(i+1),'num.m',]
      test_res[,"r",2] <- res.data[(i+1),'num.r',]
    }
    
    temp_res_nt <- t(as.matrix(temp_res[,,'n(t)']))
    temp_res_nt_1 <- t(as.matrix(temp_res[,,'n(t+1)']))

    #match seed
    #set.seed(as.integer(seedlist[i]))
    # 
    #print(test_res)
    # #match seed
    set.seed(seedlist[[i]][1])

    update_list <- update_pop(temp_res_nt,prms=prms,seedlist[[i]])


    temp_res_nt_1[,'f'] <- pmin(pmax(update_list$ppf,0),10000000)
    temp_res_nt_1[,'m'] <- pmin(pmax(update_list$ppm,0),10000000)
    temp_res_nt_1[,'r'] <- pmin(pmax(update_list$ppr,0),10000000)

    if(prms$npatch>1){

    if(prms$dispersal.case=='cloud')
      pop.after.disperse <- disperse.cloud.each(pp=temp_res_nt_1, prms=mod_prms, sid=seed_metalist[[ii]])
    
    #probably needs fixing
    if(prms$dispersal.case=='adjacent')
      pop.after.disperse <- disperse.adjacent(pp=temp_res_nt_1)
    
    temp_res_nt_1 <- pop.after.disperse
    }



    
    temp_res[,,2] <- temp_res_nt_1


    growtharr_nerfraiders[i, 'farmer_mod_growth',] <- temp_res[,"f",2] - temp_res[,"f",1]
    growtharr_nerfraiders[i, 'merc_mod_growth',] <- temp_res[,"m",2] - temp_res[,"m",1]
    growtharr_nerfraiders[i, 'raider_mod_growth',] <- temp_res[,"r",2] - temp_res[,"r",1]
    


    if(testing == TRUE){
      print("temp_res")
      print(temp_res)
      print("test_res")
      print(test_res)
      print("difference in growth: raiders - no raiders")
      print(temp_res - test_res)
    }
    
    
    if(testing == TRUE){
      print(growtharr[i,,])
      print(growtharr_nerfraiders[i,,])
      print(growtharr[i,,] - growtharr_nerfraiders[i,,])
    }
  }



   #define & populate the raider total cost array
  tcostarr_r <- growtharr - growtharr_nerfraiders
  #tcost_m_r <- array(NA, dim=c(prms$ngens-1, 3, prms$npatch))
  dimnames(tcostarr_r) <- list(1:(prms$ngens-1),
                               c('benefit_r_to_f','benefit_r_to_m','benefit_r_to_r'),
                               1:prms$npatch)
  names(dimnames(tcostarr_r)) <- c('gen', 'costs_of_r', 'patch')



  #################

  
  #save data
  total_cost.data <- list(res=res.data, growth_array = growtharr, m_totalcost_array = tcostarr_m,
                          r_totalcost_array = tcostarr_r, prms=prms)
  save(total_cost.data, file=sprintf('%s/kdm_total_costs_%s.RData',filepath, timetag))
}


#component_cost_matrix() calculates the cost of a given parameter or subset of parameters by comparing population
#growth between cases where those parameters are some non-zero value (specified) against that where those prms 
#are set to zero instead

#example call
#component_cost_matrix(res_all3, prms,prms, seedlist, fp, tt, testing =FALSE,"mu.r")

#def_prms are the default values for parameters
#curr_prms are the current parameter values being used
#an arbitrary number of prms can be evaluated for their effect on cost (by comparing to the case where they're set to zero)
#if the default and current prm sets are different, it will be noted that they're also different, and the modified values will be used
#in the comparison, so as to not change more than one prm value at once
component_cost_matrix <- function(res.data, def_prms, curr_prms, seedlist, filepath, timetag, testing, ...){
  
  #get parameters for measuring effect in a list
  get_comps <-  list(...)
  
#   if("r.f" %in% get_comps){
#     print("if r.f is zero, so will be v.m")
#     get_comps <- append(get_comps,'v.v')
#   }


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
  curr_diff_prms_names <- names(curr_prms[as.character(curr_prms)!=as.character(def_prms)])
  #print(c(curr_diff_prms_names))
  
  #names and values
  #curr_diff_prms <- as.list(c(names(Sys.getenv(curr_prms[as.character(curr_prms)!=as.character(prms)]))))
  curr_diff_prms <- as.list(c(curr_prms[as.character(curr_prms)!=as.character(def_prms)]))
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

  if(length(curr_diff_prms_names)>0){

    curr_diff_prms <- curr_diff_prms[!(names(curr_diff_prms) %in% unlist(get_comps))]
  #print(curr_diff_prms)
  

curr_diff_prms[names(curr_diff_prms)%in%curr_diff_prms_names]

as.character(curr_diff_prms[names(curr_diff_prms)%in%curr_diff_prms_names])

combine_prms <- lapply(seq_along(curr_diff_prms), function(i) {
  sprintf("%s=%s", curr_diff_prms_names[i], as.character(curr_diff_prms[i]))
})

# Combining the results into a single string
curr_input_for_make.prms <- paste(combine_prms, collapse=",")
  }else{
    curr_input_for_make.prms <- list()
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
  
  if(mod_prms$ngens!=curr_prms$ngens){
    stop("gen number must be equal")
  }
  if(mod_prms$npatch !=curr_prms$npatch){
    stop("patch number must be equal")
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
    test_res <- array(rep(0, each=curr_prms$npatch), dim = c(curr_prms$npatch, 3, 2))
    dimnames(test_res) <- list(1:curr_prms$npatch,c("f","m","r"),c("n(t)","n(t+1)"))
  }
  
  #populate modified growth array (B)
  temp_res <- array(rep(0, each=curr_prms$npatch), dim = c(curr_prms$npatch, 3, 2))
  #dimnames(temp_res) <- list(1:mod_prms$npatch,c("f","m","r"),1:2)
  dimnames(temp_res) <- list(1:curr_prms$npatch,c("f","m","r"),c('n(t)','n(t+1)'))

  for(i in seq(1,(curr_prms$ngens-1),1)){
    
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
    
    temp_res_nt <- t(as.matrix(temp_res[,,'n(t)']))
    temp_res_nt_1 <- t(as.matrix(temp_res[,,'n(t+1)']))

    
    #match seed
    #set.seed(as.integer(seedlist[i]))
    # 
    #print(test_res)
    # #match seed
    set.seed(seedlist[[i]][1])
    # 
    # ## farmer reproduction
    #HERE1
    #reprod.f.temp <- repp(pp=temp_res[,,2,drop=FALSE], prms=prms,sid=seedlist[[i]])
    #reprod.f.temp <- repp(pp=temp_res_nt_1, prms=def_prms,sid=seedlist[[i]])
#    reprod.f.temp <- repp(pp=temp_res_nt_1, prms=def_prms,sid=seedlist[[i]])

    update_list <- update_pop(temp_res_nt_1,prms=mod_prms,seedlist[[i]])

    temp_res_nt_1[,'f'] <- pmin(pmax(update_list$ppf,0),10000000)
    temp_res_nt_1[,'m'] <- pmin(pmax(update_list$ppm,0),10000000)
    temp_res_nt_1[,'r'] <- pmin(pmax(update_list$ppr,0),10000000)

    if(mod_prms$npatch>1){

    if(mod_prms$dispersal.case=='cloud')
      pop.after.disperse <- disperse.cloud.each(pp=temp_res_nt_1, prms=mod_prms, sid=seedlist[[i]])
    
    #probably needs fixing
    if(mod_prms$dispersal.case=='adjacent')
      pop.after.disperse <- disperse.adjacent(pp=temp_res_nt_1)
    
    temp_res_nt_1 <- pop.after.disperse
    }
    
  
  
    temp_res[,,2] <- temp_res_nt_1
    if(testing == TRUE){
      print("temp_res")
      print(temp_res)
      print("test_res")
      print(test_res)
      sprintf("difference in growth: %s[!=0] - %s[=0]",get_comps,get_comps)
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













