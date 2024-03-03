
makeTransparent <- function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  apply(newColor, 2, .makeTransparent, alpha=alpha)
}



plot_outcome_by_prm_values <- function(data,plot_fp,filter_yn==TRUE){

#data must be master_df
    # names(master_df)
    # head(master_df)

#get rid of a bunch of columns
    results_df <- cbind(data[,6],data[,28:ncol(data)])
    results_df <- results_df %>%
        select(-r_init,-f_init,-m_init, -dispersal.case, -raider.type, -print.pop, -npatch, -ngens)
    results_df <- cbind(results_df[,1:8],results_df[,10:ncol(results_df)])

    #head(results_df)

    if(filter_yn == TRUE){

        #identify cases where a positive feedback loop has caused population explosion in the farmers and annotate outcome by
        #whether or not these kinds of increasing returns likely happened
        results_df$filtered_outcome <- paste0(results_df[,1],"_",ifelse(as.numeric(master_df[,11]) > as.numeric(100000),"inc","cont"))

        # nrow(results_df[results_df$filtered_outcome == "coex_cont",])
        # master_df[,"F&M F avg metapop"][results_df$filtered_outcome == "coex_cont"]
        # nrow(master_df[as.numeric(master_df[,11]) < as.numeric(100000) & master_df[,"All 3"]=="coex",])
        # head(master_df[,11][as.numeric(master_df[,11]) < as.numeric(100000) & master_df[,"All 3"]=="coex"], n=1000)

        for(i in 2:ncol(results_df)){
            plot <- ggplot(results_df, aes(x = as.numeric(results_df[,i]), fill = filtered_outcome)) +
              #geom_density(alpha = 0.5) +
              geom_histogram(alpha = 0.5) +
              theme_minimal() +
              labs(title = NULL,
                   x = "prm value",
                   y = "Density")
            ggsave(sprintf("%s/plot_%s.png",plot_fp,names(results_df)[i]),plot)
        }
    }else{
        for(i in 2:ncol(results_df)){
            plot <- ggplot(results_df, aes(x = as.numeric(results_df[,i]), fill = results_df[,1])) +
              #geom_density(alpha = 0.5) +
              geom_histogram(alpha = 0.5) +
              theme_minimal() +
              labs(title = NULL,
                   x = "prm value",
                   y = "Density")
            ggsave(sprintf("%s/plot_%s.png",plot_fp,names(results_df)[i]),plot)
        }
    }


}






#plot.pop.all(test_data$rundata[[1]]$prms, test_data$rundata[[1]]$res,"/Users/keiranmaskell/Desktop/quicktest",print=TRUE)

plot.pop.all <- function(prms, res, plotfp, sumpatches=TRUE, print=FALSE){
  
  ymin <- 0
  if(dim(res)[3]==1){
    ymax <- max(res)
    #ymin <- min(res)
    numf_data <- res[,'num.f',]
    numm_data <- res[,'num.m',]
    numr_data <- res[,'num.r',]

  }else{
    if(sumpatches==TRUE){
        ymax <- max(rowSums(res[, 'num.f', ]),rowSums(res[, 'num.m', ]),rowSums(res[, 'num.r', ]))
        numf_data <- rowSums(res[,'num.f',])
        numm_data <- rowSums(res[,'num.m',])
        numr_data <- rowSums(res[,'num.r',])


      if(print==TRUE){
        quartz()
    
        plot(NA,
           xlim=c(1,prms$ngens),
           #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
           ylim=c(ymin,ymax),
           main = sprintf('Total Populations'),
           xlab='Generation',
           ylab='Number',
           las=1)
    
    
      lines(x=1:prms$ngens,
            y=numf_data,
            col='black')
      lines(x=1:prms$ngens,
            y=numm_data,
            col='dodgerblue')
      lines(x=1:prms$ngens,
            y=numr_data,
            col='red')
      legend('topright',
             legend=c('Farmers',
                      'Mercenaries',
                      'Raiders'),
             lty=rep(1,3),
             col=c('black', 'dodgerblue', 'red'),
             bty='n')
    
      }else{
        pdf(file=sprintf('%s/pop_allpatches_plot.pdf',plotfp), width=5, height=5)
    
            plot(NA,
               xlim=c(1,prms$ngens),
               #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
               ylim=c(ymin,ymax),
               main = sprintf('Total Populations'),
               xlab='Generation',
               ylab='Number',
               las=1)
            lines(x=1:prms$ngens,
                  y=numf_data,
                  col='black')
            lines(x=1:prms$ngens,
                  y=numm_data,
                  col='dodgerblue')
            lines(x=1:prms$ngens,
                  y=numr_data,
                  col='red')
            legend('topright',
                   legend=c('Farmers',
                            'Mercenaries',
                            'Raiders'),
                   lty=rep(1,3),
                   col=c('black', 'dodgerblue', 'red'),
                   bty='n')
          dev.off()
    
      }

    }else{

        if(print==TRUE){
            quartz()

            plot(NA,
               xlim=c(1,prms$ngens),
               #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
               ylim=c(ymin,ymax),
               main = sprintf('Total Populations'),
               xlab='Generation',
               ylab='Number',
               las=1)

            for(i in 1:dim(res)[3]){
            
                numf_data <- res[,'num.f',i]
                numm_data <- res[,'num.m',i]
                numr_data <- res[,'num.r',i]

                lines(x=1:prms$ngens,
                      y=numf_data,
                      col='black')
                lines(x=1:prms$ngens,
                      y=numm_data,
                      col='dodgerblue')
                lines(x=1:prms$ngens,
                      y=numr_data,
                      col='red')

            }

          legend('topright',
                 legend=c('Farmers',
                          'Mercenaries',
                          'Raiders'),
                 lty=rep(1,3),
                 col=c('black', 'dodgerblue', 'red'),
                 bty='n')
  
        }else{
          pdf(file=sprintf('%s/pop_allpatches_plot.pdf',plotfp), width=5, height=5)

             plot(NA,
                 xlim=c(1,prms$ngens),
                 #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
                 ylim=c(ymin,ymax),
                 main = sprintf('Total Populations'),
                 xlab='Generation',
                 ylab='Number',
                 las=1)
              for(i in 1:dim(res)[3]){
                
             numf_data <- res[,'num.f',i]
              numm_data <- res[,'num.m',i]
              numr_data <- res[,'num.r',i]

             lines(x=1:prms$ngens,
                    y=numf_data,
                    col='black')
              lines(x=1:prms$ngens,
                    y=numm_data,
                    col='dodgerblue')
              lines(x=1:prms$ngens,
                    y=numr_data,
                    col='red')

         }

             legend('topright',
                     legend=c('Farmers',
                              'Mercenaries',
                              'Raiders'),
                     lty=rep(1,3),
                     col=c('black', 'dodgerblue', 'red'),
                     bty='n')
            dev.off()

       }
    }
  }
}
  

plot.pop.patch(1,test_data$rundata[[1]]$prms, test_data$rundata[[1]]$res,"/Users/keiranmaskell/Desktop/quicktest",print=FALSE)


plot.pop.patch <- function(patch, prms, res, plotfp, print=FALSE){

    # if(dir.exists(sprintf('%s',plotfp))==FALSE){
    #   system(paste0('mkdir ',plotfp))
    # }

    ymin <- 0
    ymax <- max(res)
      #ymin <- min(res)
    numf_data <- res[,'num.f',patch]
    numm_data <- res[,'num.m',patch]
    numr_data <- res[,'num.r',patch]


    if(print==TRUE){
        quartz()
        plot(NA,
            xlim=c(1,prms$ngens),
            #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
            ylim=c(ymin,ymax),
            main = sprintf('Total Populations in %s', patch),
            xlab='Generation',
            ylab='Number',
            las=1)
        lines(x=1:prms$ngens,
              y=numf_data,
              col='black')
        lines(x=1:prms$ngens,
              y=numm_data,
              col='dodgerblue')
        lines(x=1:prms$ngens,
              y=numr_data,
              col='red')
        legend('topright',
               legend=c('Farmers',
                        'Mercenaries',
                        'Raiders'),
               lty=rep(1,3),
               col=c('black', 'dodgerblue', 'red'),
               bty='n')


    }else{
        pdf(file=sprintf('%s/pop_patch_%s_plot.pdf',plotfp,patch), width=5, height=5)

        plot(NA,
            xlim=c(1,prms$ngens),
            #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
            ylim=c(ymin,ymax),
            main = sprintf('Total Populations in %s', patch),
            xlab='Generation',
            ylab='Number',
            las=1)
        lines(x=1:prms$ngens,
              y=numf_data,
              col='black')
        lines(x=1:prms$ngens,
              y=numm_data,
              col='dodgerblue')
        lines(x=1:prms$ngens,
              y=numr_data,
              col='red')
        legend('topright',
               legend=c('Farmers',
                        'Mercenaries',
                        'Raiders'),
               lty=rep(1,3),
               col=c('black', 'dodgerblue', 'red'),
               bty='n')

        dev.off()



    }
}
  
  
  
  
#plot.meta 

  
  ## p1_path <- file.path("/Users/keiranmaskell/Desktop/Output",paste("p1_pops", Sys.Date(), expnumber, ".pdf", sep = ""))
  ## pdf(file=p1_path)

#tcostarr <- tcost_env$total_cost.data$m_totalcost_array
#species <- "Mercenaries"
#patch <- 1



#you are here 02/25 #################################



plot.pop.tcosts <- function(patch, prms, tcostarr, species, sep, plotfp){
  
  if(dir.exists(sprintf('%s%s/total_costs',getwd(),plotfp))==FALSE){
    system(paste0('mkdir ',getwd(),plotfp,'/total_costs'))
    }
  
  
  if(sep == 'y'){
    if(species == 'Mercenaries'){
      pdf(file=sprintf('%s%s/total_costs/cost_%s_on_farmers_patch_%s_plot.pdf',getwd(),plotfp,species,patch), width=5, height=5)
      plot(NA,
           xlim=c(1,(prms$ngens-1)),
           ylim=c(min(tcostarr[, 1, patch]),max(tcostarr[, 1, patch])),
           main = sprintf("Total Cost of %s on Farmers in Patch %s", species, patch),
           xlab='Generation',
           ylab='Net Benefit (+) or Net Cost (-)',
           las=0)
      lines(x=1:(prms$ngens-1),
            y=tcostarr[, 1, patch],
            col='black')
      lines(x=1:(prms$ngens-1),
            y=filter(tcostarr[,1,patch],rep(1/6,6),sides=2),
            col='blue')
      lines(x=1:(prms$ngens-1),
            y=rep(mean(tcostarr[,1,patch]),prms$ngens-1),
            col='red')
      legend('topright',
             legend=c('Merc ben','Merc mov_avg ben','Merc avg ben'),
              lty=rep(1,2),
              lwd=rep(2,2),
              col=c('black','blue','red'),
              bty='n')
      dev.off()
      
      pdf(file=sprintf('%s%s/total_costs/cost_%s_on_raiders_patch_%s_plot.pdf',getwd(),plotfp,species,patch), width=5, height=5)
      plot(NA,
           xlim=c(1,(prms$ngens-1)),
           ylim=c(min(tcostarr[, 3, patch]),max(tcostarr[, 3, patch])),
           main = sprintf("Total Cost of %s on Raiders in Patch %s", species, patch),
           xlab='Generation',
           ylab='Net Benefit (+) or Net Cost (-)',
           las=0)
      lines(x=1:(prms$ngens-1),
            y=tcostarr[, 3, patch],
            col='red')
      lines(x=1:(prms$ngens-1),
            y=filter(tcostarr[,3,patch],rep(1/6,6),sides=2),
            col='orange')
      lines(x=1:(prms$ngens-1),
            y=rep(mean(tcostarr[,3,patch]),prms$ngens-1),
            col='blue')
      legend('topright',
             legend=c('Merc ben','Merc mov_avg ben','Merc avg ben'),
              lty=rep(1,2),
              lwd=rep(2,2),
              col=c('red','orange','blue'),
              bty='n')
      dev.off()
    }
    if(species == 'Raiders'){
      pdf(file=sprintf('%s%s/total_costs/cost_%s_on_farmers_patch_%s_plot.pdf',getwd(),plotfp,species,patch), width=5, height=5)
      plot(NA,
           xlim=c(1,(prms$ngens-1)),
           ylim=c(min(tcostarr[,1,patch]),max(tcostarr[,1,patch])),
           main = sprintf("Total Cost of %s on Farmers in Patch %s", species, patch),
           xlab='Generation',
           ylab='Net Benefit (+) or Net Cost (-)',
           las=0)
      lines(x=1:(prms$ngens-1),
            y=tcostarr[, 1, patch],
            col='black')
      lines(x=1:(prms$ngens-1),
            y=filter(tcostarr[,1,patch],rep(1/6,6),sides=2),
            col='blue')
      lines(x=1:(prms$ngens-1),
            y=rep(mean(tcostarr[,1,patch]),prms$ngens-1),
            col='red')
      legend('topright',
             legend=c('Raiders ben','Raiders mov avg ben','Raiders avg ben'),
                      lty=rep(1,2),
                      lwd=rep(2,2),
                      col=c('black','blue','red'),
                      bty='n')
      dev.off()
      
      pdf(file=sprintf('%s%s/total_costs/cost_%s_on_mercs_patch_%s_plot.pdf',getwd(),plotfp,species,patch), width=5, height=5)
      plot(NA,
           xlim=c(1,(prms$ngens-1)),
           ylim=c(min(tcostarr[, 2, patch]),max(tcostarr[, 2, patch])),
           main = sprintf("Total Cost of %s on Mercenaries in Patch %s", species, patch),
           xlab='Generation',
           ylab='Net Benefit (+) or Net Cost (-)',
           las=0)
      lines(x=1:(prms$ngens-1),
            y=tcostarr[, 2, patch],
            col='dodgerblue')
      lines(x=1:(prms$ngens-1),
            y=filter(tcostarr[,2,patch],rep(1/6,6),sides=2),
            col='hotpink')
      lines(x=1:(prms$ngens-1),
            y=rep(mean(tcostarr[,2,patch]),prms$ngens-1),
            col='green')
      legend('topright',
             legend=c('Raiders ben',"Raiders mov avg ben", "Raiders avg ben"),
              lty=rep(1,2),
              lwd=rep(2,2),
              col=c('dodgerblue','hotpink','green'),
              bty='n')
      dev.off()
    }
  }
  if(sep == 'n'){
      
  if(species == 'Mercenaries'){
    pdf(file=sprintf('%s%s/total_costs/cost_%s_patch_%s_plot.pdf',getwd(),plotfp,species,patch), width=5, height=5)
    plot(NA,
         xlim=c(1,(prms$ngens-1)),
         ylim=c(min(tcostarr[, 1, patch],tcostarr[, 3, patch]),max(tcostarr[, 1, patch],tcostarr[, 3, patch])),
         main = sprintf("Total Cost of %s on other spp. in patch %s", species, patch),
         xlab='Generation',
         ylab='Net Benefit (+) or Net Cost (-)',
         las=0)
      lines(x=1:(prms$ngens-1),
          y=tcostarr[, 1, patch],
          col='black')
      lines(x=1:(prms$ngens-1),
          y=tcostarr[, 3, patch],
          col='red')
      lines(x=1:(prms$ngens-1),
            y=filter(tcostarr[,1,patch],rep(1/6,6),sides=2),
            col='blue')
      lines(x=1:(prms$ngens-1),
            y=filter(tcostarr[,3,patch],rep(1/6,6),sides=2),
            col='orange')
      legend('topright',
           legend=c('Farmers',
                    'Raiders','Farmers avg','Raiders avg'),
           lty=rep(1,4),
           lwd=rep(2,4),
           col=c('black', 'red','blue','orange'),
           bty='n')
  dev.off()
  }
  if(species == "Raiders"){
    pdf(file=sprintf('%s%s/total_costs/cost_%s_patch_%s_plot.pdf',getwd(),plotfp,species,patch), width=5, height=5)
    plot(NA,
         xlim=c(1,(prms$ngens-1)),
         ylim=c(min(tcostarr[, 1, patch],tcostarr[, 2, patch]),max(tcostarr[, 1, patch],tcostarr[, 2, patch])),
         main = sprintf("Total Cost of %s on other spp. in patch %s", species, patch),
         xlab='Generation',
         ylab='Net Benefit (+) or Net Cost (-)',
         las=0)
      lines(x=1:(prms$ngens-1),
          y=tcostarr[, 1, patch],
          col='black')
      lines(x=1:(prms$ngens-1),
          y=tcostarr[, 2, patch],
          col='dodgerblue')
      lines(x=1:(prms$ngens-1),
            y=filter(tcostarr[,1,patch],rep(1/6,6),sides=2),
            col='blue')
      lines(x=1:(prms$ngens-1),
            y=filter(tcostarr[,2,patch],rep(1/6,6),sides=2),
            col='hotpink')
      legend('topright',
           legend=c('Farmers',
                    'Mercenaries','Farmers avg','Mercs avg'),
           lty=rep(1,4),
           lwd=rep(2,4),
           col=c('black', 'dodgerblue','blue','hotpink'),
           bty='n')
    dev.off()
    }
  }
}





# foo <- function(fmt, ...) {
#   m <- length(gregexpr("%s", fmt, fixed = TRUE)[[1]])
#   args <- list(...)
#   n <- ceiling(m/length(args))
#   args <- rep(args, n)
#   args <- args[seq_len(m)]
#   do.call(sprintf, c(args, fmt = fmt))
# }


plot.pop.comp.costs <- function(patch, prms, ccostarr, orig.prms, mod.prms, sep, movavg, plotfp){
  
  if(dir.exists(sprintf('%s%s/component_costs',getwd(),plotfp))==FALSE){
    system(paste0('mkdir ',getwd(),plotfp,'/component_costs'))
  }
  
  modified_parameters <- orig.prms[sapply(names(mod.prms), function(x) !identical(orig.prms[[x]], mod.prms[[x]]))]
  modified_parameter_names <- names(modified_parameters)
  mod_par_str <- paste(modified_parameters,collapse=',')
  mod_par_names_str <- paste(modified_parameter_names,collapse=',')
  
  if(sep == 'y'){
    pdf(file=sprintf('%s%s/component_costs/component_cost_%s%s_on_farmers_patch_%s_plot.pdf',getwd(),plotfp,mod_par_names_str,mod_par_str,patch), width=5, height=5)
    plot(NA,
         xlim=c(1,(prms$ngens-1)),
         ylim=c(min(ccostarr[, 1, patch]),max(ccostarr[, 1, patch])),
         main = sprintf("Effect of %s = %s on Farmers in Patch %s", mod_par_names_str, mod_par_str, patch),
         xlab='Generation',
         ylab='Net Benefit (+) or Net Cost (-)',
         las=0)
    lines(x=1:(prms$ngens-1),
          y=ccostarr[, 1, patch],
          col='black')
    lines(x=1:(prms$ngens-1),
          y=filter(ccostarr[,1,patch],rep(1/movavg,movavg),sides=2),
          col='blue')
    legend('topright',
           legend=c('Farmers_raw','Farmers_av'),
           lty=rep(1,2),
           lwd=rep(1,2),
           col=c('black','blue'),
           bty='n')
    dev.off()
    
    pdf(file=sprintf('%s%s/component_costs/component_cost_%s%s_on_mercs_patch_%s_plot.pdf',getwd(),plotfp,mod_par_names_str,mod_par_str,patch), width=5, height=5)
    plot(NA,
         xlim=c(1,(prms$ngens-1)),
         ylim=c(min(ccostarr[, 2, patch]),max(ccostarr[, 2, patch])),
         main = sprintf("Effect of %s = %s on Mercenaries in Patch %s", mod_par_names_str, mod_par_str, patch),
         xlab='Generation',
         ylab='Net Benefit (+) or Net Cost (-)',
         las=0)
    lines(x=1:(prms$ngens-1),
          y=ccostarr[, 2, patch],
          col='dodgerblue')
    lines(x=1:(prms$ngens-1),
          y=filter(ccostarr[,2,patch],rep(1/movavg,movavg),sides=2),
          col='hotpink')
    legend('topright',
           legend=c('Mercenaries raw','Mercenaries avg'),
           lty=rep(1,2),
           lwd=rep(2,2),
           col=c('dodgerblue','hotpink'),
           bty='n')
    dev.off()
    
    pdf(file=sprintf('%s%s/component_costs/component_cost_%s%s_on_raiders_patch_%s_plot.pdf',getwd(),plotfp,mod_par_names_str,mod_par_str,patch), width=5, height=5)
    plot(NA,
         xlim=c(1,(prms$ngens-1)),
         ylim=c(min(ccostarr[, 3, patch]),max(ccostarr[, 3, patch])),
         main = sprintf("Effect of %s = %s on Raiders in Patch %s", mod_par_names_str, mod_par_str, patch),
         xlab='Generation',
         ylab='Net Benefit (+) or Net Cost (-)',
         las=0)
    lines(x=1:(prms$ngens-1),
          y=ccostarr[, 3, patch],
          col='red')
    lines(x=1:(prms$ngens-1),
          y=filter(ccostarr[,3,patch],rep(1/movavg,movavg),sides=2),
          col='orange')
    legend('topright',
           legend=c('Raiders raw', 'Raiders avg'),
           lty=rep(1,2),
           lwd=rep(2,2),
           col=c('red','orange'),
           bty='n')
    dev.off()
  }
    
  if(sep == 'n'){
    pdf(file=sprintf('%s%s/component_costs/component_cost_%s%s_patch_%s_plot.pdf',getwd(),plotfp,mod_par_names_str,mod_par_str,patch), width=5, height=5)
    plot(NA,
         xlim=c(1,(prms$ngens-1)),
         ylim=c(min(ccostarr[, 1, patch],ccostarr[, 2, patch],ccostarr[, 3, patch]),max(ccostarr[, 1, patch],ccostarr[, 2, patch],ccostarr[, 3, patch])),
         main = sprintf("Effect of %s = %s on other spp. in patch %s", mod_par_names_str, mod_par_str, patch),
         xlab='Generation',
         ylab='Net Benefit (+) or Net Cost (-)',
         las=0)
    lines(x=1:(prms$ngens-1),
          y=ccostarr[, 1, patch],
          col='black')
    lines(x=1:(prms$ngens-1),
          y=ccostarr[, 2, patch],
          col='dodgerblue')
    lines(x=1:(prms$ngens-1),
          y=ccostarr[, 3, patch],
          col='red')
    lines(x=1:(prms$ngens-1),
          y=filter(ccostarr[,1,patch],rep(1/movavg,movavg),sides=2),
          col='blue')
    lines(x=1:(prms$ngens-1),
          y=filter(ccostarr[,2,patch],rep(1/movavg,movavg),sides=2),
          col='hotpink')
    lines(x=1:(prms$ngens-1),
          y=filter(ccostarr[,3,patch],rep(1/movavg,movavg),sides=2),
          col='orange')
    legend('topright',
           legend=c('Farmers','Mercs',
                    'Raiders','Farmers avg','Mercs avg','Raiders avg'),
           lty=rep(1,6),
           lwd=rep(2,6),
           col=c('black','dodgerblue','red','blue','hotpink','orange'),
           bty='n')
    dev.off()
    
    
    }
  }

plot.pop.growth <- function(patch, prms, growtharr, sep, plotfp){
  
  if(dir.exists(sprintf('%s%s/pop_growth',getwd(),plotfp))==FALSE){
    system(paste0('mkdir ',getwd(),plotfp,'/pop_growth'))
  }
  
  
  if(sep=='y'){
    pdf(file=sprintf('%s%s/pop_growth/growth_farmers_patch_%s_plot.pdf',getwd(),plotfp,patch), width=5, height=5)
    plot(NA,
         xlim=c(1,(prms$ngens-1)),
         #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
         #ylim=c(-50,50),
         #ylim=c(min(growtharr[,1,patch]),max(growtharr[,1,patch])),
         ylim=c(min(growtharr[,1,patch]),max(growtharr[,1,patch])),
         main = sprintf('Step-wise growth rate of farmers in %s', patch),
         xlab='Generation',
         ylab='Step-wise growth rate',
         las=1)
    lines(x=1:(prms$ngens-1),
          y=growtharr[, 1, patch],
          col='black')
    legend('topright',
           legend='Farmer Pop. Growth',
           lty=1,
           lwd=2,
           col='black',
           bty='n')
    dev.off()
    
    pdf(file=sprintf('%s%s/pop_growth/growth_mercs_patch_%s_plot.pdf',getwd(),plotfp,patch), width=5, height=5)
    plot(NA,
         xlim=c(1,(prms$ngens-1)),
         #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
         ylim=c(min(growtharr[, 2, patch]),max(growtharr[, 2, patch])),
         #ylim=c(min(growtharr[,2,patch]),max(growtharr[,2,patch])),
         main = sprintf('Step-wise growth rate of mercenaries in %s', patch),
         xlab='Generation',
         ylab='Step-wise growth rate',
         las=1)
    lines(x=1:(prms$ngens-1),
          y=growtharr[, 2, patch],
          col='dodgerblue')
    legend('topright',
           legend='Mercenary Pop. Growth',
           lty=1,
           lwd=2,
           col='dodgerblue',
           bty='n')
    dev.off()
    
    pdf(file=sprintf('%s%s/pop_growth/growth_raiders_patch_%s_plot.pdf',getwd(),plotfp,patch), width=5, height=5)
    plot(NA,
         xlim=c(1,(prms$ngens-1)),
         #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
         #ylim=c(-50,50),
         ylim=c(min(growtharr[,3,patch]),max(growtharr[,3,patch])),
         main = sprintf('Step-wise growth rate of raiders in %s', patch),
         xlab='Generation',
         ylab='Step-wise growth rate',
         las=1)
    lines(x=1:(prms$ngens-1),
          y=growtharr[, 3, patch],
          col='red')
    legend('topright',
           legend='Raider Pop. Growth',
           lty=1,
           lwd=2,
           col='red',
           bty='n')
    dev.off()
    
  }
  if(sep=='n')
  
  pdf(file=sprintf('%s%s/pop_growth/growth_patch_%s_plot.pdf',getwd(),plotfp,patch), width=5, height=5)
  plot(NA,
       xlim=c(1,(prms$ngens-1)),
       #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
       ylim=c(min(growtharr[, , patch]),max(growtharr[, , patch])),
       main = sprintf('Step-wise growth rates in %s', patch),
       xlab='Generation',
       ylab='Step-wise growth rate',
       las=1)
  lines(x=1:(prms$ngens-1),
        y=growtharr[, 1, patch],
        col='black')
  lines(x=1:(prms$ngens-1),
        y=growtharr[, 2, patch],
        col='dodgerblue')
  lines(x=1:(prms$ngens-1),
        y=growtharr[, 3, patch],
        col='red')
  legend('topright',
         legend=c('Farmer Pop. Growth',
                  'Mercenary Pop. Growth',
                  'Raider Pop. Growth'),
         lty=rep(1,3),
         lwd=rep(2,3),
         col=c('black', 'dodgerblue', 'red'),
         bty='n')
  dev.off()
}























open_fig <- function(fig_filepath, prms=prms, type){
  if(type == 'all'){
    system(paste0('for file in "',fig_filepath,'"/*.pdf; do open "${file}"; done'))
  }
  else{
    if(type == 'patch'){
      patchnumber <- readline(prompt = "Enter patch number:")
      system(paste0('for file in ',fig_filepath,'*"_patch_',patchnumber,'"*; do
                    open ${file}
                    done'))
    }else{
      for(i in 1:prms$npatch){
        if(type == 'pops'){
          path <- sprintf('%s/pop_patch_%s_plot.pdf',fig_filepath,i)
          system(paste0('open "',path, '"'))
        }
        if(type == 'costs'){
          path <- sprintf('%s/pop_patch_%s_plot.pdf',fig_filepath,i)
          system(paste0('open "',path, '"'))
          path <- sprintf('%s/pop_patch_%s_plot.pdf',fig_filepath,i)
    }
    
      }
    }
  }
}







#psplines is a WIP
library(splines2)

comp_costplot.splines <- function(cost_fp, patch, plotfp){
  
  if(dir.exists(sprintf('%s%s/component_costs',getwd(),plotfp))==FALSE){
    system(paste0('mkdir ',getwd(),plotfp,'/component_costs'))
  }
  
  #cost_fp <- sprintf('%s/%s',fp,test_costfp)
  cost.data <- load(file.path(sprintf('%s',cost_fp)), costplot.env <- new.env() )
  
  mod.prms <- costplot.env$component_cost.data$modified_prms
  orig.prms <- costplot.env$component_cost.data$original_prms
  
  modified_parameters <- orig.prms[sapply(names(mod.prms), function(x) !identical(orig.prms[[x]], mod.prms[[x]]))]
  modified_parameter_names <- names(modified_parameters)
  mod_par_str <- paste(modified_parameters,collapse=',')
  mod_par_names_str <- paste(modified_parameter_names,collapse=',')
  
  
  
  t <- seq(2,costplot.env$component_cost.data$original_prms$ngens,1)
  y <- costplot.env$component_cost.data$component_cost_array[,'benefit_to_f',patch]
  
  yf <- costplot.env$component_cost.data$component_cost_array[,'benefit_to_f',patch]
  ym <- costplot.env$component_cost.data$component_cost_array[,'benefit_to_m',patch]
  yr <- costplot.env$component_cost.data$component_cost_array[,'benefit_to_r',patch]
  
  #length(t)
  #length(y)
  pspline_fitf <- lm(yf ~ mSpline(x = t, df = 40, periodic = TRUE, data = data.frame(costplot.env$component_cost.data$component_cost_array)))
  pspline_fitm <- lm(ym ~ mSpline(x = t, df = 40, periodic = TRUE, data = data.frame(costplot.env$component_cost.data$component_cost_array)))
  pspline_fitr <- lm(yr ~ mSpline(x = t, df = 40, periodic = TRUE, data = data.frame(costplot.env$component_cost.data$component_cost_array)))
  
  predict_fitf <- predict(pspline_fitf, data.frame(x=t))
  predict_fitm <- predict(pspline_fitm, data.frame(x=t))
  predict_fitr <- predict(pspline_fitr, data.frame(x=t))
  
  
  plot(t, y, xlim=c(1,length(t)),
       #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
       ylim=c(min(predict_fitf,predict_fitm,predict_fitr),max(predict_fitf,predict_fitm,predict_fitr)), pch = 19, cex = 0.2, col = "darkgrey", main = sprintf("ben(+) of %s = %s -> 0 p%s",mod_par_names_str,mod_par_str,patch), xlab = "generation", ylab = sprintf("component cost %s = %s",mod_par_names_str,mod_par_str))
  #linear_model1 <- lm(y~t)
  #linear_model2 <- lm(y~poly(t,2,raw=TRUE))
  #lines(t, predict(linear_model1, data.frame(x=t)), col='green')
  #lines(t, predict(linear_model2, data.frame(x=t)), col='orange')
  lines(t, predict_fitf, col='darkgreen')
  lines(t, predict_fitm, col='dodgerblue')
  lines(t, predict_fitr, col='red')
  
  
  legend('right',
         legend=c('benefit_mtof',
                  'periodic spline fit raiders','periodic spline fit mercs','periodic spline fit farmers'),
         lty=rep(1,3),
         lwd=rep(2,3),
         col=c('darkgrey', 'red','dodgerblue','darkgreen'),
         bty='n')
  
  
  
  
  
  
  pdf(file=sprintf('%s%s/component_costs/compcost_splines_%s_patch_%s_plot.pdf',getwd(),plotfp,mod_par_names_str,patch), width=5, height=5)
  plot(t, y, xlim=c(1,length(t)),
       #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
       ylim=c(min(predict_fitf,predict_fitm,predict_fitr),max(predict_fitf,predict_fitm,predict_fitr)), pch = 19, cex = 0.2, col = "darkgrey", main = sprintf("ben (+) or cost (-) of %s = %s p%s",mod_par_names_str,mod_par_str,patch), xlab = "generation", ylab = sprintf("component cost %s = %s",mod_par_names_str,mod_par_str))
  #linear_model1 <- lm(y~t)
  #linear_model2 <- lm(y~poly(t,2,raw=TRUE))
  #lines(t, predict(linear_model1, data.frame(x=t)), col='green')
  #lines(t, predict(linear_model2, data.frame(x=t)), col='orange')
  lines(t, predict_fitf, col='darkgreen')
  lines(t, predict_fitm, col='dodgerblue')
  lines(t, predict_fitr, col='red')
  
  
  legend('right',
         legend=c('step-wise component cost values',
                  'periodic spline fit raiders','periodic spline fit mercs','periodic spline fit farmers'),
         lty=rep(1,3),
         lwd=rep(2,3),
         col=c('darkgrey', 'red','dodgerblue','darkgreen'),
         bty='n')
  
  dev.off()
  #original
  # #benefit m to f 
  # t <- seq(2,ccost_env$component_cost.data$original_prms$ngens,1)
  # y <- ccost_env$component_cost.data$component_cost_array[,'benefit_m_to_f',1]
  # 
  # #benefit m to r patch 11
  # yf <- ccost_env$component_cost.data$component_cost_array[,'benefit_m_to_f',12]
  # ym <- ccost_env$component_cost.data$component_cost_array[,'benefit_m_to_m',12]
  # yr <- ccost_env$component_cost.data$component_cost_array[,'benefit_m_to_r',12]
  # 
  # 
  # #length(t)
  # #length(y)
  # pspline_fitf <- lm(yf ~ mSpline(x = t, df = 15, periodic = TRUE, data = data.frame(ccost_env$component_cost.data$component_cost_array)))
  # pspline_fitm <- lm(ym ~ mSpline(x = t, df = 15, periodic = TRUE, data = data.frame(ccost_env$component_cost.data$component_cost_array)))
  # pspline_fitr <- lm(yr ~ mSpline(x = t, df = 15, periodic = TRUE, data = data.frame(ccost_env$component_cost.data$component_cost_array)))
  # 
  # 
  # plot(t, y, xlim=c(1,length(t)),
  #      #ylim=c(0,max(res[,c('num_f.p1','num_m.p1','num_r.p1')])+5),
  #      ylim=c(-5,5), pch = 19, cex = 0.2, col = "darkgrey", main = "ben (+) or cost (-) of mu.r .00055 -> 0 p12", xlab = "generation", ylab = "component cost mu.r 0.00055")
  # #linear_model1 <- lm(y~t)
  # #linear_model2 <- lm(y~poly(t,2,raw=TRUE))
  # #lines(t, predict(linear_model1, data.frame(x=t)), col='green')
  # #lines(t, predict(linear_model2, data.frame(x=t)), col='orange')
  # lines(t, predict(pspline_fitf, data.frame(x=t)), col='darkgreen')
  # lines(t, predict(pspline_fitm, data.frame(x=t)), col='dodgerblue')
  # lines(t, predict(pspline_fitr, data.frame(x=t)), col='red')
  # 
  # 
  # legend('topright',
  #        legend=c('step-wise component cost values',
  #                 'periodic spline fit raiders','periodic spline fit mercs','periodic spline fit farmers'),
  #        lty=rep(1,3),
  #        lwd=rep(2,3),
  #        col=c('darkgrey', 'red','dodgerblue','darkgreen'),
  #        bty='n')
  # 
  # 
  # 
  # 
  #                   
  # #Boundary.knots = c(2,3))
  
}
