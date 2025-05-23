foo1 <- function(x){
  for( i in x ){
    #  require returns TRUE invisibly if it was able to load package
    if( ! require( i , character.only = TRUE ) ){
      #  If package was not able to be loaded then re-install
      install.packages( i , dependencies = TRUE )
      #  Load package after installing
      require( i , character.only = TRUE )
    }
  }
}

is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

# fixing np int print issue
foo4 <- function(df){
  gsub(pattern = "np.int64\\(|\\)", replacement = "", x = df)
}

# extract fitness values and genomes of Pareto solutions
get_pareto <- function(){
  
  ## define list of two
  pareto <- setNames(replicate(2,data.frame()),c('fitness', 'genomes'))
  
  setwd(paste0(path,'/output'))

  logfile <- rev(dir(getwd(), pattern="_log", full.names = T))
  logfile_short <- rev(dir(getwd(), pattern="_log", full.names = F))
  dates <- unlist(strsplit(logfile_short, "_optimization_log.txt"))
  datetime <- as.POSIXct(dates,format="%d-%m-%Y_%H-%M-%S",tz=Sys.timezone())
  
  logfile <- file(logfile[which.max(datetime)],"r+")
  text<-readLines(logfile)
  close(logfile)
  
  if (any(grep("The optimization process needed", text) > 0)) {
    seconds <- strsplit(text[(grep("The optimization process needed", text))], "[|]") [[1]] [2]
  }
  
  if (any(grep("Best Solutions", text) > 0)) {
    bestsol <- as.matrix(text[(grep("Best Solutions", text)+2):(length(text)-1)])
    infeasible <- grep("infeasible", bestsol)
    if (length(infeasible)>0) {
      feasible <- as.matrix(c(1:length(bestsol))[-infeasible])
      bestsol <- as.matrix(bestsol[-infeasible])
    } else (feasible <- c(1:length(bestsol)))
    if (length(bestsol)>0) {
      for(k in 1:length(bestsol)){
        bestsol[k] <- strsplit(bestsol, "[|]") [[k]] [2]
      }
    }
    bestind <- bestsol
    bestfit <- bestsol
    if (length(bestsol)>0) {
      for(k in 1:length(bestsol)){
        bestind[k] <- strsplit(bestsol, "[:]") [[k]] [1]
        bestind[k] <- substr(bestind[k], 3, nchar(bestind[k])-2) 
        bestfit[k] <- strsplit(bestsol, "[:]") [[k]] [2]
        bestfit[k] <- substr(bestfit[k], 3, nchar(bestfit[k])-1) 
      }
    }
  }
  
  
  
  
  bestfit_df <- data.frame(do.call(rbind, strsplit(bestfit, ", ", fixed=TRUE)))
  bestfit_df <- data.frame(lapply(bestfit_df,as.numeric))
  names(bestfit_df) <- paste0('fit',1:length(bestfit_df))
  
  
  bestind_df <- data.frame(do.call(cbind, strsplit(bestind, ", ", fixed=TRUE)))
  bestind_df <- data.frame(lapply(bestind_df, foo4))
  bestind_df <- data.frame(lapply(bestind_df,as.numeric))
  names(bestind_df) <- paste0('ind',1:length(bestind_df))
  
  if(length(bestsol)>0){
    setwd(paste0(path,'/output_analysis'))
    write.table(file="pareto_genomes.txt",bestind_df, col.names=F, row.names=F, quote=F)
    write.table(file="pareto_fitness.txt",bestfit_df, col.names=F, row.names=F, quote=F)
  }
  
  pareto[[1]] <- bestfit_df
  pareto[[2]] <- bestind_df
  
  return(pareto)
}

write_sq <- function(){
  
  setwd(paste0(path,'/output'))
  
  nobj = length(pareto$fitness)
  
  ind_file <- rev(dir(getwd(), pattern = 'indiv', full.names = T))
  logfile_short <- rev(dir(getwd(), pattern="_log", full.names = F))
  dates <- unlist(strsplit(logfile_short, "_optimization_log.txt"))
  datetime <- as.POSIXct(dates,format="%d-%m-%Y_%H-%M-%S",tz=Sys.timezone())
  ind_all <- read.csv(ind_file[which.max(datetime)], h=F, skip=1, as.is=T)
  
  if(nobj == 2){
    ind_all.fitness <- ind_all %>%
      mutate(V3 = as.numeric(gsub('\\[', '', V3)),
             V4 = as.numeric(gsub('\\]', '', V4))) %>%
      mutate(fit1 = V3,
             fit2 = V4) %>% 
      dplyr::select(fit1, fit2)
    sq_fitness <- ind_all.fitness[1,]
    setwd(paste0(path,'/output_analysis'))
    write.table(file="sq_fitness.txt",sq_fitness, col.names=F, row.names=F, quote=F)
  }
  
  if(nobj == 3){
    ind_all.fitness <- ind_all %>%
      mutate(V3 = as.numeric(gsub('\\[', '', V3)),
             V5 = as.numeric(gsub('\\]', '', V5))) %>%
      mutate(fit1 = V3,
             fit2 = V4,
             fit3 = V5) %>% 
      dplyr::select(fit1, fit2, fit3)
    sq_fitness <- ind_all.fitness[1,]
    setwd(paste0(path,'/output_analysis'))
    write.table(file="sq_fitness.txt",sq_fitness, col.names=F, row.names=F, quote=F)
  }
  
  if(nobj == 4){
    ind_all.fitness <- ind_all %>%
      mutate(V3 = as.numeric(gsub('\\[', '', V3)),
             V6 = as.numeric(gsub('\\]', '', V6))) %>%
      mutate(fit1 = V3,
             fit2 = V4,
             fit3 = V5,
             fit4 = V6) %>% 
      dplyr::select(fit1, fit2, fit3, fit4)
    sq_fitness <- ind_all.fitness[1,]
    setwd(paste0(path,'/output_analysis'))
    write.table(file="sq_fitness.txt",sq_fitness, col.names=F, row.names=F, quote=F)
  }
}

# calculate hypervolumes for each generation
hv_generations <- function(){
  
  setwd(paste0(path,'/output'))
  
  nobj = length(pareto$fitness)
  
  ind_file <- rev(dir(getwd(), pattern = 'indiv', full.names = T))
  logfile_short <- rev(dir(getwd(), pattern="_log", full.names = F))
  dates <- unlist(strsplit(logfile_short, "_optimization_log.txt"))
  datetime <- as.POSIXct(dates,format="%d-%m-%Y_%H-%M-%S",tz=Sys.timezone())
  ind_all <- read.csv(ind_file[which.max(datetime)], h=F, skip=1, as.is=T)
  
  if(nobj == 2){
    ind_all.fitness <- ind_all %>%
      mutate(V3 = as.numeric(gsub('\\[', '', V3)),
             V4 = as.numeric(gsub('\\]', '', V4)))
    
    if(ind_all.fitness$V3[which.max(ind_all.fitness$V3)]==0) ind_all.fitness$V3[which(ind_all.fitness$V3==ind_all.fitness$V3[which.max(ind_all.fitness$V3)])] <- -0.00001
    if(ind_all.fitness$V4[which.max(ind_all.fitness$V4)]==0) ind_all.fitness$V4[which(ind_all.fitness$V4==ind_all.fitness$V4[which.max(ind_all.fitness$V4)])] <- -0.00001

    # if values range from - to +, rearrange to - only by subtracting the max value
    min_V3 <- min(ind_all.fitness$V3)
    min_V4 <- min(ind_all.fitness$V4)
    max_V3 <- max(ind_all.fitness$V3)
    max_V4 <- max(ind_all.fitness$V4)
    
    ind_all.fitness <- ind_all.fitness %>%
      dplyr::rowwise() %>%
      mutate(fit1 = ifelse(min_V3<0 & max_V3>0, V3-max_V3-0.00001, V3),
             fit2 = ifelse(min_V4<0 & max_V4>0, V4-max_V4-0.00001, V4)) %>% 
      dplyr::select(fit1, fit2)
    
    # scale dimensions by their respective maximum value (best = -1, worst = 0)
    min_fit1 <- min(ind_all.fitness$fit1)
    min_fit2 <- min(ind_all.fitness$fit2)
    max_fit1 <- max(ind_all.fitness$fit1)
    max_fit2 <- max(ind_all.fitness$fit2)
    
    ind_all.fitness <- ind_all.fitness %>%
      dplyr::rowwise() %>%
      mutate(fit1 = ifelse(max_fit1>=0, (fit1 - min_fit1)/(max_fit1 - min_fit1)*-1, (fit1 - min_fit1)/(min_fit1 - max_fit1)),
             fit2 = ifelse(max_fit2>=0, (fit2 - min_fit2)/(max_fit2 - min_fit2)*-1, (fit2 - min_fit2)/(min_fit2 - max_fit2))) %>% 
      dplyr::select(fit1, fit2)
    
    # ref point must be vector of zeros
    ref <- matrix(data=0,nrow=1, ncol=2)
    maxgen <- max(ind_all$V1)+1
    HV <- data.frame(matrix(data=NA,nrow=maxgen,ncol=1))
    
    pF1 <- paretoFilter(as.matrix(ind_all.fitness[ind_all$V1==0,]))
    HV[1,1] <- dominatedHypervolume(pF1, ref)
    
    k=1
    
    for(l in 2:maxgen){
      pF <- paretoFilter(as.matrix(ind_all.fitness[ind_all$V1==k,]))
      pF_ <-paretoFilter(as.matrix(rbind(pF,pF1)))
      pF_ <- pF_[!duplicated(pF_),]
      HV[l,1] <- dominatedHypervolume(pF_, ref)
      pF1 <- pF_
      k=k+1
    }
  }
  
  if(nobj == 3){
    ind_all.fitness <- ind_all %>%
      mutate(V3 = as.numeric(gsub('\\[', '', V3)),
             V5 = as.numeric(gsub('\\]', '', V5)))
    
    if(ind_all.fitness$V3[which.max(ind_all.fitness$V3)]==0) ind_all.fitness$V3[which(ind_all.fitness$V3==ind_all.fitness$V3[which.max(ind_all.fitness$V3)])] <- -0.00001
    if(ind_all.fitness$V4[which.max(ind_all.fitness$V4)]==0) ind_all.fitness$V4[which(ind_all.fitness$V4==ind_all.fitness$V4[which.max(ind_all.fitness$V4)])] <- -0.00001
    if(ind_all.fitness$V5[which.max(ind_all.fitness$V5)]==0) ind_all.fitness$V5[which(ind_all.fitness$V5==ind_all.fitness$V5[which.max(ind_all.fitness$V5)])] <- -0.00001

    # if values range from - to +, rearrange to - only by subtracting the max value
    min_V3 <- min(ind_all.fitness$V3)
    min_V4 <- min(ind_all.fitness$V4)
    min_V5 <- min(ind_all.fitness$V5)
    max_V3 <- max(ind_all.fitness$V3)
    max_V4 <- max(ind_all.fitness$V4)
    max_V5 <- max(ind_all.fitness$V5)
    
    ind_all.fitness <- ind_all.fitness %>%
      dplyr::rowwise() %>%
      mutate(fit1 = ifelse(min_V3<0 & max_V3>0, V3-max_V3-0.00001, V3),
             fit2 = ifelse(min_V4<0 & max_V4>0, V4-max_V4-0.00001, V4),
             fit3 = ifelse(min_V5<0 & max_V5>0, V5-max_V5-0.00001, V5)) %>% 
      dplyr::select(fit1, fit2, fit3)
    
    # scale dimensions by their respective maximum value (best = -1, worst = 0)
    min_fit1 <- min(ind_all.fitness$fit1)
    min_fit2 <- min(ind_all.fitness$fit2)
    min_fit3 <- min(ind_all.fitness$fit3)
    max_fit1 <- max(ind_all.fitness$fit1)
    max_fit2 <- max(ind_all.fitness$fit2)
    max_fit3 <- max(ind_all.fitness$fit3)
    
    ind_all.fitness <- ind_all.fitness %>%
      dplyr::rowwise() %>%
      mutate(fit1 = ifelse(max_fit1>=0, (fit1 - min_fit1)/(max_fit1 - min_fit1)*-1, (fit1 - min_fit1)/(min_fit1 - max_fit1)),
             fit2 = ifelse(max_fit2>=0, (fit2 - min_fit2)/(max_fit2 - min_fit2)*-1, (fit2 - min_fit2)/(min_fit2 - max_fit2)),
             fit3 = ifelse(max_fit3>=0, (fit3 - min_fit3)/(max_fit3 - min_fit3)*-1, (fit3 - min_fit3)/(min_fit3 - max_fit3))) %>% 
      dplyr::select(fit1, fit2, fit3)
    
    # ref point must be vector of zeros
    ref <- matrix(data=0,nrow=1, ncol=3)
    maxgen <- max(ind_all$V1)+1
    HV <- data.frame(matrix(data=NA,nrow=maxgen,ncol=1))
    
    pF1 <- paretoFilter(as.matrix(ind_all.fitness[ind_all$V1==0,]))
    HV[1,1] <- dominatedHypervolume(pF1, ref)
    
    k=1
    
    for(l in 2:maxgen){
      pF <- paretoFilter(as.matrix(ind_all.fitness[ind_all$V1==k,]))
      pF_ <-paretoFilter(as.matrix(rbind(pF,pF1)))
      pF_ <- pF_[!duplicated(pF_),]
      HV[l,1] <- dominatedHypervolume(pF_, ref)
      pF1 <- pF_
      k=k+1
    }
  }
  
  if(nobj == 4){
    ind_all.fitness <- ind_all %>%
      mutate(V3 = as.numeric(gsub('\\[', '', V3)),
             V6 = as.numeric(gsub('\\]', '', V6)))
    
    # if maximum is 0, set to -0.00001
    if(ind_all.fitness$V3[which.max(ind_all.fitness$V3)]==0) ind_all.fitness$V3[which(ind_all.fitness$V3==ind_all.fitness$V3[which.max(ind_all.fitness$V3)])] <- -0.00001
    if(ind_all.fitness$V4[which.max(ind_all.fitness$V4)]==0) ind_all.fitness$V4[which(ind_all.fitness$V4==ind_all.fitness$V4[which.max(ind_all.fitness$V4)])] <- -0.00001
    if(ind_all.fitness$V5[which.max(ind_all.fitness$V5)]==0) ind_all.fitness$V5[which(ind_all.fitness$V5==ind_all.fitness$V5[which.max(ind_all.fitness$V5)])] <- -0.00001
    if(ind_all.fitness$V6[which.max(ind_all.fitness$V6)]==0) ind_all.fitness$V6[which(ind_all.fitness$V6==ind_all.fitness$V6[which.max(ind_all.fitness$V6)])] <- -0.00001
    
    # if values range from - to +, rearrange to - only by subtracting the max value
    min_V3 <- min(ind_all.fitness$V3)
    min_V4 <- min(ind_all.fitness$V4)
    min_V5 <- min(ind_all.fitness$V5)
    min_V6 <- min(ind_all.fitness$V6)
    max_V3 <- max(ind_all.fitness$V3)
    max_V4 <- max(ind_all.fitness$V4)
    max_V5 <- max(ind_all.fitness$V5)
    max_V6 <- max(ind_all.fitness$V6)
    
    ind_all.fitness <- ind_all.fitness %>%
      dplyr::rowwise() %>%
      mutate(fit1 = ifelse(min_V3<0 & max_V3>0, V3-max_V3-0.00001, V3),
             fit2 = ifelse(min_V4<0 & max_V4>0, V4-max_V4-0.00001, V4),
             fit3 = ifelse(min_V5<0 & max_V5>0, V5-max_V5-0.00001, V5),
             fit4 = ifelse(min_V6<0 & max_V6>0, V6-max_V6-0.00001, V6)) %>% 
      dplyr::select(fit1, fit2, fit3, fit4)
    
    # scale dimensions by their respective maximum value (best = -1, worst = 0)
    min_fit1 <- min(ind_all.fitness$fit1)
    min_fit2 <- min(ind_all.fitness$fit2)
    min_fit3 <- min(ind_all.fitness$fit3)
    min_fit4 <- min(ind_all.fitness$fit4)
    max_fit1 <- max(ind_all.fitness$fit1)
    max_fit2 <- max(ind_all.fitness$fit2)
    max_fit3 <- max(ind_all.fitness$fit3)
    max_fit4 <- max(ind_all.fitness$fit4)
    
    ind_all.fitness <- ind_all.fitness %>%
      dplyr::rowwise() %>%
      mutate(fit1 = ifelse(max_fit1>=0, (fit1 - min_fit1)/(max_fit1 - min_fit1)*-1, (fit1 - min_fit1)/(min_fit1 - max_fit1)),
             fit2 = ifelse(max_fit2>=0, (fit2 - min_fit2)/(max_fit2 - min_fit2)*-1, (fit2 - min_fit2)/(min_fit2 - max_fit2)),
             fit3 = ifelse(max_fit3>=0, (fit3 - min_fit3)/(max_fit3 - min_fit3)*-1, (fit3 - min_fit3)/(min_fit3 - max_fit3)),
             fit4 = ifelse(max_fit4>=0, (fit4 - min_fit4)/(max_fit4 - min_fit4)*-1, (fit4 - min_fit4)/(min_fit4 - max_fit4)))
    
    # ref point must be vector of zeros
    ref <- matrix(data=0,nrow=1, ncol=4)

    maxgen <- max(ind_all$V1)+1
    HV <- data.frame(matrix(data=NA,nrow=maxgen,ncol=1))
    
    pF1 <- paretoFilter(as.matrix(ind_all.fitness[ind_all$V1==0,]))
    HV[1,1] <- dominatedHypervolume(pF1, ref)
    
    k=1
    
    for(l in 2:maxgen){
      pF <- paretoFilter(as.matrix(ind_all.fitness[ind_all$V1==k,]))
      pF_ <-paretoFilter(as.matrix(rbind(pF,pF1)))
      pF_ <- pF_[!duplicated(pF_),]
      HV[l,1] <- dominatedHypervolume(pF_, ref)
      pF1 <- pF_
      k=k+1
    }
  }
  
  setwd(paste0(path,'/output_analysis'))
  names(HV) <- 'HV'
  HV$Generation <- c(1:length(HV$HV))
  return(HV)
}

# plot in 2D 
plot_2D <- function(mode=1){
  setwd(paste0(path,'/output_analysis'))
  sol <- pareto$fitness
  
  # plot for 2 objectives
  if(length(sol)==2){
    names(sol) <- c("obj1","obj2")
    if(mode==1){
      p <- ggplot(sol, aes(x=obj1, y=obj2)) +
        geom_point(aes(), shape=21, col="black") +
        xlab(fit1)+
        ylab(fit2)
    }
    if(mode==2){
      p <- ggplot(sol, aes(x=obj1, y=obj2)) +
        geom_point(aes(), shape=21, col="black") +
        xlab(fit2)+
        ylab(fit1)
    }
    p
  }
  
  # plot for 3 objectives
  if(length(sol)==3){
    names(sol) <- c("obj1","obj2","obj3")
    if(mode==1){
      p <- ggplot(sol, aes(x=obj1, y=obj2)) + 
        geom_point(aes(fill=obj3), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(fill = guide_legend(title=fit3)) +
        xlab(fit1)+
        ylab(fit2)
    }
    if(mode==2){
      p <- ggplot(sol, aes(x=obj1, y=obj3)) + 
        geom_point(aes(fill=obj2), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(fill = guide_legend(title=fit2)) +
        xlab(fit1)+
        ylab(fit3)
    }
    if(mode==3){
      p <- ggplot(sol, aes(x=obj2, y=obj1)) + 
        geom_point(aes(fill=obj3), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(fill = guide_legend(title=fit3)) +
        xlab(fit2)+
        ylab(fit1)
    }
    if(mode==4){
      p <- ggplot(sol, aes(x=obj2, y=obj3)) + 
        geom_point(aes(fill=obj1), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(fill = guide_legend(title=fit1)) +
        xlab(fit2)+
        ylab(fit3)
    }
    if(mode==5){
      p <- ggplot(sol, aes(x=obj3, y=obj1)) + 
        geom_point(aes(fill=obj2), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(fill = guide_legend(title=fit2)) +
        xlab(fit3)+
        ylab(fit1)
    }
    if(mode==6){
      p <- ggplot(sol, aes(x=obj3, y=obj2)) + 
        geom_point(aes(fill=obj1), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(fill = guide_legend(title=fit1)) +
        xlab(fit3)+
        ylab(fit2)
    }
    p
  }
  
  # plot for 4 objectives
  if(length(sol)==4){
    names(sol) <- c("obj1","obj2","obj3","obj4")
    if(mode==1){
      p <- ggplot(sol, aes(x=obj1, y=obj2)) + 
        geom_point(aes(size=obj4, fill=obj3), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit4, override.aes = list(col = "black")),
               fill = guide_legend(title=fit3)) +
        xlab(fit1)+
        ylab(fit2)
    }
    if(mode==1){
      p <- ggplot(sol, aes(x=obj1, y=obj2)) + 
        geom_point(aes(size=obj4, fill=obj3), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit4, override.aes = list(col = "black")),
               fill = guide_legend(title=fit3)) +
        xlab(fit1)+
        ylab(fit2)
    }
    if(mode==2){
      p <- ggplot(sol, aes(x=obj1, y=obj2)) + 
        geom_point(aes(size=obj3, fill=obj4), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit3, override.aes = list(col = "black")),
               fill = guide_legend(title=fit4)) +
        xlab(fit1)+
        ylab(fit2)
    }
    if(mode==3){
      p <- ggplot(sol, aes(x=obj1, y=obj3)) + 
        geom_point(aes(size=obj2, fill=obj4), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit2, override.aes = list(col = "black")),
               fill = guide_legend(title=fit4)) +
        xlab(fit1)+
        ylab(fit3)
    }
    if(mode==4){
      p <- ggplot(sol, aes(x=obj1, y=obj3)) + 
        geom_point(aes(size=obj4, fill=obj2), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit4, override.aes = list(col = "black")),
               fill = guide_legend(title=fit2)) +
        xlab(fit1)+
        ylab(fit3)
    }
    if(mode==5){
      p <- ggplot(sol, aes(x=obj1, y=obj4)) + 
        geom_point(aes(size=obj2, fill=obj3), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit2, override.aes = list(col = "black")),
               fill = guide_legend(title=fit3)) +
        xlab(fit1)+
        ylab(fit4)
    }
    if(mode==6){
      p <- ggplot(sol, aes(x=obj1, y=obj4)) + 
        geom_point(aes(size=obj3, fill=obj2), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit3, override.aes = list(col = "black")),
               fill = guide_legend(title=fit2)) +
        xlab(fit1)+
        ylab(fit4)
    }
    if(mode==7){
      p <- ggplot(sol, aes(x=obj2, y=obj1)) + 
        geom_point(aes(size=obj3, fill=obj4), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit3, override.aes = list(col = "black")),
               fill = guide_legend(title=fit4)) +
        xlab(fit2)+
        ylab(fit1)
    }
    if(mode==8){
      p <- ggplot(sol, aes(x=obj2, y=obj1)) + 
        geom_point(aes(size=obj4, fill=obj3), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit4, override.aes = list(col = "black")),
               fill = guide_legend(title=fit3)) +
        xlab(fit2)+
        ylab(fit1)
    }
    if(mode==9){
      p <- ggplot(sol, aes(x=obj2, y=obj3)) + 
        geom_point(aes(size=obj1, fill=obj4), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit1, override.aes = list(col = "black")),
               fill = guide_legend(title=fit4)) +
        xlab(fit2)+
        ylab(fit3)
    }
    if(mode==10){
      p <- ggplot(sol, aes(x=obj2, y=obj3)) + 
        geom_point(aes(size=obj4, fill=obj1), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit4, override.aes = list(col = "black")),
               fill = guide_legend(title=fit1)) +
        xlab(fit2)+
        ylab(fit3)
    }
    if(mode==11){
      p <- ggplot(sol, aes(x=obj2, y=obj4)) + 
        geom_point(aes(size=obj1, fill=obj3), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit1, override.aes = list(col = "black")),
               fill = guide_legend(title=fit3)) +
        xlab(fit2)+
        ylab(fit4)
    }
    if(mode==12){
      p <- ggplot(sol, aes(x=obj2, y=obj4)) + 
        geom_point(aes(size=obj3, fill=obj1), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit3, override.aes = list(col = "black")),
               fill = guide_legend(title=fit1)) +
        xlab(fit2)+
        ylab(fit4)
    }
    if(mode==13){
      p <- ggplot(sol, aes(x=obj3, y=obj1)) + 
        geom_point(aes(size=obj2, fill=obj4), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit2, override.aes = list(col = "black")),
               fill = guide_legend(title=fit4)) +
        xlab(fit3)+
        ylab(fit1)
    }
    if(mode==14){
      p <- ggplot(sol, aes(x=obj3, y=obj1)) + 
        geom_point(aes(size=obj4, fill=obj2), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit4, override.aes = list(col = "black")),
               fill = guide_legend(title=fit2)) +
        xlab(fit3)+
        ylab(fit1)
    }
    if(mode==15){
      p <- ggplot(sol, aes(x=obj3, y=obj2)) + 
        geom_point(aes(size=obj1, fill=obj4), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit1, override.aes = list(col = "black")),
               fill = guide_legend(title=fit4)) +
        xlab(fit3)+
        ylab(fit2)
    }
    if(mode==16){
      p <- ggplot(sol, aes(x=obj3, y=obj2)) + 
        geom_point(aes(size=obj4, fill=obj1), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit4, override.aes = list(col = "black")),
               fill = guide_legend(title=fit1)) +
        xlab(fit3)+
        ylab(fit2)
    }
    if(mode==17){
      p <- ggplot(sol, aes(x=obj3, y=obj4)) + 
        geom_point(aes(size=obj2, fill=obj1), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit2, override.aes = list(col = "black")),
               fill = guide_legend(title=fit1)) +
        xlab(fit3)+
        ylab(fit4)
    }
    if(mode==18){
      p <- ggplot(sol, aes(x=obj3, y=obj4)) + 
        geom_point(aes(size=obj2, fill=obj1), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit2, override.aes = list(col = "black")),
               fill = guide_legend(title=fit1)) +
        xlab(fit3)+
        ylab(fit4)
    }
    if(mode==19){
      p <- ggplot(sol, aes(x=obj4, y=obj1)) + 
        geom_point(aes(size=obj2, fill=obj3), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit2, override.aes = list(col = "black")),
               fill = guide_legend(title=fit3)) +
        xlab(fit4)+
        ylab(fit1)
    }
    if(mode==20){
      p <- ggplot(sol, aes(x=obj4, y=obj1)) + 
        geom_point(aes(size=obj3, fill=obj2), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit3, override.aes = list(col = "black")),
               fill = guide_legend(title=fit2)) +
        xlab(fit4)+
        ylab(fit1)
    }   
    if(mode==21){
      p <- ggplot(sol, aes(x=obj4, y=obj2)) + 
        geom_point(aes(size=obj1, fill=obj3), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit1, override.aes = list(col = "black")),
               fill = guide_legend(title=fit3)) +
        xlab(fit4)+
        ylab(fit2)
    }    
    if(mode==22){
      p <- ggplot(sol, aes(x=obj4, y=obj2)) + 
        geom_point(aes(size=obj3, fill=obj1), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit3, override.aes = list(col = "black")),
               fill = guide_legend(title=fit1)) +
        xlab(fit4)+
        ylab(fit2)
    } 
    if(mode==23){
      p <- ggplot(sol, aes(x=obj4, y=obj3)) + 
        geom_point(aes(size=obj1, fill=obj2), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit1, override.aes = list(col = "black")),
               fill = guide_legend(title=fit2)) +
        xlab(fit4)+
        ylab(fit3)
    }    
    if(mode==24){
      p <- ggplot(sol, aes(x=obj4, y=obj3)) + 
        geom_point(aes(size=obj2, fill=obj1), shape=21, col="black") +
        scale_fill_gradientn(colours=viridis(100)) +
        guides(size = guide_legend(title=fit2, override.aes = list(col = "black")),
               fill = guide_legend(title=fit1)) +
        xlab(fit4)+
        ylab(fit3)
    }
  }
  p
}

