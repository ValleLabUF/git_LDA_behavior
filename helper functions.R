df.to.list=function(dat) {  #only for id as col in dat
  id<- unique(dat$id)
  n=length(id)
  dat.list<- vector("list", n)
  names(dat.list)<- id
  
  for (i in 1:length(id)) {
    dat.list[[i]]<- dat[dat$id==id[i],]
  }
  dat.list
}
#------------------------------------------------
get.summary.stats_behav=function(dat,nbins){  #dat must have time.seg assigned; for all IDs
  
  #create list of input and to store output
  dat.list<- df.to.list(dat = dat)
  id<- unique(dat$id) %>% as.character()
  n<- length(id)
  obs.list<- vector("list", n)
  names(obs.list)<- id
  
  #calculate # of obs in each bin (per move param) by tseg
  for (i in 1:length(dat.list)) {
    dat.ind=dat.list[[i]]
    ntseg=max(dat.ind$tseg)
    
    
    #generate counts of obs per param bin by tseg for all params; 'id' and 'tseg' need to be 1st two cols
    mat.list=list()
    for (j in 1:length(nbins)) {
      mat.list[[j]]<- matrix(0, ntseg, nbins[j])
      colnames(mat.list[[j]])<- paste0("y",j,".",1:nbins[j])
      for (k in 1:ntseg){
        tmp<- dat.ind %>% filter(tseg == k) %>% dplyr::select(j+2) %>% table()
        mat.list[[j]][k,as.numeric(names(tmp))]<- tmp
      }
    }
    
    id<- rep(unique(dat.ind$id) %>% as.character(), ntseg)
    tseg<- 1:ntseg
    behav.res<- do.call(cbind.data.frame, mat.list) %>% data.frame() %>% cbind(id, tseg, .)
    behav.res$id<- as.character(behav.res$id)
    obs.list[[i]]<- behav.res
  }
  #obs<- do.call(rbind.data.frame, obs.list)
  obs<- map_dfr(obs.list, `[`)
  obs[is.na(obs)]<- 0  #replace NAs w/ zero
  obs
}

#------------------------------------------------
get_behav_hist=function(res,dat_red) {  #generate DF of bin counts for histogram per behavior
  
  #summarize cluster results by frequency and proportion
  behav.list<- list()
  for (j in 1:length(res$z.agg)) {
    behav.list[[j]]<- matrix(0, dim(res$z.agg[[j]])[2]*dim(res$z.agg[[j]])[3], 4)
    
    oo = 1
    for (i in 1:dim(res$z.agg[[j]])[3]) {
      tmp<- apply(res$z.agg[[j]][,,i], 2, sum) %>% data.frame(count = ., prop = ./sum(.)) %>%
        mutate(bin = 1:dim(res$z.agg[[j]])[2]) %>% mutate(behav = i) %>% as.matrix()
      
      behav.list[[j]][oo:(oo+dim(res$z.agg[[j]])[2] - 1),]<- tmp
      oo = oo+dim(res$z.agg[[j]])[2]
    }
    
    behav.list[[j]]<- data.frame(behav.list[[j]])
    names(behav.list[[j]])<- c("count","prop","bin","behav")
    behav.list[[j]]$param<- names(dat_red[j+2])
  }
  
  #combine params
  behav.res<- map_dfr(behav.list, `[`)
  
  behav.res
}
#------------------------------------------------
obs.per.tseg=function(dat.list) {  #count the number of observations w/in time segments
  tmp<- dat.list %>% group_by(tseg) %>% tally()
  tmp$id<- unique(dat.list$id)
  
  tmp
}
#------------------------------------------------
aug_behav_df=function(dat, theta.estim, nobs) {  #augment from time segments to observations
  
  for (i in 1:nrow(theta.estim)) {
    ind<- which(dat$id == theta.estim$id[i] & dat$tseg == theta.estim$tseg[i])
    
    if (i == 1) {
      theta.estim2<- rep(theta.estim[i,], nobs$n[i]) %>%
        unlist() %>%
        matrix(nrow = nobs$n[i], ncol = ncol(theta.estim), byrow = TRUE)
    } else {
      tmp<- rep(theta.estim[i,], nobs$n[i]) %>%
        unlist() %>%
        matrix(nrow = nobs$n[i], ncol = ncol(theta.estim), byrow = TRUE)
      
      theta.estim2<- rbind(theta.estim2, tmp)
    }
  }
  
  colnames(theta.estim2)<- names(theta.estim)
  theta.estim2<- data.frame(theta.estim2, time1 = dat$time1, date = dat$date,
                            stringsAsFactors = FALSE)
  
  theta.estim2
}
#------------------------------------------------
assign_behav=function(dat.list, theta.estim2) {  #assign dominant behavior to observations
  
  for (i in 1:length(dat.list)) {
    tmp<- matrix(NA, nrow(dat.list[[i]]), 2)
    sub<- theta.estim2[theta.estim2$id == unique(dat.list[[i]]$id),]
    
    for (j in 1:nrow(sub)) {
      k<- ncol(theta.estim2)-4  # number of behaviors
      ind<- which.max(sub[j,3:(3+k-1)])
      tmp[j,1]<- names(ind) %>% as.character()  #needs to be 'j+1' for simulations
      tmp[j,2]<- round(sub[j,(2+ind)], 3) %>% as.numeric() #needs to be 'j+1' for simulations
    }
    colnames(tmp)<- c("behav","prop")
    dat.list[[i]]<- cbind(dat.list[[i]], tmp)
    dat.list[[i]]$behav<- dat.list[[i]]$behav %>% as.character()
    dat.list[[i]]$prop<- dat.list[[i]]$prop %>% as.character()
  }
  
  #Convert to DF
  dat2<- do.call(rbind.data.frame, dat.list)
  dat2$prop<- dat2$prop %>% as.numeric() 
  
  dat2
}
