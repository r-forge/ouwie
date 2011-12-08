#Does contour plot for likelihood surface for pair of parameters

#written by Brian C. O'Meara

library(parallel)
library(akima)

OUwie.contour<-function(phy,data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"),root.station=TRUE, focal.param=NULL, clade=NULL, nrep=5000, sd.mult=3, parallel=TRUE,mc.cores=detectCores(), levels=c(0.5,1,1.5,2),...){
  #focal.param is something like c("alpha_2","sigma.sq_1"). They are then split on "_"
  if(length(focal.param)!=2) {
     stop("need a focal.param vector of length two")
  }
  if(sum(grepl("theta",focal.param))>0) {
    stop("contour mapping currently only works for alpha and sigma.sq parameters") 
  }
  globalMLE<-OUwie(phy=phy,data=data,model=model,root.station=root.station,clade=clade,plot.resid=FALSE,eigenvect=TRUE)
  focal.param.df<-data.frame(strsplit(focal.param,"_"),stringsAsFactors=FALSE)
  names(focal.param.df)<-c(1,2)
  
  focal.param.df<-rbind(focal.param.df,rep(NA,2))
  focal.param.df<-rbind(focal.param.df,rep(NA,2))
  row.names(focal.param.df)<-c("parameter","element","MLE","SE")
  for (i in 1:2) {
    focal.param.df[3,i]<-as.numeric(globalMLE$Param.est[which(row.names(globalMLE$Param.est)==focal.param.df[1,i]),focal.param.df[2,i]])
    focal.param.df[4,i]<-as.numeric(globalMLE$Param.SE[which(row.names(globalMLE$Param.SE)==focal.param.df[1,i]),focal.param.df[2,i]])
  }
  rnorm.bounded<-function(mean=0,sd=1,bound=0) {
    x<-(-Inf)
    while (x<=bound) {
      x<-rnorm(1,mean,sd) 
    }
    return(x)
  }
  
  #now get our random points, sampling most densely near the MLE
  param1.points<-replicate(n=round(nrep/2),expr=rnorm.bounded(mean=as.numeric(focal.param.df[3,1]),sd=sd.mult*as.numeric(focal.param.df[4,1])))
  param2.points<-replicate(n=round(nrep/2),expr=rnorm.bounded(mean=as.numeric(focal.param.df[3,2]),sd=sd.mult*as.numeric(focal.param.df[4,2])))
  param1.points<-c(param1.points,runif(nrep-round(nrep/2),min=min(param1.points),max=max(param1.points)))
  param2.points<-c(param2.points,runif(nrep-round(nrep/2),min=min(param2.points),max=max(param2.points)))
  param1.points<-sample(param1.points,size=length(param1.points),replace=FALSE)
  param2.points<-sample(param2.points,size=length(param2.points),replace=FALSE)
  
  dev<-function(p,phy,data,model,root.station,focal.param.vector,clade,globalMLE) { #globalMLE is just for figuring out the structure of alpha and sigma.sq
    nRegimes<-dim(globalMLE$index.matrix)[2]
    alpha<-rep(NA,nRegimes)
    sigma.sq<-rep(NA,nRegimes)
    for (freeParam in sequence(max(globalMLE$index.matrix))) {
      entries<-which(globalMLE$index.matrix==freeParam,arr.ind=TRUE)
      paramValue<-NA
      firstEntry<-entries[1,]
      paramNameRoot<-row.names(entries)[1]
      firstEntryName<-paste(paramNameRoot,firstEntry[2],sep="_")
      matchingSetParams<-which(firstEntryName==names(focal.param.vector))
      if (length(matchingSetParams)==1) {
        paramValue<-as.numeric(c(focal.param.vector[matchingSetParams]))
      }
      else {
        paramValue<-p[1]
        p<-p[-1] #pop off value we just used
      }
      for(i in sequence(dim(entries)[1])) {
        if(row.names(entries)[i]=="alpha") {
          alpha[entries[i,2]]<-paramValue
        }
        else {
          sigma.sq[entries[i,2]]<-paramValue
        }
      }
    }
    loglik<-OUwie.fixed(phy=phy,data=data, model=model,root.station=root.station, alpha=alpha, sigma.sq=sigma.sq, theta=NULL, clade=clade)$loglik
    return(loglik)
  }
  
  optimizeSemifixed<-function(X,phy,data,model,root.station,clade,globalMLE) {
    focal.param.vector<-X
     np<-max(globalMLE$index.matrix)-length(focal.param.vector)
      lower = rep(0.0000000001, np)
	    upper = rep(20, np)
	    ip<-1
	  opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5, "xtol_rel"=.Machine$double.eps^0.5)
	  out = nloptr(x0=rep(ip, length.out = np), eval_f=dev, lb=lower, ub=upper, opts=opts,phy=phy,data=data, model=model,root.station=root.station, clade=clade, globalMLE=globalMLE,focal.param.vector=focal.param.vector)
     return(-1*out$objective)
  }
  
  params.points<-data.frame(param1.points,param2.points)
  names(params.points)<-focal.param
  params.points.list<-split(params.points,row(params.points),drop=TRUE)
  if(!parallel) {
    mc.cores<-1
  }
  #likelihoods<-(-1*simplify2array(lapply(params.points.list,optimizeSemifixed,phy=phy,data=data, model=model,root.station=root.station, clade=clade, globalMLE=globalMLE)))
  likelihoods<-(-1*simplify2array(lapply(params.points.list,optimizeSemifixed,phy=phy,data=data, model=model,root.station=root.station, clade=clade, globalMLE=globalMLE, mc.cores=mc.cores)))

  #include the MLE in the set
  likelihoods<-c(likelihoods,globalMLE$loglik)
  param1.points<-c(param1.points,as.numeric(focal.param.df[3,1]))
  param2.points<-c(param2.points,as.numeric(focal.param.df[3,2]))
  
  contour(interp(x= param1.points,y=param2.points,z=likelihoods-min(likelihoods),xo=seq(min(param1.points), max(param1.points), length = 400),yo=seq(min(param1.points), max(param1.points), length = 400),duplicate="strip"),levels=levels, ...)
  if(strsplit(focal.param,"_")[[1]][1] == strsplit(focal.param,"_")[[2]][1]) {
    plot.range<-range(c(param1.points,param2.points))
    lines(x= plot.range,y= plot.range,lty="dotted")
  }
  points(x=as.numeric(focal.param.df[3,1]),y=as.numeric(focal.param.df[3,2]),pch=20,col="red")
  print(paste("MLE is at c(",as.numeric(focal.param.df[3,1]),",",as.numeric(focal.param.df[3,2]),")"))
  lines(x=c(as.numeric(focal.param.df[3,1])-1.96*as.numeric(focal.param.df[4,1]),as.numeric(focal.param.df[3,1])+1.96*as.numeric(focal.param.df[4,1])),y=rep(as.numeric(focal.param.df[3,2]),2))
  lines(x=rep(as.numeric(focal.param.df[3,1]),2),y=c(as.numeric(focal.param.df[3,2])-1.96*as.numeric(focal.param.df[4,2]),as.numeric(focal.param.df[3,2])+1.96*as.numeric(focal.param.df[4,2])))
  
  
  finalResults<-data.frame(param1.points,param2.points,-1*likelihoods)
  names(finalResults)<-c(focal.param,"loglik")
  return(finalResults)             
}
