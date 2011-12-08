#Does contour plot for likelihood surface for pair of parameters

#written by Brian C. O'Meara

OUwie.contour<-function(phy,data, model=c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA"),root.station=TRUE, focal.param=NULL, clade=NULL, nrep=5000, sd.mult=4,...){
  #focal.param is something like c("alpha_2","sigma.sq_1"). They are then split on "_"
  if(length(focal.param)!=2) {
     stop("need a focal.param vector of length two")
  }
  if(sum(grepl("theta",focal.param))>0) {
    stop("contour mapping currently only works for alpha and sigma.sq parameters") 
  }
  globalMLE<-OUwie(phy=phy,data=data,model=model,root.station=root.station,clade=clade)
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
  
  #NEED TO ADD FUNCTION HERE TO OPTIMIZE WITH PARTIALLY FIXED POINTS
}
