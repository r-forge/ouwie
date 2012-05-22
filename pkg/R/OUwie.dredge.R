OUwie.dredge<-function(phy,data, criterion=c("aicc","aic","rjmcmc"), theta.max.k=3, sigma.max.k=3, alpha.max.k=3, root.station=TRUE, ip=1, lb=0.000001, ub=1000, max.generations=100, wait.generations=10, pop.size=NULL) {
  #criterion: If aicc or aic, use rgenoud, where the fitness is the score (and try to minimize)
  #    if rjmcmc, use an rjmcmc approach
  #alpha.max.k=3: allows for three OU regimes and one regime where alpha is set to ~0 (brownian motion)
  #    if alpha.max.k=0, only BM models are investigated
  #sigma.max.k, theta.max.k: work in same way, but must be at least 1 in each case
  #note that these params can be stuck together in interesting ways. For example, a BM jump model has multiple thetas but constant everything else
  #There are as many free parameters for rgenoud as the number of nodes (not edges) * 3: params for theta, sigma, and alpha
  #so each individual for rgenoud has theta1, theta2, theta3....,sigma1, sigma2, sigma3, ....,alpha1, alpha2, alpha3, where theta1 is the mapping to the theta free parameter for node 1, and so forth
  if(criterion=="aicc") {
   if (is.null(pop.size)) {
     pop.size<-___FUNCTION TO SET POP.SIZE______ 
    }
   starting.values<-matrix(c(rep(1,Nnode(phy,internal.only=FALSE),rep(1,Nnode(phy,internal.only=FALSE),rep(0,Nnode(phy,internal.only=FALSE)),nrow=1,ncol=3*Nnode(phy,internal.only=FALSE) #BM1
   starting.values<-rbind(starting.values,matrix(c(rep(1,Nnode(phy,internal.only=FALSE),rep(1,Nnode(phy,internal.only=FALSE),rep(1,Nnode(phy,internal.only=FALSE)),nrow=1,ncol=3*Nnode(phy,internal.only=FALSE)) #OU1
   results<-genoud(fn=dredge.aicc, starting.values=starting.values, max=FALSE,nvars=3*Nnode(phy,internal.only=FALSE, data.type.int=TRUE,boundary.enforcement=2,Domains=matrix(c(rep(1,Nnode(phy,internal.only=FALSE),rep(1,Nnode(phy,internal.only=FALSE),rep(0,Nnode(phy,internal.only=FALSE),rep(theta.max.k,Nnode(phy,internal.only=FALSE),rep(sigma.max.k,Nnode(phy,internal.only=FALSE),rep(alpha.max.k,Nnode(phy,internal.only=FALSE)),ncol=2,nrow=3*Nnode(phy,internal.only=FALSE,byrow=FALSE), wait.generations=wait.generations, max.generations=max.generations, pop.size=pop.size, root.station=TRUE, ip=ip, lb=lb, ub=ub, phy, data)
  }
}

dredge.aicc<-function(rgenoud.individual, phy, data, ...) {
  #convert phy+regenoud.individual to simmap.tree (later, make it so that we directly go to the proper object)
  #call dev.dredge
  #convert likelihood to AICC
  #return AICC
}




dev.dredge<-function(p,simmap.tree,data,model,root.station,focal.param.vector,clade,globalMLE) { #globalMLE is just for figuring out the structure of alpha and sigma.sq
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
  loglik<-OUwie.fixed(phy=phy,data=data, model=model,simmap.tree=simmap.tree,root.station=root.station,alpha=alpha, sigma.sq=sigma.sq, theta=NULL, clade=clade)$loglik
  return(-loglik)
}

