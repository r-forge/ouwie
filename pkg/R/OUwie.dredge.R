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
  #first check that the rgenoud.individual is well-structured: don't have just states 0 and 3 for sigma mapping, for example
  #if it fails this, reject it.
  #convert phy+regenoud.individual to simmap.tree (later, make it so that we directly go to the proper object)
  #call dev.optimize
  #which is an nloptr wrapper to call dev.dredge
  #convert likelihood to AICC
  #return AICC
}


dev.optimize<-(rgenoud.individual,edges.ouwie,regimes.mat,data,root.station) {
  
  
}


#edges.ouwie is the edges matrix with the extra columns for regimes
dev.dredge<-function(p,rgenoud.individual,edges.ouwie,regimes.mat,data,Rate.mat,index.mat,root.station) { 
  nRegimes<-dim(regimes.mat)[1]
  Rate.mat.full <- matrix(0, 3, nRegimes) #first row is alpha, second row is sigma, third is theta. Value of zero initially
  k.theta<-max(regimes.mat[,1])
  k.sigma<-max(regimes.mat[,2])
  k.alpha<-max(regimes.mat[,3])
  p.index<-1
  for(k in sequence(k.theta)) {
    Rate.mat.full[1,which(Rate.mat.full[1,]==k)]<-p[p.index]
    p.index<-p.index+1
  }
  for(k in sequence(k.sigma)) {
    Rate.mat.full[2,which(Rate.mat.full[2,]==k)]<-p[p.index]
    p.index<-p.index+1
  }
  for(k in sequence(k.alpha)) {
    Rate.mat.full[3,which(Rate.mat.full[3,]==k)]<-p[p.index]
    p.index<-p.index+1
  }
  Rate.mat<-Rate.mat.full[1:2,]
  N<-length(x[,1])
  V<-varcov.ou(phy, edges, Rate.mat, root.state=root.state)
  W<-weight.mat(phy, edges, Rate.mat, root.state=root.state, assume.station=root.station)
  
  
  
  
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
#  loglik<-OUwie.fixed(phy=phy,data=data, model=model,simmap.tree=simmap.tree,root.station=root.station,alpha=alpha, sigma.sq=sigma.sq, theta=NULL, clade=clade)$loglik
  return(-loglik)
}

#Steps:
#Create possible regimes from regime.vector -- done
#Write tree traversal to obtain possible regimes for each edge -- scores each in edges.mat
#Output is edge.ouwie
edge.mat<-function(phy,regime.vector){
  
  obj=NULL
  #Values to be used throughout
  n=max(phy$edge[,1])
  ntips=length(phy$tip.label)
  
  ##Begins the construction of the edges matrix -- similar to the ouch format##
  #Makes a vector of absolute times in proportion of the total length of the tree
  branch.lengths=rep(0,(n-1))
  branch.lengths[(ntips+1):(n-1)]=branching.times(phy)[-1]/max(branching.times(phy))
  
  #New tree matrix to be used for subsetting regimes
  edge.ouwie=cbind(c(1:(n-1)),phy$edge,phy$edge.length)
  edge.ouwie=edge.ouwie[sort.list(edge.ouwie[,3]),]
  
  edge.ouwie[,4]=branch.lengths
  
  tot<-length(regime.vector)/3
  tmp<-matrix(,tot,1)
  #Generates all strings in the regime.vector:
  for(i in 1:tot){
    tmp[i,]<-paste(regime.vector[i],regime.vector[i+tot],regime.vector[i+(2*tot)])
  }
  #Finds the unique combinations to designate regimes:
  reg.list<-unique(tmp)
  #Convert tmp to numeric regime designation:
  for(i in 1:tot){
    tmp[i,]<-which(tmp[i,]==reg.list)
  }
  #Convert to numeric:
  tmp<-as.numeric(tmp)
  #Creates the regime.mat for internal use:
  obj$regime.mat<-matrix(as.numeric(unlist(strsplit(reg.list," "))),3,3)
  
  regime <- matrix(0,nrow=length(tmp),ncol=length(unique(tmp)))
  #Generates an indicator matrix from the regime vector
  for (i in 1:length(tmp)) {
    regime[i,tmp[i]] <- 1 
  }
  #Finishes the edges matrix
  edge.ouwie=cbind(edge.ouwie,regime)
  #Resort the edge matrix so that it looks like the original matrix order:
  obj$edge.ouwie=edge.ouwie[sort.list(edge.ouwie[,1]),]
  
  obj
}