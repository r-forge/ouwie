OUwie.dredge<-function(phy,data, criterion=c("aicc","aic","rjmcmc"), theta.max.k=3, sigma.max.k=3, alpha.max.k=3, root.station=TRUE, ip=1, lb=0.000001, ub=1000, max.generations=100, wait.generations=10, pop.size=NULL) {
  #criterion: If aicc or aic, use rgenoud, where the fitness is the score (and try to minimize)
  #    if rjmcmc, use an rjmcmc approach
  #alpha.max.k=3: allows for three OU regimes and one regime where alpha is set to ~0 (brownian motion)
  #    if alpha.max.k=0, only BM models are investigated
  #sigma.max.k, theta.max.k: work in same way, but must be at least 1 in each case
  #note that these params can be stuck together in interesting ways. For example, a BM jump model has multiple thetas but constant everything else
  #There are as many free parameters for rgenoud as the number of nodes (not edges) * 3: params for theta, sigma, and alpha
  #so each individual for rgenoud has theta1, theta2, theta3....,sigma1, sigma2, sigma3, ....,alpha1, alpha2, alpha3, where theta1 is the mapping to the theta free parameter for node 1, and so forth
  data<-data.frame(data[,2], data[,2], row.names=data[,1])
  data<-data[phy$tip.label,]
  if(criterion!="rjmcmc") {
   if (is.null(pop.size)) {
     pop.size<-5 #choose better 
    }
   starting.values<-matrix(c(rep(1,Nnode(phy,internal.only=FALSE)),rep(1,Nnode(phy,internal.only=FALSE)),rep(0,Nnode(phy,internal.only=FALSE))),nrow=1,ncol=3*Nnode(phy,internal.only=FALSE)) #BM1
   starting.values<-rbind(starting.values,matrix(c(rep(1,Nnode(phy,internal.only=FALSE)),rep(1,Nnode(phy,internal.only=FALSE)),rep(1,Nnode(phy,internal.only=FALSE))),nrow=1,ncol=3*Nnode(phy,internal.only=FALSE))) #OU1
   dput(starting.values)
   nodes<-Nnode(phy,internal.only=FALSE)
   Domains<-matrix(c(rep(1,nodes),rep(1,nodes),rep(0,nodes),rep(theta.max.k,nodes),rep(sigma.max.k,nodes),rep(alpha.max.k,nodes)),ncol=2,nrow=3*Nnode(phy,internal.only=FALSE),byrow=FALSE)
	print(dim(Domains))
 #  results<-genoud(fn=dredge.akaike, starting.values=starting.values, max=FALSE,nvars=3*Nnode(phy,internal.only=FALSE), data.type.int=TRUE,boundary.enforcement=2,Domains=matrix(c(rep(1,Nnode(phy,internal.only=FALSE)),rep(1,Nnode(phy,internal.only=FALSE)),rep(0,Nnode(phy,internal.only=FALSE)),rep(theta.max.k,Nnode(phy,internal.only=FALSE)),rep(sigma.max.k,Nnode(phy,internal.only=FALSE)),rep(alpha.max.k,Nnode(phy,internal.only=FALSE))),ncol=2,nrow=3*Nnode(phy,internal.only=FALSE),byrow=FALSE), wait.generations=wait.generations, max.generations=max.generations, pop.size=pop.size, root.station=TRUE, ip=ip, lb=lb, ub=ub, phy, data, criterion=criterion)
   results<-genoud(fn=dredge.akaike, max=FALSE,nvars=3*nodes, data.type.int=TRUE, boundary.enforcement=2, Domains=Domains, wait.generations=wait.generations, max.generations=max.generations, pop.size=pop.size, root.station=TRUE, ip=ip, lb=lb, ub=ub, phy=phy, data=data, criterion=criterion)
  }
}

dredge.akaike<-function(rgenoud.individual, phy, data, criterion="aicc",lb,ub,ip,root.station,...) {
  print("in dredge akaike")
  print(dim(data))
  #first check that the rgenoud.individual is well-structured: don't have just states 0 and 3 for sigma mapping, for example
  #if it fails this, reject it.
  #convert phy+regenoud.individual to simmap.tree (later, make it so that we directly go to the proper object)
  edge.mat.all<-edge.mat(phy,rgenoud.individual)
  #call dev.optimize
  lnL<-dev.optimize(edges.ouwie=edge.mat.all$edges.ouwie, regimes.mat=edge.mat.all$regime, data=data,root.station=root.station,lb=lb, ub=ub, ip=ip)
  #which is an nloptr wrapper to call dev.dredge
  #convert likelihood to AICC
  K<-sum(apply(edge.mat.all$regime, 2, max))
  result<-2*lnL + 2*K
  if (criterion=="aicc") {
    n=Ntip(phy)
   result<-2*lnL + 2*K * n / (n - K - 1)
  }
  #return AICC
  return(result)
}


dev.optimize<-function(edges.ouwie,regimes.mat,data,root.station,lb,ub,ip) {
print("in dev.optimize")
print(dim(data))
  np<-sum(apply(regimes.mat, 2, max))
  lower = rep(lb, np)
  upper = rep(ub, np)
  ip<-ip
  opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"="1000000", "ftol_rel"=.Machine$double.eps^0.5, "xtol_rel"=.Machine$double.eps^0.5)
  out = nloptr(x0=rep(ip, length.out = np), eval_f=dev.dredge, opts=opts, data=data, root.station=root.station, lb=lower, ub=upper, edges.ouwie=edges.ouwie, regimes.mat=regimes.mat)
  return(-1*out$objective) 
}


#edges.ouwie is the edges matrix with the extra columns for regimes
dev.dredge<-function(p,edges.ouwie,regimes.mat,data,root.station) { 
print("in dev.dredge")
print(dim(data))
  nRegimes<-dim(regimes.mat)[1]
  Rate.mat.full <- matrix(0, 3, nRegimes) #first row is alpha, second row is sigma, third is theta. Value of zero initially
  k.theta<-max(regimes.mat[,1])
  k.sigma<-max(regimes.mat[,2])
  k.alpha<-max(regimes.mat[,3])
  p<-c(1:length(p)) #CHANGE AFTER DEBUGGING
  p.index<-1
  for(k in sequence(k.theta)) {
    Rate.mat.full[1,which(regimes.mat[1,]==k)]<-p[p.index]
    p.index<-p.index+1
  }
  for(k in sequence(k.sigma)) {
    Rate.mat.full[2,which(regimes.mat[2,]==k)]<-p[p.index]
    p.index<-p.index+1
  }
  for(k in sequence(k.alpha)) {
    Rate.mat.full[3,which(regimes.mat[3,]==k)]<-p[p.index]
    p.index<-p.index+1
  }
  print(Rate.mat.full)
  Rate.mat<-Rate.mat.full[1:2,]
  print(dim(regimes.mat))
  root.state<-which(regimes.mat[dim(regimes.mat)[1] , ]==1)-4
  edges.ouwie.no.root<-edges.ouwie[1:(dim(edges.ouwie)[1]-1),]
  x<-as.matrix(data[,1])
  N<-length(x[,1])
  print(Rate.mat)
  V<-varcov.ou(phy, edges.ouwie.no.root, Rate.mat, root.state=root.state)
  W<-weight.mat(phy, edges.ouwie.no.root, Rate.mat, root.state=root.state, assume.station=root.station)
  theta<-Rate.mat.full[3,]
  
  DET<-determinant(V, logarithm=TRUE)
  
  res<-W%*%theta-x		
  q<-t(res)%*%solve(V,res)
  logl <- -.5*(N*log(2*pi)+as.numeric(DET$modulus)+q[1,1])
  
  return(-logl)
}

#Steps:
#Create possible regimes from rgenoud.individual -- done
#Write tree traversal to obtain possible regimes for each edge -- scores each in edge.mat
#Output is edges.ouwie
edge.mat<-function(phy,rgenoud.individual){
  
  obj=NULL
  #Values to be used throughout
  n=max(phy$edge[,1])
  ntips=length(phy$tip.label)
  
  ##Begins the construction of the edges matrix -- similar to the ouch format##
  #Makes a vector of absolute times in proportion of the total length of the tree
  branch.lengths=rep(0,(n-1))
  branch.lengths[(ntips+1):(n-1)]=branching.times(phy)[-1]/max(branching.times(phy))
  
  #New tree matrix to be used for subsetting regimes
  edges.ouwie=cbind(c(1:(n-1)),phy$edge,phy$edge.length)
  edges.ouwie=edges.ouwie[sort.list(edges.ouwie[,3]),]
  
  edges.ouwie[,4]=branch.lengths
  
  edges.ouwie=rbind(edges.ouwie,c(n,NA,Nnode(phy,internal.only=FALSE)+1,NA)) #now adding the root
  
  
  tot<-length(rgenoud.individual)/3
  tmp<-matrix(,tot,1)
  #Generates all strings in the rgenoud.individual:
  for(i in 1:tot){
    tmp[i,]<-paste(rgenoud.individual[i],rgenoud.individual[i+tot],rgenoud.individual[i+(2*tot)])
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
  obj$regime.mat<-matrix(as.numeric(unlist(strsplit(reg.list," "))),length(reg.list),3)
  
  regime <- matrix(0,nrow=length(tmp),ncol=length(unique(tmp)))
  #Generates an indicator matrix from the regime vector
  for (i in 1:length(tmp)) {
    regime[i,tmp[i]] <- 1 
  }
  #Finishes the edges matrix
  edges.ouwie=cbind(edges.ouwie,regime)
  #Resort the edge matrix so that it looks like the original matrix order:
  obj$edges.ouwie=edges.ouwie[sort.list(edges.ouwie[,1]),]
  
  obj
}