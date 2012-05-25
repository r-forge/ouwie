library(rgenoud)
library(RColorBrewer)

OUwie.dredge<-function(phy,data, criterion=c("aicc","aic","rjmcmc"), theta.max.k=3, sigma.max.k=3, alpha.max.k=3, root.station=TRUE, ip=1, lb=0.000001, ub=1000, max.generations=100, wait.generations=10, pop.size=NULL, print.level=0, maxeval=500, logfile=NULL) {
#criterion: If aicc or aic, use rgenoud, where the fitness is the score (and try to minimize)
#    if rjmcmc, use an rjmcmc approach
#alpha.max.k=3: allows for three OU regimes and one regime where alpha is set to ~0 (brownian motion)
#    if alpha.max.k=0, only BM models are investigated
#sigma.max.k, theta.max.k: work in same way, but must be at least 1 in each case
#note that these params can be stuck together in interesting ways. For example, a BM jump model has multiple thetas but constant everything else
#There are as many free parameters for rgenoud as the number of nodes (not edges) * 3: params for theta, sigma, and alpha
#so each individual for rgenoud has theta1, theta2, theta3....,sigma1, sigma2, sigma3, ....,alpha1, alpha2, alpha3, where theta1 is the mapping to the theta free parameter for node 1, and so forth

	#Coerce the data so that it will run in OUwie -- using values of OU1 as the starting points: 
	cat("Initializing...","\n")
	
	data2<-data.frame(as.character(data[,1]),sample(c(1:2),length(data[,1]), replace=T),data[,2],stringsAsFactors=FALSE)
	phy$node.label<-sample(c(1:2),phy$Nnode, replace=T)
	start<-OUwie(phy,data2,model=c("OU1"),plot.resid=FALSE, quiet=TRUE)
	ip<-matrix(c(rep(start$theta[,1],Nnode(phy,internal.only=FALSE)),rep(start$Param.est[2],Nnode(phy,internal.only=FALSE)),rep(start$Param.est[1],Nnode(phy,internal.only=FALSE))),nrow=1,ncol=3*Nnode(phy,internal.only=FALSE)) #OU1

	data<-data.frame(data[,2], data[,2], row.names=data[,1])
	data<-data[phy$tip.label,]
	
	if(criterion!="rjmcmc") {
		if (is.null(pop.size)) {
			pop.size<-10 #choose better 
		}
		
		cat("Begin fast optimization routine -- Starting values:", c(start$theta[,1],start$Param.est[2],start$Param.est[1]), "\n")
		nodes<-Nnode(phy,internal.only=FALSE)
		starting.individuals<-matrix(c(rep(0,Nnode(phy,internal.only=FALSE)),rep(0,Nnode(phy,internal.only=FALSE)),rep(0,Nnode(phy,internal.only=FALSE))),nrow=1,ncol=3*Nnode(phy,internal.only=FALSE)) #BM1
		starting.individuals[1,c(nodes,2*nodes,3*nodes)]<-c(1,1,-1)
		starting.individuals<-rbind(starting.individuals,matrix(c(rep(0,Nnode(phy,internal.only=FALSE)),rep(0,Nnode(phy,internal.only=FALSE)),rep(0,Nnode(phy,internal.only=FALSE))),nrow=1,ncol=3*Nnode(phy,internal.only=FALSE))) #OU1
		starting.individuals[2,c(nodes,2*nodes,3*nodes)]<-c(1,1,1)
		Domains<-matrix(c(rep(0,nodes-1),1,rep(0,nodes-1),1,rep(-1,nodes),rep(theta.max.k,nodes),rep(sigma.max.k,nodes),rep(alpha.max.k,nodes)),ncol=2,nrow=3*Nnode(phy,internal.only=FALSE),byrow=FALSE)
		if(!is.null(logfile)) {
			write.table(t(c("aicc","k.theta","k.sigma.sq","k.alpha")),file=logfile,quote=F,sep="\t",row.name=F,col.name=F)
		}
		results<-genoud(fn=dredge.akaike, starting.values=starting.individuals, max=FALSE,nvars=3*nodes, data.type.int=TRUE, logfile=logfile,print.level=print.level, boundary.enforcement=2, Domains=Domains, wait.generations=wait.generations, maxeval=maxeval, max.generations=max.generations, pop.size=pop.size, root.station=TRUE, ip=ip, lb=lb, ub=ub, phy=phy, data=data, criterion=criterion)
		cat("Finished. Begin thorough optimization routine", "\n")

		edge.mat<-edge.mat(phy,results$par)
		np<-sum(apply(edge.mat$regime.mat, 2, max))
		K<-sum(apply(edge.mat$regime, 2, max))
		n=Ntip(phy)
		
		lower = rep(lb, np)
		upper = rep(ub, np)
		
		opts <- list("algorithm"="NLOPT_LN_SBPLX", "maxeval"=as.character(maxeval), "ftol_rel"=.Machine$double.eps^0.5)
		out = nloptr(x0=rep(ip, length.out = np), eval_f=dev.dredge, opts=opts, data=data, phy=phy,root.station=root.station, lb=lower, ub=upper, edges.ouwie=edge.mat$edges.ouwie, regime.mat=edge.mat$regime.mat)
		loglik = (-1) * out$objective #since dev.dredge() returns negloglik
		obj = list(loglik = loglik, AIC = -2*loglik+2*K,AICc=-2*loglik + (2*K * n / (n - K - 1)),solution=out$solution, opts=opts, data=data, phy=phy,root.station=root.station, lb=lower, ub=upper, edges.ouwie=edge.mat$edges.ouwie, regime.mat=edge.mat$regime.mat,start=start,rgenoud.individual=results$par) 
		class(obj)<-"ouwie.dredge.result"		

		return(obj)
	}
}

print.ouwie.dredge.result <- function(x, ...) {
	K<-sum(unlist(param.count(x$rgenoud.individual,x$phy)))
	n<-Ntip(x$phy)
	output<-c(x$loglik,x$AIC,x$AICc,n,dim(x$regime.mat)[1],K,c(unlist(param.count(x$rgenoud.individual,x$phy))))
	names(output)<-c("negLnL","AIC","AICc","ntax","n_regimes","K_all","K_theta","K_sigma","K_alpha")
	print(output)
	
	cat("\nRegime matrix\n")
	colnames(x$regime.mat)<-c("theta","sigma","alpha")
	rownames(x$regime.mat)<-c(paste("regime_",sequence(dim(x$regime.mat)[1]),sep=""))
	print(x$regime.mat)
	
	regime.mat.params<-x$regime.mat*0
	p.index<-1
	for(k in sequence(param.count(x$rgenoud.individual,x$phy)[[1]])) {
		regime.mat.params[which(x$regime.mat[,1]==k),1]<-x$solution[p.index]
		p.index<-p.index+1
	}
	for(k in sequence(param.count(x$rgenoud.individual,x$phy)[[2]])) {
		regime.mat.params[which(x$regime.mat[,2]==k),2]<-x$solution[p.index]
		p.index<-p.index+1
	}
	for(k in sequence(param.count(x$rgenoud.individual,x$phy)[[3]])) {
		regime.mat.params[which(x$regime.mat[,3]==k),3]<-x$solution[p.index]
		p.index<-p.index+1
	}
	cat("\nRate matrix\n")
	print(regime.mat.params)
	cat("\n")
	
}

##   Regime vs. Rate        ##
##   gradient of color      ##

plot.ouwie.dredge.result <- function(x, type=c("regime", "rate"), col.pal=c("Set1"), ...) {
	x$rgenoud.individual<-as.full.regime(x$rgenoud.individual,x$phy)
	par(mfcol=c(1,3))
	tot<-length(x$rgenoud.individual)/3
		
	if(type=="regime"){
		
		regimes<-as.factor(x$rgenoud.individual)
		regime.lvls<-length(levels(regimes))
		if(regime.lvls<3){
			regime.lvls = 3
		}
		#
		if(length(col.pal>1)){
			co<-brewer.pal(regime.lvls, col.pal)
		}
		#If not the 
		else{
			co<-col.pal
		}
		##Plot thetas
		nb.tip <- Ntip(x$phy)
		nb.node <- Nnode(x$phy)
		comp <- numeric(Nedge(x$phy))
		comp[match(1:(Ntip(x$phy)), x$phy$edge[,2])] <- as.factor(regimes[1:nb.tip])
		comp[match((2+Ntip(x$phy)):(Nedge(x$phy)+1), x$phy$edge[,2])] <- as.factor(regimes[(nb.tip+1):(nb.tip+nb.node)][-1])
		plot.phylo(x$phy,edge.color=co[comp], ...)
		title(main="theta")
		
		##Plot sigmas
		nb.tip <- Ntip(x$phy)
		nb.node <- Nnode(x$phy)
		comp <- numeric(Nedge(x$phy))
		comp[match(1:(Ntip(x$phy)), x$phy$edge[,2])] <- as.factor(regimes[(tot+1):(tot+nb.tip)])
		comp[match((2+Ntip(x$phy)):(Nedge(x$phy)+1), x$phy$edge[,2])] <- as.factor(regimes[(2*tot+1):(2*(tot))][-1])
		plot.phylo(x$phy,edge.color=co[comp], ...)
		title(main="sigma")
		
		##Plot alphas
		nb.tip <- Ntip(x$phy)
		nb.node <- Nnode(x$phy)
		comp <- numeric(Nedge(x$phy))
		comp[match(1:(Ntip(x$phy)), x$phy$edge[,2])] <- as.factor(regimes[(2*tot+1):(2*tot+nb.tip)])
		comp[match((2+Ntip(x$phy)):(Nedge(x$phy)+1), x$phy$edge[,2])] <- as.factor(regimes[(3*tot+1):length(x$rgenoud.individual)][-1])
		plot.phylo(x$phy,edge.color=co[comp], ...)	
		title(main="alpha")
	}
}

param.count<-function(rgenoud.individual,phy) {
	rgenoud.individual<-as.full.regime(rgenoud.individual,phy,alphaZero=TRUE)	
	return(list(k.theta=max(rgenoud.individual[1:(length(rgenoud.individual)/3)]),k.sigma=max(rgenoud.individual[(1+length(rgenoud.individual)/3):(2*length(rgenoud.individual)/3)]),k.alpha=max(rgenoud.individual[(1+2*length(rgenoud.individual)/3):length(rgenoud.individual)])))
}

valid.individual<-function(rgenoud.individual,phy) {
	rgenoud.individual<-as.full.regime(rgenoud.individual,phy,alphaZero=TRUE)
	#want to make sure we do not have rates 0 and 2 rather than 0 and 1
	theta.vec<-rgenoud.individual[1:(length(rgenoud.individual)/3)]
	sigma.vec<-rgenoud.individual[(1+length(rgenoud.individual)/3):(2*length(rgenoud.individual)/3)]
	alpha.vec<-rgenoud.individual[(1+2*length(rgenoud.individual)/3):length(rgenoud.individual)]

	if(length(unique(theta.vec)) != (1 + max(theta.vec) - min(theta.vec))) {
		return(FALSE)
	}
	if (min(theta.vec)!=1) {
		return(FALSE)
	}
	if(length(unique(sigma.vec)) != (1 + max(sigma.vec) - min(sigma.vec))) {
		return(FALSE)
	}
	if (min(sigma.vec)!=1) {
		return(FALSE)
	}
	if(length(unique(alpha.vec)) != (1 + max(alpha.vec) - min(alpha.vec))) {
		return(FALSE)
	}
	if(min(alpha.vec)>1) {
		return(FALSE)
	}
	if(min(alpha.vec)<0) {
		return(FALSE)
	}
	else {
		return(TRUE)
	}
}

dredge.akaike<-function(rgenoud.individual, phy, data, criterion="aicc",lb,ub,ip,root.station,maxeval,logfile=NULL,badvalue=100000000,...) {
	#first check that the rgenoud.individual is well-structured: do not have just states 0 and 3 for sigma mapping, for example
	if (!valid.individual(rgenoud.individual,phy)) {
		return(badvalue)
	}
	#if it fails this, reject it.
	#convert phy+regenoud.individual to simmap.tree (later, make it so that we directly go to the proper object)
	edge.mat.all<-edge.mat(phy,rgenoud.individual)
	#call dev.optimize
	fit<-dev.optimize(edges.ouwie=edge.mat.all$edges.ouwie, regime.mat=edge.mat.all$regime, data=data,maxeval=maxeval, root.station=root.station,lb=lb, ub=ub, ip=ip, phy=phy) 
	#which is an nloptr wrapper to call dev.dredge
	#convert likelihood to AICC
	K<-sum(apply(edge.mat.all$regime, 2, max))
	result<-(-2*fit$loglik) + (2*K)
	if (criterion=="aicc") {
		n=Ntip(phy)
		result<-(-2)*fit$loglik + (2*K * n / (n - K - 1))
	}
	tmp<-c(result,unlist(param.count(rgenoud.individual,phy)),fit$pars)
	names(tmp)<-NULL
	if(!is.null(logfile)) {
		write.table(t(tmp),file=logfile,quote=F,sep="\t",row.name=F,col.name=F,append=T)
	}
	return(result)
}


dev.optimize<-function(edges.ouwie,regime.mat,data,root.station,maxeval,lb,ub,ip,phy=phy) {
	obj<-NULL
	
	np<-sum(apply(regime.mat, 2, max))
	lower = rep(lb, np)
	upper = rep(ub, np)
	ip<-ip
	opts <- list("algorithm"="NLOPT_LN_BOBYQA", "maxeval"=as.character(maxeval), "ftol_rel"=0.01)
	
	out = nloptr(x0=rep(ip, length.out = np), eval_f=dev.dredge, opts=opts, data=data, phy=phy,root.station=root.station, lb=lower, ub=upper, edges.ouwie=edges.ouwie, regime.mat=regime.mat) 
	obj$loglik<-(-1)*out$objective 
	obj$pars<-out$solution
	obj
}


#edges.ouwie is the edges matrix with the extra columns for regimes
dev.dredge<-function(p,edges.ouwie,regime.mat,data,root.station,phy) { 
	nRegimes<-dim(regime.mat)[1]
	Rate.mat.full <- matrix(0, 3, nRegimes) #first row is alpha, second row is sigma, third is theta. Value of zero initially
	k.alpha<-max(regime.mat[,3])
	k.sigma<-max(regime.mat[,2])
	k.theta<-max(regime.mat[,1])
	p.index<-1
	for(k in sequence(k.theta)) {
		Rate.mat.full[1,which(regime.mat[,1]==k)]<-p[p.index]
		p.index<-p.index+1
	}
	for(k in sequence(k.sigma)) {
		Rate.mat.full[2,which(regime.mat[,2]==k)]<-p[p.index]
		p.index<-p.index+1
	}
	for(k in sequence(k.alpha)) {
		Rate.mat.full[3,which(regime.mat[,3]==k)]<-p[p.index]
		p.index<-p.index+1
	}
	Rate.mat<-matrix(Rate.mat.full[1:2,],nrow=2,ncol=nRegimes) #deals with case of single column Rate.mat.full being changed into a vector rather than left as a matrix
	root.state<-NULL
	edges.ouwie.no.root<-edges.ouwie[1:(dim(edges.ouwie)[1]-1),]
	x<-as.matrix(data[,1])
	N<-length(x[,1])
	
	V<-varcov.ou(phy, edges.ouwie.no.root, Rate.mat, root.state=root.state)
	W<-weight.mat(phy, edges.ouwie.no.root, Rate.mat, root.state=root.state, assume.station=root.station)
	theta<-Rate.mat.full[3,]
	
	DET<-determinant(V, logarithm=TRUE)

	logl<--.5*(t(W%*%theta-x)%*%pseudoinverse(V)%*%(W%*%theta-x))-.5*as.numeric(DET$modulus)-.5*(N*log(2*pi)) #logL, not neg logL
	neglogl<-(-1)*logl
	return(as.numeric(neglogl)) #now neg log L.
}

#this creates a vector. The third element in this vector corresponds to the number appearing in the 
#   second column (descendant) of the edge matrix in the phy object for that edge. 
#   So if you want to know the ancestor of the i-th element in the rgenoud.individual,
#   the focal node number is phy$edge[(get.mapping(phy,rgenoud.individual))[i],2]
#   and the parent node number is phy$edge[(get.mapping(phy,rgenoud.individual))[i],1]
#   the state of the parent node in the rgenoud.individual is then
#   rgenoud.individual[ which(get.mapping(phy,rgenoud.individual) == phy$edge[(get.mapping(phy,rgenoud.individual))[i],1]) ]
get.mapping<-function(rgenoud.individual,phy) {
	mapping<-match(phy$edge[,2],1:length(rgenoud.individual))
	mapping<-append(mapping,which(!(sequence(max(mapping)) %in% mapping))) #add on the root taxon
	return(mapping) 
}


get.final.label<-function(i,rgenoud.individual,phy) {
	if (rgenoud.individual[i]!=0) {
		return(rgenoud.individual[i])
	}
	indexOffset<-(length(rgenoud.individual)/3) * floor(i/(length(rgenoud.individual)/3))
	nodeCount <- Nnode(phy,internal.only=FALSE)
	mapping <- get.mapping(rgenoud.individual,phy)
	while (rgenoud.individual[i] == 0 ) {
		indexMapping <- i %% nodeCount
		current.node <- mapping[ indexMapping ]
		parent.node <- phy$edge[ which(phy$edge[,2] == current.node), 1]
		i <- which( mapping == parent.node ) + indexOffset
		if (length(i) ==0 ) {
			return ( 1) #can't find the parent node because we ARE at the root
		}
	}
	return( rgenoud.individual[i] )
}



as.full.regime<-function(rgenoud.individual,phy,alphaZero=FALSE) {
	rgenoud.individual<-sapply(X=sequence(length(rgenoud.individual)),FUN=get.final.label,rgenoud.individual=rgenoud.individual,phy=phy)
	if(alphaZero) {
		rgenoud.individual[which(rgenoud.individual==(-1))] <-0
	}
	return(rgenoud.individual)
}

#Steps:
#Create possible regimes from rgenoud.individual -- done
#Write tree traversal to obtain possible regimes for each edge -- scores each in edge.mat
#Output is edges.ouwie
edge.mat<-function(phy,rgenoud.individual){ #requires full mapping: no 0 values to match to parent node
	rgenoud.individual<-as.full.regime(rgenoud.individual,phy)
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
	obj$regime.mat<-matrix(as.numeric(unlist(strsplit(reg.list," "))),length(reg.list),3,byrow=TRUE)
	
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

