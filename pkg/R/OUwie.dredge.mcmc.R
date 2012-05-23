#do mcmc on the discrete parameters of the model; use optimization on the other parameters. Saves us from having to do the rjmcmc, since all the parameters are always there.



#allow users to pass priors as vectors of the probability of 0, 1, 2, etc. up to the maximum number of, say, k.alpha
OUwie.dredge.mcmc<-function(phy,data, prior.k.theta, prior.k.sigma, prior.k.alpha, ngen=10000000, sample.freq=100, print.freq=10, root.station=TRUE, ip=1, lb=0.000001, ub=1000,maxeval=500, samplesfile="mcmc.samples.txt") {
	data2<-data.frame(as.character(data[,1]),sample(c(1:2),length(data[,1]), replace=T),data[,2],stringsAsFactors=FALSE)
	phy$node.label<-sample(c(1:2),phy$Nnode, replace=T)
	start<-OUwie(phy,data2,model=c("OU1"),plot.resid=FALSE, quiet=TRUE)
	ip<-matrix(c(rep(start$theta[,1],Nnode(phy,internal.only=FALSE)),rep(start$Param.est[2],Nnode(phy,internal.only=FALSE)),rep(start$Param.est[1],Nnode(phy,internal.only=FALSE))),nrow=1,ncol=3*Nnode(phy,internal.only=FALSE)) #BM1

	data<-data.frame(data[,2], data[,2], row.names=data[,1])
	data<-data[phy$tip.label,]
	bayes.individual<-matrix(c(rep(1,Nnode(phy,internal.only=FALSE)),rep(1,Nnode(phy,internal.only=FALSE)),rep(0,Nnode(phy,internal.only=FALSE))),nrow=1,ncol=3*Nnode(phy,internal.only=FALSE)) #BM1
	current.state<-measure.proposal(phy=phy, data=data, new.individual=bayes.individual, prior.k.theta=prior.k.theta, prior.k.sigma=prior.k.sigma, prior.k.alpha=prior.k.alpha, root.station=root.station, lb=lb, ub=ub, ip=ip, maxeval=maxeval)
	store.state(current.state,samplesfile,generation=0)
	labels<-c("generation","posterior","loglik", "prior", "k.theta", "k.sigma", "k.alpha", "nRegimes")
	cat(labels)
	write.table(matrix(labels,nrow=1),file=samplesfile,quote=F,sep="\t",row.name=F,col.name=F,append=F)
	print.state(current.state,samplesfile,generation=0)

	for(generation in sequence(ngen)) {
		next.state<-measure.proposal(phy, data, new.individual=transform.single(current.state$new.individual), prior.k.theta,prior.k.sigma, prior.k.alpha, maxeval, root.station, lb, ub, ip)
		if ((next.state$posterior - current.state$posterior) > runif(1,0,1)) {
			current.state<-next.state
		}
		if (generation %% sample.freq == 0) {
			store.state(current.state,samplesfile,generation)
		}
		if (generation %% print.freq == 0) {
			print.state(current.state,samplesfile,generation)
		}
	}
}

truncated.geometric.prior<-function(min=1,max=5,prob=0.5) {
	if(min<0) {
		stop("min in truncated.geometric.prior shoudl be non-negative")
	}
	params<-c(min:max)
	probs<-sapply(params,dgeom,prob)
	probs<-probs/sum(probs) #remember this could be from min=2 to max=4. We want to start at zero
	allprobs<-rep(0,max+1) #first element is zero
	allprobs[params+1]<-probs #deals with the zero offset
	return(allprobs)
}

truncated.uniform.prior<-function(min=1,max=5) {
	if(min<0) {
		stop("min in truncated.uniform.prior shoudl be non-negative")
	}
	params<-c(min:max)
	probs<-rep(1,max-min+1)
	probs<-probs/sum(probs)
	allprobs<-rep(0,max+1) #first element is zero
	allprobs[params+1]<-probs #deals with the zero offset
	return(allprobs)	
}

#transforms just a single cell
transform.single<-function(bayes.individual) {
	valid.transform<-FALSE
	while(!valid.transform) { #cool thing about doing this is that the hastings ratio is one: prob of going from valid to invalid is zero, as is the reverse rate
		new.individual<-bayes.individual
		focal.param<-sample.int(length(new.individual),1)
		new.individual[focal.param]<-new.individual[focal.param]+sign(rnorm(1))
		valid.transform<-valid.individual(new.individual)
	}
	return(new.individual)
}

prior.prob<-function(bayes.individual,prior.k.theta,prior.k.sigma, prior.k.alpha) {
	k.list<-dredge.util(bayes.individual)
	return(vector.index.safe.offset(k.list$k.theta,prior.k.theta) * vector.index.safe.offset(k.list$k.sigma,prior.k.sigma) * vector.index.safe.offset(k.list$k.alpha,prior.k.alpha) )
}

vector.index.safe<-function(x,vec) {
	if (x>length(vec)) {
		return(0)
	}
	else {
		return(vec[x])
	}
}

vector.index.safe.offset<-function(x,vec) {
	x<-x+1 #deals with the indexing by zero but R indexing by 1
	if (x>length(vec)) {
		return(0)
	}
	else {
		return(vec[x])
	}
}

measure.proposal<-function(phy, data, new.individual, prior.k.theta,prior.k.sigma, prior.k.alpha, maxeval, root.station, lb, ub, ip) {
	edge.mat.all<-edge.mat(phy,new.individual)
	loglik<-(-1)*dev.optimize(edges.ouwie=edge.mat.all$edges.ouwie, regime.mat=edge.mat.all$regime, data=data,maxeval=maxeval, root.station=root.station,lb=lb, ub=ub, ip=ip, phy=phy)$loglik
	prior<-prior.prob(new.individual,prior.k.theta,prior.k.sigma, prior.k.alpha)
	posterior<-loglik+log(prior)
	k.list<-dredge.util(new.individual)
	return(list(new.individual=new.individual, loglik=loglik, prior=prior, posterior=posterior, k.theta=k.list$k.theta, k.sigma=k.list$k.sigma, k.alpha=k.list$k.alpha,nRegimes=dim(edge.mat.all$regime)[1]))
}

store.state<-function(current.state,samplesfile,generation) {
	tmp<-c(generation,current.state$posterior,current.state$loglik, current.state$prior, current.state$k.theta, current.state$k.sigma, current.state$k.alpha, current.state$nRegimes, current.state$new.individual)
	write.table(matrix(tmp,nrow=1),file=samplesfile,quote=F,sep="\t",row.name=F,col.name=F,append=T)
}

print.state<-function(current.state,samplesfile,generation) {
	cat("\n",c(generation,current.state$posterior,current.state$loglik, current.state$prior, current.state$k.theta, current.state$k.sigma, current.state$k.alpha, current.state$nRegimes))
}