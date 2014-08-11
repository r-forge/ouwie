#Generates VCV for a network
#Uses code modified from BMhyd

#note that for edges matrix, col 4 is height above root of rootward end of edge, col 5 is height of tipward end of edge, 2 is ape edge[,1], 3 is ape edge[,2]
varcov.network<-function(phy, edges, Rate.mat, root.state, simmap.tree=FALSE, scaleHeight=FALSE, flow, scaling.by.sd=1, force.hybrid.alpha=FALSE, force.hybrid.sigma=FALSE){
	if(simmap.tree) {
		stop("simmap character mapping not supported for network topologies yet")
	}
	#last parameters of Rate.mat, in order, are bt and vh
	OU.Rate.mat<-Rate.mat[-(c(-1,0)+length(Rate.mat))] #assumptions here are that Rate.mat is a vector; we're removing the last two elements
	stop("the above assumption is untested")
	vcv <- varcov.ou(phy, edges, OU.Rate.mat, root.state, simmap.tree, scaleHeight) #
	bt <- Rate.mat[(length(Rate.mat)-1)] #penultimate parameter is beta
	vh <- Rate.mat[length(Rate.mat)]
	v.modified <- vcv
	alpha=Rate.mat[1,]
	sigma=Rate.mat[2,]
	if(is.null(root.state)) {
		root.state<-which(edges[dim(edges)[1],]==1)-5
		edges<-edges[-1*dim(edges)[1],]
	}
	pp <- prop.part(phy)
	oldregime=root.state
	nodevar1=rep(0,max(edges[,3]))
	nodevar2=rep(0,max(edges[,3]))
	n.cov1=matrix(rep(0,n), n, 1)
	n.cov2=matrix(rep(0,n), n, 1)
	mm<-dim(edges)
	k<-length(6:mm[2])

	for(flow.index in sequence(dim(flow)[1])) {
		recipient.taxa <- strsplit(flow$recipient, ",")[[1]]
		donor.taxa <- strsplit(flow.donor, ",")[[1]]
		
		recipient.indices <- which(rownames(vcv(phy)) %in% recipient.taxa)
		donor.indices <- which(rownames(vcv(phy)) %in% donor.taxa)

		if(length(recipient.indices)!=length(recipient.taxa)) {
			stop(paste("Tried to find ", flow$recipient[flow.index], " but failed; make sure the taxon names in the flow dataframe recipient match that of your tree", sep=""))
		}
		donor.index <- which(rownames(v.modified)==flow$donor[flow.index])
		if(length(donor.indices)!=length(donor.taxa)) {
			stop(paste("Tried to find ", flow$donor[flow.index], " but failed; make sure the taxon names in the flow dataframe recipient match that of your tree", sep=""))
		}
		recipient.crown.node <- getMRCA(phy, tips=recipient.taxa)
		if(is.null(recipient.crown.node) { #must be a single recipient taxon
			recipient.crown.node <- which(phy$tip.label == recipient.taxa[1])
		}	
		recipient.stem.node <- edges[which(edges[,3]==recipient.crown.node,2] #the parent with the main gene flow, not the one with flow rate of m
		parent.1 <- recipient.stem.node
		donor.crown.node <- getMRCA(phy, tips=donor.taxa)
		if(is.null(donor.crown.node) { #must be a single donor taxon
			donor.crown.node <- which(phy$tip.label == donor.taxa[1])
		}
		parent.1.regime <- which(edges[parent.1,6:(k+5)]==1)
		parent.2.tipward = donor.crown.node
		parent.2.rootward = edges[donor.crown.node, 2]
		parent.2.tipward.regime <- which(edges[donor.crown.node,6:(k+5)]==1)
		parent.2.rootward.regime <- which(edges[which(edges[,3]==donor.crown.node),6:(k+5)]==1)
		parent.2.regime <- c(parent.2.rootward.regime, parent.2.tipward.regime)[which.min(abs(flow$time.from.root.recipient[flow.index] - edges(donor.crown.node,c(2:3)))))] #finds which end is closer to the time of the hybridization event, uses regime from that end
		#parent.2 is somewhere between donor.crown.node and the node below (that is, it occurs on this edge) 
		
		alpha1 <- alpha[parent.1.regime]
		alpha2 <- alpha[parent.2.regime]
		sigma1 <- sigma[parent.1.regime]
		sigma2 <- sigma[parent.2.regime]
		weights<-GetParentWeights(alpha1, alpha2, sigma1, sigma2, flow$m[flow.index], scaling.by.sd=scaling.by.sd, time.from.root.recipient=flow$time.from.root.recipient[flow.index])

		
		for (hybrid.local.index in sequence(recipient.indices)) {
			hybrid.index <- recipient.indices[hybrid.local.index]
			hybrid.regime <- which(edges[which(edges[3]==hybrid.index),6:(k+5)]==1)
			alpha.hybrid <- alpha[hybrid.regime]
			if(force.hybrid.alpha) {
				alpha.hybrid<-sum(c(alpha1, alpha2) * weights)
			}
			sigma.hybrid <- sigma[hybrid.regime]
			if(force.hybrid.sigma) {
				sigma.hybrid<-sum(c(sigma1, sigma2) * weights)
			}
			
			for (donor.local.index in sequence(donor.indices)) {
				donor.index <- donor.indices[donor.local.index]
				v.modified[hybrid.index, donor.index] <- stop("something")
				v.modified[donor.index, hybrid.index] <- v.modified[hybrid.index, donor.index]
			}
			stop("also need to adjust covariance between hybrid and parent 1")
			#so idea here is to do weighted average of the variance for both paths of the hybrid, through parent 1 and parent 2. The first path is straightfoward. 
			#The second one: we go from hybrid down to where it currently attaches, and stop. 
			#Then we go from where its other parent is down to the root, but starting partway down that first edge
			#Add up these variances, plus add on a weighted V_h
			v.modified[hybrid.index, hybrid.index] <- stop("the scaling thing here, like at the end of varcov.ou") * (weights[1] * AddVarianceFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = hybrid.index, tipward.height = edges[hybrid.index,5], rootward.height=0, Rate.mat) + weights[2] * ( 
				AddVarianceFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = hybrid.index, rootward.node=parent.1, tipward.height = edges[hybrid.index,5], rootward.height=flow$time.from.root.recipient[flow.index], Rate.mat) + 
				AddVarianceFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = parent.2.tipward, tipward.height = flow$time.from.root.recipient[flow.index], rootward.height=0, Rate.mat)
			)) + vh*exp(-2*alpha.hybrid*flow$time.from.root.recipient[flow.index])

		}
		
		
	}
}

#by default, go to root node, assumed to be at Ntip(phy)+1, according to Jeremy Beaulieu
AddVarianceFromPartwayAlongEdgeDownToAnotherPoint <- function(phy, edges, tipward.node, rootward.node=Ntip(phy)+1, tipward.height, rootward.height=0, Rate.mat) {
	alpha=Rate.mat[1,]
	sigma=Rate.mat[2,]
	done <- FALSE
	variance.sum <- 0
	anc = edges[tipward.node,2]
	oldtime=edges[tipward.node,4]
	newtime=edges[tipward.node,5]
	if(anc%in%edges[,3]){
		start=which(edges[,3]==anc)
		oldregime=which(edges[start,6:(k+5)]==1)
	}
	else{
		#For the root:
		oldregime=root.state
	}	
	newregime=which(edges[tipward.node,6:(k+5)]==1)
	oldtime <- max(oldtime, rootward.height)
	if(oldregime==newregime){
		variance.sum <- variance.sum + sigma[oldregime]*((exp(2*alpha[oldregime]*tipward.height)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))
	} 	else{
		a.time <- tipward.height - mean(c(newtime, oldtime))
		halftime=newtime-((newtime-oldtime)/2)
		if(a.time<0) { #we're starting out rootward of the half time
			variance.sum <- variance.sum + sigma[oldregime]*((exp(2*alpha[oldregime]*tipward.height)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))
		} else {
			variance.sum <- variance.sum + sigma[oldregime]*((exp(2*alpha[oldregime]*halftime)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))			
			variance.sum <- variance.sum + sigma[newregime]*((exp(2*alpha[newregime]*tipward.height)-exp(2*alpha[newregime]*tipward.height))/(2*alpha[newregime]))			
		}
	}
	current.node <- tipward.node
	if(edges[current.node,2]==rootward.node) {
		done<-TRUE
	}
	while(!done) {
		current.node <- edges[current.node,2] #move rootward
		if(current.node == rootward.node || current.node == (Ntip(phy)+1)) { #if we have gone down to the root
			done <- TRUE #so this'll be our last run here
		}
		anc = edges[current.node,2]
		oldtime=edges[current.node,4]
		newtime=edges[current.node,5]
		if(anc%in%edges[,3]){
			start=which(edges[,3]==anc)
			oldregime=which(edges[start,6:(k+5)]==1)
		}
		else{
			#For the root:
			oldregime=root.state
		}	
		newregime=which(edges[current.node,6:(k+5)]==1)
		oldtime <- max(oldtime, rootward.height)
		if(oldregime==newregime){
			variance.sum <- variance.sum + sigma[oldregime]*((exp(2*alpha[oldregime]*tipward.height)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))
		} 	else{
			a.time <- tipward.height - mean(c(newtime, oldtime))
			halftime=newtime-((newtime-oldtime)/2)
			if(a.time<0) { #we're starting out rootward of the half time
				variance.sum <- variance.sum + sigma[oldregime]*((exp(2*alpha[oldregime]*tipward.height)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))
			} else {
				variance.sum <- variance.sum + sigma[oldregime]*((exp(2*alpha[oldregime]*halftime)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))			
				variance.sum <- variance.sum + sigma[newregime]*((exp(2*alpha[newregime]*tipward.height)-exp(2*alpha[newregime]*tipward.height))/(2*alpha[newregime]))			
			}
		}
		
	}
	return(variance.sum)
}

GetOUVariance <- function(sigma, alpha, time.from.root) {
	return( (sigma^2) * (1 - exp(-2 * alpha *time.from.root)) / (2 * alpha))	
}


GetParentWeights <- function(alpha1, alpha2, sigma1, sigma2, m, scaling.by.sd, time.from.root.recipient) { #assumes m is parameter for proportion flow from parent 2: that is, the main history is from parent 1, and m (0 - 0.5 range, generally) comes from parent 2
	sd1 <- sqrt(ou.var(sigma1, alpha1, time.from.root.recipient))
	sd2 <- sqrt(ou.var(sigma2, alpha2, time.from.root.recipient))
	sd1.wt <- scaling.by.sd * ((1-m)/sd1)/ (((m)/sd2+(1-m)/sd1)) + (1-scaling.by.sd)*(1-m)
	sd2.wt <- scaling.by.sd * ((m)/sd2)/ (((m)/sd2+(1-m)/sd1)) + (1-scaling.by.sd)*(m)
	return(c(sd1.wt, sd2.wt))
}