#Generates VCV for a network
#Uses code modified from BMhyd


EdgesRowToTipwardNumber <- function(x, edges) {
	return(edges[x,3])
}

EdgesRowToRootwardNumber <- function(x, edges) {
	return(edges[x,2])
}

EdgesTipwardNumberToRow <- function(x, edges) {
	return(which(edges[,3]==x))
}

EdgesRootwardNumberToRow <- function(x, edges) {
	return(which(edges[,2]==x))
}

EdgesGetOneNodeNumberRootward <- function(x, edges) {
	return(edges[EdgesTipwardNumberToRow(x, edges), 2])
}




#note that for edges matrix, col 4 is height above root of rootward end of edge, col 5 is height of tipward end of edge, 2 is ape edge[,1], 3 is ape edge[,2]
varcov.network<-function(phy, edges, Rate.mat, root.state, simmap.tree=FALSE, scaleHeight=FALSE, flow, scaling.by.sd=1, force.hybrid.alpha=FALSE, force.hybrid.sigma=FALSE, vh=0, bt=1){
	if(simmap.tree) {
		stop("simmap character mapping not supported for network topologies yet")
	}
	vcv <- varcov.ou(phy, edges, Rate.mat, root.state, simmap.tree, scaleHeight) #
	v.modified <- vcv
	alpha=Rate.mat[1,]
	sigma.sq=Rate.mat[2,]
	if(is.null(root.state)) {
		root.state<-which(edges[dim(edges)[1],]==1)-5
		edges<-edges[-1*dim(edges)[1],]
	}
	pp <- prop.part(phy)
	oldregime=root.state
	mm<-dim(edges)
	k<-length(6:mm[2])

	for(flow.index in sequence(dim(flow)[1])) {
		recipient.taxa <- strsplit(flow$recipient, ",")[[1]]
		donor.taxa <- strsplit(flow$donor, ",")[[1]]
		
		recipient.indices <- which(phy$tip.label %in% recipient.taxa) #these are node numbers
		donor.indices <- which(phy$tip.label %in% donor.taxa)

		if(length(recipient.indices)!=length(recipient.taxa)) {
			stop(paste("Tried to find ", flow$recipient[flow.index], " but failed; make sure the taxon names in the flow dataframe recipient match that of your tree", sep=""))
		}
		if(length(donor.indices)!=length(donor.taxa)) {
			stop(paste("Tried to find ", flow$donor[flow.index], " but failed; make sure the taxon names in the flow dataframe recipient match that of your tree", sep=""))
		}
		recipient.crown.node <- getMRCA(phy, tip=recipient.taxa)
		if(is.null(recipient.crown.node)) { #must be a single recipient taxon
			recipient.crown.node <- which(phy$tip.label == recipient.taxa[1])
		}	
	#	recipient.stem.node <- edges[which(edges[,3]==recipient.crown.node),2] #the parent with the main gene flow, not the one with flow rate of m
		recipient.stem.node <- EdgesGetOneNodeNumberRootward(recipient.crown.node, edges)
		parent.1 <- recipient.stem.node
		donor.crown.node <- getMRCA(phy, tip=donor.taxa)
		if(is.null(donor.crown.node)) { #must be a single donor taxon
			donor.crown.node <- which(phy$tip.label == donor.taxa[1])
		}
#		parent.1.regime <- which(edges[parent.1,6:(k+5)]==1)
		parent.1.regime <- which(edges[EdgesTipwardNumberToRow(parent.1, edges),6:(k+5)]==1)
		parent.2.tipward = donor.crown.node
#		parent.2.rootward = edges[donor.crown.node, 2]
		parent.2.rootward = EdgesGetOneNodeNumberRootward(donor.crown.node, edges)
#		parent.2.tipward.regime <- which(edges[donor.crown.node,6:(k+5)]==1)
		parent.2.tipward.regime <- which(edges[EdgesTipwardNumberToRow(donor.crown.node, edges),6:(k+5)]==1)
#		parent.2.rootward.regime <- which(edges[which(edges[,3]==donor.crown.node),6:(k+5)]==1)
		parent.2.rootward.regime <- which(edges[EdgesRootwardNumberToRow(donor.crown.node, edges),6:(k+5)]==1)
		parent.2.regime <- c(parent.2.rootward.regime, parent.2.tipward.regime)[which.min(abs(flow$time.from.root.recipient[flow.index] - edges[EdgesTipwardNumberToRow(donor.crown.node, edges),c(4:5)]))] #finds which end is closer to the time of the hybridization event, uses regime from that end
		#parent.2 is somewhere between donor.crown.node and the node below (that is, it occurs on this edge) 
		
		alpha1 <- alpha[parent.1.regime]
		alpha2 <- alpha[parent.2.regime]
		sigma1 <- sqrt(sigma.sq[parent.1.regime])
		sigma2 <- sqrt(sigma.sq[parent.2.regime])
		weights<-GetParentWeights(alpha1, alpha2, sigma1, sigma2, flow$m[flow.index], scaling.by.sd=scaling.by.sd, time.from.root.recipient=flow$time.from.root.recipient[flow.index])
		
		for (hybrid.local.index in sequence(length(recipient.indices))) {
			hybrid.index <- recipient.indices[hybrid.local.index]
			hybrid.regime <- which(edges[EdgesTipwardNumberToRow(hybrid.index, edges),6:(k+5)]==1)
			alpha.hybrid <- alpha[hybrid.regime]
			if(force.hybrid.alpha) {
				alpha.hybrid<-sum(c(alpha1, alpha2) * weights)
			}
			sigma.hybrid <- sqrt(sigma.sq[hybrid.regime])
			if(force.hybrid.sigma) {
				sigma.hybrid<-sum(c(sigma1, sigma2) * weights)
			}
			
			for (donor.local.index in sequence(length(donor.indices))) {
				donor.index <- donor.indices[donor.local.index]
				alpha.scaling <- exp(- AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = hybrid.index, rootward.node=Ntip(phy)+1, tipward.height = edges[EdgesTipwardNumberToRow(hybrid.index,edges),5], rootward.height=0, Rate.mat=Rate.mat, root.state=root.state)$alpha 
				- AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = donor.index, rootward.node=Ntip(phy)+1, tipward.height = edges[EdgesTipwardNumberToRow(donor.index,edges),5], rootward.height=0, Rate.mat=Rate.mat, root.state=root.state)$alpha
				)
				#so weight the vcv by original and what we get from other path
				v.modified[hybrid.index, donor.index] <- (weights[1] * vcv[hybrid.index, donor.index]) + (weights[2] * alpha.scaling * AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = parent.2.tipward, tipward.height = flow$time.from.root.recipient[flow.index], rootward.height=0, Rate.mat=Rate.mat, root.state=root.state)$variance)
				v.modified[donor.index, hybrid.index] <- v.modified[hybrid.index, donor.index]
			}
			descendants.of.parent <- getDescendants(phy, parent.1) 
			descendants.of.parent <- descendants.of.parent[!(descendants.of.parent %in% recipient.indices)] #only want to look at the covariance between the nonhybrids and the hybrids
			for (parent.descendant.index in sequence(length(descendants.of.parent))) {
				parent.descendant <- descendants.of.parent[parent.descendant.index]
				#so idea here is that you decrease the covariance of the hybrid and parent 1, but this doesn't go down to zero; it is based on what it'd be if the hybrid were at the other parent
				branch1 <- AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = parent.descendant, rootward.node=Ntip(phy)+1, tipward.height = edges[EdgesTipwardNumberToRow(parent.descendant,edges),5], rootward.height=0, Rate.mat=Rate.mat, root.state=root.state)
				alpha.sum.1 <- branch1$alpha
				second.term.1 <- branch1$second.term
				branchh <- AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = hybrid.index, rootward.node=parent.2.rootward, tipward.height = edges[EdgesTipwardNumberToRow(hybrid.index,edges),5], rootward.height=flow$time.from.root.recipient[flow.index], Rate.mat=Rate.mat, root.state=root.state)
				alpha.sum.h.to.2 <- branchh$alpha
				second.term.h.to.2 <- branchh$second.term				
				branchh2 <- AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = parent.2.tipward, tipward.height = flow$time.from.root.recipient[flow.index], rootward.height=0, Rate.mat=Rate.mat, root.state=root.state)
				alpha.sum.h2node.to.root <- branchh2$alpha
				second.term.h2node.to.root <- branchh2$second.term
				node.connecting.p1.and.h2.in.pos2 <- getMRCA(phy, tip=c(donor.crown.node, parent.descendant))
				
				second.term.h2.plus.parent1.to.root <- AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = node.connecting.p1.and.h2.in.pos2, tipward.height = edges[EdgesTipwardNumberToRow(node.connecting.p1.and.h2.in.pos2,edges),5] , rootward.height=0, Rate.mat=Rate.mat, root.state=root.state)$second.term
				covariance.parent.hybrid.on.2 <- exp(-(alpha.sum.1 + alpha.sum.h.to.2 + alpha.sum.h2node.to.root)) * (second.term.h2.plus.parent1.to.root) #Beaulieu et al equation 13
				v.modified[parent.descendant, hybrid.index] <- (weights[1] * v.modified[parent.descendant, hybrid.index]) + (weights[2] * covariance.parent.hybrid.on.2)
				v.modified[hybrid.index, parent.descendant] <- v.modified[parent.descendant, hybrid.index]
			}
			
			#so idea here is to do weighted average of the variance for both paths of the hybrid, through parent 1 and parent 2. The first path is straightfoward. 
			#The second one: we go from hybrid down to where it currently attaches, and stop. 
			#Then we go from where its other parent is down to the root, but starting partway down that first edge
			#Add up these variances, plus add on a weighted V_h
			hybrid.row <- EdgesTipwardNumberToRow(hybrid.index, edges)
			alpha.scaling <- exp(-2 * AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = hybrid.index, rootward.node=Ntip(phy)+1, tipward.height = edges[hybrid.row,5], rootward.height=0, Rate.mat=Rate.mat, root.state=root.state)$alpha)
			v.modified[hybrid.index, hybrid.index] <- alpha.scaling * (weights[1] * AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = hybrid.index, tipward.height = edges[hybrid.row,5], rootward.height=0, Rate.mat=Rate.mat, root.state=root.state)$variance + weights[2] * ( 
				AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = hybrid.index, rootward.node=parent.1, tipward.height = edges[hybrid.row,5], rootward.height=flow$time.from.root.recipient[flow.index], Rate.mat=Rate.mat, root.state=root.state)$variance + 
				AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = parent.2.tipward, tipward.height = flow$time.from.root.recipient[flow.index], rootward.height=0, Rate.mat=Rate.mat, root.state=root.state)$variance
			)) + vh*exp(-2*alpha.hybrid*flow$time.from.root.recipient[flow.index])

		}
		
		
	}
	return(v.modified)
}

GetRoot <- function(phy) {
	
}

AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint <- function(phy, edges, tipward.node, rootward.node=Ntip(phy)+1, tipward.height, rootward.height=0, Rate.mat, root.state) {
	alpha=Rate.mat[1,]
	sigma.sq=Rate.mat[2,]
	done <- FALSE
	variance.sum <- 0
	alpha.sum <- 0
	weight.vector <- rep(0, length(alpha))
	second.term.sum <- 0 #second term from Beaulieu et al. 2012 equation 13, but only for one lineage (so to get the overall second term, add one of these for each lineage i, j)
	anc = EdgesGetOneNodeNumberRootward(tipward.node, edges)
	tipward.row <- EdgesTipwardNumberToRow(tipward.node, edges)
	oldtime=edges[tipward.row,4]
	newtime=edges[tipward.row,5]
	mm<-dim(edges)
	k<-length(6:mm[2])
	if(anc != (Ntip(phy)+1)){
		start=which(edges[,3]==EdgesTipwardNumberToRow(anc, edges))
		oldregime=which(edges[start,6:(k+5)]==1)
	}
	else{
		#For the root:
		oldregime=root.state
	}	
	newregime=which(edges[tipward.row,6:(k+5)]==1)
	if(oldtime<rootward.height) {
		done<-TRUE
	}
	oldtime <- max(oldtime, rootward.height)
	if(oldregime==newregime){
		variance.sum <- variance.sum + sigma.sq[oldregime]*((exp(2*alpha[oldregime]*tipward.height)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))
		if(tipward.height<oldtime) {
			stop("problem with time sign")
		}
		alpha.sum <- alpha.sum + alpha[oldregime]*(tipward.height - oldtime)
		second.term.sum <- second.term.sum + ((sigma.sq[oldregime])) * (exp(2*alpha[oldregime] * tipward.height) - exp(2*alpha[oldregime] * oldtime)) / (2 * alpha[oldregime])
	} 	else{
		a.time <- tipward.height - mean(c(newtime, oldtime))
		halftime=newtime-((newtime-oldtime)/2)
		if(a.time<0) { #we're starting out rootward of the half time
			variance.sum <- variance.sum + sigma.sq[oldregime]*((exp(2*alpha[oldregime]*tipward.height)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))
			alpha.sum <- alpha.sum + alpha[oldregime]*(tipward.height - oldtime)
			weight.vector[oldregime] <- weight.vector[oldregime] + alpha[oldregime]*(tipward.height - oldtime)
			second.term.sum <- second.term.sum + ((sigma.sq[oldregime])) * (exp(2*alpha[oldregime] * tipward.height) - exp(2*alpha[oldregime] * oldtime)) / (2 * alpha[oldregime])
		if(tipward.height<oldtime) {
			stop("problem with time sign")
		}

		} else {
			variance.sum <- variance.sum + sigma.sq[oldregime]*((exp(2*alpha[oldregime]*halftime)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))	
			alpha.sum <- alpha.sum + alpha[oldregime]*(halftime - oldtime)
			weight.vector[oldregime] <- weight.vector[oldregime] + alpha[oldregime]*(halftime - oldtime)
			second.term.sum <- second.term.sum + ((sigma.sq[oldregime])) * (exp(2*alpha[oldregime] * halftime) - exp(2*alpha[oldregime] * oldtime)) / (2 * alpha[oldregime])
		if(halftime<oldtime) {
			stop("problem with time sign")
		}

			variance.sum <- variance.sum + sigma.sq[newregime]*((exp(2*alpha[newregime]*tipward.height)-exp(2*alpha[newregime]*tipward.height))/(2*alpha[newregime]))		
			alpha.sum <- alpha.sum + alpha[newregime]*(tipward.height - halftime)
			weight.vector[newregime] <- weight.vector[newregime] + alpha[newregime]*(tipward.height - halftime)
			second.term.sum <- second.term.sum + ((sigma.sq[newregime])) * (exp(2*alpha[newregime] * tipward.height) - exp(2*alpha[newregime] * halftime)) / (2 * alpha[newregime])
		if(halftime>tipward.height) {
			stop("problem with time sign")
		}

		}
	}
	current.node <- tipward.node
	if(EdgesGetOneNodeNumberRootward(current.node, edges)==rootward.node) {
		done<-TRUE
	}
	while(!done) {
		current.node <- EdgesGetOneNodeNumberRootward(current.node, edges) #move rootward
		current.row <- EdgesTipwardNumberToRow(current.node, edges)
		anc =  EdgesGetOneNodeNumberRootward(current.node, edges)

		if(anc == rootward.node || anc == (Ntip(phy)+1)) { #if we have gone down to the root
			done <- TRUE #so this'll be our last run here
		}
		oldtime=edges[current.row,4]
		newtime=edges[current.row,5]
		if(anc != (Ntip(phy)+1)){
			start=which(edges[,3]==EdgesTipwardNumberToRow(anc, edges))
			oldregime=which(edges[start,6:(k+5)]==1)
		}
		else{
			#For the root:
			oldregime=root.state
		}	
		newregime=which(edges[current.row,6:(k+5)]==1)
		oldtime <- max(oldtime, rootward.height)

		
		if(oldregime==newregime){
			variance.sum <- variance.sum + sigma.sq[oldregime]*((exp(2*alpha[oldregime]*newtime)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))
			alpha.sum <- alpha.sum + alpha[oldregime]*(newtime - oldtime)
			weight.vector[oldregime] <- weight.vector[oldregime] + alpha[oldregime]*(newtime - oldtime)
			second.term.sum <- second.term.sum + ((sigma.sq[oldregime])) * (exp(2*alpha[oldregime] * newtime) - exp(2*alpha[oldregime] * oldtime)) / (2 * alpha[oldregime])
		if(newtime<oldtime) {
			stop("problem with time sign 253")
		}


		} 	else{
			b.time <- oldtime - mean(c(newtime, oldtime))
			halftime=newtime-((newtime-oldtime)/2)
			if(oldtime > mean(c(newtime, edges[current.row,4]))) { #we're starting out tipward of the half time
				variance.sum <- variance.sum + sigma.sq[newregime]*((exp(2*alpha[newregime]*tipward.height)-exp(2*alpha[newregime]*oldtime))/(2*alpha[newregime]))
				alpha.sum <- alpha.sum + alpha[newregime] * (newtime - oldtime)
				weight.vector[newregime] <- weight.vector[newregime] + alpha[newregime] * (newtime - oldtime)
				second.term.sum <- second.term.sum + ((sigma.sq[newregime])) * (exp(2*alpha[newregime] * newtime) - exp(2*alpha[newregime] * oldtime)) / (2 * alpha[newregime])
			if(newtime<oldtime) {
			stop("problem with time sign 267")
		}

			} else {
				variance.sum <- variance.sum + sigma.sq[oldregime]*((exp(2*alpha[oldregime]*halftime)-exp(2*alpha[oldregime]*oldtime))/(2*alpha[oldregime]))		
				alpha.sum <- alpha.sum + alpha[oldregime]*(halftime - oldtime)	
				weight.vector[oldregime] <- weight.vector[oldregime] + alpha[oldregime]*(halftime - oldtime)
				second.term.sum <- second.term.sum + ((sigma.sq[oldregime])) * (exp(2*alpha[oldregime] * halftime) - exp(2*alpha[oldregime] * oldtime)) / (2 * alpha[oldregime])
		if(halftime<oldtime) {
			stop("problem with time sign")
		}

				variance.sum <- variance.sum + sigma.sq[newregime]*((exp(2*alpha[newregime]*tipward.height)-exp(2*alpha[newregime]*tipward.height))/(2*alpha[newregime]))	
				alpha.sum <- alpha.sum + alpha[newregime]*(newtime - halftime)		
				weight.vector[newregime] <- weight.vector[newregime] + alpha[newregime]*(newtime - halftime)	
				second.term.sum <- second.term.sum + ((sigma.sq[newregime])) * (exp(2*alpha[newregime] * newtime) - exp(2*alpha[newregime] * halftime)) / (2 * alpha[newregime])
		if(newtime<halftime) {
			stop("problem with time sign")
		}
			}
		}
		
	}
	return(list(variance=variance.sum, alpha=alpha.sum, weights=weight.vector, second.term = second.term.sum ))
}



GetOUVariance <- function(sigma, alpha, time.from.root) {
	return( (sigma^2) * (1 - exp(-2 * alpha *time.from.root)) / (2 * alpha))	
}


GetParentWeights <- function(alpha1, alpha2, sigma1, sigma2, m, scaling.by.sd, time.from.root.recipient) { #assumes m is parameter for proportion flow from parent 2: that is, the main history is from parent 1, and m (0 - 0.5 range, generally) comes from parent 2
	sd1 <- sqrt(GetOUVariance(sigma1, alpha1, time.from.root.recipient))
	sd2 <- sqrt(GetOUVariance(sigma2, alpha2, time.from.root.recipient))
	sd1.wt <- scaling.by.sd * ((1-m)/sd1)/ (((m)/sd2+(1-m)/sd1)) + (1-scaling.by.sd)*(1-m)
	sd2.wt <- scaling.by.sd * ((m)/sd2)/ (((m)/sd2+(1-m)/sd1)) + (1-scaling.by.sd)*(m)
	return(c(sd1.wt, sd2.wt))
}