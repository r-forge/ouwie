#Weight matrix generator taken from Butler and King (2004) and modified to allow multiple alpha parameters. And networks



weight.mat.network<-function(phy, edges, Rate.mat, root.state, simmap.tree=FALSE, scaleHeight=FALSE, assume.station=TRUE, flow, scaling.by.sd=1, force.hybrid.alpha=FALSE, force.hybrid.sigma=FALSE, vh=0, bt=1){
	W <- weight.mat(phy, edges, Rate.mat, root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight, assume.station=assume.station, standardizeRowSums=FALSE)
	print("W")
	print(head(W))
	#This gives us weight matrix up main path. Now have to adjust for the hybrid taxa
	mm<-dim(edges)
	k<-length(6:mm[2])
	alpha=Rate.mat[1,]
	sigma=Rate.mat[2,]


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
		sigma1 <- sigma[parent.1.regime]
		sigma2 <- sigma[parent.2.regime]
		weights<-GetParentWeights(alpha1, alpha2, sigma1, sigma2, flow$m[flow.index], scaling.by.sd=scaling.by.sd, time.from.root.recipient=flow$time.from.root.recipient[flow.index])

		
		for (hybrid.local.index in sequence(length(recipient.indices))) {
			hybrid.index <- recipient.indices[hybrid.local.index]
			W.original.row.for.debugging <- W[hybrid.index,]

			hybrid.row <- EdgesTipwardNumberToRow(hybrid.index, edges)
			hybrid.regime <- which(edges[hybrid.row,6:(k+5)]==1)
			alpha.hybrid <- alpha[hybrid.regime]
			if(force.hybrid.alpha) {
				alpha.hybrid<-sum(c(alpha1, alpha2) * weights)
			}
			sigma.hybrid <- sigma[hybrid.regime]
			if(force.hybrid.sigma) {
				sigma.hybrid<-sum(c(sigma1, sigma2) * weights)
			}
			
			full.path.1 <- AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = hybrid.index, rootward.node=Ntip(phy)+1, tipward.height = edges[hybrid.row,5], rootward.height=0, Rate.mat=Rate.mat, root.state=root.state)
			post.hybrid.path.2 <- AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = hybrid.index, rootward.node=parent.1, tipward.height = edges[hybrid.row,5], rootward.height=flow$time.from.root.recipient[flow.index], Rate.mat=Rate.mat, root.state=root.state)
			pre.hybrid.path.2 <- AddVarianceAndAlphaAndWeightsFromPartwayAlongEdgeDownToAnotherPoint(phy, edges, tipward.node = parent.2.tipward, tipward.height = flow$time.from.root.recipient[flow.index], rootward.height=0, Rate.mat=Rate.mat, root.state=root.state)
			for (regime.index in sequence(dim(W)[2])) {
				W[hybrid.index, regime.index] <- weights[1] * exp(-full.path.1$alpha) * full.path.1$weights[regime.index] + weights[2] * exp(- (post.hybrid.path.2$alpha + pre.hybrid.path.2$alpha )) * (post.hybrid.path.2$weights[regime.index] + pre.hybrid.path.2$weights[regime.index])
			}
			print("weights")		
			print(cbind(matrix(W.original.row.for.debugging, ncol=1), matrix(W[hybrid.index,], ncol=1)))

		}
		
		
	}
	print("W modified")
	print(head(W))
	return(W)	
}
