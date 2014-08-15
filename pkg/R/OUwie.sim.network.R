##OUwie Simulator##

#written by Jeremy M. Beaulieu
#modified by Brian O'Meara, Paul Blischak, and Tony Jhwueng

#Simulates the Hansen model of continuous characters evolving under discrete selective 
#regimes. The input is a tree of class "phylo" that has the regimes as internal node labels 
#and a data file the contains the regime states for each species. The trait file must be in 
#the following order: Species names then Regime. The user must specify the parameters values
#for each simulation (i.e. alpha, sigma.sq, theta0, theta). 

##The following examples assume 2 selective regimes and different models can be specified:
#single rate Brownian motion BM1: alpha=c(0,0); sigma.sq=c(0.9); theta0=0; theta=c(0,0)
#two rate Brownian motion BMS: alpha=c(0,0); sigma.sq=c(0.45,.9); theta0=0; theta=c(0,0)
#global OU (OU1): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=1; theta=c(1,1)
#normal OU (OUSM): alpha=c(0.1,0.1); sigma.sq=c(0.9,0.9); theta0=0; theta=c(1,2)
#multiple sigmas (OUSMV): alpha=c(0.1,0.1); sigma.sq=c(0.45,0.9); 
#multiple alphas (OUSMA): alpha=c(0.5,0.1); sigma.sq=c(0.9,0.9); theta0=0; theta=c(1,2)
#multiple alphas and sigmas (OUSMVA): alpha=c(0.5,0.1); sigma.sq=c(0.45,0.9); theta0=0; theta=c(1,2)

OUwie.sim.network <- function(phy, data.in=NULL, simmap.tree=FALSE, scaleHeight=FALSE, alpha, sigma.sq, theta0, theta, flow, scaling.by.sd=FALSE){

	if(simmap.tree==FALSE){
		data.original <- data.in
		#This is annoying, but the second column has to be in there twice otherwise, error.
		data.shrunk<-data.frame(data.in[,2], data.in[,2], row.names=data.in[,1])
		data.shrunk<-data.shrunk[phy$tip.label,]
				
		

		
		phy.shrunk <- phy
		recipient.subtrees <- list()
		recipient.stem.lengths <- vector()
		parent.1.regimes <- c()
		parent.2.regimes <- c()
		hybrid.regimes <- c()
		print(phy.shrunk)

		for(flow.index in sequence(dim(flow)[1])) {


			n=max(phy.shrunk$edge[,1])
			ntips=length(phy.shrunk$tip.label)
		
			int.states<-factor(phy.shrunk$node.label)
			phy.shrunk$node.label=as.numeric(int.states)
			tip.states<-factor(data.shrunk[,1])
			data.shrunk[,1]<-as.numeric(tip.states)
			k<-length(levels(int.states))
		
			regime=matrix(rep(0,(n-1)*k), n-1, k)

			#Obtain root state and internal node labels
			root.state<-phy.shrunk$node.label[1]
			int.state<-phy.shrunk$node.label[-1]
		
			#New tree matrix to be used for subsetting regimes
			edges=cbind(c(1:(n-1)),phy.shrunk$edge,nodeHeights(phy.shrunk))
			if(scaleHeight==TRUE){
				edges[,4:5]<-edges[,4:5]/max(nodeHeights(phy.shrunk))
			}
			edges=edges[sort.list(edges[,3]),]

			mm<-c(data.shrunk[,1],int.state)
		
			regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
			#Generates an indicator matrix from the regime vector
			for (i in 1:length(mm)) {
				regime[i,mm[i]] <- 1 
			}
			#Finishes the edges matrix
			edges=cbind(edges,regime)
		
			#Resort the edge matrix so that it looks like the original matrix order
			edges=edges[sort.list(edges[,1]),]
		
			oldregime=root.state
		
			alpha=alpha
			alpha[alpha==0] = 1e-10
			sigma=sqrt(sigma.sq)
			theta=theta

			x <- matrix(0, n, 1)
			TIPS <- 1:ntips
			ROOT <- ntips + 1L
			x[ROOT,] <- theta0
			
			recipient.taxa <- strsplit(flow$recipient[flow.index], ",")[[1]]
			donor.taxa <- strsplit(flow$donor[flow.index], ",")[[1]]
		
			recipient.indices <- which(phy$tip.label %in% recipient.taxa) #these are node numbers
			donor.indices <- which(phy$tip.label %in% donor.taxa)

			if(length(recipient.indices)!=length(recipient.taxa)) {
				stop(paste("Tried to find ", flow$recipient[flow.index], " but failed; make sure the taxon names in the flow dataframe recipient match that of your tree", sep=""))
			}
			if(length(donor.indices)!=length(donor.taxa)) {
				stop(paste("Tried to find ", flow$donor[flow.index], " but failed; make sure the taxon names in the flow dataframe recipient match that of your tree", sep=""))
			}
			recipient.crown.node <- getMRCA(phy.shrunk, tip=recipient.taxa)
			if(is.null(recipient.crown.node)) { #must be a single recipient taxon
				recipient.crown.node <- which(phy.shrunk$tip.label == recipient.taxa[1])
				recipient.subtrees<-append(recipient.subtrees, recipient.taxa[1])
			}	else {
				recipient.subtrees <- append(recipient.subtrees, extract.clade(phy.shrunk, recipient.crown.node))
			}
			data.shrunk <- data.shrunk[!(data.shrunk[,1] %in% recipient.taxa),]
			if(length(recipient.taxa)==1) {
				phy.shrunk$tip.label[which(phy.shrunk$tip.label==recipient.taxa[1])] <- '[1_tip]'
			} else {
				phy.shrunk <- drop.tip(phy.shrunk, recipient.taxa, subtree=TRUE)
			}
			print("postdrop")
			print(phy.shrunk)
			shrunk.tip.id <- which(grepl("_tip", phy.shrunk$tip.label))
			hybrid.regime <- which(edges[EdgesTipwardNumberToRow(recipient.crown.node, edges),6:(k+5)]==1)
			hybrid.regimes <- append(hybrid.regimes, hybrid.regime)
			parent.1.regime <- which(edges[EdgesTipwardNumberToRow(EdgesGetOneNodeNumberRootward(recipient.crown.node, edges), edges),6:(k+5)]==1)
			original.length <- phy.shrunk$edge.length[which(phy.shrunk$edge[,2]==shrunk.tip.id)]
			new.rootward.length <-  original.length - ( nodeheight(phy.shrunk, shrunk.tip.id) - flow$time.from.root.recipient[flow.index] )
			new.tipward.length <- original.length - new.rootward.length
			recipient.stem.lengths <- append(recipient.stem.lengths, new.tipward.length)
			phy.shrunk$edge.length[which(phy.shrunk$edge[,2]==shrunk.tip.id)] <- original.length - ( nodeheight(phy.shrunk, shrunk.tip.id) - flow$time.from.root.recipient[flow.index] )
			phy.shrunk$tip.label[shrunk.tip.id] <- paste("Hybrid_flow_index_recipient_", flow.index, sep="") #so this is a node where the hybrid lineage attaches to parent 1
			print(phy.shrunk)
			donor.crown.node <- getMRCA(phy.shrunk, tip=donor.taxa)
			if(is.null(donor.crown.node)) { #must be a single donor taxon
				donor.crown.node <- which(phy.shrunk$tip.label == donor.taxa[1])
			}
			parent.2.tipward = donor.crown.node
			parent.2.rootward = EdgesGetOneNodeNumberRootward(donor.crown.node, edges)
			parent.2.tipward.regime <- which(edges[EdgesTipwardNumberToRow(donor.crown.node, edges),6:(k+5)]==1)
			parent.2.rootward.regime <- which(edges[EdgesTipwardNumberToRow(parent.2.rootward, edges),6:(k+5)]==1)
			parent.1.regimes <- append(parent.1.regimes, parent.1.regime)
			parent.2.regimes <- append(parent.2.regimes, parent.2.rootward.regime)
			phy.shrunk<-bind.tip(phy.shrunk, tip.label= paste("Hybrid_flow_index_donor_", flow.index, sep=""), where=donor.crown.node, edge.length=0, position=nodeheight(phy.shrunk, donor.crown.node) -  flow$time.from.root.donor[flow.index])#this is a node where the hybrid lineage attaches to parent 2
			new.data <- data.frame(parent.1.regime , parent.1.regime , stringsAsFactors=FALSE)
			colnames(new.data) <- colnames(data.shrunk)
			rownames(new.data) <- paste("Hybrid_flow_index_recipient_", flow.index, sep="")
			data.shrunk <- rbind(data.shrunk, new.data) #row.names=paste("Hybrid_flow_index_recipient_", flow.index, sep=""), # col.names=colnames(data.shrunk)
			new.data <- data.frame(parent.2.rootward.regime , parent.2.rootward.regime , stringsAsFactors=FALSE)
			colnames(new.data) <- colnames(data.shrunk)
			rownames(new.data) <- paste("Hybrid_flow_index_donor_", flow.index, sep="")
			data.shrunk <- rbind(data.shrunk, new.data) #row.names=paste("Hybrid_flow_index_recipient_", flow.index, sep=""), # col.names=colnames(data.shrunk)
		}
		
		#whew. Now we've pruned off the hybrid taxa and added the placeholders on phy.shrunk. Now simulate on that as normal
		
		n=max(phy.shrunk$edge[,1])
		ntips=length(phy.shrunk$tip.label)
	
		int.states<-factor(phy.shrunk$node.label)
		phy.shrunk$node.label=as.numeric(int.states)
		tip.states<-factor(data.shrunk[,1])
		data.shrunk[,1]<-as.numeric(tip.states)
		k<-length(levels(int.states))
	
		regime=matrix(rep(0,(n-1)*k), n-1, k)

		#Obtain root state and internal node labels
		root.state<-phy.shrunk$node.label[1]
		int.state<-phy.shrunk$node.label[-1]
	
		#New tree matrix to be used for subsetting regimes
		edges=cbind(c(1:(n-1)),phy.shrunk$edge,nodeHeights(phy.shrunk))
		if(scaleHeight==TRUE){
			edges[,4:5]<-edges[,4:5]/max(nodeHeights(phy.shrunk))
		}
		edges=edges[sort.list(edges[,3]),]

		mm<-c(data.shrunk[,1],int.state)
	
		regime <- matrix(0,nrow=length(mm),ncol=length(unique(mm)))
		#Generates an indicator matrix from the regime vector
		for (i in 1:length(mm)) {
			regime[i,mm[i]] <- 1 
		}
		#Finishes the edges matrix
		edges=cbind(edges,regime)
	
		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]
	
		oldregime=root.state
	
		alpha=alpha
		alpha[alpha==0] = 1e-10
		sigma=sqrt(sigma.sq)
		theta=theta

		x <- matrix(0, n, 1)
		TIPS <- 1:ntips
		ROOT <- ntips + 1L
		x[ROOT,] <- theta0

		
		
		for(i in 1:length(edges[,1])){
			anc = edges[i,2]
			desc = edges[i,3]
			oldtime=edges[i,4]
			newtime=edges[i,5]
			if(anc%in%edges[,3]){
				start=which(edges[,3]==anc)
				oldregime=which(edges[start,6:(k+5)]==1)
			}
			else{
				#For the root:
				oldregime=oldregime
			}	
			newregime=which(edges[i,6:(k+5)]==1)
			if(oldregime==newregime){
				x[edges[i,3],]=x[edges[i,2],]*exp(-alpha[oldregime]*(newtime-oldtime))+(theta[oldregime])*(1-exp(-alpha[oldregime]*(newtime-oldtime)))+sigma[oldregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[oldregime]*(newtime-oldtime)))/(2*alpha[oldregime]))
			}
			else{
				halftime=newtime+((oldtime-newtime)/2)
				epoch1=x[edges[i,2],]*exp(-alpha[oldregime]*(halftime-oldtime))+(theta[oldregime])*(1-exp(-alpha[oldregime]*(halftime-oldtime)))+sigma[oldregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[oldregime]*(halftime-oldtime)))/(2*alpha[oldregime]))
				oldtime=halftime
				newtime=newtime
				x[edges[i,3],]=epoch1*exp(-alpha[newregime]*(newtime-oldtime))+(theta[newregime])*(1-exp(-alpha[newregime]*(newtime-oldtime)))+sigma[newregime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[newregime]*(newtime-oldtime)))/(2*alpha[newregime]))
			}
		}
		
		sim.dat<-matrix(,ntips,3)
		sim.dat<-data.frame(sim.dat)
		
		sim.dat[,1]<-phy$tip.label
		sim.dat[,2]<-data[,1]
		sim.dat[,3]<-x[TIPS]

		colnames(sim.dat)<-c("Genus_species","Reg","X")	
		
		#and now to do the simulation for each hybrid origin
		
		
		for(flow.index in sequence(dim(flow)[1])) {
			parent.1.regime <- parent.1.regimes[flow.index]
			parent.2.regime <- parent.2.regimes[flow.index]
			hybrid.regime <- hybrid.regimes[flow.index]
			alpha1 <- alpha[parent.1.regime]
			alpha2 <- alpha[parent.2.regime]
			sigma1 <- sqrt(sigma.sq[parent.1.regime])
			sigma2 <- sqrt(sigma.sq[parent.2.regime])
			weights<-GetParentWeights(alpha1, alpha2, sigma1, sigma2, flow$m[flow.index], scaling.by.sd=scaling.by.sd, time.from.root.recipient=flow$time.from.root.recipient[flow.index])
			parent.1.value <- sim.dat[which(sim.dat[,1]==paste("Hybrid_flow_index_recipient_", flow.index, sep="")),3]
			parent.2.value <- sim.dat[which(sim.dat[,1]==paste("Hybrid_flow_index_donor_", flow.index, sep="")),3]
			theta0 <- weights[1] * parent.1.value + weights[2] * parent.2.value
			theta0.post.root.edge <- theta0*exp(-alpha[hybrid.regime]*(recipient.stem.lengths[flow.index]))+(theta[hybrid.regime])*(1-exp(-alpha[hybrid.regime]*(recipient.stem.lengths[flow.index])))+sigma[hybrid.regime]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[hybrid.regime]*(recipient.stem.lengths[flow.index])))/(2*alpha[hybrid.regime]))
			recipient.taxa <- strsplit(flow$recipient[flow.index], ",")[[1]]
			hybrid.sim <- theta0.post.root.edge
			if(length(recipient.taxa)>1) {
				hybrid.sim <- OUwie.sim(recipient.subtrees[[flow.index]], data=data.original[which(data.original[,1] %in% recipient.taxa),], simmap.tree=simmap.tree, scaleHeight=scaleHeight, alpha, sigma.sq, theta0, theta)
			}
			sim.dat <- rbind(sim.dat, hybrid.sim)
		}
		
		
			
	}
	if(simmap.tree==TRUE){
		stop("Use of simmap trees is currently not supported in OUwie for networks")
		n=max(phy$edge[,1])
		ntips=length(phy$tip.label)
		
		k=length(colnames(phy$mapped.edge))

		regimeindex<-colnames(phy$mapped.edge)
		##Begins the construction of the edges matrix -- similar to the ouch format##
		#Makes a vector of absolute times in proportion of the total length of the tree
		branch.lengths=rep(0,(n-1))
		branch.lengths[(ntips+1):(n-1)]=branching.times(phy)[-1]/max(branching.times(phy))
		
		#Obtain root state and internal node labels
		root.state<-which(colnames(phy$mapped.edge)==names(phy$maps[[1]][1]))
		
		#New tree matrix to be used for subsetting regimes
		edges=cbind(c(1:(n-1)),phy$edge,nodeHeights(phy))
		if(scaleHeight==TRUE){
			edges[,4:5]<-edges[,4:5]/max(nodeHeights(phy))
		}
		edges=edges[sort.list(edges[,3]),]
		
		#Resort the edge matrix so that it looks like the original matrix order
		edges=edges[sort.list(edges[,1]),]
		
		oldregime=root.state
		oldtime=0
		
		alpha=alpha
		sigma=sqrt(sigma.sq)
		theta=theta
		
		n.cov=matrix(rep(0,n*n), n, n)
		nodecode=matrix(c(ntips+1,1),1,2)
		
		x <- matrix(0, n, 1)
		TIPS <- 1:ntips
		ROOT <- ntips + 1L
		x[ROOT,] <- theta0
		
		for(i in 1:length(edges[,1])){
			
			if(scaleHeight==TRUE){
				currentmap<-phy$maps[[i]]/max(nodeHeights(phy))
			}
			else{
				currentmap<-phy$maps[[i]]
			}
			oldtime=edges[i,4]
			
			if(length(phy$maps[[i]])==1){
				regimeduration<-currentmap[1]
				newtime<-oldtime+regimeduration
				regimenumber<-which(colnames(phy$mapped.edge)==names(currentmap)[1])
				x[edges[i,3],]=x[edges[i,2],]*exp(-alpha[regimenumber]*(newtime-oldtime))+(theta[regimenumber])*(1-exp(-alpha[regimenumber]*(newtime-oldtime)))+sigma[regimenumber]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[regimenumber]*(newtime-oldtime)))/(2*alpha[regimenumber]))
			}
			if(length(phy$maps[[i]])>1){
				regimeduration<-currentmap[1]
				newtime<-oldtime+regimeduration
				regimenumber<-which(colnames(phy$mapped.edge)==names(currentmap)[1])
				x[edges[i,3],]=x[edges[i,2],]*exp(-alpha[regimenumber]*(newtime-oldtime))+(theta[regimenumber])*(1-exp(-alpha[regimenumber]*(newtime-oldtime)))+sigma[regimenumber]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[regimenumber]*(newtime-oldtime)))/(2*alpha[regimenumber]))				
				oldtime<-newtime
				for (regimeindex in 2:length(currentmap)){
					regimeduration<-currentmap[regimeindex]
					newtime<-oldtime+regimeduration
					regimenumber<-which(colnames(phy$mapped.edge)==names(currentmap)[regimeindex])
					x[edges[i,3],]=x[edges[i,3],]*exp(-alpha[regimenumber]*(newtime-oldtime))+(theta[regimenumber])*(1-exp(-alpha[regimenumber]*(newtime-oldtime)))+sigma[regimenumber]*rnorm(1,0,1)*sqrt((1-exp(-2*alpha[regimenumber]*(newtime-oldtime)))/(2*alpha[regimenumber]))
					oldtime<-newtime
					newregime<-regimenumber
				}
			}
		}
		
		sim.dat<-matrix(,ntips,2)
		sim.dat<-data.frame(sim.dat)
		
		sim.dat[,1]<-phy$tip.label
		sim.dat[,2]<-x[TIPS,]
		
		colnames(sim.dat)<-c("Genus_species","X")
	}
	sim.dat
}



