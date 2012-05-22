
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

#Steps:
#Create possible regimes from regime.vector -- done
#Write tree traversal to obtain possible regimes for each edge -- scores each in edges.mat
#Output is edge.ouwie