#Weight matrix generator taken from Butler and King (2004) and modified to allow multiple alpha parameters

#written by Jeremy M. Beaulieu

weight.mat<-function(phy, edges, Rate.mat, root.state, assume.station=TRUE){
	
	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
	k=max(as.numeric(phy$node.label))
	edges=edges
	oldregime=root.state
	oldtime=0
	nodeweight1=rep(0,max(edges[,3]))
	nodeweight2=rep(0,max(edges[,3]))
	edges[,4]<-1-edges[,4]
	alpha=Rate.mat[1,]
	
	if(assume.station==TRUE){
		
		W1<-matrix(0,ntips,k)
		W2<-matrix(0,ntips,k)
		
		for(j in 1:k){
			
			n.w1=matrix(0, n, n)
			n.w2=matrix(0, n, n)
			nodecode=matrix(c(ntips+1,0,oldregime),1,3)
			#Weight calculated for the root
			#Weights for each species per regime
			for(i in 1:length(edges[,1])){
				anc = edges[i, 2]
				desc = edges[i, 3]
				newregime=which(edges[i,5:(k+4)]==1)
				current=edges[i,4]
				
				if(anc%in%nodecode[,1]){
					start=which(nodecode[,1]==anc)
					oldtime=nodecode[start,2]
					oldregime=nodecode[start,3]
				}
				else{
					newrow=c(anc,newtime,oldregime)
					nodecode=rbind(nodecode,newrow)
					oldtime=newtime
				}
				
				if(oldregime==newregime){
					newtime=current
					
					if(newregime==j){
						nodeweight1[i]=alpha[oldregime]*(newtime-oldtime)
						nodeweight2[i]=exp(alpha[oldregime]*newtime)-exp(alpha[oldregime]*oldtime)
					}
					else{
						nodeweight1[i]=0
						nodeweight2[i]=0
					}
				}
				else{
					newtime=current-((current-oldtime)/2)
					epoch1a=alpha[oldregime]*(newtime-oldtime)
					epoch1b=exp(alpha[oldregime]*newtime)-exp(alpha[oldregime]*oldtime)
					oldtime=newtime
					newtime=current
					epoch2a=alpha[newregime]*(newtime-oldtime)
					epoch2b=exp(alpha[newregime]*newtime)-exp(alpha[newregime]*oldtime)
					if(oldregime==j){
						nodeweight1[i]=epoch1a
						nodeweight2[i]=epoch1b
					}
					if(newregime==j){
						nodeweight1[i]=epoch2a
						nodeweight2[i]=epoch2b
					}
				}
				oldregime=newregime
				n.w1[edges[i,2],edges[i,3]]=nodeweight1[i]
				n.w2[edges[i,2],edges[i,3]]=nodeweight2[i]
				
			}
			W1[1:(ntips),j]<-weight.gen(n.w1, phy)
			W2[1:(ntips),j]<-weight.gen(n.w2, phy)	
		}
	}
	
	if(assume.station==FALSE){
		
		W1<-matrix(0,ntips,k+1)
		W2<-matrix(0,ntips,k+1)
		
		for(j in 1:k){
			n.w1=matrix(0, n, n)
			n.w2=matrix(0, n, n)
			nodecode=matrix(c(ntips+1,0,oldregime),1,3)
			#Weight calculated for the root
			#Weights for each species per regime
			
			for(i in 1:length(edges[,1])){
				
				anc = edges[i, 2]
				desc = edges[i, 3]
				newregime=which(edges[i,5:(k+4)]==1)
				current=edges[i,4]
				
				if(anc%in%nodecode[,1]){
					start=which(nodecode[,1]==anc)
					oldtime=nodecode[start,2]
					oldregime=nodecode[start,3]
				}
				else{
					newrow=c(anc,newtime,oldregime)
					nodecode=rbind(nodecode,newrow)
					oldtime=newtime
				}
				
				if(oldregime==newregime){
					newtime=current
					if(newregime==j){
						nodeweight1[i]=alpha[oldregime]*(newtime-oldtime)
						nodeweight2[i]=exp(alpha[oldregime]*newtime)-exp(alpha[oldregime]*oldtime)
					}
					else{
						nodeweight1[i]=0
						nodeweight2[i]=0
					}
				}
				else{
					newtime=current-((current-oldtime)/2)
					epoch1a=alpha[oldregime]*(newtime-oldtime)
					epoch1b=exp(alpha[oldregime]*newtime)-exp(alpha[oldregime]*oldtime)
					oldtime=newtime
					newtime=current
					epoch2a=alpha[newregime]*(newtime-oldtime)
					epoch2b=exp(alpha[newregime]*newtime)-exp(alpha[newregime]*oldtime)
					if(oldregime==j){
						nodeweight1[i]=epoch1a
						nodeweight2[i]=epoch1b
					}
					if(newregime==j){
						nodeweight1[i]=epoch2a
						nodeweight2[i]=epoch2b
					}
				}

				oldregime=newregime
				n.w1[edges[i,2],edges[i,3]]=nodeweight1[i]
				n.w2[edges[i,2],edges[i,3]]=nodeweight2[i]
			}
			W1[1:(ntips),j+1]<-weight.gen(n.w1, phy)
			W2[1:(ntips),j+1]<-weight.gen(n.w2, phy)
		}
	}
	if(assume.station==TRUE){
		W<-exp(-rowSums(W1))*W2
		W<-W/rowSums(W)
	}
	if(assume.station==FALSE){
		W<-exp(-rowSums(W1))*W2
		W[,1]<-exp(-rowSums(W1))
		W<-W/rowSums(W)
	}
	W
}

#Utility for building a summary matrix -- slow for now, need to work on speeding up
weight.gen<-function(mat,phy){
	
	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
	k=max(as.numeric(phy$node.label))
	
	S<-matrix(0,n,n)
	for(i in 1:n){
		nn<-Ancestors(phy,i)
		S[nn,i]<-1
	}

	#Convert to a sparse matrix
	S<-as.matrix.csr(S)
	#Create a matrix of sums for each node
	
	temp<-mat%*%S+mat
	#Remove S, nodecode, and n.cov matrices from memory
	
	rm(mat)
	rm(S)

	n.covsums=apply(as.matrix(temp), 2, sum)
	
	rm(temp)
	n.covsums
	
	n.covsums<-n.covsums[1:ntips]
	n.covsums
}

#Utility function for obtaining mrcas for each species pair
Ancestors<-function (x, node, type = c("all", "parent")){
    parents <- x$edge[, 1]
    child <- x$edge[, 2]
    pvector <- numeric(max(x$edge)) # parents
    pvector[child] <- parents    
    type <- match.arg(type)
    if (type == "parent") 
	return(pvector[node])
    res <- numeric(0)
    repeat {
        anc <- pvector[node]
        if (anc == 0) break
        res <- c(res, anc)
        node <- anc
    }
    res
}


