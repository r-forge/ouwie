#OU variance-covariance matrix generator adapted from Butler and King (2004)

#written by Jeremy M. Beaulieu and Jason Shapiro

vcv.ou<-function(phy,edges,Rate.mat, root.state){
	
	n=max(phy$edge[,1])
	ntips=length(phy$tip.label)
	k=max(as.numeric(phy$node.label))

	oldregime=root.state
	oldtime=1
	nodevar=rep(0,max(edges[,3]))
	alpha=Rate.mat[1,]
	sigma=Rate.mat[2,]
	n.cov=matrix(rep(0,n*n), n, n)
	
	nodecode=matrix(c(ntips+1,1,1),1,3)
	
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
			nodevar[i]=(sigma[oldregime]/(2*alpha[oldregime]))*(exp(-2*alpha[oldregime]*newtime)-exp(-2*alpha[oldregime]*oldtime))
		}
		else{

			newtime=current+((oldtime-current)/2)
			epoch1=(sigma[oldregime]/(2*alpha[oldregime]))*(exp(-2*alpha[oldregime]*newtime)-exp(-2*alpha[oldregime]*oldtime))
			oldtime=newtime
			newtime=current
			epoch2=(sigma[newregime]/(2*alpha[newregime]))*(exp(-2*alpha[newregime]*newtime)-exp(-2*alpha[newregime]*oldtime))
			nodevar[i]<-epoch1+epoch2
		}
	
		oldregime=newregime
		n.cov[edges[i,2],edges[i,3]]=nodevar[i]
	}

	#Remove nodecode matrix from memory
	rm(nodecode)

	#Begins building summary matrix
	A=matrix(rep(0,n*n), n, n)
	A[tree$edge]=1
	mm<-c(1:n)
	Nt<-t(t(A)*mm)
	#Convert to a sparse matrix
	Nt<-as.matrix.csr(Nt)
	rm(A)
	
	#Generates the mystical S matrix
	S<-matrix(0,n,n)
	for(i in 1:n){
		nn<-Ancestors(phy,i)
		S[nn,i]<-1
	}
	#Convert to a sparse matrix
	S<-as.matrix.csr(S)
	
	#Create a matrix of sums for each node
	temp<-n.cov%*%S+n.cov
	n.covsums=apply(as.matrix(temp), 2, sum)
	rm(temp)
	rm(n.cov)

	#Generates a matrix that lists the descendants for each ancestral node
	H.temp=Nt%*%S+Nt
	rm(Nt)
	rm(S)
	H.mat=as.matrix(H.temp[(ntips+1):n,])
	
	rm(H.temp)
	
	vcv<-matrix(0,ntips,ntips)
	#Enters variances
	diag(vcv)=n.covsums[1:ntips]
	#Enters covariances
	for(i in 1:length(H.mat[,1])){
		temp=unique(H.mat[i,])
		temp=temp[temp!=0]
		tempR=which(H.mat[i,]==temp[1])
		tempR=subset(tempR,tempR%in%c(1:ntips))
		tempL=which(H.mat[i,]==temp[2])
		tempL=subset(tempL,tempL%in%c(1:ntips))
		
		vcv[tempL,tempR]=n.covsums[i+ntips]	
		vcv[tempR,tempL]=n.covsums[i+ntips]	
	}

	vcv
	
}

#Utility function for obtaining mrcas for each species pair
Ancestors<-function (x, node, type = c("all", "parent")) 
{
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

##NOTES:
#A tree of 8500 species took 4.5 minutes to run to completion; required R64 and 10G of RAM -- JMB 12-29-10