#Weight matrix generator taken from Butler and King (2004) and modified to allow multiple alpha parameters



weight.mat.network<-function(phy, edges, Rate.mat, root.state, simmap.tree=FALSE, scaleHeight=FALSE, assume.station=TRUE){
	W <- weight.mat(phy, edges, Rate.mat, root.state, simmap.tree=simmap.tree, scaleHeight=scaleHeight, assume.station=assume.station, standardizeRowSums=FALSE)
	#This gives us weight matrix up main path. Now have to adjust for the hybrid taxa
	
}

##Matrix generating function taken from vcv.phylo in ape:
mat.gen<-function(phy,piece.wise,pp){
	phy <- reorder(phy, "pruningwise")
    n <- length(phy$tip.label)
    anc <- phy$edge[,1]
    des <- phy$edge[,2]
    ep <- piece.wise[,1]
	comp <- numeric(n + phy$Nnode)
    mat <- matrix(0, n, n)
	
	for (i in length(anc):1) {
        focal <- comp[anc[i]]
        comp[des[i]] <- focal + ep[des[i]]
        j <- i - 1L
        while (anc[j] == anc[i] && j > 0) {
            left <- if (des[j] > n) pp[[des[j] - n]] else des[j]
            right <- if (des[i] > n) pp[[des[i] - n]] else des[i]
            mat[left, right] <- mat[right, left] <- focal
            j <- j - 1L
        }
    }
    diag.elts <- 1 + 0:(n - 1)*(n + 1)
    mat[diag.elts] <- comp[1:n]
	
	mat
}


