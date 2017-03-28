## January 26, 2017
## Authors: Daniel S. Caetano and Luke J. Harmon
## Supplement for the manuscript "Estimating correlated rates of trait evolution with uncertainty"

## Here we provide the implementation of the pruning algorithm as described in the supplementary material.
## The code was made to be readable rather than efficient, so it can be used as an example for future implementations.

###########################################
## Load dependencies                     ##
###########################################

library( phytools )

## Here we use 'phytools' to map the regime to the phylogenetic tree. Then, we use the element 'phy$mapped.edge' to get the information of where each rate is distributed on the branches of the tree.

###########################################
## The pruning algorithm functions       ##
###########################################

## Bellow are functions used to compute the pruning algorithm. Note that here we are assuming just two regimes applied to the tree. But the algorithm can be applied to any number of rate regimes.

logLikPruning <- function(data, phy, root, R){
    ## 'data' is a matrix with traits in the columns and trait values for the species on the rows. Matrix need to have row.names equal to the tip names of the phylogeny.
    ## 'phy' is a phylogenetic tree of the class 'simmap' with mapped regimes.
    ## 'root' is a numeric vector with length equal to the number of columns of 'data' with values for the root.
    ## 'R' is a list of matrices. The number of elements in the list is equal to the number of regimes fitted to the tree. The dimension of each matrix is equal to the number of columns in the 'data' matrix.
    
    data <-  as.matrix(data) ## 'data' needs to be a class matrix.
    root <- as.numeric(root) ## 'root' needs to be numeric vector.
    k <- ncol(data) ## Number of traits.

    ## We need to find the ancestral and the descendent nodes of the tree, so that we can traverse the tree using the pruning algorithm. So we use 'phy$edge' to get this information. This will produce a object 'anc' with the ancestral edges and another 'des' with the descendant edges. Then, the object 'nodes' have the internal nodes of the tree that we will traverse.

    ord.id <- reorder.phylo(phy, order="postorder", index.only=TRUE) ## Order for traversal.
    anc <- phy$edge[ord.id,1] ## Ancestral edges.
    des <- phy$edge[ord.id,2] ## Descendent edges.
    nodes <- unique(anc) ## The internal nodes we will traverse.

    ## We also need to reorder the 'phy$mapped.edge' matrix, so we can use the information of the mapped regime on each of the edges of the tree.
    mapped.edge <- phy$mapped.edge[ord.id,] ## The regimes.

    ## Some computations will depend on whether: (1) the node descendants are tips, (2) the node descendants are other internal nodes or (3) the node descendant is a tip and an internal node. 
    node.to.tip <- which( tabulate( anc[which(des <= length(phy$tip.label))] ) == 2 ) ## type (1).
    node.to.node <- which( tabulate( anc[which(des > length(phy$tip.label))] ) == 2 ) ## type (2).
    node.to.tip.node <- unique( anc )[!unique( anc ) %in% c(node.to.node, node.to.tip)] ## type (3).

    ## Set the node types for the different ancestral nodes:
    names(anc) <- rep(1, times=length(anc))
    names(anc)[which(anc %in% node.to.node)] <- 2
    names(anc)[which(anc %in% node.to.tip.node)] <- 3

    ## Create vectors to store results.
    X0 <- matrix(nrow=k, ncol=length(nodes)+1) ## Stores the new trait vector for the node.
    V0 <- array(dim=c(k, k, length(nodes)+1)) ## Stores the variance of the estimate of the new trait vector for the node.
    key <- vector(mode="numeric") ## Keep track of the traversal, so we know the nodes we already visited.
    ll <- vector(mode="numeric") ## Stores the node-by-node log-likelihood. The final log-likelihood is the sum of them.

    ## Now we are ready to traverse the tree. For each node we compute the contrast value of the descendents, estimate the new trait value for the node, the variance of this estimate, and, finally, the log-likelihood for the node. Some of these opperations are dependent of the type of the node (1, 2, or 3). To deal with that we define a list of functions, one for each node type, then call the function associated with the node type we are visiting. Note that we recorded this information previously.
    node_calc <- list(calcNodeToTip, calcNodeToNode, calcNodeToTipNode)

    ## Traverse the tree with a loop of 'nodes'.
    for (i in nodes) { ## Will visit all the internal nodes.
        node.id <- which(anc == i) ## The index for the 'des', 'anc', and 'mapped.edge (lines)'.
        type <- as.numeric( names(anc[node.id[1]]) )
        cache <- node_calc[[type]](data, k, des, anc, mapped.edge, R, node.id, X0, V0, key, ll)
        ## The output 'cache' is a list with the updated objects 'X0', 'V0', 'key' and 'll' that store the results for the computation on each node.
        X0 <- cache$X0
        V0 <- cache$V0
        key <- cache$key
        ll <- cache$ll
    }

    ## Append the log-likelihood for the root value.
    last.node <- length(key) ## The index for the root node.
    ll <- c(ll, logLikNode(X0[,last.node] - root, V0[,,last.node], solve(V0[,,last.node]), k))

    ## Return the sum of the log-likelihood computed over all the nodes.
    return(sum(ll))
}

calcNodeToTip <- function(data, k, des, anc, mapped.edge, R, node.id, X0, V0, key, ll){
    ## Make calculations with nodes that have only tip descendents.
    des.node <- des[ node.id ] ## The descendents of this node.
    ss <- data[des.node[1],] - data[des.node[2],] ## The contrast.
    ## The next two lines will rescale the rate matrix to each of the two branch lengths (descendents from the node). Since we have two regimes, the scaling factor can be 0 (meaning the rate regime is not present in this branch), in this case only the other rate matrix will be relevant for the computations.
    Rs1 <- ( R[[1]] * mapped.edge[node.id[1],1] ) + ( R[[2]] * mapped.edge[node.id[1],2] )
    Rs2 <- ( R[[1]] * mapped.edge[node.id[2],1] ) + ( R[[2]] * mapped.edge[node.id[2],2] )
    Rinv <- chol2inv(chol(Rs1+Rs2)) ## This is the same as: Rinv <- solve(Rs1+Rs2) , but faster.
    ## Here 'logLikNode' is a small function to compute the log-likelihood.
    ll <- c(ll, logLikNode(ss, Rs1+Rs2, Rinv, k))
    key <- c(key, anc[node.id[1]]) ## 'key' for both X0 and V0. This is the node number.
    ## X0 is the new trait value that will substitute the node as a result of the 'pruning' part.
    X0[,length(key)] <- ((Rs2 %*% Rinv) %*% data[des.node[1],]) + ((Rs1 %*% Rinv) %*% data[des.node[2],])
    ## V0 is the variance associated with the estimate of X0.
    V0[,,length(key)] <- chol2inv(chol( chol2inv(chol(Rs1)) + chol2inv(chol(Rs2)) ))
    ## Return the values so we can move forward in the pruning algorithm.
    return( list("X0"=X0, "V0"=V0, "key"=key, "ll"=ll) )
}

calcNodeToNode <- function(data, k, des, anc, mapped.edge, R, node.id, X0, V0, key, ll){
    ## Make calculations with nodes that have only node descendents.
    ## Some steps here are the same from 'calcNodeToTip'. So we are not repeating the comments.
    des.node <- des[ node.id ]
    
    ## The next line finds the position in the vector 'key' for the descendents nodes for this node. We already visited the two descendents nodes of the present node, but we need to find which are they among the other nodes visited in the tree. Check the help for 'ape:::reorder.phylo' for help about the 'postorder' scheme.
    key.id <- sapply(des.node, function(x) which(key == x) )
    
    ss <- X0[,key.id[1]] - X0[,key.id[2]] ## The contrast for the nodes.
    Rs1 <- ( R[[1]] * mapped.edge[node.id[1],1] ) + ( R[[2]] * mapped.edge[node.id[1],2] )
    Rs2 <- ( R[[1]] * mapped.edge[node.id[2],1] ) + ( R[[2]] * mapped.edge[node.id[2],2] )

    ## When computing the variance of a node which descendents are also nodes we need to add the variance of the trait values from the nodes above. So, as we go down the tree node by node, the variance accumulates.
    Rs1 <- Rs1 + V0[,,key.id[1]] ## 'key.id[1]' is one of the descendents nodes.
    Rs2 <- Rs2 + V0[,,key.id[2]] ## 'key.id[2]' is the other descendent node.
    
    Rinv <- chol2inv(chol(Rs1+Rs2)) ## Same as: Rinv <- solve(Rs1+Rs2)
    ll <- c(ll, logLikNode(ss, Rs1+Rs2, Rinv, k))
    key <- c(key, anc[node.id[1]])
    X0[,length(key)] <- ((Rs2 %*% Rinv) %*% X0[,key.id[1]]) + ((Rs1 %*% Rinv) %*% X0[,key.id[2]])
    V0[,,length(key)] <- chol2inv(chol( chol2inv(chol(Rs1)) + chol2inv(chol(Rs2)) ))
    return( list("X0"=X0, "V0"=V0, "key"=key, "ll"=ll) )
}

calcNodeToTipNode <- function(data, k, des, anc, mapped.edge, R, node.id, X0, V0, key, ll){
    ## Make calculations with nodes that have both node and tip descendents.
    ## Some steps here are the same from 'calcNodeToTip'. So we are not repeating the comments.
    des.node <- des[ node.id ]

    ## Check which of the descendants is a node.
    nn <- which(des.node > (min(anc)-1) )
    nd <- des.node[nn] ## This is the node.
    nd.id <- node.id[nn] ## This is the 'id' position of the node.
    key.id <- which(key == nd) ## This is the 'key.id' for the node.

    ## Check which of the descendants is a tip.
    ## We do not need the 'key.id' here, because we do not need to recover quantities calculated on the tips, only on the nodes.
    tt <- which(des.node <= (min(anc)-1) )
    tip <- des.node[tt] ## This is a tip.
    tip.id <- node.id[tt] ## This is the id position for this tip.
    
    ss <- data[tip,] - X0[,key.id] ## The contrast for the nodes.
    Rs1 <- ( R[[1]] * mapped.edge[tip.id,1] ) + ( R[[2]] * mapped.edge[tip.id,2] )
    Rs2 <- ( R[[1]] * mapped.edge[nd.id,1] ) + ( R[[2]] * mapped.edge[nd.id,2] )

    ## Now only the branch associated with a node need to have an additional variance.
    Rs2 <- Rs2 + V0[,,key.id] ## The additional variance for the node descendent.
    Rinv <- chol2inv(chol(Rs1+Rs2)) ## Same as: Rinv <- solve(Rs1+Rs2)
    ll <- c(ll, logLikNode(ss, Rs1+Rs2, Rinv, k))
    key <- c(key, anc[node.id[1]])
    X0[,length(key)] <- ((Rs2 %*% Rinv) %*% data[tip,]) + ((Rs1 %*% Rinv) %*% X0[,key.id])
    V0[,,length(key)] <- chol2inv(chol( chol2inv(chol(Rs1)) + chol2inv(chol(Rs2)) ))
    return( list("X0"=X0, "V0"=V0, "key"=key, "ll"=ll) )
}

logLikNode <- function(ss, sigma_len ,sigma_len_inv, k){
    ## Simple function to calculate the log-likelihood at an internal node.
    ## Here:
    ## 'ss' is the contrast value for the node.
    ## 'sigma_len' is Sigma * branch length, which is the rate matrices scaled by the branch lengths.
    ## 'sigma_len_inv' is the inverse of 'sigma_len'. Here just to avoid computing it again.
    ll <- -0.5 * ( k*log(2*pi) + log(det(sigma_len)) + (ss %*% sigma_len_inv %*% ss) )
    return(ll)
}

###########################################
## Log-likelihood using the full inverse ##
###########################################

## The following function calculates the log-likelihood for the multivariate rate by inverting the Kronecker product between the phylogenetic variance-covariance matrix and the evolutionary rate matrix.

logLikInverse <- function(data, phy, root, R){
    ## 'data' is a matrix with traits in the columns and trait values for the species on the rows. Matrix need to have row.names equal to the tip names of the phylogeny.
    ## 'phy' is a phylogenetic tree of the class 'simmap' with mapped regimes.
    ## 'root' is a numeric vector with length equal to the number of columns of 'data' with values for the root.
    ## 'R' is a list of matrices. The number of elements in the list is equal to the number of regimes fitted to the tree. The dimension of each matrix is equal to the number of columns in the 'data' matrix.
    
    p <- length(R) ## Number of regimes.
    k <- ncol(data) ## Number of traits.
    n <- length(phy$tip.label) ## Number of species.
    C.m <- multiC(phy) ## The phylogenetic variance-covariance matrices for each regime.
    D <- matrix(0, nrow=n*k, ncol=k) ## The design D, but empty.
    for(i in 1:k){
        D[((n*(i-1))+1):(n*i),i] <- 1 ## This loop fills up the D matrix.
    }
    V.m <- list() ## The Kronecker between C and R.
    for(i in 1:p){
        V.m[[i]] <- kronecker(R[[i]], C.m[[i]]) ## Loop for each rate regime.
    }
    V <- Reduce("+", V.m) ## Makes the sum of the Kronecker for each regime.
    ntot <- n*k ## Dimension of the Kronecker matrix.
    y <- c(data)
    Vinv <- chol2inv(chol(V)) ## The inverse of the Kronecker. This is computationally expensive.
    Vdet <- det(V) ## The determinat of the Kronecker. This is also computationally expensive.
    ## Now compute the log-likelihood.
    logl <- -0.5 * ( ( (y - as.vector(D %*% root)) %*% Vinv %*% (y - as.vector(D %*% root)) ) + ( n*k*log(2*pi) ) + log( Vdet ) )
    return(logl[1,1])
}

###########################################
## Example of application                ##
###########################################

## Simulate tree and data:
tree <- pbtree(n=100)
trait <- rep(c("A","B"), each = 50)
names(trait) <- tree$tip.label

## The phylogeny of class 'simmap'
map <- make.simmap(tree, x=trait, model="ARD")

## The rates and the root value.
R1 <- rbind(c(1,0.5), c(0.5,1))
R2 <- rbind(c(2,1), c(1,2))
ll <- list("A"=R1, "B"=R2)
root <- c(1,1)

## The data.
data <- sim.corrs(tree=map, vcv=ll, anc=root)

## Compute the likelihood using the pruning algorithm for multiple evolutionary rate matrices:
logLikPruning(data=data, phy=map, root=root, R=ll)

## Compute the likelihood using the full inverse. This is the same as the equation 1 of the main text.
logLikInverse(data=data, phy=map, root=root, R=ll)

## Same result. Good!
