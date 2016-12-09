build_agony_poset = function( agony_file, dataset, bootstrap_results ) {

	# read the results of agony computation
	my_ordering = read.table(file=agony_file,header=FALSE,sep=" ",check.names=FALSE,stringsAsFactors=FALSE)

	# estimate a poset for the DAG using agony (this can still contain cycles)
	agony_poset = get.poset.v2(my_ordering,ncol(dataset))

	# make the poset acyclic using bootstrap confidence
	agony_poset = perform.consensus.estimate(Reduce("+",bootstrap_results),agony_poset)
    colnames(dataset) = as.character(1:ncol(dataset))
    colnames(agony_poset) = colnames(dataset)
    rownames(agony_poset) = colnames(dataset)

    # return the poset
    return(agony_poset)

}

get.poset.v2 = function( ordering, node_size ) {
    poset = array(0,c(node_size,node_size))
     for(i in min(ordering[,2]):max(ordering[,2])) {
         curr_parents = ordering[which(ordering[,2]==i),1]
         curr_childred = ordering[which(ordering[,2]>i),1]
         # add arcs from lower to higher ranked nodes
         if(length(curr_parents)>0 && length(curr_childred)>0) {
             for(j in curr_parents) {
                 poset[j,curr_childred] = 1
             }
         }
         # add arcs in both ways for equal ranked nodes (NOTE: with this, the poset can be cyclic!)
         if(length(curr_parents)>1) {
             for(a in 1:length(curr_parents)) {
                 for(b in a:length(curr_parents)) {
                     if(a!=b) {
                         poset[curr_parents[a],curr_parents[b]] = 1
                         poset[curr_parents[b],curr_parents[a]] = 1
                     }
                 }
             }
         }
     }
     return(poset)
}

# perform a robust estimation of the likelihood fit by non-parametric bootstrap
perform.consensus.estimate <- function( bootstrap_results, poset ) {

    # constrain the confidence by bootstrap based on the poset
    constrained_bootstrap = (bootstrap_results*poset)
    
    # estimate the poset based on bootstrap confidence
    estimated.poset.confidence = build.consensus.poset_agony(constrained_bootstrap,rep(nrow(constrained_bootstrap),nrow(constrained_bootstrap)),TRUE)
    
    return(estimated.poset.confidence)
    
}

# perform a robust estimation of the poset by likelihood fit by non-parametric bootstrap
build.consensus.poset_agony <- function( bootstrap_results, parents_cardinalities, confidence ) {
    
    # create the consensus graph
    adj.matrix = array(0,c(nrow(bootstrap_results),ncol(bootstrap_results)))
    added_edges = list()
    added_scores = NULL
    for (i in 1:ncol(bootstrap_results)) {
        if(parents_cardinalities[i]>0) {
            sorted_bootstrap = sort(bootstrap_results[,i],decreasing=TRUE)
            valid_parents = sorted_bootstrap[which(sorted_bootstrap[1:parents_cardinalities[i]]>0)]
            if(length(valid_parents)>0) {
                added_scores = c(added_scores,as.vector(valid_parents))
                for(j in 1:length(valid_parents)) {
                    added_edges[[(length(added_edges)+1)]] = list(parent=names(valid_parents)[j],child=i)
                }
                adj.matrix[as.numeric(names(valid_parents)),i] = 1
            }
        }
    }
    
    # remove any loop from the adjacency matrix
    adj.matrix = remove.loops_agony(adj.matrix,added_edges,added_scores,confidence)
    colnames(adj.matrix) = as.character(1:ncol(adj.matrix))
    rownames(adj.matrix) = as.character(1:nrow(adj.matrix))
    
    return(adj.matrix)
    
}

# remove any loop in the adjacency matrix
remove.loops_agony <- function( adj.matrix, edges, scores, confidence ) {
    
    # if I have at least one edge
    if(length(edges)>0) {
        
        if(confidence) {
            ordered.scores = sort(scores,decreasing=FALSE,index.return=TRUE)
            ordered.edges = edges[ordered.scores$ix]
        }
        else {
            ordered.pos = sample(1:length(scores))
            ordered.scores = scores[ordered.pos]
            ordered.edges = edges[ordered.pos]
        }
        
        # go through the edges in decreasing order of confidence
        for (i in 1:length(edges)) {
            
            # consider any edge i --> j
            curr.edge = ordered.edges[[i]]
            curr.edge.i = as.numeric(curr.edge$parent)
            curr.edge.j = as.numeric(curr.edge$child)
            
            # search for loops between curr.edge.i and curr.edge.j
            curr.graph = graph.adjacency(adj.matrix,mode="directed")
            is.path = suppressWarnings(get.shortest.paths(curr.graph,
                                       curr.edge.j,
                                       curr.edge.i)$vpath)
            is.path = length(unlist(is.path))

            # if there is a path between the two nodes, remove edge i --> j
            if (is.path > 0) {
                adj.matrix[curr.edge.i,curr.edge.j] = 0
            }
            
        }
        
    }
    
    return(adj.matrix)
    
}
