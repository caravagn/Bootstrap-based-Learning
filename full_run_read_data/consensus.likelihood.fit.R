# perform a robust estimation of the likelihood fit by non-parametric bootstrap
perform.consensus.likelihood.fit <- function( bootstrap_results, dataset, regularization, command = "hc" ) {
    
    # estimate the poset based on bootstrap confidence
    estimated.poset.confidence = build.consensus.poset(bootstrap_results,rep(nrow(bootstrap_results[[1]]),nrow(bootstrap_results[[1]])),TRUE)
    colnames(dataset) = as.character(1:ncol(dataset))
    colnames(estimated.poset.confidence) = colnames(dataset)
    rownames(estimated.poset.confidence) = colnames(dataset)
    
    return(estimated.poset.confidence)
    
}

# perform a robust estimation of the poset by likelihood fit by non-parametric bootstrap
build.consensus.poset <- function( bootstrap_results, parents_cardinalities, confidence ) {
    
    # estimate the confidence of each edge
    bootstrap_results = Reduce("+",bootstrap_results)
    
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
    adj.matrix = remove.loops(adj.matrix,added_edges,added_scores,confidence)
    colnames(adj.matrix) = as.character(1:ncol(adj.matrix))
    rownames(adj.matrix) = as.character(1:nrow(adj.matrix))
    
    return(adj.matrix)
    
}

# remove any loop in the adjacency matrix
remove.loops <- function( adj.matrix, edges, scores, confidence ) {
    
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


#### end of file -- consensus.likelihood.fit.R
