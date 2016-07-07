# perform a robust estimation of the likelihood fit by non-parametric bootstrap
perform.consensus.likelihood.fit <- function( bootstrap_results, dataset, regularization, command = "hc" ) {
    
    # estimate the poset based on bootstrap
    estimated.poset = build.consensus.poset(bootstrap_results,rep(nrow(bootstrap_results[[1]]),nrow(bootstrap_results[[1]])))
    colnames(dataset) = as.character(1:ncol(dataset))
    colnames(estimated.poset) = colnames(dataset)
    rownames(estimated.poset) = colnames(dataset)
    
    # perform the inference on the constrained poset
    adj.matrix = constrained.likelihood.fit(dataset,estimated.poset,regularization,command)
    
    return(adj.matrix)
    
}

# perform a robust estimation of the poset by likelihood fit by non-parametric bootstrap
build.consensus.poset <- function( bootstrap_results, parents_cardinalities ) {
    
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
    adj.matrix = remove.loops(adj.matrix,added_edges,added_scores)
    colnames(adj.matrix) = as.character(1:ncol(adj.matrix))
    rownames(adj.matrix) = as.character(1:nrow(adj.matrix))
    
    return(adj.matrix)
    
}

# remove any loop in the adjacency matrix
remove.loops <- function( adj.matrix, edges, scores ) {
    
    # if I have at least one edge
    if(length(edges)>0) {
            
        ordered.scores = sort(scores,decreasing=FALSE,index.return=TRUE)
        ordered.edges = edges[ordered.scores$ix]
        
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

# perform structural inference by likelihood fit on the constrained poset
constrained.likelihood.fit <- function( dataset, poset.adj.matrix, regularization, command = "hc" ) {

    # structure to save the adjacency matrix of the reconstructed topolgy
    adj.matrix = array(0, c(ncol(dataset), ncol(dataset)))
    rownames(adj.matrix) = colnames(dataset)
    colnames(adj.matrix) = colnames(dataset)

    # set up the data structure on which to perform the likelihood fit
    if(regularization%in%c("loglik","aic","bic")) {
        
        # discrete case, create a categorical data frame from the dataset
        data = array("missing",c(nrow(dataset),ncol(dataset)))
        for (i in 1:nrow(dataset)) {
            for (j in 1:ncol(dataset)) {
                if(dataset[i,j]==1) {
                    data[i,j] = "observed"
                }
            }
        }
        data = as.data.frame(data, stringsAsFactors = TRUE)
        for (n in names(data)) {
            levels(data[[n]]) = c('missing', 'observed')
        }

        # renaming
        colnames(data) = as.character(1:ncol(dataset))
        rownames(data) = as.character(1:nrow(dataset))
        
    }
    else if(regularization%in%c("loglik-g","aic-g","bic-g")) {
        
        # continuous case
        data = dataset

        # renaming
        colnames(data) = as.character(1:ncol(dataset))
        rownames(data) = as.character(1:nrow(dataset))
        
    }
    
    # create the blacklist based on the given poset
    cont = 0
    parent = -1
    child = -1

    # set the blacklisted nodes
    for (i in rownames(poset.adj.matrix)) {
        for (j in colnames(poset.adj.matrix)) {
            if(i != j) {
                if (poset.adj.matrix[i,j] == 0) {
                    # [i,j] refers to arc i --> j
                    cont = cont + 1
                    if (cont == 1) {
                        parent = i
                        child = j
                    }
                    else {
                        parent = c(parent, i)
                        child = c(child, j)
                    }
                }
            }
        }
    }

    # perform the reconstruction by likelihood fit
    if (cont > 0) {
        blacklist = data.frame(from = parent,to = child)
        if (command == "hc") {
            my.net = hc(data,score = regularization, blacklist = blacklist)
        }
        else if (command == "tabu") {
            my.net = tabu(data,score = regularization, blacklist = blacklist)
        }
    }
    else {
        if (command == "hc") {
            my.net = hc(data, score = regularization)
        }
        else if (command == "tabu") {
            my.net = tabu(data, score = regularization)
        }
    }
    my.arcs = my.net$arcs
    
    # create the adjacency matrix of the reconstructed topology
    if(length(nrow(my.arcs))>0 && nrow(my.arcs)>0) {
        for (i in 1:nrow(my.arcs)) {
            # [i,j] refers to the edge i --> j
            adj.matrix[as.numeric(my.arcs[i,1]),as.numeric(my.arcs[i,2])] = 1
        }
    }
    
    return(adj.matrix)
    
}


#### end of file -- consensus.likelihood.fit.R
