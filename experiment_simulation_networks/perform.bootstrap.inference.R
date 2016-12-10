# perform the boostrap based maximum likelihood structure inference
perform.bootstrap.inference <- function( dataset, regularization, nboot.first, agony.files, nboot.second ) {

    # structure where to save the results
    results = list()

    # perform the first pass of bootstrap
    cat("Performing the first step of bootstrap estimation...","\n")
    bootstrap.first.pass = bootstrap.estimation.first.pass(dataset=dataset,regularization=regularization,nboot=nboot.first)
    results[["bootstrap.first.pass"]] = bootstrap.first.pass

    # estimate the poset by bootstrap confidence
    cat("Performing the estimation of the confidence based poset...","\n")
    confidence.poset = estimate.confidence.poset(bootstrap_results=bootstrap.first.pass,dataset=dataset)
    results[["confidence.poset"]] = confidence.poset
    
    # estimate the poset by agony and bootstrap confidence
    cat("Performing the estimation of the agony based poset...","\n")
    get.agony.edges.list(bootstrap.first.pass,paste0(agony.files,"/agony_inputs.txt"))
    system(paste0("./agony ",agony.files,"/agony_inputs.txt ",agony.files,"/agony_outputs.txt"))
    agony.poset = build.agony.poset(paste0(agony.files,"/agony_outputs.txt"),dataset,bootstrap.first.pass)
    results[["agony.poset"]] = agony.poset
    
    # perform the second pass of bootstrap on both posets
    cat("Performing the second step of bootstrap estimation...","\n")
    bootstrap.second.pass = list()
    bootstrap.second.pass[["confidence.poset"]] = bootstrap.estimation.second.pass(dataset=dataset,regularization=regularization,nboot=nboot.second,contraint.poset=confidence.poset)
    bootstrap.second.pass[["agony.poset"]] = bootstrap.estimation.second.pass(dataset=dataset,regularization=regularization,nboot=nboot.second,contraint.poset=agony.poset)
    results[["bootstrap.second.pass"]] = bootstrap.second.pass

    return(results)

}

# perform a robust estimation of the likelihood fit by non-parametric bootstrap
bootstrap.estimation.first.pass <- function( dataset, regularization, command = "hc",
                                    nboot = 100, random.seed = NULL, cores.ratio = 1,
                                    verbose = FALSE ) {
    
    # structure to save the results of the bootstrap
    bootstrap.results = NULL
    
    # set the seed to be used for the random samplings
    set.seed(random.seed)

    # set the number of cores to be used in the parallelization of the bootstrap
    cores = as.integer(cores.ratio * (detectCores() - 1))
    if (cores < 1) {
        cores = 1
    }
    
    # setup the parallelization to perform the bootstrap
    cl = makeCluster(cores)
    registerDoParallel(cl)
    
    if (verbose) {
        cat('\tPerforming bootstrap with', cores, 'cores via "parallel." \n')
    }
    
    # perform nboot bootstrap resampling
    r = foreach(num = 1:nboot, .packages='bnlearn', .export='perform.likelihood.fit') %dopar% {
        
        # create the sampled dataset for the current iteration
        samples = sample(1:nrow(dataset),size=nrow(dataset),replace=TRUE)
        resampled.dataset = dataset[samples,]
        
        # perform the likelihood fit on the bootstrapped dataset
        reconstructed.adj.matrix = perform.likelihood.fit(resampled.dataset,regularization,command)
        
    }

    stopCluster(cl)
    
    if (verbose) {
        cat("\tBootstrap completed.\n")
    }
    
    bootstrap.results = r
    
    return(bootstrap.results)
    
}

# perform structural inference by maximum likelihood
perform.likelihood.fit <- function( dataset, regularization, command = "hc" ) {

    # structure to save the adjacency matrix of the reconstructed topolgy
    adj.matrix = array(0, c(ncol(dataset), ncol(dataset)))
    rownames(adj.matrix) = colnames(dataset)
    colnames(adj.matrix) = colnames(dataset)

    # perform maximum likelihood estimation either by hill climbing or tabu search
    if(command=="hc") {
        my.net = hc(dataset,score=regularization)
    }
    else if(command=="tabu") {
        my.net = tabu(dataset,score=regularization)
    }
    my.arcs = my.net$arcs
    
    # create the adjacency matrix of the reconstructed topology
    if(length(nrow(my.arcs))>0 && nrow(my.arcs)>0) {
        for (i in 1:nrow(my.arcs)) {
            # [i,j] refers to the edge i --> j
            adj.matrix[which(colnames(dataset)%in%as.character(my.arcs[i,1])),which(colnames(dataset)%in%as.character(my.arcs[i,2]))] = 1
        }
    }
    
    return(adj.matrix)
    
}

# perform a robust estimation of the likelihood fit by non-parametric bootstrap
estimate.confidence.poset <- function( bootstrap_results, dataset ) {
    
    # estimate the poset based on bootstrap confidence
    estimated.poset.confidence = build.confidence.poset(bootstrap_results,rep(nrow(bootstrap_results[[1]]),nrow(bootstrap_results[[1]])))
    colnames(estimated.poset.confidence) = colnames(dataset)
    rownames(estimated.poset.confidence) = colnames(dataset)
    
    return(estimated.poset.confidence)
    
}

# perform a robust estimation of the poset by likelihood fit and non-parametric bootstrap
build.confidence.poset <- function( bootstrap_results, parents_cardinalities ) {
    
    # estimate the confidence of each edge
    bootstrap_results = Reduce("+",bootstrap_results)
    
    # create the consensus graph
    adj.matrix = array(0,c(nrow(bootstrap_results),ncol(bootstrap_results)))
    added_edges = list()
    added_scores = NULL
    for (i in 1:ncol(bootstrap_results)) {
        if(parents_cardinalities[i]>0) {
            # consider the best parents_cardinalities arcs greater than 0
            sorted_bootstrap = sort(bootstrap_results[,i],decreasing=TRUE)
            valid_parents = sorted_bootstrap[which(sorted_bootstrap[1:parents_cardinalities[i]]>0)]
            if(length(valid_parents)>0) {
                added_scores = c(added_scores,as.vector(valid_parents))
                added_valid_parents = NULL
                for(j in 1:length(valid_parents)) {
                    new_valid_parent_pos = which(colnames(bootstrap_results)==names(valid_parents)[j])
                    added_edges[[(length(added_edges)+1)]] = list(parent=new_valid_parent_pos,child=i)
                    added_valid_parents = c(added_valid_parents,new_valid_parent_pos)
                }
                adj.matrix[added_valid_parents,i] = 1
            }
        }
    }
    
    # remove any loop from the adjacency matrix
    adj.matrix = remove.loops(adj.matrix,added_edges,added_scores)
    
    return(adj.matrix)
    
}

# remove any loop in the given adjacency matrix
remove.loops <- function( adj.matrix, edges, scores ) {
    
    # if I have at least one edge
    if(length(edges)>0) {
        
        # compute the ordered confidence scores
        duplicated_scores = unique(scores[duplicated(scores)])
        ordered.scores = sort(scores,decreasing=FALSE,index.return=TRUE)
        ordered.scores_values = ordered.scores$x
        ordered.scores_idx = ordered.scores$ix
        ordered.scores_idx_random = ordered.scores_idx
        
        # randomly order arcs with equal scores to avoid any bias
        for(my_dup_score in duplicated_scores) {
            curr_dup_scores = which(ordered.scores_values==my_dup_score)
            curr_dup_scores_random = sample(curr_dup_scores,size=length(curr_dup_scores),replace=FALSE)
            ordered.scores_idx_random[curr_dup_scores] = ordered.scores_idx_random[curr_dup_scores_random]
        }
        ordered.edges = edges[ordered.scores_idx_random]
        
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

# set the files for agony computation
get.agony.edges.list <- function( bootstrap_results, agony_file ) {
    for(i in 1:length(bootstrap_results)) {
        curr_adj.matrix = bootstrap_results[[i]]
        for(j in 1:nrow(curr_adj.matrix)) {
            for(k in 1:ncol(curr_adj.matrix)) {
                if(curr_adj.matrix[j,k]==1) {
                    cat(as.character(j),as.character(k),"\n",file=agony_file,append=TRUE)
                }
            }
        }
    }
}

# build the agony based poset
build.agony.poset <- function( agony_file, dataset, bootstrap_results ) {

    # read the results from agony computation
    my_ordering = read.table(file=agony_file,header=FALSE,sep=" ",check.names=FALSE,stringsAsFactors=FALSE)

    # estimate a poset for the DAG using agony. 
    # Note: this can still contain cycles, in fact elements at the same ordering, have arcs both ways
    agony_poset = get.agony.poset(my_ordering,ncol(dataset))

    # make the poset acyclic using bootstrap confidence
    agony_poset = perform.agony.estimate(Reduce("+",bootstrap_results),agony_poset)
    colnames(agony_poset) = colnames(dataset)
    rownames(agony_poset) = colnames(dataset)

    # return the poset
    return(agony_poset)

}

get.agony.poset <- function( ordering, node_size ) {
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

# perform a robust estimation of the poset by likelihood fit and non-parametric bootstrap
perform.agony.estimate <- function( bootstrap_results, poset ) {

    # constrain the confidence by bootstrap based on the poset
    constrained_bootstrap = (bootstrap_results*poset)
    
    # estimate the poset based on bootstrap confidence
    estimated.poset.agony = build.consensus.poset_agony(constrained_bootstrap,rep(nrow(constrained_bootstrap),nrow(constrained_bootstrap)))
    
    return(estimated.poset.agony)
    
}

# perform a robust estimation of the poset by likelihood fit by non-parametric bootstrap
build.consensus.poset_agony <- function( bootstrap_results, parents_cardinalities ) {
    
    # create the consensus graph
    adj.matrix = array(0,c(nrow(bootstrap_results),ncol(bootstrap_results)))
    added_edges = list()
    added_scores = NULL
    for (i in 1:ncol(bootstrap_results)) {
        if(parents_cardinalities[i]>0) {
            # consider the best parents_cardinalities arcs greater than 0
            sorted_bootstrap = sort(bootstrap_results[,i],decreasing=TRUE)
            valid_parents = sorted_bootstrap[which(sorted_bootstrap[1:parents_cardinalities[i]]>0)]
            if(length(valid_parents)>0) {
                added_scores = c(added_scores,as.vector(valid_parents))
                added_valid_parents = NULL
                for(j in 1:length(valid_parents)) {
                    new_valid_parent_pos = which(colnames(bootstrap_results)==names(valid_parents)[j])
                    added_edges[[(length(added_edges)+1)]] = list(parent=new_valid_parent_pos,child=i)
                    added_valid_parents = c(added_valid_parents,new_valid_parent_pos)
                }
                adj.matrix[added_valid_parents,i] = 1
            }
        }
    }
    
    # remove any loop from the adjacency matrix
    adj.matrix = remove.loops(adj.matrix,added_edges,added_scores)
    
    return(adj.matrix)
    
}



# perform a robust estimation of the likelihood fit by non-parametric bootstrap for a given poset both for model and random
bootstrap.estimation.second.pass <- function( dataset, regularization, command = "hc",
                                    nboot = 100, random.seed = NULL, cores.ratio = 1,
                                    verbose = FALSE, contraint.poset ) {
    
    # structure to save the results of the bootstrap
    bootstrap.results = NULL
    
    # set the seed to be used for the random samplings
    set.seed(random.seed)

    # set the number of cores to be used in the parallelization of the bootstrap
    cores = as.integer(cores.ratio * (detectCores() - 1))
    if (cores < 1) {
        cores = 1
    }
    
    # setup the parallelization to perform the bootstrap
    cl = makeCluster(cores)
    registerDoParallel(cl)
    
    if (verbose) {
        cat('\tPerforming bootstrap with', cores, 'cores via "parallel." \n')
    }
    
    # perform nboot bootstrap resampling
    r = foreach(num = 1:nboot, .packages='bnlearn', .export='perform.constrained.likelihood.fit') %dopar% {
        
        # create the sampled datasets for the current iteration
        samples = sample(1:nrow(dataset),size=nrow(dataset),replace=TRUE)
        resampled.dataset = dataset[samples,]
        resampled.dataset_null = resampled.dataset
        for(null_model in 1:ncol(resampled.dataset_null)) {
            resampled.dataset_null[,null_model] = resampled.dataset_null[sample(1:nrow(resampled.dataset_null),size=nrow(resampled.dataset_null),replace=FALSE),null_model]
        }
        
        # perform the likelihood fit on the bootstrapped dataset
        reconstructed.adj.matrix = list()
        reconstructed.adj.matrix[["model"]] = perform.constrained.likelihood.fit(resampled.dataset,contraint.poset,regularization,command)
        reconstructed.adj.matrix[["null"]] = perform.constrained.likelihood.fit(resampled.dataset_null,contraint.poset,regularization,command)
        reconstructed.adj.matrix
        
    }

    stopCluster(cl)
    
    if (verbose) {
        cat("\tBootstrap completed.\n")
    }
    
    bootstrap.results = r
    
    bootstrap.results_model = list()
    bootstrap.results_null = list()
    for(i in 1:length(bootstrap.results)) {
        bootstrap.results_model[[paste0("iteration_",as.character(i))]] = bootstrap.results[[i]][["model"]]
        bootstrap.results_null[[paste0("iteration_",as.character(i))]] = bootstrap.results[[i]][["null"]]
    }
    bootstrap.results = list()
    bootstrap.results[["model"]] = bootstrap.results_model
    bootstrap.results[["null"]] = bootstrap.results_null
    
    return(bootstrap.results)
    
}

# perform structural inference by maximum likelihood
perform.constrained.likelihood.fit <- function( dataset, poset.adj.matrix, regularization, command = "hc" ) {

    # structure to save the adjacency matrix of the reconstructed topolgy
    adj.matrix = array(0, c(ncol(dataset), ncol(dataset)))
    rownames(adj.matrix) = colnames(dataset)
    colnames(adj.matrix) = colnames(dataset)
    
    # create the blacklist based on the given poset
    cont = 0
    parent = NULL
    child = NULL

    # set the blacklisted nodes
    colnames(poset.adj.matrix) = colnames(dataset)
    rownames(poset.adj.matrix) = colnames(dataset)
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

    # perform maximum likelihood estimation either by hill climbing or tabu search
    if (cont > 0) {
        blacklist = data.frame(from = parent,to = child)
        if (command == "hc") {
            my.net = hc(dataset,score=regularization, blacklist = blacklist)
        }
        else if (command == "tabu") {
            my.net = tabu(dataset,score=regularization, blacklist = blacklist)
        }
    }
    else {
        if (command == "hc") {
            my.net = hc(dataset,score=regularization)
        }
        else if (command == "tabu") {
            my.net = tabu(dataset,score=regularization)
        }
    }
    my.arcs = my.net$arcs
    
    # create the adjacency matrix of the reconstructed topology
    if(length(nrow(my.arcs))>0 && nrow(my.arcs)>0) {
        for (i in 1:nrow(my.arcs)) {
            # [i,j] refers to the edge i --> j
            adj.matrix[which(colnames(dataset)%in%as.character(my.arcs[i,1])),which(colnames(dataset)%in%as.character(my.arcs[i,2]))] = 1
        }
    }
    
    return(adj.matrix)
    
}
