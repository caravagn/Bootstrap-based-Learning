# perform a robust estimation of the likelihood fit by non-parametric bootstrap
build.consensus.poset.skeleton <- function( bootstrap_results, dataset, regularization, command = "hc", algorithm = "mmpc" ) {
    
    # estimate the poset based on bootstrap
    estimated.poset = build.consensus(bootstrap_results,rep(nrow(bootstrap_results[[1]]),nrow(bootstrap_results[[1]])))
    colnames(dataset) = as.character(1:ncol(dataset))
    colnames(estimated.poset) = colnames(dataset)
    rownames(estimated.poset) = colnames(dataset)
    
    # perform the inference of the skeleton
    skeleton = skeleton.learning(dataset,algorithm,regularization)
    
    # perform the structural inference
    adj.matrix = build.constrained.skeleton(estimated.poset,skeleton)
    
    return(adj.matrix)
    
}

# perform structural inference by likelihood fit
skeleton.learning <- function( dataset, algorithm = "mmpc", regularization ) {

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
    
    if(algorithm =="mmpc") {
        my.skeleton = mmpc(data)
    }
    my.arcs = my.skeleton$arcs
    
    # create the adjacency matrix of the reconstructed topology
    if(length(nrow(my.arcs))>0 && nrow(my.arcs)>0) {
        for (i in 1:nrow(my.arcs)) {
            # [i,j] refers to the edge i --> j
            adj.matrix[as.numeric(my.arcs[i,1]),as.numeric(my.arcs[i,2])] = 1
            adj.matrix[as.numeric(my.arcs[i,2]),as.numeric(my.arcs[i,1])] = 1
        }
    }
    
    return(adj.matrix)
    
}

# combine the poset and the skeleton
build.constrained.skeleton = function( poset, skeleton ) {

    # structure to save the adjacency matrix of the reconstructed topolgy
    adj.matrix = array(0, c(ncol(poset), ncol(poset)))
    rownames(adj.matrix) = colnames(poset)
    colnames(adj.matrix) = colnames(poset)
    
    # combine poset and skeleton
    for (i in 1:nrow(adj.matrix)) {
        for (j in 1:ncol(adj.matrix)) {
            if(poset[i,j] == 1 && skeleton[i,j]) {
                adj.matrix[i,j] = 1
            }
        }
    }
    
    return(adj.matrix)
    
}


#### end of file -- build.consensus.poset.skeleton.R
