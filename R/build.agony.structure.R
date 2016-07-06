# perform a robust estimation of the likelihood fit by non-parametric bootstrap
build.agony.structure <- function( dataset, estimated.poset, regularization, command = "hc" ) {
    
    # estimate the poset based on bootstrap
    colnames(dataset) = as.character(1:ncol(dataset))
    colnames(estimated.poset) = colnames(dataset)
    rownames(estimated.poset) = colnames(dataset)
    
    # perform the inference
    adj.matrix = constrained.likelihood.fit(dataset,estimated.poset,regularization,command)
    
    return(adj.matrix)
    
}

# perform structural inference by likelihood fit
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

    for (i in rownames(poset.adj.matrix)) {
        for (j in colnames(poset.adj.matrix)) {
            if(i != j) {
                if (poset.adj.matrix[i,j] == 0) {
                    # [i,j] refers to arc i --> j
                    cont = cont + 1
                    if (cont == 1) {
                        parent = i
                        child = j
                    } else {
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
        } else if (command == "tabu") {
            my.net = tabu(data,score = regularization, blacklist = blacklist)
        }
    } else {
        if (command == "hc") {
            my.net = hc(data, score = regularization)
        } else if (command == "tabu") {
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


#### end of file -- build.consensus.poset.R
