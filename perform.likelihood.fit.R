# reconstruct the best topology by likelihood fit
# dataset: a valid dataset
# command: type of search, either hill climbing (hc) or tabu (tabu)
# regularization: regularization term to be used in the likelihood fit
# return: adj.matrix.fit: the adjacency matrix of the reconstructed topology
perform.my.fit = function( dataset, command, regularization ) {
    
    # load the bnlearn package
    require("bnlearn")
    
    # adjacency matrix of the topology reconstructed by likelihood fit
    adj.matrix.fit = array(0,c(ncol(dataset),ncol(dataset)))
    rownames(adj.matrix.fit) = colnames(dataset)
    colnames(adj.matrix.fit) = colnames(dataset)

    # set up the data structure on which to perform the likelihood fit
    if(regularization%in%c("loglik","aic","bic")) {
        for (i in 1:ncol(dataset)) {
            if (sum(dataset[, i]) == 0) {
                dataset[sample(1:nrow(dataset), size=1), i] = 1
            } else if (sum(dataset[, i]) == nrow(dataset)) {
                dataset[sample(1:nrow(dataset), size=1), i] = 0
            }
        }
        # create a categorical data frame from the dataset
        data = array("missing",c(nrow(dataset),ncol(dataset)))
        for (i in 1:nrow(dataset)) {
            for (j in 1:ncol(dataset)) {
                if(dataset[i,j]==1) {
                    data[i,j] = "observed"
                }
            }
        }
        data = as.data.frame(data)
        my.names = names(data)
        for (i in 1:length(my.names)) {
            my.names[i] = toString(i)
        }
        names(data) = my.names
    }
    else {
        data = dataset
        my.names = colnames(dataset)
        for (i in 1:length(my.names)) {
            my.names[i] = toString(i)
        }
        names(data) = my.names
    }
    
    if(command=="hc") {
        my.net = hc(data,score=regularization)
    }
    else if(command=="tabu") {
        my.net = tabu(data,score=regularization)
    }
    my.arcs = my.net$arcs
    
    # build the adjacency matrix of the reconstructed topology
    if(length(nrow(my.arcs))>0 && nrow(my.arcs)>0) {
        for (i in 1:nrow(my.arcs)) {
            # [i,j] refers to arc i --> j
            adj.matrix.fit[as.numeric(my.arcs[i,1]),as.numeric(my.arcs[i,2])] = 1
        }
    }
    
    return(adj.matrix.fit)

}
