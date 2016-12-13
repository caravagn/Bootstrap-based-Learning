# perform a robust estimation of the likelihood fit by non-parametric bootstrap
bootstrap.estimation_constrained <- function( dataset, regularization, poset.constrain = NA, 
                                command = "hc", nboot = 100, random.seed = NULL, cores.ratio = 1, 
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
    r = foreach(num = 1:nboot, .packages='bnlearn', .export='constrained.likelihood.fit') %dopar% {
        
        # create the sampled datasets for the current iteration
        samples = sample(1:nrow(dataset),size=nrow(dataset),replace=TRUE)
        resampled.dataset = dataset[samples,]
        rownames(resampled.dataset) = as.character(1:nrow(resampled.dataset))
        colnames(resampled.dataset) = as.character(1:ncol(resampled.dataset))
        resampled.dataset_null = resampled.dataset
        for(null_model in 1:ncol(resampled.dataset_null)) {
            resampled.dataset_null[,null_model] = resampled.dataset_null[sample(1:nrow(resampled.dataset_null),size=nrow(resampled.dataset_null),replace=FALSE),null_model]
        }
        rownames(resampled.dataset_null) = as.character(1:nrow(resampled.dataset_null))
        colnames(resampled.dataset_null) = as.character(1:ncol(resampled.dataset_null))
        
        # perform the likelihood fit on the bootstrapped dataset
        reconstructed.adj.matrix = list()
        reconstructed.adj.matrix[["model"]] = constrained.likelihood.fit(resampled.dataset,poset.constrain,regularization,command)
        reconstructed.adj.matrix[["null"]] = constrained.likelihood.fit(resampled.dataset_null,poset.constrain,regularization,command)
        reconstructed.adj.matrix
    }

    stopCluster(cl)
    
    if (verbose) {
        cat("\tBootstrap completed.\n")
    }
    
    bootstrap.results = r
    
    return(bootstrap.results)
    
}

# perform structural inference by likelihood fit on the constrained poset
constrained.likelihood.fit <- function( dataset, poset.adj.matrix, regularization, command = "hc" ) {

    # structure to save the adjacency matrix of the reconstructed topolgy
    adj.matrix = array(0, c(ncol(dataset), ncol(dataset)))
    rownames(adj.matrix) = colnames(dataset)
    colnames(adj.matrix) = colnames(dataset)

    # # set up the data structure on which to perform the likelihood fit
    # if(regularization%in%c("loglik","aic","bic")) {
        
    #     # discrete case, create a categorical data frame from the dataset
    #     data = array("missing",c(nrow(dataset),ncol(dataset)))
    #     for (i in 1:nrow(dataset)) {
    #         for (j in 1:ncol(dataset)) {
    #             if(dataset[i,j]==1) {
    #                 data[i,j] = "observed"
    #             }
    #         }
    #     }
    #     data = as.data.frame(data, stringsAsFactors = TRUE)
    #     for (n in names(data)) {
    #         levels(data[[n]]) = c('missing', 'observed')
    #     }

    #     # renaming
    #     colnames(data) = as.character(1:ncol(dataset))
    #     rownames(data) = as.character(1:nrow(dataset))
        
    # }
    # else if(regularization%in%c("loglik-g","aic-g","bic-g")) {
        
    #     # continuous case
    #     data = dataset

    #     # renaming
    #     colnames(data) = as.character(1:ncol(dataset))
    #     rownames(data) = as.character(1:nrow(dataset))
        
    # }


    # modify for discrete variables (not binary)
    data = as.data.frame(dataset, stringsAsFactors = TRUE)
    
    # renaming
    colnames(data) = as.character(1:ncol(dataset))
    rownames(data) = as.character(1:nrow(dataset))
    
    # create the blacklist based on the given poset
    cont = 0
    parent = -1
    child = -1

    # set the blacklisted nodes
    colnames(poset.adj.matrix) = as.character(1:ncol(poset.adj.matrix))
    rownames(poset.adj.matrix) = as.character(1:nrow(poset.adj.matrix))
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


#### end of file -- bootstrap.likelihood.fit.constrained.R
