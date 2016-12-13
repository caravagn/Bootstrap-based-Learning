# perform the boostrap based maximum likelihood structure inference
perform.bootstrap.inference <- function(dataset, regularization = "bic", nboot.first = 100, nboot.second = 100, test.pvalue = 0.01, agony_files = paste0(getwd(), "/agony_files"), cores.ratio = 0.9) {

    TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")

    # set the number of cores to be used in the parallelization of the bootstrap
    cores = as.integer(cores.ratio * (detectCores() - 1))
    if (cores < 1) {
        cores = 1
    }

    
    # setup the parallelization to perform the bootstrap
    cat(paste("[*] Registering to use", cores, "/", detectCores(), "cores via \"parallel\" ..."))
    cl = makeCluster(cores)
    registerDoParallel(cl)
    cat(paste(" OK.\n"))

    # structure where to save the results
    results = list()

    # perform the first pass of bootstrap
    time = as.POSIXct(Sys.time(), format = "%H:%M:%S")
    cat(paste("[*] Bootstrapping Pi_boot: generating", nboot.first, "estimates ..."))
    bootstrap.first.pass = bootstrap.estimation.first.pass(dataset = dataset, regularization = regularization, nboot = nboot.first)
    results[["bootstrap.first.pass"]] = bootstrap.first.pass
    results[["bootstrap.first.pass.scores"]] = Reduce("+",bootstrap.first.pass)
    cat(paste(" OK [", round(as.POSIXct(Sys.time(), format = "%H:%M:%S") - time, 3), " sec].\n", sep = ""))

    # estimate the poset by bootstrap confidence
    time = as.POSIXct(Sys.time(), format = "%H:%M:%S")
    cat(paste("[*] Estimating the confidence poset (Pi_conf) ..."))
    confidence.poset = estimate.confidence.poset(bootstrap_results = bootstrap.first.pass, dataset = dataset)
    results[["confidence.poset"]] = confidence.poset
    cat(paste(" OK [", round(as.POSIXct(Sys.time(), format = "%H:%M:%S") - time, 3), " sec].\n", sep = ""))

    # estimate the poset by agony plus bootstrap confidence
    time = as.POSIXct(Sys.time(), format = "%H:%M:%S")
    unlink(agony_files, recursive = TRUE, force = TRUE)
    dir.create(agony_files, showWarnings = FALSE)
    cat(paste("[*] Estimating the agony poset (Pi_ago) ..."))
    get.agony.edges.list(bootstrap.first.pass, paste0(agony_files, "/inputs.txt"))
    system(paste0("./agony ", agony_files, "/inputs.txt ", agony_files, "/outputs.txt"), ignore.stdout = TRUE)
    agony.poset = build.agony.poset(paste0(agony_files, "/outputs.txt"), dataset, bootstrap.first.pass)
    results[["agony.poset"]] = agony.poset
    cat(paste(" OK [", round(as.POSIXct(Sys.time(), format = "%H:%M:%S") - time, 3), " sec].\n", sep = ""))

    # perform the second pass of bootstrap on both posets
    time = as.POSIXct(Sys.time(), format = "%H:%M:%S")
    cat("[*] Bootstrapping for Multiple Hypotheses Testing with Pi_conf: generating", nboot.second, "estimates ...")
    bootstrap.second.pass = list()
    bootstrap.second.pass[["confidence.poset"]] = bootstrap.estimation.second.pass(dataset = dataset, regularization = regularization, 
        nboot = nboot.second, contraint.poset = confidence.poset)
    cat(paste(" OK [", round(as.POSIXct(Sys.time(), format = "%H:%M:%S") - time, 3), " sec].\n", sep = ""))

    time = as.POSIXct(Sys.time(), format = "%H:%M:%S")
    cat("[*] Bootstrapping for Multiple Hypotheses Testing with Pi_ago: generating", nboot.second, "estimates ...")
    bootstrap.second.pass[["agony.poset"]] = bootstrap.estimation.second.pass(dataset = dataset, regularization = regularization, 
        nboot = nboot.second, contraint.poset = agony.poset)
    results[["bootstrap.second.pass"]] = bootstrap.second.pass
    bootstrap.second.pass.scores = list()
    bootstrap.second.pass.scores[["confidence.poset"]][["model"]] = Reduce("+",results$bootstrap.second.pass$confidence.poset$model)
    bootstrap.second.pass.scores[["confidence.poset"]][["null"]] = Reduce("+",results$bootstrap.second.pass$confidence.poset$null)
    bootstrap.second.pass.scores[["agony.poset"]][["model"]] = Reduce("+",results$bootstrap.second.pass$agony.poset$model)
    bootstrap.second.pass.scores[["agony.poset"]][["null"]] = Reduce("+",results$bootstrap.second.pass$agony.poset$null)
    results[["bootstrap.second.pass.scores"]] = bootstrap.second.pass.scores
    cat(paste(" OK [", round(as.POSIXct(Sys.time(), format = "%H:%M:%S") - time, 3), " sec].\n", sep = ""))

    # estimate the final Bayesian Network with the estimated posets
    time = as.POSIXct(Sys.time(), format = "%H:%M:%S")
    cat("[*] Estimating the final Bayesian Network with Pi_conf ...")
    confidence.inference = perform.bn.inference(bootstrap.second.pass[["confidence.poset"]],confidence.poset,test.pvalue)
    results[["confidence.inference"]] = confidence.inference[["inference"]]
    results[["confidence.inference.pvalues"]] = confidence.inference[["pvalues"]]
    cat(paste(" OK [", round(as.POSIXct(Sys.time(), format = "%H:%M:%S") - time, 3), " sec].\n", sep = ""))

    time = as.POSIXct(Sys.time(), format = "%H:%M:%S")
    cat("[*] Estimating the final Bayesian Network with Pi_ago ...")
    agony.inference = perform.bn.inference(bootstrap.second.pass[["agony.poset"]],agony.poset,test.pvalue)
    results[["agony.inference"]] = agony.inference[["inference"]]
    results[["agony.inference.pvalues"]] = agony.inference[["pvalues"]]
    cat(paste(" OK [", round(as.POSIXct(Sys.time(), format = "%H:%M:%S") - time, 3), " sec].\n", sep = ""))
    
    # Matrix compression (last step)
    results$bootstrap.first.pass = lapply(results$bootstrap.first.pass, compress.matrix)
    results$bootstrap.second.pass$agony.poset$model = lapply(results$bootstrap.second.pass$agony.poset$model, compress.matrix)
    results$bootstrap.second.pass$agony.poset$null = lapply(results$bootstrap.second.pass$agony.poset$null, compress.matrix)
    results$bootstrap.second.pass$confidence.poset$model =lapply(results$bootstrap.second.pass$confidence.poset$model, compress.matrix)
    results$bootstrap.second.pass$confidence.poset$null =lapply(results$bootstrap.second.pass$confidence.poset$null, compress.matrix)

    stopCluster(cl)
    
    cat(paste("Total time ", round(as.POSIXct(Sys.time(), format = "%H:%M:%S") - TIME, 3), " sec.\n", sep = ""))

    return(results)

}

# perform a robust estimation of the likelihood fit by non-parametric bootstrap
bootstrap.estimation.first.pass <- function(dataset, regularization, command = "hc", nboot = 100, random.seed = NULL, verbose = FALSE) {

    # structure to save the results of the bootstrap
    bootstrap.results = NULL

    # set the seed to be used for the random samplings
    set.seed(random.seed)

    # perform nboot bootstrap resampling
    r = foreach(num = 1:nboot, .packages = "bnlearn", .export = c("perform.likelihood.fit", "compress.matrix")) %dopar% {

        # create the sampled dataset for the current iteration
        samples = sample(1:nrow(dataset), size = nrow(dataset), replace = TRUE)
        resampled.dataset = dataset[samples, ]

        # perform the likelihood fit on the bootstrapped dataset
        reconstructed.adj.matrix = perform.likelihood.fit(resampled.dataset, regularization, command)
    }

    bootstrap.results = r

    return(bootstrap.results)

}

# compress a matrix to bnlearn representation
compress.matrix <- function(M) {

    string = ""
    
    string = paste(
        sapply(colnames(M),
            function(x){ 
                str = paste("[", x, sep = '')
                if(any(M[,x]==1)) 
                    str = paste(str, '|', paste(rownames(M)[which(M[,x]==1)], collapse = ':'), sep ='')    
                str = paste(str, "]", sep = '')
            }
        ),
        collapse = '')
    
    return(string) 

}

# perform structural inference by maximum likelihood
perform.likelihood.fit <- function(dataset, regularization, command = "hc") {

    # create an empty network
    empty_net = empty.graph(nodes=colnames(dataset))
    curr_random_set = sample(colnames(dataset),size=2)
    curr_parent = curr_random_set[1]
    curr_child = curr_random_set[2]
    curr_start = set.arc(empty_net,from=curr_parent,to=curr_child)

    # perform maximum likelihood estimation either by hill climbing or tabu search
    if (command == "hc") 
        my.net = hc(dataset,start=curr_start,score=regularization)
    if (command == "tabu") 
        my.net = tabu(dataset,start=curr_start,score=regularization)

    return(amat(my.net))

}

# perform a robust estimation of the likelihood fit by non-parametric bootstrap
estimate.confidence.poset <- function(bootstrap_results, dataset) {

    # estimate the poset based on bootstrap confidence
    estimated.poset.confidence = build.confidence.poset(bootstrap_results, rep(nrow(bootstrap_results[[1]]), nrow(bootstrap_results[[1]])))
    colnames(estimated.poset.confidence) = colnames(dataset)
    rownames(estimated.poset.confidence) = colnames(dataset)

    return(estimated.poset.confidence)

}

# perform a robust estimation of the poset by likelihood fit and non-parametric bootstrap
build.confidence.poset <- function(bootstrap_results, parents_cardinalities) {

    # estimate the confidence of each edge
    bootstrap_results = Reduce("+", bootstrap_results)

    # create the consensus graph
    adj.matrix = array(0, c(nrow(bootstrap_results), ncol(bootstrap_results)))
    added_edges = list()
    added_scores = NULL
    for (i in 1:ncol(bootstrap_results)) {
        if (parents_cardinalities[i] > 0) {
            # consider the best parents_cardinalities arcs greater than 0
            sorted_bootstrap = sort(bootstrap_results[, i], decreasing = TRUE)
            valid_parents = sorted_bootstrap[which(sorted_bootstrap[1:parents_cardinalities[i]] > 0)]
            if (length(valid_parents) > 0) {
                added_scores = c(added_scores, as.vector(valid_parents))
                added_valid_parents = NULL
                for (j in 1:length(valid_parents)) {
                    new_valid_parent_pos = which(colnames(bootstrap_results) == names(valid_parents)[j])
                    added_edges[[(length(added_edges) + 1)]] = list(parent = new_valid_parent_pos, child = i)
                    added_valid_parents = c(added_valid_parents, new_valid_parent_pos)
                }
                adj.matrix[added_valid_parents, i] = 1
            }
        }
    }

    # remove any loop from the adjacency matrix
    adj.matrix = remove.loops(adj.matrix, added_edges, added_scores)

    return(adj.matrix)

}

# remove any loop in the given adjacency matrix
remove.loops <- function(adj.matrix, edges, scores) {

    # if I have at least one edge
    if (length(edges) > 0) {

        # compute the ordered confidence scores
        duplicated_scores = unique(scores[duplicated(scores)])
        ordered.scores = sort(scores, decreasing = FALSE, index.return = TRUE)
        ordered.scores_values = ordered.scores$x
        ordered.scores_idx = ordered.scores$ix
        ordered.scores_idx_random = ordered.scores_idx

        # randomly order arcs with equal scores to avoid any bias
        for (my_dup_score in duplicated_scores) {
            curr_dup_scores = which(ordered.scores_values == my_dup_score)
            curr_dup_scores_random = sample(curr_dup_scores, size = length(curr_dup_scores), replace = FALSE)
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
            curr.graph = graph.adjacency(adj.matrix, mode = "directed")
            is.path = suppressWarnings(get.shortest.paths(curr.graph, curr.edge.j, curr.edge.i)$vpath)
            is.path = length(unlist(is.path))

            # if there is a path between the two nodes, remove edge i --> j
            if (is.path > 0) {
                adj.matrix[curr.edge.i, curr.edge.j] = 0
            }

        }

    }

    return(adj.matrix)

}

# set the files for agony computation
get.agony.edges.list <- function(bootstrap_results, agony_file) {
    for (i in 1:length(bootstrap_results)) {
        curr_adj.matrix = bootstrap_results[[i]]
        for (j in 1:nrow(curr_adj.matrix)) {
            for (k in 1:ncol(curr_adj.matrix)) {
                if (curr_adj.matrix[j, k] == 1) {
                    cat(as.character(j), as.character(k), "\n", file = agony_file, append = TRUE)
                }
            }
        }
    }
}

# build the agony based poset
build.agony.poset <- function(agony_file, dataset, bootstrap_results) {

    # read the results from agony computation
    my_ordering = read.table(file = agony_file, header = FALSE, sep = " ", check.names = FALSE, stringsAsFactors = FALSE)

    # estimate a poset for the DAG using agony. 
    # Note: this can still contain cycles, in fact elements at the same ordering, have arcs both ways
    agony_poset = get.agony.poset(my_ordering, ncol(dataset))

    # make the poset acyclic using bootstrap confidence
    agony_poset = perform.agony.estimate(Reduce("+", bootstrap_results), agony_poset)
    colnames(agony_poset) = colnames(dataset)
    rownames(agony_poset) = colnames(dataset)

    # return the poset
    return(agony_poset)

}

# construct the poset by agony
get.agony.poset <- function(ordering, node_size) {
    poset = array(0, c(node_size, node_size))
    for (i in min(ordering[, 2]):max(ordering[, 2])) {
        curr_parents = ordering[which(ordering[, 2] == i), 1]
        curr_childred = ordering[which(ordering[, 2] > i), 1]
        # add arcs from lower to higher ranked nodes
        if (length(curr_parents) > 0 && length(curr_childred) > 0) {
            for (j in curr_parents) {
                poset[j, curr_childred] = 1
            }
        }
        # add arcs in both ways for equal ranked nodes (NOTE: with this, the poset can be cyclic!)
        if (length(curr_parents) > 1) {
            for (a in 1:length(curr_parents)) {
                for (b in a:length(curr_parents)) {
                    if (a != b) {
                        poset[curr_parents[a], curr_parents[b]] = 1
                        poset[curr_parents[b], curr_parents[a]] = 1
                    }
                }
            }
        }
    }
    return(poset)
}

# perform a robust estimation of the poset by likelihood fit and non-parametric bootstrap
perform.agony.estimate <- function(bootstrap_results, poset) {

    # constrain the confidence by bootstrap based on the poset
    constrained_bootstrap = (bootstrap_results * poset)

    # estimate the poset based on bootstrap confidence
    estimated.poset.agony = build.consensus.poset_agony(constrained_bootstrap, rep(nrow(constrained_bootstrap), nrow(constrained_bootstrap)))

    return(estimated.poset.agony)

}

# perform a robust estimation of the poset by likelihood fit by non-parametric bootstrap
build.consensus.poset_agony <- function(bootstrap_results, parents_cardinalities) {

    # create the consensus graph
    adj.matrix = array(0, c(nrow(bootstrap_results), ncol(bootstrap_results)))
    added_edges = list()
    added_scores = NULL
    for (i in 1:ncol(bootstrap_results)) {
        if (parents_cardinalities[i] > 0) {
            # consider the best parents_cardinalities arcs greater than 0
            sorted_bootstrap = sort(bootstrap_results[, i], decreasing = TRUE)
            valid_parents = sorted_bootstrap[which(sorted_bootstrap[1:parents_cardinalities[i]] > 0)]
            if (length(valid_parents) > 0) {
                added_scores = c(added_scores, as.vector(valid_parents))
                added_valid_parents = NULL
                for (j in 1:length(valid_parents)) {
                    new_valid_parent_pos = which(colnames(bootstrap_results) == names(valid_parents)[j])
                    added_edges[[(length(added_edges) + 1)]] = list(parent = new_valid_parent_pos, child = i)
                    added_valid_parents = c(added_valid_parents, new_valid_parent_pos)
                }
                adj.matrix[added_valid_parents, i] = 1
            }
        }
    }

    # remove any loop from the adjacency matrix
    adj.matrix = remove.loops(adj.matrix, added_edges, added_scores)

    return(adj.matrix)

}

# perform a robust estimation of the likelihood fit by non-parametric bootstrap for a given poset both for model and random
bootstrap.estimation.second.pass <- function(dataset, regularization, command = "hc", nboot = 100, random.seed = NULL, verbose = FALSE, contraint.poset) {

    # structure to save the results of the bootstrap
    bootstrap.results = NULL

    # set the seed to be used for the random samplings
    set.seed(random.seed)


    # perform nboot bootstrap resampling
    r = foreach(num = 1:nboot, .packages = "bnlearn", .export = "perform.constrained.likelihood.fit") %dopar% {

        # create the sampled datasets for the current iteration
        samples = sample(1:nrow(dataset), size = nrow(dataset), replace = TRUE)
        resampled.dataset = dataset[samples, ]
        resampled.dataset_null = resampled.dataset
        for (null_model in 1:ncol(resampled.dataset_null)) {
            resampled.dataset_null[, null_model] = resampled.dataset_null[sample(1:nrow(resampled.dataset_null), size = nrow(resampled.dataset_null), 
                replace = FALSE), null_model]
        }

        # perform the likelihood fit on the bootstrapped dataset
        reconstructed.adj.matrix = list()
        reconstructed.adj.matrix[["model"]] = perform.constrained.likelihood.fit(resampled.dataset, contraint.poset, regularization, 
            command)
        reconstructed.adj.matrix[["null"]] = perform.constrained.likelihood.fit(resampled.dataset_null, contraint.poset, regularization, 
            command)
        reconstructed.adj.matrix

    }

    if (verbose) {
        cat("\tBootstrap completed.\n")
    }

    bootstrap.results = r

    bootstrap.results_model = list()
    bootstrap.results_null = list()
    for (i in 1:length(bootstrap.results)) {
        bootstrap.results_model[[paste0("iteration_", as.character(i))]] = bootstrap.results[[i]][["model"]]
        bootstrap.results_null[[paste0("iteration_", as.character(i))]] = bootstrap.results[[i]][["null"]]
    }
    bootstrap.results = list()
    bootstrap.results[["model"]] = bootstrap.results_model
    bootstrap.results[["null"]] = bootstrap.results_null

    return(bootstrap.results)

}

# perform structural inference by maximum likelihood
perform.constrained.likelihood.fit <- function(dataset, poset.adj.matrix, regularization, command = "hc") {

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
            if (i != j) {
                if (poset.adj.matrix[i, j] == 0) {
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

    # create an empty network
    empty_net = empty.graph(nodes=colnames(dataset))
    curr_random_set = sample(colnames(dataset),size=2)
    curr_parent = curr_random_set[1]
    curr_child = curr_random_set[2]
    curr_start = set.arc(empty_net,from=curr_parent,to=curr_child)

    # perform maximum likelihood estimation either by hill climbing or tabu search
    if (cont > 0) {
        blacklist = data.frame(from = parent, to = child)
        if (command == "hc") {
            my.net = hc(dataset,start=curr_start,score=regularization,blacklist=blacklist)
        } else if (command == "tabu") {
            my.net = tabu(dataset,start=curr_start,score=regularization,blacklist=blacklist)
        }
    } else {
        if (command == "hc") {
            my.net = hc(dataset,start=curr_start,score=regularization)
        } else if (command == "tabu") {
            my.net = tabu(dataset,start=curr_start,score=regularization)
        }
    }

    return(amat(my.net))

}

# perform the inference of the Bayesian network with the given bootstrap estimates
perform.bn.inference <- function( bootstrap_results, curr.poset, test_pvalue ) {
    pvalues_results = list()
    pvalues_estimates = get.pvalue.estimate(bootstrap_results,curr.poset)
    pvalues_results[["pvalues"]] = pvalues_estimates
    for(vals in names(pvalues_estimates)) {
        curr_res = pvalues_estimates[[vals]]
        curr_res[which(curr_res<test_pvalue,arr.ind=TRUE)] = -1
        curr_res[which(curr_res>=test_pvalue,arr.ind=TRUE)] = 0
        curr_res = curr_res * -1
        colnames(curr_res) = colnames(curr.poset)
        rownames(curr_res) = rownames(curr.poset)
        pvalues_estimates[[vals]] = curr_res
    }
    pvalues_results[["inference"]] = pvalues_estimates
    return(pvalues_results)
}

# estimate the confidence of each arc by chi-squared test (i.e., success rate consistently lower in the null estimate by bootstrap)
get.pvalue.estimate <- function( bootstrap_results, curr.poset ) {

    bootstrap_res = get.bootstrap.second.pass.res(bootstrap_results)
    bootstrap_model = bootstrap_res[["model"]]
    bootstrap_null = bootstrap_res[["null"]]

    curr_confidence = array(0,c(nrow(bootstrap_model),ncol(bootstrap_model)))
    for(p in 1:nrow(bootstrap_model)) {
        for(q in 1:ncol(bootstrap_model)) {
            curr_model_success_estimate = unlist(bootstrap_model[p,q])
            curr_null_success_estimate = unlist(bootstrap_null[p,q])
            curr_success_rates = c(sum(curr_model_success_estimate),sum(curr_null_success_estimate))
            curr_trial_sizes = c(length(curr_model_success_estimate),length(curr_null_success_estimate))
            if(curr_success_rates[1]==0 && curr_success_rates[2]==0) {
                curr_confidence[p,q] = 1
            }
            else {
                curr_confidence[p,q] = prop.test(curr_success_rates,curr_trial_sizes,alternative="greater")$p.value
            }
        }
    }
    curr_valid_arcs = which(curr.poset==1,arr.ind=TRUE)
    # perform pvalue correction by FDR
    curr_confidence_fdr_adj = p.adjust(as.vector(curr_confidence[curr_valid_arcs]),method="fdr")
    curr_confidence_fdr = curr_confidence
    curr_confidence_fdr[curr_valid_arcs] = curr_confidence_fdr_adj
    # perform pvalue correction by Holm
    curr_confidence_holm_adj = p.adjust(as.vector(curr_confidence[curr_valid_arcs]),method="holm")
    curr_confidence_holm = curr_confidence
    curr_confidence_holm[curr_valid_arcs] = curr_confidence_holm_adj

    return(list(pvalues=curr_confidence,qvalues.fdr=curr_confidence_fdr,qvalues.holm=curr_confidence_holm))

}

# pre-process the results of the second bootstrap
get.bootstrap.second.pass.res <- function( bootstrap_second_pass ) {
    
    boot_confidence = array(list(),c(nrow(bootstrap_second_pass[[1]][[1]]),ncol(bootstrap_second_pass[[1]][[1]])))
    boot_null_confidence = array(list(),c(nrow(bootstrap_second_pass[[1]][[1]]),ncol(bootstrap_second_pass[[1]][[1]])))
    for(m in 1:length(bootstrap_second_pass[["model"]])) {
        curr_boot_confidence = bootstrap_second_pass[["model"]][[m]]
        curr_boot_null_confidence = bootstrap_second_pass[["null"]][[m]]
        for(n in 1:nrow(curr_boot_confidence)) {
            for(o in 1:ncol(curr_boot_confidence)) {
                boot_confidence[n,o] = list(c(unlist(boot_confidence[n,o]),curr_boot_confidence[n,o]))
                boot_null_confidence[n,o] = list(c(unlist(boot_null_confidence[n,o]),curr_boot_null_confidence[n,o]))
            }
        }
    }
    res = list()
    res[["model"]] = boot_confidence
    res[["null"]] = boot_null_confidence
    
    bootstrap_confidence = res
    
    return(bootstrap_confidence)
}
