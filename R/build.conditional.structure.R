# perform a robust estimation of the likelihood fit by conditional non-parametric bootstrap
build.conditional.consensus <- function( bootstrap_results, parents_cardinalities ) {
    
    # estimate the conditional adjacency matrix
    adj.matrix = algorithm.select.cmax(bootstrap_results,parents_cardinalities)$model
    
    # remove any loop
    adj.matrix = remove.loops(adj.matrix, to.list.of.edges(adj.matrix),to.list.of.scores(adj.matrix,bootstrap_results))
    
    return(adj.matrix)
    
}

# function that, for node n (node), k (numparents) and score matrix (data)
# returns the highest k scored parents. 
#
# If more than k are equally ranked, only the first k are returned -- bias towards
# the positioning in the adjacency matrix.
parent_select_cmax = function(node, numparents, data, debug = FALSE) {
    
    # no parents required, stop recursion
    if(numparents == 0) return(c());
          
    # sum up scores -- npb scores
    scores = Reduce('+', data)
     
    # sort parent set
    candidates = sort(scores[node, ], decreasing = TRUE)
    
    # this best
    parent = names(candidates[1, drop = FALSE])

    if(debug) {        
        print(paste('** Ordered parent set out of ', length(data), 'matrices'))
        print(candidates)
        print('** Best parent')
        print(parent)
    }
    
    # recursive call will pick the next maximum, among
    # those with edge node -> parent
    filter.fun = function(M){ return(M[node, parent] == 1)}
    
    # data subsetting 
    data = Filter(filter.fun, data)
    
    # data are also reshaped, we drop the column/row for parent
    prune.fun = function(M) {
        
        M = M[ , !(colnames(M) %in% parent) ]
        M = M[ !(colnames(M) %in% parent),  ]
        
        return(M)
        
    }
    data = lapply(data, prune.fun)

    if(length(data) == 0) {
        
        if(debug) print('Stopping search -- max parent set reached.')
        return(parent) 
        
    }
    
    # Recurisive call
    parent = c(parent,parent_select_cmax(node, numparents - 1, data, debug))
    
    return(parent)
    
}

# Algorithm - Reconstruction by locally maximum bootstrap consensus
algorithm.select.cmax = function( bootstrapped_models, num_parents_set, debug = FALSE ) {

        result = NULL
        
        # resulting matrix
        model = matrix(0, nrow = nrow(bootstrapped_models[[1]]), ncol = ncol(bootstrapped_models[[1]])) 
        rownames(model) = rownames(bootstrapped_models[[1]])
        colnames(model) = colnames(bootstrapped_models[[1]])
        
        parentset = rep(0, ncol(model))
        
        # for each node
        for(i in 1:nrow(model)) {
            if(debug) print(paste('** Node ', i, ', k_i = ', num_parents_set[i]))
            
            # a local procedure estimates the best parents
            parents = parent_select_cmax(i, num_parents_set[i], bootstrapped_models, debug)
            
            if(debug) print(paste('** Node ', i, ', true k = ', length(parents)))
            parentset[i] = length(parents)
            
            model[i, parents] = 1                    
        }    
        
        result$model = model
        result$parentset = parentset
        
        return(result)
        
}

to.list.of.edges = function( M ) {
    
    res = list()
    k = 1
    
    for(i in 1:nrow(M)) {
        for(j in 1:ncol(M)) {
            if(M[i,j] == 1) {
                    res[[k]] = list(parent=i,child=j)
                    k = k+1
            }
        }
    }
                
    return(res)
    
}

to.list.of.scores = function( M, data ) {
    
    # sump up scores -- npb scores
    scores = Reduce('+', data)

    res = c()
    
    for(i in 1:nrow(M)) {
        for(j in 1:ncol(M)) {
            if(M[i,j] == 1) {
                    res = c(res, scores[i,j])
            }
        }
    }
                
    return(res)
    
}
