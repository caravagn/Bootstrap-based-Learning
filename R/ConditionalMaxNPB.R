# function that, for node n (node), k (numparents) and score matrix (data)
# returns the highest k scored parents. 
#
# If more than k are equally ranked, only the first k are returned -- bias towards
# the positioning in the adjacency matrix.
parent_select_cmax = function(node, numparents, data, debug = FALSE)
{
	# no parents required, stop recursion
	if(numparents == 0) return(c());
	 	 
	# sump up scores -- npb scores
	scores = Reduce('+', data)
 	
 	# print(scores) 
 	 
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
	prune.fun = function(M){
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
	parent = c(parent,
		parent_select_cmax(node, numparents - 1, data, debug)
		)
	
	# # you should not have that k parents are required, any has 0-score
	# if(any(selected == 0)) stop('parent_select_max: You ask for ', numparents, ' the candidates are ', candidates)

	return(parent)
}

# Algorithm - Reconstruction by locally maximum bootstrap consensus
algorithm_select_cmax = function(bootstrapped_models, num_parents_set, debug = FALSE)
{
		result = NULL
		
		# resulting matrix
		model = matrix(0, nrow = nrow(bootstrapped_models[[1]]), ncol = ncol(bootstrapped_models[[1]])) 
		rownames(model) = rownames(bootstrapped_models[[1]])
		colnames(model) = colnames(bootstrapped_models[[1]])
		
		parentset = rep(0, ncol(model))
		
		# for each node
		for(i in 1:nrow(model))
		{
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

to.list.of.edges = function(M)
{
	res = list()
	k = 1
	
	for(i in 1:nrow(M))
		for(j in 1:ncol(M))
			if(M[i,j] == 1)
				{
					res[[k]] = list(parent=i,child=j)
					k = k+1
				}
				
	return(res)
}

to.list.of.scores = function(M,data)
{
	# sump up scores -- npb scores
	scores = Reduce('+', data)

	res = c()
	
	for(i in 1:nrow(M))
		for(j in 1:ncol(M))
			if(M[i,j] == 1)
					res = c(res, scores[i,j])
				
	return(res)
}

## Da qui, va nel main.R
source('main.test.R')
library(pheatmap)
library(gridExtra)
source('LocalMaxNPB.R')
library(bnlearn)

## Note: adjacency matrices are transposed
npb_local_max = algorithm_select_max(t(bootstrap_results_discrete), discrete.cardinalities, TRUE)

npb_local_max$model = t(npb_local_max$model)
npb_local_max$scores = t(npb_local_max$scores)

npb_local_max$model = remove.loops(
	npb_local_max$model,
	to.list.of.edges(npb_local_max$model),
	to.list.of.scores(npb_local_max$model, bootstrap_results_discrete)
	)

## Note: adjacency matrices are transposed
npb_local_cmax = algorithm_select_cmax(t(bootstrap_results_discrete), discrete.cardinalities, TRUE)

npb_local_cmax$model = t(npb_local_cmax$model)
npb_local_cmax$scores = Reduce('+', bootstrap_results_discrete)

npb_local_cmax$model = remove.loops(
	npb_local_cmax$model,
	to.list.of.edges(npb_local_cmax$model),
	to.list.of.scores(npb_local_cmax$model, bootstrap_results_discrete)
	)
	
parentset.matrix = function()
{
	trick = matrix(0, ncol = 3, nrow = nrow(npb_local_max$model))
	rownames(trick) = rownames(npb_local_max$model)
	colnames(trick) = c( 'k', 'LM', 'CLM')

	trick[,1] = discrete.cardinalities
	trick[,2] = discrete.cardinalities
	trick[,3] = npb_local_cmax$parentset
	
	trick = rbind(trick, rep(0,ncol(trick)))
	
	return(trick)
}
	
grid.arrange(
	pheatmap(npb_local_max $model, 
		main = 'LM bootstrap consensus',
		color = c("lightgray", "brown3"),
		legend = F,
		cluster_cols = F, cluster_rows = F, silent = T)$gtable,
	pheatmap(npb_local_cmax$model, 
		main = 'CLM bootstrap consensus',
		color = c("lightgray", "brown3"),
		legend = F,
		cluster_cols = F, cluster_rows = F, silent = T)$gtable,
	pheatmap(npb_local_max$scores, 
		main = 'Non-Parametric Bootstrap',
		color = colorRampPalette(c("white", "red"))(max(npb_local_max$scores)),
		display_numbers = T,
		number_format = "%d",
		cluster_cols = F, cluster_rows = F, silent = T)$gtable,
	pheatmap(parentset.matrix(), 
		main = 'Parents',
		color = colorRampPalette(c("white", "blue"))(max(trick)),
		number_format = "%d",
		display_numbers = T,
		cellwidth = 10,
		cluster_cols = F, cluster_rows = F, silent = T)$gtable,		
	ncol = 2
)

library(igraph)

G_local = graph.adjacency(npb_local_max$model)
G_conditional = graph.adjacency(npb_local_cmax$model)
G_results_discrete = graph.adjacency(results_discrete)

dev.new()
plot(G_local)
title(paste('LM bootstrap consensus model    n=', sum(npb_local_max$model)))
dev.new()
plot(G_conditional)
title(paste('CLM bootstrap consensus model    n=', sum(npb_local_cmax$model)))
dev.new()
plot(G_results_discrete)
title(paste('AIC model    n=', sum(results_discrete)))


