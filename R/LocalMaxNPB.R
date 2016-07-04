# function that, for node n (node), k (numparents) and score matrix (data)
# returns the highest k scored parents. 
#
# If more than k are equally ranked, only the first k are returned -- bias towards
# the positioning in the adjacency matrix.
parent_select_max = function(node, numparents, data, debug = FALSE)
{
	# no parents required
	if(numparents == 0) return(c());
	 
	# sort parent set, selecte first numparents
	candidates = sort(data[node, ], decreasing = TRUE)
	selected = candidates[1: numparents]
	
	if(debug) {		
		print('** Ordered parent set ')
		print(candidates)
		print('** Selected parents')
		print(names(selected))
	}

	# you should not have that k parents are required, any has 0-score
	if(any(selected == 0)) stop('parent_select_max: You ask for ', numparents, ' the candidates are ', candidates)

	return(names(selected))
}

# Algorithm - Reconstruction by locally maximum bootstrap consensus
algorithm_select_max = function(bootstrapped_models, num_parents_set, debug = FALSE)
{
		result = NULL
			
		# sump up scores -- npb scores
		scores = Reduce('+', bootstrapped_models)
		
		if(debug) {
				print('** Non-parametric Bootstrap Scores ');
				print(scores);
		}
		
		# resulting matrix
		model = matrix(0, nrow = nrow(scores), ncol = ncol(scores)) 
		rownames(model) = rownames(scores)
		colnames(model) = colnames(scores)
		
		parentset = rep(0, ncol(model))

		# for each node
		for(i in 1:nrow(model))
		{
			if(debug) print(paste('** Node ', i, ', k_i = ', num_parents_set[i]))
			
			# a local procedure estimates the best parents
			parents = parent_select_max(i, num_parents_set[i], scores, debug)
			model[i, parents] = 1					

			parentset[i] = length(parents)
		}	
		
		result$model = model
		result$scores = scores
		result$parentset = parentset

		return(result)
}

## Test just for this
# source('main.R')

# library(pheatmap)
# library(gridExtra)

# ## Note: adjacency matrix are transposed
# npb_local_max = algorithm_select_max(t(bootstrap_results_discrete), discrete.cardinalities, TRUE)

# grid.arrange(
# 	pheatmap(npb_local_max$model, 
# 		main = 'Locally maximum bootstrap consensus',
# 		color = c("lightgray", "brown3"),
# 		legend = F,
# 		cluster_cols = F, cluster_rows = F, silent = T)$gtable,
# 	pheatmap(npb_local_max$scores, 
# 		main = 'Non-Parametric Bootstrap Scores',
# 		color = colorRampPalette(c("white", "red"))(max(npb_local_max$scores)),
# 		display_numbers = T,
# 		cluster_cols = F, cluster_rows = F, silent = T)$gtable,
# 	ncol = 2
# )

# dev.new()
# library(igraph)
# G = graph.adjacency(npb_local_max$model)
# plot(G)
# title('Locally maximum bootstrap consensus model')
